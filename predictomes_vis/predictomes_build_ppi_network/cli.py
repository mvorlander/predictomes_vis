#!/usr/bin/env python3
"""Query SPOC AlphaFold pairwise PPI predictions and build a seed-centric network.

The CSV exported from SPOC contains confident AlphaFold protein-protein interaction
predictions. This helper lets you look up rows involving a UniProt ID or name of
interest, expand those hits to include every confident partner, and emit outputs
including CSV, JSON, and a PyVis HTML network enriched with UniProt Function annotations.
"""

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import os
import select
import difflib
import urllib.error
import urllib.request
import urllib.parse
import ssl
import certifi
from collections import defaultdict
from pathlib import Path
from textwrap import dedent
from typing import Dict, Iterable, List, Sequence, Tuple

SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{identifier}.json"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search?size=1&format=json&query={query}"


def parse_args():
    parser = argparse.ArgumentParser(
        prog="predictomes-vis network",
        description=dedent(
            """\
            Predictomes network visualization (predictomes.org data).
            Data sources:
            - FDR (default): High-confidence SPOC results (<10%% FDR).
            - peak_score: Full predictome matrix (~1.6M pairs parsed by M. Vorländer, Plaschka lab, IMP Vienna).
            """
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    src_group = parser.add_argument_group("Data source")
    input_group = parser.add_argument_group("Input POIs")
    filter_group = parser.add_argument_group("Filters / expansion")
    viz_group = parser.add_argument_group("Visualization fields")
    output_group = parser.add_argument_group("Output")

    default_csv = Path(__file__).resolve().parents[1] / "predictome_data" / "20251113_download_hsbps.csv"
    default_matrix = None  # resolved later to allow multiple fallbacks
    src_group.add_argument(
        "--metric",
        choices=["FDR", "peak_score"],
        default="FDR",
        help="Metric to drive network building: 'FDR' uses high-confidence SPOC results (<10%% FDR); 'peak_score' uses the full predictome matrix thresholded by peak score.",
    )
    src_group.add_argument(
        "--predictome-matrix",
        type=Path,
        default=default_matrix,
        help="Path to AlphaFold predictome matrix TSV (used when --metric peak_score; default: bundled lookup).",
    )
    src_group.add_argument(
        "-f",
        "--file",
        dest="csv_path",
        default=str(default_csv),
        help="Path to the SPOC CSV export (used when --metric FDR).",
    )
    src_group.add_argument(
        "--peak-threshold",
        "-k",
        dest="predictome_peak_threshold",
        type=float,
        default=0.75,
        help="Minimum normalized peak average (1-score/ceiling) to include an edge when --metric peak_score.",
    )
    src_group.add_argument(
        "--peak-ceiling",
        "-K",
        dest="predictome_peak_ceiling",
        type=float,
        default=30.0,
        help="Peak score upper bound used for normalization when --metric peak_score.",
    )
    src_group.add_argument(
        "--iptm-threshold",
        "-I",
        dest="predictome_iptm_threshold",
        type=float,
        default=0.5,
        help="Minimum ipTM max to include an edge when --metric peak_score.",
    )
    output_group.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("outputs/ppi_networks"),
        help="Base directory for generated files when paths are derived automatically (default: %(default)s).",
    )
    input_group.add_argument(
        "--poi",
        action="append",
        help="UniProt ID or UniProt entry name to match (repeatable; comma or space separated allowed).",
    )
    input_group.add_argument(
        "--poi-file",
        help="Text file with POIs (one per line or comma/space separated).",
    )
    viz_group.add_argument(
        "-l",
        "--edge-label-field",
        default="ipTM_max",
        help="CSV column to use for edge labels (default: %(default)s).",
    )
# Output toggles are hidden; defaults apply.
    viz_group.add_argument(
        "--edge-width-field",
        default="ipTM_max",
        help="CSV column controlling PyVis edge thickness (default: %(default)s).",
    )
    viz_group.add_argument(
        "--edge-color-field",
        default="SPOC_score",
        help="CSV column controlling PyVis edge color gradient (default: %(default)s).",
    )
    filter_group.add_argument(
        "--no-function-fetch",
        action="store_true",
        help="Skip retrieving UniProt Function annotations for nodes.",
    )
    filter_group.add_argument(
        "--large-network-threshold",
        type=int,
        default=100,
        help="Warn when node count exceeds this value before fetching UniProt data (default: %(default)s).",
    )
    filter_group.add_argument(
        "--expansion-depth",
        type=int,
        default=1,
        help="Extra hops beyond direct neighbors (0=only direct neighbors).",
    )
    filter_group.add_argument(
        "--query-file",
        help="Optional text file with one query per line as 'UniProt_ID[,UniProt_Name]'.",
    )
    filter_group.add_argument(
        "--per-query-output",
        action="store_true",
        help="Generate separate outputs for every provided POI (CSV written to <stem>.csv).",
    )
    viz_group.add_argument(
        "-C",
        "--node-color",
        dest="node_color_field",
        help="CSV column whose numeric average per node controls PyVis node color.",
    )
    viz_group.add_argument(
        "-S",
        "--node-size",
        dest="node_size_field",
        help="CSV column whose numeric average per node controls PyVis node size.",
    )
    parser.set_defaults(
        no_json=False,
        no_pyvis=False,
        json_out=None,
        pyvis_out=None,
        uniprot_id=[],
        uniprot_name=[],
    )
    args = parser.parse_args()
    # Adjust defaults for predictome mode to avoid missing fields
    if args.metric == "peak_score":
        # Sensible defaults for predictome: color by ipTM_max, width by peak score
        if args.edge_color_field == "SPOC_score":
            args.edge_color_field = "iptm_max"
        if args.edge_label_field == "SPOC_score":
            args.edge_label_field = "peak_avg_norm"
        if args.edge_width_field in {"ipTM_max", "SPOC_score"}:
            args.edge_width_field = "peak_avg_norm"

    def normalize_queries(raw):
        queries = []
        for item in raw or []:
            for part in re.split(r"[,\s]+", item.strip()):
                if part:
                    queries.append(part)
        return queries

    poi_file_tokens: list[str] = []
    if args.poi_file:
        poi_file_tokens = []
        qf_path = Path(args.poi_file)
        if not qf_path.exists():
            sys.exit(f"POI file not found: {args.poi_file}")
        for line in qf_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            poi_file_tokens.extend([tok for tok in re.split(r"[,\s]+", line) if tok])

    args.uniprot_id = normalize_queries(args.uniprot_id) + normalize_queries(args.poi) + poi_file_tokens
    args.uniprot_name = normalize_queries(args.uniprot_name) + normalize_queries(args.poi) + poi_file_tokens

    if not args.uniprot_id and not args.uniprot_name:
        parser.error("Provide --poi or --poi-file.")
    return args


def split_tokens(value: str) -> List[str]:
    if not value:
        return []
    return [token.strip() for token in re.split(r"\s*/\s*", value) if token.strip()]


def matches_any(value: str, queries: List[str]) -> bool:
    if not queries:
        return True
    tokens = split_tokens(value)
    lowered = {q.lower() for q in queries}
    return any(token.lower() in lowered for token in tokens)


def suggest_matches(queries: List[str], rows: List[dict]) -> List[str]:
    search_space: set[str] = set()
    for row in rows:
        search_space.update(split_tokens(row.get("UniProt_ID", "")))
        search_space.update(split_tokens(row.get("UniProt_Name", "")))
    suggestions: set[str] = set()
    for q in queries:
        suggestions.update(difflib.get_close_matches(q, list(search_space), n=5, cutoff=0.6))
    return sorted(suggestions, key=str.lower)


def load_rows(csv_path: Path) -> Tuple[Sequence[str], List[dict]]:
    print(f"[info] Loading FDR CSV: {csv_path}", file=sys.stderr)
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        cleaned_fields = [name for name in reader.fieldnames or [] if name]
        rows = []
        for row in reader:
            rows.append({field: row.get(field, "") for field in cleaned_fields})
    print(f"[info] Loaded {len(rows):,} rows from FDR CSV.", file=sys.stderr)
    return cleaned_fields, rows


def float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_peak_matrix(tsv_path: Path, peak_threshold: float, iptm_threshold: float, peak_ceiling: float) -> Tuple[Sequence[str], List[dict]]:
    print(f"[info] Loading peak matrix: {tsv_path}", file=sys.stderr)
    with tsv_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cleaned_fields = [name for name in reader.fieldnames or [] if name]
        required = {"protein_a", "protein_b", "iptm_max", "peak_score_rank_1", "peak_score_rank_2", "peak_score_rank_3"}
        missing = required - set(cleaned_fields)
        if missing:
            raise SystemExit(f"Missing required columns in peak matrix: {', '.join(sorted(missing))}")

        rows: List[dict] = []
        total_rows = 0
        for row in reader:
            total_rows += 1
            peak_cols = [row.get(f"peak_score_rank_{i}") for i in (1, 2, 3)]
            peak_vals = [float_or_none(val) for val in peak_cols if float_or_none(val) is not None]
            if len(peak_vals) != 3:
                continue
            peaks_norm = [1 - (val / peak_ceiling) for val in peak_vals]
            peak_avg_norm = sum(peaks_norm) / len(peaks_norm)
            if peak_avg_norm < peak_threshold:
                continue

            prot_a = (row.get("protein_a") or "").strip()
            prot_b = (row.get("protein_b") or "").strip()
            if not prot_a or not prot_b:
                continue

            iptm_max = float_or_none(row.get("iptm_max", ""))
            if iptm_max is not None and iptm_max < iptm_threshold:
                continue
            iptm_avg = None
            iptm_cols = [float_or_none(row.get(f"iptm_rank_{i}")) for i in (1, 2, 3)]
            iptm_vals = [v for v in iptm_cols if v is not None]
            if iptm_vals:
                iptm_avg = sum(iptm_vals) / len(iptm_vals)

            rec = {
                "UniProt_ID": f"{prot_a} / {prot_b}",
                "UniProt_Name": f"{prot_a} / {prot_b}",
                "protein_a": prot_a,
                "protein_b": prot_b,
                "ipTM_max": iptm_max if iptm_max is not None else "",
                "iptm_max": iptm_max if iptm_max is not None else "",
                "iptm_avg": iptm_avg if iptm_avg is not None else "",
                "peak_avg_norm": peak_avg_norm,
                "peak_min_norm": min(peaks_norm),
                "peak_max_norm": max(peaks_norm),
                "source": "peak_matrix",
            }
            # keep original numeric columns too
            rec.update({k: row.get(k, "") for k in cleaned_fields})
            rows.append(rec)
    field_set = set(rows[0].keys()) if rows else set(cleaned_fields)
    print(f"[info] Peak matrix rows: {total_rows:,} read, {len(rows):,} passed filters (peak>={peak_threshold}, iptm>={iptm_threshold})", file=sys.stderr)
    return list(field_set), rows


def find_peak_matrix(user_path: Path | None) -> Path:
    if user_path:
        return user_path
    env_path = os.environ.get("PREDICTOMES_PEAK_MATRIX")
    if env_path:
        return Path(env_path)
    package_root = Path(__file__).resolve().parents[1]
    repo_root = Path(__file__).resolve().parents[3]
    candidates = [
        package_root / "predictome_data" / "json_metrics_matrix.tsv",
        repo_root / "data" / "predictome_json_stats_251115" / "json_metrics_matrix.tsv",
        repo_root / "alphafold_peaks" / "data" / "predictome_json_stats_251115" / "json_metrics_matrix.tsv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    # fall back to first candidate even if missing, main will error
    return candidates[0]


def extract_interaction(row: dict) -> Tuple[str, str, str, str] | None:
    ids = split_tokens(row.get("UniProt_ID", ""))
    names = split_tokens(row.get("UniProt_Name", ""))
    if len(ids) != 2:
        return None
    name_a = names[0] if len(names) > 0 else ""
    name_b = names[1] if len(names) > 1 else ""
    return ids[0], ids[1], name_a, name_b


def collect_seed_ids(rows: Iterable[dict]) -> set[str]:
    seeds: set[str] = set()
    for row in rows:
        interaction = extract_interaction(row)
        if not interaction:
            continue
        seeds.update(interaction[:2])
    return seeds


def gather_network_rows(all_rows: Iterable[dict], seeds: set[str]) -> List[dict]:
    network_rows: List[dict] = []
    for row in all_rows:
        ids = split_tokens(row.get("UniProt_ID", ""))
        names = split_tokens(row.get("UniProt_Name", ""))
        if any(protein in seeds for protein in ids) or any(name in seeds for name in names):
            network_rows.append(row)
    return network_rows


def expand_seed_network(rows: List[dict], initial_seeds: set[str], depth: int) -> tuple[List[dict], set[str]]:
    # Interpret depth as "extra hops beyond direct neighbors".
    # depth=0 -> only direct neighbors; depth=1 -> neighbors + neighbors-of-neighbors, etc.
    collected_rows: List[dict] = []
    remaining_rows = list(rows)

    seeds = set(initial_seeds)
    direct_rows = gather_network_rows(remaining_rows, seeds)
    collected_rows.extend(direct_rows)
    for row in direct_rows:
        seeds.update(split_tokens(row.get("UniProt_ID", "")))
    remaining_rows = [row for row in remaining_rows if row not in direct_rows]

    frontier = set(seeds)
    if depth > 0:
        for _ in range(depth):
            added = gather_network_rows(remaining_rows, frontier)
            if not added:
                break
            collected_rows.extend(added)
            frontier = set()
            for row in added:
                frontier.update(split_tokens(row.get("UniProt_ID", "")))
            seeds.update(frontier)
            remaining_rows = [row for row in remaining_rows if row not in added]

    return collected_rows, seeds


def build_graph(rows: Iterable[dict]):
    adjacency: Dict[str, List[Tuple[str, dict]]] = defaultdict(list)
    node_labels: Dict[str, set[str]] = defaultdict(set)
    edges: List[Tuple[str, str, dict]] = []

    for row in rows:
        interaction = extract_interaction(row)
        if not interaction:
            continue
        prot_a, prot_b, name_a, name_b = interaction
        adjacency[prot_a].append((prot_b, row))
        adjacency[prot_b].append((prot_a, row))
        if name_a:
            node_labels[prot_a].add(name_a)
        if name_b:
            node_labels[prot_b].add(name_b)
        edges.append((prot_a, prot_b, row))
    return adjacency, node_labels, edges


def preferred_label(labels: Iterable[str], fallback: str) -> str:
    for label in labels:
        if label:
            return label
    return fallback


def edge_label(row: dict, field: str) -> str:
    value = row.get(field, "")
    if value:
        return value
    if field != "SPOC_score":
        return row.get("SPOC_score", "")
    return ""


def render_network_summary(adjacency, node_labels, seeds, label_field: str):
    if not adjacency:
        print("No PPI network could be assembled from the matches.", file=sys.stderr)
        return

    print("\nPPI network (seed nodes marked with *):", file=sys.stderr)
    for node in sorted(adjacency.keys()):
        star = "*" if node in seeds else " "
        label = preferred_label(node_labels.get(node, []), node)
        partners = ", ".join(
            f"{partner} (label={edge_label(row, label_field) or 'NA'})"
            for partner, row in sorted(adjacency[node], key=lambda item: item[0])
        )
        print(f"{star} {label} [{node}] -> {partners}", file=sys.stderr)


def load_query_file(path: str | None) -> List[dict]:
    if not path:
        return []
    entries: List[dict] = []
    query_path = Path(path)
    if not query_path.exists():
        sys.exit(f"Query file not found: {path}")
    with query_path.open() as handle:
        for line in handle:
            raw = line.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = [part.strip() for part in raw.split(",", 1)]
            ids = [parts[0]] if parts and parts[0] else []
            names = [parts[1]] if len(parts) > 1 and parts[1] else []
            if not ids and not names:
                continue
            entries.append({"ids": ids, "names": names})
    return entries


def build_queries(args) -> List[dict]:
    file_queries = load_query_file(args.query_file)
    queries: List[dict] = []
    if args.per_query_output:
        queries.extend({"ids": [uid], "names": []} for uid in args.uniprot_id)
        queries.extend({"ids": [], "names": [uname]} for uname in args.uniprot_name)
        queries.extend(file_queries)
    else:
        combined_ids = list(args.uniprot_id)
        combined_names = list(args.uniprot_name)
        for entry in file_queries:
            combined_ids.extend(entry["ids"])
            combined_names.extend(entry["names"])
        queries = [{"ids": combined_ids, "names": combined_names}]
    return queries


def collect_identifiers(queries: List[dict]) -> tuple[set[str], set[str]]:
    ids: set[str] = set()
    names: set[str] = set()
    for query in queries:
        ids.update(query.get("ids") or [])
        names.update(query.get("names") or [])
    return ids, names


def derive_output_stem(ids: List[str], names: List[str], fallback: str = "ppi") -> str:
    parts = list(ids) + list(names)
    if not parts:
        parts.append(fallback)
    slug = re.sub(r"[^A-Za-z0-9]+", "-", "-".join(parts)).strip("-")
    return (slug or fallback).lower()


def unique_stem(base: str, used: set[str]) -> str:
    stem = base
    counter = 2
    while stem in used:
        stem = f"{base}-{counter}"
        counter += 1
    used.add(stem)
    return stem


def resolve_output_path(user_value: str | None, stem: str, suffix: str, per_query: bool, output_dir: Path) -> Path:
    if user_value and not per_query:
        return Path(user_value)
    return output_dir / f"{stem}{suffix}"


def fetch_json(url: str, timeout: int = 20) -> dict:
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout, context=SSL_CONTEXT) as response:
        return json.load(response)


def fetch_uniprot_functions(nodes: Iterable[str], enabled: bool) -> Dict[str, str]:
    functions: Dict[str, str] = {}
    if not enabled:
        return functions

    for node in sorted(nodes):
        url = f"https://rest.uniprot.org/uniprotkb/{node}.json"
        print(f"Fetching UniProt Function annotation for {node}...", file=sys.stderr)
        try:
            data = fetch_json(url)
        except urllib.error.URLError as exc:
            message = getattr(exc, "reason", exc)
            print(f"  Failed to fetch {node}: {message}", file=sys.stderr)
            continue
        except Exception as exc:  # noqa: BLE001
            print(f"  Unexpected error fetching {node}: {exc}", file=sys.stderr)
            continue

        comment = ""
        for entry in data.get("comments", []):
            if entry.get("commentType") == "FUNCTION":
                texts = entry.get("texts") or []
                if texts:
                    comment = texts[0].get("value", "")
                    break
        if comment:
            functions[node] = comment.strip()
        else:
            print(f"  No Function annotation available for {node}.", file=sys.stderr)
    return functions


def collect_node_table(
    adjacency,
    node_labels,
    seeds,
    functions: Dict[str, str] | None = None,
    node_metrics: Dict[str, Dict[str, float]] | None = None,
):
    for node in sorted(adjacency.keys()):
        data = {
            "id": node,
            "label": preferred_label(node_labels.get(node, []), node),
            "seed": str(node in seeds),
            "degree": str(len(adjacency[node])),
            "function": (functions or {}).get(node, ""),
        }
        if node_metrics:
            for metric_name, values in node_metrics.items():
                metric_value = values.get(node)
                if metric_value is not None:
                    data[metric_name] = metric_value
        yield data


def parse_numeric(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def aggregate_node_metric(edges, field: str) -> Dict[str, float]:
    if not field:
        return {}
    buckets: Dict[str, List[float]] = defaultdict(list)
    for a, b, row in edges:
        value = parse_numeric(row.get(field))
        if value is None:
            continue
        buckets[a].append(value)
        buckets[b].append(value)
    return {node: sum(values) / len(values) for node, values in buckets.items() if values}


def scale_values(values: Dict[str, float], min_scale: float, max_scale: float) -> Dict[str, float]:
    if not values:
        return {}
    min_value = min(values.values())
    max_value = max(values.values())
    if max_value == min_value:
        return {node: (min_scale + max_scale) / 2.0 for node in values}

    def scale(val: float) -> float:
        return min_scale + ((val - min_value) / (max_value - min_value)) * (max_scale - min_scale)

    return {node: scale(val) for node, val in values.items()}


def build_color_map(values: Dict[str, float]):
    if not values:
        return {}, None, None
    min_value = min(values.values())
    max_value = max(values.values())
    if max_value == min_value:
        max_value += 1.0

    def lerp(a: float, b: float, t: float) -> float:
        return a + (b - a) * t

    def hex_component(value: float) -> str:
        return f"{int(round(value)):02x}"

    colors = {}
    for node, value in values.items():
        t = (value - min_value) / (max_value - min_value)
        r = lerp(31, 214, t)
        g = lerp(119, 96, t)
        b = lerp(180, 77, t)
        colors[node] = f"#{hex_component(r)}{hex_component(g)}{hex_component(b)}"
    return colors, min_value, max_value


def build_edge_color_series(edges, field: str):
    values = [parse_numeric(row.get(field)) for _, _, row in edges]
    numeric_values = [value for value in values if value is not None]
    if not numeric_values:
        return ["#c0c0c0"] * len(edges), None, None
    min_value = min(numeric_values)
    max_value = max(numeric_values)
    if max_value == min_value:
        max_value += 1.0

    silver = (192, 192, 192)
    blue = (31, 119, 180)

    def to_color(value: float | None) -> str:
        if value is None:
            return "#c0c0c0"
        t = (value - min_value) / (max_value - min_value)
        r = int(silver[0] + (blue[0] - silver[0]) * t)
        g = int(silver[1] + (blue[1] - silver[1]) * t)
        b = int(silver[2] + (blue[2] - silver[2]) * t)
        return f"#{r:02x}{g:02x}{b:02x}"

    return [to_color(value) for value in values], min_value, max_value


def edge_numeric_series(edges, field: str):
    values = [parse_numeric(row.get(field)) for _, _, row in edges]
    numeric_values = [value for value in values if value is not None]
    if not numeric_values:
        return values, None, None
    return values, min(numeric_values), max(numeric_values)


def validate_uniprot_entries(ids: Iterable[str], names: Iterable[str]) -> None:
    missing: List[str] = []

    for accession in sorted({value.strip() for value in ids if value}):
        url = UNIPROT_ENTRY_URL.format(identifier=urllib.parse.quote(accession))
        try:
            fetch_json(url)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                missing.append(f"UniProt ID not valid or not found: {accession}")
            else:
                raise

    for name in sorted({value.strip() for value in names if value}):
        query = urllib.parse.quote(f"id:{name}")
        url = UNIPROT_SEARCH_URL.format(query=query)
        try:
            data = fetch_json(url)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                missing.append(f"{name} is not a valid UniProt entry name!")
            else:
                raise
        else:
            if not data.get("results"):
                missing.append(f"{name} is not a valid UniProt entry name!")

    if missing:
        raise SystemExit("Input validation failed:\n" + "\n".join(missing))


def prompt_with_timeout(prompt: str, timeout: int = 30, default: str = "") -> str:
    try:
        print(prompt, end="", flush=True)
        rlist, _, _ = select.select([sys.stdin], [], [], timeout)
        if rlist:
            return sys.stdin.readline().strip()
    except Exception:
        pass
    print("")  # newline after timeout or error
    return default


def gradient_block(label: str, start_color: str, end_color: str, min_val: float, max_val: float) -> str:
    return f"""
    <div style="margin-top:8px;">
      <div>{label}</div>
      <div style="background: linear-gradient(90deg, {start_color}, {end_color});
                  height:12px; width:200px; border:1px solid #999; border-radius:4px;"></div>
      <div style="display:flex; justify-content:space-between; width:200px; font-size:12px;">
        <span>{min_val:.3g}</span>
        <span>{max_val:.3g}</span>
      </div>
    </div>
    """


def width_block(label: str, min_val: float, max_val: float) -> str:
    return f"""
    <div style="margin-top:8px;">
      <div>{label}</div>
      <div style="display:flex; align-items:center; gap:12px; font-size:12px;">
        <div style="height:2px; width:50px; background:#555;"></div><span>{min_val:.3g}</span>
        <div style="height:8px; width:50px; background:#555;"></div><span>{max_val:.3g}</span>
      </div>
    </div>
    """


def write_json_network(
    adjacency,
    node_labels,
    edges,
    seeds,
    label_field: str,
    path: Path,
    functions: Dict[str, str],
    node_metrics: Dict[str, Dict[str, float]],
):
    nodes = list(
        collect_node_table(adjacency, node_labels, seeds, functions, node_metrics)
    )
    json_edges = [
        {
            "source": a,
            "target": b,
            "label": edge_label(row, label_field),
            "attributes": row,
        }
        for a, b, row in edges
    ]
    payload = {"nodes": nodes, "edges": json_edges}
    path.write_text(json.dumps(payload, indent=2))
    print(f"JSON network written to {path}", file=sys.stderr)


def write_pyvis_html(
    adjacency,
    node_labels,
    edges,
    seeds,
    label_field: str,
    path: Path,
    functions: Dict[str, str],
    color_map: Dict[str, str],
    size_map: Dict[str, float],
    edge_colors: List[str],
    width_scores: List[float | None],
    width_range: tuple[float | None, float | None],
    legend_html: str,
    node_urls: Dict[str, str],
    use_peak_matrix: bool = False,
    title_text: str | None = None,
):
    try:
        from pyvis.network import Network
    except ImportError:
        print(
            "PyVis is not installed; skipping PyVis HTML export. Install via 'pip install pyvis'.",
            file=sys.stderr,
        )
        return

    net = Network(height="750px", width="100%", notebook=False)
    net.toggle_stabilization(True)
    options = {
        "nodes": {"shape": "dot", "size": 18, "font": {"size": 16}},
        "edges": {
            "smooth": True,
            "color": {
                "opacity": 0.6,
            },
        },
        "physics": {"stabilization": {"iterations": 200}},
    }
    net.set_options(json.dumps(options))

    min_width, max_width = width_range

    def edge_width(score: float | None) -> float:
        if score is None or max_width is None:
            return 2.0
        baseline = min_width or 0.0
        if max_width == baseline:
            return 6.0
        norm = (score - baseline) / (max_width - baseline)
        return 2.0 + norm * 6.0

    for node, neighbors in adjacency.items():
        label = preferred_label(node_labels.get(node, []), node)
        function_text = functions.get(node, "").replace("\n", " ")
        title_lines = [f"{label} ({node})", f"Degree: {len(neighbors)}"]
        if function_text:
            title_lines.append("")
            title_lines.append(f"Function: {function_text}")
        title = "\n".join(title_lines)
        color_value = color_map.get(node) if color_map else None
        default_color = "#ff7f0e" if node in seeds else "#1f77b4"
        node_color = color_value or default_color
        color_config = {"background": node_color}
        border_width = 1
        if node in seeds:
            color_config["border"] = "#000000"
            border_width = 3
        size = size_map.get(node, 18.0) if size_map else 18.0
        size = size_map.get(node, 18.0) if size_map else 18.0
        net.add_node(
            node,
            label=label,
            title=title,
            color=color_config,
            value=size,
            size=size,
            borderWidth=border_width,
            borderWidthSelected=border_width + 1,
        )

    for (source, target, row), score, color in zip(edges, width_scores, edge_colors):
        label = edge_label(row, label_field) or ""
        if use_peak_matrix:
            uni_name = row.get("UniProt_Name") or f"{source} / {target}"
            peak_val = row.get("peak_max_norm", row.get("peak_avg_norm", ""))
            iptm_val = row.get("ipTM_max", row.get("iptm_max", ""))
            title = "\n".join(
                [
                    str(uni_name),
                    f"peak_max_norm: {peak_val}",
                    f"ipTM_max: {iptm_val}",
                ]
            )
        else:
            title = "\n".join(f"{k}: {v}" for k, v in row.items() if v)
        net.add_edge(
            source,
            target,
            label=label,
            title=title,
            width=edge_width(score),
            color=color,
            value=score if score is not None else 0,
        )

    net.write_html(str(path), notebook=False)
    tip_html = (
        "<div style=\"padding:10px;font-size:14px;max-width:420px;font-weight:bold;\">"
        "Double-click a node to open its UniProt entry."
        "</div>"
                
        "<div style=\"padding:10px;font-size:14px;max-width:420px;font-weight:bold;\">"

        "Hover over an egde to see prediction scores."
        "</div>"

        "<div style=\"padding:10px;font-size:14px;max-width:420px;font-weight:bold;\">"

        "Drag a node to re-arrange the network layout."
        "</div>"

    )
    extra_html = tip_html
    if legend_html:
        extra_html += (
            '<div style="padding:10px;font-size:14px;max-width:240px;">'
            f"<b>Legend:</b><br>{legend_html}</div>"
        )
    script = ""
    if node_urls:
        script = (
            "<script>"
            "window.addEventListener('load', function() {"
            f"const nodeUrlMap = {json.dumps(node_urls)};"
            "if (typeof network === 'undefined') return;"
            "network.on('doubleClick', function(params) {"
            "if (params.nodes.length) {"
            "const nodeId = params.nodes[0];"
            "const url = nodeUrlMap[nodeId];"
            "if (url) window.open(url, '_blank');"
            "}"
            "});"
            "});"
            "</script>"
        )
    header_html = ""
    if title_text:
        header_html = f'<div style="padding:10px;font-size:18px;font-weight:bold;">{title_text}</div>'
    contents = path.read_text()
    contents = contents.replace("</body>", f"{header_html}{extra_html}{script}</body>")
    path.write_text(contents)
    print(f"PyVis HTML network written to {path}", file=sys.stderr)


def process_query(
    rows: List[dict],
    fieldnames: Sequence[str],
    ids: List[str],
    names: List[str],
    args,
    stem: str,
    write_stdout: bool,
):
    matches_found = []
    for row in rows:
        id_match = matches_any(row.get("UniProt_ID", ""), ids) if ids else False
        name_match = matches_any(row.get("UniProt_Name", ""), names) if names else False
        # Match if IDs or names hit (not strictly both)
        if (ids or names) and (id_match or name_match):
            matches_found.append(row)

    print(
        f"[info] Query IDs={ids or '-'} names={names or '-'} -> {len(matches_found)} match(es)",
        file=sys.stderr,
    )
    if not matches_found:
        suggested = suggest_matches(ids + names, rows)
        msg = f"No matches found for query IDs={ids or '-'} names={names or '-'}"
        if suggested:
            msg += f". Closest matches: {', '.join(suggested)}"
        print(msg, file=sys.stderr)
        return

    if write_stdout:
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(matches_found)
    else:
        csv_path = Path(f"{stem}.csv")
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(matches_found)
        print(f"CSV matches written to {csv_path}", file=sys.stderr)

    initial_seeds = set(ids)
    if not initial_seeds and names:
        target_names = {name.lower() for name in names}
        name_based: set[str] = set()
        for row in matches_found:
            id_tokens = split_tokens(row.get("UniProt_ID", ""))
            name_tokens = split_tokens(row.get("UniProt_Name", ""))
            for tok_id, tok_name in zip(id_tokens, name_tokens):
                if tok_name.lower() in target_names:
                    name_based.add(tok_id)
        if name_based:
            initial_seeds = name_based
    if not initial_seeds:
        initial_seeds = collect_seed_ids(matches_found)
    network_rows, _expanded_seeds = expand_seed_network(rows, initial_seeds, args.expansion_depth)
    adjacency, node_labels, edges = build_graph(network_rows)
    render_network_summary(adjacency, node_labels, initial_seeds, args.edge_label_field)

    node_count = len(adjacency)
    print(
        f"[info] Expanded to {len(network_rows)} interaction rows, {node_count} nodes, {len(edges)} edges",
        file=sys.stderr,
    )
    large_graph = node_count > args.large_network_threshold

    if large_graph and not args.no_function_fetch:
        print(
            f"Warning: the network contains {node_count} nodes (> {args.large_network_threshold}). "
            "Consider reducing --expansion-depth if this is unintended.",
            file=sys.stderr,
        )
        confirm = prompt_with_timeout(
            "Fetch UniProt Function annotations for this large network? [y/N] ",
            timeout=30,
            default="",
        ).lower()
        if confirm not in {"y", "yes"}:
            args.no_function_fetch = True

    functions = fetch_uniprot_functions(adjacency.keys(), not args.no_function_fetch)
    degree_map = {node: len(neighbors) for node, neighbors in adjacency.items()}
    color_values = degree_map
    size_values = aggregate_node_metric(edges, args.node_size_field)
    node_metrics: Dict[str, Dict[str, float]] = {}
    legend_blocks: List[str] = []

    node_metrics["degree"] = color_values
    color_map, color_min, color_max = build_color_map(color_values)
    legend_blocks.append(
        gradient_block(
            "Node color (degree)",
            "#cfe2f3",
            "#1c4587",
            color_min or 0,
            color_max or 0,
        )
    )

    size_map = {}
    if size_values:
        node_metrics[f"{args.node_size_field}_avg"] = size_values
        size_map = scale_values(size_values, 12.0, 28.0)
        size_min = min(size_values.values())
        size_max = max(size_values.values())
        legend_blocks.append(
            f"<div style='margin-top:8px;'>Node size ({args.node_size_field}) average range "
            f"{size_min:.3g}-{size_max:.3g}</div>"
        )

    edge_colors, edge_color_min, edge_color_max = build_edge_color_series(
        edges, args.edge_color_field
    )
    if edge_color_min is not None and edge_color_max is not None:
        legend_blocks.append(
            gradient_block(
                f"Edge color ({args.edge_color_field})",
                "#c0c0c0",
                "#1f77b4",
                edge_color_min,
                edge_color_max,
            )
        )

    width_scores, width_min, width_max = edge_numeric_series(edges, args.edge_width_field)
    if width_min is not None and width_max is not None:
        legend_blocks.append(
            width_block(f"Edge width ({args.edge_width_field})", width_min, width_max)
        )
    legend_blocks.append(
        f"<div style='margin-top:8px;'><b>Expansion depth:</b> {args.expansion_depth}</div>"
    )
    if initial_seeds:
        legend_blocks.append("<div style='margin-top:8px;'>Seed nodes outlined in black</div>")
    legend_blocks.append("<div style='margin-top:8px;'>Double-click a node to open its UniProt entry.</div>")
    legend_html = "".join(legend_blocks)
    node_urls = {node: f"https://www.uniprot.org/uniprotkb/{node}" for node in adjacency.keys()}

    id_label = ", ".join(ids) if ids else "-"
    name_label = ", ".join(names) if names else "-"
    metric_label = "Predictome (peak score)" if args.metric == "peak_score" else "SPOC (FDR)"
    title_text = f"{metric_label} network for IDs={id_label} names={name_label}"
    json_path = resolve_output_path(args.json_out, stem, ".json", args.per_query_output, args.output_dir)
    pyvis_path = resolve_output_path(args.pyvis_out, stem, "_network.html", args.per_query_output, args.output_dir)

    if not args.no_json:
        write_json_network(
            adjacency,
            node_labels,
            edges,
            initial_seeds,
            args.edge_label_field,
            json_path,
            functions,
            node_metrics,
        )

    if not args.no_pyvis:
        write_pyvis_html(
            adjacency,
            node_labels,
            edges,
            initial_seeds,
            args.edge_label_field,
            pyvis_path,
            functions,
            color_map,
            size_map,
            edge_colors,
            width_scores,
            (width_min, width_max),
            legend_html,
            node_urls,
            use_peak_matrix=args.metric == "peak_score",
            title_text=title_text,
        )


def main():
    args = parse_args()
    if args.metric == "peak_score":
        matrix_path = find_peak_matrix(args.predictome_matrix)
        if not matrix_path.exists():
            sys.exit(f"Peak matrix file not found: {matrix_path}")
        fieldnames, rows = load_peak_matrix(
            matrix_path,
            args.predictome_peak_threshold,
            args.predictome_iptm_threshold,
            args.predictome_peak_ceiling,
        )
        print(
            f"[info] Using predictome matrix: {matrix_path} "
            f"(peak>={args.predictome_peak_threshold}, iptm>={args.predictome_iptm_threshold})",
            file=sys.stderr,
        )
        source_label = "peak_score"
    else:
        csv_path = Path(args.csv_path)
        if not csv_path.exists():
            sys.exit(f"CSV file not found: {csv_path}")
        fieldnames, rows = load_rows(csv_path)
        print(f"[info] Using SPOC CSV source: {csv_path}", file=sys.stderr)
        source_label = "FDR"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    queries = build_queries(args)
    if not any(query["ids"] or query["names"] for query in queries):
        sys.exit("Provide at least one UniProt ID or name (CLI or query file).")

    all_ids, all_names = collect_identifiers(queries)
    validate_uniprot_entries(all_ids, all_names)

    used_stems: set[str] = set()
    write_stdout = not args.per_query_output
    for idx, query in enumerate(queries, start=1):
        ids = query["ids"]
        names = query["names"]
        if not ids and not names:
            continue
        base_stem = derive_output_stem(ids, names)
        stem = unique_stem(f"{source_label}-{base_stem}", used_stems)
        process_query(rows, fieldnames, ids, names, args, stem, write_stdout)
        write_stdout = False
        print(
            f"[{idx}/{len(queries)}] Processed query with expansion depth {args.expansion_depth}. "
            "Increase --expansion-depth to include more interaction layers.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
