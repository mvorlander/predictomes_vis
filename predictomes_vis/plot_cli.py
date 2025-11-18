#!/usr/bin/env python3
"""
Generate plots from the AlphaFold `json_metrics_matrix.tsv` output for a
protein of interest (POI).
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create interaction plots for a given POI "
        "from json_metrics_matrix.tsv."
    )
    parser.add_argument(
        "--matrix",
        default=None,
        type=Path,
        help="Path to json_metrics_matrix.tsv (tab-separated).",
    )
    parser.add_argument(
        "--poi",
        nargs="+",
        help="Protein(s) of interest (entry name preferred). Space-separated allowed.",
    )
    parser.add_argument(
        "--poi-file",
        type=Path,
        help="Optional text file with one POI per line (entry name or accession).",
    )
    parser.add_argument(
        "--outdir",
        default=Path("outputs/predictome_binary_plots"),
        type=Path,
        help="Directory where plots will be written.",
    )
    parser.add_argument(
        "--peak-ceiling",
        type=float,
        default=30.0,
        help="Upper bound used to normalise peak scores (1 - score/ceiling). "
        "The default mirrors the R reference script (30).",
    )
    parser.add_argument(
        "--csv-out",
        type=Path,
        default=None,
        help="Optional path to write the filtered/derived table as CSV. "
        "Defaults to <outdir>/<poi>_filtered.csv.",
    )
    parser.add_argument(
        "--iptm-cutoff",
        type=float,
        default=0.50,
        help="Horizontal guide line for IPTM beeswarm and histogram.",
    )
    parser.add_argument(
        "--peak-cutoff",
        type=float,
        default=0.75,
        help="Horizontal guide line for Peak beeswarm and histogram.",
    )
    return parser.parse_args()


def ensure_columns(df: pd.DataFrame, cols: Iterable[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing expected columns: {', '.join(missing)}")


def looks_like_uniprot_id(text: str) -> bool:
    import re

    core = text.split("-")[0]
    patterns = [
        r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$",
        r"^[A-NR-Z][0-9]{2}[A-Z][A-Z0-9]{2}[0-9]$",
        r"^[A-Z0-9]{10}$",
    ]
    return any(re.match(pat, core) for pat in patterns)


def looks_like_entry_name(text: str) -> bool:
    import re

    return bool(re.match(r"^[A-Z0-9\\-]{1,16}_[A-Z0-9\\-]{1,16}$", text))


def resolve_poi(raw_poi: str, df: pd.DataFrame) -> str:
    poi = raw_poi.strip().upper()
    proteins = pd.unique(pd.concat([df["protein_a"], df["protein_b"]]))
    if poi in proteins:
        return poi
    base = poi.split("-")[0]
    if base in proteins:
        print(f"POI '{poi}' not found; using base entry '{base}'.")
        return base
    resembles_id = looks_like_uniprot_id(poi)
    resembles_entry = looks_like_entry_name(poi)
    suggestions = []
    try:
        import difflib

        suggestions = difflib.get_close_matches(poi, proteins.tolist(), n=5, cutoff=0.4)
    except Exception:
        suggestions = []
    suggestion_txt = f" Closest matches: {', '.join(suggestions)}" if suggestions else ""
    if resembles_id and not resembles_entry:
        raise SystemExit(
            f"POI '{raw_poi}' looks like a UniProt accession but is not present in the matrix."
            " Please supply the UniProt entry name (e.g., A1AG2_HUMAN)."
            f"{suggestion_txt}"
        )
    raise SystemExit(
        f"POI '{raw_poi}' not found in protein_a/protein_b columns.{suggestion_txt}"
    )


def gather_pois(args: argparse.Namespace) -> list[str]:
    pois: list[str] = []
    if args.poi:
        pois.extend(args.poi)
    if args.poi_file:
        if not args.poi_file.exists():
            raise SystemExit(f"POI file not found: {args.poi_file}")
        with args.poi_file.open() as fh:
            pois.extend([line.strip() for line in fh if line.strip()])
    if not pois:
        raise SystemExit("No POI provided. Use --poi or --poi-file.")
    seen = set()
    uniq = []
    for p in pois:
        if p not in seen:
            seen.add(p)
            uniq.append(p)
    return uniq


def get_default_matrix() -> Path:
    env_path = os.environ.get("PREDICTOMES_DEFAULT_MATRIX")
    if env_path:
        return Path(env_path)
    candidate = Path(__file__).resolve().parents[1] / "predictome_data" / "json_metrics_matrix.tsv"
    return candidate


def prepare_dataframe(df: pd.DataFrame, poi: str, peak_ceiling: float) -> pd.DataFrame:
    keep = (df["protein_a"] == poi) | (df["protein_b"] == poi)
    filtered = df.loc[keep].copy()
    if filtered.empty:
        raise SystemExit(f"No rows found for POI '{poi}'.")

    peak_cols = [f"peak_score_rank_{i}" for i in (1, 2, 3)]
    ptm_cols = [f"ptm_rank_{i}" for i in (1, 2, 3)]
    iptm_cols = [f"iptm_rank_{i}" for i in (1, 2, 3)]

    ensure_columns(filtered, peak_cols + ptm_cols + iptm_cols)

    filtered["bait"] = poi
    filtered["partner"] = filtered.apply(
        lambda row: row["protein_b"] if row["protein_a"] == poi else row["protein_a"],
        axis=1,
    )

    peaks_norm = 1 - filtered[peak_cols].astype(float) / peak_ceiling
    filtered["peak_avg"] = peaks_norm.mean(axis=1)
    filtered["peak_max_norm"] = peaks_norm.max(axis=1)
    filtered["peak_min_norm"] = peaks_norm.min(axis=1)
    filtered["peak_values_norm"] = peaks_norm.round(3).astype(str).agg(":".join, axis=1)

    filtered["ptm_avg"] = filtered[ptm_cols].mean(axis=1)
    filtered["iptm_avg"] = filtered[iptm_cols].mean(axis=1)

    filtered["ptm_values"] = filtered[ptm_cols].round(3).astype(str).agg(":".join, axis=1)
    filtered["iptm_values"] = filtered[iptm_cols].round(3).astype(str).agg(":".join, axis=1)

    filtered["uniprot_url"] = "https://www.uniprot.org/uniprotkb/" + filtered["partner"] + "/entry"
    return filtered


def save_plot(fig: go.Figure, outdir: Path, name: str, write_html: bool = True) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if write_html:
        html_path = outdir / f"{name}.html"
        fig.write_html(html_path)
        print(f"Wrote {html_path}")

    for ext in ("png", "pdf"):
        img_path = outdir / f"{name}.{ext}"
        try:
            fig.write_image(img_path)
        except Exception as exc:  # noqa: BLE001
            print(f"Skipping {img_path} ({exc})")
            continue
        print(f"Wrote {img_path}")


def scatter_ptm_vs_iptm(df: pd.DataFrame, poi: str) -> go.Figure:
    fig = px.scatter(
        df,
        x="ptm_avg",
        y="iptm_avg",
        color="iptm_max",
        labels={"iptm_max": "max IPTM (color)", "ptm_avg": "PTMavg", "iptm_avg": "IPTMavg"},
        hover_name="partner",
        hover_data={
            "bait": True,
            "partner": True,
            "rank_count": True,
            "peak_values_norm": True,
            "ptm_values": True,
            "iptm_values": True,
            "ptm_avg": ":.3f",
            "iptm_avg": ":.3f",
        },
        range_x=[0, 1],
        range_y=[0, 1],
        color_continuous_scale="Viridis",
        title=f"{poi}: PTMavg vs IPTMavg",
    )
    fig.update_layout(coloraxis_reversescale=True)
    fig.update_traces(marker=dict(size=10, opacity=0.8, line=dict(width=0.5, color="white")))
    return fig


def scatter_peak_vs_iptm(df: pd.DataFrame, poi: str) -> go.Figure:
    max_size = float(df["iptm_max"].max()) if not df.empty else 0.0
    sizeref = 2.0 * max_size / (18.0 ** 2) if max_size > 0 else None

    fig = px.scatter(
        df,
        x="peak_avg",
        y="iptm_avg",
        color="peak_max_norm",
        size="iptm_max",
        labels={
            "peak_max_norm": "max Peak (color)",
            "iptm_max": "max IPTM (size)",
            "peak_avg": "Peakavg (norm.)",
            "iptm_avg": "IPTMavg",
        },
        hover_name="partner",
        hover_data={
            "bait": True,
            "partner": True,
            "peak_values_norm": True,
            "iptm_values": True,
            "peak_avg": ":.3f",
            "iptm_avg": ":.3f",
        },
        range_x=[0, 1],
        range_y=[0, 1],
        color_continuous_scale="Viridis",
        title=f"{poi}: Peak (norm.) vs IPTMavg",
    )
    fig.update_layout(
        coloraxis_reversescale=True,
        annotations=[
            dict(
                text="Size = max IPTM; double-click a dot to open UniProt.",
                xref="paper",
                yref="paper",
                x=0,
                y=1.12,
                showarrow=False,
                font=dict(size=12),
            )
        ],
    )
    marker_cfg = dict(
        opacity=0.7,
        line=dict(width=0.5, color="white"),
        sizemode="area",
        sizemin=6,
    )
    if sizeref:
        marker_cfg["sizeref"] = sizeref
    fig.update_traces(marker=marker_cfg)
    return fig


def beeswarm_metric(
    df: pd.DataFrame,
    poi: str,
    metric: str,
    max_metric: str,
    cutoff: float,
    title: str,
) -> go.Figure:
    tmp = df.copy()
    tmp["x_jitter"] = compute_swarm_offsets_circle(tmp[metric].to_numpy(), radius=0.025)

    fig = px.scatter(
        tmp,
        x="x_jitter",
        y=metric,
        color=max_metric,
        labels={max_metric: f"max {metric.upper()} (color)", metric: metric},
        hover_name="partner",
        hover_data={
            "bait": True,
            "partner": True,
            "peak_values_norm": True,
            "ptm_values": True,
            "iptm_values": True,
            metric: ":.3f",
        },
        color_continuous_scale="Viridis",
        title=title,
    )
    fig.update_layout(
        coloraxis_reversescale=True,
        xaxis=dict(
            title="",
            tickmode="array",
            tickvals=[0],
            ticktext=[poi],
            range=[-0.22, 0.22],
        ),
        yaxis=dict(range=[0, 1]),
        showlegend=True,
    )
    fig.add_hline(y=cutoff, line_width=1, line_dash="dot", line_color="red")
    fig.update_traces(marker=dict(size=9, opacity=0.9, line=dict(width=0.4, color="white")))
    return fig


def histogram_metric(df: pd.DataFrame, metric: str, cutoff: float, title: str) -> go.Figure:
    fig = px.histogram(
        df,
        x=metric,
        nbins=25,
        range_x=[0, 1],
        title=title,
        color_discrete_sequence=["#4c78a8"],
    )
    fig.add_vline(x=cutoff, line_width=2, line_dash="dot", line_color="red")
    fig.update_layout(bargap=0.05)
    return fig


def compute_swarm_offsets_circle(values: np.ndarray, radius: float = 0.03) -> np.ndarray:
    offsets = np.zeros_like(values, dtype=float)
    order = np.argsort(values)
    placed: list[tuple[float, float]] = []

    for idx in order:
        y = float(values[idx])
        if not placed:
            offsets[idx] = 0.0
            placed.append((0.0, y))
            continue

        step = radius * 0.9
        k = 0
        while True:
            candidates = [k * step, -k * step] if k else [0.0]
            for cand in candidates:
                if all((cand - px) ** 2 + (y - py) ** 2 >= radius**2 for px, py in placed):
                    offsets[idx] = cand
                    placed.append((cand, y))
                    k = None
                    break
            if k is None:
                break
            k += 1
    return offsets


def add_uniprot_click(fig: go.Figure, div_id: str) -> str:
    fig.update_layout(
        clickmode="event+select",
        annotations=[
            dict(
                text="Double click a dot to open the UniProt entry (opens new tab).",
                xref="paper",
                yref="paper",
                x=0,
                y=1.08,
                showarrow=False,
                font=dict(size=12),
            )
        ],
    )
    script = f"""
    const plot = document.getElementById('{div_id}');
    plot.on('plotly_click', function(data) {{
      const pt = data.points[0];
      const url = pt?.customdata?.[0];
      if (url) {{
        window.open(url, '_blank');
      }}
    }});
    """
    return script


def write_html_with_click(fig: go.Figure, outdir: Path, name: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    div_id = f"{name}-plot"
    post_script = add_uniprot_click(fig, div_id)
    html_path = outdir / f"{name}.html"
    fig.write_html(html_path, div_id=div_id, post_script=post_script)
    print(f"Wrote {html_path}")


def main(args: argparse.Namespace | None = None) -> None:
    if args is None:
        args = parse_args()

    matrix_path = args.matrix or get_default_matrix()
    if not matrix_path.exists():
        raise SystemExit(f"Matrix file not found: {matrix_path}")

    df = pd.read_csv(matrix_path, sep="\t")
    pois = gather_pois(args)
    args.outdir.mkdir(parents=True, exist_ok=True)

    for poi_in in pois:
        resolved_poi = resolve_poi(poi_in, df)
        prepared = prepare_dataframe(df, poi=resolved_poi, peak_ceiling=args.peak_ceiling)

        csv_path = args.csv_out or (args.outdir / f"{resolved_poi}_filtered.csv")
        prepared.to_csv(csv_path, index=False)
        print(f"[{resolved_poi}] Wrote {csv_path}")

        scatter_ptm = scatter_ptm_vs_iptm(prepared, poi=resolved_poi)
        scatter_ptm.update_traces(customdata=prepared[["uniprot_url"]].to_numpy())
        base = f"{resolved_poi}_ptm_vs_iptm"
        write_html_with_click(scatter_ptm, args.outdir, base)
        save_plot(scatter_ptm, args.outdir, base, write_html=False)

        scatter_peak = scatter_peak_vs_iptm(prepared, poi=resolved_poi)
        scatter_peak.update_traces(customdata=prepared[["uniprot_url"]].to_numpy())
        base = f"{resolved_poi}_peak_vs_iptm"
        write_html_with_click(scatter_peak, args.outdir, base)
        save_plot(scatter_peak, args.outdir, base, write_html=False)

        save_plot(
            beeswarm_metric(
                prepared,
                poi=resolved_poi,
                metric="peak_avg",
                max_metric="peak_max_norm",
                cutoff=args.peak_cutoff,
                title=f"{resolved_poi}: Peak (norm.) distribution",
            ),
            args.outdir,
            f"{resolved_poi}_peak_only",
        )
        save_plot(
            beeswarm_metric(
                prepared,
                poi=resolved_poi,
                metric="iptm_avg",
                max_metric="iptm_max",
                cutoff=args.iptm_cutoff,
                title=f"{resolved_poi}: IPTM distribution",
            ),
            args.outdir,
            f"{resolved_poi}_iptm_only",
        )

        save_plot(
            histogram_metric(
                prepared,
                metric="peak_avg",
                cutoff=args.peak_cutoff,
                title=f"{resolved_poi}: Peak (norm.) histogram",
            ),
            args.outdir,
            f"{resolved_poi}_peak_hist",
        )
        save_plot(
            histogram_metric(
                prepared,
                metric="iptm_avg",
                cutoff=args.iptm_cutoff,
                title=f"{resolved_poi}: IPTM histogram",
            ),
            args.outdir,
            f"{resolved_poi}_iptm_hist",
        )


if __name__ == "__main__":
    main()
