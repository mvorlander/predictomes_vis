# predictomes-vis (containers)

Toolkit for visualizing predictomes.org confidence scores:

- **network** mode: build PPI networks from either
  - **FDR** (default): high-confidence SPOC results (<10% FDR)
  - **peak_score**: full predictome matrix (~1.6M pairs parsed by M. Vorländer, Plaschka lab, IMP Vienna)
- **binary** mode: per-POI plots (ptm/iptm/peak scatter, beeswarm, histograms) from the predictome matrix.

## Running on the host (Apptainer wrapper, recommended)

The container bundles both data sources:
- `20251113_download_hsbps.csv` (high-confidence SPOC, FDR mode)
- `predictome_json_stats_251115/json_metrics_matrix.tsv` (predictome matrix, peak_score mode)

Wrapper on the host (after `build_and_copy_predictomes_vis.sh` has uploaded and built the `.sif`):
```bash
# Network, FDR (default), direct neighbors:
/groups/plaschka/shared/software/predictomes/wrappers/predictomes_vis.sh \
  network --metric FDR --poi A1AG2_HUMAN --expansion-depth 0

# Network, peak_score with thresholds:
/groups/plaschka/shared/software/predictomes/wrappers/predictomes_vis.sh \
  network --metric peak_score --poi A1AG2_HUMAN \
  --predictome-peak-threshold 0.75 --predictome-iptm-threshold 0.5 \
  --expansion-depth 1

# Binary plots (per-POI metrics):
/groups/plaschka/shared/software/predictomes/wrappers/predictomes_vis.sh \
  binary --poi A1AG2_HUMAN --outdir /tmp/predictome_binary_plots
```

Defaults:
- Output dirs: `outputs/ppi_networks` (network) and `outputs/predictome_binary_plots` (binary)
- DOT/PNG disabled by default; JSON/PyVis enabled for network outputs.
- Predictome matrix path is set via `PREDICTOMES_PEAK_MATRIX` inside the container.

## Building and deploying the container

From the repo root:
```bash
bash scripts/build_and_copy_predictomes_vis.sh
```
This builds the Docker image (linux/amd64), saves a tar, uploads it and the wrapper to the host, and builds the Apptainer `.sif` on the host (`/groups/plaschka/shared/software/predictomes/containers/predictomes_vis.sif`).

## Local usage (non-container)

Install deps and run:
```bash
python3 -m predictomes_vis.cli network --metric FDR --poi A1AG2_HUMAN --expansion-depth 0
python3 -m predictomes_vis.cli binary --poi A1AG2_HUMAN --outdir outputs/predictome_binary_plots
```
You need `certifi` (network mode), `pandas/plotly/kaleido` (binary), and Graphviz for PNG.

## Notes
- POI can be UniProt ID or entry name; unmatched POIs suggest closest matches.
- Large network function-fetch prompt times out after 30s and skips by default.
- When using peak_score metric, edges default to width = `peak_avg_norm`, color = `iptm_max`.
