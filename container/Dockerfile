FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /opt/predictomes

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
# retry install to mitigate transient apt/dpkg errors
RUN set -euo pipefail; \
    tries=3; \
    while (( tries > 0 )); do \
      apt-get update && apt-get install -y --no-install-recommends build-essential && break; \
      tries=$((tries-1)); \
      echo "apt install failed, retrying... ($tries left)"; \
      sleep 5; \
    done; \
    rm -rf /var/lib/apt/lists/*

# Copy code
COPY predictomes_vis/predictomes_vis /opt/predictomes/predictomes_vis
COPY predictomes_vis/pyproject.toml /opt/predictomes/pyproject.toml

# Copy data into package path
COPY predictomes_vis/predictomes_vis/predictome_data/20251113_download_hsbps.csv /opt/predictomes/predictomes_vis/predictome_data/20251113_download_hsbps.csv
COPY predictomes_vis/predictomes_vis/predictome_data/json_metrics_matrix.tsv /opt/predictomes/predictomes_vis/predictome_data/json_metrics_matrix.tsv

ENV PREDICTOMES_PEAK_MATRIX=/opt/predictomes/predictomes_vis/predictome_data/json_metrics_matrix.tsv
ENV PIP_DEFAULT_TIMEOUT=180

RUN pip install --no-cache-dir --progress-bar off \
    pyvis==0.3.2 \
    certifi==2024.8.30 \
    "pandas>=2.0.0" \
    "plotly>=5.0.0" \
    kaleido==0.2.1

ENV PYTHONPATH=/opt/predictomes/predictomes_vis:/opt/predictomes
ENTRYPOINT ["python", "-m", "predictomes_vis.cli"]
