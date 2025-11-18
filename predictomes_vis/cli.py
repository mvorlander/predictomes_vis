#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Alphafold2-based protein-protein interaction prediction visualization toolkit (data from predictomes.org).\n"
            "Main use cases:\n"
            "  network (FDR): build networks from the high-confidence SPOC subset (<10% FDR)\n"
            "  network (peak_score): build networks from the full predictome matrix (~1.6M pairs, thresholdable), using the "
            "Peak score as implemented in https://gitlab.com/BrenneckeLab/ht-colabfold to detect interactions\n"
            "  binary: plot per-POI metrics (ptm/iptm/peak) in the style of VBC HT-colabfold pipeline\n"
            "\nExamples:\n"
            "  predictomes-vis network --metric FDR --poi A1AG2_HUMAN --expansion-depth 0\n"
            "  predictomes-vis network --metric peak_score --poi A1AG2_HUMAN --predictome-peak-threshold 0.8 --expansion-depth 1\n"
            "  predictomes-vis binary --poi A1AG2_HUMAN --outdir outputs/predictome_binary_plots\n"
            "\nRun  \"predictomes-vis network --help\"   or \"predictomes-vis binary --help\")  to see full list of options."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True, metavar="{network,binary}")

    subparsers.add_parser(
        "network",
        help="Build PPI networks (FDR/SPOC or predictome peak_score).",
        add_help=False,
    )

    subparsers.add_parser(
        "binary",
        help="Plot per-POI metrics from predictome matrix.",
        add_help=False,
    )

    args, remainder = parser.parse_known_args()

    if args.command == "network":
        try:
            from predictomes_build_ppi_network import cli as network_cli  # type: ignore
        except ModuleNotFoundError as exc:
            if exc.name == "predictomes_build_ppi_network":
                import pathlib

                package_root = pathlib.Path(__file__).resolve().parents[1]
                sys.path.append(str(package_root))
                from predictomes_build_ppi_network import cli as network_cli  # type: ignore
            else:
                print(f"Missing dependency for network mode: {exc.name}. Try pip install certifi pyvis.", file=sys.stderr)
                sys.exit(1)

        if not remainder or any(flag in remainder for flag in ("-h", "--help")):
            remainder = ["--help"]
        sys.argv = ["predictomes-vis network"] + remainder
        network_cli.main()
    else:
        try:
            from predictomes_vis import plot_cli
        except ModuleNotFoundError as exc:
            print(f"Missing dependency for binary mode: {exc.name}. Try pip install pandas plotly kaleido.", file=sys.stderr)
            sys.exit(1)

        if not remainder or any(flag in remainder for flag in ("-h", "--help")):
            remainder = ["--help"]
        sys.argv = ["predictomes-vis binary"] + remainder
        plot_cli.main()


if __name__ == "__main__":
    main()
