# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Command-line interface."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__, cluster14q32, pipeline, synthetic
from .config import ConfigError, load_config
from .de import DeseqError
from .io import InputError


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="irbd-mirna",
        description="Plasma microRNA analysis for iRBD, de novo PD and healthy controls.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show debug messages.")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Run the whole analysis from one configuration file.")
    run.add_argument("--config", required=True, help="Path to the YAML configuration file.")
    run.add_argument("--output-dir", help="Override the output directory of the configuration.")
    run.add_argument(
        "--skip-de",
        action="store_true",
        help="Reuse the DESeq2 tables already in the output directory.",
    )
    run.add_argument(
        "--strict",
        action="store_true",
        help="Fail when a key number differs from the manuscript values.",
    )

    example = sub.add_parser("make-example-data", help="Write the synthetic example dataset.")
    example.add_argument("--outdir", default="examples/data", help="Destination directory.")
    example.add_argument("--seed", type=int, default=20261126, help="Random seed.")

    res = sub.add_parser(
        "build-cluster-resource", help="Derive the 14q32 domain table from a miRBase GFF3 file."
    )
    res.add_argument("--gff", required=True, help="miRBase hsa.gff3 (GRCh38).")
    res.add_argument("--out", required=True, help="Output TSV path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%H:%M:%S",
    )
    try:
        if args.command == "run":
            cfg = load_config(args.config)
            if args.output_dir:
                cfg["output_dir"] = args.output_dir
            pipeline.run(cfg, skip_de=args.skip_de, strict=args.strict)
        elif args.command == "make-example-data":
            paths = synthetic.write_dataset(args.outdir, args.seed)
            for name, path in paths.items():
                print(f"{name}: {path}")
        elif args.command == "build-cluster-resource":
            table = cluster14q32.extract_domain(args.gff)
            Path(args.out).parent.mkdir(parents=True, exist_ok=True)
            table.to_csv(args.out, sep="\t", index=False)
            print(f"{len(table)} mature miRNA rows written to {args.out}")
    except (ConfigError, InputError, DeseqError, pipeline.ManuscriptMismatch, ValueError) as err:
        logging.getLogger("irbdmirna").error("%s", err)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
