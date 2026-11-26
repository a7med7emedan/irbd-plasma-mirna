# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Differential expression: comparison definitions and the DESeq2 runner."""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

log = logging.getLogger(__name__)

R_DIR = Path(__file__).resolve().parent / "r"


@dataclass(frozen=True)
class Comparison:
    """A pairwise comparison. Fold changes are numerator relative to reference."""

    name: str
    numerator: tuple[str, ...]
    reference: tuple[str, ...]
    title: str


COMPARISONS: tuple[Comparison, ...] = (
    Comparison(
        "hc_high_vs_low_prs", ("hc_high_prs",), ("hc_low_prs",), "High-PRS vs low-PRS controls"
    ),
    Comparison("irbd_vs_controls", ("irbd",), ("hc_high_prs", "hc_low_prs"), "iRBD vs controls"),
    Comparison("irbd_vs_pd", ("irbd",), ("pd",), "iRBD vs de novo PD"),
    Comparison("pd_vs_controls", ("pd",), ("hc_high_prs", "hc_low_prs"), "de novo PD vs controls"),
)
MAIN_COMPARISONS = ("irbd_vs_controls", "irbd_vs_pd", "pd_vs_controls")


class DeseqError(RuntimeError):
    """Raised when the R step fails."""


def group_labels(cfg: dict[str, Any], keys: tuple[str, ...]) -> list[str]:
    """Translate internal group keys into the labels used in the sample sheet."""
    return [cfg["groups"][k] for k in keys]


def build_job(cfg: dict[str, Any], counts_path: Path, samples_path: Path, outdir: Path) -> dict:
    """Assemble the job description passed to the R script."""
    return {
        "counts": str(counts_path),
        "samples": str(samples_path),
        "outdir": str(outdir),
        "alpha": float(cfg["deseq2"]["independent_filtering_alpha"]),
        "covariates": list(cfg["deseq2"]["covariates"]),
        "min_count": int(cfg["deseq2"]["min_count"]),
        "min_samples": int(cfg["deseq2"]["min_samples"]),
        "within_group": cfg["groups"]["irbd"],
        "comparisons": [
            {
                "name": c.name,
                "numerator": group_labels(cfg, c.numerator),
                "reference": group_labels(cfg, c.reference),
                "prefilter": cfg["deseq2"]["prefilter"][c.name],
            }
            for c in COMPARISONS
        ],
    }


def run_deseq2(
    cfg: dict[str, Any],
    counts: pd.DataFrame,
    samples: pd.DataFrame,
    de_dir: Path,
    work_dir: Path,
) -> None:
    """Write the analysis inputs, call Rscript and check that every table was produced."""
    rscript = shutil.which(cfg["rscript"])
    if rscript is None:
        raise DeseqError(
            f"'{cfg['rscript']}' was not found. Install R with DESeq2, edgeR and apeglm, "
            "or run the pipeline in the Docker image."
        )
    de_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    counts_path = work_dir / "counts_analysis_set.tsv"
    samples_path = work_dir / "samples_analysis_set.tsv"
    counts.loc[:, samples["library_id"]].to_csv(counts_path, sep="\t")
    samples.to_csv(samples_path, sep="\t", index=False)
    job_path = work_dir / "deseq2_job.yaml"
    with job_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(build_job(cfg, counts_path, samples_path, de_dir), handle, sort_keys=False)
    cmd = [rscript, str(R_DIR / "run_deseq2.R"), str(job_path)]
    log.info("Running DESeq2: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    for line in proc.stderr.splitlines():
        log.info("R | %s", line)
    if proc.returncode != 0:
        raise DeseqError(
            "The DESeq2 step failed. The R messages are in the log above. "
            "Check that group labels and covariates in the sample sheet are complete."
        )
    missing = [c.name for c in COMPARISONS if not (de_dir / f"{c.name}_all.tsv").is_file()]
    if missing:
        raise DeseqError(f"DESeq2 finished without writing tables for: {missing}.")


def read_de_tables(outdir: Path) -> dict[str, pd.DataFrame]:
    """Read the per-comparison DESeq2 tables, indexed by miRNA."""
    tables = {}
    for c in COMPARISONS:
        path = outdir / f"{c.name}_all.tsv"
        if not path.is_file():
            raise DeseqError(f"Missing DESeq2 table: {path}. Run the 'de' step first.")
        tables[c.name] = pd.read_csv(path, sep="\t").set_index("miRNA")
    return tables
