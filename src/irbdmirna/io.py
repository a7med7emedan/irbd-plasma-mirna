# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Input tables: reading, validation and basic count summaries."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED_SAMPLE_COLUMNS = ("library_id", "group", "age", "sex", "batch", "qc_pass")
_TRUE = {"true", "1", "yes", "y", "ok", "pass"}
_FALSE = {"false", "0", "no", "n", "fail", ""}


class InputError(ValueError):
    """Raised when an input table does not have the expected shape or content."""


def _to_bool(series: pd.Series, column: str) -> pd.Series:
    text = series.astype(str).str.strip().str.lower().replace({"nan": ""})
    unknown = sorted(set(text) - _TRUE - _FALSE)
    if unknown:
        raise InputError(
            f"Column '{column}' has values that are not true or false: {unknown}. "
            "Use true or false."
        )
    return text.isin(_TRUE)


def read_counts(path: str | Path) -> pd.DataFrame:
    """Read a tab-separated miRNA by library count matrix. The first column holds miRNA names."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"Count table not found: {path}")
    counts = pd.read_csv(path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    counts.index.name = "miRNA"
    if counts.index.has_duplicates:
        dup = counts.index[counts.index.duplicated()].unique()[:5].tolist()
        raise InputError(f"The count table has duplicated miRNA names, for example {dup}.")
    if counts.columns.has_duplicates:
        raise InputError("The count table has duplicated library columns.")
    non_numeric = [c for c in counts.columns if not pd.api.types.is_numeric_dtype(counts[c])]
    if non_numeric:
        raise InputError(f"These count columns are not numeric: {non_numeric[:5]}.")
    if counts.isna().any().any():
        raise InputError("The count table has missing values. Counts must be complete.")
    if (counts < 0).any().any():
        raise InputError("The count table has negative values.")
    return counts.astype(float)


def read_samples(path: str | Path, group_labels: list[str]) -> pd.DataFrame:
    """Read the tab-separated sample sheet and validate its columns and values."""
    path = Path(path)
    if not path.is_file():
        raise InputError(f"Sample sheet not found: {path}")
    samples = pd.read_csv(path, sep="\t", dtype={"library_id": str, "group": str, "sex": str})
    missing = [c for c in REQUIRED_SAMPLE_COLUMNS if c not in samples.columns]
    if missing:
        raise InputError(
            f"The sample sheet lacks these columns: {missing}. "
            f"Required columns: {list(REQUIRED_SAMPLE_COLUMNS)}."
        )
    if samples["library_id"].duplicated().any():
        raise InputError("The sample sheet lists a library_id more than once.")
    unknown = sorted(set(samples["group"]) - set(group_labels))
    if unknown:
        raise InputError(
            f"The sample sheet has unknown group labels: {unknown}. Expected: {group_labels}."
        )
    samples["qc_pass"] = _to_bool(samples["qc_pass"], "qc_pass")
    if "converter" in samples.columns:
        samples["converter"] = _to_bool(samples["converter"], "converter")
    else:
        samples["converter"] = False
    samples["age"] = pd.to_numeric(samples["age"], errors="coerce")
    samples["batch"] = samples["batch"].astype(str)
    return samples


def analysis_set(samples: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    """Return the libraries that passed quality control, in sample-sheet order."""
    keep = samples[samples["qc_pass"]].reset_index(drop=True)
    if keep.empty:
        raise InputError("No library has qc_pass set to true.")
    absent = [lib for lib in keep["library_id"] if lib not in counts.columns]
    if absent:
        raise InputError(
            f"{len(absent)} libraries of the sample sheet are not in the count table, "
            f"for example {absent[:5]}."
        )
    for column in ("age", "sex"):
        if keep[column].isna().any():
            raise InputError(f"Column '{column}' has missing values among the analysed libraries.")
    return keep


def round_counts(counts: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """Round counts to integers. Return the integer table and the share of non-integer values."""
    values = counts.to_numpy(dtype=float)
    non_integer = float(np.mean(values != np.floor(values)))
    rounded = pd.DataFrame(
        np.rint(values).astype(np.int64), index=counts.index, columns=counts.columns
    )
    return rounded, non_integer


def detection_summary(counts: pd.DataFrame) -> dict[str, int]:
    """Count annotated miRNAs and those detected in at least one or at least half of the libraries."""
    detected = counts > 0
    return {
        "annotated": int(counts.shape[0]),
        "detected_in_at_least_one": int(detected.any(axis=1).sum()),
        "detected_in_at_least_half": int((detected.mean(axis=1) >= 0.5).sum()),
    }


def rpmmm(counts: pd.DataFrame) -> pd.DataFrame:
    """Reads per million mapped miRNA reads, per library."""
    totals = counts.sum(axis=0)
    if (totals <= 0).any():
        raise InputError("At least one library has no miRNA reads.")
    return counts.div(totals, axis=1) * 1e6
