# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Converters versus non-converters within iRBD: rank test and permutation calibration."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def detection_threshold(n_samples: int, fraction: float, minimum: int = 3) -> int:
    """Smallest number of libraries with a non-zero count for a miRNA to be tested.

    The threshold is floor(fraction x n), with a lower bound of `minimum` libraries.
    For 42 libraries and a fraction of 0.2 this gives 8.
    """
    return max(int(minimum), math.floor(fraction * n_samples))


def normalise(counts: pd.DataFrame, size_factors: pd.Series) -> pd.DataFrame:
    """Divide each library by its DESeq2 size factor."""
    missing = [c for c in counts.columns if c not in size_factors.index]
    if missing:
        raise ValueError(f"No size factor for libraries: {missing[:5]}.")
    return counts.div(size_factors.loc[counts.columns], axis=1)


def mannwhitney_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Two-sided Mann-Whitney U p value for each row.

    The method follows the SciPy default applied to each miRNA on its own: the exact
    distribution when one group has at most 8 libraries and the row has no ties,
    the normal approximation with tie and continuity correction otherwise.
    """
    n1, n2 = a.shape[1], b.shape[1]
    if n1 > 8 and n2 > 8:
        res = stats.mannwhitneyu(a, b, axis=1, alternative="two-sided", method="asymptotic")
        return np.asarray(res.pvalue, dtype=float)
    ordered = np.sort(np.concatenate([a, b], axis=1), axis=1)
    ties = (np.diff(ordered, axis=1) == 0).any(axis=1)
    p = np.empty(a.shape[0], dtype=float)
    for mask, method in ((ties, "asymptotic"), (~ties, "exact")):
        if mask.any():
            res = stats.mannwhitneyu(
                a[mask], b[mask], axis=1, alternative="two-sided", method=method
            )
            p[mask] = res.pvalue
    return p


def _test(values: np.ndarray, is_converter: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    a = values[:, is_converter]
    b = values[:, ~is_converter]
    p = mannwhitney_rows(a, b)
    lfc = np.log2((a.mean(axis=1) + 1.0) / (b.mean(axis=1) + 1.0))
    return p, lfc


def compare(
    normalised: pd.DataFrame, is_converter: np.ndarray, fdr: float = 0.05, abs_log2fc: float = 0.5
) -> pd.DataFrame:
    """Two-sided Mann-Whitney U test per miRNA with Benjamini-Hochberg adjustment."""
    is_converter = np.asarray(is_converter, dtype=bool)
    if is_converter.sum() < 2 or (~is_converter).sum() < 2:
        raise ValueError("The converter comparison needs at least two libraries per group.")
    p, lfc = _test(normalised.to_numpy(dtype=float), is_converter)
    out = pd.DataFrame(
        {
            "baseMean": normalised.mean(axis=1).to_numpy(),
            "log2FoldChange": lfc,
            "pvalue": p,
            "padj": multipletests(p, method="fdr_bh")[1],
        },
        index=normalised.index,
    )
    out["dysregulated"] = (out["padj"] < fdr) & (out["log2FoldChange"].abs() >= abs_log2fc)
    return out.sort_values(["pvalue", "padj"])


def permutation_calibration(
    normalised: pd.DataFrame,
    is_converter: np.ndarray,
    n_permutations: int,
    seed: int,
    fdr: float = 0.05,
    abs_log2fc: float = 0.5,
) -> dict[str, float]:
    """Repeat the full procedure on datasets with randomly permuted group labels."""
    is_converter = np.asarray(is_converter, dtype=bool)
    values = normalised.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    with_calls = 0
    for _ in range(n_permutations):
        labels = rng.permutation(is_converter)
        p, lfc = _test(values, labels)
        q = multipletests(p, method="fdr_bh")[1]
        if np.any((q < fdr) & (np.abs(lfc) >= abs_log2fc)):
            with_calls += 1
    return {
        "n_permutations": int(n_permutations),
        "permutations_with_calls": int(with_calls),
        "fraction_with_calls": with_calls / n_permutations if n_permutations else float("nan"),
    }
