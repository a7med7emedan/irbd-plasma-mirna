# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Haemolysis index, 14q32 cluster score and their relation to study group."""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats


def haemolysis_index(rpm: pd.DataFrame, numerator: str, denominator: str) -> pd.Series:
    """log2((numerator + 1) / (denominator + 1)) on reads per million mapped miRNA reads."""
    for name in (numerator, denominator):
        if name not in rpm.index:
            raise ValueError(f"{name} is not in the count table. The haemolysis index needs it.")
    index = np.log2((rpm.loc[numerator] + 1.0) / (rpm.loc[denominator] + 1.0))
    index.name = "haemolysis_index"
    return index


def cluster_score(rpm: pd.DataFrame, members: list[str]) -> pd.Series:
    """Mean over cluster miRNAs of the z-scored log2(RPMMM + 1), per library."""
    if not members:
        raise ValueError("No cluster miRNA was tested in iRBD versus controls. No score.")
    logged = np.log2(rpm.loc[members] + 1.0)
    centred = logged.sub(logged.mean(axis=1), axis=0)
    z = centred.div(logged.std(axis=1, ddof=1) + 1e-9, axis=0)
    score = z.mean(axis=0)
    score.name = "cluster_score"
    return score


def pairwise_mannwhitney(values: pd.Series, groups: pd.Series, pairs) -> pd.DataFrame:
    """Two-sided Mann-Whitney U test for each (first, second) pair of group labels."""
    rows = []
    for first, second in pairs:
        x = values[groups == first]
        y = values[groups == second]
        res = stats.mannwhitneyu(x, y, alternative="two-sided")
        rows.append(
            {
                "comparison": f"{first} vs {second}",
                "n_first": int(len(x)),
                "n_second": int(len(y)),
                "median_first": float(x.median()),
                "median_second": float(y.median()),
                "U": float(res.statistic),
                "p_value": float(res.pvalue),
            }
        )
    return pd.DataFrame(rows)


def spearman_by_group(x: pd.Series, y: pd.Series, groups: pd.Series) -> pd.DataFrame:
    """Spearman correlation of x and y within each group."""
    rows = []
    for g in pd.unique(groups):
        mask = groups == g
        rho, p = stats.spearmanr(x[mask], y[mask])
        rows.append(
            {"group": g, "n": int(mask.sum()), "spearman_rho": float(rho), "p_value": float(p)}
        )
    return pd.DataFrame(rows)


def score_models(data: pd.DataFrame, reference: str, target: str) -> pd.DataFrame:
    """Linear regressions of the cluster score on study group, then with adjustment.

    `data` needs the columns cluster_score, group3, haemolysis_index, age and sex.
    The reported coefficient is `target` relative to `reference`.
    """
    d = data.copy()
    levels = [reference] + sorted(set(d["group3"]) - {reference})
    d["group3"] = pd.Categorical(d["group3"], categories=levels)
    term = f"C(group3)[T.{target}]"
    formulas = [
        ("group only", "cluster_score ~ C(group3)"),
        ("plus haemolysis index", "cluster_score ~ C(group3) + haemolysis_index"),
        (
            "plus haemolysis index, age and sex",
            "cluster_score ~ C(group3) + haemolysis_index + age + C(sex)",
        ),
    ]
    rows = []
    base = None
    for label, formula in formulas:
        fit = smf.ols(formula, d).fit()
        coef = float(fit.params[term])
        base = coef if base is None else base
        rows.append(
            {
                "model": label,
                "formula": formula,
                "n": int(fit.nobs),
                "coefficient": coef,
                "std_error": float(fit.bse[term]),
                "p_value": float(fit.pvalues[term]),
                "percent_reduction_vs_group_only": 100.0 * (1.0 - coef / base),
            }
        )
    return pd.DataFrame(rows)
