# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Batch-effect analysis: balance across groups, principal components, per-miRNA F test."""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests


def balance(groups: pd.Series, batch: pd.Series) -> dict:
    """Chi-square test of batch by group, with Cramer's V."""
    table = pd.crosstab(groups, batch)
    if table.shape[1] < 2:
        raise ValueError("Only one batch is present. The batch analysis needs at least two.")
    chi2, p, _, _ = stats.chi2_contingency(table)
    n = table.to_numpy().sum()
    v = float(np.sqrt(chi2 / (n * (min(table.shape) - 1))))
    return {"table": table, "chi2": float(chi2), "p_value": float(p), "cramers_v": v}


def expressed_log2(rpm: pd.DataFrame, min_detection_fraction: float) -> pd.DataFrame:
    """log2(RPMMM + 1) of the miRNAs detected in at least the given share of libraries."""
    keep = (rpm > 0).mean(axis=1) >= min_detection_fraction
    return np.log2(rpm.loc[keep] + 1.0)


def principal_components(
    logged: pd.DataFrame, groups: pd.Series, batch: pd.Series, n_components: int = 3
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """PCA on standardised miRNAs. Kruskal-Wallis test of each component by batch and group."""
    x = logged.T.to_numpy()
    x = (x - x.mean(axis=0)) / (x.std(axis=0) + 1e-9)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    k = min(n_components, len(s))
    explained = s**2 / (s**2).sum() * 100.0
    scores = pd.DataFrame(
        u[:, :k] * s[:k], index=logged.columns, columns=[f"PC{i + 1}" for i in range(k)]
    )
    rows = []
    for i in range(k):
        pc = scores[f"PC{i + 1}"]
        by_batch = [pc[(batch == b).to_numpy()] for b in pd.unique(batch)]
        by_group = [pc[(groups == g).to_numpy()] for g in pd.unique(groups)]
        rows.append(
            {
                "component": f"PC{i + 1}",
                "variance_explained_percent": float(explained[i]),
                "p_batch": float(stats.kruskal(*by_batch).pvalue),
                "p_group": float(stats.kruskal(*by_group).pvalue),
            }
        )
    return scores, pd.DataFrame(rows)


def per_mirna_batch_test(logged: pd.DataFrame, groups: pd.Series, batch: pd.Series) -> pd.DataFrame:
    """F test of a model with group and batch against a model with group alone, per miRNA."""
    design = pd.DataFrame({"group": groups.to_numpy(), "batch": batch.to_numpy()})
    rows = []
    for name, values in logged.iterrows():
        d = design.assign(y=values.to_numpy())
        if d["y"].std() == 0:
            continue
        small = smf.ols("y ~ C(group)", d).fit()
        large = smf.ols("y ~ C(group) + C(batch)", d).fit()
        rows.append((name, float(anova_lm(small, large)["Pr(>F)"].iloc[1])))
    out = pd.DataFrame(rows, columns=["miRNA", "p_value"]).dropna().set_index("miRNA")
    out["padj"] = multipletests(out["p_value"], method="fdr_bh")[1] if len(out) else []
    return out.sort_values("p_value")
