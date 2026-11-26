# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Clinical group comparisons and the missing-data overview."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _describe_continuous(x: pd.Series) -> str:
    x = x.dropna()
    return f"{x.mean():.2f} ± {x.std(ddof=1):.2f}" if len(x) > 1 else "na"


def _describe_binary(x: pd.Series, positive) -> str:
    x = x.dropna()
    return f"{int((x == positive).sum())} ({100 * (x == positive).mean():.1f}%)" if len(x) else "na"


def compare_groups(
    data: pd.DataFrame,
    group_column: str,
    first: list[str],
    second: list[str],
    continuous: list[str],
    categorical: list[str],
) -> pd.DataFrame:
    """Welch's t test for continuous and Fisher's exact test for binary variables.

    Missing values are removed per variable (pairwise deletion). P values are adjusted with
    the Benjamini-Hochberg procedure across all variables of this comparison.
    """
    a = data[data[group_column].isin(first)]
    b = data[data[group_column].isin(second)]
    rows = []
    for var in continuous:
        x, y = a[var].dropna().astype(float), b[var].dropna().astype(float)
        p = np.nan
        if len(x) > 1 and len(y) > 1:
            p = float(stats.ttest_ind(x, y, equal_var=False).pvalue)
        rows.append(
            {
                "variable": var,
                "test": "Welch t",
                "n_first": len(x),
                "n_second": len(y),
                "first": _describe_continuous(x),
                "second": _describe_continuous(y),
                "p_value": p,
            }
        )
    for var in categorical:
        levels = sorted(pd.unique(data[var].dropna()).tolist(), key=str)
        if len(levels) > 2:
            raise ValueError(
                f"Categorical variable '{var}' has {len(levels)} levels. "
                "Fisher's exact test is implemented for two levels. Recode the variable."
            )
        x, y = a[var].dropna(), b[var].dropna()
        p = np.nan
        positive = levels[-1] if levels else None
        if len(levels) == 2 and len(x) and len(y):
            table = [
                [int((x == positive).sum()), int((x != positive).sum())],
                [int((y == positive).sum()), int((y != positive).sum())],
            ]
            p = float(stats.fisher_exact(table).pvalue)
        rows.append(
            {
                "variable": f"{var} ({positive})",
                "test": "Fisher exact",
                "n_first": len(x),
                "n_second": len(y),
                "first": _describe_binary(x, positive),
                "second": _describe_binary(y, positive),
                "p_value": p,
            }
        )
    out = pd.DataFrame(rows)
    out["p_adjusted"] = np.nan
    ok = out["p_value"].notna()
    if ok.any():
        out.loc[ok, "p_adjusted"] = multipletests(out.loc[ok, "p_value"], method="fdr_bh")[1]
    return out


def missing_data(data: pd.DataFrame, group_column: str, variables: list[str]) -> pd.DataFrame:
    """Number and percentage of missing values per variable and group."""
    rows = []
    for var in variables:
        row = {"variable": var}
        for g, sub in data.groupby(group_column, sort=True):
            n = int(sub[var].isna().sum())
            row[str(g)] = f"{n} ({100 * n / len(sub):.1f}%)"
        rows.append(row)
    return pd.DataFrame(rows)
