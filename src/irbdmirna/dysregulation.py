# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Dysregulation rule applied to DESeq2 tables."""

from __future__ import annotations

import numpy as np
import pandas as pd


def annotate(table: pd.DataFrame, fdr: float = 0.05, abs_log2fc: float = 0.5) -> pd.DataFrame:
    """Add the columns tested, dysregulated and direction.

    A miRNA is tested when it has an adjusted p value. It is dysregulated when the adjusted
    p value is below `fdr` and the absolute shrunken log2 fold change is at least `abs_log2fc`.
    """
    out = table.copy()
    out["tested"] = out["padj"].notna()
    out["dysregulated"] = (
        out["tested"] & (out["padj"] < fdr) & (out["log2FoldChange"].abs() >= abs_log2fc)
    )
    out["direction"] = np.where(
        out["dysregulated"], np.where(out["log2FoldChange"] > 0, "up", "down"), "none"
    )
    return out


def summarise(table: pd.DataFrame) -> dict[str, float]:
    """Counts for one annotated table."""
    tested = table[table["tested"]]
    return {
        "fitted": int(len(table)),
        "tested": int(len(tested)),
        "dysregulated": int(table["dysregulated"].sum()),
        "up": int((table["direction"] == "up").sum()),
        "down": int((table["direction"] == "down").sum()),
        "smallest_padj": float(tested["padj"].min()) if len(tested) else float("nan"),
    }


def significant(table: pd.DataFrame) -> pd.DataFrame:
    """Rows that meet the dysregulation rule, most significant first."""
    cols = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "direction"]
    return table.loc[table["dysregulated"], cols].sort_values(["padj", "pvalue"])
