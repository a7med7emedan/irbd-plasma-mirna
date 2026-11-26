# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Predefined brain-enriched miRNA panel: counts, overlap and over-representation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy import stats


def read_panel(path: str | Path) -> list[str]:
    """Read the panel table. The column miRNA holds the mature miRNA names."""
    panel = pd.read_csv(path, sep="\t")
    if "miRNA" not in panel.columns:
        raise ValueError(f"The panel table {path} needs a column named 'miRNA'.")
    names = panel["miRNA"].astype(str).tolist()
    if len(set(names)) != len(names):
        raise ValueError(f"The panel table {path} lists a miRNA more than once.")
    return names


def panel_status(table: pd.DataFrame, panel: list[str]) -> pd.DataFrame:
    """Per panel miRNA: whether it was tested, its estimates and its call."""
    rows = []
    for name in panel:
        if name in table.index:
            r = table.loc[name]
            status = "tested" if r["tested"] else "not tested (independent filtering)"
            rows.append(
                {
                    "miRNA": name,
                    "status": status,
                    "log2FoldChange": r["log2FoldChange"],
                    "padj": r["padj"],
                    "dysregulated": bool(r["dysregulated"]),
                    "direction": r["direction"],
                }
            )
        else:
            rows.append(
                {
                    "miRNA": name,
                    "status": "not tested (low counts)",
                    "log2FoldChange": float("nan"),
                    "padj": float("nan"),
                    "dysregulated": False,
                    "direction": "none",
                }
            )
    return pd.DataFrame(rows).set_index("miRNA")


def over_representation(
    table: pd.DataFrame, panel: list[str], exclude: set[str] | None = None
) -> dict[str, float]:
    """One-sided Fisher's exact test for panel miRNAs among dysregulated miRNAs.

    The background is every miRNA with an adjusted p value in the comparison. Names in
    `exclude` are removed from the panel and from the background before testing.
    """
    tested = table[table["tested"]]
    if exclude:
        tested = tested[~tested.index.isin(exclude)]
    in_panel = tested.index.isin(panel)
    dys = tested["dysregulated"].to_numpy()
    a = int((in_panel & dys).sum())
    b = int((in_panel & ~dys).sum())
    c = int((~in_panel & dys).sum())
    d = int((~in_panel & ~dys).sum())
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "panel_tested": a + b,
        "panel_dysregulated": a,
        "background_tested": c + d,
        "background_dysregulated": c,
        "odds_ratio": float(odds),
        "p_value": float(p),
    }


def shared_candidates(first: pd.DataFrame, second: pd.DataFrame) -> pd.DataFrame:
    """Panel miRNAs dysregulated in both comparisons, with direction concordance."""
    both = first.index[first["dysregulated"] & second.loc[first.index, "dysregulated"]]
    out = pd.DataFrame(
        {
            "direction_first": first.loc[both, "direction"],
            "direction_second": second.loc[both, "direction"],
        }
    )
    out["concordant"] = out["direction_first"] == out["direction_second"]
    return out
