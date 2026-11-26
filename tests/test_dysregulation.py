# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import numpy as np
import pandas as pd

from irbdmirna import dysregulation


def _table():
    return pd.DataFrame(
        {
            "baseMean": 10.0,
            "log2FoldChange": [0.5, -0.5, 0.49, 2.0, -1.0, 3.0],
            "lfcSE": 0.1,
            "stat": 1.0,
            "pvalue": 0.001,
            "padj": [0.049, 0.01, 0.001, 0.05, 0.2, np.nan],
        },
        index=["at_cut_up", "at_cut_down", "small_effect", "padj_at_cut", "not_sig", "untested"],
    )


def test_rule_boundaries():
    out = dysregulation.annotate(_table(), fdr=0.05, abs_log2fc=0.5)
    assert out["dysregulated"].to_dict() == {
        "at_cut_up": True,  # |log2FC| of exactly 0.5 counts
        "at_cut_down": True,
        "small_effect": False,
        "padj_at_cut": False,  # adjusted p must be below 0.05
        "not_sig": False,
        "untested": False,
    }
    assert out.loc["at_cut_up", "direction"] == "up"
    assert out.loc["at_cut_down", "direction"] == "down"
    assert not out.loc["untested", "tested"]


def test_summary_and_significant_table():
    out = dysregulation.annotate(_table())
    s = dysregulation.summarise(out)
    assert (s["fitted"], s["tested"], s["dysregulated"], s["up"], s["down"]) == (6, 5, 2, 1, 1)
    assert s["smallest_padj"] == 0.001
    sig = dysregulation.significant(out)
    assert sig.index.tolist() == ["at_cut_down", "at_cut_up"]
