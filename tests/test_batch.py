# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from irbdmirna import batch


def test_balance_chi_square_and_cramers_v():
    groups = pd.Series(["a"] * 20 + ["b"] * 20)
    batches = pd.Series((["x"] * 10 + ["y"] * 10) * 2)
    res = batch.balance(groups, batches)
    assert res["p_value"] == pytest.approx(1.0)
    assert res["cramers_v"] == pytest.approx(0.0)
    skewed = pd.Series(["x"] * 18 + ["y"] * 2 + ["x"] * 2 + ["y"] * 18)
    res = batch.balance(groups, skewed)
    chi2 = stats.chi2_contingency(pd.crosstab(groups, skewed))[0]
    assert res["cramers_v"] == pytest.approx(np.sqrt(chi2 / 40))
    assert res["p_value"] < 1e-4
    with pytest.raises(ValueError, match="Only one batch"):
        batch.balance(groups, pd.Series(["x"] * 40))


def _expression(seed=2, n=60, m=80):
    rng = np.random.default_rng(seed)
    libs = [f"L{i}" for i in range(n)]
    groups = pd.Series(np.repeat(["a", "b", "c"], n // 3), index=libs)
    batches = pd.Series(np.tile(["x", "y"], n // 2), index=libs)
    rpm = pd.DataFrame(
        rng.gamma(5.0, 40.0, size=(m, n)), index=[f"m{i}" for i in range(m)], columns=libs
    )
    return rpm, groups, batches


def test_expressed_filter_keeps_detected_mirnas():
    rpm, _, _ = _expression()
    rpm.iloc[0, :40] = 0.0  # detected in one third of the libraries
    rpm.iloc[1, :30] = 0.0  # detected in exactly half
    logged = batch.expressed_log2(rpm, 0.5)
    assert "m0" not in logged.index
    assert "m1" in logged.index
    assert logged.loc["m2", "L0"] == pytest.approx(np.log2(rpm.loc["m2", "L0"] + 1))


def test_per_mirna_test_finds_planted_batch_effect():
    rpm, groups, batches = _expression()
    rpm.loc["m5", (batches == "y").to_numpy()] *= 6.0
    per = batch.per_mirna_batch_test(batch.expressed_log2(rpm, 0.5), groups, batches)
    assert per.index[0] == "m5"
    assert per.loc["m5", "padj"] < 0.001
    assert int((per["padj"] < 0.05).sum()) == 1


def test_principal_components_detect_group_structure():
    rpm, groups, batches = _expression()
    rpm.loc[rpm.index[:30], (groups == "c").to_numpy()] *= 4.0
    scores, comps = batch.principal_components(batch.expressed_log2(rpm, 0.5), groups, batches, 3)
    assert scores.shape == (60, 3)
    assert comps.loc[0, "p_group"] < 1e-6
    assert comps.loc[0, "p_batch"] > 0.05
    assert comps["variance_explained_percent"].is_monotonic_decreasing
