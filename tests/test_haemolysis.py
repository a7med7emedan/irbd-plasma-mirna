# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from irbdmirna import haemolysis


def test_haemolysis_index_formula():
    rpm = pd.DataFrame(
        {"s1": [7.0, 1.0], "s2": [0.0, 3.0]}, index=["hsa-miR-451a", "hsa-miR-23a-3p"]
    )
    index = haemolysis.haemolysis_index(rpm, "hsa-miR-451a", "hsa-miR-23a-3p")
    assert index["s1"] == pytest.approx(np.log2(8 / 2))
    assert index["s2"] == pytest.approx(np.log2(1 / 4))
    with pytest.raises(ValueError, match="hsa-miR-16-5p"):
        haemolysis.haemolysis_index(rpm, "hsa-miR-16-5p", "hsa-miR-23a-3p")


def test_cluster_score_is_mean_z_score():
    rng = np.random.default_rng(1)
    rpm = pd.DataFrame(rng.gamma(2.0, 50.0, size=(5, 12)), index=[f"m{i}" for i in range(5)])
    score = haemolysis.cluster_score(rpm, ["m0", "m1", "m2"])
    logged = np.log2(rpm.loc[["m0", "m1", "m2"]] + 1)
    z = (logged.to_numpy() - logged.to_numpy().mean(axis=1, keepdims=True)) / logged.to_numpy().std(
        axis=1, ddof=1, keepdims=True
    )
    assert score.to_numpy() == pytest.approx(z.mean(axis=0), abs=1e-6)
    assert score.mean() == pytest.approx(0.0, abs=1e-9)
    with pytest.raises(ValueError, match="No cluster miRNA"):
        haemolysis.cluster_score(rpm, [])


def test_pairwise_mannwhitney_matches_scipy():
    values = pd.Series([1.0, 2.0, 3.0, 4.0, 10.0, 11.0, 12.0, 13.0, 5.0])
    groups = pd.Series(["a"] * 4 + ["b"] * 4 + ["c"])
    out = haemolysis.pairwise_mannwhitney(values, groups, [("a", "b")])
    expected = stats.mannwhitneyu(values[:4], values[4:8], alternative="two-sided").pvalue
    assert out.loc[0, "p_value"] == pytest.approx(expected)
    assert (out.loc[0, "n_first"], out.loc[0, "n_second"]) == (4, 4)


def test_score_models_recover_planted_effects():
    rng = np.random.default_rng(3)
    n = 120
    group = np.repeat(["controls", "de novo PD", "iRBD"], n // 3)
    haem = rng.normal(0, 1, n) + 1.0 * (group == "iRBD")
    score = -1.0 * (group == "iRBD") - 0.5 * haem + rng.normal(0, 0.05, n)
    data = pd.DataFrame(
        {
            "cluster_score": score,
            "group3": group,
            "haemolysis_index": haem,
            "age": rng.normal(65, 8, n),
            "sex": rng.choice(["male", "female"], n),
        }
    )
    models = haemolysis.score_models(data, reference="controls", target="iRBD")
    assert models["model"].tolist()[0] == "group only"
    assert models.loc[0, "coefficient"] == pytest.approx(-1.5, abs=0.2)
    assert models.loc[1, "coefficient"] == pytest.approx(-1.0, abs=0.05)
    assert models.loc[0, "percent_reduction_vs_group_only"] == 0.0
    assert models.loc[1, "percent_reduction_vs_group_only"] == pytest.approx(33.0, abs=6.0)
    assert models.loc[2, "p_value"] < 1e-10


def test_spearman_by_group():
    x = pd.Series([1, 2, 3, 4, 1, 2, 3, 4], dtype=float)
    y = pd.Series([1, 2, 3, 4, 4, 3, 2, 1], dtype=float)
    out = haemolysis.spearman_by_group(x, y, pd.Series(["a"] * 4 + ["b"] * 4)).set_index("group")
    assert out.loc["a", "spearman_rho"] == pytest.approx(1.0)
    assert out.loc["b", "spearman_rho"] == pytest.approx(-1.0)
