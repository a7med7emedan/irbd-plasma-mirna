# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from statsmodels.stats.multitest import multipletests

from irbdmirna import converters


def test_detection_threshold_follows_the_analysis():
    assert converters.detection_threshold(42, 0.2) == 8  # floor(8.4)
    assert converters.detection_threshold(10, 0.2) == 3  # lower bound of three libraries
    assert converters.detection_threshold(100, 0.2) == 20


def _data(seed=5, n_mirna=60, n_conv=8, n_non=34):
    rng = np.random.default_rng(seed)
    values = rng.gamma(2.0, 50.0, size=(n_mirna, n_conv + n_non))
    values[:10, :] = np.round(values[:10, :] / 40)  # rows with ties and zeros
    labels = np.array([True] * n_conv + [False] * n_non)
    return pd.DataFrame(values, index=[f"m{i}" for i in range(n_mirna)]), labels


def test_row_wise_test_equals_scipy_default_per_mirna():
    table, labels = _data()
    x = table.to_numpy()
    got = converters.mannwhitney_rows(x[:, labels], x[:, ~labels])
    expected = [stats.mannwhitneyu(r[labels], r[~labels]).pvalue for r in x]
    assert got == pytest.approx(expected, rel=1e-12)


def test_compare_matches_manual_calculation():
    table, labels = _data()
    table.iloc[20, np.where(labels)[0]] *= 30.0  # one strongly shifted miRNA
    out = converters.compare(table, labels)
    x = table.to_numpy()
    p = np.array([stats.mannwhitneyu(r[labels], r[~labels]).pvalue for r in x])
    q = multipletests(p, method="fdr_bh")[1]
    lfc = np.log2((x[:, labels].mean(axis=1) + 1) / (x[:, ~labels].mean(axis=1) + 1))
    assert out.loc[table.index, "padj"].to_numpy() == pytest.approx(q)
    assert out.loc[table.index, "log2FoldChange"].to_numpy() == pytest.approx(lfc)
    assert bool(out.loc["m20", "dysregulated"])
    assert int(out["dysregulated"].sum()) == 1


def test_normalise_divides_by_size_factor():
    counts = pd.DataFrame({"a": [10.0, 4.0], "b": [30.0, 6.0]})
    out = converters.normalise(counts, pd.Series({"a": 2.0, "b": 3.0}))
    assert out["a"].tolist() == [5.0, 2.0]
    assert out["b"].tolist() == [10.0, 2.0]
    with pytest.raises(ValueError, match="No size factor"):
        converters.normalise(counts, pd.Series({"a": 2.0}))


def test_permutation_calibration_is_reproducible_and_calibrated():
    table, labels = _data(seed=11)
    first = converters.permutation_calibration(table, labels, n_permutations=40, seed=1)
    again = converters.permutation_calibration(table, labels, n_permutations=40, seed=1)
    assert first == again
    assert first["n_permutations"] == 40
    assert first["fraction_with_calls"] <= 0.15  # null data: few permutations give a call


def test_groups_need_two_libraries():
    table, labels = _data(n_conv=1, n_non=5)
    with pytest.raises(ValueError, match="at least two"):
        converters.compare(table, labels)
