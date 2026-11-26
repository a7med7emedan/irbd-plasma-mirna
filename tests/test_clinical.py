# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from statsmodels.stats.multitest import multipletests

from irbdmirna import clinical


def _data():
    rng = np.random.default_rng(4)
    n = 40
    data = pd.DataFrame(
        {
            "group": ["a"] * n + ["b"] * n,
            "age": np.r_[rng.normal(60, 5, n), rng.normal(66, 9, n)],
            "score": rng.normal(25, 3, 2 * n),
            "cvd": ["yes"] * 12 + ["no"] * 28 + ["yes"] * 3 + ["no"] * 37,
        }
    )
    data.loc[[0, 1, 45], "score"] = np.nan
    return data


def test_welch_fisher_and_benjamini_hochberg():
    data = _data()
    out = clinical.compare_groups(data, "group", ["a"], ["b"], ["age", "score"], ["cvd"])
    a, b = data[data.group == "a"], data[data.group == "b"]
    p_age = stats.ttest_ind(a["age"], b["age"], equal_var=False).pvalue
    p_score = stats.ttest_ind(a["score"].dropna(), b["score"].dropna(), equal_var=False).pvalue
    p_cvd = stats.fisher_exact([[12, 28], [3, 37]]).pvalue
    assert out["p_value"].tolist() == pytest.approx([p_age, p_score, p_cvd])
    assert out["p_adjusted"].tolist() == pytest.approx(
        multipletests([p_age, p_score, p_cvd], method="fdr_bh")[1]
    )
    # pairwise deletion: only the variable with missing values loses individuals
    assert out.set_index("variable").loc["score", ["n_first", "n_second"]].tolist() == [38, 39]
    assert out.set_index("variable").loc["age", ["n_first", "n_second"]].tolist() == [40, 40]
    assert out.set_index("variable").loc["cvd (yes)", "first"] == "12 (30.0%)"


def test_more_than_two_levels_is_rejected():
    data = _data()
    data["stage"] = ["i", "ii", "iii", "i"] * 20
    with pytest.raises(ValueError, match="two levels"):
        clinical.compare_groups(data, "group", ["a"], ["b"], [], ["stage"])


def test_missing_data_table():
    out = clinical.missing_data(_data(), "group", ["age", "score"]).set_index("variable")
    assert out.loc["age", "a"] == "0 (0.0%)"
    assert out.loc["score", "a"] == "2 (5.0%)"
    assert out.loc["score", "b"] == "1 (2.5%)"
