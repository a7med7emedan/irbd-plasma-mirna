# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import numpy as np
import pandas as pd
import pytest

from irbdmirna import io

GROUPS = ["iRBD", "de_novo_PD", "HC_high_PRS", "HC_low_PRS"]


def test_round_counts_reports_share_of_non_integers():
    counts = pd.DataFrame({"a": [1.0, 2.5, 0.0, 7.0], "b": [3.5, 4.0, 1.0, 2.0]})
    rounded, share = io.round_counts(counts)
    assert share == pytest.approx(2 / 8)
    assert rounded["a"].tolist() == [1, 2, 0, 7]  # 2.5 rounds half to even, as in R
    assert rounded["b"].tolist() == [4, 4, 1, 2]
    assert rounded.dtypes.eq(np.int64).all()


def test_detection_summary():
    counts = pd.DataFrame(
        {"a": [0, 5, 1, 0], "b": [0, 2, 0, 0], "c": [0, 1, 0, 0], "d": [0, 3, 2, 1]},
        index=["m0", "m1", "m2", "m3"],
    )
    assert io.detection_summary(counts) == {
        "annotated": 4,
        "detected_in_at_least_one": 3,
        "detected_in_at_least_half": 2,
    }


def test_rpmmm_sums_to_one_million():
    counts = pd.DataFrame({"a": [10.0, 30.0], "b": [5.0, 5.0]})
    rpm = io.rpmmm(counts)
    assert rpm.sum().tolist() == pytest.approx([1e6, 1e6])
    assert rpm.loc[0, "a"] == pytest.approx(250000.0)


def test_read_samples_and_analysis_set(tmp_path, dataset):
    path = tmp_path / "samples.tsv"
    dataset["samples"].to_csv(path, sep="\t", index=False)
    samples = io.read_samples(path, GROUPS)
    keep = io.analysis_set(samples, dataset["counts"])
    assert keep["qc_pass"].all()
    assert len(keep) == int(dataset["samples"]["qc_pass"].sum())
    assert samples["converter"].dtype == bool


def test_unknown_group_label(tmp_path, dataset):
    bad = dataset["samples"].copy()
    bad.loc[0, "group"] = "PSP"
    path = tmp_path / "samples.tsv"
    bad.to_csv(path, sep="\t", index=False)
    with pytest.raises(io.InputError, match="unknown group labels"):
        io.read_samples(path, GROUPS)


def test_missing_column(tmp_path, dataset):
    path = tmp_path / "samples.tsv"
    dataset["samples"].drop(columns="batch").to_csv(path, sep="\t", index=False)
    with pytest.raises(io.InputError, match="batch"):
        io.read_samples(path, GROUPS)


def test_library_missing_from_counts(dataset):
    samples = dataset["samples"].copy()
    counts = dataset["counts"].drop(columns=samples.loc[0, "library_id"])
    with pytest.raises(io.InputError, match="not in the count table"):
        io.analysis_set(samples, counts)


def test_counts_validation(tmp_path):
    path = tmp_path / "counts.tsv"
    pd.DataFrame({"a": [1, -2]}, index=["m1", "m2"]).to_csv(path, sep="\t")
    with pytest.raises(io.InputError, match="negative"):
        io.read_counts(path)
    pd.DataFrame({"a": [1, 2]}, index=["m1", "m1"]).to_csv(path, sep="\t")
    with pytest.raises(io.InputError, match="duplicated miRNA"):
        io.read_counts(path)
