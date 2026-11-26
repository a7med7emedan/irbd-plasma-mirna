# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""The clinical step writes the five comparisons and five missing-data columns of the manuscript."""

import pandas as pd
import pytest

from irbdmirna import clinical, config, pipeline, synthetic


@pytest.fixture
def run(tmp_path):
    paths = synthetic.write_dataset(tmp_path / "data")
    cfg = config.DEFAULTS | {
        "inputs": {"counts": None, "samples": None, "clinical": str(paths["clinical"])},
        "clinical": {
            "continuous": ["moca", "disease_duration_years"],
            "categorical": ["cardiovascular_disease"],
        },
    }
    samples = pd.read_csv(paths["samples"], sep="\t", dtype={"library_id": str})
    samples = samples[samples["qc_pass"]].reset_index(drop=True)
    group3 = pd.Series(pipeline._group3(samples, cfg).to_numpy(), index=samples["library_id"])
    out = tmp_path / "out"
    pipeline.step_clinical(cfg, {"samples": samples, "group3": group3}, out)
    return cfg, samples, out


def test_five_comparisons_as_in_table_1(run):
    _, _, out = run
    names = sorted(p.name for p in (out / "clinical").glob("comparison_*.tsv"))
    assert names == [
        "comparison_de_novo_PD_vs_controls.tsv",
        "comparison_iRBD_vs_controls.tsv",
        "comparison_iRBD_vs_de_novo_PD.tsv",
        "comparison_iRBD_vs_high-PRS_controls.tsv",
        "comparison_iRBD_vs_low-PRS_controls.tsv",
    ]


def test_stratum_comparison_equals_direct_test(run):
    cfg, samples, out = run
    g = cfg["groups"]
    table = pd.read_csv(cfg["inputs"]["clinical"], sep="\t", dtype={"library_id": str})
    keep = samples.loc[samples["group"].isin([g["irbd"], g["hc_high_prs"]]), "library_id"]
    sub = table[table["library_id"].isin(keep)].merge(
        samples[["library_id", "group"]], on="library_id", suffixes=("_table", "")
    )
    expected = clinical.compare_groups(
        sub,
        "group",
        [g["irbd"]],
        [g["hc_high_prs"]],
        ["moca", "disease_duration_years"],
        ["cardiovascular_disease"],
    )
    got = pd.read_csv(out / "clinical" / "comparison_iRBD_vs_high-PRS_controls.tsv", sep="\t")
    assert got["p_value"].tolist() == pytest.approx(expected["p_value"].tolist(), nan_ok=True)
    assert got["p_adjusted"].tolist() == pytest.approx(expected["p_adjusted"].tolist(), nan_ok=True)


def test_missing_data_has_the_columns_of_table_s5(run):
    _, _, out = run
    missing = pd.read_csv(out / "clinical" / "missing_data.tsv", sep="\t")
    assert missing.columns.tolist() == [
        "variable",
        "iRBD",
        "de novo PD",
        "high-PRS controls",
        "low-PRS controls",
        "all controls",
    ]
    counts = missing.set_index("variable").map(lambda v: int(str(v).split(" ")[0]))
    total = counts["high-PRS controls"] + counts["low-PRS controls"]
    assert counts["all controls"].tolist() == total.tolist()
