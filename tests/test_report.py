# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
from pathlib import Path

import pytest

from irbdmirna import report

REGISTRY = Path(__file__).resolve().parents[1] / "config" / "manuscript_values.yaml"


@pytest.mark.parametrize(
    ("value", "spec", "expected"),
    [
        (250, {"value": 250}, "match"),
        (249, {"value": 250}, "differs"),
        (0.0074, {"value": 0.007, "compare": "round", "digits": 3}, "match"),
        (0.0076, {"value": 0.007, "compare": "round", "digits": 3}, "differs"),
        (-0.634, {"value": -0.63, "compare": "round", "digits": 2}, "match"),
        (31.7, {"value": 32, "compare": "round", "digits": 0}, "match"),
        (0.0004, {"value": 0.001, "compare": "below"}, "match"),
        (0.002, {"value": 0.001, "compare": "below"}, "differs"),
        ("not tested", {"value": "not tested", "compare": "text"}, "match"),
        (0.3, {"value": "not tested", "compare": "text"}, "differs"),
        ("not tested", {"value": 0.25, "compare": "round", "digits": 2}, "differs"),
        (float("nan"), {"value": 3}, "differs"),
    ],
)
def test_compare_value(value, spec, expected):
    assert report.compare_value(value, spec) == expected


def test_build_table_lists_every_key():
    numbers = {"a": 1, "b": 0.0004, "c": 5}
    manuscript = {"a": {"value": 1}, "b": {"value": 0.001, "compare": "below"}, "d": {"value": 9}}
    table = report.build_table(numbers, manuscript).set_index("key")
    assert table["status"].to_dict() == {
        "a": "match",
        "b": "match",
        "c": "not reported",
        "d": "differs",
    }
    assert table.loc["d", "pipeline_value"] == "not computed"
    assert report.count_status(table.reset_index()) == {"match": 2, "differs": 1, "not reported": 1}


def test_manuscript_registry_is_valid_and_holds_the_reported_counts():
    values = report.load_manuscript_values(REGISTRY)
    assert values["samples.analysed"]["value"] == 165
    assert values["de.irbd_vs_controls.dysregulated"]["value"] == 250
    assert values["de.irbd_vs_pd.dysregulated"]["value"] == 161
    assert values["de.pd_vs_controls.dysregulated"]["value"] == 9
    assert values["cluster.irbd_vs_controls.tested"]["value"] == 57
    assert values["cluster.irbd_vs_controls.down"]["value"] == 54
    up_down = [values[f"de.irbd_vs_controls.{k}"]["value"] for k in ("up", "down")]
    assert sum(up_down) == values["de.irbd_vs_controls.dysregulated"]["value"]
    up_down = [values[f"de.irbd_vs_pd.{k}"]["value"] for k in ("up", "down")]
    assert sum(up_down) == values["de.irbd_vs_pd.dysregulated"]["value"]
    groups = sum(
        values[f"samples.{k}"]["value"] for k in ("irbd", "pd", "hc_high_prs", "hc_low_prs")
    )
    assert groups == values["samples.analysed"]["value"]


def test_invalid_registry_entries(tmp_path):
    path = tmp_path / "values.yaml"
    path.write_text("a: {value: 1, compare: nearly}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unknown compare rule"):
        report.load_manuscript_values(path)
    path.write_text("a: {value: 1.2, compare: round}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="needs 'digits'"):
        report.load_manuscript_values(path)
