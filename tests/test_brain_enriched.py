# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import pytest
from scipy import stats

from conftest import annotated_table
from irbdmirna import brain_enriched
from irbdmirna.config import packaged_resource

NAN = float("nan")


def test_packaged_panel_has_thirty_unique_mirnas():
    panel = brain_enriched.read_panel(packaged_resource("brain_enriched_mirnas.tsv"))
    assert len(panel) == 30
    assert len(set(panel)) == 30
    assert all(name.startswith("hsa-miR-") for name in panel)
    assert {"hsa-miR-7-5p", "hsa-miR-3200-3p", "hsa-miR-127-3p"} <= set(panel)


def test_panel_status_distinguishes_untested_reasons():
    table = annotated_table(
        {"p1": (-1.0, 0.001), "p2": (0.1, 0.8), "p3": (2.0, NAN), "x": (1, 0.01)}
    )
    status = brain_enriched.panel_status(table, ["p1", "p2", "p3", "p4"])
    assert status["status"].tolist() == [
        "tested",
        "tested",
        "not tested (independent filtering)",
        "not tested (low counts)",
    ]
    assert status["dysregulated"].tolist() == [True, False, False, False]


def test_over_representation_matches_fisher():
    rows = {f"p{i}": (-1.0, 0.001) for i in range(4)}
    rows |= {f"q{i}": (0.0, 0.9) for i in range(2)}
    rows |= {f"b{i}": (1.0, 0.001) for i in range(3)}
    rows |= {f"n{i}": (0.0, 0.9) for i in range(11)}
    rows |= {"untested": (3.0, NAN)}
    table = annotated_table(rows)
    panel = [f"p{i}" for i in range(4)] + [f"q{i}" for i in range(2)]
    res = brain_enriched.over_representation(table, panel)
    assert (res["panel_tested"], res["panel_dysregulated"]) == (6, 4)
    assert (res["background_tested"], res["background_dysregulated"]) == (14, 3)
    expected = stats.fisher_exact([[4, 2], [3, 11]], alternative="greater").pvalue
    assert res["p_value"] == pytest.approx(expected)
    without = brain_enriched.over_representation(table, panel, exclude={"p0", "p1", "b0"})
    assert (without["panel_tested"], without["panel_dysregulated"]) == (4, 2)
    assert (without["background_tested"], without["background_dysregulated"]) == (13, 2)


def test_shared_candidates_and_concordance():
    first = annotated_table(
        {"a": (-1.0, 0.001), "b": (1.0, 0.001), "c": (1.0, 0.001), "d": (0, 0.9)}
    )
    second = annotated_table({"a": (-2.0, 0.01), "b": (-1.0, 0.01), "c": (0.1, 0.5), "d": (0, 0.9)})
    panel = ["a", "b", "c", "d"]
    shared = brain_enriched.shared_candidates(
        brain_enriched.panel_status(first, panel), brain_enriched.panel_status(second, panel)
    )
    assert shared.index.tolist() == ["a", "b"]
    assert shared["concordant"].tolist() == [True, False]
