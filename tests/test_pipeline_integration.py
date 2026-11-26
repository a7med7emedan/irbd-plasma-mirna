# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Full pipeline on synthetic data. Needs Rscript with DESeq2, edgeR and apeglm."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from irbdmirna import cli, cluster14q32, report, synthetic
from irbdmirna.config import packaged_resource

pytestmark = pytest.mark.integration
REGISTRY = Path(__file__).resolve().parents[1] / "config" / "manuscript_values.yaml"


def _config(tmp_path: Path, name: str, manuscript: Path | None = None) -> Path:
    paths = synthetic.write_dataset(tmp_path / "data", seed=20261126)
    cfg = {
        "inputs": {k: str(v) for k, v in paths.items()},
        "output_dir": str(tmp_path / name),
        "manuscript_values": str(manuscript) if manuscript else None,
        "converters": {"n_permutations": 20},
        "clinical": {
            "continuous": ["moca", "disease_duration_years"],
            "categorical": ["cardiovascular_disease"],
        },
    }
    path = tmp_path / f"{name}.yaml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    return path


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory, rscript) -> Path:
    tmp = tmp_path_factory.mktemp("pipeline")
    assert cli.main(["run", "--config", str(_config(tmp, "run1", REGISTRY))]) == 0
    return tmp


def _numbers(out: Path) -> pd.Series:
    table = pd.read_csv(out / "summary" / "key_numbers.tsv", sep="\t", dtype=str)
    return table.set_index("key")["pipeline_value"]


def test_every_output_is_written(run_dir):
    out = run_dir / "run1"
    expected = [
        "qc/library_depth.tsv",
        "de/comparison_summary.tsv",
        "de/R_session_info.txt",
        "de/irbd_vs_controls_all.tsv",
        "de/irbd_vs_pd_significant.tsv",
        "de/pd_vs_controls_annotated.tsv",
        "de/hc_high_vs_low_prs_all.tsv",
        "brain_enriched/over_representation.tsv",
        "brain_enriched/shared_candidates.tsv",
        "cluster_14q32/cluster_counts.tsv",
        "cluster_14q32/other_domain_mirnas.tsv",
        "haemolysis/score_models.tsv",
        "converters/permutation_calibration.tsv",
        "batch/per_mirna_batch_test.tsv",
        "clinical/missing_data.tsv",
        "figures/figure_1_volcano_irbd_vs_controls.pdf",
        "figures/figure_2_brain_enriched_panel.png",
        "figures/figure_s1_volcano_hc_high_vs_low_prs.png",
        "figures/figure_s2_converters.png",
        "figures/figure_s3_batch_pca.png",
        "summary/key_numbers.tsv",
        "summary/run_info.yaml",
    ]
    missing = [p for p in expected if not (out / p).is_file()]
    assert not missing


def test_models_follow_the_methods(run_dir):
    summary = pd.read_csv(run_dir / "run1" / "de" / "comparison_summary.tsv", sep="\t")
    summary = summary.set_index("comparison")
    assert (summary["design"] == "~age + sex + cond").all()
    assert summary.loc["irbd_vs_controls", "prefilter"] == "filterByExpr"
    assert summary.loc["irbd_vs_pd", "prefilter"] == "filterByExpr"
    assert summary.loc["pd_vs_controls", "prefilter"] == "min_count"
    assert summary.loc["irbd_vs_controls", ["n_numerator", "n_reference"]].tolist() == [14, 20]
    assert summary.loc["pd_vs_controls", "reference"] == "HC_high_PRS+HC_low_PRS"
    assert (summary["n_with_padj"] <= summary["n_after_prefilter"]).all()
    assert (summary["n_after_prefilter"] < summary["n_input"]).all()


def test_planted_effects_are_recovered(run_dir):
    out = run_dir / "run1"
    n = _numbers(out)
    resource = cluster14q32.read_resource(packaged_resource("mirbase_v22_chr14_DLK1_DIO3.tsv"))
    members = cluster14q32.cluster_members(resource)
    irbd = pd.read_csv(out / "de" / "irbd_vs_controls_annotated.tsv", sep="\t").set_index("miRNA")
    in_cluster = irbd[irbd.index.isin(members) & irbd["tested"]]
    assert len(in_cluster) >= 40
    assert (in_cluster["direction"] == "down").mean() > 0.9  # the cluster is planted lower in iRBD
    assert irbd.loc["hsa-miR-7-5p", "direction"] == "up"
    assert int(n["de.hc_high_vs_low_prs.dysregulated"]) <= 2  # no planted difference
    assert int(n["de.pd_vs_controls.dysregulated"]) >= 3  # three planted miRNAs
    assert int(n["cluster.pd_vs_controls.dysregulated"]) == 0
    assert int(n["cluster.irbd_vs_controls.down"]) == int((in_cluster["direction"] == "down").sum())
    assert int(n["panel.shared"]) == int(n["panel.shared_concordant"])
    assert int(n["converters.dysregulated"]) == 0


def test_significant_tables_obey_the_rule(run_dir):
    out = run_dir / "run1"
    for name in ("irbd_vs_controls", "irbd_vs_pd", "pd_vs_controls"):
        sig = pd.read_csv(out / "de" / f"{name}_significant.tsv", sep="\t")
        assert (sig["padj"] < 0.05).all()
        assert (sig["log2FoldChange"].abs() >= 0.5).all()
        every = pd.read_csv(out / "de" / f"{name}_annotated.tsv", sep="\t")
        assert len(sig) == int(every["dysregulated"].sum())


def test_every_manuscript_key_is_computed(run_dir):
    table = pd.read_csv(run_dir / "run1" / "summary" / "key_numbers.tsv", sep="\t", dtype=str)
    reported = set(report.load_manuscript_values(REGISTRY))
    computed = set(table.loc[table["pipeline_value"] != "not computed", "key"])
    assert reported <= computed


def test_run_is_reproducible(run_dir):
    assert cli.main(["run", "--config", str(_config(run_dir, "run2"))]) == 0
    first, second = _numbers(run_dir / "run1"), _numbers(run_dir / "run2")
    pd.testing.assert_series_equal(first, second, check_names=False)


def test_skip_de_reuses_tables_and_strict_mode_fails_on_synthetic_data(run_dir):
    cfg = str(_config(run_dir, "run1", REGISTRY))
    assert cli.main(["run", "--config", cfg, "--skip-de"]) == 0
    assert cli.main(["run", "--config", cfg, "--skip-de", "--strict"]) == 1
