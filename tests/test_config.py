# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import pytest
import yaml

from irbdmirna.config import ConfigError, load_config


def _write(tmp_path, content):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(content), encoding="utf-8")
    return path


def test_defaults_follow_the_methods(tmp_path):
    cfg = load_config(_write(tmp_path, {"inputs": {"counts": "c.tsv", "samples": "s.tsv"}}))
    assert cfg["thresholds"] == {"fdr": 0.05, "abs_log2fc": 0.5}
    assert cfg["deseq2"]["covariates"] == ["age", "sex"]
    assert cfg["deseq2"]["prefilter"]["irbd_vs_controls"] == "filterByExpr"
    assert cfg["deseq2"]["prefilter"]["irbd_vs_pd"] == "filterByExpr"
    assert cfg["deseq2"]["prefilter"]["pd_vs_controls"] == "min_count"
    assert (cfg["deseq2"]["min_count"], cfg["deseq2"]["min_samples"]) == (10, 3)
    assert cfg["converters"]["n_permutations"] == 200
    assert cfg["converters"]["min_detection_fraction"] == 0.2
    assert cfg["resources"]["brain_enriched"].endswith("brain_enriched_mirnas.tsv")


def test_missing_input_is_reported(tmp_path):
    with pytest.raises(ConfigError, match="inputs.samples"):
        load_config(_write(tmp_path, {"inputs": {"counts": "c.tsv"}}))


def test_unknown_key_is_rejected(tmp_path):
    content = {"inputs": {"counts": "c.tsv", "samples": "s.tsv"}, "treshold": 1}
    with pytest.raises(ConfigError, match="Unknown configuration key 'treshold'"):
        load_config(_write(tmp_path, content))


def test_invalid_prefilter_rule(tmp_path):
    content = {
        "inputs": {"counts": "c.tsv", "samples": "s.tsv"},
        "deseq2": {"prefilter": {"pd_vs_controls": "none"}},
    }
    with pytest.raises(ConfigError, match="filterByExpr, min_count"):
        load_config(_write(tmp_path, content))


def test_missing_file(tmp_path):
    with pytest.raises(ConfigError, match="not found"):
        load_config(tmp_path / "absent.yaml")
