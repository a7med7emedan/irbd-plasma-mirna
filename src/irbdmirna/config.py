# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Pipeline configuration: defaults, loading and validation."""

from __future__ import annotations

import copy
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

PREFILTER_RULES = ("filterByExpr", "min_count")

DEFAULTS: dict[str, Any] = {
    "inputs": {"counts": None, "samples": None, "clinical": None},
    "resources": {"brain_enriched": None, "cluster_14q32": None},
    "manuscript_values": None,
    "output_dir": "results",
    "rscript": "Rscript",
    "groups": {
        "irbd": "iRBD",
        "pd": "de_novo_PD",
        "hc_high_prs": "HC_high_PRS",
        "hc_low_prs": "HC_low_PRS",
    },
    "thresholds": {"fdr": 0.05, "abs_log2fc": 0.5},
    "deseq2": {
        "covariates": ["age", "sex"],
        "independent_filtering_alpha": 0.1,
        "min_count": 10,
        "min_samples": 3,
        "prefilter": {
            "hc_high_vs_low_prs": "filterByExpr",
            "irbd_vs_controls": "filterByExpr",
            "irbd_vs_pd": "filterByExpr",
            "pd_vs_controls": "min_count",
        },
    },
    "haemolysis": {"numerator": "hsa-miR-451a", "denominator": "hsa-miR-23a-3p"},
    "converters": {
        "min_detection_fraction": 0.2,
        "min_detected_samples": 3,
        "n_permutations": 200,
        "seed": 20260830,
    },
    "batch": {"min_detection_fraction": 0.5, "n_components": 3},
    "clinical": {"continuous": [], "categorical": []},
}


class ConfigError(ValueError):
    """Raised when the configuration is incomplete or inconsistent."""


def _merge(base: dict[str, Any], override: dict[str, Any], path: str = "") -> dict[str, Any]:
    for key, value in override.items():
        where = f"{path}.{key}" if path else key
        if key not in base:
            raise ConfigError(
                f"Unknown configuration key '{where}'. Remove it or check the spelling."
            )
        if isinstance(base[key], dict) and isinstance(value, dict):
            _merge(base[key], value, where)
        elif isinstance(base[key], dict) and value is not None:
            raise ConfigError(f"Configuration key '{where}' must be a mapping.")
        else:
            base[key] = value
    return base


def packaged_resource(name: str) -> Path:
    """Return the path of a resource file shipped with the package."""
    return Path(str(resources.files("irbdmirna").joinpath("resources", name)))


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML configuration file, apply defaults and validate it."""
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"Configuration file not found: {path}")
    with path.open(encoding="utf-8") as handle:
        user = yaml.safe_load(handle) or {}
    if not isinstance(user, dict):
        raise ConfigError("The configuration file must contain a YAML mapping at the top level.")
    cfg = _merge(copy.deepcopy(DEFAULTS), user)
    validate_config(cfg)
    return cfg


def validate_config(cfg: dict[str, Any]) -> None:
    """Check required entries and value ranges. Raise ConfigError with the fix."""
    for key in ("counts", "samples"):
        if not cfg["inputs"].get(key):
            raise ConfigError(f"Set 'inputs.{key}' to the path of the {key} table.")
    if cfg["resources"]["brain_enriched"] is None:
        cfg["resources"]["brain_enriched"] = str(packaged_resource("brain_enriched_mirnas.tsv"))
    if cfg["resources"]["cluster_14q32"] is None:
        cfg["resources"]["cluster_14q32"] = str(
            packaged_resource("mirbase_v22_chr14_DLK1_DIO3.tsv")
        )
    labels = list(cfg["groups"].values())
    if len(set(labels)) != 4:
        raise ConfigError("The four entries under 'groups' must have four distinct labels.")
    fdr = cfg["thresholds"]["fdr"]
    if not 0 < fdr < 1:
        raise ConfigError("'thresholds.fdr' must lie between 0 and 1.")
    if cfg["thresholds"]["abs_log2fc"] < 0:
        raise ConfigError("'thresholds.abs_log2fc' must not be negative.")
    alpha = cfg["deseq2"]["independent_filtering_alpha"]
    if not 0 < alpha < 1:
        raise ConfigError("'deseq2.independent_filtering_alpha' must lie between 0 and 1.")
    for name, rule in cfg["deseq2"]["prefilter"].items():
        if rule not in PREFILTER_RULES:
            raise ConfigError(
                f"'deseq2.prefilter.{name}' is '{rule}'. Use one of: {', '.join(PREFILTER_RULES)}."
            )
    if cfg["converters"]["n_permutations"] < 0:
        raise ConfigError("'converters.n_permutations' must not be negative.")
    for key in ("min_detection_fraction",):
        for section in ("converters", "batch"):
            if not 0 < cfg[section][key] <= 1:
                raise ConfigError(f"'{section}.{key}' must lie in (0, 1].")
