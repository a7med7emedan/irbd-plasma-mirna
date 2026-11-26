# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Key numbers of a run and their comparison with the values reported in the manuscript."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

COMPARE_KINDS = ("exact", "round", "below", "text")


def load_manuscript_values(path: str | Path) -> dict[str, dict[str, Any]]:
    """Read the registry of reported values. Each entry has a value and a comparison rule."""
    with Path(path).open(encoding="utf-8") as handle:
        values = yaml.safe_load(handle) or {}
    for key, spec in values.items():
        if not isinstance(spec, dict) or "value" not in spec:
            raise ValueError(f"Manuscript value '{key}' needs a mapping with a 'value' entry.")
        kind = spec.get("compare", "exact")
        if kind not in COMPARE_KINDS:
            raise ValueError(f"Manuscript value '{key}' has unknown compare rule '{kind}'.")
        if kind == "round" and "digits" not in spec:
            raise ValueError(f"Manuscript value '{key}' uses 'round' and needs 'digits'.")
    return values


def _is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def compare_value(pipeline: Any, spec: dict[str, Any]) -> str:
    """Return 'match' or 'differs' for one number under its comparison rule."""
    kind = spec.get("compare", "exact")
    reported = spec["value"]
    if _is_missing(pipeline):
        return "differs"
    if kind == "text":
        return "match" if str(pipeline) == str(reported) else "differs"
    if isinstance(pipeline, str):
        return "differs"
    if kind == "exact":
        return "match" if float(pipeline) == float(reported) else "differs"
    if kind == "round":
        digits = int(spec["digits"])
        return (
            "match"
            if round(float(pipeline), digits) == round(float(reported), digits)
            else "differs"
        )
    return "match" if float(pipeline) < float(reported) else "differs"


def build_table(
    numbers: dict[str, Any], manuscript: dict[str, dict[str, Any]] | None
) -> pd.DataFrame:
    """One row per key number: pipeline value, reported value and status."""
    rows = []
    manuscript = manuscript or {}
    for key, value in numbers.items():
        spec = manuscript.get(key)
        if spec is None:
            reported, rule, status = "", "", "not reported"
        else:
            rule = spec.get("compare", "exact")
            reported = f"< {spec['value']}" if rule == "below" else spec["value"]
            status = compare_value(value, spec)
        rows.append(
            {
                "key": key,
                "pipeline_value": value,
                "manuscript_value": reported,
                "rule": rule,
                "status": status,
            }
        )
    for key, spec in manuscript.items():
        if key not in numbers:
            rows.append(
                {
                    "key": key,
                    "pipeline_value": "not computed",
                    "manuscript_value": spec["value"],
                    "rule": spec.get("compare", "exact"),
                    "status": "differs",
                }
            )
    return pd.DataFrame(rows)


def count_status(table: pd.DataFrame) -> dict[str, int]:
    """Number of rows per status."""
    return {s: int((table["status"] == s).sum()) for s in ("match", "differs", "not reported")}
