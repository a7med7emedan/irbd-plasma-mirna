# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Shared fixtures. Every test builds its own data."""

from __future__ import annotations

import shutil
import subprocess

import pandas as pd
import pytest

from irbdmirna import synthetic


@pytest.fixture(scope="session")
def dataset() -> dict[str, pd.DataFrame]:
    return synthetic.make_dataset(seed=7)


@pytest.fixture(scope="session")
def rscript() -> str:
    """Path to Rscript with the required packages, or skip."""
    path = shutil.which("Rscript")
    if path is None:
        pytest.skip("Rscript is not installed")
    code = "suppressMessages({library(DESeq2);library(edgeR);library(apeglm);library(yaml)})"
    if subprocess.run([path, "-e", code], capture_output=True, check=False).returncode != 0:
        pytest.skip("R packages DESeq2, edgeR, apeglm or yaml are not installed")
    return path


def annotated_table(rows: dict[str, tuple[float, float]]) -> pd.DataFrame:
    """Small annotated DESeq2-like table from name -> (log2FC, padj)."""
    from irbdmirna import dysregulation

    table = pd.DataFrame(
        {
            "baseMean": 100.0,
            "log2FoldChange": [v[0] for v in rows.values()],
            "lfcSE": 0.1,
            "stat": 1.0,
            "pvalue": [v[1] / 2 if v[1] == v[1] else float("nan") for v in rows.values()],
            "padj": [v[1] for v in rows.values()],
        },
        index=pd.Index(list(rows), name="miRNA"),
    )
    return dysregulation.annotate(table)
