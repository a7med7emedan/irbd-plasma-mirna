# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Synthetic example data with the same shape as the study inputs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import brain_enriched, cluster14q32
from .config import packaged_resource

GROUP_SIZES = {"iRBD": 14, "de_novo_PD": 14, "HC_high_PRS": 10, "HC_low_PRS": 10}
N_FAILED_QC = 4
N_BACKGROUND = 300
LOW_PANEL = ("hsa-miR-132-3p", "hsa-miR-329-5p", "hsa-miR-431-3p")
UP_PANEL = ("hsa-miR-7-5p", "hsa-miR-124-3p", "hsa-miR-3200-3p")


def _negative_binomial(rng: np.random.Generator, mean: np.ndarray) -> np.ndarray:
    """Negative binomial counts with a dispersion that falls with the mean, as in sequencing data."""
    dispersion = 0.08 + 1.5 / np.maximum(mean, 0.05)
    shape = 1.0 / dispersion
    return rng.poisson(rng.gamma(shape, mean / shape))


def make_dataset(seed: int = 20261126) -> dict[str, pd.DataFrame]:
    """Create synthetic counts, a sample sheet and a clinical table."""
    rng = np.random.default_rng(seed)
    resource = cluster14q32.read_resource(packaged_resource("mirbase_v22_chr14_DLK1_DIO3.tsv"))
    cluster = sorted(cluster14q32.cluster_members(resource))
    other_domain = sorted(cluster14q32.other_domain_members(resource))
    panel = brain_enriched.read_panel(packaged_resource("brain_enriched_mirnas.tsv"))
    fixed = ["hsa-miR-451a", "hsa-miR-23a-3p", "hsa-miR-16-5p", "hsa-miR-486-5p"]
    background = [f"hsa-miR-syn{i:04d}" for i in range(1, N_BACKGROUND + 1)]
    names = list(dict.fromkeys(cluster + other_domain + panel + fixed + background))

    rows = []
    for group, n in GROUP_SIZES.items():
        rows += [group] * n
    extra = list(rng.choice(list(GROUP_SIZES), size=N_FAILED_QC))
    groups = np.array(rows + extra)
    n_lib = len(groups)
    qc_pass = np.array([True] * len(rows) + [False] * N_FAILED_QC)
    samples = pd.DataFrame(
        {
            "library_id": [f"LIB{i:03d}" for i in range(1, n_lib + 1)],
            "group": groups,
            "age": np.round(rng.normal(65, 8, n_lib), 1),
            "sex": rng.choice(["male", "female"], size=n_lib, p=[0.7, 0.3]),
            "batch": rng.choice(["batch1", "batch2"], size=n_lib, p=[0.2, 0.8]),
            "qc_pass": qc_pass,
            "converter": False,
        }
    )
    irbd_rows = samples.index[(samples["group"] == "iRBD") & samples["qc_pass"]]
    samples.loc[irbd_rows[:4], "converter"] = True

    base = pd.Series(np.exp(rng.uniform(np.log(20), np.log(3000), len(names))), index=names)
    base[background[:40]] = rng.uniform(0.2, 3.0, 40)  # low counts, removed by the pre-filter
    base[background[40:50]] = 0.0  # never detected
    base[background[110:160]] = rng.uniform(3.0, 25.0, 50)  # near the pre-filter limits
    base[cluster] = np.exp(rng.uniform(np.log(30), np.log(600), len(cluster)))
    base[cluster[::3]] = rng.uniform(0.2, 2.0, len(cluster[::3]))  # low members, not testable
    for name in panel:
        if name not in cluster:
            base[name] = rng.uniform(80, 800)
    base[list(LOW_PANEL)] = 0.5
    base[["hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p"]] = [60000.0, 20000.0, 9000.0]
    base["hsa-miR-23a-3p"] = 2500.0
    for name in other_domain:
        base[name] = 0.3
    if "hsa-miR-1247-5p" in base.index:
        base["hsa-miR-1247-5p"] = 25.0

    fold = pd.DataFrame(1.0, index=names, columns=samples["library_id"])
    irbd = (samples["group"] == "iRBD").to_numpy()
    pd_mask = (samples["group"] == "de_novo_PD").to_numpy()
    fold.loc[cluster, irbd] = 2.0**-1.3
    fold.loc[background[60:90], irbd] = 2.0**1.0
    fold.loc[[n for n in UP_PANEL if n in fold.index], irbd] = 2.0**0.9
    fold.loc[background[100:103], pd_mask] = 2.0**1.2

    haem = np.exp(rng.normal(0.0, 0.5, n_lib) + 0.3 * irbd)
    fold.loc[["hsa-miR-451a", "hsa-miR-16-5p", "hsa-miR-486-5p"]] *= haem
    fold.loc[cluster] *= haem**-0.3

    depth = np.exp(rng.normal(0.0, 0.3, n_lib))
    mean = base.to_numpy()[:, None] * fold.to_numpy() * depth[None, :]
    counts = _negative_binomial(rng, mean).astype(float)
    half = (rng.random(counts.shape) < 0.005) & (counts > 0)
    counts[half] += 0.5  # reads shared between miRNAs give non-integer values
    counts_table = pd.DataFrame(
        counts, index=pd.Index(names, name="miRNA"), columns=samples["library_id"]
    )

    passed = samples[samples["qc_pass"]]
    clinical = pd.DataFrame(
        {
            "library_id": passed["library_id"].to_numpy(),
            "group": passed["group"].to_numpy(),
            "moca": np.round(rng.normal(26, 2.5, len(passed)), 0),
            "disease_duration_years": np.round(rng.gamma(2.0, 1.5, len(passed)), 1),
            "cardiovascular_disease": rng.choice(["no", "yes"], size=len(passed), p=[0.85, 0.15]),
        }
    )
    clinical.loc[rng.choice(clinical.index, size=5, replace=False), "moca"] = np.nan
    clinical.loc[clinical["group"].str.startswith("HC"), "disease_duration_years"] = np.nan
    return {"counts": counts_table, "samples": samples, "clinical": clinical}


def write_dataset(outdir: str | Path, seed: int = 20261126) -> dict[str, Path]:
    """Write the synthetic tables as tab-separated files and return their paths."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = make_dataset(seed)
    paths = {
        "counts": outdir / "counts.tsv",
        "samples": outdir / "samples.tsv",
        "clinical": outdir / "clinical.tsv",
    }
    data["counts"].to_csv(paths["counts"], sep="\t")
    data["samples"].to_csv(paths["samples"], sep="\t", index=False)
    data["clinical"].to_csv(paths["clinical"], sep="\t", index=False)
    return paths
