# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Run every analysis step in order and collect the key numbers."""

from __future__ import annotations

import hashlib
import logging
import platform
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy
import statsmodels
import yaml

from . import (
    __version__,
    batch,
    brain_enriched,
    clinical,
    cluster14q32,
    converters,
    de,
    dysregulation,
    figures,
    haemolysis,
    io,
    report,
)

log = logging.getLogger(__name__)

CONTROLS = "controls"
HIGH_PRS, LOW_PRS = "high-PRS controls", "low-PRS controls"
ALL_CONTROLS = "all controls"


class ManuscriptMismatch(RuntimeError):
    """Raised in strict mode when a key number differs from the reported value."""


def _write(table: pd.DataFrame, path: Path, index: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=index)


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _group3(samples: pd.DataFrame, cfg: dict[str, Any]) -> pd.Series:
    g = cfg["groups"]
    mapping = {
        g["irbd"]: "iRBD",
        g["pd"]: "de novo PD",
        g["hc_high_prs"]: CONTROLS,
        g["hc_low_prs"]: CONTROLS,
    }
    return samples["group"].map(mapping)


def step_inputs(cfg: dict[str, Any], out: Path, numbers: dict[str, Any]) -> dict[str, Any]:
    """Read and validate the inputs, round the counts and summarise detection."""
    g = cfg["groups"]
    counts_all = io.read_counts(cfg["inputs"]["counts"])
    samples_all = io.read_samples(cfg["inputs"]["samples"], list(g.values()))
    samples = io.analysis_set(samples_all, counts_all)
    raw = counts_all.loc[:, samples["library_id"]]
    rounded, non_integer = io.round_counts(raw)
    detection = io.detection_summary(raw)
    depth = raw.sum(axis=0) / 1e6
    numbers["samples.sequenced"] = int(len(samples_all))
    numbers["samples.analysed"] = int(len(samples))
    for key in ("irbd", "pd", "hc_high_prs", "hc_low_prs"):
        numbers[f"samples.{key}"] = int((samples["group"] == g[key]).sum())
    numbers["counts.percent_non_integer"] = 100.0 * non_integer
    numbers["mirnas.annotated"] = detection["annotated"]
    numbers["mirnas.detected_in_at_least_one"] = detection["detected_in_at_least_one"]
    numbers["mirnas.detected_in_at_least_half"] = detection["detected_in_at_least_half"]
    numbers["reads.median_million"] = float(depth.median())
    numbers["reads.min_million"] = float(depth.min())
    numbers["reads.max_million"] = float(depth.max())
    summary = pd.DataFrame(
        {
            "library_id": samples["library_id"],
            "group": samples["group"],
            "mirna_reads": raw.sum(axis=0).to_numpy(),
        }
    )
    _write(summary, out / "qc" / "library_depth.tsv", index=False)
    log.info("Inputs: %d libraries analysed, %d miRNAs annotated", len(samples), len(raw))
    return {
        "samples": samples,
        "raw": raw,
        "rounded": rounded,
        "rpm": io.rpmmm(raw),
        "group3": pd.Series(_group3(samples, cfg).to_numpy(), index=samples["library_id"]),
    }


def step_de(
    cfg: dict[str, Any], data: dict[str, Any], out: Path, numbers: dict[str, Any], skip_de: bool
) -> dict[str, pd.DataFrame]:
    """Fit the four DESeq2 comparisons and apply the dysregulation rule."""
    de_dir = out / "de"
    if not skip_de:
        de.run_deseq2(cfg, data["rounded"], data["samples"], de_dir, out / "work")
    tables = {}
    thr = cfg["thresholds"]
    for name, table in de.read_de_tables(de_dir).items():
        annotated = dysregulation.annotate(table, thr["fdr"], thr["abs_log2fc"])
        tables[name] = annotated
        _write(annotated, de_dir / f"{name}_annotated.tsv")
        _write(dysregulation.significant(annotated), de_dir / f"{name}_significant.tsv")
        s = dysregulation.summarise(annotated)
        for key in ("tested", "dysregulated", "up", "down"):
            numbers[f"de.{name}.{key}"] = s[key]
        numbers[f"de.{name}.smallest_padj"] = s["smallest_padj"]
        log.info(
            "%s: %d tested, %d dysregulated (%d up, %d down)",
            name,
            *(s[k] for k in ("tested", "dysregulated", "up", "down")),
        )
    return tables


def step_panel(
    cfg: dict[str, Any],
    tables: dict[str, pd.DataFrame],
    cluster: set[str],
    out: Path,
    numbers: dict[str, Any],
) -> list[str]:
    """Brain-enriched panel: status per comparison, shared candidates, over-representation."""
    panel = brain_enriched.read_panel(cfg["resources"]["brain_enriched"])
    numbers["panel.size"] = len(panel)
    status, tests = {}, []
    for name in de.MAIN_COMPARISONS:
        status[name] = brain_enriched.panel_status(tables[name], panel)
        _write(status[name], out / "brain_enriched" / f"panel_{name}.tsv")
        numbers[f"panel.{name}.tested"] = int((status[name]["status"] == "tested").sum())
        numbers[f"panel.{name}.dysregulated"] = int(status[name]["dysregulated"].sum())
        for label, exclude in (("all miRNAs", None), ("without 14q32 cluster", cluster)):
            res = brain_enriched.over_representation(tables[name], panel, exclude)
            tests.append({"comparison": name, "background": label, **res})
    _write(pd.DataFrame(tests), out / "brain_enriched" / "over_representation.tsv", index=False)
    shared = brain_enriched.shared_candidates(status["irbd_vs_controls"], status["irbd_vs_pd"])
    shared["in_14q32_cluster"] = shared.index.isin(cluster)
    _write(shared, out / "brain_enriched" / "shared_candidates.tsv")
    numbers["panel.shared"] = int(len(shared))
    numbers["panel.shared_concordant"] = int(shared["concordant"].sum())
    numbers["panel.shared_in_cluster"] = int(shared["in_14q32_cluster"].sum())
    return panel


def step_cluster(
    cfg: dict[str, Any], tables: dict[str, pd.DataFrame], out: Path, numbers: dict[str, Any]
) -> set[str]:
    """14q32 cluster counts per comparison and the remaining miRNAs of the domain."""
    resource = cluster14q32.read_resource(cfg["resources"]["cluster_14q32"])
    members = cluster14q32.cluster_members(resource)
    for key, value in cluster14q32.describe(resource).items():
        numbers[f"cluster.{key}"] = value
    rows, other_rows = [], []
    for name in de.MAIN_COMPARISONS:
        c = cluster14q32.count_in_comparison(tables[name], members)
        rows.append({"comparison": name, **c})
        numbers[f"cluster.{name}.tested"] = c["tested"]
        numbers[f"cluster.{name}.dysregulated"] = c["dysregulated"]
        numbers[f"cluster.{name}.down"] = c["down"]
        for mirna in sorted(cluster14q32.other_domain_members(resource)):
            t = tables[name]
            tested = mirna in t.index and bool(t.loc[mirna, "tested"])
            other_rows.append(
                {
                    "miRNA": mirna,
                    "comparison": name,
                    "status": "tested" if tested else "not tested",
                    "log2FoldChange": t.loc[mirna, "log2FoldChange"] if tested else np.nan,
                    "padj": t.loc[mirna, "padj"] if tested else np.nan,
                }
            )
            short = mirna.replace("hsa-miR-", "mir").replace("-", "_")
            numbers[f"domain.{short}.{name}.log2fc"] = (
                float(t.loc[mirna, "log2FoldChange"]) if tested else "not tested"
            )
            numbers[f"domain.{short}.{name}.padj"] = (
                float(t.loc[mirna, "padj"]) if tested else "not tested"
            )
    _write(pd.DataFrame(rows), out / "cluster_14q32" / "cluster_counts.tsv", index=False)
    _write(pd.DataFrame(other_rows), out / "cluster_14q32" / "other_domain_mirnas.tsv", index=False)
    return members


def step_haemolysis(
    cfg: dict[str, Any],
    data: dict[str, Any],
    tables: dict[str, pd.DataFrame],
    cluster: set[str],
    out: Path,
    numbers: dict[str, Any],
) -> None:
    """Haemolysis index, cluster score, correlations and adjusted group differences."""
    g = cfg["groups"]
    samples, rpm, group3 = data["samples"], data["rpm"], data["group3"]
    index = haemolysis.haemolysis_index(
        rpm, cfg["haemolysis"]["numerator"], cfg["haemolysis"]["denominator"]
    )
    ref = tables["irbd_vs_controls"]
    members = sorted(m for m in cluster if m in ref.index and bool(ref.loc[m, "tested"]))
    score = haemolysis.cluster_score(rpm, members)
    per_library = pd.DataFrame(
        {
            "group": samples["group"].to_numpy(),
            "group3": group3.to_numpy(),
            "age": samples["age"].to_numpy(),
            "sex": samples["sex"].to_numpy(),
            "haemolysis_index": index.to_numpy(),
            "cluster_score": score.to_numpy(),
        },
        index=pd.Index(samples["library_id"], name="library_id"),
    )
    hdir = out / "haemolysis"
    _write(per_library, hdir / "per_library.tsv")
    pairs = [("iRBD", CONTROLS), ("iRBD", "de novo PD"), ("de novo PD", CONTROLS)]
    tests = haemolysis.pairwise_mannwhitney(
        per_library["haemolysis_index"], per_library["group3"], pairs
    )
    _write(tests, hdir / "index_group_tests.tsv", index=False)
    keys = ("irbd_vs_controls", "irbd_vs_pd", "pd_vs_controls")
    for key, p in zip(keys, tests["p_value"], strict=True):
        numbers[f"haemolysis.p.{key}"] = float(p)
    rho = haemolysis.spearman_by_group(
        per_library["haemolysis_index"], per_library["cluster_score"], per_library["group3"]
    )
    _write(rho, hdir / "spearman_index_vs_score.tsv", index=False)
    labels = {CONTROLS: "controls", "de novo PD": "pd", "iRBD": "irbd"}
    for _, row in rho.iterrows():
        numbers[f"haemolysis.rho.{labels[row['group']]}"] = float(row["spearman_rho"])
    models = haemolysis.score_models(per_library, reference=CONTROLS, target="iRBD")
    _write(models, hdir / "score_models.tsv", index=False)
    numbers["score.members"] = len(members)
    numbers["score.percent_reduction_haemolysis"] = float(
        models["percent_reduction_vs_group_only"].iloc[1]
    )
    numbers["score.p_adjusted_haemolysis"] = float(models["p_value"].iloc[1])
    numbers["score.p_adjusted_haemolysis_age_sex"] = float(models["p_value"].iloc[2])
    strata = haemolysis.pairwise_mannwhitney(
        per_library["cluster_score"],
        per_library["group"],
        [
            (g["irbd"], g["hc_high_prs"]),
            (g["irbd"], g["hc_low_prs"]),
            (g["hc_high_prs"], g["hc_low_prs"]),
        ],
    )
    _write(strata, hdir / "score_by_prs_stratum.tsv", index=False)
    numbers["score.p_irbd_vs_hc_high_prs"] = float(strata["p_value"].iloc[0])
    numbers["score.p_irbd_vs_hc_low_prs"] = float(strata["p_value"].iloc[1])
    controls = per_library[per_library["group3"] == CONTROLS]
    age = haemolysis.spearman_by_group(
        controls["age"], controls["cluster_score"], controls["group3"]
    )
    _write(age, hdir / "spearman_age_vs_score_controls.tsv", index=False)
    numbers["score.age_rho_controls"] = float(age["spearman_rho"].iloc[0])
    numbers["score.age_p_controls"] = float(age["p_value"].iloc[0])


def step_converters(
    cfg: dict[str, Any], data: dict[str, Any], panel: list[str], out: Path, numbers: dict[str, Any]
) -> None:
    """Converters versus non-converters within iRBD."""
    samples = data["samples"]
    irbd = samples[samples["group"] == cfg["groups"]["irbd"]]
    is_conv = irbd["converter"].to_numpy(dtype=bool)
    sf_path = out / "de" / "size_factors_within_group.tsv"
    if is_conv.sum() < 2 or (~is_conv).sum() < 2 or not sf_path.is_file():
        log.warning("Converter comparison skipped: it needs two libraries per group in iRBD.")
        return
    c = cfg["converters"]
    thr = cfg["thresholds"]
    sf = pd.read_csv(sf_path, sep="\t", dtype={"library_id": str}).set_index("library_id")[
        "size_factor"
    ]
    counts = data["raw"].loc[:, irbd["library_id"]]
    need = converters.detection_threshold(
        len(irbd), c["min_detection_fraction"], c["min_detected_samples"]
    )
    keep = (counts > 0).sum(axis=1) >= need
    norm = converters.normalise(counts.loc[keep], sf)
    result = converters.compare(norm, is_conv, thr["fdr"], thr["abs_log2fc"])
    result["brain_enriched"] = result.index.isin(panel)
    _write(result, out / "converters" / "converters_vs_non_converters.tsv")
    calib = converters.permutation_calibration(
        norm, is_conv, c["n_permutations"], c["seed"], thr["fdr"], thr["abs_log2fc"]
    )
    _write(pd.DataFrame([calib]), out / "converters" / "permutation_calibration.tsv", index=False)
    numbers["converters.n_converters"] = int(is_conv.sum())
    numbers["converters.n_non_converters"] = int((~is_conv).sum())
    numbers["converters.min_libraries_detected"] = int(need)
    numbers["converters.tested"] = int(len(result))
    numbers["converters.dysregulated"] = int(result["dysregulated"].sum())
    numbers["converters.panel_detected"] = int(result["brain_enriched"].sum())
    ratio = sf.loc[irbd["library_id"]].to_numpy()
    numbers["converters.size_factor_ratio"] = float(
        np.median(ratio[is_conv]) / np.median(ratio[~is_conv])
    )
    numbers["converters.n_permutations"] = calib["n_permutations"]
    numbers["converters.permutations_with_calls"] = calib["permutations_with_calls"]
    annotated = result.assign(
        direction=np.where(
            result["dysregulated"], np.where(result["log2FoldChange"] > 0, "up", "down"), "none"
        )
    )
    figures.volcano(
        annotated,
        panel,
        "iRBD converters vs non-converters",
        out / "figures" / "figure_s2_converters",
        thr["fdr"],
        thr["abs_log2fc"],
        label_all_panel=True,
    )


def step_batch(
    cfg: dict[str, Any], data: dict[str, Any], out: Path, numbers: dict[str, Any]
) -> None:
    """Batch balance, principal components and per-miRNA batch association."""
    samples = data["samples"]
    groups = pd.Series(samples["group"].to_numpy(), index=samples["library_id"])
    batches = pd.Series(samples["batch"].to_numpy(), index=samples["library_id"])
    if batches.nunique() < 2:
        log.warning("Batch analysis skipped: only one batch among the analysed libraries.")
        return
    bdir = out / "batch"
    bal = batch.balance(groups, batches)
    _write(bal["table"], bdir / "batch_by_group.tsv")
    logged = batch.expressed_log2(data["rpm"], cfg["batch"]["min_detection_fraction"])
    scores, comps = batch.principal_components(
        logged, groups, batches, cfg["batch"]["n_components"]
    )
    _write(scores, bdir / "pca_scores.tsv")
    _write(comps, bdir / "pca_components.tsv", index=False)
    per = batch.per_mirna_batch_test(logged, groups, batches)
    _write(per, bdir / "per_mirna_batch_test.tsv")
    sizes = batches.value_counts()
    numbers["batch.smallest_batch_n"] = int(sizes.min())
    numbers["batch.largest_batch_n"] = int(sizes.max())
    numbers["batch.chi2_p"] = bal["p_value"]
    numbers["batch.cramers_v"] = bal["cramers_v"]
    for i in range(min(2, len(comps))):
        numbers[f"batch.pc{i + 1}.p_group"] = float(comps["p_group"].iloc[i])
        numbers[f"batch.pc{i + 1}.p_batch"] = float(comps["p_batch"].iloc[i])
    numbers["batch.mirnas_tested"] = int(len(per))
    numbers["batch.mirnas_associated"] = int((per["padj"] < cfg["thresholds"]["fdr"]).sum())
    numbers["batch.nominal_p_below_0.05"] = int((per["p_value"] < 0.05).sum())
    numbers["batch.nominal_expected"] = int(round(0.05 * len(per)))
    figures.batch_pca(scores, groups, batches, comps, out / "figures" / "figure_s3_batch_pca")


def step_clinical(cfg: dict[str, Any], data: dict[str, Any], out: Path) -> None:
    """Clinical comparisons and missing-data overview, when a clinical table is configured.

    The five comparisons and the five missing-data columns follow Table 1, Supplementary
    Table S1 and Supplementary Table S5 of the manuscript.
    """
    path = cfg["inputs"].get("clinical")
    if not path:
        return
    g = cfg["groups"]
    table = pd.read_csv(path, sep="\t", dtype={"library_id": str})
    table = table[table["library_id"].isin(data["samples"]["library_id"])].copy()
    table["group3"] = table["library_id"].map(data["group3"])
    strata = {
        g["irbd"]: "iRBD",
        g["pd"]: "de novo PD",
        g["hc_high_prs"]: HIGH_PRS,
        g["hc_low_prs"]: LOW_PRS,
    }
    stratum = data["samples"].set_index("library_id")["group"].map(strata)
    table["stratum"] = table["library_id"].map(stratum)
    cont, cat = list(cfg["clinical"]["continuous"]), list(cfg["clinical"]["categorical"])
    missing = [v for v in cont + cat if v not in table.columns]
    if missing:
        raise io.InputError(f"Clinical variables not found in the clinical table: {missing}.")
    comparisons = [
        ("group3", "iRBD", "de novo PD"),
        ("stratum", "iRBD", HIGH_PRS),
        ("stratum", "iRBD", LOW_PRS),
        ("group3", "iRBD", CONTROLS),
        ("group3", "de novo PD", CONTROLS),
    ]
    for column, first, second in comparisons:
        res = clinical.compare_groups(table, column, [first], [second], cont, cat)
        name = f"{first}_vs_{second}".replace(" ", "_")
        _write(res, out / "clinical" / f"comparison_{name}.tsv", index=False)
    by_stratum = clinical.missing_data(table, "stratum", cont + cat).set_index("variable")
    merged = clinical.missing_data(table, "group3", cont + cat).set_index("variable")
    order = [c for c in ("iRBD", "de novo PD", HIGH_PRS, LOW_PRS) if c in by_stratum.columns]
    overview = by_stratum[order]
    if CONTROLS in merged.columns:
        overview = overview.assign(**{ALL_CONTROLS: merged[CONTROLS]})
    _write(overview.reset_index(), out / "clinical" / "missing_data.tsv", index=False)


def step_figures(
    cfg: dict[str, Any],
    data: dict[str, Any],
    tables: dict[str, pd.DataFrame],
    panel: list[str],
    out: Path,
) -> None:
    """Volcano plots and the brain-enriched panel figure."""
    thr = cfg["thresholds"]
    fdir = out / "figures"
    titles = {c.name: c.title for c in de.COMPARISONS}
    for name in de.MAIN_COMPARISONS:
        figures.volcano(
            tables[name],
            panel,
            titles[name],
            fdir / f"figure_1_volcano_{name}",
            thr["fdr"],
            thr["abs_log2fc"],
        )
    figures.volcano(
        tables["hc_high_vs_low_prs"],
        panel,
        titles["hc_high_vs_low_prs"],
        fdir / "figure_s1_volcano_hc_high_vs_low_prs",
        thr["fdr"],
        thr["abs_log2fc"],
    )
    norm_path = out / "de" / "normalised_counts_all_samples.tsv"
    if norm_path.is_file():
        norm = pd.read_csv(norm_path, sep="\t").set_index("miRNA")
        pair_tables = {
            ("controls", "de novo PD"): tables["pd_vs_controls"],
            ("de novo PD", "iRBD"): tables["irbd_vs_pd"],
            ("controls", "iRBD"): tables["irbd_vs_controls"],
        }
        figures.panel_figure(
            norm, data["group3"], pair_tables, panel, fdir / "figure_2_brain_enriched_panel"
        )


def run(cfg: dict[str, Any], skip_de: bool = False, strict: bool = False) -> pd.DataFrame:
    """Run the whole pipeline. Return the key-numbers table."""
    out = Path(cfg["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    numbers: dict[str, Any] = {}
    data = step_inputs(cfg, out, numbers)
    tables = step_de(cfg, data, out, numbers, skip_de)
    cluster = step_cluster(cfg, tables, out, numbers)
    panel = step_panel(cfg, tables, cluster, out, numbers)
    step_haemolysis(cfg, data, tables, cluster, out, numbers)
    step_converters(cfg, data, panel, out, numbers)
    step_batch(cfg, data, out, numbers)
    step_clinical(cfg, data, out)
    step_figures(cfg, data, tables, panel, out)

    manuscript = None
    if cfg.get("manuscript_values"):
        manuscript = report.load_manuscript_values(cfg["manuscript_values"])
    table = report.build_table(numbers, manuscript)
    _write(table, out / "summary" / "key_numbers.tsv", index=False)
    info = {
        "pipeline_version": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "statsmodels": statsmodels.__version__,
        "inputs_sha256": {key: _sha256(path) for key, path in cfg["inputs"].items() if path},
        "config": cfg,
    }
    with (out / "summary" / "run_info.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(info, handle, sort_keys=False)
    status = report.count_status(table)
    if manuscript is not None:
        log.info(
            "Key numbers against the manuscript: %d match, %d differ, %d not reported",
            status["match"],
            status["differs"],
            status["not reported"],
        )
        if strict and status["differs"]:
            keys = table.loc[table["status"] == "differs", "key"].tolist()
            raise ManuscriptMismatch(
                f"{status['differs']} key numbers differ from the manuscript values: "
                f"{', '.join(keys[:8])}{' ...' if len(keys) > 8 else ''}. "
                f"See {out / 'summary' / 'key_numbers.tsv'}."
            )
    log.info("Finished. Results are in %s", out)
    return table
