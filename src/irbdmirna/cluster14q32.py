# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""14q32 (DLK1-DIO3) miRNA cluster, defined by genomic position."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# GRCh38 gene spans (Ensembl): DLK1 chr14:100,725,705-100,738,224, DIO3 chr14:101,560,139-101,563,452.
DOMAIN_CHROMOSOME = "chr14"
DOMAIN_START = 100_725_705
DOMAIN_END = 101_563_452
RESOURCE_COLUMNS = ["miRNA", "gene", "chromosome", "gene_start", "gene_end", "strand"]


def extract_domain(
    gff_path: str | Path,
    chromosome: str = DOMAIN_CHROMOSOME,
    start: int = DOMAIN_START,
    end: int = DOMAIN_END,
) -> pd.DataFrame:
    """List the mature miRNAs of all miRNA genes inside the domain, from a miRBase GFF3 file."""
    genes: dict[str, dict] = {}
    mature: list[tuple[str, str]] = []
    with Path(gff_path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                raise ValueError(f"Malformed GFF3 line in {gff_path}: {line[:60]!r}")
            attrs = dict(item.split("=", 1) for item in f[8].split(";") if "=" in item)
            if f[2] == "miRNA_primary_transcript":
                if f[0] == chromosome and int(f[3]) >= start and int(f[4]) <= end:
                    genes[attrs["ID"]] = {
                        "gene": attrs["Name"],
                        "chromosome": f[0],
                        "gene_start": int(f[3]),
                        "gene_end": int(f[4]),
                        "strand": f[6],
                    }
            elif f[2] == "miRNA":
                mature.append((attrs["Name"], attrs["Derives_from"]))
    rows = [{"miRNA": name, **genes[parent]} for name, parent in mature if parent in genes]
    table = pd.DataFrame(rows, columns=RESOURCE_COLUMNS)
    return table.sort_values(["gene_start", "miRNA"]).reset_index(drop=True)


def read_resource(path: str | Path) -> pd.DataFrame:
    """Read the domain table produced by `extract_domain`."""
    table = pd.read_csv(path, sep="\t", comment="#")
    missing = [c for c in RESOURCE_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(f"The 14q32 resource {path} lacks the columns {missing}.")
    return table


def cluster_members(resource: pd.DataFrame) -> set[str]:
    """Mature miRNAs of forward-strand genes: the cluster."""
    return set(resource.loc[resource["strand"] == "+", "miRNA"])


def other_domain_members(resource: pd.DataFrame) -> set[str]:
    """Mature miRNAs of reverse-strand genes in the domain, reported separately."""
    return set(resource.loc[resource["strand"] != "+", "miRNA"]) - cluster_members(resource)


def describe(resource: pd.DataFrame) -> dict[str, int]:
    """Gene and mature miRNA counts of the domain."""
    fwd = resource[resource["strand"] == "+"]
    return {
        "genes_in_domain": int(resource["gene"].nunique()),
        "genes_forward_strand": int(fwd["gene"].nunique()),
        "mature_forward_strand": int(fwd["miRNA"].nunique()),
        "mature_in_domain": int(resource["miRNA"].nunique()),
    }


def count_in_comparison(table: pd.DataFrame, members: set[str]) -> dict[str, int]:
    """Tested and dysregulated cluster miRNAs in one annotated DESeq2 table."""
    sub = table[table.index.isin(members) & table["tested"]]
    return {
        "tested": int(len(sub)),
        "dysregulated": int(sub["dysregulated"].sum()),
        "down": int((sub["direction"] == "down").sum()),
        "up": int((sub["direction"] == "up").sum()),
        "negative_log2fc": int((sub["log2FoldChange"] < 0).sum()),
    }
