# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
from conftest import annotated_table
from irbdmirna import cluster14q32
from irbdmirna.config import packaged_resource

GFF = """##gff-version 3
chr14\t.\tmiRNA_primary_transcript\t100814491\t100814574\t.\t+\t.\tID=MI1;Alias=MI1;Name=hsa-mir-a
chr14\t.\tmiRNA\t100814500\t100814520\t.\t+\t.\tID=MIMAT1;Alias=MIMAT1;Name=hsa-miR-a-5p;Derives_from=MI1
chr14\t.\tmiRNA\t100814540\t100814560\t.\t+\t.\tID=MIMAT2;Alias=MIMAT2;Name=hsa-miR-a-3p;Derives_from=MI1
chr14\t.\tmiRNA_primary_transcript\t101560287\t101560422\t.\t-\t.\tID=MI2;Alias=MI2;Name=hsa-mir-b
chr14\t.\tmiRNA\t101560300\t101560320\t.\t-\t.\tID=MIMAT3;Alias=MIMAT3;Name=hsa-miR-b-5p;Derives_from=MI2
chr14\t.\tmiRNA_primary_transcript\t100307859\t100307956\t.\t+\t.\tID=MI3;Alias=MI3;Name=hsa-mir-upstream
chr14\t.\tmiRNA\t100307870\t100307890\t.\t+\t.\tID=MIMAT4;Alias=MIMAT4;Name=hsa-miR-upstream-5p;Derives_from=MI3
chr1\t.\tmiRNA_primary_transcript\t100814491\t100814574\t.\t+\t.\tID=MI4;Alias=MI4;Name=hsa-mir-other
chr1\t.\tmiRNA\t100814500\t100814520\t.\t+\t.\tID=MIMAT5;Alias=MIMAT5;Name=hsa-miR-other-5p;Derives_from=MI4
"""


def test_extract_domain_uses_position_and_strand(tmp_path):
    path = tmp_path / "toy.gff3"
    path.write_text(GFF, encoding="utf-8")
    table = cluster14q32.extract_domain(path)
    assert set(table["miRNA"]) == {"hsa-miR-a-5p", "hsa-miR-a-3p", "hsa-miR-b-5p"}
    assert cluster14q32.cluster_members(table) == {"hsa-miR-a-5p", "hsa-miR-a-3p"}
    assert cluster14q32.other_domain_members(table) == {"hsa-miR-b-5p"}
    assert cluster14q32.describe(table) == {
        "genes_in_domain": 2,
        "genes_forward_strand": 1,
        "mature_forward_strand": 2,
        "mature_in_domain": 3,
    }


def test_packaged_resource_matches_mirbase_v22():
    resource = cluster14q32.read_resource(packaged_resource("mirbase_v22_chr14_DLK1_DIO3.tsv"))
    d = cluster14q32.describe(resource)
    assert d["genes_in_domain"] == 54
    assert d["genes_forward_strand"] == 53
    assert d["mature_forward_strand"] == 93
    forward = resource[resource["strand"] == "+"]
    assert forward["gene_start"].min() == 100814491  # hsa-mir-2392
    assert forward["gene_end"].max() == 101066801  # hsa-mir-656
    assert cluster14q32.other_domain_members(resource) == {"hsa-miR-1247-5p", "hsa-miR-1247-3p"}
    members = cluster14q32.cluster_members(resource)
    twelve = "127-3p 136-3p 154-5p 323a-3p 382-5p 409-3p 410-3p 432-5p 485-5p 495-3p 487b-3p 654-3p"
    assert {f"hsa-miR-{m}" for m in twelve.split()} <= members
    assert not {"hsa-miR-7-5p", "hsa-miR-3200-3p"} & members


def test_count_in_comparison_counts_only_tested_members():
    table = annotated_table(
        {
            "c1": (-1.2, 0.001),
            "c2": (-0.9, 0.2),
            "c3": (-2.0, float("nan")),
            "c4": (0.8, 0.01),
            "x": (-3.0, 0.001),
        }
    )
    counts = cluster14q32.count_in_comparison(table, {"c1", "c2", "c3", "c4", "absent"})
    assert counts == {"tested": 3, "dysregulated": 2, "down": 1, "up": 1, "negative_log2fc": 2}
