# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import pandas as pd
import pytest
import yaml

from irbdmirna import __version__, cli


def test_version(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--version"])
    assert exit_info.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_make_example_data(tmp_path, capsys):
    assert cli.main(["make-example-data", "--outdir", str(tmp_path / "data"), "--seed", "1"]) == 0
    assert (tmp_path / "data" / "counts.tsv").is_file()
    assert "samples:" in capsys.readouterr().out


def test_build_cluster_resource(tmp_path):
    gff = tmp_path / "toy.gff3"
    gff.write_text(
        "chr14\t.\tmiRNA_primary_transcript\t100900000\t100900080\t.\t+\t.\tID=MI1;Name=hsa-mir-a\n"
        "chr14\t.\tmiRNA\t100900010\t100900030\t.\t+\t.\tID=M1;Name=hsa-miR-a-5p;Derives_from=MI1\n",
        encoding="utf-8",
    )
    out = tmp_path / "res" / "domain.tsv"
    assert cli.main(["build-cluster-resource", "--gff", str(gff), "--out", str(out)]) == 0
    assert pd.read_csv(out, sep="\t")["miRNA"].tolist() == ["hsa-miR-a-5p"]


def test_run_reports_input_problems_without_a_traceback(tmp_path, caplog):
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        yaml.safe_dump({"inputs": {"counts": "absent.tsv", "samples": "absent.tsv"}}),
        encoding="utf-8",
    )
    assert cli.main(["run", "--config", str(cfg), "--output-dir", str(tmp_path / "out")]) == 1
    assert "Count table not found" in caplog.text
