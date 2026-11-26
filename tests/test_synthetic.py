# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
import pandas as pd

from irbdmirna import io, synthetic


def test_dataset_is_deterministic_and_well_formed():
    first, second = synthetic.make_dataset(seed=3), synthetic.make_dataset(seed=3)
    pd.testing.assert_frame_equal(first["counts"], second["counts"])
    pd.testing.assert_frame_equal(first["samples"], second["samples"])
    samples, counts = first["samples"], first["counts"]
    assert list(counts.columns) == samples["library_id"].tolist()
    assert set(io.REQUIRED_SAMPLE_COLUMNS) <= set(samples.columns)
    assert {"hsa-miR-451a", "hsa-miR-23a-3p", "hsa-miR-1247-5p"} <= set(counts.index)
    assert (~samples["qc_pass"]).sum() == synthetic.N_FAILED_QC
    assert samples.loc[samples["converter"], "group"].eq("iRBD").all()
    _, share = io.round_counts(counts)
    assert 0 < share < 0.02


def test_write_dataset(tmp_path):
    paths = synthetic.write_dataset(tmp_path, seed=3)
    assert all(p.is_file() for p in paths.values())
    counts = io.read_counts(paths["counts"])
    assert counts.shape[1] == sum(synthetic.GROUP_SIZES.values()) + synthetic.N_FAILED_QC
