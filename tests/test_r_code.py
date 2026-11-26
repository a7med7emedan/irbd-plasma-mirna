# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
"""Run the testthat suite of the R code from pytest."""

import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration
TESTS = Path(__file__).resolve().parent / "testthat"


def test_testthat_suite(rscript):
    code = f'testthat::test_dir("{TESTS.as_posix()}", stop_on_failure = TRUE, reporter = "summary")'
    proc = subprocess.run([rscript, "-e", code], capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stdout + proc.stderr
