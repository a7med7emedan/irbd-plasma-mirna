# Contributing

## Report a problem

1. Open an issue and choose the bug report template.
2. Give the command, the last 30 log lines and the file `summary/run_info.yaml`.
3. Never attach participant data. Reproduce the problem with the synthetic example when possible.

## Propose a change

1. Open an issue first when the change alters a statistical method or a default value. The
   defaults follow the published Methods and are changed only with a documented reason.
2. Fork the repository and create a branch.
3. Set up the environment with `micromamba create -f environment.yml`, or build the Docker image.
4. Make the change and add a test for it.
5. Run `make lint` and `make test`. Both must pass.
6. Update `CHANGELOG.md` and the documentation.
7. Open a pull request and fill in the template.

## Style

1. Python code is formatted and linted with ruff. Run `make format`.
2. Every source file starts with the SPDX licence identifier and the copyright line.
3. Error messages say what the user should change.
4. Documentation uses plain sentences and numbered steps.

## Tests

1. Unit tests cover one function each and build their own data.
2. `tests/test_pipeline_integration.py` runs the whole pipeline on synthetic data.
3. `tests/testthat` holds the tests of the R code.
