# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan

PYTHON ?= python
CONFIG ?= config/config.example.yaml
IMAGE  ?= irbd-plasma-mirna:1.0.0

.PHONY: help check-python install lint format test test-fast test-r example-data run \
        docker-build docker-test docker-run doc-figures clean

help:
	@echo "make install       install the package with the test tools"
	@echo "make lint          run the linter and the format check"
	@echo "make test          run all tests (needs R with DESeq2, edgeR and apeglm)"
	@echo "make test-fast     run the tests that do not need R"
	@echo "make test-r        run the R unit tests"
	@echo "make example-data  rebuild the synthetic example dataset"
	@echo "make run           run the pipeline (CONFIG=$(CONFIG))"
	@echo "make docker-build  build the reference image ($(IMAGE))"
	@echo "make docker-test   run lint and all tests inside the image"
	@echo "make docker-run    run the pipeline on the example data inside the image"
	@echo "make doc-figures    redraw the figures in docs/img (needs cairosvg)"
	@echo "make clean         remove results and caches"

check-python:
	@$(PYTHON) -c "import sys; sys.exit('Python 3.10 or newer is required, found ' + sys.version.split()[0]) if sys.version_info < (3, 10) else None"

install: check-python
	$(PYTHON) -m pip install -r requirements.txt -r requirements-dev.txt
	$(PYTHON) -m pip install --no-deps -e .

lint:
	ruff check src tests
	ruff format --check src tests

format:
	ruff check --fix src tests
	ruff format src tests

test: check-python
	$(PYTHON) -m pytest

test-fast: check-python
	$(PYTHON) -m pytest -m "not integration"

test-r:
	Rscript -e 'testthat::test_dir("tests/testthat", stop_on_failure = TRUE)'

example-data:
	irbd-mirna make-example-data --outdir examples/data

run:
	irbd-mirna run --config $(CONFIG)

docker-build:
	docker build -t $(IMAGE) .

docker-test:
	docker run --rm --entrypoint make $(IMAGE) lint test

docker-run:
	mkdir -p results
	docker run --rm -v "$(CURDIR)/results:/opt/irbd-plasma-mirna/results" $(IMAGE)

doc-figures:
	$(PYTHON) tools/make_doc_figures.py

clean:
	rm -rf results .pytest_cache .ruff_cache build dist src/*.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
