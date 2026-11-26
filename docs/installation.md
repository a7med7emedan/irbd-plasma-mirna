# Installation

## Docker

1. Install Docker.
2. Clone the repository and build the image.

```bash
git clone https://github.com/a7med7emedan/irbd-plasma-mirna.git
cd irbd-plasma-mirna
make docker-build
```

3. Check the image.

```bash
make docker-test     # linter and all tests inside the image
make docker-run      # pipeline on the example data, results in ./results/example
```

The build stops when DESeq2 is not 1.52, apeglm is not 1.34, Python is not 3.10.12, SciPy is not
1.15.3 or statsmodels is not 0.14.6. These are the versions stated in the manuscript.

Files written to a mounted directory belong to root. Add `--user "$(id -u):$(id -g)"` to
`docker run` to write them as your own user.

## Conda

```bash
micromamba create -f environment.yml
micromamba activate irbd-plasma-mirna
make test
make run
```

`conda` and `mamba` accept the same file with `env create -f environment.yml`. The environment holds
R 4.5 with DESeq2 1.50.2 and apeglm 1.32.0, the nearest builds on bioconda, and the Python packages
of `requirements.txt`.

## Manual

1. Install R 4.5 or newer, then the R packages.

```r
install.packages(c("BiocManager", "yaml", "testthat"))
BiocManager::install(c("DESeq2", "edgeR", "apeglm"))
```

2. Install Python 3.10 or newer, then the package.

```bash
python -m venv .venv
. .venv/bin/activate
make install
```

3. Check the installation.

```bash
irbd-mirna --version
make test
```

`make install` stops with a clear message when Python is older than 3.10.
