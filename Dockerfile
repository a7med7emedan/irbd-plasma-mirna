# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
#
# Reference environment with the software versions stated in the manuscript:
# R 4.6, DESeq2 1.52 and apeglm 1.34 (Bioconductor 3.23), Python 3.10.12,
# SciPy 1.15.3 and statsmodels 0.14.6.

FROM ghcr.io/astral-sh/uv:0.13.10 AS uv

FROM bioconductor/bioconductor_docker:RELEASE_3_23-R-4.6.1

LABEL org.opencontainers.image.title="irbd-plasma-mirna" \
      org.opencontainers.image.description="Plasma microRNA analysis for iRBD, de novo PD and healthy controls" \
      org.opencontainers.image.authors="Ahmed Hemedan" \
      org.opencontainers.image.licenses="Apache-2.0" \
      org.opencontainers.image.version="1.0.0"

# R packages from Bioconductor 3.23, then a check of the versions.
RUN Rscript -e 'BiocManager::install(c("DESeq2", "edgeR", "apeglm"), ask = FALSE, update = FALSE)' \
 && Rscript -e 'install.packages(c("yaml", "testthat"))' \
 && Rscript -e 'v <- function(p) as.character(packageVersion(p)); \
      stopifnot(startsWith(v("DESeq2"), "1.52."), startsWith(v("apeglm"), "1.34."), \
                requireNamespace("edgeR"), requireNamespace("yaml"), requireNamespace("testthat")); \
      cat("DESeq2", v("DESeq2"), "apeglm", v("apeglm"), "edgeR", v("edgeR"), "\n")'

# Python 3.10.12 in its own environment.
COPY --from=uv /uv /usr/local/bin/uv
ENV UV_PYTHON_INSTALL_DIR=/opt/python \
    UV_LINK_MODE=copy \
    VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:$PATH \
    MPLBACKEND=Agg
RUN uv python install 3.10.12 && uv venv --python 3.10.12 /opt/venv

WORKDIR /opt/irbd-plasma-mirna
COPY requirements.txt requirements-dev.txt ./
RUN uv pip install -r requirements.txt -r requirements-dev.txt

COPY . .
RUN uv pip install --no-deps -e . \
 && python -c "import scipy, statsmodels, sys; \
      assert sys.version_info[:3] == (3, 10, 12), sys.version; \
      assert scipy.__version__ == '1.15.3' and statsmodels.__version__ == '0.14.6'" \
 && irbd-mirna --version

ENTRYPOINT ["irbd-mirna"]
CMD ["run", "--config", "config/config.example.yaml"]
