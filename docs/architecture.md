# Architecture

## Components

```mermaid
flowchart TB
    subgraph cli [Command line]
        A[cli.py]
    end
    subgraph core [Python package irbdmirna]
        B[config.py]
        C[io.py]
        D[pipeline.py]
        E[de.py]
        F[dysregulation.py]
        G[brain_enriched.py]
        H[cluster14q32.py]
        I[haemolysis.py]
        J[converters.py]
        K[batch.py]
        L[clinical.py]
        M[figures.py]
        N[report.py]
    end
    subgraph rcode [R code]
        R1[run_deseq2.R]
        R2[deseq2_functions.R]
    end
    subgraph res [Resources]
        P[brain_enriched_mirnas.tsv]
        Q[mirbase_v22_chr14_DLK1_DIO3.tsv]
    end
    A --> B
    A --> D
    D --> C
    D --> E
    E --> R1
    R1 --> R2
    D --> F
    D --> G
    D --> H
    D --> I
    D --> J
    D --> K
    D --> L
    D --> M
    D --> N
    G --> P
    H --> Q
```

## Responsibilities

| Module | Responsibility |
|---|---|
| `cli.py` | Parse the command line, set up logging, turn errors into exit code 1 |
| `config.py` | Defaults, YAML loading, validation |
| `io.py` | Read and validate the inputs, round counts, compute reads per million |
| `pipeline.py` | Run the steps in order, write tables, collect the key numbers |
| `de.py` | Define the four comparisons, write the job file, call Rscript |
| `r/deseq2_functions.R` | Pre-filter, design, DESeq2 fit, Wald test, apeglm shrinkage, size factors |
| `r/run_deseq2.R` | Read the job file, run all comparisons, write the tables |
| `dysregulation.py` | Apply the dysregulation rule |
| `brain_enriched.py` | Panel status, shared candidates, Fisher's exact test |
| `cluster14q32.py` | Derive and read the 14q32 table, count cluster miRNAs |
| `haemolysis.py` | Haemolysis index, cluster score, correlations, regressions |
| `converters.py` | Converter comparison and permutation calibration |
| `batch.py` | Batch balance, principal components, per-miRNA F test |
| `clinical.py` | Clinical comparisons and missing-data overview |
| `figures.py` | All figures |
| `report.py` | Key numbers and their comparison with the manuscript values |
| `synthetic.py` | Synthetic example data |

## Design choices

1. DESeq2 runs in R and everything else in Python, as in the manuscript. The two parts exchange
   plain tab-separated files and one YAML job file.
2. Each comparison is fitted on its own subset of libraries. The R script receives the list of
   comparisons and has no study-specific logic.
3. All thresholds live in the configuration. The defaults are the values of the manuscript.
4. Random steps use a fixed seed, so two runs on the same input give the same numbers.
5. Every run records software versions and input checksums in `summary/run_info.yaml`.

## Deployment

```mermaid
flowchart LR
    subgraph host [Host]
        D1[Data directory, read only]
        D2[Results directory]
        D3[Configuration file]
    end
    subgraph image [Docker image]
        E1[R 4.6, DESeq2 1.52, apeglm 1.34, edgeR]
        E2[Python 3.10.12, SciPy 1.15.3, statsmodels 0.14.6]
        E3[irbd-mirna]
    end
    D1 --> E3
    D3 --> E3
    E3 --> D2
    E3 --> E1
    E3 --> E2
```
