# Usage and configuration

## Commands

| Command | Purpose |
|---|---|
| `irbd-mirna run --config FILE` | Run the whole analysis |
| `irbd-mirna make-example-data --outdir DIR` | Write the synthetic example dataset |
| `irbd-mirna build-cluster-resource --gff FILE --out FILE` | Derive the 14q32 table from a miRBase GFF3 file |
| `irbd-mirna --version` | Show the version |

Options of `run`:

| Option | Meaning |
|---|---|
| `--output-dir DIR` | Write results to DIR and ignore `output_dir` of the configuration |
| `--skip-de` | Reuse the DESeq2 tables already in the output directory |
| `--strict` | Exit with code 1 when a key number differs from the manuscript value |
| `-v` | Show debug messages. Place it before `run`. |

The exit code is 0 on success and 1 on an error. Error messages say what to change.

## Configuration file

Paths are relative to the directory from which the pipeline is started. Unknown keys are rejected,
so that a misspelled setting cannot be ignored without notice.

```yaml
inputs:
  counts: data/counts.tsv
  samples: data/samples.tsv
  clinical: null                 # optional

output_dir: results/study
manuscript_values: config/manuscript_values.yaml   # optional
rscript: Rscript                 # R executable

groups:                          # labels used in the sample sheet
  irbd: iRBD
  pd: de_novo_PD
  hc_high_prs: HC_high_PRS
  hc_low_prs: HC_low_PRS

thresholds:
  fdr: 0.05
  abs_log2fc: 0.5

deseq2:
  covariates: [age, sex]
  independent_filtering_alpha: 0.1
  min_count: 10
  min_samples: 3
  prefilter:
    hc_high_vs_low_prs: filterByExpr
    irbd_vs_controls: filterByExpr
    irbd_vs_pd: filterByExpr
    pd_vs_controls: min_count

haemolysis:
  numerator: hsa-miR-451a
  denominator: hsa-miR-23a-3p

converters:
  min_detection_fraction: 0.2
  min_detected_samples: 3
  n_permutations: 200
  seed: 20260830

batch:
  min_detection_fraction: 0.5
  n_components: 3

clinical:
  continuous: []
  categorical: []

resources:
  brain_enriched: null           # null uses the table shipped with the package
  cluster_14q32: null
```

## Comparing a run with the manuscript

1. Set `manuscript_values` to `config/manuscript_values.yaml`.
2. Run with `--strict`.
3. Open `summary/key_numbers.tsv`. Each row holds the pipeline value, the reported value, the
   comparison rule and the status.

Comparison rules:

| Rule | A value matches when |
|---|---|
| `exact` | it equals the reported value |
| `round` | it equals the reported value after rounding to the stated digits |
| `below` | it is smaller than the reported bound, for values reported as p < bound |
| `text` | it equals the reported text, for example `not tested` |

## Steps that are skipped

1. The converter comparison is skipped when fewer than two converters or two non-converters are
   flagged in the sample sheet.
2. The batch analysis is skipped when all analysed libraries share one batch.
3. The clinical comparisons are skipped when no clinical table is configured.

The log names each skipped step.
