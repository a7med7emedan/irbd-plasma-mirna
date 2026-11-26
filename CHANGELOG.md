# Changelog

All notable changes to this project are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-11-26

### Added

1. DESeq2 analysis of four comparisons with pre-filter, age and sex adjustment, Wald test,
   independent filtering and apeglm shrinkage.
2. Dysregulation rule, brain-enriched panel summary and one-sided Fisher's exact test.
3. 14q32 cluster defined by genomic position from miRBase release 22.
4. Haemolysis index, cluster score, correlations and adjusted regression models.
5. Converter comparison within iRBD with permutation calibration.
6. Batch-effect analysis and clinical group comparisons.
7. Figures 1 and 2 and Supplementary Figures S1 to S3.
8. Key numbers report with comparison against the values reported in the manuscript.
9. Synthetic example dataset, unit tests, end-to-end test, Docker image and continuous integration.
10. Workflow figure and run figure, drawn by `tools/make_doc_figures.py`, with a text version of
    the workflow in `docs/workflow.md`.
11. Separate Tests, Conda and Docker workflows, each with its own status badge.
12. Clinical comparisons of iRBD against high-PRS and low-PRS controls, and missing data per
    control stratum, as in Table 1, Supplementary Table S1 and Supplementary Table S5.
