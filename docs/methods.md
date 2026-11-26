# Methods and their implementation

Each statement of the Methods section, paraphrased, with the code that implements it, the setting
that controls it and the test that checks it.

## Sequencing data and quality control

| Methods statement | Implementation | Setting | Test |
|---|---|---|---|
| The analysis uses the libraries that passed quality control. | `io.analysis_set` | column `qc_pass` | `test_io.py` |
| Counts come from miRMaster 2.0. Non-integer counts are rounded, because DESeq2 models integer counts. | `io.round_counts`, `round_counts` in `deseq2_functions.R` | none | `test_io.py`, `test_deseq2_functions.R` |
| Number of annotated miRNAs, of miRNAs detected in at least one sample and in at least half of the samples. | `io.detection_summary` | none | `test_io.py` |
| Median and range of miRNA reads per sample. | `pipeline.step_inputs` | none | `test_pipeline_integration.py` |

## Differential expression

| Methods statement | Implementation | Setting | Test |
|---|---|---|---|
| Differential expression is analysed with DESeq2. Each comparison is fitted separately, on the samples of the two groups compared. | `fit_comparison` in `deseq2_functions.R`, `de.COMPARISONS` | `groups` | `test_deseq2_functions.R` |
| High-PRS and low-PRS controls are compared first. Controls are pooled for the other comparisons. | comparison `hc_high_vs_low_prs`, then pooled reference groups | `groups` | `test_pipeline_integration.py` |
| miRNAs with low counts are removed before fitting: edgeR `filterByExpr` for the two iRBD comparisons, at least 10 counts in at least 3 samples for de novo PD versus controls. | `prefilter_keep` | `deseq2.prefilter`, `deseq2.min_count`, `deseq2.min_samples` | `test_deseq2_functions.R`, `test_config.py` |
| Library size is normalised with median-of-ratios size factors. Dispersions are shrunk towards a fitted trend. | `DESeq2::DESeq` with default settings | none | `test_deseq2_functions.R` |
| The model includes study group, age and sex. Sequencing batch is not included. | `build_design`: `~ age + sex + group` | `deseq2.covariates` | `test_deseq2_functions.R`, `test_pipeline_integration.py` |
| Group differences are tested with the Wald test. | `DESeq2::results` | none | `test_deseq2_functions.R` |
| Log2 fold changes are shrunk with apeglm. | `DESeq2::lfcShrink(type = "apeglm")` | none | `test_deseq2_functions.R` |
| Independent filtering excludes miRNAs with low mean counts from the multiple-testing correction. P values are adjusted with the Benjamini-Hochberg procedure. | `DESeq2::results` | `deseq2.independent_filtering_alpha` | `test_deseq2_functions.R` |
| A miRNA is dysregulated when its adjusted p value is below 0.05 and its absolute log2 fold change is at least 0.5. The threshold applies to the apeglm-shrunken estimate. | `dysregulation.annotate` | `thresholds.fdr`, `thresholds.abs_log2fc` | `test_dysregulation.py` |

## Brain-enriched miRNAs

| Methods statement | Implementation | Setting | Test |
|---|---|---|---|
| Thirty predefined brain-enriched miRNAs are followed in every comparison. | `brain_enriched.panel_status`, `resources/brain_enriched_mirnas.tsv` | `resources.brain_enriched` | `test_brain_enriched.py` |
| Candidates are the panel miRNAs dysregulated in both iRBD comparisons, with concordant direction. | `brain_enriched.shared_candidates` | none | `test_brain_enriched.py` |
| Over-representation of the panel among dysregulated miRNAs is tested with a one-sided Fisher's exact test. | `brain_enriched.over_representation` | none | `test_brain_enriched.py` |

## 14q32 cluster

| Methods statement | Implementation | Setting | Test |
|---|---|---|---|
| The cluster is defined by genomic position, not by miRNA name. Members are all miRNA genes annotated in miRBase v22 (GRCh38) on the forward strand of chromosome 14 between DLK1 and DIO3, with their mature miRNAs. | `cluster14q32.extract_domain`, `cluster14q32.cluster_members`, `resources/mirbase_v22_chr14_DLK1_DIO3.tsv` | `resources.cluster_14q32` | `test_cluster14q32.py` |
| For each comparison, the tested and the dysregulated cluster miRNAs are counted. A miRNA is tested when it has an adjusted p value. | `cluster14q32.count_in_comparison` | none | `test_cluster14q32.py` |
| The gene on the reverse strand of the domain, hsa-mir-1247, is reported separately. | `cluster14q32.other_domain_members`, `pipeline.step_cluster` | none | `test_cluster14q32.py` |

## Haemolysis index and cluster score

| Methods statement | Implementation | Setting | Test |
|---|---|---|---|
| The haemolysis index is the log2 ratio of miR-451a to miR-23a-3p, using reads per million mapped miRNA reads plus one. | `haemolysis.haemolysis_index`, `io.rpmmm` | `haemolysis.numerator`, `haemolysis.denominator` | `test_haemolysis.py` |
| Groups are compared on the index with the two-sided Mann-Whitney U test. | `haemolysis.pairwise_mannwhitney` | none | `test_haemolysis.py` |
| The cluster score is the mean standardised expression of the cluster miRNAs tested in iRBD versus controls. | `haemolysis.cluster_score` | none | `test_haemolysis.py` |
| The relation of score and index is assessed with Spearman correlation within each group. | `haemolysis.spearman_by_group` | none | `test_haemolysis.py` |
| The score is regressed on study group, with controls as the reference. The haemolysis index, then age and sex, are added. | `haemolysis.score_models` | none | `test_haemolysis.py` |

## Converters within iRBD

| Methods statement | Implementation | Setting | Test |
|---|---|---|---|
| Counts are normalised with DESeq2 size factors. | `size_factors` in `deseq2_functions.R`, `converters.normalise` | none | `test_deseq2_functions.R`, `test_converters.py` |
| miRNAs detected in at least 20% of iRBD samples are tested. | `converters.detection_threshold` | `converters.min_detection_fraction` | `test_converters.py` |
| Groups are compared with the two-sided Mann-Whitney U test. | `converters.mannwhitney_rows`, `converters.compare` | none | `test_converters.py` |
| Fold changes are log2 of (mean converters + 1) divided by (mean non-converters + 1). | `converters.compare` | none | `test_converters.py` |
| Multiple-testing correction and thresholds are the same as for the main comparisons. | `converters.compare` | `thresholds` | `test_converters.py` |
| The full procedure is repeated on 200 datasets with randomly permuted group labels. | `converters.permutation_calibration` | `converters.n_permutations`, `converters.seed` | `test_converters.py` |

## Batch-effect analysis (Supplementary Figure S3)

| Methods statement | Implementation | Setting | Test |
|---|---|---|---|
| Balance of batch across study groups is tested with a chi-square test. Cramer's V is reported. | `batch.balance` | none | `test_batch.py` |
| Principal component analysis uses the log2 expression of the miRNAs detected in at least half of the samples. | `batch.expressed_log2`, `batch.principal_components` | `batch.min_detection_fraction` | `test_batch.py` |
| Each leading component is tested for association with batch and with study group by the Kruskal-Wallis test. | `batch.principal_components` | `batch.n_components` | `test_batch.py` |
| For each miRNA, a linear regression on group alone is compared with one that adds batch, using an F test. P values are adjusted with the Benjamini-Hochberg procedure. | `batch.per_mirna_batch_test` | none | `test_batch.py` |

## Clinical variables

| Methods statement | Implementation | Setting | Test |
|---|---|---|---|
| Continuous variables are compared with Welch's t test, categorical variables with Fisher's exact test. | `clinical.compare_groups` | `clinical.continuous`, `clinical.categorical` | `test_clinical.py` |
| Five comparisons, as in the columns of Table 1 and Supplementary Table S1: iRBD vs de novo PD, iRBD vs high-PRS controls, iRBD vs low-PRS controls, iRBD vs all controls, de novo PD vs all controls. | `pipeline.step_clinical` | `groups` | `test_pipeline_clinical.py` |
| P values are adjusted with the Benjamini-Hochberg procedure across all clinical variables within each pairwise group comparison. | `clinical.compare_groups` | none | `test_clinical.py` |
| Missing clinical values are handled by pairwise deletion. The number of missing values is listed per variable and group, with the columns of Supplementary Table S5. | `clinical.compare_groups`, `clinical.missing_data`, `pipeline.step_clinical` | none | `test_clinical.py`, `test_pipeline_clinical.py` |

## Software versions

The manuscript names R 4.6 with DESeq2 1.52 and apeglm 1.34 (Bioconductor 3.23), and Python 3.10.12
with SciPy 1.15.3 and statsmodels 0.14.6. The `Dockerfile` installs these versions and stops the
build when one of them differs.

## Implementation details that the Methods do not spell out

1. Independent filtering uses the DESeq2 default significance level of 0.1.
2. `filterByExpr` runs with its default settings and the two comparison groups. It is also used for
   the comparison of high-PRS and low-PRS controls.
3. Age is centred and scaled before fitting. This leaves the group coefficient unchanged.
4. Fold changes are always the first-named group relative to the second-named group.
5. Reads per million mapped miRNA reads are the counts divided by the total miRNA count of the
   library, times one million.
6. The cluster score standardises each miRNA across all analysed libraries, with the sample
   standard deviation.
7. In the converter comparison, a miRNA is tested when at least floor(0.2 x n) iRBD libraries have a
   non-zero count, with a lower bound of 3. For 42 libraries this is 8.
8. The Mann-Whitney U test uses the exact distribution when one group has at most 8 libraries and
   the miRNA has no tied values. Otherwise it uses the normal approximation.
9. The batch analysis uses the four sequencing groups, with high-PRS and low-PRS controls kept apart.
10. Figure 2 shows log2 of the DESeq2-normalised count plus one, with size factors estimated over
    all analysed libraries.

## Key numbers

`config/manuscript_values.yaml` lists every number reported in the manuscript that the pipeline
computes. A run writes each computed number next to the reported value in
`summary/key_numbers.tsv`, with the status `match`, `differs` or `not reported`.
