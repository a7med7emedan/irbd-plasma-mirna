# Outputs

All paths are relative to `output_dir`. Tables are tab-separated.

## qc

| File | Content |
|---|---|
| `library_depth.tsv` | miRNA reads per analysed library |

## de

| File | Content |
|---|---|
| `<comparison>_all.tsv` | DESeq2 output for every fitted miRNA |
| `<comparison>_annotated.tsv` | The same table with the columns `tested`, `dysregulated` and `direction` |
| `<comparison>_significant.tsv` | Dysregulated miRNAs only |
| `comparison_summary.tsv` | Groups, library numbers, pre-filter, design and miRNA numbers per comparison |
| `size_factors_all_samples.tsv`, `normalised_counts_all_samples.tsv` | DESeq2 size factors and normalised counts over all analysed libraries |
| `size_factors_within_group.tsv` | DESeq2 size factors within iRBD |
| `R_session_info.txt` | R and package versions of the run |

Comparison names are `hc_high_vs_low_prs`, `irbd_vs_controls`, `irbd_vs_pd` and `pd_vs_controls`.

Columns of the DESeq2 tables:

| Column | Content |
|---|---|
| `baseMean` | Mean of normalised counts |
| `log2FoldChange`, `lfcSE` | apeglm-shrunken estimate and its posterior standard deviation |
| `log2FoldChange_mle`, `lfcSE_mle` | Maximum-likelihood estimate and its standard error |
| `stat`, `pvalue` | Wald statistic and p value |
| `padj` | Benjamini-Hochberg adjusted p value. Empty when independent filtering removed the miRNA. |

## brain_enriched

| File | Content |
|---|---|
| `panel_<comparison>.tsv` | Status, estimate and call of each panel miRNA |
| `shared_candidates.tsv` | Panel miRNAs dysregulated in both iRBD comparisons |
| `over_representation.tsv` | One-sided Fisher's exact test, with and without the 14q32 cluster |

## cluster_14q32

| File | Content |
|---|---|
| `cluster_counts.tsv` | Tested and dysregulated cluster miRNAs per comparison |
| `other_domain_mirnas.tsv` | miRNAs of the reverse-strand gene of the domain |

## haemolysis

| File | Content |
|---|---|
| `per_library.tsv` | Haemolysis index and cluster score per library |
| `index_group_tests.tsv` | Mann-Whitney U tests of the index between groups |
| `spearman_index_vs_score.tsv` | Spearman correlation of index and score within each group |
| `score_models.tsv` | Linear regressions of the score on group, with adjustment |
| `score_by_prs_stratum.tsv` | Score of iRBD against each control stratum |
| `spearman_age_vs_score_controls.tsv` | Spearman correlation of age and score in controls |

## converters

| File | Content |
|---|---|
| `converters_vs_non_converters.tsv` | Test result per miRNA |
| `permutation_calibration.tsv` | Number of permuted datasets with at least one call |

## batch

| File | Content |
|---|---|
| `batch_by_group.tsv` | Cross table of batch and group |
| `pca_scores.tsv`, `pca_components.tsv` | Component scores, explained variance and Kruskal-Wallis p values |
| `per_mirna_batch_test.tsv` | F test p value and adjusted p value per miRNA |

## clinical

| File | Content |
|---|---|
| `comparison_<groups>.tsv` | Descriptive values, p value and adjusted p value per variable |
| `missing_data.tsv` | Number and percentage of missing values per variable and group |

The five comparisons follow the columns of Table 1 and Supplementary Table S1: iRBD vs de novo PD,
iRBD vs high-PRS controls, iRBD vs low-PRS controls, iRBD vs controls and de novo PD vs controls.
`missing_data.tsv` has the columns of Supplementary Table S5: iRBD, de novo PD, high-PRS controls,
low-PRS controls and all controls.

## figures

| File | Manuscript item |
|---|---|
| `figure_1_volcano_<comparison>` | Figure 1 |
| `figure_2_brain_enriched_panel` | Figure 2 |
| `figure_s1_volcano_hc_high_vs_low_prs` | Supplementary Figure S1 |
| `figure_s2_converters` | Supplementary Figure S2 |
| `figure_s3_batch_pca` | Supplementary Figure S3 |

Each figure is written as PNG and PDF.

## summary

| File | Content |
|---|---|
| `key_numbers.tsv` | Every key number of the run, the reported value and the status |
| `run_info.yaml` | Software versions, checksums of the inputs and the configuration used |

## work

Intermediate files passed to R: the analysis set of the counts, the sample sheet and the job file.
