#!/usr/bin/env Rscript
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
#
# Run all DESeq2 comparisons described in a job file.
# Usage: Rscript run_deseq2.R <job.yaml>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_deseq2.R <job.yaml>")
}
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
source(file.path(script_dir, "deseq2_functions.R"))

job <- yaml::read_yaml(args[1])
set.seed(1)
dir.create(job$outdir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(job$counts, row.names = 1, check.names = FALSE)
samples <- read.delim(job$samples, stringsAsFactors = FALSE, check.names = FALSE,
                      colClasses = c(library_id = "character", group = "character"))
stopifnot(!anyDuplicated(samples$library_id), all(samples$library_id %in% colnames(counts)))
counts <- round_counts(counts[, samples$library_id, drop = FALSE])

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

summary_rows <- list()
for (cmp in job$comparisons) {
  fit <- fit_comparison(
    counts, samples,
    numerator = unlist(cmp$numerator), reference = unlist(cmp$reference),
    prefilter = cmp$prefilter, covariates = unlist(job$covariates),
    alpha = job$alpha, min_count = job$min_count, min_samples = job$min_samples
  )
  write_tsv(fit$table, file.path(job$outdir, paste0(cmp$name, "_all.tsv")))
  summary_rows[[cmp$name]] <- data.frame(
    comparison = cmp$name,
    numerator = paste(unlist(cmp$numerator), collapse = "+"),
    reference = paste(unlist(cmp$reference), collapse = "+"),
    n_numerator = fit$n_numerator, n_reference = fit$n_reference,
    prefilter = cmp$prefilter, design = fit$design,
    n_input = fit$n_input, n_after_prefilter = fit$n_after_prefilter,
    n_with_padj = fit$n_with_padj, stringsAsFactors = FALSE
  )
  message(sprintf("%-22s %d vs %d libraries, %d miRNAs fitted, %d with adjusted p",
                  cmp$name, fit$n_numerator, fit$n_reference,
                  fit$n_after_prefilter, fit$n_with_padj))
}
write_tsv(do.call(rbind, summary_rows), file.path(job$outdir, "comparison_summary.tsv"))

# Size factors and normalised counts over all analysed libraries (used for display).
sf_all <- size_factors(counts)
write_tsv(data.frame(library_id = names(sf_all), size_factor = sf_all),
          file.path(job$outdir, "size_factors_all_samples.tsv"))
norm <- sweep(counts, 2, sf_all, "/")
write_tsv(data.frame(miRNA = rownames(norm), norm, check.names = FALSE),
          file.path(job$outdir, "normalised_counts_all_samples.tsv"))

# Size factors within one group (used for the converter comparison).
if (!is.null(job$within_group)) {
  libs <- samples$library_id[samples$group == job$within_group]
  if (length(libs) >= 2) {
    sf_g <- size_factors(counts[, libs, drop = FALSE])
    write_tsv(data.frame(library_id = names(sf_g), size_factor = sf_g),
              file.path(job$outdir, "size_factors_within_group.tsv"))
  }
}

writeLines(capture.output(sessionInfo()), file.path(job$outdir, "R_session_info.txt"))
