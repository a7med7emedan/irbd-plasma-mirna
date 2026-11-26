# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan
#
# DESeq2 differential expression for one pairwise comparison.
# Pre-filter, median-of-ratios normalisation, dispersion estimation, Wald test,
# independent filtering with Benjamini-Hochberg adjustment, apeglm shrinkage.

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(apeglm)
})

# Round a numeric count matrix to integers, as DESeq2 models integer counts.
round_counts <- function(counts) {
  m <- round(as.matrix(counts))
  storage.mode(m) <- "integer"
  m
}

# Rows to keep before model fitting.
#   "filterByExpr": edgeR::filterByExpr with the comparison groups.
#   "min_count":    at least `min_count` reads in at least `min_samples` libraries.
prefilter_keep <- function(counts, cond, rule, min_count = 10, min_samples = 3) {
  if (rule == "filterByExpr") {
    keep <- edgeR::filterByExpr(edgeR::DGEList(counts = counts), group = cond)
  } else if (rule == "min_count") {
    keep <- rowSums(counts >= min_count) >= min_samples
  } else {
    stop("Unknown pre-filter rule '", rule, "'. Use 'filterByExpr' or 'min_count'.")
  }
  as.logical(keep)
}

# Build the column data and the design formula: covariates first, group last.
# Numeric covariates are centred and scaled, which leaves the group coefficient unchanged.
# A covariate that is constant in the subset is dropped.
build_design <- function(samples, cond, covariates) {
  coldata <- data.frame(row.names = samples$library_id, cond = cond)
  terms <- character(0)
  for (cv in covariates) {
    if (!cv %in% colnames(samples)) {
      stop("Covariate '", cv, "' is not a column of the sample sheet.")
    }
    v <- samples[[cv]]
    if (anyNA(v)) {
      stop("Covariate '", cv, "' has missing values. Complete it or remove it from the model.")
    }
    if (is.numeric(v)) {
      if (stats::sd(v) == 0) next
      coldata[[cv]] <- (v - mean(v)) / stats::sd(v)
    } else {
      f <- factor(v)
      if (nlevels(f) < 2) next
      coldata[[cv]] <- f
    }
    terms <- c(terms, cv)
  }
  design <- stats::as.formula(paste("~", paste(c(terms, "cond"), collapse = " + ")))
  list(coldata = coldata, design = design)
}

# Fit one comparison. Fold changes are numerator relative to reference.
fit_comparison <- function(counts, samples, numerator, reference, prefilter,
                           covariates = c("age", "sex"), alpha = 0.1,
                           min_count = 10, min_samples = 3) {
  s <- samples[samples$group %in% c(numerator, reference), , drop = FALSE]
  if (!any(s$group %in% numerator) || !any(s$group %in% reference)) {
    stop("A comparison group has no libraries. Check the group labels in the sample sheet.")
  }
  cond <- factor(ifelse(s$group %in% numerator, "numerator", "reference"),
                 levels = c("reference", "numerator"))
  m <- counts[, s$library_id, drop = FALSE]
  n_input <- nrow(m)
  keep <- prefilter_keep(m, cond, prefilter, min_count, min_samples)
  m <- m[keep, , drop = FALSE]
  if (nrow(m) < 2) {
    stop("Fewer than two miRNAs passed the pre-filter. The comparison cannot be fitted.")
  }
  d <- build_design(s, cond, covariates)
  dds <- DESeq2::DESeqDataSetFromMatrix(m, d$coldata, d$design)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  coef <- "cond_numerator_vs_reference"
  stopifnot(coef %in% DESeq2::resultsNames(dds))
  res <- DESeq2::results(dds, name = coef, alpha = alpha)
  shr <- DESeq2::lfcShrink(dds, coef = coef, type = "apeglm", res = res, quiet = TRUE)
  table <- data.frame(
    miRNA = rownames(res),
    baseMean = res$baseMean,
    log2FoldChange = shr$log2FoldChange,
    lfcSE = shr$lfcSE,
    log2FoldChange_mle = res$log2FoldChange,
    lfcSE_mle = res$lfcSE,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj,
    stringsAsFactors = FALSE
  )
  table <- table[order(table$padj, table$pvalue, table$miRNA, na.last = TRUE), ]
  rownames(table) <- NULL
  list(
    table = table,
    n_numerator = sum(cond == "numerator"),
    n_reference = sum(cond == "reference"),
    n_input = n_input,
    n_after_prefilter = nrow(m),
    n_with_padj = sum(!is.na(table$padj)),
    design = paste(deparse(d$design), collapse = "")
  )
}

# DESeq2 median-of-ratios size factors for a set of libraries.
size_factors <- function(counts) {
  coldata <- data.frame(row.names = colnames(counts))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, ~1)
  DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds))
}
