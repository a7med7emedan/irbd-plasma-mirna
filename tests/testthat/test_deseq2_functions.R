# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Ahmed Hemedan

r_dir <- Sys.getenv("IRBD_R_DIR", file.path("..", "..", "src", "irbdmirna", "r"))
source(file.path(r_dir, "deseq2_functions.R"))

simulate <- function(n_ref = 12, n_num = 12, n_mirna = 200, n_up = 15, seed = 1) {
  set.seed(seed)
  n <- n_ref + n_num
  mu <- exp(runif(n_mirna, log(20), log(2000)))
  mu[181:200] <- runif(20, 0.2, 2)
  fold <- matrix(1, n_mirna, n)
  fold[seq_len(n_up), (n_ref + 1):n] <- 4
  depth <- exp(rnorm(n, 0, 0.2))
  m <- matrix(rnbinom(n_mirna * n, mu = mu * fold * rep(depth, each = n_mirna), size = 8),
              n_mirna, n)
  rownames(m) <- paste0("mir", seq_len(n_mirna))
  colnames(m) <- paste0("L", seq_len(n))
  samples <- data.frame(
    library_id = colnames(m),
    group = c(rep("ref_a", n_ref / 2), rep("ref_b", n_ref / 2), rep("num", n_num)),
    age = round(rnorm(n, 65, 8), 1),
    sex = sample(c("female", "male"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  list(counts = m, samples = samples)
}

test_that("counts are rounded to integers, half to even", {
  m <- round_counts(matrix(c(0.5, 1.5, 2.4, 7), 2, 2))
  expect_identical(storage.mode(m), "integer")
  expect_equal(as.vector(m), c(0L, 2L, 2L, 7L))
})

test_that("min_count pre-filter keeps miRNAs with 10 reads in 3 libraries", {
  m <- rbind(a = c(10, 10, 10, 0, 0), b = c(10, 10, 9, 9, 9), c = c(50, 50, 0, 0, 0),
             d = c(11, 12, 13, 14, 15))
  cond <- factor(c("reference", "reference", "numerator", "numerator", "numerator"))
  expect_equal(prefilter_keep(m, cond, "min_count"), c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(prefilter_keep(m, cond, "min_count", min_count = 50, min_samples = 2),
               c(FALSE, FALSE, TRUE, FALSE))
})

test_that("filterByExpr pre-filter equals edgeR with the comparison groups", {
  sim <- simulate()
  cond <- factor(ifelse(sim$samples$group == "num", "numerator", "reference"),
                 levels = c("reference", "numerator"))
  expected <- edgeR::filterByExpr(edgeR::DGEList(counts = sim$counts), group = cond)
  expect_equal(prefilter_keep(sim$counts, cond, "filterByExpr"), as.logical(expected))
  expect_true(sum(expected) < nrow(sim$counts))
  expect_error(prefilter_keep(sim$counts, cond, "other"), "Unknown pre-filter rule")
})

test_that("design holds covariates first and the group term last", {
  sim <- simulate()
  cond <- factor(ifelse(sim$samples$group == "num", "numerator", "reference"),
                 levels = c("reference", "numerator"))
  d <- build_design(sim$samples, cond, c("age", "sex"))
  expect_equal(paste(deparse(d$design), collapse = ""), "~age + sex + cond")
  expect_equal(mean(d$coldata$age), 0, tolerance = 1e-12)
  expect_equal(sd(d$coldata$age), 1, tolerance = 1e-12)
  expect_true(is.factor(d$coldata$sex))

  constant <- sim$samples
  constant$sex <- "male"
  d2 <- build_design(constant, cond, c("age", "sex"))
  expect_equal(paste(deparse(d2$design), collapse = ""), "~age + cond")

  missing_age <- sim$samples
  missing_age$age[3] <- NA
  expect_error(build_design(missing_age, cond, c("age", "sex")), "missing values")
  expect_error(build_design(sim$samples, cond, "bmi"), "not a column")
})

test_that("a comparison returns the DESeq2 table with shrunken fold changes", {
  sim <- simulate()
  fit <- fit_comparison(sim$counts, sim$samples, numerator = "num",
                        reference = c("ref_a", "ref_b"), prefilter = "filterByExpr")
  t <- fit$table
  expect_named(t, c("miRNA", "baseMean", "log2FoldChange", "lfcSE", "log2FoldChange_mle",
                    "lfcSE_mle", "stat", "pvalue", "padj"))
  expect_equal(c(fit$n_numerator, fit$n_reference), c(12, 12))
  expect_equal(fit$design, "~age + sex + cond")
  expect_true(fit$n_after_prefilter < fit$n_input)
  expect_true(fit$n_with_padj <= fit$n_after_prefilter)

  planted <- t[t$miRNA %in% paste0("mir", 1:15), ]
  expect_true(all(planted$log2FoldChange > 1))         # numerator relative to reference
  expect_true(all(planted$padj < 0.05))
  null <- t[!t$miRNA %in% paste0("mir", 1:15) & !is.na(t$padj), ]
  expect_true(mean(null$padj < 0.05 & abs(null$log2FoldChange) >= 0.5) < 0.05)

  # Wald p values and Benjamini-Hochberg adjustment over the miRNAs kept by independent filtering
  expect_equal(t$pvalue, 2 * pnorm(-abs(t$stat)), tolerance = 1e-8)
  ok <- !is.na(t$padj)
  expect_equal(t$padj[ok], p.adjust(t$pvalue[ok], "BH"), tolerance = 1e-10)
  # apeglm shrinks estimates towards zero
  expect_true(mean(abs(t$log2FoldChange) <= abs(t$log2FoldChange_mle) + 1e-6) > 0.95)
})

test_that("swapping numerator and reference flips the sign", {
  sim <- simulate(seed = 2)
  a <- fit_comparison(sim$counts, sim$samples, "num", c("ref_a", "ref_b"), "min_count")$table
  b <- fit_comparison(sim$counts, sim$samples, c("ref_a", "ref_b"), "num", "min_count")$table
  b <- b[match(a$miRNA, b$miRNA), ]
  expect_equal(a$log2FoldChange_mle, -b$log2FoldChange_mle, tolerance = 1e-4)
  expect_equal(a$pvalue, b$pvalue, tolerance = 1e-4)
})

test_that("covariates enter the model", {
  sim <- simulate(seed = 3)
  with_cov <- fit_comparison(sim$counts, sim$samples, "num", c("ref_a", "ref_b"), "min_count")
  without <- fit_comparison(sim$counts, sim$samples, "num", c("ref_a", "ref_b"), "min_count",
                            covariates = character(0))
  expect_equal(without$design, "~cond")
  expect_false(isTRUE(all.equal(with_cov$table$pvalue, without$table$pvalue)))
})

test_that("an empty group is reported", {
  sim <- simulate()
  expect_error(fit_comparison(sim$counts, sim$samples, "absent", "num", "min_count"),
               "has no libraries")
})

test_that("size factors are DESeq2 median-of-ratios factors", {
  sim <- simulate()
  expect_equal(unname(size_factors(sim$counts)),
               unname(DESeq2::estimateSizeFactorsForMatrix(sim$counts)), tolerance = 1e-12)
})
