# End-to-end expression-shift estimation on the synthetic Cacoa object (current pair-LM engine; the
# engine is replaced in Track A.3/A.4 but the public entry point and result shape are kept).

test_that("estimateExpressionShiftMagnitudes runs end to end on a synthetic object", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  set.seed(1)
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 99, verbose = FALSE)
  expect_s3_class(res$results, "data.frame")
  expect_setequal(res$results$celltype, c("ct1", "ct2"))
  expect_true(all(res$results$pvalue > 0 & res$results$pvalue <= 1))
  # a strong planted shift is detected in both cell types
  expect_true(all(res$results$pvalue < 0.1))
  expect_identical(res, cao$test.results$expression.shifts)
  # per-call formula override with a covariate also runs and stores under its own name
  res2 <- cao$estimateExpressionShiftMagnitudes(formula = ~ group + batch, n.permutations = 49, verbose = FALSE, name = "adj")
  expect_setequal(res2$results$celltype, c("ct1", "ct2"))
  expect_true("adj" %in% names(cao$test.results))
})

test_that("a null design gives no systematic signal", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 0, seed = 7)
  set.seed(2)
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 99, verbose = FALSE)
  expect_true(all(res$results$pvalue > 0.01))
})
