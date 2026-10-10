# End-to-end expression-shift estimation on the synthetic Cacoa object (new engine, Track B.5 / B.6).

test_that("estimateExpressionShiftMagnitudes runs end to end and finds a planted shift", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 99, verbose = FALSE)
  expect_s3_class(res$results, "data.frame")
  expect_setequal(unique(res$results$celltype), c("ct1", "ct2"))
  expect_setequal(unique(res$results$effect), c("shift", "var", "total"))
  expect_true(all(c("estimate", "estimate.norm", "se.jk", "p", "padj", "p.fwer", "scheme", "p.floor", "flags") %in% names(res$results)))
  sh <- res$results[res$results$effect == "shift", ]
  expect_true(all(sh$p > 0 & sh$p <= 1))
  expect_true(all(sh$p < 0.1))                    # strong planted shift
  expect_true(all(sh$estimate > 0))
  expect_equal(nrow(res$global), 3)
  expect_true(res$global$p[res$global$effect == "shift"] < 0.1)
  expect_identical(res, cao$test.results$expression.shifts)
  expect_s3_class(res$model, "cacoaModel")
  expect_equal(res$settings$n.permutations, 99)
  expect_equal(res$settings$seed, 1)              # option default
  # distances and adjusted distances per cell type, with cell counts attached
  expect_setequal(names(res$distances), c("ct1", "ct2"))
  expect_equal(dim(res$adjusted.distances[[1]]$ct1), c(10, 10))
  expect_equal(unname(attr(res$distances$ct1, "n.cells")[1]), 20)   # 40 cells per sample, 2 cell types
  # influence: samples x cell types
  expect_equal(dim(res$influence[[1]]$shift), c(10, 2))
  # the same seed reproduces the p-values
  res2 <- cao$estimateExpressionShiftMagnitudes(n.permutations = 99, verbose = FALSE, name = "again")
  expect_equal(res2$results$p, res$results$p)
  # pseudobulk and distances are cached
  expect_false(is.null(cao$cache$pseudobulk)); expect_false(is.null(cao$cache$sample.distances))
})

test_that("per-call model override, options precedence and legacy arguments", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  cao$setOptions(n.permutations = 49, verbose = FALSE)
  expect_message(res <- cao$estimateExpressionShiftMagnitudes(formula = ~ group + batch, verbose = TRUE, name = "adj"), "temporary model")
  expect_equal(res$settings$n.permutations, 49)
  expect_true("batch" %in% all.vars(res$model$formula))
  expect_true(all(res$results$scheme == "block"))
  expect_equal(all.vars(cao$model$formula), "group")          # stored model unchanged
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 19, verbose = FALSE)
  expect_equal(res$settings$n.permutations, 19)                # explicit call wins over options
  # legacy arguments are translated / announced, unknown ones are errors
  expect_message(res <- cao$estimateExpressionShiftMagnitudes(perm.method = "freedman-lane", dist.type = "shift", n.permutations = 19, verbose = TRUE), "perm.method")
  expect_true(all(res$results$scheme == "freedman-lane"))
  expect_error(cao$estimateExpressionShiftMagnitudes(nonsense = 1, verbose = FALSE), "unknown argument")
  # numeric test and term test through the same entry point
  res <- cao$estimateExpressionShiftMagnitudes(test = "age", n.permutations = 19, verbose = FALSE, name = "age")
  expect_setequal(unique(res$results$effect), c("shift", "var", "total"))
  cao$sample.meta$stage <- factor(rep(c("I", "II", "III"), length.out = nrow(cao$sample.meta)))
  res <- cao$estimateExpressionShiftMagnitudes(test = "stage: all", n.permutations = 19, verbose = FALSE, name = "stage")
  expect_setequal(unique(res$results$effect), c("location", "dispersion"))
  expect_true(all(c("R2.adj", "F", "p.location", "F.disp", "p.dispersion") %in% names(res$wide[[1]])))
  expect_equal(dim(res$fits[[1]]$ct1$cell.table), c(3, 3))
  expect_equal(nrow(res$fits[[1]]$ct1$pair.table), 3)
  # several tests at once
  res <- cao$estimateExpressionShiftMagnitudes(formula = ~ group + stage, test = c("group", "stage: all"), n.permutations = 19, verbose = FALSE, name = "two")
  expect_equal(length(unique(res$results$test)), 2)
  expect_equal(nrow(res$global), 5)
})

test_that("a null design gives no systematic signal and plots render", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 0, seed = 7)
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 99, verbose = FALSE)
  expect_true(all(res$results$p[res$results$effect == "shift"] > 0.01))
  gg <- cao$plotExpressionShiftMagnitudes()
  expect_s3_class(gg, "ggplot")
  expect_match(gg$labels$subtitle, "B vs A")
  expect_match(gg$labels$subtitle, "permutations")
  gg2 <- cao$plotExpressionShiftMagnitudes(type = "bar", normalized = FALSE, significance = "p.fwer", show.provenance = FALSE)
  expect_s3_class(gg2, "ggplot"); expect_null(gg2$labels$subtitle)
  gg3 <- cao$plotExpressionShiftMagnitudes(effects = "shift", show.pvalues = "adjusted")   # legacy argument ignored
  expect_s3_class(gg3, "ggplot")
  expect_s3_class(cao$plotSampleInfluence(), "ggplot")
  expect_warning(cao$plotExpressionShiftResiduals(), "deprecated")
})

test_that("cell types with too few samples per level are skipped with a reason", {
  cao <- makeToyCacoa(n.per.group = c(A = 3, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1)
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 19, min.samp.per.level = 4, verbose = FALSE)
  expect_equal(nrow(res$results), 0)
  expect_equal(nrow(res$skipped), 2)
  expect_match(res$skipped$reason[1], "fewer than 4 samples")
  expect_error(cao$plotExpressionShiftMagnitudes(), "no cell type")
})

test_that("sample distance accessor and plots work on the new results", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  cao$estimateExpressionShiftMagnitudes(n.permutations = 19, verbose = FALSE)
  D <- cao$getSampleDistanceMatrix(space = "expression.shifts", cell.type = "ct1")
  expect_equal(dim(D), c(10, 10)); expect_equal(unname(diag(D)), rep(0, 10))
  Dj <- cao$getSampleDistanceMatrix(space = "expression.shifts")              # joint over cell types
  expect_equal(dim(Dj), c(10, 10))
  Da <- cao$getSampleDistanceMatrix(space = "expression.shifts", values = "adjusted", cell.type = "ct1")
  expect_equal(dim(Da), c(10, 10)); expect_true(all(Da >= -1e-10))
  expect_s3_class(cao$plotSampleDistances(space = "expression.shifts", values = "unadjusted", color.by = "group"), "ggplot")
  expect_true(inherits(cao$plotSampleDistances(space = "expression.shifts", values = "both", color.by = "group"), "ggplot"))
  sep <- cao$estimateMetadataSeparation(cao$sample.meta[, c("group", "batch")], n.permutations = 99, show.warning = FALSE, verbose = FALSE)
  expect_setequal(names(sep$pvalues), c("group", "batch"))
})
