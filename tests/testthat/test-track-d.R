# Track D: C++ permutation kernel, sensitivity, density weights, DE / CoDA on the model object.

test_that("the C++ permutation kernel reproduces the R statistics (block and Freedman-Lane)", {
  set.seed(4)
  sim <- simulateIndividualModel(n.per.group = c(A = 8, B = 8), p = 40, group.effect = 0.3, batch.levels = c("x", "y"), batch.effect = 0.3)
  des <- suppressMessages(buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = ~ group + batch))
  G <- cacoa:::gowerCenter(sim$D2); X <- des$F; cvec <- des$contrast.F
  pre <- cacoa:::inferencePrecompute(G, X, matrix(1, nrow(X), 1), cvec)
  P <- replicate(25, sample.int(nrow(X)))
  f.r <- apply(P, 2, function(p) cacoa:::permutedStats(G, X, matrix(1, nrow(X), 1), cvec, list(num = 1, den = 1), pre, p, TRUE, FALSE)["F"])
  f.c <- as.numeric(cacoa:::permuted_contrast_F(G, pre$a, pre$H, pre$cXc, pre$q, P))
  expect_equal(f.c, unname(f.r), tolerance = 1e-10)
  parts <- cacoa:::flGowerParts(G, X %*% cacoa:::contrastNullBasis(cvec))
  f.r2 <- apply(P, 2, function(p) { Gs <- cacoa:::flGowerPermute(parts, p); pre.b <- pre; pre.b$trG <- sum(diag(Gs)); pre.b$dG <- diag(Gs)
    cacoa:::permutedStats(Gs, X, matrix(1, nrow(X), 1), cvec, list(num = 1, den = 1), pre.b, NULL, TRUE, FALSE)["F"] })
  f.c2 <- as.numeric(cacoa:::permuted_contrast_F_fl(parts$K1, parts$K2, parts$K3, parts$K4, pre$a, pre$H, pre$cXc, pre$q, P))
  expect_equal(f.c2, unname(f.r2), tolerance = 1e-10)
  # identity permutation gives the observed F
  expect_equal(as.numeric(cacoa:::permuted_contrast_F(G, pre$a, pre$H, pre$cXc, pre$q, matrix(seq_len(nrow(X)), ncol = 1))), unname(cacoa:::contrastF(G, cacoa:::hatInfo(X), cvec)))
})

test_that("regression weights: balanced two-group weights, zero sums within strata, density uses them", {
  meta <- data.frame(group = factor(rep(c("A", "B"), each = 6)), batch = factor(rep(c("x", "y"), 6)), row.names = sprintf("s%02d", 1:12))
  d0 <- suppressMessages(buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group))
  w0 <- regressionWeights(d0)
  expect_equal(unname(w0[meta$group == "B"]), rep(1 / 6, 6)); expect_equal(unname(w0[meta$group == "A"]), rep(-1 / 6, 6))
  # the old weights F c differ only by a constant factor and offset here; signs agree
  old <- drop(d0$F %*% d0$contrast.F); expect_equal(sign(w0), sign(old - mean(old)))
  d1 <- buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch)
  w1 <- regressionWeights(d1)
  expect_equal(sum(w1), 0, tolerance = 1e-12)
  expect_equal(as.numeric(tapply(w1, meta$batch, sum)), c(0, 0), tolerance = 1e-12)
  expect_equal(sum(w1 * (meta$group == "B")), 1, tolerance = 1e-12)   # estimates the B - A difference
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 30, n.genes = 60)
  cao$embedding <- matrix(rnorm(2 * length(cao$cell.groups)), ncol = 2, dimnames = list(names(cao$cell.groups), c("x", "y")))
  res <- cao$estimateCellDensity(method = "kde", bins = 20, verbose = FALSE)
  expect_equal(sum(res$sample.weights), 0, tolerance = 1e-12)
  expect_s3_class(res$model, "cacoaModel")
  expect_equal(sum(res$cond.densities$num), 1, tolerance = 1e-8)
  dd <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 19, adjust = FALSE, verbose = FALSE)
  expect_true(length(dd$diff$permutation$raw) > 0)
})

test_that("sensitivity verdicts: robust under a null covariate, dependent under a confounded one", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  cao$sample.meta$nuisance <- factor(rep(c("u", "v"), 6))                      # unrelated to group
  cao$sample.meta$proxy <- factor(c(rep("p", 5), "q", "p", rep("q", 5)))        # almost identical to group
  cao$setModel(~ group + nuisance, test = "group", verbose = FALSE)
  cao$estimateExpressionShiftMagnitudes(n.permutations = 49, verbose = FALSE)
  sens <- cao$checkSensitivity(covariate.sets = list(unadjusted = ~ group, current = ~ group + nuisance, `plus proxy` = ~ group + nuisance + proxy),
                               n.permutations = 49, verbose = FALSE)
  expect_s3_class(sens, "cacoaSensitivity")
  expect_equal(levels(sens$table$model), c("unadjusted", "current", "plus proxy"))
  expect_true(all(c("celltype", "verdict", "max.rel.change", "influential.sample") %in% names(sens$summary)))
  # adjusting for the (near-collinear) proxy changes the answer: no cell type is "robust"
  sh <- sens$table[sens$table$effect == "shift", ]
  expect_equal(nrow(sh), 6)
  expect_true(all(grepl("depends on adjustment|changes by", sens$summary$verdict)))
  expect_output(print(sens), "Sensitivity of")
  expect_s3_class(cao$plotSensitivity(), "ggplot")
  # automatic covariate sets: unadjusted, current, minus nuisance (duplicates removed), plus screened candidates
  cao$screenCovariates(n.permutations = 29, verbose = FALSE)
  sens2 <- cao$checkSensitivity(n.permutations = 29, verbose = FALSE)
  expect_true(all(c("unadjusted", "current") %in% names(sens2$formulas)))
  # a planted outlier sample is flagged as influential
  cao2 <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.genes = 80, shift = 0.4, seed = 9)
  res0 <- cao2$estimateExpressionShiftMagnitudes(n.permutations = 29, verbose = FALSE)
  D <- res0$distances
  for (ct in names(D)) { D[[ct]]["B1", ] <- D[[ct]]["B1", ] * 4; D[[ct]][, "B1"] <- D[[ct]][, "B1"] * 4; diag(D[[ct]]) <- 0 }
  sens3 <- checkSensitivity(D, res0$model, cao2$sample.meta, formulas = list(current = ~ group), n.permutations = 29, seed = 1)
  expect_true(any(sens3$summary$influential.sample == "B1", na.rm = TRUE))
})

test_that("DE and CoDA run on the model object, including whole-factor tests", {
  skip_if_not_installed("limma"); skip_if_not_installed("edgeR"); skip_if_not_installed("DESeq2")
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 60, n.genes = 80, n.cell.types = 4, shift = 1.5)
  de <- cao$estimateDEPerCellType(test = "limma-voom", verbose = FALSE, min.cell.count = 5)
  expect_s3_class(attr(de, "model"), "cacoaModel")
  expect_true(all(c("ct1", "ct2") %in% names(de)))
  expect_true(any(de$ct1$res$padj < 0.1))
  cao$sample.meta$stage <- factor(rep(c("I", "II", "III"), length.out = 10))
  de3 <- suppressWarnings(cao$estimateDEPerCellType(test = "limma-voom", test.spec = "stage", formula = ~ stage, verbose = FALSE, name = "de.stage", min.cell.count = 5))
  expect_true(all(c("stat", "pvalue", "padj") %in% names(de3$ct1$res)))
  expect_equal(attr(de3, "model")$tests[[1]]$kind, "term")
  de4 <- suppressWarnings(cao$estimateDEPerCellType(test = "edgeR", test.spec = "stage", formula = ~ stage, verbose = FALSE, name = "de.stage2", min.cell.count = 5))
  expect_true(all(de4$ct1$res$pvalue >= 0 & de4$ct1$res$pvalue <= 1))
  expect_error(cao$estimateDEPerCellType(test = "Wilcoxon", test.spec = "stage", formula = ~ stage, verbose = FALSE), "two-group")
  # CoDA: contrast path unchanged, term path through the distance engine
  skip_if_not_installed("coda.base"); skip_if_not_installed("psych")
  coda <- cao$estimateCellLoadings(n.permutations = 49, verbose = FALSE)
  expect_s3_class(coda$model, "cacoaModel")
  coda3 <- cao$estimateCellLoadings(test = "stage", formula = ~ stage, n.permutations = 49, verbose = FALSE, name = "coda.stage")
  expect_equal(coda3$kind, "term")
  expect_true(all(c("R2.adj", "p.location", "p.dispersion") %in% names(coda3$results)))
})

test_that("cluster-free shifts on the engine: constant within whole-type neighbourhoods, planted type on top", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  cm <- Matrix::t(cao$data.object)                                  # genes x cells
  cg <- cao$cell.groups
  # neighbourhood of every cell = all cells of its type (0-based indices)
  idx <- split(seq_along(cg) - 1L, cg)
  nns <- lapply(as.character(cg), function(t) idx[[t]]); names(nns) <- names(cg)
  res <- clusterFreeExpressionShifts(cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, n.permutations = 49, seed = 1, min.n.obs.per.samp = 3)
  expect_equal(length(res$stat), length(cg))
  for (t in levels(cg)) expect_equal(length(unique(round(res$stat[cg == t], 10))), 1)   # identical statistic within a type
  expect_true(all(is.finite(res$z.adj)))
  expect_true(all(res$p.value[is.finite(res$p.value)] < 0.1))                            # planted shift in both types
  expect_equal(unname(res$n.samples[1]), 10)
  # the whole-type neighbourhood statistic equals the cell-type engine's F on the same mean profiles
  for (ct in levels(cg)) {
    prof <- sapply(levels(cao$sample.per.cell), function(s) Matrix::rowMeans(cm[, cg == ct & cao$sample.per.cell == s, drop = FALSE]))
    D <- 1 - cor(log10(1e3 * prof + 1))
    r <- testPairwiseEffects(list(x = D), cao$model, cao$sample.meta, dist = "cor", n.permutations = 19, seed = 1)
    expect_equal(unname(res$stat[cg == ct][1]), r$results$F, tolerance = 1e-8, info = ct)
  }
  # smoothing and adjustment can be switched off
  res2 <- clusterFreeExpressionShifts(cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, n.permutations = 19, seed = 1, adjust = FALSE, smooth = FALSE)
  expect_null(res2$z.adj); expect_null(res2$shifts.smoothed)
  # a whole-factor model is refused
  cao$sample.meta$stage <- factor(rep(c("I", "II", "III"), length.out = 10))
  m3 <- buildCacoaModel(cao$sample.meta, formula = ~ stage, test = "stage")
  expect_error(clusterFreeExpressionShifts(cm, cao$sample.per.cell, nns, m3, cao$sample.meta, n.permutations = 9), "contrast test")
})

test_that("regressions found by the notebook run: test lookup by variable, CoDA plot, volcano tables, FL fitter without residuals", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 60, n.genes = 80, n.cell.types = 4, shift = 1.5)
  cao$sample.meta$stage <- factor(rep(c("I", "II", "III"), length.out = 10))
  cao$setModel(~ group + stage, test = c("group: B vs A", "stage"), verbose = FALSE)
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 19, verbose = FALSE)
  # a test can be addressed by label, by variable name or by index
  expect_equal(cacoa:::matchModelTest(cao$model, "stage"), 2L); expect_equal(cacoa:::matchModelTest(cao$model, "stage (3 levels)"), 2L)
  expect_equal(cacoa:::matchModelTest(cao$model, "group"), 1L); expect_equal(cacoa:::matchModelTest(cao$model, 2), 2L)
  expect_error(cacoa:::matchModelTest(cao$model, "nothing"), "test not found")
  expect_s3_class(cao$plotExpressionShiftMagnitudes(test = "stage"), "ggplot")
  expect_s3_class(cao$plotExpressionShiftMagnitudes(test = 2), "ggplot")
  expect_equal(as.character(unname(cao$getSampleGroups("stage"))), as.character(cao$sample.meta$stage))
  # Freedman-Lane fitter: columns with missing values, residuals not requested
  X <- cbind(1, as.numeric(cao$sample.meta$group == "B")); Y <- matrix(rnorm(30), 10, 3); Y[c(2, 14)] <- NA
  f <- cacoa:::fl_fwl_cpp(X = X, Z = matrix(0, 0, 0), Y = Y, contrast = c(0, 1), n_randomizations = 9, return_residuals = FALSE, seed = 1L)
  expect_length(f$stat, 3); expect_true(all(is.finite(f$stat)))
  cao$embedding <- matrix(rnorm(2 * length(cao$cell.groups)), ncol = 2, dimnames = list(names(cao$cell.groups), c("x", "y")))
  cao$estimateCellDensity(method = "kde", bins = 20, verbose = FALSE)
  dd <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 19, adjust = FALSE, verbose = FALSE)   # default perm.method = freedman-lane
  expect_true(all(is.finite(dd$diff$permutation$raw)))
  # CoDA loadings plot on the lmCoda structure
  skip_if_not_installed("coda.base"); skip_if_not_installed("psych")
  coda <- cao$estimateCellLoadings(n.permutations = 49, verbose = FALSE)
  expect_equal(dim(coda$contrast$loadings$perm), c(4, 49))
  expect_s3_class(cao$plotCellLoadings(), "ggplot")
  expect_s3_class(cao$plotCellLoadings(show.pvals = FALSE, show.null = FALSE, ordering = "loading"), "ggplot")
  cao$estimateCellLoadings(test = "stage", formula = ~ stage, n.permutations = 19, verbose = FALSE, name = "coda.stage")
  expect_error(cao$plotCellLoadings(name = "coda.stage"), "no per-cell-type loadings")
  # volcano plot reads the DE tables (list(res = ...)) and finds CellFrac
  skip_if_not_installed("limma"); skip_if_not_installed("EnhancedVolcano")
  cao$setModel(~ group, test = "group", verbose = FALSE)
  de <- cao$estimateDEPerCellType(test = "limma-voom", verbose = FALSE, min.cell.count = 5)
  expect_true("CellFrac" %in% names(de$ct1$res))
  expect_s3_class(cao$plotVolcano(cell.types = "ct1"), "ggplot")
  expect_true(inherits(cao$plotVolcano(cell.types = c("ct1", "ct2")), "ggplot"))
})
