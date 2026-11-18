# intervals of the loadings, pseudocount, per-sample diagnostics, grouping terms and the corrected jackknife error

simComp <- function(n = 8, D = 6, seed = 1, paired = TRUE) {
  set.seed(seed); pat <- rep(paste0("p", seq_len(n)), 2); cond <- rep(c("a", "b"), each = n)
  u <- matrix(rnorm(n * D, 0, 0.8), n, D)[rep(seq_len(n), 2), ]; eff <- c(1.2, -1.2, rep(0, D - 2))
  lp <- matrix(rep(log(c(0.3, 0.2, 0.2, 0.15, 0.1, 0.05)[seq_len(D)]), each = 2 * n), 2 * n) + outer(cond == "b", eff) + (if (paired) u else 0) + rnorm(2 * n * D, 0, 0.3)
  cn <- t(apply(lp, 1, function(l) rmultinom(1, 2000, exp(l)))); dimnames(cn) <- list(paste0("s", seq_len(2 * n)), paste0("t", seq_len(D)))
  list(cnts = cn, meta = data.frame(row.names = rownames(cn), cond = factor(cond), patient = pat, stringsAsFactors = FALSE))
}

test_that("loadings carry model-based standard errors that equal those of lm on the log ratios", {
  s <- simComp(); m <- buildCacoaModel(s$meta, ~ cond, test = "cond")
  r <- cacoa:::lmCoda(s$cnts, m, perm.method = "freedman-lane", n.permutations = 99, seed = 1)
  tab <- codaLoadingsTable(r, min.effect = 1.5)
  expect_true(all(c("loading", "se", "ci.low", "ci.high", "beyond.min.effect", "median.cells", "n.zero", "pseudocount.shift") %in% names(tab)))
  Y <- r$ilr %*% t(r$psi); j <- tab$celltype[1]; co <- summary(lm(Y[, j] ~ s$meta$cond))$coefficients[2, ]
  expect_equal(tab$loading[1] * log(2), unname(co[1]), tolerance = 1e-8); expect_equal(tab$se[1] * log(2), unname(co[2]), tolerance = 1e-8)
  expect_true(all(tab$beyond.min.effect[tab$celltype %in% c("t1", "t2")]))
  expect_s3_class(plotCodaLoadings(r, min.effect = 1.5), "ggplot")
})

test_that("the pseudocount is half a cell in a median sample, scaled with the sample's cell number", {
  cn <- rbind(a = c(10, 0, 90), b = c(100, 0, 900), c = c(30, 5, 265)); colnames(cn) <- c("x", "y", "z")
  f <- cacoa:::computeILRMatrix(cn)$freqs; n <- rowSums(cn); pc <- 0.5 * n / median(n)
  expect_equal(unname(f[, "y"]), unname((cn[, "y"] + pc) / (n + 3 * pc)))
  expect_equal(unname(f["a", "y"]), unname(f["b", "y"]))            # no cells: the same proportion in a small and in a large sample
  f0 <- cacoa:::computeILRMatrix(cn, zero.pseudocount = 0.1, scale.pseudocount = FALSE)$freqs
  expect_equal(unname(f0["a", "y"]), 0.1 / (100 + 0.3))
})

test_that("per-sample influence equals the change of the loading on refitting without the sample", {
  s <- simComp(seed = 2, paired = FALSE); m <- buildCacoaModel(s$meta, ~ cond, test = "cond")
  r <- cacoa:::lmCoda(s$cnts, m, perm.method = "freedman-lane", n.permutations = 49, seed = 1, scale.pseudocount = FALSE)
  d <- codaSampleDiagnostics(r); expect_equal(nrow(d), prod(dim(s$cnts)))
  Y <- r$ilr %*% t(r$psi); colnames(Y) <- rownames(r$psi); x <- as.numeric(s$meta$cond == "b"); i <- 3
  full <- coef(lm(Y[, "t1"] ~ x))[2]; loo <- coef(lm(Y[-i, "t1"] ~ x[-i]))[2]; se <- r$contrast$loadings$se[["t1"]]
  expect_equal(d$influence[d$sample == rownames(s$cnts)[i] & d$celltype == "t1"], unname((loo - full) / se), tolerance = 1e-8)
  expect_s3_class(plotCodaSamples(r, n.top = 4), "ggplot"); expect_s3_class(plotCodaSamples(r, n.top = 4, pair.by = "patient"), "ggplot")
  expect_s3_class(plotSampleHeatmap(matrix(d$residual, nrow(s$cnts), dimnames = list(rownames(s$cnts), unique(d$celltype))), mark = 3), "ggplot")
})

test_that("grouping terms: the transformation is the inverse square root of the covariance, and the variance ratio is recovered", {
  g <- rep(1:4, c(1, 2, 3, 4)); th <- 1.7; W <- cacoa:::groupWhitener(g, th); U <- cacoa:::groupUnwhitener(g, th); I <- diag(length(g))
  V <- I + th * outer(g, g, "=="); expect_equal(unname(W(I) %*% V %*% t(W(I))), I, tolerance = 1e-10); expect_equal(unname(U(W(I))), I, tolerance = 1e-10)
  set.seed(3); g <- rep(1:30, each = 4); Fm <- cbind(1, rep(0:1, 60)); Y <- sapply(1:40, function(j) rnorm(30, 0, sqrt(2))[g] + rnorm(120))
  expect_lt(abs(log(cacoa:::remlVarianceRatio(Y, Fm, g) / 2)), log(1.6))
  sp <- cacoa:::splitGroupingTerms(~ cond + sex + (1 | patient)); expect_equal(sp$groups, "patient"); expect_equal(deparse(sp$formula), "~cond + sex")
})

test_that("a random intercept for the patient uses the pairing in the composition test", {
  s <- simComp(n = 8, seed = 4)
  m0 <- buildCacoaModel(s$meta, ~ cond, test = "cond"); m1 <- buildCacoaModel(s$meta, ~ cond + (1 | patient), test = "cond")
  expect_equal(m1$random.info$structure, "within")
  r0 <- cacoa:::lmCoda(s$cnts, m0, perm.method = "freedman-lane", n.permutations = 99, seed = 1); r1 <- cacoa:::lmCoda(s$cnts, m1, perm.method = "freedman-lane", random.baseline = TRUE, n.permutations = 99, seed = 1)
  # default: the grouping term only restricts the relabeling to within patients, as in the other analyses
  rd <- cacoa:::lmCoda(s$cnts, m1, perm.method = "block", n.permutations = 99, seed = 1)
  expect_null(rd$random); expect_true(all(apply(rd$fit$P, 2, function(p) all(s$meta$patient[p] == s$meta$patient))))
  expect_equal(unname(rd$contrast$loadings$se), unname(r0$contrast$loadings$se)); expect_error(cacoa:::lmCoda(s$cnts, m0, random.baseline = TRUE, n.permutations = 9), "grouping term")
  expect_equal(r1$random$group, "patient"); expect_true(median(r1$random$icc) > 0.3)
  expect_true(all(r1$contrast$loadings$se[c("t1", "t2")] < r0$contrast$loadings$se[c("t1", "t2")]))      # the patient effect is removed from the error
  expect_true(all(r1$contrast$per_cell$pval[c("t1", "t2")] <= 0.02))
  # between-group structure: the tested variable is constant within the group; the composition model still tests it
  s$meta$donor <- paste0(s$meta$cond, rep(1:4, length.out = nrow(s$meta)))
  m2 <- buildCacoaModel(s$meta, ~ cond + (1 | donor), test = "cond"); expect_equal(m2$random.info$structure, "between")
  r2 <- cacoa:::lmCoda(s$cnts, m2, perm.method = "freedman-lane", random.baseline = TRUE, n.permutations = 49, seed = 1); expect_true(all(is.finite(r2$contrast$per_cell$pval)))
})

test_that("the corrected jackknife error lies between the jackknife error and that error over root two", {
  expect_equal(cacoa:::correctedJackknifeSE(2, 0), 2 / sqrt(2)); expect_equal(cacoa:::correctedJackknifeSE(2, 100), 2); expect_equal(cacoa:::correctedJackknifeSE(2, -5), 2 / sqrt(2))
  expect_equal(cacoa:::correctedJackknifeSE(2, NA), 2)
  set.seed(5); n <- 6; g <- rep(1:2, each = n); X <- matrix(rnorm(200 * 2 * n), 200); D <- 1 - cor(X); dimnames(D) <- list(paste0("s", 1:(2 * n)), paste0("s", 1:(2 * n)))
  meta <- data.frame(row.names = rownames(D), cond = factor(g)); m <- buildCacoaModel(meta, ~ cond, test = "cond")
  f <- expressionShiftsForModel(list(ct = D), m, dist = "cor", n.permutations = 19, influence = TRUE, seed = 1)$fits[[1]][["ct"]]
  jk <- f$influence$se.jackknife[["shift"]]; expect_lte(f$influence$se[["shift"]], jk + 1e-12); expect_gte(f$influence$se[["shift"]], jk / sqrt(2) - 1e-12)
})

test_that("with a patient covariate the cell types are tested by the studentized loading, and the test holds its level", {
  # no condition effect, strong patient effect, 7 pairs + 3 unpaired samples: the unstudentized loading rejected in 21 % of such tests
  # before the residuals of the smaller model were brought to the scale of the errors (kappa); it now holds its level as well
  set.seed(7); pat <- c(rep(paste0("p", 1:7), each = 2), "p8", "p9", "p10"); cond <- c(rep(c("a", "b"), 7), "a", "b", "b"); P <- sapply(1:199, function(b) sample.int(17))
  Fm <- cbind(cond == "a", cond == "b", model.matrix(~ factor(pat))[, -1]) + 0; Fm <- Fm[, qr(Fm)$pivot[seq_len(qr(Fm)$rank)]]; cv <- c(-1, 1, rep(0, ncol(Fm) - 2))
  Y <- sapply(1:300, function(j) rnorm(10, 0, 2)[as.integer(factor(pat))] + rnorm(17))
  ct <- cacoa:::contrastPermutationTest(Y, Fm, cv, P)
  p.t <- (rowSums(abs(ct$t.perm) >= abs(ct$t)) + 1) / 200; p.e <- (rowSums(abs(ct$est.perm) >= abs(ct$est)) + 1) / 200
  expect_lt(mean(p.t < 0.05), 0.09); expect_lt(mean(p.e < 0.05), 0.09); expect_gt(ct$kappa, 1.5)
  s <- simComp(seed = 8); m <- buildCacoaModel(s$meta, ~ cond + patient, test = "cond")
  r <- cacoa:::lmCoda(s$cnts, m, perm.method = "freedman-lane", n.permutations = 99, seed = 1)
  expect_equal(r$settings$per.celltype.statistic, "t"); expect_match(r$contrast$global$statistic, "squared t")
})

test_that("when all distinct relabelings are enumerated the observed labeling is not counted twice", {
  s <- simComp(n = 5); m <- buildCacoaModel(s$meta, ~ cond, test = "cond", block.vars = "patient")   # 5 pairs: 2^5 = 32 relabelings within patients
  r <- cacoa:::lmCoda(s$cnts, m, perm.method = "block", n.permutations = 199, seed = 1)
  expect_true(isTRUE(r$fit$exhaustive)); expect_equal(ncol(r$fit$P), 32)
  p <- r$contrast$per_cell$pval
  expect_equal(p * 32, round(p * 32)); expect_true(all(p >= 2 / 32))                 # the mirror image of the observed labeling ties with it
  expect_equal(min(p), 2 / 32)                                                       # the two changed cell types reach the smallest attainable value
  expect_equal(unname(r$fit$pval * 32), round(unname(r$fit$pval * 32)))
  r2 <- cacoa:::lmCoda(s$cnts, buildCacoaModel(s$meta, ~ cond, test = "cond"), perm.method = "block", n.permutations = 199, seed = 1)   # random draws: (k + 1) / (B + 1)
  expect_false(isTRUE(r2$fit$exhaustive)); expect_equal(r2$contrast$per_cell$pval * 200, round(r2$contrast$per_cell$pval * 200))
})

test_that("exchange of residuals: the full model is refitted to every relabeled data set and the residual scale is corrected", {
  s <- simComp(n = 8, seed = 2); m <- buildCacoaModel(s$meta, ~ cond + patient, test = "cond"); n <- nrow(s$meta)
  set.seed(3); Y <- matrix(rnorm(n * 5), n, 5, dimnames = list(rownames(s$meta), paste0("y", 1:5))); P <- vapply(1:60, function(b) sample.int(n), integer(n))
  f <- cacoa:::performLMPermutations(m, Y, perm.method = "freedman-lane", P = P, return.sampled.stats = TRUE)
  ct <- cacoa:::contrastPermutationTest(Y, m$F, m$contrast.F, apply(P, 2, order))        # the R function moves the residuals, the fitter the design: inverse relabelings
  expect_equal(unname(t(f$stats.perm)), unname(ct$t.perm), tolerance = 1e-8)             # t statistic of the full model on the relabeled data
  expect_equal(unname(t(f$effects.perm)), unname(ct$est.perm), tolerance = 1e-8)         # estimates on the corrected residual scale
  expect_true(ct$kappa > 1.5 && ct$kappa < 2.4)                                          # 16 samples, 9 parameters without the tested variable: about 15 / 8
  expect_true(abs(mean(apply(ct$est.perm, 1, var) / ct$se^2) - 1) < 0.5)                # relabeled estimates as variable as the observed one
  Pw <- vapply(1:60, function(b) { p <- seq_len(n); for (g in split(seq_len(n), s$meta$patient)) if (runif(1) < 0.5) p[g] <- rev(g); p }, integer(n))
  expect_equal(cacoa:::contrastPermutationTest(Y, m$F, m$contrast.F, Pw)$kappa, 1)       # exchange within the fitted patients: nothing to correct
  fb <- cacoa:::performLMPermutations(m, Y, perm.method = "block", P = Pw, return.sampled.stats = TRUE)
  cb <- cacoa:::contrastPermutationTest(Y, m$F, m$contrast.F, Pw, scheme = "block")
  expect_equal(unname(t(fb$stats.perm)), unname(cb$t.perm), tolerance = 1e-8)            # exchange of the tested variable: the same statistic in both routes
})

test_that("distance model under the exchange of residuals: the Gower parts are those of the reduced fit plus the rescaled residuals", {
  s <- simComp(n = 8, seed = 5); m <- buildCacoaModel(s$meta, ~ cond + patient, test = "cond"); n <- nrow(s$meta)
  set.seed(6); Y <- scale(matrix(rnorm(n * 7), n, 7), scale = FALSE); G <- tcrossprod(Y); P <- vapply(1:40, function(b) sample.int(n), integer(n))
  Xr <- m$F %*% cacoa:::contrastNullBasis(m$contrast.F); a <- drop(t(cacoa:::hatInfo(m$F)$A) %*% m$contrast.F)
  kap <- cacoa:::flKappa(a, Xr, P); expect_equal(kap, cacoa:::contrastPermutationTest(Y, m$F, m$contrast.F, P)$kappa, tolerance = 1e-10); expect_gt(kap, 1.5)
  parts <- cacoa:::flGowerParts(G, Xr, kap); p <- P[, 3]; F0 <- qr.fitted(qr(Xr), Y); R0 <- Y - F0
  expect_equal(cacoa:::flGowerPermute(parts, p), tcrossprod(F0 + sqrt(kap) * R0[p, ]), tolerance = 1e-8, ignore_attr = TRUE)
  Pw <- vapply(1:20, function(b) { q <- seq_len(n); for (g in split(seq_len(n), s$meta$patient)) if (runif(1) < 0.5) q[g] <- rev(g); q }, integer(n))
  expect_equal(cacoa:::flKappa(a, Xr, Pw), 1)                                           # exchange within the fitted patients
})

test_that("samples dropped for a missing covariate, a continuous tested variable, and an exclusion that removes a level of the tested variable", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.cell.types = 4, seed = 7)
  meta <- cao$sample.meta; meta$sex <- rep(c("f", "m"), length.out = nrow(meta)); meta$sex[c(2, 9)] <- NA; cao$sample.meta <- meta
  # composition with samples dropped by the model: the counts follow the model's samples
  cao$setModel(~ group + sex, test = "group", verbose = FALSE); expect_equal(length(cao$model$samples$dropped), 2)
  r <- cao$estimateCellLoadings(verbose = FALSE); expect_equal(nrow(r$cnts), nrow(meta) - 2); expect_setequal(rownames(r$cnts), rownames(cao$model$meta))
  # cell-type expression shifts and DE with the dropped samples
  sh <- cao$estimateExpressionShiftMagnitudes(verbose = FALSE, influence = FALSE, n.permutations = 19); expect_true(all(sh$results$n <= nrow(meta) - 2, na.rm = TRUE))
  # a continuous tested variable: the per-sample plot draws the samples against it
  cao$setModel(~ age + group, test = "age", verbose = FALSE); r2 <- cao$estimateCellLoadings(name = "coda.age", verbose = FALSE)
  ga <- cao$plotCodaSamples(name = "coda.age", n.top = 3); expect_s3_class(ga, "ggplot")
  txt <- unlist(lapply(ga$data, function(col) if (is.character(col) || is.factor(col)) unique(as.character(col))))
  expect_false(any(grepl("(^|[^.0-9])1-fold", txt)))      # a change per unit near 1 keeps two decimals (1.04-fold), not "1-fold"
  # excluding every sample of one level: the stored model no longer applies and is removed with a message, not an error
  cao$setModel(~ group, test = "group", verbose = FALSE)
  expect_message(cao$excludeSamples(rownames(meta)[meta$group == "B"]), "does not apply")
  expect_null(cao$model); expect_equal(nrow(cao$sample.meta), sum(meta$group == "A"))
  cao$excludeSamples(NULL, verbose = FALSE); expect_equal(nrow(cao$sample.meta), nrow(meta))
  # a column added while samples are excluded is kept by the next exclusion (values of the excluded samples unknown: NA)
  cao$excludeSamples(rownames(meta)[1], verbose = FALSE)
  m2 <- cao$sample.meta; m2$extra <- seq_len(nrow(m2)); m2$group <- factor(as.character(m2$group)); cao$sample.meta <- m2
  cao$excludeSamples(rownames(meta)[2], verbose = FALSE)
  expect_true("extra" %in% names(cao$sample.meta)); expect_equal(cao$sample.meta[rownames(meta)[1], "extra"], NA_integer_)
  expect_equal(cao$sample.meta[rownames(meta)[3], "extra"], m2[rownames(meta)[3], "extra"])
  cao$excludeSamples(NULL, verbose = FALSE); expect_equal(nrow(cao$sample.meta), nrow(meta)); expect_true(is.na(cao$sample.meta[rownames(meta)[1], "extra"]))
})

test_that("the step of a continuous tested variable is kept when DE fits the grouping term as a fixed term and after excludeSamples()", {
  skip_if_not_installed("DESeq2")
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.genes = 60, n.cell.types = 2, seed = 8)
  lfc <- function(step) {
    cao$setModel(~ age + (1 | batch), test = "age", numeric.step = step, verbose = FALSE)
    de <- suppressWarnings(suppressMessages(cao$estimateDEPerCellType(verbose = FALSE, name = paste0("de", step))))
    r <- cacoa:::getResDf(de[[1]]); setNames(r$log2FoldChange, r$Gene)
  }
  l1 <- lfc(1); l10 <- lfc(10); g <- intersect(names(l1), names(l10)); g <- g[abs(l1[g]) > 1e-6]
  expect_equal(unname(median(l10[g] / l1[g])), 10, tolerance = 1e-3)
  cao$excludeSamples(rownames(cao$sample.meta)[1], verbose = FALSE)
  expect_equal(cao$model$numeric.step, 10)
})

test_that("composition sensitivity fits alternative models on samples the current model drops", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.cell.types = 4, shift = 1, seed = 7)
  meta <- cao$sample.meta; meta$sex <- rep(c("f", "m"), length.out = nrow(meta)); meta$sex[c(2, 9)] <- NA; cao$sample.meta <- meta
  cao$setModel(~ group + sex, test = "group", verbose = FALSE)
  suppressMessages(cao$estimateCellLoadings(n.permutations = 29, verbose = FALSE))
  s <- suppressMessages(cao$checkSensitivity(analysis = "composition", covariate.sets = list(unadjusted = ~ group), n.permutations = 29, jackknife = FALSE, verbose = FALSE))
  expect_equal(s$dropped, character(0))
  expect_equal(unique(s$table$n[s$table$model == "unadjusted"]), nrow(meta))      # all 12 samples, not the 10 of the current model
})
