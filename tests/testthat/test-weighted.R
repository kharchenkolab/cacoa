# Step 4b: weighted individual-level model (R reference): weights of one reproduce the plain fit; robust weights
# down-weight a planted outlying sample; weak imputation reproduces the dropped-sample fit.

test_that("weights of one reproduce the plain contrast and term statistics", {
  set.seed(7)
  sim <- simulateIndividualModel(n.per.group = c(A = 7, B = 7), p = 40, group.effect = 0.3, batch.levels = c("x", "y"), batch.effect = 0.3, group.sd = c(1, 1.8))
  des <- suppressMessages(buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = ~ group + batch))
  G <- cacoa:::gowerCenter(sim$D2); X <- des$F; cvec <- des$contrast.F; Z <- cbind(1, as.numeric(sim$meta$group == "B")); zend <- list(num = c(1, 1), den = c(1, 0))
  w1 <- cacoa:::weightedContrastStats(G, X, Z, cvec, zend)
  e <- cacoa:::estimatePairwiseEffects(NULL, X, cvec, Z, zend, TRUE, G = G)
  expect_equal(w1$shift, e$shift, tolerance = 1e-10); expect_equal(w1$var, e$var, tolerance = 1e-10); expect_equal(w1$total, e$total, tolerance = 1e-10)
  expect_equal(w1$F, unname(cacoa:::contrastF(G, cacoa:::hatInfo(X), cvec)), tolerance = 1e-10)
  Xf <- stats::model.matrix(~ group + batch, sim$meta); Xr <- stats::model.matrix(~ batch, sim$meta); Zf <- stats::model.matrix(~ group, sim$meta); Zr <- matrix(1, 14, 1)
  t1 <- cacoa:::weightedTermStats(G, Xf, Xr, Zf, Zr)
  expect_equal(t1$F, unname(cacoa:::termTestGower(G, Xf, Xr)["F"]), tolerance = 1e-10)
  expect_equal(t1$F.disp, unname(cacoa:::dispersionTermTest(G, Xf, Zf, Zr)["F.disp"]), tolerance = 1e-10)
  # the weighted hat matrix is the W-orthogonal projection: H' W = W H and H X = X
  w <- runif(14, 0.2, 1); hw <- cacoa:::hatInfoW(X, w)
  expect_equal(unname(t(hw$H) %*% diag(w)), unname(diag(w) %*% hw$H), tolerance = 1e-10); expect_equal(as.vector(hw$H %*% X), as.vector(X), tolerance = 1e-10)
})

test_that("robust weights down-weight a planted outlying sample and pull the estimate back", {
  set.seed(8)
  sim <- simulateIndividualModel(n.per.group = c(A = 8, B = 8), p = 40, group.effect = 0.3)
  des <- suppressMessages(buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = ~ group))
  X <- des$F; cvec <- des$contrast.F; Z <- cbind(1, as.numeric(sim$meta$group == "B")); zend <- list(num = c(1, 1), den = c(1, 0))
  G0 <- cacoa:::gowerCenter(sim$D2); clean <- cacoa:::weightedContrastStats(G0, X, Z, cvec, zend)
  D2 <- sim$D2; D2[3, ] <- D2[3, ] * 6; D2[, 3] <- D2[, 3] * 6; diag(D2) <- 0           # sample 3 far from everything
  G <- cacoa:::gowerCenter(D2)
  plain <- cacoa:::weightedContrastStats(G, X, Z, cvec, zend)
  hub <- cacoa:::weightedContrastStats(G, X, Z, cvec, zend, robust = "huber")
  win <- cacoa:::weightedContrastStats(G, X, Z, cvec, zend, robust = "winsor")
  expect_lt(hub$w[3], 0.5); expect_true(all(hub$w[-3] > 0.9)); expect_lt(win$w[3], 0.5)
  expect_lt(abs(hub$shift - clean$shift), abs(plain$shift - clean$shift))
  expect_lt(abs(hub$var - clean$var), abs(plain$var - clean$var))
  # the dropped-sample fit is the limit of very small weight
  w <- rep(1, 16); w[3] <- 1e-6
  drop3 <- cacoa:::weightedContrastStats(cacoa:::gowerCenter(D2[-3, -3]), X[-3, ], Z[-3, ], cvec, zend)
  tiny <- cacoa:::weightedContrastStats(G, X, Z, cvec, zend, w = w)
  expect_equal(tiny$shift, drop3$shift, tolerance = 1e-3); expect_equal(tiny$var, drop3$var, tolerance = 1e-3)
})

test_that("weak imputation of an absent sample reproduces the dropped-sample fit while keeping the full sample set", {
  set.seed(9)
  sim <- simulateIndividualModel(n.per.group = c(A = 6, B = 6), p = 40, group.effect = 0.4, batch.levels = c("x", "y"), batch.effect = 0.2)
  des <- suppressMessages(buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = ~ group + batch))
  all.s <- rownames(sim$meta); present <- all.s[-c(2, 9)]
  imp <- cacoa:::imputeAbsentSamples(sim$D2, present, all.s)
  expect_equal(dim(imp$D2), c(12, 12)); expect_equal(unname(imp$w[c(2, 9)]), c(1e-4, 1e-4))
  X <- des$F; cvec <- des$contrast.F; Z <- cbind(1, as.numeric(sim$meta$group == "B")); zend <- list(num = c(1, 1), den = c(1, 0))
  weak <- cacoa:::weightedContrastStats(cacoa:::gowerCenter(imp$D2), X, Z, cvec, zend, w = imp$w)
  keep <- all.s %in% present
  drop <- cacoa:::weightedContrastStats(cacoa:::gowerCenter(sim$D2[keep, keep]), X[keep, ], Z[keep, ], cvec, zend)
  expect_equal(weak$shift, drop$shift, tolerance = 1e-3); expect_equal(weak$var, drop$var, tolerance = 1e-3); expect_equal(weak$F, drop$F, tolerance = 1e-3)
  tw <- cacoa:::weightedTermStats(cacoa:::gowerCenter(imp$D2), stats::model.matrix(~ group + batch, sim$meta), stats::model.matrix(~ batch, sim$meta),
                                  stats::model.matrix(~ group, sim$meta), matrix(1, 12, 1), w = imp$w)
  td <- cacoa:::weightedTermStats(cacoa:::gowerCenter(sim$D2[keep, keep]), stats::model.matrix(~ group + batch, sim$meta[keep, ]), stats::model.matrix(~ batch, sim$meta[keep, ]),
                                  stats::model.matrix(~ group, sim$meta[keep, ]), matrix(1, 10, 1))
  expect_equal(tw$F, td$F, tolerance = 1e-3); expect_equal(tw$F.disp, td$F.disp, tolerance = 1e-3)
})

test_that("C++ weighted kernels equal the R reference under relabeling: plain, huber and winsor weights, weak imputation, block and FL", {
  set.seed(11)
  sim <- simulateIndividualModel(n.per.group = c(A = 7, B = 7), p = 40, group.effect = 0.3, batch.levels = c("x", "y"), batch.effect = 0.3, group.sd = c(1, 1.6))
  des <- suppressMessages(buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = ~ group + batch))
  D2 <- sim$D2; D2[5, ] <- D2[5, ] * 4; D2[, 5] <- D2[, 5] * 4; diag(D2) <- 0
  G <- cacoa:::gowerCenter(D2); X <- des$F; cvec <- des$contrast.F; Z <- cbind(1, as.numeric(sim$meta$group == "B")); zend <- list(num = c(1, 1), den = c(1, 0))
  n <- 14; P <- replicate(8, sample.int(n)); base <- rep(1, n); weak <- base; weak[c(2, 11)] <- 1e-4
  Xf <- stats::model.matrix(~ group + batch, sim$meta); Xr <- stats::model.matrix(~ batch, sim$meta); Zf <- stats::model.matrix(~ group, sim$meta); Zr <- matrix(1, n, 1)
  for (rob in c("none", "huber", "winsor")) for (w0 in list(base, weak)) {
    rc <- c(none = 0L, huber = 1L, winsor = 2L)[[rob]]
    obs <- cacoa:::weightedContrastStats(G, X, Z, cvec, zend, w = w0, robust = rob)
    cf <- cacoa:::weighted_contrast_fit(G, X, Z, w0, cvec, zend$num, zend$den, TRUE, rc)
    expect_equal(as.numeric(cf$stats), c(obs$F, obs$shift, obs$var, obs$total), tolerance = 1e-8, info = rob)
    expect_equal(as.numeric(cf$w), unname(obs$w), tolerance = 1e-8)
    # block: relabel the design, weights stay with the samples
    ref <- t(apply(P, 2, function(p) { o <- cacoa:::weightedContrastStats(G, X[p, ], Z[p, , drop = FALSE], cvec, zend, w = w0, robust = rob); c(o$F, o$shift, o$var, o$total) }))
    cpp <- cacoa:::permuted_contrast_stats_w(G, X, Z, w0, cvec, zend$num, zend$den, P, FALSE, w0, TRUE, TRUE, rc)
    expect_equal(unname(cpp), unname(ref), tolerance = 1e-8, info = paste("block", rob))
    # Freedman-Lane: residual kernel of the weighted reduced fit (final observed weights) re-indexed
    parts <- cacoa:::flGowerPartsW(G, X %*% cacoa:::contrastNullBasis(cvec), obs$w)
    reff <- t(apply(P, 2, function(p) { Gs <- parts$K1 + parts$K2[, p] + parts$K3[p, ] + parts$K4[p, p]; o <- cacoa:::weightedContrastStats(Gs, X, Z, cvec, zend, w = w0, robust = rob); c(o$F, o$shift, o$var, o$total) }))
    cppf <- cacoa:::permuted_contrast_stats_w(G, X, Z, w0, cvec, zend$num, zend$den, P, TRUE, obs$w, TRUE, TRUE, rc)
    expect_equal(unname(cppf), unname(reff), tolerance = 1e-8, info = paste("fl", rob))
    # terms
    tobs <- cacoa:::weightedTermStats(G, Xf, Xr, Zf, Zr, w = w0, robust = rob)
    tref <- t(apply(P, 2, function(p) { o <- cacoa:::weightedTermStats(G, Xf[p, ], Xr[p, , drop = FALSE], Zf[p, ], Zr[p, , drop = FALSE], w = w0, robust = rob); c(o$F, o$F.disp) }))
    tcpp <- cacoa:::permuted_term_stats_w(G, Xf, Xr, Zf, Zr, w0, P, FALSE, w0, TRUE, rc)
    expect_equal(unname(tcpp), unname(tref), tolerance = 1e-8, info = paste("term", rob))
    tparts <- cacoa:::flGowerPartsW(G, Xr, tobs$w)
    treff <- t(apply(P, 2, function(p) { Gs <- tparts$K1 + tparts$K2[, p] + tparts$K3[p, ] + tparts$K4[p, p]; o <- cacoa:::weightedTermStats(Gs, Xf, Xr, Zf, Zr, w = w0, robust = rob); c(o$F, o$F.disp) }))
    tcppf <- cacoa:::permuted_term_stats_w(G, Xf, Xr, Zf, Zr, w0, P, TRUE, tobs$w, TRUE, rc)
    expect_equal(unname(tcppf), unname(treff), tolerance = 1e-8, info = paste("term fl", rob))
  }
  # with unit weights and no robust step the weighted kernels equal the fast kernels
  pre <- cacoa:::inferencePrecompute(G, X, Z, cvec)
  fast <- cacoa:::permuted_contrast_stats(G, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, zend$num, zend$den, P, TRUE, TRUE)
  slow <- cacoa:::permuted_contrast_stats_w(G, X, Z, base, cvec, zend$num, zend$den, P, FALSE, base, TRUE, TRUE, 0L)
  expect_equal(unname(slow), unname(fast), tolerance = 1e-8)
})

test_that("robust and na.mode options run end to end: contrast and whole-factor tests, screen, sensitivity, cluster-free", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.genes = 80, shift = 1.2)
  res0 <- cao$estimateExpressionShiftMagnitudes(n.permutations = 39, verbose = FALSE)
  D <- res0$distances
  # a planted outlying sample in ct1 (its distances scaled up) and a sample missing from ct2
  D$ct1["A2", ] <- D$ct1["A2", ] * 4; D$ct1[, "A2"] <- D$ct1[, "A2"] * 4; diag(D$ct1) <- 0
  D$ct2 <- D$ct2[rownames(D$ct2) != "B3", colnames(D$ct2) != "B3"]
  plain <- testPairwiseEffects(D, cao$model, cao$sample.meta, n.permutations = 39, seed = 1)
  hub <- testPairwiseEffects(D, cao$model, cao$sample.meta, n.permutations = 39, seed = 1, robust = "huber")
  expect_lt(hub$fits$ct1$w[["A2"]], 0.5); expect_equal(names(which.min(hub$fits$ct1$w)), "A2"); expect_gte(sum(hub$fits$ct1$w > 0.99), 9)
  expect_true(all(is.finite(hub$results$p.shift)) && all(is.finite(hub$results$p.var)))
  expect_equal(hub$call.info$robust, "huber")
  expect_true(abs(hub$results$shift[1] - res0$fits[[1]]$ct1$shift) < abs(plain$results$shift[1] - res0$fits[[1]]$ct1$shift))   # pulled back towards the clean estimate
  # weak imputation: ct2 keeps all 12 samples, estimates match the dropped-sample fit
  weak <- testPairwiseEffects(D, cao$model, cao$sample.meta, n.permutations = 39, seed = 1, na.mode = "impute_weak")
  expect_equal(weak$results$n[weak$results$celltype == "ct2"], 12); expect_equal(plain$results$n[plain$results$celltype == "ct2"], 11)
  expect_equal(weak$results$shift[2], plain$results$shift[2], tolerance = 1e-3); expect_equal(unname(weak$results$var[2]), unname(plain$results$var[2]), tolerance = 1e-2)
  expect_equal(unname(weak$fits$ct2$base.w[["B3"]]), 1e-4)
  expect_equal(weak$fits$ct2$influence$se[["shift"]], plain$fits$ct2$influence$se[["shift"]], tolerance = 0.05)
  # the R6 entry point, options and provenance
  cao$setOptions(robust = "winsor", na.mode = "impute_weak")
  r6 <- cao$estimateExpressionShiftMagnitudes(n.permutations = 29, verbose = FALSE, name = "rob")
  expect_equal(r6$settings$robust, "winsor"); expect_equal(r6$settings$na.mode, "impute_weak")
  gg <- cao$plotExpressionShiftMagnitudes(name = "rob"); expect_true(grepl("robust: winsor", gg$labels$subtitle))
  expect_error(cao$setOptions(robust = "median"), "robust must be")
  cao$setOptions(robust = "none", na.mode = "drop")
  # whole-factor test with robust weights and weak imputation
  cao$sample.meta$stage <- factor(rep(c("I", "II", "III"), length.out = 12))
  m3 <- buildCacoaModel(cao$sample.meta, formula = ~ stage, test = "stage")
  t0 <- testTermEffects(D, m3, cao$sample.meta, n.permutations = 29, seed = 2)
  t1 <- testTermEffects(D, m3, cao$sample.meta, n.permutations = 29, seed = 2, robust = "huber", na.mode = "impute_weak")
  expect_true(all(is.finite(t1$results$p.location))); expect_equal(t1$results$n[2], 12); expect_equal(t0$results$n[2], 11)
  expect_lt(t1$fits$ct1$w[["A2"]], 0.5)
  # screen and sensitivity accept the options
  sc <- screenCovariates(D, cao$sample.meta, covariates = c("group", "stage"), n.permutations = 29, seed = 3, robust = "huber")
  expect_true(all(is.finite(sc$table$p.perm[sc$table$n.used >= 4]))); expect_equal(sc$settings$robust, "huber")
  sens <- checkSensitivity(D, cao$model, cao$sample.meta, formulas = list(current = ~ group), n.permutations = 19, seed = 1, robust = "huber", na.mode = "impute_weak")
  expect_s3_class(sens, "cacoaSensitivity")
  # cluster-free robust kernel equals a per-cell weighted R reference for whole-type neighbourhoods
  cg <- cao$cell.groups; idx <- split(seq_along(cg) - 1L, cg)
  nns <- lapply(as.character(cg), function(t) idx[[t]]); names(nns) <- names(cg)
  cm <- Matrix::t(cao$data.object)
  cf <- clusterFreeExpressionShifts(cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, n.permutations = 19, seed = 1, adjust = FALSE, smooth = FALSE, robust = "huber")
  expect_true(all(is.finite(cf$stat)))
  for (t in levels(cg)) {
    prof <- sapply(levels(cao$sample.per.cell), function(s) Matrix::rowMeans(cm[, cg == t & cao$sample.per.cell == s, drop = FALSE]))
    G <- cacoa:::gowerCenter(1 - cor(log10(1e3 * prof + 1)))
    X <- cao$model$F; cvec <- cao$model$contrast.F
    ref <- cacoa:::weightedContrastStats(G, X, matrix(1, nrow(X), 1), cvec, list(num = 1, den = 1), robust = "huber")
    expect_equal(unname(cf$stat[cg == t][1]), ref$F, tolerance = 1e-8, info = t); expect_equal(unname(cf$shifts[cg == t][1]), ref$shift, tolerance = 1e-8, info = t)
  }
})
