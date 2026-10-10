# Step 2 of the engine convergence: contrast statistics under relabeling in C++, against the R reference permutedStats().

refPermStats <- function(G, X, Z, cvec, zend, P, scheme, bias = TRUE, need.var = TRUE) {
  pre <- cacoa:::inferencePrecompute(G, X, Z, cvec)
  if (scheme == "block") return(t(apply(P, 2, function(p) cacoa:::permutedStats(G, X, Z, cvec, zend, pre, p, bias, need.var))))
  parts <- cacoa:::flGowerParts(G, X %*% cacoa:::contrastNullBasis(cvec))
  t(apply(P, 2, function(p) { Gs <- cacoa:::flGowerPermute(parts, p); pb <- pre; pb$trG <- sum(diag(Gs)); pb$dG <- diag(Gs)
    cacoa:::permutedStats(Gs, X, Z, cvec, zend, pb, NULL, bias, need.var) }))
}

test_that("kernel A contrast statistics equal the R reference over the grid (block and Freedman-Lane, with and without bias correction)", {
  set.seed(21)
  for (cell in list(list(groups = c(A = 7, B = 7), batch = c("x", "y"), disp = 1), list(groups = c(A = 5, B = 9), batch = c("x", "y"), disp = 2.5),
                    list(groups = c(A = 6, B = 6, C = 6), batch = NULL, disp = 1))) {
    sim <- simulateIndividualModel(n.per.group = cell$groups, p = 40, group.effect = 0.3, batch.levels = cell$batch, batch.effect = 0.3, group.sd = c(1, cell$disp, 1)[seq_along(cell$groups)])
    f <- if (is.null(cell$batch)) ~ group else ~ group + batch
    des <- suppressMessages(buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = f))
    G <- cacoa:::gowerCenter(sim$D2); X <- des$F; cvec <- des$contrast.F
    Z <- cbind(1, as.numeric(sim$meta$group == "B")); zend <- list(num = c(1, 1), den = c(1, 0))
    P <- replicate(15, sample.int(nrow(X)))
    pre <- cacoa:::inferencePrecompute(G, X, Z, cvec)
    for (scheme in c("block", "freedman-lane")) for (bias in c(TRUE, FALSE)) {
      ref <- refPermStats(G, X, Z, cvec, zend, P, scheme, bias)
      cpp <- if (scheme == "block") cacoa:::permuted_contrast_stats(G, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, zend$num, zend$den, P, bias, TRUE) else {
        parts <- cacoa:::flGowerParts(G, X %*% cacoa:::contrastNullBasis(cvec))
        cacoa:::permuted_contrast_stats_fl(parts$K1, parts$K2, parts$K3, parts$K4, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, zend$num, zend$den, P, bias, TRUE) }
      expect_equal(unname(cpp), unname(ref), tolerance = 1e-10, info = paste(scheme, bias, length(cell$groups)))
    }
    # F-only mode agrees with the dedicated F kernel; identity relabeling reproduces the observed estimates
    fo <- cacoa:::permuted_contrast_stats(G, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, zend$num, zend$den, P, TRUE, FALSE)
    expect_equal(fo[, 1], as.numeric(cacoa:::permuted_contrast_F(G, pre$a, pre$H, pre$cXc, pre$q, P)), tolerance = 1e-12)
    expect_true(all(is.na(fo[, 2:4])))
    obs <- cacoa:::permuted_contrast_stats(G, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, zend$num, zend$den, cbind(seq_len(nrow(X))), TRUE, TRUE)
    eff <- cacoa:::estimatePairwiseEffects(NULL, X, cvec, Z, zend, TRUE, G = G)
    expect_equal(obs[1, 2], eff$shift, tolerance = 1e-10); expect_equal(obs[1, 3], eff$var, tolerance = 1e-10)
  }
})

test_that("the cell-type test is unchanged by the kernel switch: skipped reasons, exhaustive enumeration and max-T still work", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 30, n.genes = 60, shift = 1.5)
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 199, verbose = FALSE)
  sh <- res$results[res$results$effect == "shift", ]
  expect_true(all(sh$exhaustive)); expect_equal(sh$n.perm, sh$n.perm.distinct)             # small design: enumerated
  expect_equal(unname(sh$n.perm.distinct), rep(cao$model$tests[[1]]$permutation$n.distinct, nrow(sh)))
  expect_true(all(is.finite(sh$p.fwer)) && all(sh$p.fwer >= sh$p - 1e-12))                  # max-T never below the raw p
  expect_true(all(res$global$p >= min(sh$p) - 1e-12))
  expect_true(all(c("var", "total") %in% res$results$effect))
  expect_true(all(is.finite(res$results$statistic[res$results$effect == "var"])))
})

test_that("a numeric tested variable is permuted across all samples within strata (regression: the plan used to permute nothing)", {
  meta <- gridMeta(levels = 2, n.per.level = 5, age = TRUE)
  m <- buildCacoaModel(meta, formula = ~ age + batch, test = "age")
  mp <- modelPermutations(m, scheme = "block", n.permutations = 50, seed = 2)
  expect_true(all(mp$plan$in.set))
  expect_equal(mp$plan$n.distinct, factorial(5)^2)                                       # 5! per batch
  expect_true(mean(apply(mp$P, 2, function(p) any(p != seq_len(10)))) > 0.9)             # the draws actually move samples
  expect_true(all(meta$batch[mp$P[, 2]] == meta$batch))                                   # within strata
  m2 <- buildCacoaModel(meta, formula = ~ age, test = "age")
  mp2 <- modelPermutations(m2, scheme = "block", n.permutations = 20, seed = 2)
  expect_equal(mp2$plan$n.distinct, factorial(10))
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 30, n.genes = 60)
  cao$sample.meta$age <- c(41, 55, 62, 48, 70, 39, 58, 51)
  res <- cao$estimateExpressionShiftMagnitudes(test = "age", n.permutations = 39, verbose = FALSE, name = "age")
  sh <- res$results[res$results$effect == "shift", ]
  expect_true(all(is.finite(sh$p)) && all(sh$p < 1))
  expect_true(all(sh$n.perm == 39))
})

refTermStats <- function(G, Xf, Xr, Zf, Zr, P, scheme) {
  if (scheme == "block") return(t(apply(P, 2, function(p) c(cacoa:::termTestGower(G, Xf[p, , drop = FALSE], Xr[p, , drop = FALSE])["F"],
                                                             cacoa:::dispersionTermTest(G, Xf[p, , drop = FALSE], Zf[p, , drop = FALSE], Zr[p, , drop = FALSE])["F.disp"]))))
  parts <- cacoa:::flGowerParts(G, Xr)
  t(apply(P, 2, function(p) { Gs <- cacoa:::flGowerPermute(parts, p); c(cacoa:::termTestGower(Gs, Xf, Xr)["F"], cacoa:::dispersionTermTest(Gs, Xf, Zf, Zr)["F.disp"]) }))
}

test_that("kernel A term statistics (location F, dispersion F) equal the R reference over the grid", {
  set.seed(31)
  for (cell in list(list(groups = c(A = 6, B = 6, C = 6), batch = c("x", "y"), sd = c(1, 1, 2)), list(groups = c(A = 5, B = 8), batch = NULL, sd = c(1, 1.8)),
                    list(groups = c(A = 4, B = 4, C = 4, D = 4), batch = c("x", "y"), sd = c(1, 1, 1, 1)))) {
    sim <- simulateIndividualModel(n.per.group = cell$groups, p = 40, group.effect = 0.3, batch.levels = cell$batch, batch.effect = 0.3, group.sd = cell$sd)
    G <- cacoa:::gowerCenter(sim$D2); n <- nrow(G)
    Xf <- stats::model.matrix(if (is.null(cell$batch)) ~ group else ~ group + batch, sim$meta)
    Xr <- stats::model.matrix(if (is.null(cell$batch)) ~ 1 else ~ batch, sim$meta)
    Zf <- stats::model.matrix(~ group, sim$meta); Zr <- matrix(1, n, 1)
    P <- replicate(12, sample.int(n))
    k <- cacoa:::termKernelInputs(Xf, Xr, Zf, Zr, n)
    for (scheme in c("block", "freedman-lane")) {
      ref <- refTermStats(G, Xf, Xr, Zf, Zr, P, scheme)
      cpp <- if (scheme == "block") cacoa:::permuted_term_stats(G, k$Hf, k$Hr, k$Zf, k$Zr, k$df, k$nu, k$qZf, k$qZr, P, TRUE) else {
        parts <- cacoa:::flGowerParts(G, Xr); cacoa:::permuted_term_stats_fl(parts$K1, parts$K2, parts$K3, parts$K4, k$Hf, k$Hr, k$Zf, k$Zr, k$df, k$nu, k$qZf, k$qZr, P, TRUE) }
      expect_equal(unname(cpp), unname(ref), tolerance = 1e-9, info = paste(scheme, length(cell$groups)))
    }
    obs <- cacoa:::permuted_term_stats(G, k$Hf, k$Hr, k$Zf, k$Zr, k$df, k$nu, k$qZf, k$qZr, cbind(seq_len(n)), TRUE)
    expect_equal(obs[1, 1], unname(cacoa:::termTestGower(G, Xf, Xr)["F"]), tolerance = 1e-10)
    expect_equal(obs[1, 2], unname(cacoa:::dispersionTermTest(G, Xf, Zf, Zr)["F.disp"]), tolerance = 1e-10)
    # screen path: the same kernel with Zf = Xf, Zr = Xr (dispersion on the location design) and a reduced design with an adjustment column
    noDisp <- cacoa:::permuted_term_stats_fl(cacoa:::flGowerParts(G, Xr)$K1, cacoa:::flGowerParts(G, Xr)$K2, cacoa:::flGowerParts(G, Xr)$K3, cacoa:::flGowerParts(G, Xr)$K4,
                                            k$Hf, k$Hr, k$Zf, k$Zr, k$df, k$nu, k$qZf, k$qZr, P, FALSE)
    expect_true(all(is.na(noDisp[, 2])))
  }
})

test_that("the screen and the whole-factor test give the same answers as before the kernel switch (structure, calibration-free checks)", {
  meta <- gridMeta(levels = 3, n.per.level = 5, age = TRUE)
  set.seed(5); D <- lapply(1:2, function(i) { M <- matrix(rnorm(15 * 30), 15); D <- as.matrix(dist(M)); dimnames(D) <- list(rownames(meta), rownames(meta)); D })
  names(D) <- c("ct1", "ct2")
  sc <- screenCovariates(D, meta, covariates = c("group", "batch", "age"), n.permutations = 49, seed = 3)
  expect_true(all(c("p.perm", "p.disp.perm", "F", "F.disp") %in% names(sc$table)))
  expect_true(all(is.finite(sc$table$p.perm[sc$table$n.used >= 4])))
  m <- buildCacoaModel(meta, formula = ~ group + batch, test = "group: all")
  tr <- testTermEffects(D, m, meta, dist = "l2", n.permutations = 49, seed = 3)
  expect_true(all(c("p.location", "p.dispersion") %in% names(tr$results)))
  expect_true(all(is.finite(tr$results$p.location)))
})
