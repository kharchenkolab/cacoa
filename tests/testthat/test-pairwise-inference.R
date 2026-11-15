# Inference layer: permutation plans, generators, induced permutations, permuted statistics, p-values,
# max-statistic combination, calibration on small simulations.

simMeta <- function(n.per.group = c(A = 6, B = 6, C = 4), seed = 1) {
  set.seed(seed)
  g <- factor(rep(names(n.per.group), n.per.group), levels = names(n.per.group))
  data.frame(group = g, batch = factor(rep(c("b1", "b2"), length.out = length(g))),
             age = round(rnorm(length(g), 50, 8)), row.names = sprintf("s%02d", seq_along(g)))
}

test_that("block plan permutes only the compared conditions within strata, counts relabelings, enumerates", {
  meta <- simMeta()
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch)
  plan <- permutationPlan(des, meta, scheme = "auto", n.permutations = 99)
  expect_equal(plan$scheme, "block")
  expect_equal(unname(plan$in.set), meta$group %in% c("A", "B"))
  expect_equal(nlevels(plan$strata), 2)
  # distinct relabelings: per batch, choose which of the A/B samples are labelled B
  nd <- prod(sapply(split(meta$group[meta$group != "C"], meta$batch[meta$group != "C"]), function(g) choose(length(g), sum(g == "B"))))
  expect_equal(plan$n.distinct, nd)
  expect_false(plan$exhaustive)
  P <- drawPermutations(plan, 50)
  expect_equal(dim(P), c(nrow(meta), 50))
  fixed <- which(meta$group == "C")
  expect_true(all(P[fixed, ] == fixed))
  for (b in 1:50) {
    p <- P[, b]
    expect_setequal(p, seq_len(nrow(meta)))
    expect_true(all(meta$batch[p] == meta$batch))          # strata respected
    expect_true(all(meta$group[p] %in% c("A", "B") == (meta$group %in% c("A", "B"))))
  }
  # exhaustive enumeration: every distinct relabeling exactly once, identity included
  plan2 <- permutationPlan(des, meta, scheme = "block", n.permutations = 10000)
  expect_true(plan2$exhaustive)
  P2 <- drawPermutations(plan2)
  expect_equal(ncol(P2), nd)
  keys <- apply(P2, 2, function(p) paste(meta$group[p], collapse = ""))
  expect_equal(length(unique(keys)), nd)
  expect_true(paste(meta$group, collapse = "") %in% keys)
  expect_equal(plan2$p.floor, 1 / nd)
})

test_that("plan handles interaction cells, marginal and numeric contrasts, and continuous nuisance", {
  meta <- simMeta(c(A = 8, B = 8))
  # interaction cell B:b1 vs A:b1 -> only b1 samples of A/B swap
  des <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch,
                                     contrast = list(type = "simple", term = "group:batch", num = "B:b1", den = "A:b1"))
  plan <- permutationPlan(des, meta)
  expect_equal(unname(plan$in.set), meta$batch == "b1")
  # simple contrast at batch = b2
  des2 <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch,
                                      contrast = list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b2")))
  expect_equal(unname(permutationPlan(des2, meta)$in.set), meta$batch == "b2")
  # marginal over batch -> all A/B samples, strata by batch
  des3 <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch,
                                      contrast = list(type = "marginal", term = "group", num = "B", den = "A", over = "batch"))
  p3 <- permutationPlan(des3, meta)
  expect_true(all(p3$in.set)); expect_equal(nlevels(p3$strata), 2)
  # numeric contrast with discrete nuisance: block scheme over all samples, n! per stratum
  des4 <- cacoa:::buildDesignMatrices(meta, contrast = c(age = 1), formula = ~ batch + age)
  p4 <- permutationPlan(des4, meta)
  expect_equal(p4$scheme, "block"); expect_true(all(p4$in.set))
  expect_equal(p4$log.n.distinct, sum(lgamma(table(meta$batch) + 1)))
  # continuous nuisance -> Freedman-Lane
  des5 <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + age)
  p5 <- permutationPlan(des5, meta)
  expect_equal(p5$scheme, "freedman-lane"); expect_match(p5$notes[1], "Freedman-Lane")
  expect_equal(permutationPlan(des5, meta, scheme = "block")$scheme, "block")
  # block.vars add strata
  p6 <- permutationPlan(des5, meta, scheme = "block", block.vars = "batch")
  expect_equal(nlevels(p6$strata), 2)
})

test_that("induced permutations are uniform within strata and coupled", {
  meta <- simMeta(c(A = 3, B = 3))
  meta$batch <- factor("b1")
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group)
  plan <- permutationPlan(des, meta, scheme = "block", n.permutations = 10000)
  P <- drawPermutations(plan)                      # all 20 relabelings of 6 samples (3 A, 3 B) ... as index perms
  # induce on a 4-sample subset (2 A, 2 B): each of its 6 relabelings should appear equally often over all
  # 720 full permutations
  full <- cacoa:::distinctLabelPermutations(as.character(1:6))   # all 720 permutations
  sub <- c(1, 2, 4, 5)
  sub.plan <- permutationPlan(des, meta, samples = rownames(meta)[sub], scheme = "block", n.permutations = 10000)
  keys <- sapply(full, function(p) paste(meta$group[sub][cacoa:::inducePermutation(p, sub, sub.plan)], collapse = ""))
  tab <- table(keys)
  expect_equal(length(tab), 6)
  expect_true(all(tab == 720 / 6))
})

test_that("permuted statistics equal direct recomputation on the relabelled design (block and FL)", {
  sim <- simulateIndividualModel(c(A = 7, B = 8), p = 40, group.effect = 0.4, group.sd = c(A = 1, B = 1.4),
                                 batch.levels = c("b1", "b2"), batch.effect = 0.3, seed = 4)
  meta <- sim$meta
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch)
  eff <- pairwiseEffectsFromDesign(sqrt(sim$D2), des, meta, dist = "l2")
  pre <- cacoa:::inferencePrecompute(eff$G, eff$X, eff$Z, eff$contrast)
  obs <- cacoa:::permutedStats(eff$G, eff$X, eff$Z, eff$contrast, eff$z.end, pre)
  expect_equal(unname(obs[c("shift", "var", "total")]), c(eff$shift, eff$var, eff$total))
  expect_equal(unname(obs["F"]), eff$F.shift)
  plan <- permutationPlan(des, meta, scheme = "block")
  P <- drawPermutations(plan, 5)
  for (b in 1:5) {
    p <- P[, b]
    st <- cacoa:::permutedStats(eff$G, eff$X, eff$Z, eff$contrast, eff$z.end, pre, p)
    direct <- estimatePairwiseEffects(NULL, eff$X[p, ], eff$contrast, eff$Z[p, ], eff$z.end, G = eff$G)
    expect_equal(unname(st[c("shift", "var", "total")]), c(direct$shift, direct$var, direct$total))
    expect_equal(unname(st["F"]), cacoa:::contrastF(eff$G, cacoa:::hatInfo(eff$X[p, ]), eff$contrast))
  }
  # Freedman-Lane Gower identity on Euclidean data: G* equals the Gram matrix of H_r Y + P R_r Y
  Y <- sim$Y; Yc <- sweep(Y, 2, colMeans(Y)); G <- tcrossprod(Yc)
  X <- model.matrix(~ group + batch, meta); cvec <- c(0, 1, 0)
  Xr <- X %*% cacoa:::contrastNullBasis(cvec); Hr <- cacoa:::hatInfo(Xr)$H; Rr <- diag(nrow(X)) - Hr
  parts <- cacoa:::flGowerParts(G, Xr)
  p <- sample.int(nrow(X))
  Ystar <- Hr %*% Yc + (Rr %*% Yc)[p, ]
  expect_equal(cacoa:::flGowerPermute(parts, p), tcrossprod(Ystar), tolerance = 1e-8)
})

test_that("testPairwiseEffects: results table, reproducibility, exhaustive p-values, max-T", {
  sim <- simulateIndividualModel(c(A = 6, B = 6, C = 4), p = 40, group.effect = 0.8 * sqrt(40), batch.levels = c("b1", "b2"),
                                 batch.effect = 0.3, seed = 18)
  meta <- sim$meta
  D1 <- sqrt(sim$D2)
  sim2 <- simulateIndividualModel(c(A = 6, B = 6, C = 4), p = 40, group.effect = 0, seed = 9)   # null cell type
  D2 <- sqrt(sim2$D2); dimnames(D2) <- dimnames(D1)
  D3 <- D1[-c(1, 8), -c(1, 8)]                                                                  # a sample subset
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch)
  res <- testPairwiseEffects(list(ct1 = D1, ct0 = D2, ct3 = D3), des, meta, dist = "l2", n.permutations = 99, seed = 3)
  expect_setequal(res$results$celltype, c("ct1", "ct0", "ct3"))
  expect_true(all(c("shift", "var", "total", "F", "p.shift", "p.var", "p.total", "padj.shift", "p.fwer.shift",
                    "scheme", "n.perm.distinct", "exhaustive", "p.floor", "flags") %in% names(res$results)))
  expect_true(all(res$results$scheme == "block"))
  expect_true(all(res$results$p.fwer.shift >= res$results$p.shift - 1e-12))
  expect_equal(res$global$effect, c("shift", "var", "total"))
  expect_true(all(is.finite(res$global$p)))
  expect_lte(res$global$p[1], min(res$results$p.fwer.shift) + 1e-12)
  # a planted strong shift (delta^2 = 25.6) is estimated and found in ct1, not in the null cell type
  expect_lt(abs(res$results$shift[res$results$celltype == "ct1"] - (0.8 * sqrt(40))^2), 15)
  expect_lt(res$results$p.shift[res$results$celltype == "ct1"], 0.05)
  expect_gt(res$results$p.shift[res$results$celltype == "ct0"], 0.05)
  # reproducibility: same seed -> identical; different seed -> different permutation stats
  res2 <- testPairwiseEffects(list(ct1 = D1, ct0 = D2, ct3 = D3), des, meta, dist = "l2", n.permutations = 99, seed = 3)
  expect_equal(res$results, res2$results)
  res4 <- testPairwiseEffects(list(ct1 = D1, ct0 = D2, ct3 = D3), des, meta, dist = "l2", n.permutations = 99, seed = 3, n.cores = 2)
  expect_equal(res$results, res4$results)
  res3 <- testPairwiseEffects(list(ct1 = D1), des, meta, dist = "l2", n.permutations = 99, seed = 4, return.perm.stats = TRUE)
  res3b <- testPairwiseEffects(list(ct1 = D1), des, meta, dist = "l2", n.permutations = 99, seed = 3, return.perm.stats = TRUE)
  expect_false(isTRUE(all.equal(res3$perm.stats$ct1, res3b$perm.stats$ct1)))
  # set.seed() applies when seed is NULL
  set.seed(11); a <- testPairwiseEffects(list(ct1 = D1), des, meta, dist = "l2", n.permutations = 49)
  set.seed(11); b <- testPairwiseEffects(list(ct1 = D1), des, meta, dist = "l2", n.permutations = 49)
  expect_equal(a$results, b$results)
  # exhaustive: a tiny design enumerates all relabelings and p is a multiple of 1/n.distinct
  small <- rownames(meta)[c(1:3, 7:9, 13)]
  rs <- testPairwiseEffects(list(ct = D1[small, small]), des, meta, dist = "l2", n.permutations = 999, seed = 1)
  expect_true(rs$results$exhaustive)
  ms <- meta[small, ]; ab <- ms$group != "C"
  expect_equal(rs$results$n.perm.distinct, prod(sapply(split(ms$group[ab], ms$batch[ab]), function(g) choose(length(g), sum(g == "B")))))
  expect_equal(rs$results$p.shift * rs$results$n.perm.distinct, round(rs$results$p.shift * rs$results$n.perm.distinct))
  expect_match(rs$results$flags, "few-permutations")
  # skipped cell types are reported
  sk <- testPairwiseEffects(list(ct = D1[meta$group != "B", meta$group != "B"]), des, meta, dist = "l2", n.permutations = 19)
  expect_equal(nrow(sk$results), 0); expect_equal(sk$skipped$celltype, "ct")
})

test_that("Freedman-Lane scheme runs on a shifted design with a continuous covariate", {
  sim <- simulateIndividualModel(c(A = 8, B = 8), p = 30, group.effect = 0.8 * sqrt(30), age = TRUE, age.effect = 0.2, seed = 12)
  meta <- sim$meta
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + age)
  r <- testPairwiseEffects(list(ct = sqrt(sim$D2)), des, meta, dist = "l2", n.permutations = 99, seed = 2)
  expect_equal(r$results$scheme, "freedman-lane"); expect_match(r$results$flags, "approximate")
  expect_lt(r$results$p.shift, 0.05)
  expect_error(testPairwiseEffects(list(ct = sqrt(sim$D2)), des, meta, dist = "l2", n.permutations = 9, permutation = "huh-jhun"))   # scheme removed
})

test_that("calibration: block scheme is exact under the null and detects dispersion and shift", {
  nrep <- 120
  pv <- replicate(nrep, {
    sim <- simulateIndividualModel(c(A = 6, B = 6, C = 4), p = 30, group.effect = 0, batch.levels = c("b1", "b2"), batch.effect = 0.4)
    des <- suppressMessages(cacoa:::buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = ~ group + batch))
    r <- testPairwiseEffects(list(ct = sqrt(sim$D2)), des, sim$meta, dist = "l2", n.permutations = 59)
    c(r$results$p.shift, r$results$p.var)
  })
  rej <- rowMeans(pv < 0.05)
  expect_true(rej[1] >= 0.005 && rej[1] <= 0.12)      # binomial 95% range around 0.05 at n = 120 is ~[0.01, 0.09]
  expect_true(rej[2] >= 0.005 && rej[2] <= 0.12)
  pw <- replicate(40, {
    sim <- simulateIndividualModel(c(A = 6, B = 6), p = 30, group.effect = 0, group.sd = c(A = 1, B = 2))
    des <- suppressMessages(cacoa:::buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = ~ group))
    r <- testPairwiseEffects(list(ct = sqrt(sim$D2)), des, sim$meta, dist = "l2", n.permutations = 59)
    c(r$results$p.shift, r$results$p.var)
  })
  expect_gt(mean(pw[2, ] < 0.05), 0.7)                 # var test has power against a dispersion change
  ps <- replicate(40, {
    sim <- simulateIndividualModel(c(A = 6, B = 6), p = 30, group.effect = 0.6 * sqrt(30))
    des <- suppressMessages(cacoa:::buildDesignMatrices(sim$meta, contrast = c("group", "B", "A"), formula = ~ group))
    testPairwiseEffects(list(ct = sqrt(sim$D2)), des, sim$meta, dist = "l2", n.permutations = 59)$results$p.shift
  })
  expect_gt(mean(ps < 0.05), 0.7)
})

test_that("distance model under the exchange of residuals: the shift F of every relabeling on the scale of the observed one", {
  set.seed(21); n <- 16
  meta <- data.frame(row.names = paste0("s", seq_len(n)), cond = factor(rep(c("a", "b"), n / 2)), batch = rep(c("b1", "b2", "b3", "b4"), each = 4), x = round(rnorm(n), 2), stringsAsFactors = FALSE)
  Y <- matrix(rnorm(n * 40), n) + outer(meta$x, rnorm(40)) + outer(as.integer(factor(meta$batch)), rnorm(40)); G <- tcrossprod(scale(Y, scale = FALSE))
  ref <- function(d, P, scale) {   # explicit hat matrices of every relabeled arrangement
    X <- d$F; cvec <- d$contrast.F[colnames(X)]; hi <- cacoa:::hatInfo(X); H <- hi$H; a <- drop(t(hi$A) %*% cvec); cXc <- drop(crossprod(cvec, hi$XtXi %*% cvec))
    Xr <- X %*% cacoa:::contrastNullBasis(cvec); R0 <- diag(n) - cacoa:::hatInfo(Xr)$H; K4 <- R0 %*% G %*% R0
    apply(P, 2, function(p) { ab <- numeric(n); ab[p] <- a; Hb <- matrix(0, n, n); Hb[p, p] <- H; Rb <- diag(n) - Hb; at <- drop(R0 %*% ab); cb <- sum(at^2) / sum(a^2)
      E <- Rb %*% R0; ss <- drop(t(ab) %*% K4 %*% ab) / cXc; rss <- sum(diag(Rb %*% K4 %*% Rb)); df <- sum(diag(E %*% t(E)))
      if (scale == "kernel") return(ss / (rss / (n - hi$rank)))
      if (scale == "independent") { k <- drop(Rb %*% at); L <- Rb - outer(k, ab) / sum(at^2); rss <- sum(diag(L %*% K4 %*% t(L))); df <- sum(diag(L %*% R0 %*% t(L))) }
      (ss / cb) / (rss / df) }) }
  kernelF <- function(d, P) { X <- d$F; cvec <- d$contrast.F[colnames(X)]; Xr <- X %*% cacoa:::contrastNullBasis(cvec); parts <- cacoa:::flGowerParts(G, Xr)
    apply(P, 2, function(p) { Gs <- cacoa:::flGowerPermute(parts, p); hi <- cacoa:::hatInfo(X); a <- drop(t(hi$A) %*% cvec); (drop(t(a) %*% Gs %*% a) / drop(crossprod(cvec, hi$XtXi %*% cvec))) / ((sum(diag(Gs)) - sum(hi$H * Gs)) / (n - hi$rank)) }) }
  d <- buildCacoaModel(meta, ~ cond + x + batch, test = "cond")$tests[[1]]$design; mp <- modelPermutations(d, n.permutations = 60, seed = 2); expect_equal(mp$plan$scheme, "freedman-lane")
  X <- d$F; cvec <- d$contrast.F[colnames(X)]; Xr <- X %*% cacoa:::contrastNullBasis(cvec)
  expect_equal(ref(d, mp$P, "kernel"), kernelF(d, mp$P), tolerance = 1e-10)                       # the reference reproduces the unscaled F
  for (sc in c("relabeling", "independent")) expect_equal(cacoa:::flShiftF(G, X, cvec, Xr, mp$P, sc), ref(d, mp$P, sc), tolerance = 1e-9)
  expect_gt(max(abs(cacoa:::flShiftF(G, X, cvec, Xr, mp$P, "relabeling") / kernelF(d, mp$P) - 1)), 1e-3)
  # the identity is the observed F; a fitted factor alone (exchange within its levels) is unchanged
  id <- matrix(seq_len(n), ncol = 1); for (sc in c("relabeling", "independent")) expect_equal(cacoa:::flShiftF(G, X, cvec, Xr, id, sc), kernelF(d, id), tolerance = 1e-10)
  d0 <- buildCacoaModel(meta, ~ cond + batch, test = "cond")$tests[[1]]$design; m0 <- modelPermutations(d0, n.permutations = 40, seed = 2, scheme = "freedman-lane")
  X0 <- d0$F; c0 <- d0$contrast.F[colnames(X0)]; Xr0 <- X0 %*% cacoa:::contrastNullBasis(c0)
  for (sc in c("relabeling", "independent")) expect_equal(cacoa:::flShiftF(G, X0, c0, Xr0, m0$P, sc), kernelF(d0, m0$P), tolerance = 1e-9)
  # the option reaches the test of a cell type: "average" leaves the kernel's F
  eff <- list(G = G, X = X, Z = cbind(1, as.numeric(meta$cond == "b")), contrast = cvec, z.end = list(num = c(1, 1), den = c(1, 0)))
  st <- cacoa:::permutationStatsForCellType(eff, mp$plan, mp$P); expect_equal(unname(st$perm[, "F"]), cacoa:::flShiftF(G, X, cvec, Xr, mp$P, "independent"), tolerance = 1e-9)
  op <- options(cacoa.residual.scale = "average"); on.exit(options(op))
  expect_equal(unname(cacoa:::permutationStatsForCellType(eff, mp$plan, mp$P)$perm[, "F"]), kernelF(d, mp$P), tolerance = 1e-8)
})

test_that("dispersion change: not tested when every exchange cell holding both labels is a fitted pair (mirror-image residuals)", {
  set.seed(8); nd <- 10
  # partly paired: donors 1-5 hold one sample of each label, donors 6-10 two samples of one label; x a continuous covariate
  meta <- data.frame(row.names = paste0("s", 1:(2 * nd)), donor = rep(paste0("d", 1:nd), each = 2), lab = factor(c(rep(c("A", "B"), 5), rep(c("A", "A", "B", "B"), length.out = 10))), x = round(rnorm(2 * nd), 2), stringsAsFactors = FALSE)
  Y <- matrix(rnorm(2 * nd * 30), 2 * nd) + matrix(rnorm(nd * 30, 0, 2), nd)[rep(1:nd, each = 2), ]; D <- as.matrix(dist(Y)); dimnames(D) <- list(rownames(meta), rownames(meta))
  run <- function(f) { m <- buildCacoaModel(meta, f, test = "lab"); testPairwiseEffects(list(ct = D), m$tests[[1]]$design, meta, dist = "l2", n.permutations = 99, seed = 1, block.vars = m$block.vars)$results }
  # the donor fitted: mirrored pairs only -> no dispersion test, flagged, shift still tested
  for (f in list(~ lab + donor, ~ lab + x + donor)) { r <- run(f); expect_true(is.na(r$p.var)); expect_true(is.na(r$var)); expect_true(is.na(r$total)); expect_match(r$flags, "dispersion-not-identifiable"); expect_true(is.finite(r$p.shift)) }
  # the residual projection gives mirror images for a fitted pair, with or without the covariate
  X <- buildCacoaModel(meta, ~ lab + x + donor, test = "lab")$tests[[1]]$design$F; R <- diag(nrow(X)) - cacoa:::hatInfo(X)$H
  expect_lt(max(abs(R[1, ] + R[2, ])), 1e-10)
  # the donor as a grouping term: residual distances are not mirrored, the dispersion change is tested
  for (f in list(~ lab + (1 | donor), ~ lab + x + (1 | donor))) expect_true(is.finite(run(f)$p.var))
  # a donor with three samples (two of one label, one of the other) carries information: tested
  meta3 <- rbind(meta, data.frame(row.names = "s21", donor = "d1", lab = factor("A", levels = c("A", "B")), x = 0.3))
  Y3 <- rbind(Y, Y[1, ] + rnorm(30, 0, 0.5)); D3 <- as.matrix(dist(Y3)); dimnames(D3) <- list(rownames(meta3), rownames(meta3))
  m3 <- buildCacoaModel(meta3, ~ lab + donor, test = "lab"); r3 <- testPairwiseEffects(list(ct = D3), m3$tests[[1]]$design, meta3, dist = "l2", n.permutations = 99, seed = 1)$results
  expect_true(is.finite(r3$p.var))
  # unpaired samples: unaffected
  m0 <- buildCacoaModel(meta, ~ lab + x, test = "lab"); expect_true(is.finite(testPairwiseEffects(list(ct = D), m0$tests[[1]]$design, meta, dist = "l2", n.permutations = 99, seed = 1)$results$p.var))
})
