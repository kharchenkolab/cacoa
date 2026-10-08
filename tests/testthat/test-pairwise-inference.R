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

test_that("Freedman-Lane and Huh-Jhun schemes run and agree with block on a shifted design", {
  sim <- simulateIndividualModel(c(A = 8, B = 8), p = 30, group.effect = 0.8 * sqrt(30), age = TRUE, age.effect = 0.2, seed = 12)
  meta <- sim$meta
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + age)
  r <- testPairwiseEffects(list(ct = sqrt(sim$D2)), des, meta, dist = "l2", n.permutations = 99, seed = 2)
  expect_equal(r$results$scheme, "freedman-lane"); expect_match(r$results$flags, "approximate")
  expect_lt(r$results$p.shift, 0.05)
  rh <- testPairwiseEffects(list(ct = sqrt(sim$D2)), des, meta, dist = "l2", n.permutations = 99, seed = 2, permutation = "huh-jhun")
  expect_lt(rh$results$p.shift, 0.05)
  expect_true(is.na(rh$results$p.var))
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
