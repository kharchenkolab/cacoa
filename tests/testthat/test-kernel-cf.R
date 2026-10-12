# Step 4: batched cluster-free kernel against the former per-cell R loop (helper-reference.R).

cfInputs <- function(cao, nns, permutation = "auto", B = 29, seed = 1, min.n.obs = 3) {
  cm <- Matrix::t(cao$data.object); samples <- levels(cao$sample.per.cell); n <- length(samples)
  spc <- factor(as.character(cao$sample.per.cell[colnames(cm)]), levels = samples)
  pairs <- t(utils::combn(n, 2)); storage.mode(pairs) <- "integer"
  Y <- cacoa:::estimateExpressionShiftsPairsLM(cm = cm, sample_per_cell = as.integer(spc), nn_ids = lapply(nns, as.integer), pairs_mat = pairs, min_n_obs_per_samp = min.n.obs, dist = "cor", log_vecs = TRUE)
  des <- cao$model; des$F <- cao$model$F[samples, , drop = FALSE]
  gplan <- cacoa:::permutationPlan(des, cao$sample.meta[samples, ], samples, scheme = permutation, n.permutations = B, max.enumerate = 0)
  P <- cacoa:::withSeed(seed, drawPermutations(gplan, B)); storage.mode(P) <- "integer"
  list(Y = Y, pairs = pairs, samples = samples, gplan = gplan, P = P, meta = cao$sample.meta[samples, ], cm = cm)
}

test_that("batched kernel equals the R reference: whole-type and random neighbourhoods, block and Freedman-Lane, missing samples", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.2)
  cg <- cao$cell.groups; idx <- split(seq_along(cg) - 1L, cg)
  nns.type <- lapply(as.character(cg), function(t) idx[[t]]); names(nns.type) <- names(cg)
  set.seed(3); nns.rand <- lapply(seq_along(cg), function(i) sample(seq_along(cg) - 1L, 90)); names(nns.rand) <- names(cg)   # some samples fall below min.n.obs
  for (nns in list(nns.type, nns.rand)) {
    inp <- cfInputs(cao, nns)
    expect_true(any(is.na(inp$Y)) || identical(nns, nns.type))
    ref <- refClusterFreeShifts(inp$Y, inp$pairs, inp$samples, cao$model, inp$meta, inp$gplan, inp$P)
    res <- clusterFreeExpressionShifts(inp$cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, n.permutations = 29, seed = 1, adjust = TRUE, smooth = FALSE)
    expect_equal(unname(res$stat), ref$F, tolerance = 1e-9); expect_equal(unname(res$p.value), ref$p, tolerance = 1e-12)
    expect_equal(unname(res$shifts), ref$shift, tolerance = 1e-9); expect_equal(unname(res$n.samples), ref$n)
    expect_equal(unname(res$z.score), ref$z, tolerance = 1e-9)
    expect_true(sum(is.finite(res$stat)) == sum(is.finite(ref$F)))
    expect_true(all(is.finite(res$z.adj[is.finite(res$z.score)])))
  }
  # Freedman-Lane (continuous nuisance covariate) and a design with a third level
  cao$sample.meta$age <- c(41, 55, 62, 48, 70, 39, 58, 51, 66, 44)
  cao$setModel(~ group + age, test = "group", verbose = FALSE)
  inp <- cfInputs(cao, nns.rand, permutation = "freedman-lane")
  expect_equal(inp$gplan$scheme, "freedman-lane")
  ref <- refClusterFreeShifts(inp$Y, inp$pairs, inp$samples, cao$model, inp$meta, inp$gplan, inp$P)
  res <- clusterFreeExpressionShifts(inp$cm, cao$sample.per.cell, nns.rand, cao$model, cao$sample.meta, permutation = "freedman-lane", n.permutations = 29, seed = 1, adjust = FALSE, smooth = FALSE)
  expect_equal(unname(res$stat), ref$F, tolerance = 1e-9); expect_equal(unname(res$p.value), ref$p); expect_equal(unname(res$z.score), ref$z, tolerance = 1e-9)
  # results do not depend on the number of threads
  res4 <- clusterFreeExpressionShifts(inp$cm, cao$sample.per.cell, nns.rand, cao$model, cao$sample.meta, permutation = "freedman-lane", n.permutations = 29, seed = 1, adjust = TRUE, smooth = TRUE, n.cores = 4)
  res1 <- clusterFreeExpressionShifts(inp$cm, cao$sample.per.cell, nns.rand, cao$model, cao$sample.meta, permutation = "freedman-lane", n.permutations = 29, seed = 1, adjust = TRUE, smooth = TRUE, n.cores = 1)
  expect_equal(res4$stat, res1$stat); expect_equal(res4$z.adj, res1$z.adj); expect_equal(res4$shifts.smoothed, res1$shifts.smoothed)
})

test_that("streaming kernel equals the batched kernel on a precomputed distance matrix; 0-based neighbourhoods are read as such", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.2)
  cg <- cao$cell.groups; idx <- split(seq_along(cg) - 1L, cg)
  nns <- lapply(as.character(cg), function(t) idx[[t]]); names(nns) <- names(cg)
  inp <- cfInputs(cao, nns)
  lev <- rep(-1L, length(inp$samples))
  bat <- cacoa:::cluster_free_shift_batch(inp$Y, inp$pairs - 1L, length(inp$samples), cao$model$F[inp$samples, ], cao$model$contrast.F, lev, 2L,
                                          as.integer(inp$gplan$strata), as.integer(inp$gplan$in.set), inp$P, FALSE, TRUE, 1L)
  res <- clusterFreeExpressionShifts(inp$cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, n.permutations = 29, seed = 1, adjust = FALSE, smooth = FALSE, n.cores = 3)
  expect_equal(unname(res$stat), as.numeric(bat$stat), tolerance = 1e-10); expect_equal(unname(res$p.value), as.numeric(bat$p))
  expect_equal(unname(res$shifts), as.numeric(bat$shift), tolerance = 1e-10)
  # every whole-type neighbourhood reproduces the cell-type engine's F on the same mean profiles (regression: the
  # distance builder used to shift 0-based neighbourhoods without cell 0 by one cell, so only the first type matched)
  for (t in levels(cg)) {
    prof <- sapply(levels(cao$sample.per.cell), function(s) Matrix::rowMeans(inp$cm[, cg == t & cao$sample.per.cell == s, drop = FALSE]))
    D <- 1 - cor(log10(1e3 * prof + 1))
    r <- testPairwiseEffects(list(x = D), cao$model, cao$sample.meta, dist = "cor", n.permutations = 9, seed = 1)
    expect_equal(unname(res$stat[cg == t][1]), r$results$F, tolerance = 1e-8, info = t)
  }
  # explicit index conventions of the pair-distance builder
  nb <- idx[[2]]                                                     # 0-based ids of the second type: none is 0
  Y0 <- cacoa:::estimateExpressionShiftsPairsLM(inp$cm, as.integer(cao$sample.per.cell), list(nb), inp$pairs, 1L, "cor", TRUE)
  Y1 <- cacoa:::estimateExpressionShiftsPairsLM(inp$cm, as.integer(cao$sample.per.cell), list(nb + 1L), inp$pairs, 1L, "cor", TRUE, nn_one_based = TRUE)
  expect_equal(Y0, Y1)
  expect_error(cacoa:::estimateExpressionShiftsPairsLM(inp$cm, as.integer(cao$sample.per.cell), list(nb + 1L), inp$pairs, 1L, "cor", TRUE), "out-of-range")   # no silent re-basing
})

test_that("a sample with too few cells in a neighbourhood drops itself, not the samples it is paired with", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.2)
  cg <- cao$cell.groups; spc <- cao$sample.per.cell
  thin <- levels(spc)[1]; drop <- which(spc == thin)[-(1:2)]                         # every neighbourhood keeps only two cells of the first sample
  nns <- lapply(seq_along(cg), function(i) setdiff(seq_along(cg) - 1L, drop - 1L)); names(nns) <- names(cg)
  res <- clusterFreeExpressionShifts(Matrix::t(cao$data.object), spc, nns, cao$model, cao$sample.meta, n.permutations = 19, seed = 1, adjust = FALSE, smooth = FALSE)
  expect_true(all(res$n.samples == nlevels(spc) - 1)); expect_true(all(is.finite(res$stat)))
  # the R helper applies the same rule
  n <- nlevels(spc); pairs <- t(utils::combn(n, 2)); y <- runif(nrow(pairs)); y[pairs[, 1] == 1 | pairs[, 2] == 1] <- NA
  expect_equal(rownames(cacoa:::pairVectorToMatrix(y, pairs, levels(spc))), levels(spc)[-1])
})

test_that("label permutations are global relabelings: a cell is skipped where the relabeling leaves too few samples, and pooled group means stay calibrated with absent samples", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.2)
  cg <- cao$cell.groups
  set.seed(5); nns <- lapply(seq_along(cg), function(i) sample(seq_along(cg) - 1L, 60)); names(nns) <- names(cg)   # several samples absent per cell
  res <- clusterFreeExpressionShifts(Matrix::t(cao$data.object), cao$sample.per.cell, nns, cao$model, cao$sample.meta, n.permutations = 49, seed = 1,
                                     adjust = TRUE, smooth = FALSE, groups = cg, min.samp.per.level = 3)
  ok <- is.finite(res$z.score)
  expect_true(any(res$n.samples[ok] < nlevels(cao$sample.per.cell)))
  expect_true(all(res$n.perm[ok] >= 1 & res$n.perm[ok] <= 49)); expect_true(any(res$n.perm[ok] < 49))
  pm <- attr(res$group.test, "perm.means"); expect_equal(dim(pm), c(49L, nlevels(cg)))
  for (g in seq_len(nlevels(cg))) {   # the group p is the exceedance over the relabelings under which cells of the group were tested
    pk <- is.finite(pm[, g]); mo <- res$group.test$mean.z[g]
    expect_equal(res$group.test$p.perm[g], (sum(pm[pk, g] >= mo - 1e-12) + 1) / (sum(pk) + 1))
  }
  # calibration of the group test with absent samples (batched kernel on iid distance matrices; the induced-within-cell
  # permutation of the former kernel rejected about 30 % of the groups at the 5 % level in this setting)
  set.seed(7); n <- 16; B <- 199; m <- 600
  pairs <- t(utils::combn(n, 2)); storage.mode(pairs) <- "integer"
  regions <- sample(3, m, replace = TRUE); base <- lapply(1:3, function(r) matrix(rnorm(n * 5), n))
  Y <- sapply(seq_len(m), function(i) { X <- base[[regions[i]]] + matrix(rnorm(n * 5, 0, 0.7), n); d <- as.matrix(dist(X))^2
    miss <- sample(n, 6); d[miss, ] <- NA; d[, miss] <- NA; d[cbind(pairs[, 1], pairs[, 2])] })
  meta0 <- data.frame(row.names = paste0("s", seq_len(n))); pg <- c()
  for (k in 1:20) {
    meta <- meta0; meta$cond <- factor(sample(rep(c("A", "B"), n / 2)))
    model <- buildCacoaModel(meta, formula = ~ cond, test = "cond")
    gp <- cacoa:::permutationPlan(model, model$meta, rownames(meta), scheme = "auto", n.permutations = B, max.enumerate = 0)
    P <- cacoa:::withSeed(k, drawPermutations(gp, B)); storage.mode(P) <- "integer"
    r <- cacoa:::cluster_free_shift_batch(Y, pairs - 1L, n, model$F, model$contrast.F, ifelse(meta$cond == "A", 1L, 2L), 2L, as.integer(gp$strata), as.integer(gp$in.set), P,
                                          FALSE, TRUE, 2L, 0L, 1.345, as.numeric(c()), as.integer(regions) - 1L)
    mo <- r$group_obs / r$group_n; pm <- r$group_sum / r$group_count; pm[r$group_count == 0] <- NA
    pg <- c(pg, sapply(1:3, function(g) { pk <- is.finite(pm[, g]); (sum(pm[pk, g] >= mo[g] - 1e-12) + 1) / (sum(pk) + 1) }))
  }
  expect_lt(mean(pg < 0.05), 0.2); expect_gt(mean(pg), 0.3)
})

test_that("cluster-free shifts return pooled per-cell p-values and accept p.source = 'pooled'", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  cg <- cao$cell.groups; idx <- split(seq_along(cg) - 1L, cg)
  set.seed(4); nns <- lapply(seq_along(cg), function(i) sample(idx[[as.character(cg[i])]], 60)); names(nns) <- names(cg)
  inp <- cfInputs(cao, nns)
  r <- clusterFreeExpressionShifts(inp$cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, n.permutations = 49, seed = 1, adjust = TRUE, smooth = FALSE,
                                   adjust.method = "bh", p.source = "pooled")
  ok <- is.finite(r$z.score)
  expect_true(all(is.finite(r$p.pooled[ok]))); expect_true(all(r$p.pooled[ok] > 0 & r$p.pooled[ok] <= 1))
  expect_true(all(is.na(r$p.pooled[!ok])))
  o <- order(r$z.score[ok]); expect_true(all(diff(r$p.pooled[ok][o]) <= 1e-12))    # non-increasing in z
  expect_true(min(r$p.pooled[ok]) < 1 / 50)                                        # below the per-cell floor
  expect_equal(r$settings$p.source, "pooled"); expect_true(all(is.finite(r$p.adj[ok])))
  cao$embedding <- matrix(rnorm(2 * length(cg)), ncol = 2, dimnames = list(names(cg), c("x", "y")))
  cao$estimateCellDensity(method = "kde", bins = 20, verbose = FALSE)
  d <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 49, adjust.method = "bh", p.source = "pooled", verbose = FALSE)$diff$permutation
  expect_equal(d$p.source, "pooled"); expect_true(all(is.finite(d$p.pooled))); expect_true(all(d$p.pooled > 0 & d$p.pooled <= 1))
  expect_gt(cor(d$p.pooled, d$p.perm, method = "spearman"), 0.4)        # both rank the statistic against permutation nulls (p.perm is coarse at 49 permutations)
  expect_true(min(d$p.pooled) < min(d$p.perm))                             # the pooled null resolves below the per-cell floor
})
