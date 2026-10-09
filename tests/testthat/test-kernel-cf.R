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
