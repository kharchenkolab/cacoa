# Step 1 of the engine convergence: the per-column fitter consumes R-drawn permutations (perm_matrix).

test_that("perm_matrix: OLS statistics equal brute-force refits on the permuted design (complete data)", {
  meta <- gridMeta(levels = 2, n.per.level = 6); d <- gridDesign(meta); Y <- gridY(meta, m = 4)
  mp <- modelPermutations(d, scheme = "block", n.permutations = 25, seed = 3)
  expect_equal(dim(mp$P), c(12, 25)); expect_length(mp$cells, 2)                        # swaps within each batch
  f <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, n_randomizations = 0,
                                 return_sampled_stats = TRUE, perm_matrix = mp$P)
  expect_equal(dim(f$sampled_stats), c(25, 4))
  for (j in 1:4) for (b in c(1, 7, 25)) expect_equal(f$sampled_stats[b, j], refContrastStat(d$F, Y[, j], d$contrast.F, mp$P[, b]), tolerance = 1e-10)
  # identity column reproduces the observed statistic
  P1 <- cbind(seq_len(12)); f1 <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = P1)
  expect_equal(as.numeric(f1$sampled_stats[1, ]), as.numeric(f1$stat), tolerance = 1e-12)
  # the permutations are shared by all columns: permuting the columns of P permutes the rows of sampled_stats
  f2 <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P[, 25:1])
  expect_equal(f2$sampled_stats, f$sampled_stats[25:1, ], tolerance = 1e-12)
})

test_that("perm_matrix with missing data: drop induces P on the observed rows per NA pattern; impute_weak uses all rows", {
  meta <- gridMeta(levels = 2, n.per.level = 6); d <- gridDesign(meta); Y <- gridY(meta, m = 5, na = "pattern")
  mp <- modelPermutations(d, scheme = "block", n.permutations = 30, seed = 4)
  f <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, na_mode = "drop")
  for (j in c(2, 3)) {
    obs <- which(!is.na(Y[, j]))
    sub.plan <- permutationPlan(d, meta, rownames(meta)[obs], scheme = "block", n.permutations = 30)
    for (b in c(1, 11, 30)) {
      q <- cacoa:::inducePermutation(mp$P[, b], obs, sub.plan)
      expect_equal(f$sampled_stats[b, j], refContrastStat(d$F[obs, , drop = FALSE], Y[obs, j], d$contrast.F, q), tolerance = 1e-10)
    }
  }
  fw <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, na_mode = "impute_weak", na_weight = 1e-4)
  j <- 3; w <- ifelse(is.na(Y[, j]), 1e-4, 1); yc <- Y[, j]; yc[is.na(yc)] <- mean(yc, na.rm = TRUE)
  for (b in c(2, 30)) expect_equal(fw$sampled_stats[b, j], refContrastStat(d$F, yc, d$contrast.F, mp$P[, b], w = w), tolerance = 1e-8)
})

test_that("perm_matrix: robust path and Freedman-Lane path consume the same permutations", {
  meta <- gridMeta(levels = 2, n.per.level = 6, age = TRUE); Y <- gridY(meta, m = 3, na = "pattern"); Y[2, 1] <- 8   # an outlier
  d <- gridDesign(meta, formula = ~ group + batch)
  mp <- modelPermutations(d, scheme = "block", n.permutations = 20, seed = 5)
  fh <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, robust = "huber")
  fh2 <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P[, 20:1], robust = "huber")
  expect_equal(fh2$sampled_stats, fh$sampled_stats[20:1, ], tolerance = 1e-12)
  P1 <- cbind(seq_len(12)); fh1 <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = P1, robust = "winsor")
  expect_equal(as.numeric(fh1$sampled_stats[1, ]), as.numeric(fh1$stat), tolerance = 1e-12)
  # Freedman-Lane with an empty Z (fast path) equals the block fitter on the same P
  fb <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P)
  ff <- cacoa:::fl_fwl_cpp(X = d$F, Z = matrix(0, 0, 0), Y = Y, contrast = d$contrast.F, core_perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P)
  expect_equal(ff$sampled_stats, fb$sampled_stats, tolerance = 1e-12)
  # with a continuous nuisance column: identity reproduces the observed statistic, reversed P reverses the rows
  dX <- suppressMessages(buildDesignMatrices(meta, contrast = c("group", "G2", "G1"), formula = ~ group + batch + age))
  mpf <- modelPermutations(dX, scheme = "freedman-lane", n.permutations = 15, seed = 6)
  expect_true(all(mpf$plan$in.set))
  g1 <- cacoa:::fl_fwl_cpp(X = dX$X, Z = dX$Z, Y = Y, contrast = dX$contrast.X, core_rows = dX$core.rows, core_perm_groups = cacoa:::coreCells(mpf$cells, dX$core.rows), return_sampled_stats = TRUE, perm_matrix = mpf$P)
  g2 <- cacoa:::fl_fwl_cpp(X = dX$X, Z = dX$Z, Y = Y, contrast = dX$contrast.X, core_rows = dX$core.rows, core_perm_groups = cacoa:::coreCells(mpf$cells, dX$core.rows), return_sampled_stats = TRUE, perm_matrix = mpf$P[, 15:1])
  expect_equal(g2$sampled_stats, g1$sampled_stats[15:1, ], tolerance = 1e-12)
  gid <- cacoa:::fl_fwl_cpp(X = dX$X, Z = dX$Z, Y = Y, contrast = dX$contrast.X, core_rows = dX$core.rows, return_sampled_stats = TRUE, perm_matrix = cbind(seq_len(12)))
  expect_equal(as.numeric(gid$sampled_stats[1, ]), as.numeric(gid$stat), tolerance = 1e-10)
})

test_that("performLMPermutations draws from the model's plan: exhaustive small designs, shared P across analyses", {
  meta <- gridMeta(levels = 2, n.per.level = 2, batch = FALSE); d <- gridDesign(meta); Y <- gridY(meta, m = 3)
  r <- cacoa:::performLMPermutations(d, Y, perm.method = "block", n.permutations = 99, seed = 1)
  expect_equal(ncol(r$P), 6); expect_equal(nrow(r$stats.perm), 6)                      # 2 vs 2: six distinct relabelings
  expect_true(all(r$pval >= 1 / 6 - 1e-12))
  meta2 <- gridMeta(levels = 4, n.per.level = 5); d2 <- gridDesign(meta2); Y2 <- gridY(meta2, m = 2)
  r2 <- cacoa:::performLMPermutations(d2, Y2, perm.method = "block", n.permutations = 40, seed = 7)
  mp <- modelPermutations(d2, scheme = "block", n.permutations = 40, seed = 7)
  expect_equal(r2$P, mp$P)                                                             # same seed, same P as the shift engine would use
  moved <- apply(mp$P, 2, function(p) unique(as.character(meta2$group[p != seq_len(20)])))
  expect_true(all(unlist(moved) %in% c("G1", "G2")))                                   # only the compared levels move
  # the whole toy pipeline: composition and shifts on one object share the permutations
  skip_if_not_installed("coda.base"); skip_if_not_installed("psych")
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 60, n.genes = 80, n.cell.types = 4, shift = 1.5)
  cao$setOptions(seed = 11)
  coda <- cao$estimateCellLoadings(n.permutations = 49, verbose = FALSE)
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 49, verbose = FALSE)
  expect_equal(dim(coda$fit$P), c(10, 49))
  expect_equal(unname(coda$fit$P), unname(modelPermutations(cao$model, scheme = "block", n.permutations = 49, seed = 11)$P))
  expect_true(all(c("P", "perm.cells") %in% names(coda$fit)))
})
