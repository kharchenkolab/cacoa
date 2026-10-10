# The per-column fitter (fit_and_randomize): R-drawn permutations, NA modes, robust fits, nuisance residualization
# (Freedman-Lane on all samples) and the studentized statistic.

test_that("perm_matrix: OLS statistics equal brute-force refits on the permuted design (complete data)", {
  meta <- gridMeta(levels = 2, n.per.level = 6); d <- gridDesign(meta); Y <- gridY(meta, m = 4)
  mp <- modelPermutations(d, scheme = "block", n.permutations = 25, seed = 3)
  expect_equal(dim(mp$P), c(12, 25)); expect_length(mp$cells, 2)                        # swaps within each batch
  f <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, n_randomizations = 0,
                                 return_sampled_stats = TRUE, perm_matrix = mp$P, statistic = "coef")
  expect_equal(dim(f$sampled_stats), c(25, 4))
  for (j in 1:4) for (b in c(1, 7, 25)) expect_equal(f$sampled_stats[b, j], refContrastStat(d$F, Y[, j], d$contrast.F, mp$P[, b]), tolerance = 1e-10)
  expect_equal(f$stat, f$effect)                                                          # statistic = "coef": the raw contrast
  # studentized statistic: brute-force t on the permuted design; the observed one equals lm's t
  ft <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P)
  for (j in 1:4) for (b in c(1, 7, 25)) expect_equal(ft$sampled_stats[b, j], refContrastT(d$F, Y[, j], d$contrast.F, mp$P[, b]), tolerance = 1e-10)
  expect_equal(ft$sampled_effects, f$sampled_stats, tolerance = 1e-10)                  # the raw effects travel along
  lmfit <- summary(lm(Y[, 1] ~ group + batch, meta))$coefficients["groupG2", ]
  expect_equal(unname(c(ft$effect[1], ft$se[1], ft$stat[1])), unname(lmfit[1:3]), tolerance = 1e-10); expect_equal(ft$df[1], 9L)
  # identity column reproduces the observed statistic
  P1 <- cbind(seq_len(12)); f1 <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = P1)
  expect_equal(as.numeric(f1$sampled_stats[1, ]), as.numeric(f1$stat), tolerance = 1e-12)
  # the permutations are shared by all columns: permuting the columns of P permutes the rows of sampled_stats
  f2 <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P[, 25:1])
  expect_equal(f2$sampled_stats, ft$sampled_stats[25:1, ], tolerance = 1e-12)
})

test_that("perm_matrix with missing data: drop induces P on the observed rows per NA pattern; impute_weak uses all rows", {
  meta <- gridMeta(levels = 2, n.per.level = 6); d <- gridDesign(meta); Y <- gridY(meta, m = 5, na = "pattern")
  mp <- modelPermutations(d, scheme = "block", n.permutations = 30, seed = 4)
  f <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, na_mode = "drop", statistic = "coef")
  ft <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, na_mode = "drop")
  for (j in c(2, 3)) {
    obs <- which(!is.na(Y[, j]))
    sub.plan <- permutationPlan(d, meta, rownames(meta)[obs], scheme = "block", n.permutations = 30)
    for (b in c(1, 11, 30)) {
      q <- cacoa:::inducePermutation(mp$P[, b], obs, sub.plan)
      expect_equal(f$sampled_stats[b, j], refContrastStat(d$F[obs, , drop = FALSE], Y[obs, j], d$contrast.F, q), tolerance = 1e-10)
      expect_equal(ft$sampled_stats[b, j], refContrastT(d$F[obs, , drop = FALSE], Y[obs, j], d$contrast.F, q), tolerance = 1e-10)
    }
    expect_equal(ft$df[j], length(obs) - ncol(d$F))
  }
  fw <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, na_mode = "impute_weak", na_weight = 1e-4, statistic = "coef")
  fwt <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, na_mode = "impute_weak", na_weight = 1e-4)
  j <- 3; w <- ifelse(is.na(Y[, j]), 1e-4, 1); yc <- Y[, j]; yc[is.na(yc)] <- mean(yc, na.rm = TRUE)
  for (b in c(2, 30)) {
    expect_equal(fw$sampled_stats[b, j], refContrastStat(d$F, yc, d$contrast.F, mp$P[, b], w = w), tolerance = 1e-8)
    expect_equal(fwt$sampled_stats[b, j], refContrastT(d$F, yc, d$contrast.F, mp$P[, b], w = w, df = sum(w == 1) - ncol(d$F)), tolerance = 1e-8)
  }
})

test_that("perm_matrix: robust paths and the nuisance (Freedman-Lane) path consume the same permutations", {
  meta <- gridMeta(levels = 2, n.per.level = 6, age = TRUE); Y <- gridY(meta, m = 3, na = "pattern"); Y[2, 1] <- 8   # an outlier
  d <- gridDesign(meta, formula = ~ group + batch)
  mp <- modelPermutations(d, scheme = "block", n.permutations = 20, seed = 5)
  fh <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, robust = "huber")
  fh2 <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P[, 20:1], robust = "huber")
  expect_equal(fh2$sampled_stats, fh$sampled_stats[20:1, ], tolerance = 1e-12)
  P1 <- cbind(seq_len(12)); fh1 <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = P1, robust = "winsor")
  expect_equal(as.numeric(fh1$sampled_stats[1, ]), as.numeric(fh1$stat), tolerance = 1e-12)
  # an empty Z is the same as no Z
  fb <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P)
  fz <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, Z = matrix(0, 12, 0), perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P)
  expect_equal(fz$sampled_stats, fb$sampled_stats, tolerance = 1e-12)
  # with a continuous nuisance column (Z): every sample stays in the fit; the observed fit equals lm on the full model
  dX <- suppressMessages(buildDesignMatrices(meta, contrast = c("group", "G2", "G1"), formula = ~ group + batch + age))
  mpf <- modelPermutations(dX, scheme = "freedman-lane", n.permutations = 15, seed = 6)
  expect_true(all(mpf$plan$in.set))
  g1 <- cacoa:::fit_and_randomize(X = dX$X, Y = Y, contrast = dX$contrast.X, Z = dX$Z, perm_groups = mpf$cells, return_sampled_stats = TRUE, perm_matrix = mpf$P)
  g2 <- cacoa:::fit_and_randomize(X = dX$X, Y = Y, contrast = dX$contrast.X, Z = dX$Z, perm_groups = mpf$cells, return_sampled_stats = TRUE, perm_matrix = mpf$P[, 15:1])
  expect_equal(g2$sampled_stats, g1$sampled_stats[15:1, ], tolerance = 1e-12)
  gid <- cacoa:::fit_and_randomize(X = dX$X, Y = Y, contrast = dX$contrast.X, Z = dX$Z, return_sampled_stats = TRUE, perm_matrix = cbind(seq_len(12)))
  expect_equal(as.numeric(gid$sampled_stats[1, ]), as.numeric(gid$stat), tolerance = 1e-10)
  for (j in 1:3) {
    obs <- which(!is.na(Y[, j])); lmj <- summary(lm(Y[obs, j] ~ group + batch + age, meta[obs, ]))$coefficients["groupG2", ]
    expect_equal(unname(c(g1$effect[j], g1$se[j], g1$stat[j])), unname(lmj[1:3]), tolerance = 1e-9)
    expect_equal(g1$df[j], length(obs) - 4L)
    expect_equal(sum(is.finite(g1$residuals[, j])), length(obs))                         # residuals for every observed sample
    for (b in c(1, 8, 15)) {                                                               # permuted: brute-force Freedman-Lane on the observed rows
      sub.plan <- permutationPlan(dX, meta, rownames(meta)[obs], scheme = "freedman-lane", n.permutations = 15)
      q <- cacoa:::inducePermutation(mpf$P[, b], obs, sub.plan)
      expect_equal(g1$sampled_stats[b, j], refFLT(dX$X[obs, , drop = FALSE], dX$Z[obs, , drop = FALSE], Y[obs, j], dX$contrast.X, q), tolerance = 1e-9)
    }
  }
  # weak imputation with Z: the design rows are permuted, the weights stay with the response rows
  gw <- cacoa:::fit_and_randomize(X = dX$X, Y = Y, contrast = dX$contrast.X, Z = dX$Z, perm_groups = mpf$cells, return_sampled_stats = TRUE, perm_matrix = mpf$P, na_mode = "impute_weak", statistic = "coef")
  j <- 2; miss <- is.na(Y[, j]); w <- ifelse(miss, 1e-4, 1)
  yc <- Y[, j]; yc[miss] <- mean(yc[!miss]); Xr <- qr.resid(qr(dX$Z), dX$X); yr <- qr.resid(qr(dX$Z), yc); yr[miss] <- mean(yr[!miss])
  for (b in c(3, 15)) expect_equal(gw$sampled_stats[b, j], refContrastStat(Xr, yr, dX$contrast.X, mpf$P[, b], w = w), tolerance = 1e-8)
  # robust fits with Z run on every NA mode and keep the structure
  for (rb in c("huber", "winsor")) for (na in c("drop", "impute_weak")) {
    r1 <- cacoa:::fit_and_randomize(X = dX$X, Y = Y, contrast = dX$contrast.X, Z = dX$Z, perm_groups = mpf$cells, return_sampled_stats = TRUE, perm_matrix = mpf$P, robust = rb, na_mode = na)
    r2 <- cacoa:::fit_and_randomize(X = dX$X, Y = Y, contrast = dX$contrast.X, Z = dX$Z, perm_groups = mpf$cells, return_sampled_stats = TRUE, perm_matrix = mpf$P[, 15:1], robust = rb, na_mode = na)
    rid <- cacoa:::fit_and_randomize(X = dX$X, Y = Y, contrast = dX$contrast.X, Z = dX$Z, return_sampled_stats = TRUE, perm_matrix = cbind(seq_len(12)), robust = rb, na_mode = na)
    expect_true(all(is.finite(r1$stat)) && all(is.finite(r1$se)), info = paste(rb, na))
    expect_equal(r2$sampled_stats, r1$sampled_stats[15:1, ], tolerance = 1e-12, info = paste(rb, na))
    expect_equal(as.numeric(rid$sampled_stats[1, ]), as.numeric(rid$stat), tolerance = 1e-8, info = paste(rb, na))
  }
})

test_that("Freedman-Lane fits all samples: a two-level test with a numeric covariate matches lm and has power", {
  # regression: the former fitter kept only the rows with a non-zero contrast weight, i.e. the target level alone
  set.seed(2); n <- 12
  meta <- data.frame(Group = factor(rep(c("G1", "G2"), each = n / 2)), age = round(rnorm(n, 50, 5)), row.names = paste0("s", 1:n))
  m <- buildCacoaModel(meta, ~ Group + age, test = "Group"); d <- m$tests[[1]]$design
  expect_equal(m$tests[[1]]$permutation$scheme, "freedman-lane")
  Y <- cbind(eff = rnorm(n) + 2 * (meta$Group == "G2") + 0.1 * meta$age, null = rnorm(n) + 0.1 * meta$age)
  r <- cacoa:::performLMPermutations(d, Y, perm.method = "freedman-lane", n.permutations = 999, seed = 1)
  lmc <- sapply(1:2, function(j) summary(lm(Y[, j] ~ Group + age, meta))$coefficients["GroupG2", ])
  expect_equal(unname(r$effect), unname(lmc[1, ]), tolerance = 1e-10)
  expect_equal(unname(r$se), unname(lmc[2, ]), tolerance = 1e-10)
  expect_equal(unname(r$stat.obs), unname(lmc[3, ]), tolerance = 1e-10)
  expect_equal(unname(r$df), c(9L, 9L))
  expect_lt(r$pval[["eff"]], 0.01); expect_gt(r$pval[["null"]], 0.2)
  expect_equal(nrow(r$residuals), n)                                                      # residuals for all samples
  expect_equal(rownames(r$coef), colnames(d$X))
  expect_equal(dim(r$effects.perm), dim(r$stats.perm)); expect_equal(r$statistic, "t")
  # three-level factor, contrast of two levels: all three levels contribute to the pooled residual variance
  meta3 <- data.frame(Group = factor(rep(c("G1", "G2", "G3"), each = 4)), age = round(rnorm(12, 50, 5)), row.names = paste0("s", 1:12))
  m3 <- buildCacoaModel(meta3, ~ Group + age, test = "Group: G2 vs G1"); d3 <- m3$tests[[1]]$design
  y3 <- rnorm(12) + 0.1 * meta3$age
  r3 <- cacoa:::performLMPermutations(d3, y3, perm.method = "freedman-lane", n.permutations = 99, seed = 1)
  lm3 <- summary(lm(y3 ~ Group + age, meta3))$coefficients["GroupG2", ]
  expect_equal(unname(c(r3$effect, r3$se, r3$stat.obs)), unname(lm3[1:3]), tolerance = 1e-10); expect_equal(unname(r3$df), 8L)
  expect_true(any(r3$P[meta3$Group == "G3", ] != which(meta3$Group == "G3")))              # Freedman-Lane permutes the residuals of every sample
  # statistic = "coef" keeps the raw contrast as the statistic
  rc <- cacoa:::performLMPermutations(d, Y, perm.method = "freedman-lane", n.permutations = 99, seed = 1, statistic = "coef")
  expect_equal(rc$stat.obs, rc$effect); expect_equal(rc$stats.perm, rc$effects.perm); expect_equal(rc$statistic, "coef")
})

test_that("t and coef statistics: identical p-values for two groups without covariates, t is calibrated and more powerful otherwise", {
  set.seed(3); n <- 12
  meta <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)), row.names = paste0("s", 1:n))
  d <- suppressMessages(buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group))
  mp <- modelPermutations(d, scheme = "block", n.permutations = 199, seed = 1)
  Y <- matrix(rnorm(n * 50), n, 50) + 1.2 * (meta$group == "B")
  pt <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, perm_matrix = mp$P)$p_value
  pc <- cacoa:::fit_and_randomize(X = d$F, Y = Y, contrast = d$contrast.F, perm_groups = mp$cells, perm_matrix = mp$P, statistic = "coef")$p_value
  expect_true(all(abs(pt - pc) <= 1 / 200 + 1e-12))                                        # monotone-equivalent statistics (ties at rounding level may differ by one count)
  # a numeric covariate permuted along with the labels (block scheme): the raw coefficient is not pivotal
  meta$age <- round(rnorm(n, 50, 5))
  d2 <- suppressMessages(buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + age))
  mp2 <- modelPermutations(d2, scheme = "block", n.permutations = 199, seed = 2)
  Y0 <- matrix(rnorm(n * 400), n, 400) + 0.1 * meta$age                                   # null columns
  Y1 <- Y0 + 1.5 * (meta$group == "B")                                                       # planted effect
  for (st in c("t", "coef")) {
    p0 <- cacoa:::fit_and_randomize(X = d2$F, Y = Y0, contrast = d2$contrast.F, perm_groups = mp2$cells, perm_matrix = mp2$P, statistic = st)$p_value
    expect_true(mean(p0 < 0.05) > 0.02 && mean(p0 < 0.05) < 0.09, info = st)              # calibrated (400 columns)
  }
  pow <- sapply(c("t", "coef"), function(st) mean(cacoa:::fit_and_randomize(X = d2$F, Y = Y1, contrast = d2$contrast.F, perm_groups = mp2$cells, perm_matrix = mp2$P, statistic = st)$p_value < 0.05))
  expect_gt(pow[["t"]], pow[["coef"]] + 0.1)
  # Freedman-Lane with Z, robust fits and weak imputation stay calibrated under the null
  mpf <- modelPermutations(d2, scheme = "freedman-lane", n.permutations = 199, seed = 3)
  Yna <- Y0; Yna[1, 1:200] <- NA
  for (rb in c("none", "huber", "winsor")) for (na in c("drop", "impute_weak")) {
    p0 <- cacoa:::fit_and_randomize(X = d2$X, Y = Yna, contrast = d2$contrast.X, Z = d2$Z, perm_groups = mpf$cells, perm_matrix = mpf$P, robust = rb, na_mode = na)$p_value
    expect_true(mean(p0 < 0.05) > 0.02 && mean(p0 < 0.05) < 0.09, info = paste(rb, na))
  }
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
