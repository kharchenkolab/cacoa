# Normalized, noise-weighted profiles of the cluster-free shift map (profile = "normalized").

toyCurves <- function(scale, phi = 0.2, tau2 = 1e-3) {
  mu.grid <- seq(-7, -1, by = 0.05)
  list(mu.grid = mu.grid, phi = rep(phi, length(mu.grid)), tau2 = rep(tau2, length(mu.grid)), scale = scale)
}
profOpts <- function(scale, dist, curves, correct = TRUE, tau = TRUE) {
  m1 <- cacoa:::shiftPoissonTable()
  list(scale = scale, dist = dist, noise = TRUE, tau = tau, correct = correct, mu_grid = curves$mu.grid, phi = curves$phi, tau2 = curves$tau2,
       lam_grid = m1$lam.grid, s_grid = m1$s.grid, V1 = m1$V1)
}

test_that("Poisson lookup table: exact variance of the transform, delta method at large lambda", {
  m1 <- cacoa:::shiftPoissonTable()
  i <- which.min(abs(m1$lam.grid - 0)); j <- which.min(abs(m1$s.grid - 0))   # lambda = 1, s = 1
  y <- 0:60; p <- dpois(y, 1); t <- log10(y + 1)
  expect_equal(m1$V1[i, j], sum(p * t^2) - sum(p * t)^2, tolerance = 1e-10)
  i <- length(m1$lam.grid); j <- which.min(abs(m1$s.grid + 3)); L <- 1e3; s <- 1e-3
  expect_equal(m1$V1[i, j], (s / log(10))^2 * L / (s * L + 1)^2, tolerance = 0.01)
})

test_that("unit profile: transform, M2 noise, weights and both distances match their formulas", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 30, n.genes = 60, shift = 0.8)
  cm <- Matrix::t(cao$data.object); spc <- factor(as.character(cao$sample.per.cell[colnames(cm)])); n.s <- nlevels(spc)
  lib <- Matrix::colSums(cm) * 3   # library sizes over "all genes" (more than the profile genes)
  ids <- which(cao$cell.groups[colnames(cm)] == "ct1") - 1L
  pairs <- t(utils::combn(n.s, 2)) - 1L; storage.mode(pairs) <- "integer"
  scale <- median(lib); cv <- toyCurves(scale)
  for (dist in 1:2) for (correct in c(FALSE, TRUE)) {
    r <- cacoa:::cf_unit_profiles(cm, as.integer(spc), list(as.integer(ids)), as.numeric(lib), n.s, 3L, pairs, profOpts(scale, dist, cv, correct))[[1]]
    sa <- as.integer(spc)[ids + 1L]
    Y <- sapply(seq_len(n.s), function(a) Matrix::rowSums(cm[, ids[sa == a] + 1L, drop = FALSE])); N <- sapply(seq_len(n.s), function(a) sum(lib[ids[sa == a] + 1L]))
    keep <- rowSums(Y) > 0; Y <- Y[keep, ]
    x <- log10(scale * sweep(Y, 2, N, "/") + 1)
    expect_equal(r$x, unname(x), tolerance = 1e-12)
    # noise: Poisson part (exact for small lambda) + heterogeneity part; weights 1 / (mean noise + tau2)
    pbar <- rowSums(Y) / sum(N); tp <- scale / ((scale * pbar + 1) * log(10))
    S2 <- sapply(seq_len(n.s), function(a) sum(lib[ids[sa == a] + 1L]^2)); neff <- N^2 / S2
    het <- outer(tp^2 * 0.2 * pbar^2, 1 / neff)
    pois <- r$v - het
    expect_true(all(pois > 0))
    k <- which.min(abs(pbar * N[1] - 1)); lam <- pbar[k] * N[1]; yy <- 0:200; pp <- dpois(yy, lam); tt <- log10(scale * yy / N[1] + 1)
    expect_equal(unname(pois[k, 1]), sum(pp * tt^2) - sum(pp * tt)^2, tolerance = 0.02)   # bilinear interpolation of the table
    w <- 1 / (rowMeans(r$v) + 1e-3)
    expect_equal(as.numeric(r$w), w, tolerance = 1e-12)
    cc <- colSums(w * r$v); expect_equal(as.numeric(r$c), cc, tolerance = 1e-10)
    Z <- sqrt(w) * (if (dist == 1) x - rowMeans(x) else x)
    D <- if (dist == 2) { d <- as.matrix(dist(t(Z)))^2; if (correct) d <- d - outer(cc, cc, "+"); d } else {
      n2 <- colSums(Z^2); if (correct) n2 <- pmax(n2 - cc, 0.05 * n2); R <- crossprod(Z) / sqrt(outer(n2, n2)); R[] <- pmin(1, pmax(-1, R)); 1 - R }   # clamped to [-1, 1]
    expect_equal(as.numeric(r$y), unname(D[pairs + 1L]), tolerance = 1e-10)
  }
})

test_that("normalized shift map: runs for both distances, thread-independent, scale and curve checks", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.2)
  cm <- Matrix::t(cao$data.object); cg <- cao$cell.groups; idx <- split(seq_along(cg) - 1L, cg)
  set.seed(3); nns <- lapply(seq_along(cg), function(i) sample(seq_along(cg) - 1L, 120)); names(nns) <- names(cg)
  lib <- Matrix::colSums(cm)
  cv <- cacoa:::fitShiftNoiseCurves(cm, lib, cao$sample.per.cell[colnames(cm)], nns, scale = median(lib), n.units = 50)
  expect_true(all(is.finite(cv$phi)) && all(cv$phi >= 0) && all(is.finite(cv$tau2)) && all(cv$tau2 >= 0))
  for (dist in c("cor", "euclidean")) {
    r1 <- clusterFreeExpressionShifts(statistic = "distance", cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, dist = dist, lib = lib, noise.curves = cv,
                                      n.permutations = 29, seed = 1, smooth = FALSE, n.cores = 1)
    r3 <- clusterFreeExpressionShifts(statistic = "distance", cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, dist = dist, lib = lib, noise.curves = cv,
                                      n.permutations = 29, seed = 1, smooth = FALSE, n.cores = 3)
    expect_equal(r1$stat, r3$stat, tolerance = 1e-12); expect_equal(r1$p.value, r3$p.value)
    expect_true(sum(is.finite(r1$stat)) > 0.9 * length(nns))
    expect_equal(r1$settings$profile, "normalized"); expect_equal(r1$settings$scale, median(lib))
    expect_false(r1$settings$noise.correction)
    # the shift is detected in this toy (a quarter of the genes shifted in group B)
    expect_lt(median(r1$p.value, na.rm = TRUE), 0.2)
  }
  # curves fitted at another scale are refused; "median" equals the explicit median
  expect_error(clusterFreeExpressionShifts(statistic = "distance", cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, lib = lib, profile.scale = 100, noise.curves = cv, n.permutations = 9), "different scale")
  ra <- clusterFreeExpressionShifts(statistic = "distance", cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, lib = lib, noise.curves = cv, n.permutations = 9, seed = 1, smooth = FALSE)
  rb <- clusterFreeExpressionShifts(statistic = "distance", cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, lib = lib, profile.scale = median(lib), noise.curves = cv, n.permutations = 9, seed = 1, smooth = FALSE)
  expect_equal(ra$stat, rb$stat)
  # unweighted normalized profiles, and the guard on distance / profile combinations
  rn <- clusterFreeExpressionShifts(statistic = "distance", cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, lib = lib, weighting = "none", n.permutations = 9, seed = 1, smooth = FALSE)
  expect_true(sum(is.finite(rn$stat)) > 0.9 * length(nns))
  expect_error(clusterFreeExpressionShifts(statistic = "distance", cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, dist = "js"), "normalized")
  expect_error(clusterFreeExpressionShifts(statistic = "distance", cm, cao$sample.per.cell, nns, cao$model, cao$sample.meta, dist = "euclidean", profile = "raw"), "normalized")
})
