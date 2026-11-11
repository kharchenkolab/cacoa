# Per-gene combination statistic of the cluster-free shift map (statistic = "gsd"): the kernel's observed and permuted unit
# statistics against an R reimplementation from the unit profiles, for two groups (subset-table kernel) and for general
# designs (continuous variable, covariates, block relabeling and Freedman-Lane), with and without winsorization; the
# tested centres (centres.per.cell) against a direct average of the per-neighbourhood statistics; argument checks.

gsdToy <- function(seed = 42, shift = 1) {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = shift, seed = seed)
  cm <- Matrix::t(cao$data.object); cg <- cao$cell.groups; lib <- Matrix::colSums(cm) * 2
  set.seed(seed); units <- lapply(seq_along(cg), function(i) sort(sample(seq_along(cg) - 1L, 120))); names(units) <- names(cg)
  cv <- list(mu.grid = seq(-7, -1, by = 0.05), phi = rep(0.2, 121), tau2 = rep(1e-3, 121), scale = median(lib))
  list(cao = cao, cm = cm, lib = lib, units = units, cv = cv)
}
# reference statistic of one unit from its profile r (cf_unit_profiles): rows of P are the samples' relabelings (1-based)
refGsd <- function(r, X, cvec, P, fl, strata, mode = "moderated", gamma = 0, d0 = 3, winsor.k = 0, lev = NULL, minl = 2) {
  pres <- r$present + 1L; m <- length(pres); Y <- r$x[, pres, drop = FALSE]; V <- r$v[, pres, drop = FALSE] + as.numeric(r$tau2)
  if (winsor.k > 0) Y <- t(sapply(seq_len(nrow(Y)), function(g) { x <- Y[g, ]; md <- median(x); s2 <- (1.4826 * median(abs(x - md)))^2; sd <- sqrt(pmax(V[g, ], s2))
    pmin(pmax(x, md - winsor.k * sd), md + winsor.k * sd) }))
  Y <- Y - rowMeans(Y)
  stat1 <- function(Xp, Yp, cp, Vp = V, mv.scale = 1) {
    XtXi <- MASS::ginv(crossprod(Xp)); a <- Xp %*% XtXi %*% cp; beta <- as.numeric(Yp %*% a); mv <- mv.scale * as.numeric(Vp %*% a^2)
    H <- Xp %*% XtXi %*% t(Xp); rss <- rowSums((Yp - Yp %*% H)^2); cxc <- as.numeric(t(cp) %*% XtXi %*% cp); df <- m - qr(Xp)$rank
    va <- if (mode == "moderated") (d0 * mv + cxc * rss) / (d0 + df) else mv; z2 <- beta^2 / va; c(sum(sqrt(z2)), sum(z2), sum(z2^2), sum(z2^4), max(z2)) }
  Xs <- X[pres, , drop = FALSE]
  if (!fl) Tm <- cbind(stat1(Xs, Y, cvec), sapply(seq_len(ncol(P)), function(b) stat1(X[P[pres, b], , drop = FALSE], Y, cvec)))
  else {
    Qc <- qr.Q(qr(cvec), complete = TRUE); Xr <- Xs %*% Qc[, -1, drop = FALSE]; H0 <- Xr %*% MASS::ginv(crossprod(Xr)) %*% t(Xr); E <- Y - Y %*% H0; F0 <- Y %*% H0
    blocks <- split(seq_len(m), strata[pres])
    induced <- function(p) { q <- seq_len(m); for (b in blocks) { if (length(b) < 2) next; blk <- sort(b); q[blk[order(p[pres[blk]], method = "radix")]] <- blk }; q }
    # the residuals of the model without the tested contrast are smaller than the errors: the model variance of the relabeled data sets
    # is multiplied by 1 / kappa = (|a|^2 - mean_b |Q0' a_b|^2) / |a|^2 (a_b: the contrast's weights moved with the residuals)
    a0 <- as.numeric(Xs %*% MASS::ginv(crossprod(Xs)) %*% cvec); aa <- sum(a0^2); sv <- svd(Xr); Q0 <- sv$u[, sv$d > 1e-7 * max(sv$d), drop = FALSE]
    qs <- lapply(seq_len(ncol(P)), function(b) induced(P[, b])); gs <- mean(vapply(qs, function(q) { ab <- numeric(m); ab[q] <- a0; sum(crossprod(Q0, ab)^2) }, numeric(1)))
    mvs <- if (aa > 0 && aa - gs > 1e-8 * aa) (aa - gs) / aa else 1
    Tm <- cbind(stat1(Xs, Y, cvec), sapply(seq_len(ncol(P)), function(b) { q <- qs[[b]]; stat1(Xs, F0 + E[, q, drop = FALSE], cvec, V[, q, drop = FALSE], mvs) }))
  }
  adm <- rep(TRUE, ncol(Tm))
  if (!is.null(lev)) { idx <- cbind(pres, P[pres, , drop = FALSE]); adm <- apply(idx, 2, function(ii) sum(lev[ii] == 1) >= minl && sum(lev[ii] == 2) >= minl); Tm[, !adm] <- NA }
  if (!adm[1]) return(rep(NA_real_, ncol(Tm)))
  if (gamma > 0) return(Tm[match(gamma, c(1, 2, 4, 8, 99)), ])
  f <- rep(-Inf, ncol(Tm)); okp <- c(FALSE, adm[-1])
  for (k in 1:5) { lt <- log(pmax(Tm[k, ], 1e-300)); mu <- mean(lt[okp]); s <- sd(lt[okp]); if (s > 0) f <- pmax(f, (lt - mu) / s) }
  f[!adm] <- NA; f
}
checkAgainstReference <- function(ty, model, scheme = "auto", n.units = 12, B = 29, ...) {
  a <- list(...); cao <- ty$cao
  run <- function() clusterFreeExpressionShifts(ty$cm, cao$sample.per.cell, ty$units[seq_len(n.units)], model, model$meta, lib = ty$lib, noise.curves = ty$cv,
                                                n.permutations = B, seed = 7, permutation = scheme, adjust = FALSE, smooth = FALSE, post = "none", return.perm = TRUE,
                                                gsd.variance = a$variance %||% "moderated", gsd.gamma = a$gamma %||% 0, robust = if ((a$winsor.k %||% 0) > 0) "winsor" else "none", winsor.k = a$winsor.k %||% 3)
  k <- run()
  samples <- levels(droplevels(factor(cao$sample.per.cell))); n <- length(samples); des <- model; des$F <- model$F[samples, , drop = FALSE]
  gplan <- permutationPlan(des, model$meta[samples, , drop = FALSE], samples, scheme = scheme, block.vars = model$block.vars, n.permutations = B, max.enumerate = 0)
  P <- withSeed(7, drawPermutations(gplan, B)); fl <- gplan$scheme == "freedman-lane"; cF <- model$contrast.F[colnames(des$F)]
  spec <- model$contrast_spec; lev <- NULL
  if (!is.null(spec) && spec$type %in% c("simple", "marginal") && spec$term %in% names(model$meta) && !is.numeric(model$meta[[spec$term]])) {
    g <- as.character(model$meta[samples, spec$term]); lev <- ifelse(g == spec$den, 1, ifelse(g == spec$num, 2, 0)) }
  m1 <- shiftPoissonTable(); pairs <- t(utils::combn(n, 2)) - 1L; storage.mode(pairs) <- "integer"
  opts <- list(scale = ty$cv$scale, dist = 2L, noise = TRUE, tau = TRUE, correct = FALSE, unit_detection = 0.05, mu_grid = ty$cv$mu.grid, phi = ty$cv$phi, tau2 = ty$cv$tau2,
               lam_grid = m1$lam.grid, s_grid = m1$s.grid, V1 = m1$V1)
  spc <- factor(as.character(cao$sample.per.cell[colnames(ty$cm)]), levels = samples)
  prof <- cf_unit_profiles(ty$cm, as.integer(spc), lapply(ty$units[seq_len(n.units)], as.integer), as.numeric(ty$lib), n, 3L, pairs, opts)
  err <- sapply(seq_len(n.units), function(i) {
    f <- refGsd(prof[[i]], des$F, cF, P, fl, as.integer(gplan$strata), mode = a$variance %||% "moderated", gamma = a$gamma %||% 0, winsor.k = a$winsor.k %||% 0, lev = lev)
    if (!is.finite(f[1])) return(c(if (is.finite(k$stat[i])) Inf else 0, 0))
    c(abs(f[1] - k$stat[i]) / max(1, abs(f[1])), max(abs(f[-1] - k$stats.perm[, i]) / pmax(1, abs(f[-1])), na.rm = TRUE)) })
  list(err = max(err), tested = sum(is.finite(k$stat)), scheme = gplan$scheme)
}

test_that("gsd, two groups: observed and permuted unit statistics equal the R reference", {
  ty <- gsdToy()
  for (variance in c("moderated", "model")) for (gamma in c(0, 2, 99)) {
    r <- checkAgainstReference(ty, ty$cao$model, variance = variance, gamma = gamma)
    expect_lt(r$err, 1e-10); expect_gt(r$tested, 5)
  }
  expect_lt(checkAgainstReference(ty, ty$cao$model, winsor.k = 1.5)$err, 1e-10)
})

test_that("gsd, general designs: continuous variable, covariates, Freedman-Lane", {
  ty <- gsdToy(); meta <- ty$cao$sample.meta
  m.cont <- buildCacoaModel(meta, formula = ~ age, test = "age")
  m.adj <- buildCacoaModel(meta, formula = ~ group + age, test = "age")
  r1 <- checkAgainstReference(ty, m.cont); expect_lt(r1$err, 1e-10)
  r2 <- checkAgainstReference(ty, m.adj); expect_lt(r2$err, 1e-10)
  r3 <- checkAgainstReference(ty, m.adj, scheme = "freedman-lane"); expect_equal(r3$scheme, "freedman-lane"); expect_lt(r3$err, 1e-10)
  expect_lt(checkAgainstReference(ty, m.cont, gamma = 99, variance = "model")$err, 1e-10)
  expect_lt(checkAgainstReference(ty, m.adj, scheme = "freedman-lane", winsor.k = 2)$err, 1e-10)
})

test_that("tested centres: each cell's statistic is the mean over the tested neighbourhoods containing it", {
  ty <- gsdToy(); cao <- ty$cao; m <- length(ty$units)
  set.seed(5); graph <- lapply(seq_len(m), function(i) sort(unique(c(i - 1L, sample(m, 8) - 1L)))); names(graph) <- names(ty$units)   # a sparse "graph"
  args <- list(ty$cm, cao$sample.per.cell, ty$units, cao$model, cao$sample.meta, lib = ty$lib, noise.curves = ty$cv, n.permutations = 19, seed = 3, smooth = FALSE, nns.graph = graph)
  full <- do.call(clusterFreeExpressionShifts, c(args, list(post = "none", adjust = FALSE, return.perm = TRUE)))   # every unit's own statistic
  thin <- do.call(clusterFreeExpressionShifts, c(args, list(centres.per.cell = 2, return.perm = TRUE)))
  cover <- integer(m); idx <- logical(m); ord <- withSeed(3, sample.int(m))
  for (i in ord) if (cover[i] < 2) { idx[i] <- TRUE; j <- graph[[i]] + 1L; cover[j] <- cover[j] + 1L }
  expect_equal(thin$settings$n.neighbourhoods.tested, sum(idx)); expect_lt(sum(idx), m)
  memb <- lapply(which(idx), function(i) unique(c(i, ty$units[[i]] + 1L)))
  contrib <- split(rep(which(idx), lengths(memb)), unlist(memb))
  expect_equal(sum(is.finite(thin$stat)), length(contrib))   # every covered cell is tested
  cells <- as.integer(names(contrib))
  ref.obs <- vapply(contrib, function(u) mean(full$stat[u][is.finite(full$stat[u])]), 0)
  expect_equal(unname(thin$stat[cells]), unname(ref.obs), tolerance = 1e-5)   # float32 averaging buffer
  ref.perm <- sapply(contrib, function(u) rowMeans(full$stats.perm[, u, drop = FALSE]))
  expect_equal(unname(thin$stats.perm[, cells]), unname(ref.perm), tolerance = 1e-5)
})

test_that("gsd defaults: thread-independent, calls by BH on the pooled p, a shift is found and a null is quiet", {
  ty <- gsdToy(shift = 1.2); cao <- ty$cao
  r1 <- clusterFreeExpressionShifts(ty$cm, cao$sample.per.cell, ty$units, cao$model, cao$sample.meta, lib = ty$lib, noise.curves = ty$cv, n.permutations = 29, seed = 1, smooth = FALSE, n.cores = 1)
  r2 <- clusterFreeExpressionShifts(ty$cm, cao$sample.per.cell, ty$units, cao$model, cao$sample.meta, lib = ty$lib, noise.curves = ty$cv, n.permutations = 29, seed = 1, smooth = FALSE, n.cores = 3)
  expect_equal(r1$stat, r2$stat, tolerance = 1e-6); expect_equal(r1$p.pooled, r2$p.pooled, tolerance = 1e-6)
  expect_equal(r1$settings$statistic, "gsd"); expect_equal(r1$settings$adjust.method, "bh"); expect_equal(r1$settings$p.source, "pooled"); expect_equal(r1$post, "averaged")
  expect_equal(unname(r1$p.adj), unname(stats::p.adjust(r1$p.pooled, "BH")))
  expect_true(all(is.na(r1$p.param)))
  expect_gt(mean(r1$p.adj < 0.1, na.rm = TRUE), 0.5)
  t0 <- gsdToy(shift = 0); c0 <- t0$cao
  r0 <- clusterFreeExpressionShifts(t0$cm, c0$sample.per.cell, t0$units, c0$model, c0$sample.meta, lib = t0$lib, noise.curves = t0$cv, n.permutations = 29, seed = 1, smooth = FALSE)
  expect_lt(mean(r0$p.adj < 0.1, na.rm = TRUE), 0.1)
})

test_that("gsd argument checks", {
  ty <- gsdToy(); cao <- ty$cao
  f <- function(...) clusterFreeExpressionShifts(ty$cm, cao$sample.per.cell, ty$units, cao$model, cao$sample.meta, lib = ty$lib, noise.curves = ty$cv, n.permutations = 9, ...)
  expect_error(f(profile = "raw"), "normalized")
  expect_error(f(weighting = "none"), "weighting")
  expect_error(f(robust = "huber"), "huber")
  expect_error(f(centres.per.cell = 1.5), "centres.per.cell")
  expect_error(f(centres.per.cell = 0), "centres.per.cell")
  expect_error(f(neighbourhood.gene.frac = 1), "neighbourhood.gene.frac")
  expect_error(f(gsd.gamma = 3), "gamma")
  expect_silent(f(centres.per.cell = Inf)); expect_silent(f(centres.per.cell = NULL))
})
