## Sampling-noise model of the normalized neighbourhood profiles of the cluster-free shift map (profile = "normalized").
##
## Per neighbourhood unit and sample a, the profile is x_ga = log10(scale p_ga + 1), p_ga = Y_ga / N_a (summed counts over
## the summed library sizes of the sample's cells in the unit). The noise of x_ga is modelled without the labels, with
## pbar_g pooled over the unit's samples: the exact variance of the transform of a Poisson count (shiftPoissonTable()) plus
## the within-unit heterogeneity t'(pbar)^2 phi(mu) pbar^2 / n_eff,a (mu = log10 pbar; n_eff = N^2 / sum of squared library
## sizes). The gene weights are 1 / (mean noise + tau2(mu)), with tau2 the between-sample variance floor. phi(mu) and
## tau2(mu) are smooth curves fitted on a random subset of the analysed units (fitShiftNoiseCurves()).

shift.noise.env <- new.env(parent = emptyenv())

#' Variance of log10(s Y + 1) for a Poisson count Y
#'
#' Lookup table used by the noise model of the normalized cluster-free shift profiles: the exact variance over the
#' Poisson distribution on a grid of log10(lambda) in \[-3, 3\] and log10(s) in \[-5, 4\] (s = scale / N, with N the
#' sample's summed library size in the unit). Above the lambda range the kernel uses the delta method. Computed once
#' per session.
#' @return list: `lam.grid`, `s.grid` (log10 axes) and `V1` (lambda x s matrix of variances)
#' @keywords internal
shiftPoissonTable <- function() {
  if (!is.null(shift.noise.env$m1)) return(shift.noise.env$m1)
  lam <- seq(-3, 3, by = 0.05); sg <- seq(-5, 4, by = 0.1)
  V <- matrix(NA_real_, length(lam), length(sg))
  for (i in seq_along(lam)) {
    L <- 10^lam[i]; y <- 0:stats::qpois(1 - 1e-13, L); p <- stats::dpois(y, L); p <- p / sum(p)
    tt <- log10(outer(y, 10^sg) + 1)                       # counts x s
    m <- colSums(p * tt); V[i, ] <- colSums(p * tt^2) - m^2
  }
  shift.noise.env$m1 <- list(lam.grid = lam, s.grid = sg, V1 = pmax(V, 0))
  shift.noise.env$m1
}

# per bin of mu (width 0.1) the mean of the values capped at the bin's 99.9 % quantile (the moment estimates of low genes
# are zero-inflated and heavy-tailed), loess over the bin centres (weights sqrt(bin size)), evaluated on mu.grid and held
# flat outside the fitted range
binnedNoiseCurve <- function(mu, val, mu.grid, cap = 0.999, min.bin = 20) {
  ok <- is.finite(mu) & is.finite(val); key <- floor(mu[ok] / 0.1); v <- val[ok]
  # values grouped by bin through integer codes and a stable order (split() on the numeric bin would convert every value to a
  # string); bins, their order and the centres (via the same character round-trip as split's names) are those of split(v, b)
  uk <- sort(unique(key)); code <- match(key, uk); o <- order(code, method = "radix"); vs <- v[o]
  cnt <- tabulate(code, length(uk)); ends <- cumsum(cnt); starts <- ends - cnt + 1L
  m <- vapply(seq_along(uk), function(i) { z <- vs[starts[i]:ends[i]]; mean(pmin(z, stats::quantile(z, cap, names = FALSE))) }, 0)
  d <- data.frame(b = as.numeric(as.character(uk * 0.1 + 0.05)), m = m, n = cnt)
  d <- d[d$n >= min.bin, , drop = FALSE]
  if (nrow(d) < 4) stop("too few expressed genes to fit the noise curves")
  lo <- withCallingHandlers(stats::loess(m ~ b, data = d, weights = sqrt(d$n), span = 0.5),   # few bins: loess's residual scale is NaN (harmless), its fit is not
    warning = function(w) if (grepl("NaNs produced", conditionMessage(w))) invokeRestart("muffleWarning"))
  pr <- stats::predict(lo, data.frame(b = pmin(pmax(mu.grid, min(d$b)), max(d$b))))
  list(curve = pmax(as.numeric(pr), 0), bins = d)
}

#' Fit the noise curves of the normalized cluster-free shift profiles
#'
#' Label-free: uses the counts, library sizes and sample of each cell in a random subset of the neighbourhood units,
#' never the sample groups. Per unit, sample a (at least `min.n.obs.per.samp` cells) and gene g:
#' - phi (within-unit heterogeneity, the squared CV of the per-cell rates): the moment estimate
#'   \[sum_i (y_i - p m_i)^2 - Y\] / (p^2 sum_i m_i^2), p = Y / N;
#' - tau2 (between-sample floor): the variance of x_ga across the unit's samples minus their mean predicted sampling
#'   variance (Poisson + phi).
#' Each is binned on mu = log10 pbar (pbar pooled over the unit's samples) and smoothed (see the source for details).
#' @param cm genes x cells count matrix (the genes of the shift profiles)
#' @param lib library size per cell over all genes (aligned to the columns of `cm`)
#' @param sample.per.cell factor of samples (aligned to the columns of `cm`)
#' @param nns.per.cell list of 0-based cell indices per unit
#' @param scale transform scale (x = log10(scale p + 1))
#' @param n.units units used for the fit (default 200)
#' @param min.n.obs.per.samp minimum cells of a sample in a unit (default 3)
#' @param seed seed of the unit draw
#' @return list: `mu.grid`, `phi`, `tau2`, `scale`, `n.units`, `bins`
#' @keywords internal
fitShiftNoiseCurves <- function(cm, lib, sample.per.cell, nns.per.cell, scale, n.units = 200, min.n.obs.per.samp = 3, seed = 1) {
  spc <- as.factor(sample.per.cell); spcI <- as.integer(spc); n.s <- nlevels(spc)
  mu.grid <- seq(-7, -1, by = 0.05)
  sizes <- lengths(nns.per.cell)
  cand <- which(sizes >= 3 * min.n.obs.per.samp)
  if (!length(cand)) stop("no neighborhood unit is large enough to fit the noise curves")
  units <- withSeed(seed, cand[sample.int(length(cand), min(n.units, length(cand)))])
  # only the cells of the drawn units: the C++ calls below convert the matrix they get
  cells <- sort(unique(unlist(lapply(nns.per.cell[units], as.integer)))) + 1L
  remap <- integer(ncol(cm)); remap[cells] <- seq_along(cells) - 1L
  cm <- methods::as(cm[, cells, drop = FALSE], "CsparseMatrix"); lib <- as.numeric(lib)[cells]; spcI <- spcI[cells]
  nns.per.cell <- lapply(nns.per.cell[units], function(v) remap[as.integer(v) + 1L]); units <- seq_along(units)
  ph <- cf_phi_moments(cm, spcI, lapply(nns.per.cell[units], as.integer), as.numeric(lib), n.s, as.integer(min.n.obs.per.samp))   # one C++ pass per unit
  if (ph$n_units == 0) stop("no unit has 3 samples with enough cells to fit the noise curves")
  fphi <- binnedNoiseCurve(ph$mu, ph$v, mu.grid)
  m1 <- shiftPoissonTable()
  opts <- list(scale = scale, dist = 2L, noise = TRUE, tau = FALSE, correct = FALSE, mu_grid = mu.grid, phi = fphi$curve, tau2 = rep(0, length(mu.grid)),
               lam_grid = m1$lam.grid, s_grid = m1$s.grid, V1 = m1$V1)
  pairs <- matrix(0L, 0, 2)
  td <- list()
  for (r in cf_unit_profiles(cm, spcI, lapply(nns.per.cell[units], as.integer), lib, n.s, as.integer(min.n.obs.per.samp), pairs, opts)) {
    pres <- as.integer(r$present) + 1L; if (length(pres) < 3) next
    x <- r$x[, pres, drop = FALSE]; v <- r$v[, pres, drop = FALSE]
    td[[length(td) + 1]] <- cbind(as.numeric(r$mu), matrixStats::rowVars(x) - rowMeans(v))
  }
  td <- do.call(rbind, td); ftau <- binnedNoiseCurve(td[, 1], td[, 2], mu.grid)
  list(mu.grid = mu.grid, phi = fphi$curve, tau2 = ftau$curve, scale = scale, n.units = length(units), bins = list(phi = fphi$bins, tau2 = ftau$bins))
}
