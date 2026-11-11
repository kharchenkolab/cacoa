## Random intercept for one grouping variable in a per-feature linear model, with inference by relabeling.
##
## Model per feature: y = F beta + u[group] + e, u ~ N(0, theta sigma^2), e ~ N(0, sigma^2). theta (the between-group variance relative to
## the within-group variance) is estimated by restricted maximum likelihood, profiled over a grid. Multiplying the data and the design by
## V^(-1/2), V = I + theta J within a group, gives uncorrelated errors of equal variance, so that least squares on the transformed model
## is the generalized least-squares fit and its residuals can be exchanged between samples (Freedman-Lane). The estimate of theta is
## uncertain in small cohorts; treating it as known makes the test anticonservative (6-8 % at a nominal 5 % with 17-24 samples), so the
## whole procedure, including the estimation of theta, is repeated on every relabeled data set (4-6 % in the same simulations, for
## paired, partly confounded and between-group designs).

# W = V^(-1/2) as a function of a matrix (rows = samples): W = I - k_b J_b / m_b per group b, k_b = 1 - 1 / sqrt(1 + m_b theta)
groupWhitener <- function(g, theta) {
  g <- as.integer(factor(g)); m <- tabulate(g)[g]; k <- 1 - 1 / sqrt(1 + m * theta)
  function(M) { M <- as.matrix(M); R <- rowsum(M, g)[g, , drop = FALSE]; dimnames(R) <- NULL; M - (R / m) * k }
}
# V^(1/2) = I + k_b J_b / m_b, k_b = sqrt(1 + m_b theta) - 1
groupUnwhitener <- function(g, theta) {
  g <- as.integer(factor(g)); m <- tabulate(g)[g]; k <- sqrt(1 + m * theta) - 1
  function(M) { M <- as.matrix(M); R <- rowsum(M, g)[g, , drop = FALSE]; dimnames(R) <- NULL; M + (R / m) * k }
}

# restricted maximum likelihood estimate of theta for the columns of Y (one common value), profiled over `grid`
remlVarianceRatio <- function(Y, Fm, g, grid = c(0, 10^seq(-2, 2, length.out = 25))) {
  Y <- as.matrix(Y); n <- nrow(Fm); mb <- tabulate(as.integer(factor(g)))
  crit <- vapply(grid, function(th) {
    W <- groupWhitener(g, th); q <- qr(W(Fm)); p <- q$rank; rss <- colSums(qr.resid(q, W(Y))^2)
    if (n - p < 1 || any(rss <= 0)) return(Inf)
    sum((n - p) * log(rss / (n - p))) + ncol(Y) * (sum(log(1 + mb * th)) + 2 * sum(log(abs(diag(qr.R(q))[seq_len(p)]))))
  }, numeric(1))
  grid[which.min(crit)]
}

#' Contrast test with a random intercept per group, by relabeling
#'
#' For every column of `Y` (a feature observed on the samples): the generalized least-squares estimate of the contrast under a random
#' intercept per level of `g` (variance ratio by restricted maximum likelihood, per feature), its standard error and t statistic, and
#' the same quantities on `ncol(P)` relabeled data sets. A relabeled data set is the fit of the model without the tested contrast plus
#' the permuted residuals of that fit, both taken on the transformed (uncorrelated) scale and transformed back; the variance ratio is
#' re-estimated on each.
#'
#' @param Y samples x features matrix
#' @param Fm design matrix (samples x coefficients)
#' @param cvec contrast over the columns of `Fm`
#' @param g grouping variable (one value per sample)
#' @param P samples x B matrix of permutation indices (free permutations of the samples)
#' @return list: `est`, `se`, `t`, `theta` (per feature), `est.perm`, `t.perm` (features x B), `df`
#' @keywords internal
randomInterceptContrast <- function(Y, Fm, cvec, g, P) {
  Y <- as.matrix(Y); n <- nrow(Fm); N0 <- Fm %*% MASS::Null(matrix(cvec, ncol = 1)); B <- ncol(P)
  stat <- function(v) {
    th <- remlVarianceRatio(v, Fm, g); W <- groupWhitener(g, th); vw <- drop(W(v)); Fw <- W(Fm); q <- qr(Fw); df <- n - q$rank
    XtXi <- chol2inv(qr.R(q)); a <- drop(Fw %*% XtXi %*% cvec); est <- sum(a * vw); se <- sqrt(sum(cvec * (XtXi %*% cvec)) * sum(qr.resid(q, vw)^2) / df)
    c(est, se, th, df)
  }
  out <- list(est = numeric(ncol(Y)), se = numeric(ncol(Y)), theta = numeric(ncol(Y)), est.perm = matrix(NA_real_, ncol(Y), B), t.perm = matrix(NA_real_, ncol(Y), B))
  for (j in seq_len(ncol(Y))) {
    y <- Y[, j]; s0 <- stat(y); out$est[j] <- s0[1]; out$se[j] <- s0[2]; out$theta[j] <- s0[3]; out$df <- s0[4]
    th0 <- if (ncol(N0)) remlVarianceRatio(y, N0, g) else 0; W0 <- groupWhitener(g, th0); U0 <- groupUnwhitener(g, th0); yw <- drop(W0(y))
    fit0 <- if (ncol(N0)) qr.fitted(qr(W0(N0)), yw) else rep(0, n); r0 <- yw - fit0
    for (b in seq_len(B)) { sp <- stat(drop(U0(fit0 + r0[P[, b]]))); out$est.perm[j, b] <- sp[1]; out$t.perm[j, b] <- sp[1] / sp[2] }
  }
  out$t <- out$est / out$se; names(out$est) <- names(out$se) <- names(out$t) <- names(out$theta) <- colnames(Y)
  rownames(out$est.perm) <- rownames(out$t.perm) <- colnames(Y)
  out
}

#' Studentized permutation test of a contrast for every column of a response matrix
#'
#' Relabeling in R for a small number of features. `scheme = "freedman-lane"`: the fit of the model without the tested contrast plus its
#' exchanged residuals gives a relabeled data set, on which the contrast and its t statistic are recomputed. `scheme = "block"`: the rows
#' of the design are exchanged among the samples (within the cells of the model's relabeling plan) and the model is refitted. The t statistic (the estimate over
#' its standard error in the same relabeled data set) is used for the p-value: with covariates in the model, the estimate alone is not
#' comparable between the observed and the relabeled data, because the residuals of the smaller model are spread over fewer degrees
#' of freedom (false-positive rate 21 % at a nominal 5 % for the composition model with a patient effect, 17 samples, 11 parameters;
#' 4-5 % with the t statistic).
#'
#' @param Y samples x features matrix
#' @param Fm design matrix
#' @param cvec contrast over the columns of `Fm`
#' @param P samples x B matrix of permutation indices (as drawn by the model's relabeling plan)
#' @param scheme `"freedman-lane"` (default) or `"block"`
#' @return list: `est`, `se`, `t` (per feature), `est.perm`, `t.perm` (features x B), `df`
#' @keywords internal
contrastPermutationTest <- function(Y, Fm, cvec, P, scheme = c("freedman-lane", "block")) {
  scheme <- match.arg(scheme)
  Y <- as.matrix(Y); n <- nrow(Fm); hi <- hatInfo(Fm); a <- drop(t(hi$A) %*% cvec); cc <- sum(a^2); q <- qr(Fm); df <- n - hi$rank
  if (df < 1) stop("no residual degrees of freedom")
  N0 <- Fm %*% MASS::Null(matrix(cvec, ncol = 1)); q0 <- if (ncol(N0)) qr(N0) else NULL
  F0 <- if (is.null(q0)) Y * 0 else qr.fitted(q0, Y); R0 <- Y - F0
  # the residuals of the smaller model are smaller than the errors: they are multiplied by sqrt(kappa), kappa = the variance of the
  # observed over the mean variance of the relabeled estimate, |a|^2 / mean_b(|a|^2 - |Q0' a_b|^2) (a_b: the contrast's weights moved
  # with the residuals), so that the relabeled estimates are on the scale of the observed one (the t statistic does not depend on it)
  kappa <- 1
  if (!is.null(q0) && q0$rank > 0) { Q0 <- qr.Q(q0)[, seq_len(q0$rank), drop = FALSE]
    g <- mean(apply(P, 2, function(p) { ab <- numeric(n); ab[p] <- a; sum(crossprod(Q0, ab)^2) })); if (cc - g > 1e-8 * cc) kappa <- cc / (cc - g) }
  R0 <- R0 * sqrt(kappa)
  one <- function(Z) { est <- drop(crossprod(a, Z)); se <- sqrt(cc * colSums(qr.resid(q, Z)^2) / df); rbind(est, se) }
  o <- one(Y); B <- ncol(P); ep <- tp <- matrix(NA_real_, ncol(Y), B, dimnames = list(colnames(Y), NULL))
  # block: the design F[p, ] fitted to Y is the design F fitted to the responses moved the other way
  for (b in seq_len(B)) { z <- one(if (scheme == "block") Y[order(P[, b]), , drop = FALSE] else F0 + R0[P[, b], , drop = FALSE]); ep[, b] <- z[1, ]; tp[, b] <- z[1, ] / z[2, ] }
  tp[is.nan(tp)] <- 0     # no estimate and no residual variance in a relabeled data set (a feature constant over the samples)
  list(est = setNames(o[1, ], colnames(Y)), se = setNames(o[2, ], colnames(Y)), t = setNames(o[1, ] / o[2, ], colnames(Y)), est.perm = ep, t.perm = tp, df = df, kappa = kappa)
}
