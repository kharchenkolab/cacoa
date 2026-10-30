# Inverse of the trigamma function by Newton's method (for the moments fit of the prior degrees of freedom)
trigammaInverse <- function(y) {
  y <- as.numeric(y); x <- rep(NA_real_, length(y)); ok <- is.finite(y) & y > 0
  if (any(ok)) { x[ok] <- 0.5 + 1 / y[ok]
    for (k in seq_len(50)) { tri <- trigamma(x[ok]); dif <- tri * (1 - tri / y[ok]) / psigamma(x[ok], deriv = 2); x[ok] <- x[ok] + dif
      if (max(abs(dif / x[ok]), na.rm = TRUE) < 1e-8) break } }
  x[is.finite(y) & y <= 0] <- Inf   # trigamma(x) -> 0 as x -> Inf
  x
}

#' Scaled-F prior of a set of variances
#'
#' The variances `v` (each on `df1` degrees of freedom) are modeled as `s0^2 F(df1, d0)`, as in limma's `fitFDist`, by the
#' moments of `log v`: `E log v = log s0^2 + digamma(df1/2) - log(df1/2) - digamma(d0/2) + log(d0/2)` and
#' `Var log v = trigamma(df1/2) + trigamma(d0/2)`. With a `covariate`, `log s0^2` follows a lowess trend of `log v` on it.
#' When the spread of the variances is no wider than their sampling noise the fit gives `d0 = Inf`, a moderated variance
#' that is the prior alone whatever a column's own variance says; `d0.max` caps the prior degrees of freedom so that a
#' column's own variance keeps a share of at least `df1 / (df1 + d0.max)`.
#' @param v variances (one per column; `NA`, zero or infinite values are left out of the fit)
#' @param df1 their degrees of freedom (one value or one per variance)
#' @param covariate optional per-column covariate of the trend (e.g. the column mean)
#' @param d0.max cap on the prior degrees of freedom (default `Inf`, none)
#' @param span lowess span of the trend (default 0.3)
#' @return list: `var` (prior variance per column), `df` (prior degrees of freedom, capped), `df.fit` (as fitted),
#'   `trend` (the fitted trend of `log v`), `n` (variances used)
#' @export
fitFPrior <- function(v, df1, covariate = NULL, d0.max = Inf, span = 0.3) {
  v <- as.numeric(v); ok <- is.finite(v) & v > 0
  if (sum(ok) < 10) stop("too few variances for the prior")
  df1 <- rep_len(as.numeric(df1), length(v)); ok <- ok & is.finite(df1) & df1 > 0
  z <- log(v[ok]); d1 <- df1[ok]
  trend <- rep(mean(z), length(v))
  if (!is.null(covariate)) {
    cv <- as.numeric(covariate); fin <- ok & is.finite(cv)
    if (sum(fin) >= 10 && stats::sd(cv[fin]) > 0) {
      o <- order(cv[fin]); lo <- stats::lowess(cv[fin][o], log(v[fin])[o], f = span, iter = 3)
      tr <- stats::approx(lo$x, lo$y, xout = cv, rule = 2, ties = mean)$y
      trend[is.finite(tr)] <- tr[is.finite(tr)]
    }
  }
  e <- z - trend[ok] - (digamma(d1 / 2) - log(d1 / 2))   # log v less its expectation under a common variance
  m <- mean(e); ex <- stats::var(e) - mean(trigamma(d1 / 2))
  d0.fit <- if (is.finite(ex) && ex > 0) 2 * trigammaInverse(ex) else Inf
  shift <- m + (if (is.finite(d0.fit)) digamma(d0.fit / 2) - log(d0.fit / 2) else 0)   # the scale from the fitted spread; the cap only limits the shrinkage
  list(var = exp(trend + shift), df = as.numeric(min(d0.fit, d0.max)), df.fit = as.numeric(d0.fit), trend = trend, n = sum(ok))
}

# Label-free prior of the moderated permutation statistic: per column, the total variance of the response after removing
# only the nuisance part of the design (Z with an intercept; the intercept alone without Z), which no relabeling changes.
# "fixed": prior variance = the median of these variances, prior df = d0; "squeeze": the scaled-F prior of fitFPrior()
# (prior variance and df estimated); "trend": the same with the prior variance a trend on the column mean; the estimated
# prior df is capped at d0.max. Returns list(var, df) for fit_and_randomize (df = 0: no moderation).
moderationPrior <- function(Y, Z, mode = "none", d0 = 3, d0.max = Inf) {
  if (identical(mode, "none")) return(list(var = numeric(0), df = 0))
  Zr <- if (is.null(Z)) matrix(1, nrow(Y), 1) else cbind(1, Z); qz <- qr(Zr); dft <- nrow(Y) - qz$rank
  if (dft < 1) stop("no residual degrees of freedom for the moderation prior")
  ok <- colSums(!is.finite(Y)) == 0
  v <- rep(NA_real_, ncol(Y)); v[ok] <- colSums(qr.resid(qz, Y[, ok, drop = FALSE])^2) / dft
  good <- is.finite(v) & v > 0
  if (sum(good) < 10) stop("too few columns with variance for the moderation prior")
  if (mode == "fixed") return(list(var = rep(stats::median(v[good]), ncol(Y)), df = as.numeric(d0)))
  pr <- fitFPrior(v, df1 = dft, covariate = if (mode == "trend") colMeans(Y) else NULL, d0.max = d0.max)
  s0 <- pr$var; s0[!is.finite(s0)] <- stats::median(s0[is.finite(s0)])
  list(var = as.numeric(s0), df = as.numeric(min(pr$df, 1e6)))
}


#' Permutation tests for linear models (block or Freedman-Lane)
#'
#' @description
#' One call to the C++ fitter `fit_and_randomize()` for a matrix of responses: every column is fit on the model's
#' design, a contrast of the coefficients is estimated and tested against the model's permutations
#' ([modelPermutations()]), which are shared with the distance-based tests of the same model.
#'
#' - `perm.method = "block"`: fit on the full design `F`; the relabelings permute rows of `F` within the plan's cells.
#' - `perm.method = "freedman-lane"`: residualize `X` and the responses on the nuisance design `Z` using all
#'   samples, fit on the residualized design and permute the residualized responses (Freedman-Lane). Every sample
#'   stays in the fit; the estimate equals the full-model coefficient (Frisch-Waugh-Lovell).
#'
#' The test statistic is the studentized contrast (`statistic = "t"`, the default: the estimate divided by its
#' standard error from each fit's own residual variance, so the permutation null does not inherit the variance of
#' the permuted fits) or the raw contrast estimate (`"coef"`). The raw estimate (`effect`) and its standard error
#' (`se`) are returned in either case as the effect size.
#'
#' @param x design from [buildDesignMatrices()] or a test design of [buildCacoaModel()] (`F`, `X`, `Z`, `contrast.F`,
#'   `contrast.X`, `meta`)
#' @param y numeric vector (`n`) or matrix (`n x m`) of responses; `NA` allowed
#' @param perm.method `"block"` or `"freedman-lane"`
#' @param robust.method `"none"`, `"huber"` or `"winsor"` (robust fits, re-estimated under every relabeling)
#' @param na.mode `"drop"` (rows with `NA` are dropped for that column) or `"impute_weak"` (kept with a near-zero weight)
#' @param na.center fill value of the weakly imputed rows: the observed mean (`"mean"`) or `"zero"`
#' @param statistic `"t"` (default) or `"coef"`
#' @param alternative `"two-sided"`, `"greater"` or `"less"`
#' @param cross.rows optional row index to subset the residual matrices after fitting
#' @param return.residuals,return.sampled.stats,return.sampled.fits optional outputs
#' @param n.permutations number of permutations (ignored when `P` is given)
#' @param n.cores threads (results do not depend on it)
#' @param seed seed of the permutation draw (default: taken from R's RNG, so `set.seed()` applies)
#' @param P,perm.cells permutations and cells from [modelPermutations()]; drawn from the design when `NULL`
#' @param block.vars additional strata for the permutation plan
#' @param moderate.perm moderated statistic inside the permutation test: `"none"` (default), `"fixed"` (residual variances
#'   shrunk towards the median label-free variance with `moderate.d0` prior degrees of freedom), `"squeeze"` (prior fitted by
#'   limma's `fitFDist` on the label-free variances) or `"trend"` (the same, prior variance a trend on the response mean);
#'   label-free = after removing only the nuisance design, so the test stays exact (closed-form path only)
#' @param moderate.d0 prior degrees of freedom of `"fixed"` (default 3)
#' @param moderate.d0.max cap on the estimated prior degrees of freedom of `"squeeze"` and `"trend"` (default `Inf`; see [fitFPrior()])
#' @return list: `coef` (`p x m`), `effect` and `se` (contrast estimate and standard error), `df`, `stat.obs` (the
#'   statistic), `pval`, `z.score` (normal quantile of `pval` with the sign of the effect), `stats.perm` and
#'   `effects.perm` (`B x m`, when requested), `sampled.fits`, `residuals`, `residuals.pearson`, `P`, `perm.cells`,
#'   `statistic`, `perm.method`
#' @seealso [buildDesignMatrices()], [modelPermutations()]
#' @keywords internal
performLMPermutations <- function(x, y,
                  perm.method = c("block","freedman-lane"),
                  robust.method = c("none", "huber", "winsor"),
                  na.mode = c("drop", "impute_weak"), na.center = c("mean", "zero"),
                  statistic = c("t", "coef"),
                  alternative = c("two-sided","greater","less"),
                  cross.rows = NULL,
                  return.residuals = TRUE,
                  return.sampled.stats = TRUE,
                  return.sampled.fits = FALSE,
                  n.permutations = 1000,
                  n.cores=1, seed = NULL, P = NULL, perm.cells = NULL, block.vars = NULL,
                  pooled = FALSE, z.grid = NULL, use.abs = TRUE, group = NULL, average = FALSE, nn.graph = NULL, closed.form = TRUE, pooled.standardize = FALSE,
                  moderate.perm = c("none", "fixed", "squeeze", "trend"), moderate.d0 = 3, moderate.d0.max = Inf) {

  robust.method <- match.arg(robust.method); moderate.perm <- match.arg(moderate.perm)
  perm.method   <- match.arg(perm.method)
  na.mode       <- match.arg(na.mode)
  na.center     <- match.arg(na.center)
  statistic     <- match.arg(statistic)
  alternative   <- match.arg(alternative)
  # responses are matched to the design's samples by name whenever both are named (positional pairing is unsafe)
  if (is.matrix(y) && !is.null(rownames(y)) && !is.null(rownames(x$F))) {
    miss <- setdiff(rownames(x$F), rownames(y))
    if (length(miss)) stop("response rows missing for sample(s): ", paste(utils::head(miss, 5), collapse = ", "))
    y <- y[rownames(x$F), , drop = FALSE]
  }
  # seed for the R-drawn permutations: taken from R's RNG when not given, so set.seed() makes runs reproducible
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  seed <- as.integer(seed)
  # permutations from the model's plan (shared with the distance-based tests) unless given explicitly
  if (is.null(P) && n.permutations > 0) {
    if (is.null(x$meta) || is.null(x$F)) stop("the design carries no sample metadata: pass P (see modelPermutations()) or build it with buildDesignMatrices() / buildCacoaModel()")
    mp <- modelPermutations(x, scheme = perm.method, n.permutations = n.permutations, block.vars = block.vars, seed = seed)
    P <- mp$P; perm.cells <- mp$cells
  }
  if (!is.null(P)) {
    storage.mode(P) <- "integer"
    if (nrow(P) != nrow(x$F)) stop("P must have one row per sample of the design (", nrow(x$F), ")")
    n.permutations <- ncol(P)
    if (is.null(perm.cells)) perm.cells <- list(seq_len(nrow(P)))
  }
  Y <- if (is.matrix(y)) y else matrix(y, ncol = 1L)
  if (!is.numeric(Y)) stop("y must be a numeric vector or matrix.")
  n.F <- nrow(x$F)
  if (nrow(Y) != n.F) stop("nrow(Y)=", nrow(Y), " != nrow(x$F)=", n.F, ".")
  y.names <- colnames(Y)
  if (is.null(y.names)) y.names <- paste0("Y", seq_len(ncol(Y)))

  # design: the full design for block relabeling; X with Z residualized out for Freedman-Lane
  fl <- perm.method == "freedman-lane"
  D <- as.matrix(if (fl) x$X else x$F)
  Z <- if (fl && !is.null(x$Z) && ncol(x$Z)) as.matrix(x$Z) else NULL
  cvec <- as.numeric(if (fl) x$contrast.X else x$contrast.F)
  mod <- moderationPrior(Y, Z, moderate.perm, moderate.d0, moderate.d0.max)
  fit <- fit_and_randomize(
    X = D, Y = Y, contrast = cvec, Z = Z,
    perm_groups = perm.cells, perm_matrix = P, n_randomizations = n.permutations,
    alternative = alternative, statistic = statistic,
    return_residuals = return.residuals, return_sampled_fits = return.sampled.fits, return_sampled_stats = return.sampled.stats,
    robust = robust.method, huber_k = 1.345, huber_maxit = 8, huber_tol = 1e-6,
    na_mode = na.mode, na_weight = 1e-4, na_center = na.center,
    illcond_rcond = 1e-12, pinv_tol = 0.0, n_cores = n.cores,
    pooled = pooled, z_grid = if (is.null(z.grid)) numeric(0) else as.numeric(z.grid), use_abs = use.abs,
    group = if (is.null(group)) integer(0) else as.integer(group), average = average, nn_graph = if (is.null(nn.graph)) NULL else unname(nn.graph), closed_form = closed.form, pooled_standardize = pooled.standardize,
    prior_var = mod$var, prior_df = mod$df
  )

  coef <- fit$coef
  stat <- as.vector(fit$stat); p <- as.vector(fit$p_value)
  effect <- as.vector(fit$effect); se <- as.vector(fit$se); df <- as.vector(fit$df)
  # z-like score from the permutation p-value, signed by the effect
  p.cl <- pmin(pmax(p, 1e-16), 1 - 1e-16)
  z <- switch(alternative, "two-sided" = qnorm(1 - p.cl / 2) * sign(effect), greater = qnorm(1 - p.cl), less = qnorm(p.cl))
  names(stat) <- names(z) <- names(p) <- names(effect) <- names(se) <- names(df) <- y.names
  if (is.matrix(coef)) dimnames(coef) <- list(colnames(D), y.names)

  ## --- permutation statistics and effects ---
  stats.perm <- effects.perm <- NULL
  if (pooled && !is.null(fit$pooled) && return.sampled.stats && length(fit$pooled$perm)) {   # the accumulator's matrix (averaged statistics under `average`)
    stats.perm <- fit$pooled$perm; dimnames(stats.perm) <- list(paste0("perm", seq_len(nrow(stats.perm))), y.names); stats.perm[!is.finite(stats.perm)] <- NA_real_
  }
  if (return.sampled.stats && !is.null(fit$sampled_stats)) {
    S <- fit$sampled_stats; dimnames(S) <- list(paste0("perm", seq_len(nrow(S))), y.names); S[!is.finite(S)] <- NA_real_
    E <- if (is.null(fit$sampled_effects)) S else fit$sampled_effects
    dimnames(E) <- dimnames(S); E[!is.finite(E)] <- NA_real_
    stats.perm <- S; effects.perm <- E
  }

  ## --- residuals (all samples; NA where the response was missing) and Pearson residuals ---
  residuals <- residuals.pearson <- NULL
  if (return.residuals && !is.null(fit$residuals)) {
    residuals <- fit$residuals
    rn <- rownames(Y); if (is.null(rn)) rn <- rownames(x$F); if (is.null(rn)) rn <- as.character(seq_len(nrow(Y)))
    dimnames(residuals) <- list(rn, y.names)
    sigma.hat <- sqrt(colSums(residuals^2, na.rm = TRUE) / pmax(1, df))
    sigma.hat[!is.finite(sigma.hat) | sigma.hat <= 0] <- NA_real_
    if (!is.null(cross.rows)) residuals <- residuals[cross.rows, , drop = FALSE]
    residuals.pearson <- sweep(residuals, 2, sigma.hat, "/")
  }

  ## --- failed columns ---
  failed <- !is.finite(stat) | (if (n.permutations > 0) !is.finite(p) else FALSE)
  failed <- failed | (if (is.matrix(coef)) apply(!is.finite(coef), 2L, any) else !is.finite(coef))
  if (any(failed)) {
    stat[failed] <- z[failed] <- p[failed] <- NA_real_
    if (is.matrix(coef)) coef[, failed] <- NA_real_ else coef[failed] <- NA_real_
    warning("Model fit failed for ", sum(failed), " out of ", length(stat), " columns of Y (too many missing values, ",
            "a constant response or a singular design); their statistics are NA. Consider na.mode = 'impute_weak'.", call. = FALSE)
  }

  list(
    P = P, perm.cells = perm.cells, statistic = statistic, perm.method = perm.method,
    coef = coef, effect = effect, se = se, df = df,
    stat.obs = stat, stats.perm = stats.perm, effects.perm = effects.perm,
    sampled.fits = fit$sampled_fits, pooled = fit$pooled, effect.perm.mean = setNames(as.numeric(fit$effect_perm_mean), y.names), closed.form = isTRUE(fit$closed_form),
    pval = p, z.score = z,
    residuals = residuals, residuals.pearson = residuals.pearson
  )
}

# helper to sanitize groups against F
sanitizeGroups <- function(groups, F) {
    if (is.null(groups)) {
      g <- as.list(colnames(F))
      names(g) <- colnames(F)
      return(g)
    } else {
      g2 <- lapply(groups, function(v) intersect(v, colnames(F)))
      ok <- lengths(g2) > 0
      if (!all(ok)) warning("Some groups had no matching columns in F and were dropped.")
      g2 <- g2[ok]
      if (!length(g2)) stop("No valid groups after intersecting with F's columns.")
      return(g2)
    }
  }


# Fit and get sum of squares (SSE) 
#' @keywords internal
getSSE <- function(X, y) {
    fit <- lm.fit(X, y)
    bhat <- fit$coefficients; bhat[!is.finite(bhat)] <- 0
    res  <- y - as.numeric(X %*% bhat)
    sum(res^2)
}


