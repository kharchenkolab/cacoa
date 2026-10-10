
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
                  n.cores=1, seed = NULL, P = NULL, perm.cells = NULL, block.vars = NULL) {

  robust.method <- match.arg(robust.method)
  perm.method   <- match.arg(perm.method)
  na.mode       <- match.arg(na.mode)
  na.center     <- match.arg(na.center)
  statistic     <- match.arg(statistic)
  alternative   <- match.arg(alternative)
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
  fit <- fit_and_randomize(
    X = D, Y = Y, contrast = cvec, Z = Z,
    perm_groups = perm.cells, perm_matrix = P, n_randomizations = n.permutations,
    alternative = alternative, statistic = statistic,
    return_residuals = return.residuals, return_sampled_fits = return.sampled.fits, return_sampled_stats = return.sampled.stats,
    robust = robust.method, huber_k = 1.345, huber_maxit = 8, huber_tol = 1e-6,
    na_mode = na.mode, na_weight = 1e-4, na_center = na.center,
    illcond_rcond = 1e-12, pinv_tol = 0.0, n_cores = n.cores
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
    sampled.fits = fit$sampled_fits,
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


