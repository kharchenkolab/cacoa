
#' Permutation tests for linear models (block or Freedman–Lane)
#'
#' @description
#' Wrapper to run permutation-based inference for linear models using either
#' **block permutations** (shuffle within nuisance-defined blocks on the full design `F`)
#' or the **Freedman–Lane** scheme (residualize out nuisance `Z`, then permute).
#' Supports robust fitting and NA handling, and returns neatly named outputs.
#'
#' @param x A design bundle from [buildDesignMatrices()], containing at least:
#'   `F`, `X`, `Z`, `contrast.F`, `contrast.X`, `perm.groups.full`, `perm.groups.core`,
#'   and optionally `core.rows` and `pairs`.
#' @param y Numeric vector (`n`) or matrix (`n × m`) of responses.
#' @param n.permutations Integer, number of randomizations.
#' @param perm.method One of `"block"` or `"freedman-lane"`.
#' @param robust.method One of `"none"`, `"huber"`, or `"winsor"`. Passed through to the fitter.
#' @param na.mode One of `"drop"` or `"impute.weak"`.
#'   - `"drop"`: for each response column, rows with `NA` are dropped **for that column's fit**.
#'   - `"impute.weak"`: rows with `NA` are **kept**; nuisance residualization uses
#'     weighted LS (tiny weight on missing rows) to get Xr (falls back to NAs), and the fitter treats missing rows with
#'     a tiny weights similarly (observed rows are the ones permuted).
#' @param cross.rows Optional integer/logical index to subset residual rows to cross-group sample pairs 
#'   for visualization **after** fitting.
#'   Indices are with respect to the rows actually used: all rows for `"block"`,
#'   or `which(x$core.rows)` when `perm.method = "freedman-lane"` 
#' @param return.residuals Logical; return residual matrix.
#' @param return.sampled.stats Logical; return the matrix of sampled permutation statistics.
#'
#' @details
#' **Block**: calls the underlying C++ `fit_and_randomize()` on the full design `F`,
#' shuffling within `x$perm.groups.full`.
#'
#' **Freedman–Lane**: calls C++ helper `fl_fwl_cpp()`
#' Column names from `colnames(y)` are propagated to `stat.obs`, `pval`, `z.score`,
#' and (when requested) to permutation statistics and residual matrices.
#'
#' @return A list with:
#' \itemize{
#'   \item `stat.obs` — numeric length-`m` vector of observed contrast statistics.
#'   \item `pval` — numeric length-`m` vector of permutation p-values.
#'   \item `z.score` — numeric length-`m` vector of z-like scores derived from permutation p.
#'   \item `stats.perm` — if requested, a matrix `n.permutations × m` (when available), otherwise a list.
#'     Columns are named by `colnames(y)`; rows are `"perm1"`, `"perm2"`, ….
#'   \item `residual` — if requested, residual matrix:
#'     - `"block"`: `n × m` (all rows);
#'     - `"freedman-lane"`: `|core rows| × m`, otherwise `n × m`.
#'     Row names reflect the used rows; `cross.rows` may further subset rows.
#'   \item `y.resid` — for `"freedman-lane"`, the `|core rows| × m` matrix of residualized responses (if requested).
#' }
#'
#' @section Notes:
#' * For `"freedman-lane"`, if some columns are not estimable under a subset (e.g., too few rows),
#'   their outputs may be `NA`. 
#'
#' @seealso [buildDesignMatrices()], [modelPermutations()], `fit_and_randomize()`
#'
#' @examples
#' \dontrun{
#' # Assume 'x' from buildDesignMatrices(sample.meta, contrast = "group=B") or buildPairDesignMatrices()
#' set.seed(1)
#' Y <- matrix(rnorm(nrow(x$F) * 3L), ncol = 3L)
#' colnames(Y) <- c("g1","g2","g3")
#'
#' # Block permutations
#' out.block <- performLMPermutations(
#'   x, Y, n.permutations = 999, perm.method = "block",
#'   robust.method = "none", na.mode = "drop", return.sampled.stats = TRUE
#' )
#'
#' # Freedman–Lane on core rows with weak imputation
#' out.fl <- performLMPermutations(
#'   x, Y, n.permutations = 999, perm.method = "freedman-lane",
#'   na.mode = "impute_weak", return.residuals = TRUE
#' )
#' }
#' @keywords internal
performLMPermutations <- function(x, y,
                  perm.method = c("block","freedman-lane"),
                  robust.method = c("none", "huber", "winsor"),
                  na.mode = c("drop", "impute_weak"), na.center = c("mean","median"),
                  alternative = c("two-sided","greater","less"),
                  cross.rows = NULL,
                  return.residuals = TRUE,
                  return.sampled.stats = TRUE,
                  return.sampled.fits = FALSE,
                  return.y.resid = TRUE,
                  n.permutations = 1000,
                  n.cores=1, seed = NULL, P = NULL, perm.cells = NULL, block.vars = NULL) {

  robust.method <- match.arg(robust.method)
  perm.method   <- match.arg(perm.method)
  na.mode       <- match.arg(na.mode)
  alternative   <- match.arg(alternative)
  na.center     <- match.arg(na.center)
  # seed for the R-drawn permutations: taken from R's RNG when not given, so set.seed() makes runs reproducible
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  seed <- as.integer(seed)
  # permutations from the model's plan (shared with the distance-based tests) unless given explicitly;
  # the legacy internal generator is used only for designs without metadata
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
  # response vector/matrix
  if (is.numeric(y) && !is.matrix(y)) {
    Y <- matrix(y, ncol = 1L)
  } else if (is.matrix(y) && is.numeric(y)) {
    Y <- y
  } else {
    stop("y must be a numeric vector or matrix.")
  }
  # sanity checks
  n.F <- nrow(x$F)
  if (nrow(Y) != n.F) stop("nrow(Y)=", nrow(Y), " != nrow(x$F)=", n.F, ".")
  if (!is.null(x$X) && nrow(x$X) != n.F) stop("nrow(x$X) must equal nrow(x$F).")
  if (!is.null(x$Z) && nrow(x$Z) != n.F) stop("nrow(x$Z) must equal nrow(x$F).")
  
  y.names <- colnames(Y)
  if (is.null(y.names)) y.names <- paste0("Y", seq_len(ncol(Y)))
  
  y.resid <- NULL

  if (perm.method == "block") {
    # ----------------------------------------------------------
    # BLOCK: fit on F; randomize within blocks
    # ----------------------------------------------------------
    fit <- fit_and_randomize(
      X = x$F,
      Y = Y,
      contrast = x$contrast.F,
      perm_groups = perm.cells,
      perm_matrix = P,
      n_randomizations = n.permutations,
      alternative = alternative,
      return_residuals = return.residuals,
      return_sampled_fits = return.sampled.fits,
      return_sampled_stats = return.sampled.stats,
      robust = robust.method, huber_k = 1.345, huber_maxit = 8, huber_tol = 1e-6,
      na_mode = na.mode, na_weight = 1e-4, na_center = na.center,
      illcond_rcond = 1e-12, pinv_tol = 0.0, n_cores = n.cores
    )
  } else if (perm.method == "freedman-lane") {
    # ------------------------------------------------------------
    # FREEDMAN–LANE: NA-aware, per-column residualization
    # ------------------------------------------------------------
    fit <- fl_fwl_cpp(
      X = x$X,
      Z = if (is.null(x$Z)) matrix(0,0,0) else x$Z,
      Y = Y,
      contrast = x$contrast.X,
      core_rows = x$core.rows, 
      core_perm_groups = if (!is.null(P)) coreCells(perm.cells, x$core.rows %||% rep(TRUE, nrow(P))) else NULL,
      perm_matrix = P,
      n_randomizations = n.permutations,
      alternative = alternative,
      robust = robust.method,
      huber_k = 1.345, huber_maxit = 8, huber_tol = 1e-6,
      na_mode = na.mode, na_weight = 1e-4, na_center = na.center,
      illcond_rcond = 1e-12, pinv_tol = 0.0, n_cores = n.cores, return_residuals = return.residuals,
      return_sampled_fits=return.sampled.fits, return_sampled_stats = return.sampled.stats)

    y.resid <- if(return.y.resid) fit$partial_core else NULL
  }
  
  # extract and format results
  coef <- fit$coef
  stat <- as.vector(fit$stat)    
  z    <- as.vector(fit$z_score) 
  p    <- as.vector(fit$p_value) 

  
  # calculate z scores
  p_clamped <- pmin(pmax(p, 1e-16), 1 - 1e-16)
  if (alternative == "two-sided") {
    # For two-sided: p = 2 * (1 - pnorm(|Z|))  =>  |Z| = qnorm(1 - p/2)
    # We restore the sign from the original statistic (assuming symmetry under null)
    # Note: If stat is 0, sign is 0.
    z <- qnorm(1 - p_clamped / 2) * sign(stat)
  } else if (alternative == "greater") {
    # For greater: p = 1 - pnorm(Z)  =>  Z = qnorm(1 - p)
    z <- qnorm(1 - p_clamped)
  } else { 
    # For less: p = pnorm(Z)  =>  Z = qnorm(p)
    z <- qnorm(p_clamped)
  }
  
  names(stat) <- names(z) <- names(p) <- y.names
  
  ## dimnames for coef (p × m)
  if (is.matrix(coef)) {
    colnames(coef) <- y.names
    rownames(coef) <- if (perm.method == "block") colnames(x$F) else colnames(x$X)
  }
  
  ## --- permutation statistics ---
  stats.perm <- NULL
  max.perm <- min.perm <- NULL
  if (return.sampled.stats && !is.null(fit$sampled_stats)) {
    S <- fit$sampled_stats
    if (is.matrix(S)) {
      colnames(S) <- y.names
      rownames(S) <- paste0("perm", seq_len(nrow(S)))
      S[!is.finite(S)] <- NA_real_
      #perm.rng <- apply(S, 1, range, na.rm=TRUE)
      #min.perm <- perm.rng[1,]
      #max.perm <- perm.rng[2,]
    }
    stats.perm <- S
  }

  ## --- residuals: row/colnames & optional subsetting ---
  residuals <- NULL
  if (return.residuals && !is.null(fit$residuals)) {
    residuals <- fit$residuals
    colnames(residuals) <- y.names
    
    rn.all <- rownames(Y)
    if (is.null(rn.all)) rn.all <- as.character(seq_len(nrow(Y)))
    
    used.idx <- if (perm.method == "freedman-lane" && !is.null(x$core.rows)) {
      which(x$core.rows)
    } else {
      seq_len(nrow(Y))
    }
    rownames(residuals) <- rn.all[used.idx]
    
    if (!is.null(cross.rows)) {
      residuals <- residuals[cross.rows, , drop = FALSE]
    }
  }

  ## -- Pearson's residuals ---
  residuals.pearson <- NULL
  sigma.hat <- NULL
  df.resid  <- NULL

  if (return.residuals && !is.null(residuals)) {
    D <- if (perm.method == "block") as.matrix(x$F) else as.matrix(x$X)
    p.eff <- qr(D)$rank # effective parameter count

    # Per-response-column n used 
    n.obs <- colSums(is.finite(residuals))
    df.resid <- pmax(1, n.obs - p.eff) # Residual df per column
    rss <- colSums(residuals^2, na.rm = TRUE) # Sigma-hat per column
    sigma.hat <- sqrt(rss / df.resid)
    sigma.hat[!is.finite(sigma.hat) | sigma.hat <= 0] <- NA_real_

    # Pearson residuals: r / sigma
    residuals.pearson <- sweep(residuals, 2, sigma.hat, "/")
  }
  
  ## --- y.resid (partial_core) naming for FL ---
  if (!is.null(y.resid)) {
    colnames(y.resid) <- y.names
    rn.all <- rownames(Y)
    if (is.null(rn.all)) rn.all <- as.character(seq_len(nrow(Y)))
    used.idx <- if (!is.null(x$core.rows)) which(x$core.rows) else seq_len(nrow(Y))
    rownames(y.resid) <- rn.all[used.idx]
  }
  
  ## --- detect failed columns and warn ---
  if (n.permutations > 0) {
    failed_stat <- !is.finite(stat) | !is.finite(z) | !is.finite(p)
  } else {
    failed_stat <- !is.finite(stat)
  }
  
  if (is.matrix(coef)) {
    failed_coef <- apply(!is.finite(coef), 2L, any)
    failed      <- failed_stat | failed_coef
  } else {
    failed      <- failed_stat | !is.finite(coef)
  }
  
  n.failed <- sum(failed)
  
  if (n.failed > 0L) {
    stat[failed] <- NA_real_
    z[failed]    <- NA_real_
    p[failed]    <- NA_real_
    
    if (is.matrix(coef)) {
      coef[, failed] <- NA_real_
    } else {
      coef[failed] <- NA_real_
    }
    
    warning(
      "Model fit failed for ", n.failed, " out of ", length(stat),
      " columns of Y; statistics, p-values, and coefficients for these columns are NA/NaN. ",
      "This is usually due to too many missing values or a singular design. ",
      "Consider using perm.method = 'freedman-lane' and/or na.mode = 'impute_weak'.",
      call. = FALSE
    )
  }
  
  
  ## --- final result ---
  list(
    P = P, perm.cells = perm.cells,
    y.resid    = if(perm.method == "freedman-lane") y.resid else NULL,
    coef       = coef,
    stat.obs   = stat,
    stats.perm = stats.perm,
    sampled.fits = fit$sampled_fits,
    pval       = p,
    z.score    = z,
    residuals  = residuals,
    residuals.pearson = residuals.pearson
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


# cells of row indices restricted to the core rows, re-indexed to core positions (for the Freedman-Lane fitter)
coreCells <- function(cells, core.rows) {
  pos <- integer(length(core.rows)); pos[core.rows] <- seq_len(sum(core.rows))
  out <- lapply(cells, function(i) pos[i[core.rows[i]]])
  out[lengths(out) >= 2]
}
