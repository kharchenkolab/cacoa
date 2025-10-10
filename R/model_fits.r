
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
#'   and optionally `core.rows`.
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
#' @seealso [buildDesignMatrices()], `cpp_fl()`, `fit_and_randomize()`, [permutationGroups()]
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
                  n.cores=1) {

  robust.method <- match.arg(robust.method)
  perm.method   <- match.arg(perm.method)
  na.mode       <- match.arg(na.mode)
  alternative   <- match.arg(alternative)
  na.center     <- match.arg(na.center)
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

  if (perm.method == "block") {
  # ----------------------------------------------------------
  # BLOCK: fit on F; randomize within blocks
  # ----------------------------------------------------------
  fit <- fit_and_randomize(
      X = x$F,
      Y = Y,
      contrast = x$contrast.F,
      perm_groups = x$perm.groups.full,
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
  core_rows = if (is.null(x$core.rows)) NULL else x$core.rows,  # logical 
  n_randomizations = n.permutations,
  alternative = alternative,
  robust = robust.method,
  huber_k = 1.345, huber_maxit = 8, huber_tol = 1e-6,
  na_mode = na.mode, na_weight = 1e-4, na_center = na.center,
  illcond_rcond = 1e-12, pinv_tol = 0.0, n_cores = n.cores, return_residuals = return.residuals,
  return_sampled_fits=return.sampled.fits, return_sampled_stats = return.sampled.stats)
  y.resid <- if(return.y.resid) fit$partial_core else NULL
  } else {
    stop("Unknown permutation method: ", perm.method)
  }
   
  # extract and format results
  coef <- as.numeric(fit$coef)
  stat <- as.numeric(fit$stat)
  z    <- as.numeric(if (!is.null(fit$z.score)) fit$z.score else fit$z_score)
  p    <- as.numeric(if (!is.null(fit$p.value)) fit$p.value else fit$p_value)
  #names(stat) <- names(z) <- names(p) <- y.names

  stats.perm <- NULL
  if (return.sampled.stats) {
    S <- fit$sampled_stats
    if (!is.null(S)) {
      if (is.matrix(S) && nrow(S) == n.permutations) {
        colnames(S) <- y.names
        rownames(S) <- paste0("perm", seq_len(nrow(S)))
        stats.perm <- S
      } else {
        stats.perm <- S
        if (!is.matrix(stats.perm)) names(stats.perm) <- y.names
      }
    }
  }
  #block: residuals are for all rows. freedman-lane: residuals are for core rows only (i.e., which(x$core.rows)).
  residuals <- NULL
  if (return.residuals && !is.null(fit$residuals)) {
    residuals <- fit$residuals
    colnames(residuals) <- y.names
    # set rownames depending on method
    rn.all <- rownames(Y); if (is.null(rn.all)) rn.all <- as.character(seq_len(nrow(Y)))
    used.idx <- if (perm.method == "freedman-lane" && !is.null(x$core.rows)) {
      which(x$core.rows)
    } else {
      seq_len(nrow(Y))
    }
    rownames(residuals) <- rn.all[used.idx]
  }
  
  list(
    y.resid    = if(perm.method == "freedman-lane") y.resid else NULL,
    coef       = coef,
    stat.obs   = stat,
    stats.perm = stats.perm,
    pval       = p,
    z.score    = z,
    residuals  = residuals
  )
}

## Compute unique R^2 for each term/group in a linear model
# - F: full model matrix (n x p)
# - y: numeric response vector (length n) OR matrix Y (n x K) with columns = cell types
# - groups: list of character vectors, each naming columns of F that form a term/group.
# - idx, sample.meta, sample.id, pair.sep: optional helpers to align Y->F row order when y is a matrix.
#   If provided, expected rownames of Y are built as paste(sample.meta[[sample.id]][idx$i], sample.meta[[sample.id]][idx$j], sep=pair.sep)
#   If not provided, function assumes Y is already in the same row order as F.
#' @keywords internal
estimateR2PerTerm <- function(F, y, groups = NULL,
                              idx = NULL,           # data.frame(i,j) for design order
                              pairs.Y = NULL,       # OPTIONAL: data.frame(i,j) giving Y row order
                              pair.sep = "__", verbose=FALSE) {
  stopifnot(is.matrix(F) || is.data.frame(F))
  F <- as.matrix(F)

  # --- original scalar-y path (unchanged math) ---
  if (is.vector(y) || (is.numeric(y) && is.null(dim(y)))) {
    y <- as.numeric(y)
    has.intercept <- "(Intercept)" %in% colnames(F)
    TSS <- if (has.intercept) sum((y - mean(y))^2) else sum(y^2)
    if (TSS <= 0) stop("TSS is zero; cannot compute R^2.")

    SSE.full  <- getSSE(F, y)
    R2.full   <- 1 - SSE.full / TSS
    groups2   <- sanitizeGroups(groups, F)

    out <- lapply(names(groups2), function(g) {
      cols.g     <- groups2[[g]]
      cols.other <- setdiff(colnames(F), cols.g)
      if (!length(cols.other)) {
        SSE.other <- if (has.intercept) sum((y - mean(y))^2) else sum(y^2)
      } else {
        SSE.other <- getSSE(F[, cols.other, drop = FALSE], y)
      }
      R2.unique <- max(0, (SSE.other - SSE.full) / TSS)
      data.frame(term = g, R2.unique = R2.unique, R2.full = R2.full,
                 stringsAsFactors = FALSE)
    })
    df <- do.call(rbind, out)
    # Drop intercept 
    df <- df[df$term != "Intercept", , drop = FALSE]
    rownames(df) <- NULL
    return(df)
  }

  # --- matrix-Y path (pairs x celltypes) ---
  stopifnot(is.matrix(y) || is.data.frame(y))
  Y <- as.matrix(y)

  # 1) Align Y to idx by indices
  if (!is.null(idx)) {
    stopifnot(is.data.frame(idx), all(c("i","j") %in% names(idx)))
    key.idx <- paste(idx$i, idx$j, sep=":")

    if (!is.null(pairs.Y)) {
      stopifnot(is.data.frame(pairs.Y), all(c("i","j") %in% names(pairs.Y)))
      key.y  <- paste(pairs.Y$i, pairs.Y$j, sep=":")
      map <- match(key.idx, key.y)
      if (anyNA(map)) stop("pairs.Y does not cover all rows in idx; cannot align Y.")
      Y <- Y[map, , drop = FALSE]
    } else if (!is.null(rownames(Y))) {
      # try to parse numeric i,j from rownames(Y)
      rn <- rownames(Y)
      # split on any non-digit (e.g., "__", "___", ":", ",", "-")
      ij <- lapply(strsplit(rn, "[^0-9]+"), function(v) as.integer(v[nzchar(v)]))
      ok <- vapply(ij, function(v) length(v) >= 2L && all(is.finite(v[1:2])), logical(1))
      if (all(ok)) {
        iY <- vapply(ij, function(v) v[1], integer(1))
        jY <- vapply(ij, function(v) v[2], integer(1))
        key.y <- paste(iY, jY, sep=":")
        map <- match(key.idx, key.y)
        if (anyNA(map)) {
          if (nrow(Y) != nrow(F)) stop("Failed to align Y by parsing rownames; and nrow(Y) != nrow(F).")
          if (verbose) warning("Could not fully align Y by parsed rownames; assuming Y already in F/idx order.")
        } else {
          Y <- Y[map, , drop = FALSE]
        }
      } else {
        if (nrow(Y) != nrow(F)) stop("Y rownames not parseable and nrow(Y) != nrow(F).")
        if (verbose) warning("Y rownames not parseable to (i,j); assuming Y already in F/idx order.")
      }
    } else {
      if (nrow(Y) != nrow(F)) stop("Y has no rownames and no pairs.Y; nrow(Y) must equal nrow(F).")
      # assume already aligned
    }
  } else {
    # No idx given: require Y already aligned to F
    if (nrow(Y) != nrow(F))
      stop("idx not provided; ensure Y has same row order and nrow as F.")
  }

  groups2 <- sanitizeGroups(groups, F)
  has.intercept <- "(Intercept)" %in% colnames(F)

  # 2) Fit per cell type (drop NA rows per column)
  ct.names <- colnames(Y)
  res <- lapply(ct.names, function(ct) {
    y.ct <- Y[, ct]
    keep <- which(is.finite(y.ct))
    if (length(keep) < 2L) {
      return(data.frame(term = names(groups2), celltype = ct,
                        R2.unique = NA_real_, R2.full = NA_real_,
                        stringsAsFactors = FALSE))
    }
    F.sub <- F[keep, , drop = FALSE]
    y.sub <- y.ct[keep]

    TSS <- if (has.intercept) sum((y.sub - mean(y.sub))^2) else sum(y.sub^2)
    if (TSS <= 0) {
      return(data.frame(term = names(groups2), celltype = ct,
                        R2.unique = NA_real_, R2.full = NA_real_,
                        stringsAsFactors = FALSE))
    }

    SSE.full <- getSSE(F.sub, y.sub)
    R2.full  <- 1 - SSE.full / TSS

    out.ct <- lapply(names(groups2), function(g) {
      cols.g     <- groups2[[g]]
      cols.other <- setdiff(colnames(F.sub), cols.g)
      if (!length(cols.other)) {
        SSE.other <- if (has.intercept) sum((y.sub - mean(y.sub))^2) else sum(y.sub^2)
      } else {
        SSE.other <- getSSE(F.sub[, cols.other, drop = FALSE], y.sub)
      }
      R2.unique <- max(0, (SSE.other - SSE.full) / TSS)
      data.frame(term = g, celltype = ct,
                 R2.unique = R2.unique, R2.full = R2.full,
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, out.ct)
  })

  df <- do.call(rbind, res)
  # Drop intercept
  df <- df[df$term != "Intercept", , drop = FALSE]
  rownames(df) <- NULL
  df
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