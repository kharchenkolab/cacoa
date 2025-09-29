
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
#' @param core.only Logical. For `"freedman-lane"`, fit is restricted to `x$core.rows` if `TRUE`.
#'   Ignored for `"block"`.
#' @param cross.rows Optional integer/logical index to subset residual rows to cross-group sample pairs 
#'   for visualization **after** fitting.
#'   Indices are with respect to the rows actually used: all rows for `"block"`,
#'   or `which(x$core.rows)` when `perm.method = "freedman-lane"` and `core.only = TRUE`.
#' @param return.residuals Logical; return residual matrix.
#' @param return.sampled.stats Logical; return the matrix of sampled permutation statistics.
#'
#' @details
#' **Block**: calls the underlying C++ `fit_and_randomize()` on the full design `F`,
#' shuffling within `x$perm.groups.full`.
#'
#' **Freedman–Lane**: calls C++ helper `cpp_fl()` which:
#' 1) residualizes `X` against `Z` on the **full** data once,
#' 2) handles `NA` per response column:
#'    - batches no-`NA` columns together,
#'    - for `NA` columns with `na.mode = "drop"`, groups by identical row mask and refits on that subset,
#'    - for `na.mode = "impute.weak"`, performs weighted residualization and preserves `NA` flags so the fitter
#'      can apply tiny-weight handling and permute only observed rows,
#' 3) then restricts to `x$core.rows` if `core.only = TRUE`, remapping `x$perm.groups.core` accordingly.
#' 4) fits and randomizes each column using `fit_and_randomize()`.
#'
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
#'   \item `y.resid` — if requested, residual matrix:
#'     - `"block"`: `n × m` (all rows);
#'     - `"freedman-lane"`: `|core rows| × m` if `core.only = TRUE`, otherwise `n × m`.
#'     Row names reflect the used rows; `cross.rows` may further subset rows.
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
#'   core.only = TRUE, na.mode = "impute_weak", return.residuals = TRUE
#' )
#' }
#' @keywords internal
performLMPermutations <- function(x, y,
                  n.permutations = 1000,
                  perm.method = c("block","freedman-lane"),
                  robust.method = c("none", "huber", "winsor"),
                  na.mode = c("drop", "impute_weak"),
                  core.only = TRUE,
                  cross.rows = NULL,
                  return.residuals = FALSE,
                  return.sampled.stats = TRUE) {

  robust.method <- match.arg(robust.method)
  perm.method   <- match.arg(perm.method)
  na.mode       <- match.arg(na.mode)
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
      alternative = "two-sided",
      return_residuals = return.residuals,
      return_sampled_fits = FALSE,
      return_sampled_stats = return.sampled.stats,
      robust = robust.method, huber_k = 1.345, huber_maxit = 8, huber_tol = 1e-6,
      na_mode = na.mode, na_weight = 1e-4, na_center = "mean",
      illcond_rcond = 1e-12, pinv_tol = 0.0, n_cores = 1
    )
  } else if (perm.method == "freedman-lane") {
  # ------------------------------------------------------------
  # FREEDMAN–LANE: NA-aware, per-column residualization
  # ------------------------------------------------------------
  fit <- cpp_fl(
  X = x$X,
  Z = if (is.null(x$Z)) matrix(0,0,0) else x$Z,
  Y = Y,
  contrast = x$contrast.X,
  core_rows_opt = if (is.null(x$core.rows)) NULL else x$core.rows,  # logical
  perm_groups_core_opt = x$perm.groups.core,                        # 1-based core space
  n_randomizations = n.permutations,
  alternative = "two-sided",
  return_residuals = TRUE,
  return_sampled_stats = TRUE,
  robust = robust.method,
  huber_k = 1.345, huber_maxit = 8, huber_tol = 1e-6,
  na_mode = na.mode, na_weight = 1e-4, na_center = "mean",
  illcond_rcond = 1e-12, pinv_tol = 0.0, n_cores = 1)
  } else {
    stop("Unknown permutation method: ", perm.method)
  }
   
  # extract and format results
  stat <- as.numeric(fit$stat)
  z    <- as.numeric(if (!is.null(fit$z.score)) fit$z.score else fit$z_score)
  p    <- as.numeric(if (!is.null(fit$p.value)) fit$p.value else fit$p_value)
  #names(stat) <- names(z) <- names(p) <- y.names

  stats.perm <- NULL
  if (return.sampled.stats) {
    S <- if (!is.null(fit$sampled.stats)) fit$sampled.stats else fit$sampled_stats
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
  y.resid <- NULL
  if (return.residuals && !is.null(fit$residuals)) {
    y.resid <- fit$residuals
    colnames(y.resid) <- y.names
    # set rownames depending on method 
    rn.all <- rownames(Y); if (is.null(rn.all)) rn.all <- as.character(seq_len(nrow(Y)))
    used.idx <- if (perm.method == "freedman-lane" && core.only && !is.null(x$core.rows)) {
      which(x$core.rows)
    } else {
      seq_len(nrow(Y))
    }
    rownames(y.resid) <- rn.all[used.idx]

    # optional cross.rows slice (indices relative to the used set)
    #if (!is.null(cross.rows)) y.resid <- y.resid[cross.rows, , drop = FALSE]
  }

  list(
    stat.obs   = stat,
    stats.perm = stats.perm,
    pval       = p,
    z.score    = z,
    y.resid    = y.resid
  )
}

## Compute unique R^2 for each term/group in a linear model
# - F: full model matrix (n x p)
# - y: numeric response vector (length n)
# - groups: list of character vectors, each naming columns of F that form a term/group. compute using .makeGroupsPair() or similar.
#' @keywords internal
estimateR2PerTerm <- function(F, y, groups = NULL) {
    stopifnot(is.matrix(F) || is.data.frame(F))
    F <- as.matrix(F)
    y <- as.numeric(y)
    has.intercept <- "(Intercept)" %in% colnames(F)
    
    # Total Sum of Squares (center if intercept, else uncentered)
    TSS <- if (has.intercept) sum((y - mean(y))^2) else sum(y^2)
    if (TSS <= 0) stop("TSS is zero; cannot compute R^2.")

    # Full model SSE and rank (for reference)
    SSE.full <- getSSE(F, y)
    R2.full  <- 1 - SSE.full / TSS
    rank.full <- qr(F)$rank

    # Default groups = each column as its own term
    if (is.null(groups)) {
        groups <- as.list(colnames(F))
        names(groups) <- colnames(F)
    } else {
        # sanity: drop any names not present
        groups <- lapply(groups, function(v) intersect(v, colnames(F)))
        ok <- lengths(groups) > 0
        if (!all(ok)) warning("Some groups had no matching columns in F and were dropped.")
        groups <- groups[ok]
    }
    
    # Compute unique R^2 by refitting without each group's columns
    out <- lapply(names(groups), function(g) {
        cols.g <- groups[[g]]
        cols.other <- setdiff(colnames(F), cols.g)
        if (!length(cols.other)) {
            # Removing the only columns leaves empty model: SSE_other = TSS (if no intercept),
            # or SSE w.r.t intercept-only. Handle gracefully:
            F.other <- matrix(1, nrow(F), 1); colnames(F.other) <- "(Intercept*)"
            SSE.other <- if (has.intercept) sum((y - mean(y))^2) else sum(y^2)
            rank.other <- 1L
        } else {
            F.other <- F[, cols.other, drop = FALSE]
            SSE.other <- getSSE(F.other, y)
            rank.other <- qr(F.other)$rank
        }
        R2.unique <- max(0, (SSE.other - SSE.full) / TSS)  # clip tiny negatives from numeric noise
        data.frame(term = g,
                   #k = length(cols.g),
                   #rank.minus = rank.other,
                   R2.unique = R2.unique,
                   R2.full = R2.full,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

# Fit and get sum of squares (SSE) 
#' @keywords internal
getSSE <- function(X, y) {
    fit <- lm.fit(X, y)
    bhat <- fit$coefficients; bhat[!is.finite(bhat)] <- 0
    res  <- y - as.numeric(X %*% bhat)
    sum(res^2)
}

# Summarize a list of LM permutation test results (from performLMPermutations)
# into a data.frame with p-values, adjusted p-values, z-scores, and observed
# statistics, plus lists of p-value distribution info, distance matrices per type,
# and R^2 values per type.
# If residuals were returned, also include a list of residual matrices per type.
# If permutation stats were returned, also include a list of permutation stats matrices per type.
#' @param res.list List of results from performLMPermutations(), named by cell type
#' @param adj.method Method for p-value adjustment (default: "BH")
#' @section Notes:
#' * Written only for cluster-based so far. 
#' @keywords internal
summarizeLMResults <- function(res.list, adj.method="BH") {
    ct   <- names(res.list)
    eff  <- sapply(res.list, \(x) unname(x$res$stat.obs[1]))
    pvalues <- sapply(res.list, \(x) x$res$pval)
    zscores <- sapply(res.list, \(x) x$res$z.score)
    if(!is.null(res.list[[1]]$res$y.resid)) { residuals <- sapply(res.list, \(x) x$res$y.resid, simplify=FALSE) }
    out  <- data.frame(celltype=ct, obs.stat=eff, pvalue=pvalues, zscore=zscores, stringsAsFactors=FALSE)
    padjust <- p.adjust(pvalues, method=adj.method)
    out$p.adj <- padjust
    out <- out[order(out$obs.stat, decreasing=TRUE), ]
    p.dist.info    <- lapply(res.list, function(x) x$dist.mat)
    dists.per.type <- lapply(res.list, function(x) x$dists)
    r2 <- lapply(res.list, function(x) x$r2)
    #if(!is.null(res.list[[1]]$res$y.resid)) { 
    #  residuals <- lapply(res.list, function(x) x$res$y.resid)
    #}
    if(!is.null(res.list[[1]]$res$stats.perm)) { 
      stats.perm <- lapply(res.list, function(x) x$res$stats.perm)
    }
    ret <- list(
      results = out,
      p.dist.info = p.dist.info,
      dists.per.type = dists.per.type,
      r2 = r2,
      pvalues = pvalues,
      padjust = padjust,
      zscores = zscores
    )
    if (exists("residuals")) ret$residuals <- residuals
    if (exists("stats.perm")) ret$stats.perm <- stats.perm
    ret
}

