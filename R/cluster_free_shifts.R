## Cluster-free expression shifts on the pairwise-effects engine (Track D.3).
##
## For every cell, the samples' mean expression over the cell's neighbourhood gives a sample x sample distance
## matrix (computed in C++ by estimateExpressionShiftsPairsLM over all sample pairs). The contrast's shift is
## estimated per neighbourhood with the individual-level model (common dispersion), tested by permuting sample
## labels (one set of permutations shared by all cells, induced on each neighbourhood's sample subset), and the
## per-cell z-scores are adjusted by the max-statistic over cells (adjustedZScoresMaxStat) and smoothed over
## the graph (applyMedianFilterES).

# square sample distance matrix from one column of the pair layout
pairVectorToMatrix <- function(y, pairs, samples) {
  n <- length(samples); D <- matrix(NA_real_, n, n, dimnames = list(samples, samples)); diag(D) <- 0
  D[pairs] <- y; D[pairs[, 2:1, drop = FALSE]] <- y
  ok <- rowSums(is.na(D)) == 0
  D[ok, ok, drop = FALSE]
}

#' Cluster-free expression shifts for one contrast
#'
#' @param cm genes x cells sparse count matrix
#' @param sample.per.cell factor of sample per cell (named by cell, aligned to `colnames(cm)`)
#' @param nns.per.cell list (per cell) of neighbour indices (0-based, as returned by the graph adjacency)
#' @param design sample-level design of a contrast test (`buildDesignMatrices()` output or a `cacoaModel`)
#' @param meta sample metadata
#' @param dist `"cor"` (1 - Pearson correlation between the two samples' mean profiles; default), `"cosine"` or `"js"`
#' @param min.n.obs.per.samp minimum cells of a sample inside the neighbourhood (default 3)
#' @param min.samp.per.level minimum samples per compared level (default 2)
#' @param permutation,n.permutations,block.vars,seed permutation settings (see [permutationPlan()])
#' @param adjust max-statistic adjustment of the z-scores across cells (default TRUE)
#' @param smooth median-filter the shifts and adjusted z-scores over the graph (default TRUE)
#' @param wins winsorizing fraction for the adjustment (default 0.025)
#' @param log.vectors log10(1e3 x + 1) transform of the mean profiles (default TRUE)
#' @param n.cores cores for the per-cell loop
#' @param verbose progress
#' @return list: `stat` (shift F per cell), `p.value`, `z.score` (and `z.scores`), `z.adj`, `shifts` (bias-corrected shift
#'   estimate), `shifts.smoothed`, `n.samples`, `settings`
#' @export
clusterFreeExpressionShifts <- function(cm, sample.per.cell, nns.per.cell, design, meta, dist = c("cor", "cosine", "js"), min.n.obs.per.samp = 3,
                                        min.samp.per.level = 2, permutation = "auto", n.permutations = 199, block.vars = NULL, seed = 1,
                                        adjust = TRUE, smooth = TRUE, wins = 0.025, log.vectors = TRUE, n.cores = 1, verbose = FALSE) {
  dist <- match.arg(dist)
  if (identical(design$contrast_spec$type, "term")) stop("cluster-free shifts need a contrast test (two groups or a numeric step)")
  samples <- intersect(rownames(design$F), levels(droplevels(factor(sample.per.cell))))
  if (length(samples) < 4) stop("fewer than 4 samples shared by the design and the cells")
  spc <- factor(as.character(sample.per.cell[colnames(cm)]), levels = samples)
  keep.cells <- !is.na(spc)
  if (!all(keep.cells)) { cm <- cm[, keep.cells, drop = FALSE]; spc <- spc[keep.cells] }
  n <- length(samples)
  pairs <- t(utils::combn(n, 2)); storage.mode(pairs) <- "integer"
  m <- length(nns.per.cell)
  nn.list <- lapply(nns.per.cell, as.integer)
  if (verbose) message(sprintf("Computing neighbourhood sample distances for %d cells (%d samples)...", m, n))
  Y <- estimateExpressionShiftsPairsLM(cm = cm, sample_per_cell = as.integer(spc), nn_ids = nn.list, pairs_mat = pairs,
                                       min_n_obs_per_samp = min.n.obs.per.samp, dist = dist, log_vecs = log.vectors)
  # global permutations over all samples
  meta <- meta[samples, , drop = FALSE]
  des <- design; des$F <- design$F[samples, , drop = FALSE]
  gplan <- permutationPlan(des, meta, samples, scheme = permutation, block.vars = block.vars, n.permutations = n.permutations, max.enumerate = 0)
  P.global <- withSeed(seed, drawPermutations(gplan, n.permutations))
  B <- ncol(P.global)
  cF <- design$contrast.F[colnames(des$F)]
  spec <- design$contrast_spec
  lev.var <- if (!is.null(spec) && spec$type %in% c("simple", "marginal") && !grepl(":", spec$term, fixed = TRUE) && spec$term %in% names(meta) && !is.numeric(meta[[spec$term]])) spec$term else NULL
  plan.cache <- new.env()

  oneCell <- function(k) {
    out <- list(F = NA_real_, p = NA_real_, shift = NA_real_, n = 0L, perm = NULL)
    y <- Y[, k]
    if (all(is.na(y))) return(out)
    D <- pairVectorToMatrix(y, pairs, samples)
    s <- rownames(D)
    if (length(s) < 3) return(out)
    if (!is.null(lev.var)) { g <- as.character(meta[s, lev.var]); tb <- table(g)[c(spec$den, spec$num)]; if (any(is.na(tb)) || min(tb) < min.samp.per.level) return(out) }
    X <- des$F[s, , drop = FALSE]; keep <- colSums(abs(X)) > 1e-12
    if (any(abs(cF[!keep]) > 1e-12)) return(out)
    X <- X[, keep, drop = FALSE]; cvec <- cF[keep]
    if (ncol(X) >= nrow(X) || !isEstimable(X, cvec)) return(out)
    G <- gowerCenter(asSquaredDistance(D, "cor"))
    pre <- inferencePrecompute(G, X, matrix(1, nrow(X), 1), cvec)
    eff <- tryCatch(estimatePairwiseEffects(NULL, X, cvec, NULL, NULL, TRUE, G = G), error = function(e) NULL)
    key <- paste(s, collapse = ";")
    plan <- plan.cache[[key]]
    if (is.null(plan)) { plan <- permutationPlan(des, meta, s, scheme = gplan$scheme, block.vars = block.vars, n.permutations = n.permutations, max.enumerate = 0); assign(key, plan, envir = plan.cache) }
    sub <- match(s, samples)
    P <- apply(P.global, 2, inducePermutation, sub = sub, sub.plan = plan); storage.mode(P) <- "integer"
    Fp <- if (gplan$scheme == "freedman-lane") {
      parts <- flGowerParts(G, X %*% contrastNullBasis(cvec))
      permuted_contrast_F_fl(parts$K1, parts$K2, parts$K3, parts$K4, pre$a, pre$H, pre$cXc, pre$q, P)
    } else permuted_contrast_F(G, pre$a, pre$H, pre$cXc, pre$q, P)
    Fo <- contrastF(G, hatInfo(X), cvec)
    list(F = Fo, p = (sum(Fp >= Fo - 1e-12) + 1) / (B + 1), shift = if (is.null(eff)) NA_real_ else eff$shift, n = length(s), perm = as.numeric(Fp))
  }
  if (verbose) message(sprintf("Testing %d cells with %d permutations (%s)...", m, B, gplan$scheme))
  cells <- seq_len(m)
  res <- if (n.cores > 1) sccore::plapply(cells, oneCell, n.cores = n.cores, progress = verbose, fail.on.error = TRUE) else lapply(cells, oneCell)
  Fobs <- vapply(res, `[[`, numeric(1), "F"); pval <- vapply(res, `[[`, numeric(1), "p"); shifts <- vapply(res, `[[`, numeric(1), "shift")
  n.samp <- vapply(res, `[[`, integer(1), "n")
  valid <- which(is.finite(Fobs))
  # z-scores per cell and max-statistic extremes per permutation (same scale as z)
  z <- rep(NA_real_, m); mx <- rep(-Inf, B); mn <- rep(Inf, B)
  for (k in valid) {
    Fp <- res[[k]]$perm; mu <- mean(Fp); sdv <- stats::sd(Fp)
    if (!is.finite(sdv) || sdv == 0) { z[k] <- NA_real_; next }
    z[k] <- (Fobs[k] - mu) / sdv
    zp <- (Fp - mu) / sdv
    mx <- pmax(mx, zp); mn <- pmin(mn, zp)
  }
  valid <- which(is.finite(z))
  cell.names <- names(nns.per.cell) %||% colnames(cm)
  names(Fobs) <- names(pval) <- names(shifts) <- names(z) <- names(n.samp) <- cell.names
  z.adj <- NULL
  if (adjust && length(valid)) {
    z0 <- z; z0[!is.finite(z0)] <- 0
    z.adj <- adjustedZScoresMaxStat(z_obs = z0, alt = 1L, max_vals_in = mx, min_vals_in = mn, wins = wins, smooth = smooth,
                                    nn_ids = nn.list, non_zero_ids = as.integer(valid - 1L))
    z.adj[!is.finite(z)] <- NA_real_
    names(z.adj) <- cell.names
  }
  shifts.smoothed <- NULL
  if (smooth) {
    s0 <- shifts; s0[!is.finite(s0)] <- 0
    shifts.smoothed <- applyMedianFilterES(s0, nn_ids = nn.list, non_zero_ids = as.integer(which(is.finite(shifts)) - 1L))
    shifts.smoothed[!is.finite(shifts)] <- NA_real_
    names(shifts.smoothed) <- cell.names
  }
  list(stat = Fobs, p.value = pval, z.score = z, z.scores = z, z.adj = z.adj, shifts = shifts, shifts.smoothed = shifts.smoothed, n.samples = n.samp,
       settings = list(dist = dist, permutation = gplan$scheme, n.permutations = B, seed = seed, min.n.obs.per.samp = min.n.obs.per.samp,
                       adjust = adjust, smooth = smooth, n.cells.tested = length(valid)))
}
