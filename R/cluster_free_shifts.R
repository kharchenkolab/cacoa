## Cluster-free expression shifts on the pairwise-effects engine (Track D.3).
##
## For every cell, the samples' mean expression over the cell's neighbourhood gives a sample x sample distance
## matrix (computed in C++ by estimateExpressionShiftsPairsLM over all sample pairs). The contrast's shift is
## estimated per neighbourhood with the individual-level model (common dispersion), tested by permuting sample
## labels (one set of permutations shared by all cells, induced on each neighbourhood's sample subset), and the
## per-cell z-scores are adjusted by the max-statistic over cells (adjustedZScoresMaxStat) and smoothed over
## the graph (applyMedianFilterES).

# square sample distance matrix from one column of the pair layout. Samples with missing pair distances (too few
# cells in the neighbourhood) are dropped one at a time, the one with the most missing pairs first, until no missing
# pair is left: a thin sample removes itself, not the samples it is paired with (same rule as the C++ kernel).
pairVectorToMatrix <- function(y, pairs, samples) {
  n <- length(samples); D <- matrix(NA_real_, n, n, dimnames = list(samples, samples)); diag(D) <- 0
  D[pairs] <- y; D[pairs[, 2:1, drop = FALSE]] <- y
  ok <- rep(TRUE, n)
  repeat {
    miss <- rowSums(is.na(D[ok, ok, drop = FALSE]))
    if (!length(miss) || max(miss) == 0) break
    ok[which(ok)[which.max(miss)]] <- FALSE
  }
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
#' @param n.cores threads for the per-cell loop (C++)
#' @param verbose progress
#' @param robust `"none"`, `"huber"` or `"winsor"`: down-weight samples with outlying residual distances (re-estimated under every relabeling)
#' @param robust.k robust tuning constant
#' @return list: `stat` (shift F per cell), `p.value`, `z.score` (and `z.scores`), `z.adj`, `shifts` (bias-corrected shift
#'   estimate), `shifts.smoothed`, `n.samples`, `settings`
#' @export
clusterFreeExpressionShifts <- function(cm, sample.per.cell, nns.per.cell, design, meta, dist = c("cor", "cosine", "js"), min.n.obs.per.samp = 3,
                                        min.samp.per.level = 2, permutation = "auto", n.permutations = 199, block.vars = NULL, seed = 1,
                                        adjust = TRUE, smooth = TRUE, wins = 0.025, log.vectors = TRUE, n.cores = 1, verbose = FALSE,
                                        robust = c("none", "huber", "winsor"), robust.k = 1.345) {
  dist <- match.arg(dist); robust <- match.arg(robust)
  if (identical(design$contrast_spec$type, "term")) stop("cluster-free shifts need a contrast test (two groups or a numeric step)")
  samples <- intersect(rownames(design$F), levels(droplevels(factor(sample.per.cell))))
  if (length(samples) < 4) stop("fewer than 4 samples shared by the design and the cells")
  spc <- factor(as.character(sample.per.cell[colnames(cm)]), levels = samples)
  keep.cells <- !is.na(spc)
  nn.map <- NULL
  if (!all(keep.cells)) {                      # cells of samples outside the design are dropped: re-index the neighbourhoods
    cm <- cm[, keep.cells, drop = FALSE]; spc <- spc[keep.cells]
    nn.map <- cumsum(keep.cells) - 1L; nn.map[!keep.cells] <- NA_integer_
  }
  n <- length(samples)
  pairs <- t(utils::combn(n, 2)); storage.mode(pairs) <- "integer"
  m <- length(nns.per.cell)
  nn.list <- lapply(nns.per.cell, as.integer)
  if (!is.null(nn.map)) nn.list <- lapply(nn.list, function(v) { w <- nn.map[v + 1L]; w[!is.na(w)] })
  # global permutations over all samples; neighbourhood profiles, distances and the test per cell run in C++
  meta <- meta[samples, , drop = FALSE]
  des <- design; des$F <- design$F[samples, , drop = FALSE]
  gplan <- permutationPlan(des, meta, samples, scheme = permutation, block.vars = block.vars, n.permutations = n.permutations, max.enumerate = 0)
  P.global <- withSeed(seed, drawPermutations(gplan, n.permutations)); storage.mode(P.global) <- "integer"
  B <- ncol(P.global)
  cF <- design$contrast.F[colnames(des$F)]
  spec <- design$contrast_spec
  lev.var <- if (!is.null(spec) && spec$type %in% c("simple", "marginal") && !grepl(":", spec$term, fixed = TRUE) && spec$term %in% names(meta) && !is.numeric(meta[[spec$term]])) spec$term else NULL
  level.code <- rep(-1L, n)
  if (!is.null(lev.var)) { g <- as.character(meta[[lev.var]]); level.code <- ifelse(is.na(g), -1L, ifelse(g == spec$den, 1L, ifelse(g == spec$num, 2L, 0L))) }
  if (verbose) message(sprintf("Testing %d cells (%d samples) with %d permutations (%s)...", m, n, B, gplan$scheme))
  kr <- cluster_free_shift_stream(cm, as.integer(spc), nn.list, FALSE, as.integer(min.n.obs.per.samp), dist, log.vectors,
                                  pairs - 1L, n, des$F, cF, as.integer(level.code), as.integer(min.samp.per.level), as.integer(gplan$strata), as.integer(gplan$in.set),
                                  P.global, gplan$scheme == "freedman-lane", TRUE, as.integer(n.cores), c(none = 0L, huber = 1L, winsor = 2L)[[robust]], robust.k)
  Fobs <- as.numeric(kr$stat); pval <- as.numeric(kr$p); shifts <- as.numeric(kr$shift); n.samp <- as.integer(kr$n); z <- as.numeric(kr$z)
  mx <- as.numeric(kr$max); mn <- as.numeric(kr$min)
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
       settings = list(dist = dist, permutation = gplan$scheme, n.permutations = B, seed = seed, min.n.obs.per.samp = min.n.obs.per.samp, robust = robust,
                       adjust = adjust, smooth = smooth, n.cells.tested = length(valid)))
}
