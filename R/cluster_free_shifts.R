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
#' @param permutation,n.permutations,block.vars,seed permutation settings (see [permutationPlan()]); `block.vars` defaults to the model's
#' @param adjust max-statistic adjustment of the z-scores across cells (default TRUE)
#' @param smooth median-filter the shifts and adjusted z-scores over the graph (default TRUE)
#' @param wins winsorizing fraction for the adjustment (default 0.025)
#' @param adjust.method `"perm-fdr"` (default: pooled permutation FDR over cells, see [permutationFDR()]; no per-cell
#'   floor), `"bh"`, `"spatial-fdr"`, `"neighbourhood"` or `"max-t"`, see [adjustCellScores()]
#' @param p.source per-cell p-values for the BH-type adjustments: `"permutation"` (default; calibrated but floored at
#'   1 / (n.permutations + 1), so with many cells and sparse signal BH may call nothing), `"pooled"` (the cell's
#'   standardized z against the pooled null of all cells' permuted z-scores, see [pooledPermutationP()]: no floor,
#'   calibrated with enough permutations; about 2,000 keep the far tail within a factor of 2, 10,000 make it exact),
#'   `"z"` (the normal tail of the permutation-standardized z-score: no floor, but anti-conservative in the tail by
#'   orders of magnitude) or `"parametric"` (F statistic; not calibrated for this statistic)
#' @param z.grid thresholds of the permutation-standardized z at which the pooled null exceedance counts are
#'   accumulated for `"perm-fdr"` and the pooled p-values (default -6 to 25 in steps of 0.02)
#' @param groups optional factor of cell groups (named by cell): the kernel then also accumulates the sum of the
#'   permuted z-scores per group, from which `group.test` gives the exact group-level permutation test (the mean z of
#'   the group's cells against its permutation distribution; one-sided) with BH across groups; the permutation
#'   distribution of the group means is stored as `attr(group.test, "perm.means")`
#' @param return.perm return the permutations x cells matrix of the permuted statistics as `stats.perm` (default
#'   FALSE; `n.permutations` x cells doubles, NaN where a cell is not tested under a relabeling)
#' @param post post-processing on the cell graph: `"none"` (default), `"averaged"` (the shift statistic averaged over
#'   each cell and its graph neighbours, under every permutation alike, before the pooled FDR: an exact test of the
#'   one-hop-larger unit with calibrated q-values, [averageOverNeighbours()]; needs `adjust.method = "perm-fdr"` and
#'   keeps the permutation matrix in memory) or `"inherit"` (every cell takes the smallest q among itself and its
#'   graph neighbours after the test, [inheritNeighbourQ()])
#' @param nns.graph the 1-hop graph neighbourhoods used by `post` (list of 0-based indices per cell, the cell
#'   itself included; default: `nns.per.cell`, which is the 1-hop graph unless a kNN unit was passed)
#' @details Label permutations (the block scheme) are applied as global relabelings of the samples: in a neighbourhood
#'   where some samples are absent, a present sample takes the design row of the sample it is mapped to, so that one
#'   permutation means the same relabeling in every cell and statistics pooled over cells (the group test, the
#'   pooled permutation FDR, the max statistic) are compared with a coherent null. A cell is not tested under a
#'   relabeling that leaves fewer than `min.samp.per.level` present samples at a compared level (`n.perm` counts the
#'   relabelings under which each cell was tested). Residual permutations (Freedman-Lane) act on the present samples'
#'   residuals and are induced within each cell's sample set.
#' @param log.vectors log10(1e3 x + 1) transform of the mean profiles (default TRUE)
#' @param n.cores threads for the per-cell loop (C++)
#' @param verbose progress
#' @param robust `"none"`, `"huber"` or `"winsor"`: down-weight samples with outlying residual distances (re-estimated under every relabeling)
#' @param robust.k robust tuning constant
#' @return list: `stat` (shift F per cell), `p.value` (exact permutation p), `p.pooled` (pooled-null p, no floor), `p.param`,
#'   `z.pooled` (the tail z of `p.pooled`, `qnorm(1 - p)`), `z.std` (the permutation-standardized
#'   statistic: observed minus permutation mean over permutation sd; the input of the max-t adjustment and of the
#'   pooled null), `z.adj`, `shifts` (bias-corrected shift
#'   estimate), `shifts.smoothed`, `n.samples`, `settings`
#' @export
clusterFreeExpressionShifts <- function(cm, sample.per.cell, nns.per.cell, design, meta, dist = c("cor", "cosine", "js"), min.n.obs.per.samp = 3,
                                        min.samp.per.level = 2, permutation = "auto", n.permutations = 199, block.vars = NULL, seed = 1,
                                        adjust = TRUE, smooth = TRUE, wins = 0.025, log.vectors = TRUE, n.cores = 1, verbose = FALSE,
                                        robust = c("none", "huber", "winsor"), robust.k = 1.345,
                                        adjust.method = c("perm-fdr", "bh", "spatial-fdr", "neighbourhood", "max-t"), p.source = c("permutation", "pooled", "z", "parametric"),
                                        z.grid = seq(-6, 25, by = 0.02), groups = NULL, return.perm = FALSE, post = c("none", "averaged", "inherit"), nns.graph = NULL) {
  dist <- match.arg(dist); robust <- match.arg(robust); adjust.method <- match.arg(adjust.method); p.source <- match.arg(p.source); post <- match.arg(post)
  if (post == "averaged" && adjust.method != "perm-fdr") { if (verbose) message("post = 'averaged' needs adjust.method = 'perm-fdr': post = 'none'"); post <- "none" }
  if (post == "averaged" && !adjust) post <- "none"
  keep.perm <- return.perm || post == "averaged"
  block.vars <- block.vars %||% design$block.vars   # the model's permutation strata unless overridden
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
  cell.names <- names(nns.per.cell) %||% colnames(cm)
  grp.code <- integer(0); grp.levels <- NULL
  if (!is.null(groups)) { gf <- as.factor(groups[cell.names]); grp.levels <- levels(gf); grp.code <- as.integer(gf) - 1L; grp.code[is.na(grp.code)] <- -1L }
  if (verbose) message(sprintf("Testing %d cells (%d samples) with %d permutations (%s)...", m, n, B, gplan$scheme))
  kr <- cluster_free_shift_stream(cm, as.integer(spc), nn.list, FALSE, as.integer(min.n.obs.per.samp), dist, log.vectors,
                                  pairs - 1L, n, des$F, cF, as.integer(level.code), as.integer(min.samp.per.level), as.integer(gplan$strata), as.integer(gplan$in.set),
                                  P.global, gplan$scheme == "freedman-lane", TRUE, as.integer(n.cores), c(none = 0L, huber = 1L, winsor = 2L)[[robust]], robust.k,
                                  as.numeric(z.grid), as.integer(grp.code), keep.perm)
  Fobs <- as.numeric(kr$stat); pval <- as.numeric(kr$p); shifts <- as.numeric(kr$shift); n.samp <- as.integer(kr$n); n.perm <- as.integer(kr$n_perm); z <- as.numeric(kr$z)
  mx <- as.numeric(kr$max); mn <- as.numeric(kr$min)
  valid <- which(is.finite(z))
  names(Fobs) <- names(pval) <- names(shifts) <- names(z) <- names(n.samp) <- names(n.perm) <- cell.names
  z.adj <- NULL; p.adj <- NULL; index.cells <- NULL; fdr.curve <- NULL
  perm.hist <- list(grid = as.numeric(kr$grid), hist = kr$hist)
  q <- ncol(des$F)
  p.param <- stats::pf(Fobs, 1, pmax(n.samp - q, 1), lower.tail = FALSE); p.param[!is.finite(Fobs)] <- NA_real_
  names(p.param) <- cell.names
  p.pooled <- pooledPermutationP(z, perm.hist, alternative = "greater"); names(p.pooled) <- cell.names   # no per-cell floor
  # the reported z-score (z.pooled) is the tail z of the pooled p (one-sided: z = qnorm(1 - p)), calibrated and not floored;
  # the permutation-standardized statistic (observed minus permutation mean over permutation sd) is kept as z.std
  z.std <- z
  z <- stats::qnorm(pmin(pmax(p.pooled, 1e-300), 1 - 1e-16), lower.tail = FALSE); z[!is.finite(p.pooled)] <- NA_real_; names(z) <- cell.names
  stats.perm <- NULL
  if (keep.perm) { stats.perm <- kr$perm; colnames(stats.perm) <- cell.names }
  nn.graph <- if (is.null(nns.graph)) nn.list else { g <- lapply(nns.graph[cell.names], as.integer); if (!is.null(nn.map)) g <- lapply(g, function(v) { w <- nn.map[v + 1L]; w[!is.na(w)] }); g }
  if (post == "averaged") {
    # the one-hop-larger unit: the shift statistic averaged over the cell and its graph neighbours, observed and permuted
    An <- neighbourAveragingList(nn.graph, cell.names)
    s0 <- Fobs; s0[!is.finite(s0)] <- 0; P0 <- stats.perm; P0[!is.finite(P0)] <- 0
    Fobs <- stats::setNames(as.numeric(An %*% s0), cell.names); stats.perm <- as.matrix(Matrix::t(An %*% Matrix::t(P0))); colnames(stats.perm) <- cell.names
    tested <- is.finite(z); Fobs[!tested] <- NA_real_
    e0 <- shifts; e0[!is.finite(e0)] <- 0; shifts <- stats::setNames(as.numeric(An %*% e0), cell.names); shifts[!tested] <- NA_real_
    pval <- (colSums(stats.perm >= rep(Fobs, each = nrow(stats.perm)), na.rm = TRUE) + 1) / (nrow(stats.perm) + 1); pval[!tested] <- NA_real_; names(pval) <- cell.names
    mu <- colMeans(stats.perm); sdv <- apply(stats.perm, 2, stats::sd); z.std.av <- (Fobs - mu) / sdv; z.std.av[!tested | !is.finite(z.std.av)] <- NA_real_
    p.pooled <- pooledPermutationP(Fobs, stats.perm, alternative = "greater"); p.pooled[!tested] <- NA_real_; names(p.pooled) <- cell.names
    z <- stats::qnorm(pmin(pmax(p.pooled, 1e-300), 1 - 1e-16), lower.tail = FALSE); z[!is.finite(p.pooled)] <- NA_real_; names(z) <- cell.names
    p.param[] <- NA_real_
  }
  if (adjust && length(valid)) {
    if (adjust.method == "max-t") {
      z0 <- z.std; z0[!is.finite(z0)] <- 0
      z.adj <- adjustedZScoresMaxStat(z_obs = z0, alt = 1L, max_vals_in = mx, min_vals_in = mn, wins = wins, smooth = smooth,
                                      nn_ids = nn.list, non_zero_ids = as.integer(valid - 1L))
      z.adj[!is.finite(z.std)] <- NA_real_
    } else if (adjust.method == "perm-fdr") {
      ad <- if (post == "averaged") adjustCellScores(Fobs, NULL, method = "perm-fdr", alternative = "greater", perm = stats.perm) else adjustCellScores(z.std, NULL, method = "perm-fdr", alternative = "greater", perm = perm.hist)
      z.adj <- ad$z.adj; p.adj <- ad$p.adj; names(p.adj) <- cell.names; fdr.curve <- ad$fdr.curve
    } else {
      p.in <- switch(p.source, parametric = p.param, permutation = pval, pooled = p.pooled,
                     z = { pz <- stats::pnorm(z.std, lower.tail = FALSE); pz[!is.finite(z.std)] <- NA_real_; names(pz) <- cell.names; pz })   # normal tail of the permutation-standardized z: no floor
      z.in <- stats::qnorm(pmin(pmax(p.in, 1e-300), 1 - 1e-16), lower.tail = FALSE); z.in[!is.finite(p.in)] <- NA_real_
      ad <- adjustCellScores(z.in, p.in, method = adjust.method, nn.list = nn.list, alternative = "greater", seed = seed)
      z.adj <- ad$z.adj; p.adj <- ad$p.adj; names(p.adj) <- cell.names
      if (!is.null(ad$index)) index.cells <- cell.names[ad$index]
    }
    names(z.adj) <- cell.names
    if (post == "inherit" && !is.null(p.adj)) {
      ih <- inheritNeighbourQList(p.adj, shifts, nn.graph, cell.names)
      p.adj.cell <- p.adj; p.adj <- ih$q
      z.adj <- pmax(stats::qnorm(pmin(pmax(p.adj, 1e-300), 1 - 1e-16), lower.tail = FALSE), 0); z.adj[!is.finite(p.adj)] <- NA_real_; names(z.adj) <- cell.names
    }
  }
  group.test <- NULL
  if (length(grp.levels)) {
    gn <- as.integer(kr$group_n); gobs <- as.numeric(kr$group_obs)
    pm <- kr$group_sum / kr$group_count; pm[kr$group_count == 0] <- NA_real_   # permutations x groups: the mean z of the cells tested under each relabeling
    mean.obs <- ifelse(gn > 0, gobs / pmax(gn, 1), NA_real_)
    p.g <- vapply(seq_along(grp.levels), function(g) { ok <- is.finite(pm[, g]); if (gn[g] > 0 && any(ok)) (sum(pm[ok, g] >= mean.obs[g] - 1e-12) + 1) / (sum(ok) + 1) else NA_real_ }, numeric(1))
    group.test <- data.frame(group = grp.levels, n.tested = gn, mean.z = mean.obs, p.perm = p.g, padj = stats::p.adjust(p.g, "BH"), stringsAsFactors = FALSE)
    attr(group.test, "perm.means") <- `colnames<-`(pm, grp.levels)
  }
  shifts.smoothed <- NULL
  if (smooth) {
    s0 <- shifts; s0[!is.finite(s0)] <- 0
    shifts.smoothed <- applyMedianFilterES(s0, nn_ids = nn.list, non_zero_ids = as.integer(which(is.finite(shifts)) - 1L))
    shifts.smoothed[!is.finite(shifts)] <- NA_real_
    names(shifts.smoothed) <- cell.names
  }
  if (!return.perm) stats.perm <- NULL
  if (post == "averaged") z.std <- z.std.av
  out.post <- list(post = post); if (post == "inherit" && exists("p.adj.cell", inherits = FALSE)) out.post$p.adj.cell <- p.adj.cell
  list(stat = Fobs, p.value = pval, p.param = p.param, p.pooled = p.pooled, z.pooled = z, z.std = z.std, z.adj = z.adj, p.adj = p.adj, index.cells = index.cells,
       fdr.curve = fdr.curve, perm.hist = perm.hist, group.test = group.test, stats.perm = stats.perm, post = post, p.adj.cell = out.post$p.adj.cell,
       shifts = shifts, shifts.smoothed = shifts.smoothed, n.samples = n.samp, n.perm = n.perm,
       settings = list(dist = dist, permutation = gplan$scheme, block.vars = block.vars, n.permutations = B, seed = seed, min.n.obs.per.samp = min.n.obs.per.samp, robust = robust,
                       adjust = adjust, adjust.method = adjust.method, p.source = p.source, smooth = smooth, n.cells.tested = length(valid)))
}
