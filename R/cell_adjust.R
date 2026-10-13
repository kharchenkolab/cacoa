## Multiplicity adjustment of per-cell statistics for the cluster-free tests (density, expression shifts).
## Besides the max-statistic adjustment (familywise, over all cells) the package offers false-discovery-rate
## alternatives: BH over the tested cells, a spatially weighted BH (Milo's spatial FDR, cells in dense regions carry
## less weight because their neighbourhoods overlap), and a neighbourhood-level test in which a covering set of index
## cells is tested and every cell inherits the smallest adjusted p-value of the neighbourhoods it belongs to.

# Weighted Benjamini-Hochberg (as in miloR::graphSpatialFDR)
weightedBH <- function(p, w) {
  o <- order(p); ps <- p[o]; ws <- w[o]
  adj <- rev(cummin(rev(sum(ws) * ps / cumsum(ws))))
  out <- numeric(length(p)); out[o] <- pmin(adj, 1); out
}

# A covering set of index cells: random order, a cell becomes an index cell unless it is already inside an index
# neighbourhood (itself or a neighbour of an earlier index cell)
coveringIndexCells <- function(nn.list, ok, seed = 1) {
  covered <- !ok; idx <- integer(0)
  ord <- withSeed(seed, sample(which(ok)))
  for (i in ord) {
    if (covered[i]) next
    idx <- c(idx, i); covered[i] <- TRUE; covered[nn.list[[i]] + 1L] <- TRUE
  }
  idx
}

# 0-based neighbour lists of a cell graph, restricted and ordered to `cells`
graphNeighbourList <- function(graph, cells) {
  adj <- igraph::as_adj(graph)
  cells <- intersect(cells, rownames(adj))
  adj <- methods::as(adj[cells, cells, drop = FALSE], "TsparseMatrix")
  nn <- split(adj@j, factor(adj@i + 1L, levels = seq_along(cells)))
  nn <- lapply(nn, as.integer); names(nn) <- cells
  nn
}

#' Adjust per-cell test statistics for multiplicity across cells
#'
#' @param z named numeric vector of signed per-cell z-scores (`NA` for cells that were not tested)
#' @param p per-cell p-values (default: two-sided normal p-values from `z`; pass one-sided p-values with
#'   `alternative = "greater"` for a one-sided statistic such as the expression shift)
#' @param method `"bh"` (Benjamini-Hochberg over the tested cells), `"spatial-fdr"` (weighted BH with weights
#'   `1 / (1 + degree)`, so that cells in dense regions, whose neighbourhoods overlap most, count less; the spatial FDR
#'   of Milo), `"neighbourhood"` (a covering set of index cells is tested, weighted BH over those neighbourhoods
#'   with weights inverse to their overlap, and every cell inherits the smallest adjusted p-value of the index
#'   neighbourhoods it belongs to; the neighbourhood-level testing of Milo without its negative-binomial model) or
#'   `"perm-fdr"` (pooled permutation FDR, see [permutationFDR()]: the expected number of cells exceeding a threshold
#'   under the permutations, pooled over cells, against the observed number; no per-cell p-value floor; needs `perm`)
#' @param nn.list list of 0-based neighbour indices per cell (needed for `"spatial-fdr"` and `"neighbourhood"`)
#' @param weights optional per-cell weights for `"spatial-fdr"`
#' @param alternative `"two.sided"` (default) or `"greater"`: how `p` relates to `z` and how `z.adj` is derived
#' @param seed seed of the index-cell sampling
#' @param perm for `"perm-fdr"`: the permuted statistics, either a permutations x cells matrix on the scale of `z`, or a
#'   list `(grid, hist)` of pooled exceedance histograms (see [permutationFDR()])
#' @return list: `p.adj` (per cell), `z.adj` (signed z-scores of the adjusted p-values; 0 when `p.adj` is 1),
#'   `index` (index cells of the neighbourhood method), `method`, and for `"perm-fdr"` the `fdr.curve`
#' @export
adjustCellScores <- function(z, p = NULL, method = c("bh", "spatial-fdr", "neighbourhood", "perm-fdr"), nn.list = NULL, weights = NULL,
                             alternative = c("two.sided", "greater"), seed = 1, perm = NULL) {
  method <- match.arg(method); alternative <- match.arg(alternative)
  if (is.null(p)) p <- if (alternative == "two.sided") 2 * stats::pnorm(-abs(z)) else stats::pnorm(z, lower.tail = FALSE)
  ok <- is.finite(z) & is.finite(p)
  p.adj <- rep(NA_real_, length(z)); names(p.adj) <- names(z); idx <- NULL; curve <- NULL
  if (!any(ok)) return(list(p.adj = p.adj, z.adj = p.adj, index = idx, method = method))
  if (method == "bh") {
    p.adj[ok] <- stats::p.adjust(p[ok], "BH")
  } else if (method == "perm-fdr") {
    if (is.null(perm)) stop("perm (the permuted statistics) is needed for method 'perm-fdr'")
    pf <- permutationFDR(z, perm, alternative = alternative)
    p.adj <- pf$q; names(p.adj) <- names(z); curve <- pf$curve
  } else {
    if (is.null(nn.list) || length(nn.list) != length(z)) stop("nn.list (one entry per cell) is needed for method '", method, "'")
    deg <- vapply(nn.list, length, integer(1))
    if (method == "spatial-fdr") {
      w <- weights %||% (1 / (1 + deg))
      p.adj[ok] <- weightedBH(p[ok], w[ok])
    } else {
      idx <- coveringIndexCells(nn.list, ok, seed = seed)
      nb <- lapply(idx, function(i) unique(c(i, nn.list[[i]] + 1L)))
      overlap <- tabulate(unlist(nb), nbins = length(z))                     # index neighbourhoods containing each cell
      w.n <- vapply(nb, function(m) 1 / mean(overlap[m]), numeric(1))         # overlapping neighbourhoods weigh less
      pa <- weightedBH(p[idx], w.n)
      best <- rep(NA_real_, length(z))
      for (j in seq_along(idx)) { m <- nb[[j]]; best[m] <- pmin(best[m], pa[j], na.rm = TRUE) }
      p.adj[ok] <- best[ok]
    }
  }
  pa <- pmin(pmax(p.adj, 1e-300), 1)
  z.adj <- if (alternative == "two.sided") sign(z) * stats::qnorm(pa / 2, lower.tail = FALSE) else pmax(stats::qnorm(pa, lower.tail = FALSE), 0)   # one-sided: 0 at p >= 0.5
  z.adj[!is.finite(z.adj)] <- NA_real_
  list(p.adj = p.adj, z.adj = z.adj, index = idx, method = method, fdr.curve = curve)
}

#' Pooled permutation false discovery rate for per-cell statistics
#'
#' The SAM-type estimate (Tusher et al. 2001; Storey & Tibshirani 2003): for a threshold t the expected number of
#' cells with a null statistic at or above t is the average, over the permutations, of the number of cells whose
#' permuted statistic is at or above t, pooled over all cells; the FDR at t is that count divided by the observed
#' number of cells at or above t, and a cell's q-value is the smallest FDR over the thresholds at or below its own
#' statistic. Pooling over cells removes the per-cell permutation floor (a cell can only be as significant as
#' 1 / (B + 1) on its own, but the pooled null count resolves far smaller rates), and the permuted maps carry the
#' dependence between cells, so no independence assumption is made. Cells must be exchangeable with their own
#' permuted statistics under the null, which holds for the permutation-standardized z-scores (shift F) and for the
#' studentized t of the density contrast; their null distributions need not be identical across cells.
#'
#' @param stat observed per-cell statistics (`NA` for untested cells)
#' @param perm either a permutations x cells numeric matrix of the permuted statistics (`NA` allowed), or a list with
#'   `grid` (increasing thresholds) and `hist` (a permutations x (length(grid) + 1) count matrix: column k + 1 holds the
#'   number of cells whose permuted statistic fell in `[grid[k], grid[k + 1])`, column 1 those below `grid[1]`, the last
#'   column those at or above the last threshold), as returned by the streaming cluster-free kernel
#' @param alternative `"two.sided"` (thresholds on `abs(stat)`) or `"greater"`
#' @param pi0 proportion of null cells assumed (1, conservative, by default)
#' @param n.grid number of thresholds when `perm` is a matrix (quantiles of the observed statistics)
#' @return list: `q` (per cell), `curve` (data.frame: `threshold`, `n.observed`, `n.null` (expected), `fdr`)
#' @export
permutationFDR <- function(stat, perm, alternative = c("two.sided", "greater"), pi0 = 1, n.grid = 2000) {
  alternative <- match.arg(alternative)
  s <- if (alternative == "two.sided") abs(stat) else stat
  ok <- is.finite(s)
  q <- rep(NA_real_, length(stat)); names(q) <- names(stat)
  if (!any(ok)) return(list(q = q, curve = NULL))
  if (is.list(perm)) {
    grid <- as.numeric(perm$grid); H <- perm$hist
    if (ncol(H) != length(grid) + 1L) stop("perm$hist must have length(perm$grid) + 1 columns")
    if (alternative == "two.sided") stop("histogram input is one-sided (thresholds on the statistic itself); use alternative = 'greater'")
    # expected null count at or above grid[k]: mean over permutations of the counts in columns k + 1 .. end
    cum <- t(apply(H[, ncol(H):1, drop = FALSE], 1, cumsum))[, ncol(H):1, drop = FALSE]   # reverse cumulative sums
    n.null <- colMeans(cum[, -1, drop = FALSE])                                              # columns 2..: at or above grid[k]
  } else {
    P <- if (alternative == "two.sided") abs(perm) else perm
    if (ncol(P) != length(stat)) stop("perm must have one column per cell")
    P <- P[, ok, drop = FALSE]
    so <- sort(s[ok])
    grid <- if (length(so) > n.grid) unique(stats::quantile(so, probs = seq(0, 1, length.out = n.grid), names = FALSE, type = 1)) else unique(so)
    n.null <- numeric(length(grid))
    for (b in seq_len(nrow(P))) {
      v <- P[b, ]; v <- v[is.finite(v)]
      h <- tabulate(findInterval(v, grid) + 1L, nbins = length(grid) + 1L)
      n.null <- n.null + rev(cumsum(rev(h)))[-1]
    }
    n.null <- n.null / nrow(P)
  }
  so <- sort(s[ok]); n.obs <- length(so) - findInterval(grid, so, left.open = TRUE)     # cells with s >= grid[k]
  fdr <- pmin(1, pi0 * n.null / pmax(n.obs, 1)); fdr[n.obs == 0] <- NA_real_
  # q-value: the smallest FDR over the thresholds at or below the cell's statistic (the most favourable rejection set
  # that still contains the cell), which makes q non-increasing in the statistic
  fdr.run <- fdr; fdr.run[is.na(fdr.run)] <- 1
  qmin <- cummin(fdr.run)
  k <- findInterval(s[ok], grid)
  q[ok] <- ifelse(k == 0, 1, qmin[pmax(k, 1)])
  list(q = q, curve = data.frame(threshold = grid, n.observed = n.obs, n.null = n.null, fdr = fdr))
}

#' Per-cell p-values from the pooled null of the permutation-standardized statistics
#'
#' Every cell's statistic is standardized by the mean and standard deviation of its own permutation distribution,
#' and the standardized permuted statistics of all cells are pooled into one null sample, against which each cell's
#' standardized observed statistic is ranked. The resolution is set by the pooled count (cells x permutations) rather
#' than by `1 / (B + 1)`, so the p-values have no per-cell floor. No parametric form is assumed: the normal tail of the
#' standardized statistic and moment-matched gamma tails were both found to be anti-conservative by orders of magnitude
#' in the far tail (`misc/benchmark/README.md`), whereas the pooled null is calibrated as long as the permutations
#' sample the extreme relabelings well. The far tail of the pooled null is formed by few relabelings (those that
#' isolate a deviant sample are extreme for many cells at once), so its accuracy below about p = 1e-3 is governed by
#' the number of permutations, not by the number of cells: about 2,000 permutations keep the false-positive rate within
#' a factor of about 2 of the nominal level down to 1e-4, and 10,000 or more make it exact there.
#'
#' @param stat observed per-cell statistics (`NA` for untested cells); with a histogram `perm` these must be the
#'   per-cell standardized z-scores themselves
#' @param perm either a permutations x cells matrix of the permuted statistics (`NA` allowed; standardized here per
#'   column), or the kernel histogram list (`grid`, `hist`) of already standardized z-scores as returned by the
#'   cluster-free shift kernel (see [permutationFDR()])
#' @param alternative `"two.sided"` (ranks `abs(z)`) or `"greater"`; a histogram input is one-sided
#' @param resolution bin width of the pooled standardized null when `perm` is a matrix (default 0.005; a cell's
#'   p-value counts the pooled values at or above the lower edge of its bin)
#' @return named numeric vector of p-values
#' @export
pooledPermutationP <- function(stat, perm, alternative = c("two.sided", "greater"), resolution = 0.005) {
  alternative <- match.arg(alternative)
  p <- rep(NA_real_, length(stat)); names(p) <- names(stat)
  if (is.list(perm)) {
    if (alternative != "greater") stop("histogram input is one-sided; use alternative = 'greater'")
    grid <- as.numeric(perm$grid); H <- perm$hist
    if (ncol(H) != length(grid) + 1L) stop("perm$hist must have length(perm$grid) + 1 columns")
    counts <- colSums(H); z <- as.numeric(stat)
  } else {
    if (ncol(perm) != length(stat)) stop("perm must have one column per cell")
    n <- colSums(is.finite(perm)); mu <- colMeans(perm, na.rm = TRUE)
    v <- (colSums(perm^2, na.rm = TRUE) - n * mu^2) / pmax(n - 1, 1); sdv <- sqrt(pmax(v, 0)); sdv[n < 3 | sdv <= 0] <- NA_real_
    z <- (as.numeric(stat) - mu) / sdv
    if (alternative == "two.sided") z <- abs(z)
    zf <- z[is.finite(z)]; lo <- if (length(zf)) min(zf) else 0; hi <- if (length(zf)) max(zf) else 0
    grid <- seq(floor(min(lo, -1) / resolution) * resolution, max(hi, 1) + resolution, by = resolution)
    counts <- numeric(length(grid) + 1L)
    for (cols in split(seq_len(ncol(perm)), ceiling(seq_len(ncol(perm)) / 2000))) {   # chunked histogram of the pooled z
      Z <- sweep(sweep(perm[, cols, drop = FALSE], 2, mu[cols]), 2, sdv[cols], "/")
      if (alternative == "two.sided") Z <- abs(Z)
      Z <- Z[is.finite(Z)]
      counts <- counts + tabulate(findInterval(Z, grid) + 1L, nbins = length(grid) + 1L)
    }
  }
  N <- sum(counts)
  if (N == 0) return(p)
  cum <- rev(cumsum(rev(counts)))                       # cum[j]: pooled values in bins j..end (bin 1: below grid[1])
  ok <- is.finite(z)
  k <- findInterval(z[ok], grid) + 1L                   # the bin holding the cell's z
  p[ok] <- (cum[k] + 1) / (N + 1)
  p
}

#' Cauchy combination of dependent p-values (Liu & Xie 2020)
#'
#' The weighted mean of `tan((0.5 - p) * pi)` is standard Cauchy under the global null for any dependence between the
#' p-values; the combined p-value is its upper tail. Useful for summarizing per-cell tests over a region (a cluster, a
#' neighbourhood of the embedding). Very small p-values dominate, so the combination behaves like a minimum-p test
#' that needs no multiplicity correction within the region.
#'
#' @param p p-values (`NA` dropped)
#' @param weights non-negative weights (default equal)
#' @param p.min,p.max p-values are clamped into `[p.min, p.max]` before the transform (permutation p-values of exactly
#'   1 would otherwise cancel any signal)
#' @return the combined p-value (`NA` when no finite p-value is given)
#' @export
cauchyCombination <- function(p, weights = NULL, p.min = 1e-15, p.max = 1 - 1e-3) {
  ok <- is.finite(p)
  if (!any(ok)) return(NA_real_)
  p <- pmin(pmax(p[ok], p.min), p.max)
  w <- if (is.null(weights)) rep(1, length(p)) else weights[ok]
  w <- w / sum(w)
  t.stat <- sum(w * tan((0.5 - p) * pi))
  stats::pcauchy(t.stat, lower.tail = FALSE)
}

#' Exact group-level permutation test of per-cell statistics
#'
#' For every group, the mean of the per-cell statistics over the group's tested cells is compared with the same mean
#' under every permutation of the sample labels (the permutations are shared by all cells, so the dependence between
#' neighbouring cells is carried along exactly). Two-sided for signed statistics (the density contrast), one-sided
#' for the shift statistic.
#'
#' @param stat observed per-cell statistics named by cell (`NA` for untested cells)
#' @param stats.perm permutations x cells matrix of the permuted statistics (columns named by cell)
#' @param groups factor of cell groups named by cell
#' @param alternative `"two.sided"` or `"greater"`
#' @param min.cells groups with fewer tested cells get `NA`
#' @return data.frame: `group`, `n.tested`, `mean.stat`, `p.perm`, `padj` (BH across groups)
#' @export
groupPermutationTest <- function(stat, stats.perm, groups, alternative = c("two.sided", "greater"), min.cells = 10) {
  alternative <- match.arg(alternative); groups <- as.factor(groups)
  ok <- names(stat)[is.finite(stat)]; ok <- intersect(ok, intersect(names(groups), colnames(stats.perm)))
  B <- nrow(stats.perm)
  res <- lapply(levels(groups), function(g) {
    idx <- ok[groups[ok] == g]
    if (length(idx) < min.cells) return(data.frame(group = g, n.tested = length(idx), mean.stat = NA_real_, p.perm = NA_real_, stringsAsFactors = FALSE))
    obs <- mean(stat[idx]); perm <- rowMeans(stats.perm[, idx, drop = FALSE], na.rm = TRUE)
    p <- if (alternative == "two.sided") (sum(abs(perm) >= abs(obs) - 1e-12, na.rm = TRUE) + 1) / (B + 1) else (sum(perm >= obs - 1e-12, na.rm = TRUE) + 1) / (B + 1)
    data.frame(group = g, n.tested = length(idx), mean.stat = obs, p.perm = p, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res); res$padj <- stats::p.adjust(res$p.perm, "BH"); res
}

#' Summarize per-cell test results by cell group
#'
#' One p-value per group with BH adjustment across groups, next to the median effect, the median z-score and the
#' fraction of the group's cells significant after the per-cell adjustment. Two group tests are offered: the exact
#' permutation test of the group's mean statistic (see [groupPermutationTest()]; the default whenever the permuted
#' statistics or a precomputed `group.test` are given) and the Cauchy combination of the cells' p-values (see
#' [cauchyCombination()]; valid under dependence but driven by the smallest p-values, so it needs per-cell p-values
#' that are calibrated in the far tail: permutation p-values rather than parametric ones). Either gives a group-level
#' verdict where the per-cell adjustment calls nothing.
#'
#' @param p per-cell p-values (named by cell; untested cells `NA`)
#' @param groups factor of cell groups named by cell
#' @param effect optional per-cell effects (signed; for the direction and the median)
#' @param z optional per-cell z-scores
#' @param p.adj optional per-cell adjusted p-values (for the fraction of significant cells)
#' @param alpha significance level for the fraction of significant cells and the `significant` flag (default 0.05)
#' @param min.cells groups with fewer tested cells are reported with `NA` (default 10)
#' @param stat,stats.perm observed and permuted per-cell statistics for the permutation group test
#' @param group.test a precomputed permutation group test (data.frame with `group` and `p.perm`)
#' @param alternative for the permutation group test
#' @param method which p-value fills `p` and `padj`: `"permutation"` (default when available) or `"cauchy"`
#' @return data.frame: `group`, `n.cells`, `n.tested`, `p`, `padj`, `p.perm`, `p.cauchy`, `median.effect`, `median.z`,
#'   `frac.sig`, `direction` (sign of the median effect: "up", "down" or "none"), `significant`, `method`
#' @export
summarizeCellScoresByGroup <- function(p, groups, effect = NULL, z = NULL, p.adj = NULL, alpha = 0.05, min.cells = 10,
                                       stat = NULL, stats.perm = NULL, group.test = NULL, alternative = c("two.sided", "greater"),
                                       method = c("permutation", "cauchy")) {
  method <- match.arg(method); alternative <- match.arg(alternative)
  groups <- as.factor(groups)
  cells <- names(groups)
  if (is.null(cells) || is.null(names(p))) stop("p and groups must be named by cell")
  pc <- p[cells]; ec <- if (!is.null(effect)) effect[cells] else NULL; zc <- if (!is.null(z)) z[cells] else NULL; ac <- if (!is.null(p.adj)) p.adj[cells] else NULL
  res <- lapply(levels(groups), function(g) {
    i <- which(groups == g); tested <- i[is.finite(pc[i])]
    med.e <- if (!is.null(ec) && length(tested)) stats::median(ec[tested], na.rm = TRUE) else NA_real_
    data.frame(group = g, n.cells = length(i), n.tested = length(tested),
               p.cauchy = if (length(tested) >= min.cells) cauchyCombination(pc[tested]) else NA_real_,
               median.effect = med.e, median.z = if (!is.null(zc) && length(tested)) stats::median(zc[tested], na.rm = TRUE) else NA_real_,
               frac.sig = if (!is.null(ac) && length(tested)) mean(ac[tested] < alpha, na.rm = TRUE) else NA_real_,
               direction = if (!is.finite(med.e) || med.e == 0) "none" else if (med.e > 0) "up" else "down", stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  if (is.null(group.test) && !is.null(stat) && !is.null(stats.perm)) group.test <- groupPermutationTest(stat, stats.perm, groups, alternative = alternative, min.cells = min.cells)
  res$p.perm <- if (!is.null(group.test)) group.test$p.perm[match(res$group, group.test$group)] else NA_real_
  if (method == "permutation" && is.null(group.test)) method <- "cauchy"
  res$p <- if (method == "permutation") res$p.perm else res$p.cauchy
  res$padj <- stats::p.adjust(res$p, "BH")
  res$significant <- is.finite(res$padj) & res$padj < alpha
  res$method <- method
  res[, c("group", "n.cells", "n.tested", "p", "padj", "p.perm", "p.cauchy", "median.effect", "median.z", "frac.sig", "direction", "significant", "method")]
}

#' Group-level permutation test from the accumulated per-group sums of the C++ kernels
#'
#' The kernels accumulate, per permutation and group, the sum and the count of the permuted standardized
#' statistics of the group's tested cells (and the observed sum): the mean per group is compared with its own
#' permutation distribution, one-sided or two-sided, with BH across groups (see [groupPermutationTest()]).
#' @param acc the `pooled` list of a kernel result (`group_sum`, `group_count`, `group_obs`, `group_n`)
#' @param levels group levels in the order of the codes
#' @param alternative `"greater"` or `"two.sided"`
#' @param min.cells groups with fewer tested cells get `NA`
#' @param negate test `-z` (the one-sided "less" alternative)
#' @keywords internal
groupTestFromAccum <- function(acc, levels, alternative = c("greater", "two.sided"), min.cells = 10, negate = FALSE) {
  alternative <- match.arg(alternative)
  if (!length(acc$group_n) || length(acc$group_n) != length(levels)) return(NULL)
  gn <- as.integer(acc$group_n); gobs <- as.numeric(acc$group_obs)
  pm <- acc$group_sum / acc$group_count; pm[acc$group_count == 0] <- NA_real_
  mean.obs <- ifelse(gn > 0, gobs / pmax(gn, 1), NA_real_)
  if (negate) { pm <- -pm; mean.obs <- -mean.obs }   # one-sided "less": the test of -z
  p.g <- vapply(seq_along(levels), function(g) { ok <- is.finite(pm[, g]); if (gn[g] < min.cells || !any(ok)) return(NA_real_)
    if (alternative == "greater") (sum(pm[ok, g] >= mean.obs[g] - 1e-12) + 1) / (sum(ok) + 1) else (sum(abs(pm[ok, g]) >= abs(mean.obs[g]) - 1e-12) + 1) / (sum(ok) + 1) }, numeric(1))
  gt <- data.frame(group = levels, n.tested = gn, mean.z = mean.obs, p.perm = p.g, padj = stats::p.adjust(p.g, "BH"), stringsAsFactors = FALSE)
  attr(gt, "perm.means") <- `colnames<-`(pm, levels); gt
}

#' Fixed grids for the pooled histograms of the standardized statistic
#' @param abs nonnegative grid for `|z|` (two-sided maps) or a signed grid (one-sided)
#' @keywords internal
pooledZGrid <- function(abs = FALSE) if (abs) seq(0, 25, by = 0.01) else seq(-6, 25, by = 0.02)
