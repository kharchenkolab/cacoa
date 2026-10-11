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
#'   of Milo), or `"neighbourhood"` (a covering set of index cells is tested, weighted BH over those neighbourhoods
#'   with weights inverse to their overlap, and every cell inherits the smallest adjusted p-value of the index
#'   neighbourhoods it belongs to; the neighbourhood-level testing of Milo without its negative-binomial model)
#' @param nn.list list of 0-based neighbour indices per cell (needed for `"spatial-fdr"` and `"neighbourhood"`)
#' @param weights optional per-cell weights for `"spatial-fdr"`
#' @param alternative `"two.sided"` (default) or `"greater"`: how `p` relates to `z` and how `z.adj` is derived
#' @param seed seed of the index-cell sampling
#' @return list: `p.adj` (per cell), `z.adj` (signed z-scores of the adjusted p-values; 0 when `p.adj` is 1),
#'   `index` (index cells of the neighbourhood method), `method`
#' @export
adjustCellScores <- function(z, p = NULL, method = c("bh", "spatial-fdr", "neighbourhood"), nn.list = NULL, weights = NULL,
                             alternative = c("two.sided", "greater"), seed = 1) {
  method <- match.arg(method); alternative <- match.arg(alternative)
  if (is.null(p)) p <- if (alternative == "two.sided") 2 * stats::pnorm(-abs(z)) else stats::pnorm(z, lower.tail = FALSE)
  ok <- is.finite(z) & is.finite(p)
  p.adj <- rep(NA_real_, length(z)); names(p.adj) <- names(z); idx <- NULL
  if (!any(ok)) return(list(p.adj = p.adj, z.adj = p.adj, index = idx, method = method))
  if (method == "bh") {
    p.adj[ok] <- stats::p.adjust(p[ok], "BH")
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
  list(p.adj = p.adj, z.adj = z.adj, index = idx, method = method)
}
