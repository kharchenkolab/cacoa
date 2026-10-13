## Neighbourhood units for the cluster-free maps: kNN neighbourhoods on the cell graph, the effective size of the heat
## kernel, and the two post-processing steps (neighbour-averaged statistic, neighbourhood inheritance).
## Benchmark behind the defaults: misc/benchmark/README.md ("unit x post-processing grid").

#' Default neighbourhood size for the kNN density unit
#'
#' The power of the neighbourhood test is set by the cells each sample contributes to a neighbourhood (about
#' `k / n.samples`), its resolution by how far the neighbourhood reaches: the default takes the larger of 100 cells
#' and 10 cells per sample, capped at 1 % of the cells (never below 20).
#' @param n.cells number of cells in the graph
#' @param n.samples number of samples
#' @param per.sample cells per sample (default 10)
#' @param min.k lower bound (default 100)
#' @param max.frac upper bound as a fraction of the cells (default 0.01)
#' @return integer neighbourhood size
#' @export
defaultNeighbourhoodSize <- function(n.cells, n.samples, per.sample = 10, min.k = 100, max.frac = 0.01) {
  k <- max(min.k, per.sample * n.samples)
  as.integer(max(20, min(k, floor(max.frac * n.cells))))
}

#' kNN neighbourhoods on the cell graph
#'
#' Every cell is the centre of its own neighbourhood of `k` cells (itself included): its graph neighbours ranked by
#' edge weight, then, if those are fewer than `k - 1`, two-hop neighbours ranked by the number of shared neighbours.
#' @param graph igraph cell graph (edge attribute `weight` used for ranking when present)
#' @param k neighbourhood size including the cell itself
#' @param cells cells to use (default: all vertices)
#' @param chunk rows per block of the two-hop product (default 4000)
#' @return a sparse `cells x cells` 0/1 matrix, row i marking the neighbourhood of cell i
#' @export
knnNeighbourhoodMatrix <- function(graph, k, cells = NULL, chunk = 4000L) {
  cells <- cells %||% igraph::V(graph)$name; n <- length(cells); k <- as.integer(k)
  if (k < 2) stop("k must be at least 2")
  has.w <- !is.null(igraph::edge_attr(graph, "weight"))
  W <- igraph::as_adjacency_matrix(graph, attr = if (has.w) "weight" else NULL, sparse = TRUE)[cells, cells, drop = FALSE]
  A <- W; A@x[] <- 1; Matrix::diag(A) <- 1; A <- methods::as(A, "CsparseMatrix")
  out.i <- vector("list", 0); out.j <- vector("list", 0)
  for (s in seq(1, n, by = chunk)) {
    rows <- s:min(n, s + chunk - 1L)
    W1 <- methods::as(W[rows, , drop = FALSE], "TsparseMatrix"); i1 <- W1@i + 1L; j1 <- W1@j + 1L; x1 <- W1@x
    o <- order(i1, -x1); i1 <- i1[o]; j1 <- j1[o]; r1 <- stats::ave(seq_along(i1), i1, FUN = seq_along)
    keep1 <- r1 <= k - 1L
    sel.i <- c(seq_along(rows), i1[keep1]); sel.j <- c(rows, j1[keep1])
    need <- pmax(0L, k - 1L - as.integer(table(factor(i1[keep1], levels = seq_along(rows)))))
    if (any(need > 0)) {
      A1 <- A[rows, , drop = FALSE]; A2 <- A1 %*% A
      M <- Matrix::drop0(A2 - A2 * A1); M <- methods::as(M, "TsparseMatrix")      # two-hop candidates scored by shared neighbours
      i2 <- M@i + 1L; j2 <- M@j + 1L; x2 <- M@x; o <- order(i2, -x2); i2 <- i2[o]; j2 <- j2[o]
      r2 <- stats::ave(seq_along(i2), i2, FUN = seq_along); keep2 <- r2 <= need[i2]
      sel.i <- c(sel.i, i2[keep2]); sel.j <- c(sel.j, j2[keep2])
    }
    out.i[[length(out.i) + 1]] <- rows[sel.i]; out.j[[length(out.j) + 1]] <- sel.j
  }
  Matrix::sparseMatrix(i = unlist(out.i), j = unlist(out.j), x = 1, dims = c(n, n), dimnames = list(cells, cells))
}

#' Effective neighbourhood size of the heat kernel
#'
#' For the kernel K used by the graph density (`sccore::heatFilter`), the effective number of cells a cell averages
#' over is `r_i^2 / sum_k K_ik^2` with `r_i = sum_k K_ik` (the size of an equally weighted neighbourhood with the same
#' variance). Both sums are estimated by a Hutchinson probe: `sum_k K_ik^2 = E[(K v)_i^2]` for random `v` of +-1 and
#' `r_i = (K 1)_i`, which costs one smoothing of `n.vec + 1` vectors.
#' @param graph igraph cell graph
#' @param beta kernel strength (as in `estimateCellDensity`)
#' @param n.vec number of probe vectors (default 24)
#' @param m Chebyshev order (default 50)
#' @param n.cores cores for the smoothing
#' @param seed seed of the probe
#' @return list with `size` (per cell), `row.sum` (per cell) and `beta`
#' @export
kernelEffectiveSize <- function(graph, beta, n.vec = 24, m = 50, n.cores = 1, seed = 7) {
  cells <- igraph::V(graph)$name
  R <- withSeed(seed, cbind(1, matrix(sample(c(-1, 1), length(cells) * n.vec, replace = TRUE), length(cells), n.vec)))
  rownames(R) <- cells
  KR <- as.matrix(sccore::smoothSignalOnGraph(R, filter = function(...) sccore::heatFilter(..., beta = beta), graph = graph, m = m, n.cores = n.cores, progress = FALSE))
  r <- KR[, 1]
  list(size = stats::setNames(r^2 / rowMeans(KR[, -1, drop = FALSE]^2), cells), row.sum = stats::setNames(r, cells), beta = beta)
}

#' Kernel strength for a target effective neighbourhood size
#'
#' Bisection on `log(beta)` so that the median effective size (`kernelEffectiveSize`) matches `target.size`. The
#' same `beta` gives very different widths on different graphs (a 20-fold spread between datasets in the benchmark),
#' so the width is better stated in cells.
#' @inheritParams kernelEffectiveSize
#' @param target.size target median effective size in cells
#' @param range search interval for beta (default c(0.5, 500))
#' @param tol relative tolerance on the size (default 0.1)
#' @param max.iter bisection steps (default 12)
#' @return list with `beta`, `size` (per cell at that beta) and the search `trace`
#' @export
solveKernelBeta <- function(graph, target.size, range = c(0.5, 500), tol = 0.1, max.iter = 12, n.vec = 24, m = 50, n.cores = 1, verbose = FALSE) {
  lo <- log(range[1]); hi <- log(range[2]); trace <- NULL; best <- NULL
  f <- function(lb) { e <- kernelEffectiveSize(graph, exp(lb), n.vec = n.vec, m = m, n.cores = n.cores); e$median <- stats::median(e$size); e }
  e.lo <- f(lo); e.hi <- f(hi); trace <- data.frame(beta = c(exp(lo), exp(hi)), median.size = c(e.lo$median, e.hi$median))
  if (target.size <= e.lo$median) return(list(beta = exp(lo), size = e.lo$size, trace = trace))
  if (target.size >= e.hi$median) return(list(beta = exp(hi), size = e.hi$size, trace = trace))
  for (it in seq_len(max.iter)) {
    mid <- (lo + hi) / 2; e <- f(mid); trace <- rbind(trace, data.frame(beta = exp(mid), median.size = e$median))
    if (verbose) message(sprintf("beta %.3g: median effective size %.0f", exp(mid), e$median))
    best <- list(beta = exp(mid), size = e$size)
    if (abs(log(e$median / target.size)) <= log1p(tol)) break
    if (e$median < target.size) lo <- mid else hi <- mid
  }
  best$trace <- trace; best
}

#' Row-normalized self-plus-neighbours operator of a cell graph
#' @param graph igraph cell graph
#' @param cells cells (rows and columns)
#' @return sparse matrix averaging over a cell and its graph neighbours
#' @keywords internal
neighbourAveragingMatrix <- function(graph, cells) {
  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)[cells, cells, drop = FALSE]; A@x[] <- 1; Matrix::diag(A) <- 1
  Matrix::Diagonal(x = 1 / Matrix::rowSums(A)) %*% A
}

#' Neighbour-averaged statistics
#'
#' Averages the observed and the permuted per-cell statistics over each cell and its graph neighbours, under every
#' permutation alike, so the result is an exact test of the one-hop-larger unit. Non-finite statistics count as 0.
#' @param stat named vector of observed statistics
#' @param stats.perm permutations x cells matrix of permuted statistics (columns in the order of `stat`)
#' @param graph igraph cell graph
#' @param effect optional per-cell effect to average alongside
#' @return list with `stat`, `stats.perm` and `effect` (averaged)
#' @export
averageOverNeighbours <- function(stat, stats.perm, graph, effect = NULL) {
  cells <- names(stat); An <- neighbourAveragingMatrix(graph, cells)
  s <- stat; s[!is.finite(s)] <- 0; P <- stats.perm; P[!is.finite(P)] <- 0
  out <- list(stat = stats::setNames(as.numeric(An %*% s), cells), stats.perm = as.matrix(Matrix::t(An %*% Matrix::t(P))))
  colnames(out$stats.perm) <- cells
  if (!is.null(effect)) { e <- effect[cells]; e[!is.finite(e)] <- 0; out$effect <- stats::setNames(as.numeric(An %*% e), cells) }
  out
}

#' Neighbourhood inheritance of q-values
#'
#' Every cell takes the smallest q among itself and its graph neighbours (and that cell's effect), as Milo lets a
#' cell inherit the best spatial FDR of the neighbourhoods containing it. A selection step, not a test: the realized
#' false-call rate at a nominal threshold is about the nominal rate against a localized truth, where the averaged
#' statistic stays well below it (benchmark). Calibrated under the null (nothing to select).
#' @param q named vector of q-values
#' @param effect named vector of effects
#' @param graph igraph cell graph
#' @return list with `q` and `effect`
#' @export
inheritNeighbourQ <- function(q, effect, graph) {
  cells <- names(q)
  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)[cells, cells, drop = FALSE]; A@x[] <- 1; Matrix::diag(A) <- 1; A <- methods::as(A, "TsparseMatrix")
  qq <- q; qq[!is.finite(qq)] <- Inf
  v <- qq[A@j + 1L]; o <- order(A@i, v); i <- A@i[o] + 1L; j <- A@j[o] + 1L; first <- !duplicated(i)
  out <- rep(NA_real_, length(cells)); out[i[first]] <- qq[j[first]]; out[!is.finite(out)] <- NA_real_
  e2 <- effect[cells]; e2[i[first]] <- effect[cells][j[first]]
  names(out) <- names(e2) <- cells
  list(q = out, effect = e2)
}

#' Row-normalized averaging operator from a list of 0-based neighbour indices (cell itself included)
#' @param nn.list list of integer vectors (0-based indices into `cells`)
#' @param cells cell names
#' @keywords internal
neighbourAveragingList <- function(nn.list, cells) {
  i <- rep(seq_along(cells), lengths(nn.list)); j <- unlist(nn.list) + 1L
  A <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(length(cells), length(cells)), dimnames = list(cells, cells))
  A@x[] <- 1; Matrix::diag(A) <- 1
  Matrix::Diagonal(x = 1 / Matrix::rowSums(A)) %*% A
}

#' Neighbourhood inheritance of q-values from a list of 0-based neighbour indices
#' @inheritParams neighbourAveragingList
#' @param q named vector of q-values
#' @param effect named vector of effects
#' @keywords internal
inheritNeighbourQList <- function(q, effect, nn.list, cells) {
  qq <- q[cells]; qq[!is.finite(qq)] <- Inf; ee <- effect[cells]
  out <- rep(NA_real_, length(cells)); e2 <- ee
  for (i in seq_along(cells)) { nb <- unique(c(i, nn.list[[i]] + 1L)); j <- nb[which.min(qq[nb])]; if (is.finite(qq[j])) { out[i] <- qq[j]; e2[i] <- ee[j] } }
  names(out) <- names(e2) <- cells
  list(q = out, effect = e2)
}
