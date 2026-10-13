## Neighbourhood units for the cluster-free maps: kNN neighbourhoods on the cell graph, the effective size of the heat
## kernel, and the two post-processing steps (neighbour-averaged statistic, neighbourhood inheritance).
## Benchmark behind the defaults: misc/benchmark/README.md ("unit x post-processing grid").

#' Default neighbourhood size for the kNN density unit
#'
#' The power of the neighbourhood test is set by the cells each sample contributes to a neighbourhood (about
#' `k / n.samples`), its resolution by how far the neighbourhood reaches. The default takes the larger of `min.k`
#' cells and `per.sample` cells per sample, capped at `max.frac` of the cells, but the cap never goes below
#' `min.per.sample` cells per sample (a count test on one cell per sample has no power, whatever the dataset size).
#' @param n.cells number of cells in the graph
#' @param n.samples number of samples
#' @param per.sample cells per sample (default 10)
#' @param min.k lower bound (default 100)
#' @param max.frac upper bound as a fraction of the cells (default 0.01)
#' @param min.per.sample the cap's own floor in cells per sample (default 5)
#' @return integer neighbourhood size
#' @export
defaultNeighbourhoodSize <- function(n.cells, n.samples, per.sample = 10, min.k = 100, max.frac = 0.01, min.per.sample = 5) {
  k.power <- max(min.k, per.sample * n.samples)
  k.cap <- max(floor(max.frac * n.cells), min.per.sample * n.samples, 20)
  as.integer(min(k.power, k.cap))
}

#' kNN neighbourhoods on the cell graph
#'
#' Every cell is the centre of its own neighbourhood of at least `k` cells (itself included): all of its graph
#' neighbours, then, if those are fewer than `k - 1`, two-hop neighbours ranked by the number of shared neighbours.
#' The 1-hop neighbourhood is never cut down to the strongest edges: when expression differs between the compared
#' arms, the strongest edges join cells of the same arm, and a unit built from them hides the composition change
#' it is meant to test (on a dense simulated graph the own-arm share of such a unit was 0.60 against 0.50 for the
#' full neighbourhood, and the planted changes vanished). Cells with more than `k - 1` neighbours therefore get a
#' larger unit; `lengths()` of the lists or `Matrix::rowSums()` of the matrix give the sizes. Computed in C++.
#' @param graph igraph cell graph
#' @param k minimum neighbourhood size including the cell itself
#' @param cells cells to use (default: all vertices)
#' @param n.cores threads
#' @return `knnNeighbourhoodLists()`: a list named by cell of 0-based integer vectors into `cells` (the cell itself
#'   first); `knnNeighbourhoodMatrix()`: the same as a sparse `cells x cells` 0/1 matrix, row i marking the
#'   neighbourhood of cell i
#' @export
knnNeighbourhoodLists <- function(graph, k, cells = NULL, n.cores = 1) {
  cells <- cells %||% igraph::V(graph)$name; k <- as.integer(k)
  if (k < 2) stop("k must be at least 2")
  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)[cells, cells, drop = FALSE]; A@x[] <- 1
  nn <- knn_neighbourhoods(methods::as(A, "dgCMatrix"), k, as.integer(n.cores)); names(nn) <- cells
  nn
}

#' @rdname knnNeighbourhoodLists
#' @export
knnNeighbourhoodMatrix <- function(graph, k, cells = NULL, n.cores = 1) {
  nn <- knnNeighbourhoodLists(graph, k, cells = cells, n.cores = n.cores); cells <- names(nn)
  Matrix::sparseMatrix(i = rep(seq_along(nn), lengths(nn)), j = unlist(nn, use.names = FALSE) + 1L, x = 1, dims = c(length(cells), length(cells)), dimnames = list(cells, cells))
}

#' Per-sample counts of the cells of each neighbourhood
#' @param nn neighbourhood lists (0-based, see [knnNeighbourhoodLists()])
#' @param sample.per.cell factor of samples named by cell, in the order of the cells the lists index
#' @param samples column order (default: the sorted samples present)
#' @return numeric cells x samples matrix
#' @export
neighbourhoodCounts <- function(nn, sample.per.cell, samples = NULL) {
  spc <- as.character(sample.per.cell[names(nn)]); samples <- samples %||% sort(unique(spc)); spc <- factor(spc, levels = samples)
  cnt <- neighbourhood_counts(unname(nn), as.integer(spc), length(samples), FALSE)
  dimnames(cnt) <- list(names(nn), samples); cnt
}

#' 1-hop graph neighbourhoods as 0-based lists (cell itself first)
#' @param graph igraph cell graph
#' @param cells cells (the index space of the lists)
#' @keywords internal
graphNeighbourLists0 <- function(graph, cells) {
  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)[cells, cells, drop = FALSE]; A@x[] <- 1; Matrix::diag(A) <- 1
  A <- methods::as(A, "TsparseMatrix"); o <- order(A@i, A@j != A@i, A@j)
  nn <- split(A@j[o], factor(A@i[o] + 1L, levels = seq_along(cells))); nn <- lapply(nn, as.integer); names(nn) <- cells; nn
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

#' Average a per-cell vector over 0-based neighbour lists (cell itself included in the lists)
#' @param x named numeric vector (non-finite as 0)
#' @param nn neighbour lists (0-based into `names(x)`)
#' @keywords internal
averageVectorOverLists <- function(x, nn) {
  v <- as.numeric(x); v[!is.finite(v)] <- 0
  out <- vapply(nn, function(ids) mean(v[ids + 1L]), numeric(1)); names(out) <- names(x); out
}

#' Diagnostics of the kNN neighbourhood units
#'
#' What a neighbourhood size `k` means on a dataset, without running any test: how many cells each sample
#' contributes to a unit (the power side: a sample with fewer than `min.cells` cells in a unit is dropped from
#' that unit's test), how large the units are, and how pure they are with respect to an annotation (the
#' resolution side: the fraction of a unit's cells that share the centre cell's label).
#' @param nn neighbourhood lists (see [knnNeighbourhoodLists()])
#' @param sample.per.cell factor of samples named by cell
#' @param cell.groups optional factor of cell groups named by cell
#' @param min.cells cells per sample below which a sample is dropped from a unit's test (default 3)
#' @param k the requested minimum size (reported; default: the smallest unit)
#' @return list: `size` (quantiles of the unit size), `per.sample` (data.frame: sample, median cells per unit,
#'   fraction of units with fewer than `min.cells` of the sample's cells), `purity` (data.frame per cell group:
#'   number of cells, median purity, fraction of cells with purity below 0.5; `NULL` without groups), `summary`
#'   (one row: k, median size, median cells per sample per unit, fraction of unit-sample pairs below `min.cells`,
#'   median purity) and `text` (one sentence)
#' @export
neighbourhoodDiagnostics <- function(nn, sample.per.cell, cell.groups = NULL, min.cells = 3, k = NULL) {
  cells <- names(nn); spc <- as.character(sample.per.cell[cells]); samples <- sort(unique(spc))
  cnt <- neighbourhoodCounts(nn, sample.per.cell, samples = samples)
  size <- stats::quantile(lengths(nn), c(0.05, 0.25, 0.5, 0.75, 0.95))
  per.sample <- data.frame(sample = samples, median.cells = apply(cnt, 2, stats::median), frac.units.below = colMeans(cnt < min.cells), stringsAsFactors = FALSE)
  purity <- NULL; med.purity <- NA_real_
  if (!is.null(cell.groups)) {
    cg <- as.character(cell.groups[cells]); pur <- vapply(seq_along(nn), function(i) mean(cg[nn[[i]] + 1L] == cg[i], na.rm = TRUE), numeric(1))
    purity <- do.call(rbind, lapply(split(seq_along(cells), cg), function(idx) data.frame(n.cells = length(idx), median.purity = stats::median(pur[idx]), frac.below.0.5 = mean(pur[idx] < 0.5))))
    purity <- data.frame(group = rownames(purity), purity, row.names = NULL, stringsAsFactors = FALSE); med.purity <- stats::median(pur)
  }
  summary <- data.frame(k = as.integer(k %||% min(lengths(nn))), median.size = unname(size["50%"]), median.cells.per.sample = stats::median(cnt), frac.below.min = mean(cnt < min.cells), median.purity = med.purity)
  text <- sprintf("kNN units of %d cells (median size %d): a sample contributes %d cells to a unit (median), %.0f%% of the unit-sample pairs have fewer than %d cells%s",
                  summary$k, as.integer(summary$median.size), as.integer(summary$median.cells.per.sample), 100 * summary$frac.below.min, min.cells,
                  if (is.na(med.purity)) "" else sprintf("; median purity %.2f against the annotation (lowest: %s)", med.purity, paste(utils::head(purity$group[order(purity$median.purity)], 3), collapse = ", ")))
  list(size = size, per.sample = per.sample, purity = purity, summary = summary, text = text)
}
