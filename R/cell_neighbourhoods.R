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

#' Heat-kernel smoothing on the cell graph
#'
#' Applies `exp(-beta L / l.max)` (`sccore::heatFilter`, `L` the graph Laplacian, `l.max` its largest eigenvalue) to the
#' columns of `signal` through a Chebyshev approximation of order `m`. A strong kernel is applied in `ceiling(beta /
#' step.max)` equal steps (the kernel composes exactly, `exp(-a L) exp(-b L) = exp(-(a + b) L)`), since a single
#' polynomial of moderate order approximates `exp(-beta x)` poorly for large `beta`; up to `step.max` it is one step,
#' as before. `l.max` is estimated directly, so the graph need not be connected.
#' @param signal cells x columns matrix (rows named by the graph's cells)
#' @param beta kernel strength
#' @param graph igraph cell graph (or `lap`)
#' @param lap graph Laplacian (default: from `graph`)
#' @param l.max its largest eigenvalue (default: estimated)
#' @param m Chebyshev order (default 50)
#' @param n.cores cores (blocks of columns are smoothed in parallel)
#' @param step.max largest strength applied in one step (default 120)
#' @param symmetric whether the Laplacian is symmetric (default: from the graph, or checked)
#' @return smoothed matrix (rows in the order of the Laplacian)
#' @keywords internal
heatSmooth <- function(signal, beta, graph = NULL, lap = NULL, l.max = NULL, m = 50, n.cores = 1, step.max = 120, symmetric = NULL) {
  if (is.null(lap)) { gl <- cachedLaplacian(graph); lap <- gl$lap; if (is.null(l.max)) l.max <- gl$l.max; if (is.null(symmetric)) symmetric <- gl$symmetric }
  if (is.null(l.max)) { checkPackageInstalled("irlba", cran = TRUE); l.max <- irlba::partial_eigen(lap, n = 1)$values }
  vec <- is.null(dim(signal)); X <- if (vec) matrix(signal, ncol = 1, dimnames = list(names(signal), NULL)) else as.matrix(signal)
  if (!is.null(colnames(lap)) && !is.null(rownames(X))) X <- X[colnames(lap), , drop = FALSE] else if (nrow(X) != ncol(lap)) stop("signal must have one value per graph vertex")
  n.steps <- max(1L, as.integer(ceiling(beta / step.max))); b1 <- beta / n.steps
  co <- chebyshevCoefficients(function(x) exp(-b1 * x / l.max), l.max, m)
  # compiled path: the recurrence in C++ on threads (rows of the symmetric Laplacian split over n.cores), no forked copies
  if (isTRUE(getOption("cacoa.cpp.smooth", TRUE)) && (if (!is.null(symmetric)) symmetric else if (!is.null(graph)) !igraph::is_directed(graph) else Matrix::isSymmetric(lap))) {
    lapc <- methods::as(methods::as(methods::as(lap, "CsparseMatrix"), "generalMatrix"), "dMatrix")
    storage.mode(X) <- "double"
    Y <- chebyshev_smooth_cpp(lapc, X, l.max, as.numeric(co), n.steps, as.integer(max(1L, n.cores)))
    dimnames(Y) <- dimnames(X)
    return(if (vec) stats::setNames(Y[, 1], rownames(Y)) else Y)
  }
  smooth1 <- function(Xb) { for (s in seq_len(n.steps)) Xb <- chebyshevSmooth(Xb, lap, l.max, co); Xb }
  nb <- min(as.integer(n.cores), ncol(X))   # the columns are independent: blocks of them are smoothed in parallel (forked)
  X <- if (nb > 1) { res <- parallel::mclapply(parallel::splitIndices(ncol(X), nb), function(b) smooth1(X[, b, drop = FALSE]), mc.cores = nb)
    bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("heatSmooth: ", conditionMessage(attr(res[[which(bad)[1]]], "condition")))
    do.call(cbind, res) } else smooth1(X)
  if (vec) stats::setNames(X[, 1], rownames(X)) else X
}

#' Laplacian of a cell graph and its largest eigenvalue, kept for the last graph seen
#'
#' The kernel-strength search ([solveKernelBeta()]) and the smoothing of densities and neighborhoods all start from the same graph;
#' the Laplacian and its largest eigenvalue are computed once and kept for the last graph (identified by its numbers of
#' vertices and edges, the sum of its edge weights and its first and last vertex names).
#' @param graph igraph cell graph
#' @return list with `lap` (sparse Laplacian) and `l.max`
#' @keywords internal
cachedLaplacian <- function(graph) {
  nv <- igraph::vcount(graph); id <- if (exists("graph_id", envir = asNamespace("igraph"))) igraph::graph_id(graph) else ""
  key <- paste(id, nv, igraph::ecount(graph), if (nv) paste(igraph::V(graph)$name[c(1, nv)], collapse = "|") else "")   # reading the weights costs ~1 s on a large graph
  if (!nzchar(id) || !identical(.laplacianCache$key, key)) {
    lap <- igraph::laplacian_matrix(graph, sparse = TRUE); checkPackageInstalled("irlba", cran = TRUE)
    .laplacianCache$lap <- lap; .laplacianCache$l.max <- irlba::partial_eigen(lap, n = 1)$values; .laplacianCache$symmetric <- !igraph::is_directed(graph); .laplacianCache$key <- key
  }
  list(lap = .laplacianCache$lap, l.max = .laplacianCache$l.max, symmetric = .laplacianCache$symmetric)
}

#' Empty the Laplacian cache ([cachedLaplacian()]); called at the start of every map, so that a graph whose edge weights were
#' changed in place (igraph keeps its id) is never served a stale Laplacian across calls
#' @keywords internal
clearLaplacianCache <- function() { rm(list = ls(.laplacianCache), envir = .laplacianCache); invisible(NULL) }
.laplacianCache <- new.env(parent = emptyenv())

#' Chebyshev coefficients of a spectral filter on [0, l.max] (as sccore::smoothSignalOnGraph computes them)
#' @param filt filter of the eigenvalue
#' @param l.max largest eigenvalue of the Laplacian
#' @param m order
#' @return numeric vector of `m + 1` coefficients
#' @keywords internal
chebyshevCoefficients <- function(filt, l.max, m = 50) {
  n <- m + 1; a <- l.max / 2; tn <- 0:(n - 1); num <- cos(pi * (tn + 0.5) / n)
  sapply(0:m, function(i) 2 / n * sum(filt(a * (num + 1)) * cos(pi * i * (tn + 0.5) / n)))
}

#' Chebyshev polynomial of the Laplacian applied to a dense signal
#'
#' The recurrence of sccore::smoothSignalOnGraph with sparse x dense products on base matrices (the same result; sccore's
#' Matrix arithmetic is far slower on wide signals: a 48,000-cell graph and 100 columns take 24 s here).
#' @param X dense matrix (vertices x signals, in the Laplacian's order)
#' @param lap sparse Laplacian
#' @param l.max its largest eigenvalue
#' @param coeffs Chebyshev coefficients ([chebyshevCoefficients()])
#' @return dense matrix
#' @keywords internal
chebyshevSmooth <- function(X, lap, l.max, coeffs) {
  lap.norm <- lap / (l.max / 2); fac <- 2 * (lap.norm - Matrix::Diagonal(ncol(lap)))
  cur <- as.matrix(lap.norm %*% X) - X; old <- X; r <- 0.5 * coeffs[1] * old + coeffs[2] * cur
  for (k in 3:length(coeffs)) { new <- as.matrix(fac %*% cur) - old; r <- r + coeffs[k] * new; old <- cur; cur <- new }
  dimnames(r) <- dimnames(X); r
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
#' The same probes give the kernel's diagonal, `K_ii = E[v_i (K v)_i]` (`self`: the share of a cell's heat that stays on
#' the cell), and the effective number of cells the heat arriving from the other cells represents,
#' `(r_i - K_ii)^2 / sum_{k != i} K_ik^2` (`off.size`), which the density map uses to count the cell as one cell of its
#' neighborhood ([heatLogDensity()]). The diagonal estimate is noisy where it is small (a well-connected cell keeps about
#' one cell's share, below the probe's standard error `sqrt(sum_{k != i} K_ik^2 / n.vec)`), so `self` is the estimate
#' less two standard errors, floored at one cell's share `1 / size`: a weakly connected cell's self-mass is taken as
#' measured, a well-connected cell's as the one cell it is.
#' @param graph igraph cell graph
#' @param beta kernel strength (as in `estimateCellDensity`)
#' @param n.vec number of probe vectors (default 24)
#' @param m Chebyshev order (default 50)
#' @param n.cores cores for the smoothing
#' @param seed seed of the probe
#' @return list with `size`, `row.sum`, `self` and `off.size` (per cell) and `beta`
#' @export
kernelEffectiveSize <- function(graph, beta, n.vec = 24, m = 50, n.cores = 1, seed = 7) {
  cells <- igraph::V(graph)$name
  R <- withSeed(seed, cbind(1, matrix(sample(c(-1, 1), length(cells) * n.vec, replace = TRUE), length(cells), n.vec)))
  rownames(R) <- cells
  KR <- as.matrix(heatSmooth(R, beta, graph = graph, m = m, n.cores = n.cores))
  c(probeSummaries(R, KR, cells), list(beta = beta))
}

#' Summaries of the heat kernel from its probe vectors
#'
#' From the smoothed probes `KR = K R` (`R` a column of ones and `n.vec` columns of random +-1): the row sums, the
#' effective neighborhood size, the kernel's diagonal and the effective size of the heat arriving from the other cells.
#' @param R probe matrix (cells x (n.vec + 1), the first column ones)
#' @param KR smoothed probes (same shape)
#' @param cells cell names of the rows
#' @return list with `size`, `row.sum`, `self` and `off.size` (per cell, named)
#' @keywords internal
probeSummaries <- function(R, KR, cells) {
  r <- KR[, 1]; P <- KR[, -1, drop = FALSE]; Rp <- R[, -1, drop = FALSE]; size <- r^2 / rowMeans(P^2)
  self <- pmin(pmax(rowMeans(Rp * P), 0), pmax(r, 0))          # K_ii = E[v_i (K v)_i]
  off2 <- pmax(rowMeans((P - self * Rp)^2), 1e-12)             # sum_{k != i} K_ik^2 = E[((K v)_i - K_ii v_i)^2]
  self <- pmin(pmax(self - 2 * sqrt(off2 / ncol(Rp)), pmax(r, 0) / pmax(size, 1)), pmax(r, 0))   # less two standard errors, at least one cell's share
  list(size = stats::setNames(size, cells), row.sum = stats::setNames(r, cells),
       self = stats::setNames(self, cells), off.size = stats::setNames(pmax(r - self, 1e-8)^2 / off2, cells))
}

#' Kernel strength for a target effective neighbourhood size
#'
#' Brackets the target from `start` in steps of 4x, then bisection on `log(beta)` so that the median effective size (`kernelEffectiveSize`) matches `target.size`. The
#' same `beta` gives very different widths on different graphs (a 20-fold spread between datasets in the benchmark),
#' so the width is better stated in cells.
#' @inheritParams kernelEffectiveSize
#' @param target.size target median effective size in cells
#' @param range search interval for beta (default c(0.5, 5000))
#' @param start first strength evaluated; the target is bracketed from there in steps of 4x (default 30)
#' @param tol relative tolerance on the size (default 0.1)
#' @param max.iter bisection steps (default 12)
#' @return list with `beta`, `size`, `self` and `off.size` (per cell at that beta, [kernelEffectiveSize()]) and the search `trace`
#' @export
solveKernelBeta <- function(graph, target.size, range = c(0.5, 5000), tol = 0.1, max.iter = 12, n.vec = 24, m = 50, n.cores = 1, verbose = FALSE, start = 30) {
  trace <- NULL; best <- NULL; L <- log(range)
  f <- function(lb) { e <- kernelEffectiveSize(graph, exp(lb), n.vec = n.vec, m = m, n.cores = n.cores); e$median <- stats::median(e$size)
    trace <<- rbind(trace, data.frame(beta = exp(lb), median.size = e$median)); e }
  # bracket the target from `start` in steps of 4x (a strong kernel costs more smoothing steps, so it is evaluated only when needed)
  lb <- min(max(log(start), L[1]), L[2]); e <- f(lb)
  if (e$median < target.size) { lo <- lb; e.lo <- e; repeat { hi <- min(lo + log(4), L[2]); e.hi <- f(hi); if (e.hi$median >= target.size || hi >= L[2]) break; lo <- hi; e.lo <- e.hi } }
  else { hi <- lb; e.hi <- e; repeat { lo <- max(hi - log(4), L[1]); e.lo <- f(lo); if (e.lo$median <= target.size || lo <= L[1]) break; hi <- lo; e.hi <- e.lo } }
  if (target.size <= e.lo$median) return(c(list(beta = exp(lo)), e.lo[c("size", "self", "off.size")], list(trace = trace)))
  if (target.size >= e.hi$median) return(c(list(beta = exp(hi)), e.hi[c("size", "self", "off.size")], list(trace = trace)))
  for (it in seq_len(max.iter)) {
    mid <- (lo + hi) / 2; e <- f(mid)
    if (verbose) message(sprintf("beta %.3g: median effective size %.0f", exp(mid), e$median))
    best <- c(list(beta = exp(mid)), e[c("size", "self", "off.size")])
    if (abs(log(e$median / target.size)) <= log1p(tol)) break
    if (e$median < target.size) lo <- mid else hi <- mid
  }
  best$trace <- trace; best
}

#' The candidate scales of `scale = "auto"`
#'
#' @param auto.scales numeric vector of shares of the cells (at least two distinct values in (0, 1))
#' @return the sorted distinct scales
#' @keywords internal
checkAutoScales <- function(auto.scales) {
  if (!is.numeric(auto.scales) || !all(is.finite(auto.scales)) || any(auto.scales <= 0 | auto.scales >= 1)) stop("auto.scales must be shares of the cells in (0, 1)")
  sc <- sort(unique(auto.scales)); if (length(sc) < 2) stop("auto.scales needs at least two distinct scales to choose from (one scale: scale = ", auto.scales[1], ")")
  sc
}

#' Scales with a floor on the neighbourhood size
#'
#' A scale (share of the cells) whose neighbourhood would hold fewer than `floor.cells` cells is raised to the floor and
#' duplicates are dropped: the shift map's per-neighbourhood test needs enough cells per sample
#' ([defaultNeighbourhoodSize()]: 10 per sample, at least 100, capped at 1 % of the cells), the density map a handful of
#' cells. On a large dataset the default scales are unchanged (1 % of 35,000 cells is 350).
#' @param scale scales (shares of the cells)
#' @param n.cells number of cells
#' @param floor.cells smallest neighbourhood in cells
#' @param verbose say when a scale was raised
#' @return sorted unique scales
#' @keywords internal
resolveScales <- function(scale, n.cells, floor.cells, verbose = FALSE) {
  if (!is.numeric(scale) || !length(scale) || any(!is.finite(scale) | scale <= 0 | scale >= 1)) stop("scale must be shares of the cells in (0, 1), or \"auto\"")
  fl <- min(0.5, floor.cells / n.cells); raised <- scale < fl; out <- sort(unique(pmax(scale, fl)))
  if (any(raised) && verbose) message(sprintf("scale%s %s below the floor of %d cells: %s", if (sum(raised) > 1) "s" else "", paste(sprintf("%g", scale[raised]), collapse = ", "), as.integer(floor.cells),
                                              if (length(out) < length(scale)) sprintf("the scales are now %s", paste(sprintf("%g", out), collapse = ", ")) else sprintf("raised to %g", fl)))
  out
}

#' Diffusion neighbourhoods at a scale
#'
#' Neighbourhoods of the expression-shift map on a diffusion kernel of the cell graph. The kernel strength is solved so
#' that the median cell's effective neighbourhood ([kernelEffectiveSize()]) holds `scale` of the cells
#' ([solveKernelBeta()]), which also gives every cell's reach, the effective size of its own kernel. Cells are taken in a
#' random order (seeded); a cell in no neighbourhood yet whose reach is at least `min.reach` of the scale's cells becomes
#' a centre, and its neighbourhood is the scale's number of cells (fewer if its kernel reaches fewer, at least `min.size`)
#' that its heat reaches most. Only the centres' kernels are computed (`pool` candidates per smoothing). Afterwards every
#' cell joins, up to `cover` neighbourhoods in all, the centres whose kernel weight at the cell is at least half the
#' neighbourhood's boundary weight (the kernel is symmetric, so that is the cell's own weight at the centre), each
#' neighbourhood taking at most its own size of such cells, the best reached first; a cell in no neighbourhood (a weakly
#' connected cell, typically one without between-sample neighbours whose heat stays on itself, or one of a population far
#' smaller than the scale) joins the centre whose kernel reaches it most, and only a cell no kernel reaches stays out. So
#' every neighbourhood holds about the scale's number of cells and nearly every cell lies in one or more.
#' Indices refer to the graph's vertex order (`cells`). [scanShiftScales()] sets `min.size` to the kNN unit size
#' ([defaultNeighbourhoodSize()]) so that every neighbourhood holds enough cells per sample to be tested.
#' @param graph igraph cell graph
#' @param scale share of the cells a neighbourhood holds
#' @param cover neighbourhoods per cell the joins aim at (default 2)
#' @param seed seed of the centre order
#' @param beta kernel strength (default NULL: solved for the scale)
#' @param pool candidate centres per smoothing (default 100)
#' @param min.size smallest neighbourhood (default 10; the shift scan passes the kNN unit size)
#' @param min.reach smallest reach of a centre's kernel as a share of the scale's cells (default 1/3)
#' @param m Chebyshev order (default 50)
#' @param n.cores cores
#' @return list: `units` (per centre, 0-based indices of its cells in the graph's order, the centre included; named
#'   by the centre's 1-based index), `centres` (1-based indices), `beta`, `size` (cells of each centre's own
#'   neighbourhood before the joins), `cover` (neighbourhoods per cell), `cells`, `n.smooth` (smoothings run)
#' @export
diffusionNeighbourhoodLists <- function(graph, scale, cover = 2L, seed = 1, beta = NULL, pool = 100L, min.size = 10L, min.reach = 1/3, m = 50, n.cores = 1) {
  cells <- igraph::V(graph)$name; N <- length(cells); target <- round(scale * N)
  if (is.null(beta)) { sb <- solveKernelBeta(graph, target.size = target, n.vec = 16, m = m, n.cores = n.cores); beta <- sb$beta; reach <- sb$size }
  else reach <- kernelEffectiveSize(graph, beta, n.vec = 16, m = m, n.cores = n.cores)$size
  gl <- cachedLaplacian(graph); lap <- gl$lap; l.max <- gl$l.max; sym <- gl$symmetric
  wide <- reach >= min.reach * target   # may become a centre
  cov <- integer(N); units <- list(); centres <- integer(0); sizes <- numeric(0); thr <- numeric(N); todo <- withSeed(seed, sample.int(N))
  best.w <- matrix(0, N, 2); best.c <- matrix(0L, N, 2); n.smooth <- 0L   # the two centres reaching each cell most, and their weights there
  while (length(todo)) {
    cand <- todo[seq_len(min(pool, length(todo)))]; todo <- todo[-seq_len(length(cand))]
    cand <- cand[cov[cand] < 1L & wide[cand]]; if (!length(cand)) next
    I <- matrix(0, N, length(cand)); I[cbind(cand, seq_along(cand))] <- 1; rownames(I) <- cells
    K <- as.matrix(heatSmooth(I, beta, lap = lap, l.max = l.max, m = m, n.cores = n.cores, symmetric = sym)); K[K < 0] <- 0; n.smooth <- n.smooth + 1L
    for (j in seq_along(cand)) { c0 <- cand[j]; if (cov[c0] >= 1L) next   # covered by a centre earlier in the pool
      w0 <- K[, j]; w0[c0] <- 0; s0 <- sum(w0); ne <- if (s0 > 0) round(s0^2 / sum(w0^2)) else 0
      if (ne == 0) next
      top <- order(w0, decreasing = TRUE)[seq_len(min(max(min.size, min(ne, target)), N - 1L))]; thr[c0] <- w0[top[length(top)]] / 2
      u <- sort(unique(c(c0, top))); units[[as.character(c0)]] <- u - 1L; centres <- c(centres, c0); sizes <- c(sizes, length(top)); cov[u] <- cov[u] + 1L
      b1 <- w0 > best.w[, 1]; b2 <- !b1 & w0 > best.w[, 2]
      best.w[b1, 2] <- best.w[b1, 1]; best.c[b1, 2] <- best.c[b1, 1]; best.w[b1, 1] <- w0[b1]; best.c[b1, 1] <- c0
      best.w[b2, 2] <- w0[b2]; best.c[b2, 2] <- c0
    }
  }
  room <- stats::setNames(sizes, centres)   # joins per neighbourhood: at most its own size
  for (col in seq_len(min(2L, cover))) { cc <- best.c[, col]; ww <- best.w[, col]
    left <- which(cov < cover & ww > 0 & (ww >= thr[pmax(cc, 1L)] | cov == 0L))
    for (c0 in unique(cc[left])) { key <- as.character(c0); i <- left[cc[left] == c0]; i <- i[!((i - 1L) %in% units[[key]])]; i <- i[order(ww[i], decreasing = TRUE)]
      take <- i[cov[i] == 0L | seq_along(i) <= room[key]]; room[key] <- room[key] - length(take)   # a cell in no neighbourhood always joins
      if (length(take)) { units[[key]] <- sort(c(units[[key]], take - 1L)); cov[take] <- cov[take] + 1L } } }
  list(units = units, centres = centres, beta = beta, size = stats::setNames(sizes, centres), cover = stats::setNames(cov, cells), cells = cells, n.smooth = n.smooth)
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
