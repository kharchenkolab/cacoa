#' @importFrom sccore checkPackageInstalled
NULL

#' Estimate cell density in giving embedding, density will be estimated for individual samples
#'
#' @param emb cell embedding matrix
#' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
#' @param sample.weights Named numeric vector: per-sample signed weights.
#' @param bins number of bins for density estimation, default 400
#' @keywords internal
estimateCellDensityKde <- function(emb, sample.per.cell, sample.weights, bins, bandwidth=0.05, expansion.mult=0.05) {
  if (is.null(bandwidth)) {
    bandwidth <- apply(emb, 2, MASS::bandwidth.nrd)
  } else {
    bandwidth <- apply(emb, 2, quantile, c(0.1, 0.9)) %>% diff() %>% {. * bandwidth} %>% .[1,]
  }

  expand_limits_continuous <- utils::getFromNamespace("expand_limits_continuous", "ggplot2")

  lims <- emb %>%
    apply(2, function(x) expand_limits_continuous(range(x), expansion(mult=expansion.mult))) %>%
    as.numeric()

  cname <- intersect(names(sample.per.cell), rownames(emb))
  sample.per.cell <- sample.per.cell[cname]
  emb <- emb[cname, ]
  cells.per.samp <- split(names(sample.per.cell), sample.per.cell)
  density.mat <-lapply(cells.per.samp, function(nn) {
      MASS::kde2d(emb[nn, 1], emb[nn, 2], h=bandwidth, n=bins, lims=lims)$z %>%
      as.numeric() %>% {. / sum(.)}
    }) %>% do.call("cbind", .) %>%
    set_colnames(names(cells.per.samp)) %>%
    set_rownames(1:nrow(.)) # needed for indexing in diffCellDensity

  # ---------- condition densities from weights ----------
  # align weights to density.mat columns
  w <- sample.weights[colnames(density.mat)]
  if (any(is.na(w))) {
    warning("Some samples in score.mat have no weights; treating them as 0.")
    w[is.na(w)] <- 0
  }

  w.num <- pmax(w,  0)
  w.den <- pmax(-w, 0)

  if (sum(w.num) > 0) w.num <- w.num / sum(w.num)
  if (sum(w.den) > 0) w.den <- w.den / sum(w.den)
  dens.num <- as.numeric(density.mat %*% w.num)
  dens.den <- as.numeric(density.mat %*% w.den)

  cond.densities <- list(num = dens.num, den = dens.den)
  # ------------------------------------------------------

  # coordinate embedding space
  mat <- matrix(cond.densities[[1]], ncol = bins, byrow = FALSE)
  x1 <- seq(lims[1], lims[2], length.out=bins) %>% setNames(seq(bins))
  y1 <- seq(lims[3], lims[4], length.out=bins) %>% setNames(seq(bins))
  d1 <- setNames(melt(mat), c('x', 'y', 'z'))
  emb2 <- data.frame(x=x1[d1$x], y=y1[d1$y])

  #count cell number in each bin
  s1 <- seq(from = lims[1], to = lims[2], length.out=bins + 1)
  s2 <- seq(from = lims[3], to = lims[4], length.out=bins + 1)
  dcounts <- table(cut(emb[,1], breaks = s1), cut(emb[,2], breaks = s2)) #%>% as.matrix.data.frame
  emb2$counts <- as.numeric(dcounts)

  return(list(density.mat=density.mat, cond.densities=cond.densities, density.emb=emb2, bins=bins,
              method='kde', cell.emb=emb))
}


#' Estimate graph smooth based cell density
#'
#' @param graph input graph
#' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
#' @param sample.weights Named numeric vector: per-sample signed weights.
#' @param n.cores number of cores
#' @param m numeric Maximum order of Chebyshev coeff to compute (default=50)
#' @param smoothing `"knn"` (default): the number of each sample's cells in the cell's kNN neighbourhood of `k`
#'   cells (itself, its graph neighbours by edge weight, then two-hop neighbours by shared neighbours; see
#'   [knnNeighbourhoodMatrix()]); `"heat"`: the heat-kernel smooth of each sample's cell indicator (a soft
#'   neighbourhood whose width is set by `beta` or, better, by `target.size` cells; the cell itself counted as one
#'   cell of its neighbourhood, [ownCellDensity()]); `"neighbourhood"`: the counts
#'   within `n.hops` graph hops. The raw counts of the count units are kept in `counts.mat`
#' @param k neighbourhood size for `smoothing = "knn"` (default: [defaultNeighbourhoodSize()])
#' @param target.size heat kernel only: median effective neighbourhood size in cells; when given, `beta` is solved
#'   by [solveKernelBeta()]. The effective size at the `beta` used is returned as `effective.size`
#' @param n.hops neighbourhood radius for `smoothing = "neighbourhood"` (default 1)
#' @param neighbourhoods precomputed kNN lists ([knnNeighbourhoodLists()]) for `smoothing = "knn"` (cached by the class)
#' @return list with `density.mat` (cells x samples: the fraction of each sample's cells in the unit), `cond.densities`,
#'   `counts.mat` (count units), `sample.sizes`, `unit.size` (cells per unit: exact for the count units, the
#'   effective size for the kernel), `beta` and `smoothing`
#' @keywords internal
estimateCellDensityGraph <- function(graph, sample.per.cell, sample.weights, n.cores=1, beta=30, m=50, verbose=TRUE,
                                     smoothing=c("knn", "heat", "neighbourhood"), k=NULL, target.size=NULL, n.hops=1, neighbourhoods=NULL) {
  smoothing <- match.arg(smoothing)
  sig.mat <- unique(sample.per.cell) %>% sapply(function(s) as.numeric(sample.per.cell == s)) %>%
    set_rownames(names(sample.per.cell)) %>% set_colnames(unique(sample.per.cell))
  counts.mat <- NULL; sample.sizes <- colSums(sig.mat); unit.size <- NULL

  if (smoothing == "knn") {
    cells <- intersect(igraph::V(graph)$name, rownames(sig.mat))
    k <- k %||% defaultNeighbourhoodSize(length(cells), ncol(sig.mat))
    if (verbose) message(sprintf("kNN neighborhoods of at least %d cells (%.1f per sample)", k, k / ncol(sig.mat)))
    nn <- neighbourhoods %||% knnNeighbourhoodLists(graph, k, cells = cells, n.cores = n.cores)
    if (!identical(names(nn), cells)) stop("the neighborhood lists do not match the graph's cells")
    counts.mat <- neighbourhoodCounts(nn, sample.per.cell[cells], samples = colnames(sig.mat))
    unit.size <- stats::setNames(as.numeric(lengths(nn)), cells)
    score.mat <- t(t(counts.mat) / sample.sizes[colnames(counts.mat)])
  } else if (smoothing == "neighbourhood") {
    cells <- intersect(igraph::V(graph)$name, rownames(sig.mat))
    A <- igraph::as_adj(graph, sparse = TRUE)[cells, cells, drop = FALSE]
    A@x[] <- 1; Matrix::diag(A) <- 1
    M <- A
    if (n.hops > 1) for (h in seq_len(n.hops - 1)) { M <- M %*% A; M@x[] <- 1 }
    counts.mat <- as.matrix(M %*% sig.mat[cells, , drop = FALSE])
    rownames(counts.mat) <- cells
    unit.size <- stats::setNames(as.numeric(Matrix::rowSums(M)), cells)
    score.mat <- t(t(counts.mat) / sample.sizes[colnames(counts.mat)])   # fraction of the sample's cells in the neighbourhood
  } else {
    if (!is.null(target.size)) {
      sb <- solveKernelBeta(graph, target.size, n.cores = n.cores, verbose = verbose); beta <- sb$beta; unit.size <- sb$size
      if (verbose) message(sprintf("heat kernel beta = %.3g for a median effective size of %.0f cells (target %d)", beta, stats::median(unit.size), as.integer(target.size)))
      kes <- sb
    } else {
      kes <- kernelEffectiveSize(graph, beta, n.cores = n.cores); unit.size <- kes$size
      if (verbose) message(sprintf("heat kernel beta = %.3g: median effective neighborhood size %.0f cells (5-95 %%: %.0f-%.0f)", beta, stats::median(unit.size), stats::quantile(unit.size, 0.05), stats::quantile(unit.size, 0.95)))
    }
    score.mat <- as.matrix(heatSmooth(sig.mat, beta, graph = graph, m = m, n.cores = n.cores)); score.mat[!is.finite(score.mat) | score.mat < 0] <- 0
    cells <- rownames(score.mat); si <- match(as.character(sample.per.cell[cells]), colnames(score.mat))
    score.mat <- ownCellDensity(score.mat, si, kes$self[cells], kes$off.size[cells])   # the cell as one cell of its neighborhood
    score.mat <- t(t(score.mat) / sample.sizes[colnames(score.mat)])   # density per cell of the sample
  }


  # ---------- condition densities from weights ----------
  w <- sample.weights[colnames(score.mat)]
  if (any(is.na(w))) {
    warning("Some samples in score.mat have no weights; treating them as 0.")
    w[is.na(w)] <- 0
  }

  w.num <- pmax(w,  0)
  w.den <- pmax(-w, 0)

  if (sum(w.num) > 0) w.num <- w.num / sum(w.num)
  if (sum(w.den) > 0) w.den <- w.den / sum(w.den)
  dens.num <- as.numeric(score.mat %*% w.num)  # cells x samples * samples
  dens.den <- as.numeric(score.mat %*% w.den)

  names(dens.num) <- names(dens.den) <- rownames(score.mat)

  cond.densities <- list(num = dens.num, den = dens.den)
  # ------------------------------------------------------

  return(list(density.mat=score.mat, cond.densities=cond.densities, method='graph', smoothing=smoothing, counts.mat=counts.mat, sample.sizes=sample.sizes,
              unit.size=unit.size, beta=if (smoothing == "heat") beta else NULL, k=if (smoothing == "knn") k else NULL, neighbourhoods=if (smoothing == "knn") nn else NULL))
}

#' The cell as one cell of its heat-kernel neighborhood
#'
#' The smoothed per-sample indicators give, at each cell, the heat arriving from every sample's cells; the cell's own
#' heat (`K_ii`, `self`) counts toward its own sample. On a graph of uneven degree a weakly connected cell keeps most of
#' its heat, so its own sample would stand alone at the cell whatever the labels say, while the heat that does arrive
#' carries the composition of the neighborhood it is matched to. The cell is therefore counted as one cell of a
#' neighborhood of `off.size + 1` cells: the arriving heat (the row without `K_ii`) rescaled to `off.size` cells, the
#' effective number of cells it represents ([kernelEffectiveSize()]), plus one cell of the cell's own sample. For a
#' well-connected cell `K_ii` is about one cell's share already and nothing changes.
#' @param d cells x samples matrix of smoothed sample indicators
#' @param si the column of each cell's own sample
#' @param self the kernel's diagonal per cell
#' @param off.size the effective size of the arriving heat per cell
#' @return cells x samples matrix, rows summing to 1
#' @keywords internal
ownCellDensity <- function(d, si, self, off.size) {
  ii <- cbind(seq_len(nrow(d)), si); a <- d; a[ii] <- pmax(a[ii] - self, 0)
  asum <- rowSums(a); h <- a * ifelse(asum > 0, off.size / asum, 0); h[ii] <- h[ii] + 1
  h / (off.size + 1)
}

#' Broad heat-kernel unit of the density test
#'
#' Per-sample cell indicators smoothed on the cell graph with the heat kernel (`sccore::heatFilter`, strength `beta`),
#' the cell counted as one cell of its neighborhood ([ownCellDensity()]), normalized per sample to a density per cell of
#' the sample (mean `1 / n.cells`) and put on the log scale with a relative pseudo-count, `log2(d + alpha / n.cells)`: a
#' kernel density is far below one cell, so the count pseudo-count of the kNN units (half a cell) would swamp it. The
#' kernel's diagonal and the effective size of the arriving heat (`self`, `off.size`) come from the probe vectors of the
#' strength search ([solveKernelBeta()]) when given; otherwise `n.vec` probes are smoothed together with the indicators
#' in the same call. The largest eigenvalue of the Laplacian is estimated directly, so the graph need not be connected
#' (the kernel stays within each component); a strong kernel is applied in steps ([heatSmooth()]). Chebyshev ripples
#' below zero are set to zero.
#' @param graph igraph cell graph
#' @param sample.per.cell factor of samples named by cell
#' @param beta kernel strength (default 120)
#' @param alpha pseudo-count as a fraction of the mean density (default 0.3)
#' @param m Chebyshev order (default 50)
#' @param n.cores cores for the smoothing
#' @param self,off.size per cell, named ([kernelEffectiveSize()]); default NULL: estimated here
#' @param n.vec probe vectors when `self` and `off.size` are estimated here (default 16, the strength search's, so that a given strength gives the same map as a solved one)
#' @return cells x samples matrix of `log2` densities (the graph's cells that have a sample)
#' @keywords internal
heatLogDensity <- function(graph, sample.per.cell, beta = 120, alpha = 0.3, m = 50, n.cores = 1, self = NULL, off.size = NULL, n.vec = 16) {
  cells <- intersect(igraph::V(graph)$name, names(sample.per.cell)[!is.na(sample.per.cell)])
  if (length(cells) < 2) stop("fewer than two graph cells have a sample")
  if (length(cells) < igraph::vcount(graph)) graph <- igraph::induced_subgraph(graph, cells)
  cells <- igraph::V(graph)$name; spc <- droplevels(as.factor(sample.per.cell[cells])); si <- as.integer(spc)
  sig <- matrix(0, length(cells), nlevels(spc), dimnames = list(cells, levels(spc))); sig[cbind(seq_along(cells), si)] <- 1
  if (is.null(self) || is.null(off.size) || !all(cells %in% names(self)) || !all(cells %in% names(off.size))) {   # the probes and the indicators in one smoothing
    R <- withSeed(7, cbind(1, matrix(sample(c(-1, 1), length(cells) * n.vec, replace = TRUE), length(cells), n.vec))); rownames(R) <- cells
    KR <- as.matrix(heatSmooth(cbind(R, sig), beta, graph = graph, m = m, n.cores = n.cores))[cells, , drop = FALSE]
    ps <- probeSummaries(R, KR[, seq_len(n.vec + 1), drop = FALSE], cells); self <- ps$self; off.size <- ps$off.size
    d <- KR[, -seq_len(n.vec + 1), drop = FALSE]
  } else { d <- as.matrix(heatSmooth(sig, beta, graph = graph, m = m, n.cores = n.cores))[cells, , drop = FALSE]; self <- self[cells]; off.size <- off.size[cells] }
  d[!is.finite(d) | d < 0] <- 0
  d <- ownCellDensity(d, si, self, off.size)
  ns <- colSums(sig); ns[ns <= 0] <- 1; d <- t(t(d) / ns)
  log2(d + alpha / nrow(d))
}

#' Scale scan of the density map
#'
#' Diffusion neighborhoods on the cell graph at one or several scales, the scale being the share of the cells a typical
#' neighborhood covers (the kernel strength is solved per graph, [solveKernelBeta()]). At each scale the per-sample
#' kernel densities ([heatLogDensity()]) are tested with [diffCellDensityPermutations()] under the same relabelings, the
#' statistic is averaged over the cell and its graph neighbours (`post = "averaged"`, [averageOverNeighbours()]) and
#' turned into a p-value against that scale's pooled permutation null (every cell under every relabeling). A cell's
#' scan statistic is its smallest p over the scales; its q-value is the pooled permutation false discovery rate of
#' that statistic, i.e. the expected number of cells at or below it under relabeling (the smallest p over scales of
#' every cell under every relabeling) over the number observed, so the choice among scales costs no Bonferroni
#' factor. Each scale's own permutation FDR is kept as well (`Q`): the single-scale map at that scale.
#'
#' Memory: the permuted statistics of one scale (`n.permutations` x cells) and the running minimum over scales are held
#' at once.
#' @param graph igraph cell graph
#' @param sample.per.cell named factor of samples per cell (cells without a sample are left out)
#' @param sample.model sample model (see [buildCacoaModel()])
#' @param scale scales: shares of the cells (one or several)
#' @param beta kernel strength per scale (default NULL: solved for the scale by [solveKernelBeta()])
#' @param self,off.size lists per scale of the kernel's diagonal and the arriving heat's effective size per cell
#'   ([kernelEffectiveSize()]), for given strengths; default NULL: from the strength search, or estimated per scale
#' @param alpha pseudo-count of the kernel density on the log scale, as a fraction of the mean density ([heatLogDensity()])
#' @param post `"averaged"` (default) or `"none"`
#' @param alternative `"two-sided"` (default), `"greater"` or `"less"`
#' @param n.permutations number of relabelings
#' @param seed seed of the relabelings (the same for every scale)
#' @param groups optional cell groups: the group-level permutation test ([groupPermutationTest()]) per scale
#' @param n.cores cores
#' @param verbose print the strength per scale
#' @param ... further arguments of [diffCellDensityPermutations()] (`perm.method`, `robust.method`, `na.mode`, `moderate.perm`, ...)
#' @return list: `q` and `p` (the scan: pooled permutation FDR and p of the smallest p over scales), `best` (index of the
#'   scale of the smallest p per cell), `P`, `E`, `Q` (cells x scales: the per-scale pooled p, effect (log2 fold change)
#'   and permutation FDR), `scales` (data.frame: `scale`, `beta`, `size` (median effective size when solved)),
#'   `group.test` (list per scale), `cells`, `B` (relabelings)
#' @export
scanDensityScales <- function(graph, sample.per.cell, sample.model, scale, beta = NULL, self = NULL, off.size = NULL, alpha = 0.3, post = c("averaged", "none"), alternative = c("two-sided", "greater", "less"),
                              n.permutations = 499, seed = 1, groups = NULL, n.cores = 1, verbose = FALSE, ...) {
  post <- match.arg(post); alternative <- match.arg(alternative)
  if (!is.numeric(scale) || !length(scale) || any(!is.finite(scale) | scale <= 0 | scale >= 1)) stop("scale must be shares of the cells in (0, 1)")
  cells <- intersect(igraph::V(graph)$name, names(sample.per.cell)[!is.na(sample.per.cell)])
  if (length(cells) < 2) stop("fewer than two graph cells have a sample")
  if (length(cells) < igraph::vcount(graph)) graph <- igraph::induced_subgraph(graph, cells)
  cells <- igraph::V(graph)$name; N <- length(cells); K <- length(scale); spc <- droplevels(as.factor(sample.per.cell[cells])); names(spc) <- cells
  if (is.null(beta)) beta <- rep(NA_real_, K) else if (length(beta) != K) stop("beta must have one strength per scale")
  tail.stat <- switch(alternative, "two-sided" = abs, greater = identity, less = function(x) -x)
  P <- matrix(NA_real_, N, K, dimnames = list(cells, NULL)); E <- P; Q <- P; m.perm <- NULL; gt <- vector("list", K); B <- NULL; n.all <- NULL
  info <- data.frame(scale = scale, beta = beta, size = NA_real_); t0 <- Sys.time()
  for (j in seq_len(K)) {
    if (!is.finite(info$beta[j])) { sb <- solveKernelBeta(graph, target.size = round(scale[j] * N), n.vec = 16, n.cores = n.cores); info$beta[j] <- sb$beta; info$size[j] <- stats::median(sb$size); sj <- sb$self; oj <- sb$off.size }
    else { sj <- if (is.null(self)) NULL else self[[j]]; oj <- if (is.null(off.size)) NULL else off.size[[j]] }
    y <- heatLogDensity(graph, spc, beta = info$beta[j], alpha = alpha, n.cores = n.cores, self = sj, off.size = oj)[cells, , drop = FALSE]
    if (verbose) message(sprintf("scale %g: heat kernel beta %.3g%s", scale[j], info$beta[j], if (is.finite(info$size[j])) sprintf(" (median %.0f cells)", info$size[j]) else ""))
    d <- diffCellDensityPermutations(y, sample.model, n.permutations = n.permutations, seed = seed, alternative = alternative, pooled = FALSE, average = FALSE,
                                     return.sampled.stats = TRUE, n.cores = n.cores, verbose = FALSE, ...)
    sp <- d$stats.perm; if (ncol(sp) != N) sp <- t(sp); if (ncol(sp) != N) stop("the permuted statistics do not cover the cells"); B <- nrow(sp)
    st <- stats::setNames(as.numeric(d$stat.obs), cells); ef <- stats::setNames(as.numeric(d$effect), cells); rm(d)
    if (post == "averaged") { av <- averageOverNeighbours(st, sp, graph, effect = ef); st <- av$stat; sp <- av$stats.perm; ef <- av$effect; rm(av) }
    if (!is.null(groups)) gt[j] <- list(tryCatch(groupPermutationTest(st, sp, groups, alternative = if (alternative == "two-sided") "two.sided" else "greater"), error = function(e) NULL))
    if (is.null(n.all)) n.all <- matrix(0L, nrow(sp) + 1L, K); n.all[, j] <- scanCallCounts(st, sp, tail = tail.stat, n.cores = n.cores)   # calls of every data set, for the test of the chosen scale
    ps <- pooledScanStep(st, sp, tail.stat, n.cores = n.cores); rm(sp)
    P[, j] <- ps$p; Q[, j] <- ps$q; E[, j] <- ef; m.perm <- if (is.null(m.perm)) ps$pp else pmin(m.perm, ps$pp); rm(ps); invisible(gc())
  }
  m.obs <- apply(P, 1, min); best <- max.col(-P, ties.method = "first")
  sc <- pooledMinPFDR(m.obs, m.perm, B, n.cores = n.cores); rm(m.perm)
  auto <- autoScaleTest(n.all[1, ], n.all[-1, , drop = FALSE], sum(rowSums(is.finite(Q)) > 0))   # the test of the automatically chosen scale
  list(q = stats::setNames(sc$q, cells), p = stats::setNames(sc$p, cells), best = stats::setNames(best, cells), P = P, E = E, Q = Q, scales = info, group.test = gt, cells = cells, B = B, auto = auto,
       time = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

#' The scale of a scan to report as a single map
#'
#' The scale whose own map calls the most cells at level `alpha`. When some scale calls at least `saturated` of the
#' tested cells, the finest such scale (the scales are in increasing order): the coarser ones add no calls and only lose
#' resolution. When no scale calls a cell, the scale that is the per-cell best scale of the scan's strongest `top` share
#' of the cells.
#' @param sc result of [scanDensityScales()] or [scanShiftScales()] (any list with `Q`, `p`, `best`)
#' @param alpha level of the calls (default 0.1)
#' @param saturated share of the tested cells from which a scale counts as calling everything (default 0.9)
#' @param top share of the cells used when nothing is called (default 0.01)
#' @return index of the scale
#' @keywords internal
chooseScanScale <- function(sc, alpha = 0.1, saturated = 0.9, top = 0.01) {
  Q <- as.matrix(sc$Q); K <- ncol(Q); n <- colSums(Q < alpha, na.rm = TRUE); tested <- sum(rowSums(is.finite(Q)) > 0)
  if (!any(n > 0)) { called <- rank(sc$p, ties.method = "first") <= ceiling(top * length(sc$p)); return(as.integer(names(which.max(table(factor(sc$best[called], levels = seq_len(K))))))) }
  sat <- which(n >= saturated * tested); if (length(sat)) sat[1] else as.integer(which.max(n))
}

#' Number of cells called at one scale by the observed and by every relabeled data set, all treated alike
#'
#' For the test of the automatically chosen scale ([autoScaleTest()]) the observed data set must be handled exactly as the relabeled
#' ones. Each of the B + 1 data sets (observed, relabeled) is put through the scale's permutation FDR against the other B, and the
#' cells below `alpha` are counted. With `standardize`, the statistic of every cell is first standardized by its mean and standard
#' deviation over all B + 1 data sets (the map itself standardizes the observed statistic by the relabeled ones only, and each
#' relabeled one by a mean and standard deviation that include it, which makes the relabeled data sets slightly less extreme).
#' @param s.obs observed per-cell statistics (large values are extreme, after `tail`)
#' @param sp relabeled statistics (relabelings x cells)
#' @param alpha level of the calls (default 0.1)
#' @param tail function mapping the statistic to its extreme tail
#' @param standardize standardize every cell over the B + 1 data sets (default FALSE)
#' @param n.cores threads of the compiled path (option `cacoa.cpp.pooled`, default TRUE; the R path gives the same result)
#' @return integer vector of length B + 1: the observed data set first
#' @keywords internal
scanCallCounts <- function(s.obs, sp, alpha = 0.1, tail = identity, standardize = FALSE, n.cores = 1) {
  if (isTRUE(getOption("cacoa.cpp.pooled", TRUE)) && !is.na(tc <- tailCode(tail)) && (!standardize || tc == 0L)) {   # built and counted in C++ without R copies
    if (!is.double(sp)) storage.mode(sp) <- "double"
    return(if (standardize) scan_call_counts_std_cpp(as.numeric(s.obs), sp, alpha, as.integer(max(1L, n.cores)))
           else scan_call_counts_raw_cpp(as.numeric(s.obs), sp, tc, alpha, as.integer(max(1L, n.cores)))) }
  B <- nrow(sp); S <- rbind(tail(s.obs), tail(sp)); bad <- !is.finite(S[1, ])
  if (standardize) { mu <- colMeans(S); sdv <- sqrt(pmax(colMeans(S^2) - mu^2, 0) * (B + 1) / B); bad <- bad | !is.finite(sdv) | sdv <= 0; S <- sweep(sweep(S, 2, mu), 2, sdv, "/") }
  S[, bad] <- -Inf; S[!is.finite(S)] <- -Inf
  if (isTRUE(getOption("cacoa.cpp.pooled", TRUE))) { storage.mode(S) <- "double"; return(scan_call_counts_cpp(S, alpha, as.integer(max(1L, n.cores)))) }   # the same counts, sorted once and counted on threads
  sv <- as.vector(S); oa <- order(sv, method = "radix"); ps <- sv[oa]; Nn <- length(ps)
  # number of pooled values strictly below each value (findInterval(x, ps, left.open = TRUE)), from the one sort of the pool:
  # a value's count is the position of the first copy of it in the sorted pool, minus one. findInterval() would check that
  # the pool is sorted at every call, once per relabeling, which dominated the run time on many cells.
  below <- integer(Nn); below[oa] <- match(ps, ps) - 1L; below <- matrix(below, nrow = nrow(S)); rm(sv, oa)
  vapply(seq_len(B + 1), function(r) { o <- order(S[r, ], decreasing = TRUE, method = "radix"); m <- S[r, o]; cnt <- below[r, o]
    fin <- is.finite(m); m <- m[fin]; cnt <- cnt[fin]; if (!length(m)) return(0L)
    R <- length(m) - findInterval(m, rev(m), left.open = TRUE)                       # own values at or above each value
    fdr <- (Nn - cnt - R) / B / R; k <- which(fdr < alpha); if (length(k)) as.integer(R[k[length(k)]]) else 0L }, integer(1))
}

#' Code of a tail function for the compiled pooled helpers (0 identity, 1 abs, 2 negation; NA for any other function)
#' @param tail function
#' @return integer code or NA
#' @keywords internal
tailCode <- function(tail) {
  x <- c(-2, -0.5, 0, 1.5, Inf, -Inf, NA); y <- suppressWarnings(tail(x))
  if (identical(y, x)) 0L else if (identical(y, abs(x))) 1L else if (identical(y, -x)) 2L else NA_integer_
}

#' Test of the automatically chosen scale of a scan
#'
#' `scale = "auto"` reports the scale that calls the most cells ([chooseScanScale()]), so it reports calls whenever any scale has them,
#' and each scale's own false-discovery rate does not cover that choice. Here the choice itself is tested: the number of cells called at
#' the chosen scale is compared with the same number for every relabeled data set (each put through all scales and the same rule).
#' @param n.obs number of cells called at each scale by the observed data
#' @param n.perm relabelings x scales matrix of the numbers called by the relabeled data sets (both from [scanCallCounts()])
#' @param tested number of tested cells
#' @param saturated share of the tested cells from which a scale counts as calling everything (as in [chooseScanScale()])
#' @return list: `p` (share of the relabeled data sets, the observed one included, whose chosen scale calls at least as many cells; 1 when
#'   the observed data call none), `n` (cells called at the chosen scale), `n.null` (the same per relabeling)
#' @keywords internal
autoScaleTest <- function(n.obs, n.perm, tested, saturated = 0.9) {
  pick <- function(n) { sat <- which(n >= saturated * tested); if (length(sat)) n[sat[1]] else max(n) }
  n <- pick(n.obs); nb <- apply(n.perm, 1, pick)
  list(p = if (n > 0) (1 + sum(nb >= n)) / (length(nb) + 1) else 1, n = n, n.null = nb)
}

#' Pooled permutation p-values of one scale
#'
#' The pooled permutation p of the observed per-cell statistics and of every permuted value (the share of permuted
#' values of all cells and relabelings at or above it, the value itself included) and the scale's own permutation
#' FDR ([pooledMinPFDR()]). Non-finite statistics rank last.
#' @param s.obs observed per-cell statistics
#' @param sp permuted statistics (relabelings x cells)
#' @param tail function mapping the statistic to its extreme tail (`abs` for two-sided, `identity` for greater)
#' @param n.cores threads of the compiled path (option `cacoa.cpp.pooled`, default TRUE; the R path gives the same result)
#' @return list with `p` (per cell), `pp` (relabelings x cells) and `q`
#' @keywords internal
pooledScanStep <- function(s.obs, sp, tail = abs, n.cores = 1) {
  B <- nrow(sp); N <- ncol(sp)
  if (isTRUE(getOption("cacoa.cpp.pooled", TRUE)) && !is.na(tc <- tailCode(tail))) {   # the same pooled p-values from one threaded sort and binary searches
    if (!is.double(sp)) storage.mode(sp) <- "double"
    r <- pooled_scan_step_cpp(as.numeric(s.obs), sp, tc, as.integer(max(1L, n.cores))); dim(r$pp) <- c(B, N)
    R <- findInterval(r$p, sort(r$p)); fdr <- pmin(1, (r$ge / B) / R); o <- order(r$p, decreasing = TRUE); q <- fdr; q[o] <- cummin(fdr[o])   # pooledMinPFDR(p, pp, B)$q without the second sort
    return(list(p = r$p, pp = r$pp, q = q))
  }
  a <- tail(as.vector(sp)); a[!is.finite(a)] <- -Inf; Nn <- length(a); so <- sort(a, method = "radix", index.return = TRUE); rm(a)
  first <- cummax(ifelse(c(TRUE, so$x[-1] != so$x[-Nn]), seq_len(Nn), 0L))      # first position of each value's tie block (-Inf ties included)
  pp <- numeric(Nn); pp[so$ix] <- (Nn - first + 2) / (Nn + 1); dim(pp) <- c(B, N)
  s <- tail(s.obs); s[!is.finite(s)] <- -Inf; p <- (Nn - findInterval(s, so$x, left.open = TRUE) + 1) / (Nn + 1)
  list(p = p, pp = pp, q = pooledMinPFDR(p, pp, B)$q)
}

#' Pooled permutation FDR of a per-cell statistic against its values under relabeling
#'
#' For a statistic whose small values are extreme (a p-value, or the smallest p over scales): the expected number of
#' cells at or below a cell's value under relabeling (its count among the permuted values over the number of
#' relabelings) over the number of cells observed at or below it, made monotone from the top.
#' @param m.obs observed per-cell values
#' @param m.perm permuted values (relabelings x cells, or any vector)
#' @param B number of relabelings
#' @param n.cores threads of the compiled path (option `cacoa.cpp.pooled`, default TRUE; the R path gives the same result)
#' @return list with `q` and `p` (the pooled permutation p of the value)
#' @keywords internal
pooledMinPFDR <- function(m.obs, m.perm, B, n.cores = 1) {
  if (isTRUE(getOption("cacoa.cpp.pooled", TRUE)) && !anyNA(m.obs)) {   # the same counts from a threaded sort (NA permuted values dropped, as sort() does)
    mp <- as.vector(m.perm); if (anyNA(mp)) mp <- mp[!is.na(mp)]
    r <- pooled_min_counts_cpp(as.numeric(m.obs), as.numeric(mp), as.integer(max(1L, n.cores))); rm(mp)
    fdr <- pmin(1, (r$cnt / B) / r$R); o <- order(m.obs, decreasing = TRUE); q <- fdr; q[o] <- cummin(fdr[o])
    return(list(q = q, p = (r$cnt + 1) / (r$Nn + 1)))
  }
  ps <- sort(as.vector(m.perm), method = "radix"); Nn <- length(ps)
  cnt <- findInterval(m.obs, ps); V <- cnt / B; R <- findInterval(m.obs, sort(m.obs))
  fdr <- pmin(1, V / R); o <- order(m.obs, decreasing = TRUE); q <- fdr; q[o] <- cummin(fdr[o])
  list(q = q, p = (cnt + 1) / (Nn + 1))
}


#' Extract contour from embedding
#'
#' @param emb cell embedding matrix
#' @param cell.groups vector of cell type annotation
#' @param group specify cell types for contour, multiple cell types are also supported
#' @param conf confidence interval of contour
#' @keywords internal
getDensityContour <- function(emb, cell.groups, group,  color='black', linetype=2, conf="10%", bandwidth=NULL, ...) {
  emb %<>% .[rownames(.) %in% names(cell.groups)[cell.groups %in% group], ]

  if (is.null(bandwidth)) {
    kd <- ks::kde(emb, compute.cont=TRUE, ...)
  } else {
    h <- matrix(c(bandwidth, 0, 0, bandwidth), ncol=2)
    kd <- ks::kde(emb, compute.cont=TRUE, H=h, ...)
  }

  lcn <- kd %$% contourLines(x=eval.points[[1]], y=eval.points[[2]], z=estimate, levels=cont[conf]) %>%
    .[[1]] %>% data.frame() %>% cbind(z=1)
  cn <- geom_path(aes(x, y), data=lcn, linetype=linetype, color=color);
  return(cn)
}


#' Plot cell density
#'
#' @param bins number of bins for density estimation, should keep consistent with bins in estimateCellDensity
#' @param palette color palette function. Default: `YlOrRd`
#' @keywords internal
plotDensityKde <- function(mat, bins, cell.emb, show.grid=TRUE, lims=NULL, show.labels=FALSE, show.ticks=FALSE,
                           palette=NULL, legend.title=NULL, ...) {
  if (is.null(lims)) {
    lims <- c(min(mat$z), max(mat$z)*1.1)
  }

  breaks <- lapply(mat[c('x', 'y')], function(m) {
    seq(quantile(m, 0.1), quantile(m, 0.9), length.out=6) %>%
      signif(digits=3)
  })

  p <- ggplot(mat, aes(x, y, fill = z)) +
    geom_raster() +
    scale_x_continuous(breaks=breaks$x, expand = c(0,0)) +
    scale_y_continuous(breaks=breaks$y, expand = c(0,0)) +
    val2ggcol(mat$z, palette=palette, color.range=lims, return.fill=TRUE)

  if (!is.null(legend.title)) p <- p + guides(fill=guide_colorbar(title=legend.title))

  p %<>% sccore::styleEmbeddingPlot(show.labels=show.labels, show.ticks=show.ticks, ...)

  if (show.grid) { #  add grid manually
    p <- p +
      geom_vline(xintercept=breaks$x, col='grey', alpha=0.1) +
      geom_hline(yintercept=breaks$y, col='grey', alpha=0.1)
  }

  return(p)
}

#' Estimate differential cell density
#'
#' @param density.mat estimated cell density matrix with estimateCellDensity
#' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
#' @param target.level target/disease level for sample.group vector
#' @param type method to calculate differential cell density of each bin; subtract: target density minus ref density; entropy: estimated kl divergence entropy between sample groups ; t.test: zscore of t-test,
#'   global variance is setting for t.test
#' @keywords internal
diffCellDensity <- function(density.mat, sample.groups, ref.level, target.level, type='subtract'){
  nt <- names(sample.groups[sample.groups == target.level]) # sample name of target
  nr <- names(sample.groups[sample.groups == ref.level]) # sample name of reference

  if (type == 'subtract') {
    rmt <- rowMeans(density.mat[, nt])
    rmr <- rowMeans(density.mat[, nr])
    score <- (rmt - rmr) / max(max(rmt), max(rmr))
  } else if (type=='t.test'){
    score <- matrixTests::row_t_welch(density.mat[,nt], density.mat[,nr])$statistic %>%
      setNames(rownames(density.mat))
  } else if (type == 'wilcox') {
    pvalue <- matrixTests::row_wilcoxon_twosample(density.mat[,nt], density.mat[,nr])$pvalue
    zstat <- abs(qnorm(pvalue / 2))
    fc <- rowMeans(density.mat[,nt]) - rowMeans(density.mat[,nr])
    score <- zstat * sign(fc)
  } else stop("Unknown method: ", type)

  return(score)
}

#' Estimate differential cell density with permutation testing
#'
#' @param density.mat estimated cell density matrix with estimateCellDensity
#' @param sample.meta sample metadata data.frame
#' @param contrast contrast matrix for the design
#' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
#' @param target.level target/disease level for sample.group vector
#' @param type method to calculate differential cell density of each bin;
#' @param block.id block variable for permutation test
#' Pooled p-values, tail z-scores and exact permutation p-values of per-cell statistics
#'
#' The unadjusted z-score of the density map: the signed tail z of the per-cell p-value against the pooled
#' permutation null of all cells (`pooledPermutationP`, no per-cell floor). Also recomputes the exact per-cell
#' permutation p, which is needed when the statistics have been post-processed (`averageOverNeighbours`).
#' @param stat.obs named vector of observed statistics
#' @param stats.perm permutations x cells matrix
#' @param effect named vector of effects (sign of the two-sided z)
#' @param alternative `"two-sided"`, `"greater"` or `"less"`
#' @return list with `p.pooled`, `z.score` and `p.perm`
#' @keywords internal
pooledPermutationSummaries <- function(stat.obs, stats.perm, effect, alternative = "two-sided") {
  cells <- names(stat.obs); sgn <- if (alternative == "less") -1 else 1
  p.pooled <- pooledPermutationP(sgn * stat.obs, sgn * stats.perm, alternative = if (alternative == "two-sided") "two.sided" else "greater")
  names(p.pooled) <- cells
  pp <- pmin(pmax(p.pooled, 1e-300), 1 - 1e-16)
  z.score <- switch(alternative, "two-sided" = sign(effect[cells]) * stats::qnorm(pp / 2, lower.tail = FALSE),
                    greater = stats::qnorm(pp, lower.tail = FALSE), less = stats::qnorm(pp))
  z.score[!is.finite(p.pooled)] <- NA_real_; names(z.score) <- cells
  so <- sgn * stat.obs; P <- sgn * stats.perm
  p.perm <- if (alternative == "two-sided") (colSums(abs(P) >= rep(abs(so), each = nrow(P)), na.rm = TRUE) + 1) / (nrow(P) + 1) else (colSums(P >= rep(so, each = nrow(P)), na.rm = TRUE) + 1) / (nrow(P) + 1)
  p.perm[!is.finite(so)] <- NA_real_; names(p.perm) <- cells
  list(p.pooled = p.pooled, z.score = z.score, p.perm = p.perm)
}

#' @keywords internal
diffCellDensityPermutations <- function(density.mat, sample.model, perm.method=c("block", "freedman-lane"), 
                                        n.permutations=999, robust.method=c("none","huber","winsor"), na.mode=c("drop", "impute_weak"),
                                        alternative=c("two-sided","greater","less"), return.residuals=FALSE, return.sampled.stats=TRUE,
                                        n.cores=1, verbose=TRUE, moderate=TRUE, pooled=FALSE, group=NULL, average=FALSE, nn.graph=NULL, closed.form=TRUE,
                                        pooled.statistic=c("raw", "standardized"), ...) {
  pooled.statistic <- match.arg(pooled.statistic); standardize <- pooled.statistic == "standardized"
  if (pooled) {   # the pooled histograms: of |t| (raw; t is studentized already) on 0..50 or of the standardized statistic on the z grid
    z.grid <- if (standardize) pooledZGrid(abs = alternative[1] == "two-sided") else if (alternative[1] == "two-sided") seq(0, 50, by = 0.01) else seq(-30, 50, by = 0.01)
  }
  perm.method <- match.arg(perm.method)
  robust.method <- match.arg(robust.method)
  na.mode <- match.arg(na.mode)
  alternative <- match.arg(alternative)
  two.sided <- alternative == "two-sided"
  if (average && is.null(nn.graph)) stop("average needs nn.graph (0-based neighbour lists, the cell itself included)")
  
    if (is.null(sample.model)) {
      stop("sample.model must be provided for testing.")
    }
    res.diff <- list()
    Y <- t(density.mat)
    if (is.null(rownames(Y))) stop("density.mat must have sample names as column names")
    
    ## ---- fit & permutations ----
    res <- performLMPermutations(x = sample.model, y = Y, n.permutations = n.permutations, perm.method = perm.method,
                                 robust.method = robust.method, na.mode = na.mode, alternative = alternative,
                                 return.sampled.stats = return.sampled.stats, return.residuals = return.residuals,
                                 n.cores = n.cores, pooled = pooled, z.grid = if (pooled) z.grid else NULL, use.abs = two.sided,
                                 group = group, average = average, nn.graph = nn.graph, closed.form = closed.form, pooled.standardize = standardize, ...)
  cells <- colnames(Y)
  # effect (contrast estimate) centred by its permutation mean; the z-scores are the signed tail z of the per-cell
  # p-values: of the pooled permutation null when the permuted statistics are kept (calibrated, no floor), otherwise
  # of the exact permutation p (floored at 1 / (n.permutations + 1))
  em <- if (!is.null(res$effect.perm.mean) && any(is.finite(res$effect.perm.mean))) res$effect.perm.mean else if (!is.null(res$effects.perm)) colMeans(res$effects.perm, na.rm = TRUE) else NULL
  score.c <- as.numeric(if (is.null(em)) res$effect else res$effect - ifelse(is.finite(em), em, 0))
  names(score.c) <- cells
  z.score <- res$z.score; names(z.score) <- cells
  effect <- stats::setNames(res$effect, cells); p.perm <- stats::setNames(res$pval, cells); stat.obs <- stats::setNames(res$stat.obs, cells)
  p.pooled <- NULL; z.std <- NULL; pooled.hist <- NULL; group.acc <- NULL
  if (!is.null(res$pooled)) {
    acc <- res$pooled
    stat.obs <- stats::setNames(as.numeric(acc$stat), cells); p.perm <- stats::setNames(as.numeric(acc$p), cells); z.std <- stats::setNames(as.numeric(acc$z), cells)
    if (average) { effect <- averageVectorOverLists(effect, nn.graph); score.c <- averageVectorOverLists(score.c, nn.graph); effect[!is.finite(stat.obs)] <- NA_real_ }
    base <- if (standardize) z.std else stat.obs   # the statistic the histograms were accumulated on
    s.obs <- if (two.sided) abs(base) else if (alternative == "greater") base else -base
    pooled.hist <- list(grid = as.numeric(acc$grid), hist = acc$hist, abs = two.sided, negate = alternative == "less", statistic = pooled.statistic, s.obs = s.obs)
    p.pooled <- pooledPermutationP(s.obs, pooled.hist, alternative = "greater"); names(p.pooled) <- cells
    pp <- pmin(pmax(p.pooled, 1e-300), 1 - 1e-16)
    z.score <- switch(alternative, "two-sided" = sign(effect) * stats::qnorm(pp / 2, lower.tail = FALSE), greater = stats::qnorm(pp, lower.tail = FALSE), less = stats::qnorm(pp))
    z.score[!is.finite(p.pooled)] <- NA_real_; names(z.score) <- cells
    group.acc <- acc[c("group_sum", "group_count", "group_obs", "group_n")]
  } else if (!is.null(res$stats.perm)) {
    ps <- pooledPermutationSummaries(stat.obs, res$stats.perm, effect, alternative)
    p.pooled <- ps$p.pooled; z.score <- ps$z.score
  }

  # moderated t: the per-cell residual variances are shrunk towards a common prior estimated across cells (limma's
  # empirical Bayes), which borrows strength across the many cells the way Milo's count model does across neighbourhoods
  dfr <- res$df; se2 <- res$se^2; ok <- is.finite(se2) & se2 > 0
  mod.perm <- list(...)$moderate.perm
  if (moderate && !(length(mod.perm) && mod.perm[1] != "none") && sum(ok) > 10) {   # (already moderated with moderate.perm)
    d0.max <- list(...)$moderate.d0.max %||% Inf
    pr <- fitFPrior(se2[ok], df1 = dfr, d0.max = d0.max); d0 <- if (is.finite(pr$df)) pr$df else 1e6
    dfo <- rep_len(dfr, length(se2))[ok]
    se.mod <- res$se; se.mod[ok] <- sqrt((d0 * pr$var + dfo * se2[ok]) / (d0 + dfo)); dfr <- dfr + d0
    tstat <- res$effect / se.mod
  } else tstat <- res$effect / res$se
  p.param <- if (alternative == "two-sided") 2 * stats::pt(-abs(tstat), dfr) else if (alternative == "greater") stats::pt(tstat, dfr, lower.tail = FALSE) else stats::pt(tstat, dfr)
  if (average) { tstat[] <- NA_real_; p.param[] <- NA_real_ }   # the parametric tail does not apply to the averaged statistic
  out <- list(score = score.c, permut.scores = if (is.null(res$effects.perm)) NULL else t(res$effects.perm), stats.perm = res$stats.perm, stat.obs = stat.obs, z.score = z.score,
              effect = effect, se = setNames(res$se, cells), statistic = res$statistic,
              t = setNames(tstat, cells), df = dfr, p.param = setNames(p.param, cells), p.perm = p.perm,
              p.pooled = p.pooled, z.std = z.std, pooled.hist = pooled.hist, group.acc = group.acc, max.z = if (!is.null(res$pooled)) as.numeric(res$pooled$max) else NULL, min.z = if (!is.null(res$pooled)) as.numeric(res$pooled$min) else NULL,
              averaged = average, closed.form = isTRUE(res$closed.form), pooled.statistic = if (pooled) pooled.statistic else NULL)

  if (return.residuals) {
    out$residuals <- res$residuals
    out$residuals.pearson <- res$residuals.pearson
  }
  return(out)
}

#' @keywords internal
adjustZScoresByPermutations <- function(score, scores.shuffled, wins=0.01, smooth=FALSE, graph=NULL, beta=30,
                                        l.max=NULL, n.cores=1, verbose=TRUE) {
  checkPackageInstalled("matrixStats", details="for adjusting p-values", cran=TRUE)
  if (smooth) {
    if (is.null(graph)) stop("graph has to be provided if smooth is TRUE")

    g.filt <- function(...) sccore::heatFilter(..., beta=beta)
    score %<>% smoothSignalOnGraph(filter=g.filt, graph=graph, l.max=l.max, progress=FALSE)
    scores.shuffled %<>%
      smoothSignalOnGraph(filter=g.filt, graph=graph, n.cores=n.cores, l.max=l.max, progress=verbose) %>%
      as.matrix()
  }

  if (wins > 0) {
    uq <- 1 - wins
    lq <- wins
    score %<>% pmin(quantile(., uq, na.rm=TRUE)) %>% pmax(quantile(., lq, na.rm=TRUE))
    scores.shuffled.ranges <- matrixStats::colQuantiles(scores.shuffled, probs=c(lq, uq), na.rm=TRUE)
  } else {
    scores.shuffled.ranges <- matrixStats::colRanges(scores.shuffled, na.rm=TRUE)
  }

  p.vals <- c()
  sign.mask <- (score >= 0)
  if (any(sign.mask)) {
    p.vals %<>% c((sapply(score[sign.mask], function(x) sum((x - 1e-10) < scores.shuffled.ranges[,2])) + 1) / (ncol(scores.shuffled) + 1))
  }

  if (any(!sign.mask)) {
    p.vals %<>% c((sapply(score[!sign.mask], function(x) sum((x + 1e-10) > scores.shuffled.ranges[,1])) + 1) / (ncol(scores.shuffled) + 1))
  }

  z.scores <- pmax(1 - p.vals[names(score)], 0.5) %>% qnorm(lower.tail=TRUE)
  z.scores[!sign.mask] %<>% {. * -1}

  return(z.scores)
}

#' @keywords internal
findScoreGroupsGraph <- function(scores, min.score, graph) { # TODO: deprecated?
  graph %>% igraph::induced_subgraph(names(scores)[scores > min.score]) %>%
    igraph::components() %>% .$membership
}

#' @keywords internal
graphFromGrid <- function(grid.ids, n.bins) {
  graph <- lapply(grid.ids, function (i) c(i+1, i-1, i+n.bins, i-n.bins)) %>%
    mapIds(grid.ids) %>% igraph::graph_from_adj_list() %>% igraph::as.undirected()
  igraph::V(graph)$name <- grid.ids
  return(graph)
}
