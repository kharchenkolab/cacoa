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
#'   neighbourhood whose width is set by `beta` or, better, by `target.size` cells); `"neighbourhood"`: the counts
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
    if (verbose) message(sprintf("kNN neighbourhoods of at least %d cells (%.1f per sample)", k, k / ncol(sig.mat)))
    nn <- neighbourhoods %||% knnNeighbourhoodLists(graph, k, cells = cells, n.cores = n.cores)
    if (!identical(names(nn), cells)) stop("the neighbourhood lists do not match the graph's cells")
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
    } else {
      unit.size <- kernelEffectiveSize(graph, beta, n.cores = n.cores)$size
      if (verbose) message(sprintf("heat kernel beta = %.3g: median effective neighbourhood size %.0f cells (5-95 %%: %.0f-%.0f)", beta, stats::median(unit.size), stats::quantile(unit.size, 0.05), stats::quantile(unit.size, 0.95)))
    }
    score.mat <- sccore::smoothSignalOnGraph(
      sig.mat, filter=function(...) sccore::heatFilter(..., beta=beta), graph=graph,
      m=m, n.cores=n.cores, progress.chunk=(verbose + 1), progress=verbose,
    ) %>% as.matrix()

    score.mat %<>% {t(.) / colSums(.)} %>% t() # Normalize by columns to adjust on the number of cells per sample
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
  if (moderate && sum(ok) > 10 && requireNamespace("limma", quietly = TRUE)) {
    sq <- limma::squeezeVar(se2[ok], df = dfr)
    se.mod <- res$se; se.mod[ok] <- sqrt(sq$var.post); dfr <- dfr + sq$df.prior
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
