##' @description Estimate cell density in giving embedding, Density will estimated for indivisual sample
##' @param emb cell embedding matrix
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param bins number of bins for density estimation, default 400
estimateCellDensityKde <- function(emb, sample.per.cell, sample.groups, bins, bandwidth=0.05, expansion.mult=0.05){
  if (is.null(bandwidth)) {
    bandwidth <- apply(emb, 2, MASS::bandwidth.nrd)
  } else {
    bandwidth <- apply(emb, 2, quantile, c(0.1, 0.9)) %>% diff() %>% {. * bandwidth} %>% .[1,]
  }

  lims <- as.numeric(apply(emb,2,function(x) ggplot2:::expand_limits_continuous(range(x), expansion(mult=expansion.mult))))

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

  cond.densities <- split(names(sample.groups), sample.groups) %>%
    lapply(function(ns) rowMeans(density.mat[,ns]))

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


##' @description estimate graph smooth based cell density
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param n.cores number of cores
##' @param m numeric Maximum order of Chebyshev coeff to compute (default=50)
estimateCellDensityGraph <- function(graph, sample.per.cell, sample.groups, n.cores=1, beta=30, m=50, verbose = TRUE) {
  sig.mat <- unique(sample.per.cell) %>% sapply(function(s) as.numeric(sample.per.cell == s)) %>%
    set_rownames(names(sample.per.cell)) %>% set_colnames(unique(sample.per.cell))

  score.mat <- sccore:::smoothSignalOnGraph(
    sig.mat, filter=function(...) sccore:::heatFilter(..., beta=beta), graph=graph,
    m=m, n.cores=n.cores, progress.chunk=(verbose + 1), progress=verbose,
  ) %>% as.matrix()

  score.mat %<>% {t(.) / colSums(.)} %>% t() # Normalize by columns to adjust on the number of cells per sample


  cond.densities <- split(names(sample.groups), sample.groups) %>%
    sapply(function(samps) apply(score.mat[,samps,drop=FALSE], 1, mean, trim=0.2)) # Robust estimator of sum
  cond.densities %<>% {lapply(1:ncol(.), function(i) .[,i])} %>% setNames(colnames(cond.densities))

  return(list(density.mat=score.mat, cond.densities=cond.densities, method='graph'))
}


##' @description extract contour from embedding
##' @param emb cell embedding matrix
##' @param cell.groups vector of cell type annotation
##' @param group specify cell types for contour, multiple cell types are also supported
##' @param conf confidence interval of contour
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


##' @description Plot cell density
##' @param bins number of bins for density estimation, should keep consistent with bins in estimateCellDensity
##' @param palette color palette function. Default: `YlOrRd`
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

##' @description estimate differential cell density
##' @param density.mat estimated cell density matrix with estimateCellDensity
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @param type method to calculate differential cell density of each bin; subtract: target density minus ref density; entropy: estimated kl divergence entropy between sample groups ; t.test: zscore of t-test,
##' global variance is setting for t.test;
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

diffCellDensityPermutations <- function(density.mat, sample.groups, ref.level, target.level, type='permutation',
                                        verbose=TRUE, n.permutations=200, n.cores=1) {
  nt <- names(sample.groups[sample.groups == target.level]) # sample name of target
  nr <- names(sample.groups[sample.groups == ref.level]) # sample name of reference

  if (type %in% c('subtract', 't.test', 'wilcox')) {
    score <- diffCellDensity(density.mat, sample.groups=sample.groups, ref.level=ref.level, target.level=target.level, type=type)
    permut.scores <- plapply(1:n.permutations, function(i) {
      sg.shuff <- setNames(sample(sample.groups), names(sample.groups))
      diffCellDensity(density.mat, sample.groups=sg.shuff, ref.level=ref.level, target.level=target.level, type=type)
    }, progress=verbose, n.cores=n.cores, mc.preschedule=TRUE, fail.on.error=TRUE) %>% do.call(cbind, .)

    return(list(score=score, permut.scores=permut.scores))
  }

  if (type == 'permutation') {
    checkPackageInstalled("matrixStats", details="for `type='permutation'`", cran=TRUE)
    colRed <- matrixStats::colMedians
  } else if (type == 'permutation.mean') {
    colRed <- Matrix::colMeans
  } else stop("Unknown method: ", type)

  density.mat <- t(density.mat)
  dm.shuffled <- density.mat
  permut.diffs <- plapply(1:n.permutations, function(i) { # Null distribution looks normal, so we don't need a lot of samples
    rownames(dm.shuffled) %<>% sample()
    colRed(dm.shuffled[nt,]) - colRed(dm.shuffled[nr,])
  }, progress=verbose, n.cores=1, fail.on.error=TRUE) %>%  # according to tests, n.cores>1 here does not speed up calculations
    do.call(cbind, .) %>% set_rownames(colnames(density.mat))

  sds <- apply(permut.diffs, 1, sd)
  score <- (colRed(density.mat[nt,]) - colRed(density.mat[nr,])) / sds
  permut.diffs <- apply(permut.diffs, 2, `/`, sds)
  score[sds < 1e-10] <- 0
  permut.diffs[sds < 1e-10,] <- 0

  return(list(score=score, permut.scores=permut.diffs))
}

adjustZScoresByPermutations <- function(score, scores.shuffled, wins=0.01, smooth=FALSE, graph=NULL, beta=30, n.cores=1, verbose=TRUE) {
  checkPackageInstalled("matrixStats", details="for adjusting p-values", cran=TRUE)
  if (smooth) {
    if (is.null(graph)) stop("graph has to be provided if smooth is TRUE")

    g.filt <- function(...) heatFilter(..., beta=beta)
    score %<>% smoothSignalOnGraph(filter=g.filt, graph=graph)
    scores.shuffled %<>% smoothSignalOnGraph(filter=g.filt, graph=graph, n.cores=n.cores,
                                             progress=verbose) %>%
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
  sign_mask <- (score >= 0)
  if (any(sign_mask)) {
    p.vals %<>% c((sapply(score[sign_mask], function(x) sum((x - 1e-10) < scores.shuffled.ranges[,2])) + 1) / (ncol(scores.shuffled) + 1))
  }

  if (any(!sign_mask)) {
    p.vals %<>% c((sapply(score[!sign_mask], function(x) sum((x + 1e-10) > scores.shuffled.ranges[,1])) + 1) / (ncol(scores.shuffled) + 1))
  }

  z.scores <- pmax(1 - p.vals[names(score)], 0.5) %>% qnorm(lower.tail=TRUE)
  z.scores[!sign_mask] %<>% {. * -1}

  return(z.scores)
}

findScoreGroupsGraph <- function(scores, min.score, graph) {
  graph %>% igraph::induced_subgraph(names(scores)[scores > min.score]) %>%
    igraph::components() %>% .$membership
}


