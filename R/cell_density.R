#' @importFrom sccore checkPackageInstalled
NULL

#' Estimate cell density in giving embedding, density will be estimated for individual samples
#'
#' @param emb cell embedding matrix
#' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
#' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param bins number of bins for density estimation, default 400
#' @keywords internal
estimateCellDensityKde <- function(emb, sample.per.cell, sample.groups, bins, bandwidth=0.05, expansion.mult=0.05) {
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


#' Estimate graph smooth based cell density
#'
#' @param graph input graph
#' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
#' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param n.cores number of cores
#' @param m numeric Maximum order of Chebyshev coeff to compute (default=50)
#' @keywords internal
estimateCellDensityGraph <- function(graph, sample.per.cell, sample.groups, n.cores=1, beta=30, m=50, verbose=TRUE) {
  sig.mat <- unique(sample.per.cell) %>% sapply(function(s) as.numeric(sample.per.cell == s)) %>%
    set_rownames(names(sample.per.cell)) %>% set_colnames(unique(sample.per.cell))

  score.mat <- sccore::smoothSignalOnGraph(
    sig.mat, filter=function(...) sccore::heatFilter(..., beta=beta), graph=graph,
    m=m, n.cores=n.cores, progress.chunk=(verbose + 1), progress=verbose,
  ) %>% as.matrix()

  score.mat %<>% {t(.) / colSums(.)} %>% t() # Normalize by columns to adjust on the number of cells per sample


  cond.densities <- split(names(sample.groups), sample.groups) %>%
    sapply(function(samps) apply(score.mat[,samps,drop=FALSE], 1, mean, trim=0.2)) # Robust estimator of sum
  cond.densities %<>% {lapply(1:ncol(.), function(i) .[,i])} %>% setNames(colnames(cond.densities))

  return(list(density.mat=score.mat, cond.densities=cond.densities, method='graph'))
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
#' @param X model matrix used for fitting (preferably with intercept kept)
#' @keywords internal
diffCellDensityPermutations <- function(density.mat, X, sample.groups, ref.level, target.level, type='permutation',
                                        verbose=TRUE, smooth=FALSE, graph=NULL, l.max=NULL, beta=30, n.permutations=200, n.cores=1) {
  if (type %in% c('subtract', 't.test', 'wilcox')) { # no covariate implementation
    nt <- names(sample.groups[sample.groups == target.level]) 
   nr <- names(sample.groups[sample.groups == ref.level]) 
    score <- diffCellDensity(
      density.mat, sample.groups=sample.groups, ref.level=ref.level, target.level=target.level, type=type
    )
    permut.scores <- plapply(1:n.permutations, function(i) {
      sg.shuff <- setNames(sample(sample.groups), names(sample.groups))
      diffCellDensity(density.mat, sample.groups=sg.shuff, ref.level=ref.level, target.level=target.level, type=type)
    }, progress=verbose, n.cores=n.cores, mc.preschedule=TRUE, fail.on.error=TRUE) %>% do.call(cbind, .)

    return(list(score=score, permut.scores=permut.scores))
  }

  # linear model
  res.diff <- fit_density_lm(X, t(density.mat), n.permutations)
    rownames(res.diff$KP) <- colnames(X)
    rownames(res.diff$KPge) <- colnames(X)
    dimnames(res.diff$KP_perm)[[1]] <- colnames(X) # coef, bins, permutations

    intercept.idx <- which(rownames(res.diff$KP) == "(Intercept)") # remove intercept from output
    if (length(intercept.idx) > 0) {
        res.diff$KP <- res.diff$KP[-intercept.idx, , drop = FALSE]
        res.diff$KPge <- res.diff$KPge[-intercept.idx, , drop = FALSE]
        res.diff$KP_perm <- res.diff$KP_perm[-intercept.idx, , , drop = FALSE]
    }

    colnames(res.diff$KP) <- rownames(density.mat)
    colnames(res.diff$KPge) <- rownames(density.mat)
    dimnames(res.diff$KP_perm)[[2]] <- rownames(density.mat)

    coef.names <- rownames(res.diff$KP)
    res.diff$Z <- qnorm((res.diff$KPge + 1) / (res.diff$n_randomizations + 1))

    n.coefs <- length(coef.names)
    n.bins <- ncol(res.diff$KP)
    adjusted.scores <- matrix(NA_real_, nrow = n.coefs, ncol = n.bins,
                                                        dimnames = list(coef.names, colnames(res.diff$KP)))
    score.c <- matrix(NA_real_, nrow = n.coefs, ncol = n.bins,
                                        dimnames = list(coef.names, colnames(res.diff$KP)))
    for (i in seq_len(n.coefs)) {
        coef <- coef.names[i]
        score <- res.diff$KP[coef, ]
        scores.shuffled <- aperm(res.diff$KP_perm[coef, , , drop = FALSE], c(2, 3, 1))[, , 1, drop = FALSE]
        scores.shuffled <- matrix(scores.shuffled, nrow = n.bins, ncol = dim(res.diff$KP_perm)[3], byrow = FALSE)
        rownames(scores.shuffled) <- colnames(res.diff$KP)
        score.c[i, ] <- score - rowMeans(scores.shuffled, na.rm = TRUE) 
        adjusted.scores[i, ] <- adjustZScoresByPermutations(
            score = score.c[i, ],
            scores.shuffled = scores.shuffled,
            smooth = smooth,
            graph = graph,
            n.cores = n.cores,
            verbose = verbose,
            l.max = l.max,
            beta = beta
        )
    }
    res.diff$Z_adj <- adjusted.scores
    res.diff$score <- score.c
    return(res.diff)
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
