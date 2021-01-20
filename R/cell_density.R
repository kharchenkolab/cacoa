
##' @description Estimate cell density in giving embedding, Density will estimated for indivisual sample
##' @param emb cell embedding matrix
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param bins number of bins for density estimation, default 400
estimateCellDensityKde <- function(emb, sample.per.cell, sample.groups, bins, expansion.mult=0.05){
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("You have to install preprocessCore package from Bioconductor to do quantile normlization ")
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("You have to install MASS package to estimate density ")
  }
  lims <- as.numeric(apply(emb,2,function(x) ggplot2:::expand_limits_continuous(range(x),expansion(mult=expansion.mult))))

  cname <- intersect(names(sample.per.cell), rownames(emb))
  sample.per.cell <- sample.per.cell[cname]
  emb <- emb[cname, ]
  cells.per.samp <- split(names(sample.per.cell), sample.per.cell)
  density.mat <-lapply(cells.per.samp, function(nn) {
      MASS::kde2d(emb[nn, 1], emb[nn, 2], n = bins, lims = lims)$z %>% as.numeric()
    }) %>% do.call("cbind", .) %>%
    preprocessCore::normalize.quantiles() %>%
    set_colnames(names(cells.per.samp)) %>%
    set_rownames(1:nrow(.)) # needed for indexing in diffCellDensity

  density.fraction <- split(names(sample.groups), sample.groups) %>%
    lapply(function(ns) rowMeans(density.mat[,ns]))

  # coordinate embedding space
  mat <- matrix(density.fraction[[1]], ncol = bins, byrow = FALSE)
  x1 <- seq(lims[1], lims[2], length.out=bins) %>% setNames(seq(bins))
  y1 <- seq(lims[3], lims[4], length.out=bins) %>% setNames(seq(bins))
  d1 <- setNames(melt(mat), c('x', 'y', 'z'))
  emb2 <- data.frame(x=x1[d1$x], y=y1[d1$y])

  #count cell number in each bin
  s1 <- seq(from = lims[1], to = lims[2], length.out=bins + 1)
  s2 <- seq(from = lims[3], to = lims[4], length.out=bins + 1)
  dcounts <- table(cut(emb[,1], breaks = s1), cut(emb[,2], breaks = s2)) #%>% as.matrix.data.frame
  emb2$counts <- as.numeric(dcounts)

  return(list(density.mat=density.mat, density.fraction=density.fraction, density.emb=emb2, bins=bins,
              method='kde', cell.emb=emb))
}


##' @description estimate graph smooth based cell density
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param n.cores number of cores
##' @param m numeric Maximum order of Chebyshev coeff to compute (default=50)
estimateCellDensityGraph <- function(graph, sample.per.cell, sample.groups, n.cores=1, beta=30, m=50, verbose = TRUE) {
  tmp <-  setNames(as.numeric(sample.per.cell), names(sample.per.cell))
  scores.smoothed <- sccore:::plapply(sccore::sn(unique(tmp)), function(x) {
    tryCatch({
      x1 <-  tmp
      x1[x1 != x] <-  0
      x1[x1 == x] <-  1
      sccore:::smoothSignalOnGraph(x1, graph, function(...) sccore:::heatFilter(..., beta=beta), m = m)
    }, error = function(err) {
      return(NA)
    })
  }, n.cores = n.cores, mc.preschedule=TRUE, progress=verbose)
  score.mat <- do.call(cbind, scores.smoothed)
  colnames(score.mat) <-  unique(sample.per.cell)

  density.fraction <- lapply(sccore:::sn(as.character(unique(sample.groups))),
                             function(x) rowMeans(score.mat[, names(sample.groups[sample.groups == x])]))

  return(list(density.mat=score.mat, density.fraction=density.fraction, method='graph'))
}


##' @description extract contour from embedding
##' @param emb cell embedding matrix
##' @param cell.groups vector of cell type annotation
##' @param group specify cell types for contour, multiple cell types are also supported
##' @param conf confidence interval of contour
getDensityContour <- function(emb, cell.groups, group,  color='white', linetype = 2, conf = "10%"){
  emb %<>% .[rownames(.) %in% names(cell.groups)[cell.groups %in% group], ]
  kd <- ks::kde(emb, compute.cont=TRUE)
  lcn <- kd %$% contourLines(x=eval.points[[1]], y=eval.points[[2]], z=estimate, levels=cont[conf]) %>%
    .[[1]] %>% data.frame() %>% cbind(z=1)
  cn <- geom_path(aes(x, y), data=lcn, linetype=linetype , color=color);
  return(cn)
}


##' @description Plot cell density
##' @param bins number of bins for density estimation, should keep consistent with bins in estimateCellDensity
##' @param palette color palette function. Default: `viridis::viridis_pal(option="B")`
plotDensityKde <- function(mat, bins, cell.emb, show.grid=TRUE, lims=NULL, show.labels=FALSE, show.ticks=FALSE, palette=viridis::viridis_pal(option="B"), ...){
  if (is.null(lims)){
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

  p %<>% sccore::styleEmbeddingPlot(show.labels=show.labels, show.ticks=show.ticks, ...)

  if(show.grid){ #  add grid manually
    p <- p +
      geom_vline(xintercept=breaks$x, col='grey', alpha=0.1) +
      geom_hline(yintercept=breaks$y, col='grey', alpha=0.1)
  }

  return(p)
}

adjustPvalueScores <- function(scores) {
  scores %<>% abs() %>% pnorm(lower.tail=FALSE) %>% p.adjust(method='BH') %>%
    qnorm(lower.tail=FALSE) %>% {. * sign(scores)}
  return(scores)
}

##' @description estimate differential cell density
##' @param density.mat estimated cell density matrix with estimateCellDensity
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @param type method to calculate differential cell density of each bin; subtract: target density minus ref density; entropy: estimated kl divergence entropy between sample groups ; t.test: zscore of t-test,
##' global variance is setting for t.test;
diffCellDensity <- function(density.mat, sample.groups, ref.level, target.level, type = 'subtract',
                            z.cutoff = NULL, adjust.pvalues=TRUE){
  nt <- names(sample.groups[sample.groups == target.level]) # sample name of target
  nr <- names(sample.groups[sample.groups == ref.level]) # sample name of reference

  if (type == 'subtract') {
    score <- rowMeans(density.mat[, nt]) - rowMeans(density.mat[, nr])
  #} else if (type == 'subtract.norm'){
  #  score <- (rowMeans(density.mat[, nt]) - rowMeans(density.mat[, nr])) / rowMeans(density.mat[, nr])
  } else if (type=='t.test'){
    score <- matrixTests::row_t_welch(density.mat[,nt], density.mat[,nr])$statistic %>%
      setNames(rownames(density.mat))
    if(adjust.pvalues) score %<>% adjustPvalueScores()
  } else if (type == 'wilcox') {
    pvalue <- matrixTests::row_wilcoxon_twosample(density.mat[,nt], density.mat[,nr])$pvalue
    zstat <- abs(qnorm(pvalue / 2))
    fc <- rowMeans(density.mat[,nt]) - rowMeans(density.mat[,nr])
    score <- zstat * sign(fc)
    if(adjust.pvalues) score %<>% adjustPvalueScores()
  } else stop("Unknown method: ", type)

  if (!is.null(z.cutoff)) {
    score[abs(score) < z.cutoff] <- 0
  }

  return(score)
}



