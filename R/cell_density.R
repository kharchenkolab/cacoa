
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
  list.den <- lapply(sccore:::sn(as.character(unique(sample.per.cell))), function(x) {
    nname <- names(sample.per.cell[sample.per.cell == x])
    tmp <- emb[nname, ]
    f2 <- MASS::kde2d(tmp[, 1], tmp[, 2], n = bins, lims = lims)
    f2
  })
  den.mat <- do.call("cbind", lapply(list.den, function(x) as.numeric(x$z)))
  density.mat <- preprocessCore::normalize.quantiles(den.mat)    #quantile normalization
  colnames(density.mat) <- colnames(den.mat)

  density.fraction <- lapply(sccore:::sn(as.character(unique(sample.groups))),
                             function(x) rowMeans(density.mat[, names(sample.groups[sample.groups == x])]))
  
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

  return(list(density.mat=density.mat,density.fraction=density.fraction, density.emb=emb2, bins=bins,'method'='kde'))
}


##' @description estimate graph smooth based cell density
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param n.cores number of cores
##' @param m numeric Maximum order of Chebyshev coeff to compute (default=50)
estimateCellDensityGraph <- function(graph, sample.per.cell, sample.groups, n.cores = 1, m = 50, verbose = TRUE) {
  tmp <-  setNames(as.numeric(sample.per.cell), names(sample.per.cell))
  scoreL <- sccore:::plapply(sccore::sn(unique(tmp)), function(x) {
    tryCatch({
      x1 <-  tmp
      x1[x1 != x] <-  0
      x1[x1 == x] <-  1
      sccore:::smoothSignalOnGraph(x1, graph, sccore:::heatFilter, m = m)
    }, error = function(err) {
      return(NA)
    })
  }, n.cores = n.cores, mc.preschedule=TRUE, progress=verbose)
  scM <- do.call(cbind, scoreL)
  colnames(scM) <-  unique(sample.per.cell)
 
  density.fraction <- lapply(sccore:::sn(as.character(unique(sample.groups))),
                             function(x) rowMeans(scM[, names(sample.groups[sample.groups == x])]))
  
  return(list(density.mat=scM,density.fraction=density.fraction,'method'='graph'))
}



##' @description extract contour from embedding
##' @param emb cell embedding matrix
##' @param cell.type vector of cell type annotation
##' @param cell specify cell types for contour, mutiple cell types are also suported
##' @param conf confidence interval of contour
getContour <- function(emb, cell.type, cell,  color = 'white', linetype = 2, conf = "10%"){
  linetype <- 2
  tmp <- emb[rownames(emb) %in% names(cell.type)[cell.type %in% cell], ]
  kd <- ks::kde(tmp, compute.cont = TRUE)
  lcn <- with(kd, contourLines(x = eval.points[[1]], y = eval.points[[2]], z = estimate, levels = cont[conf])[[1]])
  #name1 <- point.in.polygon(tmp[,1], tmp[,2], cn$x, cn$y)
  dd <- data.frame(lcn)
  dd$z <- 1
  cn <- geom_path(aes(x, y), data = dd, linetype = linetype , color = color);
  return(cn)
}


##' @description Plot cell density
##' @param bins number of bins for density estimation, should keep consistent with bins in estimateCellDensity
##' @param col color palettes, default is c('blue','white','red')
plotDensity <- function(mat, bins, show.legend = FALSE, legend.position = NULL, title = NULL, show.grid = TRUE, mi=NULL, ma=NULL, method = NULL){
  if (is.null(mi)){
    mi <- min(mat$z)
  }
  if (is.null(ma)){
    ma <- max(mat$z)*1.1
  }
  
  p <- ggplot(mat, aes(x, y, fill = z)) +
    geom_raster() +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), panel.border = element_blank(),
                       panel.background = element_blank())+ #, plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank()) +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    viridis::scale_fill_viridis(option = 'B', alpha = 1, direction = 1, limits = c(mi, ma))
    if(show.grid){ #  add grid manually
      p <- p + geom_vline(xintercept=seq(quantile(mat$x,0.1),quantile(mat$x,0.9), length.out=6), col='grey', alpha=0.1)
      p <- p + geom_hline(yintercept=seq(quantile(mat$y,0.1),quantile(mat$y,0.9),, length.out=6), col='grey', alpha=0.1)
    }
    if (!show.legend){
      p <- p + theme(legend.position = "none")
    }
    if (!is.null(legend.position)){
      p <- p + theme(legend.position = legend.position)
    }
    if(!is.null(title)){
      p <- p + ggtitle(title)
    }
  return(p)
}






##' @description estimate differential cell density
##' @param density.mat estimated cell density matrix with estimateCellDensity
##' @param bins number of bins for density estimation, should keep consistent with bins in estimateCellDensity
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @param type method to calculate differential cell density of each bin; subtract: target density minus ref density; entropy: estimated kl divergence entropy between sample groups ; t.test: zscore of t-test,
##' global variance is setting for t.test;
diffCellDensity <- function(density.emb, density.mat, sample.groups, bins, ref.level, target.level, type = 'subtract',
                             z.cutoff = NULL, adjust.pvalues=TRUE){
  nt <- names(sample.groups[sample.groups == target.level]) # sample name of target
  nr <- names(sample.groups[sample.groups == ref.level]) # sample name of reference

  if (type == 'subtract') {
    score <- rowMeans(density.mat[, nt]) - rowMeans(density.mat[, nr])
  #} else if (type == 'subtract.norm'){
  #  score <- (rowMeans(density.mat[, nt]) - rowMeans(density.mat[, nr])) / rowMeans(density.mat[, nr])
  } else if (type=='t.test'){
    score <- matrixTests::row_t_welch(density.mat[,nt], density.mat[,nr])$statistic 
    if(adjust.pvalues) score <- sign(score) * qnorm(p.adjust(pnorm(abs(score),lower.tail=F),method='BH'),lower.tail=F)
  } else if (type == 'wilcox') {
    pvalue = matrixTests::row_wilcoxon_twosample(density.mat[,nt], density.mat[,nr])$pvalue
    zstat <- abs(qnorm(pvalue / 2))
    fc <- rowMeans(density.mat[,nt]) - rowMeans(density.mat[,nr])
    score <- zstat * sign(fc)
    if(adjust.pvalues) score <- sign(score) * qnorm(p.adjust(pnorm(abs(score),lower.tail=F),method='BH'),lower.tail=F)
  } else stop("Unknown method: ", type)

  mat <-  data.frame(density.emb, 'z' = score)
  if (!is.null(z.cutoff)) {
    mat[abs(mat$z) < z.cutoff, 'z'] = 0
  }
  return(mat)
}



