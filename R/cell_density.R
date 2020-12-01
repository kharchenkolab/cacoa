
##' @description Estimate cell density in giving embedding
##' @param emb cell embedding matrix
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @param bins number of bins for density estimation, default 400
##' @param by.sample  if TRUE, density will estimated by sample and quantiles normalization will applied to individual sample. If FALSE, cell condition.per.cell need to be provided and density will simply esitmated by condition.per.cell.
estimateCellDensity <- function(emb, sample.per.cell, sample.groups, bins, ref.level, target.level, condition.per.cell = NULL, by.sample = TRUE){
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("You have to install preprocessCore package from Bioconductor to do quantile normlization ")
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("You have to install MASS package to estimate density ")
  }


  cname <- intersect(names(sample.per.cell), rownames(emb))
  sample.per.cell <- sample.per.cell[cname]
  condition.per.cell <- condition.per.cell[cname]
  emb <- emb[cname, ]
  list.den <- lapply(sccore:::sn(as.character(unique(sample.per.cell))), function(x) {
    nname <- names(sample.per.cell[sample.per.cell == x])
    tmp <- emb[nname, ]
    f2 <- MASS::kde2d(tmp[, 1], tmp[, 2], n = bins, lims = c(range(emb[, 1]), range(emb[, 2])))
    f2
  })
  den.mat <- do.call("cbind", lapply(list.den, function(x) as.numeric(x$z)))
  density.mat <- preprocessCore::normalize.quantiles(den.mat)    #quantiles normlization
  colnames(density.mat) <- colnames(den.mat)

  if (by.sample){
    density.fraction <- lapply(sccore:::sn(as.character(unique(sample.groups))),
                              function(x) {
                                tmp  <-  density.mat[, names(sample.groups[sample.groups == x])]
                                rowMeans(tmp)
                                #mmatrix(rowMeans(tmp), ncol = bins, byrow = FALSE)
                              })
  }else{
    if (is.null(condition.per.cell)) { stop("'condition.per.cell' must be provided") }
    list.den <- lapply(sccore:::sn(as.character(unique(condition.per.cell))), function(x) {
      nname <- names(condition.per.cell[condition.per.cell == x])
      tmp <- emb[nname, ]
      f2 <- kde2d(tmp[, 1], tmp[, 2], n = bins, lims = c(range(emb[, 1]), range(emb[, 2])))
      f2
    })
    denMatrix <- do.call("cbind", lapply(list.den, function(x) as.numeric(x$z)))
    density.fraction <- lapply(sccore:::sn(as.character(unique(sample.groups))),
                              function(x) {
                                denMatrix[, x]
                                #matrix(denMatrix[, x], ncol = bins, byrow = FALSE)
                              })
  }

  # coordinate embedding space
  target.density = density.fraction[[target.level]]
  mat <- matrix(target.density, ncol = bins, byrow = FALSE)
  x <- emb[, 1]
  y <- emb[, 2]
  x1=seq(min(x),max(x),length.out = bins)
  y1=seq(min(y),max(y),length.out = bins)
  names(x1)=seq(bins)
  names(y1)=seq(bins)
  d1=setNames(melt(mat), c('x', 'y', 'z'))
  d1$x1=x1[d1$x]
  d1$y1=y1[d1$y]
  emb2 <- data.frame(x = d1$x1, y = d1$y1)

  #count cell number in each bin
  x <- emb[,1]
  y <- emb[,2]
  s1 <- seq(from = min(x),
           to = max(x),
           length.out = bins + 1)
  s2 <- seq(from = min(y),
           to = max(y),
           length.out = bins + 1)
  dcounts <- table(cut(x, breaks = s1), cut(y, breaks = s2)) #%>% as.matrix.data.frame
  emb2$counts <- as.numeric(dcounts)

  return(list('density.mat' = density.mat, 'density.fraction' = density.fraction, 'density.emb' = emb2))
}


##' @description estimate graph  smooth based cell density
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @param n.cores number of cores
##' @param m numeric Maximum order of Chebyshev coeff to compute (default=50)

estimateGraphDensity <- function(sample.per.cell, sample.groups, ref.level, target.level, n.cores = 1, m = 50, verbose = TRUE) {
  tmp <-  setNames(as.numeric(sample.per.cell), names(sample.per.cell))
  scoreL <- sccore:::plapply(sn(unique(tmp)), function(x) {
    tryCatch({
      x1 <-  tmp
      x1[x1 != x] <-  0
      x1[x1 == x] <-  1
      sccore:::smoothSignalOnGraph(x1, con$graph, sccore:::heatFilter, m = m)
    }, error = function(err) {
      return(NA)
    })
  }, n.cores = n.cores, mc.preschedule=T, progress=verbose)
  scM <- do.call(cbind, scoreL)
  colnames(scM) <-  unique(sample.per.cell)
  NT <-  sample.groups[sample.groups == target.level] %>% names()
  NR <- sample.groups[sample.groups == ref.level] %>% names()
  score <- apply(scM, 1, function(x) {
    x1 <- x[NT]
    x2 <- x[NR]
    t.test(x1, x2)$statistic
  })
  names(score) <- rownames(scM)
  # trim outliner
  c1 <- quantile(score, 0.9999)
  c2 <- quantile(score, 0.0001)
  score[score > c1] <- c1
  score[score < c2] <- c2
  return(score)
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
##' @param bins number of bins for density esitmation, should keep consistent with bins in estimateCellDensity
##' @param col color palettes, default is c('blue','white','red')
plotDensity <- function(mat, bins, col = c('blue','white','red'), show.legend = NULL, legend.position = NULL, title = NULL, show.grid = NULL, mi=NULL, ma=NULL, diffDensity = NULL){
  #  p  <-  mat %>% as_tibble() %>% rowid_to_column(var = "X") %>%
  #    gather(key = "Y", value = "Z", -1) %>% mutate(Y = as.numeric(gsub("V", "", Y))) %>%
  #
  if (is.null(mi)){
    mi <- min(mat$z)
  }
  if (is.null(ma)){
    ma <- max(mat$z)*1.1
  }

  if (is.null(diffDensity)){
    p <- ggplot(mat, aes(x, y, fill = z)) +
      geom_raster() +
      theme_bw() + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), panel.border = element_blank(),
                         panel.background = element_blank(), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.title.y = element_blank(), axis.text.y = element_blank()) +
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
      viridis::scale_fill_viridis(option = 'B', alpha = 1, direction = 1, limits = c(mi, ma))
      if(!is.null(show.grid)){ #  add grid manually
        p <- p + geom_vline(xintercept=seq(quantile(mat$x,0.1),quantile(mat$x,0.9), length.out=6), col='grey', alpha=0.1)
        p <- p + geom_hline(yintercept=seq(quantile(mat$y,0.1),quantile(mat$y,0.9),, length.out=6), col='grey', alpha=0.1)
      }
  }else{ # using geom_tile and keep the same
    p <- ggplot(mat, aes(x, y, fill = z)) +
      geom_tile() +
      theme_bw() +
      ggplot2::lims(x = range(mat[, 'x']), y = range(mat[, 'y']))+
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks = element_blank())+
            scale_fill_gradient2(low = col[1], high = col[3], mid = col[2], midpoint = 0, limits = c(mi, ma))

  }


  p <- p + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


  if (is.null(show.legend)){
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






##' @description esitmate differential cell density
##' @param density.mat esitmated cell density matrix with estimateCellDensity
##' @param bins number of bins for density esitmation, should keep consistent with bins in estimateCellDensity
##' @param col color palettes, default is c('blue','white','red')
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param condition.per.cell A two-level factor on the cell names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @param method method to cacuated differential cell density of each bin; substract: target density minus ref density; entropy: estimated kl divergence entropy betwwen sample grapups ; t.test: zscore of t-test,global variacen is setting for t.test;
diffCellDensity <- function(density.emb, density.mat, condition.per.cell, sample.groups, bins, ref.level, target.level, method = 'substract', show.legend = NULL,legend.position = NULL, show.grid = TRUE, col = c('blue','white','red'), title = NULL, dcount.cutoff = 0, z.cutoff = NULL){
  nt <- names(sample.groups[sample.groups == target.level]) # sample name of target
  nr <- names(sample.groups[sample.groups == ref.level]) # sample name of reference

  if (method == 'substract'){
    score = rowMeans(density.mat[, nt]) - rowMeans(density.mat[, nr])
  }else if (method == 'entropy'){
    sudo <- mean(as.numeric(density.mat)) # add sudo counts
    density.mat2 <- density.mat + sudo
    s1 <- rowSums(density.mat2[, nr])
    s2 <- rowSums(density.mat2[, nt])
    r1 <- s1 / (s1 + s2)
    r2 <- s2 / (s1 + s2)
    weight.sum.per.fac.cell <- data.frame(r1, r2)
    xt <- table(condition.per.cell)
    max.ent <- (if (xt[1] > xt[2]) c(0, 1) else c(1, 0)) %>% entropy::KL.empirical(xt, unit='log2')
    entropy.per.cell <- apply(weight.sum.per.fac.cell, 1, entropy::KL.empirical, xt, unit = 'log2') / max.ent
    score <- entropy.per.cell * sign(r2 - r1)
  }else if (method=='t.test'){
    vel <- rowMeans(density.mat)
    density.mat2 <- density.mat + quantile(vel, 0.05) # add sudo counts at 5%
    score <- apply(density.mat2, 1, function(x) {
      x1 <- x[nt]
      x2 <- x[nr]
      tryCatch({
        t.test(x1, x2)$statistic
      }, error = function(e) {
        0
      })
    })
  } else if (method == 'willcox') {
    vel <- rowMeans(density.mat)
    density.mat2 <- density.mat + quantile(vel, 0.05) # add sudo counts at 5%
    score <- apply(density.mat2, 1, function(x) {
      x1 <- x[nt]
      x2 <- x[nr]
      mw = wilcox.test(x1, x2, exact = FALSE)
      zstat <- abs(qnorm(mw$p.value / 2))
      fc <- mean(x1) - mean(x2)
      zscore <- zstat * sign(fc)
      zscore
    })
  }

  if (is.null(title)){
    title <- method
  }
  #density.score <- matrix(score, ncol = bins, byrow = FALSE)
  #density.score[dcounts < dcount.cutoff] <- 0
  mat <-  data.frame(density.emb, 'z' = score)
  mat <-  mat[mat$counts > dcount.cutoff, ]

  if (!is.null(z.cutoff))
    mat[abs(mat$z) < z.cutoff, 'z'] = 0

  p <- plotDensity(mat, bins, col = col, title = title, legend.position = legend.position, show.legend = show.legend, show.grid = show.grid, diffDensity = TRUE)

  return(list('fig' = p,'score' = mat))
}



