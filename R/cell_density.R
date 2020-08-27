
##' @description Estimate cell density in giving embedding 
##' @param emb cell embedding matrix
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @param bins number of bins for density esitmation, default 400
##' @param by.sample  if TRUE, density will esitmated by sample and quantiles normlization will applied to indivisual sample. If FALSE, cell fraction need to be provided and density will simply esitmated by fraction. 
##' @add.ponits add.ponits  show cells in density plot   
estimateCellDensity <- function(emb, anoSample, sample.groups, bins, ref.level, target.level, fraction = NULL, by.sample = TRUE){
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("You have to install preprocessCore package to do quantile normlization ")
  }
  
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("You have to install MASS package to estimate density ")
  }
  
  
  cname <- intersect(names(anoSample), rownames(emb)) 
  anoSample <- anoSample[cname]
  fraction <- fraction[cname]
  emb <- emb[cname, ]
  list.den <- lapply(sccore:::sn(as.character(unique(anoSample))), function(x) {
    nname <- names(anoSample[anoSample == x])
    tmp <- emb[nname, ]
    f2 <- MASS::kde2d(tmp[, 1], tmp[, 2], n = bins, lims = c(range(emb[, 1]), range(emb[, 2])))
    f2
  })
  denMatrix <- do.call("cbind", lapply(list.den, function(x) as.numeric(x$z)))
  denMatrix.nor <- preprocessCore::normalize.quantiles(denMatrix)    #quantiles normlization
  colnames(denMatrix.nor) <- colnames(denMatrix)
  
  if (by.sample){
    density.fraction <- lapply(sccore:::sn(as.character(unique(sample.groups))), 
                              function(x) {
                                tmp  <-  denMatrix.nor[, names(sample.groups[sample.groups == x])]
                                matrix(rowMeans(tmp), ncol = bins, byrow = FALSE)
                              })
  }else{
    if (is.null(fraction)) { stop("'fraction' must be provided") }
    list.den <- lapply(sccore:::sn(as.character(unique(fraction))), function(x) {
      nname <- names(fraction[fraction == x])
      tmp <- emb[nname, ]
      f2 <- kde2d(tmp[, 1], tmp[, 2], n = bins, lims = c(range(emb[, 1]), range(emb[, 2])))
      f2
    })
    denMatrix <- do.call("cbind", lapply(list.den, function(x) as.numeric(x$z)))
    density.fraction <- lapply(sccore:::sn(as.character(unique(sample.groups))), 
                              function(x) {
                                matrix(denMatrix[, x], ncol = bins, byrow = FALSE)
                              })
    }
  return(list('denMatrix.nor' = denMatrix.nor, 'density.fraction' = density.fraction))
}





##' @description extract Counter from embedding 
##' @param emb cell embedding matrix
##' @param cell.type specify cell types for counter, mutiple cell types are also suported 
##' @param conf confidence interval of counter
##' @param bins number of bins for density esitmation, should keep consistent with bins in estimateCellDensity
getContour <- function(emb, cell.type, bins, cell, color = 'white', linetype = 2, conf = "10%"){
  x <- emb[, 1]
  y <- emb[, 2]
  x <- (x - range(x)[1])
  x <- (x / max(x)) * bins
  y <- (y - range(y)[1])
  y <- (y / max(y)) * bins
  emb2 <- data.frame(x = x, y = y)
  linetype <- 2
  emb2 <- data.frame(x = x, y = y)
  linetype <- 2;
  tmp <- emb2[rownames(emb2) %in% names(cell.type)[cell.type %in% cell], ]
  kd <- ks::kde(tmp, compute.cont = TRUE)
  lcn <- with(kd, contourLines(x = eval.points[[1]], y = eval.points[[2]], z = estimate, levels = cont[conf])[[1]])
  #name1 <- point.in.polygon(tmp[,1], tmp[,2], cn$x, cn$y)
  dd <- data.frame(lcn)
  dd$Z <- 1
  cn <- geom_path(aes(x, y), data = dd, linetype = linetype , color = color);
  return(cn)
}


##' @description Plot cell density 
##' @param bins number of bins for density esitmation, should keep consistent with bins in estimateCellDensity
##' @param col color palettes, 4 different color palettes are supported; default is yellow-black-magenta; BWR: blue-white-red;  WR: white-read; B: magma in viridi;
plotDensity <- function(mat, bins, col = 'BWR', legend = NULL, title = NULL, grid = NULL){
  p  <-  mat %>% as_tibble() %>% rowid_to_column(var = "X") %>% 
    gather(key = "Y", value = "Z", -1) %>% mutate(Y = as.numeric(gsub("V", "", Y))) %>% ggplot(aes(X, Y, fill = Z)) + 
    geom_raster() +
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), panel.border = element_blank(), 
                       panel.background = element_blank(), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_blank(), axis.text.y = element_blank()) +  
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
  
    if (col=='BWR'){
      p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(min(mat), max(mat)))
    }else if(col=='WR'){
      p <- p + scale_fill_gradient2(low = "white", high = "red", limits = c(min(mat), max(mat)))
    }else if(col=='B'){
      p <- p + scale_fill_viridis(option = 'B', alpha = 1, direction = 1, limits = c(min(mat), max(mat)))
    }else{
      p <- p + scale_fill_gradient2(low = "yellow", high = "magenta", mid = "black", midpoint = 0, limits = c(min(mat), max(mat)))
    }
  
    if (is.null(legend)){
      p <- p + theme(legend.position = "none")
    }
  
    if(!is.null(title)){
      p <- p + ggtitle(title)
    }
  
    if(!is.null(grid)){
      p <- p + geom_vline(xintercept=seq(30, bins, length.out=6), col='grey', alpha=0.1) 
      p <- p + geom_hline(yintercept=seq(30, bins, length.out=6), col='grey', alpha=0.1) 
    }
  
  return(p)
}






##' @description esitmate differential cell density
##' @param denMatrix.nor esitmated cell density matrix with estimateCellDensity
##' @param bins number of bins for density esitmation, should keep consistent with bins in estimateCellDensity
##' @param col color palettes, 4 different color palettes are supported; default is yellow-black-magenta; BWR: blue-white-red;  WR: white-read; B: magma in viridi;
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param fraction A two-level factor on the cell names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @method method to cacuated differential cell density of each bin; substract: target density minus ref density; entropy: estimated kl divergence entropy betwwen sample grapups ; t.test: zscore of t-test,global variacen is setting for t.test; 
diffCellDensity <- function(denMatrix.nor, fraction, sample.groups, bins, ref.level, target.level, method = 'substract', legend = NULL, grid = TRUE, col = 'YBM', title = NULL){
  NT <- names(sample.groups[sample.groups == target.level])
  NR <- names(sample.groups[sample.groups == ref.level])
  
  if (method == 'substract'){
    score = rowMeans(denMatrix.nor[, NT]) - rowMeans(denMatrix.nor[, NR])
  }else if (method == 'entropy'){
    sudo <- mean(as.numeric(denMatrix.nor)) # add sudo counts
    denMatrix.nor2 <- denMatrix.nor + sudo
    s1 <- rowSums(denMatrix.nor2[, NR])
    s2 <- rowSums(denMatrix.nor2[, NT])
    #s1=rowMeans(denMatrix.nor2[,NR])
    #s2=rowMeans(denMatrix.nor2[,NT])
    r1 <- s1 / (s1 + s2)
    r2 <- s2 / (s1 + s2)
    weight.sum.per.fac.cell <- data.frame(r1, r2)
    xt <- table(fraction)
    max.ent <- (if (xt[1] > xt[2]) c(0, 1) else c(1, 0)) %>% entropy::KL.empirical(xt, unit='log2')
    entropy.per.cell <- apply(weight.sum.per.fac.cell, 1, entropy::KL.empirical, xt, unit = 'log2') / max.ent
    score <- entropy.per.cell * sign(r2 - r1)
  }else if (method=='t.test'){
    vel <- rowMeans(denMatrix.nor)
    denMatrix.nor2 <- denMatrix.nor + sudo + quantile(vel, 0.05) # add sudo counts at 5%
    N1 <- denMatrix.nor2[, NR]
    T1 <- denMatrix.nor2[, NT]
    n1 <- length(as.numeric(N1))
    n2 <- length(as.numeric(T1))
    var.pooled <- weighted.mean(x = c(var(x1), var(x2)), w = c(n1 - 1, n2 - 1)) # caculate global variance 
    score <- apply(denMatrix.nor2, 1, function(x) {
      x1 <- x[NT]
      x2 <- x[NR]
      n1 <- length(x1)
      n2 <- length(x2)
      (mean(x1) - mean(x2)) / sqrt(var.pooled / n1 + var.pooled / n2)
    })
  }

  if (is.null(title)){
    title <- method
  }
  DensitScore <- matrix(score, ncol = bins, byrow = FALSE)
  
  p <- plotDensity(DensitScore, bins, col = col, title = title, legend = legend, grid = grid)
  return(p)
}




EntropySamples <- function(denMatrix.nor, samples, bins) {
  vel <- rowMeans(denMatrix.nor)
  sudo <- mean(as.numeric(denMatrix.nor)) # add sudo counts at 5%
  denMatrix.nor <- denMatrix.nor + sudo
  Z <- apply(denMatrix.nor[, samples], 1 , function(x) entropy(x) / log(length(x))) # empirical estimate near theoretical maximum
  Z <- 1 - Z
  density.score <- matrix(Z, ncol = bins, byrow = FALSE)
  return(density.score)
}


