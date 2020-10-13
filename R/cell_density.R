
##' @description Estimate cell density in giving embedding 
##' @param emb cell embedding matrix
##' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
##' @param sample.groups @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @param bins number of bins for density esitmation, default 400
##' @param by.sample  if TRUE, density will esitmated by sample and quantiles normlization will applied to indivisual sample. If FALSE, cell condition.per.cell need to be provided and density will simply esitmated by condition.per.cell. 
##' @add.ponits add.ponits  show cells in density plot   
estimateCellDensity <- function(emb, sample.per.cell, sample.groups, bins, ref.level, target.level, condition.per.cell = NULL, by.sample = TRUE){
  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("You have to install preprocessCore package to do quantile normlization ")
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
                                matrix(rowMeans(tmp), ncol = bins, byrow = FALSE)
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
                                matrix(denMatrix[, x], ncol = bins, byrow = FALSE)
                              })
    }
  return(list('density.mat' = density.mat, 'density.fraction' = density.fraction))
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
plotDensity <- function(mat, bins, col = 'BWR', legend = NULL, title = NULL, grid = NULL, mi=NULL, ma=NULL){
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
    
    if (is.null(mi)){
      mi <- min(mat)*1.1
    }
    if (is.null(ma)){
      ma <- max(mat)*1.1
    }  
  
  
    if (col=='BWR'){ # 
      p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(mi, ma))
    }else if(col=='WR'){
      p <- p + scale_fill_gradient2(low = "white", high = "red", limits = c(mi, ma))
    }else if(col=='B'){
      p <- p + scale_fill_viridis(option = 'B', alpha = 1, direction = 1, limits = c(mi, ma))
    }else{ #purple-black-yellow
      p <- p + scale_fill_gradient2(low = "purple", high = "yellow", mid = "black", midpoint = 0, limits = c(mi, ma))
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
##' @param density.mat esitmated cell density matrix with estimateCellDensity
##' @param bins number of bins for density esitmation, should keep consistent with bins in estimateCellDensity
##' @param col color palettes, 4 different color palettes are supported; default is yellow-black-magenta; BWR: blue-white-red;  WR: white-read; B: magma in viridi;
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param condition.per.cell A two-level factor on the cell names describing the conditions being compared (default: stored vector)
##' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
##' @param target.level target/disease level for sample.group vector
##' @method method to cacuated differential cell density of each bin; substract: target density minus ref density; entropy: estimated kl divergence entropy betwwen sample grapups ; t.test: zscore of t-test,global variacen is setting for t.test; 
diffCellDensity <- function(density.mat, dcounts, condition.per.cell, sample.groups, bins, ref.level, target.level, method = 'substract', legend = NULL, grid = TRUE, col = 'YBM', title = NULL, dcount.cutoff = 2){
  nt <- names(sample.groups[sample.groups == target.level]) # sample name of target
  nr <- names(sample.groups[sample.groups == ref.level]) # sample name of reference 
  
  if (method == 'substract'){
    score = rowMeans(density.mat[, nt]) - rowMeans(density.mat[, nr])
  }else if (method == 'entropy'){
    sudo <- mean(as.numeric(density.mat)) # add sudo counts
    density.mat2 <- density.mat + sudo
    s1 <- rowSums(density.mat2[, nr])
    s2 <- rowSums(density.mat2[, nt])
    #s1=rowMeans(density.mat2[,NR])
    #s2=rowMeans(density.mat2[,NT])
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
      t.test(x1,x2)$statistic
    })
    
  } else if (method == 'willcox') {
    vel <- rowMeans(density.mat)
    density.mat2 <- density.mat + quantile(vel, 0.05) # add sudo counts at 5%
    score <- apply(density.mat2, 1, function(x) {
      mw = wilcox.test(x[nt], x[nr], exact = FALSE)
      zstat <- abs(qnorm(mw$p.value / 2))
      fc <- mean(x1) - mean(x2)
      zscore <- zstat * sign(fc)
      zscore
    })
  }

  if (is.null(title)){
    title <- method
  }
  

  density.score <- matrix(score, ncol = bins, byrow = FALSE)
  density.score[dcounts < dcount.cutoff] <- 0
  
  p <- plotDensity(density.score, bins, col = col, title = title, legend = legend, grid = grid)
  
  return(list('fig'=p,'score'=score))
}



