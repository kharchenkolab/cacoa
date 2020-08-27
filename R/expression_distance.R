


##' @description   extract aggregated counts matrix 
##' @param count.matrices list object of count matrix
getClusterCountMatrices<-function (count.matrices, groups = NULL, common.genes = TRUE, omit.na.cells = TRUE) {
  groups <- as.factor(groups)
  matl <- lapply(count.matrices, function(m) {
    cl <- factor(groups[match(rownames(m), names(groups))], 
                 levels = levels(groups))
    tc <- conos:::colSumByFactor(m, cl)
    if (omit.na.cells) {
      tc <- tc[-1, , drop = F]
    }
    t(tc)
  })
  if (common.genes) {
    gs <- unique(unlist(lapply(matl, rownames)))
    matl <- lapply(matl, function(m) {
      nm <- matrix(0, nrow = length(gs), ncol = ncol(m))
      colnames(nm) <- colnames(m)
      rownames(nm) <- gs
      mi <- match(rownames(m), gs)
      nm[mi, ] <- m
      nm
    })
  }
  return(matl)
}



##'  @description esitmate expression distance between samples of each cell type  
##'  @param dist what distance measure to use: 'JS' - Jensen-Shannon divergence, 'cor' - Pearson's linear correlation on log transformed values
##'  @param min.cluster.size minimum number of cells in a cluster (in a sample) for the distance to be estimated. default: 10
##'  @param samplef  Named sample factor with cell names (default: stored vector)
##'  @param cell.type Vector indicating cell groups with cell names (default: stored vector)
esitmateExpressionDsiatnce <- function(count.matrices,samplef,cell.type,dist='JS',min.cells =10){
  
  cm <- getClusterCountMatrices(count.matrices,groups=cell.type)
  cct <- table(cell.type,samplef[names(cell.type)])
  ctdm <- lapply(sccore:::sn(colnames(cm[[1]])),function(ct) {
    tcm <- do.call(rbind,lapply(cm,function(x) x[,ct]))
    tcm <- t(tcm/pmax(1,rowSums(tcm)))
    if(dist=='JS') {
      tcd <- pagoda2:::jsDist(tcm); dimnames(tcd) <- list(colnames(tcm),colnames(tcm));
    }else{ # correlation distance
      tc <- log10(t(tc/pmax(1,rowSums(tc)))*1e3+1)
      tcd <- 1-cor(tc)
    }
    # calculate how many cells there are
    attr(tcd,'cc') <- cct[ct,colnames(tcm)]
    tcd
  })
  
  cells <- apply(cct,1,function(x) length(x[x>min.cells])) %>% .[.>2] %>% names()
  
  xlist <- lapply(sccore:::sn(cells),function(ct) {
    #print(ct)
    xd <- ctdm[[ct]]
    nc <- attr(ctdm[[ct]],'cc');
    vi <- nc[rownames(xd)]>=min.cells;
    xd <- xd[vi,vi]
    xd
  })
  
  
   # a cube across all cell types and sample pairs
   # weights of individual cell types determined by the minimal number of cells on each side of the pairwise comparison
   x <- abind(lapply(ctdm,function(x) {
     nc <- attr(x,'cc');
     #wm <- (outer(nc,nc,FUN='pmin'))
     wm <- sqrt(outer(nc,nc,FUN='pmin'))
     return( x*wm )
   }),along=3)
   # just the weights (for total sum of weights normalization)
   y <- abind(lapply(ctdm,function(x) {
     nc <- attr(x,'cc');
     sqrt(outer(nc,nc,FUN='pmin'))
   }),along=3)
   
   # normalize by total weight sums
   xd <- apply(x,c(1,2),sum)/apply(y,c(1,2),sum)
  
   return(list('xd'=xd,'xdlist'=xlist,'cct'=cct))
}


##' @description Plot expression distance 
##' @param disData esitmated cell expression distance with esitmate.expression.dsiatnce
##' @param type figure output format to measure inter distance or intra distance;
##' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
##' @param cell.Type default is null for weigeted expression distance across mutiple cell types, if setting, draw plot for specific cell type.
##' @return A ggplot2 object 
plotExpressionDistance <- function(dist,sample.groups,cell.type=NULL,type='inter'){
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("You have to install Rtsne package to use this function")
  }
  
  cct <- dist[['cct']]
  if(is.null(cell.type)){
    x <- dist[['xd']]
    cell.type <- 'combined'
    nc <- colSums(cct)
  }else{
    x <- dist[['xdlist']][[cell.type]]
    nc <- cct[cell.type, ]
  }
  
  if (type == 'inter'){
    xde <- Rtsne::Rtsne(x, is_distance = TRUE, perplexity = 4, max_iter = 1e4)$Y
    df <- data.frame(xde)
    rownames(df) <- rownames(x)
    colnames(df) <- c("x", "y")
    df$fraction <- sample.groups[rownames(df)]
    df$ncells <- nc[rownames(df)]
    p <- ggplot(df, aes(x, y, color=fraction, shape=fraction, size=log10(ncells))) + geom_point() +
      theme_bw() + ggtitle(cell.type) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  }else if (type == 'intra'){
      x[upper.tri(x)] <- NA; diag(x) <- NA;
      df2 <- na.omit(melt(x))
      df2$fraction1 <- sample.groups[df2$Var1]
      df2$fraction2 <- sample.groups[df2$Var2]
      df2$samePatient <- df2$Var1==df2$Var2;
      df2$sameFraction <- df2$fraction1==df2$fraction2;
      
      df3 <- df2[df2$sameFraction,]
      df3$type <- df3$fraction1
      p <- ggplot(na.omit(df3), aes(x = type, y = value))+
        geom_boxplot(notch = TRUE, outlier.shape = NA, aes(fill = type))+ ggtitle(cell.type) +
        geom_jitter(position = position_jitter(0.2), color = adjustcolor('black', alpha = 0.2))+ theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) +  
        guides(fill = FALSE) + xlab('') + ylab('expression distance')
    }else{
      stop("'type' must be provided either during the object initialization or during this function call")
    }
    return(p)
  }



