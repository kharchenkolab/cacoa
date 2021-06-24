##' @title Expression shift magnitudes per cluster between conditions
##' @param count.matrices List of count matrices
##' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
##' @param cell.groups Named clustering/annotation factor with cell names
##' @param dist what distance measure to use: 'JS' - Jensen-Shannon divergence (default), 'cor' - Pearson's linear correlation on log transformed values
##' @param normalize.both (default=T)
##' @param n.cores number of cores (default=1)
##' @param verbose (default=F)
##' @param transposed.matrices (default=F)
##' @export
estimateExpressionShiftMagnitudes <- function(count.matrices, sample.groups, cell.groups, dist=c('cor', 'js'), normalize.both=TRUE,
                                              verbose=FALSE, transposed.matrices=FALSE, ref.level=NULL) {
  dist <- match.arg(dist)
  sample.groups <- as.factor(sample.groups) %>% na.omit() %>% droplevels()
  if(length(levels(sample.groups))!=2) stop("'sample.groups' must be a 2-level factor describing which samples are being contrasted")

  if (!transposed.matrices) {
    count.matrices %<>% lapply(Matrix::t)
  }

  common.genes <- Reduce(intersect, lapply(count.matrices, colnames))
  count.matrices %<>% lapply(`[`, , common.genes)

  # get a cell sample factor, restricted to the samples being contrasted
  sample.per.cell <- lapply(count.matrices[names(sample.groups)], rownames)
  sample.per.cell <- sample.per.cell %>% {rep(names(.), sapply(., length))} %>%
    setNames(unlist(sample.per.cell)) %>%  as.factor()
  cell.groups %<>% .[names(.) %in% names(sample.per.cell)]

  p.dist.info <- count.matrices %>%
    estimatePairwiseExpressionDistances(sample.per.cell=sample.per.cell, cell.groups=cell.groups,
                                        sample.groups=sample.groups, dist=dist)

  if(verbose) cat('calculating distances ... ')
  norm.type <- ifelse(normalize.both, "both", "ref")
  dist.df <- estimateExpressionShiftDf(p.dist.info, sample.groups, norm.type=norm.type, ref.level=ref.level) %>%
    mutate(Type=factor(Type, levels=names(sort(tapply(value, Type, median))))) # sort cell types
  if(verbose) cat('done!\n')

  return(list(dist.df=dist.df, p.dist.info=p.dist.info, sample.groups=sample.groups, cell.groups=cell.groups))
}

estimatePairwiseExpressionDistances <- function(count.matrices, sample.per.cell, cell.groups, sample.groups, dist) {
  # table of sample types and cells
  cell.groups %<>% as.factor()
  cct <- cell.groups %>% table(sample.per.cell[names(.)])
  caggr <- lapply(count.matrices, collapseCellsByType, groups=cell.groups, min.cell.count=1) %>%
    .[names(sample.groups)]

  # note: this is not efficient, as it will compare all samples on the two sides of the sample.groups
  #       would be faster to go only through the valid comparisons
  ctdm <- levels(cell.groups) %>% sccore:::sn() %>% lapply(function(ct) {
    tcm <- na.omit(do.call(rbind,lapply(caggr,function(x) x[match(ct,rownames(x)),])))

    if (dist=='js') {
      tcm <- t(tcm/pmax(1,rowSums(tcm)))
      dist.mat <- jsDist(tcm) %>% set_rownames(colnames(tcm)) %>% set_colnames(colnames(tcm))
    } else if (dist=='cor') {
      tcm <- log10(t(tcm / pmax(1, rowSums(tcm))) * 1e3 + 1)
      dist.mat <- 1 - cor(tcm)
      dist.mat[is.na(dist.mat)] <- 1;
    } else stop("Unknown distance: ", dist)
    # calculate how many cells there are
    attr(dist.mat, 'n.cells') <- cct[ct, colnames(tcm)]
    dist.mat
  })

  return(ctdm)
}

subsetDistanceMatrix <- function(dist.mat, sample.groups, cross.factor) {
  comp.selector <- if (cross.factor) "!=" else "=="
  selection.mask <- outer(sample.groups[rownames(dist.mat)], sample.groups[colnames(dist.mat)], comp.selector);
  diag(dist.mat) <- NA;
  dist.mat[!selection.mask] <- NA;

  if(all(is.na(dist.mat))) return(NULL);
  dist.df <- reshape2::melt(dist.mat) %>% na.omit()
  return(dist.df);
}

estimateCellTypeExpressionShiftDf <- function(dist.mat, sample.groups, norm.type=c("both", "ref", "none"), ref.level=NULL) {
  norm.type <- match.arg(norm.type)

  if ((norm.type == "ref") && is.null(ref.level))
    stop("ref.level has to be provided for norm.type='ref'")

  sample.groups %<>% .[rownames(dist.mat)]
  if (norm.type == "both") {
    sg1 <- levels(sample.groups)[1]
    m1 <- outer(sample.groups, sample.groups, function(a, b) (a == sg1) & (b == sg1))
    m2 <- outer(sample.groups, sample.groups, function(a, b) (a != sg1) & (b != sg1))
    diag(m1) <- NA; diag(m2) <- NA
    norm.const <- (median(dist.mat[m1], na.rm=TRUE) + median(dist.mat[m2], na.rm=TRUE)) / 2
  } else if (norm.type == "ref") {
    mr <- outer(sample.groups, sample.groups, function(a, b) (a == ref.level) & (b == ref.level))
    diag(mr) <- NA
    norm.const <- median(dist.mat[mr], na.rm=TRUE)
  } else {
    norm.const <- 0
  }

  dist.mat <- dist.mat - norm.const
  dist.df <- subsetDistanceMatrix(dist.mat, sample.groups=sample.groups, cross.factor=TRUE)

  return(dist.df)
}

estimateExpressionShiftDf <- function(p.dist.per.type, sample.groups, return.dists.within=FALSE, norm.type="both", ...) {
  if (return.dists.within) {
    dists.per.type <- lapply(p.dist.per.type, subsetDistanceMatrix, sample.groups, cross.factor=FALSE)
  } else {
    dists.per.type <- p.dist.per.type %>%
      lapply(estimateCellTypeExpressionShiftDf, sample.groups, norm.type=norm.type, ...)
  }

  dists.per.type %<>% .[!sapply(., is.null)]
  df <- names(dists.per.type) %>% lapply(function(n) cbind(dists.per.type[[n]], Type=n)) %>% do.call(rbind, .) %>%
    mutate(Condition=sample.groups[as.character(Var1)]) %>%
    na.omit()

  return(df)
}

prepareJointExpressionDistance <- function(p.dist.per.type, sample.groups=NULL, return.dists=TRUE) {
  # bring to a common set of cell types
  common.types <- lapply(p.dist.per.type, colnames) %>% unlist() %>% unique()

  p.dist.per.type %<>% lapply(function(x) {
    y <- matrix(0,nrow=length(common.types),ncol=length(common.types));  # can set the missing entries to zero, as they will carry zero weights
    rownames(y) <- colnames(y) <- common.types;
    y[rownames(x),colnames(x)] <- x;
    ycct <- setNames(rep(0,length(common.types)), common.types);
    ycct[colnames(x)] <- attr(x, 'n.cells')
    attr(y, 'n.cells') <- ycct
    y
  }) # reform the matrix to make sure all cell type have the same dimensions

  x <- abind::abind(lapply(p.dist.per.type, function(x) {
    nc <- attr(x, 'n.cells')
    #wm <- (outer(nc,nc,FUN='pmin'))
    wm <- sqrt(outer(nc, nc, FUN = 'pmin'))
    return(x * wm)
  }), along = 3)

  # just the weights (for total sum of weights normalization)
  y <- abind::abind(lapply(p.dist.per.type, function(x) {
    nc <- attr(x, 'n.cells')
    sqrt(outer(nc, nc, FUN = 'pmin'))
  }), along = 3)

  # normalize by total weight sums
  xd <- apply(x, c(1, 2), sum) / apply(y, c(1, 2), sum)

  if (return.dists)
    return(xd)

  cross.factor <- outer(sample.groups[rownames(xd)], sample.groups[colnames(xd)], '==')
  diag(xd) <- NA # remove self pairs
  # restrict
  xd[!cross.factor] <- NA
  if (!any(!is.na(xd)))
    return(NULL)
  xmd2 <- na.omit(reshape2::melt(xd))
  xmd2 <- na.omit(xmd2)
  xmd2$type1 <- sample.groups[as.character(xmd2$Var1)]
  xmd2$type2 <- sample.groups[as.character(xmd2$Var2)]

  return(xmd2)
}

filterExpressionDistanceInput <- function(cms, cell.groups, sample.per.cell, sample.groups,
                                          min.cells.per.sample=10, min.samp.per.type=2, min.gene.frac=0.01) {
  cell.names <- lapply(cms, rownames) %>% unlist()
  freq.table <- table(cell.groups[cell.names], sample.per.cell[cell.names]) %>%
    as.data.frame() %>%
    mutate(Condition=sample.groups[as.character(Var2)]) %>%
    filter(Freq >= min.cells.per.sample)

  filt.types <- freq.table %>% split(.$Var1) %>% sapply(function(df) {
    df %$% split(Var2, Condition) %>% sapply(length) %>% {all(. >= min.samp.per.type)}
  }) %>% which() %>% names()

  freq.table %<>% filter(Var1 %in% filt.types)
  filt.types.per.samp <- freq.table %$% split(Var1, Var2)

  cms.filt <- names(filt.types.per.samp) %>% sn() %>% lapply(function(n) {
    cms[[n]] %>% .[cell.groups[rownames(.)] %in% filt.types.per.samp[[n]],]
  })

  filt.genes <- lapply(cms.filt, function(cm) names(which(colMeans(cm > 1) > min.gene.frac))) %>%
    unlist() %>% table() %>% {. / length(cms.filt) > 0.1} %>% which() %>% names()

  cms.filt %<>% lapply(sccore:::extendMatrix, filt.genes) %>% lapply(`[`,,filt.genes)
  cell.names <- lapply(cms.filt, rownames) %>% unlist()

  return(list(cms=cms.filt, cell.groups=droplevels(cell.groups[cell.names]), sample.groups=sample.groups[names(cms)]))
}
