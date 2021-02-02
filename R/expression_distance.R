##' @title Expression shift magnitudes per cluster between conditions
##' @param count.matrices List of count matrices
##' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
##' @param cell.groups Named clustering/annotation factor with cell names
##' @param dist what distance measure to use: 'JS' - Jensen-Shannon divergence (default), 'cor' - Pearson's linear correlation on log transformed values
##' @param within.group.normalization (default=T)
##' @param valid.comparisons (default=NULL)
##' @param n.cells Number of cells per group (default=NULL)
##' @param n.top.genes (default=Inf)
##' @param n.subsamples (default=100)
##' @param min.cells (default=10)
##' @param n.cores number of cores (default=1)
##' @param verbose (default=F)
##' @param transposed.matrices (default=F)
##' @export
estimateExpressionShiftMagnitudes <- function(count.matrices, sample.groups, cell.groups, dist='JS', within.group.normalization=TRUE,
                                              valid.comparisons=NULL, n.cells=NULL, n.top.genes=Inf, n.subsamples=100, min.cells=10,
                                              n.cores=1, verbose=FALSE, transposed.matrices=FALSE) {
  sample.groups <- as.factor(sample.groups) %>% na.omit() %>% droplevels()
  if(length(levels(sample.groups))!=2) stop("'sample.groups' must be a 2-level factor describing which samples are being contrasted")

  if (!transposed.matrices) {
    count.matrices %<>% lapply(Matrix::t)
  }

  common.genes <- Reduce(intersect, lapply(count.matrices, colnames))
  count.matrices %<>% lapply(`[`, , common.genes)

  # set up comparison mask
  vc.res <- extractValidComparisons(valid.comparisons, sample.groups, within.group.normalization=within.group.normalization)
  valid.comparisons <- vc.res$valid.comparisons; sample.groups <- vc.res$sample.groups

  # get a cell sample factor, restricted to the samples being contrasted
  cl <- lapply(count.matrices[names(sample.groups)], rownames)
  cl <- rep(names(cl), sapply(cl, length)) %>% setNames(unlist(cl)) %>%  as.factor()
  cell.groups %<>% .[names(.) %in% names(cl)]

  if(is.null(n.cells)) {
    n.cells <- min(table(cf)) # use the size of the smallest group
    if(verbose) cat('setting group size of ',n.cells,' cells for comparisons\n')
  }

  if(verbose) cat('running',n.subsamples,'subsamples using ',n.cores,'cores ...\n')

  n.cells.scaled <- max(min.cells, ceiling(n.cells / length(sample.groups)))
  p.dist.info <- plapply(1:n.subsamples, function(i) {
    subsamplePairwiseExpressionDistances(count.matrices, cl=cl, cell.groups=cell.groups, n.top.genes=n.top.genes,
                                         sample.groups=sample.groups, n.cells.scaled=n.cells.scaled, dist=dist)
  },n.cores=n.cores, mc.preschedule=TRUE, progress=verbose)

  if(verbose) cat('calculating distances ... ')
  df <- aggregateExpressionShiftMagnitudes(p.dist.info, valid.comparisons, sample.groups, min.cells=min.cells,
                                           within.group.normalization=within.group.normalization, comp.filter='!=') %>%
    mutate(Type=factor(Type, levels=names(sort(tapply(value, Type, median))))) # sort cell types
  if(verbose) cat('done!\n')

  return(list(df=df, p.dist.info=p.dist.info, sample.groups=sample.groups, valid.comparisons=valid.comparisons))
}

extractValidComparisons <- function(valid.comparisons, sample.groups, within.group.normalization) {
  # set up comparison mask
  if(is.null(valid.comparisons)) {
    # all cross-level pairs will be compared
    valid.comparisons <- outer(sample.groups, sample.groups, '!=')
    diag(valid.comparisons) <- FALSE
  } else {
    # clean up valid.comparisons
    if(!all(rownames(valid.comparisons) == colnames(valid.comparisons))) stop('valid.comparisons must have the same row and column names')
    valid.comparisons %<>% {. | t(.)} %>% .[rowSums(.) > 0, colSums(.) > 0]
    # ensure that only valid.comp groups are in the sample.groups
    sample.groups %<>% .[names(.) %in% c(rownames(valid.comparisons), colnames(valid.comparisons))] %>% droplevels()
    if(length(levels(sample.groups))!=2) stop("insufficient number of levels in sample.groups after intersecting with valid.comparisons")

    # intersect with the cross-level pairs
    comp.matrix <- sample.groups %>% outer(.[rownames(valid.comparisons)], .[colnames(valid.comparisons)], '!=')
    diag(comp.matrix) <- FALSE
    valid.comparisons <- valid.comparisons & comp.matrix;
    # reduce and check again
    valid.comparisons %<>% .[rowSums(.) > 0, colSums(.) > 0]
    sample.groups %<>% .[names(.) %in% c(rownames(valid.comparisons), colnames(valid.comparisons))] %>% droplevels()
    if(length(levels(sample.groups))!=2) stop("insufficient number of levels in sample.groups after intersecting with valid.comparisons and sample.groups pairs")
    if(verbose) cat('a total of', (nrow(which(valid.comparisons, arr.ind=TRUE)) / 2), 'comparisons left after intersecting with valid.comparisons and sample.group pairs\n')
  }

  if(within.group.normalization) {
    control.matrix <- outer(sample.groups, sample.groups, '==');
    valid.comparisons <- valid.comparisons | control.matrix[rownames(valid.comparisons), colnames(valid.comparisons)]
  }

  return(list(valid.comparisons=valid.comparisons, sample.groups=sample.groups))
}

subsamplePairwiseExpressionDistances <- function(count.matrices, cl, cell.groups, sample.groups, n.cells.scaled,
                                                 n.top.genes, dist) {
  ## # draw cells without sample stratification - this can drop certain samples, particularly those with lower total cell numbers
  ## cf <- tapply(names(cf),cf,function(x) {
  ##   if(length(x)<=n.cells) { return(cf[x]) } else { setNames(rep(cf[x[1]],n.cells), sample(x,n.cells)) }
  ## })

  # calculate expected mean number of cells per sample and aim to sample that
  cf <- tapply(names(cell.groups), list(cell.groups, cl[names(cell.groups)]),function(x) {
    if(length(x)<=n.cells.scaled) { return(cell.groups[x]) } else { setNames(rep(cell.groups[x[1]], n.cells.scaled), sample(x, n.cells.scaled)) }
  })

  cf <- as.factor(setNames(unlist(lapply(cf,as.character)),unlist(lapply(cf,names))))

  # table of sample types and cells
  cct <- table(cf, cl[names(cf)])
  caggr <- lapply(count.matrices, collapseCellsByType, groups=as.factor(cf), min.cell.count=1) %>%
    .[names(sample.groups)]

  # note: this is not efficient, as it will compare all samples on the two sides of the sample.groups
  #       would be faster to go only through the valid comparisons
  ctdm <- lapply(sccore:::sn(levels(cf)),function(ct) {
    tcm <- na.omit(do.call(rbind,lapply(caggr,function(x) x[match(ct,rownames(x)),])))

    # restrict to top expressed genes
    if(n.top.genes<ncol(tcm)) tcm <- tcm[,rank(-colSums(tcm))>=n.top.genes]

    if(dist=='JS') {
      tcm <- t(tcm/pmax(1,rowSums(tcm)))
      tcd <- pagoda2:::jsDist(tcm, ncores = 1); dimnames(tcd) <- list(colnames(tcm),colnames(tcm));
    } else {
      tcm <- log10(t(tcm/pmax(1,rowSums(tcm)))*1e3+1)
      tcd <- 1-cor(tcm)
      tcd[is.na(tcd)] <- 1;
    }
    # calculate how many cells there are
    attr(tcd, 'n.cells') <- cct[ct, colnames(tcm)]
    tcd
  })

  return(ctdm)
}

aggregateExpressionShiftMagnitudes <- function(p.dist.info, valid.comparisons, sample.groups, min.cells,
                                               within.group.normalization=FALSE, comp.filter='!=') {
  df <- do.call(rbind, lapply(p.dist.info, function(ctdm) {
    x <- lapply(ctdm, function(xm) {
      n.cells <- attr(xm, 'n.cells');
      wm <- outer(n.cells, n.cells, FUN='pmin')

      cross.factor <- outer(sample.groups[rownames(xm)], sample.groups[colnames(xm)], comp.filter);
      frm <- valid.comparisons[rownames(xm),colnames(xm)] & cross.factor

      if (within.group.normalization) {
        frm.cont <- valid.comparisons[rownames(xm),colnames(xm)] & !cross.factor
        med.cont <- median(na.omit(xm[frm.cont]))
        xm <- xm / med.cont
      }

      diag(xm) <- NA;

      # remove self pairs
      xm[!frm] <- NA;
      xm[wm < min.cells] <- NA;
      if(!any(!is.na(xm))) return(NULL);
      xmd <- na.omit(reshape2::melt(xm))
      wm[is.na(xm)] <- NA;
      xmd$n <- na.omit(reshape2::melt(wm))$value
      return(xmd);
    })

    x <- x[!sapply(x, is.null)]
    df <- names(x) %>% lapply(function(n) cbind(x[[n]], Type=n)) %>% do.call(rbind, .)
    df
  }))

  df %<>% group_by(Var1, Var2, Type) %>%
    summarize(value=median(value), n=median(n)) %>%
    mutate(Condition=sample.groups[as.character(Var1)]) %>%
    na.omit()

  return(df)
}

prepareJointExpressionDistance <- function(p.dist.info, valid.comparisons=NULL, sample.groups=NULL) {
  df <- lapply(p.dist.info, function(ctdm) {
    # bring to a common set of cell types
    commoncell <- unique( unlist( lapply(ctdm, function(x) colnames(x)) ))

    ctdm <-  lapply(ctdm, function(x) {
      y <- matrix(0,nrow=length(commoncell),ncol=length(commoncell)); rownames(y) <- colnames(y) <- commoncell; # can set the missing entries to zero, as they will carry zero weights
      y[rownames(x),colnames(x)] <- x;
      ycct <- setNames(rep(0,length(commoncell)), commoncell);
      ycct[colnames(x)] <- attr(x, 'n.cells')
      attr(y, 'n.cells') <- ycct
      y
    }) # reform the matrix to make sure all cell type have the same dimensions

    x <- abind::abind(lapply(ctdm, function(x) {
      nc <- attr(x, 'n.cells')
      #wm <- (outer(nc,nc,FUN='pmin'))
      wm <- sqrt(outer(nc, nc, FUN = 'pmin'))
      return(x * wm)
    }), along = 3)

    # just the weights (for total sum of weights normalization)
    y <- abind::abind(lapply(ctdm, function(x) {
      nc <- attr(x, 'n.cells')
      sqrt(outer(nc, nc, FUN = 'pmin'))
    }), along = 3)

    # normalize by total weight sums
    xd <- apply(x, c(1, 2), sum) / apply(y, c(1, 2), sum)

    if (is.null(valid.comparisons))
      return(xd)

    cross.factor <- outer(sample.groups[rownames(xd)], sample.groups[colnames(xd)], '==')
    frm <- valid.comparisons[rownames(xd), colnames(xd)] & cross.factor

    diag(xd) <- NA # remove self pairs
    # restrict
    xd[!frm] <- NA
    if (!any(!is.na(xd)))
      return(NULL)
    xmd2 <- na.omit(reshape2::melt(xd))
    xmd2 <- na.omit(xmd2)
    xmd2$type1 <- sample.groups[xmd2$Var1]
    xmd2$type2 <- sample.groups[xmd2$Var2]
    xmd2
  })

  return(df)
}
