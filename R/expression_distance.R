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
  if(!is.factor(sample.groups)) sample.groups <- as.factor(sample.groups)
  sample.groups <- droplevels(na.omit(sample.groups))
  if(length(levels(sample.groups))!=2) stop("'sample.groups' must be a 2-level factor describing which samples are being contrasted")

  if (!transposed.matrices) {
    count.matrices %<>% lapply(Matrix::t)
  }

  common.genes <- Reduce(intersect, lapply(count.matrices, colnames))
  count.matrices %<>% lapply(`[`, , common.genes)

  comp.matrix <- outer(sample.groups,sample.groups,'!='); diag(comp.matrix) <- FALSE

  # set up comparison mask
  if(is.null(valid.comparisons)) {
    # all cross-level pairs will be compared
    valid.comparisons <- comp.matrix;
  } else {
    # clean up valid.comparisons
    if(!all(rownames(valid.comparisons)==colnames(valid.comparisons))) stop('valid.comparisons must have the same row and column names')
    valid.comparisons %<>% {. | t(.)} %>% .[rowSums(.) > 0, colSums(.) > 0]
    # ensure that only valid.comp groups are in the sample.groups
    sample.groups %<>% .[names(.) %in% c(rownames(valid.comparisons), colnames(valid.comparisons))] %>% droplevels()
    if(length(levels(sample.groups))!=2) stop("insufficient number of levels in sample.groups after intersecting with valid.comparisons")

    # intersect with the cross-level pairs
    comp.matrix <- sample.groups %>% outer(.[rownames(valid.comparisons)], .[colnames(valid.comparisons)],'!=')
    diag(comp.matrix) <- FALSE
    valid.comparisons <- valid.comparisons & comp.matrix;
    # reduce and check again
    valid.comparisons <- valid.comparisons[rowSums(valid.comparisons)>0,colSums(valid.comparisons)>0]
    sample.groups %<>% .[names(.) %in% c(rownames(valid.comparisons), colnames(valid.comparisons))] %>% droplevels()
    if(length(levels(sample.groups))!=2) stop("insufficient number of levels in sample.groups after intersecting with valid.comparisons and sample.groups pairs")
    if(verbose) cat('a total of',(nrow(which(valid.comparisons,arr.ind=T))/2),'comparisons left after intersecting with valid.comparisons and sample.group pairs\n')
  }

  if(within.group.normalization) {
    control.matrix <- outer(sample.groups,sample.groups,'==');
    valid.comparisons <- valid.comparisons | control.matrix[rownames(valid.comparisons),colnames(valid.comparisons)]
  }

  # get a cell sample factor, restricted to the samples being contrasted
  cl <- lapply(count.matrices[names(sample.groups)], rownames)
  cl <- rep(names(cl), sapply(cl, length)) %>% setNames(unlist(cl)) %>%  as.factor()

  # cell factor
  cf <- cell.groups
  cf <- cf[names(cf) %in% names(cl)]

  if(is.null(n.cells)) {
    n.cells <- min(table(cf)) # use the size of the smallest group
    if(verbose) cat('setting group size of ',n.cells,' cells for comparisons\n')
  }

  if(verbose) cat('running',n.subsamples,'subsamples using ',n.cores,'cores ...\n')
  ctdml <- plapply(1:n.subsamples,function(i) {
    # subsample cells

    ## # draw cells without sample stratification - this can drop certain samples, particularly those with lower total cell numbers
    ## cf <- tapply(names(cf),cf,function(x) {
    ##   if(length(x)<=n.cells) { return(cf[x]) } else { setNames(rep(cf[x[1]],n.cells), sample(x,n.cells)) }
    ## })

    # calculate expected mean number of cells per sample and aim to sample that
    n.cells.scaled <- max(min.cells,ceiling(n.cells/length(sample.groups)));
    cf <- tapply(names(cf),list(cf,cl[names(cf)]),function(x) {
      if(length(x)<=n.cells.scaled) { return(cf[x]) } else { setNames(rep(cf[x[1]],n.cells.scaled), sample(x,n.cells.scaled)) }
    })

    cf <- as.factor(setNames(unlist(lapply(cf,as.character)),unlist(lapply(cf,names))))

    # table of sample types and cells
    cct <- table(cf,cl[names(cf)])
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
      attr(tcd,'cc') <- cct[ct,colnames(tcm)]
      tcd
    })

  },n.cores=n.cores, mc.preschedule=TRUE, progress=verbose)

  if(verbose) cat('calculating distances ... ')
  df <- do.call(rbind,lapply(ctdml,function(ctdm) {

    x <- lapply(ctdm,function(xm) {
      nc <- attr(xm,'cc');
      wm <- outer(nc,nc,FUN='pmin')

      cross.factor <- outer(sample.groups[rownames(xm)],sample.groups[colnames(xm)],'!=');
      frm <- valid.comparisons[rownames(xm),colnames(xm)] & cross.factor

      if(within.group.normalization) {
        frm.cont <- valid.comparisons[rownames(xm),colnames(xm)] & !cross.factor
        med.cont <- median(na.omit(xm[frm.cont]))
        xm <- xm/med.cont
      }

      diag(xm) <- NA;

      # restrict
      xm[!frm] <- NA;
      xm[wm<min.cells] <- NA;
      if(!any(!is.na(xm))) return(NULL);
      xmd <- na.omit(reshape2::melt(xm))
      wm[is.na(xm)] <- NA;
      xmd$n <- na.omit(reshape2::melt(wm))$value
      return(xmd);
    })

    x <- x[!unlist(lapply(x,is.null))]
    df <- do.call(rbind,lapply(sccore:::sn(names(x)),function(n) { z <- x[[n]]; z$Type <- n; z }))
    df$patient <- df$Var1
    df
  }))

  # median across pairs
  df <- do.call(rbind,tapply(1:nrow(df),paste(df$Var1,df$Var2,df$Type,sep='!!'),function(ii) {
    ndf <- data.frame(df[ii[1],,drop=F]);
    ndf$value <- median(df$value[ii])
    ndf$n <- median(df$n[ii])
    ndf
  }))

  # sort cell types
  df$Type <- factor(df$Type,levels=names(sort(tapply(df$value,as.factor(df$Type),median))))

  if(verbose) cat('done!\n')
  return(list(df=df, ctdml=ctdml, sample.groups=sample.groups, valid.comparisons=valid.comparisons))
}


plotExpressionDistanceIndividual <- function(ctdml, valid.comparisons, sample.groups=NULL, min.cells=10, ...) {
  df <- do.call(rbind,lapply(ctdml,function(ctdm) {
    x <- lapply(ctdm, function(xm) {
      nc <- attr(xm, 'cc')
      wm <- outer(nc, nc, FUN = 'pmin')
      cross.factor <- outer(sample.groups[rownames(xm)], sample.groups[colnames(xm)], '==')
      frm <- valid.comparisons[rownames(xm), colnames(xm)] & cross.factor
      diag(xm) <- NA
      # remove self pairs
      xm[!frm] <- NA
      xm[wm < min.cells] <- NA
      if (!any(!is.na(xm)))
        return(NULL)
      xmd <- na.omit(reshape2::melt(xm))
      wm[is.na(xm)] <- NA
      xmd$n <- na.omit(reshape2::melt(wm))$value
      return(xmd)
    })
    x <- x[!unlist(lapply(x, is.null))]
    df <- do.call(rbind, lapply(sccore:::sn(names(x)), function(n) {
      z <- x[[n]]
      z$type <- n
      z
    }))
    df$patient <- df$Var1
    df$type1 <- sample.groups[df$Var1]
    df$type2 <- sample.groups[df$Var2]
    df
  }))

  # median across pairs
  df <- tapply(1:nrow(df), paste(df$Var1, df$Var2, df$type, sep = '!!'), function(ii) {
    data.frame(df[ii[1],, drop = F]) %>%
      mutate(value=median(df$value[ii]), n=median(df$n[ii]))
  }) %>%
    do.call(rbind, .) %>% rename(group=type1, variable=type) %>% na.omit()

  gg <- plotCountBoxplotsPerType(df, y.lab="expression distance", y.expand=c(0, max(df$value) * 0.1), ...)

  return(gg)
}

prepareJointExpressionDistance <- function(ctdml, valid.comparisons=NULL, sample.groups=NULL) {
  df <- lapply(ctdml, function(ctdm) {
    # bring to a common set of cell types
    commoncell <- unique( unlist( lapply(ctdm, function(x) colnames(x)) ))

    ctdm <-  lapply(ctdm, function(x) {
      y <- matrix(0,nrow=length(commoncell),ncol=length(commoncell)); rownames(y) <- colnames(y) <- commoncell; # can set the missing entries to zero, as they will carry zero weights
      y[rownames(x),colnames(x)] <- x;
      ycct <- setNames(rep(0,length(commoncell)), commoncell);
      ycct[colnames(x)] <- attr(x,'cc')
      attr(y, 'cc') <- ycct
      y
    }) # reform the matrix to make sure all cell type have the same dimensions

    x <- abind::abind(lapply(ctdm, function(x) {
      nc <- attr(x, 'cc')
      #wm <- (outer(nc,nc,FUN='pmin'))
      wm <- sqrt(outer(nc, nc, FUN = 'pmin'))
      return(x * wm)
    }), along = 3)

    # just the weights (for total sum of weights normalization)
    y <- abind::abind(lapply(ctdm, function(x) {
      nc <- attr(x, 'cc')
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

plotExpressionDistanceJoint <- function(ctdml, valid.comparisons, sample.groups=NULL, ...) {
  df <- prepareJointExpressionDistance(ctdml, valid.comparisons=valid.comparisons, sample.groups=sample.groups) %>%
    do.call(rbind, .)
  df2 <- do.call(rbind, tapply(1:nrow(df), paste(df$Var1, df$Var2, sep = '!!'), function(ii) {
    ndf <- data.frame(df[ii[1],, drop = FALSE])
    ndf$value <- median(df$value[ii])
    ndf$n <- median(df$n[ii])
    ndf
  }))

  df2$group <- df2$type1
  df2$variable <- ""

  gg <- plotCountBoxplotsPerType(df2, y.lab="expression distance", y.expand=c(0, max(df$value) * 0.1), ...) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())

  return(gg)
}
