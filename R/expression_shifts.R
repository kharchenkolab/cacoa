##' Evaluate expression distances between clusters, separately within each sample, combining, optionally controling for cluster sizes.
##'
##' When using per-sample distance estimation (default, unless usr.aggregated.matrices=T), the distances for the cluster pairs that
##' are not seen sufficient number of times together (see min.samples) are set to NA.
##'
##' @title expression distance between clusters
##' @param con conos object
##' @param cell.groups Named factor with cell names and clustering/annotation
##' @param dist what distance measure to use: 'JS' - Jensen-Shannon divergence, 'cor' - Pearson's linear correlation on log transformed values
##' @param use.aggregated.matrices whether to simply aggregate all the molecules from all of the samples in which cluster appears (i.e. don't do aggregation of per-sample distance estimates, but just add all the data and measure the distance once)
##' @param use.single.cell.comparisons whether instead of adding up all molecules from a cluster in a given sample, instead compare n.cell draws of single cells (note: doesn't seem to work that well)
##' @param min.cluster.size minimum number of cells in a cluster (in a sample) for the distance to be estimated. default: 1
##' @param min.samples minimum number of samples in which a given pair of clusters has been compared (i.e. the two clusters were present in sufficient size). default:1
##' @param max.n.cells maximum number of cells to take from a cluster (to avoid very large clusters,default: Inf)
##' @param aggr aggregation function (default: median)
##' @param n.cores number of cores
##' @param return.details whether to return a list of individual sample distances
##' @param n.cells number of single cells to sample when use.single.cell.comparisons=T
##' @return distance matrix (with possible NAs), or if return.details=T a list with a distance matrix ($mdist) and a list of distance matrices computed on the individual samples ($dc)
##' @export
cluster.expression.distances <- function(con,cell.groups=NULL,dist='JS',n.cores=con$n.cores,min.cluster.size=1,min.samples=1,max.n.cells=Inf,aggr=median,return.details=FALSE,use.aggregated.matrices=FALSE,use.single.cell.comparisons=FALSE,n.cells=100) {
  # TODO: switch to abstracted accessor methods for con access
  if(is.null(cell.groups)) {
    if(is.null(con$clusters)) stop('no "cell.groups" specified and no clusterings found')
    cell.groups <- as.factor(con$clusters[[1]]$groups)
  } else {
    cell.groups <- as.factor(cell.groups)
  }

  valid.dists <- c('JS','cor');
  if(!dist %in% valid.dists) stop(paste('only the following distance types are supported:',paste(valid.dists,collapse=', ')))


  if(use.aggregated.matrices) {
    # in this case, we simply add all the counts from all the clusters into a common matrix and simply calculate distances on that
    tc <- conos:::rawMatricesWithCommonGenes(con) %>%
      lapply(conos:::collapseCellsByType, groups=cell.groups, min.cell.count=0) %>%
      abind::abind(along=3) %>%
      apply(c(1,2),sum,na.rm=T)

    if(dist=='JS') {
      tc <- t(tc/pmax(1,rowSums(tc)))
      tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
    } else { # correlation distance
      tc <- log10(t(tc/pmax(1,rowSums(tc)))*1e3+1)
      tcd <- 1-cor(tc)
    }
    diag(tcd) <- 0;
    if(return.details) {
      return(list(dc=tc,mdist=tcd))
    } else {
      return(tcd)
    }
  }


  # determine distance matrices for each sample
  dc <- abind::abind(conos:::papply(con$samples,function(s) {
    m <- s$misc$rawCounts
    cl <- factor(cell.groups[match(rownames(m),names(cell.groups))],levels=levels(cell.groups));
    tt <- table(cl);

    empty <- tt<min.cluster.size;

    if(use.single.cell.comparisons) { # sample individual cells and compare

      tcd <- lapply(1:n.cells,function(i) {
        scn <- unlist(tapply(names(cl),cl,function(x) sample(x,1)))
        tc <- as.matrix(m[na.omit(as.character(scn)),,drop=F])
        rownames(tc) <- names(scn)[!is.na(scn)]
        tc <- tc[match(levels(cl),rownames(tc)),]
        rownames(tc) <- levels(cl)
        if(dist=='JS') {
          tc <- t(tc/pmax(1,rowSums(tc)))
          tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
        } else { # correlation distance
          tc <- log10(t(tc/pmax(1,rowSums(tc)))*1e3+1)
          tcd <- 1-cor(tc)
        }
        tcd[empty,] <- tcd[,empty] <- NA;
        tcd
      }) %>% abind::abind(along=3) %>% apply(c(1,2),median,na.rm=T)

    } else { # aggregated clusters
      if(any(tt>max.n.cells)) { # need subsampling
        scn <- unlist(tapply(names(cl),cl,function(x) sample(x,min(max.n.cells,length(x)))))
        cl[!(names(cl) %in% scn)] <- NA; tt <- table(cl);
      }
      tc <- conos:::colSumByFactor(m,cl);
      tc <- tc[-1,,drop=F]  # omit NA cells
      if(dist=='JS') {
        tc <- t(tc/pmax(1,rowSums(tc)))
        tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
      } else { # correlation distance
        tc <- log10(t(tc/pmax(1,rowSums(tc)))*1e3+1)
        tcd <- 1-cor(tc)
      }
    }

    tcd[empty,] <- tcd[,empty] <- NA;
    diag(tcd) <- 0;
    # calculate how many cells there are
    attr(tcd,'cc') <- table(cl)
    tcd
  },n.cores=n.cores,mc.preschedule=T),along=3)

  # summarize across samples
  n.valid.obs <- apply(dc,c(1,2),function(x) sum(!is.na(x)))
  mdist <- apply(dc,c(1,2),aggr,na.rm=T)

  mdist[n.valid.obs<min.samples] <- NA;
  if(return.details) {
    return(list(dc=dc,mdist=mdist))
  } else {
    return(mdist)
  }
}

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
    valid.comparisons <- valid.comparisons | t(valid.comparisons)
    valid.comparisons <- valid.comparisons[rowSums(valid.comparisons)>0,colSums(valid.comparisons)>0]
    # ensure that only valid.comp groups are in the sample.groups
    sample.groups <- droplevels(sample.groups[names(sample.groups) %in% c(rownames(valid.comparisons),colnames(valid.comparisons))])
    if(length(levels(sample.groups))!=2) stop("insufficient number of levels in sample.groups after intersecting with valid.comparisons")

    # intersect with the cross-level pairs
    comp.matrix <- outer(sample.groups[rownames(valid.comparisons)],sample.groups[colnames(valid.comparisons)],'!='); diag(comp.matrix) <- FALSE
    valid.comparisons <- valid.comparisons & comp.matrix;
    # reduce and check again
    valid.comparisons <- valid.comparisons[rowSums(valid.comparisons)>0,colSums(valid.comparisons)>0]
    sample.groups <- droplevels(sample.groups[names(sample.groups) %in% c(rownames(valid.comparisons),colnames(valid.comparisons))])
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

  if(verbose) cat('running',n.subsamples,'subsamples ... \n')
  ctdml <- sccore:::plapply(1:n.subsamples,function(i) {
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
    caggr <- lapply(count.matrices, conos:::collapseCellsByType, groups=as.factor(cf), min.cell.count=1) %>%
      .[names(sample.groups)]

    # note: this is not efficient, as it will compare all samples on the two sides of the sample.groups
    #       would be faster to go only through the valid comparisons
    ctdm <- lapply(sn(levels(cf)),function(ct) {
      tcm <- na.omit(do.call(rbind,lapply(caggr,function(x) x[match(ct,rownames(x)),])))

      # restrict to top expressed genes
      if(n.top.genes<ncol(tcm)) tcm <- tcm[,rank(-colSums(tcm))>=n.top.genes]

      if(dist=='JS') {
        tcm <- t(tcm/pmax(1,rowSums(tcm)))
        tcd <- pagoda2:::jsDist(tcm); dimnames(tcd) <- list(colnames(tcm),colnames(tcm));
      } else {
        tcm <- log10(t(tcm/pmax(1,rowSums(tcm)))*1e3+1)
        tcd <- 1-cor(tcm)
        tcd[is.na(tcd)] <- 1;
      }
      # calculate how many cells there are
      attr(tcd,'cc') <- cct[ct,colnames(tcm)]
      tcd
    })

  },n.cores=n.cores, mc.preschedule=TRUE, mc.allow.recursive=FALSE, progress=verbose)

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
    df <- do.call(rbind,lapply(sn(names(x)),function(n) { z <- x[[n]]; z$Type <- n; z }))
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
##' @title Expression shift Z scores per cluster between conditions
##' @param pca Principal components
##' @param sample.per.cell Named factor with cell names indicating sample origin per cell
##' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
##' @param cell.groups Named factor with cell names and clustering/annotation
##' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt
estimateExpressionShiftZScores <- function(pca, sample.per.cell, sample.groups, cell.groups, ref.level) {
  sample.groups %<>% as.character() %>% setNames(names(sample.groups))
  mean.pc.per.samp.per.type <- split(names(sample.per.cell), sample.per.cell) %>%
    lapply(function(nsa) split(nsa, cell.groups[nsa]) %>% .[sapply(., length) > 0] %>%
             sapply(function(nss) colMeans(pca[nss,,drop=F])))

  res.df <- combn(names(mean.pc.per.samp.per.type), 2) %>% apply(2, function(ns)
    data.frame(S1=ns[1], S2=ns[2], value=getColumnwiseCorrelations(mean.pc.per.samp.per.type[[ns[1]]], mean.pc.per.samp.per.type[[ns[2]]]), stringsAsFactors=F) %>%
      tibble::as_tibble(rownames="Type")) %>% Reduce(rbind, .) %>%
    dplyr::mutate(SameCondition=(sample.groups[S1] == sample.groups[S2]),
                  Condition=ifelse(!SameCondition, "Between", sample.groups[S1])) %>%
    split(.$Type) %>% lapply(function(df)
      dplyr::mutate(df, distance=(mean(value[Condition == ref.level], trim=0.4) - value) / mad(value[Condition == ref.level]))) %>%
    Reduce(rbind, .) %>% dplyr::rename(correlation=value)

  return(res.df)
}
