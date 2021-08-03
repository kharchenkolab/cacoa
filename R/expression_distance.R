##' @title Expression shift magnitudes per cluster between conditions
##' @param cm.per.type List of count matrices per cell type
##' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
##' @param cell.groups Named clustering/annotation factor with cell names
##' @param dist what distance measure to use: 'JS' - Jensen-Shannon divergence (default), 'cor' - Pearson's linear correlation on log transformed values
##' @param normalize.both (default=T)
##' @param n.cores number of cores (default=1)
##' @param verbose (default=F)
##' @param transposed.matrices (default=F)
##' @export
estimateExpressionShiftMagnitudes <- function(cm.per.type, sample.groups, cell.groups, sample.per.cell,
                                              dist=c('cor', 'js', 'l2'), normalize.both=TRUE, verbose=FALSE,
                                              ref.level=NULL, ...) {
  dist <- match.arg(dist)
  sample.groups <- as.factor(sample.groups) %>% na.omit() %>% droplevels()
  if(length(levels(sample.groups)) != 2)
    stop("'sample.groups' must be a 2-level factor describing which samples are being contrasted")

  cell.groups %<>% .[names(.) %in% names(sample.per.cell)]

  if(verbose) cat('Calculating distances ... ')
  p.dist.info <- cm.per.type %>%
    estimatePairwiseExpressionDistances(sample.per.cell=sample.per.cell, cell.groups=cell.groups,
                                        sample.groups=sample.groups, dist=dist, ...)

  norm.type <- ifelse(normalize.both, "both", "ref")
  dist.df <- estimateExpressionShiftDf(p.dist.info, sample.groups, norm.type=norm.type, ref.level=ref.level) %>%
    mutate(Type=factor(Type, levels=names(sort(tapply(value, Type, median))))) # sort cell types
  if(verbose) cat('done!\n')

  return(list(dist.df=dist.df, p.dist.info=p.dist.info, sample.groups=sample.groups, cell.groups=cell.groups))
}

estimatePairwiseExpressionDistances <- function(cm.per.type, sample.per.cell, cell.groups, sample.groups, dist,
                                                top.n.genes=NULL, n.pcs=NULL, ...) {
  cell.groups %<>% as.factor()
  # table of sample types and cells
  cct <- cell.groups %>% table(sample.per.cell[names(.)])

  ctdm <- levels(cell.groups) %>% sccore:::sn() %>% lapply(function(ct) {
    tcm <- cm.per.type[[ct]]
    cm.norm <- tcm / pmax(1, rowSums(tcm))

    if (!is.null(top.n.genes)) {
      sel.genes <- filterGenesForCellType(tcm, sample.groups=sample.groups, top.n.genes=top.n.genes, ...)
      cm.norm <- cm.norm[,sel.genes,drop=FALSE]
    }

    if (dist=='js') {
      dist.mat <- t(cm.norm) %>% jsDist() %>%
        set_rownames(rownames(cm.norm)) %>% set_colnames(rownames(cm.norm))
    } else if (dist %in% c('cor', 'l2')) {
      cm.norm <- log10(cm.norm * 1e3 + 1)
      if (!is.null(n.pcs)) {
        min.dim <- min(dim(cm.norm)) - 1
        if (n.pcs > min.dim) {
          warning("n.pcs is too large for cell type ", ct, ". Setting it to maximal allowed value ", min.dim)
          n.pcs <- min.dim
        }
        pcs <- irlba::irlba(cm.norm, nv=n.pcs, nu=0, right_only=FALSE, fastpath=TRUE,
                            maxit=3000, reorth=TRUE, verbose=FALSE)
        cm.norm <- as.matrix(cm.norm %*% pcs$v)
      }

      if (dist == 'cor') {
        dist.mat <- 1 - cor(t(cm.norm))
      } else {
        dist.mat <- dist(cm.norm) %>% as.matrix()
      }

      dist.mat[is.na(dist.mat)] <- 1;
    } else stop("Unknown distance: ", dist)
    # calculate how many cells there are
    attr(dist.mat, 'n.cells') <- cct[ct, rownames(cm.norm)]
    dist.mat
  })

  return(ctdm)
}

subsetDistanceMatrix <- function(dist.mat, sample.groups, cross.factor, return.raw=FALSE) {
  comp.selector <- if (cross.factor) "!=" else "=="
  selection.mask <- outer(sample.groups[rownames(dist.mat)], sample.groups[colnames(dist.mat)], comp.selector);
  diag(dist.mat) <- NA;
  if (return.raw)
    return(na.omit(dist.mat[selection.mask]))

  dist.mat[!selection.mask] <- NA;

  if(all(is.na(dist.mat))) return(NULL);
  dist.df <- reshape2::melt(dist.mat) %>% na.omit()
  return(dist.df);
}

estimateCellTypeExpressionShiftDf <- function(dist.mat, sample.groups, norm.type=c("both", "ref", "none"),
                                              ref.level=NULL, return.raw=FALSE) {
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
  dist.df <- subsetDistanceMatrix(dist.mat, sample.groups=sample.groups, cross.factor=TRUE, return.raw=return.raw)

  return(dist.df)
}

estimateExpressionShiftDf <- function(p.dist.per.type, sample.groups, return.dists.within=FALSE, norm.type="both",
                                      return.raw=FALSE, ...) {
  if (return.dists.within) {
    dists.per.type <- lapply(p.dist.per.type, subsetDistanceMatrix, sample.groups,
                             return.raw=return.raw, cross.factor=FALSE)
  } else {
    dists.per.type <- p.dist.per.type %>%
      lapply(estimateCellTypeExpressionShiftDf, sample.groups, norm.type=norm.type, return.raw=return.raw, ...)
  }

  if (return.raw)
    return(dists.per.type)

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
                                          min.cells.per.sample=10, min.samp.per.type=2, min.gene.frac=0.01,
                                          genes=NULL) {
  # Filter rare samples per cell type
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

  cell.names <- lapply(cms.filt, rownames) %>% unlist()

  # Filter low-expressed genes
  if (is.null(genes)) {
    genes <- lapply(cms.filt, function(cm) {
      cm@x <- 1 * (cm@x > 1)
      names(which(colMeans(cm) > min.gene.frac))
    }) %>% unlist() %>% table() %>% {. / length(cms.filt) > 0.1} %>% which() %>% names()
  }

  # Collapse matrices and extend to the same genes
  cms.filt %<>% lapply(collapseCellsByType, groups=cell.groups, min.cell.count=1)

  cms.filt %<>% lapply(sccore:::extendMatrix, genes) %>% lapply(`[`,,genes, drop=FALSE)

  # Group matrices by cell type
  cell.groups <- droplevels(cell.groups[cell.names])

  cm.per.type <- levels(cell.groups) %>% sccore:::sn() %>% lapply(function(ct) {
    lapply(cms.filt, function(x) x[match(ct, rownames(x)),]) %>%
      do.call(rbind, .) %>% na.omit()
  })

  return(list(cm.per.type=cm.per.type, cell.groups=cell.groups, sample.groups=sample.groups[names(cms)]))
}

estimateExpressionShiftPValues <- function(p.dist.info, sample.groups, n.permutations=1000, verbose=TRUE, n.cores=1, trim=0.1, ...) {
  dist.per.type <- estimateExpressionShiftDf(p.dist.info, sample.groups, return.raw=TRUE, ...)
  obs.diffs <- sapply(dist.per.type, mean, trim=trim)
  comp.res <- plapply(1:n.permutations, function(i) {
    sg.shuff <- sample.groups %>% {setNames(sample(.), names(.))}
    res <- estimateExpressionShiftDf(p.dist.info, sg.shuff, return.raw=TRUE)
    res[sapply(res, is.null)] <- NA

    sapply(res, mean, trim=trim) %>% .[names(obs.diffs)] %>%
      {. >= obs.diffs}
  }, progress=verbose, n.cores=n.cores, mc.preschedule=TRUE, fail.on.error=TRUE) %>%
    do.call(rbind, .)

  pvalues <- (colSums(comp.res, na.rm=TRUE) + 1) / (colSums(!is.na(comp.res)) + 1)

  return(list(dists=dist.per.type, pvalues=pvalues, obs.diffs=obs.diffs))
}

## Common shifts

##' @description calculate consensus change direction and distances between samples along this axis
consensusShiftDistances <- function(tcm, sample.groups, use.median=FALSE, mean.trim=0, use.cpp=TRUE) {
  if(min(table(sample.groups[colnames(tcm)])) < 1) return(NA); # not enough samples
  g1 <- which(sample.groups[colnames(tcm)]==levels(sample.groups)[1])
  g2 <- which(sample.groups[colnames(tcm)]==levels(sample.groups)[2])
  if (use.cpp)
    return(as.numeric(projdiff(tcm, g1 - 1, g2 - 1)))

  dm <- do.call(rbind, lapply(g1, function(n1) { # R
    do.call(rbind, lapply(g2, function(n2) {
      tcm[,n1] - tcm[,n2]
    }))
  }))

  if (use.median) {
    checkPackageInstalled("matrixStats", details="when `use.median=TRUE`", cran=TRUE)
    dmm <- matrixStats::colMedians(dm)
  } else {
    if (mean.trim > 0) {
      dmm <- apply(dm, 2, mean, trim=mean.trim)
    } else {
      dmm <- colMeans(dm)
    }
  }

  dmm <- dmm / sqrt(sum(dmm^2)) # normalize

  # project samples and calculate distances
  return(as.numeric(dm %*% dmm))
}

estimateExplainedVariance <- function(cm, sample.groups) {
  spg <- rownames(cm) %>% split(droplevels(as.factor(sample.groups[.])))
  if (length(spg) == 1) return(NULL)

  sapply(spg, function(spc) apply(cm[spc,,drop=FALSE], 2, var)) %>%
    rowMeans(na.rm=TRUE) %>% # TODO: add weights here
    {1 - . / apply(cm, 2, var)}
}

filterGenesForCellType <- function(cm, sample.groups, top.n.genes=500, selection=c("wilcox", "var"),
                                   exclude.genes=NULL) {
  selection <- match.arg(selection)

  if (selection == "var") {
    sel.genes <- estimateExplainedVariance(cm, sample.groups=sample.groups) %>%
      sort(decreasing=TRUE) %>% names()
  } else {
    spg <- rownames(cm) %>% split(sample.groups[.])
    cm <- cm / pmax(1, rowSums(cm))
    test.res <- matrixTests::col_wilcoxon_twosample(cm[spg[[1]],], cm[spg[[2]],])$pvalue
    sel.genes <- test.res %>% setNames(colnames(cm)) %>% sort() %>% names()
  }

  sel.genes %<>% setdiff(exclude.genes) %>% head(top.n.genes)
  return(sel.genes)
}
