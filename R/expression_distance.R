## Modified version of expression_distances.R to support blocking by co-variates

#' Expression shift magnitudes per cluster between conditions
#'
#' @param cm.per.type List of normalized count matrices per cell type
#'        (rows = samples, cols = genes/features; rownames are sample IDs)
#' @param sample.groups Named factor with sample IDs indicating condition (e.g., ctrl/disease)
#' @param cell.groups Named clustering/annotation factor with cell names
#' @param sample.per.cell Named vector mapping cell IDs -> sample ID
#' @param dist Distance measure: 'cor' (1 - Pearson r), 'l1' (Manhattan), 'l2' (Euclidean)
#' @param dist.type One of c("shift","total","var"); governs baseline & which pairs are used
#' @param pair.on Optional *named* factor (or data.frame of factors) keyed by sample IDs;
#'        only pairs matching all covariates are considered; permutations are blocked within these strata.
#' @param n.cores number of cores (default=1)
#' @param verbose (default=FALSE)
#' @param ref.level Reference level of `sample.groups` (required for dist.type="total" or "var")
#' @param n.permutations number of label permutations for p-value (default 1000)
#' @param p.adjust.method method for p.adjust across cell types (default "BH")
#' @param top.n.genes optional integer: select top-N genes before distances
#' @param gene.selection "wilcox" (default), "var", or "od"
#'        If `pair.on` is provided **and** gene.selection == "wilcox",
#'        a blocked Wilcoxon is run within each `pair.on` stratum and per-gene
#'        p-values are combined across strata by Fisher's method.
#' @param n.pcs optional dimensionality reduction by SVD before distances
#' @param trim trimmed-mean trim fraction (default 0.2)
#' @return List: dists.per.type (vector per cell type), p.dist.info (distance matrices),
#'         sample.groups, cell.groups, pvalues, padjust
#' @export
estimateExpressionChange <- function(cm.per.type, sample.groups, cell.groups, sample.per.cell,
                                     dist=NULL, dist.type=c("shift", "total", "var"), verbose=FALSE,
                                     ref.level=NULL, n.permutations=1000, p.adjust.method="BH",
                                     top.n.genes=NULL, gene.selection="wilcox", n.pcs=NULL,
                                     trim=0.2, n.cores=1, pair.on=NULL, ...) {
  dist.type <- match.arg(dist.type)
  dist <- parseDistance(dist, top.n.genes=top.n.genes, n.pcs=n.pcs)
  
  norm.type <- ifelse(dist.type == "shift", "both", "ref")
  r.type <- ifelse(dist.type == "var", "target", "cross")
  
  cell.groups %<>% as.factor() %>% droplevels()
  
  sample.type.table <- cell.groups %>% table(sample.per.cell[names(.)]) # table of sample types and cells
  
  if (verbose) message("Calculating pairwise distances using dist='", dist, "'...\n", sep="")
  
  n.cores.inner <- max(n.cores %/% length(levels(cell.groups)), 1)
  res.per.type <- levels(cell.groups) %>% sccore::sn() %>% plapply(function(ct) {
    cm.norm <- cm.per.type[[ct]]
    dist.mat <- estimateExpressionShiftsForCellType(
      cm.norm,
      sample.groups=sample.groups,
      dist=dist, n.pcs=n.pcs,
      top.n.genes=top.n.genes, gene.selection=gene.selection,
      pair.on=pair.on, ...
    )
    attr(dist.mat, 'n.cells') <- sample.type.table[ct, rownames(cm.norm)] # calculate how many cells there are
    
    # Build covariate-based pair mask (TRUE if pair matches all covariates)
    pair.mask <- .build_pair_mask(rownames(cm.norm), pair.on)
    
    dists <- estimateExpressionShiftsByDistMat(
      dist.mat, sample.groups,
      norm.type=norm.type, return.type=r.type,
      ref.level=ref.level, pair.mask=pair.mask
    )
    
    # If no eligible pairs exist under the mask, skip permutations
    if (length(dists) == 0) {
      warning("No eligible pairs under 'pair.on' for cell type ", ct)
      return(list(dists=numeric(0), dist.mat=dist.mat, pvalue=NA_real_))
    }
    
    obs.diff <- mean(dists, trim=trim)
    
    randomized.dists <- plapply(1:n.permutations, function(i) {
      sg.orig <- sample.groups[rownames(cm.norm)] %>% as.factor() %>% droplevels()
      
      # Permute within covariate strata (if provided)
      if (is.null(pair.on)) {
        sg.shuff <- setNames(sample(sg.orig), names(sg.orig))
      } else {
        block <- if (is.data.frame(pair.on)) interaction(pair.on[rownames(cm.norm)], drop=TRUE)
        else pair.on[rownames(cm.norm)]
        block <- as.factor(block)
        sg.shuff <- sg.orig
        for (b in levels(block)) {
          idx <- which(block == b)
          if (length(idx) > 1) sg.shuff[idx] <- sample(sg.shuff[idx])
        }
      }
      sg.shuff %<>% factor(levels=levels(sg.orig))
      
      dm <- dist.mat
      # If gene selection depends on labels, recompute under permutation (still respecting pair.on)
      if (!is.null(top.n.genes) && (gene.selection != "od")) {
        dm <- estimateExpressionShiftsForCellType(
          cm.norm, sample.groups=sg.shuff, dist=dist, n.pcs=n.pcs,
          top.n.genes=top.n.genes, gene.selection=gene.selection,
          pair.on=pair.on, ...
        )
      }
      
      estimateExpressionShiftsByDistMat(
        dm, sg.shuff, norm.type=norm.type,
        return.type=r.type, ref.level=ref.level, pair.mask=pair.mask
      ) %>% mean(trim=trim)
    }, progress=FALSE, n.cores=n.cores.inner, mc.preschedule=TRUE, fail.on.error=TRUE) %>% unlist()
    
    pvalue <- (sum(randomized.dists >= obs.diff) + 1) / (sum(!is.na(randomized.dists)) + 1)
    
    # Null-centering by the median null statistic
    dists <- dists - median(randomized.dists, na.rm=TRUE)
    dists.scaled <- dists/mad(randomized.dists, center = median(randomized.dists, na.rm=TRUE), constant = 1, na.rm=TRUE)
    
    list(dists=dists, dist.mat=dist.mat, pvalue=pvalue, dists.scaled=dists.scaled)
  }, progress=verbose, n.cores=n.cores, mc.preschedule=TRUE, mc.allow.recursive=TRUE, fail.on.error=TRUE)
  
  if (verbose) message("Done!\n")
  
  pvalues <- sapply(res.per.type, `[[`, "pvalue")
  dists.per.type <- lapply(res.per.type, `[[`, "dists")
  dists.per.type.scaled <- lapply(res.per.type, `[[`, "dists.scaled")
  p.dist.info <- lapply(res.per.type, `[[`, "dist.mat")
  
  padjust <- p.adjust(pvalues, method=p.adjust.method)
  
  return(list(dists.per.type=dists.per.type, p.dist.info=p.dist.info, sample.groups=sample.groups,
              cell.groups=cell.groups, pvalues=pvalues, padjust=padjust, dists.per.type.scaled=dists.per.type.scaled))
}


#' @keywords internal
estimateExpressionShiftsForCellType <- function(cm.norm, sample.groups, dist,
                                                top.n.genes=NULL, n.pcs=NULL,
                                                gene.selection="wilcox",
                                                exclude.genes=NULL,
                                                pair.on=NULL) {
  if (!is.null(top.n.genes)) {
    sel.genes <- filterGenesForCellType(
      cm.norm, sample.groups=sample.groups, top.n.genes=top.n.genes,
      gene.selection=gene.selection, exclude.genes=exclude.genes,
      pair.on=pair.on
    )
    cm.norm <- cm.norm[, sel.genes, drop=FALSE]
  }
  
  if (!is.null(n.pcs)) {
    min.dim <- min(dim(cm.norm)) - 1
    if (n.pcs > min.dim) {
      n.pcs <- min.dim
      warning("n.pcs is too large. Setting it to maximal allowed value ", min.dim)
    }
    pcs <- svd(cm.norm, nv=n.pcs, nu=0)
    cm.norm <- as.matrix(cm.norm %*% pcs$v)
  }
  
  if (dist == 'cor') {
    dist.mat <- 1 - cor(t(cm.norm))
  } else if (dist == 'l2') {
    dist.mat <- dist(cm.norm, method="euclidean") %>% as.matrix()
  } else if (dist == 'l1') {
    dist.mat <- dist(cm.norm, method="manhattan") %>% as.matrix()
  } else {
    stop("Unknown distance: ", dist)
  }
  
  dist.mat[is.na(dist.mat)] <- 1
  return(dist.mat)
}

#' @keywords internal
.build_pair_mask <- function(ids, pair.on) {
  if (is.null(pair.on)) return(NULL)
  cols <- if (is.data.frame(pair.on)) as.list(pair.on) else list(pair.on)
  cols <- lapply(cols, function(v) as.factor(v[ids]))
  mask <- Reduce("&", lapply(cols, function(v) outer(v, v, "==")))
  diag(mask) <- FALSE
  mask
}

#' @keywords internal
subsetDistanceMatrix <- function(dist.mat, sample.groups, cross.factor,
                                 pair.mask=NULL, build.df=FALSE) {
  comp.selector <- if (cross.factor) "!=" else "=="
  selection.mask <- outer(sample.groups[rownames(dist.mat)], sample.groups[colnames(dist.mat)], comp.selector)
  if (!is.null(pair.mask)) selection.mask <- selection.mask & pair.mask
  diag(dist.mat) <- NA
  if (!build.df) {
    return(na.omit(dist.mat[selection.mask]))
  }
  dist.mat[!selection.mask] <- NA
  if (all(is.na(dist.mat))) return(NULL)
  dist.df <- reshape2::melt(dist.mat) %>% na.omit()
  return(dist.df)
}

#' @keywords internal
estimateExpressionShiftsByDistMat <- function(dist.mat, sample.groups,
                                              norm.type=c("both", "ref", "none"),
                                              ref.level=NULL,
                                              return.type=c("cross", "target"),
                                              pair.mask=NULL) {
  norm.type <- match.arg(norm.type)
  return.type <- match.arg(return.type)
  
  if (((norm.type == "ref") || (return.type == "target")) && is.null(ref.level))
    stop("ref.level has to be provided for norm.type='ref' or return.type='target'")
  
  sample.groups %<>% .[rownames(dist.mat)]
  
  # --- Baseline normalization (respecting pair.mask) ---
  if (norm.type == "both") {
    sg1 <- levels(sample.groups)[1]
    m1 <- outer(sample.groups, sample.groups, function(a, b) (a == sg1) & (b == sg1))
    m2 <- outer(sample.groups, sample.groups, function(a, b) (a != sg1) & (b != sg1))
    if (!is.null(pair.mask)) { m1 <- m1 & pair.mask; m2 <- m2 & pair.mask }
    diag(m1) <- NA; diag(m2) <- NA
    norm.const <- (median(dist.mat[m1], na.rm=TRUE) + median(dist.mat[m2], na.rm=TRUE)) / 2
  } else if (norm.type == "ref") {
    mr <- outer(sample.groups, sample.groups, function(a, b) (a == ref.level) & (b == ref.level))
    if (!is.null(pair.mask)) mr <- mr & pair.mask
    diag(mr) <- NA
    norm.const <- median(dist.mat[mr], na.rm=TRUE)
  } else {
    norm.const <- 0
  }
  
  dist.mat <- dist.mat - norm.const
  
  # --- Select which pairs to return (respecting pair.mask) ---
  if (return.type == "cross") {
    dists <- subsetDistanceMatrix(dist.mat, sample.groups=sample.groups,
                                  cross.factor=TRUE, pair.mask=pair.mask)
  } else {
    # target-target, upper triangle only
    is.target <- sample.groups != ref.level
    tt <- outer(is.target, is.target, "&")
    if (!is.null(pair.mask)) tt <- tt & pair.mask
    diag(dist.mat) <- NA
    sel <- dist.mat
    sel[!tt] <- NA
    if (!any(!is.na(sel))) return(numeric(0))
    dists <- sel[upper.tri(sel)]
    dists <- dists[!is.na(dists)]
  }
  
  return(dists)
}

#' @keywords internal
joinExpressionShiftDfs <- function(dist.df.per.type, sample.groups) {
  dist.df.per.type %<>% .[!sapply(., is.null)]
  df <- names(dist.df.per.type) %>%
    lapply(function(n) cbind(dist.df.per.type[[n]], Type=n)) %>%
    do.call(rbind, .) %>%
    mutate(Condition=sample.groups[as.character(Var1)]) %>%
    na.omit()
  return(df)
}

#' @keywords internal
prepareJointExpressionDistance <- function(p.dist.per.type, sample.groups=NULL, return.dists=TRUE) {
  checkPackageInstalled(c("abind"), cran=TRUE)
  # bring to a common set of cell types
  common.types <- lapply(p.dist.per.type, colnames) %>% unlist() %>% unique()
  
  p.dist.per.type %<>% lapply(function(x) {
    y <- matrix(0,nrow=length(common.types),ncol=length(common.types))  # missing entries -> zero, carry zero weights
    rownames(y) <- colnames(y) <- common.types
    y[rownames(x),colnames(x)] <- x
    ycct <- setNames(rep(0,length(common.types)), common.types)
    ycct[colnames(x)] <- attr(x, 'n.cells')
    attr(y, 'n.cells') <- ycct
    y
  })
  
  x <- abind::abind(lapply(p.dist.per.type, function(x) {
    nc <- attr(x, 'n.cells')
    wm <- sqrt(outer(nc, nc, FUN='pmin'))
    return(x * wm)
  }), along = 3)
  
  # just the weights
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
  xd[!cross.factor] <- NA
  if (!any(!is.na(xd)))
    return(NULL)
  xmd2 <- na.omit(reshape2::melt(xd))
  xmd2 <- na.omit(xmd2)
  xmd2$type1 <- sample.groups[as.character(xmd2$Var1)]
  xmd2$type2 <- sample.groups[as.character(xmd2$Var2)]
  return(xmd2)
}

#' @keywords internal
filterCellTypesByNSamples <- function(cell.groups, sample.per.cell, sample.groups,
                                      min.cells.per.sample, min.samp.per.type, verbose=TRUE) {
  freq.table <- table(Type=cell.groups, Sample=sample.per.cell) %>% as.data.frame() %>%
    mutate(Condition=sample.groups[as.character(Sample)]) %>%
    filter(Freq >= min.cells.per.sample)
  
  if (length(unique(freq.table$Condition)) != 2)
    stop("'sample.groups' must be a 2-level factor describing which samples are being contrasted")
  
  removed.types <- freq.table %>% split(.$Type) %>% sapply(function(df) {
    df %$% split(Sample, Condition) %>% sapply(length) %>% {any(. < min.samp.per.type)}
  }) %>% which() %>% names()
  
  if (verbose && (length(removed.types) > 0)) {
    message("Excluding cell types ", paste(removed.types, collapse=", "), " that don't have enough samples\n")
  }
  
  freq.table %<>% filter(!(Type %in% removed.types))
  return(freq.table)
}

#' @keywords internal
filterExpressionDistanceInput <- function(cms, cell.groups, sample.per.cell, sample.groups,
                                          min.cells.per.sample=10, min.samp.per.type=2, min.gene.frac=0.01,
                                          genes=NULL, verbose=FALSE) {
  # Filter rare samples per cell type
  cell.names <- lapply(cms, rownames) %>% unlist()
  freq.table <- filterCellTypesByNSamples(
    cell.groups[cell.names], sample.per.cell[cell.names], sample.groups=sample.groups,
    min.cells.per.sample=min.cells.per.sample, min.samp.per.type=min.samp.per.type, verbose=verbose
  )
  
  filt.types.per.samp <- freq.table %$% split(Type, Sample)
  cms.filt <- names(filt.types.per.samp) %>% sn() %>% lapply(function(n) {
    cms[[n]] %>% .[cell.groups[rownames(.)] %in% filt.types.per.samp[[n]],, drop=FALSE]
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
  cms.filt %<>% lapply(sccore::extendMatrix, genes) %>% lapply(`[`,,genes, drop=FALSE)
  
  # Group matrices by cell type
  cell.groups <- droplevels(cell.groups[cell.names])
  
  cm.per.type <- levels(cell.groups) %>% sccore::sn() %>% lapply(function(ct) {
    lapply(cms.filt, function(x) if (ct %in% rownames(x)) x[ct,] else NULL) %>%
      do.call(rbind, .) %>% {. / pmax(1, rowSums(.))} %>% {log10(. * 1e3 + 1)}
  })
  
  return(list(cm.per.type=cm.per.type, cell.groups=cell.groups, sample.groups=sample.groups[names(cms)]))
}

#' @keywords internal
estimateExplainedVariance <- function(cm, sample.groups) {
  checkPackageInstalled("matrixStats", cran=TRUE)
  spg <- rownames(cm) %>% split(droplevels(as.factor(sample.groups[.])))
  if (length(spg) == 1){
    return(NULL)
  }
  sapply(spg, function(spc) matrixStats::colVars(cm[spc,,drop=FALSE]) * (length(spc) - 1)) %>%
    rowSums(na.rm=TRUE) %>%
    {1 - (. / (matrixStats::colVars(cm) * (nrow(cm) - 1)))} %>%
    setNames(colnames(cm))
}

#' @keywords internal
#' Gene selection with optional blocking by 'pair.on' using Fisher's method
filterGenesForCellType <- function(cm.norm, sample.groups, top.n.genes=500,
                                   gene.selection=c("wilcox", "var", "od"),
                                   exclude.genes=NULL, pair.on=NULL) {
  gene.selection <- match.arg(gene.selection)
  
  if (gene.selection == "var") {
    sel.genes <- estimateExplainedVariance(cm.norm, sample.groups=sample.groups) %>%
      sort(decreasing=TRUE) %>% names()
    
  } else if (gene.selection == "wilcox") {
    # If pair.on provided, do blocked Wilcoxon and Fisher combine across blocks
    if (!is.null(pair.on)) {
      pvals <- .blocked_wilcox_fisher(cm.norm, sample.groups, pair.on=pair.on)
      sel.genes <- setNames(pvals, colnames(cm.norm)) %>% sort() %>% names()
    } else {
      spg <- rownames(cm.norm) %>% split(sample.groups[.])
      test.res <- matrixTests::col_wilcoxon_twosample(
        cm.norm[spg[[1]],,drop=FALSE], cm.norm[spg[[2]],,drop=FALSE], exact=FALSE
      )$pvalue
      sel.genes <- test.res %>% setNames(colnames(cm.norm)) %>% sort() %>% names()
    }
    
  } else { # "od"
    checkPackageInstalled("pagoda2", details="for gene.selection='od'", cran=TRUE)
    p2 <- pagoda2::Pagoda2$new(t(cm.norm), modelType="raw", verbose=FALSE, n.cores=1)
    p2$adjustVariance(verbose=FALSE)
    sel.genes <- p2$getOdGenes(Inf)
  }
  
  sel.genes %<>% setdiff(exclude.genes) %>% head(top.n.genes)
  return(sel.genes)
}

#' @keywords internal
#' Blocked Wilcoxon per-gene with Fisher p-value combination across blocks
.blocked_wilcox_fisher <- function(cm.norm, sample.groups, pair.on) {
  checkPackageInstalled("matrixTests", cran=TRUE)
  # Align inputs
  ids <- rownames(cm.norm)
  groups <- as.factor(sample.groups[ids])
  if (length(levels(groups)) != 2)
    stop("Blocked Wilcoxon requires 'sample.groups' to have exactly 2 levels.")
  g1 <- levels(groups)[1]; g2 <- levels(groups)[2]
  
  block <- if (is.data.frame(pair.on)) interaction(pair.on[ids], drop=TRUE) else pair.on[ids]
  block <- as.factor(block)
  
  # Per-block p-values (vector per block; length = ncol(cm.norm))
  p_per_block <- lapply(levels(block), function(b) {
    i1 <- which(block == b & groups == g1)
    i2 <- which(block == b & groups == g2)
    if (length(i1) == 0 || length(i2) == 0) {
      rep(NA_real_, ncol(cm.norm))
    } else {
      matrixTests::col_wilcoxon_twosample(
        cm.norm[i1,,drop=FALSE], cm.norm[i2,,drop=FALSE], exact=FALSE
      )$pvalue %>% setNames(colnames(cm.norm))
    }
  })
  # Combine by Fisher across blocks, gene-wise
  p_mat <- do.call(cbind, p_per_block) # genes x blocks
  combine_fisher <- function(pv) {
    pv <- pv[is.finite(pv) & !is.na(pv) & pv > 0]
    if (length(pv) == 0) return(1)
    stat <- -2 * sum(log(pv))
    stats::pchisq(stat, df = 2 * length(pv), lower.tail = FALSE)
  }
  p_comb <- apply(p_mat, 1, combine_fisher)
  return(p_comb)
}

#' @keywords internal
parseDistance <- function(dist, top.n.genes, n.pcs) {
  n.comps <- min(top.n.genes, n.pcs, Inf)
  if (is.null(dist)) {
    dist <- ifelse(n.comps < 20, 'l1', 'cor')
    return(dist)
  }
  
  dist %<>% tolower()
  if (dist == 'l2') {
    warning("Using dist='l2' is not recommended, as it may introduce unwanted dependency ",
            "on the number of cells per cluster. Please, consider using 'l1' instead.")
  } else if (dist == 'cor') {
    if (n.comps < 20) {
      warning("dist='cor' is not recommended for data with dimensionality < 20. ",
              "Please, consider using 'l1' instead.")
    }
  } else if (dist == 'l1') {
    if (n.comps > 30) {
      warning("dist='l1' is not recommended for data with dimensionality > 30. ",
              "Please, consider using 'cor' instead.")
    }
  } else {
    stop("Unknown dist: ", dist)
  }
  return(dist)
}

#' @keywords internal
checkPackageInstalled <- function(pkgs, cran=FALSE, details=NULL) {
  pkgs <- as.character(pkgs)
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    msg <- paste0("Missing required package(s): ", paste(missing, collapse=", "))
    if (!is.null(details)) msg <- paste0(msg, ". ", details)
    stop(msg, call. = FALSE)
  }
}
