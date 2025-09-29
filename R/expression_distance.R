#' @importFrom sccore checkPackageInstalled
NULL

#' Expression shift magnitudes per cluster between conditions
#'
#' @param cm.per.type List of normalized count matrices per cell type
#' @param sample.meta Data frame with sample-level covariates, row names are sample IDs
#' @param sample.per.cell Named vector with cell names indicating sample/condition
#' @param formula Formula specifying the model to fit
#' @param contrast Character vector of length 3 specifying the contrast to test: c(var, alt, ref)
#' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
#' @param cell.groups Named clustering/annotation factor with cell names
#' @param dist what distance measure to use: 'JS' - Jensen-Shannon divergence (default), 'cor' - Pearson's linear correlation on log transformed values
#' @param dist.type one of "shift" (default), "total", or "var"
#' @param perm.method permutation method: "freedman-lane" (default) or "block"
#' @param block.id Optional column in `sample.meta` specifying blocks for restricted randomizations if perm.method="block"

#' @param n.cores number of cores (default=1)
#' @param verbose (default=FALSE)
#' @return List of distances per cell type, distance matrices, sample groups, cell types, pvalues, and adjusted p-values.
#' @export
estimateExpressionChange <- function(cm.per.type, sample.meta, cell.groups, sample.groups, sample.id = NULL,
                                     sample.per.cell, formula = NULL, contrast = NULL, block.id = NULL,
                                     dist = "cor", dist.type = c("shift", "total", "var"), robust.method = c("none", "huber", "winsor"),
                                     gene.selection = "wilcox", perm.method = c("freedman-lane", "block"), na.mode = c("drop", "impute_weak"),
                                     n.permutations = 1000, p.adjust.method = "BH", trim = 0.2, return.residuals = FALSE, return.sampled.stats = TRUE,
                                     top.n.genes = NULL, n.pcs = NULL, n.cores = 1, verbose = TRUE, ...) {
  dist.type <- match.arg(dist.type)
  dist <- parseDistance(dist, top.n.genes = top.n.genes, n.pcs = n.pcs)
  perm.method <- match.arg(perm.method)
  robust.method <- match.arg(robust.method)
  na.mode <- match.arg(na.mode)

  cell.groups <- droplevels(factor(cell.groups))
  sample.type.table <- table(cell.groups, sample.per.cell[names(cell.groups)])

  contrast.var <- contrast[1]

  covariates <- all.vars(formula)
  sample.meta.df <- sample.meta[, covariates, drop = FALSE]
  sample.meta.df[[contrast.var]] <- factor(sample.meta.df[[contrast.var]])
  if (ncol(sample.meta.df) == 0) stop("Covariate frame is empty after selecting terms from the formula.")

  if (verbose) message("Fitting LM with formula: ", deparse(formula))

  n.cores.inner <- max(floor(n.cores / max(1L, length(levels(cell.groups)))), 1)

  res.per.type <- levels(cell.groups) %>%
                          sccore::sn() %>%
                            plapply(function(ct) {
                              if (verbose) message("Processing cell type '", ct, "'...")
                                cm.norm <- cm.per.type[[ct]]
                                

                                ## distance matrix
                                dist.mat <- estimateExpressionShiftsForCellType(
                                 cm.norm, sample.groups = sample.groups, dist = dist, n.pcs = n.pcs,
                                 top.n.genes = top.n.genes, gene.selection = gene.selection, meta= sample.meta.df[rownames(cm.norm), , drop = FALSE],
                                 formula = formula, sample.id = sample.id, contrast = contrast)
                                attr(dist.mat, "n.cells") <- sample.type.table[ct, rownames(cm.norm)]

                                dists = subsetDistanceMatrix(dist.mat, sample.groups, cross.factor=TRUE, build.df=FALSE)
                                cross.rows <- which(outer(rownames(dist.mat), rownames(dist.mat), function(a,b) sample.groups[a] != sample.groups[b]))

                                # model matrices and pairwise distances
                                sample.meta.df.ct <- sample.meta.df[rownames(dist.mat), , drop = FALSE]
                                px <- buildPairDesignMatrices(sample.meta = sample.meta.df.ct,
                                       triplet = c(contrast[1],contrast[2],contrast[3]),  # var, alt, ref
                                       dist.type = "shift", block.vars = ifelse(!is.null(block.id), paste0(block.id, "_pair"), NULL))

                                y <- vectorizeLowerTri(dist.mat, px$pairs)

                                # Linear models and Permutations
                                res <- performLMPermutations(y=y, x=px$model, n.permutations=n.permutations, perm.method=perm.method, robust.method = robust.method,
                                                             na.mode = na.mode, core.only=TRUE, cross.rows= cross.rows, return.residuals = return.residuals, 
                                                             return.sampled.stats = return.sampled.stats)

                                # R^2 for each term
                                groups <- makeGroupsPair(px$model$F)
                                r2.df <- estimateR2PerTerm(px$model$F, y, groups)
                                r2.df <- r2.df[-grepl("Intercept", r2.df$term), ]

                                if(return.residuals){ residuals <- res$y.resid } else { residuals <- NULL }

                                list(res = res, dist.mat = dist.mat, dists = dists, r2 = r2.df, model.diag=px$model.diag, residuals = residuals)
                 }, progress = verbose, n.cores = n.cores.inner, mc.preschedule = TRUE, mc.allow.recursive = TRUE, fail.on.error = TRUE)

  if (verbose) message("Done!\n")
  
  summary <- summarizeLMResults(res.per.type, adj.method=p.adjust.method)

  ret <- list(results= summary$results, pvalues= summary$pvalues, padjust= summary$padjust, zscores= summary$zscores, perm.method = perm.method, 
              dist.type = dist.type, dists.per.type = summary$dists.per.type, p.dist.info = summary$p.dist.info, r2 = summary$r2, cell.groups=cell.groups, 
              model.diagnostics = lapply(res.per.type, `[[`, "model.diag"), robust.method= robust.method)
  if (return.residuals) ret$residuals <- summary$residuals
  if (return.sampled.stats) ret$permutation.stats <- summary$stats.perm
  ret
}


#' @keywords internal
estimateExpressionShiftsForCellType <- function(cm.norm, sample.groups, dist, top.n.genes=NULL, n.pcs=NULL,
                                                gene.selection="wilcox", exclude.genes=NULL, meta= sample.meta, 
                                                formula = formula, sample.id = sample.id, contrast = contrast) {
  if (!is.null(top.n.genes)) {
    sel.genes <- filterGenesForCellType(
      cm.norm, sample.groups=sample.groups, top.n.genes=top.n.genes, gene.selection=gene.selection,
      exclude.genes=exclude.genes, meta= meta, formula = formula, sample.id = sample.id, contrast = contrast
    )
    cm.norm <- cm.norm[,sel.genes,drop=FALSE]
  }

  if (!is.null(n.pcs)) {
    min.dim <- min(dim(cm.norm)) - 1
    if (n.pcs > min.dim) {
      n.pcs <- min.dim
      warning("n.pcs is too large. Setting it to maximal allowed value ", min.dim)
    }

  cm.norm <- getTopPCs(cm.norm, n.pcs = n.pcs)

  if (any(!is.finite(cm.norm))) {
    stop("PCA output contains non-finite values (NA/NaN/Inf).")
  }
  zero.var.cols <- which(apply(cm.norm, 2, var) < .Machine$double.eps)
  if (length(zero.var.cols) > 0) {
    warning("PCA output contains zero-variance components; removing them.")
    cm.norm <- cm.norm[, -zero.var.cols, drop = FALSE]
  }
  if (ncol(cm.norm) == 0) {
    stop("No valid principal components remain after filtering zero-variance columns.")
    }
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

  dist.mat[is.na(dist.mat)] <- 1;
  return(dist.mat)
}



#' @keywords internal
subsetDistanceMatrix <- function(dist.mat, sample.groups, cross.factor, build.df=FALSE) {
  comp.selector <- if (cross.factor) "!=" else "=="
  selection.mask <- outer(sample.groups[rownames(dist.mat)], sample.groups[colnames(dist.mat)], comp.selector);
  diag(dist.mat) <- NA;
  if (!build.df)
    return(na.omit(dist.mat[selection.mask]))

  dist.mat[!selection.mask] <- NA;

  if(all(is.na(dist.mat))) return(NULL);
  dist.df <- reshape2::melt(dist.mat) %>% na.omit()
  return(dist.df);
  return(na.omit(dist.mat[selection.mask]))
}


###################### Expression distance helper functions #########################

# Vectorize the lower triangle of a square matrix using a precomputed index order
# If pairs is NULL, use standard lower.tri order for that matrix.
#' @keywords internal
vectorizeLowerTri <- function(d, pairs = NULL) {
  stopifnot(is.matrix(d), nrow(d) == ncol(d))
  if (is.null(pairs)) {
    ij <- which(lower.tri(d), arr.ind = TRUE)
    return(d[cbind(ij[,1], ij[,2])])
  } else {
    stopifnot(all(c("i","j") %in% names(pairs)),
              nrow(d) >= max(pairs$i, pairs$j))
    return(d[cbind(pairs$i, pairs$j)])
  }
}

###################### Expression distance helper functions #########################

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
filterExpressionDistanceInput <- function(cms, cell.groups, sample.per.cell, sample.groups, keep.all=FALSE,
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
filterGenesForCellType <- function(cm.norm, meta, sample.groups, top.n.genes=500, gene.selection=c("wilcox", "var", "od", "deseq2"),
                                   exclude.genes=NULL, formula = NULL, sample.id = NULL, contrast = NULL) {
  gene.selection <- match.arg(gene.selection)

  if (gene.selection == "var") {
    sel.genes <- estimateExplainedVariance(cm.norm, sample.groups=sample.groups) %>%
      sort(decreasing=TRUE) %>% names()
  } else if (gene.selection == "wilcox") {
    spg <- rownames(cm.norm) %>% split(sample.groups[.])
    test.res <- matrixTests::col_wilcoxon_twosample(cm.norm[spg[[1]],,drop=FALSE], cm.norm[spg[[2]],,drop=FALSE], exact=FALSE)$pvalue
    sel.genes <- test.res %>% setNames(colnames(cm.norm)) %>% sort() %>% names()
  } else if (gene.selection == "deseq2") {
    checkPackageInstalled("DESeq2", details="for gene.selection='deseq2'", cran=TRUE)
    # Sample size filtering 
    n.samples <- table(meta[[contrast[1]]])
    if (any(n.samples < 2)) {
      stop("Each group must be present in at least two samples. Change gene.selection method or filter cell types with too few samples.")
    }
    # drop formula terms with insufficient samples
    terms.in.formula <- all.vars(stats::terms(formula))
    terms.in.formula <- setdiff(terms.in.formula, sample.id)
    for (v in terms.in.formula) {
      if (is.character(meta[[v]]) || is.logical(meta[[v]])) {
        meta[[v]] <- factor(meta[[v]])
      }
    }
    valid.terms <- Filter(function(v) {
      x <- meta[[v]]
      if (is.factor(x)) nlevels(droplevels(x)) >= 2 else TRUE
    }, terms.in.formula)
    
    formula.inner <- stats::reformulate(valid.terms)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = cm.norm, colData = meta, design = formula.inner)
    DESeq2::DESeq(dds, quiet = TRUE, test = 'Wald')
    sel.genes <- rownames(dds)[order(DESeq2::results(dds)$padj)]
  } else {
    checkPackageInstalled("pagoda2", details="for gene.selection='od'", cran=TRUE)
    # TODO: we need to extract the OD function from Pagoda and scITD into sccore
    # Pagoda2 should not be in DESCRIPTION
    p2 <- pagoda2::Pagoda2$new(t(cm.norm), modelType="raw", verbose=FALSE, n.cores=1)
    p2$adjustVariance(verbose=FALSE)
    sel.genes <- p2$getOdGenes(Inf)
  }



  sel.genes %<>% setdiff(exclude.genes) %>% head(top.n.genes)
  return(sel.genes)
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
getTopPCs <- function(cm.norm, n.pcs) {
  samp.names <- rownames(cm.norm)

  if (!is.matrix(cm.norm)) {
    cm.norm <- as.matrix(cm.norm)
  }
  if (!is.numeric(cm.norm)) {
    stop("Input matrix must be numeric.")
  }
  cm.norm.p <- pca_project(cm.norm, n.pcs)
  rownames(cm.norm.p) <- samp.names
  return(cm.norm.p)
}
