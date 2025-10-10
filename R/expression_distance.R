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
estimateExpressionChange <- function(cm.per.type, sample.meta, cell.groups, design.mat, sample.per.cell,
                                     top.n.genes = NULL, cm.raw.per.type = NULL, n.pcs = NULL, sample.id = NULL, formula = NULL, 
                                     dist = "cor", dist.type = c("shift", "total", "var"), contrast = NULL,
                                     perm.method = c("freedman-lane", "block"), robust.method = c("none", "huber", "winsor"),
                                     gene.selection = "deseq2", na.mode = c("drop", "impute_weak"), alternative = c("two-sided", "greater", "less"),
                                     n.permutations = 1000, p.adjust.method = "BH", trim = 0.2, return.residuals = FALSE, 
                                     return.sampled.stats = TRUE, return.sampled.fits = FALSE, 
                                     n.cores = 1, verbose = TRUE, ...) {
  dist.type <- match.arg(dist.type)
  dist <- parseDistance(dist, top.n.genes = top.n.genes, n.pcs = n.pcs)
  perm.method <- match.arg(perm.method)
  robust.method <- match.arg(robust.method)
  na.mode <- match.arg(na.mode)

  cell.groups <- droplevels(factor(cell.groups))
  sample.type.table <- table(cell.groups, sample.per.cell[names(cell.groups)])
  #n.cores.inner <- max(floor(n.cores / max(1L, length(levels(cell.groups)))), 1)

  # Pairwise Distances
  p.dist <- estimateExpressionShiftsForCellType(cm.per.type, dist = dist, top.n.genes = top.n.genes, cm.raw = cm.raw.per.type,
                                                gene.selection = gene.selection, n.pcs = n.pcs, meta = sample.meta,
                                                sample.id = sample.id, formula=formula, contrast = contrast, cross.only=TRUE, 
                                                core.rows=which(design.mat$model$core.rows), ...)
  # Fitting and randomization
  res <- performLMPermutations(y=p.dist$Y, x=design.mat$model, n.permutations=n.permutations, perm.method=perm.method, 
                               robust.method = robust.method, na.mode = na.mode, cross.rows= cross.rows, alternative = alternative,
                               return.residuals = return.residuals,return.sampled.stats = return.sampled.stats,
                               return.sampled.fits = return.sampled.fits, n.cores = n.cores, ...)
  # R2 estimation
  r2 <- estimateR2PerTerm(design.mat$model$F, p.dist$Y, groups = makeGroupsPair(design.mat$model$F), idx = design.mat$pairs)

  #model.diagnostics <- lapply(res.per.type, `[[`, "model.diag") // TODO

  if (verbose) message("Done!\n")
  ct <- levels(cell.groups)
  names(res$stat.obs) <- names(res$pval) <- names(res$z.score) <- ct
  summary  <- data.frame(celltype=ct, obs.stat=res$stat.obs, pvalue=res$pval, zscore=res$z.score, stringsAsFactors=FALSE)
  summary$padjust <- p.adjust(summary$pvalue, method=p.adjust.method)
  summary <- summary[order(summary$obs.stat, decreasing=TRUE), ]

  out <- list(results= summary, res = res, perm.method = perm.method, robust.method= robust.method,
              dist.type = dist.type, dists.per.type = p.dist$dists, p.dist = p.dist, r2 = r2, sample.table = sample.type.table
              )
  out
}

#' @keywords internal
estimateExpressionShiftsForCellType <- function(
  cm.norm, dist, top.n.genes=NULL, cm.raw=NULL, n.pcs=NULL,
  gene.selection="wilcox", exclude.genes=NULL, meta = NULL, 
  formula = NULL, sample.id = NULL, contrast = NULL,
  core.rows = NULL, cross.only = FALSE # keep only alt/ref cross pairs if TRUE
) {
  if (is.list(cm.norm)) {
    stopifnot(is.data.frame(meta))
    n <- nrow(meta)
    idx <- lowerTriIndices(n)
    samp.names <- if (!is.null(sample.id) && sample.id %in% names(meta)) {
      as.character(meta[[sample.id]])
    } else if (!is.null(rownames(meta)) && all(nzchar(rownames(meta)))) {
      rownames(meta)
    } else {
      NULL
    }
    pair.names <- if (!is.null(samp.names)) {
      paste(samp.names[idx$i], samp.names[idx$j], sep="__")
    } else {
      paste(idx$i, idx$j, sep="__")
    }

    ctnames <- names(cm.norm)
    if (is.null(ctnames)) ctnames <- paste0("CT", seq_along(cm.norm))

    # Reorder rows to match meta row order 
    rn <- rownames(meta)
    for (k in seq_along(cm.norm)) {
      if (!is.null(rownames(cm.norm[[k]])) && !is.null(rn) &&
          setequal(rownames(cm.norm[[k]]), rn)) {
        cm.norm[[k]] <- cm.norm[[k]][rn, , drop=FALSE]
      }
    }

    # Optional: keep only cross-group (alt/ref) pairs
    keep.rows <- rep(TRUE, nrow(idx))
    if (cross.only) {
      tr <- parseTriplet(contrast) 
      g  <- factor(meta[[tr$var]])
      gi <- as.character(g[idx$i]); gj <- as.character(g[idx$j])
      keep.rows <- ( (gi == tr$ref & gj == tr$alt) | (gi == tr$alt & gj == tr$ref) )
    }

    # Allocate Y: rows = sample pairs (design order), cols = cell types
    Y <- matrix(NA_real_, nrow = nrow(idx), ncol = length(cm.norm),
                dimnames = list(pair.names, ctnames))

    for (t in seq_along(cm.norm)) {
      X <- cm.norm[[t]]

      # (1) Optional gene selection
      if (!is.null(top.n.genes)) {
        if(gene.selection == "deseq2"){
          if (is.null(formula) || is.null(cm.raw)) {
          stop("Design formula and raw counts must be provided for gene.selection='deseq2'")
        }
        X.use <- cm.raw[[t]]  # use raw counts for DESeq2
        } else {
          X.use <- X
        }
        sel.genes <- filterGenesForCellType(
          X.use, top.n.genes=top.n.genes, gene.selection=gene.selection,
          exclude.genes=exclude.genes, meta=meta, formula=formula,
          sample.id=sample.id, contrast=contrast
        )
        X <- X[, sel.genes, drop=FALSE]
      }

      # (2) PCA (skip if any NA to preserve NA semantics)
      Xp <- X
      if (!is.null(n.pcs) && !anyNA(X)) {
        min.dim <- min(dim(X)) - 1
        if (n.pcs > min.dim) {
          n.pcs <- min.dim
          warning("n.pcs is too large for type '", ctnames[t], "'. Setting to ", min.dim)
        }
        Xp <- getTopPCs(X, n.pcs = n.pcs)

        if (any(!is.finite(Xp))) stop("PCA output contains non-finite values (NA/NaN/Inf).")
        zero.var.cols <- which(apply(Xp, 2, var) < .Machine$double.eps)
        if (length(zero.var.cols) > 0) {
          warning("PCA output has zero-variance components for '", ctnames[t], "'; removing.")
          Xp <- Xp[, -zero.var.cols, drop = FALSE]
        }
        if (ncol(Xp) == 0) stop("No valid principal components remain after filtering zero-variance columns.")
      } else if (!is.null(n.pcs) && anyNA(X)) {
        warning("Type '", ctnames[t], "': data contain NA; skipping PCA to preserve NA-pair semantics.")
      }

      # (3) Distances — produce a square *base* matrix D
      D <- NULL
      if (dist == 'cor') {
        Xc <- as.matrix(Xp)                 # ensure base matrix
        n_rows <- nrow(Xc)
        if (is.null(n_rows) || n_rows < 2L) {
          D <- matrix(0, n_rows, n_rows, dimnames = list(rownames(Xc), rownames(Xc)))
        } else {
          R <- stats::cor(t(Xc), use = "pairwise.complete.obs", method = "pearson")
          D <- 1 - R
        }
      } else if (dist %in% c('l2','l1')) {
        Xd <- as.matrix(Xp)                 # ensure base matrix
        n_rows <- nrow(Xd)
        if (n_rows < 2L) {
          D <- matrix(0, n_rows, n_rows, dimnames = list(rownames(Xd), rownames(Xd)))
        } else if (anyNA(Xd)) {
          D <- matrix(NA_real_, n_rows, n_rows, dimnames = list(rownames(Xd), rownames(Xd)))
          diag(D) <- 0
          for (i in seq_len(n_rows-1L)) {
            xi <- Xd[i,]
            for (j in (i+1L):n_rows) {
              xj <- Xd[j,]
              if (anyNA(xi) || anyNA(xj)) {
                d <- NA_real_
              } else if (dist == 'l2') {
                d <- sqrt(sum((xi - xj)^2))
              } else {
                d <- sum(abs(xi - xj))
              }
              D[i,j] <- D[j,i] <- d
            }
          }
        } else {
          method <- if (dist == 'l2') "euclidean" else "manhattan"
          D <- as.matrix(stats::dist(Xd, method = method))
          # set dimnames if lost
          if (is.null(rownames(D)) && !is.null(rownames(Xd))) {
            rownames(D) <- colnames(D) <- rownames(Xd)
          }
        }
      } else {
        stop("Unknown distance: ", dist)
      }

      # Harden: ensure D is a square base matrix
      if (!is.matrix(D)) D <- as.matrix(D)
      dimD <- dim(D)
      if (length(dimD) != 2L || dimD[1] != dimD[2]) {
        stop("Distance is not a square matrix for cell type '", ctnames[t], "': got ", paste(dimD, collapse = "x"))
      }

      # (4) Vectorize using the same pair order as the design (idx)
      Y[, t] <- vectorizeLowerTri(D, pairs = idx)
    }
    
    # rows to show in dists (visualization)
    if (is.null(core.rows)) {
      core.rows <- seq_len(nrow(idx))
    } else {
      core.rows <- as.integer(core.rows)
      stopifnot(all(core.rows >= 1L & core.rows <= nrow(idx)))
    }
    rows.use <- if (cross.only) keep.rows else core.rows
    dists <- Y[rows.use, , drop = FALSE]

    return(list(Y = Y, dists = dists))

  } else {
    # ---- single-matrix path (unchanged, but coerce for cor) ----
    if (!is.null(top.n.genes)) {
      if(gene.selection == "deseq2"){
        if(is.null(formula) || is.null(cm.raw)) {
          stop("Design formula and raw counts must be provided for gene.selection='deseq2'")
        }
        X.use <- cm.raw  # use raw counts for DESeq2
      } else {
        X.use <- cm.norm
      }
      sel.genes <- filterGenesForCellType(
        X.use, top.n.genes=top.n.genes, gene.selection=gene.selection,
        exclude.genes=exclude.genes, meta= meta, formula = formula,
        sample.id = sample.id, contrast = contrast
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
      dist.mat <- 1 - stats::cor(t(as.matrix(cm.norm)), use = "pairwise.complete.obs")
    } else if (dist == 'l2') {
      dist.mat <- as.matrix(stats::dist(cm.norm, method="euclidean"))
    } else if (dist == 'l1') {
      dist.mat <- as.matrix(stats::dist(cm.norm, method="manhattan"))
    } else {
      stop("Unknown distance: ", dist)
    }
    return(dist.mat)
  }
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
filterExpressionDistanceInput <- function(
  cms, cell.groups, sample.per.cell, sample.meta, sample.id, keep.all=FALSE,
  min.cells.per.sample=10, min.samp.per.type=2, min.gene.frac=0.01,
  genes=NULL, verbose=FALSE,
  gene.selection=c("wilcox","var","od","deseq2")
) {
  gene.selection <- match.arg(gene.selection)

  stopifnot(is.list(cms), length(cms) > 0)
  all.samples <- names(cms)

  if (!keep.all) {
    sample.groups <- getSampleGroups(sample.meta, sample.id=sample.id)
    cell.names <- lapply(cms, rownames) %>% unlist()
    freq.table <- filterCellTypesByNSamples(
      cell.groups[cell.names], sample.per.cell[cell.names], sample.groups=sample.groups,
      min.cells.per.sample=min.cells.per.sample, min.samp.per.type=min.samp.per.type, verbose=verbose
    )
    filt.types.per.samp <- freq.table %$% split(Type, Sample)
    cms.filt <- names(filt.types.per.samp) %>% sn() %>% lapply(function(n) {
      cms[[n]] %>% .[cell.groups[rownames(.)] %in% filt.types.per.samp[[n]], , drop=FALSE]
    })
    names(cms.filt) <- names(filt.types.per.samp)
  } else {
    cms.filt <- cms
  }

  cell.names <- lapply(cms.filt, rownames) %>% unlist()

  # ---- Gene filtering across samples ----
  if (is.null(genes)) {
    genes <- lapply(cms.filt, function(cm) {
      cm@x <- 1 * (cm@x > 1)
      names(which(colMeans(cm) > min.gene.frac))
    }) %>%
      unlist() %>% table() %>% {. / length(cms.filt) > 0.1} %>% which() %>% names()
  }

  # ---- Collapse by cell type and align genes ----
  cms.filt %<>% lapply(collapseCellsByType, groups=cell.groups, min.cell.count=1)
  cms.filt %<>% lapply(sccore::extendMatrix, genes) %>% lapply(`[`, , genes, drop=FALSE)

  # ---- Build per-cell-type matrices (always make RAW + NORM) ----
  cell.groups <- droplevels(cell.groups[cell.names])
  all.types <- levels(cell.groups)

  cms.filt <- cms.filt[all.samples]

  cm.raw.per.type <- sccore::sn(all.types) %>% lapply(function(ct) {
    # missing -> zero row (sparse)
    rows_raw <- lapply(cms.filt, function(x) {
      if (ct %in% rownames(x)) {
        x[ct, , drop=FALSE]
      } else {
        Matrix::Matrix(0, nrow=1, ncol=ncol(x),
                       dimnames=list(ct, colnames(x)), sparse=TRUE)  # keep 'double' slot
      }
    })
    mat_raw <- do.call(rbind, rows_raw)                # samples x genes
    # ensure integer-valued while keeping dgCMatrix 'double' storage
    if (inherits(mat_raw, "dgCMatrix")) {
      if (!all(abs(mat_raw@x - round(mat_raw@x)) < 1e-8)) {
        stop("Counts for DESeq2 are not whole numbers; ensure cms are raw counts.")
      }
      mat_raw@x <- round(mat_raw@x)
    } else {
      if (!all(abs(mat_raw - round(mat_raw)) < 1e-8, na.rm = TRUE)) {
        stop("Counts for DESeq2 are not whole numbers; ensure cms are raw counts.")
      }
      mat_raw <- round(mat_raw)
    }
    rownames(mat_raw) <- names(cms.filt)
    mat_raw
  })
  names(cm.raw.per.type) <- all.types

  cm.per.type <- sccore::sn(all.types) %>% lapply(function(ct) {
    # norm: row present -> take it; missing -> NA row (as before)
    rows_norm <- lapply(cms.filt, function(x) {
      if (ct %in% rownames(x)) {
        x[ct, , drop=FALSE]
      } else {
        m <- matrix(NA_real_, nrow=1, ncol=ncol(x), dimnames=list(ct, colnames(x)))
        as(m, "dgCMatrix")
      }
    })
    mat <- do.call(rbind, rows_norm)                   # samples x genes
    denom <- pmax(1, rowSums(mat, na.rm=TRUE))
    mat <- mat / denom
    mat <- log10(mat * 1e3 + 1)
    rownames(mat) <- names(cms.filt)
    mat
  })
  names(cm.per.type) <- all.types

  return(list(
    cm.per.type      = cm.per.type,       # normalized (use everywhere downstream)
    cm.raw.per.type  = cm.raw.per.type,   # raw (use only for DESeq2 gene selection)
    cell.groups      = cell.groups
  ))
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
filterGenesForCellType <- function(cm.norm, meta, top.n.genes=500, gene.selection=c("wilcox", "var", "od", "deseq2"),
                                   exclude.genes=NULL, formula = NULL, sample.id = NULL, contrast = NULL) {
  gene.selection <- match.arg(gene.selection)
  sample.groups <- getSampleGroups(meta, contrast = contrast, sample.id = sample.id)

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

    # align samples 
    if (is.null(rownames(meta))) stop("meta must have rownames as sample IDs.")
    common <- intersect(rownames(meta), rownames(cm.norm))
    if (length(common) < 2L) stop("Too few overlapping samples between meta and counts.")
    #cm.norm  <- cm.norm[common, , drop=FALSE]
    meta <- meta    [common, , drop=FALSE]
    
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
    counts <- t(cm.norm)       # genes x samples (currently numeric)
    counts <- round(as.matrix(counts))
    storage.mode(counts) <- "integer"

    meta <- meta[colnames(counts), , drop = FALSE] # reorder samples
    stopifnot(identical(colnames(counts), rownames(meta)))
    # run deseq2
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta, design = formula.inner)
    dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE, test = 'Wald'))
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

#' Extract Partial or Fitted Expression Shifts for Visualization. "partial" corresponds to the distance explained
#'   by the contrast of interest only (i.e., after regressing out nuisance covariates).
#' @param p.dist Output of \code{estimateExpressionShiftsForCellType} function.
#' @param res Output of \code{performLMPermutations} function.
#' @param design.mat Output of \code{buildDesignMatrix} function. (or pairwise design matrices)
#' @param contrast Contrast triplet (variable, ref, alt) used to define the groups.
#' @param sample.meta Sample metadata data frame used to build the design matrix.
#' @param sample.id Column name in \code{sample.meta} that contains sample IDs
#'   matching the row names of the input cell matrices.
#' @param cov.plot.keys Character vector of column names in \code{design.mat$pair.meta}
#'   to use for coding covariate patterns in the output data frame.
#' @param block.vars Optional character vector of column names in \code{design.mat$pair.meta}
#'   to use for block-wise summary of shifts (e.g., batch or other grouping variable).
#'   If provided, the function will compute mean shifts within each block and
#'   return an additional data frame with this summary.
#' @return A list with two data frames:
#'   \item{df.cov.keys}{Data frame suitable for plotting shifts with covariate patterns.}
#'   \item{df.blocks}{(optional) Data frame with block-wise summary of shifts.}
#'   \item{df.perm}{(optional) Data frame with permutation statistics for background distribution.}
#' @details This function extracts expression shifts from the results of
#'   \code{performLMPermutations} and prepares data frames for visualization.
#'   It computes partial shifts (explained by the contrast of interest only), centers them
#'   by the mean of ref-ref (AA) pairs and/or alt-alt (BB) pairs, and codes covariate patterns
#'   for plotting. If \code{block.vars} is provided, it also computes mean shifts
#'   within each block defined by the specified variables.
#' @keywords internal
extractPairwiseShifts <- function(res, p.dist, design.mat, dist.type, perm.method= c("freedman-lane", "block"),
                           sample.meta, sample.id, contrast, block.vars = NULL, cov.plot.keys = NULL) {
  perm.method <- match.arg(perm.method)
  dist.type <- match.arg(dist.type, c("shift", "total"))
  if (perm.method == "block") {
    blk <- extractFitsBlock(res, p.dist, design.model = design.mat)
    partial.fit <- blk$partial
  } else if (perm.method == "freedman-lane") {
    partial.fit <- extractFitsFL(res, p.dist, design.mat)
  } else {
    stop("Unknown model type: ", perm.method)
  }

  delta <- extractContrastDeltas(yhat.all = partial.fit, pair.meta = design.mat$pair.meta, pairs = design.mat$pairs,
                                 sample.meta = sample.meta, sample.id = sample.id, contrast = contrast,
                                 center.by = dist.type, block.vars = block.vars)

  ## --- get covariate patterns ---
  ab.pos <- alignABrows(p.dist$dists, design.mat$pairs, sample.meta, sample.id = sample.id)
  pm.ab <- design.mat$pair.meta[ab.pos, , drop = FALSE]
  pattern.keys <- paste0(cov.plot.keys, "_pair")
  keys <- intersect(pattern.keys, colnames(pm.ab))
  if (!length(keys)) warning("No key columns found in pair.meta. Check cov.plot.keys.")
  patt.coded <- if (length(keys)) do.call(paste, c(pm.ab[, keys, drop = FALSE], sep = " | "))
                else rep("(no key cols in pair.meta)", nrow(pm.ab))

  ## --- panel I ---
  df.shifts <- as.data.frame(delta$shifts) |>
    tibble::rownames_to_column("pair") |>
    tidyr::pivot_longer(-pair, names_to = "celltype", values_to = "shifts") |>
    dplyr::filter(is.finite(shifts))

  pattern.by.pair <- setNames(patt.coded, rownames(delta$shifts))
  df.cov.keys <- df.shifts
  df.cov.keys$pattern_coded <- factor(pattern.by.pair[df.cov.keys$pair])

  ## --- panel II: block-wise summary (optional) ---
  df.blocks <- delta$delta.block

  ## --- panel III: permutation background + observed ---
  stats.perm <- if (!is.null(res$stats.perm)) res$stats.perm else NULL
  stat.obs   <- if (!is.null(res$stat.obs))   res$stat.obs   else NULL
  df.perm     <- NULL
  df.perm.obs <- NULL
  if (!is.null(stats.perm)) {
    ct <- colnames(stats.perm)
    df.perm <- as.data.frame(stats.perm) |>
      tidyr::pivot_longer(tidyselect::all_of(ct), names_to = "celltype", values_to = "perm_stat") |>
      dplyr::filter(is.finite(perm_stat))
    if (!is.null(stat.obs)) {
      df.perm.obs <- data.frame(celltype = ct, obs = as.numeric(stat.obs), stringsAsFactors = FALSE)
    }
  }
  out <- list(df.shifts = df.shifts, df.cov.keys = df.cov.keys)
  if (!is.null(df.blocks))    out$df.blocks    <- df.blocks
  if (!is.null(df.perm))      out$df.perm      <- df.perm
  if (!is.null(df.perm.obs))  out$df.perm.obs  <- df.perm.obs
  return(out)
}

#' Extract Contrast-Specific Expression Shifts. // TODO modify to use the contrast vector instead
#' @keywords internal
extractContrastDeltas <- function(yhat.all,      # (all.pairs × celltypes),
                 pair.meta,                        # design.mat$pair.meta (rows in design order)
                 pairs,                            # design.mat$pairs with integer i,j
                 sample.meta, sample.id,           # to classify AA vs AB via contrast var
                 contrast,                         # triplet for ref/alt
                 center.by = c("shift","total"),  # how to center baselines
                 block.vars = NULL) {
  center.by <- match.arg(center.by)
  
  # classify rows
  tr <- parseTriplet(contrast)
  g  <- factor(sample.meta[[tr$var]])
  gi <- as.character(g[pairs$i])
  gj <- as.character(g[pairs$j])
  is.AA <- (gi == tr$ref & gj == tr$ref)
  is.AB <- ((gi == tr$alt & gj == tr$ref) | (gi == tr$ref & gj == tr$alt))
  is.BB <- (gi == tr$alt & gj == tr$alt)
  if (!any(is.AB)) stop("No AB rows found.")
  ct <- colnames(yhat.all)
  
  # --- centering baselines ---
  if (center.by == "total") {
    mu.AA <- colMeans(yhat.all[is.AA, , drop = FALSE], na.rm = TRUE)
    base   <- matrix(rep(mu.AA, each = nrow(yhat.all)), nrow(yhat.all), byrow = FALSE,
             dimnames = list(rownames(yhat.all), ct))

  } else if (center.by == "shift") {
    mu.AA <- colMeans(yhat.all[is.AA, , drop = FALSE], na.rm = TRUE)
    mu.BB <- colMeans(yhat.all[is.BB, , drop = FALSE], na.rm = TRUE)
    base   <- matrix(rep((mu.AA + mu.BB) / 2, each = nrow(yhat.all)), nrow(yhat.all), byrow = FALSE,
             dimnames = list(rownames(yhat.all), ct))
  }
  # Δ for AB rows only
  delta.df <- yhat.all[is.AB, , drop = FALSE] - base[is.AB, , drop = FALSE]

  # --- block-wise summary (optional) ---
  delta.block <- NULL
  if(!is.null(block.vars)) {
  blk.cols <- intersect(block.vars, colnames(pair.meta))
  if (!length(blk.cols)) return(NULL)
  block.all <- do.call(paste, c(pair.meta[, blk.cols, drop = FALSE], sep = " | "))
  ct <- colnames(yhat.all)

  # per block: centering by AA (or (AA+BB)/2), then mean shift of AB
  delta.block <- do.call(rbind, lapply(unique(block.all), function(b) {
    idx.ab <- which(is.AB & block.all == b)
    idx.aa <- which(is.AA & block.all == b)
    idx.bb <- which(is.BB & block.all == b)
    if (!length(idx.ab)) return(NULL)

    mean.ab <- colMeans(yhat.all[idx.ab, , drop = FALSE], na.rm = TRUE)
    mu.aa   <- if (length(idx.aa)) colMeans(yhat.all[idx.aa, , drop = FALSE], na.rm = TRUE)
               else rep(NA_real_, length(ct))
    mu.bb  <- if (length(idx.bb)) colMeans(yhat.all[idx.bb, , drop = FALSE], na.rm = TRUE)
               else rep(NA_real_, length(ct))
    if (center.by == "total") {
      shift <- mean.ab - mu.aa
    } else if (center.by == "shift") {
      shift <- mean.ab - (mu.aa + mu.bb) / 2
    }
    data.frame(block = b, celltype = ct, shifts = as.numeric(shift),
               stringsAsFactors = FALSE)
  }))
  }
    return(list(shifts = delta.df, delta.block = delta.block))
}

extractFitsBlock <- function(res, p.dist, design.model) {
  F <- as.matrix(design.model$model$F)
  X <- as.matrix(design.model$model$X)

  # Coefs are on F's columns
  B_F <- matrix(res$coef, nrow = ncol(F), dimnames = list(colnames(F), colnames(p.dist$Y)))

  # Full fitted values: yhat = F %*% B_F
  y.fitted.all <- as.matrix(F %*% B_F)

  # Partial *without noise*: X %*% beta_X
  xnames <- intersect(colnames(X), rownames(B_F))
  partial.fit.all <- as.matrix(X[, xnames, drop = FALSE] %*% B_F[xnames, , drop = FALSE])

  # Partial *with noise*: (X %*% beta_X) + residuals from F-fit
  partial.all <- partial.fit.all + as.matrix(res$residuals)

  list(partial.fit = partial.fit.all,  # X * beta_X
       partial     = partial.all,      # X * beta_X + e
       fitted      = y.fitted.all)     # F * beta
}


extractFitsFL <- function(res, p.dist, design.model) {
    X    <- as.matrix(design.model$model$X)
    core <- design.model$model$core.rows  # logical length n, or NULL
    rn   <- rownames(X); cn <- colnames(p.dist$Y)
    y.r     <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(p.dist$Y),
                              dimnames = list(rn, cn))
    if (!is.null(res$y.resid)) {
        idx <- if (is.null(core)) seq_len(nrow(X)) else which(core)
        stopifnot(nrow(res$y.resid) == length(idx),
                  ncol(res$y.resid) == ncol(p.dist$Y))
        y.r[idx, ] <- as.matrix(res$y.resid)           # X_core %*% beta_X
    } else { stop("Residualized values not provided.")}
    return(y.r)   # X * beta_X (on core rows if partial_core provided; NA elsewhere)
}


## Helpers for visualization

# Align AB rows in shifts to design.pairs using sample.meta
# Returns integer vector of row indices into design.pairs
# Arguments:
## - shifts: AB-only matrix with rownames "SampleA__SampleB"
## - x.Pair$pairs: data.frame with columns i,j in the design/lower-tri order
## - sample.meta: data.frame in the same sample order used to build x.Pair
## - sample.id: the column in sample.meta that produced those sample labels
#' @keywords internal
alignABrows <- function(shifts, pairs, sample.meta, sample.id = NULL) {
    if (!is.null(sample.id) && sample.id %in% names(sample.meta)) {
        s <- as.character(sample.meta[[sample.id]])
    } else if (!is.null(rownames(sample.meta)) && all(nzchar(rownames(sample.meta)))) {
        s <- rownames(sample.meta)
    } else {
        stop("Need sample.id (column in sample.meta) OR rownames(sample.meta) to make pair labels.")
    }
    s <- trimws(s)
    pair_labels_all <- paste(s[pairs$i], s[pairs$j], sep="__")
    ab_rows <- rownames(shifts)
    m1 <- match(ab_rows, pair_labels_all)
    if (anyNA(m1)) {
        pair_labels_swapped <- paste(s[pairs$j], s[pairs$i], sep="__")
        m2 <- match(ab_rows, pair_labels_swapped)
        use_swapped <- sum(!is.na(m2)) > sum(!is.na(m1))
        if (use_swapped) {
            ab_pos <- m2
            which_unmatched <- which(is.na(m2))
            if (length(which_unmatched)) {
                warning("Some AB rows did not match even after swapping: ",
                        paste(head(ab_rows[which_unmatched], 5), collapse=", "), " ...")
            }
        } else {
            ab_pos <- m1
            which_unmatched <- which(is.na(m1))
            if (length(which_unmatched)) {
                warning("Some AB rows did not match: ",
                        paste(head(ab_rows[which_unmatched], 5), collapse=", "), " ...")
            }
        }
    } else {
        ab_pos <- m1
    }
    if (anyNA(ab_pos)) stop("Could not align all AB rows. Check that sample.id matches AB labels.")
    
    ab_pos
}
  
# readable label from pair.meta keys, fallback to design hash
#' @keywords internal
designHash <- function(X) apply(round(X, 12), 1, paste, collapse="|")