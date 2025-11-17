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
estimateExpressionChange <- function(cm.per.type, cell.groups, pair.model, sample.per.cell,
                                     sample.ids = NULL, dist = "cor", dist.type = c("shift", "total", "var"), 
                                     perm.method = c("freedman-lane", "block"), robust.method = c("none", "huber", "winsor"),
                                     na.mode = c("drop", "impute_weak"), alternative = c("two-sided", "greater", "less"),
                                     n.permutations = 1000, p.adjust.method = "BH", trim = 0.2, return.residuals = FALSE, 
                                     return.sampled.stats = TRUE, return.sampled.fits = FALSE, 
                                     n.cores = 1, verbose = TRUE, ...) {
  dist.type <- match.arg(dist.type)
  perm.method <- match.arg(perm.method)
  robust.method <- match.arg(robust.method)
  na.mode <- match.arg(na.mode)

  cell.groups <- droplevels(factor(cell.groups))
  sample.type.table <- table(cell.groups, sample.per.cell[names(cell.groups)])
  #n.cores.inner <- max(floor(n.cores / max(1L, length(levels(cell.groups)))), 1)

  # Pairwise Distances
  p.dist <- estimateExpressionShiftsForCellType(cm.per.type, dist = dist, pair.model = pair.model, sample.ids = sample.ids)

  # Fitting and randomization
  res <- performLMPermutations(y=p.dist$Y, x=pair.model, n.permutations=n.permutations, perm.method=perm.method, 
                               robust.method = robust.method, na.mode = na.mode, alternative = alternative,
                               return.residuals = return.residuals,return.sampled.stats = return.sampled.stats,
                               return.sampled.fits = return.sampled.fits, n.cores = n.cores, ...)
  # R2 estimation
  r2 <- estimateR2PerTerm(pair.model$F, p.dist$Y, groups = makeGroupsPair(pair.model$F))

  #model.diagnostics <- lapply(res.per.type, `[[`, "model.diag") // TODO

  if (verbose) message("Done!\n")
  ct <- levels(cell.groups)
  names(res$stat.obs) <- names(res$pval) <- names(res$z.score) <- ct
  summary  <- data.frame(celltype=ct, obs.stat=res$stat.obs, pvalue=res$pval, zscore=res$z.score, stringsAsFactors=FALSE)
  summary$padjust <- p.adjust(summary$pvalue, method=p.adjust.method)
  summary <- summary[order(summary$obs.stat, decreasing=TRUE), ]

  out <- list(results= summary, res = res, perm.method = perm.method, robust.method= robust.method,
              dist.type = dist.type, p.dist = p.dist, r2 = r2, sample.table = sample.type.table, design.mat = pair.model
              )
  out
}

#' Estimate expression-based pairwise distances for all cell types
#' @param cm.norm List of normalized count matrices per cell type
#' @param dist what distance measure to use: 'cor' - Pearson's correlation, 'l2' - Euclidean distance, 'l1' - Manhattan distance
#' @param pair.model result of buildPairDesignMatrices()
#' @param sample.ids global sample IDs in the same order as rows of sample.metadata
#' @keywords internal
estimateExpressionShiftsForCellType <- function( cm.norm, dist, pair.model, sample.ids) {
    # ---- basic checks ----
    if (!is.list(cm.norm))
        stop("cm.norm must be a list of count matrices per cell type.")
    if (missing(pair.model) || is.null(pair.model$pairs))
        stop("pair.model (from buildPairDesignMatrices) with $pairs is required.")
    if (missing(sample.ids))
        stop("sample.ids must be provided and must align with rows of sample.metadata.")
    
    idx <- pair.model$pairs
    if (!(is.matrix(idx) && ncol(idx) == 2L))
        stop("pair.model$pairs must be an n_pairs x 2 matrix of sample indices (i,j).")
    
    if (length(sample.ids) < max(idx))
        stop("sample.ids length is smaller than max(pair.model$pairs).")
    
    # pair names for rows of Y / dists (for printing)
    pair.names <- paste(sample.ids[idx[, "i"]], sample.ids[idx[, "j"]], sep = "__")
    
    # cell type names
    ctnames <- names(cm.norm)
    if (is.null(ctnames)) ctnames <- paste0("CT", seq_along(cm.norm))
    
    # ---- allocate Y: rows = pairs, cols = cell types ----
    Y <- matrix(NA_real_, nrow = nrow(idx), ncol = length(cm.norm),
                dimnames = list(pair.names, ctnames))
    
    # ---- loop over cell types ----
    for (t in seq_along(cm.norm)) {
        X <- cm.norm[[t]]    # samples x genes for this type (subset of samples)
        Xp <- X              # removed any PCA / gene selection
        
        # ensure rownames are sample IDs; needed for mapping
        if (is.null(rownames(Xp))) {
            stop("cm.norm[[", t, "]] has no rownames; cannot map to sample.ids.")
        }
        
        # (1) Distances — square matrix D over *local* rows of this type
        D <- NULL
        if (dist == "cor") {
            Xc <- as.matrix(Xp)
            n_rows <- nrow(Xc)
            if (is.null(n_rows) || n_rows < 2L) {
                D <- matrix(0, n_rows, n_rows,
                            dimnames = list(rownames(Xc), rownames(Xc)))
            } else {
                R <- stats::cor(t(Xc), use = "pairwise.complete.obs", method = "pearson")
                D <- 1 - R
                if (is.null(rownames(R))) {
                    rownames(D) <- colnames(D) <- rownames(Xc)
                }
            }
        } else if (dist %in% c("l2", "l1")) {
            Xd <- as.matrix(Xp)
            n_rows <- nrow(Xd)
            if (n_rows < 2L) {
                D <- matrix(0, n_rows, n_rows,
                            dimnames = list(rownames(Xd), rownames(Xd)))
            } else if (anyNA(Xd)) {
                D <- matrix(NA_real_, n_rows, n_rows,
                            dimnames = list(rownames(Xd), rownames(Xd)))
                diag(D) <- 0
                for (i in seq_len(n_rows - 1L)) {
                    xi <- Xd[i, ]
                    for (j in (i + 1L):n_rows) {
                        xj <- Xd[j, ]
                        if (anyNA(xi) || anyNA(xj)) {
                            d <- NA_real_
                        } else if (dist == "l2") {
                            d <- sqrt(sum((xi - xj)^2))
                        } else {
                            d <- sum(abs(xi - xj))
                        }
                        D[i, j] <- D[j, i] <- d
                    }
                }
            } else {
                method <- if (dist == "l2") "euclidean" else "manhattan"
                D <- as.matrix(stats::dist(Xd, method = method))
                if (is.null(rownames(D)) && !is.null(rownames(Xd))) {
                    rownames(D) <- colnames(D) <- rownames(Xd)
                }
            }
        } else {
            stop("Unknown distance: ", dist)
        }
        
        # ensure D is square and has rownames
        if (!is.matrix(D)) D <- as.matrix(D)
        dimD <- dim(D)
        if (length(dimD) != 2L || dimD[1] != dimD[2]) {
            stop("Distance is not a square matrix for cell type '", ctnames[t],
                 "': got ", paste(dimD, collapse = "x"))
        }
        if (is.null(rownames(D))) {
            stop("Distance matrix D for cell type '", ctnames[t],
                 "' has no rownames; cannot map to sample.ids.")
        }
        
        # (2) Vectorize using global pairs, mapping to local D indices via sample IDs
        Y[, t] <- vectorizeLowerTri(D, pairs = idx, sample.ids = sample.ids, row.ids = rownames(D))
    }
    
    list(Y = Y)        # all pairs x cell types (NA where sample missing in type)
    
}

#' @keywords internal
vectorizeLowerTri <- function(D, pairs, sample.ids, row.ids = rownames(D)) {
    # D: square distance matrix over *local* rows
    # pairs: n_pairs x 2 integer matrix with columns i,j (global sample indices)
    # sample.ids: length n_global, IDs in the same order as rows of sample.meta
    # row.ids: IDs corresponding to rows of D (subset of sample.ids)
    
    if (is.null(row.ids)) {
        stop("Distance matrix D must have rownames (sample IDs) to map pairs correctly.")
    }
    
    # Map global sample indices -> local row indices in D
    # global k -> sample.ids[k] -> match in row.ids
    global2local <- match(sample.ids, row.ids)  # length n_global, can be NA
    
    vapply(seq_len(nrow(pairs)), function(k) {
        gi <- pairs[k, 1]
        gj <- pairs[k, 2]
        li <- global2local[gi]
        lj <- global2local[gj]
        
        if (is.na(li) || is.na(lj)) {
            NA_real_              # this pair involves a sample missing in this cell type
        } else {
            D[li, lj]
        }
    }, numeric(1))
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

#' Filter cell types for paired analysis using paired design coverage
#'
#' This replaces the old `filterCellTypesByNSamples()` logic in a *paired*
#' differential analysis context.
#'
#' It keeps only those cell types that:
#'   1. Have enough cells per sample (`min.cells.per.sample`) AND
#'   2. Are observed in at least `min.samp.per.type` *distinct samples* that
#'      actually participate in the paired contrast of interest.
#'
#' How "participate in the paired contrast" is determined:
#'   - We look at `pairDesign$model$core.rows`, which is a logical vector
#'     (length = number of pair rows in the paired design) telling us which
#'     *pairs of samples (i,j)* contribute to the core contrast `X`.
#'   - We then mark any sample that appears in at least one such "core" pair
#'     as a "usable" sample.
#'
#' Steps performed internally:
#'   1. Extract the set of usable samples (based on core pairs).
#'   2. Build a table of how many single cells of each `cell.type` each `sample.id`
#'      contributed (`Freq`).
#'   3. Keep only those (Type,Sample) combos where:
#'        - `Freq >= min.cells.per.sample`, and
#'        - `Sample` is among usable samples.
#'   4. For each Type, count how many distinct usable samples survive. If this
#'      count is at least `min.samp.per.type`, we keep that Type.
#'
#' The output can be fed into downstream preprocessing to drop unsupported
#' cell types before building per-type matrices.
#'
#' @param cell.groups factor/character vector of length = total cells, giving
#'        the cell type (cluster / annotation) for each cell.
#' @param sample.per.cell character or factor vector of same length as
#'        `cell.groups`, giving which biological sample each cell came from.
#' @param pairDesign result of `buildPairDesignMatrices()`. We assume:
#'        - `pairDesign$pairs`: integer matrix (n_pairs x 2), each row is the
#'          indices `(i,j)` of the two samples in that pair, 1-based and
#'          referring to rows of the original `sample.meta`.
#'        - `pairDesign$model$core.rows`: logical vector of length n_pairs;
#'          TRUE where that pair contributes to the tested contrast (i.e. has
#'          non-negligible X signal after contrast splitting).
#' @param sample.ids character vector of length `nrow(sample.meta)` in the same
#'        row order that was given to `buildPairDesignMatrices()`. So if
#'        `sample.meta` had rownames, you usually pass `rownames(sample.meta)`.
#'        These IDs must match the values used in `sample.per.cell`.
#' @param min.cells.per.sample integer. Minimum number of cells of a given type
#'        that must appear in a given sample for that (Type,Sample) combo to count.
#' @param min.samp.per.type integer. Minimum number of *distinct usable samples*
#'        in which a type must appear (above the per-sample threshold) in order
#'        for that type to be kept.
#' @param use.core.rows logical (default TRUE). If TRUE, only pairs with
#'        `core.rows==TRUE` are used to establish which samples are "usable".
#'        If FALSE, *all* pairs in `pairDesign$pairs` are considered usable.
#' @param verbose logical. If TRUE, emit messages about dropped cell types.
#'
#' @return A list with:
#'   - `freq.table`: data.frame with columns
#'        `Type`, `Sample`, `Freq`, `Usable`
#'     after applying the per-sample coverage threshold.
#'   - `kept.types`: character vector of cell types that survive filtering.
#'   - `kept.by.sample`: named list. For each surviving sample ID, a unique
#'        character vector of the kept types that sample supports.
#'
#' Typical usage:
#' \preformatted{
#' out <- filterCellTypesByCoveragePairs(
#'   cell.groups         = cell.groups,
#'   sample.per.cell     = sample.per.cell,
#'   pairDesign          = pairDesign,
#'   sample.ids          = rownames(sample.meta),
#'   min.cells.per.sample= 10,
#'   min.samp.per.type   = 2,
#'   use.core.rows       = TRUE,
#'   verbose             = TRUE
#' )
#'
#' ## Keep only out$kept.types when building per-type expression matrices.
#' }
#'
#' Rationale:
#'   - We no longer assume a 2-condition factor like "Condition".
#'   - We only keep types with enough replication across biologically distinct
#'     samples that actually drive the paired contrast (i.e. show up in core pairs).
#'
#' @export
filterCellTypesByCoveragePairs <- function(cell.groups, sample.per.cell, pairDesign, sample.ids,
                                           min.cells.per.sample, min.samp.per.type, use.core.rows = TRUE,
                                           verbose = TRUE) {
  # ---- basic checks ----
  if (length(cell.groups) != length(sample.per.cell)) {
    stop("cell.groups and sample.per.cell must have the same length (one entry per cell).")
  }
  if (!is.list(pairDesign) ||
      is.null(pairDesign$pairs) ||
      is.null(pairDesign$core.rows)) {
    stop("pairDesign must be the full list returned by buildPairDesignMatrices() ",
         "and must include $pairs and $core.rows.")
  }
  pairs.mat <- pairDesign$pairs
  if (!(is.matrix(pairs.mat) && ncol(pairs.mat) == 2L)) {
    stop("pairDesign$pairs must be an n_pairs x 2 integer matrix of sample indices.")
  }
  n.pairs <- nrow(pairs.mat)
  core.rows <- pairDesign$core.rows
  if (!is.null(core.rows) && length(core.rows) != n.pairs) {
    stop("pairDesign$core.rows length does not match nrow(pairDesign$pairs).")
  }
  if (length(sample.ids) < max(pairs.mat)) {
    stop("sample.ids does not cover all indices mentioned in pairDesign$pairs.")
  }
  
  # ---- 1. determine usable samples from the paired design ----
  if (!n.pairs) {
    if (verbose) {
      message("No pairs in pairDesign$pairs; returning empty filter result.")
    }
    return(list(
      freq.table     = data.frame(Type=character(0), Sample=character(0),
                                  Freq=integer(0), Usable=logical(0),
                                  stringsAsFactors = FALSE),
      kept.types     = character(0),
      kept.by.sample = list()
    ))
  }
  
  if (use.core.rows && !is.null(core.rows)) {
    keep.pair.rows <- which(core.rows)
  } else {
    keep.pair.rows <- seq_len(n.pairs)
  }
  
  if (!length(keep.pair.rows)) {
    if (verbose) {
      message("No informative pairs (core.rows is empty TRUE set); all types dropped.")
    }
    return(list(
      freq.table     = data.frame(Type=character(0), Sample=character(0),
                                  Freq=integer(0), Usable=logical(0),
                                  stringsAsFactors = FALSE),
      kept.types     = character(0),
      kept.by.sample = list()
    ))
  }
  
  sample.idx.used <- unique(as.integer(pairs.mat[keep.pair.rows, , drop = FALSE]))
  usable.samples  <- unique(sample.ids[sample.idx.used])
  
  # ---- 2. build per-(Type,Sample) cell counts ----
  # coerce to plain factors/chars
  cell.types.vec     <- droplevels(factor(cell.groups))
  sample.per.cell.vec<- as.character(sample.per.cell)
  
  # table(Type, Sample)
  freq.table <- table(Type = cell.types.vec,
                      Sample = sample.per.cell.vec)
  freq.table <- as.data.frame(freq.table, stringsAsFactors = FALSE)
  
  # ---- 3. enforce per-sample minimum coverage and "sample is usable" ----
  # keep only Type/Sample combos with enough cells
  freq.table <- freq.table[freq.table$Freq >= min.cells.per.sample, , drop = FALSE]
  if (!nrow(freq.table)) {
    if (verbose) {
      message("No (Type,Sample) pairs pass min.cells.per.sample in paired context.")
    }
    return(list(
      freq.table     = transform(freq.table, Usable = logical(0)),
      kept.types     = character(0),
      kept.by.sample = list()
    ))
  }
  
  # mark which Sample entries are actually usable for the contrast
  freq.table$Usable <- freq.table$Sample %in% usable.samples
  
  # ---- 4. count usable samples per Type ----
  type.sample.split <- split(freq.table$Sample[freq.table$Usable],
                             freq.table$Type[freq.table$Usable])
  usable.counts <- vapply(type.sample.split,
                          function(sv) length(unique(sv)),
                          integer(1))
  
  kept.types <- names(usable.counts)[usable.counts >= min.samp.per.type]
  kept.types <- kept.types[!is.na(kept.types)]
  
  # produce kept.by.sample for convenience
  keep.mask <- freq.table$Usable & freq.table$Type %in% kept.types
  kept.by.sample <- split(freq.table$Type[keep.mask],
                          freq.table$Sample[keep.mask])
  kept.by.sample <- lapply(kept.by.sample, function(v) unique(as.character(v)))
  
  # ---- 5. messaging ----
  if (verbose) {
    removed.types <- setdiff(levels(cell.types.vec), kept.types)
    if (length(removed.types)) {
      message(
        "Excluding cell types in paired context: ",
        paste(removed.types, collapse = ", "),
        " (not enough distinct usable samples with >= ",
        min.cells.per.sample, " cells)"
      )
    }
  }
  
  # ---- 6. return ----
  list(
    freq.table     = freq.table,
    kept.types     = kept.types,
    kept.by.sample = kept.by.sample
  )
}



#' @keywords internal
filterExpressionDistanceInput <- function(
  cms, cell.groups, sample.per.cell, pair.model, sample.ids, keep.all=FALSE,
  min.cells.per.sample=10, min.samp.per.type=2, min.gene.frac=0.01,
  genes=NULL, verbose=FALSE
) {
  stopifnot(is.list(cms), length(cms) > 0)
  all.samples <- names(cms)

  if (!keep.all) {
    cell.names <- lapply(cms, rownames) %>% unlist()
    freq.table <- filterCellTypesByCoveragePairs(cell.groups[cell.names], sample.per.cell[cell.names], pair.model, sample.ids, min.cells.per.sample, min.samp.per.type, verbose=verbose)
    filt.types.per.samp <- freq.table$kept.by.sample #freq.table %$% split(Type, Sample)
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

  # ---- Build per-cell-type matrices (kept NORM only) ----
  cell.groups <- droplevels(cell.groups[cell.names])
  all.types <- levels(cell.groups)

  cms.filt <- cms.filt[all.samples]

  cm.per.type <- sccore::sn(all.types) %>% lapply(function(ct) {
    # norm: row present -> take it; missing -> NA row (as before)
    rows.norm <- lapply(cms.filt, function(x) {
      # skip NULL/non-matrix entries
      if (is.null(x) || is.null(ncol(x))) {
        return(NULL)
      }

      if (ct %in% rownames(x)) {
        x[ct, , drop=FALSE]
      } else {
        m <- matrix(NA_real_, nrow=1, ncol=ncol(x), dimnames=list(ct, colnames(x)))
        as(m, "dgCMatrix")
      }
    })
    rows.norm <- Filter(Negate(is.null), rows.norm)
    mat <- do.call(rbind, rows.norm)                   # samples x genes
    denom <- pmax(1, rowSums(mat, na.rm=TRUE))
    mat <- mat / denom
    mat <- log10(mat * 1e3 + 1)
    # set rownames only for the non-NULL entries
    valid <- !vapply(cms.filt, function(x) is.null(x) || is.null(ncol(x)), logical(1))
    rownames(mat) <- names(cms.filt)[valid]
    mat
  })
  names(cm.per.type) <- all.types
  return(list("cm.per.type" = cm.per.type))       # normalized (use everywhere downstream)
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


## Generic grouping for pair-level design matrices, used in R2 estimation per variable across cell types (updated for new pair design)
## - collapses all "<var>_pair*" columns into one group per <var>
##   (and prefers any "<var>_pair_contrast" column if present)
## - groups numeric pair columns "pair_<var>_(mean|diff)" into one group "pair_<var>"
## - keeps "(Intercept)" if present
## - anything else stays as-is (unless include.other=TRUE, then grouped under "Other")
makeGroupsPair <- function(M, include.other = FALSE) {
    stopifnot(!is.null(colnames(M)))
    cn <- colnames(M)
    groups <- list()
    
    ## 1) Intercept
    if ("(Intercept)" %in% cn) {
        groups$Intercept <- "(Intercept)"
    }
    
    ## 2) Factor pair terms: "<var>_pair*"
    pair.terms <- grep("_pair", cn, value = TRUE)
    if (length(pair.terms)) {
        # extract variable name before "_pair"
        varnames <- sub("_pair.*$", "", pair.terms)
        for (v in unique(varnames)) {
            cols <- pair.terms[varnames == v]
            # if there's a collapsed contrast, prefer that
            contrast.cols <- grep("_pair_contrast$", cols, value = TRUE)
            if (length(contrast.cols)) {
                groups[[v]] <- contrast.cols
            } else {
                groups[[v]] <- cols
            }
        }
    }
    
    ## 3) Numeric pair terms: "pair_<var>_(mean|diff)"
    num.pair <- grep("^pair_[A-Za-z0-9.]+_(mean|diff)$", cn, value = TRUE)
    if (length(num.pair)) {
        base <- sub("^pair_([A-Za-z0-9.]+)_(mean|diff)$", "\\1", num.pair)
        for (v in unique(base)) {
            cols <- num.pair[base == v]
            gname <- paste0("pair_", v)
            if (gname %in% names(groups)) {
                groups[[gname]] <- unique(c(groups[[gname]], cols))
            } else {
                groups[[gname]] <- cols
            }
        }
    }
    
    ## 4) Leftovers
    used <- unlist(groups, use.names = FALSE)
    leftovers <- setdiff(cn, used)
    if (length(leftovers)) {
        if (include.other) {
            groups$Other <- leftovers
        } else {
            for (lf in leftovers) groups[[lf]] <- lf
        }
    }
    
    groups
}




#' Extract Partial or Fitted Expression Shifts for Visualization. "partial" corresponds to the distance explained
#'   by the contrast of interest only (i.e., after regressing out nuisance covariates).
#' @param p.dist Output of \code{estimateExpressionShiftsForCellType} function for pairwise distances.
#' @param res Output of \code{performLMPermutations} function.
#' @param design.mat Output of \code{buildDesignMatrix} function. (or pairwise design matrices)
#' @param block.vars Optional character vector of column names in \code{design.mat$pair.meta}
#'   to use for block-wise summary of shifts (e.g., batch or other grouping variable).
#'   If provided, the function will compute mean shifts within each block and
#'   return an additional data frame with this summary.
#' @return A list with two data frames:
#'   \item{pair.covariates}{Data frame suitable for plotting shifts with covariate patterns.}
#'   \item{df.blocks}{(optional) Data frame with block-wise summary of shifts.}
#'   \item{df.perm}{(optional) Data frame with permutation statistics for background distribution.}
#' @details This function extracts expression shifts from the results of
#'   \code{performLMPermutations} and prepares data frames for visualization.
#'   It computes partial shifts (explained by the contrast of interest only), centers them
#'   by the mean of ref-ref (AA) pairs and/or alt-alt (BB) pairs, and codes covariate patterns
#'   for plotting. If \code{block.vars} is provided, it also computes mean shifts
#'   within each block defined by the specified variables.
#' @keywords internal
extractPairwiseShifts <- function(res, p.dist, design.mat,
                                  perm.method = c("freedman-lane", "block"),
                                  block.vars = NULL) {
    perm.method <- match.arg(perm.method)
    
    if (perm.method == "block") {
        blk <- extractFitsBlock(res, p.dist, design.model = design.mat)
        partial.fit  <- blk$partial
        contrast.vec <- design.mat$contrast.F
        F <- as.matrix(design.mat$F)
        
        coef.mat <- NULL
        if (!is.null(res$coef)) {
            expect <- ncol(F) * ncol(p.dist$Y)
            if (length(res$coef) == expect) {
                coef.mat <- matrix(res$coef,
                                   nrow = ncol(F),
                                   dimnames = list(colnames(F), colnames(p.dist$Y)))
            }
        }
        
    } else if (perm.method == "freedman-lane") {
        partial.fit  <- extractFitsFL(res, p.dist, design.model = design.mat)
        contrast.vec <- design.mat$contrast.X
        X <- as.matrix(design.mat$X)
        
        coef.mat <- NULL
        if (!is.null(res$coef)) {
            expect <- ncol(X) * ncol(p.dist$Y)
            if (length(res$coef) == expect) {
                coef.mat <- matrix(res$coef,
                                   nrow = ncol(X),
                                   dimnames = list(colnames(X), colnames(p.dist$Y)))
            }
        }
        
    } else {
        stop("Unknown model type: ", perm.method)
    }
    
    ci <- getContrastInfo(contrast.vec)
    center.contrast <- ci$center
    x0.override <- NULL   # <- new: numeric baseline for continuous contrasts
    
    ## --- continuous contrasts: pick x0 using numeric_ref_used when possible ---
    if (ci$kind == "continuous") {
        var <- ci$var  # e.g. "pair_age_diff" or "pair_age_mean"
        
        # Heuristic for pair numerics: pair_<v>_(mean|diff)
        m  <- regexec("^pair_([A-Za-z0-9.]+)_(mean|diff)$", var)
        mt <- regmatches(var, m)[[1]]
        
        if (length(mt) == 3L) {
            base <- mt[2]     # "age"
            kind <- mt[3]     # "mean" or "diff"
            
            if (kind == "diff") {
                # natural: baseline is zero difference
                x0.override <- 0
            } else if (kind == "mean") {
                # natural: baseline is numeric_ref_used[[base]] if present
                if (!is.null(design.mat$numeric_ref_used) &&
                    base %in% names(design.mat$numeric_ref_used)) {
                    x0.override <- design.mat$numeric_ref_used[[base]]
                } else {
                    # leave NULL here; we'll fall back to mean(x) inside extractContrastDeltas()
                }
            }
        }
        
        # if we still don't have a numeric baseline, fall back to the old
        # factor-centering mechanism via center.contrast
        if (is.null(center.contrast)) {
            for (v in list(design.mat$contrast.F, design.mat$contrast.X)) {
                if (is.null(v)) next
                v.center <- attr(v, "center")
                if (!is.null(v.center)) { center.contrast <- v.center; break }
                if (any(grepl("=", names(v), fixed = TRUE))) {
                    center.contrast <- v
                    break
                }
            }
            # note: if center.contrast is still NULL and x0.override is NULL,
            # extractContrastDeltas() will error with a clear message.
        }
    }
    
    delta <- extractContrastDeltas(
        yhat.all        = partial.fit,      # contrast-specific fits
        pair.meta       = design.mat$pair.meta,
        contrast.vec    = contrast.vec,
        block.vars      = block.vars,
        coef.mat        = coef.mat,
        center.contrast = center.contrast,
        x0.override     = x0.override       # <- NEW
    )
    
    df.shifts <- as.data.frame(delta$shifts) |>
        tibble::rownames_to_column("pair") |>
        tidyr::pivot_longer(-pair, names_to = "celltype", values_to = "shifts") |>
        dplyr::filter(is.finite(shifts))
    
    stats.perm <- if (!is.null(res$stats.perm)) res$stats.perm else NULL
    stat.obs   <- if (!is.null(res$stat.obs))   res$stat.obs   else NULL
    df.perm     <- NULL
    df.perm.obs <- NULL
    
    if (!is.null(stats.perm)) {
        ct <- colnames(stats.perm)
        df.perm <- as.data.frame(stats.perm) |>
            tidyr::pivot_longer(tidyselect::all_of(ct),
                                names_to = "celltype",
                                values_to = "perm_stat") |>
            dplyr::filter(is.finite(perm_stat))
        if (!is.null(stat.obs)) {
            df.perm.obs <- data.frame(celltype = ct,
                                      obs      = as.numeric(stat.obs),
                                      stringsAsFactors = FALSE)
        }
    }
    
    out <- list(
        df.shifts       = df.shifts,
        pair.covariates = design.mat$pair.meta[rownames(delta$shifts), , drop = FALSE],
        dist.type       = delta$inferred,
        contrast        = contrast.vec
    )
    if (!is.null(delta$delta.block)) out$df.blocks   <- delta$delta.block
    if (!is.null(df.perm))           out$df.perm     <- df.perm
    if (!is.null(df.perm.obs))       out$df.perm.obs <- df.perm.obs
    if (!is.null(delta$x.center))    out$x.center    <- delta$x.center
    
    out
}


#' Extract Contrast-Specific Expression Shifts
#' @keywords internal
extractContrastDeltas <- function(yhat.all, pair.meta, contrast.vec,
                                  block.vars = NULL, tol = 1e-8,
                                  coef.mat = NULL, center.contrast = NULL,
                                  x0.override = NULL) {
    ci <- getContrastInfo(contrast.vec)
    
    # Align pair.meta to yhat rows
    if (is.null(rownames(yhat.all)))
        stop("yhat.all must have rownames to align with pair.meta.")
    if (is.null(rownames(pair.meta)))
        stop("pair.meta must have rownames to align with yhat.all.")
    pos.pm <- match(rownames(yhat.all), rownames(pair.meta))
    if (any(is.na(pos.pm)))
        stop("Could not align pair.meta to yhat rows; check rownames.")
    pm.aligned <- pair.meta[pos.pm, , drop = FALSE]
    
    ## ----- FACTOR CONTRAST PATH -----
    if (ci$kind == "factor") {
        w        <- ci$w
        pair.var <- ci$var
        inferred <- inferCenteringType(w, tol)
        
        lab.aligned <- as.character(pm.aligned[[pair.var]])
        if (is.null(lab.aligned))
            stop("pair.meta[['", pair.var, "']] not found.")
        
        is.pos <- lab.aligned %in% names(w)[w > tol]
        if (!any(is.pos))
            stop("No rows match positive-weight levels in pair.meta[['", pair.var, "']].")
        
        neg.levels <- names(w)[w < -tol]
        denom <- sum(w[w > tol]); if (denom <= 0) stop("Sum of positive weights must be > 0.")
        
        mu.neg <- lapply(neg.levels, function(L) {
            idx <- (lab.aligned == L)
            if (!any(idx)) return(rep(NA_real_, ncol(yhat.all)))
            colMeans(yhat.all[idx, , drop = FALSE], na.rm = TRUE)
        })
        names(mu.neg) <- neg.levels
        
        neg.sum <- Reduce(`+`, Map(function(L) { (-w[L]) * mu.neg[[L]] }, neg.levels),
                          init = rep(0, ncol(yhat.all)))
        baseline <- neg.sum / denom
        
        base.mat <- matrix(rep(baseline, each = nrow(yhat.all)),
                           nrow = nrow(yhat.all), byrow = FALSE,
                           dimnames = list(rownames(yhat.all), colnames(yhat.all)))
        
        delta.df <- yhat.all[is.pos, , drop = FALSE] - base.mat[is.pos, , drop = FALSE]
        
        # block-wise summary
        delta.block <- NULL
        if (!is.null(block.vars)) {
            blk.cols <- intersect(block.vars, colnames(pm.aligned))
            if (length(blk.cols)) {
                block.aligned <- do.call(paste, c(pm.aligned[, blk.cols, drop = FALSE], sep = " | "))
                ublk <- unique(block.aligned)
                delta.block <- do.call(rbind, lapply(ublk, function(b) {
                    inb <- (block.aligned == b)
                    inb.pos <- which(inb & is.pos); if (!length(inb.pos)) return(NULL)
                    
                    mu.neg.b <- lapply(neg.levels, function(L) {
                        idx <- (lab.aligned == L) & inb
                        if (!any(idx)) return(rep(NA_real_, ncol(yhat.all)))
                        colMeans(yhat.all[idx, , drop = FALSE], na.rm = TRUE)
                    })
                    names(mu.neg.b) <- neg.levels
                    
                    neg.sum.b  <- Reduce(`+`, Map(function(L) { (-w[L]) * mu.neg.b[[L]] }, neg.levels),
                                         init = rep(0, ncol(yhat.all)))
                    baseline.b <- neg.sum.b / denom
                    mean.pos.b <- colMeans(yhat.all[inb.pos, , drop = FALSE], na.rm = TRUE)
                    shift.b    <- mean.pos.b - baseline.b
                    
                    data.frame(block   = b,
                               celltype = colnames(yhat.all),
                               shifts   = as.numeric(shift.b),
                               stringsAsFactors = FALSE)
                }))
            }
        }
        return(list(shifts = delta.df,
                    delta.block = delta.block,
                    inferred = inferred))
    }
    
    ## ----- CONTINUOUS CONTRAST PATH -----
    var <- ci$var
    if (!var %in% names(pm.aligned))
        stop("Continuous contrast variable '", var, "' not found in pair.meta.")
    x <- pm.aligned[[var]]
    if (!is.numeric(x))
        stop("pair.meta[['", var, "']] must be numeric for a continuous contrast.")
    
    if (is.null(coef.mat) || !(var %in% rownames(coef.mat)))
        stop("Provide coef.mat with a row named '", var, "' (rows = term names; cols = cell types).")
    
    # choose numeric baseline x0 (scalar) for continuous contrast
    x0.scalar <- NA_real_
    
    if (!is.null(x0.override)) {
        # If user supplies an anchor (from numeric_ref_used), use it directly
        x0.scalar <- as.numeric(x0.override)[1]
    } else {
        # Heuristics for pair_<var>_(mean|diff)
        m  <- regexec("^pair_([A-Za-z0-9.]+)_(mean|diff)$", var)
        mt <- regmatches(var, m)[[1]]
        if (length(mt) == 3L) {
            kind <- mt[3]
            if (kind == "diff") {
                # zero difference is natural baseline
                x0.scalar <- 0
            } else if (kind == "mean") {
                # fallback: mean of x when no explicit anchor
                x0.scalar <- mean(x, na.rm = TRUE)
            }
        }
    }
    
    # If still NA, use factor-based centering rule (as before)
    if (!is.finite(x0.scalar)) {
        if (!is.null(ci$center)) center.contrast <- ci$center
        if (is.null(center.contrast))
            stop("Provide centering via attr(contrast,'center'), 'center.contrast', ",
                 "or x0.override for continuous contrasts.")
        
        cf <- getContrastInfo(center.contrast)
        if (cf$kind != "factor")
            stop("'center.contrast' must be factor-like with names '<pair_var>=<level>'.")
        pair.var <- cf$var
        if (!(pair.var %in% names(pm.aligned)))
            stop("pair.meta lacks the factor variable '", pair.var, "' for centering.")
        
        x0.scalar <- contrastCenterX0(x,
                                      pair.level = as.character(pm.aligned[[pair.var]]),
                                      w = cf$w, tol = tol)
    }
    
    beta <- coef.mat[var, , drop = TRUE]
    x.c  <- x - x0.scalar
    delta.df <- x.c %*% t(beta)
    rownames(delta.df) <- rownames(yhat.all)
    colnames(delta.df) <- colnames(yhat.all)
    
    delta.block <- NULL
    if (!is.null(block.vars)) {
        blk.cols <- intersect(block.vars, colnames(pm.aligned))
        if (length(blk.cols)) {
            block.aligned <- do.call(paste, c(pm.aligned[, blk.cols, drop = FALSE], sep = " | "))
            ublk <- unique(block.aligned)
            delta.block <- do.call(rbind, lapply(ublk, function(b) {
                idx <- which(block.aligned == b); if (!length(idx)) return(NULL)
                shift.b <- colMeans(delta.df[idx, , drop = FALSE], na.rm = TRUE)
                data.frame(block   = b,
                           celltype = colnames(delta.df),
                           shifts   = as.numeric(shift.b),
                           stringsAsFactors = FALSE)
            }))
        }
    }
    
    list(
        shifts      = delta.df,
        delta.block = delta.block,
        inferred    = "continuous",
        x.center    = x0.scalar   # << scalar baseline used
    )
}



extractFitsBlock <- function(res, p.dist, design.model) {
  F <- as.matrix(design.model$F)
  X <- as.matrix(design.model$X)
  na.rows <- apply(is.na(p.dist$Y), 1, all)
  # Coefs are on F's columns
  B.F <- matrix(res$coef, nrow = ncol(F), dimnames = list(colnames(F), colnames(p.dist$Y)))

  # Full fitted values: yhat = F %*% B.F
  y.fitted.all <- as.matrix(F %*% B.F)

  # Partial *without noise*: X %*% beta_X
  xnames <- intersect(colnames(X), rownames(B.F))
  partial.fit.all <- as.matrix(X[, xnames, drop = FALSE] %*% B.F[xnames, , drop = FALSE])

  # Partial *with noise*: (X %*% beta_X) + residuals from F-fit
  partial.all <- partial.fit.all + as.matrix(res$residuals)
  
  ## Drop rows that had NA in the original response
  y.fitted.all    <- y.fitted.all[which(!na.rows), , drop = FALSE]
  partial.fit.all <- partial.fit.all[which(!na.rows), , drop = FALSE]
  partial.all     <- partial.all[which(!na.rows), , drop = FALSE]

  list(partial.fit = partial.fit.all,  # X * beta_X
       partial     = partial.all,      # X * beta_X + e
       fitted      = y.fitted.all)     # F * beta
}


extractFitsFL <- function(res, p.dist, design.model) {
    X    <- as.matrix(design.model$X)
    core <- design.model$core.rows  # logical length n, or NULL
    rn   <- rownames(X); cn <- colnames(p.dist$Y)
    na.rows <- apply(is.na(p.dist$Y), 1, all)
    y.r     <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(p.dist$Y),
                              dimnames = list(rn, cn))
    if (!is.null(res$y.resid)) {
        idx <- if (is.null(core)) seq_len(nrow(X)) else which(core)
        stopifnot(nrow(res$y.resid) == length(idx),
                  ncol(res$y.resid) == ncol(p.dist$Y))
        y.r[idx, ] <- as.matrix(res$y.resid)           # X_core %*% beta_X
    } else { stop("Residualized values not provided.")}
    y.r <- y.r[which(!na.rows), , drop = FALSE]
    return(y.r)   # X * beta_X (on core rows if partial_core provided; NA elsewhere)
}


## Helpers for visualization

# infer centering & relevant levels from weights
# Return one of: "shift", "total", "var", "custom"
#' @keywords internal
inferCenteringType <- function(w, tol = 1e-8) {
  stopifnot(!is.null(names(w)))
  pos <- names(w)[w >  tol]
  neg <- names(w)[w < -tol]
  nz  <- names(w)[abs(w) > tol]

  is.same <- function(lab) grepl("^(.+)\\1$", lab)   # "AA" pattern (e.g., Group1Group1)
  all.same <- function(labs) length(labs) > 0 && all(vapply(labs, is.same, logical(1)))
  any.cross <- function(labs) any(!vapply(labs, is.same, logical(1)))

  # VAR: contrast between two same levels (AA vs BB), no cross terms involved
  if (length(nz) == 2 && sign(w[nz[1]]) != sign(w[nz[2]]) && all.same(nz)) {
    return("var")
  }

  # TOTAL: exactly one positive (must be cross) and one negative (must be same)
  if (length(pos) == 1 && length(neg) == 1 &&
      !is.same(pos) && is.same(neg)) {
    return("total")
  }

  # SHIFT: one positive (cross) and two negatives (both same)
  if (length(pos) == 1 && length(neg) == 2 &&
      !is.same(pos) && all.same(neg)) {
    return("shift")
  }
  "custom"
}

# names are either "<pair_var>=<level>" (factor) or "<numeric_var>" (continuous)
#' @keywords internal
getContrastInfo <- function(contrast.vec) {
    stopifnot(length(contrast.vec) > 0, !is.null(names(contrast.vec)))
    nm <- names(contrast.vec)
    ## ---------- FACTOR CASE: "<var>_pair<level>" (new pair design style) ----------
    # e.g., "group_pairA_and_A", "group_pairB_and_B", "group_pairA_and_B"
    # We want:
    #   var = "group_pair"
    #   w   = weights with names = "A_and_A", "B_and_B", "A_and_B"
    #
    # Try to split at the "_pair" boundary:
    m <- regexec("^(.+?_pair)(.+)$", nm)
    parts <- regmatches(nm, m)
    
    if (all(lengths(parts) == 3L)) {
        var.prefix <- parts[[1]][2]  # "group_pair"
        lev <- vapply(parts, function(p) p[3], character(1))
        w   <- setNames(as.numeric(contrast.vec), lev)
        return(list(kind = "factor", var = var.prefix, w = w, center = NULL))
    }
    
    ## ---------- CONTINUOUS CASE ----------
    nz <- nm[abs(contrast.vec) > 0]
    stopifnot(length(nz) == 1)
    ctr.center <- attr(contrast.vec, "center")  # for continuous contrasts
    list(
        kind   = "continuous",
        var    = nz,
        w      = setNames(as.numeric(contrast.vec[nz]), nz),
        center = ctr.center
    )
}

# x0 from factor-contrast weights applied to the predictor x (mirrors total/shift)
#' @keywords internal
contrastCenterX0 <- function(x, pair.level, w, tol = 1e-8) {
  pos.lv <- names(w)[w >  tol]
  neg.lv <- names(w)[w < -tol]
  if (!length(pos.lv)) stop("Factor centering contrast has no positive weights.")
  denom <- sum(w[pos.lv]); if (denom <= 0) stop("Sum of positive weights must be > 0.")

  mu.x.neg <- lapply(neg.lv, function(L) {
    idx <- (pair.level == L)
    if (!any(idx)) return(NA_real_)
    mean(x[idx], na.rm = TRUE)
  })
  names(mu.x.neg) <- neg.lv

  num <- sum(mapply(function(L, mu) { (-w[L]) * mu }, neg.lv, mu.x.neg))
  num / denom
}


# helper: detect factor vs numeric for a single column name in a data.frame
#' @keywords internal
isContinuousCol <- function(df, col) {
  isTRUE(col %in% names(df)) && is.numeric(df[[col]])
}
