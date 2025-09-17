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
#' @param perm.method permutation method: "freedman-lane" (default) or "full"
#' @param block.id Optional column in `sample.meta` specifying blocks for restricted randomizations
#' @param n.cores number of cores (default=1)
#' @param verbose (default=FALSE)
#' @return List of distances per cell type, distance matrices, sample groups, cell types, pvalues, and adjusted p-values.
#' @export
estimateExpressionChange <- function(cm.per.type, sample.meta, cell.groups, sample.groups, sample.id = NULL,
                                     sample.per.cell, formula = NULL, contrast = NULL, block.id = NULL,
                                     dist = "cor", dist.type = c("shift", "total", "var"), 
                                     gene.selection = "wilcox", perm.method = c("freedman-lane", "full"), 
                                     n.permutations = 1000, p.adjust.method = "BH", trim = 0.2, 
                                     top.n.genes = NULL, n.pcs = NULL, n.cores = 1, verbose = TRUE, ...) {
  dist.type <- match.arg(dist.type)
  dist <- parseDistance(dist, top.n.genes = top.n.genes, n.pcs = n.pcs)
  perm.method <- match.arg(perm.method)

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

                                # model matrices and pairwise distances
                                sample.meta.df.ct <- sample.meta.df[rownames(dist.mat), , drop = FALSE]
                                px <- buildPairDesignMatrices(sample.meta = sample.meta.df.ct,
                                       triplet = c(contrast[1],contrast[3],contrast[2]),  # var, alt, ref
                                       dist.type = "shift", block.vars = ifelse(!is.null(block.id), paste0(block.id, "_pair"), NULL))

                                y <- vectorizeLowerTri(dist.mat, px$pairs)

                                # Linear models and Permutations
                                res <- performLMPermutations(y=y, x=px, n.permutations=n.permutations, perm.method=perm.method, 
                                                             block=TRUE, core.only=TRUE)

                                # R^2 for each term
                                groups <- makeGroupsPair(px$model$F)
                                r2.df <- estimateR2PerTerm(px$model$F, y, groups)
                                r2.df <- r2.df[-grepl("Intercept", r2.df$term), ]

                                list(res = res, dist.mat = dist.mat, dists = res$yadj, r2 = r2.df)
                 }, progress = verbose, n.cores = n.cores.inner, mc.preschedule = TRUE, mc.allow.recursive = TRUE, fail.on.error = TRUE)

  if (verbose) message("Done!\n")
  
  summary <- summarizeLMResults(res.per.type, adj.method=p.adjust.method)

  list(results= summary$results, pvalues= summary$pvalues, padjust= summary$padjust, perm.method = perm.method, dist.type = dist.type, dists.per.type = summary$dists.per.type, 
       p.dist.info = summary$p.dist.info, r2 = summary$r2, cell.groups=cell.groups)
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
performLMPermutations <- function(y, x, n.permutations = 1000, perm.method = "freedman-lane", block = TRUE, core.only = TRUE, ids = NULL) {
  # blocks
  rows <- rep(TRUE, nrow(x$model$F))
  blocks.use <- NULL
  if (block && core.only) {
    rows <- if (!is.null(x$model$core.rows)) x$model$core.rows else rows
    blocks.use <- if (is.null(x$model$blocks)) NULL else droplevels(x$model$blocks[rows])
  }

  B <- n.permutations
  stats.perm <- numeric(B)
  adj <- list()

  # fitting within permutation loops
  if (perm.method == "full") {
    fit.full <- lm.fit(x$model$F, y[rows])
    betaF <- fit.full$coefficients; betaF[!is.finite(betaF)] <- 0
    stat.obs <- filterContrastFit(betaF, x$model$contrast.F)

    yhat.full <- as.numeric(x$model$F %*% betaF)  # fitted values in full space: distances with nuisance effects removed
    y.resid <- y[rows] - yhat.full # unexplained (full space) for qc/model check
    yadj <- yhat.full # still has contrast non-participating signal + random noise

    for (b in seq_len(B)) {
      fit.perm <- lm.fit(x = x$model$F, y = permuteWithinBlocks(y, blocks.use))
      stats.perm[b] <- filterContrastFit(fit.perm$coefficients, x$model$contrast.F)
    }
    stats.perm <- stats.perm[is.finite(stats.perm)]
    stopifnot(is.finite(stat.obs))
    pval <- (1 + sum(abs(stats.perm) >= abs(stat.obs))) / (length(stats.perm) + 1)

  } else if (perm.method == "freedman-lane") { # (residualization)
    rz <- residualizeForFL(y, x$model$qrZ, x$model$X)
    fit.core <- lm.fit(rz$X.r, rz$y.r)
    betaX <- fit.core$coefficients; betaX[!is.finite(betaX)] <- 0
    stat.obs <- filterContrastFit(betaX, x$model$contrast.X)

    X.r.use <- rz$X.r[rows, , drop = FALSE]
    y.r.use <- rz$y.r[rows]

    yadj <- y.r.use # distances with all nuisance effects removed; still contain the contrast signal and random noise
    y.resid <- y.r.use - as.numeric(X.r.use %*% fit.core$coefficients)  # unexplained (FL space) for qc/model check

    for (b in seq_len(B)) {
      y.perm <- permuteWithinBlocks(y.r.use, blocks.use)
      fit.perm <- lm.fit(X.r.use, y.perm)
      stats.perm[b] <- as.numeric(crossprod(fit.perm$coefficients, x$model$contrast.X))
    }

    stats.perm <- stats.perm[is.finite(stats.perm)]
    stopifnot(is.finite(stat.obs))
    pval <- (1 + sum(abs(stats.perm) >= abs(stat.obs))) / (length(stats.perm) + 1)

  } else {
    stop("Unknown permutation method: ", perm.method)
  }

  stat.obs <- stat.obs - mean(stats.perm)

  list(stat.obs = stat.obs, stats.perm = stats.perm, pval = pval, yadj = yadj, y.resid = y.resid)
}

#' To remove NA and non-overlapping names from coef and contrast vectors
#' @keywords internal
filterContrastFit <- function(coef.vec, contrast.vec) {
  common <- intersect(names(coef.vec), names(contrast.vec))
  if (!length(common)) return(NA_real_)
  a <- coef.vec[common]; a[!is.finite(a)] <- 0
  c <- contrast.vec[common]
  sum(a * c)
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


###################### Pair Metadata functions #########################

# Prepare pair-level design for modeling lower-tri distances
#
# Builds per-pair metadata (factor pairs "AA","AB",...; numeric diffs |xi-xj|),
# constructs a pair-level contrast from a sample-level triplet (var, alt, ref)
# according to dist.type, and calls `prepare_design()` to produce F, X, Z, qrZ,
# contrast vectors, core.rows, and blocks. Returns the (i,j) index table so you
# can vectorize ANY distance matrix later in the same order.
#
# @param meta data.frame of per-sample covariates; row order defines sample order
# @param triplet character triplet c(var, alt, ref) (or named var/alt/ref)
# @param dist.type one of "shift","total","var" (pair-level contrast form)
# @param block.vars optional character vector to override permutation blocks
# @param na.action NA handler for model.matrix
# @return list with:
#   - pairs: data.frame of (i,j) indices (1-based) in lower-tri order
#   - pair.meta: data.frame of per-pair covariates
#   - prep: design list from buildDesignMatrices() (F, X, Z, qrZ, contrast.F/.X, ...)
#' @keywords internal
buildPairDesignMatrices <- function(sample.meta, triplet,
                                dist.type = c("shift","total","var"),
                                block.vars = NULL,
                                na.action = stats::na.pass) {
  stopifnot(is.data.frame(sample.meta))
  dist.type <- match.arg(dist.type)
  
  # 1) Indices for lower triangle, derived from nrow(meta)
  n   <- nrow(sample.meta)
  idx <- lowerTriIndices(n)

  # 2) Pair-level metadata from sample-level meta
  pair.meta <- pairifyMeta(sample.meta, idx)

  # 3) Triplet (var, alt, ref) → pair-level contrast weights over AA/AB/BB
  tr <- parseTriplet(triplet)                # var, alt, ref
  if (!(tr$var %in% names(sample.meta)))
    stop("Triplet 'var' not found in meta: ", tr$var)

  f.levels <- levels(factor(sample.meta[[tr$var]]))
  if (!all(c(tr$ref, tr$alt) %in% f.levels))
    stop("alt/ref levels not found in meta[['", tr$var, "']]: ",
         tr$alt, ", ", tr$ref)

  weights <- pairContrastWeights(tr$ref, tr$alt, dist.type, level.order = f.levels)

  # 4) Build pair-level variable name for the contrasted factor (e.g., "group_pair")
  pair.var <- paste0(tr$var, "_pair")
  if (!(pair.var %in% names(pair.meta)))
    stop("Pair variable not found in pair_meta: ", pair.var,
         " (did '", tr$var, "' exist and was it factor-like?)")

  # 5) Call buildDesignMatrices() on pair-level metadata
  prep <- buildDesignMatrices(
    sample.meta       = pair.meta,
    contrast   = list(var = pair.var, weights = weights),
    #nuisance   = setdiff(names(pair.meta), pair.var),  # core is implicitly pair.var
    core.extra = NULL,
    block.vars = block.vars,
    na.action  = na.action
  )

  list(pairs = idx, pair.meta = pair.meta, model = prep)
}

#' Return lower-triangular (row, col) index pairs for an n x n matrix
#' @param n integer number of samples
#' @return data.frame with columns i, j (1-based), in row-major lower-tri order
#' @keywords internal
lowerTriIndices <- function(n) {
  ij <- which(lower.tri(matrix(NA_real_, n, n)), arr.ind = TRUE)
  data.frame(i = ij[,1], j = ij[,2])
}

#' Build pair-level metadata from sample-level metadata in lower-tri order
#' @param meta data.frame of per-sample covariates (rows correspond to d's order)
#' @param idx data.frame with columns i, j (lower-tri indices)
#' @return data.frame of per-pair covariates
#' @keywords internal
pairifyMeta <- function(meta, idx) {
  stopifnot(is.data.frame(meta), all(c("i","j") %in% names(idx)))
  n <- nrow(meta)
  stopifnot(all(idx$i >= 1 & idx$i <= n & idx$j >= 1 & idx$j <= n & idx$i > idx$j))
  
  out <- list()
  for (v in names(meta)) {
    x <- meta[[v]]
    if (is.factor(x) || is.character(x) || is.logical(x)) {
      f <- factor(x)
      out[[paste0(v, "_pair")]] <- factor(pairCodeFactor(f, idx$i, idx$j))
    } else if (is.numeric(x)) {
      out[[paste0(v, "_diff")]] <- pairDiffNumeric(x, idx$i, idx$j)
    } else {
      # unsupported type → coerce to factor then pair-code
      f <- factor(as.character(x))
      out[[paste0(v, "_pair")]] <- factor(pairCodeFactor(f, idx$i, idx$j))
    }
  }
  as.data.frame(out, stringsAsFactors = TRUE)
}

#' Symmetric pair-code for a factor: "AA", "AB", "BB" (AB ≡ BA)
#' Uses the factor's level order to define A<B for concatenation.
#' @param f factor vector (length n)
#' @param i integer indices for first sample in each pair
#' @param j integer indices for second sample in each pair
#' @return character vector of codes ("AA","AB",...)
#' @keywords internal
pairCodeFactor <- function(f, i, j) {
  lev <- levels(f); fi <- as.integer(f[i]); fj <- as.integer(f[j])
  check.na <- is.na(fi) | is.na(fj)
  lo <- pmin(fi, fj); hi <- pmax(fi, fj)
  code <- paste0(lev[lo], lev[hi])
  code[check.na] <- NA_character_
  code
}

#' Symmetric numeric pair difference |x_i - x_j|
#' @param x numeric vector (length n)
#' @param i, j integer index vectors
#' @return numeric vector of absolute differences
#' @keywords internal
pairDiffNumeric <- function(x, i, j) {
  xi <- x[i]; xj <- x[j]
  out <- abs(xi - xj)
  out
}

#' Parse a sample-level "triplet" contrast specification
#' @param triplet list(var="<factor>", ref="<A>", alt="<B>")
#' @return list(var, ref, alt)
#' @keywords internal
parseTriplet <- function(triplet) {
  stopifnot(length(triplet)==3)
  triplet <- as.list(triplet)
  names(triplet) <- c('var','alt','ref')
  if (triplet$ref == triplet$alt) stop("ref and alt must be different.")
  list(var = as.character(triplet$var),
       alt = as.character(triplet$alt),
       ref = as.character(triplet$ref))
}

# Translate (ref, alt) + dist.type into pair-level contrast weights
# - level.order: optional character vector giving the factor's level order
#   (used to decide the label for the cross pair, e.g., "AB" vs "BA")
# - Returns a named numeric vector, e.g. c("AA"=-0.5, "AB"=1, "BB"=-0.5)
#' @keywords internal
pairContrastWeights <- function(ref, alt, type = c("shift","total","var"), level.order = NULL) {
  type <- match.arg(type)
  
  # Cross label uses factor-level order if provided; otherwise alphabetical
  two <- c(ref, alt)
  if (!is.null(level.order)) {
    ord <- two[order(match(two, level.order))]
  } else {
    ord <- sort(two)  # fallback: alphabetical
  }
  AB <- paste0(ord[1], ord[2])
  
  AA <- paste0(ref, ref)
  BB <- paste0(alt, alt)
  
  if (type == "shift") {
    w <- c(-0.5, 1, -0.5)
    names(w) <- c(AA, AB, BB)
    return(w)
  } else if (type == "total") {
    w <- c(-1, 1)
    names(w) <- c(AA, AB)
    return(w)
  } else { # type == "var"
    w <- c(-1, 1)
    names(w) <- c(AA, BB)
    return(w)
  }
}

## Generic grouping for pair-level design matrices
## - collapses all "Var_pair*" columns into one group per Var
## - maps "Var_pair_contrast" to Var
## - keeps "(Intercept)" if present
## - anything else stays as-is (unless include.other=TRUE, then grouped under "Other")
#' @keywords internal
makeGroupsPair <- function(M, include.other = FALSE) {
  stopifnot(!is.null(colnames(M)))
  cn <- colnames(M)
  groups <- list()
  if ("(Intercept)" %in% cn) {
    groups$Intercept <- "(Intercept)"
  }
  # Handle all _pair terms
  pair.terms <- grep("_pair", cn, value = TRUE)
  if (length(pair.terms)) {
    # Extract the variable name before "_pair"
    varnames <- sub("_pair.*$", "", pair.terms)
    for (v in unique(varnames)) {
      cols <- pair.terms[varnames == v]
      # if there's a collapsed contrast, rename group nicely
      if (any(grepl("_pair_contrast$", cols))) {
        groups[[v]] <- cols[grepl("_pair_contrast$", cols)]
      } else {
        groups[[v]] <- cols
      }
    }
  }

  # Leftovers (e.g. numeric pair variables not matching _pair)
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
###################### Additional fitting functions #########################

## Compute unique R^2 for each term/group in a linear model
# - F: full model matrix (n x p)
# - y: numeric response vector (length n)
# - groups: list of character vectors, each naming columns of F that form a term/group. compute using .makeGroupsPair() or similar.
#' @keywords internal
estimateR2PerTerm <- function(F, y, groups = NULL) {
    stopifnot(is.matrix(F) || is.data.frame(F))
    F <- as.matrix(F)
    y <- as.numeric(y)
    has.intercept <- "(Intercept)" %in% colnames(F)
    
    # Total Sum of Squares (center if intercept, else uncentered)
    TSS <- if (has.intercept) sum((y - mean(y))^2) else sum(y^2)
    if (TSS <= 0) stop("TSS is zero; cannot compute R^2.")

    # Full model SSE and rank (for reference)
    SSE.full <- getSSE(F, y)
    R2.full  <- 1 - SSE.full / TSS
    rank.full <- qr(F)$rank

    # Default groups = each column as its own term
    if (is.null(groups)) {
        groups <- as.list(colnames(F))
        names(groups) <- colnames(F)
    } else {
        # sanity: drop any names not present
        groups <- lapply(groups, function(v) intersect(v, colnames(F)))
        ok <- lengths(groups) > 0
        if (!all(ok)) warning("Some groups had no matching columns in F and were dropped.")
        groups <- groups[ok]
    }
    
    # Compute unique R^2 by refitting without each group's columns
    out <- lapply(names(groups), function(g) {
        cols.g <- groups[[g]]
        cols.other <- setdiff(colnames(F), cols.g)
        if (!length(cols.other)) {
            # Removing the only columns leaves empty model: SSE_other = TSS (if no intercept),
            # or SSE w.r.t intercept-only. Handle gracefully:
            F.other <- matrix(1, nrow(F), 1); colnames(F.other) <- "(Intercept*)"
            SSE.other <- if (has.intercept) sum((y - mean(y))^2) else sum(y^2)
            rank.other <- 1L
        } else {
            F.other <- F[, cols.other, drop = FALSE]
            SSE.other <- getSSE(F.other, y)
            rank.other <- qr(F.other)$rank
        }
        R2.unique <- max(0, (SSE.other - SSE.full) / TSS)  # clip tiny negatives from numeric noise
        data.frame(term = g,
                   k = length(cols.g),
                   rank.minus = rank.other,
                   R2.unique = R2.unique,
                   R2.full = R2.full,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

# Fit and get sum of squares (SSE) robustly (treat NAs in coef as 0 via lm.fit behavior)
#' @keywords internal
getSSE <- function(X, y) {
    fit <- lm.fit(X, y)
    bhat <- fit$coefficients; bhat[!is.finite(bhat)] <- 0
    res  <- y - as.numeric(X %*% bhat)
    sum(res^2)
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
