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
                                     dist = "cor", dist.type = c("shift", "total", "var"), 
                                     gene.selection = "wilcox", perm.method = c("freedman-lane", "block"), 
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

                                dists = subsetDistanceMatrix(dist.mat, sample.groups, cross.factor=TRUE, build.df=FALSE)
                                cross.rows <- which(outer(rownames(dist.mat), rownames(dist.mat), function(a,b) sample.groups[a] != sample.groups[b]))

                                # model matrices and pairwise distances
                                sample.meta.df.ct <- sample.meta.df[rownames(dist.mat), , drop = FALSE]
                                px <- buildPairDesignMatrices(sample.meta = sample.meta.df.ct,
                                       triplet = c(contrast[1],contrast[3],contrast[2]),  # var, alt, ref
                                       dist.type = "shift", block.vars = ifelse(!is.null(block.id), paste0(block.id, "_pair"), NULL))

                                y <- vectorizeLowerTri(dist.mat, px$pairs)

                                # Linear models and Permutations
                                res <- performLMPermutations(y=y, x=px, n.permutations=n.permutations, perm.method=perm.method, 
                                                             block=TRUE, core.only=TRUE, cross.rows= cross.rows)

                                # R^2 for each term
                                groups <- makeGroupsPair(px$model$F)
                                r2.df <- estimateR2PerTerm(px$model$F, y, groups)
                                r2.df <- r2.df[-grepl("Intercept", r2.df$term), ]

                                list(res = res, dist.mat = dist.mat, dists = res$yadj, r2 = r2.df, model.diag=px$model.diag)
                 }, progress = verbose, n.cores = n.cores.inner, mc.preschedule = TRUE, mc.allow.recursive = TRUE, fail.on.error = TRUE)

  if (verbose) message("Done!\n")
  
  summary <- summarizeLMResults(res.per.type, adj.method=p.adjust.method)

  list(results= summary$results, pvalues= summary$pvalues, padjust= summary$padjust, perm.method = perm.method, dist.type = dist.type, dists.per.type = summary$dists.per.type, 
       p.dist.info = summary$p.dist.info, r2 = summary$r2, cell.groups=cell.groups, model.diagnostics = lapply(res.per.type, `[[`, "model.diag"))

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
performLMPermutations <- function(y, x, n.permutations = 1000, perm.method = "freedman-lane", block = TRUE, core.only = TRUE, ids = NULL, cross.rows = NULL) {
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
  if (perm.method == "block") {
    fit.full <- lm.fit(x$model$F, y[rows])
    betaF <- fit.full$coefficients; betaF[!is.finite(betaF)] <- 0
    stat.obs <- filterContrastFit(betaF, x$model$contrast.F)

    # add fit for yhat values. 

    yhat.full <- as.numeric(x$model$F %*% betaF)  # fitted values in full space: distances with nuisance effects removed
    y.resid <- y[rows] - yhat.full # unexplained (full space) for qc/model check
    yadj <- yhat.full[cross.rows] # still has contrast non-participating signal + random noise

    for (b in seq_len(B)) {
      y.perm <- permuteWithinBlocks(y, blocks.use)
      fit.perm <- lm.fit(x = x$model$F, y = y.perm)
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

    yadj <- y.r.use[cross.rows] # distances with all nuisance effects removed; still contain the contrast signal and random noise
    y.resid <- y.r.use - as.numeric(X.r.use %*% fit.core$coefficients)  # unexplained (FL space) for qc/model check

    for (b in seq_len(B)) {
      y.perm <- sample(y.r.use)
      fit.perm <- lm.fit(X.r.use, y.perm)
      stats.perm[b] <- filterContrastFit(fit.perm$coefficients, x$model$contrast.X)
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
                   #k = length(cols.g),
                   #rank.minus = rank.other,
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
  as.data.frame(out, stringsAsFactors = TRUE)
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
  list(dists = dists)
}


#' @keywords internal
fitCellTypePairwiseDistances <- function(dist.df, ct=NULL, formula = NULL, contrast = NULL, sample.meta.df = NULL, diff.term.map = NULL, cm.norm, 
                                         sample.groups, top.n.genes = NULL, gene.selection = "wilcox", norm.type = NULL, r.type = NULL, n.pcs = NULL,
                                         dist = NULL, ref.level = NULL, target.level = NULL, n.permutations = 1000, n.cores.inner = 1, 
                                         return.all.cov = FALSE, verbose = TRUE, ...) {
  # Model formula from *_diff columns generated by computePairwiseDistances 
  diff.cols <- grep("_diff$", names(dist.df), value = TRUE)
  if (!length(diff.cols)) stop("No *_diff columns found in dist.df for cell type: ", ct)
  model.formula <- stats::as.formula(paste("Distance ~", paste(diff.cols, collapse = " + ")))

  # Design matrix & response 
  tr <- terms(model.formula)
  X  <- stats::model.matrix(tr, data = dist.df)  # includes intercept
  y  <- dist.df$Distance
  fit.obs <- lm_fit(X, y) 

  # normalize shapes & strip intercept in outputs 
  coef.names <- colnames(X)
  if (is.null(dim(fit.obs$coefficients))) {
    fit.obs$coefficients <- matrix(fit.obs$coefficients, ncol = 1,
                                   dimnames = list(coef.names, "Estimate"))
  } else rownames(fit.obs$coefficients) <- coef.names

  if (is.null(dim(fit.obs$t_values))) {
    fit.obs$t_values <- matrix(fit.obs$t_values, ncol = 1,
                               dimnames = list(coef.names, "t_value"))
  } else rownames(fit.obs$t_values) <- coef.names

  names(fit.obs$p_values) <- coef.names

  non.intercept <- setdiff(coef.names, "(Intercept)") # remove intercept 
  fit.obs$coefficients <- fit.obs$coefficients[non.intercept, , drop = FALSE]
  fit.obs$t_values     <- fit.obs$t_values[non.intercept, , drop = FALSE]
  fit.obs$p_values     <- fit.obs$p_values[non.intercept]
  coef.names <- non.intercept

  contrast.var <- contrast[1]
  cont <- grep(paste0("^", contrast.var, ".*_diff$"), coef.names, value = TRUE)
  coef.name <- if (length(cont)) cont[1] else coef.names[1]
  obs.stat <- if (!return.all.cov) fit.obs$coefficients[coef.name, ] else fit.obs$coefficients

  # Group columns into terms for R2 
  non.intercept <- setdiff(colnames(X), "(Intercept)")
  if (is.null(names(diff.term.map))) { # needed in case of interaction terms
    stop("diff.term.map must be a named character vector mapping *_diff columns to terms")
  }
  names(diff.term.map) <- make.names(names(diff.term.map)) # needed if there are spaces or special characters
  diff.term.map <- diff.term.map[non.intercept]
  groups <- split(match(names(diff.term.map), colnames(X)), diff.term.map)
  rss.full <- fit.obs$rss

  partial.r2 <- setNames(numeric(length(groups)), names(groups))
  for (nm in names(groups)) {
    drop.idx <- groups[[nm]]
    X.red <- X[, setdiff(seq_len(ncol(X)), drop.idx), drop = FALSE]
    if (!"(Intercept)" %in% colnames(X.red)) X.red <- cbind("(Intercept)" = 1, X.red)
    fit.red <- try(lm_fit(X.red, y), silent = TRUE)
    if (!inherits(fit.red, "try-error")) {
      rss.red <- fit.red$rss
      partial.r2[nm] <- pmax(0, (rss.red - rss.full) / pmax(rss.red, .Machine$double.eps))
    } else {
      partial.r2[nm] <- NA_real_
    }
  }
  ord <- c(sort(setdiff(names(partial.r2), grep(":", names(partial.r2), value = TRUE))),
           sort(grep(":", names(partial.r2), value = TRUE)))
  partial.r2.df <- as.data.frame(t(partial.r2[ord]), check.names = FALSE)

  # Permutations (shuffle labels, recompute distances & refit)
  samples <- unique(c(dist.df$Sample1, dist.df$Sample2))
  contrast.var <- contrast[1]

  perm.results <- plapply(seq_len(n.permutations), function(i) {
    valid.samples <- intersect(samples, rownames(sample.meta.df))
    if (!length(valid.samples)) return(NULL)

    # 1) shuffle labels (preserve factor levels)
    sg.shuff <- sample(sample.meta.df[valid.samples, contrast.var])
    names(sg.shuff) <- valid.samples
    smp.perm <- sample.meta.df[valid.samples, , drop = FALSE]
    smp.perm[[contrast.var]] <- factor(sg.shuff, levels = levels(sample.meta.df[[contrast.var]]))

    # 2) recompute distances under permutation
    dm.perm <- estimateExpressionShiftsForCellType(cm.norm, sample.groups = sg.shuff, dist = dist, n.pcs = n.pcs, top.n.genes = top.n.genes,
                                                       gene.selection = gene.selection, ...)

    # 3) normalize by same scheme
    if (norm.type == "both") {
      sg1 <- levels(sg.shuff)[1]
      m1 <- outer(sg.shuff, sg.shuff, function(a, b) (a == sg1) & (b == sg1))
      m2 <- outer(sg.shuff, sg.shuff, function(a, b) (a != sg1) & (b != sg1))
      diag(m1) <- NA; diag(m2) <- NA
      norm.const <- (stats::median(dm.perm[m1], na.rm = TRUE) + stats::median(dm.perm[m2], na.rm = TRUE)) / 2
      dm.perm.norm <- dm.perm - norm.const
    } else if (norm.type == "ref") { # TODO: check normalization method
    if (is.null(ref.level))  ref.level  <- levels(sample.groups)[1]
    if (is.null(target.level)) { target.level <- setdiff(levels(sample.groups), ref.level)[1] }

    rr <- outer(sample.groups, sample.groups, function(a, b) a == ref.level    & b == ref.level)
    tt <- outer(sample.groups, sample.groups, function(a, b) a == target.level & b == target.level)
    diag(rr) <- FALSE
    diag(tt) <- FALSE

    med <- function(x) stats::median(x, na.rm = TRUE)
    mR  <- if (any(rr)) med(dm.perm[rr]) else 0
    # ref group scale (MAD/SD, then to 1)
    sR <- if (any(rr)) stats::mad(dm.perm[rr], center = 0, constant = 1, na.rm = TRUE) else 0
    if (!is.finite(sR) || sR <= 0) sR <- stats::sd(dm.perm[rr], na.rm = TRUE)
    if (!is.finite(sR) || sR <= 0) sR <- 1

    if (any(tt)) { # Target group scale
    sT <- stats::mad(dm.perm[tt], center = 0, constant = 1, na.rm = TRUE)
    if (!is.finite(sT) || sT <= 0) sT <- stats::sd(dm.perm[tt], na.rm = TRUE)
    if (!is.finite(sT) || sT <= 0) sT <- sR
    } else {
    sT <- sR
    }
    s.vec <- ifelse(sg == ref.level, sR, sT)
    eps <- 1e-8
    scale.mat  <- outer(s.vec, s.vec, function(a, b) sqrt(pmax(a, eps) * pmax(b, eps)))
    center.mat <- matrix(mR, nrow(dm.perm), ncol(dm.perm))
    dm.perm.norm <- (dm.perm - center.mat) / scale.mat
  } else {
    dm.perm.norm <- dm.perm
  }
    
    # 4) pairwise design + fit
      model.mat.perm <- buildModelMatrix(sample.meta = smp.perm, formula = formula, contrast = contrast)
      dist.df.perm <- computePairwiseDistances(dm.perm.norm, model.mat.perm, valid.samples)

    tt.perm <- terms(model.formula)
    X.perm  <- stats::model.matrix(tt.perm, data = dist.df.perm)
    y.perm  <- dist.df.perm$Distance
    fit.perm <- try(lm_fit(X.perm, y.perm), silent = TRUE)
    if (inherits(fit.perm, "try-error")) {
      return(list(coef = setNames(rep(NA_real_, length(coef.names)), coef.names),
                  r2   = setNames(rep(NA_real_, length(groups)), names(groups))))
    }
    # Remove intercept
    non.intercept.perm <- setdiff(colnames(X.perm), "(Intercept)")
    if (is.null(dim(fit.perm$coefficients))) {
    fit.perm$coefficients <- matrix(fit.perm$coefficients, ncol = 1,
                                   dimnames = list(colnames(X.perm), "Estimate"))
     } else rownames(fit.perm$coefficients) <- colnames(X.perm)
    coef.mat.perm <- fit.perm$coefficients[non.intercept.perm, , drop = FALSE]
    if (is.null(dim(coef.mat.perm)))
      coef.mat.perm <- matrix(coef.mat.perm, ncol = 1,
                               dimnames = list(non.intercept.perm, "Estimate"))
    # Align group indices using observed column names
    rss.full.p <- fit.perm$rss
    Xn <- colnames(X.perm)
    groups.perm <- lapply(groups, function(idx) match(colnames(X)[idx], Xn))
    groups.perm <- lapply(groups.perm, function(v) v[!is.na(v)])

    r2.perm <- setNames(numeric(length(groups.perm)), names(groups.perm))
    for (nm in names(groups.perm)) {
      idx.drop <- groups.perm[[nm]]
      if (!length(idx.drop)) { r2.perm[nm] <- NA_real_; next }
      Xr <- X.perm[, setdiff(seq_len(ncol(X.perm)), idx.drop), drop = FALSE]
      if (!"(Intercept)" %in% colnames(Xr)) Xr <- cbind("(Intercept)" = 1, Xr)
      fr <- try(lm_fit(Xr, y.perm), silent = TRUE)
      if (inherits(fr, "try-error")) {
        r2.perm[nm] <- NA_real_
      } else {
        rssr <- fr$rss
        r2.perm[nm] <- pmax(0, (rssr - rss.full.p) / pmax(rssr, .Machine$double.eps))
      }
    }

    list(
      coef = setNames(as.vector(coef.mat.perm[, 1]), rownames(coef.mat.perm)),
      r2   = r2.perm
    )
  }, n.cores = n.cores.inner, progress = verbose, fail.on.error = FALSE)
  

  perm.results <- perm.results[!vapply(perm.results, is.null, FALSE)]
  # Combine permutation results
  if (!return.all.cov) {
    perm.coefs <- unlist(lapply(perm.results, function(x) x$coef[coef.name]), use.names = FALSE)
    expected.stat <- mean(perm.coefs, na.rm = TRUE)
    centered.coefficients <- obs.stat - expected.stat
    pvalue <- (sum(abs(perm.coefs) >= abs(centered.coefficients), na.rm = TRUE) + 1) /
              (sum(!is.na(perm.coefs)) + 1)
    #names(pvalue) <- coef.name
    perm.stats <- perm.coefs
  } else {
    # matrix (#perms x #coefs)
    perm.coef.mat <- do.call(rbind, lapply(perm.results, `[[`, "coef"))
    miss <- setdiff(coef.names, colnames(perm.coef.mat))
    if (length(miss)) perm.coef.mat[, miss] <- NA_real_
    perm.coef.mat <- perm.coef.mat[, coef.names, drop = FALSE]

    expected.stat <- colMeans(perm.coef.mat, na.rm = TRUE)
    centered.coefficients <- obs.stat - expected.stat

    pvalue <- sapply(seq_along(coef.names), function(j) {
      (sum(abs(perm.coef.mat[, j]) >= abs(centered.coefficients[j, 1]), na.rm = TRUE) + 1) /
      (sum(!is.na(perm.coef.mat[, j])) + 1)
    })
    names(pvalue) <- coef.names
    perm.stats <- perm.coef.mat
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
