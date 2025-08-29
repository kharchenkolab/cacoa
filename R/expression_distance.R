#' @importFrom sccore checkPackageInstalled
NULL

#' Expression shift magnitudes per cluster between conditions
#'
#' @param cm.per.type List of normalized count matrices per cell type
#' @param sample.meta Data frame with sample-level covariates, row names are sample IDs
#' @param sample.per.cell Named vector with cell names indicating sample/condition
#' @param formula Formula specifying the model to fit
#' @param contrast Character vector of length 3 specifying the contrast to test: c("condition_variable", "reference_level", "target_level")
#' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
#' @param cell.groups Named clustering/annotation factor with cell names
#' @param dist what distance measure to use: 'JS' - Jensen-Shannon divergence (default), 'cor' - Pearson's linear correlation on log transformed values
#' @param return.all.cov Logical indicating whether to return results for all covariates or only contrast variable (default=FALSE)
#' @param n.cores number of cores (default=1)
#' @param verbose (default=FALSE)
#' @return List of distances per cell type, distance matrices, sample groups, cell types, pvalues, and adjusted p-values.
#' @export
estimateExpressionChange <- function(cm.per.type, sample.groups, cell.groups, 
                                     sample.meta, sample.per.cell, sample.id = NULL, 
                                     formula = NULL, contrast = NULL, dist = "cor", dist.type = c("shift", "total", "var"), 
                                     ref.level = NULL, target.level = NULL, gene.selection = "wilcox", 
                                     n.permutations = 1000, return.all.cov = FALSE, p.adjust.method = "BH", trim = 0.2, 
                                     top.n.genes = NULL, n.pcs = NULL, n.cores = 1, verbose = TRUE, ...) {
  dist.type <- match.arg(dist.type)
  dist <- parseDistance(dist, top.n.genes = top.n.genes, n.pcs = n.pcs)
  norm.type <- ifelse(dist.type == "shift", "both", "ref")
  #r.type    <- ifelse(dist.type == "var", "target", "cross")

  cell.groups <- droplevels(factor(cell.groups))
  sample.type.table <- table(cell.groups, sample.per.cell[names(cell.groups)])

  contrast.var <- contrast[1]
  ref.level    <- contrast[2]
  target.level <- contrast[3]

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
                                # distances per sample
                                dist.mat <- estimateExpressionShiftsForCellType(cm.norm, sample.groups = sample.groups, dist = dist, n.pcs = n.pcs,
                                                                                top.n.genes = top.n.genes, gene.selection = gene.selection, ...)
                                attr(dist.mat, "n.cells") <- sample.type.table[ct, rownames(cm.norm)]

                                d <- estimateExpressionShiftsByDistMat(dist.mat = dist.mat, sample.groups = sample.groups, sample.meta.df = sample.meta.df, 
                                                                       formula = formula, contrast = contrast, norm.type = norm.type, 
                                                                       return.type = ifelse(dist.type == "var", "target", "cross"),
                                                                       ref.level = ref.level, target.level = target.level, verbose = verbose)

                                fit <- fitCellTypePairwiseDistances(dist.df = d$dist.df, ct = ct, formula = formula, contrast = contrast, sample.meta.df = sample.meta.df, 
                                                                    diff.term.map = d$diff.term.map %||% NULL, cm.norm = cm.norm, sample.groups = sample.groups, 
                                                                    top.n.genes = top.n.genes, norm.type = norm.type, r.type = ifelse(dist.type == "var", "target", "cross"), 
                                                                    n.pcs = n.pcs, dist = dist, ref.level = ref.level, target.level = target.level, return.all.cov = return.all.cov, 
                                                                    n.permutations = n.permutations, n.cores.inner = n.cores.inner, verbose = verbose, ...)

                                list(fit = fit, dist.mat = d$dist.mat.norm, dists = d$dist.df$Distance)
                            }, progress = verbose, n.cores = n.cores, mc.preschedule = TRUE, mc.allow.recursive = TRUE, fail.on.error = TRUE)

  if (verbose) message("Done!\n")

  fits           <- lapply(res.per.type, `[[`, "fit")
  p.dist.info    <- lapply(res.per.type, function(x) x$dist.mat)
  dists.per.type <- lapply(res.per.type, function(x) x$dists)

  if (!return.all.cov) {
    coefs.per.type <- do.call(rbind, lapply(fits, `[[`, "coefficients"))[, 1]
    names(coefs.per.type) <- gsub("_diff$", "", names(coefs.per.type))
    pvalue <- sapply(fits, `[[`, "pvalue")
    names(pvalue) <- gsub("_diff$", "", names(pvalue))
    obs.stat <- do.call(rbind, lapply(fits, `[[`, "obs.stat"))[, 1]
    names(obs.stat) <- gsub("_diff$", "", names(obs.stat))
    exp.stat <- sapply(fits, `[[`, "expected.stat")
    names(exp.stat) <- gsub("_diff$", "", names(exp.stat))
    padjust <- stats::p.adjust(pvalue, method = p.adjust.method)
    names(padjust) <- gsub("_diff$", "", names(padjust))
    perm.stats <- lapply(fits, `[[`, "perm.stats")
    partial.r2.df <- lapply(fits, `[[`, "partial.r2.df")
    partial.r2.perm.mean <- lapply(fits, `[[`, "partial.r2.perm.mean")
    partial.r2.pvalue <- lapply(fits, `[[`, "partial.r2.pvalue")
  } else {
    coefs.per.type <- lapply(fits, `[[`, "coefficients")
    all.coefs <- unique(unlist(lapply(coefs.per.type, rownames)))
    coefs.per.type <- do.call(rbind, lapply(names(coefs.per.type), function(nm) {
      x <- coefs.per.type[[nm]]
      coefs <- setNames(rep(NA_real_, length(all.coefs)), all.coefs)
      coefs[rownames(x)] <- as.vector(x)
      data.frame(t(coefs), row.names = nm, check.names = FALSE)
    }))
    colnames(coefs.per.type) <- gsub("_diff$", "", colnames(coefs.per.type))
    pvalue <- do.call(rbind, lapply(fits, `[[`, "pvalue"))
    colnames(pvalue) <- gsub("_diff$", "", colnames(pvalue))
    obs.stat <- lapply(fits, `[[`, "obs.stat")
    all.obs <- unique(unlist(lapply(obs.stat, rownames)))
    obs.stat <- do.call(rbind, lapply(names(obs.stat), function(nm) {
      x <- obs.stat[[nm]]
      obs <- setNames(rep(NA_real_, length(all.obs)), all.obs)
      obs[rownames(x)] <- as.vector(x)
      data.frame(t(obs), row.names = nm, check.names = FALSE)
    }))
    colnames(obs.stat) <- gsub("_diff$", "", colnames(obs.stat))
    exp.stat <- do.call(rbind, lapply(fits, `[[`, "expected.stat"))
    colnames(exp.stat) <- gsub("_diff$", "", colnames(exp.stat))
    padjust <- apply(pvalue, 2, function(x) stats::p.adjust(x, method = p.adjust.method))
    colnames(padjust) <- gsub("_diff$", "", colnames(padjust))
    perm.stats <- lapply(fits, `[[`, "perm.stats")

    partial.r2.df <- do.call(rbind,lapply(fits, `[[`, "partial.r2.df"))
    partial.r2.perm.mean <- lapply(fits, `[[`, "partial.r2.perm.mean")
    partial.r2.pvalue <- do.call(rbind, lapply(fits, `[[`, "partial.r2.pvalue"))
  }

  list(dists.per.type = dists.per.type, p.dist.info = p.dist.info, sample.groups = sample.groups, coefs.per.type = coefs.per.type, 
       partial.r2.df = partial.r2.df, obs.stat = obs.stat, exp.stat = exp.stat, cell.groups = cell.groups, pvalues = pvalue, 
       padjust = padjust, perm.stat = perm.stats, partial.r2.perm.mean = partial.r2.perm.mean, partial.r2.pvalue = partial.r2.pvalue, 
       return.all.cov = return.all.cov, contrast = contrast, formula = formula)
}


#' @keywords internal
estimateExpressionShiftsForCellType <- function(cm.norm, sample.groups, dist, top.n.genes=NULL, n.pcs=NULL,
                                                gene.selection="wilcox", exclude.genes=NULL) {
  if (!is.null(top.n.genes)) {
    sel.genes <- filterGenesForCellType(
      cm.norm, sample.groups=sample.groups, top.n.genes=top.n.genes, gene.selection=gene.selection,
      exclude.genes=exclude.genes
    )
    cm.norm <- cm.norm[,sel.genes,drop=FALSE]
  }

  if (!is.null(n.pcs)) {
    min.dim <- min(dim(cm.norm)) - 1
    if (n.pcs > min.dim) {
      n.pcs <- min.dim
      warning("n.pcs is too large. Setting it to maximal allowed value ", min.dim)
    }

  cm.norm <- getTopPCs(cm.norm, n_pcs = n.pcs)
  rownames(cm.norm) <- n

  if (any(!is.finite(cm.norm))) {
    stop("PCA output contains non-finite values (NA/NaN/Inf).")
  }
  zero_var_cols <- which(apply(cm.norm, 2, var) < .Machine$double.eps)
  if (length(zero_var_cols) > 0) {
    warning("PCA output contains zero-variance components; removing them.")
    cm.norm <- cm.norm[, -zero_var_cols, drop = FALSE]
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

#' @keywords internal
estimateExpressionShiftsByDistMat <- function(dist.mat, sample.groups, formula, sample.meta.df = NULL, contrast = NULL, norm.type = c("both", "ref", "none"),
                                              ref.level = NULL, target.level = NULL, return.type = c("cross", "target"), verbose = FALSE) {
  norm.type <- match.arg(norm.type)
  return.type <- match.arg(return.type)

  if (((norm.type == "ref") || (return.type == "target")) && is.null(ref.level))
    stop("ref.level has to be provided for norm.type='ref' or return.type='target'")

  sample.groups %<>% .[rownames(dist.mat)]
    if (norm.type == "both") {
    sg1 <- levels(sample.groups)[1]
    m1 <- outer(sample.groups, sample.groups, function(a, b) (a == sg1) & (b == sg1))
    m2 <- outer(sample.groups, sample.groups, function(a, b) (a != sg1) & (b != sg1))
    diag(m1) <- NA; diag(m2) <- NA
    norm.const <- (stats::median(dist.mat[m1], na.rm = TRUE) + stats::median(dist.mat[m2], na.rm = TRUE)) / 2
    dist.mat.norm <- dist.mat - norm.const

  } else if (norm.type == "ref") {
    if (is.null(ref.level))  ref.level  <- levels(sample.groups)[1]
    if (is.null(target.level)) { target.level <- setdiff(levels(sample.groups), ref.level)[1] }

    rr <- outer(sample.groups, sample.groups, function(a, b) a == ref.level    & b == ref.level)
    tt <- outer(sample.groups, sample.groups, function(a, b) a == target.level & b == target.level)
    diag(rr) <- FALSE
    diag(tt) <- FALSE

    med <- function(x) stats::median(x, na.rm = TRUE)
    mR  <- if (any(rr)) med(dist.mat[rr]) else 0
    # ref group scale (MAD/SD, then to 1)
    sR <- if (any(rr)) stats::mad(dist.mat[rr], center = 0, constant = 1, na.rm = TRUE) else 0
    if (!is.finite(sR) || sR <= 0) sR <- stats::sd(dist.mat[rr], na.rm = TRUE)
    if (!is.finite(sR) || sR <= 0) sR <- 1

    if (any(tt)) { # Target group scale
    sT <- stats::mad(dist.mat[tt], center = 0, constant = 1, na.rm = TRUE)
    if (!is.finite(sT) || sT <= 0) sT <- stats::sd(dist.mat[tt], na.rm = TRUE)
    if (!is.finite(sT) || sT <= 0) sT <- sR
    } else {
    sT <- sR
    }
    s.vec <- ifelse(sg == ref.level, sR, sT)
    eps <- 1e-8
    scale.mat  <- outer(s.vec, s.vec, function(a, b) sqrt(pmax(a, eps) * pmax(b, eps)))
    center.mat <- matrix(mR, nrow(dist.mat), ncol(dist.mat))
    dist.mat.norm <- (dist.mat - center.mat) / scale.mat
  } else {
    dist.mat.norm <- dist.mat
  }

  samples <- rownames(dist.mat)

  if (!is.null(sample.meta.df)) {
    if (!is.data.frame(sample.meta.df) || is.null(rownames(sample.meta.df))) {
      stop("sample.meta.df must be a data.frame with samples as row names.")
    }
    sample.meta.df.ct <- sample.meta.df[samples, , drop = FALSE]
    model.mat <- buildModelMatrix(sample.meta = sample.meta.df.ct, formula = formula, contrast = contrast)
    # pairwise distances
    dist.df <- computePairwiseDistances(dist.mat, model.mat, samples)

    # backtrace *_diff columns to original terms (no intercept in model.mat)
    assign.vec  <- attr(model.mat, "assign")
    term.labels <- attr(terms(formula), "term.labels")
    col2term    <- setNames(term.labels[assign.vec], colnames(model.mat))
    diff.term.map <- setNames(unname(col2term), paste0(names(col2term), "_diff"))

    return(list(dist.df = dist.df, dist.mat.norm = dist.mat.norm, diff.term.map = diff.term.map))
  }

  # Fallback (no covariates)
  if (return.type == "cross") {
    dists <- subsetDistanceMatrix(dist.mat.norm, sample.groups = sample.groups, cross.factor = TRUE)
  } else {
    dists <- dist.mat.norm %>%
      .[(sample.groups[rownames(.)] != ref.level), (sample.groups[colnames(.)] != ref.level)] %>%
      .[upper.tri(.)]
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

  # Partial R² permutation aggregation
  perm.r2.mat <- do.call(rbind, lapply(perm.results, `[[`, "r2"))
  miss.r2 <- setdiff(names(groups), colnames(perm.r2.mat))
  if (length(miss.r2)) perm.r2.mat[, miss.r2] <- NA_real_
  perm.r2.mat <- perm.r2.mat[, names(groups), drop = FALSE]
  partial.r2.perm.mean <- colMeans(perm.r2.mat, na.rm = TRUE)
  partial.r2.pvalue <- sapply(names(groups), function(nm) {
    (sum(perm.r2.mat[, nm] >= partial.r2[nm], na.rm = TRUE) + 1) /
    (sum(!is.na(perm.r2.mat[, nm])) + 1)
  })

  list(model = fit.obs, pvalue = pvalue, obs.stat = obs.stat, coefficients = centered.coefficients, perm.stats = perm.stats,
       expected.stat = expected.stat, partial.r2.df = partial.r2.df, partial.r2.perm.mean = partial.r2.perm.mean, partial.r2.pvalue = partial.r2.pvalue)
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
filterGenesForCellType <- function(cm.norm, sample.groups, top.n.genes=500, gene.selection=c("wilcox", "var", "od"),
                                   exclude.genes=NULL) {
  gene.selection <- match.arg(gene.selection)

  if (gene.selection == "var") {
    sel.genes <- estimateExplainedVariance(cm.norm, sample.groups=sample.groups) %>%
      sort(decreasing=TRUE) %>% names()
  } else if (gene.selection == "wilcox") {
    spg <- rownames(cm.norm) %>% split(sample.groups[.])
    test.res <- matrixTests::col_wilcoxon_twosample(cm.norm[spg[[1]],,drop=FALSE], cm.norm[spg[[2]],,drop=FALSE], exact=FALSE)$pvalue
    sel.genes <- test.res %>% setNames(colnames(cm.norm)) %>% sort() %>% names()
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
getTopPCs <- function(cm.norm, n_pcs) {
  if (!is.matrix(cm.norm)) {
    cm.norm <- as.matrix(cm.norm)
  }
  if (!is.numeric(cm.norm)) {
    stop("Input matrix must be numeric.")
  }
  pca_project(cm.norm, n_pcs)
}
