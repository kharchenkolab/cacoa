#' @rawNamespace import(ape, except = where)
#' @importFrom sccore checkPackageInstalled
NULL


## ------------------------------------------------------------
## lmCoda: CoDA-aware LM + permutations
## ------------------------------------------------------------
#' Fits a multivariate linear model to sample-by-cell-type compositions using an
#' isometric log-ratio (ILR) transform, and maps covariate/contrast effects back to
#' cell-type loadings in clr (log-ratio) space.
#'
#' Given counts per sample and cell type, a pseudocount is added to every count (default 0.5 cells in a sample of median size,
#' scaled with the sample's cell number so that it is the same proportion in every sample; half a count is the classical choice
#' that makes the logarithm of a small count nearly unbiased) and sample-wise frequencies are formed. An ILR basis matrix \eqn{\Psi \in \mathbb{R}^{D \times (D-1)}}
#' with orthonormal columns is used to compute ILR coordinates
#' \deqn{Z = \log(X)\Psi,}
#' where \eqn{X} stacks compositions row-wise. A multivariate linear model is then fit in
#' ILR space,
#' \deqn{Z = F B + E,}
#' with design matrix \eqn{F} (or its core subspace for Freedman--Lane),
#' coefficient matrix \eqn{B}, and residuals \eqn{E}.
#'
#' For each model coefficient \eqn{j} (a column of the design matrix) the ILR-space effect is
#' \eqn{\beta^{(j)} = B_{j,\cdot}}, and for a user-defined contrast vector \eqn{\theta} the ILR-space
#' contrast effect is \eqn{\beta^{(\theta)} = \theta^\top B}. Effects are mapped to cell-type
#' (clr) loadings via
#' \deqn{\ell = \Psi \beta,}
#' yielding a length-\eqn{D} loading vector that satisfies the zero-sum constraint and is
#' interpretable only in relative terms. Differences \eqn{\ell_a - \ell_b} correspond to
#' changes in log-ratios of cell-type proportions.
#'
#' Statistical significance is assessed by sample randomization (block permutations or the
#' Freedman--Lane procedure) applied to the fitted ILR-space model. For each coefficient and
#' the specified contrast, the function reports (i) a global composition-level statistic based on
#' the fraction of ILR variance explained by the corresponding fitted component, and (ii)
#' cell-type–specific empirical p-values for clr loadings, with FDR adjustment across cell types.
#'
#' @param cnts Matrix of counts (samples × cell types). Rows are samples; columns are parts.
#' @param model Design specification returned by \code{buildDesignMatrices()}, containing
#'   design matrices (e.g. \code{F} and/or \code{X}) and contrast information aligned to the
#'   coefficient ordering.
#' @param perm.method Permutation scheme used for inference. \code{"block"} permutes in the
#'   full design space; \code{"freedman-lane"} uses the core design space.
#' @param random.baseline fit the grouping term `(1 | g)` of the model as a random baseline per level of `g` (generalized least
#'   squares with a variance ratio estimated per cell type; the residuals, freed of the correlation within levels, are exchanged
#'   among all samples). Default `FALSE`: the grouping term only sets the relabeling, as in every other analysis (exchange within
#'   the levels of `g`; levels that hold one value of the tested variable are exchanged as wholes, see [permutationPlan()]).
#' @param zero.pseudocount Non-negative pseudocount added to every count prior to forming sample-wise frequencies (default 0.5:
#'   the logarithm of a small count plus one half is nearly unbiased for the logarithm of its expectation; in simulations from
#'   the example cohorts the loadings of rare cell types were least biased at 0.25-0.5 and overstated by 0.2-0.4 log2 at 0.1).
#'   The result records how much each loading moves when the pseudocount is halved or doubled (`contrast$loadings$pseudocount.shift`,
#'   in standard errors; see [codaLoadingsTable()])
#' @param scale.pseudocount scale the pseudocount with the sample's cell number relative to the median sample (default TRUE: the
#'   same proportion is added in every sample)
#' @param basis.type ILR basis type passed to \code{coda.base::ilr_basis()}.
#' @param ref.method reference cell types for the presentation of the loadings: `"block"` (default) the largest set of
#'   cell types whose proportions do not change relative to each other ([referenceBlock()]: pairwise permutation tests of the
#'   loading differences, merged while p >= `ref.block.thresh`), or `"pvalue"` (the previous rule: the cell types with the
#'   weakest evidence of a change relative to the geometric mean, `ref.p.thresh`, `ref.min.size`, `ref.max.size`)
#' @param ref.block.thresh p-value below which two groups of cell types are kept apart in the block search (default 0.05)
#' @param ref.size.frac blocks at least this fraction of the largest block's size compete for the reference, the one closest to the
#'   center of the composition winning (default 0.5; [referenceBlock()])
#' @param ref.p.thresh P-value threshold used to select a "least-affected" reference set for
#'   centering loadings (`ref.method = "pvalue"`).
#' @param ci.level level of the intervals of the predicted endpoint compositions (default 0.95; `predicted$baseline.lo` etc.)
#' @param ref.min.size Minimum size of the reference set (`ref.method = "pvalue"`).
#' @param ref.max.size Maximum size of the reference set (`ref.method = "pvalue"`).
#' @param ... Additional arguments forwarded to the permutation/LM routine
#'   (e.g. \code{n.permutations}, na.mode, robust, alternative, etc.).
#' 
#' @return A list containing ILR coordinates and basis, fitted model/permutation output, and
#'   coefficient- and contrast-level results including global variance-explained statistics,
#'   clr-space cell-type loadings, empirical p-values/FDR (relative to the geometric mean of all cell types, `pval` /
#'   `padj`, and relative to the reference cell types, `pval.ref` / `padj.ref`), and (if available) predicted endpoint
#'   compositions for the specified contrast; `reference` (the reference cell types, [referenceBlock()]) and `tree`
#'   (the contrast tree, [contrastTree()]).
#' @keywords internal
lmCoda <- function(cnts, model, perm.method = c("block", "freedman-lane"), random.baseline = FALSE, zero.pseudocount = 0.5, scale.pseudocount = TRUE,
                   basis.type = c("default"), ref.method = c("block", "pvalue"), ref.block.thresh = 0.05, ref.size.frac = 0.5, ref.p.thresh = 0.3, ref.min.size = 1,
                   ref.max.size = 3, ci.level = 0.95, ...) {
  
  basis.type  <- match.arg(basis.type); ref.method <- match.arg(ref.method)
  perm.method <- match.arg(perm.method)
  # a grouping term `(1 | g)` restricts the relabeling to within the levels of g (the model's block.vars), as in every other analysis;
  # with `random.baseline = TRUE` it is fitted as a random baseline per level instead
  use.rnd <- isTRUE(random.baseline)
  if (use.rnd && !length(model$random.vars)) stop("random.baseline = TRUE needs a grouping term `(1 | g)` in the model formula")
  
  cnts <- as.matrix(cnts)
  if (is.null(rownames(cnts)))
    rownames(cnts) <- paste0("sample_", seq_len(nrow(cnts)))
  
  ## 1) ILR transform
  ilr.res <- computeILRMatrix(cnts = cnts, zero.pseudocount = zero.pseudocount, basis.type = basis.type, scale.pseudocount = scale.pseudocount)
  B.ilr <- ilr.res$ilr     # n_samples × K
  psi   <- ilr.res$psi     # D × K
  freqs <- ilr.res$freqs # 
  
  K <- ncol(B.ilr)
  D <- nrow(psi)
  cell.types <- rownames(psi)
  
  ## 1b) Design space: the full design F (level means for the tested factor); the contrast scheme only affects the permutations
  design <- model$F
  contrast.vec <- model$contrast.F
  endpoints <- model$contrast_endpoints_F
  if (is.null(design)) stop("Design matrix F is missing in model.")
  if (nrow(design) != nrow(B.ilr))
    stop("Row counts of design matrix and ILR matrix do not match.")
  
  ## Total variance in ILR space (common denominator)
  B.centered <- scale(B.ilr, center = TRUE, scale = FALSE)
  total.ss   <- sum(B.centered^2)
  
  ## 2) Linear model + permutations in ILR space
  fit <- performLMPermutations(model, B.ilr, perm.method = perm.method,
                               return.sampled.stats = TRUE,
                               return.sampled.fits  = FALSE, ...)
  
  # coefficients in F space (level means of the tested factor in ILR coordinates) for predictions and adjusted ILR
  coef.mat <- qr.coef(qr(design), B.ilr); coef.mat[is.na(coef.mat)] <- 0
  dimnames(coef.mat) <- list(colnames(design), colnames(B.ilr))
  fit$coef <- coef.mat
  stat.obs     <- as.numeric(fit$effect)        # length K: contrast effect in ILR space (the loadings are built from effects, not from the test statistic)
  stats.perm   <- as.matrix(fit$effects.perm)   # n_perm × K: the same under relabeling
  
  n.coef     <- nrow(coef.mat)
  n.perm     <- nrow(stats.perm)
  
  ## --------------------------------------------------------
  ## 3) USER-DEFINED CONTRAST
  ##    - uses contrast in the *same space* as coef_mat
  ## --------------------------------------------------------
  if (is.null(contrast.vec))
    stop("Contrast vector for perm.method='", perm.method, "' is missing in model.")
  
  beta.contrast.obs  <- stat.obs              # length K (ILR-space contrast effect)
  beta.contrast.perm <- stats.perm            # n_perm × K
  
  # Cell-type loadings for contrast
  ell.contrast.obs  <- drop(psi %*% beta.contrast.obs)          # D
  ell.contrast.perm <- psi %*% t(beta.contrast.perm)            # D × n_perm
  rownames(ell.contrast.perm) <- cell.types
  
  # Per-cell p-values (reference-free, two-sided)
  ex <- isTRUE(fit$exhaustive)   # all distinct relabelings enumerated: the observed labeling is one of them (relabelingPValue())
  p.contrast.cell <- relabelingPValue(rowSums(abs(ell.contrast.perm) >= abs(ell.contrast.obs) - 1e-12), n.perm, ex)
  names(p.contrast.cell) <- cell.types
  # With `statistic = "t"` (the default) the cell types are tested by the studentized loading, recomputed on every relabeled data set
  # (contrastPermutationTest(); the same relabelings as the fit). The unstudentized loading is not a valid statistic when the model
  # has covariates. Robust fits keep the fitter's own relabeled effects.
  p.statistic <- "coef"
  if (identical((list(...)$statistic %||% "t")[1], "t") && identical((list(...)$robust.method %||% "none")[1], "none") && !is.null(fit$P) && !use.rnd) {
    ct <- tryCatch({ Yc <- B.ilr %*% t(psi); colnames(Yc) <- cell.types; contrastPermutationTest(Yc, design, contrast.vec, fit$P, scheme = perm.method) }, error = function(e) NULL)
    if (!is.null(ct) && all(is.finite(ct$t)) && max(abs(ct$est - ell.contrast.obs)) < 1e-6 * max(1, max(abs(ell.contrast.obs)))) {
      # the null loadings shown as a band and used by the reference block and the contrast tree: the relabeled t statistics at the
      # observed standard errors (what the loading would be without a change, at the precision of the observed fit)
      ell.contrast.perm <- ct$t.perm * ct$se; bad <- !is.finite(ell.contrast.perm); ell.contrast.perm[bad] <- ct$est.perm[bad]; rownames(ell.contrast.perm) <- cell.types; n.perm <- ncol(ct$t.perm)
      p.contrast.cell <- setNames(relabelingPValue(rowSums(abs(ct$t.perm) >= abs(ct$t) - 1e-12), n.perm, ex), cell.types); p.statistic <- "t"
    }
  }

  # Grouping term `(1 | g)` of the model: a random intercept per level of g for every cell type's log ratio. The loadings are the
  # generalized least-squares estimates (variance ratio per cell type), tested by their t statistics against relabeled data sets on
  # which the whole fit is repeated (randomInterceptContrast()); the global test is the sum of the squared t statistics.
  rnd <- NULL
  if (use.rnd) {
    if (length(model$random.vars) > 1) stop("one grouping term is supported in the composition model")
    gv <- as.character(model$meta[rownames(design) %||% rownames(cnts), model$random.vars]); if (anyNA(gv)) stop("missing values in the grouping variable")
    n.perm <- as.integer(list(...)$n.permutations %||% max(n.perm, 499L))   # free relabelings of the transformed residuals: the requested number, not the number of distinct within-level relabelings
    Pm <- withSeed(list(...)$seed %||% 1L, vapply(seq_len(n.perm), function(b) sample.int(nrow(design)), integer(nrow(design))))
    Yc <- B.ilr %*% t(psi); colnames(Yc) <- cell.types
    rnd <- randomInterceptContrast(Yc, design, contrast.vec, gv, Pm)
    ell.contrast.obs <- rnd$est; ell.contrast.perm <- rnd$t.perm * rnd$se; rownames(ell.contrast.perm) <- cell.types
    p.contrast.cell <- setNames((rowSums(abs(rnd$t.perm) >= abs(rnd$t) - 1e-12) + 1) / (n.perm + 1), cell.types)
    p.statistic <- "t"
    rnd$global <- list(stat = sum(rnd$t^2), stat_perm = colSums(rnd$t.perm^2)); rnd$global$p <- (sum(rnd$global$stat_perm >= rnd$global$stat) + 1) / (n.perm + 1)
    rnd$group <- model$random.vars; rnd$icc <- rnd$theta / (1 + rnd$theta)
  }
  padj.contrast.cell <- p.adjust(p.contrast.cell, method = "fdr")
  
  ## Global stat for contrast = variance explained in ILR space
  # NOTE: contrast.vec must live in same space (F or X) as 'design'
  z      <- as.numeric(design %*% contrast.vec)  # n_samples
  sum.z2 <- sum(z^2)
  scale.contrast <- sum.z2 / total.ss           # scaling factor
  
  # variance explained (observed + permuted)
  T.contrast.obs  <- scale.contrast * sum(beta.contrast.obs^2)
  T.contrast.perm <- scale.contrast * rowSums(beta.contrast.perm^2)
  
  # p-value (same ordering, now in var-expl units)
  p.contrast.global <- relabelingPValue(sum(T.contrast.perm >= T.contrast.obs * (1 - 1e-12)), length(T.contrast.perm), ex)
  var.expl <- T.contrast.obs; global.statistic <- "variance explained in ILR space"
  if (!is.null(rnd)) { T.contrast.obs <- rnd$global$stat; T.contrast.perm <- rnd$global$stat_perm; p.contrast.global <- rnd$global$p; global.statistic <- "sum of squared t statistics of the cell types" }
  else if (p.statistic == "t") {   # the global test from the studentized loadings as well: the unstudentized variance explained is not a valid statistic with covariates
    T.contrast.obs <- sum(ct$t^2); T.contrast.perm <- colSums(ct$t.perm^2); T.contrast.perm[!is.finite(T.contrast.perm)] <- Inf
    p.contrast.global <- relabelingPValue(sum(T.contrast.perm >= T.contrast.obs - 1e-12), length(T.contrast.perm), ex); global.statistic <- "sum of squared t statistics of the cell types"
  }
  
  # Reference cell types for the presentation of the loadings: the largest block of cell types unchanged relative to each
  # other (the default), or the cell types with the weakest evidence of a change relative to the geometric mean
  ref <- if (ref.method == "block") referenceBlock(ell.contrast.obs, ell.contrast.perm, p.thresh = ref.block.thresh, size.frac = ref.size.frac)
         else c(pickReferenceCluster(loadings = ell.contrast.obs, pval = p.contrast.cell, p.thresh = ref.p.thresh, min.size = ref.min.size, max.size = ref.max.size), list(method = "pvalue"))
  ell.contrast.ref <- ell.contrast.obs - ref$mean
  # model-based standard errors of the loadings (least squares on the log ratios): se_j^2 = c'(F'F)^-1 c x residual variance of the
  # cell type's clr value; `se.ref` for the loading relative to the reference cell types (conditional on that data-chosen set)
  lse <- tryCatch({
    qd <- qr(design); df.res <- nrow(design) - qd$rank
    if (df.res < 1) stop("no residual degrees of freedom")
    XtXi <- chol2inv(qr.R(qd)); pv <- qd$pivot[seq_len(qd$rank)]; cvp <- contrast.vec[pv]
    if (any(abs(contrast.vec[-pv]) > 1e-12) && qd$rank < ncol(design)) stop("contrast involves an aliased column")
    cc <- sum(cvp * (XtXi[seq_len(qd$rank), seq_len(qd$rank), drop = FALSE] %*% cvp))
    Yc <- B.ilr %*% t(psi); colnames(Yc) <- cell.types; rs <- qr.resid(qd, Yc)
    rr <- if (length(ref$celltypes)) rs - rowMeans(rs[, ref$celltypes, drop = FALSE]) else rs
    list(se = sqrt(cc * colSums(rs^2) / df.res), se.ref = sqrt(cc * colSums(rr^2) / df.res), df = df.res)
  }, error = function(e) NULL)
  # dependence on the pseudocount: the largest change of each loading when the pseudocount is halved or doubled, in standard errors
  pc.shift <- tryCatch({
    lo <- function(pc) { z <- computeILRMatrix(cnts, zero.pseudocount = pc, basis.type = basis.type, scale.pseudocount = scale.pseudocount)$ilr
      b <- qr.coef(qr(design), z); b[is.na(b)] <- 0; drop(psi %*% drop(contrast.vec %*% b)) }
    l0 <- lo(zero.pseudocount); d <- pmax(abs(lo(zero.pseudocount / 2) - l0), abs(lo(zero.pseudocount * 2) - l0)); setNames(d / lse$se, cell.types)
  }, error = function(e) NULL)
  if (!is.null(rnd) && !is.null(lse)) {   # standard errors of the generalized least-squares loadings; relative to the reference, from the fit of that log ratio
    lse$se <- rnd$se; lse$df <- rnd$df
    if (length(ref$celltypes)) lse$se.ref <- tryCatch({ Yr <- (B.ilr %*% t(psi)); colnames(Yr) <- cell.types; Yr <- Yr - rowMeans(Yr[, ref$celltypes, drop = FALSE])
      gv <- as.character(model$meta[rownames(design) %||% rownames(cnts), model$random.vars])
      setNames(vapply(seq_len(ncol(Yr)), function(j) { th <- remlVarianceRatio(Yr[, j], design, gv); W <- groupWhitener(gv, th); Fw <- W(design); q <- qr(Fw); XtXi <- chol2inv(qr.R(q))
        sqrt(sum(contrast.vec * (XtXi %*% contrast.vec)) * sum(qr.resid(q, drop(W(Yr[, j])))^2) / (nrow(Fw) - q$rank)) }, numeric(1)), cell.types) }, error = function(e) lse$se)
    names(lse$se) <- cell.types
  }
  cnt.summary <- data.frame(celltype = cell.types, median.cells = apply(cnts, 2, stats::median), n.zero = colSums(cnts == 0), n.samples = nrow(cnts), row.names = NULL, stringsAsFactors = FALSE)
  ref.rel <- referenceRelativeP(ell.contrast.obs, ell.contrast.perm, ref$celltypes)   # each cell type against the reference, same relabelings
  
  ## --------------------------------------------------------
  ## 3b) PREDICTED BASELINE vs TARGET COMPOSITIONS (if endpoints present)
  ##      endpoints must be in the *same coefficient space* as coef.mat
  ## --------------------------------------------------------
  predicted <- NULL
  if (!is.null(endpoints)) {
    if (!all(c("num", "den") %in% names(endpoints))) {
      warning("contrast endpoints present but missing 'num'/'den'; skipping predicted compositions.")
    } else {
      x.den <- endpoints$den
      x.num <- endpoints$num
      
      if (length(x.den) != n.coef || length(x.num) != n.coef) {
        stop("Length of contrast endpoints (", length(x.den),
             ") does not match number of coefficients (", n.coef, ").")
      }
      
      # baseline/target ILR: 1×K = (1×n_coef) %*% (n_coef×K)
      b.den <- as.numeric(x.den %*% coef.mat)
      b.num <- as.numeric(x.num %*% coef.mat)
      
      ilr.den <- matrix(b.den, nrow = 1)
      ilr.num <- matrix(b.num, nrow = 1)
      
      # back-transform to compositions
      f.den <- ilrInverseFromBasis(ilr.den, psi)[1, ]
      f.num <- ilrInverseFromBasis(ilr.num, psi)[1, ]
      
      # normalize & name (defensive, though inverse already normalizes)
      f.den <- f.den / sum(f.den)
      f.num <- f.num / sum(f.num)
      names(f.den) <- cell.types
      names(f.num) <- cell.types
      
      # choose human-readable labels from model, with safe defaults
      endpoint.labels <- list(baseline = "baseline", target   = "target")
      if (!is.null(model$contrast_endpoint_labels)) {
        if (!is.null(model$contrast_endpoint_labels$baseline))
          endpoint.labels$baseline <- model$contrast_endpoint_labels$baseline
        if (!is.null(model$contrast_endpoint_labels$target))
          endpoint.labels$target <- model$contrast_endpoint_labels$target
      }
      
      predicted <- list(baseline = f.den, target = f.num, delta = f.num - f.den,
                        labels = endpoint.labels)
      # 95 % intervals of the fitted endpoint compositions, from the ILR fit (see predictedCompositionCI())
      ci <- tryCatch(predictedCompositionCI(design, coef.mat, B.ilr, psi, list(baseline = x.den, target = x.num), level = ci.level), error = function(e) NULL)
      if (!is.null(ci)) predicted <- c(predicted, list(baseline.lo = ci$baseline$lo, baseline.hi = ci$baseline$hi, target.lo = ci$target$lo, target.hi = ci$target$hi,
                                                       ci.level = ci.level, ci.method = "ILR fit, multivariate t"))
    }
  }
  
  ## --------------------------------------------------------
  ## 5) Final result structure
  ## --------------------------------------------------------
  out <- list(ilr = B.ilr, psi = psi, freqs = freqs,
       fit = fit, perm.method = perm.method, # Original LM/permutation output
       contrast = list(beta_ilr = beta.contrast.obs, # ILR effect vector
                       global = list(stat = T.contrast.obs, # the global statistic (see `statistic`)
                                     stat_perm = T.contrast.perm, # the same under relabeling
                                     p = p.contrast.global, statistic = global.statistic, var.explained = var.expl),
                       loadings = list(obs = ell.contrast.obs, # cell-type loadings (ref-free)
                                       ref_centered = ell.contrast.ref,
                                       se = lse$se, se.ref = lse$se.ref, df = lse$df,   # model-based standard errors (natural log, as the loadings)
                                       pseudocount.shift = pc.shift, counts = cnt.summary,
                                       perm = ell.contrast.perm), # D x n_perm null loadings
                       per_cell = list(pval = p.contrast.cell, padj = padj.contrast.cell, pval.ref = ref.rel$pval, padj.ref = ref.rel$padj),   # relative to the geometric mean; relative to the reference cell types
                       predicted = predicted,# baseline / target / delta compositions (if available)
                       label = if (!is.null(model$contrast_label)) model$contrast_label else NULL),
       coef = coef.mat,              # F-space coefficients (level means in ILR space)
       random = if (!is.null(rnd)) list(group = rnd$group, theta = rnd$theta, icc = rnd$icc, global = rnd$global) else NULL,   # random intercept: variance ratio and share of variance between groups, per cell type
       reference = ref, model = model,
       settings = list(zero.pseudocount = zero.pseudocount, scale.pseudocount = scale.pseudocount, per.celltype.statistic = p.statistic, basis.type = basis.type, perm.method = perm.method, ref.method = ref.method,
                       ref.block.thresh = ref.block.thresh, ref.size.frac = ref.size.frac, ref.p.thresh = ref.p.thresh, ref.min.size = ref.min.size,
                       ref.max.size = ref.max.size, n.permutations = n.perm))
  out$cnts <- cnts
  out$tree <- tryCatch(contrastTree(out), error = function(e) NULL)   # the contrast tree ([contrastTree()]) from the same loadings and null
  out
}

## ILR transform helper
computeILRMatrix <- function(cnts, zero.pseudocount = 0.5, basis.type = c("default"), scale.pseudocount = TRUE) {
  basis.type <- match.arg(basis.type)
  cnts <- as.matrix(cnts)
  
  if (any(rowSums(cnts) <= 0))
    stop("Some samples have zero total count; cannot compute frequencies.")
  
  # the pseudocount is added to every count; scaled with the sample's cell number it is the same proportion in every sample
  n.s <- rowSums(cnts)
  cnts <- cnts + if (scale.pseudocount) zero.pseudocount * n.s / stats::median(n.s) else zero.pseudocount
  
  freqs <- sweep(cnts, 1, rowSums(cnts), "/")
  
  # ILR basis: rows = cell types, cols = ILR dimensions
  psi <- coda.base::ilr_basis(ncol(freqs), type = basis.type)
  rownames(psi) <- colnames(freqs)
  
  # ILR coordinates: samples × ILR dims
  ilr <- log(freqs) %*% psi
  
  list(
    ilr   = ilr,    # n_samples × (D-1)
    psi   = psi,    # D × (D-1)
    freqs = freqs   # n_samples × D
  )
}

#' Reference cluster helper: pick "least affected" cell types
#' Chooses a small set of cell types with minimal evidence of association and uses it to
#' compute reference-centered (presentation-only) loadings. The reference set is intended
#' to anchor interpretation of clr-space loadings by selecting cell types that appear
#' approximately unchanged under the tested contrast.
#' @param loadings Named numeric vector of clr-space cell-type loadings \eqn{\ell}.
#' @param pval Named numeric vector of empirical p-values for each cell type.
#' @param p.thresh P-value threshold used to define candidate reference cell types.
#' @param min.size Minimum number of cell types to include in the reference set.
#' @param max.size Maximum number of cell types to include in the reference set.
#'
#' @return A list with indices and names of the selected reference cell types and the mean
#'   loading over the reference set (useful for centering).
#' @keywords internal
pickReferenceCluster <- function(loadings, pval, p.thresh = 0.3, min.size = 1,
                                 max.size = 3) {
  stopifnot(length(loadings) == length(pval))
  ct <- names(loadings)
  # candidates: high p-value (weak evidence)
  idx <- which(pval > p.thresh)
  
  # if too few candidates, fall back to smallest absolute loadings
  if (length(idx) < min.size) {
    idx <- order(abs(loadings))[seq_len(min.size)]
  }
  # if too many candidates, keep those closest to zero
  if (length(idx) > max.size) {
    idx <- idx[order(abs(loadings[idx]))[seq_len(max.size)]]
  }

  list(idx = idx, celltypes = ct[idx], mean = mean(loadings[idx]))
}

## Inverse ILR given coordinates (rows) and basis psi (D × (D-1))
ilrInverseFromBasis <- function(z, psi) {
  z   <- as.matrix(z)          # n_samples × (D-1)
  logf <- z %*% t(psi)         # n_samples × D
  f    <- exp(logf)
  f / rowSums(f)
}

# Intervals of the fitted compositions at the contrast endpoints. The endpoint's ILR mean x_e B has covariance v_e Sigma, with
# v_e = x_e (F'F)^- x_e' and Sigma the residual covariance of the ILR fit (K x K, residual df = n - rank F); draws from the
# multivariate t around the fitted mean are mapped to proportions, and the per-cell-type quantiles are the interval.
# `endpoints`: named list of design rows (length ncol(design)); returns, per endpoint, list(lo, hi) named by cell type.
predictedCompositionCI <- function(design, coef, Y, psi, endpoints, level = 0.95, n.draws = 4000, seed = 1) {
  qd <- qr(design); r <- qd$rank; n <- nrow(design); df <- n - r
  if (df < 2) return(NULL)
  R <- Y - design %*% coef; Sigma <- crossprod(R) / df; K <- ncol(Y)
  piv <- qd$pivot[seq_len(r)]; XtX.inv <- matrix(0, ncol(design), ncol(design))
  XtX.inv[piv, piv] <- chol2inv(qr.R(qd)[seq_len(r), seq_len(r), drop = FALSE])
  L <- tryCatch(t(chol(Sigma)), error = function(e) { e <- eigen(Sigma, symmetric = TRUE); e$vectors %*% diag(sqrt(pmax(e$values, 0)), K) })
  if (exists(".Random.seed", envir = .GlobalEnv)) { old <- get(".Random.seed", envir = .GlobalEnv); on.exit(assign(".Random.seed", old, envir = .GlobalEnv)) }
  set.seed(seed)
  lapply(endpoints, function(x) {
    x <- as.numeric(x); v <- as.numeric(t(x) %*% XtX.inv %*% x); z0 <- as.numeric(x %*% coef)
    U <- (matrix(stats::rnorm(n.draws * K), n.draws, K) %*% t(L)) * sqrt(v)       # N(0, v Sigma)
    W <- sqrt(stats::rchisq(n.draws, df) / df)
    P <- ilrInverseFromBasis(sweep(U / W, 2, z0, "+"), psi)                        # multivariate t around the fitted endpoint, as proportions
    q <- apply(P, 2, stats::quantile, probs = c((1 - level) / 2, 1 - (1 - level) / 2), names = FALSE)
    list(lo = stats::setNames(q[1, ], rownames(psi)), hi = stats::setNames(q[2, ], rownames(psi)))
  })
}

#' Get loadings by different methods (old two-group comparison)
#'
#' @param cnts Counts of cell typer in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param method Method to get loadings
#' @param ref.cell.type Reference cell type
#' @return Updated data frame with Z scores
#' @keywords internal
getLoadings <- function(cnts, groups, method = c('lda', 'svm', 'cda', 'cda.std'), ref.cell.type = NULL) {
  # Checks
  method <- match.arg(method)
  if (method == "svm") checkPackageInstalled("e1071", cran=TRUE)
  if (method == "cda") checkPackageInstalled("candisc", cran=TRUE)
  if(!is.null(ref.cell.type) && !(ref.cell.type %in% colnames(cnts)))
     stop(paste('Reference cell type', ref.cell.type, 'is not correct. Correct cell types are:',
                paste0(colnames(cnts), collapse = ', ') ))
  
  #Get freqs
  cnts[cnts == 0] <- 0.1
  freqs <- cnts/rowSums(cnts)

  # Get ilr
  psi <- coda.base::ilr_basis(ncol(freqs), type = "default") %>% 
    `rownames<-`(cnts %>% colnames())
  b <- log(freqs) %*% psi

  # PCA
  pca.res <- apply(b, 2, function(y) y - mean(y)) %>% 
    prcomp()
  pca.loadings <- psi %*% pca.res$rotation

  psi <- psi %*% pca.res$rotation
  b <- (log(freqs) %*% psi) %>% 
    apply(2, scale)

  idx.na = which(colSums(is.nan(b)) == 0)
  b <- b[,idx.na]
  psi <- psi[,idx.na]

  # Create dataframe
  b.df <- data.frame(b) %>% 
    mutate(groups = 1*groups)

  # Optimization
  if(method == 'lda') {
    if(is.null(ref.cell.type)) {  # pure linear regression which is proportional LDA when number of classes = 2
      res.lm <- lm(groups ~ ., data = b.df)
      w <- res.lm$coefficients[-1]
      w[is.na(w)] <- 0

      res.sum = summary(res.lm)
    } else {  # Quadratic optimization with the linear condition corresponding to the reference level
      X <- b
      Y <- groups * 1
      X2 <- cbind(X, 1)
      Rinv <- solve(chol(t(X2) %*% X2))
      dvec <- t(Y) %*% X2
      Amat <- t(psi[ref.cell.type,,drop=FALSE])
      Amat <- rbind(Amat, 0)
      res.qp <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)

      n.qp <- ncol(X2)
      w <- res.qp$solution[-n.qp]
    }
  } else if(method == 'svm') {
    b.model <- e1071::svm(groups~., data = b.df, kernel = "linear", scale = FALSE)
    # Get hyperplane parameters
    w <- t(b.model$SV) %*% b.model$coefs  # slope
  } else if(method == 'cda') {  # Canonical discriminant analysis
    model <- lm(b ~ groups)
    cda <- candisc::candisc(model, ndim=1)
    w <- cda$coeffs.raw
  } else if(method == 'cda.std') {
    pca.res <- apply(b, 2, function(y) y - mean(y)) %>% 
      prcomp()
    
    w <- pca.res$rotation %*% t(cor(groups, pca.res$x))
  }

  w <- w / sqrt(sum(w ^ 2))
  scores <- b %*% w
  
  if(mean(scores[groups]) < mean(scores[!groups])) w = -w

  loadings <- (psi %*% w) %>% 
    `names<-`(cnts %>% colnames())

  return(loadings)
}

#' This function produces the resampling dataset (old two-group comparison)
#'
#' @param cnts Counts of cell types in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param n.iter Number of permutations
#' @param remain.groups TRUE - if the labels of groups remain, FALSE - if to produce the null distribution
#' @param replace.samples TRUE - is bootstrap on samples, FALAE - if to remain samples
#' @param seed random seed
#' @return Updated data frame with Z scores
#' @keywords internal
produceResampling <- function(cnts, groups, n.perm = 1000, seed = 239) {
  cnts.perm <- list()
  groups.perm <- list()
  set.seed(seed)

  for (i in 1:n.perm) {
    samples.tmp <- rownames(cnts) %>% split(groups[.]) %>%
      lapply(sample, replace=TRUE) %>% do.call(c, .)
    groups.tmp <- groups[samples.tmp]

    # sample cells
    cnts.resampling <- apply(cnts[samples.tmp,], 1, function(v) {
      if (sum(v) == 0) return(setNames(rep(0, length(v)), names(v)))

      sample(names(v), size=sum(v), replace=TRUE, prob = v / sum(v)) %>%
        {c(names(v), .)} %>% table() %>% {. - 1}
    }) %>% t()

    cnts.perm[[length(cnts.perm) + 1]] <- cnts.resampling
    groups.perm[[length(groups.perm) + 1]] <- groups.tmp
  }

  return(list(cnts = cnts.perm, groups = groups.perm))
}

#' This function runs the coda analysis (old two-group comparison)
#' @keywords internal 
runCoda <- function(cnts, groups, n.seed=239, n.boot=1000, ref.cell.type=NULL, null.distr=FALSE, method="lda", n.cores=1, verbose=TRUE) {
  # Create datasets as
  samples.init <- produceResampling(cnts = cnts, groups = groups, n.perm = n.boot, seed = n.seed)
  loadings <- do.call(cbind, lapply(1:length(samples.init$cnts), function(ib) {
    getLoadings(samples.init$cnts[[ib]], samples.init$groups[[ib]], method=method)
  })) %>% 
    data.frame()

  # Calculate p-values of confidence interval by bootstrap
  tmp <- referenceSet(cnts, groups, p.thresh = 0.1)
  cell.list <- tmp$cell.list

  mean.list <- c()  # mean value of loadings in a list
  for (i.list in 1:length(cell.list)) {
    loadings.list <- loadings[cell.list[[i.list]],]
    mean.list <- c(mean.list, mean(unlist(loadings.list)))

    if (length(cell.list[[i.list]]) == 1) {  # if a list contains only one sample - it cannot be a reference group
      mean.list[i.list] <- 10
    }
  }

  # Define a cluster with reference cell type
  if (is.null(ref.cell.type)) {
    id.ref.cluster <- which(abs(mean.list) == min(abs(mean.list)))  # <- working version
  } else {
    id.ref.cluster <- -1
    for (i.list in 1:length(cell.list)) {
      if (ref.cell.type[1] %in% cell.list[[i.list]]) {
        id.ref.cluster <- i.list
      }
    }
    if (id.ref.cluster == -1) stop('Wrong name of reference cell type')
  }
  
  ref.cell.type <- cell.list[[id.ref.cluster]]

  # sorting of cell types and reference level
  cell.types.order <- c()
  for (i.list in order(-abs(mean.list - mean.list[id.ref.cluster]))) {
    cell.types.tmp <- cell.list[[i.list]] %>%
      .[order(-abs(rowMeans(loadings[.,, drop=FALSE]) - mean.list[id.ref.cluster]))]
    cell.types.order <- c(cell.types.order, cell.types.tmp)
  }
  
  loadings <- loadings[cell.types.order,]
  ref.load.level <- loadings[ref.cell.type,] %>% 
    unlist() %>% 
    mean()

  if (!null.distr) {
    tmp <- rowSums(ref.load.level > loadings)
    pval <- (pmin(tmp, ncol(loadings) - tmp) + 1) / (ncol(loadings) + 1) * 2
    names(pval) <- rownames(loadings)
  } else {
    loadings.perm <- do.call(cbind, plapply(1:n.boot, function(ib) {
      groups.perm <- sample(groups)
      names(groups.perm) <- names(groups)
      samples.perm <- produceResampling(cnts=cnts, groups=groups.perm, n.perm = 100, seed = n.seed + ib)

      do.call(cbind, lapply(1:length(samples.perm$cnts), function(ib) {
        getLoadings(samples.perm$cnts[[ib]], samples.perm$groups[[ib]])
      })) %>% rowMeans()
    }, n.cores=n.cores, progress=verbose, mc.preschedule=TRUE))

    loadings.stat <- rowMeans(loadings) - ref.load.level
    pval <- sapply(names(loadings.stat), function(s) {
      tmp <- sum(loadings.stat[s] > loadings.perm[s,])
      p <- (pmin(tmp, ncol(loadings) - tmp) + 1) / (ncol(loadings) + 1) * 2
      return(p)
    })
  }

  padj <- p.adjust(pval, method = "fdr")
  loadings <- loadings[names(pval),]
  return(list(padj=padj,
              loadings=loadings,
              pval=pval,
              ref.load.level=ref.load.level,
              ref.cell.type=ref.cell.type,
              cell.list=cell.list))
}

#' Plot compositional loadings of a contrast
#'
#' One row per cell type, in the style of [plotShiftMagnitudes()]: the point is the loading of the contrast on the cell type (log2;
#' the change of the cell type's log proportion relative to the geometric mean of all cell types), centered on the reference cell
#' types chosen by `lmCoda()` ([referenceBlock()], marked "(ref)") or not (`center`), in the cell-type color, filled when the cell
#' type is significant by the primary test (relative to the geometric mean; `significance`). The grey band is the 95 % range of
#' the tested quantity over the permutations on the same axis: what the loading is without a composition change.
#'
#' @param res result of `lmCoda()` (as stored by `cao$estimateCellLoadings()`)
#' @param ref.level,target.level labels for the axis annotation
#' @param alpha significance threshold on the adjusted per-cell-type p-values (default 0.05)
#' @param palette named colours per cell type (optional)
#' @param ordering `"loading"` (default: by the loading) or `"pvalue"` (significant cell types first, by p-value)
#' @param show.pvals print the adjusted p-values next to the points (default TRUE)
#' @param show.null show the permutation band (default TRUE)
#' @param plot.theme ggplot2 theme
#' @param legend.position legend position (default "bottom"; "none" hides it)
#' @param significance the test behind the filled points, the p-values and the band: `"clr"` (default: the primary per-cell-type
#'   test, each cell type relative to the geometric mean of all cell types, which does not depend on the reference) or
#'   `"reference"` (relative to the reference cell types; conditional on that data-chosen set)
#' @param center `"reference"` (default): the loadings centred on the reference cell types (zero = no change relative to them; the
#'   geometric mean of all cell types is marked by a dashed line), or `"none"`: the loadings as they are (zero = the geometric
#'   mean; the reference block's mean is marked). The band is the null of the test behind the points placed on the same axis
#' @param title plot title (optional)
#' @param interval draw the 95 % confidence interval of each loading (default TRUE): the model-based interval of the least-squares fit
#'   on the log ratios (t quantile with the residual degrees of freedom; `interval.level`). With `center = "reference"` it is the
#'   interval of the loading relative to the reference cell types, given that set. It is not the test: the filled points and the
#'   p-values come from the permutations
#' @param interval.level level of the interval (default 0.95)
#' @param min.effect optional fold change (e.g. 1.5) marked by dotted lines at +- log2(`min.effect`): a cell type whose interval lies
#'   beyond a line changes by at least that factor relative to the zero of the axis (see [codaLoadingsTable()])
#' @return ggplot object
#' @export
plotCodaLoadings <- function(res, ref.level = "reference", target.level = "target", alpha = 0.05, palette = NULL,
                             ordering = c("loading", "pvalue"), show.pvals = TRUE, show.null = TRUE, plot.theme = ggplot2::theme_bw(),
                             legend.position = "bottom", significance = c("clr", "reference"), center = c("reference", "none"), title = NULL,
                             interval = TRUE, interval.level = 0.95, min.effect = NULL) {
  ordering <- match.arg(ordering); significance <- match.arg(significance); center <- match.arg(center)
  if (is.null(res$contrast)) stop("not a contrast result of lmCoda()")
  l <- res$contrast$loadings$obs; cts <- names(l); ref <- res$reference$celltypes; ref <- ref[ref %in% cts]
  off <- if (length(ref)) mean(l[ref]) else 0                      # the reference block's mean loading (natural log)
  use.ref <- significance == "reference"
  if (use.ref && is.null(res$contrast$per_cell$padj.ref)) stop("no reference-relative p-values in the result")
  pc <- res$contrast$per_cell; pval <- (if (use.ref) pc$pval.ref else pc$pval)[cts]; padj <- (if (use.ref) pc$padj.ref else pc$padj)[cts]
  shift <- if (center == "reference") -off else 0                 # displayed value = (loading + shift) / log(2)
  d <- data.frame(celltype = cts, value = unname(l + shift) / log(2), pval = unname(pval), padj = unname(padj), nlo = NA_real_, nhi = NA_real_,
                  sig = is.finite(padj) & unname(padj) <= alpha, reference = cts %in% ref, stringsAsFactors = FALSE)
  d$lo <- d$hi <- NA_real_
  se <- if (center == "reference") res$contrast$loadings$se.ref else res$contrast$loadings$se
  if (interval && !is.null(se) && is.finite(res$contrast$loadings$df %||% NA)) {
    tq <- stats::qt(1 - (1 - interval.level) / 2, res$contrast$loadings$df); d$lo <- d$value - tq * unname(se[cts]) / log(2); d$hi <- d$value + tq * unname(se[cts]) / log(2)
  }
  perm <- res$contrast$loadings$perm
  if (show.null && !is.null(perm)) {
    # the null of the test behind the points (the loading, or its difference from the reference block), placed on the axis of the plot:
    # with `significance = "clr"` the band is centred where a cell type that followed the geometric mean would sit
    pr <- perm[cts, , drop = FALSE]; if (use.ref) pr <- sweep(pr, 2, colMeans(perm[ref, , drop = FALSE])) + off   # back on the loading axis
    q <- apply(pr, 1, stats::quantile, c(0.025, 0.975), na.rm = TRUE); d$nlo <- (q[1, ] + shift) / log(2); d$nhi <- (q[2, ] + shift) / log(2)
  }
  ord <- if (ordering == "loading") order(d$value) else order(d$sig, -d$padj, abs(d$value))
  d$label <- ifelse(d$reference, paste0(d$celltype, " (ref)"), d$celltype)
  d$label <- factor(d$label, levels = d$label[ord])
  sig.lab <- sprintf("padj < %.2g", alpha)
  d$sig.f <- factor(ifelse(d$sig, sig.lab, "not significant"), levels = c(sig.lab, "not significant"))
  pal <- if (is.null(palette)) stats::setNames(rep("grey25", nrow(d)), d$celltype) else palette[d$celltype]
  pal[is.na(pal)] <- "grey40"; names(pal) <- d$celltype
  # the other frame's zero: the geometric mean of all cell types on the reference-centred axis, the reference block on the clr axis
  other <- if (center == "reference") -off / log(2) else off / log(2); other.lab <- if (center == "reference") "geometric mean" else "reference"
  gg <- ggplot2::ggplot(d, ggplot2::aes(y = .data$label)) + ggplot2::geom_vline(xintercept = 0, colour = "grey70")
  if (length(ref) && abs(other) > 1e-8) gg <- gg + ggplot2::geom_vline(xintercept = other, colour = "grey55", linetype = 2)
  if (show.null && any(is.finite(d$nlo)))
    gg <- gg + ggplot2::geom_linerange(ggplot2::aes(xmin = .data$nlo, xmax = .data$nhi, alpha = "null: 95 % of permutations"), linewidth = 4.5, colour = "grey88", na.rm = TRUE) +
      ggplot2::scale_alpha_manual(values = 1, name = NULL, guide = ggplot2::guide_legend(override.aes = list(linewidth = 4.5, colour = "grey88"), order = 2))
  if (!is.null(min.effect) && is.finite(min.effect) && min.effect > 1)
    gg <- gg + ggplot2::geom_vline(xintercept = c(-1, 1) * log2(min.effect), colour = "grey45", linetype = 3)
  if (any(is.finite(d$lo)))
    gg <- gg + ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lo, xmax = .data$hi, colour = .data$celltype, linewidth = sprintf("%g %% confidence interval", 100 * interval.level)), na.rm = TRUE, show.legend = c(linewidth = TRUE, colour = FALSE)) +
      ggplot2::scale_linewidth_manual(values = 0.5, name = NULL, guide = ggplot2::guide_legend(override.aes = list(colour = "grey25"), order = 3))
  gg <- gg + ggplot2::geom_point(ggplot2::aes(x = .data$value, colour = .data$celltype, fill = .data$celltype, shape = .data$sig.f), size = 2.8, stroke = 0.9,
                                 show.legend = c(shape = TRUE, alpha = FALSE, colour = FALSE, fill = FALSE)) +
    ggplot2::geom_point(data = d[!d$sig, , drop = FALSE], ggplot2::aes(x = .data$value, colour = .data$celltype), shape = 21, fill = "white", size = 2.8, stroke = 0.9,
                        show.legend = FALSE) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(21, 21), c(sig.lab, "not significant")), drop = FALSE, name = NULL,
                                guide = ggplot2::guide_legend(override.aes = list(colour = "grey25", fill = c("grey25", "white")), order = 1)) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") + ggplot2::scale_fill_manual(values = pal, guide = "none")
  # the marker's label goes after the layers that set the discrete y scale (an annotation with a numeric y placed first would make the scale continuous)
  if (length(ref) && abs(other) > 1e-8) gg <- gg + ggplot2::annotate("text", x = other, y = 0.45, label = other.lab, size = 2.4, colour = "grey40", hjust = if (other < 0) 1.05 else -0.05)
  xr <- range(c(0, other, d$value, if (show.null) c(d$nlo, d$nhi), d$lo, d$hi, if (!is.null(min.effect) && is.finite(min.effect) && min.effect > 1) c(-1, 1) * log2(min.effect)), na.rm = TRUE)
  if (show.pvals) {
    xt <- xr[2] + 0.07 * diff(xr)   # clear of the largest point
    gg <- gg + ggplot2::geom_text(ggplot2::aes(x = xt, label = ifelse(is.finite(.data$padj), format(signifBelow(.data$padj, alpha), drop0trailing = TRUE), "")),
                                  hjust = 0, size = 2.6, colour = "grey35") +
      ggplot2::annotate("text", x = xt, y = nrow(d) + 0.7, label = "padj", hjust = 0, size = 2.6, colour = "grey35", fontface = "italic")
    xr[2] <- xr[2] + 0.18 * diff(xr)
  }
  xlab <- sprintf("loading (log2)   (<- %s   |   %s ->)", ref.level, target.level)   # the frame is named by the dashed marker and, with center = "both", the panel title
  gg + ggplot2::coord_cartesian(xlim = xr, ylim = c(0.3, nrow(d) + if (show.pvals) 0.7 else 0), clip = "off") +
    plot.theme + ggplot2::labs(x = xlab, y = NULL, title = title) +
    ggplot2::theme(legend.position = legend.position, panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(),
                   legend.box = "vertical", legend.spacing.y = ggplot2::unit(0, "pt"), legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.margin = ggplot2::margin(0, 0, 0, 0))
}

#' Table of the compositional loadings of a contrast
#'
#' One row per cell type: the loading (log2; the change of the cell type's log proportion relative to the geometric mean of all cell
#' types, or relative to the reference cell types with `center = "reference"`), its model-based standard error and confidence
#' interval (least-squares fit on the log ratios), the permutation p-values, and, with `min.effect`, whether the interval lies
#' entirely beyond that fold change (`beyond.min.effect`: the change is at least `min.effect`-fold at the given confidence).
#' The minimal-effect column is a statement about the size of the change; the test of any change is the permutation p-value.
#' `median.cells` and `n.zero` (samples without the cell type) show how many cells a loading rests on, and `pseudocount.shift`
#' how far it moves, in standard errors, when the pseudocount is halved or doubled: values near or above 1 mark loadings whose
#' size (not their significance, which simulations show to be insensitive) depends on the treatment of small counts.
#'
#' @param res result of `lmCoda()` (as stored by `cao$estimateCellLoadings()`)
#' @param center `"none"` (default; relative to the geometric mean) or `"reference"`
#' @param level confidence level (default 0.95)
#' @param min.effect optional fold change (> 1), e.g. 1.5
#' @return data.frame ordered by the p-value
#' @export
codaLoadingsTable <- function(res, center = c("none", "reference"), level = 0.95, min.effect = NULL) {
  center <- match.arg(center); ld <- res$contrast$loadings; pc <- res$contrast$per_cell; cts <- names(ld$obs)
  l <- (if (center == "reference") ld$ref_centered else ld$obs)[cts] / log(2); se <- (if (center == "reference") ld$se.ref else ld$se)
  se <- if (is.null(se)) rep(NA_real_, length(cts)) else unname(se[cts]) / log(2); tq <- if (is.null(ld$df)) NA_real_ else stats::qt(1 - (1 - level) / 2, ld$df)
  pv <- if (center == "reference") pc$pval.ref else pc$pval; pa <- if (center == "reference") pc$padj.ref else pc$padj
  d <- data.frame(celltype = cts, loading = unname(l), fold.change = 2^unname(l), se = se, ci.low = unname(l) - tq * se, ci.high = unname(l) + tq * se,
                  pval = unname(pv[cts]), padj = unname(pa[cts]), reference = cts %in% res$reference$celltypes, stringsAsFactors = FALSE)
  if (!is.null(ld$counts)) { k <- match(cts, ld$counts$celltype); d$median.cells <- ld$counts$median.cells[k]; d$n.zero <- ld$counts$n.zero[k] }
  if (!is.null(ld$pseudocount.shift)) d$pseudocount.shift <- unname(ld$pseudocount.shift[cts])
  if (!is.null(min.effect)) { th <- log2(min.effect); d$beyond.min.effect <- is.finite(d$ci.low) & (d$ci.low > th | d$ci.high < -th) }
  d <- d[order(d$pval, -abs(d$loading)), ]; rownames(d) <- NULL
  attr(d, "level") <- level; attr(d, "min.effect") <- min.effect; d
}

#' Per-sample diagnostics of the composition model
#'
#' For every sample and cell type: the count, the observed proportion, the log ratio (clr value) the model was fitted to and its
#' fitted value, the studentized residual, and the influence of the sample on the cell type's loading (the change of the loading
#' when the sample is left out, in standard errors of the loading; exact for the least-squares fit, no refit needed).
#'
#' @param res result of `lmCoda()` (as stored by `cao$estimateCellLoadings()`), with the counts (`res$cnts`)
#' @return long data.frame: `sample`, `celltype`, `count`, `n.cells` (of the sample), `proportion` (count / n.cells), `value`
#'   (proportion with the pseudocount, as fitted), `clr`, `fitted.clr`, `fitted.proportion`, `residual` (studentized),
#'   `influence` (change of the loading without the sample, in standard errors), `leverage`
#' @export
codaSampleDiagnostics <- function(res) {
  if (is.null(res$cnts) || is.null(res$ilr)) stop("the result holds no counts; run estimateCellLoadings() again")
  cn <- as.matrix(res$cnts); design <- res$model$F; cvec <- res$model$contrast.F; cts <- rownames(res$psi)
  Y <- res$ilr %*% t(res$psi); colnames(Y) <- cts; rownames(Y) <- rownames(cn)
  q <- qr(design); rk <- q$rank; df <- nrow(design) - rk; fit <- qr.fitted(q, Y); rs <- Y - fit
  h <- pmin(rowSums(qr.Q(q)[, seq_len(rk), drop = FALSE]^2), 1); s <- sqrt(colSums(rs^2) / max(df, 1))
  rstud <- sweep(rs, 2, s, "/") / sqrt(pmax(1 - h, 1e-8))
  pv <- q$pivot[seq_len(rk)]; XtXi <- chol2inv(qr.R(q))[seq_len(rk), seq_len(rk), drop = FALSE]
  a <- drop(design[, pv, drop = FALSE] %*% XtXi %*% cvec[pv]); se <- sqrt(sum(cvec[pv] * (XtXi %*% cvec[pv]))) * s
  infl <- sweep(-(a / pmax(1 - h, 1e-8)) * rs, 2, se, "/")            # loading without sample i minus the loading, in standard errors
  fp <- exp(fit); fp <- fp / rowSums(fp); n <- rowSums(cn)
  data.frame(sample = rep(rownames(cn), ncol(cn)), celltype = rep(cts, each = nrow(cn)), count = as.vector(cn[, cts]), n.cells = rep(n, ncol(cn)),
             proportion = as.vector(cn[, cts] / n), value = as.vector(res$freqs[, cts]), clr = as.vector(Y) / log(2), fitted.clr = as.vector(fit) / log(2),
             fitted.proportion = as.vector(fp), residual = as.vector(rstud), influence = as.vector(infl), leverage = rep(h, ncol(cn)), stringsAsFactors = FALSE)
}

#' Samples behind the composition result, per cell type
#'
#' One panel per cell type: the samples' proportions by level of the tested variable (logarithmic axis), with the model's fitted
#' proportion at the two endpoints of the contrast and its confidence interval behind them (box and line; from the fit, adjusted
#' for the model's other terms). Samples in which the cell type was not observed are drawn as open points at the value the
#' pseudocount gives them. Samples that carry the result are named: those whose omission changes the cell type's loading by more
#' than `influence.se` standard errors, and those with a studentized residual beyond `resid.thresh`. With `pair.by`, samples of
#' the same unit (e.g. patient) are joined.
#'
#' @param res result of `lmCoda()` with the counts (as stored by `cao$estimateCellLoadings()`)
#' @param cell.types cell types to show (default: the `n.top` with the smallest p-values)
#' @param n.top number of cell types shown by default (12)
#' @param pair.by optional column of the sample metadata; samples sharing a value are joined by a line
#' @param influence.se,resid.thresh thresholds for naming a sample (defaults 1 standard error of the loading; residual 3)
#' @param palette named colours of the levels of the tested variable (optional)
#' @param label.samples name the flagged samples (default TRUE)
#' @param ncol number of panel columns (default: chosen from the number of cell types)
#' @param plot.theme ggplot2 theme
#' @return ggplot object
#' @export
plotCodaSamples <- function(res, cell.types = NULL, n.top = 12, pair.by = NULL, influence.se = 1, resid.thresh = 3, palette = NULL, label.samples = TRUE,
                            ncol = NULL, plot.theme = ggplot2::theme_bw()) {
  d <- codaSampleDiagnostics(res); meta <- res$model$meta; tv <- res$model$tested.variable
  if (is.null(tv) || !tv %in% names(meta)) stop("the tested variable is not a column of the model's sample metadata")
  pc <- res$contrast$per_cell; pval <- pc$pval; padj <- pc$padj
  if (is.null(cell.types)) cell.types <- names(sort(pval))[seq_len(min(n.top, length(pval)))]
  cell.types <- intersect(cell.types, unique(d$celltype)); if (!length(cell.types)) stop("no such cell types in the result")
  d <- d[d$celltype %in% cell.types, , drop = FALSE]
  g <- meta[d$sample, tv]; num <- is.numeric(g); spec <- res$model$contrast_spec; lb <- res$contrast$predicted$labels
  lab <- list(baseline = as.character(spec$den %||% sub("^.* = ", "", lb$baseline %||% "")), target = as.character(spec$num %||% sub("^.* = ", "", lb$target %||% "")))   # the levels of the contrast endpoints
  if (!num) { lv <- unique(c(lab$baseline, lab$target, if (is.factor(g)) levels(droplevels(g)) else sort(unique(as.character(g))))); lv <- lv[lv %in% as.character(g)]; d$x <- factor(as.character(g), levels = lv) } else d$x <- g
  d$observed <- factor(ifelse(d$count > 0, "observed", "no cells (pseudocount)"), levels = c("observed", "no cells (pseudocount)"))
  d$flag <- abs(d$influence) > influence.se | abs(d$residual) > resid.thresh
  pr0 <- res$contrast$predicted; fc <- if (!is.null(pr0)) pr0$target[cell.types] / pr0$baseline[cell.types] else rep(NA_real_, length(cell.types))
  fc.lab <- ifelse(!is.finite(fc), "", ifelse(fc >= 1, sprintf("%s-fold higher, ", format(signif(fc, 2), drop0trailing = TRUE)), sprintf("%s-fold lower, ", format(signif(1 / fc, 2), drop0trailing = TRUE))))
  strip <- setNames(sprintf("%s\n%spadj %s", cell.types, fc.lab, format(signif(padj[cell.types], 2), drop0trailing = TRUE)), cell.types)
  d$panel <- factor(strip[d$celltype], levels = strip[cell.types[order(pval[cell.types])]])
  gg <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$value))
  pr <- res$contrast$predicted
  if (!num && !is.null(pr) && all(c(lab$baseline, lab$target) %in% levels(d$x))) {      # fitted endpoint proportions and their intervals
    f <- rbind(data.frame(celltype = cell.types, x = lab$baseline, y = pr$baseline[cell.types], lo = (pr$baseline.lo %||% pr$baseline)[cell.types], hi = (pr$baseline.hi %||% pr$baseline)[cell.types]),
               data.frame(celltype = cell.types, x = lab$target, y = pr$target[cell.types], lo = (pr$target.lo %||% pr$target)[cell.types], hi = (pr$target.hi %||% pr$target)[cell.types]))
    f$x <- factor(f$x, levels = levels(d$x)); f$panel <- factor(strip[f$celltype], levels = levels(d$panel)); f$lo <- pmax(f$lo, min(d$value) / 3)
    gg <- gg + ggplot2::geom_crossbar(data = f, ggplot2::aes(x = .data$x, y = .data$y, ymin = .data$lo, ymax = .data$hi, fill = "fitted proportion and its confidence interval"), colour = "grey55", width = 0.55, linewidth = 0.25, middle.linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = "grey94", name = NULL)
  }
  if (!is.null(pair.by)) {
    if (!pair.by %in% names(meta)) stop("pair.by is not a column of the sample metadata")
    d$pair <- as.character(meta[d$sample, pair.by]); gg <- gg + ggplot2::geom_line(data = d, ggplot2::aes(group = .data$pair), colour = "grey70", linewidth = 0.3)
  }
  pos <- if (is.null(pair.by) && !num) ggplot2::position_jitter(width = 0.12, height = 0, seed = 1) else ggplot2::position_identity()
  gg <- gg + ggplot2::geom_point(data = d, ggplot2::aes(colour = if (num) NULL else .data$x, shape = .data$observed), size = 1.8, stroke = 0.7, position = pos) +
    ggplot2::scale_shape_manual(values = c(observed = 16, `no cells (pseudocount)` = 1), drop = FALSE, name = NULL)
  if (!num) gg <- gg + (if (is.null(palette)) ggplot2::scale_colour_brewer(palette = "Set1", guide = "none") else ggplot2::scale_colour_manual(values = palette, guide = "none"))
  if (label.samples && any(d$flag)) {
    checkPackageInstalled("ggrepel", cran = TRUE)
    gg <- gg + ggrepel::geom_text_repel(data = d[d$flag, , drop = FALSE], ggplot2::aes(label = .data$sample), size = 2.3, colour = "grey20", min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf, position = pos, seed = 1)
  }
  pct <- function(v) paste0(format(100 * v, drop0trailing = TRUE, scientific = FALSE, trim = TRUE), " %")
  gg + ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = ncol %||% min(6, ceiling(sqrt(1.6 * length(cell.types))))) + ggplot2::scale_y_log10(labels = pct) + plot.theme +
    ggplot2::labs(x = if (num) tv else NULL, y = "proportion of the sample's cells",
                  caption = sprintf("named samples: omission changes the loading by more than %g standard error%s, or studentized residual beyond %g", influence.se, if (influence.se == 1) "" else "s", resid.thresh)) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 8, hjust = 0, lineheight = 0.95), strip.background = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 7), plot.caption = ggplot2::element_text(size = 7, colour = "grey35"))
}

#' old two-group comparison: reference set selection
#' @keywords internal
referenceSet <- function(freqs, groups, p.thresh=0.05) {
  checkPackageInstalled("psych", cran=TRUE)
  freqs[freqs == 0] <- min(freqs[freqs != 0])/2
  cell.types <- colnames(freqs)
  cell.list <- lapply(cell.types, function(x) x)

  for (it in 1:(length(cell.list) - 2)) {

    mx <- matrix(0, nrow=length(cell.list), ncol=length(cell.list))
    for(i in 1:length(cell.list)){
      for(j in 1:length(cell.list)){
        if (j <= i) next

        ratio <- log(apply(freqs[,cell.list[[i]],drop=FALSE], 1, psych::geometric.mean )) -
          log(apply(freqs[,cell.list[[j]],drop=FALSE], 1, psych::geometric.mean ))
      
        mod <- lm(groups ~ ratio)
        res <- summary(mod)
        pval <- res$coefficients[2,4]
        mx[i, j] <- pval
      }
    }
    if(max(max(mx)) < p.thresh) break
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)
    cell.list[[length(cell.list) + 1]] <- c(cell.list[[i]], cell.list[[j]])
    cell.list <- cell.list[-c(i,j)]

    if(it == 1){
      mx.first <- mx
      rownames(mx.first) <- cell.types
      colnames(mx.first) <- cell.types
    }
  }
  len.cell.list <- sapply(cell.list, length)
  id.ref <- which(len.cell.list == max(len.cell.list))

  # Make mx.first symmetric
  if (!exists("mx.first")) {
    mx.first <- mx
  }
  
  for(i in 1:nrow(mx.first)){
    for(j in 1:ncol(mx.first)){
      if (j <= i) next
      mx.first[j, i] <- mx.first[i, j]
    }
  }

  return(list(ref.set = cell.list[[id.ref[1]]],
              cell.list = cell.list,
              mx.first = mx.first,
              mx.final = mx))
}

#' @keywords internal
sbpDiff <- function(freqs, groups){
  freqs[freqs == 0] <- min(freqs[freqs != 0])/2
  cell.types <- colnames(freqs)
  cell.list <- lapply(cell.types, function(x) x)

  sbp <- matrix(nrow = 0, ncol = length(cell.types), dimnames = list(c(), cell.types))
  for(it in 1:(length(cell.list) - 1)){

    mx <- matrix(0, nrow = length(cell.list), ncol = length(cell.list))
    for(i in 1:length(cell.list)){
      for(j in 1:length(cell.list)){
        if (j <= i) next

        ratio <- log(apply(freqs[,cell.list[[i]],drop=FALSE], 1, psych::geometric.mean )) -
          log(apply(freqs[,cell.list[[j]],drop=FALSE], 1, psych::geometric.mean ))

        mod <- lm(groups ~ ratio)
        res <- summary(mod)
        pval <- res$coefficients[2,4]
        mx[i, j] <- pval
      }
    }
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)

    sbp.tmp <- rep(0, length(cell.types))
    sbp.tmp[which(cell.types %in% cell.list[[i]])] <- 1
    sbp.tmp[which(cell.types %in% cell.list[[j]])] <- -1
    sbp <- rbind(sbp.tmp, sbp)
    cell.list[[length(cell.list) + 1]] <- c(cell.list[[i]], cell.list[[j]])
    cell.list <- cell.list[-c(i,j)]
  }

  return(t(sbp))
}

#' Check cell count data
#' 
#' @param d.counts Cell count table
#' @keywords internal
checkData <- function(d.counts){
  if(is.null(d.counts)) 
    stop('Cell count matrix is not provided')
  if(nrow(d.counts) == 0) 
    stop('Sample size is zero')
  if(ncol(d.counts) == 0) 
    stop('Cell count information is missed')
}

#' Check groups
#' 
#' @param d.groups Groups variable for samples
#' @keywords internal
checkGroups <- function(d.groups){
  if(is.null(d.groups)) 
    stop('Groups are not provided')
  if(length(d.groups) == 0) 
    stop('Group information is missed')
}

#' Check cell count data and groups
#' 
#' @param d.counts Cell count table
#' @param d.groups Groups variable for samples
#' @keywords internal
checkDataGroups <- function(d.counts, d.groups){
  checkData(d.counts)
  checkGroups(d.groups)
  if(nrow(d.counts) != length(d.groups)) 
    stop('Sample size in cell count matrix and in group variable is not the same')
}

#' Check sequential binary partitions(sbp)
#' 
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkSbp <- function(sbp){
  if(is.null(sbp)) 
    stop('Sequential binary partitions are not provided')
  if(nrow(sbp) == 0)
    stop('No partitions are provided')
}

#' Check sequential binary partitions(sbp) for the whole set of balances
#' 
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkSbpWhole <- function(sbp){
  checkSbp(sbp)
  if(nrow(sbp) < 2) 
    stop('Balances require at least 2 cell types')
  if((nrow(sbp) - ncol(sbp)) != 1) 
    stop('Wrong number of balances provided')
}

#' Check the agreement between sbp and cell count table
#' 
#' @param d.counts Cell count table
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkDataSbp <- function(d.counts, sbp){
  checkData(d.counts)
  checkSbp(sbp)
  if(ncol(sbp) != ncol(d.counts))
    stop('Number of cell types in data and sbp should be equal')
}

#' Check data and list of cells
#' 
#' @param d.counts Cell count table
#' @param cells List of cell types of interest
#' @keywords internal
checkDataAndCells <- function(d.counts, cells){
  if(is.null(cells))
    stop('Cell types are not provided')
  if(length(cells) == 0)
    stop('Cell types are not provided')
  if(sum(cells %in% colnames(d.counts)) != length(cells))
    stop('Some cell types do not exist in the data')
}

#' Check Tree
#' 
#' @param tree phylo tree
#' @keywords internal
checkTree <- function(tree){
  if(is.null(tree))
    stop('Tree is not provided')
}

#' Check agreement between Tree and the Cell count table
#' 
#' @param d.counts Cell count table
#' @param tree phylo tree
#' @keywords internal
checkDataTree <- function(d.counts, tree){
  checkTree(tree)
  checkData(d.counts)
  checkDataAndCells(d.counts, tree$tip.label)
}

#' Construct the canonical tree (old two-group comparison)
#'
#' @param cnts Table with cell type counts
#' @param groups Groups variable for samples
#' @return phylo tree
#' @keywords internal
constructTree <- function(cnts, groups, partition.thresh = 0){
  checkDataGroups(cnts, groups)
  n.cells <- ncol(cnts)
  sbp.cda <- matrix(0, nrow = n.cells, ncol = n.cells-1, dimnames = list(colnames(cnts), c()))
  
  unsolved.cells <- list(rownames(sbp.cda))  # List of cell types to separate
  for(id.bal in 1:(n.cells-1)){
    # If list for separation contains two cell types: perform random separation
    if(length(unsolved.cells[[id.bal]]) == 2){
      sbp.cda[unsolved.cells[[id.bal]][1],id.bal] <- 1
      sbp.cda[unsolved.cells[[id.bal]][2],id.bal] <- -1
      next
    }
    #  Data for the current subset of cells
    d.tmp <- cnts[, colnames(cnts) %in% unsolved.cells[[id.bal]]]
    
    # Define the most contrast balance
    can.loadings <- getLoadings(d.tmp, groups)
    cells.tmp <- rownames(can.loadings)
    
    # Get cell types from opposite sides of the principal balance
    cells.plus <- cells.tmp[can.loadings > partition.thresh]
    cells.minus <- cells.tmp[can.loadings < partition.thresh]
    
    sbp.cda[cells.minus, id.bal] <- -1
    sbp.cda[cells.plus, id.bal] <- 1
    
    if(length(cells.minus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.minus
    }
    
    if(length(cells.plus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.plus
    }
  }
  
  tree <- sbp2tree(sbp.cda)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)
  
  return(list(tree = tree, sbp = sbp.cda, dendro = d.cur))
}

#' (old two-group comparison)
#' @keywords internal
constructTreeUp <- function(freqs, groups) {
  
  sbp <- sbpDiff(freqs, groups)
  
  tree <- sbp2tree(sbp)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)
  return(list(tree = tree, sbp = sbp, dendro = d.cur))
}

#' (old two-group comparison)
#' @keywords internal
constructTreeUpDown <- function(cnts, groups){
  cnts[cnts <= 0] <- 0.5
  
  ref.set <- referenceSet(cnts, groups)
  cell.lists <- ref.set$cell.list
  cell.types <- colnames(cnts)
  
  # Constrtuct sbp from lists
  freqs <- (cnts)/rowSums(cnts)
  freqs.lists <- c()
  for (i in 1:length(cell.lists)) {
    freqs.lists <- cbind(freqs.lists, apply(freqs[,cell.lists[[i]],drop=FALSE], 1, psych::geometric.mean))
  }
  colnames(freqs.lists) <- paste('tmp', 1:length(cell.lists), sep = '')
  t.list <- constructBestPartitionTree(freqs.lists, groups)
  sbp.list <- t.list$sbp
  
  sbp.all <- matrix(ncol = 0, nrow = length(cell.types), dimnames = list(cell.types, c()))
  for(k in 1:ncol(sbp.list)){
    p <- sbp.list[,k]
    i.list.plus <- which(p > 0)
    i.list.minus <- which(p < 0)
    
    sbp.tmp <- rep(0, length(cell.types))
    for(i.plus in i.list.plus){
      sbp.tmp[cell.types %in% cell.lists[[i.plus]]] = 1
    }
    for(i.plus in i.list.minus){
      sbp.tmp[cell.types %in% cell.lists[[i.plus]]] = -1
    }
    sbp.all <- cbind(sbp.all, sbp.tmp)
  }
  
  # Construct sbps from sublists
  for(i in 1:length(cell.lists)){
    freqs.tmp <- freqs[,cell.lists[[i]], drop = FALSE]
    if(ncol(freqs.tmp) == 1) next
    
    tmp <- constructBestPartitionTree(freqs.tmp, groups)
    sbp.tmp <- tmp$sbp
    
    sbp.add <- matrix(0, nrow = nrow(sbp.all), ncol = ncol(sbp.tmp))
    rownames(sbp.add) <- rownames(sbp.all)
    sbp.add[rownames(sbp.tmp),] <- sbp.tmp
    sbp.all <- cbind(sbp.all, sbp.add)
  }
  
  tree <- sbp2tree(sbp.all)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)
  return(list(tree = tree, sbp = sbp.all, dendro = d.cur))
}


#' Get tree from balances
#'
#' @param sbpart A contrast matrix for balances, sequential binary partition
#' @return phylo tree
#' @details We use ellipsis because a previous version of this function
#'          contains additional and insignificant parameters.
#' @keywords internal
sbp2tree_old <- function(sbpart){
  checkSbpWhole(sbpart)
  
  n.cells <- nrow(sbpart)
  edges <- c(n.cells+1, n.cells+1)
  id <- 2*n.cells - 1
  
  collapse <- sbpart
  rownames(collapse) <- 1:n.cells
  
  for(i in rev(1:(n.cells - 1)) ){
    ids <- which(collapse[,i] != 0)
    
    edges <- rbind(edges, as.integer(c(id, rownames(collapse)[ids[1]])))
    edges <- rbind(edges, as.integer(c(id, rownames(collapse)[ids[2]])))
    rownames(collapse)[ids[2]] <- id
    id <- id - 1
    
    collapse[ids[1],] <- 0
  }
  
  # Create face tree, because we need a tidytree object
  x <- tidytree::as_tibble(rtree(n = n.cells, br = NULL))
  x[,1:2] <- edges
  
  t.tmp <- as.phylo(x)
  t.tmp$tip.label <- rownames(sbpart)
  return(t.tmp)
}

#' Get tree from balances (safe phylo constructor)
#' @keywords internal
sbp2tree <- function(sbpart){
  checkSbpWhole(sbpart)

  n.cells <- nrow(sbpart)
  edges <- c(n.cells+1, n.cells+1)
  id <- 2*n.cells - 1

  collapse <- sbpart
  rownames(collapse) <- 1:n.cells

  #  map each balance column -> internal node id
  node.id.by.balance <- integer(n.cells - 1)

  for(i in rev(1:(n.cells - 1)) ){
    ids <- which(collapse[,i] != 0)

    node.id.by.balance[i] <- id  # NEW

    edges <- rbind(edges, as.integer(c(id, rownames(collapse)[ids[1]])))
    edges <- rbind(edges, as.integer(c(id, rownames(collapse)[ids[2]])))
    rownames(collapse)[ids[2]] <- id
    id <- id - 1

    collapse[ids[1],] <- 0
  }

  edge <- edges[-1, , drop=FALSE]
  phy <- list(edge = edge, tip.label = rownames(sbpart), Nnode = n.cells - 1)
  class(phy) <- "phylo"

  phy <- ape::reorder.phylo(phy)

  attr(phy, "balance_node_id") <- node.id.by.balance  # NEW
  phy
}



#' Get balances for all inner nodes of the tree
#'
#' @param t A phylo object
#' @param log.f Logarithms of frequencies
#' @return Balances, contrast matrix and names on cell types
#' @keywords internal
getNodeBalances <- function(log.f, sbp.my){
  
  # get psi matrix
  psi <- sbp.my * 0
  for(k in 1:ncol(sbp.my))
  {
    p <- sbp.my[, k]  # k-th partition
    n1 <- sum(p > 0)
    n2 <- sum(p < 0)
    p[p > 0] <- 1/n1
    p[p < 0] <- -1/n2
    psi[, k] <- p * sqrt(n1*n2/(n1+n2))
  }
  
  balances <- log.f %*% psi
  return(balances)
}

#' @keywords internal
bestPartition <- function(freqs.tmp, groups){
  
  
  n.cells <- ncol(freqs.tmp)
  sbp <- 1
  for(i in 2:n.cells){
    sbp <- rbind(cbind(sbp, 1),
                 cbind(sbp, -1))
  }
  sbp <- sbp[-c(1, nrow(sbp)),]
  log.f <- log(freqs.tmp)
  
  for(i in 1:nrow(sbp)){
    p <- sbp[i,]
    p[p < 0] <- p[p < 0] / sum(p < 0)
    p[p > 0] <- p[p > 0] / sum(p > 0)
    sbp[i,] <- p
  }
  colnames(sbp) <- colnames(freqs.tmp)
  
  
  bals <- sbp %*% t(log.f)
  p <- c()
  for(i in 1:nrow(bals)){
    b <- bals[i,]
    mod <- lm(groups ~ b)
    res <- summary(mod)
    pval <- res$coefficients[2,4]
    p <- c(p, pval)
  }
  
  sbp.best <- (sbp[p == min(p),] > 0) * 2 - 1
  names(sbp.best) <- colnames(freqs.tmp)
  
  return(sbp.best)
}



#' Construct the canonical tree (old two-group comparison)
#'
#' @param cnts Table with cell type counts
#' @param groups Groups variable for samples
#' @return phylo tree
#' @keywords internal
constructBestPartitionTree <- function(cnts, groups, partition.thresh = 0){
  checkDataGroups(cnts, groups)
  n.cells <- ncol(cnts)
  sbp.cda <- matrix(0, nrow = n.cells, ncol = n.cells-1, dimnames = list(colnames(cnts), c()))
  
  unsolved.cells <- list(rownames(sbp.cda))  # List of cell types to separate
  for(id.bal in 1:(n.cells-1)){
    # If list for separation contains two cell types: perform random separation
    if(length(unsolved.cells[[id.bal]]) == 2){
      sbp.cda[unsolved.cells[[id.bal]][1],id.bal] <- 1
      sbp.cda[unsolved.cells[[id.bal]][2],id.bal] <- -1
      next
    }
    #  Data for the current subset of cells
    d.tmp <- cnts[, colnames(cnts) %in% unsolved.cells[[id.bal]]]
    
    # Find the best partition
    sbp.best <- bestPartition(d.tmp, groups)
    
    # Get cell types from opposite sides of the principal balance
    cells.tmp <- colnames(d.tmp)
    
    cells.plus <- cells.tmp[sbp.best > 0]
    cells.minus <- cells.tmp[sbp.best < 0]
    
    sbp.cda[cells.minus, id.bal] <- -1
    sbp.cda[cells.plus, id.bal] <- 1
    
    if(length(cells.minus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.minus
    }
    
    if(length(cells.plus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.plus
    }
    
  }
  
  tree <- sbp2tree(sbp.cda)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)
  
  return(list(tree = tree, sbp = sbp.cda, dendro = d.cur))
}

#' @keywords internal
sbpInNodes <- function(tree){
  edges <- tree$edge
  n.nodes <- max(edges)
  n.leaves <- length(tree$tip.label)
  edge.lists <- list()
  for(i in 1:n.leaves){
    edge.lists[[i]] <- c(tree$tip.label[i])
  }
  
  sbp <- matrix(0, nrow = n.nodes, ncol = n.leaves, dimnames = list(c(), tree$tip.label))
  for(i in rev((n.leaves+1):n.nodes)){
    idx <- which(edges[,1] == i)
    cells.plus <- edge.lists[[edges[idx[1],2] ]]
    cells.minus <- edge.lists[[edges[idx[2],2] ]]
    sbp[i, cells.plus] <- 1
    sbp[i, cells.minus] <- -1
    
    edge.lists[[i]] <- c(cells.plus, cells.minus)
  }
  sbp <- sbp[-c(1:n.leaves),]
  sbp <- t(sbp)
  
  return(sbp)
}

#' Generate one random partition of given length
#'
#' @param bal.len length of partition
#' @return random partition - vector with values from {1, -1}
#' @keywords internal
getRndPartition <- function(bal.len){
  if(bal.len < 2)
    stop('Balance cannot be generated')
  
  bal.values <- c(0)
  while(length(unique(bal.values)) == 1){
    bal.values <- sample(c(-1, 1), bal.len, replace = TRUE)
  }
  return(bal.values)
}


#' Generate random sbp
#'
#' @param n.cell.types number of cell types
#' @return list with binary and normalized sbp matrices (partitions are in rows)
#' @keywords internal
getRndSbp <- function(n.cell.types, n.seed=239){
  if(n.cell.types < 2)
    stop('For a partition at least 2 cell types should be provided')
  
  set.seed(n.seed)
  sbp.bin <- matrix(1, 1, n.cell.types)
  bal.stack <- c(1) # stack of untreated balances
  while(length(bal.stack) > 0){
    
    parent.bal <- sbp.bin[bal.stack[1],]  # parental balance
    bal.stack <- bal.stack[-1]
    
    for(bal.direction in c(-1, 1)){
      if(sum(parent.bal == bal.direction) < 2) next
      tmp.bal <- parent.bal * 0
      tmp.bal[parent.bal == bal.direction] <- getRndPartition(sum(parent.bal == bal.direction))
      sbp.bin <- rbind(sbp.bin, tmp.bal)
      bal.stack <- c(bal.stack, nrow(sbp.bin))
      
    }
  }
  sbp.bin <- sbp.bin[-1,] # remove initialization
  
  sbp.norm <- sbp.bin
  for(i in 1:nrow(sbp.norm)){
    n1 <- sum(sbp.norm[i,] == 1)
    n2 <- sum(sbp.norm[i,] == -1)
    sbp.norm[i,sbp.norm[i,] == 1] <- 1/n1
    sbp.norm[i,sbp.norm[i,] == -1] <- -1/n2
    sbp.norm[i,] <- sbp.norm[i,] * sqrt(n1*n2/(n1+n2))
  }
  
  return(list(bin = sbp.bin, norm = sbp.norm))
}


#' Get balances based on the random ilr basis
#'
#' @param d.counts Cell count table
#' @param n.seed Random seed
#' @return Matrix with balances in columns
#' @keywords internal
getRndBalances <- function(d.counts, n.seed=239){
  checkData(d.counts)
  
  n.cell.types <- ncol(d.counts)
  sbp <- getRndSbp(n.cell.types, n.seed)
  
  log.freq <- getLogFreq(d.counts)
  
  psi <- as.matrix(t(sbp$norm))
  rownames(psi) <- colnames(d.counts)
  
  balances <- as.matrix(log.freq) %*% psi
  balances.norm <-  apply(balances, 2, function(y) y - mean(y))
  
  return(list(values = balances, norm = balances.norm, psi = psi))
}


#' Get log frequencies from the cell count table
#' @keywords internal
getLogFreq <- function(d.counts){
  checkData(d.counts)
  
  data.norm <- d.counts
  data.norm[data.norm == 0] <- 1
  
  log.freq <- t(apply(data.norm, 1, function(y) log(y/sum(y))))
  return(log.freq)
}

#' Remove the strongest group effect from balances
#'
#' @param d.used Correct values of balances
#' @param d.groups if provided then resampling controls presence of both groups in a new dataset
#' @param n.seed Random seed
#' @param thresh.pc.var percentage of variance which should be characterized by PSc
#' @return Balances without group effect
#' @keywords internal
removeGroupEffect <- function(d.used, d.groups, thresh.pc.var = 0.95){
  checkDataGroups(d.used, d.groups)
  
  pca.res <- prcomp(d.used, center = FALSE, scale. = FALSE)
  
  # Calculate cummularive variance explained by PCs
  expl.var <- cumsum(pca.res$sdev^2/sum(pca.res$sdev^2))
  n.pc.var <- sum(expl.var < thresh.pc.var) + 1
  n.pc.var <- ncol(pca.res$x) - 1
  
  for (i in 0:n.pc.var) {
    stop.loop <- TRUE
    d.working <- pca.res$x[,1:(n.pc.var - i)] # data to cda
    
    model<-lm(d.working ~ d.groups)
    tryCatch(suppressWarnings(cda <- candisc::candisc(model)), error = function(e) { stop.loop <<- FALSE})
    if(stop.loop) { break }
  }
  n.pc.var <- n.pc.var - i
  
  model<-lm(d.working ~ d.groups)
  cda <- candisc::candisc(model)
  
  cda.rotation <- cda$structure # already normalized
  
  d.scores <- d.working %*% cda.rotation  # Projection values
  d.part <- d.scores %*% t(cda.rotation)  # Projection vectors
  df.remain <- d.working - d.part  # Orthogonal part
  
  d.used.part <- d.part  %*% t(pca.res$rotation[,1:n.pc.var])
  d.used.remain <- df.remain %*% t(pca.res$rotation[,1:n.pc.var])
  
  cumm.rotation <- pca.res$rotation[,1:n.pc.var] %*% cda.rotation
  return(list(rotation = cumm.rotation,
              used.part = d.used.part,
              remain = d.used.remain,
              scores = d.scores,
              res = cda))
}
