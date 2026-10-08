## Pairwise expression-shift engine (individual-level model).
##
## Location model   : y_i = B' x_i + e_i            x_i: row of the sample-level design X
## Dispersion model : s_i = E||e_i||^2 = z_i' gamma  z_i: row of the dispersion design Z
## Implied pairwise model (never fitted): E[d^2_ij] = (x_i - x_j)' M (x_i - x_j) + s_i + s_j,  M = B B'.
##
## For a contrast with endpoints x_T (target/num) and x_R (reference/den), c = x_T - x_R:
##   shift = c' M c            var = 2 (s_T - s_R)         total = shift + var / 2
## In the one-factor case these equal the classical pair-cell quantities
##   m_RT - (m_RR + m_TT)/2,   m_TT - m_RR,   m_RT - m_RR.
## Everything is computed from the Gower-centred matrix G = -1/2 J D^2 J, so it works for any
## squared-Euclidean-type distance matrix (see asSquaredDistance()).

#' Map a distance matrix to the matrix used as a squared Euclidean distance
#'
#' @param D numeric sample x sample distance matrix (symmetric, zero diagonal)
#' @param dist distance type: `"cor"` (gene-centred cosine distance `1 - cos`, which already equals half a
#'   squared Euclidean distance of unit-norm centred profiles), `"l2"` (Euclidean; squared here) or `"l1"`
#'   (Manhattan; of negative type, used as is: location/dispersion separation is approximate)
#' @return numeric matrix playing the role of `D^2`
#' @export
asSquaredDistance <- function(D, dist = c("cor", "l2", "l1")) {
  dist <- match.arg(dist)
  D <- as.matrix(D)
  if (dist == "l2") D^2 else D
}

# Gower centring: G = -1/2 J D2 J with J = I - 11'/n
gowerCenter <- function(D2) {
  D2 <- as.matrix(D2); n <- nrow(D2)
  if (n == 0) return(D2)
  rm <- rowMeans(D2); gm <- mean(rm)
  G <- -0.5 * (D2 - outer(rm, rep(1, n)) - outer(rep(1, n), rm) + gm)
  dimnames(G) <- dimnames(D2)
  (G + t(G)) / 2
}

# Hat-matrix information for a (possibly rank-deficient) design
hatInfo <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) == 0) {
    n <- nrow(X)
    return(list(XtXi = matrix(0, 0, 0), A = matrix(0, 0, n), H = matrix(0, n, n), h = rep(0, n), rank = 0L))
  }
  XtXi <- MASS::ginv(crossprod(X))
  A <- XtXi %*% t(X)                       # q x n: maps responses to coefficients
  H <- X %*% A
  list(XtXi = XtXi, A = A, H = H, h = pmin(diag(H), 1), rank = qr(X)$rank)
}

# Orthonormal basis of {b : c'b = 0}; X %*% contrastNullBasis(c) is the reduced design with c'beta = 0 imposed
contrastNullBasis <- function(cvec) {
  Q <- qr.Q(qr(matrix(cvec, ncol = 1)), complete = TRUE)
  Q[, -1, drop = FALSE]
}

# Is the contrast c estimable for design X (c in the row space of X)?
isEstimable <- function(X, cvec, tol = 1e-8) {
  if (!length(cvec) || all(abs(cvec) < tol)) return(FALSE)
  XtX <- crossprod(as.matrix(X))
  P <- MASS::ginv(XtX) %*% XtX
  sqrt(sum((drop(P %*% cvec) - cvec)^2)) <= tol * max(1, sqrt(sum(cvec^2)))
}

#' Pairwise effects (shift / var / total) from a squared-distance matrix and individual-level designs
#'
#' Core estimator of the individual-level model. `X` is the location design, `Z` the dispersion design;
#' `contrast` is a numeric contrast over `colnames(X)`; `z.end` gives the rows of `Z` at the contrast
#' endpoints (`num` = target, `den` = reference). Returns the three effects, their normalized versions,
#' the estimated effect matrix and per-sample dispersions.
#'
#' @param D2 squared-type distance matrix (see [asSquaredDistance()]); ignored when `G` is given
#' @param X location design matrix (n x q)
#' @param contrast numeric contrast vector over the columns of `X`
#' @param Z dispersion design matrix (n x r); default: intercept only
#' @param z.end list(num =, den =) rows of `Z` at the contrast endpoints; default: both equal to the
#'   intercept-only row (var and total are then 0)
#' @param bias.correct logical; subtract the model-based dispersion term `A diag(s) A'` from `M` (default TRUE)
#' @param G optional precomputed Gower matrix (`gowerCenter(D2)`)
#' @details The dispersion model is fitted without bias under arbitrary (heteroscedastic) per-sample
#'   dispersions: with `r_i = (R G R)_ii` the squared residual norm and `R` the residual projection,
#'   `E[r] = (R o R) s` exactly, so `gamma` is the least-squares solution of `r ~ (R o R) Z`. In the one-way
#'   layout this equals the classical leverage-corrected estimator `r_i / (1 - h_i)` averaged per group.
#'   `v` reports the per-sample leverage-corrected values (used for diagnostics and influence).
#' @return list with `shift`, `var`, `total`, `ratio`, `s.alt`, `s.ref`, `shift.raw` (uncorrected),
#'   `shift.norm`, `var.norm`, `total.norm`, `M` (q x q), `gamma`, `s` (model-based per-sample dispersion),
#'   `v` (leverage-corrected per-sample dispersion), `h` (leverage), `n`, `rank`
#' @export
estimatePairwiseEffects <- function(D2, X, contrast, Z = NULL, z.end = NULL, bias.correct = TRUE, G = NULL) {
  X <- as.matrix(X); n <- nrow(X)
  if (is.null(G)) G <- gowerCenter(D2)
  if (nrow(G) != n) stop("D2/G and X must have the same number of rows.")
  if (is.null(Z)) Z <- matrix(1, n, 1, dimnames = list(rownames(X), "(Intercept)"))
  Z <- as.matrix(Z)
  if (is.null(z.end)) z.end <- list(num = rep(1, ncol(Z)), den = rep(1, ncol(Z)))
  cvec <- as.numeric(contrast)
  if (length(cvec) != ncol(X)) stop("contrast length must equal ncol(X).")

  hi <- hatInfo(X); h <- hi$h
  R <- diag(n) - hi$H
  RG <- R %*% G
  r <- rowSums(RG * R)                           # diag(R G R): squared distance of i from its fitted position
  ok <- h < 1 - 1e-8
  v <- rep(NA_real_, n); v[ok] <- r[ok] / (1 - h[ok])   # leverage-corrected per-sample dispersion (diagnostic)
  M <- hi$A %*% G %*% t(hi$A)                    # = Bhat Bhat' for Euclidean data
  shift.raw <- drop(crossprod(cvec, M %*% cvec))
  # dispersion model: E[r] = (R o R) s with s = Z gamma  ->  unbiased gamma for any dispersion pattern
  gamma <- fitDispersion(R, Z, r)
  names(gamma) <- colnames(Z)
  s <- drop(Z %*% gamma); names(s) <- rownames(X)
  if (bias.correct) M <- M - hi$A %*% (s * t(hi$A))     # E[A G A'] = M + A diag(s) A'
  s.alt <- sum(as.numeric(z.end$num) * gamma); s.ref <- sum(as.numeric(z.end$den) * gamma)
  shift <- drop(crossprod(cvec, M %*% cvec))
  var   <- 2 * (s.alt - s.ref)
  total <- shift + 0.5 * var
  dimnames(M) <- list(colnames(X), colnames(X))
  names(v) <- names(h) <- rownames(X)
  list(shift = shift, var = var, total = total,
       ratio = 1 + shift / (s.alt + s.ref),
       s.alt = s.alt, s.ref = s.ref, shift.raw = shift.raw,
       shift.norm = shift / (s.alt + s.ref), var.norm = log2(s.alt / s.ref), total.norm = total / (2 * s.ref),
       M = M, gamma = gamma, s = s, v = v, h = h, n = n, rank = hi$rank)
}

# Least-squares fit of the dispersion model through E[r] = (R o R) Z gamma. Falls back to the intercept-only
# solution when the regressor matrix is rank deficient.
fitDispersion <- function(R, Z, r) {
  W <- (R * R) %*% Z
  g <- tryCatch(drop(MASS::ginv(crossprod(W)) %*% crossprod(W, r)), error = function(e) rep(NA_real_, ncol(Z)))
  g
}

# Model-implied pair cell m(a, b) = (x_a - x_b)' M (x_a - x_b) + s_a + s_b for rows of a location design
# `Xrows` (k x q) and dispersion design `Zrows` (k x r). Returns a k x k matrix.
impliedCellTable <- function(M, gamma, Xrows, Zrows) {
  Xrows <- as.matrix(Xrows); Zrows <- as.matrix(Zrows); k <- nrow(Xrows)
  s <- drop(Zrows %*% gamma)
  out <- matrix(NA_real_, k, k, dimnames = list(rownames(Xrows), rownames(Xrows)))
  for (a in seq_len(k)) for (b in seq_len(k)) {
    d <- Xrows[a, ] - Xrows[b, ]
    out[a, b] <- drop(crossprod(d, M %*% d)) + s[a] + s[b]
  }
  out
}

# ---- term-level statistics (multi-df) ------------------------------------------------------------

# Effective dimension r = tr(S)^2 / tr(S^2) of the residual covariance, estimated without bias from the
# residual Gram matrix W = R G R with nu residual degrees of freedom (Wishart moments).
effectiveDimension <- function(W, nu) {
  if (nu < 3) return(NA_real_)
  a <- sum(W * W); b <- sum(diag(W))^2
  T2 <- (a - b / nu) / ((nu - 1) * (nu + 2))
  T1sq <- (b - 2 * nu * T2) / nu^2
  r <- T1sq / T2
  if (!is.finite(r) || r < 1) r <- 1
  r
}

# Location test of the columns of Xf not in Xr (Xr nested in Xf): PERMANOVA-type F, raw / chance-corrected /
# partial R^2 and effective dimension. `G` is the Gower matrix.
termTestGower <- function(G, Xf, Xr) {
  n <- nrow(G); hf <- hatInfo(Xf); hr <- hatInfo(Xr)
  df <- hf$rank - hr$rank; nu <- n - hf$rank
  if (df <= 0 || nu <= 0) return(c(R2 = NA, R2.adj = NA, R2.partial = NA, F = NA, df = df, nu = nu, r.eff = NA, p.analytic = NA))
  Rf <- diag(n) - hf$H; W <- Rf %*% G %*% Rf
  ss <- sum((hf$H - hr$H) * G); rss <- sum(diag(W)); tss <- sum(diag(G))
  Fobs <- (ss / df) / (rss / nu)
  r <- effectiveDimension(W, nu)
  c(R2 = ss / tss,
    R2.adj = (ss - df * rss / nu) / (tss + rss / nu),   # omega^2-type; chance level removed, may be < 0
    R2.partial = ss / (ss + rss),
    F = Fobs, df = df, nu = nu, r.eff = r,
    p.analytic = if (is.finite(r)) stats::pf(Fobs, df * r, nu * r, lower.tail = FALSE) else NA_real_)
}

# Single-contrast pivotal F for shift: [a' G a / c'(X'X)^- c] / [tr(R G) / (n - q)], a = X (X'X)^- c.
# Equals termTestGower()$F for the 1-df term X vs X %*% contrastNullBasis(c).
contrastF <- function(G, hi, cvec) {
  n <- nrow(G)
  a <- drop(t(hi$A) %*% cvec)
  ss.c <- drop(crossprod(a, G %*% a)) / drop(crossprod(cvec, hi$XtXi %*% cvec))
  rss <- sum(diag(G)) - sum(hi$H * G)
  (ss.c) / (rss / (n - hi$rank))
}

# Dispersion association of a term: per-sample leverage-corrected dispersion v_i from the location model
# Xloc, variance-stabilized as sqrt(v), regressed on Zf vs the nested Zr; parametric F (used as a screen
# statistic; permutation p-values are produced by the inference layer).
dispersionTermTest <- function(G, Xloc, Zf, Zr) {
  n <- nrow(G); hi <- hatInfo(Xloc); R <- diag(n) - hi$H
  r <- rowSums((R %*% G) * R); ok <- hi$h < 0.99
  v <- sqrt(pmax(r[ok], 0) / (1 - hi$h[ok])); m <- sum(ok)
  Zf <- as.matrix(Zf)[ok, , drop = FALSE]; Zr <- as.matrix(Zr)[ok, , drop = FALSE]
  rssOf <- function(Z) { H <- hatInfo(Z)$H; sum(((diag(m) - H) %*% v)^2) }
  df <- qr(Zf)$rank - qr(Zr)$rank; nu <- m - qr(Zf)$rank
  if (df <= 0 || nu <= 0) return(c(F.disp = NA_real_, p.disp = NA_real_, df.disp = df, nu.disp = nu))
  rss.f <- rssOf(Zf); rss.r <- rssOf(Zr)
  Fd <- ((rss.r - rss.f) / df) / (rss.f / nu)
  c(F.disp = Fd, p.disp = stats::pf(Fd, df, nu, lower.tail = FALSE), df.disp = df, nu.disp = nu)
}

# Partial R^2 per term of a location design: tr((H_full - H_red) G) / tr((I - H_red) G), where H_red
# drops the term's columns. `assign` maps columns to term indices (attr "assign" of model.matrix).
partialR2PerTerm <- function(G, X, assign, term.labels) {
  X <- as.matrix(X); hf <- hatInfo(X)
  out <- setNames(rep(NA_real_, length(term.labels)), term.labels)
  for (k in seq_along(term.labels)) {
    keep <- assign != k
    hr <- hatInfo(X[, keep, drop = FALSE])
    num <- sum((hf$H - hr$H) * G); den <- sum(diag(G)) - sum(hr$H * G)
    out[k] <- if (den > 0) num / den else NA_real_
  }
  out
}

# ---- leave-one-sample-out ------------------------------------------------------------------------

# Refit without each sample; returns an n x 3 matrix of (shift, var, total) and jackknife SEs.
looPairwiseEffects <- function(G, X, contrast, Z, z.end, bias.correct = TRUE, min.n = 4) {
  n <- nrow(X)
  eff <- matrix(NA_real_, n, 3, dimnames = list(rownames(X), c("shift", "var", "total")))
  if (n <= min.n) return(list(effects = eff, se = c(shift = NA, var = NA, total = NA)))
  for (i in seq_len(n)) {
    Xi <- X[-i, , drop = FALSE]
    keep <- colSums(abs(Xi)) > 0
    if (!isEstimable(Xi[, keep, drop = FALSE], contrast[keep]) || any(abs(contrast[!keep]) > 1e-12)) next
    Zi <- Z[-i, , drop = FALSE]
    e <- tryCatch(estimatePairwiseEffects(NULL, Xi[, keep, drop = FALSE], contrast[keep], Zi, z.end, bias.correct,
                                          G = gowerCenter(uncenterGower(G)[-i, -i])),
                  error = function(err) NULL)
    if (!is.null(e)) eff[i, ] <- c(e$shift, e$var, e$total)
  }
  m <- colSums(!is.na(eff))
  se <- sqrt((m - 1) / m * colSums(sweep(eff, 2, colMeans(eff, na.rm = TRUE))^2, na.rm = TRUE))
  se[m < 2] <- NA_real_
  list(effects = eff, se = se)
}

# squared-distance matrix from a Gower matrix: D2_ij = G_ii + G_jj - 2 G_ij
uncenterGower <- function(G) {
  g <- diag(G)
  D2 <- outer(g, rep(1, length(g))) + outer(rep(1, length(g)), g) - 2 * G
  D2[D2 < 0] <- 0; diag(D2) <- 0
  dimnames(D2) <- dimnames(G)
  D2
}

# ---- per-cell-type sample distance matrices --------------------------------------------------------

#' Sample-sample distance matrices per cell type
#'
#' @param cm.per.type list (per cell type) of normalized sample x gene matrices (samples lacking the cell
#'   type may be absent or all-NA; they are dropped)
#' @param dist `"cor"` (gene-centred cosine), `"l2"` or `"l1"`
#' @param n.pcs optional number of principal components to reduce each matrix to before computing distances
#' @param center.genes logical; centre every gene across the samples present for the cell type before the
#'   cosine distance (default TRUE; this is what makes `"cor"` a squared-Euclidean-type distance on unit
#'   profiles). Ignored for `"l1"`/`"l2"`.
#' @return named list of square distance matrices (samples present for that cell type)
#' @export
sampleDistanceMatrices <- function(cm.per.type, dist = c("cor", "l2", "l1"), n.pcs = NULL, center.genes = TRUE) {
  dist <- match.arg(dist)
  lapply(cm.per.type, function(M) {
    if (is.null(M)) return(NULL)
    M <- as.matrix(M)
    M <- M[rowSums(!is.na(M)) > 0 & !apply(M, 1, anyNA), , drop = FALSE]
    n <- nrow(M)
    if (n < 2) return(matrix(0, n, n, dimnames = list(rownames(M), rownames(M))))
    if (!is.null(n.pcs) && n.pcs > 0) {
      k <- min(n.pcs, n - 1, ncol(M))
      pca <- stats::prcomp(M, center = TRUE, scale. = FALSE, rank. = k)
      M <- pca$x[, seq_len(k), drop = FALSE]
    }
    if (dist == "cor") {
      Xc <- if (center.genes) sweep(M, 2, colMeans(M), "-") else M
      nrm <- sqrt(rowSums(Xc^2)); nrm[nrm < 1e-12] <- NA
      U <- Xc / nrm
      D <- 1 - tcrossprod(U); diag(D) <- 0
      D[D < 0 & D > -1e-12] <- 0
    } else {
      D <- as.matrix(stats::dist(M, method = if (dist == "l2") "euclidean" else "manhattan"))
    }
    dimnames(D) <- list(rownames(M), rownames(M))
    D
  })
}

# ---- design-level wrapper -------------------------------------------------------------------------

# Default dispersion formula (D1): main effects of the contrast variables; `~ 1` for coefficient contrasts.
defaultDispersionFormula <- function(spec) {
  if (is.null(spec) || inherits(spec, "try-error")) return(~ 1)
  vars <- switch(spec$type,
                 simple   = parseTermVars(spec$term),
                 marginal = spec$term,
                 lincomb  = parseTermVars(spec$term),
                 coef     = character(0))
  if (!length(vars)) return(~ 1)
  stats::as.formula(paste("~", paste(unique(vars), collapse = " + ")))
}

#' Pairwise effects for one cell type from a sample-level design
#'
#' Subsets the design to the samples present in `D`, drops empty columns, checks that the contrast is
#' estimable, builds the dispersion design from `dispersion.formula` and evaluates it at the contrast
#' endpoints, then calls [estimatePairwiseEffects()].
#'
#' @param D distance matrix over (a subset of) the samples in `design`
#' @param design output of `buildDesignMatrices()` (sample level)
#' @param meta sample metadata used to build `design` (rows named by sample)
#' @param dispersion.formula formula for the dispersion model (default: main effects of the contrast
#'   variables)
#' @param dist distance type of `D` (see [asSquaredDistance()])
#' @param bias.correct see [estimatePairwiseEffects()]
#' @param influence logical; also compute leave-one-sample-out effects and jackknife SEs
#' @param min.samp.per.level minimum number of samples at each contrasted level (for factor contrasts)
#' @return list with `ok` (logical), `reason` (when skipped), the elements of [estimatePairwiseEffects()],
#'   `G`, `X`, `Z`, `contrast`, `z.end`, `samples`, `r2` (partial R^2 per location term), `cell.table`
#'   (model-implied table over the contrasted factor's levels, when applicable), `n.ref`, `n.alt`,
#'   `influence`
#' @export
pairwiseEffectsFromDesign <- function(D, design, meta, dispersion.formula = NULL, dist = c("cor", "l2", "l1"),
                                      bias.correct = TRUE, influence = FALSE, min.samp.per.level = 3) {
  dist <- match.arg(dist)
  skip <- function(reason) list(ok = FALSE, reason = reason)
  samples <- intersect(rownames(design$F), rownames(D))
  if (length(samples) < 3) return(skip("fewer than 3 samples with this cell type"))
  D2 <- asSquaredDistance(D[samples, samples], dist)
  G <- gowerCenter(D2)

  F <- design$F[samples, , drop = FALSE]
  cF <- design$contrast.F[colnames(F)]
  keep <- colSums(abs(F)) > 1e-12
  if (any(abs(cF[!keep]) > 1e-12)) return(skip("contrast involves a level absent in this cell type"))
  X <- F[, keep, drop = FALSE]; cvec <- cF[keep]
  if (!isEstimable(X, cvec)) return(skip("contrast not estimable for this cell type's design"))
  if (ncol(X) >= nrow(X)) return(skip("no residual degrees of freedom"))

  # samples per contrasted level (factor contrasts only)
  spec <- design$contrast_spec; n.ref <- n.alt <- NA_integer_
  if (!is.null(spec) && spec$type %in% c("simple", "marginal") && !grepl(":", spec$term, fixed = TRUE)) {
    g <- as.character(meta[samples, spec$term])
    n.ref <- sum(g == spec$den); n.alt <- sum(g == spec$num)
    if (min(n.ref, n.alt) < min.samp.per.level)
      return(skip(sprintf("fewer than %d samples in a contrasted level (%s: %d, %s: %d)", min.samp.per.level,
                          spec$den, n.ref, spec$num, n.alt)))
  }

  # dispersion design and its endpoint rows
  if (is.null(dispersion.formula)) dispersion.formula <- defaultDispersionFormula(spec)
  dispersion.formula <- stats::as.formula(checkFormula(dispersion.formula))
  Zfull <- buildFullDesign(dispersion.formula, meta[samples, , drop = FALSE])
  zkeep <- colSums(abs(Zfull)) > 1e-12
  Z <- Zfull[, zkeep, drop = FALSE]
  z.end <- NULL
  if (!is.null(design$contrast_endpoints_at)) {
    ze <- endpointRowsFromDesign(Zfull, meta[samples, , drop = FALSE], design$contrast_endpoints_at,
                                 numericRef = design$numeric_ref_used)
    z.end <- list(num = ze$num[zkeep], den = ze$den[zkeep])
    if (!isEstimable(Z, z.end$num - z.end$den) && any(abs(z.end$num - z.end$den) > 1e-12))
      return(skip("dispersion endpoint difference not estimable"))
  } else {
    z.end <- list(num = rep(NA_real_, ncol(Z)), den = rep(NA_real_, ncol(Z)))   # var/total undefined
  }

  eff <- estimatePairwiseEffects(NULL, X, cvec, Z, z.end, bias.correct = bias.correct, G = G)
  eff$ok <- TRUE; eff$reason <- NULL
  eff$G <- G; eff$X <- X; eff$Z <- Z; eff$contrast <- cvec; eff$z.end <- z.end; eff$samples <- samples
  eff$n.ref <- n.ref; eff$n.alt <- n.alt
  eff$dispersion.formula <- dispersion.formula
  eff$F.shift <- contrastF(G, hatInfo(X), cvec)

  # partial R^2 per location term
  assign <- attr(design$F, "assign")
  if (!is.null(assign)) {
    tl <- attr(attr(design$F, "terms"), "term.labels")
    a <- assign[keep]
    if (length(tl)) eff$r2 <- partialR2PerTerm(G, X, a, tl)
  }
  # model-implied cell table over the levels of a single contrasted factor
  if (!is.null(spec) && spec$type %in% c("simple", "marginal") && !grepl(":", spec$term, fixed = TRUE)) {
    levs <- levels(droplevels(factor(meta[samples, spec$term])))
    at0 <- if (spec$type == "simple") spec$at else list()
    one <- function(Fm, at) .oneRowFromFormula(attr(Fm, "terms"), attr(Fm, "xlevels"), attr(Fm, "contrasts"),
                                               colnames(Fm), meta[samples, , drop = FALSE], at, design$numeric_ref_used)
    Xrows <- t(sapply(levs, function(l) one(design$F, utils::modifyList(at0, setNames(list(l), spec$term)))[colnames(F)][keep]))
    Zrows <- t(sapply(levs, function(l) one(Zfull, utils::modifyList(at0, setNames(list(l), spec$term)))[zkeep]))
    if (length(levs) == 1) { Xrows <- matrix(Xrows, 1); Zrows <- matrix(Zrows, 1) }
    rownames(Xrows) <- rownames(Zrows) <- levs
    eff$cell.table <- impliedCellTable(eff$M, eff$gamma, Xrows, Zrows)
  }
  if (influence) eff$influence <- looPairwiseEffects(G, X, cvec, Z, z.end, bias.correct)
  eff
}
