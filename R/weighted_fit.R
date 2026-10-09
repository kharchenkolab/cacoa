## Weighted individual-level model on the Gower matrix (step 4b of the engine convergence): sample weights carry
## robust fits (Huber / winsorized down-weighting of samples with outlying residual distances) and weak imputation
## of samples absent from a cell type. These R functions are the reference for the C++ kernels and are used by the
## observed fits; permutation statistics come from the kernels.

# weighted hat information: A = (X'WX)^- X'W (q x n), H = X A (W-orthogonal projection, not symmetric)
hatInfoW <- function(X, w) {
  X <- as.matrix(X)
  XtWXi <- MASS::ginv(crossprod(X, w * X))
  A <- XtWXi %*% t(w * X)
  H <- X %*% A
  list(XtWXi = XtWXi, A = A, H = H, h = pmin(diag(H), 1), rank = qr(X)$rank)
}

# weighted least squares of the dispersion model E[r] = (R o R) Z gamma
fitDispersionW <- function(R, Z, r, w) {
  Wm <- (R * R) %*% Z; sw <- sqrt(w)
  tryCatch(drop(MASS::ginv(crossprod(sw * Wm)) %*% crossprod(sw * Wm, sw * r)), error = function(e) rep(NA_real_, ncol(Z)))
}

# robust sample weights from the leverage-corrected residual sizes (score_i = sqrt(r_i / (1 - h_i))): samples whose
# score exceeds the median by more than k robust standard deviations are down-weighted; `base` carries the
# weak-imputation weights and defines the reference set of the scale
robustWeights <- function(score, base, robust = c("huber", "winsor"), k = 1.345) {
  robust <- match.arg(robust)
  ref <- score[base >= 0.5 & is.finite(score)]
  med <- stats::median(ref); s <- 1.4826 * stats::median(abs(ref - med))
  if (!is.finite(s) || s <= 1e-12) return(base)
  f <- rep(1, length(score)); big <- is.finite(score) & score > med + k * s
  f[big] <- if (robust == "huber") (k * s) / (score[big] - med) else ((med + k * s) / score[big])^2
  base * f
}

#' Weighted contrast statistics of the individual-level model (reference implementation)
#'
#' @param G Gower matrix; `X` location design; `Z` dispersion design; `cvec` contrast; `z.end` endpoint rows of `Z`
#' @param w sample weights (1 = full weight; near zero for weakly imputed samples)
#' @param bias.correct bias-correct the shift
#' @param robust `"none"`, `"huber"` (iterated) or `"winsor"` (one step)
#' @param k robust tuning constant (in robust standard deviations of the residual sizes)
#' @param maxit Huber iterations
#' @return list: `F`, `shift`, `var`, `total`, `w` (final weights), `gamma`, `s`, `score`, `df.res`, `h`
#' @keywords internal
weightedContrastStats <- function(G, X, Z, cvec, z.end, w = rep(1, nrow(X)), bias.correct = TRUE, robust = c("none", "huber", "winsor"), k = 1.345, maxit = 5) {
  robust <- match.arg(robust); n <- nrow(X); base <- w
  one <- function(w) {
    hi <- hatInfoW(X, w); R <- diag(n) - hi$H; r <- rowSums((R %*% G) * R)
    ok <- hi$h < 1 - 1e-8; score <- rep(NA_real_, n); score[ok] <- sqrt(pmax(r[ok], 0) / (1 - hi$h[ok]))
    list(hi = hi, R = R, r = r, score = score)
  }
  f <- one(w)
  if (robust != "none") for (it in seq_len(if (robust == "huber") maxit else 1)) {
    w.new <- robustWeights(f$score, base, robust, k)
    if (max(abs(w.new - w)) < 1e-8) break
    w <- w.new; f <- one(w)
  }
  hi <- f$hi; R <- f$R; r <- f$r
  a <- drop(t(hi$A) %*% cvec); cXc <- drop(crossprod(cvec, hi$XtWXi %*% cvec))
  ss <- drop(crossprod(a, G %*% a)) / cXc; rss <- sum(w * r); df.res <- sum(w) - hi$rank
  Fs <- ss / (rss / df.res)
  M <- hi$A %*% G %*% t(hi$A); shift.raw <- drop(crossprod(cvec, M %*% cvec))
  gamma <- fitDispersionW(R, Z, r, w); s <- drop(Z %*% gamma)
  if (bias.correct) M <- M - hi$A %*% (s * t(hi$A))
  shift <- drop(crossprod(cvec, M %*% cvec))
  s.alt <- sum(as.numeric(z.end$num) * gamma); s.ref <- sum(as.numeric(z.end$den) * gamma)
  var <- 2 * (s.alt - s.ref)
  list(F = Fs, shift = shift, var = var, total = shift + var / 2, w = w, gamma = gamma, s = s, score = f$score, df.res = df.res, h = hi$h,
       s.alt = s.alt, s.ref = s.ref, M = M, shift.raw = shift.raw, rank = hi$rank, v = f$score^2)
}

#' Weighted term statistics (location F, dispersion F) of the individual-level model (reference implementation)
#' @inheritParams weightedContrastStats
#' @param Xf,Xr full and reduced location designs; `Zf`,`Zr` full and reduced dispersion designs
#' @return list: `F`, `F.disp`, `w`, `score`, `df`, `nu`
#' @keywords internal
weightedTermStats <- function(G, Xf, Xr, Zf, Zr, w = rep(1, nrow(G)), robust = c("none", "huber", "winsor"), k = 1.345, maxit = 5) {
  robust <- match.arg(robust); n <- nrow(G); base <- w
  one <- function(w) {
    hf <- hatInfoW(Xf, w); hr <- hatInfoW(Xr, w); R <- diag(n) - hf$H; r <- rowSums((R %*% G) * R)
    ok <- hf$h < 0.99; score <- rep(NA_real_, n); score[ok] <- sqrt(pmax(r[ok], 0) / (1 - hf$h[ok]))
    list(hf = hf, hr = hr, r = r, score = score, ok = ok)
  }
  f <- one(w)
  if (robust != "none") for (it in seq_len(if (robust == "huber") maxit else 1)) {
    w.new <- robustWeights(f$score, base, robust, k)
    if (max(abs(w.new - w)) < 1e-8) break
    w <- w.new; f <- one(w)
  }
  ss <- sum((w * (f$hf$H - f$hr$H)) * G); rss <- sum(w * f$r)
  df <- f$hf$rank - f$hr$rank; nu <- sum(w) - f$hf$rank
  Fl <- if (df > 0 && nu > 0) (ss / df) / (rss / nu) else NA_real_
  ok <- f$ok; v <- f$score[ok]; ww <- w[ok]
  rssOf <- function(Z) { Z <- as.matrix(Z)[ok, , drop = FALSE]; b <- MASS::ginv(crossprod(Z, ww * Z)) %*% crossprod(Z, ww * v); sum(ww * (v - Z %*% b)^2) }
  dfd <- qr(as.matrix(Zf))$rank - qr(as.matrix(Zr))$rank; nud <- sum(ww) - qr(as.matrix(Zf))$rank
  Fd <- if (dfd > 0 && nud > 0) ((rssOf(Zr) - rssOf(Zf)) / dfd) / (rssOf(Zf) / nud) else NA_real_
  list(F = Fl, F.disp = Fd, w = w, score = f$score, df = df, nu = nu)
}

# weak imputation: samples absent from a unit are kept with the mean squared distance of the present samples and a
# near-zero weight, so the design keeps its rows (and its levels) while the estimates are driven by the present ones
imputeAbsentSamples <- function(D2, present, all.samples, weak.weight = 1e-4) {
  n <- length(all.samples); out <- matrix(NA_real_, n, n, dimnames = list(all.samples, all.samples))
  out[present, present] <- D2[present, present]
  absent <- setdiff(all.samples, present)
  if (length(absent)) {
    fill <- mean(D2[present, present][upper.tri(D2[present, present])])
    out[absent, ] <- fill; out[, absent] <- fill; diag(out) <- 0
  }
  w <- stats::setNames(rep(1, n), all.samples); w[absent] <- weak.weight
  list(D2 = out, w = w)
}

# Freedman-Lane kernel parts under weights: the reduced model is fitted with the weights (H_r = X_r (X_r'WX_r)^- X_r'W);
# G* = K1 + K2[, p] + K3[p, ] + K4[p, p] re-indexes the residual part under the relabeling p
flGowerPartsW <- function(G, Xr, w) {
  n <- nrow(G); Hr <- if (ncol(Xr)) hatInfoW(Xr, w)$H else matrix(0, n, n); Rr <- diag(n) - Hr
  list(K1 = Hr %*% G %*% t(Hr), K2 = Hr %*% G %*% t(Rr), K3 = Rr %*% G %*% t(Hr), K4 = Rr %*% G %*% t(Rr))
}
