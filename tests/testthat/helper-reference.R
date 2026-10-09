# R reference implementations kept only to validate the C++ kernels (engine_convergence.md §11, class 1).

# former per-cell R loop of clusterFreeExpressionShifts() (before step 4): returns stat, p, shift, n, z, max, min
refClusterFreeShifts <- function(Y, pairs, samples, design, meta, gplan, P.global, min.samp.per.level = 2, block.vars = NULL) {
  n <- length(samples); B <- ncol(P.global)
  des <- design; des$F <- design$F[samples, , drop = FALSE]
  cF <- design$contrast.F[colnames(des$F)]; spec <- design$contrast_spec
  lev.var <- if (!is.null(spec) && spec$type %in% c("simple", "marginal") && !grepl(":", spec$term, fixed = TRUE) && spec$term %in% names(meta) && !is.numeric(meta[[spec$term]])) spec$term else NULL
  m <- ncol(Y); out <- list(F = rep(NA_real_, m), p = rep(NA_real_, m), shift = rep(NA_real_, m), n = rep(0L, m), z = rep(NA_real_, m))
  mx <- rep(-Inf, B); mn <- rep(Inf, B)
  for (k in seq_len(m)) {
    y <- Y[, k]; if (all(is.na(y))) next
    D <- cacoa:::pairVectorToMatrix(y, pairs, samples); s <- rownames(D)
    if (length(s) < 3) next
    if (!is.null(lev.var)) { g <- as.character(meta[s, lev.var]); tb <- table(g)[c(spec$den, spec$num)]; if (any(is.na(tb)) || min(tb) < min.samp.per.level) next }
    X <- des$F[s, , drop = FALSE]; keep <- colSums(abs(X)) > 1e-12
    if (any(abs(cF[!keep]) > 1e-12)) next
    X <- X[, keep, drop = FALSE]; cvec <- cF[keep]
    if (ncol(X) >= nrow(X) || !cacoa:::isEstimable(X, cvec)) next
    G <- cacoa:::gowerCenter(cacoa:::asSquaredDistance(D, "cor"))
    pre <- cacoa:::inferencePrecompute(G, X, matrix(1, nrow(X), 1), cvec)
    eff <- cacoa:::estimatePairwiseEffects(NULL, X, cvec, NULL, NULL, TRUE, G = G)
    plan <- cacoa:::permutationPlan(des, meta, s, scheme = gplan$scheme, block.vars = block.vars, n.permutations = B, max.enumerate = 0)
    sub <- match(s, samples)
    P <- apply(P.global, 2, cacoa:::inducePermutation, sub = sub, sub.plan = plan); storage.mode(P) <- "integer"
    Fp <- if (gplan$scheme == "freedman-lane") { parts <- cacoa:::flGowerParts(G, X %*% cacoa:::contrastNullBasis(cvec))
      cacoa:::permuted_contrast_F_fl(parts$K1, parts$K2, parts$K3, parts$K4, pre$a, pre$H, pre$cXc, pre$q, P) } else cacoa:::permuted_contrast_F(G, pre$a, pre$H, pre$cXc, pre$q, P)
    Fp <- as.numeric(Fp); Fo <- cacoa:::contrastF(G, cacoa:::hatInfo(X), cvec)
    out$F[k] <- Fo; out$p[k] <- (sum(Fp >= Fo - 1e-12) + 1) / (B + 1); out$shift[k] <- eff$shift; out$n[k] <- length(s)
    mu <- mean(Fp); sdv <- sd(Fp)
    if (is.finite(sdv) && sdv > 0) { out$z[k] <- (Fo - mu) / sdv; zp <- (Fp - mu) / sdv; mx <- pmax(mx, zp); mn <- pmin(mn, zp) }
  }
  out$max <- mx; out$min <- mn; out
}
