## Tail approximations of per-cell permutation p-values, evaluated on held-out null draws.
## For every cell the first B.small permutations are the "available" null sample (as in a routine run) and a random set
## of the remaining permutations are held-out null draws of the same statistic; an approximation is calibrated when the
## pooled held-out p-values are uniform, and its accuracy is measured against the rank of each draw among all B.
## Usage: Rscript tail_analyze.R <collected.rds> <out.prefix> [n.heldout] [n.cores]
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
res <- readRDS(args[1]); out <- args[2]; n.held <- if (length(args) > 2) as.integer(args[3]) else 4000L; n.cores <- if (length(args) > 3) as.integer(args[4]) else 32L
B.small.set <- as.integer(strsplit(Sys.getenv("CACOA_BSMALL", "199,499,1999"), ",")[[1]])   # training sizes (permutations available to a routine run)
alphas <- c(0.05, 0.01, 1e-3, 1e-4, 1e-5)
brk <- seq(0, 12, by = 0.05)   # histogram of -log10 p for the QQ curves
set.seed(7)

# generalized Pareto tail by probability-weighted moments (Hosking & Wallis 1987); returns the survival function of the
# exceedances y >= 0 (xi > 0 heavy, xi < 0 bounded)
gpdPWM <- function(y) {
  y <- sort(y); n <- length(y)
  b0 <- mean(y); b1 <- sum((seq_len(n) - 1) / (n - 1) * y) / n
  k <- b0 / b1 - 2; alpha <- b0 * (b0 / b1 - 1)   # from b0 = alpha / (1 + k), b1 = alpha / ((1 + k) (2 + k))
  xi <- -k; sigma <- alpha
  function(x) { s <- 1 + xi * x / sigma; s[s < 0] <- 0; if (abs(xi) < 1e-8) exp(-x / sigma) else s^(-1 / xi) }
}
# generalized Pareto tail by maximum likelihood (xi, sigma), exceedances y >= 0
gpdMLE <- function(y) {
  nll <- function(par) { xi <- par[1]; sg <- exp(par[2]); s <- 1 + xi * y / sg; if (any(s <= 0)) return(1e10); length(y) * log(sg) + (1 + 1 / xi) * sum(log(s)) }
  o <- try(stats::optim(c(0.1, log(mean(y))), nll, method = "Nelder-Mead", control = list(maxit = 500)), silent = TRUE)
  if (inherits(o, "try-error") || !is.finite(o$value) || o$value >= 1e10) return(NULL)
  xi <- o$par[1]; sigma <- exp(o$par[2])
  function(x) { s <- 1 + xi * x / sigma; s[s < 0] <- 0; if (abs(xi) < 1e-8) exp(-x / sigma) else s^(-1 / xi) }
}
# number of elements of v that are >= each x (sorted lookup)
countGE <- function(v, x) { sv <- sort(v); length(sv) - findInterval(x - 1e-12, sv) }
# per-cell approximations on the non-negative scale S; signed.* give the signed statistic for the two-sided normal tail
approxP <- function(train, held, signed.train = NULL, signed.held = NULL, pool = NULL, pool.loo = NULL, pool.iqr = NULL, pool.shrink = NULL, sd.shrink = NULL) {
  train <- train[is.finite(train)]; n <- length(train)
  out <- list()
  out$exact <- (countGE(train, held) + 1) / (n + 1)
  mu <- mean(train); v <- var(train)
  st <- if (!is.null(signed.train)) signed.train[is.finite(signed.train)] else train; sh <- if (!is.null(signed.train)) signed.held else held
  zs <- (sh - mean(st)) / sd(st); if (!is.null(signed.train)) zs <- abs(zs)
  out$pooled <- (countGE(pool, zs) + 1) / (length(pool) + 1)
  out$pooled.loo <- (countGE(pool.loo, zs) + 1) / (length(pool.loo) + 1)
  zi <- (sh - median(st)) / stats::IQR(st); if (!is.null(signed.train)) zi <- abs(zi)   # robust location and scale
  out$pooled.iqr <- (countGE(pool.iqr, zi) + 1) / (length(pool.iqr) + 1)
  zk <- (sh - mean(st)) / sd.shrink; if (!is.null(signed.train)) zk <- abs(zk)          # empirical-Bayes shrunk sd
  out$pooled.shrink <- (countGE(pool.shrink, zk) + 1) / (length(pool.shrink) + 1)
  out$gamma <- stats::pgamma(held, shape = mu^2 / v, rate = mu / v, lower.tail = FALSE)
  if (!is.null(signed.train)) { z <- (signed.held - mean(signed.train)) / sd(signed.train); out$normal <- 2 * stats::pnorm(-abs(z)) }
  else { z <- (held - mu) / sqrt(v); out$normal <- stats::pnorm(z, lower.tail = FALSE) }
  out
}
evalTest <- function(name, perm, signed = FALSE, strata = NULL) {
  keep <- colSums(is.finite(perm)) > 0.9 * nrow(perm)
  perm <- perm[, keep, drop = FALSE]; if (!is.null(strata)) strata <- strata[keep]
  B <- nrow(perm); m <- ncol(perm)
  held.rows <- sort(sample(seq(max(B.small.set) + 1L, B), n.held))
  S <- if (signed) perm^2 else perm
  # oracle: standardization by the moments of all non-held permutations, strata by the kurtosis of those
  full <- if (signed) perm[-held.rows, , drop = FALSE] else S[-held.rows, , drop = FALSE]
  mu.full <- colMeans(full, na.rm = TRUE); sd.full <- apply(full, 2, sd, na.rm = TRUE)
  zfull <- scale(full, center = mu.full, scale = sd.full); if (signed) zfull <- abs(zfull)
  kurt.full <- apply(zfull, 2, function(v) { v <- v[is.finite(v)]; mean(v^4) - 3 })
  kbin.full <- cut(kurt.full, unique(quantile(kurt.full, seq(0, 1, 0.1), na.rm = TRUE)), include.lowest = TRUE)
  oracle <- list(all = sort(as.numeric(zfull[is.finite(zfull)])))
  # combined: log S (shift) or |t| (density), leave-one-out standardized, stratified by the data covariate
  # pooled standardized training statistics per B.small: |z| of the signed statistic, or z of S
  pools <- lapply(B.small.set, function(Bs) {
    tr <- if (signed) perm[seq_len(Bs), , drop = FALSE] else S[seq_len(Bs), , drop = FALSE]
    z <- scale(tr); if (signed) z <- abs(z)
    # robust standardization: median and IQR
    med <- apply(tr, 2, median, na.rm = TRUE); iqr <- apply(tr, 2, stats::IQR, na.rm = TRUE)
    zi <- scale(tr, center = med, scale = iqr); if (signed) zi <- abs(zi)
    # empirical-Bayes shrinkage of the per-cell variance of the permutation statistic (limma's squeezeVar) with the
    # effective degrees of freedom of a heavy-tailed variance estimate, 2 B / (kurtosis + 2), kurtosis pooled by stratum
    v.tr <- apply(tr, 2, var, na.rm = TRUE)
    kurt.pool <- if (is.null(strata)) rep(mean(z^4, na.rm = TRUE) - 3, m) else { kk <- tapply(seq_len(m), strata, function(idx) mean(z[, idx]^4, na.rm = TRUE) - 3); as.numeric(kk[as.character(strata)]) }
    df.eff <- pmax(2 * Bs / (pmax(kurt.pool, 0) + 2), 4)
    sq <- if (is.null(strata)) limma::squeezeVar(v.tr, df = df.eff) else { vp <- v.tr; for (lv in unique(strata)) { idx <- which(strata == lv); if (length(idx) > 5) vp[idx] <- limma::squeezeVar(v.tr[idx], df = df.eff[idx])$var.post }; list(var.post = vp) }
    sd.shr <- sqrt(sq$var.post)
    zk <- scale(tr, center = colMeans(tr, na.rm = TRUE), scale = sd.shr); if (signed) zk <- abs(zk)
    # leave-one-out standardization of the training values: each value against the mean and sd of the other B - 1
    mu <- colMeans(tr, na.rm = TRUE); ss <- colSums(tr^2, na.rm = TRUE); nB <- colSums(is.finite(tr))
    mu.loo <- sweep(-tr, 2, nB * mu, "+"); mu.loo <- sweep(mu.loo, 2, nB - 1, "/")
    var.loo <- sweep(-tr^2, 2, ss, "+"); var.loo <- sweep(var.loo, 2, nB - 1, "/") - mu.loo^2; var.loo <- var.loo * sweep(matrix(1, nrow(tr), ncol(tr)), 2, (nB - 1) / (nB - 2), "*")
    zloo <- (tr - mu.loo) / sqrt(var.loo); if (signed) zloo <- abs(zloo)
    list(all = sort(as.numeric(z[is.finite(z)])), loo = sort(as.numeric(zloo[is.finite(zloo)])), iqr = sort(as.numeric(zi[is.finite(zi)])),
         shrink = sort(as.numeric(zk[is.finite(zk)])), sd.shrink = sd.shr, sd.train = sqrt(v.tr))
  }); names(pools) <- B.small.set
  shape <- t(apply(S, 2, function(s) { s <- s[is.finite(s)]; c(mean = mean(s), sd = sd(s), skew = mean((s - mean(s))^3) / sd(s)^3, kurt = mean((s - mean(s))^4) / sd(s)^4 - 3) }))
  per.cell <- mclapply(seq_len(m), function(j) {
    s <- S[, j]; ok <- is.finite(s)
    all.s <- s[ok]; h <- s[held.rows]; hok <- is.finite(h); h <- h[hok]
    p.ref <- countGE(all.s, h) / length(all.s)   # the draw itself counts: >= 1/B
    lo <- p.ref <= 0.01
    cnt <- list(); hist <- list(); acc <- list()
    for (Bs in B.small.set) {
      tr <- s[seq_len(Bs)]
      pl <- pools[[as.character(Bs)]]
      a <- if (signed) approxP(tr, h, signed.train = perm[seq_len(Bs), j], signed.held = perm[held.rows, j][hok], pool = pl$all, pool.loo = pl$loo, pool.iqr = pl$iqr, pool.shrink = pl$shrink, sd.shrink = pl$sd.shrink[j])
           else approxP(tr, h, pool = pl$all, pool.loo = pl$loo, pool.iqr = pl$iqr, pool.shrink = pl$shrink, sd.shrink = pl$sd.shrink[j])
      zo <- if (signed) abs((perm[held.rows, j][hok] - mu.full[j]) / sd.full[j]) else (h - mu.full[j]) / sd.full[j]
      a$oracle <- (countGE(oracle$all, zo) + 1) / (length(oracle$all) + 1)
      for (mth in names(a)) {
        p <- a[[mth]]
        cnt[[length(cnt) + 1]] <- data.frame(B.small = Bs, method = mth, n = length(p), t(vapply(alphas, function(al) sum(p < al), numeric(1))))
        hist[[length(hist) + 1]] <- data.frame(B.small = Bs, method = mth, bin = seq_along(brk), count = tabulate(findInterval(pmin(-log10(pmax(p, 1e-300)), 11.99), brk), nbins = length(brk)))
        if (any(lo)) acc[[length(acc) + 1]] <- data.frame(B.small = Bs, method = mth, p.ref = p.ref[lo], p = p[lo])
      }
    }
    list(cnt = bind_rows(cnt), hist = bind_rows(hist), acc = bind_rows(acc))
  }, mc.cores = n.cores)
  sdiag <- bind_rows(lapply(B.small.set, function(Bs) data.frame(B.small = Bs, kbin = kbin.full, ratio.train = pools[[as.character(Bs)]]$sd.train / sd.full, ratio.shrink = pools[[as.character(Bs)]]$sd.shrink / sd.full))) %>%
    group_by(B.small, kbin) %>% summarise(n = n(), med.train = median(ratio.train), frac.train.low = mean(ratio.train < 1), med.shrink = median(ratio.shrink), frac.shrink.low = mean(ratio.shrink < 1), .groups = "drop")
  cat("\n", name, ": sd of the training permutations relative to the true sd, by kurtosis decile\n"); print(as.data.frame(sdiag), digits = 3)
  cnt <- bind_rows(lapply(per.cell, `[[`, "cnt")) %>% group_by(B.small, method) %>% summarise(across(everything(), sum), .groups = "drop")
  names(cnt)[4:8] <- paste0("n_", alphas)
  hist <- bind_rows(lapply(per.cell, `[[`, "hist")) %>% group_by(B.small, method, bin) %>% summarise(count = sum(count), .groups = "drop")
  acc <- bind_rows(lapply(per.cell, `[[`, "acc"))
  list(cnt = mutate(cnt, test = name), hist = mutate(hist, test = name), acc = mutate(acc, test = name),
       shape = as.data.frame(shape) %>% mutate(test = name, cell = colnames(perm)))
}
tests <- list(
  evalTest("expression shift (F)", res$shifts$perm, signed = FALSE, strata = res$n.samples[colnames(res$shifts$perm)]),
  evalTest("density: heat kernel (t)", res$heat$perm, signed = TRUE, strata = if (is.null(res$heat$ystats)) NULL else cut(res$heat$ystats[colnames(res$heat$perm), "n.eff"], unique(quantile(res$heat$ystats$n.eff, seq(0, 1, 0.1))), include.lowest = TRUE)),
  evalTest("density: 1-hop log counts (t)", res$nb$perm, signed = TRUE, strata = if (is.null(res$nb$ystats)) NULL else cut(res$nb$ystats[colnames(res$nb$perm), "n.eff"], unique(quantile(res$nb$ystats$n.eff, seq(0, 1, 0.1))), include.lowest = TRUE)))
# relation of the per-cell null shape with the data covariate
if (!is.null(res$heat$ystats)) { k.heat <- tests[[2]]$shape; k.heat$n.eff <- res$heat$ystats[k.heat$cell, "n.eff"]; k.heat$max.share <- res$heat$ystats[k.heat$cell, "max.share"]
  cat("heat kernel: Spearman cor of the null kurtosis with n.eff", cor(k.heat$kurt, k.heat$n.eff, method = "spearman"), "with max.share", cor(k.heat$kurt, k.heat$max.share, method = "spearman"), "\n")
  print(k.heat %>% mutate(bin = cut(n.eff, quantile(n.eff, seq(0, 1, 0.2)), include.lowest = TRUE)) %>% group_by(bin) %>% summarise(n = n(), kurt.med = median(kurt), kurt.q90 = quantile(kurt, .9), skew.med = median(skew))) }
k.sh <- tests[[1]]$shape; k.sh$n.samples <- res$n.samples[k.sh$cell]
print(k.sh %>% group_by(n.samples) %>% summarise(n = n(), kurt.med = median(kurt), kurt.q90 = quantile(kurt, .9), skew.med = median(skew)))
cnt <- bind_rows(lapply(tests, `[[`, "cnt")); hist <- bind_rows(lapply(tests, `[[`, "hist")); acc <- bind_rows(lapply(tests, `[[`, "acc")); shape <- bind_rows(lapply(tests, `[[`, "shape"))
calib <- cnt; for (al in alphas) calib[[paste0("ratio_", al)]] <- cnt[[paste0("n_", al)]] / (cnt$n * al)
calib <- calib %>% select(test, method, B.small, n, starts_with("ratio_"))
write.csv(calib, paste0(out, "_calibration.csv"), row.names = FALSE)
accs <- acc %>% mutate(bin = cut(p.ref, c(0, 1e-4, 1e-3, 1e-2), labels = c("<1e-4", "1e-4..1e-3", "1e-3..1e-2")), lr = log10(pmax(p, 1e-300) / p.ref)) %>%
  group_by(test, method, B.small, bin) %>% summarise(n = n(), median.log10.ratio = median(lr), q10 = quantile(lr, 0.1), q90 = quantile(lr, 0.9), .groups = "drop")
write.csv(accs, paste0(out, "_accuracy.csv"), row.names = FALSE)
write.csv(shape %>% group_by(test) %>% summarise(skew.med = median(skew), skew.q10 = quantile(skew, .1), skew.q90 = quantile(skew, .9), kurt.med = median(kurt), kurt.q90 = quantile(kurt, .9)), paste0(out, "_shape.csv"), row.names = FALSE)
# QQ from the pooled histograms: observed -log10 p (bin lower edge) against the expected quantile of the pooled sample
qq <- hist %>% group_by(test, method, B.small) %>% arrange(desc(bin), .by_group = TRUE) %>% mutate(cum = cumsum(count), N = sum(count), exp = -log10(pmax(cum, 0.5) / N), obs = brk[bin]) %>% ungroup()
g <- ggplot(qq, aes(exp, obs, colour = method)) + geom_abline(linetype = 2, colour = "grey50") + geom_step(direction = "vh") +
  facet_grid(test ~ factor(paste("B =", B.small), levels = paste("B =", sort(unique(B.small))))) + coord_cartesian(xlim = c(0, 7), ylim = c(0, 7)) + theme_bw() +
  labs(x = "expected -log10 p (uniform)", y = "observed -log10 p", title = sprintf("Held-out null draws, %s relabeling %d (%d cells, B = %d)", res$ds, res$k, length(res$cells), res$B))
ggsave(paste0(out, "_qq.png"), g, width = 11, height = 8, dpi = 110)
options(width = 200); print(as.data.frame(calib), digits = 3); print(as.data.frame(accs), digits = 3); print(read.csv(paste0(out, "_shape.csv")), digits = 3)
