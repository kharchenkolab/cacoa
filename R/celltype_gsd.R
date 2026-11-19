## Cell-type expression shifts with the per-gene statistic gsd: every cell type is tested as one neighbourhood by the
## cluster-free kernel (clusterFreeExpressionShifts(statistic = "gsd")), and the shift magnitudes, their jackknife errors and
## the per-gene z-scores are computed from the same profiles, noise model and samples.

# per-gene fit of one cell type: x, v (genes x present samples), tau2 (per gene), F (present samples x columns), cF (contrast)
gsdGeneFit <- function(x, v, tau2, F, cF, d0 = 3) {
  keep <- colSums(abs(F)) > 1e-12; F <- F[, keep, drop = FALSE]; cF <- cF[keep]
  q <- qr(F)
  if (q$rank < ncol(F)) {   # aliased nuisance columns are dropped; an aliased contrast leaves nothing to estimate
    keep <- sort(q$pivot[seq_len(q$rank)]); if (any(abs(cF[-keep]) > 0)) return(NULL)
    F <- F[, keep, drop = FALSE]; cF <- cF[keep]; q <- qr(F)
  }
  if (nrow(F) - q$rank < 1) return(NULL)
  XtXi <- chol2inv(qr.R(q))[order(q$pivot), order(q$pivot), drop = FALSE]
  a <- as.numeric(F %*% XtXi %*% cF); cc <- sum(cF * (XtXi %*% cF)); df <- nrow(F) - q$rank
  beta <- as.numeric(x %*% a)
  rss <- colSums(qr.resid(q, t(x))^2)
  m <- as.numeric((v + tau2) %*% a^2)
  s2 <- (d0 * m + cc * rss) / (d0 + df)
  list(beta = beta, se = sqrt(s2), z = beta / sqrt(s2), rss = rss, df = df, c = cc)
}

# bias-corrected magnitudes from a fit, per gene num = beta^2 - c RSS / df (unbiased for the squared true change):
#   total: sum over genes of num, the squared distance between the condition means (log10 units);
#   std: the same sum over the summed donor-level variance (RSS / df minus the mean sampling noise; the denominator is floored at
#        a tenth of the summed residual variance, so that it stays positive when the donors barely differ);
#   n.eff: breadth, the effective number of changed genes (sum b2)^2 / sum b4 with b2, b4 unbiased for the squared and fourth
#        powers of the true change (normal errors with the gene's variance se^2)
gsdMagnitudes <- function(f, vbar, tau2, breadth = TRUE) {
  if (is.null(f)) return(c(total = NA_real_, std = NA_real_, n.eff = NA_real_))
  s2 <- f$rss / f$df; num <- f$beta^2 - f$c * s2
  ne <- NA_real_
  if (breadth) {
    v <- f$se^2; b2 <- sum(f$beta^2 - v); b4 <- sum(f$beta^4 - 6 * f$beta^2 * v + 3 * v^2)
    if (b2 > 0 && b4 > 0) ne <- min(b2^2 / b4, length(num))
  }
  c(total = sum(num), std = sum(num) / max(sum(s2 - vbar), 0.1 * sum(s2)), n.eff = ne)
}

# first-order variance of the total shift (sum over genes of the squared change) from the residuals of the fit: the genes x samples
# residual matrix enters only through its samples x samples products (squaredDistanceLinearVar())
gsdLinearVar <- function(x, F, cF, beta) {
  keep <- colSums(abs(F)) > 1e-12; F <- F[, keep, drop = FALSE]; cF <- cF[keep]; q <- qr(F)
  if (q$rank < ncol(F)) { kp <- sort(q$pivot[seq_len(q$rank)]); F <- F[, kp, drop = FALSE]; cF <- cF[kp]; q <- qr(F) }
  nu <- nrow(F) - q$rank; if (nu < 3) return(NA_real_)
  XtXi <- chol2inv(qr.R(q))[order(q$pivot), order(q$pivot), drop = FALSE]; a <- as.numeric(F %*% XtXi %*% cF)
  E <- t(qr.resid(q, t(x)))                     # genes x samples
  squaredDistanceLinearVar(crossprod(E), drop(crossprod(E, beta)), a, rowSums(qr.Q(q)^2), nu)
}

# label-free winsorization of a gene's sample values (as in the kernel): median +- k max(sqrt(v + tau2), MAD sd)
gsdWinsorize <- function(x, v, tau2, k) {
  med <- apply(x, 1, stats::median); mad <- 1.4826 * apply(abs(x - med), 1, stats::median)
  sd <- pmax(sqrt(v + tau2), mad)
  pmin(pmax(x, med - k * sd), med + k * sd)
}

#' Cell-type expression shifts with the per-gene statistic gsd
#'
#' Each cell type is tested as one group of cells with the statistic of the cluster-free map (see
#' [clusterFreeExpressionShifts()], `statistic = "gsd"`): per gene the contrast of the profile values over the samples,
#' standardized by a moderated standard error (sampling-noise model and residual variance), and the z-scores of the
#' genes combined adaptively over powers. P-values are exact permutation p-values (relabelings of the samples), adjusted
#' by BH across cell types (`padj`) and by the max-statistic across cell types (`p.fwer`).
#'
#' The magnitude of the shift is reported separately from its significance, in two bias-corrected forms that do not depend on
#' the numbers of cells (the expected contribution of sampling and donor noise is removed): `estimate.norm`, the standardized
#' shift, is the summed squared change of the genes relative to their summed donor-level variance (the change in units of
#' within-group variation, comparable across cell types); `estimate`, the total shift, is the summed squared change in log10
#' expression (the squared distance between the condition means). `n.eff` is the breadth of the change: the effective number of
#' changed genes. All are fitted with the test's design plus the permutation strata (`block.vars`) as nuisance factors; the
#' magnitudes can be slightly negative when there is no change.
#'
#' @param cm gene x cell count matrix (sparse) of the tested genes
#' @param sample.per.cell named factor of sample per cell
#' @param cell.groups named factor of cell type per cell
#' @param model a `cacoaModel` (see [buildCacoaModel()]); whole-factor (term) tests are skipped (use the distance test)
#' @param lib library size per cell over all genes (named by cell)
#' @param min.cells.per.sample minimum cells of a sample in a cell type for the sample to enter its test (default 5)
#' @param min.samp.per.level minimum samples per compared level (default 3)
#' @param neighbourhood.gene.frac test, per cell type, the genes detected in at least this fraction of its cells (default 0.05)
#' @param permutation,n.permutations,block.vars,seed permutation settings (see [permutationPlan()])
#' @param robust `"none"` or `"winsor"` (label-free winsorization of each gene's sample values at `winsor.k` robust sd)
#' @param winsor.k winsorization threshold (default 3)
#' @param gsd.d0 prior degrees of freedom of the moderated variance (default 3)
#' @param influence jackknife standard errors and leave-one-sample-out changes of the magnitude (default TRUE)
#' @param n.null relabelings at which the magnitudes are also computed, for the null band of the plots (default 99; 0 skips)
#' @param n.cores threads
#' @param verbose progress
#' @param return.profiles also return, per test, the profiles the test was computed from (`profiles`: per cell type the genes x samples
#'   matrix `x` of normalized expression, its sampling variance `v`, the gene-wise extra variance `tau2`, the samples and genes; the
#'   design `F` and the contrast), for evaluations of the statistic outside the compiled kernel (default FALSE)
#' @return list: `results` (long table as for the distance test, `effect = "shift"`, plus `n.genes`, `n.eff`, `se.jk.norm`),
#'   `global` (per test: max-statistic p across cell types), `genes` (per test, per cell type: data.frame of `gene`,
#'   `beta` (change in log10 expression), `se`, `z`, ordered by |z|), `influence` (per test: list with the samples x cell
#'   types matrices `shift` (standardized) and `shift.total`: the leave-one-out change), `null.magnitudes` (per test: relabelings x
#'   cell types matrices `std` and `total`), `skipped`, `noise.curves`, `model`, `settings`
#' @export
gsdShiftsPerCellType <- function(cm, sample.per.cell, cell.groups, model, lib, min.cells.per.sample = 5, min.samp.per.level = 2,
                                 neighbourhood.gene.frac = 0.05, permutation = "auto", n.permutations = 999, block.vars = model$block.vars,
                                 seed = NULL, robust = c("none", "winsor"), winsor.k = 3, gsd.d0 = 3, influence = TRUE, n.null = 99, n.cores = 1, verbose = FALSE, return.profiles = FALSE) {
  robust <- match.arg(robust); prof.out <- list()
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  meta <- model$meta
  samples <- intersect(rownames(meta), levels(droplevels(as.factor(sample.per.cell))))
  cells <- colnames(cm)[as.character(sample.per.cell[colnames(cm)]) %in% samples & !is.na(cell.groups[colnames(cm)])]
  cm <- cm[, cells, drop = FALSE]; lib <- as.numeric(lib[cells])
  spc <- setNames(factor(as.character(sample.per.cell[cells]), levels = samples), cells)
  cg <- droplevels(as.factor(cell.groups[cells]))
  units <- split(seq_along(cells) - 1L, cg); units <- units[lengths(units) > 0]
  types <- names(units); B <- n.permutations
  scale <- stats::median(lib[lib > 0])
  skipped <- data.frame(test = character(0), celltype = character(0), reason = character(0))
  out <- list(); genes.out <- list(); infl.out <- list(); null.out <- list(); global <- list(); curves <- NULL
  for (t in model$tests) {
    if (t$kind != "contrast") {
      skipped <- rbind(skipped, data.frame(test = t$label, celltype = types, reason = "whole-factor test: use statistic = 'distance'"))
      next
    }
    if (verbose) message(sprintf("Testing %d cell types for '%s' (gsd, %d permutations)...", length(units), t$label, B))
    r <- clusterFreeExpressionShifts(cm, spc, units, t$design, meta, statistic = "gsd", lib = lib, profile.scale = scale, noise.curves = curves,
                                     n.fit.units = length(units), min.n.obs.per.samp = min.cells.per.sample, min.samp.per.level = min.samp.per.level,
                                     permutation = permutation, n.permutations = B, block.vars = block.vars, seed = seed, adjust = FALSE, smooth = FALSE,
                                     post = "none", return.perm = TRUE, n.cores = n.cores, neighbourhood.gene.frac = neighbourhood.gene.frac,
                                     robust = robust, winsor.k = winsor.k, gsd.d0 = gsd.d0, max.enumerate = B)
    curves <- r$noise.curves; ex <- isTRUE(r$settings$exhaustive)
    stat <- setNames(as.numeric(r$stat), types); p <- setNames(as.numeric(r$p.value), types)
    tested <- is.finite(stat) & is.finite(p)
    # max-statistic across cell types on the permutation-standardized statistics
    P <- r$stats.perm; mu <- colMeans(P, na.rm = TRUE); sdv <- apply(P, 2, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv <= 0] <- NA
    Zp <- sweep(sweep(P, 2, mu, "-"), 2, sdv, "/"); zo <- (stat - mu) / sdv
    mx <- apply(Zp[, tested, drop = FALSE], 1, function(z) if (all(is.na(z))) NA else max(z, na.rm = TRUE)); mx <- mx[is.finite(mx)]
    p.fwer <- vapply(zo, function(z) if (!is.finite(z)) NA_real_ else if (ex) max(sum(mx >= z), 1) / length(mx) else (1 + sum(mx >= z)) / (1 + length(mx)), numeric(1))   # enumerated: identity included
    padj <- rep(NA_real_, length(types)); padj[tested] <- stats::p.adjust(p[tested], "BH")
    # magnitudes, jackknife and per-gene z from the same profiles
    n.s <- length(samples); pairs <- matrix(c(0L, 1L), 1)
    m1 <- shiftPoissonTable()
    opts <- list(scale = scale, dist = 2L, noise = TRUE, tau = TRUE, correct = FALSE, unit_detection = neighbourhood.gene.frac,
                 mu_grid = curves$mu.grid, phi = curves$phi, tau2 = curves$tau2, lam_grid = m1$lam.grid, s_grid = m1$s.grid, V1 = m1$V1)
    prof <- cf_unit_profiles(cm, as.integer(spc), lapply(units, as.integer), lib, n.s, as.integer(min.cells.per.sample), pairs, opts)
    if (return.profiles) prof.out[[t$label]] <- list(profiles = setNames(lapply(seq_along(types), function(k) list(x = prof[[k]]$x, v = prof[[k]]$v, tau2 = as.numeric(prof[[k]]$tau2),
                                                     samples = samples[prof[[k]]$present + 1L], genes = rownames(cm)[prof[[k]]$genes + 1L])), types), F = t$design$F[samples, , drop = FALSE], contrast = t$design$contrast.F, samples = samples)
    Fd <- t$design$F[samples, , drop = FALSE]; cF <- t$design$contrast.F[colnames(Fd)]
    # magnitudes: the test's design plus the permutation strata (block.vars) not in it, as nuisance factors; otherwise differences
    # between strata inflate the residual variance and the bias correction overcorrects
    bv <- setdiff(intersect(block.vars, names(meta)), all.vars(t$design$formula %||% model$formula %||% ~ 1))
    Fm <- Fd; cFm <- cF
    if (length(bv)) {
      mm <- stats::model.matrix(~ ., as.data.frame(lapply(meta[samples, bv, drop = FALSE], function(z) factor(as.character(z))), row.names = samples))[, -1, drop = FALSE]
      colnames(mm) <- paste0(".block.", colnames(mm)); Fm <- cbind(Fd, mm); cFm <- c(cF, setNames(rep(0, ncol(mm)), colnames(mm)))
      # a tested variable with one value in every level of the grouping variable (whole levels are exchanged) cannot be separated from
      # the level baselines: the magnitudes are then those of the test's design
      if (!isEstimable(Fm, cFm)) { Fm <- Fd; cFm <- cF }
    }
    mag <- matrix(NA_real_, length(types), 5, dimnames = list(types, c("total", "std", "n.eff", "se.total", "se.std"))); ng <- setNames(rep(0L, length(types)), types)
    nsamp <- setNames(as.integer(r$n.samples), types)
    infl <- infl.tot <- matrix(NA_real_, n.s, length(types), dimnames = list(samples, types)); gt <- list()
    # relabelings for the null band of the magnitudes (labels permuted within the strata, as the test)
    P0 <- NULL
    if (n.null > 0) {
      des0 <- t$design; des0$F <- Fd
      pl0 <- permutationPlan(des0, meta[samples, , drop = FALSE], samples, scheme = "block", block.vars = block.vars, n.permutations = n.null, max.enumerate = 0)
      P0 <- withSeed(as.integer(seed %% 1e9) + 1L, drawPermutations(pl0, n.null))
    }
    null.mag <- list(total = matrix(NA_real_, ncol(P0) %||% 0, length(types), dimnames = list(NULL, types))); null.mag$std <- null.mag$total
    # distinct relabelings of each type's own samples (the floor of its p-value)
    # when every relabeling was enumerated, the p-value of a cell type is a count among the global relabelings: under label relabeling a present
    # sample takes the design row of the sample the global relabeling maps it to, which may be absent from the cell type (the kernel's rule,
    # cluster_free_shift_kernel.cpp); under residual relabeling the relabeling is induced onto the present samples. The floor is the share of
    # the relabelings that reproduce the observed arrangement, and n.perm.distinct the number of distinct arrangements among them (the cell
    # type's own plan counts differently when samples are missing, e.g. a pair with one sample left)
    if (ex) { des.g <- t$design; des.g$F <- Fd
      gpl <- permutationPlan(des.g, meta[samples, , drop = FALSE], samples, scheme = permutation, block.vars = block.vars, n.permutations = B, max.enumerate = B)
      Pg <- withSeed(seed, drawPermutations(gpl, B)) }
    floors <- rep(NA_real_, length(types)); swaps <- rep(FALSE, length(types))
    n.dist <- vapply(seq_along(types), function(k) { pres <- prof[[k]]$present + 1L; if (!tested[k] || length(pres) < 2) return(NA_real_)
      dk <- t$design; dk$F <- Fd[pres, , drop = FALSE]
      sub.plan <- tryCatch(permutationPlan(dk, meta[samples[pres], , drop = FALSE], samples[pres], scheme = permutation, block.vars = block.vars, n.permutations = B, max.enumerate = 0), error = function(e) NULL)
      if (is.null(sub.plan)) return(NA_real_)
      swaps[k] <<- isTRUE(sub.plan$label.swap)
      if (!ex) return(sub.plan$n.distinct)
      rows <- apply(Fd, 1, paste, collapse = "\r")
      fl <- identical(gpl$scheme, "freedman-lane")   # residuals are moved: every distinct induced permutation counts
      keys <- if (fl) apply(Pg, 2, function(p) paste(inducePermutation(p, sub = pres, sub.plan = sub.plan, plan = gpl), collapse = ","))
              else apply(Pg, 2, function(p) paste(rows[p[pres]], collapse = "\n"))
      obs <- if (fl) paste(seq_along(pres), collapse = ",") else paste(rows[pres], collapse = "\n")
      # the gsd statistic does not change when the two labels are exchanged: when that exchange is one of the global relabelings, as
      # many relabelings reproduce the exchanged arrangement as the observed one
      floors[k] <<- (1 + isTRUE(gpl$label.swap)) * sum(keys == obs) / length(keys)
      length(unique(keys)) }, numeric(1))
    # a cell type whose samples allow a single arrangement under the relabeling cannot be tested (p = 1 by construction): not tested
    one <- tested & is.finite(n.dist) & n.dist <= 1
    if (any(one)) { tested[one] <- FALSE; padj[] <- NA_real_; padj[tested] <- stats::p.adjust(p[tested], "BH") }
    for (k in seq_along(types)) {
      pr <- prof[[k]]; pres <- pr$present + 1L
      if (!tested[k] || length(pres) < 3) next
      x <- pr$x[, pres, drop = FALSE]; v <- pr$v[, pres, drop = FALSE]; t2 <- as.numeric(pr$tau2)
      if (robust == "winsor") x <- gsdWinsorize(x, v, t2, winsor.k)
      f <- gsdGeneFit(x, v, t2, Fd[pres, , drop = FALSE], cF, gsd.d0)   # the test's model: gene changes and z
      Fk <- Fm[pres, , drop = FALSE]; fm <- gsdGeneFit(x, v, t2, Fk, cFm, gsd.d0)   # magnitudes
      if (is.null(f) || is.null(fm)) next
      vbar <- rowMeans(v); mag[k, 1:3] <- gsdMagnitudes(fm, vbar, t2); ng[k] <- nrow(x)
      gn <- rownames(cm)[pr$genes + 1L]
      gt[[types[k]]] <- data.frame(gene = gn, beta = f$beta, se = f$se, z = f$z, stringsAsFactors = FALSE)[order(-abs(f$z)), ]
      rownames(gt[[types[k]]]) <- NULL
      if (influence) {
        jk <- t(vapply(seq_along(pres), function(i) gsdMagnitudes(gsdGeneFit(x[, -i, drop = FALSE], v[, -i, drop = FALSE], t2, Fk[-i, , drop = FALSE], cFm, gsd.d0),
                                                                 rowMeans(v[, -i, drop = FALSE]), t2, breadth = FALSE)[1:2], numeric(2)))
        ok <- stats::complete.cases(jk); nj <- sum(ok)
        if (nj >= 3) {
          mag[k, 4:5] <- sqrt((nj - 1) / nj * colSums(sweep(jk[ok, , drop = FALSE], 2, colMeans(jk[ok, , drop = FALSE]))^2))
          # the jackknife counts the second-order part of the variance of a squared distance twice: count it once (correctedJackknifeSE())
          v1 <- gsdLinearVar(x, Fk, cFm, fm$beta)
          den <- max(sum(fm$rss / fm$df - vbar), 0.1 * sum(fm$rss / fm$df))
          mag[k, 4] <- correctedJackknifeSE(mag[k, 4], v1); mag[k, 5] <- correctedJackknifeSE(mag[k, 5], v1 / den^2)
          infl.tot[samples[pres], k] <- jk[, 1]; infl[samples[pres], k] <- jk[, 2]
        }
      }
      if (!is.null(P0)) for (bb in seq_len(ncol(P0))) {
        fp <- gsdGeneFit(x, v, t2, Fm[P0[, bb], , drop = FALSE][pres, , drop = FALSE], cFm, gsd.d0)
        if (!is.null(fp)) { mm0 <- gsdMagnitudes(fp, vbar, t2, breadth = FALSE); null.mag$total[bb, k] <- mm0[1]; null.mag$std[bb, k] <- mm0[2] }
      }
    }
    z <- stats::qnorm(0.975)
    rows <- data.frame(test = t$label, test.id = t$id, celltype = types, effect = "shift", estimate = mag[, "total"], estimate.norm = mag[, "std"],
                       se.jk = mag[, "se.total"], ci.low = mag[, "total"] - z * mag[, "se.total"], ci.high = mag[, "total"] + z * mag[, "se.total"],
                       statistic = stat, p = p, padj = padj, p.fwer = p.fwer, n = nsamp, n.ref = NA_integer_, n.alt = NA_integer_,
                       scheme = r$settings$permutation, n.perm = as.integer(r$settings$n.permutations), n.perm.distinct = round(pmin(n.dist, if (ex) Inf else B)),
                       p.floor = ifelse(is.finite(floors), floors, ifelse(is.finite(n.dist) & n.dist <= B, (1 + swaps) / n.dist, 1 / (B + 1))), flags = "", n.genes = ng, n.eff = mag[, "n.eff"], se.jk.norm = mag[, "se.std"], stringsAsFactors = FALSE)
    rows$flags <- ifelse(rows$p.floor > 0.05 / 10, "few-permutations", "")
    rows <- rows[tested, , drop = FALSE]; rownames(rows) <- NULL
    if (any(!tested)) skipped <- rbind(skipped, data.frame(test = t$label, celltype = types[!tested],
                                                          reason = ifelse(one[!tested], "a single arrangement of its samples under the relabeling", "too few samples with enough cells")))
    out[[t$label]] <- rows; genes.out[[t$label]] <- gt
    infl.out[[t$label]] <- list(shift = sweep(infl[, tested, drop = FALSE], 2, mag[tested, "std"]), shift.total = sweep(infl.tot[, tested, drop = FALSE], 2, mag[tested, "total"]))
    null.out[[t$label]] <- lapply(null.mag, function(m) m[, tested, drop = FALSE])
    global[[t$label]] <- data.frame(test = t$label, test.id = t$id, effect = "shift", p = if (any(tested) && any(is.finite(p.fwer[tested]))) min(p.fwer[tested], na.rm = TRUE) else NA_real_)
  }
  results <- do.call(rbind, unname(out)); if (is.null(results)) results <- data.frame()
  global <- do.call(rbind, unname(global)); rownames(skipped) <- NULL
  list(results = results, global = global, genes = genes.out, influence = infl.out, null.magnitudes = null.out, skipped = skipped, noise.curves = curves, model = model, profiles = if (return.profiles) prof.out else NULL,
       settings = list(statistic = "gsd", permutation = permutation, n.permutations = B, seed = seed, block.vars = block.vars,
                       min.cells.per.sample = min.cells.per.sample, min.samp.per.level = min.samp.per.level,
                       neighbourhood.gene.frac = neighbourhood.gene.frac, robust = robust, winsor.k = winsor.k, gsd.d0 = gsd.d0, influence = influence, scale = scale))
}

# rbind of data.frames with different columns (missing ones filled with NA)
rbindFill <- function(dfs) {
  dfs <- Filter(function(d) !is.null(d) && nrow(d) > 0, dfs)
  if (!length(dfs)) return(data.frame())
  cols <- unique(unlist(lapply(dfs, names)))
  do.call(rbind, lapply(dfs, function(d) { for (cl in setdiff(cols, names(d))) d[[cl]] <- NA; d[, cols, drop = FALSE] }))
}

# results of estimateExpressionShiftMagnitudes() with statistic = "gsd": the distance tests (expressionShiftsForModel(), all tests of the
# model) give the effect sizes (normalized shift, var, total), their jackknife errors, relabeling bands and influence; for the contrast
# tests the significance of the shift is replaced by the gsd test (gsdShiftsPerCellType()), the distance test's own kept as
# p.shift.distance / padj.shift.distance (wide) and p.distance (long)
mergeShiftResults <- function(res.gsd, res.dist, model, statistic) {
  if (is.null(res.gsd)) { res.dist$genes <- list(); return(res.dist) }
  g <- res.gsd$results; key <- function(d) paste(d$test, d$celltype, sep = "\r")
  res <- res.dist
  # long table: shift rows take the gsd p-values
  L <- res$results; L$p.distance <- NA_real_; L$padj.distance <- NA_real_
  i <- which(L$effect == "shift" & key(L) %in% key(g)); j <- match(key(L)[i], key(g))
  L$p.distance[i] <- L$p[i]; L$padj.distance[i] <- L$padj[i]
  L$n.gsd <- NA_real_; L$n.gsd[i] <- g$n[j]   # samples of the gsd test (n: samples of the distance model, min.cells.magnitude)
  for (cl in c("p", "padj", "p.fwer", "statistic", "n.perm", "n.perm.distinct", "p.floor", "flags")) if (cl %in% names(L)) L[[cl]][i] <- g[[cl]][j]
  L$scheme[i] <- g$scheme[j]
  # gsd-tested types without a distance fit (e.g. a degenerate distance matrix): significance only
  miss <- g[!(key(g) %in% key(L[L$effect == "shift", ])), , drop = FALSE]
  if (nrow(miss)) L <- rbindFill(list(L, data.frame(test = miss$test, test.id = miss$test.id, celltype = miss$celltype, effect = "shift", statistic = miss$statistic,
                                                     p = miss$p, padj = miss$padj, p.fwer = miss$p.fwer, n = NA_real_, n.gsd = miss$n, scheme = miss$scheme,
                                                     n.perm = miss$n.perm, n.perm.distinct = miss$n.perm.distinct, p.floor = miss$p.floor, flags = miss$flags, stringsAsFactors = FALSE)))
  L$n.genes <- g$n.genes[match(key(L), key(g))]; L$n.genes[L$effect != "shift"] <- NA
  res$results <- L
  # wide tables
  for (tl in intersect(names(res$wide), unique(g$test))) {
    w <- res$wide[[tl]]; gt <- g[g$test == tl, ]; k <- match(w$celltype, gt$celltype)
    w$p.shift.distance <- w$p.shift; w$padj.shift.distance <- w$padj.shift
    w$p.shift <- gt$p[k]; w$padj.shift <- gt$padj[k]; w$p.fwer.shift <- gt$p.fwer[k]; w$statistic <- gt$statistic[k]; w$n.genes <- gt$n.genes[k]
    w$n.perm.distinct.shift <- gt$n.perm.distinct[k]; w$p.floor.shift <- gt$p.floor[k]; w$n.gsd <- gt$n[k]
    extra <- gt[!(gt$celltype %in% w$celltype), , drop = FALSE]
    if (nrow(extra)) w <- rbindFill(list(w, data.frame(celltype = extra$celltype, p.shift = extra$p, padj.shift = extra$padj, p.fwer.shift = extra$p.fwer,
                                                       statistic = extra$statistic, n = NA_real_, n.gsd = extra$n, n.genes = extra$n.genes,
                                                       n.perm.distinct.shift = extra$n.perm.distinct, p.floor.shift = extra$p.floor, stringsAsFactors = FALSE)))
    res$wide[[tl]] <- w
  }
  # global p of the shift: the gsd max-statistic over cell types
  gg <- res.gsd$global
  if (!is.null(gg) && nrow(gg)) for (r in seq_len(nrow(gg))) {
    h <- which(res$global$test == gg$test[r] & res$global$effect == "shift"); if (length(h)) res$global$p[h] <- gg$p[r]
  }
  res$genes <- res.gsd$genes; res$noise.curves <- res.gsd$noise.curves
  # skipped: not tested by gsd (types tested by gsd without a distance fit appear in the results with NA magnitudes)
  sk.d <- res.dist$skipped; if (!is.null(sk.d) && nrow(sk.d)) sk.d <- sk.d[!(key(sk.d) %in% key(g)), , drop = FALSE]
  res$skipped <- rbindFill(list(sk.d, res.gsd$skipped[!(key(res.gsd$skipped) %in% key(sk.d)), , drop = FALSE]))
  if (!nrow(res$skipped)) res$skipped <- data.frame(test = character(0), celltype = character(0), reason = character(0))
  res$settings$gsd <- res.gsd$settings[c("permutation", "n.permutations", "seed", "block.vars", "min.samp.per.level", "neighbourhood.gene.frac",
                                         "robust", "winsor.k", "gsd.d0", "scale")]
  res
}
