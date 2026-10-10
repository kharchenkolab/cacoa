## Covariate screen (Track C): which metadata variables are associated with sample-level variation, in which
## cell types? Exploratory, term-level, marginal and partial, FDR over the whole grid; no contrasts.
## Statistics (api note §4.5, D26): location R2 / R2.adj (chance-corrected) / R2.partial / F with Freedman-Lane
## permutation p-values on the Gower matrix (or an analytic effective-dimension preview), dispersion F / R2;
## global per-covariate p by shared-permutation max-T across cell types.

# term columns of one covariate (no intercept); numeric covariates per SD (D32); transform strings allowed
covariateColumns <- function(meta, cov) {
  f <- stats::as.formula(paste("~", cov))
  mf <- stats::model.frame(f, meta, na.action = stats::na.pass)
  X <- stats::model.matrix(f, mf)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  for (j in seq_len(ncol(X))) if (length(unique(stats::na.omit(X[, j]))) > 2) { s <- stats::sd(X[, j], na.rm = TRUE); if (is.finite(s) && s > 0) X[, j] <- X[, j] / s }
  rownames(X) <- rownames(meta)
  X
}

# marginal / partial term tests for one Gower matrix and a set of covariates
screenOneMatrix <- function(G, meta, covariates, mode, adjust.cols = NULL, n.permutations = 0, P = NULL, dispersion = TRUE, robust = "none", robust.k = 1.345) {
  samples <- rownames(G)
  cols <- lapply(covariates, function(v) covariateColumns(meta[samples, , drop = FALSE], v)); names(cols) <- covariates
  base <- cbind(`(Intercept)` = rep(1, length(samples)), if (!is.null(adjust.cols)) adjust.cols[samples, , drop = FALSE])
  rows <- list(); perm.F <- list()
  for (v in covariates) {
    Xv <- cols[[v]]
    others <- if (mode == "partial") do.call(cbind, cols[setdiff(covariates, v)]) else NULL
    Xr <- cbind(base, others)
    ok <- stats::complete.cases(cbind(Xr, Xv))
    n.used <- sum(ok)
    if (n.used < 4) { rows[[v]] <- data.frame(covariate = v, mode = mode, n.used = n.used, stringsAsFactors = FALSE); next }
    Gk <- G[ok, ok]; Xr.k <- Xr[ok, , drop = FALSE]; Xv.k <- Xv[ok, , drop = FALSE]
    Xr.k <- Xr.k[, colSums(abs(Xr.k)) > 0 & c(TRUE, apply(Xr.k[, -1, drop = FALSE], 2, stats::sd) > 0), drop = FALSE]
    Xv.k <- Xv.k[, apply(Xv.k, 2, stats::sd) > 0, drop = FALSE]
    if (!ncol(Xv.k)) { rows[[v]] <- data.frame(covariate = v, mode = mode, n.used = n.used, stringsAsFactors = FALSE); next }
    Xf.k <- cbind(Xr.k, Xv.k)
    if (!all(ok)) Gk <- gowerCenter(uncenterGower(Gk))
    tt <- termTestGower(Gk, Xf.k, Xr.k)
    dd <- if (dispersion) dispersionTermTest(Gk, Xf.k, Xf.k, Xr.k) else c(F.disp = NA, p.disp = NA, df.disp = NA, nu.disp = NA)
    w.obs <- rep(1, n.used)
    if (!identical(robust, "none")) {                         # robust statistics (R2 stays unweighted)
      wt <- weightedTermStats(Gk, Xf.k, Xr.k, Xf.k, Xr.k, robust = robust, k = robust.k)
      tt["F"] <- wt$F; tt["p.analytic"] <- NA_real_; if (dispersion) dd["F.disp"] <- wt$F.disp; w.obs <- wt$w
    }
    r2d <- if (dispersion) dispersionR2(Gk, Xf.k, Xf.k, Xr.k) else NA_real_
    p.perm <- NA_real_; p.disp.perm <- NA_real_; Fp <- NULL; Fdp <- NULL
    if (n.permutations > 0 && is.finite(tt["F"])) {
      parts <- flGowerParts(Gk, Xr.k)
      Pk <- if (!is.null(P)) inducePermutationSimple(P, which(ok)) else replicate(n.permutations, sample.int(n.used))
      if (!is.matrix(Pk)) Pk <- matrix(Pk, ncol = 1)
      storage.mode(Pk) <- "integer"; B <- ncol(Pk)
      st <- if (identical(robust, "none")) {                        # C++ kernels; R references: termTestGower / dispersionTermTest / weightedTermStats loops
        k <- termKernelInputs(Xf.k, Xr.k, Xf.k, Xr.k, n.used)
        permuted_term_stats_fl(parts$K1, parts$K2, parts$K3, parts$K4, k$Hf, k$Hr, k$Zf, k$Zr, k$df, k$nu, k$qZf, k$qZr, Pk, dispersion && is.finite(dd["F.disp"]))
      } else permuted_term_stats_w(Gk, Xf.k, Xr.k, Xf.k, Xr.k, rep(1, n.used), Pk, TRUE, w.obs, dispersion && is.finite(dd["F.disp"]),
                                   c(huber = 1L, winsor = 2L)[[robust]], robust.k, 5L)
      Fp <- st[, 1]; Fdp <- st[, 2]
      p.perm <- (sum(Fp >= tt["F"] - 1e-12) + 1) / (B + 1)
      if (dispersion && is.finite(dd["F.disp"])) p.disp.perm <- (sum(Fdp >= dd["F.disp"] - 1e-12, na.rm = TRUE) + 1) / (sum(is.finite(Fdp)) + 1)
    }
    rows[[v]] <- data.frame(covariate = v, mode = mode, n.used = n.used, R2 = unname(tt["R2"]), R2.adj = unname(tt["R2.adj"]),
                            R2.partial = unname(tt["R2.partial"]), F = unname(tt["F"]), df = unname(tt["df"]), r.eff = unname(tt["r.eff"]),
                            p.analytic = unname(tt["p.analytic"]), p.perm = p.perm, F.disp = unname(dd["F.disp"]), R2.disp.adj = r2d,
                            p.disp.analytic = unname(dd["p.disp"]), p.disp.perm = p.disp.perm, stringsAsFactors = FALSE)
    perm.F[[v]] <- list(F = Fp, F.disp = Fdp)
  }
  list(table = do.call(plyr::rbind.fill, rows), perm = perm.F)
}

# restrict a global permutation matrix (over all samples) to the positions `idx`, as a uniform permutation of them
inducePermutationSimple <- function(P, idx) {
  apply(P, 2, function(p) rank(p[idx], ties.method = "first"))
}

#' Screen covariates against per-cell-type sample distances
#'
#' For every covariate and cell type: how much of the between-sample variation the covariate explains alone
#' (`marginal`) and after adjusting for the other screened covariates (`partial`), with chance-corrected R2,
#' Freedman-Lane permutation p-values (shared permutations across cell types), BH adjustment over the whole
#' grid, and a per-covariate global p-value (max-statistic across cell types). Also an association of the
#' covariate with per-sample dispersion.
#'
#' @param D.list named list of sample x sample distance matrices (one per cell type)
#' @param meta sample metadata (rows named by sample)
#' @param covariates covariates to screen (default: usable columns of `meta`, see [describeMetadata()]);
#'   transform strings such as `"log(pmi)"` are allowed
#' @param mode `"both"` (default), `"partial"` or `"marginal"`
#' @param adjust.for optional formula of covariates always adjusted for (both modes)
#' @param dist distance type of the matrices
#' @param test.variable optional variable the covariates are related to (over-adjustment annotation)
#' @param n.permutations number of permutations (default 199)
#' @param min.samples.per.type cell types with fewer samples are skipped (default 6)
#' @param max.partial.df when the partial adjustment set would leave fewer than this many residual degrees of
#'   freedom in a cell type, partial tests fall back to `adjust.for` only (default 5)
#' @param alpha FDR level for the summary (default 0.05)
#' @param seed integer seed
#' @param n.cores cores for the per-cell-type loop
#' @param verbose print the summary
#' @return list of class `cacoaCovariateScreen`: `table` (long: celltype, covariate, mode, n.used, R2, R2.adj,
#'   R2.partial, F, df, r.eff, p, padj, F.disp, R2.disp.adj, p.disp, padj.disp, p.source, fallback),
#'   `global` (per covariate x mode: max-T p, number of significant cell types, median R2.adj), `suggestion`
#'   (formula), `settings`, `notes`
#' @export
#' @param robust,robust.k robust location / dispersion statistics (`"none"`, `"huber"`, `"winsor"`; tuning constant), weights
#'   re-estimated under every relabeling
screenCovariates <- function(D.list, meta, covariates = NULL, mode = c("both", "partial", "marginal"), adjust.for = NULL,
                             dist = c("cor", "l2", "l1"), test.variable = NULL,
                             n.permutations = 199, min.samples.per.type = 6, max.partial.df = 5, alpha = 0.05, seed = 1,
                             n.cores = 1, verbose = FALSE, robust = c("none", "huber", "winsor"), robust.k = 1.345) {
  mode <- match.arg(mode); dist <- match.arg(dist); robust <- match.arg(robust); p.values <- "permutation"
  desc <- describeMetadata(meta)
  if (is.null(covariates)) covariates <- desc$column[desc$role == "usable"]
  covariates <- unique(covariates)
  bad <- intersect(covariates, desc$column[desc$role != "usable"])
  notes <- character(0)
  if (length(bad)) { notes <- c(notes, sprintf("excluded: %s", paste(sprintf("%s (%s)", bad, desc$role[match(bad, desc$column)]), collapse = ", "))); covariates <- setdiff(covariates, bad) }
  if (!length(covariates)) stop("no covariates to screen")
  adj.vars <- if (!is.null(adjust.for)) intersect(all.vars(stats::as.formula(adjust.for)), names(meta)) else character(0)
  covariates <- setdiff(covariates, adj.vars)
  modes <- if (mode == "both") c("marginal", "partial") else mode
  all.samples <- rownames(meta)
  D.list <- Filter(function(D) !is.null(D) && nrow(D) >= min.samples.per.type, D.list)
  if (!length(D.list)) stop("no cell type has at least ", min.samples.per.type, " samples")
  if (any(vapply(D.list, function(D) is.null(rownames(D)) || !all(rownames(D) %in% all.samples), logical(1))))
    stop("every distance matrix must be named by samples present in `meta`")
  G.list <- lapply(D.list, function(D) gowerCenter(asSquaredDistance(D, dist)))
  adjust.cols <- if (length(adj.vars)) do.call(cbind, lapply(adj.vars, function(v) covariateColumns(meta, v))) else NULL
  nperm <- if (p.values == "permutation") n.permutations else 0
  P <- if (nperm > 0) withSeed(seed, replicate(nperm, sample.int(length(all.samples)))) else NULL
  if (!is.null(P)) rownames(P) <- all.samples

  runOne <- function(ct) {
    G <- G.list[[ct]]; samples <- rownames(G)
    Pk <- if (!is.null(P)) P[samples, , drop = FALSE] else NULL
    if (!is.null(Pk)) Pk <- apply(Pk, 2, rank, ties.method = "first")
    out <- list()
    for (m in modes) {
      fallback <- FALSE
      covs.m <- covariates
      if (m == "partial") {
        q.all <- 1 + (if (!is.null(adjust.cols)) ncol(adjust.cols) else 0) + sum(vapply(covariates, function(v) ncol(covariateColumns(meta[samples, , drop = FALSE], v)), integer(1)))
        if (length(samples) - q.all < max.partial.df) { fallback <- TRUE }
      }
      r <- if (m == "partial" && fallback) screenOneMatrix(G, meta, covariates, "marginal", adjust.cols, nperm, Pk, robust = robust, robust.k = robust.k) else
        screenOneMatrix(G, meta, covariates, m, adjust.cols, nperm, Pk, robust = robust, robust.k = robust.k)
      tb <- r$table; tb$mode <- m; tb$fallback <- fallback; tb$celltype <- ct
      out[[m]] <- list(table = tb, perm = r$perm)
    }
    out
  }
  cts <- names(G.list)
  runs <- if (n.cores > 1 && length(cts) > 1) sccore::plapply(cts, runOne, n.cores = n.cores, progress = verbose, fail.on.error = TRUE) else lapply(cts, runOne)
  names(runs) <- cts
  tab <- do.call(plyr::rbind.fill, unlist(lapply(runs, function(r) lapply(r, `[[`, "table")), recursive = FALSE))
  tab$p <- if (p.values == "permutation") tab$p.perm else tab$p.analytic
  tab$p.disp <- if (p.values == "permutation") tab$p.disp.perm else tab$p.disp.analytic
  tab$p.source <- p.values
  tab$padj <- NA_real_; tab$padj.disp <- NA_real_
  for (m in modes) { i <- tab$mode == m; tab$padj[i] <- stats::p.adjust(tab$p[i], "BH"); tab$padj.disp[i] <- stats::p.adjust(tab$p.disp[i], "BH") }
  # global per covariate x mode: max-T over cell types from the shared permutations
  global <- list()
  for (m in modes) for (v in covariates) {
    rows <- tab[tab$mode == m & tab$covariate == v, ]
    gp <- NA_real_
    if (nperm > 0) {
      obs <- rows$F; perm <- lapply(rows$celltype, function(ct) runs[[ct]][[m]]$perm[[v]]$F)
      okc <- is.finite(obs) & vapply(perm, function(x) !is.null(x) && length(x) == nperm, logical(1))
      if (any(okc)) {
        Pm <- do.call(cbind, perm[okc]); mu <- colMeans(Pm); sdv <- apply(Pm, 2, stats::sd); sdv[!is.finite(sdv) | sdv == 0] <- NA
        z.obs <- (obs[okc] - mu) / sdv; z.perm <- sweep(sweep(Pm, 2, mu), 2, sdv, "/")
        u <- is.finite(z.obs)
        if (any(u)) { mx <- apply(z.perm[, u, drop = FALSE], 1, max, na.rm = TRUE); gp <- (sum(mx >= max(z.obs[u])) + 1) / (length(mx) + 1) }
      }
    } else if (any(is.finite(rows$p))) gp <- min(1, min(rows$p, na.rm = TRUE) * sum(is.finite(rows$p)))   # Bonferroni-min for the analytic preview
    global[[length(global) + 1]] <- data.frame(covariate = v, mode = m, p.global = gp, n.sig = sum(rows$padj < alpha, na.rm = TRUE),
                                               n.sig.disp = sum(rows$padj.disp < alpha, na.rm = TRUE), n.types = nrow(rows),
                                               median.R2.adj = stats::median(rows$R2.adj, na.rm = TRUE), stringsAsFactors = FALSE)
  }
  global <- do.call(rbind, global)
  if (!is.null(test.variable) && test.variable %in% names(meta)) {
    global$assoc.test <- vapply(global$covariate, function(v) if (v == test.variable) NA_real_ else
      tryCatch(covariateAssociation(meta[[test.variable]], stats::model.frame(stats::as.formula(paste("~", v)), meta, na.action = stats::na.pass)[[1]]), error = function(e) NA_real_), numeric(1))
  }
  # suggestion: covariates associated in >= min(3, n.types/4) cell types in the partial (or only) mode
  sm <- if ("partial" %in% modes) "partial" else modes[1]
  g <- global[global$mode == sm, ]
  thr <- max(1, min(3, floor(length(cts) / 4)))
  sel <- g$covariate[g$n.sig >= thr]
  if (!is.null(test.variable)) sel <- setdiff(sel, test.variable)
  over <- if (!is.null(global$assoc.test)) g$covariate[is.finite(g$assoc.test) & g$assoc.test > 0.5] else character(0)
  terms <- unique(c(test.variable, adj.vars, setdiff(sel, over)))
  suggestion <- if (length(terms)) stats::as.formula(paste("~", paste(terms, collapse = " + "))) else NULL
  res <- structure(list(table = tab, global = global, suggestion = suggestion, covariates = covariates, modes = modes,
                        settings = list(dist = dist, p.values = p.values, robust = robust, robust.k = robust.k, n.permutations = nperm, adjust.for = adjust.for, alpha = alpha,
                                        test.variable = test.variable, threshold.types = thr, seed = seed),
                        notes = notes, celltypes = cts, over.adjustment = over),
                   class = c("cacoaCovariateScreen", "list"))
  if (verbose) print(res)
  res
}

#' @exportS3Method base::print
print.cacoaCovariateScreen <- function(x, ...) {
  a <- x$settings$alpha; sm <- if ("partial" %in% x$modes) "partial" else x$modes[1]
  g <- x$global[x$global$mode == sm, ]; g <- g[order(-g$n.sig, g$p.global), ]
  cat(sprintf("Covariate screen: %d covariates x %d cell types, %s mode, %s p-values%s\n", length(x$covariates), length(x$celltypes), paste(x$modes, collapse = " + "),
              x$settings$p.values, if (x$settings$p.values == "analytic") " (approximate)" else sprintf(" (%d permutations)", x$settings$n.permutations)))
  hit <- g[g$n.sig >= x$settings$threshold.types, ]
  if (nrow(hit)) cat(sprintf("Associated with expression in >= %d cell types (%s, FDR %.0f%%): %s\n", x$settings$threshold.types, sm, 100 * a,
                             paste(sprintf("%s (%d types, global p %.3g)", hit$covariate, hit$n.sig, hit$p.global), collapse = ", ")))
  else cat(sprintf("No covariate is associated with expression in >= %d cell types (%s, FDR %.0f%%).\n", x$settings$threshold.types, sm, 100 * a))
  if ("marginal" %in% x$modes && "partial" %in% x$modes) {
    gm <- x$global[x$global$mode == "marginal", ]; gp <- x$global[x$global$mode == "partial", ]
    only <- gm$covariate[gm$n.sig >= x$settings$threshold.types & gp$n.sig[match(gm$covariate, gp$covariate)] < x$settings$threshold.types]
    if (length(only)) cat(sprintf("Marginal only (explained by other covariates): %s\n", paste(only, collapse = ", ")))
  }
  dsp <- g[g$n.sig.disp >= x$settings$threshold.types, ]
  if (nrow(dsp)) cat(sprintf("Associated with sample dispersion: %s\n", paste(sprintf("%s (%d types)", dsp$covariate, dsp$n.sig.disp), collapse = ", ")))
  if (length(x$over.adjustment)) cat(sprintf("Strongly tied to '%s' (possible over-adjustment): %s\n", x$settings$test.variable,
                                             paste(sprintf("%s (%.2f)", x$over.adjustment, g$assoc.test[match(x$over.adjustment, g$covariate)]), collapse = ", ")))
  if (!is.null(x$suggestion)) cat(sprintf("Suggested model: %s   (not applied; see cao$setModel())\n", paste(deparse(x$suggestion), collapse = "")))
  for (nt in x$notes) cat("Note: ", nt, "\n", sep = "")
  invisible(x)
}

#' Variance partition of a small set of covariates
#'
#' Chance-corrected fractions of the between-sample variation explained uniquely by each covariate, shared among
#' them, and residual (api note E2).
#'
#' @param G Gower matrix (or a distance matrix with `dist`)
#' @param meta sample metadata
#' @param covariates 2-4 covariates
#' @param dist when `G` is a distance matrix, its type; `NULL` means `G` is already a Gower matrix
#' @return named numeric: one entry per covariate (unique), `shared`, `residual`; attribute `raw` holds the
#'   uncorrected fractions
#' @export
variancePartition <- function(G, meta, covariates, dist = NULL) {
  if (!is.null(dist)) G <- gowerCenter(asSquaredDistance(G, dist))
  samples <- rownames(G); meta <- meta[samples, , drop = FALSE]
  cols <- lapply(covariates, function(v) covariateColumns(meta, v)); names(cols) <- covariates
  ok <- stats::complete.cases(do.call(cbind, cols))
  G <- gowerCenter(uncenterGower(G[ok, ok])); n <- sum(ok); one <- rep(1, n); tss <- sum(diag(G))
  r2 <- function(vs) {
    if (!length(vs)) return(c(r2 = 0, q = 0))
    X <- cbind(one, do.call(cbind, lapply(cols[vs], function(x) x[ok, , drop = FALSE])))
    hi <- hatInfo(X); ss <- sum(hi$H * G) - sum(G) / n
    c(r2 = ss / tss, q = hi$rank - 1)
  }
  full <- r2(covariates)
  adj <- function(ss.frac, q) { # omega^2-type correction with the full model's residual mean square
    rss <- tss * (1 - full[["r2"]]); nu <- n - 1 - full[["q"]]
    (ss.frac * tss - q * rss / nu) / (tss + rss / nu)
  }
  uniq.raw <- vapply(covariates, function(v) { red <- r2(setdiff(covariates, v)); full[["r2"]] - red[["r2"]] }, numeric(1))
  uniq.q <- vapply(covariates, function(v) full[["q"]] - r2(setdiff(covariates, v))[["q"]], numeric(1))
  uniq <- mapply(adj, uniq.raw, uniq.q)
  full.adj <- adj(full[["r2"]], full[["q"]])
  shared <- full.adj - sum(uniq)
  out <- c(uniq, shared = shared, residual = 1 - full.adj)
  attr(out, "raw") <- c(uniq.raw, shared = full[["r2"]] - sum(uniq.raw), residual = 1 - full[["r2"]])
  out
}

# ---- plots ---------------------------------------------------------------------------------------

#' Heatmap of a covariate screen
#'
#' @param x `cacoaCovariateScreen`
#' @param effect `"location"` (default) or `"dispersion"`
#' @param mode `"partial"` (default when available), `"marginal"` or `"both"` (two panels)
#' @param value fill: `"R2.adj"` (default), `"R2"` or `"neglog10p"`
#' @param cell.types,covariates optional subsets
#' @param cluster order rows and columns by clustering (default TRUE)
#' @param alpha significance level for the marks (default: the screen's)
#' @param plot.theme ggplot2 theme
#' @return ggplot2 object. Dot = `padj < alpha`; ring = significant marginally but not partially (when both modes
#'   are available); a "global" column shows the per-covariate global p-value.
#' @export
plotCovariateScreen <- function(x, effect = c("location", "dispersion"), mode = NULL, value = c("R2.adj", "R2", "neglog10p"),
                                cell.types = NULL, covariates = NULL, cluster = TRUE, alpha = NULL, plot.theme = ggplot2::theme_bw()) {
  effect <- match.arg(effect); value <- match.arg(value); alpha <- alpha %||% x$settings$alpha
  if (is.null(mode)) mode <- if ("partial" %in% x$modes) "partial" else x$modes[1]
  modes <- if (mode == "both") x$modes else mode
  tb <- x$table[x$table$mode %in% modes, ]
  if (!is.null(cell.types)) tb <- tb[tb$celltype %in% cell.types, ]
  if (!is.null(covariates)) tb <- tb[tb$covariate %in% covariates, ]
  if (effect == "dispersion") { tb$est <- tb$R2.disp.adj; tb$pv <- tb$p.disp; tb$pa <- tb$padj.disp } else { tb$est <- if (value == "R2") tb$R2 else tb$R2.adj; tb$pv <- tb$p; tb$pa <- tb$padj }
  tb$fill <- if (value == "neglog10p") -log10(pmax(tb$pv, 1e-6)) else pmax(tb$est, 0)
  tb$neg <- !is.na(tb$est) & tb$est < 0 & value != "neglog10p"
  tb$sig <- !is.na(tb$pa) & tb$pa < alpha
  # ring: significant marginally but not partially
  tb$ring <- FALSE
  if (all(c("marginal", "partial") %in% x$modes)) {
    full <- x$table
    key <- paste(full$celltype, full$covariate)
    pm <- full[full$mode == "marginal", ]; pp <- full[full$mode == "partial", ]
    sig.m <- setNames(!is.na(if (effect == "dispersion") pm$padj.disp else pm$padj) & (if (effect == "dispersion") pm$padj.disp else pm$padj) < alpha, paste(pm$celltype, pm$covariate))
    sig.p <- setNames(!is.na(if (effect == "dispersion") pp$padj.disp else pp$padj) & (if (effect == "dispersion") pp$padj.disp else pp$padj) < alpha, paste(pp$celltype, pp$covariate))
    k <- paste(tb$celltype, tb$covariate)
    tb$ring <- tb$mode == "partial" & !tb$sig & (sig.m[k] %in% TRUE)
  }
  ord.ct <- unique(tb$celltype); ord.cv <- unique(tb$covariate)
  if (cluster && length(ord.ct) > 2 && length(ord.cv) > 1) {
    m <- reshape2::acast(tb[tb$mode == modes[1], ], covariate ~ celltype, value.var = "fill", fun.aggregate = mean)
    m[!is.finite(m)] <- 0
    if (nrow(m) > 2) ord.cv <- rownames(m)[stats::hclust(stats::dist(m))$order]
    if (ncol(m) > 2) ord.ct <- colnames(m)[stats::hclust(stats::dist(t(m)))$order]
  }
  # global column
  g <- x$global[x$global$mode %in% modes & x$global$covariate %in% ord.cv, ]
  gl <- data.frame(celltype = "global", covariate = g$covariate, mode = g$mode, fill = NA_real_, neg = FALSE, sig = !is.na(g$p.global) & g$p.global < alpha, ring = FALSE,
                   label = ifelse(is.na(g$p.global), "", formatC(g$p.global, format = "g", digits = 2)), stringsAsFactors = FALSE)
  tb$label <- ""
  cols <- c("celltype", "covariate", "mode", "fill", "neg", "sig", "ring", "label")
  df <- rbind(tb[, cols], gl[, cols])
  df$celltype <- factor(df$celltype, levels = c(ord.ct, "global")); df$covariate <- factor(df$covariate, levels = ord.cv)
  df$mode <- factor(df$mode, levels = modes)
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$celltype, y = .data$covariate)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white") +
    ggplot2::geom_tile(data = df[df$neg, ], fill = "grey85", colour = "white") +
    ggplot2::geom_point(data = df[df$sig & df$celltype != "global", ], size = 1.8, colour = "black") +
    ggplot2::geom_point(data = df[df$ring, ], size = 2.6, shape = 1, colour = "black", stroke = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.6) +
    ggplot2::geom_point(data = df[df$sig & df$celltype == "global", ], ggplot2::aes(y = .data$covariate), x = length(ord.ct) + 1.4, size = 1.8, colour = "black") +
    ggplot2::scale_fill_gradient(low = "white", high = "#2c7fb8", na.value = "grey95", name = if (value == "neglog10p") "-log10 p" else if (effect == "dispersion") "adj. R2 (disp.)" else value) +
    plot.theme + ggplot2::labs(x = NULL, y = NULL, title = sprintf("Covariate screen: %s (%s)", effect, paste(modes, collapse = " / ")),
                               subtitle = sprintf("dot: FDR < %.2g; ring: marginal only; grey: negative adj. R2; 'global' column: max-statistic p across cell types", alpha)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"))
  if (length(modes) > 1) gg <- gg + ggplot2::facet_wrap(~ mode, nrow = 1)
  gg
}

#' Summary bars of a covariate screen
#'
#' One bar per covariate: median chance-corrected R2 across cell types for location and dispersion, with the
#' number of significant cell types as a label.
#' @param x `cacoaCovariateScreen`
#' @param mode `"partial"` (default when available) or `"marginal"`
#' @param plot.theme ggplot2 theme
#' @return ggplot2 object
#' @export
plotCovariateSummary <- function(x, mode = NULL, plot.theme = ggplot2::theme_bw()) {
  if (is.null(mode)) mode <- if ("partial" %in% x$modes) "partial" else x$modes[1]
  tb <- x$table[x$table$mode == mode, ]
  agg <- do.call(rbind, lapply(split(tb, tb$covariate), function(d) data.frame(
    covariate = d$covariate[1], effect = c("location", "dispersion"),
    value = c(stats::median(d$R2.adj, na.rm = TRUE), stats::median(d$R2.disp.adj, na.rm = TRUE)),
    n.sig = c(sum(d$padj < x$settings$alpha, na.rm = TRUE), sum(d$padj.disp < x$settings$alpha, na.rm = TRUE)), stringsAsFactors = FALSE)))
  agg$value[!is.finite(agg$value)] <- 0
  ord <- agg$covariate[agg$effect == "location"][order(agg$value[agg$effect == "location"])]
  agg$covariate <- factor(agg$covariate, levels = ord); agg$effect <- factor(agg$effect, levels = c("location", "dispersion"))
  ggplot2::ggplot(agg, ggplot2::aes(x = .data$covariate, y = pmax(.data$value, 0), fill = .data$effect)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$n.sig > 0, .data$n.sig, "")), position = ggplot2::position_dodge(width = 0.8), hjust = -0.2, size = 3) +
    ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = c(location = "#2c7fb8", dispersion = "#f03b20")) + plot.theme +
    ggplot2::labs(x = NULL, y = "median chance-corrected R2 across cell types", fill = NULL,
                  title = sprintf("Covariate summary (%s); labels: cell types with FDR < %.2g", mode, x$settings$alpha))
}

#' Stacked bars of a variance partition per cell type
#' @param parts named list (per cell type) of [variancePartition()] results
#' @param plot.theme ggplot2 theme
#' @return ggplot2 object
#' @export
plotVariancePartition <- function(parts, plot.theme = ggplot2::theme_bw()) {
  df <- do.call(rbind, lapply(names(parts), function(ct) data.frame(celltype = ct, component = names(parts[[ct]]), value = pmax(as.numeric(parts[[ct]]), 0), stringsAsFactors = FALSE)))
  comps <- unique(df$component); comps <- c(setdiff(comps, c("shared", "residual")), "shared", "residual")
  df$component <- factor(df$component, levels = rev(comps))
  ord <- names(parts)[order(vapply(parts, function(p) 1 - p[["residual"]], numeric(1)))]
  df$celltype <- factor(df$celltype, levels = ord)
  pal <- c(setNames(grDevices::hcl.colors(length(comps) - 2, "Dark 3"), comps[seq_len(length(comps) - 2)]), shared = "grey60", residual = "grey90")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$celltype, y = .data$value, fill = .data$component)) + ggplot2::geom_col(width = 0.75, colour = "grey30", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = pal, breaks = comps) + ggplot2::coord_flip() + plot.theme +
    ggplot2::labs(x = NULL, y = "fraction of between-sample variation (chance-corrected)", fill = NULL, title = "Variance partition")
}


#' Screen covariates against the residuals of a per-column test
#'
#' Diagnostic for omitted covariates on the cell-density and cluster-free DE results: the residuals of the fitted
#' model (samples x bins, or one samples x cells matrix per gene) are turned into sample distance matrices and
#' screened marginally with [screenCovariates()], so a covariate that still structures the residuals shows up with
#' the same table, global max-statistic p-value and plots as the covariate screen of the expression shifts. The
#' covariates already in the model carry no residual structure by construction and are excluded.
#'
#' @param residuals samples x features residual matrix (rows named by sample), or a named list of such matrices
#'   (one per gene for cluster-free DE)
#' @param meta sample metadata (rows named by sample)
#' @param covariates covariates to screen (default: all usable columns minus `exclude`)
#' @param exclude covariates not to screen (the model's variables)
#' @param n.permutations,seed,n.cores,test.variable,robust,robust.k,verbose passed to [screenCovariates()]
#' @param min.samples matrices with fewer samples carrying finite residuals are skipped (default 4)
#' @return object of class `cacoaCovariateScreen` (settings `space = "residuals"`)
#' @export
screenResidualCovariates <- function(residuals, meta, covariates = NULL, exclude = NULL, n.permutations = 199, seed = 1, n.cores = 1,
                                     test.variable = NULL, robust = c("none", "huber", "winsor"), robust.k = 1.345, min.samples = 4, verbose = FALSE) {
  robust <- match.arg(robust)
  mats <- if (is.matrix(residuals) || is.data.frame(residuals)) list(residuals = as.matrix(residuals)) else residuals
  if (!length(mats)) stop("no residual matrices given")
  D.list <- lapply(mats, function(R) {
    R <- as.matrix(R); if (is.null(rownames(R))) stop("residual matrices must have the samples as row names")
    R <- R[rowSums(is.finite(R)) > 0, , drop = FALSE]
    if (nrow(R) < min.samples) return(NULL)
    as.matrix(stats::dist(R))
  })
  D.list <- Filter(Negate(is.null), D.list)
  if (!length(D.list)) stop("no residual matrix has at least ", min.samples, " samples")
  desc <- describeMetadata(meta)
  if (is.null(covariates)) covariates <- desc$column[desc$role == "usable"]
  covariates <- setdiff(covariates, exclude)
  if (!length(covariates)) stop("no covariate left to screen (all are in the model)")
  res <- screenCovariates(D.list, meta, covariates = covariates, mode = "marginal", dist = "l2", test.variable = test.variable,
                          n.permutations = n.permutations, min.samples.per.type = min.samples, seed = seed, n.cores = n.cores, verbose = verbose,
                          robust = robust, robust.k = robust.k)
  res$settings$space <- "residuals"; res$settings$excluded <- exclude
  res
}
