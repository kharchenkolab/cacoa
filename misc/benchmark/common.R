## Shared code of the cluster-free adjustment benchmark (misc/benchmark/README.md).
## Inputs are light-weight lists (graph, embedding, counts [cells x genes, raw], cell.groups, sample.per.cell, sample.meta)
## cached under $S/bench/inputs_<ds>.rds (the Milo comparison's extract_inputs.R format).
suppressPackageStartupMessages({library(cacoa); library(Matrix); library(igraph)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
BENCH <- file.path(S, "bench"); dir.create(BENCH, showWarnings = FALSE, recursive = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
tm <- function(label, expr) { t0 <- Sys.time(); r <- expr; message(sprintf("[%s] %s: %.1f s", format(Sys.time(), "%H:%M:%S"), label, as.numeric(difftime(Sys.time(), t0, units = "secs")))); invisible(r) }

loadInputs <- function(ds) {
  f <- file.path(BENCH, paste0("inputs_", ds, ".rds"))
  if (!file.exists(f)) {
    alt <- file.path(S, "milo", paste0("inputs_", ds, ".rds"))
    if (file.exists(alt)) f <- alt else stop("no inputs for ", ds, ": run extract_inputs.R ", ds)
  }
  inp <- readRDS(f)
  cells <- V(inp$graph)$name
  inp$cell.groups <- droplevels(factor(inp$cell.groups[cells])); names(inp$cell.groups) <- cells
  inp$sample.per.cell <- factor(as.character(inp$sample.per.cell[cells])); names(inp$sample.per.cell) <- cells
  inp$counts <- inp$counts[cells, , drop = FALSE]
  inp$embedding <- inp$embedding[cells, , drop = FALSE]
  inp
}
# the condition column and reference level per dataset
condSpec <- function(ds) switch(ds, scc = c("condition", "Normal"), ms = c("diagnosis", "Control"), pf = c("Diagnosis", "Control"), sim = c("treatment", "control"), stop(ds))

# restrict an input list to a subset of cells (induced subgraph) and/or samples
subsetInputs <- function(inp, cells = NULL, samples = NULL) {
  if (!is.null(samples)) cells <- intersect(cells %||% V(inp$graph)$name, names(inp$sample.per.cell)[inp$sample.per.cell %in% samples])
  g <- induced_subgraph(inp$graph, cells)
  comp <- components(g)                                  # the heat kernel needs a connected graph: keep the largest component
  if (comp$no > 1) { g <- induced_subgraph(g, V(g)$name[comp$membership == which.max(comp$csize)]) }
  cells <- V(g)$name
  list(graph = g, embedding = inp$embedding[cells, , drop = FALSE], counts = inp$counts[cells, , drop = FALSE],
       cell.groups = droplevels(inp$cell.groups[cells]), sample.per.cell = droplevels(inp$sample.per.cell[cells]),
       sample.meta = inp$sample.meta[if (is.null(samples)) rownames(inp$sample.meta) else samples, , drop = FALSE])
}

# per-sample density matrices on the graph: heat kernel (cacoa default) and 1-hop neighbourhood counts
densityMats <- function(graph, spc, n.cores = 8, beta = 30) {
  w <- setNames(rep(1, nlevels(spc)), levels(spc))
  heat <- cacoa:::estimateCellDensityGraph(graph, spc, w, n.cores = n.cores, beta = beta, verbose = FALSE)
  nb <- cacoa:::estimateCellDensityGraph(graph, spc, w, smoothing = "neighbourhood", n.hops = 1)
  list(heat = heat$density.mat[V(graph)$name, ], nb = nb$density.mat[V(graph)$name, ], sizes = nb$sample.sizes)
}
logCounts <- function(dmat, sizes) { sizes <- as.numeric(sizes[colnames(dmat)]); log2(t((t(dmat) * sizes + 0.5) / (sizes + 1)) * 1e4) }

# 0-based neighbour lists (self included), as the cluster-free kernel expects
neighbourLists <- function(graph) {
  cells <- V(graph)$name
  adj <- as_adj(graph, sparse = TRUE); diag(adj) <- 1
  adj <- as(adj[cells, cells], "TsparseMatrix")
  nn <- split(adj@j, factor(adj@i + 1L, levels = seq_along(cells)))
  nn <- lapply(nn, as.integer); names(nn) <- cells
  nn
}
topGenes <- function(counts, n = 1000) { cs <- colSums(counts); names(sort(cs, decreasing = TRUE))[seq_len(min(n, length(cs)))] }

# a two-group model on a sample metadata frame with a `cond` column (levels A = reference, B)
condModel <- function(meta) buildCacoaModel(meta, formula = ~ cond, test = "cond")

# ter Braak null scale: residuals of the FULL model are permuted (Y* = F B + P E) and the statistic refitted, so the
# permuted t* are centred at the observed t and their spread is the cell's null scale free of the effect (ter Braak
# 1992). Returns the null sd per cell and the centred standardized draws (B x m) for a pooled null.
terBraakScale <- function(Fm, cvec, Y, P, chunk = 200L) {
  Fm <- as.matrix(Fm); n <- nrow(Fm); p <- ncol(Fm); cvec <- as.numeric(cvec)
  XtXi <- solve(crossprod(Fm)); A <- XtXi %*% t(Fm); cc <- as.numeric(t(cvec) %*% XtXi %*% cvec)
  Bh <- A %*% Y; fit <- Fm %*% Bh; E <- Y - fit
  t.fun <- function(Ystar) { Bs <- A %*% Ystar; rss <- colSums((Ystar - Fm %*% Bs)^2); (as.numeric(cvec %*% Bs)) / sqrt(rss / (n - p) * cc) }
  t.obs <- t.fun(Y)
  Z <- matrix(NA_real_, ncol(P), ncol(Y))
  for (b in seq_len(ncol(P))) Z[b, ] <- t.fun(fit + E[P[, b], , drop = FALSE])
  Z <- sweep(Z, 2, t.obs)                       # centred at the observed statistic: effect-free null draws
  sd.tb <- apply(Z, 2, sd)
  list(sd = sd.tb, Z = sweep(Z, 2, sd.tb, "/"), t.obs = t.obs)
}

# Weighted-least-squares permutation test for a two-level contrast (design: intercept + indicator x), closed form per
# unit, vectorized over units: Y (units x samples), w (units x samples, e.g. voom precision weights), x (0/1 per
# sample), P (samples x B index permutations of the design rows, as the "block" scheme of the engine). Prototype of
# observation weights in the permutation LM (the engine has none yet). Returns the observed t and the B x m permuted t.
permWLS <- function(Y, w, x, P) {
  n <- ncol(Y); WY <- w * Y; Sw <- rowSums(w); Swy <- rowSums(WY); Swyy <- rowSums(WY * Y)
  tfun <- function(xb) {
    Swx <- as.numeric(w %*% xb); Swxy <- as.numeric(WY %*% xb)          # x is 0/1: sum w x^2 = sum w x
    det <- Sw * Swx - Swx^2; b1 <- (Sw * Swxy - Swx * Swy) / det; b0 <- (Swy - b1 * Swx) / Sw
    rss <- Swyy - 2 * b0 * Swy - 2 * b1 * Swxy + b0^2 * Sw + 2 * b0 * b1 * Swx + b1^2 * Swx
    rss <- pmax(rss, 0); se <- sqrt(rss / (n - 2) * Sw / det); b1 / se
  }
  t.obs <- tfun(x); Tp <- matrix(NA_real_, ncol(P), nrow(Y))
  for (b in seq_len(ncol(P))) Tp[b, ] <- tfun(x[P[, b]])
  list(t.obs = setNames(t.obs, rownames(Y)), t.perm = Tp, effect = setNames({ Swx <- as.numeric(w %*% x); Swxy <- as.numeric(WY %*% x); (Sw * Swxy - Swx * Swy) / (Sw * Swx - Swx^2) }, rownames(Y)))
}
# voom precision weights and log-CPM of a count matrix (units x samples) under ~ condition
voomWeights <- function(counts, meta) {
  design <- model.matrix(~ condition, data = meta[colnames(counts), , drop = FALSE])
  y <- edgeR::DGEList(counts = counts, lib.size = colSums(counts)); y <- edgeR::calcNormFactors(y, method = "TMM")
  v <- limma::voom(y, design); list(E = v$E, w = v$weights, x = as.numeric(design[, 2]))
}
# pooled FDR of a permWLS result: raw t and per-unit standardized
permWLSq <- function(r) {
  mu <- colMeans(r$t.perm); sdv <- apply(r$t.perm, 2, sd)
  list(perm.fdr = permutationFDR(r$t.obs, r$t.perm)$q, fdr.sd = permutationFDR((r$t.obs - mu) / sdv, sweep(sweep(r$t.perm, 2, mu), 2, sdv, "/"))$q)
}

# pooled-null variants of the density test from the permuted statistics: raw t (the package's perm-fdr), per-cell
# standardized by sd / by mad, and the ter Braak scale (needs the design and the permutation matrix)
pooledVariants <- function(r, model = NULL, Y = NULL, P = NULL) {
  st <- r$stat.obs; Pm <- r$stats.perm
  mu <- colMeans(Pm, na.rm = TRUE); sdv <- apply(Pm, 2, sd, na.rm = TRUE)
  med <- apply(Pm, 2, median, na.rm = TRUE); madv <- apply(Pm, 2, mad, na.rm = TRUE)
  out <- list(scale = data.frame(sd = sdv, mad = madv, row.names = names(st)))
  q <- list(fdr.raw = permutationFDR(st, Pm)$q)
  Z <- sweep(sweep(Pm, 2, mu), 2, sdv, "/"); q$fdr.sd <- permutationFDR((st - mu) / sdv, Z)$q
  Z <- sweep(sweep(Pm, 2, med), 2, madv, "/"); q$fdr.mad <- permutationFDR((st - med) / madv, Z)$q
  if (!is.null(model) && !is.null(Y) && !is.null(P)) {
    tb <- terBraakScale(model$F, model$contrast.F[colnames(model$F)], Y[rownames(model$F), , drop = FALSE], P)
    out$scale$sd.tb <- tb$sd
    q$fdr.tb <- permutationFDR(st / tb$sd, tb$Z)$q          # observed t against the pooled effect-free null
    q$fdr.tbscale <- permutationFDR(st / tb$sd, sweep(sweep(Pm, 2, mu), 2, tb$sd, "/"))$q   # tb scale, label-permutation draws
  }
  out$q <- q
  out
}

# density test with every adjustment of interest; returns per-cell p-values / q-values and the fraction called
runDensity <- function(dmat, model, nn, B = 1999, n.cores = 8, seed = 1, label = "heat", transform = "none", sizes = NULL, max.t = FALSE, groups = NULL, variants = FALSE) {
  Y <- if (transform == "log") logCounts(dmat, sizes) else dmat
  r <- cacoa:::diffCellDensityPermutations(Y, model, perm.method = "freedman-lane", statistic = "t", n.permutations = B, n.cores = n.cores, seed = seed, moderate = TRUE)
  cells <- names(r$p.param)
  zq <- function(p, eff) { z <- sign(eff) * qnorm(pmin(pmax(p, 1e-300), 1 - 1e-16) / 2, lower.tail = FALSE); z[!is.finite(p)] <- NA; names(z) <- cells; z }
  zperm <- zq(r$p.perm, r$effect)                                       # the floored exact permutation p as a signed z
  out <- list(p.param = r$p.param, p.perm = r$p.perm, p.pooled = r$p.pooled, t = r$t, effect = r$effect, z = r$z.score)
  q <- list()
  q[["bh.param"]] <- adjustCellScores(r$t, r$p.param, method = "bh")$p.adj
  q[["bh.perm"]] <- adjustCellScores(zperm, r$p.perm, method = "bh")$p.adj
  q[["bh.pooled"]] <- adjustCellScores(r$z.score, r$p.pooled, method = "bh")$p.adj
  q[["spatial.param"]] <- adjustCellScores(r$t, r$p.param, method = "spatial-fdr", nn.list = nn[cells])$p.adj
  q[["spatial.pooled"]] <- adjustCellScores(r$z.score, r$p.pooled, method = "spatial-fdr", nn.list = nn[cells])$p.adj
  q[["nbhd.param"]] <- adjustCellScores(r$t, r$p.param, method = "neighbourhood", nn.list = nn[cells], seed = seed)$p.adj
  q[["nbhd.pooled"]] <- adjustCellScores(r$z.score, r$p.pooled, method = "neighbourhood", nn.list = nn[cells], seed = seed)$p.adj
  pf <- adjustCellScores(r$stat.obs, NULL, method = "perm-fdr", perm = r$stats.perm)
  q[["perm.fdr"]] <- pf$p.adj; out$fdr.curve <- pf$fdr.curve
  if (max.t) {
    zm <- cacoa:::adjustZScoresByPermutations(r$score, r$permut.scores, smooth = FALSE)
    q[["max.t"]] <- 2 * pnorm(-abs(zm))      # familywise-adjusted p from the adjusted z
  }
  if (variants) {
    P <- modelPermutations(model, scheme = "freedman-lane", n.permutations = B, seed = seed)$P
    v <- pooledVariants(r, model, t(Y), P)
    for (m in names(v$q)) q[[m]] <- v$q[[m]]
    out$scale <- v$scale
  }
  out$q <- as.data.frame(q); rownames(out$q) <- cells
  out$label <- label
  if (!is.null(groups)) {
    gs <- summarizeCellScoresByGroup(r$p.perm, groups, effect = r$effect, p.adj = q$bh.param, stat = r$stat.obs, stats.perm = r$stats.perm)
    gs$p.cauchy.param <- summarizeCellScoresByGroup(r$p.param, groups, method = "cauchy")$p.cauchy
    gs$p.cauchy.pooled <- summarizeCellScoresByGroup(r$p.pooled, groups, method = "cauchy")$p.cauchy
    out$groups <- gs
  }
  out
}

# cluster-free expression shifts with every adjustment
runShifts <- function(cm.t, spc, nn, model, B = 1999, n.cores = 8, seed = 1, groups = NULL) {
  r <- clusterFreeExpressionShifts(cm.t, spc, nn, model, model$meta, dist = "cor", n.permutations = B, seed = seed, adjust = TRUE, smooth = FALSE,
                                   n.cores = n.cores, verbose = FALSE, adjust.method = "perm-fdr", groups = groups)
  cells <- names(r$p.value)
  pz <- pnorm(r$z.std, lower.tail = FALSE); pz[!is.finite(r$z.std)] <- NA
  q <- list(perm.fdr = r$p.adj)
  zin <- function(p) { z <- qnorm(pmin(pmax(p, 1e-300), 1 - 1e-16), lower.tail = FALSE); z[!is.finite(p)] <- NA; names(z) <- cells; z }
  q[["bh.perm"]] <- adjustCellScores(zin(r$p.value), r$p.value, method = "bh", alternative = "greater")$p.adj
  q[["bh.z"]] <- adjustCellScores(zin(pz), pz, method = "bh", alternative = "greater")$p.adj
  q[["bh.F"]] <- adjustCellScores(zin(r$p.param), r$p.param, method = "bh", alternative = "greater")$p.adj
  q[["bh.pooled"]] <- adjustCellScores(zin(r$p.pooled), r$p.pooled, method = "bh", alternative = "greater")$p.adj
  q[["spatial.perm"]] <- adjustCellScores(zin(r$p.value), r$p.value, method = "spatial-fdr", nn.list = nn[cells], alternative = "greater")$p.adj
  q[["spatial.pooled"]] <- adjustCellScores(zin(r$p.pooled), r$p.pooled, method = "spatial-fdr", nn.list = nn[cells], alternative = "greater")$p.adj
  q[["nbhd.pooled"]] <- adjustCellScores(zin(r$p.pooled), r$p.pooled, method = "neighbourhood", nn.list = nn[cells], alternative = "greater", seed = seed)$p.adj
  out <- list(p.perm = r$p.value, p.z = pz, p.F = r$p.param, p.pooled = r$p.pooled, z = r$z.std, shifts = r$shifts, n.samples = r$n.samples, q = as.data.frame(q), fdr.curve = r$fdr.curve, n.tested = r$settings$n.cells.tested)
  if (!is.null(groups)) {
    out$groups <- summarizeCellScoresByGroup(r$p.value, groups, effect = r$shifts, p.adj = r$p.adj, group.test = r$group.test, alternative = "greater")
    out$groups$p.cauchy.pooled <- summarizeCellScoresByGroup(r$p.pooled, groups, method = "cauchy", alternative = "greater")$p.cauchy
  }
  out
}

callFractions <- function(q, alphas = c(0.05, 0.1)) {
  do.call(rbind, lapply(names(q), function(m) do.call(rbind, lapply(alphas, function(a) data.frame(method = m, alpha = a, n.tested = sum(is.finite(q[[m]])), n.called = sum(q[[m]] < a, na.rm = TRUE))))))
}
