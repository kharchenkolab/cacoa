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

# density test with every adjustment of interest; returns per-cell p-values / q-values and the fraction called
runDensity <- function(dmat, model, nn, B = 1999, n.cores = 8, seed = 1, label = "heat", transform = "none", sizes = NULL, max.t = FALSE, groups = NULL) {
  Y <- if (transform == "log") logCounts(dmat, sizes) else dmat
  r <- cacoa:::diffCellDensityPermutations(Y, model, perm.method = "freedman-lane", statistic = "t", n.permutations = B, n.cores = n.cores, seed = seed, moderate = TRUE)
  cells <- names(r$p.param)
  out <- list(p.param = r$p.param, p.perm = r$p.perm, t = r$t, effect = r$effect, z = r$z.score)
  q <- list()
  q[["bh.param"]] <- adjustCellScores(r$t, r$p.param, method = "bh")$p.adj
  q[["bh.perm"]] <- adjustCellScores(r$z.score, r$p.perm, method = "bh")$p.adj
  q[["spatial.param"]] <- adjustCellScores(r$t, r$p.param, method = "spatial-fdr", nn.list = nn[cells])$p.adj
  q[["nbhd.param"]] <- adjustCellScores(r$t, r$p.param, method = "neighbourhood", nn.list = nn[cells], seed = seed)$p.adj
  pf <- adjustCellScores(r$stat.obs, NULL, method = "perm-fdr", perm = r$stats.perm)
  q[["perm.fdr"]] <- pf$p.adj; out$fdr.curve <- pf$fdr.curve
  if (max.t) {
    zm <- cacoa:::adjustZScoresByPermutations(r$score, r$permut.scores, smooth = FALSE)
    q[["max.t"]] <- 2 * pnorm(-abs(zm))      # familywise-adjusted p from the adjusted z
  }
  out$q <- as.data.frame(q); rownames(out$q) <- cells
  out$label <- label
  if (!is.null(groups)) {
    gs <- summarizeCellScoresByGroup(r$p.perm, groups, effect = r$effect, p.adj = q$bh.param, stat = r$stat.obs, stats.perm = r$stats.perm)
    gs$p.cauchy.param <- summarizeCellScoresByGroup(r$p.param, groups, method = "cauchy")$p.cauchy
    out$groups <- gs
  }
  out
}

# cluster-free expression shifts with every adjustment
runShifts <- function(cm.t, spc, nn, model, B = 1999, n.cores = 8, seed = 1, groups = NULL) {
  r <- clusterFreeExpressionShifts(cm.t, spc, nn, model, model$meta, dist = "cor", n.permutations = B, seed = seed, adjust = TRUE, smooth = FALSE,
                                   n.cores = n.cores, verbose = FALSE, adjust.method = "perm-fdr", groups = groups)
  cells <- names(r$p.value)
  pz <- pnorm(r$z.score, lower.tail = FALSE); pz[!is.finite(r$z.score)] <- NA
  q <- list(perm.fdr = r$p.adj)
  zin <- function(p) { z <- qnorm(pmin(pmax(p, 1e-300), 1 - 1e-16), lower.tail = FALSE); z[!is.finite(p)] <- NA; names(z) <- cells; z }
  q[["bh.perm"]] <- adjustCellScores(zin(r$p.value), r$p.value, method = "bh", alternative = "greater")$p.adj
  q[["bh.z"]] <- adjustCellScores(zin(pz), pz, method = "bh", alternative = "greater")$p.adj
  q[["bh.F"]] <- adjustCellScores(zin(r$p.param), r$p.param, method = "bh", alternative = "greater")$p.adj
  q[["spatial.perm"]] <- adjustCellScores(zin(r$p.value), r$p.value, method = "spatial-fdr", nn.list = nn[cells], alternative = "greater")$p.adj
  out <- list(p.perm = r$p.value, p.z = pz, p.F = r$p.param, z = r$z.score, shifts = r$shifts, n.samples = r$n.samples, q = as.data.frame(q), fdr.curve = r$fdr.curve, n.tested = r$settings$n.cells.tested)
  if (!is.null(groups)) out$groups <- summarizeCellScoresByGroup(r$p.value, groups, effect = r$shifts, p.adj = r$p.adj, group.test = r$group.test, alternative = "greater")
  out
}

callFractions <- function(q, alphas = c(0.05, 0.1)) {
  do.call(rbind, lapply(names(q), function(m) do.call(rbind, lapply(alphas, function(a) data.frame(method = m, alpha = a, n.tested = sum(is.finite(q[[m]])), n.called = sum(q[[m]] < a, na.rm = TRUE))))))
}
