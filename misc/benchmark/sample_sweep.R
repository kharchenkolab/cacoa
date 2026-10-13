## Sample-number sweep on the final density configuration (package code path): the whole-type planted sets (100 %
## removal, medium and large types, 3 reps) with each arm subsampled to n = 4 / 6 / 9 samples (balanced within the
## true-condition strata), tested with the kNN count unit at the default size rule and at a fixed k = 200, with the
## averaged post-processing and without. Usage: Rscript sample_sweep.R <ds> <B> <n.cores>  -> $S/bench/sweep/<ds>/<id>_n<n>.rds
args <- commandArgs(trailingOnly = TRUE); ds <- args[1]; B <- as.integer(args[2]); n.cores <- as.integer(args[3])
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
out.dir <- file.path(BENCH, "sweep", ds); dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
arm.sizes <- c(4, 6, 9)
runOne <- function(graph, spc, meta2, k, B, seed, n.cores, post) {
  model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
  w <- setNames(rep(1, nlevels(spc)), levels(spc))
  d <- cacoa:::estimateCellDensityGraph(graph, spc, w, smoothing = "knn", k = k, verbose = FALSE)
  Y <- logCounts(d$density.mat, d$sample.sizes)
  r <- cacoa:::diffCellDensityPermutations(Y, model, perm.method = "block", statistic = "t", n.permutations = B, n.cores = n.cores, seed = seed, moderate = FALSE)
  res <- list()
  res[["none"]] <- list(q = permutationFDR(r$stat.obs, r$stats.perm)$q, effect = r$effect)
  av <- averageOverNeighbours(r$stat.obs, r$stats.perm, graph, effect = r$effect)
  res[["averaged"]] <- list(q = permutationFDR(av$stat, av$stats.perm)$q, effect = av$effect)
  res
}
for (d in list.dirs(file.path(BENCH, "planted", ds), recursive = FALSE)) {
  tr <- readRDS(file.path(d, "truth.rds")); if (is.na(tr$type) || tr$frac < 1 || tr$size.class == "small") next
  for (n.arm in arm.sizes) {
    out <- file.path(out.dir, sprintf("%s_n%d.rds", tr$id, n.arm)); lock <- paste0(out, ".lock"); if (file.exists(out) || file.exists(lock)) next
    file.create(lock); t0 <- Sys.time()
    sub <- readRDS(file.path(d, "inputs.rds")); meta <- sub$sample.meta
    set.seed(300 + 10 * n.arm + tr$rep)
    keep <- unlist(lapply(levels(meta$condition), function(a) {   # n.arm samples per arm, allocated across the strata in proportion
      s <- rownames(meta)[meta$condition == a]; st <- as.character(meta[s, "stratum"]); tab <- table(st); n.take <- min(n.arm, length(s))
      alloc <- floor(n.take * tab / length(s)); rem <- n.take - sum(alloc); if (rem > 0) { o <- order(-(n.take * tab / length(s) - alloc)); alloc[o[seq_len(rem)]] <- alloc[o[seq_len(rem)]] + 1 }
      unlist(lapply(names(tab), function(l) { ss <- s[st == l]; if (alloc[[l]] > 0) sample(ss, alloc[[l]]) else character(0) })) }))
    if (length(keep) < 2 * n.arm) { message(tr$id, " n", n.arm, ": fewer samples than requested (", length(keep), ")") }
    inp <- subsetInputs(sub, samples = keep); spc <- droplevels(inp$sample.per.cell[V(inp$graph)$name]); cells <- V(inp$graph)$name
    meta2 <- data.frame(row.names = keep, cond = meta[keep, "condition"], stratum = meta[keep, "stratum"])
    ks <- c(rule = defaultNeighbourhoodSize(length(cells), length(keep)), fixed = 200L)
    res <- list()
    for (kn in names(ks)) { r <- runOne(inp$graph, spc, meta2, ks[[kn]], B, tr$rep, n.cores); for (p in names(r)) res[[sprintf("knn %s (k=%d) | %s", kn, ks[[kn]], p)]] <- r[[p]] }
    saveRDS(list(id = tr$id, n.arm = n.arm, samples = keep, res = res, n.cells = length(cells), time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), out); unlink(lock)
    message(sprintf("%s n=%d: %d cells, %d samples, done in %.0f s", tr$id, n.arm, length(cells), length(keep), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
message("done")
