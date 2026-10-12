## Per-cell null distributions for the tail-approximation study: one stratified null relabeling, a random subset of
## cells, B permutations with the full permutations x cells matrices kept (shift F, heat-kernel t, 1-hop log-count t).
## Usage: Rscript tail_collect.R <scc|ms|pf> <n.cells> <B> <n.cores> [relabel.seed]
args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]; n.cells <- as.integer(args[2]); B <- as.integer(args[3]); n.cores <- as.integer(args[4]); k <- if (length(args) > 4) as.integer(args[5]) else 1L
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
out.dir <- file.path(BENCH, "tail"); dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
inp <- tm("load inputs", loadInputs(ds))
cs <- condSpec(ds); meta <- inp$sample.meta
cond0 <- factor(as.character(meta[[cs[1]]])); cond0 <- relevel(cond0, cs[2])
spc <- inp$sample.per.cell; cells.all <- V(inp$graph)$name
dm <- tm("density matrices", densityMats(inp$graph, spc, n.cores = n.cores))
nn <- tm("neighbour lists", neighbourLists(inp$graph))
genes <- topGenes(inp$counts, 1000)
cm.t <- Matrix::t(inp$counts[cells.all, genes])
set.seed(1000 + k)
lab <- rep(NA_character_, nrow(meta))
for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
meta2 <- data.frame(row.names = rownames(meta), cond = factor(lab, levels = c("A", "B")), stratum = cond0)
model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
# cells stratified by type (at least 20 per type where available)
cg <- inp$cell.groups
per <- pmax(20L, round(n.cells * as.numeric(table(cg)) / length(cg)))
cells <- unlist(lapply(levels(cg), function(l) { i <- names(cg)[cg == l]; sample(i, min(length(i), per[match(l, levels(cg))])) }))
message(sprintf("%d cells from %d types", length(cells), nlevels(cg)))
sh <- tm("shifts", clusterFreeExpressionShifts(cm.t, spc, nn[cells], model, model$meta, dist = "cor", n.permutations = B, seed = k, adjust = FALSE, smooth = FALSE,
                                               n.cores = n.cores, verbose = FALSE, return.perm = TRUE))
dens <- function(Y) {
  r <- cacoa:::diffCellDensityPermutations(Y, model, perm.method = "freedman-lane", statistic = "t", n.permutations = B, n.cores = n.cores, seed = k, moderate = TRUE)
  Yc <- Y - apply(Y, 1, min)   # per-cell concentration of the density over samples (log counts are shifted to min 0)
  ystats <- data.frame(row.names = rownames(Y), n.eff = rowSums(Yc)^2 / rowSums(Yc^2), max.share = apply(Yc, 1, max) / rowSums(Yc), n.nonzero = rowSums(Y > min(Y) + 1e-12))
  list(stat = r$stat.obs, perm = r$stats.perm, p.perm = r$p.perm, p.param = r$p.param, t = r$t, df = r$df, z = r$z.score, ystats = ystats)
}
heat <- tm("heat", dens(dm$heat[cells, ]))
nb <- tm("nb", dens(logCounts(dm$nb, dm$sizes)[cells, ]))
res <- list(ds = ds, k = k, B = B, cells = cells, groups = cg[cells], labels = meta2$cond, n.samples = sh$n.samples,
            shifts = list(stat = sh$stat, perm = sh$stats.perm, p.perm = sh$p.value, p.param = sh$p.param, z = sh$z.std, n.perm = sh$n.perm),
            heat = heat, nb = nb)
saveRDS(res, file.path(out.dir, sprintf("%s_r%d_B%d.rds", ds, k, B)))
message("done")
