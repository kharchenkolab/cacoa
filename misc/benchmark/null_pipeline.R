## Calibration of the harmonized package pipeline on the balanced relabelings: the kNN density map (raw vs standardized
## pooled statistic, averaged post-processing) and the kNN shift map (standardized vs raw pooled statistic, averaged),
## through the class-level code path (estimateCellDensityGraph + diffCellDensityPermutations, clusterFreeExpressionShifts).
## Usage: Rscript null_pipeline.R <scc|pf|ms> <n.relabel> <B> <n.cores>  -> $S/bench/null_pipeline/<ds>.rds
args <- commandArgs(trailingOnly = TRUE); ds <- args[1]; n.relabel <- as.integer(args[2]); B <- as.integer(args[3]); n.cores <- as.integer(args[4])
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
out.dir <- file.path(BENCH, "null_pipeline"); dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
inp <- tm("load inputs", loadInputs(ds)); cs <- condSpec(ds); meta <- inp$sample.meta; g <- inp$graph; cells <- V(g)$name
cond0 <- relevel(factor(as.character(meta[[cs[1]]])), cs[2]); spc <- droplevels(inp$sample.per.cell[cells]); cg <- inp$cell.groups[cells]
k <- defaultNeighbourhoodSize(length(cells), nlevels(spc)); nn <- tm("kNN units", knnNeighbourhoodLists(g, k, n.cores = n.cores)); nn1 <- cacoa:::graphNeighbourLists0(g, cells)
w <- setNames(rep(1, nlevels(spc)), levels(spc)); dm <- cacoa:::estimateCellDensityGraph(g, spc, w, smoothing = "knn", k = k, verbose = FALSE, neighbourhoods = nn)
Y <- logCounts(dm$density.mat, dm$sample.sizes); codes <- as.integer(as.factor(cg)) - 1L
genes <- topGenes(inp$counts[, setdiff(colnames(inp$counts), c("Patient", "Tissue: 0=Normal, 1=Tumor"))], 1000); cm.t <- Matrix::t(inp$counts[cells, genes])
rows <- list()
for (r in seq_len(n.relabel)) {
  set.seed(1000 + r); lab <- rep(NA_character_, nrow(meta))
  for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
  meta2 <- data.frame(row.names = rownames(meta), cond = factor(lab, levels = c("A", "B")), stratum = cond0)
  model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum"); t0 <- Sys.time()
  for (ps in c("raw", "standardized")) for (av in c(FALSE, TRUE)) {
    d <- cacoa:::diffCellDensityPermutations(Y, model, perm.method = "freedman-lane", n.permutations = B, n.cores = n.cores, seed = r, pooled = TRUE, group = codes, average = av, nn.graph = nn1, return.sampled.stats = FALSE, pooled.statistic = ps)
    q <- permutationFDR(d$pooled.hist$s.obs, d$pooled.hist, alternative = "greater")$q
    rows[[length(rows) + 1]] <- data.frame(dataset = ds, relabel = r, map = "density", pooled = ps, averaged = av, n.tested = sum(is.finite(q)), called05 = sum(q < 0.05, na.rm = TRUE), called10 = sum(q < 0.1, na.rm = TRUE),
                                           p.pooled.q01 = mean(d$p.pooled < 0.01, na.rm = TRUE), p.pooled.q001 = mean(d$p.pooled < 0.001, na.rm = TRUE))
  }
  for (ps in c("standardized", "raw")) for (post in c("none", "averaged")) {
    s <- clusterFreeExpressionShifts(cm.t, spc, nn, model, model$meta, n.permutations = B, seed = r, adjust = TRUE, smooth = FALSE, n.cores = n.cores, verbose = FALSE, nns.graph = nn1, post = post, pooled.statistic = ps)
    rows[[length(rows) + 1]] <- data.frame(dataset = ds, relabel = r, map = "shift", pooled = ps, averaged = post == "averaged", n.tested = sum(is.finite(s$p.adj)), called05 = sum(s$p.adj < 0.05, na.rm = TRUE), called10 = sum(s$p.adj < 0.1, na.rm = TRUE),
                                           p.pooled.q01 = mean(s$p.pooled < 0.01, na.rm = TRUE), p.pooled.q001 = mean(s$p.pooled < 0.001, na.rm = TRUE))
  }
  message(sprintf("[%s] relabeling %d done in %.0f s", ds, r, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  saveRDS(do.call(rbind, rows), file.path(out.dir, paste0(ds, ".rds")))
}
message("done")
