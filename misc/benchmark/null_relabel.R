## Null calibration of the per-cell tests by relabeling the condition across samples.
## Usage: Rscript null_relabel.R <scc|ms|pf> <n.relabel> <B> <n.cores> [n.genes]
## Per relabeling: heat-kernel density test, 1-hop neighbourhood log-count density test, cluster-free expression shifts,
## every adjustment, and the per-group Cauchy summaries; results in $S/bench/null/<ds>/r<k>.rds
args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]; n.relabel <- as.integer(args[2]); B <- as.integer(args[3]); n.cores <- as.integer(args[4]); n.genes <- if (length(args) > 4) as.integer(args[5]) else 1000
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
out.dir <- file.path(BENCH, Sys.getenv("CACOA_NULL_DIR", "null"), ds); dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
inp <- tm("load inputs", loadInputs(ds))
cs <- condSpec(ds); meta <- inp$sample.meta
cond0 <- factor(as.character(meta[[cs[1]]])); cond0 <- relevel(cond0, cs[2])
spc <- inp$sample.per.cell; cells <- V(inp$graph)$name
dm <- tm("density matrices", densityMats(inp$graph, spc, n.cores = n.cores))
nn <- tm("neighbour lists", neighbourLists(inp$graph))
genes <- topGenes(inp$counts, n.genes)
cm.t <- Matrix::t(inp$counts[cells, genes])
saveRDS(list(cells = cells, cell.groups = inp$cell.groups, sample.per.cell = spc, cond = cond0, sizes = dm$sizes, degree = degree(inp$graph)), file.path(out.dir, "meta.rds"))
for (k in seq_len(n.relabel)) {
  f <- file.path(out.dir, sprintf("r%03d.rds", k))
  if (file.exists(f)) next
  set.seed(1000 + k)
  # arms balanced within the true condition: every stratum is split in half (the odd sample at random), so that the
  # real condition effect is the same in both arms and only the arm labels are random; the model blocks the
  # permutations on the true condition for the same reason (the observed labeling is then exchangeable with them)
  lab <- rep(NA_character_, nrow(meta))
  for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
  meta2 <- data.frame(row.names = rownames(meta), cond = factor(lab, levels = c("A", "B")), stratum = cond0)
  model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
  t0 <- Sys.time()
  groups <- inp$cell.groups
  heat <- runDensity(dm$heat, model, nn, B = B, n.cores = n.cores, seed = k, label = "heat", groups = groups)
  nb <- runDensity(dm$nb, model, nn, B = B, n.cores = n.cores, seed = k, label = "nb.log", transform = "log", sizes = dm$sizes, groups = groups)
  sh <- runShifts(cm.t, spc, nn, model, B = B, n.cores = n.cores, seed = k, groups = groups)
  gs <- list(heat = heat$groups, nb = nb$groups, shifts = sh$groups)
  calls <- rbind(cbind(test = "heat", callFractions(heat$q)), cbind(test = "nb.log", callFractions(nb$q)), cbind(test = "shifts", callFractions(sh$q)))
  saveRDS(list(k = k, labels = meta2$cond, heat = heat, nb = nb, shifts = sh, groups = gs, calls = calls, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), f)
  message(sprintf("[%s] relabeling %d done in %.0f s; calls at 5%%: %s", ds, k, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  paste(sprintf("%s/%s %d", calls$test[calls$alpha == 0.05], calls$method[calls$alpha == 0.05], calls$n.called[calls$alpha == 0.05]), collapse = ", ")))
}
message("done")
