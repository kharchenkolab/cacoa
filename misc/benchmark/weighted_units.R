## Precision-weighted permutation tests on cacoa's per-cell units (heat-kernel soft counts and 1-hop counts), with and
## without neighbourhood inheritance: the prototype of observation weights in the engine, on the planted sets (power)
## and on the balanced relabelings (calibration).
## Usage: Rscript weighted_units.R <scc|pf|ms> planted|null <B> <n.cores> [n.relabel]
## planted: <set>/hw.rds (q and effect per method); null: $S/bench/null_w/<ds>/r<k>.rds (calls per method, pooled p)
args <- commandArgs(trailingOnly = TRUE); ds <- args[1]; mode <- args[2]; B <- as.integer(args[3]); n.cores <- as.integer(args[4]); n.relabel <- if (length(args) > 4) as.integer(args[5]) else 15
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
suppressPackageStartupMessages(library(edgeR))
runUnits <- function(dm, graph, cond, model, B, seed, n.cores) {
  cells <- V(graph)$name; A <- as_adjacency_matrix(graph, sparse = TRUE); A@x[] <- 1
  P <- modelPermutations(model, scheme = "block", n.permutations = B, seed = seed)$P
  res <- list()
  for (u in c("heat", "nb")) {
    v <- if (u == "heat") linearWeights(dm$heat, cond) else softVoom(dm$nb, dm$sizes, cond)   # kernel density: linear-scale weights; counts: voom
    pw <- permWLS(v$E, v$w, v$x, P); qw <- permWLSq(pw)
    pu <- permWLS(v$E, v$w * 0 + 1, v$x, P); qu <- permWLSq(pu)
    lab <- if (u == "heat") "heat" else "1-hop"
    res[[paste0(lab, " + weighted perm: perm.fdr")]] <- list(q = qw$perm.fdr, effect = pw$effect, p.pooled = pooledPermutationP(pw$t.obs, pw$t.perm))
    res[[paste0(lab, " + weighted perm: fdr.sd")]] <- list(q = qw$fdr.sd, effect = pw$effect)
    res[[paste0(lab, " + unweighted perm: perm.fdr")]] <- list(q = qu$perm.fdr, effect = pu$effect, p.pooled = pooledPermutationP(pu$t.obs, pu$t.perm))
    for (m in c("perm.fdr", "fdr.sd")) { h <- inheritQ(res[[paste0(lab, " + weighted perm: ", m)]]$q, pw$effect, A, cells); res[[paste0(lab, " + weighted perm: ", m, " + inheritance")]] <- h }
    h <- inheritQ(qu$perm.fdr, pu$effect, A, cells); res[[paste0(lab, " + unweighted perm: perm.fdr + inheritance")]] <- h
    rm(pw, pu); gc()
  }
  res
}
if (mode == "planted") {
  root <- file.path(BENCH, "planted", ds)
  for (d in list.dirs(root, recursive = FALSE)) {
    out <- file.path(d, "hw.rds"); lock <- file.path(d, "hw.lock"); if (file.exists(out) || file.exists(lock)) next
    file.create(lock); t0 <- Sys.time()
    tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds")); spc <- sub$sample.per.cell
    dm <- densityMats(sub$graph, spc, n.cores = n.cores)
    meta2 <- data.frame(row.names = rownames(sub$sample.meta), cond = sub$sample.meta$condition, stratum = sub$sample.meta$stratum)
    model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
    res <- runUnits(dm, sub$graph, setNames(sub$sample.meta$condition, rownames(sub$sample.meta)), model, B, tr$rep, n.cores)
    saveRDS(list(id = tr$id, res = res, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), out); unlink(lock)
    message(tr$id, " done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  }
} else {
  out.dir <- file.path(BENCH, "null_w", ds); dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  inp <- tm("load inputs", loadInputs(ds)); cs <- condSpec(ds); meta <- inp$sample.meta
  cond0 <- relevel(factor(as.character(meta[[cs[1]]])), cs[2]); spc <- inp$sample.per.cell
  dm <- tm("density matrices", densityMats(inp$graph, spc, n.cores = n.cores))
  for (k in seq_len(n.relabel)) {
    f <- file.path(out.dir, sprintf("r%03d.rds", k)); lock <- paste0(f, ".lock"); if (file.exists(f) || file.exists(lock)) next
    file.create(lock); t0 <- Sys.time(); set.seed(1000 + k)
    lab <- rep(NA_character_, nrow(meta))
    for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
    meta2 <- data.frame(row.names = rownames(meta), cond = factor(lab, levels = c("A", "B")), stratum = cond0)
    model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
    res <- runUnits(dm, inp$graph, setNames(meta2$cond, rownames(meta2)), model, B, k, n.cores)
    calls <- do.call(rbind, lapply(names(res), function(m) do.call(rbind, lapply(c(0.05, 0.1), function(a) data.frame(method = m, alpha = a, n.tested = sum(is.finite(res[[m]]$q)), n.called = sum(res[[m]]$q < a, na.rm = TRUE))))))
    pq <- lapply(res[!sapply(lapply(res, `[[`, "p.pooled"), is.null)], function(r) { p <- r$p.pooled[is.finite(r$p.pooled)]; list(n = length(p), q = quantile(p, c(1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.25, 0.5))) })
    saveRDS(list(k = k, labels = meta2$cond, calls = calls, p.quantiles = pq, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), f); unlink(lock)
    message(sprintf("[%s] relabeling %d done in %.0f s; calls at 5%%: %s", ds, k, as.numeric(difftime(Sys.time(), t0, units = "secs")), paste(sprintf("%s %d", calls$method[calls$alpha == 0.05], calls$n.called[calls$alpha == 0.05]), collapse = ", ")))
  }
}
message("done")
