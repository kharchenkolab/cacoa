## cacoa's cluster-free density test (BH default) per cell, for the head-to-head with Milo
suppressPackageStartupMessages(library(cacoa))
S <- "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad"
run <- function(con, cell.groups, sample.per.cell, meta, test, out, n.cores = 16) {
  cao <- Cacoa$new(con, sample.metadata = meta, cell.groups = cell.groups, sample.per.cell = sample.per.cell, n.cores = n.cores, verbose = FALSE)
  cao$setOptions(seed = 1); cao$setModel(stats::as.formula(paste("~", test)), test = test, verbose = FALSE)
  t0 <- Sys.time(); cao$estimateCellDensity(method = "graph", verbose = FALSE); t1 <- Sys.time()
  res <- list()
  for (m in c("bh", "max-t")) {
    r <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 999, adjust.method = m, verbose = FALSE)$diff$permutation
    res[[m]] <- data.frame(cell = names(r$raw), celltype = as.character(cao$cell.groups[names(r$raw)]), sample = as.character(cao$sample.per.cell[names(r$raw)]),
                           z.raw = unname(r$raw), z.adj = unname(r$adj[names(r$raw)]), p.adj = if (!is.null(r$p.adj)) unname(r$p.adj[names(r$raw)]) else NA, effect = unname(r$effect[names(r$raw)]),
                           stringsAsFactors = FALSE)
  }
  t2 <- Sys.time()
  cl <- cao$estimateCellLoadings(verbose = FALSE)
  saveRDS(list(bh = res[["bh"]], maxt = res[["max-t"]], loadings = cl$contrast$loadings$ref_centered, loadings.padj = cl$contrast$per_cell$padj,
               time = c(density = as.numeric(difftime(t1, t0, units = "secs")), diff.x2 = as.numeric(difftime(t2, t1, units = "secs")))), out)
  cat("saved", out, "\n")
}
d <- readRDS("/home/pkharchenko/cacoa-dev/test/sim_density.rds")
run(d$con, d$cell.groups, d$sample.per.cell, d$sample.meta, "treatment", file.path(S, "milo/cacoa_sim.rds"))
sc <- readRDS(file.path(S, "ex_scc/scc_refreshed.rds")); sm <- data.frame(row.names = names(sc$con$samples)); sm$condition <- factor(sub(".*_", "", rownames(sm)), levels = c("Normal", "Tumor")); sm$patient <- sub("_.*", "", rownames(sm))
run(sc$con, sc$cell.groups, sc$sample.per.cell, sm, "condition", file.path(S, "milo/cacoa_scc.rds"))
ms <- readRDS(file.path(S, "ex_ms/cao_step1.rds"))
run(ms$data.object, ms$cell.groups, ms$sample.per.cell, ms$sample.meta[, "diagnosis", drop = FALSE], "diagnosis", file.path(S, "milo/cacoa_ms.rds"))
