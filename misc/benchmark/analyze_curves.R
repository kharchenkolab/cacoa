## Recall against the REALIZED false-call rate (other cells called down beyond two hops of any planted cell), sweeping
## the q-value threshold of every method on the planted sets, so that methods are compared at matched achieved error.
## Usage: Rscript analyze_curves.R [out.dir]   -> figures/planted_curves.png, planted_curves.csv, planted_recall_at_fdp.csv
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"
thr <- c(1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
rows <- list()
for (ds in c("scc", "pf", "ms")) for (d in list.dirs(file.path(S, "bench/planted", ds), recursive = FALSE)) {
  tr <- readRDS(file.path(d, "truth.rds")); if (is.na(tr$type)) next
  g <- readRDS(file.path(d, "inputs.rds"))$graph; cells.g <- igraph::V(g)$name
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE); A@x[] <- 1
  v <- as.numeric(cells.g %in% tr$true.cells); h1 <- as.numeric(A %*% v) > 0; h2 <- as.numeric(A %*% h1) > 0
  far <- setNames(!h2 & v == 0, cells.g); tp <- setNames(v == 1, cells.g)
  methods <- list()
  if (file.exists(file.path(d, "cacoa.rds"))) { ca <- readRDS(file.path(d, "cacoa.rds"))
    for (tst in c("heat", "nb")) { x <- ca[[tst]]; for (m in names(x$q)) methods[[paste(if (tst == "heat") "heat" else "1-hop log", m, sep = ": ")]] <- list(q = setNames(x$q[[m]], rownames(x$q)), effect = x$effect[rownames(x$q)]) } }
  if (file.exists(file.path(d, "milo.rds"))) { mi <- readRDS(file.path(d, "milo.rds"))$cells; methods[["Milo (graph, spatial FDR)"]] <- list(q = setNames(mi$milo.fdr, mi$cell), effect = setNames(mi$milo.logFC, mi$cell)) }
  if (file.exists(file.path(d, "decomp.rds"))) { dc <- readRDS(file.path(d, "decomp.rds"))$res; names(dc) <- sub("^Milo nhoods", "kNN nhoods", names(dc)); for (m in names(dc)) methods[[paste0("decomp: ", m)]] <- dc[[m]] }
  if (file.exists(file.path(d, "hw.rds"))) { hw <- readRDS(file.path(d, "hw.rds"))$res; for (m in names(hw)) methods[[paste0("weighted: ", m)]] <- hw[[m]][c("q", "effect")] }
  for (m in names(methods)) {
    q <- methods[[m]]$q; e <- methods[[m]]$effect; cells <- names(q)[is.finite(q)]; q <- q[cells]; e <- e[cells]
    tpc <- tp[cells]; tpc[is.na(tpc)] <- FALSE; farc <- far[cells]; farc[is.na(farc)] <- FALSE
    for (t in thr) { called <- q < t
      rows[[length(rows) + 1]] <- data.frame(dataset = ds, id = tr$id, type = tr$type, size.class = tr$size.class, frac = tr$frac, method = m, thr = t, n.called = sum(called),
                                             recall = mean(called[tpc] & e[tpc] < 0), fdp2 = if (any(called)) sum(called & farc & e < 0) / sum(called) else 0) }
  }
  message(ds, " ", tr$id)
}
cv <- bind_rows(rows); write.csv(cv, file.path(out.dir, "planted_curves.csv"), row.names = FALSE)
# recall at a realized false-call rate: the largest recall over thresholds whose realized fdp2 is at or below the target (per set), then averaged
rat <- cv %>% group_by(dataset, id, type, size.class, frac, method) %>% summarise(recall.fdp05 = max(c(0, recall[fdp2 <= 0.05])), recall.fdp10 = max(c(0, recall[fdp2 <= 0.10])), recall.nominal10 = recall[thr == 0.1], fdp2.nominal10 = fdp2[thr == 0.1], .groups = "drop") %>%
  group_by(dataset, method, size.class, type, frac) %>% summarise(across(c(recall.fdp05, recall.fdp10, recall.nominal10, fdp2.nominal10), mean), n = n(), .groups = "drop")
write.csv(rat, file.path(out.dir, "planted_recall_at_fdp.csv"), row.names = FALSE)
main <- c("heat: perm.fdr" = "heat kernel, pooled FDR (default)", "weighted: heat + weighted perm: perm.fdr" = "heat kernel, weighted permutation", "weighted: heat + weighted perm: perm.fdr + inheritance" = "heat kernel, weighted permutation + inheritance",
          "decomp: heat + inheritance: perm.fdr" = "heat kernel, pooled FDR + inheritance", "weighted: 1-hop + weighted perm: perm.fdr" = "1-hop counts, weighted permutation", "weighted: 1-hop + weighted perm: perm.fdr + inheritance" = "1-hop counts, weighted permutation + inheritance",
          "decomp: kNN nhoods + voom-weighted perm (libsize): perm.fdr" = "kNN nhoods (~95 cells), weighted permutation + inheritance", "decomp: kNN nhoods + edgeR QL (libsize): bh" = "kNN nhoods, edgeR QL, BH + inheritance", "Milo (graph, spatial FDR)" = "Milo")
cm <- cv %>% filter(frac == 1, method %in% names(main)) %>% mutate(method = factor(main[method], levels = main)) %>% group_by(dataset, size.class, type, method, thr) %>% summarise(recall = mean(recall), fdp2 = mean(fdp2), .groups = "drop") %>%
  mutate(panel = paste0(size.class, ": ", type))
p <- ggplot(cm, aes(fdp2, recall, colour = method)) + geom_path() + geom_point(size = 1) + facet_grid(dataset ~ panel) + geom_vline(xintercept = 0.1, linetype = 2, colour = "grey40") + coord_cartesian(xlim = c(0, 0.6)) +
  labs(x = "realized false-call rate (other cells called down beyond 2 hops of the planted cells)", y = "recall of the planted cells (100 % removal)", colour = NULL, title = "Recall at matched realized error (threshold sweep)") + theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(out.dir, "planted_curves.png"), p, width = 12, height = 9, dpi = 110)
cat("\n== recall at realized fdp2 <= 0.10 (and <= 0.05), 100 % removal, mean over types\n")
print(as.data.frame(rat %>% filter(frac == 1) %>% group_by(dataset, method) %>% summarise(r10 = mean(recall.fdp10), r05 = mean(recall.fdp05), r.nom = mean(recall.nominal10), fdp.nom = mean(fdp2.nominal10), .groups = "drop") %>% arrange(dataset, -r10)), digits = 2)
