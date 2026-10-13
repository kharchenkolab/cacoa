## Sample-number sweep (sample_sweep.R): recall of the planted type at a realized false-call rate <= 10 % (far negatives
## beyond two hops, called down) against the number of samples per arm, kNN unit sized by the rule vs fixed k = 200,
## with and without the averaged statistic. Usage: Rscript analyze_sweep.R [out.dir] -> figures/sweep_summary.csv, sweep.png
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"
thr <- c(1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
rows <- list()
for (ds in c("scc", "pf", "ms")) for (f in list.files(file.path(S, "bench/sweep", ds), pattern = "\\.rds$", full.names = TRUE)) {
  sw <- readRDS(f); tr <- readRDS(file.path(S, "bench/planted", ds, sw$id, "truth.rds")); g <- readRDS(file.path(S, "bench/planted", ds, sw$id, "inputs.rds"))$graph
  cells.g <- igraph::V(g)$name; A <- igraph::as_adjacency_matrix(g, sparse = TRUE); A@x[] <- 1
  v <- as.numeric(cells.g %in% tr$true.cells); h1 <- as.numeric(A %*% v) > 0; h2 <- as.numeric(A %*% h1) > 0; far <- setNames(!h2 & v == 0, cells.g); tp <- setNames(v == 1, cells.g)
  for (m in names(sw$res)) {
    q <- sw$res[[m]]$q; e <- sw$res[[m]]$effect; cells <- names(q)[is.finite(q)]; q <- q[cells]; e <- e[cells]; tpc <- tp[cells]; farc <- far[cells]
    for (t in thr) { called <- q < t
      rows[[length(rows) + 1]] <- data.frame(dataset = ds, id = sw$id, type = tr$type, size.class = tr$size.class, n.arm = sw$n.arm, n.samples = length(sw$samples), method = m, thr = t,
        recall = mean(called[tpc] & e[tpc] < 0), fdp2 = if (any(called)) sum(called & farc & e < 0) / sum(called) else 0) }
  }
}
cv <- bind_rows(rows)
rat <- cv %>% group_by(dataset, id, size.class, n.arm, n.samples, method) %>% summarise(recall.fdp10 = max(c(0, recall[fdp2 <= 0.10])), recall.nominal10 = recall[thr == 0.1], fdp2.nominal10 = fdp2[thr == 0.1], .groups = "drop") %>%
  tidyr::separate(method, c("unit", "post"), sep = " \\| ", remove = FALSE)
sm <- rat %>% group_by(dataset, size.class, n.arm, unit, post) %>% summarise(across(c(recall.fdp10, recall.nominal10, fdp2.nominal10, n.samples), mean), n = n(), .groups = "drop")
write.csv(sm, file.path(out.dir, "sweep_summary.csv"), row.names = FALSE)
options(width = 180); cat("== recall at realized fdp2 <= 0.10 by samples per arm (mean over 3 reps); unit = kNN sized by the rule (k) or fixed 200\n")
print(as.data.frame(sm %>% mutate(v = sprintf("%.2f", recall.fdp10)) %>% tidyr::pivot_wider(id_cols = c(dataset, size.class, unit, post), names_from = n.arm, values_from = v) %>% arrange(dataset, size.class, post, unit)))
cat("\n== realized fdp2 at the nominal 0.1\n"); print(as.data.frame(sm %>% group_by(dataset, n.arm, post) %>% summarise(fdp = round(mean(fdp2.nominal10), 3), .groups = "drop") %>% tidyr::pivot_wider(names_from = n.arm, values_from = fdp)))
sm$unit.kind <- ifelse(grepl("rule", sm$unit), "k by rule", "k = 200")
pg <- ggplot(sm, aes(n.arm, recall.fdp10, colour = post, linetype = unit.kind)) + geom_line() + geom_point(size = 2) + facet_grid(size.class ~ dataset) + theme_bw() + theme(legend.position = "bottom") +
  labs(x = "samples per arm", y = "recall at realized false-call rate <= 10 % (100 % removal)", colour = "post-processing", linetype = NULL, title = "Sample-number sweep, kNN count unit")
ggsave(file.path(out.dir, "sweep.png"), pg, width = 9, height = 6, dpi = 110)
