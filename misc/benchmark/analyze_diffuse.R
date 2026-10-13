## Diffuse / gradient planted changes (planted_diffuse.R + grid_units.R diffuse mode): recall of the target cells
## (removal probability >= 0.4) at a realized false-call rate <= 10 % among the clean negatives (p < 0.05, beyond two
## hops of any target), and the rank correlation between the estimated effect and the planted removal probability.
## Usage: Rscript analyze_diffuse.R [out.dir]  -> figures/diffuse_summary.csv, diffuse_grid.png
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"
thr <- c(1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
rows <- list()
for (ds in c("scc", "pf", "ms")) for (d in list.dirs(file.path(S, "bench/diffuse", ds), recursive = FALSE)) {
  if (!file.exists(file.path(d, "grid.rds"))) next
  tr <- readRDS(file.path(d, "truth.rds")); gr <- readRDS(file.path(d, "grid.rds"))
  for (m in names(gr$res)) {
    q <- gr$res[[m]]$q; e <- gr$res[[m]]$effect; cells <- names(q)[is.finite(q)]; q <- q[cells]; e <- e[cells]
    tp <- cells %in% tr$true.cells; ng <- cells %in% tr$negatives; p <- tr$p.remove[cells]
    rho <- suppressWarnings(cor(-e, p, method = "spearman", use = "complete.obs"))
    for (t in thr) { called <- q < t
      rows[[length(rows) + 1]] <- data.frame(dataset = ds, id = tr$id, field = tr$field, rep = tr$rep, method = m, thr = t, n.called = sum(called),
        recall = mean(called[tp] & e[tp] < 0), fpr = mean(called[ng] & e[ng] < 0), fdp.neg = if (any(called)) sum(called & ng & e < 0) / sum(called) else 0, rho = rho) }
  }
  message(ds, " ", tr$id)
}
cv <- bind_rows(rows); write.csv(cv, file.path(out.dir, "diffuse_curves.csv"), row.names = FALSE)
rat <- cv %>% group_by(dataset, id, field, method) %>% summarise(recall.fdp10 = max(c(0, recall[fdp.neg <= 0.10])), recall.nominal10 = recall[thr == 0.1], fdp.nominal10 = fdp.neg[thr == 0.1], fpr.nominal10 = fpr[thr == 0.1], rho = rho[1], .groups = "drop") %>%
  tidyr::separate(method, c("unit", "weighting", "post"), sep = " \\| ", remove = FALSE)
sm <- rat %>% group_by(dataset, field, unit, post) %>% summarise(across(c(recall.fdp10, recall.nominal10, fdp.nominal10, fpr.nominal10, rho), mean), n = n(), .groups = "drop")
write.csv(sm, file.path(out.dir, "diffuse_summary.csv"), row.names = FALSE)
options(width = 180)
for (fl in c("gradient", "bump")) { cat("\n== ", fl, ": recall at realized fdp <= 0.10 among clean negatives (mean of 3 reps) | recall at nominal 0.1 | realized fdp at nominal | Spearman(effect, p)\n")
  print(as.data.frame(sm %>% filter(field == fl) %>% arrange(dataset, post, unit) %>% select(dataset, unit, post, recall.fdp10, recall.nominal10, fdp.nominal10, rho)), digits = 2) }
sm$family <- ifelse(grepl("^kernel", sm$unit), "kernel", "kNN"); sm$unit <- factor(sm$unit, levels = c("kernel b=10", "kernel b=30", "kernel b=100", "knn k=100", "knn k=200", "knn k=400"))
pg <- ggplot(sm, aes(unit, recall.fdp10, colour = post, group = post)) + geom_line() + geom_point(size = 2) + facet_grid(field ~ dataset) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") + labs(x = NULL, y = "recall at realized false-call rate <= 10 %", colour = "post-processing", title = "Diffuse and gradient changes")
ggsave(file.path(out.dir, "diffuse_grid.png"), pg, width = 10, height = 6, dpi = 110)
