## Planted expression shifts (planted_shift.R): recall of the shifted region at a realized false-call rate <= 10 % among
## the clean negatives (beyond two hops of the region), by unit and post-processing; calibration from null_shift.
## Usage: Rscript analyze_shift.R [out.dir] -> figures/shift_summary.csv, shift_grid.png, null_shift_calls.csv
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"
thr <- c(1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5)
rows <- list()
for (ds in c("scc", "pf", "ms")) for (d in list.dirs(file.path(S, "bench/shift", ds), recursive = FALSE)) {
  if (!file.exists(file.path(d, "shift.rds"))) next
  tr <- readRDS(file.path(d, "truth.rds")); sr <- readRDS(file.path(d, "shift.rds"))
  for (m in names(sr$res)) {
    q <- sr$res[[m]]$q; cells <- names(q)[is.finite(q)]; q <- q[cells]
    tp <- cells %in% tr$true.cells; ng <- cells %in% tr$negatives
    for (t in thr) { called <- q < t
      rows[[length(rows) + 1]] <- data.frame(dataset = ds, id = tr$id, what = tr$what, strength = tr$strength, rep = tr$rep, method = m, thr = t, n.called = sum(called), n.tested = length(cells),
        recall = mean(called[tp]), fdp.neg = if (any(called)) sum(called & ng) / sum(called) else 0) }
  }
  message(ds, " ", tr$id)
}
cv <- bind_rows(rows); write.csv(cv, file.path(out.dir, "shift_curves.csv"), row.names = FALSE)
rat <- cv %>% group_by(dataset, id, what, strength, method) %>% summarise(recall.fdp10 = max(c(0, recall[fdp.neg <= 0.10])), recall.nominal10 = recall[thr == 0.1], fdp.nominal10 = fdp.neg[thr == 0.1], .groups = "drop") %>%
  tidyr::separate(method, c("unit", "post"), sep = " \\| ", remove = FALSE)
sm <- rat %>% group_by(dataset, what, strength, unit, post) %>% summarise(across(c(recall.fdp10, recall.nominal10, fdp.nominal10), mean), n = n(), .groups = "drop")
write.csv(sm, file.path(out.dir, "shift_summary.csv"), row.names = FALSE)
options(width = 180)
for (st in sort(unique(sm$strength))) { cat("\n== strength", st, ": recall at realized fdp <= 0.10 (mean over reps), rows = unit x post, columns = planted region\n")
  w <- sm %>% filter(strength == st) %>% mutate(v = sprintf("%.2f", recall.fdp10)) %>% tidyr::pivot_wider(id_cols = c(dataset, unit, post), names_from = what, values_from = v)
  print(as.data.frame(w %>% arrange(dataset, post, unit))) }
cat("\n== realized fdp at the nominal 0.1 (mean over everything)\n"); print(as.data.frame(sm %>% group_by(dataset, unit, post) %>% summarise(fdp = mean(fdp.nominal10), rec = mean(recall.nominal10), .groups = "drop")), digits = 2)
sm$unit <- factor(sm$unit, levels = c("1-hop", "knn k=100", "knn k=200", "knn k=400"))
pg <- ggplot(sm, aes(unit, recall.fdp10, colour = post, group = post)) + geom_line() + geom_point(size = 2) + facet_grid(what + strength ~ dataset) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") + labs(x = NULL, y = "recall at realized false-call rate <= 10 %", colour = "post-processing", title = "Planted expression shifts: unit x post-processing")
ggsave(file.path(out.dir, "shift_grid.png"), pg, width = 10, height = 12, dpi = 110)
# calibration
calls <- list()
for (ds in c("scc", "pf", "ms")) for (f in list.files(file.path(S, "bench/null_shift", ds), pattern = "^r[0-9]+\\.rds$", full.names = TRUE)) { r <- readRDS(f); calls[[length(calls) + 1]] <- cbind(dataset = ds, k = r$k, r$calls) }
if (length(calls)) {
  calls <- bind_rows(calls); calls$frac <- calls$n.called / pmax(calls$n.tested, 1)
  cs <- calls %>% group_by(dataset, method, alpha) %>% summarise(mean.frac = mean(frac), max.frac = max(frac), relabelings.with.calls = mean(n.called > 0), n = n(), .groups = "drop")
  write.csv(cs, file.path(out.dir, "null_shift_calls.csv"), row.names = FALSE)
  cat("\n== shifts under the null: fraction of cells called at FDR 5 % (mean / max over relabelings)\n")
  print(as.data.frame(cs %>% filter(alpha == 0.05) %>% arrange(dataset, -max.frac)), digits = 2)
}
