## Calibration of the unit x post-processing grid on the balanced relabelings (grid_units.R null mode). Usage: Rscript analyze_null_grid.R [out.dir]
suppressPackageStartupMessages(library(dplyr))
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"; dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
calls <- list()
for (ds in c("scc", "pf", "ms")) for (f in list.files(file.path(S, "bench/null_grid", ds), pattern = "^r[0-9]+\\.rds$", full.names = TRUE)) { r <- readRDS(f); calls[[length(calls) + 1]] <- cbind(dataset = ds, k = r$k, r$calls) }
calls <- bind_rows(calls); calls$frac <- calls$n.called / pmax(calls$n.tested, 1)
cs <- calls %>% group_by(dataset, method, alpha) %>% summarise(mean.frac = mean(frac), max.frac = max(frac), relabelings.with.calls = mean(n.called > 0), n = n(), .groups = "drop")
write.csv(cs, file.path(out.dir, "null_grid_calls.csv"), row.names = FALSE)
options(width = 170); cat("== grid under the null: fraction of cells called at FDR 5 % (mean / max over relabelings); only rows with any call\n")
print(as.data.frame(cs %>% filter(alpha == 0.05, max.frac > 0) %>% arrange(dataset, -mean.frac)), digits = 2)
cat("\nrelabelings per dataset:", paste(sapply(split(calls$k, calls$dataset), function(k) length(unique(k))), collapse = " / "), "\n")
