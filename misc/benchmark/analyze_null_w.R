## Calibration of the weighted tests and the inheritance on the balanced relabelings (weighted_units.R null mode).
## Usage: Rscript analyze_null_w.R [out.dir] -> figures/null_w_calls.csv, null_w_pooled_p.csv
suppressPackageStartupMessages(library(dplyr))
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
calls <- list(); pq <- list()
for (ds in c("scc", "pf", "ms")) for (f in list.files(file.path(S, "bench/null_w", ds), pattern = "^r[0-9]+\\.rds$", full.names = TRUE)) {
  r <- readRDS(f); calls[[length(calls) + 1]] <- cbind(dataset = ds, k = r$k, r$calls)
  for (m in names(r$p.quantiles)) pq[[length(pq) + 1]] <- data.frame(dataset = ds, k = r$k, method = m, n = r$p.quantiles[[m]]$n, prob = as.numeric(sub("%", "", names(r$p.quantiles[[m]]$q))) / 100, q = as.numeric(r$p.quantiles[[m]]$q))
}
calls <- bind_rows(calls); calls$frac <- calls$n.called / pmax(calls$n.tested, 1)
cs <- calls %>% group_by(dataset, method, alpha) %>% summarise(mean.frac = mean(frac), max.frac = max(frac), relabelings.with.calls = mean(n.called > 0), n = n(), .groups = "drop")
write.csv(cs, file.path(out.dir, "null_w_calls.csv"), row.names = FALSE)
pq <- bind_rows(pq); ps <- pq %>% group_by(dataset, method, prob) %>% summarise(observed = mean(q), .groups = "drop") %>% mutate(ratio = prob / observed)   # observed quantile of the pooled p at the nominal probability: ratio > 1 = anti-conservative
write.csv(ps, file.path(out.dir, "null_w_pooled_p.csv"), row.names = FALSE)
options(width = 170)
cat("== fraction of cells called at FDR 5 % under the null (mean / max over relabelings)\n"); print(as.data.frame(cs %>% filter(alpha == 0.05) %>% arrange(dataset, method)), digits = 2)
cat("\n== pooled p of the weighted tests: nominal probability / observed quantile (1 = calibrated, > 1 = anti-conservative)\n")
print(as.data.frame(ps %>% filter(prob %in% c(1e-4, 1e-3, 1e-2, 0.05)) %>% select(-observed) %>% tidyr::pivot_wider(names_from = prob, values_from = ratio)), digits = 2)
