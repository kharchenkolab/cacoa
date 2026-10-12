## Figures and tables of the null-relabeling study. Usage: Rscript analyze_null.R [out.dir]
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
dsets <- c(scc = "SCC (18 samples, 35.6k cells)", ms = "MS (21 samples, 48.9k cells)", pf = "PF (22 samples, 89.3k cells)")
# QQ points: expected vs observed -log10 p on a quantile grid (pooled over cells and relabelings)
qq.points <- function(p, n = 400) {
  p <- p[is.finite(p)]; if (!length(p)) return(NULL)
  m <- length(p); ranks <- unique(round(exp(seq(log(1), log(m), length.out = n))))
  ps <- sort(p)[ranks]
  data.frame(expected = -log10(ranks / (m + 1)), observed = -log10(pmax(ps, 1e-300)), n = m)
}
all.calls <- list(); all.qq <- list(); all.grp <- list(); all.grp.calls <- list(); times <- list()
for (ds in names(dsets)) {
  files <- list.files(file.path(S, Sys.getenv("CACOA_NULL_DIR", "bench/null3"), ds), pattern = "^r[0-9]+\\.rds$", full.names = TRUE)
  if (!length(files)) next
  pools <- list(); grp <- list()
  for (f in files) {
    r <- readRDS(f)
    pools[["density, heat kernel: parametric (moderated t)"]] <- c(pools[["density, heat kernel: parametric (moderated t)"]], r$heat$p.param)
    pools[["density, heat kernel: permutation"]] <- c(pools[["density, heat kernel: permutation"]], r$heat$p.perm)
    pools[["density, 1-hop counts (log): parametric (moderated t)"]] <- c(pools[["density, 1-hop counts (log): parametric (moderated t)"]], r$nb$p.param)
    pools[["density, 1-hop counts (log): permutation"]] <- c(pools[["density, 1-hop counts (log): permutation"]], r$nb$p.perm)
    if (!is.null(r$heat$p.pooled)) {
      pools[["density, heat kernel: pooled null"]] <- c(pools[["density, heat kernel: pooled null"]], r$heat$p.pooled)
      pools[["density, 1-hop counts (log): pooled null"]] <- c(pools[["density, 1-hop counts (log): pooled null"]], r$nb$p.pooled)
      pools[["expression shift: pooled null"]] <- c(pools[["expression shift: pooled null"]], r$shifts$p.pooled)
    }
    pools[["expression shift: permutation"]] <- c(pools[["expression shift: permutation"]], r$shifts$p.perm)
    pools[["expression shift: normal tail of z"]] <- c(pools[["expression shift: normal tail of z"]], r$shifts$p.z)
    pools[["expression shift: parametric F"]] <- c(pools[["expression shift: parametric F"]], r$shifts$p.F)
    all.calls[[length(all.calls) + 1]] <- cbind(dataset = ds, k = r$k, r$calls)
    for (g in names(r$groups)) { x <- r$groups[[g]]; x$k <- r$k; x$test <- g; x$dataset <- ds; grp[[length(grp) + 1]] <- x }
    times[[length(times) + 1]] <- data.frame(dataset = ds, k = r$k, seconds = r$time)
  }
  for (nm in names(pools)) { q <- qq.points(pools[[nm]]); if (!is.null(q)) all.qq[[length(all.qq) + 1]] <- cbind(dataset = ds, source = nm, q) }
  all.grp[[length(all.grp) + 1]] <- bind_rows(grp)
  message(ds, ": ", length(files), " relabelings")
}
qq <- do.call(rbind, all.qq); calls <- do.call(rbind, all.calls); grp <- bind_rows(all.grp); tms <- do.call(rbind, times)
qq$dataset <- factor(dsets[qq$dataset], levels = dsets); calls$dataset <- factor(dsets[calls$dataset], levels = dsets)
qq$kind <- ifelse(grepl("^density, heat", qq$source), "density: heat kernel", ifelse(grepl("^density, 1-hop", qq$source), "density: 1-hop counts", "expression shift"))
qq$source <- sub("^[^:]*: ", "", qq$source)
# --- figure 1: QQ of the per-cell p-values
p1 <- ggplot(qq, aes(expected, observed, colour = source)) + geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) +
  geom_line(linewidth = 0.8) + facet_grid(kind ~ dataset) + coord_cartesian(xlim = c(0, 7), ylim = c(0, 7)) +
  labs(x = "expected -log10 p (uniform)", y = "observed -log10 p", colour = NULL, title = "Per-cell p-values under condition relabeling (pooled over cells and relabelings)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(out.dir, "null_qq_cells.png"), p1, width = 11, height = 9, dpi = 110)
# --- figure 2: fraction of cells called per relabeling
calls$frac <- calls$n.called / pmax(calls$n.tested, 1)
calls$test <- factor(calls$test, levels = c("heat", "nb.log", "shifts"), labels = c("density: heat kernel", "density: 1-hop counts (log)", "expression shift"))
lv <- c(bh.param = "BH (parametric p)", bh.perm = "BH (permutation p)", bh.z = "BH (normal-tail z)", bh.F = "BH (parametric F)", bh.pooled = "BH (pooled p)",
        spatial.param = "spatial FDR (parametric)", spatial.perm = "spatial FDR (permutation)", spatial.pooled = "spatial FDR (pooled p)", nbhd.param = "neighbourhood (parametric)", nbhd.pooled = "neighbourhood (pooled p)",
        perm.fdr = "pooled permutation FDR", fdr.raw = "pooled FDR, raw t", fdr.sd = "pooled FDR, standardized (sd)", fdr.mad = "pooled FDR, standardized (mad)", fdr.tb = "pooled FDR, ter Braak scale", fdr.tbscale = "pooled FDR, ter Braak scale (label draws)", max.t = "max-t")
calls <- calls %>% filter(method %in% names(lv)); calls$method <- factor(lv[calls$method], levels = lv)
cs <- calls %>% group_by(dataset, test, method, alpha) %>% summarise(mean.frac = mean(frac), max.frac = max(frac), relabelings.with.calls = mean(n.called > 0), n = n(), .groups = "drop")
write.csv(cs, file.path(out.dir, "null_calls_summary.csv"), row.names = FALSE)
p2 <- ggplot(calls %>% filter(alpha == 0.05), aes(method, frac, colour = test)) + geom_jitter(width = 0.2, height = 0, size = 1, alpha = 0.6) +
  facet_grid(test ~ dataset, scales = "free_y") + coord_flip() + geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
  labs(x = NULL, y = "fraction of tested cells called at FDR 5 % (one point per relabeling)", title = "Calls under the null") + theme_bw(base_size = 10) + theme(legend.position = "none")
ggsave(file.path(out.dir, "null_calls_cells.png"), p2, width = 11, height = 7, dpi = 110)
# --- figure 3: per-group p-values (exact group permutation test; Cauchy combination of the permutation p; Cauchy of the parametric p)
grp$dataset <- factor(dsets[grp$dataset], levels = dsets)
grp$test <- factor(grp$test, levels = c("heat", "nb", "shifts"), labels = c("density: heat kernel", "density: 1-hop counts (log)", "expression shift"))
if (is.null(grp$p.cauchy.pooled)) grp$p.cauchy.pooled <- NA_real_
gl <- grp %>% select(dataset, test, k, group, p.perm, p.cauchy, p.cauchy.param, p.cauchy.pooled) %>%
  pivot_longer(c(p.perm, p.cauchy, p.cauchy.param, p.cauchy.pooled), names_to = "source", values_to = "p") %>% filter(is.finite(p)) %>%
  mutate(source = factor(source, levels = c("p.perm", "p.cauchy", "p.cauchy.param", "p.cauchy.pooled"), labels = c("group permutation test", "Cauchy (permutation p)", "Cauchy (parametric p)", "Cauchy (pooled p)")))
gq <- gl %>% group_by(dataset, test, source) %>% group_modify(~ qq.points(.x$p, n = 200)) %>% ungroup()
p3 <- ggplot(gq, aes(expected, observed, colour = source)) + geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = 2) + geom_point(size = 0.8) + geom_line() +
  facet_grid(test ~ dataset) + coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
  labs(x = "expected -log10 p", y = "observed -log10 p (one value per cell type and relabeling)", colour = NULL, title = "Per-cell-type p-values under relabeling") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out.dir, "null_qq_groups.png"), p3, width = 10, height = 8, dpi = 110)
gc.sum <- gl %>% group_by(dataset, test, source, k) %>% summarise(n.groups = n(), n.raw05 = sum(p < 0.05), n.bh05 = sum(p.adjust(p, "BH") < 0.05), .groups = "drop") %>%
  group_by(dataset, test, source) %>% summarise(relabelings = n(), groups = max(n.groups), frac.raw05 = sum(n.raw05) / sum(n.groups), relabelings.with.bh.call = mean(n.bh05 > 0), .groups = "drop")
write.csv(gc.sum, file.path(out.dir, "null_groups_summary.csv"), row.names = FALSE)
print(as.data.frame(cs %>% filter(alpha == 0.05) %>% select(dataset, test, method, mean.frac, relabelings.with.calls, n)), digits = 3)
print(as.data.frame(gc.sum), digits = 3)
print(tms %>% group_by(dataset) %>% summarise(n = n(), median.seconds = median(seconds)))
