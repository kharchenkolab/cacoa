## Figures and tables of the semi-synthetic power study. Usage: Rscript analyze_semi.R [out.dir]
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
alphas <- c(0.05, 0.1)
rows <- list(); grows <- list()
for (ds in c("scc", "pf", "scc_paired", "pf_paired")) for (d in list.dirs(file.path(S, "bench/semi", ds), recursive = FALSE)) {
  if (!file.exists(file.path(d, "cacoa.rds"))) next
  tr <- readRDS(file.path(d, "truth.rds")); ca <- readRDS(file.path(d, "cacoa.rds"))
  regime <- if (grepl("paired", ds)) "paired" else "unpaired"; base <- sub("_paired", "", ds)
  true.cells <- tr$true.cells; adjacent <- setdiff(ca$adjacent, true.cells)
  score <- function(q, effect, method, cells) {
    ok <- is.finite(q); cells <- cells[ok]; q <- q[ok]; effect <- effect[ok]
    tp.set <- cells %in% true.cells; near <- cells %in% adjacent
    do.call(rbind, lapply(alphas, function(a) {
      called <- q < a
      data.frame(method = method, alpha = a, n.tested = length(q), n.true.tested = sum(tp.set),
                 recall = mean(called[tp.set] & effect[tp.set] < 0), n.called = sum(called),
                 fdp.strict = if (any(called)) mean(!tp.set[called]) else NA_real_,
                 fdp.lenient = if (any(called)) mean(!tp.set[called] & !near[called]) else NA_real_,
                 frac.other.called = mean(called[!tp.set & !near]))
    }))
  }
  res <- list()
  for (tst in c("heat", "nb")) { x <- ca[[tst]]; cells <- rownames(x$q)
    for (m in names(x$q)) res[[length(res) + 1]] <- score(x$q[[m]], x$effect[cells], paste(if (tst == "heat") "heat" else "1-hop log", m, sep = ": "), cells) }
  if (file.exists(file.path(d, "milo.rds"))) {
    mi <- readRDS(file.path(d, "milo.rds"))$cells
    res[[length(res) + 1]] <- score(mi$milo.fdr, mi$milo.logFC, "Milo (graph, spatial FDR)", mi$cell)
  }
  r <- do.call(rbind, res); r$dataset <- base; r$regime <- regime; r$id <- tr$id; r$type <- tr$type; r$size.class <- tr$size.class; r$frac <- tr$frac; r$rep <- tr$rep
  rows[[length(rows) + 1]] <- r
  for (tst in c("heat", "nb")) { g <- ca$groups[[tst]]; if (is.null(g)) next
    isp <- !is.na(tr$type) & g$group == tr$type
    grows[[length(grows) + 1]] <- data.frame(dataset = base, regime = regime, id = tr$id, type = tr$type, size.class = tr$size.class, frac = tr$frac, rep = tr$rep, test = tst,
      planted.sig.perm = if (any(isp)) g$significant[isp] & g$direction[isp] == "down" else NA,
      planted.p.perm = if (any(isp)) g$p.perm[isp] else NA_real_, planted.p.cauchy = if (any(isp)) g$p.cauchy[isp] else NA_real_,
      other.sig.perm = sum(g$significant[g$group != tr$type], na.rm = TRUE),
      other.sig.cauchy = sum(p.adjust(g$p.cauchy, "BH")[g$group != tr$type] < 0.05, na.rm = TRUE),
      other.sig.cauchy.param = sum(p.adjust(g$p.cauchy.param, "BH")[g$group != tr$type] < 0.05, na.rm = TRUE), n.groups = sum(is.finite(g$p.perm))) }
}
res <- do.call(rbind, rows); gres <- do.call(rbind, grows)
write.csv(res, file.path(out.dir, "semi_cells_all.csv"), row.names = FALSE); write.csv(gres, file.path(out.dir, "semi_groups_all.csv"), row.names = FALSE)
res$dataset <- factor(paste(res$dataset, res$regime), levels = c("scc unpaired", "scc paired", "pf unpaired", "pf paired"),
                      labels = c("SCC normals, unpaired (4 + 5)", "SCC normals, paired (9 donors)", "PF controls, unpaired (5 + 5)", "PF controls, paired (10 donors)"))
res$size.class <- factor(res$size.class, levels = c("null", "small", "medium", "large"))
nulls <- res %>% filter(size.class == "null") %>% group_by(dataset, method, alpha) %>% summarise(mean.frac.called = mean(n.called / n.tested), max.frac.called = max(n.called / n.tested), n = n(), .groups = "drop")
write.csv(nulls, file.path(out.dir, "semi_null_calls.csv"), row.names = FALSE)
res <- res %>% filter(size.class != "null")
sm <- res %>% group_by(dataset, method, alpha, size.class, type, frac) %>% summarise(recall = mean(recall), fdp.strict = mean(fdp.strict, na.rm = TRUE), fdp.lenient = mean(fdp.lenient, na.rm = TRUE), frac.other = mean(frac.other.called), n = n(), .groups = "drop")
write.csv(sm, file.path(out.dir, "semi_cells_summary.csv"), row.names = FALSE)
keep <- c("heat: bh.param", "heat: perm.fdr", "heat: max.t", "1-hop log: bh.param", "1-hop log: perm.fdr", "1-hop log: nbhd.param", "Milo (graph, spatial FDR)")
sm2 <- sm %>% filter(alpha == 0.1, method %in% keep)
sm2$method <- factor(sm2$method, levels = keep, labels = c("heat density, BH", "heat density, pooled perm. FDR", "heat density, max-t", "1-hop log counts, BH", "1-hop log counts, pooled perm. FDR", "1-hop log counts, index-cell BH", "Milo (same graph, spatial FDR)"))
sm2$panel <- sprintf("%s: %s", sm2$size.class, sm2$type)
p1 <- ggplot(sm2, aes(frac, recall, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) + scale_x_continuous(breaks = c(0.5, 0.8, 1)) +
  labs(x = "fraction of the type's cells removed from the treated arm", y = "recall (true cells called down, FDR 10 %)", colour = NULL, title = "Power on the planted depletion") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 3))
ggsave(file.path(out.dir, "semi_recall.png"), p1, width = 11, height = 10, dpi = 110)
p2 <- ggplot(sm2, aes(frac, fdp.lenient, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) + geom_hline(yintercept = 0.1, linetype = 2, colour = "grey40") +
  labs(x = "fraction removed", y = "false-discovery proportion (called cells of other types, not adjacent to the planted type), FDR 10 %", colour = NULL, title = "False discoveries") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 3))
ggsave(file.path(out.dir, "semi_fdp.png"), p2, width = 11, height = 10, dpi = 110)
gres$regime <- factor(gres$regime, levels = c("unpaired", "paired"))
gnull <- gres %>% filter(size.class == "null") %>% group_by(dataset, regime, test) %>% summarise(groups.sig.perm = mean(other.sig.perm), groups.sig.cauchy = mean(other.sig.cauchy), groups.sig.cauchy.param = mean(other.sig.cauchy.param), n.groups = max(n.groups), n = n(), .groups = "drop")
write.csv(gnull, file.path(out.dir, "semi_null_groups.csv"), row.names = FALSE); print(as.data.frame(gnull), digits = 2)
gres <- gres %>% filter(size.class != "null")
gs <- gres %>% group_by(dataset, regime, test, size.class, frac) %>% summarise(planted.called = mean(planted.sig.perm), median.p.perm = median(planted.p.perm), other.perm = mean(other.sig.perm), other.cauchy = mean(other.sig.cauchy), other.cauchy.param = mean(other.sig.cauchy.param), n = n(), .groups = "drop")
write.csv(gs, file.path(out.dir, "semi_groups_summary.csv"), row.names = FALSE)
print(as.data.frame(nulls %>% filter(alpha == 0.1)), digits = 2)
print(as.data.frame(sm2 %>% select(dataset, panel, frac, method, recall, fdp.lenient) %>% arrange(dataset, panel, frac, method)), digits = 2)
print(as.data.frame(gs), digits = 2)
