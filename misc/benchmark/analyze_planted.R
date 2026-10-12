## Figures and tables of the planted-relabeling power study (planted_relabel.R). Usage: Rscript analyze_planted.R [out.dir]
## Truth with compositional closure: planted cells should be called DOWN; other cells are expected UP by the closure
## factor; a call of another cell DOWN is a false discovery (fdp.down); the graph neighbours of planted cells are ambiguous.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
args <- commandArgs(trailingOnly = TRUE); out.dir <- if (length(args)) args[1] else "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
alphas <- c(0.05, 0.1)
rows <- list(); srows <- list(); grows <- list()
for (ds in c("scc", "pf", "ms")) for (d in list.dirs(file.path(S, "bench/planted", ds), recursive = FALSE)) {
  has.ca <- file.exists(file.path(d, "cacoa.rds")); if (!has.ca && !file.exists(file.path(d, "decomp.rds"))) next
  tr <- readRDS(file.path(d, "truth.rds")); ca <- if (has.ca) readRDS(file.path(d, "cacoa.rds")) else NULL
  true.cells <- tr$true.cells
  adjacent <- if (has.ca) ca$adjacent else { g <- readRDS(file.path(d, "inputs.rds"))$graph; unique(names(unlist(igraph::adjacent_vertices(g, intersect(true.cells, igraph::V(g)$name))))) }
  adjacent <- setdiff(adjacent, true.cells)
  score <- function(q, effect, method, cells) {
    ok <- is.finite(q); cells <- cells[ok]; q <- q[ok]; effect <- effect[ok]
    tp.set <- cells %in% true.cells; near <- cells %in% adjacent; other <- !tp.set & !near
    do.call(rbind, lapply(alphas, function(a) {
      called <- q < a
      data.frame(method = method, alpha = a, n.tested = length(q), n.true = sum(tp.set), n.called = sum(called),
                 recall = if (any(tp.set)) mean(called[tp.set] & effect[tp.set] < 0) else NA_real_,
                 other.up = mean(called[other] & effect[other] > 0),                      # expected compositional response, detected
                 other.down = mean(called[other] & effect[other] < 0),                    # false (wrong direction)
                 fdp.down = if (any(called)) sum(called & other & effect < 0) / sum(called) else NA_real_,
                 fdp.lenient = if (any(called)) mean(!tp.set[called] & !near[called]) else NA_real_)   # the old accounting (any other call is false)
    }))
  }
  res <- list()
  for (tst in (if (has.ca) c("heat", "nb") else character(0))) { x <- ca[[tst]]; cells <- rownames(x$q)
    for (m in names(x$q)) res[[length(res) + 1]] <- score(x$q[[m]], x$effect[cells], paste(if (tst == "heat") "heat" else "1-hop log", m, sep = ": "), cells)
    if (!is.null(x$scale) && length(true.cells)) {   # null-scale estimates by cell class
      cls <- ifelse(cells %in% true.cells, "planted", ifelse(cells %in% adjacent, "adjacent", "other"))
      for (sc in names(x$scale)) srows[[length(srows) + 1]] <- data.frame(dataset = ds, id = tr$id, type = tr$type, size.class = tr$size.class, frac = tr$frac, test = tst, scale = sc,
                                                                            class = c("planted", "other"), median = c(median(x$scale[[sc]][cls == "planted"], na.rm = TRUE), median(x$scale[[sc]][cls == "other"], na.rm = TRUE)))
    }
  }
  if (file.exists(file.path(d, "milo.rds"))) {
    mi <- readRDS(file.path(d, "milo.rds"))$cells
    res[[length(res) + 1]] <- score(mi$milo.fdr, mi$milo.logFC, "Milo (graph, spatial FDR)", mi$cell)
  }
  if (file.exists(file.path(d, "decomp.rds"))) {   # unit x test x adjustment decomposition (milo_decompose.R)
    dc <- readRDS(file.path(d, "decomp.rds"))$res
    for (m in names(dc)) res[[length(res) + 1]] <- score(dc[[m]]$q, dc[[m]]$effect, paste0("decomp: ", m), names(dc[[m]]$q))
  }
  res <- do.call(rbind, res); res$dataset <- ds; res$id <- tr$id; res$type <- tr$type; res$size.class <- tr$size.class; res$frac <- tr$frac; res$rep <- tr$rep; res$closure <- tr$closure
  rows[[length(rows) + 1]] <- res
  for (tst in (if (has.ca) c("heat", "nb") else character(0))) { g <- ca$groups[[tst]]; if (is.null(g)) next
    i <- if (is.na(tr$type)) NA else which(g$group == tr$type)
    grows[[length(grows) + 1]] <- data.frame(dataset = ds, id = tr$id, type = tr$type, size.class = tr$size.class, frac = tr$frac, test = tst,
      p.perm = if (is.na(i)) NA else g$p.perm[i], p.cauchy.perm = if (is.na(i)) NA else g$p.cauchy[i], p.cauchy.pooled = if (is.na(i)) NA else g$p.cauchy.pooled[i],
      n.sig.perm = sum(p.adjust(g$p.perm, "BH") < 0.05, na.rm = TRUE), n.sig.cauchy.perm = sum(p.adjust(g$p.cauchy, "BH") < 0.05, na.rm = TRUE), n.sig.cauchy.pooled = sum(p.adjust(g$p.cauchy.pooled, "BH") < 0.05, na.rm = TRUE)) }
}
res <- bind_rows(rows); write.csv(res, file.path(out.dir, "planted_cells_all.csv"), row.names = FALSE)
if (length(srows)) { sc <- bind_rows(srows); write.csv(sc, file.path(out.dir, "planted_scales.csv"), row.names = FALSE) }
if (length(grows)) { gr <- bind_rows(grows); write.csv(gr, file.path(out.dir, "planted_groups.csv"), row.names = FALSE) }
null <- res %>% filter(is.na(type)) %>% group_by(dataset, method, alpha) %>% summarise(mean.frac = mean(n.called / n.tested), max.frac = max(n.called / n.tested), n = n(), .groups = "drop")
write.csv(null, file.path(out.dir, "planted_null_calls.csv"), row.names = FALSE)
sm <- res %>% filter(!is.na(type)) %>% group_by(dataset, method, alpha, size.class, type, frac) %>%
  summarise(recall = mean(recall), other.up = mean(other.up), other.down = mean(other.down), fdp.down = mean(fdp.down, na.rm = TRUE), fdp.lenient = mean(fdp.lenient, na.rm = TRUE), n = n(), .groups = "drop")
write.csv(sm, file.path(out.dir, "planted_cells_summary.csv"), row.names = FALSE)
sm2 <- sm %>% filter(alpha == 0.1) %>% mutate(panel = paste0(size.class, ": ", type), size.class = factor(size.class, c("small", "medium", "large"))) %>% arrange(size.class)
sm2$panel <- factor(sm2$panel, unique(sm2$panel))
main <- c("heat: perm.fdr" = "heat kernel, pooled FDR (default)", "heat: fdr.sd" = "heat kernel, pooled FDR, standardized", "heat: fdr.tb" = "heat kernel, pooled FDR, ter Braak scale",
          "heat: bh.pooled" = "heat kernel, BH on pooled p", "heat: spatial.pooled" = "heat kernel, spatial FDR on pooled p", "heat: nbhd.pooled" = "heat kernel, neighbourhood BH on pooled p",
          "1-hop log: perm.fdr" = "1-hop counts, pooled FDR", "1-hop log: fdr.tb" = "1-hop counts, pooled FDR, ter Braak scale",
          "decomp: Milo nhoods + edgeR QL: bh" = "Milo neighbourhoods, edgeR QL, BH", "decomp: Milo nhoods + voom-weighted perm: fdr.sd" = "Milo neighbourhoods, voom-weighted permutation, pooled FDR",
          "decomp: 1-hop + voom-weighted perm: fdr.sd" = "1-hop counts, voom-weighted permutation, pooled FDR", "Milo (graph, spatial FDR)" = "Milo (edgeR QL, spatial FDR with inheritance)")
sm2 <- sm2 %>% filter(method %in% names(main)) %>% mutate(method = factor(main[method], levels = main))
p1 <- ggplot(sm2, aes(frac, recall, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) + scale_x_continuous(breaks = c(0.5, 0.8, 1)) +
  labs(x = "fraction of the type's cells removed from arm B", y = "recall (planted cells called down, FDR 10 %)", colour = NULL, title = "Power on the planted depletion (balanced relabeling, whole samples per arm)") + theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(out.dir, "planted_recall.png"), p1, width = 12, height = 9.5, dpi = 110)
p2 <- ggplot(sm2, aes(frac, fdp.down, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) + geom_hline(yintercept = 0.1, linetype = 2, colour = "grey40") +
  labs(x = "fraction removed", y = "false-discovery proportion (other cells called down, FDR 10 %)", colour = NULL, title = "False discoveries under the compositional truth") + theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(out.dir, "planted_fdp.png"), p2, width = 12, height = 9.5, dpi = 110)
p3 <- ggplot(sm2, aes(frac, other.up, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) +
  labs(x = "fraction removed", y = "fraction of the other cells called up (FDR 10 %)", colour = NULL, title = "Detection of the compositional response") + theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(out.dir, "planted_up.png"), p3, width = 12, height = 9.5, dpi = 110)
cat("\n== recall / fdp.down / other.up at FDR 10 % (mean over replicates)\n")
print(as.data.frame(sm2 %>% select(dataset, panel, frac, method, recall, fdp.down, other.up) %>% arrange(dataset, panel, frac, method)), digits = 2)
dc <- sm2 %>% filter(grepl("^decomp: |^Milo|1-hop log: perm.fdr|1-hop log: fdr.tb", method)) %>% group_by(dataset, method, frac) %>% summarise(recall = mean(recall), fdp.down = mean(fdp.down, na.rm = TRUE), other.up = mean(other.up), .groups = "drop")
write.csv(dc, file.path(out.dir, "planted_decomp_summary.csv"), row.names = FALSE)
cat("\n== Milo decomposition (unit x test x adjustment), mean over planted types and replicates, FDR 10 %\n"); print(as.data.frame(dc %>% arrange(dataset, frac, method)), digits = 2)
cat("\n== null sets: fraction called\n"); print(as.data.frame(null %>% filter(alpha == 0.1)), digits = 2)
if (length(srows)) { cat("\n== null-scale estimates by class (median over cells, mean over sets)\n")
  print(as.data.frame(sc %>% group_by(dataset, test, scale, size.class, frac, class) %>% summarise(median = mean(median), .groups = "drop") %>% pivot_wider(names_from = class, values_from = median)), digits = 3) }
