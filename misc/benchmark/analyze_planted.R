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
  g <- readRDS(file.path(d, "inputs.rds"))$graph; cells.g <- igraph::V(g)$name
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE); A@x[] <- 1
  v <- as.numeric(cells.g %in% true.cells); h1 <- as.numeric(A %*% v) > 0; h2 <- as.numeric(A %*% h1) > 0
  adjacent <- setdiff(cells.g[h1], true.cells)                 # 1-hop neighbours of planted cells (ambiguous)
  adjacent2 <- setdiff(cells.g[h2], true.cells)                # within 2 hops: the spill-over radius of a ~100-cell neighbourhood with inheritance
  score <- function(q, effect, method, cells) {
    ok <- is.finite(q); cells <- cells[ok]; q <- q[ok]; effect <- effect[ok]
    tp.set <- cells %in% true.cells; near <- cells %in% adjacent; other <- !tp.set & !near; far <- !tp.set & !(cells %in% adjacent2)
    do.call(rbind, lapply(alphas, function(a) {
      called <- q < a
      data.frame(method = method, alpha = a, n.tested = length(q), n.true = sum(tp.set), n.called = sum(called),
                 recall = if (any(tp.set)) mean(called[tp.set] & effect[tp.set] < 0) else NA_real_,
                 other.up = mean(called[other] & effect[other] > 0),                      # expected compositional response, detected
                 other.down = mean(called[other] & effect[other] < 0),                    # false (wrong direction)
                 fdp.down = if (any(called)) sum(called & other & effect < 0) / sum(called) else NA_real_,
                 fdp.down2 = if (any(called)) sum(called & far & effect < 0) / sum(called) else NA_real_,   # beyond 2 hops: false, not spill-over
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
  if (file.exists(file.path(d, "hw.rds"))) {   # weighted permutation tests on cacoa's units (weighted_units.R)
    hw <- readRDS(file.path(d, "hw.rds"))$res
    for (m in names(hw)) res[[length(res) + 1]] <- score(hw[[m]]$q, hw[[m]]$effect, paste0("weighted: ", m), names(hw[[m]]$q))
  }
  if (file.exists(file.path(d, "decomp.rds"))) {   # unit x test x adjustment decomposition (milo_decompose.R)
    dc <- readRDS(file.path(d, "decomp.rds"))$res; names(dc) <- sub("^Milo nhoods", "kNN nhoods", names(dc))   # earlier label of the same unit
    for (m in names(dc)) res[[length(res) + 1]] <- score(dc[[m]]$q, dc[[m]]$effect, paste0("decomp: ", m), names(dc[[m]]$q))
    # neighbourhood inheritance (as in Milo) applied to cacoa's per-cell units: a cell takes the smallest q among itself and its graph
    # neighbours (the effect of that neighbour); isolates what the inheritance step contributes on a per-cell unit
    adj <- A; Matrix::diag(adj) <- 1; adj <- as(adj, "TsparseMatrix")
    inherit <- function(q, eff) { q <- q[cells.g]; eff <- eff[cells.g]; qq <- q; qq[!is.finite(qq)] <- Inf
      v <- qq[adj@j + 1L]; o <- order(adj@i, v); i <- adj@i[o] + 1L; j <- adj@j[o] + 1L; first <- !duplicated(i)
      out <- rep(NA_real_, length(q)); out[i[first]] <- qq[j[first]]; out[!is.finite(out)] <- NA; e <- eff; e[i[first]] <- eff[j[first]]
      names(out) <- names(e) <- cells.g; list(q = out, effect = e) }
    if (has.ca) { h <- inherit(setNames(ca$heat$q$perm.fdr, rownames(ca$heat$q)), ca$heat$effect); res[[length(res) + 1]] <- score(h$q, h$effect, "decomp: heat + inheritance: perm.fdr", cells.g) }
    m1 <- "1-hop + voom-weighted perm (libsize): fdr.sd"; if (!is.null(dc[[m1]])) { h <- inherit(dc[[m1]]$q, dc[[m1]]$effect); res[[length(res) + 1]] <- score(h$q, h$effect, paste0("decomp: ", m1, " + inheritance"), cells.g) }
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
  summarise(recall = mean(recall), other.up = mean(other.up), other.down = mean(other.down), fdp.down = mean(fdp.down, na.rm = TRUE), fdp.down2 = mean(fdp.down2, na.rm = TRUE), fdp.lenient = mean(fdp.lenient, na.rm = TRUE), n = n(), .groups = "drop")
write.csv(sm, file.path(out.dir, "planted_cells_summary.csv"), row.names = FALSE)
sm2 <- sm %>% filter(alpha == 0.1) %>% mutate(panel = paste0(size.class, ": ", type), size.class = factor(size.class, c("small", "medium", "large"))) %>% arrange(size.class)
sm2$panel <- factor(sm2$panel, unique(sm2$panel))
main <- c("heat: perm.fdr" = "heat kernel, pooled FDR (default)", "weighted: heat + weighted perm: perm.fdr" = "heat kernel, weighted permutation, pooled FDR", "weighted: heat + weighted perm: perm.fdr + inheritance" = "heat kernel, weighted permutation + inheritance",
          "weighted: 1-hop + weighted perm: perm.fdr" = "1-hop counts, weighted permutation, pooled FDR", "weighted: 1-hop + weighted perm: perm.fdr + inheritance" = "1-hop counts, weighted permutation + inheritance",
          "heat: fdr.sd" = "heat kernel, pooled FDR, standardized", "heat: fdr.tb" = "heat kernel, pooled FDR, ter Braak scale",
          "heat: bh.pooled" = "heat kernel, BH on pooled p", "heat: spatial.pooled" = "heat kernel, spatial FDR on pooled p", "heat: nbhd.pooled" = "heat kernel, neighbourhood BH on pooled p",
          "1-hop log: perm.fdr" = "1-hop counts, pooled FDR", "1-hop log: fdr.tb" = "1-hop counts, pooled FDR, ter Braak scale",
          "decomp: kNN nhoods + edgeR QL: bh" = "refined kNN neighbourhoods (~95 cells), edgeR QL, BH", "decomp: kNN nhoods + voom-weighted perm (libsize): fdr.sd" = "refined kNN neighbourhoods (~95 cells), voom-weighted permutation, pooled FDR",
          "decomp: 1-hop + voom-weighted perm (libsize): fdr.sd" = "1-hop counts, voom-weighted permutation, pooled FDR", "Milo (graph, spatial FDR)" = "Milo (edgeR QL, spatial FDR with inheritance)")
sm2 <- sm2 %>% filter(method %in% names(main)) %>% mutate(method = factor(main[method], levels = main))
p1 <- ggplot(sm2, aes(frac, recall, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) + scale_x_continuous(breaks = c(0.5, 0.8, 1)) +
  labs(x = "fraction of the type's cells removed from arm B", y = "recall (planted cells called down, FDR 10 %)", colour = NULL, title = "Power on the planted depletion (balanced relabeling, whole samples per arm)") + theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(out.dir, "planted_recall.png"), p1, width = 12, height = 9.5, dpi = 110)
p2 <- ggplot(sm2, aes(frac, fdp.down, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) + geom_hline(yintercept = 0.1, linetype = 2, colour = "grey40") +
  labs(x = "fraction removed", y = "false-discovery proportion (other cells called down, FDR 10 %)", colour = NULL, title = "False discoveries under the compositional truth") + theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(out.dir, "planted_fdp.png"), p2, width = 12, height = 9.5, dpi = 110)
p2b <- ggplot(sm2, aes(frac, fdp.down2, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) + geom_hline(yintercept = 0.1, linetype = 2, colour = "grey40") +
  labs(x = "fraction removed", y = "false-discovery proportion beyond 2 hops of the planted cells (FDR 10 %)", colour = NULL, title = "False discoveries excluding the spill-over radius") + theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(out.dir, "planted_fdp2.png"), p2b, width = 12, height = 9.5, dpi = 110)
p3 <- ggplot(sm2, aes(frac, other.up, colour = method)) + geom_line() + geom_point() + facet_grid(dataset ~ panel) +
  labs(x = "fraction removed", y = "fraction of the other cells called up (FDR 10 %)", colour = NULL, title = "Detection of the compositional response") + theme_bw() + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 3))
ggsave(file.path(out.dir, "planted_up.png"), p3, width = 12, height = 9.5, dpi = 110)
cat("\n== recall / fdp.down / other.up at FDR 10 % (mean over replicates)\n")
print(as.data.frame(sm2 %>% select(dataset, panel, frac, method, recall, fdp.down, other.up) %>% arrange(dataset, panel, frac, method)), digits = 2)
dc <- sm %>% filter(alpha == 0.1, grepl("^decomp: |^weighted: |^Milo|1-hop log: perm.fdr|1-hop log: fdr.tb|heat: perm.fdr", method)) %>% group_by(dataset, method, frac) %>% summarise(recall = mean(recall), fdp.down = mean(fdp.down, na.rm = TRUE), fdp.down2 = mean(fdp.down2, na.rm = TRUE), other.up = mean(other.up), .groups = "drop")
write.csv(dc, file.path(out.dir, "planted_decomp_summary.csv"), row.names = FALSE)
cat("\n== Milo decomposition (unit x test x adjustment), mean over planted types and replicates, FDR 10 %\n"); print(as.data.frame(dc %>% arrange(dataset, frac, method)), digits = 2)
cat("\n== null sets: fraction called\n"); print(as.data.frame(null %>% filter(alpha == 0.1)), digits = 2)
if (length(srows)) { cat("\n== null-scale estimates by class (median over cells, mean over sets)\n")
  print(as.data.frame(sc %>% group_by(dataset, test, scale, size.class, frac, class) %>% summarise(median = mean(median), .groups = "drop") %>% pivot_wider(names_from = class, values_from = median)), digits = 3) }
