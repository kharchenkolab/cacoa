## Why is the pooled FDR conservative against the localized truth? From stored results: the realized false-call rate at
## the nominal q thresholds under three accountings (strict: every non-planted call is false; 1-hop: the graph
## neighbours of planted cells are excluded; 2-hop), for the default heat-kernel map; and the pooled FDR against BH on
## the exact per-cell permutation p and on the pooled p. Usage: Rscript analyze_fdr_cons.R
suppressPackageStartupMessages(library(dplyr))
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
thr <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5); rows <- list()
for (ds in c("scc", "pf", "ms")) for (d in list.dirs(file.path(S, "bench/planted", ds), recursive = FALSE)) {
  tr <- readRDS(file.path(d, "truth.rds")); if (is.na(tr$type) || tr$frac < 1) next
  ca <- readRDS(file.path(d, "cacoa.rds")); g <- readRDS(file.path(d, "inputs.rds"))$graph; cells <- igraph::V(g)$name
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE); A@x[] <- 1; v <- as.numeric(cells %in% tr$true.cells); h1 <- as.numeric(A %*% v) > 0; h2 <- as.numeric(A %*% h1) > 0
  x <- ca$heat; e <- x$effect[cells]; p.perm <- x$p.perm[cells]; p.pooled <- x$p.pooled[cells]
  qs <- list(`pooled FDR (default)` = x$q$perm.fdr[match(cells, rownames(x$q))], `BH on exact permutation p` = p.adjust(p.perm, "BH"), `BH on pooled p` = p.adjust(p.pooled, "BH"))
  for (m in names(qs)) for (t in thr) { q <- qs[[m]]; c1 <- is.finite(q) & q < t & e < 0; tp <- v == 1
    rows[[length(rows) + 1]] <- data.frame(dataset = ds, id = tr$id, type = tr$type, method = m, thr = t, n.called = sum(c1), recall = mean(c1[tp]),
      fdp.strict = if (any(c1)) mean(!tp[c1]) else NA, fdp.1hop = if (any(c1)) mean(!tp[c1] & !h1[c1]) else NA, fdp.2hop = if (any(c1)) mean(!tp[c1] & !h2[c1]) else NA) }
}
r <- bind_rows(rows) %>% group_by(dataset, method, thr) %>% summarise(across(c(n.called, recall, fdp.strict, fdp.1hop, fdp.2hop), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
options(width = 160); cat("== realized false-call rate of DOWN calls at nominal q thresholds (100 % removal, mean over types and relabelings)\n"); print(as.data.frame(r), digits = 2)
write.csv(r, "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark/figures/planted_fdr_accounting.csv", row.names = FALSE)
