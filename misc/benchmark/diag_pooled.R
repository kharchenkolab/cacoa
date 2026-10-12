source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
S <- "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad"
args <- commandArgs(trailingOnly = TRUE); id <- args[1]; B <- as.integer(args[2]); n.cores <- as.integer(args[3])
d <- file.path(S, "bench/semi/scc_paired", id)
tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds"))
cells <- V(sub$graph)$name; spc <- sub$sample.per.cell
dm <- densityMats(sub$graph, spc, n.cores = n.cores); nn <- neighbourLists(sub$graph)
meta2 <- data.frame(row.names = rownames(sub$sample.meta), cond = sub$sample.meta$condition, patient = sub$sample.meta$patient)
model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "patient")
true.cells <- tr$true.cells
adjacent <- setdiff(unique(unlist(lapply(which(cells %in% true.cells), function(i) cells[nn[[i]] + 1L]))), true.cells)
grp <- ifelse(cells %in% true.cells, "planted", ifelse(cells %in% adjacent, "adjacent", "other"))
models <- list(unpaired.model = model, donor.covariate = buildCacoaModel(meta2, formula = ~ patient + cond, test = "cond", block.vars = "patient"))
for (pm in names(models)) {
  r <- cacoa:::diffCellDensityPermutations(dm$heat, models[[pm]], perm.method = "freedman-lane", statistic = "t", n.permutations = B, n.cores = n.cores, seed = 1, moderate = TRUE)
  P <- r$stats.perm; st <- r$stat.obs
  sdp <- apply(P, 2, sd, na.rm = TRUE); mup <- colMeans(P, na.rm = TRUE)
  cat(sprintf("\n=== %s: permutation sd of t by cell class (median): %s\n", pm, paste(sprintf("%s %.2f", names(tapply(sdp, grp, median)), tapply(sdp, grp, median)), collapse = ", ")))
  cat(sprintf("    permutation mean of t (median): %s\n", paste(sprintf("%s %.2f", names(tapply(mup, grp, median)), tapply(mup, grp, median)), collapse = ", ")))
  cat(sprintf("    observed t (median): %s\n", paste(sprintf("%s %.2f", names(tapply(st, grp, median)), tapply(st, grp, median)), collapse = ", ")))
  cat(sprintf("    effect sign among 'other' cells: %.2f positive\n", mean(r$effect[grp == "other"] > 0)))
  # pooled p with and without per-cell standardization; pooled FDR with and without standardization
  p.std <- pooledPermutationP(st, P, alternative = "two.sided")
  # raw pooling: rank |t| in the pool of |permuted t| (no standardization)
  ap <- abs(as.numeric(P)); ap <- ap[is.finite(ap)]; N <- length(ap)
  p.raw <- (N - findInterval(abs(st), sort(ap)) + 1) / (N + 1)
  Z <- sweep(sweep(P, 2, mup), 2, sdp, "/"); z <- (st - mup) / sdp
  q.fdr.raw <- permutationFDR(st, P)$q; q.fdr.std <- permutationFDR(z, Z)$q
  q.bh.std <- p.adjust(p.std, "BH"); q.bh.raw <- p.adjust(p.raw, "BH")
  tab <- function(q, nm) { c1 <- q < 0.1; data.frame(method = nm, n.called = sum(c1, na.rm = TRUE), recall = mean(c1[grp == "planted"] & r$effect[grp == "planted"] < 0, na.rm = TRUE),
    other.called = sum(c1[grp == "other"], na.rm = TRUE), other.up = sum(c1[grp == "other"] & r$effect[grp == "other"] > 0, na.rm = TRUE)) }
  print(rbind(tab(q.fdr.raw, "perm-fdr, raw t"), tab(q.fdr.std, "perm-fdr, standardized"), tab(q.bh.raw, "BH on pooled p, raw t"), tab(q.bh.std, "BH on pooled p, standardized")), row.names = FALSE, digits = 3)
  cat("    p-value quantiles among 'other' cells (1%, 5%, 25%, 50%): standardized", sprintf("%.3f", quantile(p.std[grp == "other"], c(.01, .05, .25, .5), na.rm = TRUE)), " raw", sprintf("%.3f", quantile(p.raw[grp == "other"], c(.01, .05, .25, .5), na.rm = TRUE)), " exact perm", sprintf("%.3f", quantile(r$p.perm[grp == "other"], c(.01, .05, .25, .5), na.rm = TRUE)), "\n")
}
