## Oracle normalization on PF: per-sample sizes computed WITHOUT the planted type (removes the sample-specific closure
## shift of a large depletion) for the kNN unit (weighted permutation) and for the heat kernel (weighted permutation),
## on a medium (7 %) and a large (20 %) removal. Usage: R_LIBS=$S/rlib:$S/rlib_milo Rscript diag_pf_oracle.R [B] [n.cores]
args <- commandArgs(trailingOnly = TRUE); B <- if (length(args)) as.integer(args[1]) else 1999; n.cores <- if (length(args) > 1) as.integer(args[2]) else 8
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R")); suppressPackageStartupMessages(library(edgeR))
mapCells <- function(nh, q, eff) { tr <- summary(nh); tr <- tr[tr$x > 0, ]; o <- order(tr$i, q[tr$j]); tr <- tr[o, ]; first <- !duplicated(tr$i); qc <- rep(NA_real_, nrow(nh)); ec <- qc; qc[tr$i[first]] <- q[tr$j[first]]; ec[tr$i[first]] <- eff[tr$j[first]]; names(qc) <- names(ec) <- rownames(nh); list(q = qc, effect = ec) }
for (id in c("medium_f100_r1", "large_f100_r1", "large_f100_r2")) {
  d <- file.path(BENCH, "planted/pf", id); tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds")); g <- sub$graph; cells <- V(g)$name; spc <- sub$sample.per.cell[cells]
  A <- as_adjacency_matrix(g, sparse = TRUE); A@x[] <- 1; v <- as.numeric(cells %in% tr$true.cells); h2 <- as.numeric(A %*% (as.numeric(A %*% v) > 0)) > 0; far <- !h2 & v == 0
  score <- function(r, nm) { q <- r$q[cells]; e <- r$effect[cells]; c1 <- is.finite(q) & q < 0.1; cat(sprintf("   %-58s recall %.2f  far-down %5d (%.3f of calls)  up-calls %5d\n", nm, mean(c1[v == 1] & e[v == 1] < 0), sum(c1 & far & e < 0), sum(c1 & far & e < 0) / max(1, sum(c1)), sum(c1 & e > 0))) }
  meta <- sub$sample.meta; meta$condition <- relevel(factor(as.character(meta$condition)), "A"); cond <- setNames(meta$condition, rownames(meta))
  meta2 <- data.frame(row.names = rownames(meta), cond = meta$condition, stratum = meta$stratum)
  model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum"); P <- modelPermutations(model, scheme = "block", n.permutations = B, seed = 1)$P
  ind <- Matrix::sparse.model.matrix(~ 0 + spc); colnames(ind) <- levels(spc); rownames(ind) <- cells
  sizes <- colSums(ind); sizes.o <- colSums(ind[!(cells %in% tr$true.cells), ])
  cat(sprintf("\n== %s (%s, closure %.3f): arm-B sizes lost %s %%\n", id, tr$type, tr$closure, paste(round(100 * (1 - sizes.o[tr$arm.b] / sizes[tr$arm.b])), collapse = "/")))
  mi <- readRDS(file.path(d, "milo.rds")); nh <- mi$nhood.matrix[cells, , drop = FALSE]; nh.counts <- as.matrix(Matrix::t(nh) %*% ind); rownames(nh.counts) <- paste0("nh", seq_len(nrow(nh.counts)))
  for (sz in list(full = colSums(nh.counts), oracle = colSums(nh.counts) * as.numeric(sizes.o / sizes)[match(colnames(nh.counts), names(sizes))])) {
    y <- DGEList(counts = nh.counts, lib.size = as.numeric(sz)); design <- model.matrix(~ condition, meta[colnames(nh.counts), ]); vo <- limma::voom(y, design)
    po <- permWLS(vo$E, vo$weights, as.numeric(design[, 2]), P); qo <- permWLSq(po)
    score(mapCells(nh, qo$perm.fdr, po$effect), paste0("kNN unit, weighted perm, raw pooling, sizes ", if (identical(sz, colSums(nh.counts))) "full" else "WITHOUT planted type"))
  }
  dm <- densityMats(g, spc, n.cores = n.cores)
  for (nm in c("full", "oracle")) {
    sz <- if (nm == "full") dm$sizes else sizes.o[names(dm$sizes)]
    v1 <- softVoom(dm$heat, sz, cond); pw <- permWLS(v1$E, v1$w, v1$x, P); qw <- permWLSq(pw)
    score(list(q = qw$perm.fdr, effect = pw$effect), paste0("heat kernel, weighted perm, raw pooling, sizes ", nm))
    h <- inheritQ(qw$perm.fdr, pw$effect, A, cells); score(h, paste0("heat kernel, weighted perm + inheritance, sizes ", nm))
    r <- cacoa:::diffCellDensityPermutations(t(t(dm$heat) * as.numeric(sz[colnames(dm$heat)]) / as.numeric(dm$sizes[colnames(dm$heat)])), model, perm.method = "freedman-lane", statistic = "t", n.permutations = B, n.cores = n.cores, seed = 1, moderate = TRUE)
    score(list(q = permutationFDR(r$stat.obs, r$stats.perm)$q, effect = r$effect), paste0("heat kernel, cacoa default (pooled FDR), sizes ", nm))
  }
}
message("done")
