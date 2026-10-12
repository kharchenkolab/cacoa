## Where are the far-away "down" calls of the kNN-unit rows on PF? By cell type and processing site, their consistency
## across relabelings, and two candidate fixes: library sizes that exclude the planted type (oracle normalization) and
## the processing site as a model covariate. Usage: R_LIBS=$S/rlib:$S/rlib_milo Rscript diag_pf_knn.R [B] [n.cores]
args <- commandArgs(trailingOnly = TRUE); B <- if (length(args)) as.integer(args[1]) else 1999; n.cores <- if (length(args) > 1) as.integer(args[2]) else 8
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R")); suppressPackageStartupMessages(library(edgeR))
site <- readRDS(file.path(BENCH, "inputs_pf.rds"))$sample.meta$Processing_site; names(site) <- rownames(readRDS(file.path(BENCH, "inputs_pf.rds"))$sample.meta)
mname <- "kNN nhoods + voom-weighted perm (libsize): fdr.sd"; far.sets <- list()
for (id in c("large_f100_r1", "large_f100_r2", "large_f100_r3", "medium_f100_r1")) {
  d <- file.path(BENCH, "planted/pf", id); tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds")); g <- sub$graph; cells <- V(g)$name
  A <- as_adjacency_matrix(g, sparse = TRUE); A@x[] <- 1; v <- as.numeric(cells %in% tr$true.cells); h2 <- as.numeric(A %*% (as.numeric(A %*% v) > 0)) > 0
  dc <- readRDS(file.path(d, "decomp.rds"))$res; names(dc) <- sub("^Milo nhoods", "kNN nhoods", names(dc)); r <- dc[[mname]]
  fd <- cells[!h2 & v == 0 & r$q[cells] < 0.1 & r$effect[cells] < 0]; far.sets[[id]] <- fd
  cat(sprintf("\n== %s: %d far-down calls of %d calls; by cell type (share of the type's cells):\n", id, length(fd), sum(r$q[cells] < 0.1, na.rm = TRUE)))
  tb <- table(sub$cell.groups[fd]); sh <- tb / table(sub$cell.groups)[names(tb)]; print(round(sort(sh[tb > 20], decreasing = TRUE), 2))
  cat("   by processing site of the cell's sample (share of the site's cells):\n"); st <- site[as.character(sub$sample.per.cell[fd])]; print(round(table(st) / table(site[as.character(sub$sample.per.cell[cells])]), 3))
  cat("   arm B samples by site:", paste(sprintf("%s=%s", tr$arm.b, site[tr$arm.b]), collapse = " "), "\n")
}
j <- function(a, b) length(intersect(a, b)) / length(union(a, b)); cat("\nJaccard of far-down sets r1/r2/r3:", round(c(j(far.sets[[1]], far.sets[[2]]), j(far.sets[[1]], far.sets[[3]]), j(far.sets[[2]], far.sets[[3]])), 2), "\n")
## candidate fixes on large_f100_r1
id <- "large_f100_r1"; d <- file.path(BENCH, "planted/pf", id); tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds")); g <- sub$graph; cells <- V(g)$name; spc <- sub$sample.per.cell[cells]
A <- as_adjacency_matrix(g, sparse = TRUE); A@x[] <- 1; v <- as.numeric(cells %in% tr$true.cells); h2 <- as.numeric(A %*% (as.numeric(A %*% v) > 0)) > 0; far <- !h2 & v == 0
mi <- readRDS(file.path(d, "milo.rds")); nh <- mi$nhood.matrix[cells, , drop = FALSE]
ind <- Matrix::sparse.model.matrix(~ 0 + spc); colnames(ind) <- levels(spc); rownames(ind) <- cells
nh.counts <- as.matrix(Matrix::t(nh) %*% ind); rownames(nh.counts) <- paste0("nh", seq_len(nrow(nh.counts)))
meta <- sub$sample.meta; meta$condition <- relevel(factor(as.character(meta$condition)), "A"); meta$site <- factor(site[rownames(meta)])
meta2 <- data.frame(row.names = rownames(meta), cond = meta$condition, stratum = meta$stratum, site = meta$site)
model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum"); P <- modelPermutations(model, scheme = "block", n.permutations = B, seed = 1)$P
mapCells <- function(nh, q, eff) { tr <- summary(nh); tr <- tr[tr$x > 0, ]; o <- order(tr$i, q[tr$j]); tr <- tr[o, ]; first <- !duplicated(tr$i); qc <- rep(NA_real_, nrow(nh)); ec <- qc; qc[tr$i[first]] <- q[tr$j[first]]; ec[tr$i[first]] <- eff[tr$j[first]]; names(qc) <- names(ec) <- rownames(nh); list(q = qc, effect = ec) }
score <- function(r, nm) { q <- r$q[cells]; e <- r$effect[cells]; c1 <- is.finite(q) & q < 0.1; cat(sprintf("   %-62s recall %.2f  far-down %5d (%.3f of calls)  up-calls %5d\n", nm, mean(c1[v == 1] & e[v == 1] < 0), sum(c1 & far & e < 0), sum(c1 & far & e < 0) / max(1, sum(c1)), sum(c1 & e > 0))) }
cat("\n== large_f100_r1, kNN unit, variants:\n")
sizes <- colSums(ind); vw <- voomWeights(nh.counts, meta, norm = "libsize"); pw <- permWLS(vw$E, vw$w, vw$x, P); qw <- permWLSq(pw)
score(mapCells(nh, qw$fdr.sd, pw$effect), "weighted perm, fdr.sd, library size (as run)"); score(mapCells(nh, qw$perm.fdr, pw$effect), "weighted perm, raw pooling, library size (as run)")
## (a) oracle normalization: library sizes computed without the planted type (removes the closure shift of arm B)
keep <- !(cells %in% c(tr$true.cells)); sizes.o <- colSums(ind[keep, ]); y <- DGEList(counts = nh.counts, lib.size = as.numeric(sizes.o[colnames(nh.counts)])); design <- model.matrix(~ condition, meta[colnames(nh.counts), ]); vo <- limma::voom(y, design)
po <- permWLS(vo$E, vo$weights, as.numeric(design[, 2]), P); qo <- permWLSq(po); score(mapCells(nh, qo$fdr.sd, po$effect), "weighted perm, fdr.sd, sizes WITHOUT the planted type"); score(mapCells(nh, qo$perm.fdr, po$effect), "weighted perm, raw, sizes WITHOUT the planted type")
## (b) processing site as a covariate (cacoa's engine, unweighted, on the voom log-CPM), against the same without site
for (f in c("~ cond", "~ site + cond")) { m <- buildCacoaModel(meta2, formula = as.formula(f), test = "cond", block.vars = "stratum")
  r <- cacoa:::diffCellDensityPermutations(vw$E, m, perm.method = "freedman-lane", statistic = "t", n.permutations = B, n.cores = n.cores, seed = 1, moderate = FALSE)
  score(mapCells(nh, permutationFDR(r$stat.obs, r$stats.perm)$q, r$effect), paste("unweighted cacoa perm on voom log-CPM,", f)) }
## (c) no inheritance: index cell only
idx <- mi$nhoods$index.cell; qi <- setNames(rep(NA_real_, length(cells)), cells); ei <- qi; qi[idx] <- qw$perm.fdr; ei[idx] <- pw$effect; score(list(q = qi, effect = ei), "weighted perm, raw, index cells only (no inheritance)")
message("done")
