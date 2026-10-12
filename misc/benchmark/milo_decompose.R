## Decomposition of the Milo comparison on the planted-relabeling sets: unit (Milo's refined neighbourhoods vs cacoa's
## 1-hop neighbourhoods) x test (edgeR quasi-likelihood F vs cacoa's permutation LM) x adjustment (BH / spatial FDR /
## pooled permutation FDR). Usage: R_LIBS=$S/rlib_milo:$S/rlib Rscript milo_decompose.R <scc|pf> [B] [n.cores]
## Needs <set>/milo.rds (nhood.matrix) and inputs.rds; writes <set>/decomp.rds (per-cell q-values per combination).
args <- commandArgs(trailingOnly = TRUE); ds <- args[1]; B <- if (length(args) > 1) as.integer(args[2]) else 1999; n.cores <- if (length(args) > 2) as.integer(args[3]) else 8
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
suppressPackageStartupMessages(library(edgeR))
root <- file.path(BENCH, "planted", ds)
edgerQL <- function(counts, meta) {   # counts: units x samples; Milo's test: TMM + QL F on ~ condition
  y <- DGEList(counts = counts, lib.size = colSums(counts)); y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ condition, data = meta[colnames(counts), , drop = FALSE])
  y <- estimateDisp(y, design); fit <- glmQLFit(y, design, robust = TRUE); res <- glmQLFTest(fit, coef = ncol(design))
  data.frame(logFC = res$table$logFC, p = res$table$PValue, row.names = rownames(counts))
}
voomBH <- function(counts, meta) {    # limma-voom (precision weights from the mean-variance trend) + moderated t, BH
  design <- model.matrix(~ condition, data = meta[colnames(counts), , drop = FALSE])
  y <- DGEList(counts = counts, lib.size = colSums(counts)); y <- calcNormFactors(y, method = "TMM")
  v <- limma::voom(y, design); fit <- limma::eBayes(limma::lmFit(v, design))
  data.frame(logFC = fit$coefficients[, ncol(design)], p = fit$p.value[, ncol(design)], row.names = rownames(counts))
}
mapCells <- function(nh, q, eff) {   # Milo's inheritance: per cell the smallest q over the neighbourhoods containing it, with that neighbourhood's effect
  tr <- summary(nh); tr <- tr[tr$x > 0, ]; o <- order(tr$i, q[tr$j]); tr <- tr[o, ]; first <- !duplicated(tr$i)
  qc <- rep(NA_real_, nrow(nh)); ec <- qc; qc[tr$i[first]] <- q[tr$j[first]]; ec[tr$i[first]] <- eff[tr$j[first]]
  names(qc) <- names(ec) <- rownames(nh); list(q = qc, effect = ec)
}
for (d in list.dirs(root, recursive = FALSE)) {
  out <- file.path(d, "decomp.rds"); if (file.exists(out) || !file.exists(file.path(d, "milo.rds"))) next
  t0 <- Sys.time()
  sub <- readRDS(file.path(d, "inputs.rds")); mi <- readRDS(file.path(d, "milo.rds"))
  cells <- V(sub$graph)$name; spc <- sub$sample.per.cell[cells]
  meta <- sub$sample.meta; meta$condition <- relevel(factor(as.character(meta$condition)), "A")
  meta2 <- data.frame(row.names = rownames(meta), cond = meta$condition, stratum = meta$stratum)
  model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
  ind <- Matrix::sparse.model.matrix(~ 0 + spc); colnames(ind) <- levels(spc); rownames(ind) <- cells
  sizes <- colSums(ind)
  ## (a) Milo's neighbourhoods, cacoa's permutation test on the log counts
  nh <- mi$nhood.matrix[cells, , drop = FALSE]; nh.counts <- as.matrix(Matrix::t(nh) %*% ind); rownames(nh.counts) <- paste0("nh", seq_len(nrow(nh.counts)))
  r <- cacoa:::diffCellDensityPermutations(logCounts(nh.counts, sizes), model, perm.method = "freedman-lane", statistic = "t", n.permutations = B, n.cores = n.cores, seed = 1, moderate = TRUE)
  q.nh <- list(perm.fdr = permutationFDR(r$stat.obs, r$stats.perm)$q, bh.pooled = p.adjust(r$p.pooled, "BH"), bh.perm = p.adjust(r$p.perm, "BH"))
  q.nh$fdr.sd <- { mu <- colMeans(r$stats.perm); sdv <- apply(r$stats.perm, 2, sd); permutationFDR((r$stat.obs - mu) / sdv, sweep(sweep(r$stats.perm, 2, mu), 2, sdv, "/"))$q }
  res <- list()
  for (m in names(q.nh)) res[[paste0("Milo nhoods + cacoa perm: ", m)]] <- mapCells(nh, q.nh[[m]], r$effect)
  ## (b) Milo's neighbourhoods, edgeR QL, BH (Milo without the spatial weighting) and Milo itself (from milo.rds)
  e <- edgerQL(nh.counts, meta); res[["Milo nhoods + edgeR QL: bh"]] <- mapCells(nh, p.adjust(e$p, "BH"), e$logFC)
  vm <- voomBH(nh.counts, meta); res[["Milo nhoods + voom: bh"]] <- mapCells(nh, p.adjust(vm$p, "BH"), vm$logFC)
  ## cacoa-style permutation test WITH voom precision weights (prototype of observation weights in the engine)
  P <- modelPermutations(model, scheme = "block", n.permutations = B, seed = 1)$P
  vw <- voomWeights(nh.counts, meta); pw <- permWLS(vw$E, vw$w, vw$x, P); qw <- permWLSq(pw)
  for (m in names(qw)) res[[paste0("Milo nhoods + voom-weighted perm: ", m)]] <- mapCells(nh, qw[[m]], pw$effect)
  pu <- permWLS(vw$E, vw$w * 0 + 1, vw$x, P); qu <- permWLSq(pu)        # the same closed-form test without weights (control)
  res[["Milo nhoods + unweighted perm (voom E): perm.fdr"]] <- mapCells(nh, qu$perm.fdr, pu$effect)
  res[["Milo (edgeR QL + spatial FDR)"]] <- list(q = setNames(mi$cells$milo.fdr, mi$cells$cell)[cells], effect = setNames(mi$cells$milo.logFC, mi$cells$cell)[cells])
  ## (c) cacoa's 1-hop neighbourhoods, edgeR QL with BH and with cacoa's spatial FDR
  dm <- densityMats(sub$graph, spc, n.cores = n.cores); nn <- neighbourLists(sub$graph)
  cnt <- round(t(t(dm$nb) * as.numeric(dm$sizes[colnames(dm$nb)]))); e1 <- edgerQL(cnt, meta)
  z1 <- sign(e1$logFC) * qnorm(pmin(pmax(e1$p, 1e-300), 1 - 1e-16) / 2, lower.tail = FALSE); names(z1) <- cells
  res[["1-hop + edgeR QL: bh"]] <- list(q = setNames(p.adjust(e1$p, "BH"), cells), effect = setNames(e1$logFC, cells))
  v1 <- voomBH(cnt, meta); res[["1-hop + voom: bh"]] <- list(q = setNames(p.adjust(v1$p, "BH"), cells), effect = setNames(v1$logFC, cells))
  vw1 <- voomWeights(cnt, meta); pw1 <- permWLS(vw1$E, vw1$w, vw1$x, P); qw1 <- permWLSq(pw1)
  for (m in names(qw1)) res[[paste0("1-hop + voom-weighted perm: ", m)]] <- list(q = qw1[[m]], effect = pw1$effect)
  ## a larger cacoa unit: 2-hop neighbourhood counts
  w2 <- setNames(rep(1, nlevels(spc)), levels(spc))
  nb2 <- cacoa:::estimateCellDensityGraph(sub$graph, spc, w2, smoothing = "neighbourhood", n.hops = 2)
  cnt2 <- round(t(t(nb2$density.mat[cells, ]) * as.numeric(nb2$sample.sizes[colnames(nb2$density.mat)])))
  v2 <- voomBH(cnt2, meta); res[["2-hop + voom: bh"]] <- list(q = setNames(p.adjust(v2$p, "BH"), cells), effect = setNames(v2$logFC, cells))
  vw2 <- voomWeights(cnt2, meta); pw2 <- permWLS(vw2$E, vw2$w, vw2$x, P); qw2 <- permWLSq(pw2)
  for (m in names(qw2)) res[[paste0("2-hop + voom-weighted perm: ", m)]] <- list(q = qw2[[m]], effect = pw2$effect)
  nsize <- c(`1-hop` = median(rowSums(cnt)), `2-hop` = median(rowSums(cnt2)), milo = median(colSums(nh)))
  res[["1-hop + edgeR QL: spatial-fdr"]] <- list(q = adjustCellScores(z1, setNames(e1$p, cells), method = "spatial-fdr", nn.list = nn[cells])$p.adj, effect = setNames(e1$logFC, cells))
  saveRDS(list(id = basename(d), res = res, n.nhoods = ncol(nh), nsize = nsize, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), out)
  message(basename(d), ": ", ncol(nh), " neighbourhoods; done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
}
message("done")
