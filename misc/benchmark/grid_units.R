## Unit x post-processing grid for the density map, at matched effective neighbourhood size.
## Units: the diffusion kernel at several strengths (beta) with its per-cell effective neighbourhood size n_eff
## (1 / sum_k K_ik^2, Hutchinson estimate; the kernel is symmetric with unit row sums, so the density is the kernel-
## weighted fraction p_ij of the cells around i that belong to sample j, divided by the sample size), and a cacoa-native
## kNN count unit (every cell its own centre: all 1-hop neighbours ranked by edge weight, then 2-hop neighbours ranked
## by shared-neighbour count, up to k cells). Tests: the closed-form permutation LM (permWLS), unweighted and weighted
## (kernel: binomial weights size_j^2 / (p(1-p)) of the smoothed fraction; counts: voom). Post-processing: none, the
## neighbour-averaged statistic recomputed under every permutation (an exact test of the 1-hop-larger unit), and the
## neighbourhood inheritance (min q over the cell and its neighbours). Raw pooled FDR throughout.
## Usage: Rscript grid_units.R <scc|pf|ms> planted|null <B> <n.cores> [n.relabel]   -> <set>/grid.rds | $S/bench/null_grid/<ds>/r<k>.rds
args <- commandArgs(trailingOnly = TRUE); ds <- args[1]; mode <- args[2]; B <- as.integer(args[3]); n.cores <- as.integer(args[4]); n.relabel <- if (length(args) > 4) as.integer(args[5]) else 15
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
suppressPackageStartupMessages(library(edgeR))
betas <- c(1, 3, 10, 30, 100); ks <- c(50, 100, 200, 400)   # default beta = 30 has a median effective size of ~1,300 cells on SCC: the sweep goes down from it

# effective neighbourhood size of the heat kernel per cell: 1 / sum_k K_ik^2 with sum_k K_ik^2 = E[(K r)_i^2] for r = +-1
# The Chebyshev-approximated kernel is not exactly row-stochastic (row sums 0.4-2 at beta = 30), so the row sum r_i is
# estimated too (smoothing the all-ones vector): n_eff_i = r_i^2 / sum_k K_ik^2, and the smoothed fraction is c_ij / r_i.
kernelNeff <- function(graph, beta, n.vec = 24, n.cores = 1) {
  cells <- V(graph)$name; set.seed(7); R <- cbind(1, matrix(sample(c(-1, 1), length(cells) * n.vec, replace = TRUE), length(cells), n.vec)); rownames(R) <- cells
  KR <- as.matrix(sccore::smoothSignalOnGraph(R, filter = function(...) sccore::heatFilter(..., beta = beta), graph = graph, m = 50, n.cores = n.cores, progress = FALSE))
  r <- KR[, 1]; list(neff = setNames(r^2 / rowMeans(KR[, -1, drop = FALSE]^2), cells), rowsum = setNames(r, cells))
}
# kNN neighbourhoods on the graph: 1-hop by edge weight, then 2-hop by shared neighbours, up to k cells (self included)
knnNeighbourhoods <- function(graph, k, chunk = 4000L) {
  cells <- V(graph)$name; n <- length(cells)
  W <- igraph::as_adjacency_matrix(graph, attr = if (!is.null(igraph::edge_attr(graph, "weight"))) "weight" else NULL, sparse = TRUE)[cells, cells]
  A <- W; A@x[] <- 1; Matrix::diag(A) <- 1; A <- methods::as(A, "CsparseMatrix")
  deg <- Matrix::rowSums(A)
  out.i <- list(); out.j <- list()
  for (s in seq(1, n, by = chunk)) {
    rows <- s:min(n, s + chunk - 1L)
    W1 <- methods::as(W[rows, , drop = FALSE], "TsparseMatrix"); i1 <- W1@i + 1L; j1 <- W1@j + 1L; x1 <- W1@x
    o <- order(i1, -x1); i1 <- i1[o]; j1 <- j1[o]; r1 <- ave(seq_along(i1), i1, FUN = seq_along)
    keep1 <- r1 <= k - 1L                                                     # 1-hop neighbours kept (self is added below)
    sel.i <- c(seq_along(rows), i1[keep1]); sel.j <- c(rows, j1[keep1])
    need <- pmax(0L, k - 1L - as.integer(table(factor(i1[keep1], levels = seq_along(rows)))))
    if (any(need > 0)) {
      A2 <- A[rows, , drop = FALSE] %*% A; A1 <- A[rows, , drop = FALSE]
      M <- A2 - A2 * A1; M <- Matrix::drop0(M); M <- methods::as(M, "TsparseMatrix")       # 2-hop candidates with shared-neighbour counts
      i2 <- M@i + 1L; j2 <- M@j + 1L; x2 <- M@x; o <- order(i2, -x2); i2 <- i2[o]; j2 <- j2[o]
      r2 <- ave(seq_along(i2), i2, FUN = seq_along); keep2 <- r2 <= need[i2]
      sel.i <- c(sel.i, i2[keep2]); sel.j <- c(sel.j, j2[keep2])
    }
    out.i[[length(out.i) + 1]] <- rows[sel.i]; out.j[[length(out.j) + 1]] <- sel.j
  }
  N <- Matrix::sparseMatrix(i = unlist(out.i), j = unlist(out.j), x = 1, dims = c(n, n), dimnames = list(cells, cells))
  N
}
# binomial precision weights of a kernel-smoothed fraction: var(d_ij) ~ p(1-p) / (n_eff_i size_j^2) -> w_ij = size_j^2 / (p(1-p)), p from the group means
kernelWeights <- function(dmat, sizes, cond, rowsum, p.floor = 0.01) {   # p floored at 1 %: with tiny effective sizes p(1-p) -> 0 and one sample would carry the whole fit
  x <- as.numeric(factor(as.character(cond[colnames(dmat)]))) - 1
  sz <- as.numeric(sizes[colnames(dmat)]); P <- t(t(dmat) * sz) / as.numeric(rowsum[rownames(dmat)])   # smoothed fraction p_ij = c_ij / r_i
  P <- pmin(pmax(P, 0), 1)
  mA <- rowMeans(P[, x == 0, drop = FALSE]); mB <- rowMeans(P[, x == 1, drop = FALSE]); fit <- outer(mA, 1 - x) + outer(mB, x)
  fit <- pmin(pmax(fit, p.floor), 1 - p.floor); w <- t(t(1 / (fit * (1 - fit))) * sz^2)
  list(E = as.matrix(dmat), w = w, x = x)
}
runGrid <- function(graph, spc, cond, model, B, seed, n.cores, keep.q = TRUE) {
  cells <- V(graph)$name; sizes <- table(spc)[levels(spc)]
  A <- as_adjacency_matrix(graph, sparse = TRUE)[cells, cells]; A@x[] <- 1; An <- A; Matrix::diag(An) <- 1; An <- Matrix::Diagonal(x = 1 / Matrix::rowSums(An)) %*% An   # row-normalized self + neighbours
  P <- modelPermutations(model, scheme = "block", n.permutations = B, seed = seed)$P
  res <- list(); neff <- list()
  post <- function(prefix, pw, q) {
    res[[paste0(prefix, " | none")]] <<- list(q = q$perm.fdr, effect = pw$effect)
    tav <- as.numeric(An %*% pw$t.obs); Tav <- as.matrix(Matrix::t(An %*% Matrix::t(pw$t.perm)))   # neighbour-averaged statistic, observed and permuted
    res[[paste0(prefix, " | averaged")]] <<- list(q = permutationFDR(setNames(tav, cells), Tav)$q, effect = setNames(as.numeric(An %*% pw$effect), cells))
    h <- inheritQ(q$perm.fdr, pw$effect, A, cells); res[[paste0(prefix, " | inherit")]] <<- h
  }
  w1 <- setNames(rep(1, nlevels(spc)), levels(spc))
  for (b in betas) {
    dm <- cacoa:::estimateCellDensityGraph(graph, spc, w1, n.cores = n.cores, beta = b, verbose = FALSE)$density.mat[cells, ]
    kn <- kernelNeff(graph, b, n.cores = n.cores); neff[[paste0("kernel b=", b)]] <- kn$neff
    for (wt in c("unweighted", "weighted")) {
      v <- kernelWeights(dm, sizes, cond, kn$rowsum); if (wt == "unweighted") v$w <- v$w * 0 + 1
      pw <- permWLS(v$E, v$w, v$x, P); q <- list(perm.fdr = permutationFDR(pw$t.obs, pw$t.perm)$q)
      post(paste0("kernel b=", b, " | ", wt), pw, q); rm(pw); gc()
    }
  }
  ind <- Matrix::sparse.model.matrix(~ 0 + spc); colnames(ind) <- levels(spc); rownames(ind) <- cells
  for (k in ks) {
    N <- knnNeighbourhoods(graph, k); cnt <- as.matrix(N %*% ind); neff[[paste0("knn k=", k)]] <- setNames(Matrix::rowSums(N), cells)
    for (wt in c("unweighted", "weighted")) {
      v <- softVoom(cnt, setNames(rep(1, ncol(cnt)), colnames(cnt)), cond); if (wt == "unweighted") v$w <- v$w * 0 + 1   # counts: voom log-CPM and weights
      pw <- permWLS(v$E, v$w, v$x, P); q <- list(perm.fdr = permutationFDR(pw$t.obs, pw$t.perm)$q)
      post(paste0("knn k=", k, " | ", wt), pw, q); rm(pw); gc()
    }
  }
  list(res = res, neff = lapply(neff, function(v) quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95))))
}
if (mode == "planted") {
  root <- file.path(BENCH, "planted", ds)
  for (d in list.dirs(root, recursive = FALSE)) {
    out <- file.path(d, "grid.rds"); lock <- file.path(d, "grid.lock"); if (file.exists(out) || file.exists(lock)) next
    file.create(lock); t0 <- Sys.time()
    tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds")); spc <- droplevels(sub$sample.per.cell[V(sub$graph)$name])
    meta2 <- data.frame(row.names = rownames(sub$sample.meta), cond = sub$sample.meta$condition, stratum = sub$sample.meta$stratum)
    model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
    g <- runGrid(sub$graph, spc, setNames(sub$sample.meta$condition, rownames(sub$sample.meta)), model, B, tr$rep, n.cores)
    saveRDS(list(id = tr$id, res = g$res, neff = g$neff, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), out); unlink(lock)
    message(tr$id, " done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  }
} else {
  out.dir <- file.path(BENCH, "null_grid", ds); dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  inp <- tm("load inputs", loadInputs(ds)); cs <- condSpec(ds); meta <- inp$sample.meta
  cond0 <- relevel(factor(as.character(meta[[cs[1]]])), cs[2]); spc <- inp$sample.per.cell
  for (k in seq_len(n.relabel)) {
    f <- file.path(out.dir, sprintf("r%03d.rds", k)); lock <- paste0(f, ".lock"); if (file.exists(f) || file.exists(lock)) next
    file.create(lock); t0 <- Sys.time(); set.seed(1000 + k)
    lab <- rep(NA_character_, nrow(meta))
    for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
    meta2 <- data.frame(row.names = rownames(meta), cond = factor(lab, levels = c("A", "B")), stratum = cond0)
    model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
    g <- runGrid(inp$graph, spc, setNames(meta2$cond, rownames(meta2)), model, B, k, n.cores)
    calls <- do.call(rbind, lapply(names(g$res), function(m) do.call(rbind, lapply(c(0.05, 0.1), function(a) data.frame(method = m, alpha = a, n.tested = sum(is.finite(g$res[[m]]$q)), n.called = sum(g$res[[m]]$q < a, na.rm = TRUE))))))
    saveRDS(list(k = k, calls = calls, neff = g$neff, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), f); unlink(lock)
    message(sprintf("[%s] relabeling %d done in %.0f s; max calls at 5%%: %d", ds, k, as.numeric(difftime(Sys.time(), t0, units = "secs")), max(calls$n.called[calls$alpha == 0.05])))
  }
}
message("done")
