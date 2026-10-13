# Neighbourhood units and post-processing of the cluster-free maps (R/cell_neighbourhoods.R)

toyGraph <- function(cells, cg, per.cell = 8, seed = 5, cross = 1) {   # per.cell edges within the type plus `cross` edges to any cell (connected)
  set.seed(seed); idx <- split(cells, cg)
  el <- do.call(rbind, lapply(cells, function(c) cbind(c, c(sample(setdiff(idx[[as.character(cg[c])]], c), per.cell), sample(setdiff(cells, c), cross)))))
  g <- igraph::simplify(igraph::graph_from_edgelist(el, directed = FALSE)); igraph::E(g)$weight <- runif(igraph::ecount(g)); g
}

test_that("default neighbourhood size: 10 cells per sample, at least 100, capped at 1 % of the cells but never below 5 per sample", {
  expect_equal(defaultNeighbourhoodSize(35000, 18), 180L)
  expect_equal(defaultNeighbourhoodSize(90000, 8), 100L)
  expect_equal(defaultNeighbourhoodSize(49000, 40), 400L)
  expect_equal(defaultNeighbourhoodSize(5000, 18), 90L)        # the 1 % cap (50) is lifted to 5 per sample
  expect_equal(defaultNeighbourhoodSize(2557, 20), 100L)       # the planted simulation: 20 samples, 2.6k cells
  expect_equal(defaultNeighbourhoodSize(1000, 4), 20L)
})

test_that("kNN neighbourhoods: k is a minimum; every graph neighbour is kept, two-hop neighbours fill up by shared neighbours", {
  cao <- makeToyCacoa(n.per.group = c(A = 3, B = 3), cells.per.sample = 30, n.genes = 40)
  cg <- cao$cell.groups; cells <- names(cg); g <- toyGraph(cells, cg)
  deg <- igraph::degree(g)[cells]
  N5 <- knnNeighbourhoodMatrix(g, 5, cells = cells)
  expect_equal(dim(N5), c(length(cells), length(cells))); expect_equal(unname(Matrix::rowSums(N5)), unname(pmax(5, deg + 1)))   # no truncation of the 1-hop neighbourhood
  expect_true(all(Matrix::diag(N5) == 1))
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)[cells, cells]; A@x[] <- 1
  expect_true(all(N5[A > 0] == 1))                                                                   # every graph neighbour is in
  N40 <- knnNeighbourhoodMatrix(g, 40, cells = cells)
  expect_equal(unname(Matrix::rowSums(N40)), unname(pmax(40, deg + 1)))
  d <- igraph::distances(g, v = cells[1], to = cells, weights = NA); expect_true(all(d[1, names(which(N40[cells[1], ] > 0))] <= 2))   # nothing beyond two hops
  # the two-hop fill takes the candidates with the most shared neighbours
  two <- setdiff(names(which(N40[cells[1], ] > 0)), c(cells[1], names(which(A[cells[1], ] > 0))))
  shared <- as.numeric(A[cells[1], ] %*% A); names(shared) <- cells; cand <- setdiff(names(which(shared > 0)), c(cells[1], names(which(A[cells[1], ] > 0))))
  expect_true(min(shared[two]) >= max(shared[setdiff(cand, two)]))
})

test_that("kernel effective size: Hutchinson estimate agrees with the exact kernel; beta solved for a target size", {
  cao <- makeToyCacoa(n.per.group = c(A = 3, B = 3), cells.per.sample = 25, n.genes = 40, n.cell.types = 1)
  cells <- names(cao$cell.groups); g <- toyGraph(cells, cao$cell.groups, per.cell = 5)
  I <- diag(length(cells)); dimnames(I) <- list(cells, cells)
  K <- as.matrix(sccore::smoothSignalOnGraph(I, filter = function(...) sccore::heatFilter(..., beta = 3), graph = g, m = 50, n.cores = 1, progress = FALSE))
  exact <- rowSums(K)^2 / rowSums(K^2)
  est <- kernelEffectiveSize(g, 3, n.vec = 400)
  expect_setequal(names(est$size), cells); exact <- exact[names(est$size)]
  expect_lt(median(abs(est$size / exact - 1)), 0.1); expect_gt(cor(est$size, exact), 0.9)
  expect_equal(unname(est$row.sum), unname(rowSums(K)[names(est$size)]), tolerance = 1e-6)
  sb <- solveKernelBeta(g, target.size = 20, range = c(0.3, 100), n.vec = 100)
  expect_true(all(diff(sb$trace$median.size[order(sb$trace$beta)]) >= -1e-8))                       # wider kernel, larger size
  expect_lt(abs(log(median(sb$size) / 20)), log(1.15))
})

test_that("neighbour averaging and inheritance do what they say", {
  cao <- makeToyCacoa(n.per.group = c(A = 3, B = 3), cells.per.sample = 20, n.genes = 40)
  cells <- names(cao$cell.groups); g <- toyGraph(cells, cao$cell.groups, per.cell = 4)
  set.seed(1); st <- setNames(rnorm(length(cells)), cells); st[3] <- NA; P <- matrix(rnorm(7 * length(cells)), 7, length(cells), dimnames = list(NULL, cells))
  eff <- setNames(rnorm(length(cells)), cells)
  av <- averageOverNeighbours(st, P, g, effect = eff)
  A <- igraph::as_adjacency_matrix(g, sparse = FALSE)[cells, cells]; A[A > 0] <- 1; diag(A) <- 1; An <- A / rowSums(A)
  s0 <- st; s0[is.na(s0)] <- 0
  expect_equal(unname(av$stat), unname(as.numeric(An %*% s0)), tolerance = 1e-10); expect_equal(unname(av$stats.perm[2, ]), unname(as.numeric(An %*% P[2, ])), tolerance = 1e-10)
  expect_equal(unname(av$effect), unname(as.numeric(An %*% eff)), tolerance = 1e-10)
  q <- setNames(runif(length(cells)), cells); q[5] <- NA
  ih <- inheritNeighbourQ(q, eff, g)
  qq <- q; qq[is.na(qq)] <- Inf
  for (c in cells[1:10]) { nb <- c(c, names(which(A[c, ] > 0))); j <- nb[which.min(qq[nb])]; expect_equal(unname(ih$q[c]), unname(qq[j])); expect_equal(unname(ih$effect[c]), unname(eff[j])) }
  expect_true(all(ih$q <= qq, na.rm = TRUE))
})

test_that("graph density with the kNN unit and the post-processing options run at the class level", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 30, n.genes = 60)
  cg <- cao$cell.groups; cells <- names(cg); spc <- cao$sample.per.cell; g <- toyGraph(cells, cg, per.cell = 6)
  w <- stats::setNames(rep(1, nlevels(spc)), levels(spc))
  d <- cacoa:::estimateCellDensityGraph(g, spc, w, smoothing = "knn", k = 12, verbose = FALSE)
  deg <- igraph::degree(g)
  expect_equal(d$smoothing, "knn"); expect_equal(d$k, 12); expect_equal(unname(rowSums(d$counts.mat)), unname(pmax(12, deg[rownames(d$counts.mat)] + 1))); expect_equal(unname(d$unit.size), unname(pmax(12, deg[names(d$unit.size)] + 1)))
  expect_equal(unname(d$density.mat[cells[1], ]), unname(d$counts.mat[cells[1], ] / d$sample.sizes[colnames(d$density.mat)]))
  dh <- cacoa:::estimateCellDensityGraph(g, spc, w, smoothing = "heat", beta = 3, verbose = FALSE)
  expect_equal(dh$beta, 3); expect_equal(length(dh$unit.size), length(cells)); expect_true(all(dh$unit.size > 1))
  # class level on a Conos-like object holding the graph: default unit is kNN with the size rule, log transform, averaged post-processing
  cao$data.object <- structure(list(graph = g), class = "Conos")
  expect_message(cao$estimateCellDensity(method = "graph", verbose = TRUE), "kNN neighbourhoods of 40 cells")
  r0 <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, verbose = FALSE)$diff$permutation
  expect_equal(r0$post, "averaged"); expect_equal(r0$transform, "log"); expect_equal(r0$adjust.method, "perm-fdr")
  expect_setequal(names(r0$pooled), cells); expect_true(all(is.finite(r0$pooled))); expect_true(all(r0$p.adj >= 0 & r0$p.adj <= 1)); expect_equal(unname(r0$unit.size), unname(pmax(40, deg[names(r0$unit.size)] + 1)))
  r1 <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, post = "none", verbose = FALSE)$diff$permutation
  expect_equal(r1$post, "none"); expect_false(isTRUE(all.equal(r0$pooled, r1$pooled)))
  expect_equal(unname(r0$effect), unname(as.numeric(cacoa:::neighbourAveragingMatrix(g, names(r1$effect)) %*% r1$effect)), tolerance = 1e-10)   # averaged effect
  r2 <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, post = "inherit", verbose = FALSE)$diff$permutation
  expect_equal(r2$post, "inherit"); expect_true(all(r2$p.adj <= r2$p.adj.cell, na.rm = TRUE)); expect_equal(unname(r2$pooled), unname(r1$pooled))
  cao$embedding <- matrix(rnorm(2 * length(cells)), ncol = 2, dimnames = list(cells, c("x", "y")))
  expect_s3_class(cao$plotDiffCellDensity(), "ggplot"); expect_s3_class(cao$plotDiffCellDensity(what = "effect"), "ggplot")
  # the heat kernel sized in cells
  expect_message(cao$estimateCellDensity(method = "graph", smoothing = "heat", target.size = 15, verbose = TRUE), "median effective size")
  rh <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, verbose = FALSE)$diff$permutation
  expect_equal(rh$transform, "none"); expect_lt(abs(log(median(rh$unit.size) / 15)), log(1.2))
})

test_that("cluster-free shifts: kNN units and the post-processing options", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 30, n.genes = 60, shift = 1.5)
  cg <- cao$cell.groups; cells <- names(cg); spc <- cao$sample.per.cell; g <- toyGraph(cells, cg, per.cell = 6)
  adj <- igraph::as_adjacency_matrix(g, sparse = TRUE)[cells, cells]; adj@x[] <- 1; Matrix::diag(adj) <- 1; adj <- as(adj, "TsparseMatrix")
  nn1 <- split(adj@j, factor(adj@i + 1L, levels = seq_along(cells))); nn1 <- lapply(nn1, as.integer); names(nn1) <- cells
  N <- as(knnNeighbourhoodMatrix(g, 20, cells = cells), "TsparseMatrix"); nnk <- split(N@j, factor(N@i + 1L, levels = seq_along(cells))); nnk <- lapply(nnk, as.integer); names(nnk) <- cells
  cm.t <- Matrix::t(cao$data.object)
  s0 <- clusterFreeExpressionShifts(cm.t, spc, nn1, cao$model, cao$sample.meta, n.permutations = 29, seed = 1, adjust = TRUE, smooth = FALSE)
  sk <- clusterFreeExpressionShifts(cm.t, spc, nnk, cao$model, cao$sample.meta, n.permutations = 29, seed = 1, adjust = TRUE, smooth = FALSE, nns.graph = nn1)
  expect_equal(s0$post, "none"); expect_equal(names(sk$stat), cells); expect_false(isTRUE(all.equal(s0$stat, sk$stat)))
  sa <- clusterFreeExpressionShifts(cm.t, spc, nnk, cao$model, cao$sample.meta, n.permutations = 29, seed = 1, adjust = TRUE, smooth = FALSE, nns.graph = nn1, post = "averaged", return.perm = TRUE)
  expect_equal(sa$post, "averaged"); expect_equal(dim(sa$stats.perm), c(29, length(cells)))
  An <- cacoa:::neighbourAveragingList(nn1, cells); f0 <- sk$stat; f0[!is.finite(f0)] <- 0
  ok <- is.finite(sk$stat); expect_true(any(ok) && any(!ok))                                              # the small unit leaves some cells untested
  expect_equal(unname(sa$stat[ok]), unname(as.numeric(An %*% f0)[ok]), tolerance = 1e-8); expect_true(all(is.na(sa$stat[!ok])))   # the averaged statistic; untested cells stay NA
  expect_true(all(sa$p.adj >= 0 & sa$p.adj <= 1, na.rm = TRUE)); expect_true(all(is.na(sa$p.param)))
  expect_equal(unname(sa$z.pooled), unname(qnorm(sa$p.pooled, lower.tail = FALSE)), tolerance = 1e-8)
  si <- clusterFreeExpressionShifts(cm.t, spc, nnk, cao$model, cao$sample.meta, n.permutations = 29, seed = 1, adjust = TRUE, smooth = FALSE, nns.graph = nn1, post = "inherit")
  expect_equal(si$post, "inherit"); expect_true(all(si$p.adj <= si$p.adj.cell, na.rm = TRUE)); expect_equal(si$p.adj.cell, sk$p.adj)
  expect_message(clusterFreeExpressionShifts(cm.t, spc, nn1, cao$model, cao$sample.meta, n.permutations = 9, seed = 1, adjust = TRUE, smooth = FALSE, post = "averaged", adjust.method = "bh", verbose = TRUE), "needs adjust.method")
})
