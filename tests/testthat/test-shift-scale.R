# Diffusion neighbourhoods of the expression-shift map and the scale scan (diffusionNeighbourhoodLists, scanShiftScales,
# clusterFreeExpressionShifts(centres), estimateClusterFreeExpressionShifts(scale))

toyGraphS <- function(cells, cg, per.cell = 8, seed = 5, cross = 1) {   # per.cell edges within the type plus `cross` edges to any cell
  set.seed(seed); idx <- split(cells, cg)
  el <- do.call(rbind, lapply(cells, function(c) cbind(c, c(sample(setdiff(idx[[as.character(cg[c])]], c), per.cell), if (cross > 0) sample(setdiff(cells, c), cross)))))
  g <- igraph::simplify(igraph::graph_from_edgelist(el, directed = FALSE)); igraph::E(g)$weight <- runif(igraph::ecount(g)); g
}

test_that("diffusionNeighbourhoodLists: every cell in `cover` neighbourhoods, centres in their own, neighbourhoods of the scale's size", {
  cells <- sprintf("c%03d", 1:300); cg <- factor(rep(c("a", "b", "c"), c(150, 100, 50))); names(cg) <- cells; g <- toyGraphS(cells, cg, per.cell = 6, cross = 0)
  g <- igraph::permute(g, match(igraph::V(g)$name, cells))                                        # vertices in the cell order
  nb <- cacoa:::diffusionNeighbourhoodLists(g, scale = 0.1, cover = 2L, seed = 3)
  expect_true(all(nb$cover >= 1)); expect_gt(mean(nb$cover >= 2), 0.8); expect_equal(names(nb$units), as.character(nb$centres)); expect_equal(nb$cells, cells); expect_equal(names(nb$cover), cells)
  for (i in seq_along(nb$centres)) { u <- nb$units[[i]]; expect_true((nb$centres[i] - 1L) %in% u); expect_equal(u, sort(unique(u))); expect_true(all(u >= 0 & u < 300))
    expect_equal(length(unique(cg[u + 1L])), 1L) }                                                  # no edges between the types here: a neighbourhood stays in its population
  expect_true(all(lengths(nb$units) >= 10)); expect_equal(unname(nb$size), unname(pmax(10, round(nb$size))))
  expect_lte(max(nb$size), 30); expect_true(all(lengths(nb$units) >= 20 & lengths(nb$units) <= 60))   # neighbourhoods of about the scale's 30 cells (joins at most double them)
  expect_gt(median(lengths(nb$units)), 20); expect_lt(median(lengths(nb$units)), 60)
  nb2 <- cacoa:::diffusionNeighbourhoodLists(g, scale = 0.1, cover = 2L, seed = 3, beta = nb$beta); expect_equal(nb2$units, nb$units)   # seeded and reproducible
  nb3 <- cacoa:::diffusionNeighbourhoodLists(g, scale = 0.1, cover = 1L, seed = 3, beta = nb$beta); expect_lte(length(nb3$centres), length(nb$centres)); expect_true(all(nb3$cover >= 1))
})

test_that("clusterFreeExpressionShifts(centres): only the given neighbourhoods are tested and scattered to their cells", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 60, n.genes = 60, shift = 1.5)
  cg <- cao$cell.groups; cells <- names(cg); spc <- cao$sample.per.cell; g <- toyGraphS(cells, cg, per.cell = 6); g <- igraph::permute(g, match(igraph::V(g)$name, cells))
  nn1 <- cacoa:::graphNeighbourLists0(g, cells); cm.t <- Matrix::t(cao$data.object)
  nb <- cacoa:::diffusionNeighbourhoodLists(g, scale = 0.15, cover = 2L, seed = 3)
  nn <- rep(list(integer(0)), length(cells)); names(nn) <- cells; nn[nb$centres] <- unname(nb$units)
  s <- clusterFreeExpressionShifts(cm.t, spc, nn, cao$model, cao$sample.meta, statistic = "gsd", centres = nb$centres, post = "averaged", nns.graph = nn1, n.permutations = 29, seed = 1,
                                   adjust = TRUE, smooth = FALSE, return.perm = TRUE, pooled.statistic = "standardized")
  expect_equal(names(s$p.adj), cells); expect_equal(dim(s$stats.perm), c(29L, length(cells))); expect_gt(mean(is.finite(s$stat)), 0.9)   # cells in a tested neighbourhood (small units can lack samples)
  # the centres' units tested one by one (no averaging, the empty units untested): the scattered statistic is the mean over the centres whose neighbourhood holds the cell
  s2 <- clusterFreeExpressionShifts(cm.t, spc, nn, cao$model, cao$sample.meta, statistic = "gsd", post = "none", nns.graph = nn1, n.permutations = 29, seed = 1,
                                    adjust = TRUE, smooth = FALSE, centres.per.cell = NULL, pooled.statistic = "standardized")
  st <- s2$stat[cells][nb$centres]; ok <- is.finite(st); expect_gt(sum(ok), 0.7 * length(nb$centres)); expect_equal(sum(is.finite(s2$stat)), sum(ok))
  M <- Matrix::sparseMatrix(i = unlist(nb$units[ok]) + 1L, j = rep(seq_len(sum(ok)), lengths(nb$units[ok])), x = 1, dims = c(length(cells), sum(ok)))
  ref <- as.numeric(M %*% st[ok]) / Matrix::rowSums(M); i <- is.finite(ref) & is.finite(s$stat[cells])
  expect_equal(unname(s$stat[cells])[i], ref[i], tolerance = 1e-10); expect_equal(is.finite(unname(s$stat[cells])), Matrix::rowSums(M) > 0)
  expect_error(clusterFreeExpressionShifts(cm.t, spc, nn, cao$model, cao$sample.meta, statistic = "gsd", centres = c(0L, 1L), post = "averaged", nns.graph = nn1, n.permutations = 9), "index nns.per.cell")
  expect_error(clusterFreeExpressionShifts(cm.t, spc, nn, cao$model, cao$sample.meta, statistic = "distance", centres = nb$centres, post = "averaged", nns.graph = nn1, n.permutations = 9, profile = "raw"), "gsd")
})

test_that("scanShiftScales: per-scale maps and the scan; a planted shift is found, a null is quiet; a single scale is its own map", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 60, n.genes = 60, shift = 1.5)
  cg <- cao$cell.groups; cells <- names(cg); spc <- cao$sample.per.cell; g <- toyGraphS(cells, cg, per.cell = 6); cm.t <- Matrix::t(cao$data.object)   # (the graph in its own vertex order: reordered inside)
  sc <- cacoa:::scanShiftScales(cm.t, spc, g, cao$model, cao$sample.meta, scale = c(0.1, 0.3), n.permutations = 49, seed = 2, lib = Matrix::colSums(cm.t) * 2)
  N <- length(cells); expect_equal(sc$cells, cells); expect_equal(dim(sc$P), c(N, 2L)); expect_equal(dim(sc$Q), c(N, 2L)); expect_equal(dim(sc$E), c(N, 2L)); expect_equal(sc$B, 49L)
  expect_true(all(sc$q >= 0 & sc$q <= 1)); expect_true(all(sc$p > 0 & sc$p <= 1)); expect_true(all(sc$best %in% 1:2)); expect_equal(sc$scales$scale, c(0.1, 0.3))
  expect_true(sc$scales$size[2] > sc$scales$size[1]); expect_true(all(sc$scales$n.centres > 0)); expect_length(sc$group.test, 2)
  m <- apply(sc$P, 1, min); expect_true(all(sc$p >= m - 1e-12)); expect_true(all(diff(sc$p[order(m)]) >= -1e-12))   # the scan p follows the smallest p over scales
  expect_gt(mean(sc$q < 0.1), 0.5)                                                                  # the shift is planted in every cell
  one <- cacoa:::scanShiftScales(cm.t, spc, g, cao$model, cao$sample.meta, scale = 0.1, beta = sc$scales$beta[1], n.permutations = 49, seed = 2, lib = Matrix::colSums(cm.t) * 2)
  expect_equal(unname(one$Q[, 1]), unname(sc$Q[, 1]), tolerance = 1e-10); expect_equal(unname(one$P[, 1]), unname(sc$P[, 1]), tolerance = 1e-10); expect_equal(unname(one$E[, 1]), unname(sc$E[, 1]), tolerance = 1e-8)
  expect_equal(unname(one$q), unname(one$Q[, 1]), tolerance = 1e-12)
  cao0 <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 60, n.genes = 60, shift = 0, seed = 7)
  sc0 <- cacoa:::scanShiftScales(Matrix::t(cao0$data.object), cao0$sample.per.cell, g, cao0$model, cao0$sample.meta, scale = c(0.1, 0.3), n.permutations = 49, seed = 2)
  expect_lt(mean(sc0$q < 0.1), 0.2)
  expect_error(cacoa:::scanShiftScales(cm.t, spc, g, cao$model, cao$sample.meta, scale = 1.2, n.permutations = 9), "shares of the cells")
  expect_error(cacoa:::scanShiftScales(cm.t[, -1], spc, g, cao$model, cao$sample.meta, scale = 0.1, n.permutations = 9), "columns of cm")
})

test_that("class level: scale = a vector runs the scan, one value a single scale, 'auto' chooses one; NULL keeps the kNN units; profile.scale", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 60, n.genes = 60, shift = 1.5)
  cg <- cao$cell.groups; cells <- names(cg); g <- toyGraphS(cells, cg, per.cell = 6); cm <- cao$data.object
  cao$data.object <- structure(list(graph = g, getJointCountMatrix = function(raw = TRUE) cm), class = "Conos")   # a Conos-like object holding the graph and the counts
  k <- cao$estimateClusterFreeExpressionShifts(scale = NULL, n.permutations = 29, verbose = FALSE)
  expect_null(k$scales); expect_equal(k$post, "averaged"); cl <- names(k$p.adj); expect_setequal(cl, cells)
  s <- cao$estimateClusterFreeExpressionShifts(scale = c(0.1, 0.3), n.permutations = 29, verbose = FALSE)
  expect_equal(s$scale.mode, "scan"); expect_equal(s$scale, c(0.1, 0.3)); expect_equal(s$unit, "diffusion"); expect_equal(s$effect.scale, "finest"); expect_equal(s$statistic, "gsd")
  cl <- names(s$p.adj); expect_setequal(cl, cells); expect_equal(s$p.adj, s$scales$scan.q[cl]); expect_equal(s$p.pooled, s$scales$scan.p[cl]); expect_equal(s$p.value, s$p.pooled)
  expect_equal(unname(s$shifts), unname(s$scales$shifts[cl, 1])); expect_true(all(s$scales$best %in% c(0.1, 0.3))); expect_equal(dim(s$scales$p), c(length(cl), 2L))
  expect_equal(s$z.adj, s$adj); expect_equal(s$z.pooled, s$pooled); expect_true(all(s$z.adj >= 0, na.rm = TRUE)); expect_true(all(s$scales$n.centres > 0)); expect_equal(names(s$scales$median.size), NULL)
  expect_gt(mean(s$p.adj < 0.1, na.rm = TRUE), 0.5)                                                # the shift is planted in every cell
  b <- cao$estimateClusterFreeExpressionShifts(scale = c(0.1, 0.3), effect.scale = "best", n.permutations = 29, verbose = FALSE)
  expect_equal(unname(b$shifts), unname(b$scales$shifts[cbind(seq_along(cl), match(b$scales$best[cl], c(0.1, 0.3)))]))
  one <- cao$estimateClusterFreeExpressionShifts(scale = 0.1, n.permutations = 29, verbose = FALSE)
  expect_equal(one$scale.mode, "single"); expect_equal(one$scale, 0.1); expect_null(one$effect.scale)
  expect_equal(unname(one$p.adj[cl]), unname(s$scales$q[cl, 1]), tolerance = 1e-10); expect_equal(unname(one$shifts[cl]), unname(s$scales$shifts[cl, 1]), tolerance = 1e-8)
  a <- cao$estimateClusterFreeExpressionShifts(scale = "auto", n.permutations = 29, verbose = FALSE)
  d0 <- cao$estimateClusterFreeExpressionShifts(n.permutations = 29, verbose = FALSE)   # the default is "auto" over 1, 3 and 10 %
  expect_equal(d0$scale.mode, "auto"); expect_equal(d0$scales$scale, a$scales$scale); expect_equal(d0$p.adj, a$p.adj)
  a2 <- cao$estimateClusterFreeExpressionShifts(scale = "auto", auto.scales = c(0.3, 0.1), n.permutations = 29, verbose = FALSE)
  expect_equal(a2$scale.mode, "auto"); expect_equal(a2$scales$scale, c(0.1, 0.3)); expect_true(a2$scale %in% c(0.1, 0.3)); expect_equal(a2$scales$scan.q, s$p.adj)
  expect_error(cao$estimateClusterFreeExpressionShifts(scale = "auto", auto.scales = 0.1, n.permutations = 9, verbose = FALSE), "at least two")
  expect_equal(a$scale.mode, "auto"); expect_true(a$scale %in% a$scales$scale); j <- match(a$scale, a$scales$scale)   # (480 cells: the default scales are raised to the floor of 100 cells)
  expect_equal(unname(a$p.adj[cl]), unname(a$scales$q[cl, j])); expect_equal(unname(a$shifts[cl]), unname(a$scales$shifts[cl, j]))
  expect_error(cao$estimateClusterFreeExpressionShifts(scale = 2, n.permutations = 9, verbose = FALSE), "shares of the cells")
  expect_error(cao$estimateClusterFreeExpressionShifts(statistic = "distance", scale = 0.1, n.permutations = 9, verbose = FALSE), "gsd")
  expect_message(cao$estimateClusterFreeExpressionShifts(scale = 0.1, n.permutations = 9, verbose = TRUE), "diffusion neighborhoods at 10 %")
  p2 <- cao$estimateClusterFreeExpressionShifts(scale = 0.1, profile.scale = 1000, n.permutations = 9, verbose = FALSE); expect_equal(p2$settings$scale, 1000)
  expect_true(inherits(cao$estimateClusterFreeSummary("cluster.free.expr.shifts", method = "cauchy"), "data.frame"))
})
