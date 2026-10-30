# Diffusion neighborhoods of the density map and the scale scan (heatLogDensity, scanDensityScales, estimateDiffCellDensity(scale))

toyGraphB <- function(cells, cg, per.cell = 8, seed = 5, cross = 1) {   # per.cell edges within the type plus `cross` edges to any cell
  set.seed(seed); idx <- split(cells, cg)
  el <- do.call(rbind, lapply(cells, function(c) cbind(c, c(sample(setdiff(idx[[as.character(cg[c])]], c), per.cell), if (cross > 0) sample(setdiff(cells, c), cross)))))
  g <- igraph::simplify(igraph::graph_from_edgelist(el, directed = FALSE)); igraph::E(g)$weight <- runif(igraph::ecount(g)); g
}

test_that("heatLogDensity: per-sample densities on the log scale with a relative pseudo-count, as the benchmark computes them", {
  cells <- sprintf("c%03d", 1:240); spc <- factor(rep(sprintf("s%d", 1:6), each = 40)); names(spc) <- cells
  cg <- factor(rep(c("a", "b"), 120)); names(cg) <- cells; g <- toyGraphB(cells, cg, per.cell = 6)
  y <- cacoa:::heatLogDensity(g, spc, beta = 20, alpha = 0.3)
  expect_equal(dim(y), c(240L, 6L)); expect_setequal(rownames(y), cells); expect_equal(colnames(y), levels(spc))
  d <- 2^y - 0.3 / 240
  expect_equal(unname(colSums(d)), rep(1, 6), tolerance = 0.02)                                   # densities per cell of the sample (the cell counted as one cell of its neighbourhood: sums within the probe noise of 1)
  expect_true(all(d >= -1e-12))
  # the benchmark's computation (sccore on the connected graph, ripples clipped), the cell counted as one cell of its neighbourhood, per cell of the sample
  cells.g <- igraph::V(g)$name; sig <- sapply(levels(spc), function(s) as.numeric(spc[cells.g] == s)); rownames(sig) <- cells.g
  sm <- as.matrix(sccore::smoothSignalOnGraph(sig, filter = function(...) sccore::heatFilter(..., beta = 20), graph = g, m = 50, n.cores = 1, progress = FALSE))
  sm[sm < 0] <- 0; ks <- cacoa:::kernelEffectiveSize(g, 20, n.vec = 16)   # the map's probes
  sm <- cacoa:::ownCellDensity(sm, as.integer(spc[cells.g]), ks$self[cells.g], ks$off.size[cells.g]); sm <- t(t(sm) / colSums(sig))
  expect_equal(y[rownames(sm), ], log2(sm + 0.3 / 240), tolerance = 1e-6)
  # on a graph where every cell is well connected the change is within a few percent of the plain kernel density
  sm0 <- as.matrix(sccore::smoothSignalOnGraph(sig, filter = function(...) sccore::heatFilter(..., beta = 20), graph = g, m = 50, n.cores = 1, progress = FALSE)); sm0[sm0 < 0] <- 0; sm0 <- t(t(sm0) / colSums(sm0))
  expect_lt(median(abs(y[rownames(sm0), ] - log2(sm0 + 0.3 / 240))), 0.05)
})

test_that("ownCellDensity: the cell as one cell of its neighbourhood", {
  d <- matrix(c(0.5, 0.3, 0.2, 0.9, 0.05, 0.05), 2, byrow = TRUE); si <- c(1L, 1L)
  h <- cacoa:::ownCellDensity(d, si, self = c(0.02, 0.88), off.size = c(49, 20))
  expect_equal(rowSums(h), c(1, 1))
  expect_equal(h[1, ], (49 * c(0.48, 0.3, 0.2) / 0.98 + c(1, 0, 0)) / 50)       # well connected: about the plain row
  expect_equal(h[2, ], (20 * c(0.02, 0.05, 0.05) / 0.12 + c(1, 0, 0)) / 21)     # weakly connected: the arriving heat's composition, plus the cell itself
  expect_lt(h[2, 1], 0.5)
})

test_that("heatLogDensity works on a disconnected graph and keeps each sample's density within its component", {
  cells <- sprintf("c%03d", 1:200); comp <- factor(rep(c("A", "B"), each = 100)); names(comp) <- cells
  g <- toyGraphB(cells, comp, per.cell = 6, cross = 0); expect_gt(igraph::components(g)$no, 1)
  spc <- factor(c(rep(c("s1", "s2"), 50), rep(c("s3", "s4"), 50))); names(spc) <- cells   # s1, s2 only in A; s3, s4 only in B
  y <- cacoa:::heatLogDensity(g, spc, beta = 10, alpha = 0.3)
  floor <- log2(0.3 / 200)
  expect_true(all(abs(y[cells[comp == "B"], c("s1", "s2")] - floor) < 1e-6))                     # no mass leaks to the other component
  expect_true(all(y[cells[comp == "A"], c("s1", "s2")] > floor + 1))
  spc2 <- spc; spc2[1:3] <- NA                                                                     # cells without a sample are left out
  expect_equal(nrow(cacoa:::heatLogDensity(g, spc2, beta = 10)), 197L)
})

test_that("heatSmooth: one step up to step.max as before, a strong kernel in equal steps", {
  cells <- sprintf("c%03d", 1:200); cg <- factor(rep(c("a", "b"), 100)); names(cg) <- cells; g <- toyGraphB(cells, cg, per.cell = 6)
  set.seed(3); x <- matrix(runif(400), 200, 2, dimnames = list(cells, NULL))
  ref <- as.matrix(sccore::smoothSignalOnGraph(x, filter = function(...) sccore::heatFilter(..., beta = 20), graph = g, m = 50, n.cores = 1, progress = FALSE))
  expect_equal(as.matrix(cacoa:::heatSmooth(x, 20, graph = g)), ref, tolerance = 1e-10)
  two <- as.matrix(cacoa:::heatSmooth(x, 240, graph = g))                                          # two steps of 120
  fine <- as.matrix(sccore::smoothSignalOnGraph(x, filter = function(...) sccore::heatFilter(..., beta = 240), graph = g, m = 300, n.cores = 1, progress = FALSE))
  expect_lt(max(abs(two - fine)) / max(abs(fine)), 1e-3)
  expect_equal(unname(colSums(two)), unname(colSums(x)), tolerance = 1e-6)                          # the kernel keeps the mass
})

test_that("scanDensityScales: per-scale maps and the scan on a planted depletion; a single scale is its own map", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 40, n.genes = 40)
  cg <- cao$cell.groups; spc <- cao$sample.per.cell; cells <- names(cg)
  set.seed(8); drop <- sample(cells[cg == "ct1" & grepl("^B", as.character(spc))], 60)                # 60 of the 80 ct1 cells of group B removed
  kept <- setdiff(cells, drop); g <- toyGraphB(kept, cg[kept], per.cell = 6)
  model <- cao$.__enclos_env__$private$resolveModel(formula = ~ group, test = "group", verbose = FALSE, what = "x")
  sc <- cacoa:::scanDensityScales(g, spc[kept], model, scale = c(0.1, 0.3), n.permutations = 49, seed = 2, perm.method = "freedman-lane", groups = cg[kept])
  expect_equal(sc$cells, igraph::V(g)$name); N <- length(kept)
  expect_equal(dim(sc$P), c(N, 2L)); expect_equal(dim(sc$E), c(N, 2L)); expect_equal(dim(sc$Q), c(N, 2L)); expect_equal(sc$B, 49L)
  expect_true(all(sc$q >= 0 & sc$q <= 1)); expect_true(all(sc$p > 0 & sc$p <= 1)); expect_true(all(sc$best %in% 1:2))
  m <- apply(sc$P, 1, min); expect_true(all(sc$p >= m - 1e-12)); expect_gt(cor(rank(sc$p), rank(m)), 0.999)   # the scan p: the smallest p over scales against the permuted minima (larger, same order)
  expect_true(all(sc$scales$size[2] > sc$scales$size[1])); expect_equal(sc$scales$scale, c(0.1, 0.3)); expect_length(sc$group.test, 2)
  ct1 <- sc$cells %in% kept[cg[kept] == "ct1"]
  expect_lt(median(sc$E[ct1, 1]), -0.5); expect_gt(median(sc$E[!ct1, 1]), -0.2)                        # the depleted type has a negative log2 fold change
  expect_gt(mean(sc$q[ct1] < 0.2 & sc$E[ct1, 1] < 0), 0.5); expect_lt(mean(sc$q[!ct1] < 0.2 & sc$E[!ct1, 1] < 0), 0.1)   # (the other type rises in relative density: per-sample closure)
  one <- cacoa:::scanDensityScales(g, spc[kept], model, scale = 0.1, beta = sc$scales$beta[1], n.permutations = 49, seed = 2, perm.method = "freedman-lane")
  expect_equal(unname(one$Q[, 1]), unname(sc$Q[, 1]), tolerance = 1e-12); expect_equal(unname(one$P[, 1]), unname(sc$P[, 1]), tolerance = 1e-12)   # same relabelings: the per-scale map is the single-scale map
  expect_equal(unname(one$q), unname(one$Q[, 1]), tolerance = 1e-12); expect_true(all(diff(sc$Q[order(sc$P[, 1]), 1]) >= -1e-12))   # q non-decreasing in p
  expect_error(cacoa:::scanDensityScales(g, spc[kept], model, scale = 1.5), "shares of the cells")
  k <- cacoa:::chooseScanScale(sc); expect_true(k %in% 1:2)
  expect_equal(cacoa:::chooseScanScale(list(Q = cbind(c(0.5, 0.5, 0.5), c(0.01, 0.02, 0.5)), p = c(0.001, 0.002, 0.3), best = c(2L, 2L, 1L))), 2L)   # the scale calling the most cells
  expect_equal(cacoa:::chooseScanScale(list(Q = cbind(c(0.05, 0.05, 0.05, 0.5), c(0.01, 0.01, 0.01, 0.01)), p = rep(0.001, 4), best = rep(2L, 4))), 2L)   # 3 of 4 vs 4 of 4: only the coarse scale is saturated
  expect_equal(cacoa:::chooseScanScale(list(Q = cbind(c(0.05, 0.05, 0.05, 0.5), c(0.01, 0.01, 0.01, 0.01)), p = rep(0.001, 4), best = rep(2L, 4)), saturated = 0.7), 1L)   # both saturated: the finest
  expect_equal(cacoa:::chooseScanScale(list(Q = cbind(rep(1, 300), rep(1, 300)), p = c(0.001, 0.002, 0.003, runif(297, 0.1, 1)), best = c(1L, 1L, 2L, rep(2L, 297)))), 1L)   # nothing called: the strongest 1 %
})

test_that("class level: scale = a vector runs the scan, one value a single scale, 'auto' chooses one; NULL keeps the kNN units", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 30, n.genes = 40)
  cg <- cao$cell.groups; cells <- names(cg); g <- toyGraphB(cells, cg, per.cell = 6)
  cao$data.object <- structure(list(graph = g), class = "Conos")                                   # a Conos-like object holding the graph
  cao$estimateCellDensity(method = "graph", verbose = FALSE)
  k <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = NULL, verbose = FALSE)$diff$permutation
  expect_null(k$scales); expect_null(k$broad); expect_null(k$knn); expect_setequal(names(k$p.adj), cells)
  s <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = c(0.1, 0.3), verbose = FALSE)$diff$permutation
  expect_equal(s$scale.mode, "scan"); expect_equal(s$scale, c(0.1, 0.3)); expect_equal(s$unit, "diffusion"); expect_equal(s$moderate, "trend"); expect_equal(s$effect.scale, "finest")
  cl <- names(s$p.adj); expect_setequal(cl, cells); expect_equal(s$p.adj, s$scales$scan.q[cl]); expect_equal(s$p.pooled, s$scales$scan.p[cl])
  expect_equal(unname(s$effect), unname(s$scales$effect[cl, 1]))                                   # the finest scale's log2 fold change
  expect_true(all(s$scales$best %in% c(0.1, 0.3))); expect_equal(dim(s$scales$p), c(length(cl), 2L)); expect_equal(s$z.adj, s$adj); expect_equal(s$z.pooled, s$pooled)
  expect_lt(abs(log(s$scales$median.size[1] / 24)), log(1.2)); expect_equal(names(s$unit.size), c("0.1", "0.3"))
  b <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = c(0.1, 0.3), effect.scale = "best", verbose = FALSE)$diff$permutation
  expect_equal(unname(b$effect), unname(b$scales$effect[cbind(seq_along(cl), match(b$scales$best[cl], c(0.1, 0.3)))]))
  one <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = 0.1, verbose = FALSE)$diff$permutation
  expect_equal(one$scale.mode, "single"); expect_equal(one$scale, 0.1); expect_null(one$effect.scale)
  expect_equal(unname(one$p.adj[cl]), unname(s$scales$q[cl, 1]), tolerance = 1e-12)                  # the single-scale map is the scan's per-scale map (same relabelings)
  expect_equal(unname(one$effect[cl]), unname(s$scales$effect[cl, 1]), tolerance = 1e-8)
  a <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = "auto", verbose = FALSE)$diff$permutation
  d0 <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, verbose = FALSE)$diff$permutation   # the default is "auto" over 1, 3 and 10 %
  expect_equal(d0$scale.mode, "auto"); expect_equal(d0$scales$scale, a$scales$scale); expect_equal(d0$p.adj, a$p.adj)
  a2 <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = "auto", auto.scales = c(0.3, 0.1), verbose = FALSE)$diff$permutation
  expect_equal(a2$scale.mode, "auto"); expect_equal(a2$scales$scale, c(0.1, 0.3)); expect_true(a2$scale %in% c(0.1, 0.3)); expect_equal(a2$scales$scan.q, s$p.adj)
  expect_error(cao$estimateDiffCellDensity(type = "permutation", n.permutations = 9, scale = "auto", auto.scales = 0.1, verbose = FALSE), "at least two")
  expect_error(cao$estimateDiffCellDensity(type = "permutation", n.permutations = 9, scale = "auto", auto.scales = c(0.1, 2), verbose = FALSE), "auto.scales must")
  expect_equal(a$scale.mode, "auto"); expect_true(a$scale %in% a$scales$scale); j <- match(a$scale, a$scales$scale)   # (on 240 cells the default scales are raised to the floor)
  expect_equal(unname(a$p.adj[cl]), unname(a$scales$q[cl, j])); expect_equal(unname(a$effect[cl]), unname(a$scales$effect[cl, j]))
  expect_true(all(a$scales$scale * 240 >= 10 - 1e-9)); expect_message(cao$estimateDiffCellDensity(type = "permutation", n.permutations = 9, scale = 0.01, verbose = TRUE), "below the floor")
  expect_error(cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = c(0.1, 0.3), adjust.method = "bh", verbose = FALSE), "perm-fdr")
  expect_error(cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = 2, verbose = FALSE), "shares of the cells")
  expect_message(cao$estimateDiffCellDensity(type = "permutation", n.permutations = 29, scale = 0.1, verbose = TRUE), "diffusion neighborhoods at 10 %")
})

test_that("planted density simulation: the scan finds the changed types at the finest and the broad scale alike", {
  skip_if(isTRUE(as.logical(Sys.getenv("CACOA_SKIP_SLOW", "true"))), "slow")
  f <- "../../../test/sim_density.rds"; skip_if_not(file.exists(f), "simulated density dataset not available")
  d <- readRDS(f)
  cao <- Cacoa$new(d$con, sample.metadata = d$sample.meta, cell.groups = d$cell.groups, sample.per.cell = d$sample.per.cell, n.cores = 4, verbose = FALSE)
  cao$setOptions(seed = 1); cao$setModel(~ treatment, test = "treatment", verbose = FALSE)
  cao$estimateCellDensity(method = "graph", verbose = FALSE)
  s <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 199, scale = c(0.01, 0.03, 0.1), verbose = FALSE)$diff$permutation
  expect_equal(s$scale.mode, "scan"); expect_equal(s$scale, c(0.01, 0.03, 0.1))
  ct <- cao$cell.groups; fr <- tapply(s$p.adj[names(ct)] < 0.1, ct, mean, na.rm = TRUE)
  for (ct.i in c("IN-SST", "L2/3", "L5/6-CC", "IN-PV", "OPC")) expect_gt(fr[[ct.i]], 0.5)
  expect_lt(max(fr[c("AST-PP", "Microglia", "Neu-NRGN")]), 0.1)
  a <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 199, scale = "auto", verbose = FALSE)$diff$permutation
  expect_equal(a$scale.mode, "auto"); expect_true(a$scale %in% c(0.01, 0.03, 0.1))
})
