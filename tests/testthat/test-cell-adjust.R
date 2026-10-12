# FDR-type adjustments of per-cell statistics (cluster-free tests)

test_that("weighted BH with equal weights is BH, and the neighbourhood method covers every tested cell", {
  set.seed(2); p <- c(runif(200), rep(1e-4, 10))
  expect_equal(cacoa:::weightedBH(p, rep(1, length(p))), p.adjust(p, "BH"), tolerance = 1e-12)
  n <- 300; nn <- lapply(seq_len(n), function(i) { k <- sample(setdiff(seq_len(n), i), 10); as.integer(k - 1L) })  # 0-based
  z <- rnorm(n); z[1:30] <- 5; z[c(5, 50)] <- NA
  a <- adjustCellScores(z, method = "bh")
  expect_true(all(is.na(a$p.adj[c(5, 50)]))); expect_true(all(a$p.adj[1:30][-5] < 0.01)); expect_equal(sign(a$z.adj[1:4]), rep(1, 4))
  b <- adjustCellScores(z, method = "spatial-fdr", nn.list = nn)
  expect_equal(length(b$p.adj), n); expect_true(all(b$p.adj[is.finite(z)] <= 1))
  d <- adjustCellScores(z, method = "neighbourhood", nn.list = nn, seed = 3)
  idx <- d$index; expect_true(length(idx) > 0 && length(idx) < n)
  covered <- unique(c(idx, unlist(lapply(idx, function(i) nn[[i]] + 1L))))
  expect_true(all(which(is.finite(z)) %in% covered)); expect_true(all(is.finite(d$p.adj[is.finite(z)])))
  g <- adjustCellScores(abs(z), p = pnorm(abs(z), lower.tail = FALSE), method = "bh", alternative = "greater")
  expect_true(all(g$z.adj[is.finite(z)] >= 0))
})

test_that("density and cluster-free shift tests accept the FDR adjustments", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  cao$embedding <- matrix(rnorm(2 * length(cao$cell.groups)), ncol = 2, dimnames = list(names(cao$cell.groups), c("x", "y")))
  cao$estimateCellDensity(method = "kde", bins = 20, verbose = FALSE)
  for (m in c("bh", "spatial-fdr", "neighbourhood")) {
    r <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 19, adjust.method = m, verbose = FALSE)
    d <- r$diff$permutation
    expect_equal(d$adjust.method, m); expect_equal(length(d$adj), length(d$pooled)); expect_true(all(d$p.adj >= 0 & d$p.adj <= 1, na.rm = TRUE))
    if (m == "neighbourhood") expect_true(length(d$index.cells) > 0)
  }
  expect_s3_class(cao$plotDiffCellDensity(), "ggplot")
})

test_that("very small p-values give finite adjusted z-scores", {
  z <- c(a = 40, b = 2, c = -3, d = NA); p <- c(1e-320, 0.02, 0.001, NA)
  a <- adjustCellScores(z, p, method = "bh")
  expect_true(all(is.finite(a$z.adj[1:3]))); expect_gt(a$z.adj[["a"]], 30); expect_lt(a$z.adj[["c"]], 0)
  g <- adjustCellScores(abs(z), p, method = "bh", alternative = "greater"); expect_true(all(g$z.adj[1:3] >= 0))
})

test_that("pooled permutation FDR: calibrated under the null, finds a planted signal, matrix and histogram inputs agree", {
  set.seed(7); B <- 400; m <- 2000
  perm <- matrix(rnorm(B * m), B, m); stat <- rnorm(m); names(stat) <- paste0("c", seq_len(m))
  q0 <- permutationFDR(stat, perm)$q
  expect_equal(length(q0), m); expect_true(all(q0 >= 0 & q0 <= 1))
  expect_lt(mean(q0 < 0.1), 0.02)                                        # nothing (or almost) under the null
  stat[1:200] <- stat[1:200] + 4                                           # planted cells
  pf <- permutationFDR(stat, perm); q1 <- pf$q
  expect_gt(mean(q1[1:200] < 0.05), 0.75); expect_lt(mean(q1[-(1:200)] < 0.05), 0.05)
  expect_true(is.data.frame(pf$curve) && all(c("threshold", "n.observed", "n.null", "fdr") %in% names(pf$curve)))
  # q is monotone in |stat|
  o <- order(abs(stat)); expect_true(all(diff(q1[o]) <= 1e-12))
  # histogram input (one-sided) reproduces the matrix computation on a grid
  grid <- seq(0, 8, by = 0.01)
  hist <- t(apply(perm, 1, function(v) tabulate(findInterval(v, grid) + 1L, nbins = length(grid) + 1L)))
  qa <- permutationFDR(stat, perm, alternative = "greater")$q
  qb <- permutationFDR(stat, list(grid = grid, hist = hist), alternative = "greater")$q
  i <- which(qa < 0.2 | qb < 0.2); expect_lt(max(abs(qa[i] - qb[i])), 0.03)        # the fixed grid resolves the tail to 0.01
  a <- adjustCellScores(stat, method = "perm-fdr", perm = perm); expect_equal(a$p.adj, q1); expect_true(!is.null(a$fdr.curve))
  expect_error(adjustCellScores(stat, method = "perm-fdr"), "perm")
})

test_that("Cauchy combination and the group summaries behave", {
  set.seed(3)
  expect_equal(cauchyCombination(rep(0.5, 10)), 0.5, tolerance = 1e-8)
  expect_lt(cauchyCombination(c(1e-6, runif(50))), 1e-4)
  p.null <- replicate(2000, cauchyCombination(runif(20)))
  expect_lt(abs(mean(p.null < 0.05) - 0.05), 0.02)
  expect_true(is.na(cauchyCombination(c(NA, NA))))
  expect_false(is.na(cauchyCombination(c(0.01, 1))))                     # a p of exactly 1 is clamped, not fatal
  n <- 300; cells <- paste0("c", seq_len(n)); groups <- factor(rep(c("g1", "g2", "g3"), each = 100)); names(groups) <- cells
  stat <- rnorm(n); stat[1:100] <- stat[1:100] + 1.5; names(stat) <- cells
  perm <- matrix(rnorm(200 * n), 200, n, dimnames = list(NULL, cells))
  p <- 2 * pnorm(-abs(stat)); names(p) <- cells
  gt <- groupPermutationTest(stat, perm, groups)
  expect_equal(gt$group, c("g1", "g2", "g3")); expect_lt(gt$p.perm[1], 0.01); expect_gt(min(gt$p.perm[2:3]), 0.05)
  s <- summarizeCellScoresByGroup(p, groups, effect = stat, p.adj = p.adjust(p, "BH"), stat = stat, stats.perm = perm)
  expect_equal(s$method[1], "permutation"); expect_true(s$significant[1]); expect_equal(s$direction[1], "up"); expect_equal(s$p, s$p.perm)
  s2 <- summarizeCellScoresByGroup(p, groups, method = "cauchy"); expect_equal(s2$p, s2$p.cauchy); expect_lt(s2$p[1], 0.01)
  s3 <- summarizeCellScoresByGroup(p, groups, min.cells = 500); expect_true(all(is.na(s3$p.cauchy)))
})

test_that("neighbourhood-count density, log transform, perm-fdr and the group summaries run end to end", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  cg <- cao$cell.groups; spc <- cao$sample.per.cell; cells <- names(cg)
  # a cell graph: each cell linked to 8 random cells of its own type
  set.seed(5); idx <- split(cells, cg)
  el <- do.call(rbind, lapply(cells, function(c) cbind(c, sample(setdiff(idx[[as.character(cg[c])]], c), 8))))
  g <- igraph::simplify(igraph::graph_from_edgelist(el, directed = FALSE)); igraph::E(g)$weight <- 1
  w <- stats::setNames(rep(1, nlevels(spc)), levels(spc))
  d <- cacoa:::estimateCellDensityGraph(g, spc, w, smoothing = "neighbourhood", n.hops = 1)
  expect_equal(d$smoothing, "neighbourhood"); expect_equal(dim(d$counts.mat), c(length(cells), nlevels(spc)))
  expect_equal(unname(rowSums(d$counts.mat[cells, ])), unname(igraph::degree(g)[cells] + 1))                # self plus neighbours
  expect_equal(unname(colSums(d$counts.mat) / d$sample.sizes[colnames(d$counts.mat)] / colSums(d$density.mat)), rep(1, nlevels(spc)), tolerance = 1e-10)
  sizes <- d$sample.sizes[colnames(d$density.mat)]
  Y <- log2(t((t(d$density.mat[cells, ]) * sizes + 0.5) / (sizes + 1)) * 1e4)
  r <- cacoa:::diffCellDensityPermutations(Y, cao$model, n.permutations = 49, n.cores = 1, seed = 1)
  expect_equal(dim(r$stats.perm), c(49, length(cells))); expect_equal(names(r$stat.obs), cells)
  a <- adjustCellScores(r$stat.obs, method = "perm-fdr", perm = r$stats.perm)
  expect_true(all(a$p.adj >= 0 & a$p.adj <= 1, na.rm = TRUE)); expect_false(is.null(a$fdr.curve))
  gt <- groupPermutationTest(r$stat.obs, r$stats.perm, cg); expect_equal(nrow(gt), nlevels(cg)); expect_true(all(gt$p.perm > 0 & gt$p.perm <= 1))
  # the shift kernel accumulates the per-group sums and the exceedance histogram
  nns <- lapply(cells, function(c) match(idx[[as.character(cg[c])]], cells) - 1L); names(nns) <- cells
  sh <- clusterFreeExpressionShifts(Matrix::t(cao$data.object), spc, nns, cao$model, cao$sample.meta, n.permutations = 29, seed = 1, adjust = TRUE, smooth = FALSE, groups = cg)
  expect_equal(sh$settings$adjust.method, "perm-fdr"); expect_equal(ncol(sh$perm.hist$hist), length(sh$perm.hist$grid) + 1L)
  expect_equal(unname(rowSums(sh$perm.hist$hist)), rep(sum(is.finite(sh$z.pooled)), 29))                     # every tested cell in every permutation
  expect_equal(nrow(sh$group.test), nlevels(cg)); expect_equal(sh$group.test$n.tested, as.integer(table(cg[is.finite(sh$z.pooled)])))
  expect_true(all(sh$group.test$p.perm > 0 & sh$group.test$p.perm <= 1))
  s <- summarizeCellScoresByGroup(sh$p.value, cg, effect = sh$shifts, p.adj = sh$p.adj, group.test = sh$group.test, alternative = "greater")
  expect_equal(s$method[1], "permutation"); expect_equal(s$p, sh$group.test$p.perm)
  # class level (kde density on the toy): the group test is stored at test time and summarized / plotted
  cao$embedding <- matrix(rnorm(2 * length(cg)), ncol = 2, dimnames = list(cells, c("x", "y")))
  cao$estimateCellDensity(method = "kde", bins = 20, verbose = FALSE)
  rr <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 19, adjust.method = "perm-fdr", verbose = FALSE)$diff$permutation
  expect_equal(rr$adjust.method, "perm-fdr"); expect_false(is.null(rr$fdr.curve))
  # the unadjusted z-scores are the signed tail z of the pooled p (two-sided), the adjusted ones of the q-values
  ok <- is.finite(rr$pooled)
  expect_equal(unname(rr$pooled[ok]), unname(sign(rr$effect[ok]) * qnorm(rr$p.pooled[ok] / 2, lower.tail = FALSE)), tolerance = 1e-9)
  expect_true(all(is.finite(rr$adj[ok])))
  expect_s3_class(cao$plotDiffCellDensity(), "ggplot")
  expect_s3_class(cao$plotDiffCellDensity(adjust = FALSE), "ggplot")
  expect_s3_class(cao$plotDiffCellDensity(what = "effect"), "ggplot")
  expect_s3_class(cao$plotDiffCellDensity(color.range = c(-2, 2)), "ggplot")
  expect_equal(cacoa:::symmetricColorRange(c(-1, 0, 5, 100), "75%"), c(-1, 1) * unname(quantile(c(1, 0, 5, 100), 0.75)))   # quantile of |x|
  expect_equal(cacoa:::symmetricColorRange(c(-1, 1), NULL, floor = 1.96), c(-1.96, 1.96))                 # the floor
  expect_equal(cacoa:::symmetricColorRange(c(-3, 7), c(0, 4)), c(-4, 4))                                 # numeric range
  expect_equal(cacoa:::symmetricColorRange(c(-3, 37), NULL, cap = 8), c(-8, 8))                               # the cap
})

test_that("pooled permutation p-values: uniform under an iid null, no per-cell floor, matrix and histogram inputs agree", {
  set.seed(11); B <- 299; m <- 400
  # cells with different null scales and skews: standardization per cell makes them poolable
  scl <- runif(m, 0.5, 3); perm <- matrix(rgamma(B * m, shape = 4), B, m) * rep(scl, each = B)
  stat <- rgamma(m, shape = 4) * scl                                   # null draws from the same distributions
  p <- pooledPermutationP(stat, perm, alternative = "greater")
  expect_equal(length(p), m); expect_true(all(p > 0 & p <= 1))
  expect_gt(suppressWarnings(ks.test(p, "punif"))$p.value, 0.01)        # calibrated
  stat2 <- stat; stat2[7] <- 60 * scl[7]                                 # an extreme cell gets a p-value far below 1 / (B + 1)
  p2 <- pooledPermutationP(stat2, perm, alternative = "greater")
  expect_lt(p2[7], 1 / (B * m)); expect_equal(p2[-7], p[-7])
  # direct computation for one cell
  z <- (stat - colMeans(perm)) / apply(perm, 2, sd); Z <- scale(perm)
  expect_equal(unname(p[3]), (sum(Z >= z[3] - 0.005) + 1) / (length(Z) + 1), tolerance = 0.02)
  # two-sided on signed statistics; NA cells stay NA
  stat.na <- stat; stat.na[c(2, 9)] <- NA
  pt <- pooledPermutationP(stat.na - 4 * scl, perm - 4 * scl, alternative = "two.sided")
  expect_true(all(is.na(pt[c(2, 9)]))); expect_true(all(is.finite(pt[-c(2, 9)])))
  # histogram input (standardized z, one-sided) agrees with the matrix input up to the bin resolution
  grid <- seq(-6, 25, by = 0.02)
  H <- t(apply(Z, 1, function(v) tabulate(findInterval(v, grid) + 1L, nbins = length(grid) + 1L)))
  ph <- pooledPermutationP(z, list(grid = grid, hist = H), alternative = "greater")
  expect_equal(unname(ph), unname(p), tolerance = 0.05)
  expect_error(pooledPermutationP(z, list(grid = grid, hist = H), alternative = "two.sided"), "one-sided")
})
