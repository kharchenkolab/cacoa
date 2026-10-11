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
    expect_equal(d$adjust.method, m); expect_equal(length(d$adj), length(d$raw)); expect_true(all(d$p.adj >= 0 & d$p.adj <= 1, na.rm = TRUE))
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
