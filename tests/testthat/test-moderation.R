test_that("fitFPrior recovers the scaled-F prior and caps its degrees of freedom", {
  set.seed(3); n <- 20000; df1 <- 12; d0 <- 8; s0 <- 2
  v <- s0 * stats::rf(n, df1, d0)                       # variances with a scaled-F spread
  pr <- fitFPrior(v, df1)
  expect_equal(pr$df, d0, tolerance = 0.25)
  expect_equal(median(pr$var), s0, tolerance = 0.1)
  expect_equal(pr$n, n)
  # a trend on a covariate: the prior variance follows it
  x <- runif(n, -1, 1); v2 <- exp(x) * stats::rf(n, df1, d0)
  pr2 <- fitFPrior(v2, df1, covariate = x)
  expect_gt(cor(log(pr2$var), x), 0.98)
  expect_equal(pr2$df, d0, tolerance = 0.3)
  # homogeneous variances: no spread beyond the sampling noise gives an infinite prior df, which the cap bounds
  v3 <- s0 * stats::rchisq(n, df1) / df1
  expect_true(is.infinite(fitFPrior(v3, df1)$df))
  expect_equal(fitFPrior(v3, df1, d0.max = 4)$df, 4)
  expect_equal(median(fitFPrior(v3, df1, d0.max = 4)$var), s0, tolerance = 0.1)
  # a vector of degrees of freedom and unusable values are handled
  pr4 <- fitFPrior(c(v[1:5000], NA, 0, Inf), c(rep(df1, 5000), 5, 5, 5))
  expect_equal(pr4$n, 5000); expect_equal(length(pr4$var), 5003)
})

test_that("fitFPrior agrees with limma's fitFDist", {
  skip_if_not_installed("limma")
  set.seed(5); n <- 10000; df1 <- 10; v <- 1.5 * stats::rf(n, df1, 6); x <- runif(n)
  a <- fitFPrior(v, df1); b <- limma::fitFDist(v, df1)
  expect_equal(a$df, b$df2, tolerance = 0.1); expect_equal(median(a$var), b$scale, tolerance = 0.05)
  v2 <- exp(2 * x) * stats::rf(n, df1, 6)
  a2 <- fitFPrior(v2, df1, covariate = x); b2 <- limma::fitFDist(v2, df1, covariate = x)
  expect_equal(a2$df, b2$df2, tolerance = 0.15); expect_gt(cor(log(a2$var), log(b2$scale)), 0.99)
})

test_that("the moderated permutation statistic runs with the own prior", {
  set.seed(7); n <- 12; m <- 400
  meta <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)), row.names = paste0("s", 1:n))
  model <- cacoa:::buildCacoaModel(meta, ~ group, test = "group", n.permutations = 29)
  Y <- matrix(rnorm(n * m), n, m, dimnames = list(rownames(meta), paste0("c", 1:m))); Y[, 1:20] <- Y[, 1:20] + rep(c(0, 2), each = n / 2)
  for (mode in c("squeeze", "trend")) {
    r <- cacoa:::performLMPermutations(model, Y, n.permutations = 29, moderate.perm = mode, moderate.d0.max = 6, return.sampled.stats = FALSE, return.residuals = FALSE)
    expect_equal(length(r$pval), m); expect_true(all(is.finite(r$stat.obs)))
    expect_lt(median(r$pval[1:20]), median(r$pval[-(1:20)]))
  }
  pr <- cacoa:::moderationPrior(Y, NULL, "trend", d0.max = 6)
  expect_lte(pr$df, 6); expect_equal(length(pr$var), m)
})
