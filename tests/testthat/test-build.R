# Package-level sanity: compiled entry points are registered and permutation code is deterministic.

test_that("compiled entry points are registered", {
  ns <- asNamespace("cacoa")
  for (f in c("fit_and_randomize", "fl_fwl_cpp", "estimateExpressionShiftsPairsLM", "clusterFreeGeneMat",
              "mapIds", "applyMedianFilterES", "adjustedZScoresMaxStat", "estimateCorrelationDistance",
              "clusterFreeZScoreMat", "estimateClusterFreeExpressionShiftsC")) {
    expect_true(exists(f, envir = ns, inherits = FALSE), info = f)
  }
  expect_false(exists("fit_with_focusing", envir = ns, inherits = FALSE))
  expect_false(exists("pca_project", envir = ns, inherits = FALSE))
})

test_that("estimateCorrelationDistance matches 1 - cor", {
  set.seed(1); a <- rnorm(30); b <- rnorm(30)
  expect_equal(cacoa:::estimateCorrelationDistance(a, b, TRUE), 1 - cor(a, b))
  expect_equal(cacoa:::estimateCorrelationDistance(a, b, FALSE), 1 - sum(a * b) / sqrt(sum(a^2) * sum(b^2)))
})

test_that("fit_and_randomize is reproducible for a given seed and independent of thread count", {
  set.seed(2)
  n <- 24; X <- cbind(1, rep(0:1, each = n / 2), rnorm(n)); Y <- matrix(rnorm(n * 6), n, 6)
  ctr <- c(0, 1, 0)
  run <- function(seed, cores) cacoa:::fit_and_randomize(X, Y, ctr, n_randomizations = 99, return_sampled_stats = TRUE,
                                                        n_cores = cores, seed = seed)
  r1 <- run(11L, 1L); r4 <- run(11L, 4L); r2 <- run(12L, 1L)
  expect_equal(r1$sampled_stats, r4$sampled_stats)
  expect_equal(r1$p_value, r4$p_value)
  expect_false(isTRUE(all.equal(r1$sampled_stats, r2$sampled_stats)))
  # OLS coefficients agree with lm
  expect_equal(unname(r1$coef[, 1]), unname(coef(lm(Y[, 1] ~ X - 1))), tolerance = 1e-10)
})

test_that("fl_fwl_cpp works with zero permutations and a non-empty nuisance matrix (B6)", {
  set.seed(3)
  n <- 20; X <- cbind(rep(0:1, each = n / 2)); Z <- cbind(1, rnorm(n)); Y <- matrix(rnorm(n * 3), n, 3)
  res <- cacoa:::fl_fwl_cpp(X, Z, Y, contrast = 1, n_randomizations = 0, return_sampled_stats = TRUE)
  expect_length(res$stat, 3)
  expect_true(all(is.finite(res$stat)))
  expect_null(res$sampled_stats)
})

test_that("performLMPermutations honours set.seed through the seed argument", {
  set.seed(4)
  n <- 16; meta <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)), row.names = paste0("s", 1:n))
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group)
  Y <- matrix(rnorm(n * 4), n, 4)
  set.seed(5); r1 <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 49)
  set.seed(5); r2 <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 49)
  expect_equal(r1$pval, r2$pval)
  r3 <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 49, seed = 7)
  r4 <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 49, seed = 7)
  expect_equal(r3$stats.perm, r4$stats.perm)
})
