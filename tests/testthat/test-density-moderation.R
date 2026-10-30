# Moderated statistic inside the permutation test (moderate.perm): a label-free prior shrinks the residual variance of the
# observed and of every permuted fit the same way, so the test stays exact.

test_that("moderated permutation statistic: formula, exactness and the unmoderated default", {
  set.seed(11)
  n <- 16; m <- 300
  meta <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)), cov = rnorm(n), row.names = paste0("s", 1:n))
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group)
  sdj <- exp(rnorm(m, 0, 0.6)); Y <- matrix(rnorm(n * m), n, m) * rep(sdj, each = n); rownames(Y) <- rownames(meta)
  r0 <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 99, seed = 3)
  rn <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 99, seed = 3, moderate.perm = "none")
  expect_equal(r0$stat.obs, rn$stat.obs); expect_equal(r0$pval, rn$pval)
  # fixed prior: t = effect / sqrt(cXc (rss + d0 s0^2) / (df + d0)), s0^2 = median label-free (centred) variance
  rf <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 99, seed = 3, moderate.perm = "fixed", moderate.d0 = 4)
  X <- des$F; cvec <- des$contrast.F; fitb <- lm.fit(X, Y); rss <- colSums(fitb$residuals^2); df <- n - ncol(X)
  cXc <- as.numeric(t(cvec) %*% solve(crossprod(X)) %*% cvec); s0 <- median(colSums(scale(Y, scale = FALSE)^2) / (n - 1))
  eff <- as.numeric(t(cvec) %*% fitb$coefficients)
  expect_equal(unname(rf$stat.obs), eff / sqrt(cXc * (rss + 4 * s0) / (df + 4)), tolerance = 1e-8)
  expect_equal(unname(rf$df), rep(df + 4, m))
  # two groups, block relabeling: rss falls with |effect| under every relabeling, so a cell's own exact p is unchanged by
  # moderation (it acts across cells: the pooled p and the ranking of cells change, the statistics are rescaled per cell)
  expect_equal(r0$pval, rf$pval); expect_false(isTRUE(all.equal(unname(r0$stat.obs), unname(rf$stat.obs))))
  # squeeze and trend: priors from limma's fitFDist on the label-free variances
  rs <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 99, seed = 3, moderate.perm = "squeeze")
  rt <- cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 99, seed = 3, moderate.perm = "trend")
  expect_true(all(is.finite(rs$pval)) && all(is.finite(rt$pval)))
  pr <- cacoa:::moderationPrior(Y, NULL, "squeeze"); expect_length(pr$var, m); expect_gt(pr$df, 0)
  # exactness: under the null the permutation p-values are uniform (no planted effect; 300 columns)
  expect_lt(abs(mean(rs$pval < 0.05) - 0.05), 0.04); expect_lt(abs(mean(rf$pval < 0.2) - 0.2), 0.07)
  # pooled null (the density map's p): moderated statistics are pooled across cells; a planted effect stays detectable
  Y2 <- Y; Y2[meta$group == "B", 1:60] <- Y2[meta$group == "B", 1:60] + 1.5 * rep(sdj[1:60], each = n / 2)
  pp <- cacoa:::performLMPermutations(des, Y2, perm.method = "block", n.permutations = 199, seed = 3, moderate.perm = "squeeze",
                                      pooled = TRUE, z.grid = seq(0, 50, by = 0.01), use.abs = TRUE, return.sampled.stats = FALSE)
  expect_false(is.null(pp$pooled))
  expect_gt(mean(pp$pval[1:60] < 0.05), 0.5); expect_lt(mean(pp$pval[-(1:60)] < 0.05), 0.12)
  # Freedman-Lane with a nuisance covariate: prior from the residuals on the nuisance design
  desc <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + cov)
  rfl <- cacoa:::performLMPermutations(desc, Y, perm.method = "freedman-lane", n.permutations = 99, seed = 3, moderate.perm = "fixed")
  expect_true(all(is.finite(rfl$pval)))
  # robust fits do not use the closed form: moderation is refused
  expect_error(cacoa:::performLMPermutations(des, Y, perm.method = "block", n.permutations = 9, robust.method = "huber", moderate.perm = "fixed"), "closed-form")
})
