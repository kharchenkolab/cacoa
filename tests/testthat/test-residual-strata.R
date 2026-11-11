# exchange of residuals (Freedman-Lane) within the levels of the model's discrete covariates

pairedMeta <- function(n.pairs = 8, seed = 1) {
  set.seed(seed); pat <- rep(paste0("p", seq_len(n.pairs)), 2)
  data.frame(row.names = paste0("s", seq_len(2 * n.pairs)), cond = factor(rep(c("a", "b"), each = n.pairs)), patient = pat,
             sex = c("f", "m")[1 + (as.integer(factor(pat)) %% 2)], x = round(rnorm(2 * n.pairs), 2), stringsAsFactors = FALSE)
}

test_that("residuals are exchanged within the levels of a fitted factor; the former behavior is an option", {
  meta <- pairedMeta(); m <- buildCacoaModel(meta, ~ cond + x + patient, test = "cond"); d <- m$tests[[1]]$design
  plan <- permutationPlan(d, meta, scheme = "auto", n.permutations = 99)
  expect_equal(plan$scheme, "freedman-lane")
  expect_true(all(plan$in.set))
  expect_equal(as.integer(table(plan$strata)), rep(2L, 8))
  expect_equal(plan$n.distinct, 2^8)
  expect_true(any(grepl("within the levels of patient", plan$notes)))
  P <- drawPermutations(plan, 50)
  expect_true(all(apply(P, 2, function(p) all(meta$patient[p] == meta$patient))))            # no residual leaves its pair
  old <- permutationPlan(d, meta, scheme = "auto", n.permutations = 99, residuals.within = FALSE)
  expect_equal(nlevels(old$strata), 1)
  op <- options(cacoa.residuals.within = FALSE); on.exit(options(op))
  expect_equal(nlevels(permutationPlan(d, meta, scheme = "auto", n.permutations = 99)$strata), 1)
  options(op)
  # the tested variable never defines the cells; a model without discrete covariates is unchanged
  m0 <- buildCacoaModel(meta, ~ cond + x, test = "cond")
  expect_equal(nlevels(permutationPlan(m0$tests[[1]]$design, meta, scheme = "auto", n.permutations = 99)$strata), 1)
  # the factor of scale of the residuals is 1 when the exchange stays within the levels of the fitted terms
  Y <- matrix(rnorm(16 * 3), 16, dimnames = list(rownames(meta), paste0("f", 1:3)))
  r <- cacoa:::contrastPermutationTest(Y, d$F, d$contrast.F[colnames(d$F)], modelPermutations(m$tests[[1]]$design, n.permutations = 99, seed = 1)$P)
  expect_lt(r$kappa, 16 / (16 - 9))                    # among all samples it would be about n / (n - q)
})

test_that("a redundant nuisance term (sex is constant within patient) does not change the test", {
  meta <- pairedMeta(); set.seed(5)
  Y <- matrix(rnorm(16 * 4), 16, dimnames = list(rownames(meta), paste0("f", 1:4))) + outer(as.integer(factor(meta$patient)), 1:4) / 3
  Y[, 1] <- Y[, 1] + 1.5 * (meta$cond == "b")
  m1 <- buildCacoaModel(meta, ~ cond + x + patient, test = "cond")
  m2 <- suppressWarnings(buildCacoaModel(meta, ~ cond + x + sex + patient, test = "cond"))
  d1 <- m1$tests[[1]]$design; d2 <- m2$tests[[1]]$design
  expect_lt(qr(d2$F)$rank, ncol(d2$F) + (ncol(d2$F) == ncol(d1$F)))        # redundant columns, or dropped by the builder
  p1 <- modelPermutations(d1, n.permutations = 199, seed = 3); p2 <- modelPermutations(d2, n.permutations = 199, seed = 3)
  expect_equal(p1$plan$n.distinct, p2$plan$n.distinct)                     # cells: patient, and sex within patient = patient
  expect_true(all(tapply(as.integer(p1$plan$strata), p2$plan$strata, function(v) length(unique(v))) == 1))   # the same cells
  p2$P <- p1$P; p2$cells <- p1$cells                                       # the same relabelings for both models
  r1 <- cacoa:::contrastPermutationTest(Y, d1$F, d1$contrast.F[colnames(d1$F)], p1$P)
  r2 <- cacoa:::contrastPermutationTest(Y, d2$F, d2$contrast.F[colnames(d2$F)], p2$P)
  expect_equal(r1$est, r2$est, tolerance = 1e-8); expect_equal(r1$t, r2$t, tolerance = 1e-8)
  expect_equal(r1$t.perm, r2$t.perm, tolerance = 1e-6)
  f1 <- cacoa:::performLMPermutations(d1, Y, perm.method = "freedman-lane", P = p1$P, perm.cells = p1$cells)
  f2 <- cacoa:::performLMPermutations(d2, Y, perm.method = "freedman-lane", P = p2$P, perm.cells = p2$cells)
  expect_equal(unname(f1$stat.obs), unname(f2$stat.obs), tolerance = 1e-8)
  expect_equal(unname(f1$pval), unname(f2$pval), tolerance = 1e-8)
  expect_lt(f1$pval[1], 0.05)
})
