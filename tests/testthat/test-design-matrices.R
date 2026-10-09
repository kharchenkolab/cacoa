# Smoke tests for the sample-level design builder and its contrast types.

makeMeta <- function(n = 24, seed = 10) {
  set.seed(seed)
  data.frame(group = factor(rep(c("A", "B", "C"), length.out = n)),
             batch = factor(rep(c("b1", "b2"), each = n / 2)),
             age = round(rnorm(n, 50, 10)),
             row.names = sprintf("s%02d", seq_len(n)))
}

test_that("triple contrast builds X/Z split with endpoints", {
  meta <- makeMeta()
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch + age)
  expect_equal(nrow(des$F), nrow(meta))
  expect_setequal(names(des$contrast.F), colnames(des$F))
  expect_true(all(c("num", "den") %in% names(des$contrast_endpoints_F)))
  expect_equal(unname(des$contrast_endpoints_F$num - des$contrast_endpoints_F$den), unname(des$contrast.F))
  # nuisance columns carry zero contrast weight
  expect_true(all(abs(des$contrast.F[colnames(des$Z)]) < 1e-12))
  # core rows are the samples of the contrasted levels
  expect_equal(unname(des$core.rows), meta$group %in% c("A", "B"))
})

test_that("simple contrast with at=, marginal contrast, and interaction cell contrast build", {
  meta <- makeMeta()
  s <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch + age,
                                   contrast = list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b2")))
  m <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch + age,
                                   contrast = list(type = "marginal", term = "group", num = "B", den = "A", over = "batch", weights = "equal"))
  i <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch + age,
                                   contrast = list(type = "simple", term = "group:batch", num = "B:b1", den = "A:b1"))
  for (d in list(s, m, i)) {
    expect_equal(unname(d$contrast_endpoints_F$num - d$contrast_endpoints_F$den), unname(d$contrast.F))
    expect_true(sum(abs(d$contrast.F)) > 0)
  }
  # marginal with equal weights averages the two simple contrasts
  s1 <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch + age,
                                    contrast = list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b1")))
  expect_equal(unname(m$contrast.F[names(s$contrast.F)]), unname((s$contrast.F + s1$contrast.F[names(s$contrast.F)]) / 2))
})

test_that("numeric coefficient contrast builds without endpoints", {
  meta <- makeMeta()
  d <- cacoa:::buildDesignMatrices(meta, contrast = c(age = 1), formula = ~ group + age)
  expect_equal(unname(d$contrast.F["age"]), 1)
  expect_null(d$contrast_endpoints_F)
})

test_that("permutation strata come from the factor nuisance and the cells hold only the contrasted samples", {
  meta <- makeMeta()
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch)
  plan <- modelPermutations(des, scheme = "block", n.permutations = 5)$plan
  expect_true(is.factor(plan$strata)); expect_equal(nlevels(plan$strata), 2)
  expect_equal(as.character(plan$strata), as.character(meta$batch))
  # permutation groups only contain samples of the contrasted levels
  expect_setequal(unlist(modelPermutations(des, scheme = "block", n.permutations = 5)$cells), which(meta$group %in% c("A", "B")))
})
