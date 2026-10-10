# Smoke tests for the sample-level design builder and its contrast types.

makeMeta <- function(n = 24, seed = 10) {
  set.seed(seed)
  data.frame(group = factor(rep(c("A", "B", "C"), length.out = n)),
             batch = factor(rep(c("b1", "b2"), each = n / 2)),
             age = round(rnorm(n, 50, 10)),
             row.names = sprintf("s%02d", seq_len(n)))
}

test_that("triple contrast: level-means coding, contrast = endpoint difference, split by direction", {
  meta <- makeMeta()
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch + age)
  expect_equal(nrow(des$F), nrow(meta))
  expect_equal(colnames(des$F), c("groupA", "groupB", "groupC", "batchb2", "age"))
  expect_equal(names(des$contrast.F), colnames(des$F)); expect_equal(unname(des$contrast.F), c(-1, 1, 0, 0, 0))
  expect_true(all(c("num", "den") %in% names(des$contrast_endpoints_F)))
  expect_equal(unname(des$contrast_endpoints_F$num - des$contrast_endpoints_F$den), unname(des$contrast.F))
  # X = F c / c'c carries the contrast, Z the rest; together they span F
  expect_equal(colnames(des$X), "contrast"); expect_equal(ncol(des$Z), 4)
  set.seed(1); y <- rnorm(nrow(meta))
  b <- qr.coef(qr(cbind(des$X, des$Z)), y)
  expect_equal(unname(b[1]), unname(coef(lm(y ~ group + batch + age, meta))["groupB"]), tolerance = 1e-10)
  expect_equal(unname(qr.coef(qr(des$F), y)), unname(coef(lm(y ~ 0 + group + batch + age, meta))), tolerance = 1e-10)
  expect_null(des$core.rows)
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
