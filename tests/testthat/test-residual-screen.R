# Residual-versus-covariate diagnostic for the per-column tests (misc/issues.md item 1).

test_that("screenResidualCovariates finds an omitted covariate in residual matrices and ignores an unrelated one", {
  set.seed(1); n <- 16
  meta <- data.frame(group = rep(c("A", "B"), each = n / 2), batch = rep(c("b1", "b2"), n / 2), age = round(rnorm(n, 50, 6)), row.names = sprintf("s%02d", 1:n))
  # residuals of a model without batch: batch structure remains; age is unrelated
  R <- matrix(rnorm(n * 60), n, 60, dimnames = list(rownames(meta), NULL)); R[meta$batch == "b2", ] <- R[meta$batch == "b2", ] + 0.8
  sc <- screenResidualCovariates(R, meta, exclude = "group", n.permutations = 99, seed = 2)
  expect_s3_class(sc, "cacoaCovariateScreen"); expect_equal(sc$settings$space, "residuals")
  expect_setequal(sc$table$covariate, c("batch", "age")); expect_equal(unique(sc$table$mode), "marginal")
  expect_lt(sc$table$p[sc$table$covariate == "batch"], 0.05); expect_gt(sc$table$p[sc$table$covariate == "age"], 0.1)
  # a list of matrices (one per gene): one row per gene x covariate, global max-statistic across genes
  Rl <- list(g1 = R, g2 = R[, 1:30] * 0 + rnorm(n * 30), g3 = R[, 1:20])
  scl <- screenResidualCovariates(Rl, meta, exclude = "group", n.permutations = 99, seed = 2)
  expect_equal(sort(unique(scl$table$celltype)), c("g1", "g2", "g3")); expect_equal(nrow(scl$table), 6)
  expect_lt(scl$global$p.global[scl$global$covariate == "batch"], 0.05)
  expect_error(screenResidualCovariates(R, meta, exclude = c("group", "batch", "age")), "no covariate left")
  expect_error(screenResidualCovariates(unname(R), meta), "row names")
})

test_that("Cacoa$screenResiduals works on density and cluster-free DE results and needs stored residuals", {
  skip_if_not_installed("conos")
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 40, n.genes = 80, n.cell.types = 2, shift = 1)
  cao$sample.meta$batch <- factor(rep(c("b1", "b2"), 5)); cao$sample.meta$age <- round(rnorm(10, 50, 5))
  cao$embedding <- matrix(rnorm(2 * length(cao$cell.groups)), ncol = 2, dimnames = list(names(cao$cell.groups), c("x", "y")))
  cao$estimateCellDensity(method = "kde", bins = 15, verbose = FALSE)
  cao$estimateDiffCellDensity(type = "permutation", n.permutations = 19, adjust = FALSE, verbose = FALSE)
  expect_error(cao$screenResiduals("cell.density"), "return.residuals = TRUE")
  cao$estimateDiffCellDensity(type = "permutation", n.permutations = 19, adjust = FALSE, return.residuals = TRUE, verbose = FALSE)
  sc <- cao$screenResiduals("cell.density", n.permutations = 29, verbose = FALSE)
  expect_s3_class(sc, "cacoaCovariateScreen"); expect_setequal(sc$table$covariate, c("batch", "age"))   # group is in the model
  expect_true(!is.null(cao$test.results[["cell.density.residual.screen"]]))
  expect_s3_class(cao$plotCovariateScreen(name = "cell.density.residual.screen"), "ggplot")
  # cluster-free DE keeps one samples x cells residual matrix per gene (the toy object has no cell graph: stored directly)
  set.seed(2); mk <- function() { R <- matrix(rnorm(10 * 30), 10, 30, dimnames = list(rownames(cao$sample.meta), NULL)); R[cao$sample.meta$batch == "b2", ] <- R[cao$sample.meta$batch == "b2", ] + 1; R }
  cao$test.results[["cluster.free.de"]] <- list(residuals = list(g1 = mk(), g2 = mk()))
  sc2 <- cao$screenResiduals("cluster.free.de", n.permutations = 49, verbose = FALSE)
  expect_equal(sc2$settings$source, "cluster.free.de"); expect_setequal(sc2$table$celltype, c("g1", "g2"))
  expect_lt(sc2$global$p.global[sc2$global$covariate == "batch"], 0.1)
})
