# Track C: covariate screen, variance partition, plots, metadata-separation alias, shift detail.

gramOf <- function(Y) { Yc <- sweep(Y, 2, colMeans(Y)); tcrossprod(Yc) }
distOf <- function(Y, names = sprintf("s%02d", seq_len(nrow(Y)))) { D <- as.matrix(dist(Y)); dimnames(D) <- list(names, names); D }

test_that("E2 pattern: the real driver is flagged partially, its proxy only marginally", {
  set.seed(3); n <- 40; p <- 150
  meta <- data.frame(disease = factor(rep(c("ctrl", "dis"), each = n / 2)), row.names = sprintf("s%02d", 1:n))
  meta$medication <- factor(ifelse(meta$disease == "dis", rep(c("yes", "yes", "yes", "yes", "no"), 4), rep(c("no", "no", "no", "no", "yes"), 4)))  # 80% / 20% medicated
  expect_equal(round(covariateAssociation(meta$disease, meta$medication), 1), 0.6)
  meta$noise <- rnorm(n)
  L <- diag(sqrt((1:p)^-1)) * sqrt(p / sum((1:p)^-1))
  Y <- outer(meta$disease == "dis", rnorm(p) * .8) + matrix(rnorm(n * p), n) %*% L
  D.list <- list(ct1 = distOf(Y), ct2 = distOf(Y + matrix(rnorm(n * p, 0, 0.3), n) %*% L))
  sc <- screenCovariates(D.list, meta, mode = "both", dist = "l2", n.permutations = 199, seed = 1)
  expect_s3_class(sc, "cacoaCovariateScreen")
  tb <- sc$table
  dis.p <- tb[tb$covariate == "disease" & tb$mode == "partial", ]; med.p <- tb[tb$covariate == "medication" & tb$mode == "partial", ]
  med.m <- tb[tb$covariate == "medication" & tb$mode == "marginal", ]
  expect_true(all(dis.p$padj < 0.05))
  expect_true(all(med.p$p > 0.05))                        # partial: medication explains nothing beyond disease
  expect_true(all(med.m$p < 0.1))                         # marginal: medication looks associated (proxy)
  expect_true(all(med.m$R2.adj > med.p$R2.adj))           # ... but its partial contribution is smaller
  expect_true(all(abs(tb$R2.adj[tb$covariate == "noise"]) < 0.05))
  expect_true(all(tb$R2.adj < tb$R2 | is.na(tb$R2.adj)))  # chance correction lowers R2
  g <- sc$global
  expect_lt(g$p.global[g$covariate == "disease" & g$mode == "partial"], 0.05)
  expect_gt(g$p.global[g$covariate == "noise" & g$mode == "partial"], 0.05)
  expect_output(print(sc), "Suggested model")
  expect_true("disease" %in% all.vars(sc$suggestion))
  expect_false("medication" %in% all.vars(sc$suggestion))
  # plots
  gg <- plotCovariateScreen(sc); expect_s3_class(gg, "ggplot")
  expect_s3_class(plotCovariateScreen(sc, mode = "both", value = "neglog10p"), "ggplot")
  expect_s3_class(plotCovariateScreen(sc, effect = "dispersion"), "ggplot")
  expect_s3_class(plotCovariateSummary(sc), "ggplot")
  # ring marks: significant marginally but not partially (count agrees with the table)
  n.ring <- sum(med.m$padj < 0.05 & med.p$padj >= 0.05) + sum(dis.p$padj >= 0.05 & tb$padj[tb$covariate == "disease" & tb$mode == "marginal"] < 0.05)
  d <- ggplot2::layer_data(gg, 4)
  expect_equal(nrow(d), n.ring)
  # the analytic preview is gone: permutation p-values only (the effective-dimension F stays in the table as p.analytic)
  expect_error(screenCovariates(D.list, meta, mode = "partial", dist = "l2", p.values = "analytic"), "unused argument")
  expect_true(all(sc$table$p.source == "permutation"))
})

test_that("E1 pattern: a dispersion-only difference is picked up by the dispersion score, not the location score", {
  set.seed(11); n <- 40; p <- 120
  meta <- data.frame(cond = factor(rep(c("a", "b"), each = n / 2)), batch = factor(rep(c("x", "y"), n / 2)), row.names = sprintf("s%02d", 1:n))
  Y <- matrix(rnorm(n * p), n) * ifelse(meta$cond == "b", sqrt(3), 1)
  sc <- screenCovariates(list(ct = distOf(Y)), meta, mode = "partial", dist = "l2", n.permutations = 199, seed = 2)
  r <- sc$table[sc$table$covariate == "cond", ]
  expect_lt(r$p.disp, 0.01)
  expect_gt(r$p, 0.05)
  expect_gt(r$R2.disp.adj, 0.2)
  b <- sc$table[sc$table$covariate == "batch", ]
  expect_gt(b$p, 0.05); expect_gt(b$p.disp, 0.05)
})

test_that("R2.adj is centred near zero under the null and the variance partition sums to one", {
  set.seed(5); n <- 30; p <- 100
  meta <- data.frame(a = factor(sample(c("u", "v"), n, TRUE)), b = rnorm(n), c = factor(sample(c("p", "q", "r"), n, TRUE)), row.names = sprintf("s%02d", 1:n))
  vals <- replicate(20, { Y <- matrix(rnorm(n * p), n); screenCovariates(list(ct = distOf(Y)), meta, mode = "marginal", dist = "l2", n.permutations = 0)$table$R2.adj })
  expect_lt(abs(mean(vals)), 0.02)
  Y <- outer(meta$a == "v", rnorm(p) * .5) + matrix(rnorm(n * p), n)
  vp <- variancePartition(distOf(Y), meta, c("a", "b"), dist = "l2")
  expect_named(vp, c("a", "b", "shared", "residual"))
  expect_equal(sum(vp), 1, tolerance = 1e-8)
  expect_gt(vp[["a"]], vp[["b"]])
  expect_equal(sum(attr(vp, "raw")), 1, tolerance = 1e-8)
  expect_s3_class(plotVariancePartition(list(ct1 = vp, ct2 = vp)), "ggplot")
})

test_that("screen and partition through the Cacoa object; metadata separation is an alias", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  sc <- cao$screenCovariates(n.permutations = 49, verbose = FALSE)
  expect_s3_class(sc, "cacoaCovariateScreen")
  expect_setequal(sc$covariates, c("group", "batch", "age"))
  expect_equal(sc$settings$test.variable, "group")
  expect_identical(sc, cao$test.results$covariate.screen)
  g <- sc$table[sc$table$covariate == "group" & sc$table$mode == "partial", ]
  expect_true(all(g$padj < 0.1))
  expect_s3_class(cao$plotCovariateScreen(), "ggplot")
  expect_s3_class(cao$plotCovariateSummary(), "ggplot")
  expect_s3_class(cao$plotVariancePartition(c("group", "batch")), "ggplot")
  tb <- cao$plotVariancePartition(c("group", "batch"), return.table = TRUE)
  expect_equal(dim(tb), c(2, 4))
  # composition space: one row
  sc2 <- cao$screenCovariates(space = "composition", n.permutations = 49, verbose = FALSE, name = "comp")
  expect_equal(unique(sc2$table$celltype), "composition")
  # alias of the old exploratory tool
  cao$estimateExpressionShiftMagnitudes(n.permutations = 19, verbose = FALSE)
  sep <- cao$estimateMetadataSeparation(cao$sample.meta[, c("group", "batch")], n.permutations = 49, show.warning = FALSE, verbose = FALSE)
  expect_setequal(names(sep$pvalues), c("group", "batch"))
  expect_s3_class(sep$screen, "cacoaCovariateScreen")
  expect_s3_class(cao$plotMetadataSeparation(), "ggplot")
  # adjusted sample distances through adjust.for
  Da <- cao$getSampleDistanceMatrix(space = "expression.shifts", cell.type = "ct1", adjust.for = ~ batch)
  expect_equal(dim(Da), c(12, 12)); expect_true(all(Da >= -1e-10))
  expect_s3_class(cao$plotSampleDistances(space = "expression.shifts", values = "unadjusted", adjust.for = ~ batch, color.by = "group"), "ggplot")
  expect_true(inherits(cao$plotShiftDetail("ct1"), "ggplot"))
})
