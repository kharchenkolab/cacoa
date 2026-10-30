# Sample spaces (expression per cell type / all cell types, pooled pseudobulk, composition): distances, maps with the
# covariate R2, and the summary rows of the covariate screen.

# toy object whose composition varies between samples (cell-type shares depend on the batch)
makeSpaceToy <- function(n.per.group = c(A = 6, B = 6), cells.per.sample = 60, n.genes = 80, shift = 1.5, seed = 3, contrast = c("group", "B", "A")) {
  set.seed(seed)
  groups <- rep(names(n.per.group), n.per.group)
  samples <- sprintf("%s%d", groups, unlist(lapply(n.per.group, seq_len)))
  meta <- data.frame(group = factor(groups), batch = factor(rep(c("b1", "b2"), length.out = length(samples))),
                     age = round(rnorm(length(samples), 50, 8)), row.names = samples)
  sample.per.cell <- factor(rep(samples, each = cells.per.sample))
  share <- ifelse(meta[as.character(sample.per.cell), "batch"] == "b2", 0.75, 0.4)   # ct1 share by batch
  cell.groups <- factor(ifelse(stats::runif(length(sample.per.cell)) < share, "ct1", "ct2"))
  names(sample.per.cell) <- names(cell.groups) <- sprintf("cell%04d", seq_along(sample.per.cell))
  base <- exp(rnorm(n.genes, 1, 0.7)); type.eff <- exp(rnorm(n.genes, 0, 0.6))
  mu <- outer(base, rep(1, length(cell.groups))); mu[, cell.groups == "ct2"] <- mu[, cell.groups == "ct2"] * type.eff
  isB <- meta[as.character(sample.per.cell), "group"] == "B"; sh <- seq_len(n.genes) <= n.genes / 4
  mu[sh, isB] <- mu[sh, isB] * exp(shift)
  cm <- Matrix::Matrix(matrix(rpois(length(mu), mu), n.genes, dimnames = list(sprintf("g%03d", seq_len(n.genes)), names(sample.per.cell))), sparse = TRUE)
  suppressWarnings(suppressMessages(Cacoa$new(data.object = cm, sample.metadata = meta, sample.per.cell = sample.per.cell, cell.groups = cell.groups,
                                              contrast = contrast, n.cores = 1, verbose = FALSE)))
}

test_that("space names: new names and the earlier aliases", {
  expect_equal(normalizeSampleSpace("expression.shifts"), "expression")
  expect_equal(normalizeSampleSpace("coda"), "composition")
  expect_equal(normalizeSampleSpace("pseudo.bulk"), "pseudobulk")
  expect_equal(normalizeSampleSpace("pseudobulk"), "pseudobulk")
  expect_error(normalizeSampleSpace("umap"), "space must be one of")
})

test_that("combined expression distance: cell-weighted vs equal weights", {
  s <- c("a", "b", "c")
  D1 <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3, dimnames = list(s, s)); attr(D1, "n.cells") <- c(a = 100, b = 100, c = 100)
  D2 <- matrix(c(0, 5, 5, 5, 0, 5, 5, 5, 0), 3, dimnames = list(s, s)); attr(D2, "n.cells") <- c(a = 4, b = 4, c = 4)
  De <- prepareJointExpressionDistance(list(t1 = D1, t2 = D2), weighting = "equal")
  expect_equal(De, (D1 + D2) / 2, ignore_attr = TRUE)
  Dc <- prepareJointExpressionDistance(list(t1 = D1, t2 = D2), weighting = "cells")
  expect_equal(Dc["a", "b"], (1 * 10 + 5 * 2) / 12)                 # sqrt(min(n.cells)) weights 10 and 2
  # a type absent in one sample contributes nothing to that sample's pairs
  D3 <- D2[c("a", "b"), c("a", "b")]; attr(D3, "n.cells") <- c(a = 4, b = 4)
  De3 <- prepareJointExpressionDistance(list(t1 = D1, t3 = D3), weighting = "equal")
  expect_equal(De3["a", "c"], D1["a", "c"]); expect_equal(De3["a", "b"], (1 + 5) / 2)
})

test_that("every space is available before a model and before any shift results", {
  cao <- makeSpaceToy(contrast = NULL)
  expect_null(cao$model)
  D <- cao$getSampleDistanceMatrix(space = "expression", cell.type = "ct1")
  expect_equal(dim(D), c(12, 12)); expect_equal(attr(D, "space"), "expression"); expect_false(is.null(attr(D, "n.cells")))
  Dj <- cao$getSampleDistanceMatrix(); De <- cao$getSampleDistanceMatrix(joint.weighting = "equal")
  expect_equal(dim(Dj), c(12, 12)); expect_false(isTRUE(all.equal(Dj, De, check.attributes = FALSE)))
  Dp <- cao$getSampleDistanceMatrix(space = "pseudobulk"); Dc <- cao$getSampleDistanceMatrix(space = "coda")
  expect_equal(dim(Dp), c(12, 12)); expect_equal(attr(Dc, "space"), "composition"); expect_equal(attr(Dc, "dist.type"), "l2")
  # the pooled pseudobulk is the one-group pseudobulk of all cells
  cms <- cacoa:::extractRawCountMatrices(cao$data.object, transposed = TRUE)
  cg1 <- setNames(factor(rep("all", length(cao$sample.per.cell))), names(cao$sample.per.cell))
  pb <- pseudobulkPerCellType(cms, cg1, cao$sample.per.cell)
  ref <- sampleDistanceMatrices(pb$cm.per.type, dist = attr(Dp, "dist.type"))$all
  expect_equal(unclass(Dp)[rownames(ref), colnames(ref)], ref, ignore_attr = TRUE)
  expect_error(cao$getSampleDistanceMatrix(space = "composition", cell.type = "ct1"), "cell.type applies")
  expect_error(cao$plotSampleDistances(values = "adjusted"), "needs a model")
  for (sp in c("expression", "pseudobulk", "composition")) expect_s3_class(cao$plotSampleDistances(space = sp, color.by = "batch"), "ggplot")
})

test_that("the map subtitle reports the screen's R2 on the shown distances", {
  cao <- makeSpaceToy()
  gg <- cao$plotSampleDistances(space = "composition", color.by = "batch", shape.by = "group", n.permutations = 99)
  expect_match(gg$labels$subtitle, "^batch: adj. R\u00b2 [0-9.]+, p [0-9.e-]+; group: adj. R\u00b2")
  D <- cao$getSampleDistanceMatrix(space = "composition")
  sc <- screenCovariates(list(m = D), cao$sample.meta, covariates = "batch", mode = "marginal", dist = "l2", n.permutations = 99,
                         min.samples.per.type = 4, seed = 1)
  expect_match(gg$labels$subtitle, sprintf("batch: adj. R\u00b2 %.2f", sc$table$R2.adj[1]), fixed = TRUE)
  expect_gt(sc$table$R2.adj[1], 0.3)                                    # composition is driven by the batch
  expect_null(cao$plotSampleDistances(color.by = "batch", show.r2 = FALSE)$labels$subtitle)
  # the map embeds the geometry of the test: sqrt of the correlation distance (a squared distance there), l2 as is
  Dj <- cao$getSampleDistanceMatrix(); expect_equal(attr(Dj, "dist.type"), "cor")
  xy <- cao$plotSampleDistances(color.by = "batch", show.r2 = FALSE)$data
  ref <- cmdscale(sqrt(pmax(unclass(Dj), 0)), k = 2)[xy$sample, ]
  expect_equal(abs(unname(cbind(xy$x, xy$y))), abs(unname(ref)), tolerance = 1e-8)
  xyc <- cao$plotSampleDistances(space = "composition", show.r2 = FALSE)$data
  refc <- cmdscale(unclass(D), k = 2)[xyc$sample, ]
  expect_equal(abs(unname(cbind(xyc$x, xyc$y))), abs(unname(refc)), tolerance = 1e-8)
  expect_match(cao$plotSampleDistances(color.by = "test")$labels$subtitle, "^group: adj. R")
})

test_that("stored shift results provide the expression distances; adjusted views need a model", {
  cao <- makeSpaceToy()
  suppressWarnings(res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 19, verbose = FALSE))
  D <- cao$getSampleDistanceMatrix(cell.type = "ct1")
  expect_equal(unclass(D), unclass(res$distances$ct1), ignore_attr = TRUE)
  Da <- cao$getSampleDistanceMatrix(values = "adjusted", cell.type = "ct1")
  expect_equal(dim(Da), dim(D))
  expect_true(inherits(cao$plotSampleDistances(values = "both", color.by = "group"), "ggplot"))
  # composition: adjusted for the model's covariates (none here) -> unadjusted with a message
  expect_equal(cao$getSampleDistanceMatrix(space = "composition", values = "adjusted"), cao$getSampleDistanceMatrix(space = "composition"))
  cao$setModel(formula = ~ group + batch, test = "group", verbose = FALSE)
  Dca <- cao$getSampleDistanceMatrix(space = "composition", values = "adjusted")
  Dcb <- cao$getSampleDistanceMatrix(space = "composition", adjust.for = ~ batch)
  expect_equal(unclass(Dca), unclass(Dcb), ignore_attr = TRUE)
  gg <- cao$plotSampleDistances(space = "composition", values = "adjusted", color.by = "batch", n.permutations = 99)
  expect_match(gg$labels$title, "model-adjusted")
  expect_match(gg$labels$subtitle, "batch: removed by the adjustment")   # no R2 for a covariate that was regressed out
  gg2 <- cao$plotSampleDistances(space = "composition", values = "adjusted", color.by = "group", n.permutations = 99)
  expect_match(gg2$labels$subtitle, "^group: adj. R")
})

test_that("the screen carries summary rows that stay out of the cell-type grid", {
  cao <- makeSpaceToy()
  sc <- cao$screenCovariates(n.permutations = 199, verbose = FALSE)
  expect_equal(sc$summaries, c("all cell types", "composition"))
  expect_setequal(sc$celltypes, c("ct1", "ct2"))
  tb <- sc$table
  expect_true(all(tb$summary[tb$celltype %in% sc$summaries])); expect_false(any(tb$summary[tb$celltype %in% sc$celltypes]))
  expect_equal(unique(sc$global$n.types), 2)                           # global p and counts over cell types only
  # BH among the summary rows only
  m <- tb[tb$mode == "marginal" & tb$summary, ]
  expect_equal(m$padj, p.adjust(m$p, "BH"))
  # the "all cell types" row is the screen of the combined distance (same as the deprecated metadata separation)
  Dj <- cao$getSampleDistanceMatrix()
  one <- screenCovariates(list(j = Dj), cao$sample.meta, covariates = sc$covariates, mode = "both", dist = attr(Dj, "dist.type"), n.permutations = 199, seed = sc$settings$seed)
  a <- tb[tb$celltype == "all cell types" & tb$mode == "marginal", ]; b <- one$table[one$table$mode == "marginal", ]
  expect_equal(a$R2.adj[order(a$covariate)], b$R2.adj[order(b$covariate)], tolerance = 1e-10)
  # composition row: batch drives the composition
  cb <- tb[tb$celltype == "composition" & tb$mode == "marginal" & tb$covariate == "batch", ]
  expect_lt(cb$p, 0.05)
  expect_output(print(sc), "composition (partial, FDR 5%): batch", fixed = TRUE)
  gg <- cao$plotCovariateScreen()
  expect_equal(levels(gg$data$celltype)[1:2], c("all cell types", "composition"))
  expect_true(any(vapply(gg$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))))
  expect_false(any(c("all cell types", "composition") %in% as.character(ggplot2::layer_data(cao$plotCovariateSummary())$x)))
  sc0 <- cao$screenCovariates(n.permutations = 49, verbose = FALSE, summary.rows = character(0), name = "s0")
  expect_length(sc0$summaries, 0); expect_false(any(sc0$table$summary))
  sc2 <- cao$screenCovariates(n.permutations = 19, verbose = FALSE, summary.rows = c("expression.equal", "pseudobulk"), name = "s2")
  expect_equal(sc2$summaries, c("all cell types (equal weights)", "pooled pseudobulk"))
  expect_error(cao$screenCovariates(summary.rows = "umap"), "unknown summary.rows")
})

test_that("a summary row whose partial mode fell back to marginal says so", {
  cao <- makeSpaceToy()
  Dj <- cao$getSampleDistanceMatrix(); attr(Dj, "summary") <- TRUE
  D1 <- cao$getSampleDistanceMatrix(cell.type = "ct1")
  sc <- screenCovariates(list("all cell types" = Dj, ct1 = D1), cao$sample.meta, covariates = c("group", "batch", "age"), mode = "both",
                         dist = "cor", n.permutations = 49, max.partial.df = 50, min.samples.per.type = 4, seed = 1)
  expect_true(all(sc$table$fallback[sc$table$mode == "partial"]))
  expect_output(print(sc), "all cell types (marginal, partial fell back, FDR 5%)", fixed = TRUE)
  expect_match(plotCovariateScreen(sc)$labels$subtitle, "1 of 1 cell types and 1 of 1 summary columns")
  sc2 <- screenCovariates(list("all cell types" = Dj, ct1 = D1), cao$sample.meta, covariates = c("group", "batch"), mode = "both",
                          dist = "cor", n.permutations = 49, min.samples.per.type = 4, seed = 1)
  expect_false(any(sc2$table$fallback)); expect_output(print(sc2), "all cell types (partial, FDR 5%)", fixed = TRUE)
})
