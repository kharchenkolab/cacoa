## The contrast tree and the reference block of the compositional model, on a composition with one planted change
## (one cell type doubled in group B, the others unchanged among themselves), and the CoDA plots on a toy object.

codaToy <- function(n = 12, seed = 1, numeric = FALSE) {
  set.seed(seed); ct <- paste0("T", 1:6)
  meta <- data.frame(row.names = paste0("s", seq_len(n)), cond = factor(rep(c("A", "B"), each = n / 2)), sex = factor(rep(c("f", "m"), n / 2)), age = round(seq(40, 70, length.out = n)))
  base <- c(400, 300, 150, 80, 50, 20)
  cnts <- t(sapply(seq_len(n), function(i) { p <- base; if (numeric) p[2] <- p[2] * exp(0.03 * (meta$age[i] - 55)) else if (meta$cond[i] == "B") p[2] <- p[2] * 2; stats::rmultinom(1, 1500, p) }))
  dimnames(cnts) <- list(rownames(meta), ct); list(cnts = cnts, meta = meta)
}

test_that("the contrast tree isolates the changed cell type at the root and its shares add up", {
  skip_if_not_installed("coda.base")
  d <- codaToy(); model <- buildCacoaModel(d$meta, formula = ~ cond, test = "cond")
  x <- cacoa:::lmCoda(d$cnts, model, n.permutations = 499, seed = 1)
  tr <- contrastTree(x)
  expect_equal(nrow(tr$nodes), 5); expect_equal(dim(tr$sbp), c(6, 5))
  expect_equal(tr$sides[[1]]$plus, "T2")                               # the root separates the doubled type from the rest
  expect_true(tr$nodes$share[1] > 0.6); expect_true(tr$nodes$padj[1] < 0.1); expect_true(tr$nodes$p[1] == min(tr$nodes$p))
  expect_equal(sum(tr$nodes$share), 1, tolerance = 1e-8)             # the balances of a sequential binary partition decompose the shift
  expect_equal(max(tr$nodes$cum.share), 1, tolerance = 1e-8)
  expect_true(all(abs(colSums(tr$sbp)) <= 6) && all(colSums(tr$sbp != 0) >= 2))
  expect_equal(tr$order, names(sort(x$contrast$loadings$obs)))
  expect_true(tr$nodes$log2.ratio[1] > 0.8)                          # about one doubling of T2 against the others
  expect_s3_class(tr$tree, "phylo"); expect_equal(length(tr$tree$tip.label), 6)
  # every node's effect is the balance of the loadings and its p-value comes from the same relabelings
  v <- cacoa:::balanceVector(tr$sides[[1]]$plus, tr$sides[[1]]$minus, names(x$contrast$loadings$obs))
  expect_equal(tr$nodes$effect[1], sum(v * x$contrast$loadings$obs))
  tr2 <- contrastTree(x, split = "evidence"); expect_equal(tr2$sides[[1]]$plus, "T2")
  expect_identical(x$tree$nodes$effect, tr$nodes$effect)             # lmCoda stores the same tree
})

test_that("the reference block is the set of unchanged cell types and the reference-relative p-values follow", {
  skip_if_not_installed("coda.base")
  d <- codaToy(); model <- buildCacoaModel(d$meta, formula = ~ cond, test = "cond")
  x <- cacoa:::lmCoda(d$cnts, model, n.permutations = 499, seed = 1)
  expect_equal(x$reference$method, "block")
  expect_setequal(x$reference$celltypes, c("T1", "T3", "T4", "T5", "T6"))
  expect_true(x$contrast$per_cell$padj.ref[["T2"]] < 0.05)
  expect_true(all(x$contrast$per_cell$pval.ref[c("T1", "T3", "T5")] > 0.2))
  expect_true(abs(mean(x$contrast$loadings$ref_centered[c("T1", "T3", "T4", "T5", "T6")])) < 1e-10)
  expect_true(x$contrast$loadings$ref_centered[["T2"]] > 0.5)
  # the previous rule is still available
  y <- cacoa:::lmCoda(d$cnts, model, n.permutations = 99, seed = 1, ref.method = "pvalue")
  expect_equal(y$reference$method, "pvalue"); expect_true(length(y$reference$celltypes) <= 3)
  # the pairwise p-values: within the block large, against T2 small
  pp <- x$reference$p.pairs; expect_true(all(pp["T2", c("T1", "T3", "T4", "T5", "T6")] < 0.05)); expect_true(pp["T1", "T3"] > 0.2)
})

test_that("the CoDA plots render on the lmCoda result, for a two-level and a numeric contrast", {
  skip_if_not_installed("coda.base")
  d <- codaToy(); model <- buildCacoaModel(d$meta, formula = ~ cond + sex, test = "cond")
  x <- cacoa:::lmCoda(d$cnts, model, n.permutations = 99, seed = 1)
  expect_s3_class(plotContrastTree(x), "ggplot")
  expect_s3_class(plotContrastTree(x, show.samples = FALSE, leaf.values = "clr", show.text = FALSE), "ggplot")
  expect_s3_class(plotCodaLoadings(x), "ggplot"); expect_s3_class(plotCodaLoadings(x, significance = "reference"), "ggplot")
  expect_s3_class(plotCodaLoadings(x, center = "none"), "ggplot")
  g <- ggplot2::ggplot_build(plotCodaLoadings(x, center = "none")); expect_true(any(abs(g$data[[1]]$xintercept) < 1e-8))   # zero = the geometric mean
  expect_s3_class(cacoa:::plotCodaCellTypes(x, value = "predicted"), "ggplot")
  pr <- x$contrast$predicted   # intervals of the fitted endpoint compositions
  expect_true(all(c("baseline.lo", "baseline.hi", "target.lo", "target.hi") %in% names(pr)))
  expect_true(all(pr$baseline.lo <= pr$baseline & pr$baseline <= pr$baseline.hi))
  expect_true(all(pr$target.lo <= pr$target & pr$target <= pr$target.hi))
  expect_s3_class(cacoa:::plotCodaCellTypes(x, value = "predicted", show.ci = FALSE), "ggplot")
  expect_s3_class(cacoa:::plotCodaCellTypes(x, value = "reference-centered"), "ggplot")
  expect_s3_class(cacoa:::plotCodaCellTypes(x, value = "volcano"), "ggplot")
  expect_s3_class(cacoa:::plotCodaGlobal(x), "ggplot"); expect_s3_class(cacoa:::plotCodaHeatmap(x), "ggplot")
  # the covariate adjustment of the samples' log-frequencies moves the samples to the endpoints' covariate values
  sl <- cacoa:::codaSampleLogFreq(x); expect_equal(dim(sl$logf), c(12, 6)); expect_equal(levels(sl$group), c("baseline", "target", "other"))
  expect_equal(as.character(sl$group), rep(c("baseline", "target"), each = 6))
  dn <- codaToy(numeric = TRUE); mn <- buildCacoaModel(dn$meta, formula = ~ age, test = "age")
  xn <- cacoa:::lmCoda(dn$cnts, mn, n.permutations = 99, seed = 1)
  expect_equal(xn$tree$sides[[1]]$plus, "T2")
  expect_s3_class(plotContrastTree(xn), "ggplot")
  expect_true(cacoa:::codaSampleLogFreq(xn)$numeric)
})

test_that("the class methods draw the CoDA plots", {
  skip_if_not_installed("coda.base"); skip_if_not_installed("psych")
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 60, n.genes = 40, n.cell.types = 4)
  res <- cao$estimateCellLoadings(n.permutations = 49, verbose = FALSE)
  expect_true(!is.null(res$tree)); expect_equal(nrow(res$tree$nodes), 3)
  expect_s3_class(cao$plotContrastTree(), "ggplot")
  expect_s3_class(cao$plotContrastTree(split = "evidence", show.samples = FALSE), "ggplot")
  expect_s3_class(cao$plotCodaContrastTree(), "ggplot")
  expect_s3_class(cao$plotCodaCellTypes(value = "predicted"), "ggplot")
  expect_s3_class(cao$plotCodaCellTypes(value = "reference-centered"), "ggplot")
  expect_s3_class(cao$plotCodaGlobal(), "ggplot"); expect_s3_class(cao$plotCodaHeatmap(), "ggplot")
  expect_s3_class(cao$plotCellLoadings(), "ggplot"); expect_s3_class(cao$plotCellLoadings(center = "none"), "ggplot")
  skip_if_not_installed("cowplot"); expect_true(inherits(cao$plotCellLoadings(center = "both"), "ggplot"))
  cao$estimateCellLoadings(test = "group: all", formula = ~ group, n.permutations = 19, verbose = FALSE, name = "coda.term")
  expect_error(cao$plotContrastTree(name = "coda.term"), "no contrast")
})

test_that("significanceContours outlines the region of the significant cells", {
  set.seed(2); n <- 4000
  emb <- cbind(x = c(rnorm(n / 2, -3), rnorm(n / 2, 3)), y = rnorm(n)); rownames(emb) <- paste0("c", seq_len(n))
  sig <- setNames(c(rep(TRUE, n / 2), rep(FALSE, n / 2)), rownames(emb))      # the left cloud is significant
  d <- significanceContours(emb, sig)
  expect_true(nrow(d) > 0); expect_true(all(d$x < 1))                           # the outline lies around the left cloud
  expect_equal(nrow(significanceContours(emb, setNames(rep(FALSE, n), rownames(emb)))), 0)
  expect_s3_class(cacoa:::significanceContourLayer(emb, sig), "LayerInstance")
})
