# One gene set per cell type (typeGeneSets(): the common list's genes detected in >= type.gene.frac of the type's cells), used by
# the cell-type gsd test and the cell-type sample distances; the cluster-free maps apply the same rule per neighbourhood.

# toy whose first `n.absent` genes are rare in ct1 (detected in about 2 % of its cells: sporadic counts, not exact zeros, which the
# gene-centred distance would ignore anyway)
makeTypeGeneToy <- function(n.absent = 15, n.genes = 100, seed = 4) {
  set.seed(seed)
  groups <- rep(c("A", "B"), each = 6); samples <- sprintf("%s%d", groups, rep(1:6, 2))
  meta <- data.frame(group = factor(groups), batch = factor(rep(c("b1", "b2"), 6)), row.names = samples)
  spc <- factor(rep(samples, each = 60)); cg <- factor(ifelse(stats::runif(length(spc)) < 0.5, "ct1", "ct2"))
  names(spc) <- names(cg) <- sprintf("cell%04d", seq_along(spc))
  mu <- outer(exp(rnorm(n.genes, 1, 0.7)), rep(1, length(cg)))
  mu[seq_len(n.absent), cg == "ct1"] <- 0.02
  isB <- meta[as.character(spc), "group"] == "B"; sh <- n.absent + seq_len(20); mu[sh, isB] <- mu[sh, isB] * exp(1)
  cm <- Matrix::Matrix(matrix(rpois(length(mu), mu), n.genes, dimnames = list(sprintf("g%03d", seq_len(n.genes)), names(spc))), sparse = TRUE)
  suppressWarnings(suppressMessages(Cacoa$new(data.object = cm, sample.metadata = meta, sample.per.cell = spc, cell.groups = cg,
                                              contrast = c("group", "B", "A"), n.cores = 1, verbose = FALSE)))
}

test_that("pseudobulk detection fractions and the per-type gene sets", {
  cao <- makeTypeGeneToy()
  cms <- cacoa:::extractRawCountMatrices(cao$data.object, transposed = TRUE)
  pb <- pseudobulkPerCellType(cms, cao$cell.groups, cao$sample.per.cell)
  cm <- Matrix::t(cao$data.object); g <- pb$genes[c(1, 20, 50)]
  for (ct in c("ct1", "ct2")) {
    cells <- names(cao$cell.groups)[cao$cell.groups == ct]
    expect_equal(unname(pb$det.frac[ct, g]), unname(Matrix::rowMeans(cm[g, cells] > 0)), tolerance = 1e-12)
  }
  gs <- typeGeneSets(pb, 0.05)
  expect_false(any(sprintf("g%03d", 1:15) %in% gs$ct1)); expect_true(all(sprintf("g%03d", 1:15) %in% gs$ct2))
  expect_equal(typeGeneSets(pb, 0)$ct1, pb$genes)
})

test_that("cell-type distances and the gsd test use the same per-type genes", {
  cao <- makeTypeGeneToy()
  cms <- cacoa:::extractRawCountMatrices(cao$data.object, transposed = TRUE)
  pb <- pseudobulkPerCellType(cms, cao$cell.groups, cao$sample.per.cell); gs <- typeGeneSets(pb, 0.05)
  D <- cao$getSampleDistanceMatrix(cell.type = "ct1")
  ref <- sampleDistanceMatrices(list(ct1 = pb$cm.per.type$ct1[, gs$ct1]), dist = "cor")$ct1
  expect_equal(unclass(D)[rownames(ref), colnames(ref)], ref, ignore_attr = TRUE)
  expect_equal(attr(D, "n.genes"), length(gs$ct1))
  # the option switches the rule off (common list for every type)
  cao$setOptions(type.gene.frac = 0)
  D0 <- cao$getSampleDistanceMatrix(cell.type = "ct1")
  ref0 <- sampleDistanceMatrices(list(ct1 = pb$cm.per.type$ct1), dist = "cor")$ct1
  expect_equal(unclass(D0)[rownames(ref0), colnames(ref0)], ref0, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(unclass(D), unclass(D0), check.attributes = FALSE)))
  cao$setOptions(type.gene.frac = 0.05)
  # the shift result: gsd genes of a type are within its gene set, the distances are the per-type ones
  r <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(n.permutations = 49, verbose = FALSE))
  expect_equal(r$settings$type.gene.frac, 0.05)
  if (identical(r$settings$statistic, "gsd")) {
    g1 <- r$genes[[1]]$ct1$gene
    expect_true(all(g1 %in% gs$ct1)); expect_false(any(sprintf("g%03d", 1:15) %in% g1))
  }
  expect_equal(unclass(r$distances$ct1), unclass(D), ignore_attr = TRUE)
  # the earlier argument name is accepted
  r2 <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(n.permutations = 19, neighbourhood.gene.frac = 0, verbose = FALSE, name = "r2"))
  expect_equal(r2$settings$type.gene.frac, 0)
  expect_error(cao$setOptions(type.gene.frac = 1.5), "type.gene.frac")
})
