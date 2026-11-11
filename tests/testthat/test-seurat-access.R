## The counts read from a Seurat object are the raw counts, whichever version of SeuratObject is installed (versions before 5 name the
## argument `slot` and return the normalized data when given `layer`).

test_that("raw counts of a Seurat object are read from the counts layer", {
  skip_if_not_installed("Seurat")
  set.seed(1)
  cnt <- Matrix::Matrix(matrix(rpois(40 * 30, 2), 40, 30, dimnames = list(paste0("g", 1:40), paste0("c", 1:30))), sparse = TRUE)
  so <- suppressWarnings(Seurat::NormalizeData(Seurat::CreateSeuratObject(cnt), verbose = FALSE))
  so$sample.per.cell <- rep(c("s1", "s2", "s3"), each = 10); so@misc$assay.name <- "RNA"; so@misc$data.layer <- "data"
  cms <- expect_silent(cacoa:::extractRawCountMatrices.Seurat(so, transposed = FALSE))
  expect_equal(sum(vapply(cms, sum, numeric(1))), sum(cnt))
  expect_equal(as.matrix(cms$s2), as.matrix(cnt[, 11:20]))
  jc <- cacoa:::extractJointCountMatrix.Seurat(so, raw = TRUE, transposed = FALSE)
  expect_equal(as.matrix(jc), as.matrix(cnt))
  nd <- cacoa:::extractJointCountMatrix.Seurat(so, raw = FALSE, transposed = FALSE)
  expect_false(isTRUE(all.equal(as.matrix(nd), as.matrix(cnt))))        # the normalized layer is a different matrix
  expect_lt(max(nd), max(cnt))
  # values that are not whole numbers in the counts layer are reported
  expect_warning(cacoa:::seuratAssayData(suppressWarnings(Seurat::CreateSeuratObject(cnt * 0.5)), "counts", assay = "RNA"), "not whole numbers")
})
