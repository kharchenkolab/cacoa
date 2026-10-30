# Gene filters count a gene as detected in a cell when its count is above 0 (single-UMI detections included).

test_that("a gene detected only as single UMIs is kept by the DE filter and by the pseudobulk gene list", {
  set.seed(1)
  n.s <- 4; n.c <- 50; cells <- sprintf("c%03d", seq_len(n.s * n.c))
  spc <- setNames(factor(rep(sprintf("s%d", seq_len(n.s)), each = n.c)), cells)
  cg <- setNames(factor(rep(c("A", "B"), length.out = length(cells))), cells)
  m <- matrix(0, 3, length(cells), dimnames = list(c("single", "rare", "absent"), cells))
  m["single", sample(length(cells), 0.4 * length(cells))] <- 1   # 40 % of cells, always one UMI
  m["rare", 1:2] <- 5                                             # two cells only
  M <- Matrix::Matrix(m, sparse = TRUE)
  f <- cacoa:::getPerCellTypeGeneFilter(M, cg, threshold = 0.05)
  expect_true(all(vapply(f, function(x) x[["single"]], logical(1))))
  expect_false(any(vapply(f, function(x) x[["absent"]] || x[["rare"]], logical(1))))
  cms <- lapply(split(cells, spc), function(cc) Matrix::t(M[, cc, drop = FALSE]))
  pb <- pseudobulkPerCellType(cms, cell.groups = cg, sample.per.cell = spc, min.cells.per.sample = 5, min.gene.frac = 0.01, scale = 1e3)
  expect_true("single" %in% pb$genes); expect_false("absent" %in% pb$genes)
})
