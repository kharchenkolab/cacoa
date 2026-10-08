# End-to-end smoke test on the small preprocessed panel shipped with cacoa (two CTRL, two DISEASE
# pagoda2 objects). Slow and needs conos; skipped on CRAN or without conos.
test_that("Cacoa object can be built from a conos panel", {
  skip_if_not_installed("conos")
  skip_on_cran()
  utils::data("panel.preprocessed", package = "cacoa", envir = environment())
  con <- conos::Conos$new(panel.preprocessed, n.cores = 1)
  con$buildGraph(n.odgenes = 500)
  con$findCommunities()
  meta <- data.frame(condition = factor(sub("[0-9]+$", "", names(con$samples))),
                     row.names = names(con$samples))
  cao <- Cacoa$new(data.object = con, sample.metadata = meta,
                   contrast = c("condition", "DISEASE", "CTRL"),
                   cell.groups = con$clusters$leiden$groups, n.cores = 1, verbose = FALSE)
  expect_s3_class(cao, "Cacoa")
  expect_equal(nrow(cao$sample.meta), 4)
  expect_equal(colnames(cao$model$F), c("conditionCTRL", "conditionDISEASE"))
})
