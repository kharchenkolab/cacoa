test_that("DESeq2 path: dispersions are estimated around the model fit, as in DESeq2::DESeq() with the same design", {
  skip_if_not_installed("DESeq2")
  set.seed(3); n.genes <- 300; n <- 12
  meta <- data.frame(row.names = paste0("s", seq_len(n)), lab = factor(rep(c("A", "B"), each = n / 2)), batch = factor(rep(c("b1", "b2", "b3"), times = n / 3)))
  mu <- matrix(rgamma(n.genes, shape = 2, scale = 60), n.genes, n)
  mu[1:60, meta$lab == "B"] <- mu[1:60, meta$lab == "B"] * 3                      # label effect in a fifth of the genes
  mu[, meta$batch == "b2"] <- mu[, meta$batch == "b2"] * rep(c(1, 2), length.out = n.genes)
  cm <- matrix(rnbinom(n.genes * n, mu = mu, size = 5), n.genes, n, dimnames = list(paste0("g", seq_len(n.genes)), rownames(meta))); storage.mode(cm) <- "integer"
  for (f in list(~ lab, ~ lab + batch)) {
    X <- model.matrix(f, meta); cF <- setNames(as.numeric(colnames(X) == "labB"), colnames(X))
    res <- cacoa:::estimateDEForTypeDESeq(cm, meta, list(F = X, contrast.F = cF), test.type = "wald", independentFiltering = FALSE)
    ref <- suppressMessages(DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(cm, meta, f), quiet = TRUE))
    ref <- as.data.frame(DESeq2::results(ref, name = "lab_B_vs_A", independentFiltering = FALSE, cooksCutoff = FALSE))
    expect_equal(res$pvalue, ref$pvalue, tolerance = 1e-6)
    expect_equal(res$log2FoldChange, ref$log2FoldChange, tolerance = 1e-6)
    expect_gt(sum(res$padj[1:60] < 0.05, na.rm = TRUE), 30)
  }
})
