# Dispersion design: sample properties (nuisance covariates, the per-sample 1 / n.cells sampling-noise term) stay with the sample under
# relabeling; the dispersion change is taken at equal cell numbers.

simCellDistances <- function(n.per = 10, dims = 200, tau2 = 4, sigma2 = 0.05, cells.A = 150, cells.B = 50, seed = 1) {
  set.seed(seed)   # cell numbers vary several-fold within a condition (as in real data), B with fewer on average
  g <- rep(c("A", "B"), each = n.per); nc <- round(ifelse(g == "A", cells.A, cells.B) * exp(stats::rnorm(2 * n.per, 0, 0.8))) + 5
  Y <- t(vapply(seq_along(g), function(i) stats::rnorm(dims, 0, sqrt(sigma2 + tau2 / nc[i])), numeric(dims)))   # no condition effect
  D <- as.matrix(stats::dist(Y))^2 / dims   # squared distances (the "cor" route uses them as given)
  rownames(D) <- colnames(D) <- names(nc) <- paste0("s", seq_along(g))
  attr(D, "n.cells") <- nc
  list(D = D, meta = data.frame(row.names = rownames(D), group = factor(g)), n.cells = nc)
}

test_that("fixed dispersion columns: C++ kernels match the R reference and keep the sample's values", {
  s <- simCellDistances(seed = 2); m <- buildCacoaModel(s$meta, formula = ~ group, test = "group")
  e <- cacoa:::pairwiseEffectsFromDesign(s$D, m$tests[[1]]$design, s$meta, dist = "cor")
  expect_true(e$cell.number.adjust); expect_equal(colnames(e$Z)[ncol(e$Z)], ".inv.n.cells"); expect_equal(e$z.fixed, ncol(e$Z))
  pre <- cacoa:::inferencePrecompute(e$G, e$X, e$Z, e$contrast)
  set.seed(3); P <- replicate(20, sample.int(nrow(e$X)))
  ref <- t(apply(P, 2, function(p) cacoa:::permutedStats(e$G, e$X, e$Z, e$contrast, e$z.end, pre, p, TRUE, TRUE, z.fixed = e$z.fixed)))
  cpp <- cacoa:::permuted_contrast_stats(e$G, e$X, e$Z, pre$A, pre$H, pre$a, e$contrast, pre$cXc, pre$q, e$z.end$num, e$z.end$den, P,
                                         as.integer(e$z.fixed) - 1L, TRUE, TRUE)
  expect_equal(unname(cpp), unname(ref), tolerance = 1e-8)   # including the studentized dispersion change (column 5)
  expect_true(all(is.finite(cpp[, 5])))
  cw <- cacoa:::permuted_contrast_stats_w(e$G, e$X, e$Z, rep(1, nrow(e$X)), e$contrast, e$z.end$num, e$z.end$den, P, FALSE, rep(1, nrow(e$X)),
                                          as.integer(e$z.fixed) - 1L, TRUE, TRUE, 0L)
  expect_equal(unname(cw), unname(ref), tolerance = 1e-6)
  # the fixed column is not permuted: the same result as giving the kernel a Z whose fixed column follows the permutation back
  noFix <- cacoa:::permuted_contrast_stats(e$G, e$X, e$Z, pre$A, pre$H, pre$a, e$contrast, pre$cXc, pre$q, e$z.end$num, e$z.end$den, P, integer(0), TRUE, TRUE)
  expect_false(isTRUE(all.equal(as.numeric(noFix[, 3]), as.numeric(cpp[, 3]))))
  # nuisance covariates of the dispersion formula are fixed too; the tested variable is not
  s$meta$batch <- factor(rep(c("b1", "b2"), 10))
  e2 <- cacoa:::pairwiseEffectsFromDesign(s$D, m$tests[[1]]$design, s$meta, dispersion.formula = ~ group + batch, dist = "cor", cell.number.adjust = FALSE)
  expect_equal(colnames(e2$Z)[e2$z.fixed], "batchb2"); expect_false(e2$cell.number.adjust)
})

test_that("var at equal cell numbers: no spurious dispersion change when one condition has fewer cells", {
  v <- t(vapply(1:25, function(r) {
    s <- simCellDistances(seed = 100 + r); m <- buildCacoaModel(s$meta, formula = ~ group, test = "group")
    a <- cacoa:::pairwiseEffectsFromDesign(s$D, m$tests[[1]]$design, s$meta, dist = "cor")
    b <- cacoa:::pairwiseEffectsFromDesign(s$D, m$tests[[1]]$design, s$meta, dist = "cor", cell.number.adjust = FALSE)
    c(adjusted = a$var.norm, unadjusted = b$var.norm)
  }, numeric(2)))
  expect_gt(mean(v[, "unadjusted"]), 0.5)               # B (fewer cells) looks more dispersed without the term
  expect_lt(abs(mean(v[, "adjusted"], na.rm = TRUE)), 0.2)   # ... and not with it
  expect_lt(mean(is.na(v[, "adjusted"])), 0.1)
  # through the test: the var p-values are calibrated with the term, not without it
  p <- t(vapply(1:40, function(r) {
    s <- simCellDistances(seed = 300 + r); m <- buildCacoaModel(s$meta, formula = ~ group, test = "group")
    D.list <- list(ct = s$D)
    c(adj = testPairwiseEffects(D.list, m$tests[[1]]$design, s$meta, dist = "cor", n.permutations = 99, seed = r)$results$p.var,
      raw = testPairwiseEffects(D.list, m$tests[[1]]$design, s$meta, dist = "cor", n.permutations = 99, seed = r, cell.number.adjust = FALSE)$results$p.var)
  }, numeric(2)))
  expect_lt(mean(p[, "adj"] < 0.05, na.rm = TRUE), 0.15)
  expect_gt(mean(p[, "raw"] < 0.05), 0.5)
})

test_that("cell numbers tied to the contrast are measured and flag the dispersion test", {
  s <- simCellDistances(seed = 11, cells.A = 300, cells.B = 30)   # cell numbers strongly associated with the condition
  m <- buildCacoaModel(s$meta, formula = ~ group, test = "group")
  e <- cacoa:::pairwiseEffectsFromDesign(s$D, m$tests[[1]]$design, s$meta, dist = "cor")
  expect_equal(e$cell.assoc, abs(stats::cor(as.numeric(s$meta$group == "B"), log(s$n.cells))), tolerance = 1e-8)
  r <- testPairwiseEffects(list(ct = s$D), m$tests[[1]]$design, s$meta, dist = "cor", n.permutations = 49, seed = 1)
  expect_gt(r$results$cell.assoc, 0.5); expect_true(grepl("cells-tied-to-contrast", r$results$flags))
  D0 <- s$D; attr(D0, "n.cells") <- NULL
  r0 <- testPairwiseEffects(list(ct = D0), m$tests[[1]]$design, s$meta, dist = "cor", n.permutations = 49, seed = 1)
  expect_true(is.na(r0$results$cell.assoc)); expect_false(grepl("cells-tied", r0$results$flags))
})
