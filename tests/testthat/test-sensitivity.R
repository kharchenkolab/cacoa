test_that("composition sensitivity: alternative covariate sets, jackknife and leave-one-out ranges", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.genes = 60, n.cell.types = 4, shift = 1.5)
  res <- suppressMessages(cao$estimateCellLoadings(n.permutations = 49, verbose = FALSE))
  # the leave-one-out refits use the F-space contrast of the ILR fit: at the full data it is the stored loading
  m <- res$model; co <- qr.coef(qr(m$F), res$ilr); co[is.na(co)] <- 0
  l.full <- as.numeric(res$psi %*% as.numeric(as.numeric(m$contrast.F) %*% co))
  expect_equal(l.full, unname(res$contrast$loadings$obs), tolerance = 1e-8)
  s <- suppressMessages(cao$checkSensitivity(analysis = "composition", covariate.sets = list(unadjusted = ~ group, `plus batch` = ~ group + batch),
                                             n.permutations = 49, verbose = FALSE))
  expect_s3_class(s, "cacoaSensitivity")
  expect_identical(s$analysis, "composition"); expect_identical(s$effect, "loading")
  expect_equal(levels(s$table$model), c("current", "plus batch"))      # the stored model is "current" whatever it is called
  expect_true(all(c("estimate", "se.jk", "loo.min", "loo.max", "p", "padj") %in% names(s$table)))
  expect_equal(nrow(s$table), 2 * 4)
  cur <- s$table[s$table$model == "current", ]
  expect_equal(cur$estimate, unname(res$contrast$loadings$obs[cur$celltype]) / log(2), tolerance = 1e-8)   # loadings in log2
  expect_equal(cur$padj, unname(res$contrast$per_cell$padj[cur$celltype]))                                 # the stored test, same seed
  expect_true(all(is.finite(s$table$se.jk) & s$table$se.jk > 0))                                           # jackknife under every model
  expect_true(all(s$table$loo.min <= s$table$estimate + 1e-8 & s$table$loo.max >= s$table$estimate - 1e-8))
  expect_true(all(c("celltype", "verdict.magnitude", "verdict.significance", "band.lo", "band.hi") %in% names(s$summary)))
  expect_equal(nrow(s$summary), 4)
  expect_true(all(is.finite(s$summary$band.lo)))                                                            # the stored result's permutation band
  expect_identical(s, cao$test.results$composition.sensitivity)
  expect_output(print(s), "loading")
  expect_s3_class(cao$plotSensitivity(analysis = "composition"), "ggplot")
  expect_s3_class(cao$plotSensitivity(analysis = "composition", interval = "none"), "ggplot")
  # without the jackknife: no intervals, same estimates
  s0 <- suppressMessages(cao$checkSensitivity(analysis = "composition", covariate.sets = list(unadjusted = ~ group, `plus batch` = ~ group + batch),
                                              n.permutations = 49, jackknife = FALSE, verbose = FALSE, name = "cs0"))
  expect_true(all(is.na(s0$table$se.jk)))
  expect_equal(s0$table$estimate, s$table$estimate)
  # direct call
  s2 <- checkCompositionSensitivity(res$cnts, res$model, cao$sample.meta, formulas = list(current = ~ group), n.permutations = 29, seed = 1)
  expect_equal(nrow(s2$table), 4)
})

test_that("expression-shift sensitivity: jackknife intervals and leave-one-out ranges for every model", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 40, n.genes = 80, shift = 1.5)
  cao$estimateExpressionShiftMagnitudes(n.permutations = 29, verbose = FALSE)
  s <- suppressWarnings(cao$checkSensitivity(covariate.sets = list(unadjusted = ~ group, `plus batch` = ~ group + batch), n.permutations = 29, verbose = FALSE))
  expect_identical(s$analysis, "expression.shifts"); expect_true(isTRUE(s$settings$jackknife))
  d <- s$table[s$table$effect == "shift", ]
  expect_equal(nrow(d), 4)
  expect_true(all(is.finite(d$se.jk)))
  expect_true(all(is.finite(d$loo.min) & is.finite(d$loo.max) & d$loo.min <= d$loo.max))
  expect_s3_class(cao$plotSensitivity(), "ggplot")
  expect_s3_class(cao$plotSensitivity(interval = "none"), "ggplot")
  s0 <- suppressWarnings(cao$checkSensitivity(covariate.sets = list(unadjusted = ~ group, `plus batch` = ~ group + batch), n.permutations = 29,
                                              jackknife = FALSE, verbose = FALSE, name = "s0"))
  d0 <- s0$table[s0$table$effect == "shift", ]
  expect_true(all(is.finite(d0$loo.min[d0$model == "current"])))                                   # the stored influence of the current model
  expect_equal(d0$estimate[d0$model == "current"], d$estimate[d$model == "current"])
})
