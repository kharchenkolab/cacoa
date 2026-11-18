# Cell-type expression shifts with the per-gene statistic gsd (gsdShiftsPerCellType(), estimateExpressionShiftMagnitudes(statistic = "gsd")).

test_that("cell-type gsd: gene z-scores reproduce the kernel, magnitudes follow their formulas, thread-independent", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 50, n.genes = 100, shift = 1.0)
  cm <- Matrix::t(cao$data.object); lib <- Matrix::colSums(cm) * 2
  r1 <- suppressWarnings(gsdShiftsPerCellType(cm, cao$sample.per.cell, cao$cell.groups, cao$model, lib, n.permutations = 99, seed = 1, n.cores = 1))
  r2 <- suppressWarnings(gsdShiftsPerCellType(cm, cao$sample.per.cell, cao$cell.groups, cao$model, lib, n.permutations = 99, seed = 1, n.cores = 2))
  expect_equal(r1$results, r2$results)
  res <- r1$results
  expect_setequal(res$celltype, c("ct1", "ct2")); expect_true(all(res$effect == "shift"))
  expect_true(all(res$p < 0.05)); expect_true(all(res$estimate > 0)); expect_true(all(res$se.jk > 0 & res$se.jk.norm > 0))
  expect_true(all(res$ci.low < res$estimate & res$estimate < res$ci.high))
  # the kernel's per-type shift for gsd is the mean squared gene z: the R gene table reproduces it
  units <- split(seq_len(ncol(cm)) - 1L, droplevels(as.factor(cao$cell.groups[colnames(cm)])))
  k <- suppressWarnings(clusterFreeExpressionShifts(cm, cao$sample.per.cell, units, cao$model$tests[[1]]$design, cao$model$meta, lib = lib, n.permutations = 9,
                                                    seed = 1, post = "none", adjust = FALSE, smooth = FALSE, min.n.obs.per.samp = 5, min.samp.per.level = 3, n.fit.units = 2, gsd.d0 = 3))
  for (ct in names(units)) expect_equal(mean(r1$genes[[1]][[ct]]$z^2), unname(k$shifts[ct]), tolerance = 1e-8)
  # two groups: beta is the difference of the group means and the magnitude sum(beta^2 - (1/nA + 1/nB) s2) over the type's genes
  g <- r1$genes[[1]]$ct1; grp <- cao$sample.meta[levels(droplevels(as.factor(cao$sample.per.cell))), "group"]
  prof <- cacoa:::cf_unit_profiles(cm[g$gene, , drop = FALSE], as.integer(factor(as.character(cao$sample.per.cell[colnames(cm)]), levels = rownames(cao$sample.meta))),
                                   list(units$ct1), as.numeric(lib), nrow(cao$sample.meta), 5L, matrix(c(0L, 1L), 1),
                                   c(list(scale = r1$settings$scale, dist = 2L, noise = TRUE, tau = TRUE, correct = FALSE, unit_detection = 0),
                                     with(r1$noise.curves, list(mu_grid = mu.grid, phi = phi, tau2 = tau2)),
                                     with(cacoa:::shiftPoissonTable(), list(lam_grid = lam.grid, s_grid = s.grid, V1 = V1))))[[1]]
  x <- prof$x; gb <- cao$sample.meta$group == "B"
  beta <- rowMeans(x[, gb, drop = FALSE]) - rowMeans(x[, !gb, drop = FALSE])
  expect_equal(unname(beta[match(g$gene, rownames(cm[g$gene, ]))]), g$beta, tolerance = 1e-10)
  s2 <- (rowSums((x[, gb] - rowMeans(x[, gb]))^2) + rowSums((x[, !gb] - rowMeans(x[, !gb]))^2)) / (ncol(x) - 2)
  expect_equal(res$estimate[res$celltype == "ct1"], sum(beta^2 - (1 / sum(gb) + 1 / sum(!gb)) * s2), tolerance = 1e-10)
  # without a change the magnitude is near zero
  cao0 <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 50, n.genes = 100, shift = 0)
  cm0 <- Matrix::t(cao0$data.object)
  r0 <- suppressWarnings(gsdShiftsPerCellType(cm0, cao0$sample.per.cell, cao0$cell.groups, cao0$model, Matrix::colSums(cm0) * 2, n.permutations = 99, seed = 1))
  expect_true(all(abs(r0$results$estimate) < 0.1 * min(res$estimate)))
})

test_that("estimateExpressionShiftMagnitudes(statistic = 'gsd'): gsd significance, distance-model magnitudes, plots, whole-factor fallback", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 50, n.genes = 100, shift = 1.0)
  r <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(n.permutations = 99, verbose = FALSE))   # gsd is the default
  expect_s3_class(r, "cacoaExpressionShifts")
  expect_equal(r$settings$statistic, "gsd"); expect_equal(r$settings$min.cells.per.sample, 5)
  # the shift's significance is the gsd test's; magnitudes, var and total are the distance model's (its own shift p kept)
  g <- suppressWarnings(gsdShiftsPerCellType(Matrix::t(cao$data.object)[r$genes.tested, ], cao$sample.per.cell, cao$cell.groups, cao$model,
                                             Matrix::rowSums(cao$data.object)[colnames(Matrix::t(cao$data.object))], n.permutations = 99,
                                             seed = r$settings$gsd$seed, influence = FALSE, n.null = 0))
  w <- r$wide[[1]]
  expect_true(all(c("shift.norm", "var.norm", "p.shift", "p.shift.distance", "p.var", "null.shift.lo", "null.shift.hi", "null.var.lo", "null.var.hi", "n.genes") %in% names(w)))
  expect_equal(w$p.shift, g$results$p[match(w$celltype, g$results$celltype)])
  sh <- r$results[r$results$effect == "shift", ]
  expect_equal(sh$p, w$p.shift[match(sh$celltype, w$celltype)]); expect_true(all(is.finite(sh$p.distance)))
  expect_equal(sh$estimate.norm, w$shift.norm[match(sh$celltype, w$celltype)])
  expect_true(all(r$results$padj[r$results$effect == "shift"] < 0.05))
  expect_true(all(w$null.shift.lo < w$null.shift.hi & w$shift.norm > w$null.shift.hi))   # a strong shift lies above the relabeling band
  expect_true(all(c("var", "total") %in% r$results$effect))
  expect_true(is.data.frame(r$genes[[1]]$ct1) && all(c("gene", "beta", "se", "z") %in% names(r$genes[[1]]$ct1)))
  expect_false(is.unsorted(-abs(r$genes[[1]]$ct1$z)))
  expect_output(print(r), "statistic gsd")
  # plots: one effect per plot, options
  expect_s3_class(cao$plotExpressionShiftMagnitudes(), "ggplot")
  expect_s3_class(cao$plotExpressionShiftMagnitudes(effect = "var", show.ci = TRUE, show.pvalues = "raw"), "ggplot")
  expect_s3_class(cao$plotExpressionShiftMagnitudes(show.genes = TRUE, gene.fdr = 0.05, show.null = FALSE, legend.position = "none"), "ggplot")
  expect_s3_class(cao$plotSampleInfluence(), "ggplot")
  expect_s3_class(cao$plotSampleDistances(space = "expression.shifts", values = "unadjusted", color.by = "group"), "ggplot")
  expect_s3_class(cao$plotShiftDetail("ct1"), "ggplot")
  expect_s3_class(cao$plotExpressionShiftGenes(), "ggplot"); expect_s3_class(cao$plotExpressionShiftGenes(cell.types = "ct2", n.genes = 3), "ggplot")
  s <- suppressMessages(cao$checkSensitivity(covariate.sets = list(base = ~ group), n.permutations = 49, verbose = FALSE))
  expect_s3_class(s, "cacoaSensitivity")
  expect_error(cao$estimateExpressionShiftMagnitudes(statistic = "gsd", robust = "huber", n.permutations = 9, verbose = FALSE), "winsor")
  rw <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(statistic = "gsd", robust = "winsor", n.permutations = 99, verbose = FALSE, name = "w"))
  expect_equal(rw$settings$gsd$robust, "winsor"); expect_true(all(is.finite(rw$results$p)))
  # two contrasts of a three-level factor and a numeric contrast with a covariate
  cao$sample.meta$f3 <- factor(rep(c("x", "y", "z"), length.out = nrow(cao$sample.meta)))
  r3 <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(formula = ~ f3, test = c("f3: y vs x", "f3: z vs x"), n.permutations = 99, verbose = FALSE, name = "t"))
  tb <- table(r3$results$test, r3$results$effect)
  expect_true(all(tb[, "shift"] == 2) && all(is.finite(r3$results$p.distance[r3$results$effect == "shift"])))
  expect_s3_class(cao$plotExpressionShiftMagnitudes(name = "t", test = "f3: z vs x"), "ggplot")
  set.seed(1); cao$sample.meta$age <- ifelse(cao$sample.meta$group == "B", 60, 40) + stats::rnorm(nrow(cao$sample.meta), 0, 5)
  rn <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(formula = ~ age + f3, test = "age", n.permutations = 99, verbose = FALSE, name = "num"))
  expect_true(all(rn$results$padj[rn$results$effect == "shift"] < 0.05))
  expect_s3_class(cao$plotExpressionShiftMagnitudes(name = "num", effect = "var"), "ggplot")
})

test_that("cell-type gsd: relabelings are enumerated when there are fewer than n.permutations (exact p, per-type floor)", {
  cao <- makeToyCacoa(n.per.group = c(A = 4, B = 4), cells.per.sample = 50, n.genes = 100, shift = 1.0)
  cm <- Matrix::t(cao$data.object); lib <- Matrix::colSums(cm) * 2
  r <- suppressWarnings(gsdShiftsPerCellType(cm, cao$sample.per.cell, cao$cell.groups, cao$model, lib, n.permutations = 999, seed = 1, influence = FALSE, n.null = 0))
  res <- r$results
  expect_true(all(res$n.perm.distinct == choose(8, 4)))       # 70 relabelings of 4 + 4 samples
  expect_equal(res$p.floor, rep(1 / 70, nrow(res)))
  # exact p = (relabelings at least as extreme, the identity and its mirror included) / 70: a multiple of 1/70
  expect_equal(res$p * 70, round(res$p * 70), tolerance = 1e-8); expect_true(all(res$p >= 1 / 70))
  expect_true(all(res$flags == "few-permutations"))
})

test_that("printed p-values below alpha are never rounded up to it, and the count uses the unrounded values", {
  expect_equal(cacoa:::signifBelow(c(0.0497, 0.0503, 0.0012, NA), 0.05), c(0.0497, 0.05, 0.0012, NA))
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 50, n.genes = 100, shift = 1.0)
  r <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(n.permutations = 99, verbose = FALSE))
  r$results$padj[r$results$effect == "shift"][1] <- 0.0497
  out <- capture.output(print(r))
  expect_true(any(grepl("0.0497", out))); expect_true(any(grepl("2 of 2 cell types with a significant shift", out)))
})

test_that("magnitudes use samples with >= min.cells.magnitude cells; thinner types are tested without a magnitude", {
  cao <- makeToyCacoa(n.per.group = c(A = 5, B = 5), cells.per.sample = 60, n.genes = 400, n.cell.types = 3, shift = 0.6)
  cg <- cao$cell.groups; spc <- cao$sample.per.cell
  for (s in levels(spc)) { i <- which(spc == s & cg == "ct3"); cg[i[-(1:7)]] <- "ct1" }   # ct3: 7 cells per sample
  cao$cell.groups <- droplevels(cg)
  r <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(n.permutations = 99, seed = 1, verbose = FALSE))
  expect_equal(r$settings$min.cells.per.sample, 5); expect_equal(r$settings$min.cells.magnitude, 10)
  w <- r$wide[[1]]
  expect_true(is.finite(w$padj.shift[w$celltype == "ct3"]) && is.na(w$shift.norm[w$celltype == "ct3"]))
  expect_equal(w$n.gsd[w$celltype == "ct3"], 10)
  expect_false("ct3" %in% r$skipped$celltype)
  out <- capture.output(print(r)); expect_true(any(grepl("tested without a magnitude.*ct3", out)))
  expect_s3_class(cao$plotExpressionShiftMagnitudes(), "ggplot")
  # min.cells.magnitude = 5: the magnitude is estimated on the same samples as the test
  r5 <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(n.permutations = 99, seed = 1, verbose = FALSE, min.cells.magnitude = 5, name = "m5"))
  expect_true(is.finite(r5$wide[[1]]$shift.norm[r5$wide[[1]]$celltype == "ct3"]))
})

test_that("sensitivity of a gsd result: distance-model magnitudes, significance by gsd refitted under every model", {
  cao <- makeToyCacoa(n.per.group = c(A = 6, B = 6), cells.per.sample = 50, n.genes = 100, shift = 1.0)
  cao$sample.meta$nuisance <- factor(rep(c("u", "v"), 6))
  r <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(n.permutations = 99, verbose = FALSE))
  expect_equal(r$settings$statistic, "gsd")
  s <- suppressWarnings(cao$checkSensitivity(covariate.sets = list(`plus nuisance` = ~ group + nuisance), n.permutations = 49, verbose = FALSE))
  expect_equal(s$settings$significance, "gsd")
  tb <- s$table[s$table$effect == "shift", ]; sh <- r$results[r$results$effect == "shift", ]
  cur <- tb[tb$model == "current", ]
  expect_equal(cur$padj, sh$padj[match(cur$celltype, sh$celltype)])            # the current row is the main (gsd) result
  expect_equal(cur$estimate, sh$estimate[match(cur$celltype, sh$celltype)]); expect_equal(cur$p.distance, sh$p.distance[match(cur$celltype, sh$celltype)])
  # the alternative model: its distance magnitude, its gsd p (the distance test's own kept), as a gsd run under that model gives
  alt <- tb[tb$model == "plus nuisance", ]
  expect_true(all(is.finite(alt$p) & is.finite(alt$p.distance)))
  ra <- suppressWarnings(cao$estimateExpressionShiftMagnitudes(formula = ~ group + nuisance, test = "group", n.permutations = 99,   # the main result's count
                                                              seed = r$settings$gsd$seed, verbose = FALSE, name = "alt"))
  sa <- ra$results[ra$results$effect == "shift", ]
  expect_equal(alt$p, sa$p[match(alt$celltype, sa$celltype)]); expect_equal(alt$padj, sa$padj[match(alt$celltype, sa$celltype)])
  # two verdicts, both in the printout and the plot
  expect_true(all(c("verdict.magnitude", "verdict.significance", "verdict") %in% names(s$summary)))
  expect_true(all(grepl("^(stable|near zero|estimate changes|sign depends|driven by sample|one model)", s$summary$verdict.magnitude)))
  expect_true(all(grepl("^(significant in|not significant in|one model)", s$summary$verdict.significance)))
  expect_output(print(s), "significance: gsd")
  gg <- cao$plotSensitivity()
  expect_match(gg$labels$subtitle, "position: distance-model magnitude"); expect_match(gg$labels$subtitle, "filled: gsd adjusted p")
  expect_true(all(grepl("\n(significant in|not significant)", levels(gg$data$panel))))   # celltype, magnitude, significance: one line each
  # the distance model's own test on request (no gsd refits)
  s2 <- suppressWarnings(cao$checkSensitivity(covariate.sets = list(`plus nuisance` = ~ group + nuisance), n.permutations = 49,
                                              significance = "distance", verbose = FALSE, name = "s2"))
  expect_equal(s2$settings$significance, "distance test")
  cur2 <- s2$table[s2$table$effect == "shift" & s2$table$model == "current", ]
  expect_equal(cur2$p, sh$p.distance[match(cur2$celltype, sh$celltype)])
})

test_that("cell-type gsd, all relabelings enumerated: the smallest attainable p is that of the relabelings used, also with missing samples", {
  # 5 patients with a sample of each condition, 3 with one sample (partly paired); one cell type is missing from two samples
  cao <- makeToyCacoa(n.per.group = c(A = 8, B = 8), cells.per.sample = 60, n.genes = 200, shift = 0.6, seed = 9)
  meta <- cao$sample.meta; s <- rownames(meta); meta$patient <- NA
  a <- s[meta$group == "A"]; b <- s[meta$group == "B"]; meta[a, "patient"] <- paste0("p", 1:8); meta[b[1:5], "patient"] <- paste0("p", 1:5); meta[b[6:8], "patient"] <- paste0("q", 1:3)
  cao$sample.meta <- meta; cm <- Matrix::t(cao$data.object); lib <- Matrix::colSums(cm) * 2
  spc <- cao$sample.per.cell; drop <- names(spc)[as.character(spc) %in% c(a[1], b[2]) & cao$cell.groups[names(spc)] == "ct2"]
  keep <- setdiff(colnames(cm), drop); m <- buildCacoaModel(meta, ~ group + (1 | patient), test = "group")
  r <- suppressWarnings(gsdShiftsPerCellType(cm[, keep], spc[keep], cao$cell.groups[keep], m, lib[keep], n.permutations = 2000, seed = 1, influence = FALSE, n.null = 0))
  expect_true(all(r$results$p >= r$results$p.floor - 1e-12))
  expect_true(all(r$results$n.perm.distinct >= 1))
  # 2^5 pair flips x 20 arrangements of the 6 single samples = 640 relabelings; the remaining sample of a broken pair still takes the label
  # of its absent partner under a flip, so the cell type missing two samples has 640 distinct arrangements too (its own plan has 560)
  expect_equal(r$results$p.floor, rep(1 / 640, 2))
  expect_equal(r$results$n.perm.distinct, rep(640, 2))
})
