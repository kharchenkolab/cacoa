# Tests on the example datasets found in the development folder (skipped when the files are absent).
#  - SCC: per-cell-type pseudobulk profiles of 18 paired tumour / normal skin samples (examples/scc/x.rds)
#  - simulated Cacoa objects with known per-cell-type DE strength (test/shifts_sim_objects.rds; slow)
#  - the small conos panel shipped with the package (data/panel.preprocessed)

devRoot <- function() {
  cands <- c(Sys.getenv("CACOA_DEV_ROOT"), "../../..", "../..", normalizePath("~/cacoa-dev", mustWork = FALSE))
  for (p in cands) if (nzchar(p) && dir.exists(file.path(p, "examples"))) return(normalizePath(p))
  NULL
}

test_that("SCC tumour vs normal: paired design, block permutations within patient", {
  root <- devRoot()
  skip_if(is.null(root) || !file.exists(file.path(root, "examples/scc/x.rds")), "SCC example object not available")
  x <- readRDS(file.path(root, "examples/scc/x.rds"))
  samples <- rownames(x$sample.model$F)
  meta <- data.frame(condition = factor(sub(".*_", "", samples), levels = c("Normal", "Tumor")),
                     patient = factor(sub("_.*", "", samples)), row.names = samples)
  expect_equal(nlevels(meta$patient), 10)
  D.list <- sampleDistanceMatrices(x$cm.per.type, dist = "cor")
  D.list <- Filter(function(D) nrow(D) >= 6, D.list)
  expect_true(length(D.list) >= 10)
  # unpaired model
  m1 <- buildCacoaModel(meta, test = "condition")
  expect_equal(m1$tests[[1]]$contrast, c("condition", "Tumor", "Normal"))     # "Normal" is control-like
  r1 <- expressionShiftsForModel(D.list, m1, dist = "cor", n.permutations = 199, seed = 1, influence = TRUE)
  sh1 <- r1$results[r1$results$effect == "shift", ]
  expect_true(all(sh1$p > 0 & sh1$p <= 1))
  expect_true(sum(sh1$padj < 0.05) >= 3)                                       # tumour vs normal is a strong contrast
  expect_true(r1$global$p[r1$global$effect == "shift"] < 0.05)
  # paired model: patient as a blocking covariate -> block permutations within patient
  m2 <- buildCacoaModel(meta, formula = ~ condition + patient, test = "condition")
  expect_equal(m2$tests[[1]]$permutation$scheme, "block")
  expect_equal(m2$tests[[1]]$permutation$n.strata, 10)
  expect_equal(round(m2$tests[[1]]$permutation$n.distinct), 2^8)             # 8 patients with both samples
  r2 <- expressionShiftsForModel(D.list, m2, dist = "cor", n.permutations = 299, seed = 1)   # 299 > 256: enumerated exactly
  sh2 <- r2$results[r2$results$effect == "shift", ]
  expect_true(all(sh2$n.perm.distinct <= 256))
  expect_true(all(sh2$p >= 1 / 256 - 1e-12))
  expect_true(all(r2$wide[[1]]$exhaustive))
  # shifts agree in ranking between the two models (same cell types on top)
  common <- intersect(sh1$celltype, sh2$celltype)
  expect_gt(cor(sh1$estimate[match(common, sh1$celltype)], sh2$estimate[match(common, sh2$celltype)], method = "spearman"), 0.6)
  # adjusted distances are proper distance matrices
  A <- r2$adjusted.distances[[1]][[1]]
  expect_true(isSymmetric(unname(A))); expect_true(all(A >= -1e-10))
})

test_that("conos panel shipped with the package runs through the new entry point", {
  skip_if_not_installed("conos")
  data("panel.preprocessed", package = "cacoa")
  con <- conos::Conos$new(panel.preprocessed, n.cores = 1)
  suppressWarnings(con$buildGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 500, verbose = FALSE))
  meta <- data.frame(condition = factor(c("CTRL", "CTRL", "DISEASE", "DISEASE")), row.names = names(con$samples))
  cell.groups <- con$getDatasetPerCell() %>% as.factor()                      # 4 "cell types" = samples (smoke only)
  set.seed(1); cell.groups <- factor(sample(c("t1", "t2"), length(cell.groups), replace = TRUE)) %>% setNames(names(cell.groups))
  cao <- suppressMessages(Cacoa$new(con, sample.metadata = meta, test = "condition", cell.groups = cell.groups, n.cores = 1, verbose = FALSE))
  res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 19, min.samp.per.level = 2, verbose = FALSE)
  expect_true(nrow(res$results) > 0)
  expect_true(all(res$results$exhaustive))                                    # 2 vs 2: 6 distinct relabelings
  expect_s3_class(cao$plotExpressionShiftMagnitudes(), "ggplot")
  # cluster-free shifts through the object (graph from conos)
  cf <- cao$estimateClusterFreeExpressionShifts(n.top.genes = 300, n.permutations = 19, min.samp.per.level = 2, verbose = FALSE)
  expect_equal(length(cf$stat), length(cao$cell.groups))
  expect_true(sum(is.finite(cf$z.adj)) > 0)
  expect_s3_class(cf$model, "cacoaModel")
})

test_that("simulated objects: planted DE strength orders the shifts (slow)", {
  skip_on_cran()
  skip_if(Sys.getenv("CACOA_SKIP_SLOW") == "true", "slow test skipped")
  root <- devRoot()
  f <- if (!is.null(root)) file.path(root, "test/shifts_sim_objects.rds") else ""
  skip_if(!file.exists(f), "simulated objects not available")
  skip_if_not_installed("conos")
  strength <- c("AST-FB" = "weak", "OPC" = "moderate", "IN-SST" = "strong", "Microglia" = "none", "L5/6-CC" = "moderate",
                "IN-PV" = "strong", "Neu-NRGN" = "none", "AST-PP" = "moderate", "L2/3" = "strong")
  sims <- readRDS(f)
  for (nm in c("no_batch", "with_batch")) {
    old <- sims[[nm]]
    meta <- old$data.object$misc$sample_meta
    expect_true(all(c("Group", "Batch") %in% names(meta)))
    if (is.null(rownames(meta)) || !all(names(old$data.object$samples) %in% rownames(meta))) rownames(meta) <- meta$Individual
    cao <- suppressMessages(Cacoa$new(old$data.object, sample.metadata = meta, formula = ~ Group + Batch,
                                      test = "Group: Group2 vs Group1", cell.groups = old$cell.groups,
                                      sample.per.cell = old$sample.per.cell, n.cores = 4, verbose = FALSE))
    expect_equal(cao$model$tests[[1]]$permutation$scheme, "block")
    res <- cao$estimateExpressionShiftMagnitudes(n.permutations = 199, verbose = FALSE)
    sh <- res$results[res$results$effect == "shift", ]
    s <- strength[sh$celltype]
    # every strong cell type is significant; none-types are not; strong > none in estimate
    expect_true(all(sh$padj[s == "strong"] < 0.05), info = nm)
    expect_true(all(sh$p[s == "none"] > 0.05), info = nm)
    expect_gt(min(sh$estimate[s == "strong"]), max(sh$estimate[s == "none"]))
    expect_gt(mean(sh$estimate[s == "moderate"]), mean(sh$estimate[s == "none"]))
  }
})

test_that("simulated objects: the whole workflow runs end to end on a conos-backed object (slow)", {
  skip_on_cran()
  skip_if(Sys.getenv("CACOA_SKIP_SLOW") == "true", "slow test skipped")
  root <- devRoot()
  f <- if (!is.null(root)) file.path(root, "test/shifts_sim_objects.rds") else ""
  skip_if(!file.exists(f), "simulated objects not available")
  skip_if_not_installed("conos"); skip_if_not_installed("limma"); skip_if_not_installed("coda.base")
  old <- readRDS(f)$with_batch
  meta <- old$data.object$misc$sample_meta
  if (is.null(rownames(meta)) || !all(names(old$data.object$samples) %in% rownames(meta))) rownames(meta) <- meta$Individual
  set.seed(1); meta$age <- round(rnorm(nrow(meta), 55, 8))                    # a numeric covariate unrelated to the groups
  # 1. construct without a model: metadata audit
  expect_message(cao <- Cacoa$new(old$data.object, sample.metadata = meta, cell.groups = old$cell.groups, sample.per.cell = old$sample.per.cell,
                                  n.cores = 4, verbose = TRUE), "Sample metadata")
  expect_null(cao$model)
  cao$setOptions(n.permutations = 199, seed = 3, verbose = FALSE)
  # 2. design check and screen before any model
  chk <- cao$checkDesign(test = "Group")
  expect_s3_class(chk, "cacoaDesignCheck"); expect_false(any(chk$issues$severity == "error"))
  sc <- cao$screenCovariates(n.permutations = 99)
  expect_true("Group" %in% sc$covariates)
  g <- sc$global[sc$global$mode == "partial", ]
  expect_lt(g$p.global[g$covariate == "Group"], 0.05)                          # the planted groups are found
  expect_gt(g$p.global[g$covariate == "age"], 0.05)                            # the unrelated covariate is not
  expect_s3_class(cao$plotCovariateScreen(), "ggplot")
  # 3. model: two-level contrast adjusted for batch; whole-factor test of Group as a second test
  cao$setModel(~ Group + Batch, test = c("Group: Group2 vs Group1", "Group: all"))
  expect_length(cao$model$tests, 2); expect_equal(cao$model$tests[[2]]$kind, "term")
  expect_equal(cao$ref.level, "Group1")
  res <- cao$estimateExpressionShiftMagnitudes()
  expect_setequal(unique(res$results$effect), c("shift", "var", "total", "location", "dispersion"))
  sh <- res$results[res$results$effect == "shift", ]
  expect_true(all(sh$padj[sh$celltype %in% c("IN-SST", "IN-PV", "L2/3")] < 0.05))
  loc <- res$results[res$results$effect == "location", ]
  expect_true(all(loc$padj[loc$celltype %in% c("IN-SST", "IN-PV", "L2/3")] < 0.05))
  expect_s3_class(cao$plotExpressionShiftMagnitudes(), "ggplot")
  expect_true(inherits(cao$plotShiftDetail("IN-SST"), "ggplot"))
  expect_s3_class(cao$plotSampleInfluence(), "ggplot")
  expect_s3_class(cao$plotSampleDistances(space = "expression.shifts", cell.type = "IN-SST", values = "unadjusted", color.by = "test"), "ggplot")
  # 4. sensitivity: planted effects are robust
  sens <- cao$checkSensitivity(n.permutations = 99)
  expect_true(all(sens$summary$verdict[sens$summary$celltype %in% c("IN-SST", "IN-PV", "L2/3")] == "robust"))
  expect_s3_class(cao$plotSensitivity(), "ggplot")
  # 5. other analyses on the same model
  de <- cao$estimateDEPerCellType(test = "limma-voom", min.cell.count = 5)
  expect_s3_class(attr(de, "model"), "cacoaModel")
  expect_gt(sum(de[["IN-SST"]]$res$padj < 0.05, na.rm = TRUE), sum(de[["Microglia"]]$res$padj < 0.05, na.rm = TRUE))
  coda <- cao$estimateCellLoadings(n.permutations = 199)
  expect_s3_class(coda$model, "cacoaModel")
  dens <- cao$estimateCellDensity(method = "graph")
  expect_equal(sum(dens$sample.weights), 0, tolerance = 1e-10)
  cf <- cao$estimateClusterFreeExpressionShifts(n.top.genes = 500, n.permutations = 49)
  z <- tapply(cf$z.adj, old$cell.groups[names(cf$z.adj)], median, na.rm = TRUE); sh <- tapply(cf$shifts, old$cell.groups[names(cf$shifts)], median, na.rm = TRUE)
  expect_gte(z[["IN-SST"]], z[["Microglia"]]); expect_gt(sh[["IN-SST"]], sh[["Microglia"]])   # both adjusted medians may sit at the z cap and the pooled p at its floor
  expect_s3_class(cao$plotClusterFreeExpressionShifts(), "ggplot")                       # shift, pooled z, adjusted z
  ggs <- cao$plotClusterFreeExpressionShifts(build.panel = FALSE)
  expect_length(ggs, 3); expect_length(cao$plotClusterFreeExpressionShifts(adjusted = TRUE, build.panel = FALSE), 2)
  expect_equal(unname(cf$z.pooled[!is.na(cf$z.pooled)]), unname(qnorm(cf$p.pooled[!is.na(cf$z.pooled)], lower.tail = FALSE)), tolerance = 1e-9)
})

test_that("planted density simulation: the BH adjustment calls the changed cell types and max-t does not", {
  skip_if(isTRUE(as.logical(Sys.getenv("CACOA_SKIP_SLOW", "true"))), "slow")
  f <- "../../../test/sim_density.rds"; skip_if_not(file.exists(f), "simulated density dataset not available")
  d <- readRDS(f)
  cao <- Cacoa$new(d$con, sample.metadata = d$sample.meta, cell.groups = d$cell.groups, sample.per.cell = d$sample.per.cell, n.cores = 8, verbose = FALSE)
  cao$setOptions(seed = 1); cao$setModel(~ treatment, test = "treatment", verbose = FALSE)
  cao$estimateCellDensity(method = "graph", verbose = FALSE)
  ct <- cao$cell.groups
  frac <- function(z) tapply(abs(z[names(ct)]) > qnorm(0.975), ct, mean, na.rm = TRUE)
  bh <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 199, verbose = FALSE)$diff$permutation
  expect_equal(bh$adjust.method, "perm-fdr"); expect_false(is.null(bh$group.test))
  gs <- cao$estimateClusterFreeSummary("cell.density")
  expect_true(all(gs$significant[gs$group %in% c("IN-SST", "L2/3", "L5/6-CC", "IN-PV", "OPC")]))
  expect_false(any(gs$significant[gs$group %in% c("AST-PP", "Microglia", "Neu-NRGN")]))
  expect_equal(gs$direction[gs$group == "IN-SST"], "down"); expect_equal(gs$direction[gs$group == "OPC"], "up")
  f.bh <- frac(bh$adj)
  # changed between the arms (planted or inherited from the simulation): IN-SST, L2/3, L5/6-CC down; IN-PV, OPC up
  # (kNN unit of 100 cells with the averaged statistic: IN-PV is called in about three quarters of its cells)
  for (ct.i in c("IN-SST", "L2/3", "L5/6-CC", "IN-PV", "OPC")) expect_gt(f.bh[[ct.i]], 0.7)
  expect_lt(median(bh$adj[names(ct)[ct == "IN-SST"]], na.rm = TRUE), 0); expect_gt(median(bh$adj[names(ct)[ct == "OPC"]], na.rm = TRUE), 0)
  expect_lt(max(f.bh[c("AST-PP", "Microglia", "Neu-NRGN")]), 0.1)
  mt <- cao$estimateDiffCellDensity(type = "permutation", n.permutations = 199, adjust.method = "max-t", verbose = FALSE)$diff$permutation
  expect_lt(sum(frac(mt$adj) > 0.8), sum(f.bh > 0.8))                                   # the familywise adjustment finds fewer
  s <- cao$estimateClusterFreeExpressionShifts(n.top.genes = 300, n.permutations = 99, verbose = FALSE)
  expect_equal(s$settings$adjust.method, "perm-fdr"); expect_false(is.null(s$group.test))
  ss <- cao$estimateClusterFreeSummary("cluster.free.expr.shifts"); expect_true(all(ss$significant[ss$group %in% c("L2/3", "IN-PV", "IN-SST")]))
  fs <- tapply(s$z.adj[names(ct)] > qnorm(0.95), ct, mean, na.rm = TRUE)
  # on this graph every cell has 340+ neighbours, so the kNN unit is the whole cluster and the averaged statistic gives a
  # near per-type verdict; Neu-NRGN (nothing planted, but a fifth of its unit cells are types whose abundance changes
  # between the arms, a composition-driven neighbourhood shift) is the least called type; without the averaging it is below 0.2
  expect_gt(fs[["L2/3"]], 0.5); expect_gt(fs[["IN-PV"]], 0.5); expect_equal(names(which.min(fs)), "Neu-NRGN"); expect_lt(fs[["Neu-NRGN"]], min(fs[names(fs) != "Neu-NRGN"]))
  s0 <- cao$estimateClusterFreeExpressionShifts(n.top.genes = 300, n.permutations = 99, post = "none", verbose = FALSE)
  expect_lt(mean(s0$z.adj[names(ct)[ct == "Neu-NRGN"]] > qnorm(0.95), na.rm = TRUE), 0.2)
  expect_true(all(s$z.adj >= 0, na.rm = TRUE))
})
