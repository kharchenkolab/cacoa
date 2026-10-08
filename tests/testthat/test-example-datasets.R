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
