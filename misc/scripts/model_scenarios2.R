suppressPackageStartupMessages(library(cacoa))
sc <- function(name, expr) { r <- tryCatch(withCallingHandlers(expr, message = function(m) invokeRestart("muffleMessage"), warning = function(w) { cat(sprintf("  [warn] %s\n", conditionMessage(w))); invokeRestart("muffleWarning") }), error = function(e) paste("ERROR:", conditionMessage(e))); cat(sprintf("## %s\n   %s\n", name, paste(capture.output(print(r)), collapse = "\n   "))) }
set.seed(1); n <- 12
meta <- data.frame(group = rep(c("A", "B"), each = 6), batch = rep(c("b1", "b2"), 6), age = round(rnorm(n, 50, 6)), row.names = sprintf("s%02d", 1:n), stringsAsFactors = FALSE)
sc("A logical test variable: cause (works once converted to factor)", { mm <- meta; mm$treated <- rep(c(FALSE, TRUE), each = 6); mm$treated <- factor(mm$treated); colnames(buildCacoaModel(mm, test = "treated")$tests[[1]]$design$F) })
sc("B constant covariate in the formula is pruned silently", { mm <- meta; mm$site <- "x"; m <- buildCacoaModel(mm, ~ group + site + batch, test = "group"); list(formula = deparse(m$formula), used = deparse(m$tests[[1]]$design$formula_used), issues = m$issues$message) })
sc("C term test with interaction formula: term columns", { d <- buildCacoaModel(meta, ~ group * batch, test = list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b1"))); d2 <- buildCacoaModel(meta, ~ group * batch, test = "batch", contrast = NULL); NULL })
sc("C2 term test on a 3-level factor with interaction", { m3 <- data.frame(Group = rep(c("G1", "G2", "G3"), each = 4), Batch = rep(c("b1", "b2"), 6), row.names = sprintf("s%02d", 1:12)); buildCacoaModel(m3, ~ Group * Batch, test = "Group")$tests[[1]]$design$term.cols })
# D: kernel B Freedman-Lane with a multi-column X: the reduced model lacks the non-tested direction of the factor
sc("D kernel-B FL null rate, 3-level factor, G2 vs G1 null, strong G3 effect (should be ~0.05)", {
  set.seed(7); n <- 18
  m3 <- data.frame(Group = rep(c("G1", "G2", "G3"), each = 6), Batch = rep(c("b1", "b2"), 9), age = round(rnorm(n, 50, 6)), row.names = sprintf("s%02d", 1:n))
  m <- buildCacoaModel(m3, ~ Group + age, test = "Group: G2 vs G1"); d <- m$tests[[1]]$design
  Y <- matrix(rnorm(n * 400), n, 400) + 3 * (m3$Group == "G3") + 0.1 * m3$age
  r <- cacoa:::performLMPermutations(d, Y, perm.method = "freedman-lane", n.permutations = 199, seed = 1)
  dist.engine <- testPairwiseEffects(list(ct = as.matrix(dist(Y[, 1:30]))), d, m3, dist = "l2", n.permutations = 199, seed = 1, min.samp.per.level = 2)
  list(kernelB.FL.reject = mean(r$pval < 0.05), kernelB.block.reject = mean(cacoa:::performLMPermutations(d, Y, perm.method = "block", n.permutations = 199, seed = 1)$pval < 0.05),
       dist.engine.p.shift = dist.engine$results$p.shift, scheme = dist.engine$results$scheme, X.cols = colnames(d$X), Z.cols = colnames(d$Z))
})
sc("E same with the G3 effect removed (reference)", {
  set.seed(7); n <- 18
  m3 <- data.frame(Group = rep(c("G1", "G2", "G3"), each = 6), Batch = rep(c("b1", "b2"), 9), age = round(rnorm(n, 50, 6)), row.names = sprintf("s%02d", 1:n))
  m <- buildCacoaModel(m3, ~ Group + age, test = "Group: G2 vs G1"); d <- m$tests[[1]]$design
  Y <- matrix(rnorm(n * 400), n, 400) + 0.1 * m3$age
  mean(cacoa:::performLMPermutations(d, Y, perm.method = "freedman-lane", n.permutations = 199, seed = 1)$pval < 0.05)
})
sc("F buildDesignMatrices with structured contrast and no formula", buildDesignMatrices(meta, contrast = list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b2")))$formula_used)
sc("G checkFormula random-effect demotion", cacoa:::checkFormula(~ group + (1 | batch) + age))
sc("H numeric covariate interacting with the tested factor: what the anchor does", { d <- buildDesignMatrices(meta, contrast = list(type = "simple", term = "group", num = "B", den = "A"), formula = ~ group * age); list(cF = round(d$contrast.F, 3), anchors = d$numeric_ref_used) })
