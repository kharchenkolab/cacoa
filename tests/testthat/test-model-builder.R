# Model builder (misc/model_builder_plan.md). Step 0 freezes the behaviour that is right today and provides the
# oracle for the rewrite: hand-made treatment-coded designs that carry exactly the fields the engines read.

mbMeta <- function(n = 12, seed = 1) {
  set.seed(seed)
  data.frame(group = rep(c("A", "B"), each = n / 2), batch = rep(c("b1", "b2"), n / 2), age = round(rnorm(n, 50, 6)),
             row.names = sprintf("s%02d", seq_len(n)), stringsAsFactors = FALSE)
}
mbMeta3 <- function(n = 18, seed = 2) {
  set.seed(seed)
  data.frame(Group = rep(c("G1", "G2", "G3"), each = n / 3), Batch = rep(c("b1", "b2"), n / 2), age = round(rnorm(n, 50, 6)),
             row.names = sprintf("s%02d", seq_len(n)), stringsAsFactors = FALSE)
}
# squared-distance list for the engines: a few "cell types" with the same planted group effect
mbDist <- function(meta, var = "group", alt = "B", n.ct = 2, p = 25, effect = 1, seed = 3, drop = NULL) {
  set.seed(seed)
  lapply(seq_len(n.ct), function(k) {
    Y <- matrix(rnorm(nrow(meta) * p), nrow(meta), p, dimnames = list(rownames(meta), NULL))
    if (!is.null(alt)) Y[meta[[var]] == alt, ] <- Y[meta[[var]] == alt, ] + effect / sqrt(p)
    if (!is.null(drop) && k == 1) Y <- Y[setdiff(rownames(Y), drop), , drop = FALSE]
    as.matrix(dist(Y))
  }) |> setNames(paste0("ct", seq_len(n.ct)))
}
# hand-made design: treatment coding with `ref` as the baseline, c = row(alt) - row(ref); endpoints as the engines expect
handDesign <- function(meta, formula, var, alt, ref) {
  mm <- meta
  for (v in names(mm)) if (is.character(mm[[v]])) mm[[v]] <- factor(mm[[v]])
  numeric.test <- is.numeric(mm[[var]])
  if (!numeric.test) mm[[var]] <- stats::relevel(mm[[var]], ref = ref)
  F <- cacoa:::buildFullDesign(formula, mm)                    # model.matrix plus the terms / xlevels attributes the engines read
  row.at <- function(val) {                                   # design row at var = val, other factors at baseline, numerics at their mean
    nd <- mm
    for (v in names(nd)) nd[[v]] <- if (is.factor(nd[[v]])) factor(levels(nd[[v]])[1], levels = levels(nd[[v]])) else mean(nd[[v]])
    nd[[var]] <- if (numeric.test) val else factor(val, levels = levels(mm[[var]]))
    drop(stats::model.matrix(formula, nd)[1, colnames(F)])
  }
  cF <- row.at(alt) - row.at(ref); cF[abs(cF) < 1e-12] <- 0
  at <- function(val) list(list(at = setNames(list(val), var), w = 1))
  list(F = F, contrast.F = cF, X = F, Z = NULL, contrast.X = cF, meta = meta, formula_used = formula,
       contrast_spec = list(type = "simple", term = var, num = alt, den = ref, at = list()),
       contrast_endpoints_at = list(num = at(alt), den = at(ref)), numeric_ref_used = list())
}
lmContrast <- function(y, meta, formula, cF) {                 # lm estimate, se and t of a contrast on the hand-coded design
  mm <- meta; for (v in names(mm)) if (is.character(mm[[v]])) mm[[v]] <- factor(mm[[v]])
  fit <- stats::lm(stats::as.formula(paste("y", paste(deparse(formula), collapse = ""))), data = cbind(y = y, mm))
  b <- stats::coef(fit); V <- stats::vcov(fit); cF <- cF[names(b)]
  est <- sum(cF * b); se <- sqrt(drop(t(cF) %*% V %*% cF))
  c(effect = est, se = se, t = est / se, df = fit$df.residual)
}

test_that("grammar: the frozen rows resolve as specified", {
  meta <- mbMeta()
  m <- buildCacoaModel(meta, ~ group + batch, test = "group"); t <- m$tests[[1]]
  expect_equal(t$kind, "contrast"); expect_equal(unname(t$levels), c("A", "B")); expect_equal(t$reference.reason, "most-frequent")
  mm <- meta; mm$cond <- rep(c("disease", "control"), each = 6)
  t <- buildCacoaModel(mm, test = "cond")$tests[[1]]
  expect_equal(unname(t$levels[["ref"]]), "control"); expect_equal(t$reference.reason, "control-like"); expect_equal(t$label, "cond: disease vs control")
  mm$cond <- factor(mm$cond, levels = c("disease", "control"))
  t <- buildCacoaModel(mm, test = "cond")$tests[[1]]
  expect_equal(unname(t$levels[["ref"]]), "disease"); expect_equal(t$reference.reason, "explicit-order")
  mm$g <- c(rep("X", 8), rep("Y", 4))
  expect_equal(buildCacoaModel(mm, test = "g")$tests[[1]]$label, "g: Y vs X")
  t <- buildCacoaModel(meta, ~ age + batch, test = "age")$tests[[1]]
  expect_equal(t$step, 1); expect_equal(t$label, "age: per 1 unit")
  t <- buildCacoaModel(meta, ~ age + batch, test = c("age", "60", "50"))$tests[[1]]
  expect_equal(t$step, 10); expect_equal(unname(t$contrast$den), 50)
  t <- buildCacoaModel(meta, ~ group + batch, test = "group: A vs B")$tests[[1]]
  expect_equal(unname(t$levels), c("B", "A")); expect_equal(t$reference.reason, "user")
  labs <- vapply(buildCacoaModel(meta, ~ group + age + batch, test = "all")$tests, `[[`, character(1), "label")
  expect_equal(labs, c("group: B vs A", "age: per 1 unit", "batch: b2 vs b1"))
  expect_error(buildCacoaModel(meta, ~ group + sex, test = "group"), "not found in sample metadata")
  expect_error(buildCacoaModel(meta, test = "group: C vs A"), "not found")
  mm$k <- "x"; expect_error(buildCacoaModel(mm, test = "k"), "fewer than two")
  mm$sample_id <- rownames(mm); expect_error(buildCacoaModel(mm, ~ group + sample_id, test = "group"), "unique value per sample")
})

test_that("contrast vectors equal the endpoint difference of a hand-coded design and give the lm contrast", {
  meta <- mbMeta(); meta3 <- mbMeta3()
  set.seed(4); y <- rnorm(12); y3 <- rnorm(18)
  cases <- list(
    list(meta = meta, f = ~ group + batch, test = "group", var = "group", alt = "B", ref = "A"),
    list(meta = meta3, f = ~ Group + age, test = "Group: G2 vs G1", var = "Group", alt = "G2", ref = "G1"),
    list(meta = meta3, f = ~ Group + Batch, test = "Group: G3 vs G1", var = "Group", alt = "G3", ref = "G1"),
    list(meta = meta, f = ~ age + batch, test = "age", var = "age", alt = mean(meta$age) + 1, ref = mean(meta$age)),
    list(meta = meta, f = ~ log(age) + batch, test = "age", var = "age", alt = mean(meta$age) + 1, ref = mean(meta$age)),
    list(meta = meta, f = ~ group * batch, test = list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b2")), var = "group", alt = "B", ref = "A", at = list(batch = "b2")),
    list(meta = meta, f = ~ group * batch, test = list(type = "marginal", term = "group", num = "B", den = "A", over = "batch"), var = "group", alt = "B", ref = "A", marginal = "batch"))
  for (cs in cases) {
    m <- buildCacoaModel(cs$meta, cs$f, test = cs$test); d <- m$tests[[1]]$design
    yy <- if (nrow(cs$meta) == 12) y else y3
    hd <- handDesign(cs$meta, cs$f, cs$var, cs$alt, cs$ref)
    cF <- hd$contrast.F
    if (!is.null(cs$at)) {                                          # evaluate the hand contrast at batch = b2 instead of the baseline
      mm <- cs$meta; mm$group <- factor(mm$group); mm$batch <- factor(mm$batch)
      row <- function(g) drop(stats::model.matrix(cs$f, data.frame(group = factor(g, levels = c("A", "B")), batch = factor("b2", levels = c("b1", "b2"))))[1, ])
      cF <- row(cs$alt) - row(cs$ref)
    }
    if (!is.null(cs$marginal)) {
      row <- function(g, b) drop(stats::model.matrix(cs$f, data.frame(group = factor(g, levels = c("A", "B")), batch = factor(b, levels = c("b1", "b2"))))[1, ])
      cF <- 0.5 * (row("B", "b1") - row("A", "b1")) + 0.5 * (row("B", "b2") - row("A", "b2"))
    }
    ref <- lmContrast(yy, cs$meta, cs$f, cF)
    # the builder's contrast gives the same lm contrast whatever its coding (regression weights are coding-free)
    expect_equal(unname(drop(regressionWeights(d) %*% yy)), unname(ref[["effect"]]), tolerance = 1e-9, info = m$tests[[1]]$label)
  }
})

test_that("per-column fitter agrees with lm on the builder's designs (both schemes)", {
  meta <- mbMeta(); meta3 <- mbMeta3(); set.seed(5)
  Y <- matrix(rnorm(12 * 3), 12, 3); Y3 <- matrix(rnorm(18 * 3), 18, 3)
  chk <- function(m, Y, var, alt, ref, schemes = c("block", "freedman-lane")) {
    d <- m$tests[[1]]$design; hd <- handDesign(m$meta, m$formula, var, alt, ref)
    for (sch in schemes) {
      r <- cacoa:::performLMPermutations(d, Y, perm.method = sch, n.permutations = 19, seed = 1)
      for (j in seq_len(ncol(Y))) {
        ref.j <- lmContrast(Y[, j], m$meta, m$formula, hd$contrast.F)
        expect_equal(unname(c(r$effect[j], r$se[j], r$stat.obs[j])), unname(ref.j[c("effect", "se", "t")]), tolerance = 1e-9, info = paste(m$tests[[1]]$label, sch))
        expect_equal(unname(r$df[j]), unname(ref.j[["df"]]), info = paste(m$tests[[1]]$label, sch))
      }
    }
  }
  chk(buildCacoaModel(meta, ~ group + batch, test = "group"), Y, "group", "B", "A")
  chk(buildCacoaModel(meta3, ~ Group + age, test = "Group: G2 vs G1"), Y3, "Group", "G2", "G1")
  chk(buildCacoaModel(meta, ~ age + batch, test = "age"), Y, "age", mean(meta$age) + 1, mean(meta$age))
  chk(buildCacoaModel(meta, ~ log(age) + batch, test = "age"), Y, "age", mean(meta$age) + 1, mean(meta$age))
})

test_that("engine invariance: the builder's design and a hand-coded design give identical shift results and fitter output", {
  meta <- mbMeta(); meta3 <- mbMeta3()
  runs <- list(
    list(meta = meta, f = ~ group + batch, test = "group", var = "group", alt = "B", ref = "A", D = mbDist(meta)),
    list(meta = meta3, f = ~ Group + age, test = "Group: G2 vs G1", var = "Group", alt = "G2", ref = "G1", D = mbDist(meta3, "Group", "G2")),
    list(meta = meta, f = ~ age + batch, test = "age", var = "age", alt = mean(meta$age) + 1, ref = mean(meta$age), D = mbDist(meta, alt = NULL)),
    list(meta = meta, f = ~ group + batch, test = "group", var = "group", alt = "B", ref = "A", D = mbDist(meta, drop = c("s01", "s07")), robust = "huber", na.mode = "impute_weak"))
  for (r in runs) {
    m <- buildCacoaModel(r$meta, r$f, test = r$test); d <- m$tests[[1]]$design
    hd <- handDesign(r$meta, r$f, r$var, r$alt, r$ref)
    rob <- if (is.null(r$robust)) "none" else r$robust; nam <- if (is.null(r$na.mode)) "drop" else r$na.mode
    a <- testPairwiseEffects(r$D, d, r$meta, dist = "l2", n.permutations = 39, seed = 7, min.samp.per.level = 2, robust = rob, na.mode = nam)
    b <- testPairwiseEffects(r$D, hd, r$meta, dist = "l2", n.permutations = 39, seed = 7, min.samp.per.level = 2, robust = rob, na.mode = nam)
    cols <- c("celltype", "shift", "var", "total", "F", "p.shift", "p.var", "p.total", "n")
    expect_equal(a$results[, cols], b$results[, cols], tolerance = 1e-8, info = m$tests[[1]]$label)
    expect_equal(a$global$p, b$global$p, info = m$tests[[1]]$label)
    set.seed(8); Y <- matrix(rnorm(nrow(r$meta) * 3), nrow(r$meta), 3)
    fa <- cacoa:::performLMPermutations(d, Y, perm.method = "block", n.permutations = 39, seed = 7, robust.method = rob, na.mode = nam)
    fb <- cacoa:::performLMPermutations(hd, Y, perm.method = "block", n.permutations = 39, seed = 7, robust.method = rob, na.mode = nam)
    expect_equal(unname(fa$effect), unname(fb$effect), tolerance = 1e-9); expect_equal(unname(fa$se), unname(fb$se), tolerance = 1e-9)
    expect_equal(unname(fa$pval), unname(fb$pval)); expect_equal(unname(fa$P), unname(fb$P))
  }
})

test_that("term test: location F equals the partial ANOVA F for a one-dimensional response", {
  meta3 <- mbMeta3(); set.seed(9); y <- rnorm(18) + 0.8 * (meta3$Group == "G2")
  m <- buildCacoaModel(meta3, ~ Group + Batch, test = "Group"); d <- m$tests[[1]]$design
  expect_equal(m$tests[[1]]$kind, "term")
  D <- as.matrix(dist(y)); rownames(D) <- colnames(D) <- rownames(meta3)
  r <- testTermEffects(list(ct = D), d, meta3, dist = "l2", n.permutations = 19, seed = 1)
  mm <- meta3; mm$Group <- factor(mm$Group); mm$Batch <- factor(mm$Batch)
  Fref <- anova(lm(y ~ Batch, mm), lm(y ~ Batch + Group, mm))$F[2]
  expect_equal(r$results$F[1], Fref, tolerance = 1e-8)
})

test_that("listwise deletion and design issues are recorded", {
  meta <- mbMeta(); meta$age[3] <- NA
  m <- buildCacoaModel(meta, ~ group + age, test = "group")
  expect_equal(m$samples$dropped, "s03"); expect_equal(length(m$samples$used), 11)
  expect_true(any(grepl("dropped for missing", m$issues$message)))
  meta <- mbMeta(); meta$group[3] <- NA
  expect_equal(buildCacoaModel(meta, ~ group + batch, test = "group")$samples$dropped, "s03")
  mm <- mbMeta(); mm$copy <- mm$group                                      # a covariate that determines the test variable
  m <- buildCacoaModel(mm, ~ group + copy, test = "group")
  expect_true(any(m$issues$severity == "error")); expect_match(m$issues$message[m$issues$severity == "error"], "fully determined")
})
