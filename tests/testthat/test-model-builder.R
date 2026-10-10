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
  m <- buildCacoaModel(meta3, ~ Group + Batch, test = "Group: all"); d <- m$tests[[1]]$design
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

test_that("coding: the tested factor is coded as level means (reference first) and names are stable across tests", {
  meta <- mbMeta(); meta3 <- mbMeta3()
  d <- buildCacoaModel(meta, ~ group + batch, test = "group")$tests[[1]]$design
  expect_equal(colnames(d$F), c("groupA", "groupB", "batchb2"))
  expect_equal(unname(d$contrast.F), c(-1, 1, 0)); expect_equal(colnames(d$X), "contrast"); expect_equal(ncol(d$Z), 2)
  set.seed(1); y <- rnorm(12); mm <- meta; mm$group <- factor(mm$group); mm$batch <- factor(mm$batch)
  expect_equal(unname(qr.coef(qr(d$F), y)), unname(coef(lm(y ~ 0 + group + batch, mm))), tolerance = 1e-10)
  b <- qr.coef(qr(cbind(d$X, d$Z)), y)
  expect_equal(unname(b[1]), unname(coef(lm(y ~ group + batch, mm))["groupB"]), tolerance = 1e-10)
  # three-level factor: the same columns for every test of one model; level differences for the term test
  m <- buildCacoaModel(meta3, ~ Group + Batch, test = c("Group: G2 vs G1", "Group: G3 vs G1", "Group: all"))
  cn <- lapply(m$tests, function(t) colnames(t$design$F))
  expect_true(all(vapply(cn, identical, logical(1), cn[[1]]))); expect_equal(cn[[1]], c("GroupG1", "GroupG2", "GroupG3", "Batchb2"))
  expect_equal(unname(m$tests[[1]]$design$contrast.F), c(-1, 1, 0, 0)); expect_equal(unname(m$tests[[2]]$design$contrast.F), c(-1, 0, 1, 0))
  expect_equal(colnames(m$tests[[3]]$design$term.contrast), c("G2 vs G1", "G3 vs G1"))
  expect_equal(unname(m$tests[[3]]$design$term.contrast[, 1]), c(-1, 1, 0, 0))
  # the reference level comes first even when it is not the first level alphabetically
  mm2 <- meta; mm2$cond <- rep(c("disease", "control"), each = 6)
  expect_equal(colnames(buildCacoaModel(mm2, test = "cond")$tests[[1]]$design$F), c("condcontrol", "conddisease"))
  # logical test variable
  mm2$treated <- rep(c(FALSE, TRUE), each = 6)
  dl <- buildCacoaModel(mm2, test = "treated")$tests[[1]]$design
  expect_equal(colnames(dl$F), c("treatedFALSE", "treatedTRUE")); expect_equal(unname(dl$contrast.F), c(-1, 1))
  # numeric test variable keeps the intercept; interaction model stays full rank; explicit no-intercept formula; factor() wrapper
  expect_equal(colnames(buildCacoaModel(meta, ~ age + batch, test = "age")$tests[[1]]$design$F), c("(Intercept)", "age", "batchb2"))
  di <- buildCacoaModel(meta, ~ group * batch, test = list(type = "marginal", term = "group", num = "B", den = "A", over = "batch"))$tests[[1]]$design
  expect_equal(colnames(di$F), c("groupA", "groupB", "batchb2", "groupB:batchb2")); expect_equal(qr(di$F)$rank, 4)
  expect_equal(colnames(buildCacoaModel(meta, ~ 0 + batch + group, test = "group")$tests[[1]]$design$F), c("groupA", "groupB", "batchb2"))
  mm3 <- meta; mm3$batch <- rep(1:2, 6)
  expect_equal(colnames(buildCacoaModel(mm3, ~ group + factor(batch), test = "group")$tests[[1]]$design$F), c("groupA", "groupB", "factor(batch)2"))
  # the default formula of the expert constructor gives the same coding
  expect_equal(colnames(suppressMessages(buildDesignMatrices(meta, contrast = c("group", "B", "A")))$F), c("groupA", "groupB"))
})

test_that("grammar: a multi-level factor must name its comparison; 'var: all' asks for the whole-factor test", {
  meta3 <- mbMeta3()
  expect_error(buildCacoaModel(meta3, ~ Group + Batch, test = "Group"), "has 3 levels \\(G1, G2, G3\\)")
  expect_error(buildCacoaModel(meta3, ~ Group + Batch, test = "Group"), 'test = "Group: G2 vs G1"')
  expect_error(buildCacoaModel(meta3, ~ Group + Batch, test = "Group"), 'test = "Group: all"')
  m <- buildCacoaModel(meta3, ~ Group + Batch, test = "Group: all")
  expect_equal(m$tests[[1]]$kind, "term"); expect_equal(m$tests[[1]]$label, "Group (3 levels)")
  expect_equal(buildCacoaModel(meta3, ~ Group + Batch, test = "Group: any")$tests[[1]]$kind, "term")
  expect_error(buildCacoaModel(meta3, ~ Group + age, test = "age: all"), "is numeric")
  # test = "all" takes multi-level factors as whole-factor tests
  labs <- vapply(buildCacoaModel(meta3, ~ Group + Batch + age, test = "all")$tests, `[[`, character(1), "label")
  expect_equal(labs, c("Group (3 levels)", "Batch: b2 vs b1", "age: per 1 unit"))
})

test_that("interaction models: the plain grammar gets a documented default and a note", {
  meta <- mbMeta()
  m <- buildCacoaModel(meta, ~ group * batch, test = "group"); d <- m$tests[[1]]$design
  ms <- buildCacoaModel(meta, ~ group * batch, test = list(type = "marginal", term = "group", num = "B", den = "A", over = "batch"))
  expect_equal(d$contrast.F, ms$tests[[1]]$design$contrast.F)
  expect_true(any(grepl("compared marginally", m$issues$message[m$issues$severity == "note"])))
  expect_output(print(m), "note: group interacts with batch")
  # factor x numeric: the comparison at the covariate's mean equals the structured simple contrast (anchors at the mean)
  m2 <- buildCacoaModel(meta, ~ group * age, test = "group"); d2 <- m2$tests[[1]]$design
  s2 <- buildCacoaModel(meta, ~ group * age, test = list(type = "simple", term = "group", num = "B", den = "A"))$tests[[1]]$design
  expect_equal(d2$contrast.F, s2$contrast.F); expect_true(any(grepl("compared at age = ", m2$issues$message[m2$issues$severity == "note"])))
  # numeric x factor: the slope averaged over the factor's levels
  m3 <- buildCacoaModel(meta, ~ age * group, test = "age"); d3 <- m3$tests[[1]]$design
  mm <- meta; mm$group <- factor(mm$group)
  set.seed(2); y <- rnorm(12); fit <- lm(y ~ age * group, mm)
  slope.avg <- unname(coef(fit)["age"] + 0.5 * coef(fit)["age:groupB"])
  expect_equal(unname(drop(regressionWeights(d3) %*% y)), slope.avg, tolerance = 1e-9)
  # explicit at = / over = override the default and leave no note
  m4 <- buildCacoaModel(meta, ~ group * batch, test = list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b2")))
  expect_false(any(grepl("interacts", m4$issues$message)))
  # a term test under an interaction model is the main effect at the reference setting
  m5 <- buildCacoaModel(mbMeta3(), ~ Group * Batch, test = "Group: all"); C <- m5$tests[[1]]$design$term.contrast
  expect_equal(rownames(C)[1:3], c("GroupG1", "GroupG2", "GroupG3")); expect_true(all(C[grepl(":", rownames(C)), ] == 0))
})

test_that("formula hygiene: dropped constant terms are notes, random effects and lincomb are errors", {
  meta <- mbMeta(); meta$site <- "x"
  m <- buildCacoaModel(meta, ~ group + site + batch, test = "group")
  expect_true(any(grepl("site dropped", m$issues$message[m$issues$severity == "note"])))
  expect_equal(deparse(m$tests[[1]]$design$formula_used), "~group + batch")
  expect_output(print(m), "note: term site dropped")
  expect_error(buildCacoaModel(meta, ~ group + (1 | batch), test = "group"), "random-effect terms")
  expect_error(buildCacoaModel(meta, ~ group + batch, test = list(type = "lincomb", term = "group", cells = c(A = 1))), "no longer supported")
})

test_that("model printout snapshots", {
  meta <- mbMeta(); meta3 <- mbMeta3()
  expect_snapshot(print(buildCacoaModel(meta, ~ group + batch, test = "group")))
  expect_snapshot(print(buildCacoaModel(meta3, ~ Group + Batch, test = "Group: G2 vs G1")))
  expect_snapshot(print(buildCacoaModel(meta3, ~ Group + Batch, test = "Group: all")))
  expect_snapshot(print(buildCacoaModel(meta, ~ group * batch, test = "group")))
  expect_snapshot(print(buildCacoaModel(meta, ~ age + batch, test = "age")))
  mm <- meta; mm$site <- "x"; expect_snapshot(print(buildCacoaModel(mm, ~ group + site + batch, test = "group")))
})

test_that("subsetDesign keeps the design full rank and the contrast equivalent; refuses non-estimable contrasts", {
  meta3 <- mbMeta3(); set.seed(3); y <- rnorm(18)
  d <- buildCacoaModel(meta3, ~ Group + Batch, test = "Group: G2 vs G1")$tests[[1]]$design
  rows <- rownames(meta3)[meta3$Group != "G3"]                                 # a cell type without G3 samples
  sd <- cacoa:::subsetDesign(d$F, d$contrast.F, rows)
  expect_equal(colnames(sd$F), c("GroupG1", "GroupG2", "Batchb2")); expect_equal(unname(sd$C[, 1]), c(-1, 1, 0))
  expect_null(cacoa:::subsetDesign(d$F, buildCacoaModel(meta3, ~ Group + Batch, test = "Group: G3 vs G1")$tests[[1]]$design$contrast.F, rows))
  # a dependent column is dropped and the contrast re-expressed: same estimate for any response
  rows2 <- rownames(meta3)[meta3$Batch == "b2"]                                # batch constant: Batchb2 = sum of the level columns
  sd2 <- cacoa:::subsetDesign(d$F, d$contrast.F, rows2)
  expect_equal(ncol(sd2$F), 3); expect_equal(qr(sd2$F)$rank, 3)
  est.full <- drop(t(d$contrast.F) %*% MASS::ginv(d$F[rows2, ]) %*% y[match(rows2, rownames(meta3))])
  est.sub <- drop(t(sd2$C[, 1]) %*% MASS::ginv(sd2$F) %*% y[match(rows2, rownames(meta3))])
  expect_equal(est.sub, est.full, tolerance = 1e-10)
  # confounded within the subset (group = batch): not estimable
  meta <- mbMeta(); meta$batch <- ifelse(meta$group == "A", "b1", "b2"); meta$batch[c(1, 12)] <- c("b2", "b1")
  d2 <- buildCacoaModel(meta, ~ group + batch, test = "group")$tests[[1]]$design
  expect_null(cacoa:::subsetDesign(d2$F, d2$contrast.F, rownames(meta)[-c(1, 12)]))
  expect_false(is.null(cacoa:::subsetDesign(d2$F, d2$contrast.F, rownames(meta))))
  # term contrasts on a subset lacking a level
  dt <- buildCacoaModel(meta3, ~ Group + Batch, test = "Group: all")$tests[[1]]$design
  expect_null(cacoa:::subsetDesign(dt$F, dt$term.contrast, rows))             # "G3 vs G1" is not estimable without G3
  m.sub <- cacoa:::restrictModelToSamples(buildCacoaModel(meta3, ~ Group + Batch, test = "Group: G2 vs G1")$tests[[1]]$design, rows)
  expect_equal(rownames(m.sub$F), rows); expect_equal(ncol(m.sub$Z), 2)
})

test_that("per-column fitter: a rank-deficient subset design is fitted when the contrast is estimable", {
  meta <- mbMeta(); set.seed(4); n <- 12
  d <- buildCacoaModel(meta, ~ group + batch, test = "group")$tests[[1]]$design
  Y <- cbind(full = rnorm(n), nob2 = rnorm(n), noB = rnorm(n))
  Y[meta$batch == "b2", "nob2"] <- NA                                           # batch level absent: estimable
  Y[meta$group == "B", "noB"] <- NA                                             # contrasted level absent: not estimable
  for (sch in c("block", "freedman-lane")) {
    r <- suppressWarnings(cacoa:::performLMPermutations(d, Y, perm.method = sch, n.permutations = 29, seed = 1))
    obs <- meta$batch == "b1"; mm <- meta[obs, ]; mm$group <- factor(mm$group)
    ref <- summary(lm(Y[obs, "nob2"] ~ group, mm))$coefficients["groupB", ]
    expect_equal(unname(c(r$effect[["nob2"]], r$se[["nob2"]], r$stat.obs[["nob2"]])), unname(ref[1:3]), tolerance = 1e-8, info = sch)
    expect_equal(unname(r$df[["nob2"]]), 4L, info = sch)
    expect_true(is.finite(r$pval[["nob2"]]), info = sch)
    expect_true(is.na(r$pval[["noB"]]) && is.na(r$stat.obs[["noB"]]), info = sch)
    expect_true(is.finite(r$pval[["full"]]), info = sch)
  }
})
