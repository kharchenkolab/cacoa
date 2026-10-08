# Engine tests: closed-form estimator against raw pair-cell formulas, invariances, metric mapping,
# term statistics, leave-one-out, design-level wrapper with estimability and dispersion endpoints.

oneFactorSetup <- function(n.per.group, group.sd, p = 60, seed = 3) {
  sim <- simulateIndividualModel(n.per.group, p = p, group.effect = 0.4, group.sd = group.sd, seed = seed)
  X <- model.matrix(~ group, sim$meta); Z <- X
  list(sim = sim, X = X, Z = Z)
}

test_that("one-factor shift/var/total equal the raw pair-cell formulas (K = 2 and K = 3)", {
  for (setup in list(list(n = c(A = 7, B = 10), sd = c(A = 1, B = 1.6)),
                     list(n = c(A = 7, B = 10, C = 5), sd = c(A = 1, B = 1.5, C = 0.8)))) {
    s <- oneFactorSetup(setup$n, setup$sd)
    k <- ncol(s$X)
    cvec <- c(0, 1, rep(0, k - 2)); zn <- c(1, 1, rep(0, k - 2)); zd <- c(1, rep(0, k - 1))
    e <- estimatePairwiseEffects(s$sim$D2, s$X, cvec, s$Z, list(num = zn, den = zd))
    m <- pairCellMeans(s$sim$D2, s$sim$meta$group, "A", "B")
    expect_equal(e$shift, unname(m["RT"] - (m["RR"] + m["TT"]) / 2), tolerance = 1e-10)
    expect_equal(e$var,   unname(m["TT"] - m["RR"]), tolerance = 1e-10)
    expect_equal(e$total, unname(m["RT"] - m["RR"]), tolerance = 1e-10)
    expect_equal(e$ratio, unname(m["RT"] / mean(c(m["RR"], m["TT"]))), tolerance = 1e-10)
    expect_true(e$shift.raw > e$shift)      # uncorrected c'Mc is biased upwards
    expect_equal(e$total, e$shift + e$var / 2)
  }
})

test_that("estimates are invariant to sample order, coding and redundant columns", {
  s <- oneFactorSetup(c(A = 7, B = 10, C = 5), c(A = 1, B = 1.5, C = 0.8))
  cvec <- c(0, 1, 0); ze <- list(num = c(1, 1, 0), den = c(1, 0, 0))
  e <- estimatePairwiseEffects(s$sim$D2, s$X, cvec, s$Z, ze)
  o <- sample.int(nrow(s$X))
  e2 <- estimatePairwiseEffects(s$sim$D2[o, o], s$X[o, ], cvec, s$Z[o, ], ze)
  expect_equal(e2[c("shift", "var", "total")], e[c("shift", "var", "total")])
  X0 <- model.matrix(~ 0 + group, s$sim$meta)
  e3 <- estimatePairwiseEffects(s$sim$D2, X0, c(-1, 1, 0), s$Z, ze)
  expect_equal(e3$shift, e$shift)
  e4 <- estimatePairwiseEffects(s$sim$D2, cbind(s$X, 1), c(cvec, 0), s$Z, ze)   # duplicated intercept
  expect_equal(e4$shift, e$shift)
  # reference level is irrelevant
  meta2 <- s$sim$meta; meta2$group <- relevel(meta2$group, "C")
  X5 <- model.matrix(~ group, meta2); cn <- colnames(X5)
  c5 <- setNames(numeric(3), cn); c5["groupB"] <- 1; c5["groupA"] <- -1
  e5 <- estimatePairwiseEffects(s$sim$D2, X5, c5, s$Z, ze)
  expect_equal(e5$shift, e$shift)
})

test_that("metric mapping: cor is half a squared Euclidean distance of centred unit profiles; l1 Gower is PSD", {
  set.seed(5); Y <- matrix(rnorm(12 * 40), 12, 40, dimnames = list(paste0("s", 1:12), NULL))
  D <- sampleDistanceMatrices(list(ct = Y), dist = "cor")$ct
  Yc <- sweep(Y, 2, colMeans(Y)); U <- Yc / sqrt(rowSums(Yc^2))
  expect_equal(unname(D), unname(0.5 * as.matrix(dist(U))^2), tolerance = 1e-12)
  expect_equal(asSquaredDistance(D, "cor"), D)
  Dl2 <- sampleDistanceMatrices(list(ct = Y), dist = "l2")$ct
  expect_equal(asSquaredDistance(Dl2, "l2"), as.matrix(dist(Y))^2)
  Dl1 <- sampleDistanceMatrices(list(ct = Y), dist = "l1")$ct
  expect_gt(min(eigen(cacoa:::gowerCenter(Dl1), symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  # uncentred cosine differs from the centred one
  Dnc <- sampleDistanceMatrices(list(ct = abs(Y)), dist = "cor", center.genes = FALSE)$ct
  Dc  <- sampleDistanceMatrices(list(ct = abs(Y)), dist = "cor", center.genes = TRUE)$ct
  expect_false(isTRUE(all.equal(Dnc, Dc)))
  # Euclidean data: Gower of the squared distance equals the centred Gram matrix
  expect_equal(unname(cacoa:::gowerCenter(as.matrix(dist(Y))^2)), unname(tcrossprod(Yc)), tolerance = 1e-10)
  expect_equal(cacoa:::uncenterGower(cacoa:::gowerCenter(as.matrix(dist(Y))^2)), as.matrix(dist(Y))^2, tolerance = 1e-10)
})

test_that("1-df term F equals the contrast F, and R2.adj is centred at zero under the null", {
  s <- oneFactorSetup(c(A = 8, B = 8), c(A = 1, B = 1))
  sim <- simulateIndividualModel(c(A = 8, B = 8), p = 60, group.effect = 0.3, batch.levels = c("b1", "b2"), batch.effect = 0.3, seed = 9)
  X <- model.matrix(~ group + batch, sim$meta); G <- cacoa:::gowerCenter(sim$D2)
  cvec <- c(0, 1, 0)
  Ft <- cacoa:::termTestGower(G, X, X %*% cacoa:::contrastNullBasis(cvec))
  Fc <- cacoa:::contrastF(G, cacoa:::hatInfo(X), cvec)
  expect_equal(unname(Ft["F"]), Fc)
  expect_equal(unname(Ft["df"]), 1)
  r2 <- replicate(200, {
    simn <- simulateIndividualModel(c(A = 8, B = 8), p = 40, group.effect = 0)
    Xn <- model.matrix(~ group, simn$meta)
    cacoa:::termTestGower(cacoa:::gowerCenter(simn$D2), Xn, Xn[, 1, drop = FALSE])["R2.adj"]
  })
  expect_lt(abs(mean(r2)), 0.02)
  expect_gt(mean(replicate(50, {
    simn <- simulateIndividualModel(c(A = 8, B = 8), p = 40, group.effect = 0)
    Xn <- model.matrix(~ group, simn$meta)
    cacoa:::termTestGower(cacoa:::gowerCenter(simn$D2), Xn, Xn[, 1, drop = FALSE])["R2"]
  })), 0.04)   # raw R2 has a chance level of about 1/(n-1)
})

test_that("dispersion term test separates dispersion from location", {
  pd <- replicate(40, {
    sim <- simulateIndividualModel(c(A = 8, B = 8), p = 60, group.effect = 0, group.sd = c(A = 1, B = 2))
    X <- model.matrix(~ group, sim$meta)
    cacoa:::dispersionTermTest(cacoa:::gowerCenter(sim$D2), X, X, X[, 1, drop = FALSE])["p.disp"]
  })
  expect_gt(mean(pd < 0.05), 0.8)
  pl <- replicate(40, {
    sim <- simulateIndividualModel(c(A = 8, B = 8), p = 60, group.effect = 0.6)
    X <- model.matrix(~ group, sim$meta)
    cacoa:::dispersionTermTest(cacoa:::gowerCenter(sim$D2), X, X, X[, 1, drop = FALSE])["p.disp"]
  })
  expect_lt(mean(pl < 0.05), 0.3)
})

test_that("leave-one-out refits equal direct refits on n - 1 samples", {
  s <- oneFactorSetup(c(A = 6, B = 7), c(A = 1, B = 1.3))
  cvec <- c(0, 1); ze <- list(num = c(1, 1), den = c(1, 0))
  G <- cacoa:::gowerCenter(s$sim$D2)
  loo <- cacoa:::looPairwiseEffects(G, s$X, cvec, s$Z, ze)
  for (i in c(1, 5, 13)) {
    e <- estimatePairwiseEffects(s$sim$D2[-i, -i], s$X[-i, ], cvec, s$Z[-i, ], ze)
    expect_equal(unname(loo$effects[i, ]), c(e$shift, e$var, e$total))
  }
  expect_true(all(is.finite(loo$se)))
})

test_that("design-level wrapper: endpoints, dispersion formula, cell table and skips", {
  sim <- simulateIndividualModel(c(A = 7, B = 9, C = 6), p = 60, group.effect = 0.4, group.sd = c(A = 1, B = 1.5, C = 1),
                                 batch.levels = c("b1", "b2"), batch.effect = 0.3, seed = 21)
  meta <- sim$meta; D <- sqrt(sim$D2)   # l2 distance
  des <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group + batch)
  eff <- pairwiseEffectsFromDesign(D, des, meta, dist = "l2")
  expect_true(eff$ok)
  expect_equal(eff$n.ref, 7); expect_equal(eff$n.alt, 9)
  expect_equal(deparse(eff$dispersion.formula), "~group")
  expect_equal(sort(names(eff$r2)), sort(c("group", "batch")))
  expect_true(all(eff$r2 >= 0 & eff$r2 <= 1))
  # the cell table is a symmetric K x K table with s_a + s_a on the diagonal
  expect_equal(dim(eff$cell.table), c(3, 3))
  expect_equal(eff$cell.table, t(eff$cell.table))
  expect_equal(unname(eff$cell.table["B", "A"] - (eff$cell.table["A", "A"] + eff$cell.table["B", "B"]) / 2), eff$shift)
  expect_equal(unname(eff$cell.table["B", "B"] - eff$cell.table["A", "A"]), eff$var)
  # saturated one-way model: implied cells equal raw pair-cell means
  des1 <- cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group)
  eff1 <- pairwiseEffectsFromDesign(D, des1, meta, dist = "l2")
  for (a in c("A", "B", "C")) for (b in c("A", "B", "C")) {
    m <- pairCellMeans(sim$D2, meta$group, a, b)
    expect_equal(unname(eff1$cell.table[a, b]), unname(if (a == b) m["RR"] else m["RT"]), tolerance = 1e-8)
  }
  # dispersion formula with the batch covariate: var endpoints only differ through group
  eff2 <- pairwiseEffectsFromDesign(D, des, meta, dist = "l2", dispersion.formula = ~ group + batch)
  dz <- eff2$z.end$num - eff2$z.end$den
  expect_equal(eff2$s.alt - eff2$s.ref, sum(dz * eff2$gamma))
  expect_true(all(abs(dz[grepl("batch", names(eff2$gamma))]) < 1e-12))
  expect_equal(unname(eff2$s.alt - eff2$s.ref), unname(eff2$gamma["groupB"]))
  # interaction-cell and marginal contrasts carry endpoints into the dispersion model
  desi <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch,
                                      contrast = list(type = "simple", term = "group:batch", num = "B:b1", den = "A:b1"))
  effi <- pairwiseEffectsFromDesign(D, desi, meta, dist = "l2", dispersion.formula = ~ group)
  expect_true(effi$ok); expect_true(is.finite(effi$var))
  desm <- cacoa:::buildDesignMatrices(meta, formula = ~ group * batch,
                                      contrast = list(type = "marginal", term = "group", num = "B", den = "A", over = "batch"))
  effm <- pairwiseEffectsFromDesign(D, desm, meta, dist = "l2")
  expect_true(effm$ok); expect_true(is.finite(effm$shift))
  # numeric coefficient contrast: shift defined, var/total NA
  meta.age <- cbind(meta, age = rnorm(nrow(meta), 50, 10))
  desn <- cacoa:::buildDesignMatrices(meta.age, contrast = c(age = 1), formula = ~ group + age)
  effn <- pairwiseEffectsFromDesign(D, desn, meta.age, dist = "l2")
  expect_true(effn$ok); expect_true(is.finite(effn$shift)); expect_true(is.na(effn$var))
  # skips: a contrasted level missing, too few samples, confounded design
  sub <- rownames(meta)[meta$group != "B"]
  sk <- pairwiseEffectsFromDesign(D[sub, sub], des, meta, dist = "l2")
  expect_false(sk$ok); expect_match(sk$reason, "absent")
  sub2 <- rownames(meta)[meta$group != "B" | seq_len(nrow(meta)) %in% which(meta$group == "B")[1:2]]   # two B samples left
  sk2 <- pairwiseEffectsFromDesign(D[sub2, sub2], des, meta, dist = "l2")
  expect_false(sk2$ok); expect_match(sk2$reason, "fewer than")
  meta3 <- meta; meta3$batch <- factor(ifelse(meta3$group == "B", "b1", "b2"))
  des3 <- suppressWarnings(cacoa:::buildDesignMatrices(meta3, contrast = c("group", "B", "A"), formula = ~ group + batch))
  sk3 <- pairwiseEffectsFromDesign(D, des3, meta3, dist = "l2")
  expect_false(sk3$ok)
  # influence
  effI <- pairwiseEffectsFromDesign(D, des, meta, dist = "l2", influence = TRUE)
  expect_equal(dim(effI$influence$effects), c(nrow(meta), 3))
  expect_true(all(is.finite(effI$influence$se)))
})

test_that("dispersion fit is unbiased with a continuous covariate and unequal dispersions", {
  # the leverage-corrected per-sample estimator is biased towards the other group when a covariate couples the
  # groups; the model-based fit through E[r] = (R o R) Z gamma is not
  set.seed(31); reps <- 150; p <- 60
  est <- replicate(reps, {
    sim <- simulateIndividualModel(c(A = 8, B = 8), p = p, group.effect = 0, group.sd = c(A = 1, B = sqrt(2)), age = TRUE, age.effect = 0.3)
    sim$meta$age <- sim$meta$age + 8 * (sim$meta$group == "B")          # age correlated with group
    X <- model.matrix(~ group + age, sim$meta); Z <- model.matrix(~ group, sim$meta)
    e <- estimatePairwiseEffects(sim$D2, X, c(0, 1, 0), Z, list(num = c(1, 1), den = c(1, 0)))
    hc2 <- tapply(e$v, sim$meta$group, mean)
    c(var = e$var, var.hc2 = unname(2 * (hc2["B"] - hc2["A"])), s.ref = e$s.ref, s.alt = e$s.alt)
  })
  truth.var <- 2 * (2 - 1) * p
  expect_lt(abs(mean(est["var", ]) - truth.var) / truth.var, 0.06)       # unbiased within MC error
  expect_lt(mean(est["var.hc2", ]), mean(est["var", ]) - 0.05 * truth.var) # the old per-sample estimator is biased low
  expect_lt(abs(mean(est["s.ref", ]) - p) / p, 0.06)
  expect_lt(abs(mean(est["s.alt", ]) - 2 * p) / (2 * p), 0.06)
})
