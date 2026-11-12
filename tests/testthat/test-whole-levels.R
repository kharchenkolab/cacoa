# exchange of whole levels of a grouping variable: a condition per donor, several samples per donor

donorMeta <- function(n.donors = 6, per.donor = 2, seed = 1) {
  set.seed(seed); d <- rep(paste0("d", seq_len(n.donors)), each = per.donor)
  data.frame(row.names = sprintf("s%02d", seq_along(d)), donor = d, cond = factor(ifelse(as.integer(sub("d", "", d)) %% 2 == 1, "a", "b")),
             region = rep(c("r1", "r2"), length.out = length(d)), x = round(rnorm(length(d)), 2), stringsAsFactors = FALSE)
}
wholeOK <- function(P, meta, by = "donor") all(apply(P, 2, function(p) {
  all(sort(p) == seq_along(p)) &&                                                          # a permutation of the samples
    all(tapply(meta[[by]][p], meta[[by]], function(v) length(unique(v))) == 1) }))         # every level receives the rows of one level

test_that("levels that hold one value are exchanged as wholes; the relabelings are counted and enumerated over the levels", {
  meta <- donorMeta(); m <- buildCacoaModel(meta, ~ cond + (1 | donor), test = "cond"); d <- m$tests[[1]]$design
  expect_equal(m$random.info$structure, "between"); expect_true("donor" %in% m$block.vars)
  plan <- permutationPlan(d, meta, scheme = "auto", block.vars = "donor", n.permutations = 999)
  expect_equal(plan$scheme, "block"); expect_true(all(plan$unit > 0)); expect_equal(nlevels(plan$strata), 1)
  expect_equal(plan$n.distinct, choose(6, 3)); expect_true(plan$exhaustive)
  P <- drawPermutations(plan); expect_equal(ncol(P), 20); expect_true(wholeOK(P, meta))
  lab <- apply(P, 2, function(p) paste(meta$cond[p], collapse = "")); expect_equal(length(unique(lab)), 20)   # every arrangement once
  expect_true(all(apply(P, 2, function(p) all(table(meta$cond[p]) == 6))))
  Pr <- drawPermutations(permutationPlan(d, meta, scheme = "auto", block.vars = "donor", n.permutations = 10), 200)
  expect_true(wholeOK(Pr, meta)); expect_gt(length(unique(apply(Pr, 2, function(p) paste(meta$cond[p], collapse = "")))), 15)
  # switched off: nothing can be exchanged
  expect_equal(permutationPlan(d, meta, scheme = "auto", block.vars = "donor", whole.levels = FALSE)$n.distinct, 1)
  # one sample per level: the grouping term changes nothing
  one <- meta[!duplicated(meta$donor), ]; d1 <- buildCacoaModel(one, ~ cond + (1 | donor), test = "cond")$tests[[1]]
  expect_equal(d1$permutation$n.distinct, choose(6, 3))
})

test_that("levels are exchanged only with levels of the same make-up; mixed levels are relabeled within", {
  # donors with one and with two samples: 4 donors of two samples (2 a, 2 b), 2 donors of one sample (1 a, 1 b)
  meta <- donorMeta()[-c(10, 12), ]; d <- buildCacoaModel(meta, ~ cond + (1 | donor), test = "cond")$tests[[1]]$design
  plan <- permutationPlan(d, meta, block.vars = "donor", n.permutations = 999)
  expect_equal(plan$n.distinct, choose(4, 2) * 2); expect_equal(nlevels(droplevels(plan$strata[plan$unit > 0])), 2)
  P <- drawPermutations(plan); expect_true(wholeOK(P, meta))
  n.per <- table(meta$donor); expect_true(all(apply(P, 2, function(p) all(n.per[meta$donor[p]] == n.per[meta$donor]))))
  # a donor with a make-up of its own keeps its value
  meta2 <- donorMeta()[-12, ]; p2 <- permutationPlan(buildCacoaModel(meta2, ~ cond + (1 | donor), test = "cond")$tests[[1]]$design, meta2, block.vars = "donor", n.permutations = 999)
  expect_true(all(p2$unit[meta2$donor == "d6"] == 0)); expect_equal(p2$n.distinct, choose(5, 2))
  expect_true(any(grepl("5 of the 6 levels.*the others keep their value", p2$notes)))
  # a discrete covariate that varies within the donors: samples are matched by its level, so a row keeps its region
  meta3 <- donorMeta(); d3 <- buildCacoaModel(meta3, ~ cond + region + (1 | donor), test = "cond")$tests[[1]]$design
  p3 <- permutationPlan(d3, meta3, block.vars = "donor", n.permutations = 50); P3 <- drawPermutations(p3, 50)
  expect_true(wholeOK(P3, meta3)); expect_true(all(apply(P3, 2, function(p) all(meta3$region[p] == meta3$region))))
  # a covariate with one value per donor: only donors with the same value are exchanged
  meta4 <- donorMeta(8); meta4$sex <- ifelse(meta4$donor %in% c("d1", "d2", "d3", "d4"), "f", "m")
  d4 <- buildCacoaModel(meta4, ~ cond + sex + (1 | donor), test = "cond")$tests[[1]]$design
  p4 <- permutationPlan(d4, meta4, block.vars = "donor", n.permutations = 999); expect_equal(p4$n.distinct, choose(4, 2)^2)
  P4 <- drawPermutations(p4); expect_true(all(apply(P4, 2, function(p) all(meta4$sex[p] == meta4$sex))))
  # partly: two donors hold both values (exchange within), four hold one (exchanged as wholes)
  meta5 <- donorMeta(); meta5$cond[c(1, 3)] <- "b"; meta5$cond <- factor(meta5$cond)       # d1: b a, d2: b b, d3: a a, d4: b b, d5: a a, d6: b b
  m5 <- buildCacoaModel(meta5, ~ cond + (1 | donor), test = "cond"); expect_equal(m5$random.info$structure, "partly")
  p5 <- m5$tests[[1]]$permutation; expect_equal(p5$n.distinct, 2 * choose(5, 2)); expect_equal(p5$n.whole, 5)
  P5 <- drawPermutations(permutationPlan(m5$tests[[1]]$design, meta5, block.vars = "donor", n.permutations = 999))
  expect_true(all(apply(P5, 2, function(p) all(meta5$donor[p][1:2] == "d1"))))               # the mixed donor stays within itself
})

test_that("induction onto present samples keeps the levels whole, in R and in the compiled fitter", {
  meta <- donorMeta(8); m <- buildCacoaModel(meta, ~ cond + (1 | donor), test = "cond"); d <- m$tests[[1]]$design
  mp <- modelPermutations(d, n.permutations = 40, seed = 3, block.vars = "donor"); gplan <- mp$plan
  expect_true(all(mp$unit > 0)); expect_length(mp$cells, 1)
  set.seed(4); Y <- matrix(rnorm(16 * 3), 16, dimnames = list(rownames(meta), paste0("y", 1:3))); Y[c(2, 7, 12), 2] <- NA; Y[c(1, 2, 5), 3] <- NA
  Fd <- d$F; cF <- d$contrast.F[colnames(Fd)]
  f <- cacoa:::fit_and_randomize(X = Fd, Y = Y, contrast = cF, perm_groups = mp$cells, return_sampled_stats = TRUE, perm_matrix = mp$P, na_mode = "drop",
                                 perm_units = mp$unit, perm_slots = mp$slot)
  for (j in 1:3) {
    obs <- which(!is.na(Y[, j])); sub.plan <- permutationPlan(d, meta, rownames(meta)[obs], block.vars = "donor", n.permutations = 40)
    Q <- apply(mp$P, 2, cacoa:::inducePermutation, sub = obs, sub.plan = sub.plan, plan = gplan)
    expect_true(wholeOK(Q, meta[obs, ]))
    # a donor takes the rows of a donor with as many present samples
    npres <- table(meta$donor[obs]); expect_true(all(apply(Q, 2, function(q) all(npres[meta$donor[obs][q]] == npres[meta$donor[obs]]))))
    if (j == 1) expect_equal(Q, mp$P, ignore_attr = TRUE)                                       # nothing missing: the relabeling itself
    for (b in c(1, 9, 40)) expect_equal(f$sampled_stats[b, j], refContrastT(Fd[obs, , drop = FALSE], Y[obs, j], cF, Q[, b]), tolerance = 1e-10)
  }
  # without the level information the compiled induction would mix the rows of several donors: the plan passes it on
  r <- cacoa:::performLMPermutations(d, Y, perm.method = "block", n.permutations = 40, seed = 3, block.vars = "donor", return.sampled.stats = TRUE)
  expect_equal(unname(r$stats.perm %||% r$sampled.stats), unname(f$sampled_stats), tolerance = 1e-10)
})

test_that("a condition per donor: the test with the grouping term is valid where relabeling single samples is not", {
  skip_on_cran()
  # 8 donors with two samples each and a strong donor effect, no effect of the condition; share of p <= 0.05 over 300 data sets
  meta <- donorMeta(8); mg <- buildCacoaModel(meta, ~ cond + (1 | donor), test = "cond")$tests[[1]]$design
  m0 <- buildCacoaModel(meta, ~ cond, test = "cond")$tests[[1]]$design; set.seed(11)
  Y <- matrix(rnorm(16 * 300), 16) + matrix(rnorm(8 * 300, 0, 2), 8)[rep(1:8, each = 2), ]; rownames(Y) <- rownames(meta)
  pg <- cacoa:::performLMPermutations(mg, Y, perm.method = "block", n.permutations = 199, seed = 1, block.vars = "donor")$pval
  p0 <- cacoa:::performLMPermutations(m0, Y, perm.method = "block", n.permutations = 199, seed = 1)$pval
  expect_lt(mean(pg <= 0.05), 0.09); expect_gt(mean(p0 <= 0.05), 0.09)
  # and it detects a difference between the conditions
  Ye <- Y[, 1:100] + 4 * (meta$cond == "b")
  expect_gt(mean(cacoa:::performLMPermutations(mg, Ye, perm.method = "block", n.permutations = 199, seed = 1, block.vars = "donor")$pval <= 0.05), 0.5)
})
