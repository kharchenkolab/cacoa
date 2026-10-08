# Simulation checks of calibration, bias and power across design scenarios (compact subset of the validation
# grid in misc/validation). Slow; skipped on CRAN. Rejection-rate bands are binomial 99% limits at the given
# replicate count so that the test is stable.

skip_on_cran()
skip_if(identical(Sys.getenv("CACOA_SKIP_SLOW"), "true"))

simScenario <- function(sizes, p = 60, shift = 0, disp.ratio = 1, batch = FALSE, batch.conf = FALSE, batch.eff = 0,
                        batch.cos = 0, age = FALSE, age.eff = 0, decay = 1) {
  groups <- names(sizes); g <- factor(rep(groups, sizes), levels = groups); n <- length(g)
  meta <- data.frame(group = g, row.names = sprintf("s%02d", 1:n))
  u.g <- unitVector(p); mu <- matrix(0, n, p)
  delta <- shift * sqrt(p)
  for (k in seq_along(groups)[-1]) mu[g == groups[k], ] <- mu[g == groups[k], ] + rep((k - 1) * delta * u.g, each = sizes[k])
  if (batch) {
    pb <- if (batch.conf) ifelse(g == "B", 0.75, 0.25) else 0.5
    # redraw until every group is seen in both batches (a fully confounded design is correctly refused as not estimable)
    repeat {
      meta$batch <- factor(ifelse(runif(n) < pb, "b2", "b1"), levels = c("b1", "b2"))
      if (all(table(meta$group, meta$batch) >= 1)) break
    }
    mu <- mu + outer(meta$batch == "b2", batch.eff * sqrt(p) * mixedVector(u.g, batch.cos, p))
  }
  if (age) {
    meta$age <- round(50 + 8 * (g == "B") + rnorm(n, 0, 10))
    mu <- mu + outer((meta$age - 50) / 10, age.eff * sqrt(p) * mixedVector(u.g, 0, p))
  }
  sdg <- rep(1, n); sdg[g == "B"] <- sqrt(disp.ratio)
  spec <- (1:p)^(-decay); spec <- sqrt(spec / mean(spec))
  Y <- mu + matrix(rnorm(n * p), n, p) * rep(spec, each = n) * sdg
  rownames(Y) <- rownames(meta)
  list(meta = meta, Y = Y, truth = delta^2, truth.var = 2 * (disp.ratio - 1) * p)
}

runScenario <- function(R, n.perm = 99, formula, contrast = c("group", "B", "A"), permutation = "auto", ...) {
  args <- list(...)                      # `...` is not visible inside replicate()'s wrapper function
  skipped <- character(0)
  out <- do.call(rbind, lapply(seq_len(R), function(i) {
    sim <- do.call(simScenario, args)
    D <- sampleDistanceMatrices(list(ct = sim$Y), dist = "l2")$ct
    des <- suppressMessages(suppressWarnings(cacoa:::buildDesignMatrices(sim$meta, contrast = contrast, formula = formula)))
    res <- testPairwiseEffects(list(ct = D), des, sim$meta, dist = "l2", n.permutations = n.perm, permutation = permutation)
    r <- res$results
    if (!nrow(r)) { skipped <<- c(skipped, res$skipped$reason); r <- data.frame(shift = NA, var = NA, p.shift = NA, p.var = NA, p.total = NA) }
    c(shift = r$shift, var = r$var, p.shift = r$p.shift, p.var = r$p.var, p.total = r$p.total, truth = sim$truth, truth.var = sim$truth.var)
  }))
  attr(out, "skipped") <- skipped
  out
}
band <- function(R, alpha = 0.05) c(qbinom(0.005, R, alpha), qbinom(0.995, R, alpha)) / R
inBand <- function(p, R) { p <- p[!is.na(p)]; r <- mean(p < 0.05); b <- band(R); r >= b[1] & r <= b[2] }

test_that("null calibration: two groups, confounded batch (block), correlated age (Freedman-Lane)", {
  set.seed(101); R <- 120
  a <- runScenario(R, formula = ~ group, sizes = c(A = 8, B = 8))
  expect_true(inBand(a[, "p.shift"], R)); expect_true(inBand(a[, "p.var"], R)); expect_true(inBand(a[, "p.total"], R))
  b <- runScenario(R, formula = ~ group + batch, sizes = c(A = 8, B = 8), batch = TRUE, batch.conf = TRUE, batch.eff = 0.6)
  expect_length(attr(b, "skipped"), 0)
  expect_true(inBand(b[, "p.shift"], R)); expect_true(inBand(b[, "p.var"], R))
  c1 <- runScenario(R, formula = ~ group + age, sizes = c(A = 8, B = 8), age = TRUE, age.eff = 0.3)
  expect_lte(mean(c1[, "p.shift"] < 0.05), band(R)[2])          # Freedman-Lane: at or below the nominal level
  expect_lte(mean(c1[, "p.var"] < 0.05), 0.12)                  # mildly liberal at n = 16 (documented)
})

test_that("estimates are unbiased under confounding, for both signs of the effect-vector angle", {
  set.seed(102); R <- 100
  for (cs in c(0.6, -0.6)) {
    r <- runScenario(R, formula = ~ group + batch, sizes = c(A = 8, B = 8), shift = 0.3, batch = TRUE, batch.conf = TRUE,
                     batch.eff = 0.4, batch.cos = cs)
    truth <- r[1, "truth"]
    expect_lt(abs(mean(r[, "shift"]) - truth), 3 * sd(r[, "shift"]) / sqrt(R) + 0.1 * truth)
  }
})

test_that("dispersion change is estimated without bias and detected; shift test stays at level", {
  set.seed(103); R <- 100
  r <- runScenario(R, formula = ~ group, sizes = c(A = 8, B = 8), disp.ratio = 2)
  expect_lt(abs(mean(r[, "var"]) - r[1, "truth.var"]) / r[1, "truth.var"], 0.08)
  expect_gt(mean(r[, "p.var"] < 0.05), 0.75)
  expect_lte(mean(r[, "p.shift"] < 0.05), band(R)[2])
  # with a continuous covariate coupling the groups the dispersion fit must still be unbiased
  r2 <- runScenario(R, formula = ~ group + age, sizes = c(A = 8, B = 8), disp.ratio = 2, age = TRUE, age.eff = 0.3)
  expect_lt(abs(mean(r2[, "var"]) - r2[1, "truth.var"]) / r2[1, "truth.var"], 0.10)
})

test_that("power: planted shifts are detected under block and Freedman-Lane schemes", {
  set.seed(104); R <- 60
  r <- runScenario(R, formula = ~ group + batch, sizes = c(A = 8, B = 8), shift = 0.6, batch = TRUE, batch.eff = 0.4)
  expect_gt(mean(r[, "p.shift"] < 0.05), 0.6)
  r2 <- runScenario(R, formula = ~ group + age, sizes = c(A = 8, B = 8), shift = 0.6, age = TRUE, age.eff = 0.3)
  expect_gt(mean(r2[, "p.shift"] < 0.05), 0.5)
})

test_that("max-statistic combination across cell types is calibrated under shared sample variation", {
  set.seed(105); R <- 100; n <- 16; n.ct <- 6
  meta <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)), row.names = sprintf("s%02d", 1:n))
  des <- suppressMessages(cacoa:::buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group))
  pg <- replicate(R, {
    U <- matrix(rnorm(n * 20), n, 20)
    Ds <- lapply(seq_len(n.ct), function(k) { Y <- sqrt(.7) * U %*% matrix(rnorm(20 * 50), 20) / sqrt(20) + sqrt(.3) * matrix(rnorm(n * 50), n)
      rownames(Y) <- rownames(meta); sampleDistanceMatrices(list(ct = Y), dist = "l2")$ct })
    names(Ds) <- paste0("ct", seq_len(n.ct))
    r <- testPairwiseEffects(Ds, des, meta, dist = "l2", n.permutations = 99)
    c(global = r$global$p[1], any.fwer = any(r$results$p.fwer.shift < .05))
  })
  expect_true(inBand(pg["global", ], R))
  expect_lte(mean(pg["any.fwer", ] == 1), band(R)[2])
})

# ---- API-level simulations (api note E1 / E3) ---------------------------------------------------

test_that("E1: partial screen with seven covariates is calibrated and separates dispersion from location (slow)", {
  skip_on_cran(); skip_if(Sys.getenv("CACOA_SKIP_SLOW") == "true", "slow test skipped")
  set.seed(2024); n <- 40; p <- 200; R <- 120
  L <- diag(sqrt((1:p)^-1)) * sqrt(p / sum((1:p)^-1))
  mkMeta <- function(n) {
    cond <- factor(rep(c("ctrl", "dis"), length.out = n))
    data.frame(cond = cond, batch = factor(sample(c("b1", "b2", "b3"), n, TRUE)), sex = factor(sample(c("F", "M"), n, TRUE)),
               age = round(50 + 8 * (cond == "dis") + rnorm(n, 0, 10)), rin = rnorm(n, 7, 1), pmi = rexp(n, 1 / 12),
               site = factor(sample(c("s1", "s2"), n, TRUE)), ancestry = factor(sample(c("a1", "a2", "a3"), n, TRUE, prob = c(.6, .3, .1))),
               row.names = sprintf("s%02d", seq_len(n)))
  }
  simY <- function(meta, condVar = 1) {
    X <- model.matrix(~ batch + sex + scale(age) + scale(rin) + ancestry, meta)[, -1]
    B <- matrix(rnorm(ncol(X) * p), ncol(X)) * 0.35
    X %*% B + matrix(rnorm(nrow(meta) * p), nrow(meta)) %*% L * ifelse(meta$cond == "dis", sqrt(condVar), 1)
  }
  one <- function(condVar) {
    meta <- mkMeta(n); Y <- simY(meta, condVar); D <- as.matrix(dist(Y)); dimnames(D) <- list(rownames(meta), rownames(meta))
    sc <- screenCovariates(list(ct = D), meta, mode = "partial", dist = "l2", n.permutations = 99, seed = sample.int(1e6, 1))
    r <- sc$table[sc$table$covariate == "cond", ]
    c(p = r$p, p.disp = r$p.disp, r.eff = r$r.eff)
  }
  null <- t(replicate(R, one(1)))
  disp <- t(replicate(R, one(2)))
  band <- function(k) c(qbinom(0.005, k, 0.05), qbinom(0.995, k, 0.05)) / k
  rej.null <- mean(null[, "p"] < 0.05); rej.disp.null <- mean(null[, "p.disp"] < 0.05)
  expect_lte(rej.null, band(R)[2])                                   # FL partial test: nominal or conservative
  expect_true(rej.disp.null >= band(R)[1] - 0.02 && rej.disp.null <= band(R)[2] + 0.02)
  expect_gt(mean(disp[, "p.disp"] < 0.05), 0.8)                       # a 2x variance change is found by the dispersion score
  expect_lte(mean(disp[, "p"] < 0.05), band(R)[2] + 0.03)             # ... and does not masquerade as a location effect
  expect_true(median(null[, "r.eff"]) > 10 && median(null[, "r.eff"]) < 40)
})

test_that("E3: the global max-statistic p-value across cell types is calibrated under shared sample variation (slow)", {
  skip_on_cran(); skip_if(Sys.getenv("CACOA_SKIP_SLOW") == "true", "slow test skipped")
  set.seed(77); R <- 150; nct <- 12; n <- 24
  meta <- data.frame(x = factor(rep(c("a", "b"), each = n / 2)), row.names = sprintf("s%02d", 1:n))
  des <- suppressMessages(buildDesignMatrices(meta, contrast = c("x", "b", "a")))
  one <- function(shared) {
    U <- matrix(rnorm(n * 30), n)
    D <- lapply(seq_len(nct), function(k) { Y <- sqrt(shared) * U %*% matrix(rnorm(30 * 100), 30) / sqrt(30) + sqrt(1 - shared) * matrix(rnorm(n * 100), n)
      M <- as.matrix(dist(Y)); dimnames(M) <- list(rownames(meta), rownames(meta)); M })
    names(D) <- paste0("ct", seq_len(nct))
    r <- testPairwiseEffects(D, des, meta, dist = "l2", n.permutations = 99, seed = sample.int(1e6, 1), influence = FALSE)
    c(global = r$global$p[r$global$effect == "shift"], any.raw = any(r$results$p.shift < 0.05), any.fwer = any(r$results$p.fwer.shift < 0.05))
  }
  for (shared in c(0, 0.9)) {
    out <- t(replicate(R, one(shared)))
    rej <- mean(out[, "global"] < 0.05); fw <- mean(out[, "any.fwer"] == 1)
    band <- c(qbinom(0.005, R, 0.05), qbinom(0.995, R, 0.05)) / R
    expect_true(rej >= band[1] && rej <= band[2], info = sprintf("shared %.1f: global rejection %.3f", shared, rej))
    expect_true(fw >= band[1] && fw <= band[2], info = sprintf("shared %.1f: FWER %.3f", shared, fw))
    expect_gt(mean(out[, "any.raw"] == 1), 0.1)                       # unadjusted "any cell type" is far above 5%
  }
})
