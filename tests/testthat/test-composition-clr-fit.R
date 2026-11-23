## The composition model fitted on the clr coordinates in the compiled fitter (fit_and_randomize): the plain fit equals the
## least-squares reference on the same randomizations, and robust fits (Huber, winsorized) are tested by their studentized loadings.

codaClrToy <- function(n = 24, D = 8, seed = 1, eff = 0, outlier = 0) {
  set.seed(seed); ct <- sprintf("T%02d", seq_len(D))
  meta <- data.frame(row.names = sprintf("s%02d", seq_len(n)), cond = factor(rep(c("A", "B"), each = n / 2)), age = round(runif(n, 30, 80)),
                     batch = factor(rep(c("b1", "b2", "b3"), length.out = n)))
  base <- exp(seq(log(300), log(10), length.out = D)); w <- seq(-1, 1, length.out = D)
  cnts <- t(sapply(seq_len(n), function(i) { lp <- log(base) + rnorm(D, 0, 0.35) + 0.02 * (meta$age[i] - 55) * w
    if (meta$cond[i] == "B") lp[2] <- lp[2] + eff
    stats::rmultinom(1, round(runif(1, 1000, 3000)), exp(lp)) }))
  if (outlier > 0) cnts[n / 2 + 1, 3] <- cnts[n / 2 + 1, 3] * outlier
  dimnames(cnts) <- list(rownames(meta), ct); list(cnts = cnts, meta = meta)
}

test_that("the plain clr fit equals the least-squares reference on the same randomizations (block and Freedman-Lane)", {
  skip_if_not_installed("coda.base")
  d <- codaClrToy(eff = log(2))
  for (sc in list(list(f = ~ cond + batch, pm = "block"), list(f = ~ cond + age, pm = "freedman-lane"))) {
    m <- suppressMessages(buildCacoaModel(d$meta, formula = sc$f, test = "cond"))
    x <- suppressMessages(cacoa:::lmCoda(d$cnts, m, perm.method = sc$pm, n.permutations = 199, seed = 1))
    expect_identical(x$settings$per.celltype.statistic, "t")
    Yc <- x$ilr %*% t(x$psi); colnames(Yc) <- rownames(x$psi)
    # the compiled fitter gives sample i the residual of sample P^-1[i] under Freedman-Lane; block randomization permutes the design rows
    P <- x$fit$P; if (sc$pm == "freedman-lane") P <- apply(P, 2, order)
    r <- cacoa:::contrastPermutationTest(Yc, m$F, m$contrast.F, P, scheme = sc$pm)
    expect_lt(max(abs(x$contrast$loadings$obs - r$est)), 1e-10)
    expect_lt(max(abs(x$contrast$loadings$se - r$se)), 1e-10)
    p.ref <- (rowSums(abs(r$t.perm) >= abs(r$t) - 1e-12) + 1) / (ncol(r$t.perm) + 1)
    expect_lt(max(abs(x$contrast$per_cell$pval - p.ref)), 1e-12, label = sc$pm)
    expect_lt(abs(sum(x$contrast$loadings$obs)), 1e-10)                            # the plain loadings sum to zero
    expect_lt(abs(x$contrast$global$stat - sum(r$t^2)), 1e-8)
  }
})

test_that("robust composition fits are studentized, and on clean data stay close to the plain fit", {
  skip_if_not_installed("coda.base")
  d <- codaClrToy(eff = log(2)); m <- suppressMessages(buildCacoaModel(d$meta, formula = ~ cond + age, test = "cond"))
  x0 <- suppressMessages(cacoa:::lmCoda(d$cnts, m, perm.method = "freedman-lane", n.permutations = 199, seed = 1))
  for (rb in c("huber", "winsor")) {
    x <- suppressMessages(cacoa:::lmCoda(d$cnts, m, perm.method = "freedman-lane", n.permutations = 199, seed = 1, robust.method = rb))
    expect_identical(x$settings$per.celltype.statistic, "t", label = rb)
    expect_true(all(is.finite(x$contrast$per_cell$pval)) && all(x$contrast$per_cell$pval > 0 & x$contrast$per_cell$pval <= 1))
    expect_lt(max(abs(x$contrast$loadings$obs - x0$contrast$loadings$obs)), 0.25, label = rb)
    expect_lt(x$contrast$per_cell$pval[["T02"]], 0.05, label = paste(rb, "planted change"))
    expect_true(all(is.finite(x$contrast$loadings$se)))
  }
})

test_that("a robust fit pulls the loading of a cell type with one outlying sample back towards the clean fit", {
  skip_if_not_installed("coda.base")
  clean <- codaClrToy(seed = 3); dirty <- codaClrToy(seed = 3, outlier = 10)
  m <- suppressMessages(buildCacoaModel(clean$meta, formula = ~ cond + age, test = "cond"))
  l <- function(d, rb) suppressMessages(cacoa:::lmCoda(d$cnts, m, perm.method = "freedman-lane", n.permutations = 99, seed = 1, robust.method = rb))$contrast$loadings$obs[["T03"]]
  ref <- l(clean, "none"); dev.plain <- abs(l(dirty, "none") - ref)
  for (rb in c("huber", "winsor")) expect_lt(abs(l(dirty, rb) - ref), dev.plain, label = rb)
})

test_that("robust fits leave out cell types observed in very few samples, which get the plain fit on the same randomizations", {
  skip_if_not_installed("coda.base")
  d <- codaClrToy(eff = log(2)); d$cnts[, "T08"] <- 0L; d$cnts[d$meta$cond == "B", "T08"][1] <- 200L   # T08 present in one sample only
  m <- suppressMessages(buildCacoaModel(d$meta, formula = ~ cond + age, test = "cond"))
  x0 <- suppressMessages(cacoa:::lmCoda(d$cnts, m, perm.method = "freedman-lane", n.permutations = 199, seed = 1))
  expect_message(x <- cacoa:::lmCoda(d$cnts, m, perm.method = "freedman-lane", n.permutations = 199, seed = 1, robust.method = "huber"), "too few samples")
  expect_identical(x$settings$robust.skipped, "T08")
  expect_equal(x$contrast$loadings$obs[["T08"]], x0$contrast$loadings$obs[["T08"]], tolerance = 1e-10)
  expect_equal(x$contrast$per_cell$pval[["T08"]], x0$contrast$per_cell$pval[["T08"]])
  expect_gt(x$contrast$per_cell$pval[["T08"]], 0.05)
})
