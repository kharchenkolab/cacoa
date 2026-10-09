## Validation of the Track A engine (estimation + permutation inference) across design scenarios.
## Usage: Rscript validate_trackA.R [n.rep] [n.perm] [n.cores]   (defaults 200, 199, 32)
## Writes results/scenarios.rds (per-replicate rows), results/summary.csv and results/REPORT.md.
## The installed cacoa must be the current build (set R_LIBS to the scratch library if needed).
suppressPackageStartupMessages({ library(cacoa); library(parallel) })
args <- commandArgs(TRUE)
N.REP  <- if (length(args) >= 1) as.integer(args[1]) else 200L
N.PERM <- if (length(args) >= 2) as.integer(args[2]) else 199L
N.CORES <- if (length(args) >= 3) as.integer(args[3]) else 32L
OUT <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "results")

## ---------------------------------------------------------------- simulator
unitv <- function(p) { v <- rnorm(p); v / sqrt(sum(v^2)) }
mixv  <- function(u, cosv, p) { w <- unitv(p); w <- w - sum(w * u) * u; w <- w / sqrt(sum(w^2)); cosv * u + sqrt(1 - cosv^2) * w }
noiseSd <- function(p, decay) { s <- (1:p)^(-decay); sqrt(s / mean(s)) }            # per-gene sd, trace = p

## sc: a scenario row (list). Returns list(meta, Y, truth)
simulate <- function(sc) {
  sizes <- sc$sizes[[1]]; K <- length(sizes); groups <- LETTERS[1:K]
  g <- factor(rep(groups, sizes), levels = groups); n <- length(g); p <- sc$p
  meta <- data.frame(group = g, row.names = sprintf("s%02d", 1:n))
  u.g <- unitv(p); mu <- matrix(0, n, p)
  delta <- sc$shift * sqrt(p)                                   # ||delta|| for B vs A (sc$shift = per-gene shift)
  for (k in 2:K) mu[g == groups[k], ] <- mu[g == groups[k], ] + rep((k - 1) * delta * u.g, each = sizes[k])
  truth <- delta^2
  if (sc$nuisance %in% c("batch", "batch.conf", "batch+age")) {
    pb <- if (sc$nuisance == "batch.conf") ifelse(g == "B", 0.75, 0.25) else 0.5
    meta$batch <- factor(ifelse(runif(n) < pb, "b2", "b1"), levels = c("b1", "b2"))
    if (any(table(meta$batch) < 2)) meta$batch <- factor(rep(c("b1", "b2"), length.out = n))
    u.b <- mixv(u.g, sc$nuis.cos, p)
    mu <- mu + outer(meta$batch == "b2", sc$nuis.eff * sqrt(p) * u.b)
  }
  if (sc$nuisance %in% c("age", "batch+age")) {
    meta$age <- round(50 + 8 * (g == "B") + rnorm(n, 0, 10))     # age correlated with group
    u.a <- mixv(u.g, sc$nuis.cos, p)
    mu <- mu + outer((meta$age - 50) / 10, sc$nuis.eff * sqrt(p) * u.a)
  }
  if (sc$contrast == "interaction" || sc$contrast == "marginal") {   # group x batch interaction effect
    u.i <- mixv(u.g, 0, p)
    mu <- mu + outer(g == "B" & meta$batch == "b2", sc$inter.eff * sqrt(p) * u.i)
  }
  sdg <- rep(1, n); sdg[g == "B"] <- sqrt(sc$disp.ratio)
  E <- matrix(rnorm(n * p), n, p) * rep(noiseSd(p, sc$decay), each = n) * sdg
  Y <- mu + E; rownames(Y) <- rownames(meta)
  list(meta = meta, Y = Y, truth = truth)
}

distanceFor <- function(Y, metric) sampleDistanceMatrices(list(ct = Y), dist = metric)$ct

formulaFor <- function(sc) switch(sc$nuisance,
  none = ~ group, batch = ~ group + batch, batch.conf = ~ group + batch, age = ~ group + age, `batch+age` = ~ group + batch + age)

contrastFor <- function(sc) switch(sc$contrast,
  triple = c("group", "B", "A"),
  marginal = list(type = "marginal", term = "group", num = "B", den = "A", over = "batch", weights = "equal"),
  interaction = list(type = "simple", term = "group:batch", num = "B:b1", den = "A:b1"),
  numeric = c(age = 1))

formulaForContrast <- function(sc) {
  f <- formulaFor(sc)
  if (sc$contrast %in% c("marginal", "interaction")) f <- update(f, ~ . - batch + group * batch)
  f
}

## the old pair-level LM (dev_lm path) for single-factor contrasts: response = squared distance (as the handoff's
## comparison) or the distance itself (as implemented for l2)
oldPairShift <- function(D, meta, des, response = c("d2", "d"), n.perm) {
  response <- match.arg(response)
  pm <- suppressWarnings(suppressMessages(cacoa:::buildPairDesignMatrices(meta, des, dist.type = "shift", verbosity = "none")))
  Dm <- if (response == "d2") D^2 else D
  y <- matrix(Dm[pm$pairs], ncol = 1)
  r <- cacoa:::performLMPermutations(pm, y, perm.method = "freedman-lane", n.permutations = n.perm, alternative = "greater",
                                     return.residuals = TRUE, return.sampled.stats = FALSE)
  c(old.shift = unname(r$stat.obs), old.p = unname(r$pval))
}

oneRep <- function(sc, rep, n.perm) {
  set.seed(1000 * sc$id + rep)
  sim <- simulate(sc)
  D <- distanceFor(sim$Y, sc$metric)
  f <- formulaForContrast(sc)
  des <- tryCatch(suppressMessages(suppressWarnings(buildDesignMatrices(sim$meta, contrast = contrastFor(sc), formula = f))), error = function(e) NULL)
  if (is.null(des)) return(NULL)
  res <- tryCatch(testPairwiseEffects(list(ct = D), des, sim$meta, dist = sc$metric, n.permutations = n.perm, seed = rep,
                                      permutation = sc$scheme), error = function(e) NULL)
  if (is.null(res) || !nrow(res$results)) return(data.frame(id = sc$id, rep = rep, skipped = TRUE))
  r <- res$results
  out <- data.frame(id = sc$id, rep = rep, skipped = FALSE, n = r$n, shift = r$shift, var = r$var, total = r$total,
                    s.ref = r$s.ref, s.alt = r$s.alt, F = r$F, p.shift = r$p.shift, p.var = r$p.var, p.total = r$p.total,
                    scheme.used = r$scheme, truth = sim$truth, truth.var = 2 * (sc$disp.ratio - 1) * sc$p)
  if (sc$contrast == "triple" && sc$metric == "l2") {
    o2 <- tryCatch(oldPairShift(D, sim$meta, des, "d2", n.perm), error = function(e) c(old.shift = NA, old.p = NA))
    o1 <- tryCatch(oldPairShift(D, sim$meta, des, "d", n.perm), error = function(e) c(old.shift = NA, old.p = NA))
    out$old.d2.shift <- o2[["old.shift"]]; out$old.d2.p <- o2[["old.p"]]; out$old.d.shift <- o1[["old.shift"]]; out$old.d.p <- o1[["old.p"]]
  }
  out
}

## ---------------------------------------------------------------- scenario grid
S <- function(label, sizes, nuisance = "none", contrast = "triple", metric = "l2", shift = 0, nuis.eff = 0, nuis.cos = 0,
              disp.ratio = 1, inter.eff = 0, p = 80, decay = 1, scheme = "auto")
  data.frame(label = label, sizes = I(list(sizes)), nuisance = nuisance, contrast = contrast, metric = metric, shift = shift,
             nuis.eff = nuis.eff, nuis.cos = nuis.cos, disp.ratio = disp.ratio, inter.eff = inter.eff, p = p, decay = decay,
             scheme = scheme, stringsAsFactors = FALSE)
two <- c(A = 8, B = 8); three <- c(A = 7, B = 7, C = 5); uneq <- c(A = 4, B = 12); big <- c(A = 20, B = 20)
grid <- rbind(
  ## --- one factor, two groups: null, shift, dispersion, both; metrics
  S("2grp null",                 two),
  S("2grp shift .3",             two, shift = .3),
  S("2grp disp 2x",              two, disp.ratio = 2),
  S("2grp shift .3 + disp 2x",   two, shift = .3, disp.ratio = 2),
  S("2grp null, cor",            two, metric = "cor"),
  S("2grp shift .3, cor",        two, shift = .3, metric = "cor"),
  S("2grp disp 2x, cor",         two, disp.ratio = 2, metric = "cor"),
  S("2grp null, l1",             two, metric = "l1"),
  S("2grp shift .3, l1",         two, shift = .3, metric = "l1"),
  S("2grp disp 2x, l1",          two, disp.ratio = 2, metric = "l1"),
  S("2grp null, isotropic p200", two, p = 200, decay = 0),
  S("2grp shift .2, iso p200",   two, p = 200, decay = 0, shift = .2),
  ## --- unequal sizes
  S("uneq null",                 uneq),
  S("uneq disp 2x (B large)",    uneq, disp.ratio = 2),
  S("uneq shift .3",             uneq, shift = .3),
  ## --- three groups, contrast B vs A (C differs)
  S("3grp null (C shifted)",     three, shift = 0),
  S("3grp shift .3",             three, shift = .3),
  S("3grp null, C disp",         three, shift = 0),
  ## --- discrete nuisance, balanced and confounded, with effect-vector angles
  S("batch null, eff .4",        two, "batch", nuis.eff = .4),
  S("batch shift .3, orth",      two, "batch", shift = .3, nuis.eff = .4, nuis.cos = 0),
  S("batch.conf null, eff .4",   two, "batch.conf", nuis.eff = .4),
  S("batch.conf null, eff .8",   two, "batch.conf", nuis.eff = .8),
  S("batch.conf shift .3, cos+.6", two, "batch.conf", shift = .3, nuis.eff = .4, nuis.cos = .6),
  S("batch.conf shift .3, cos-.6", two, "batch.conf", shift = .3, nuis.eff = .4, nuis.cos = -.6),
  S("batch.conf shift .3, cor",  two, "batch.conf", shift = .3, nuis.eff = .4, nuis.cos = .6, metric = "cor"),
  ## --- continuous nuisance (Freedman-Lane) and its block-scheme misuse
  S("age null, eff .3",          two, "age", nuis.eff = .3),
  S("age shift .3, eff .3",      two, "age", shift = .3, nuis.eff = .3),
  S("age null, eff .3, iso p200", two, "age", nuis.eff = .3, p = 200, decay = 0),
  S("age null, eff .3, n40",     big, "age", nuis.eff = .3),
  S("age shift .2, n40",         big, "age", shift = .2, nuis.eff = .3),
  S("batch+age null",            two, "batch+age", nuis.eff = .3),
  S("batch+age shift .3",        two, "batch+age", shift = .3, nuis.eff = .3),
  ## --- interaction designs: marginal and cell contrasts
  S("marginal null, inter .4",   two, "batch", contrast = "marginal", inter.eff = .4),
  S("marginal shift .3",         two, "batch", contrast = "marginal", shift = .3, inter.eff = .4),
  S("interaction cell null",     c(A = 10, B = 10), "batch", contrast = "interaction", inter.eff = .4),
  S("interaction cell shift .4", c(A = 10, B = 10), "batch", contrast = "interaction", shift = .4, inter.eff = .4),
  ## --- numeric contrast (age slope) with batch nuisance
  S("age slope null",            two, "batch+age", contrast = "numeric", nuis.eff = 0),
  S("age slope eff .3",          two, "batch+age", contrast = "numeric", nuis.eff = .3),
  ## --- forced schemes
  S("batch null, FL forced",     two, "batch", nuis.eff = .4, scheme = "freedman-lane"),
  S("batch shift .3, FL forced", two, "batch", shift = .3, nuis.eff = .4, scheme = "freedman-lane"),
  S("batch null, HJ forced",     two, "batch", nuis.eff = .4, scheme = "huh-jhun"),
  S("2grp disp 2x, HJ forced",   two, disp.ratio = 2, scheme = "huh-jhun")
)
grid$id <- seq_len(nrow(grid))
## special truth adjustments
grid$shift[grid$label == "3grp null (C shifted)"] <- 0; grid$cshift <- 0
grid$cshift[grid$label %in% c("3grp null (C shifted)", "3grp shift .3")] <- .4   # C gets its own shift (handled below)

## C-level shift and C-specific dispersion are handled by post-processing the simulator:
simulateWithC <- function(sc) {
  sim <- simulate_base(sc)
  if (sc$cshift > 0 && "C" %in% levels(sim$meta$group)) {
    p <- ncol(sim$Y); sim$Y[sim$meta$group == "C", ] <- sim$Y[sim$meta$group == "C", ] + rep(sc$cshift * sqrt(p) * unitv(p), each = sum(sim$meta$group == "C"))
  }
  if (sc$label == "3grp null, C disp") sim$Y[sim$meta$group == "C", ] <- sim$Y[sim$meta$group == "C", ] * sqrt(3)
  sim
}
simulate_base <- simulate; simulate <- simulateWithC

## ---------------------------------------------------------------- max-T calibration (several cell types sharing samples)
maxTRep <- function(rep, shared, n.perm, n.ct = 8) {
  set.seed(50000 + rep)
  n <- 16; meta <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)), row.names = sprintf("s%02d", 1:n))
  U <- matrix(rnorm(n * 20), n, 20)
  Ds <- lapply(seq_len(n.ct), function(k) {
    Y <- sqrt(shared) * U %*% matrix(rnorm(20 * 60), 20) / sqrt(20) + sqrt(1 - shared) * matrix(rnorm(n * 60), n)
    rownames(Y) <- rownames(meta); distanceFor(Y, "l2") })
  names(Ds) <- paste0("ct", seq_len(n.ct))
  des <- suppressMessages(buildDesignMatrices(meta, contrast = c("group", "B", "A"), formula = ~ group))
  r <- testPairwiseEffects(Ds, des, meta, dist = "l2", n.permutations = n.perm, seed = rep)
  data.frame(shared = shared, rep = rep, p.global = r$global$p[1], any.fwer = any(r$results$p.fwer.shift < .05),
             any.raw = any(r$results$p.shift < .05), any.bh = any(r$results$padj.shift < .05),
             any.bonf = any(r$results$p.shift < .05 / n.ct))
}

## ---------------------------------------------------------------- run
message(sprintf("Running %d scenarios x %d reps, %d permutations, %d cores", nrow(grid), N.REP, N.PERM, N.CORES))
t0 <- Sys.time()
rows <- mclapply(seq_len(nrow(grid) * N.REP), function(k) {
  i <- (k - 1) %/% N.REP + 1; rep <- (k - 1) %% N.REP + 1
  tryCatch(oneRep(grid[i, ], rep, N.PERM), error = function(e) data.frame(id = i, rep = rep, skipped = TRUE, error = conditionMessage(e)))
}, mc.cores = N.CORES, mc.preschedule = FALSE)
res <- do.call(plyr::rbind.fill, rows)
res <- merge(grid[, c("id", "label", "metric", "nuisance", "contrast", "scheme")], res, by = "id")
maxt <- do.call(rbind, mclapply(seq_len(3 * N.REP), function(k) {
  shared <- c(0, .5, .9)[(k - 1) %/% N.REP + 1]; maxTRep((k - 1) %% N.REP + 1, shared, N.PERM) }, mc.cores = N.CORES))
message(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(list(grid = grid, res = res, maxt = maxt, n.rep = N.REP, n.perm = N.PERM), file.path(OUT, "scenarios.rds"))

## ---------------------------------------------------------------- summary
rej <- function(p) mean(p < .05, na.rm = TRUE)
summ <- do.call(rbind, lapply(split(res, res$id), function(d) {
  sc <- grid[grid$id == d$id[1], ]
  d <- d[!d$skipped, ]
  if (!nrow(d)) return(data.frame(id = sc$id, scenario = sc$label, metric = sc$metric, contrast = sc$contrast, scheme = sc$scheme, n.ok = 0L,
                                  truth = NA, shift.mean = NA, shift.sd = NA, var.mean = NA, truth.var = NA, rej.shift = NA, rej.var = NA,
                                  rej.total = NA, old.d2.shift = NA, old.d2.rej = NA, old.d.rej = NA, stringsAsFactors = FALSE))
  data.frame(id = sc$id, scenario = sc$label, metric = sc$metric, contrast = sc$contrast, scheme = d$scheme.used[1], n.ok = nrow(d),
             # truth for l2: ||delta_g||^2, plus the interaction's share for the marginal contrast (equal weights over batch,
             # interaction only in b2 -> half the interaction vector, orthogonal to the group vector)
             truth = if (sc$metric != "l2" || sc$contrast == "numeric") NA else
                     sc$shift^2 * sc$p + if (sc$contrast == "marginal") 0.25 * sc$inter.eff^2 * sc$p else 0,
             shift.mean = mean(d$shift), shift.sd = sd(d$shift),
             var.mean = mean(d$var), truth.var = if (sc$metric == "l2") 2 * (sc$disp.ratio - 1) * sc$p else NA,
             rej.shift = rej(d$p.shift), rej.var = rej(d$p.var), rej.total = rej(d$p.total),
             old.d2.shift = if ("old.d2.shift" %in% names(d)) mean(d$old.d2.shift, na.rm = TRUE) else NA,
             old.d2.rej = if ("old.d2.p" %in% names(d)) rej(d$old.d2.p) else NA,
             old.d.rej = if ("old.d.p" %in% names(d)) rej(d$old.d.p) else NA, stringsAsFactors = FALSE)
}))
write.csv(summ, file.path(OUT, "summary.csv"), row.names = FALSE)
msum <- do.call(rbind, lapply(split(maxt, maxt$shared), function(d) data.frame(shared = d$shared[1], global.rej = rej(d$p.global),
  any.fwer = mean(d$any.fwer), any.bh = mean(d$any.bh), any.bonf = mean(d$any.bonf), any.raw = mean(d$any.raw))))
write.csv(msum, file.path(OUT, "maxt.csv"), row.names = FALSE)

fmt <- function(x, d = 3) ifelse(is.na(x), "", formatC(x, digits = d, format = "f"))
lines <- c(sprintf("# Track A validation (%d reps, %d permutations, alpha = 0.05; exact-test reference at %d perms = %.3f)",
                   N.REP, N.PERM, N.PERM, floor(.05 * (N.PERM + 1)) / (N.PERM + 1)), "",
           "Binomial 95% band for a 0.05 rejection rate with this many reps: " ,
           sprintf("[%.3f, %.3f]", qbinom(.025, N.REP, .05) / N.REP, qbinom(.975, N.REP, .05) / N.REP), "",
           "| id | scenario | metric | contrast | scheme | truth shift | est. shift (sd) | truth var | est. var | rej shift | rej var | rej total | old pair-LM shift (d2) | old rej (d2) | old rej (d) |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(summ))) { s <- summ[i, ]
  lines <- c(lines, sprintf("| %d | %s | %s | %s | %s | %s | %s (%s) | %s | %s | %s | %s | %s | %s | %s | %s |", s$id, s$scenario, s$metric, s$contrast, s$scheme,
    fmt(s$truth, 2), fmt(s$shift.mean, 2), fmt(s$shift.sd, 2), fmt(s$truth.var, 1), fmt(s$var.mean, 1), fmt(s$rej.shift), fmt(s$rej.var), fmt(s$rej.total),
    fmt(s$old.d2.shift, 2), fmt(s$old.d2.rej), fmt(s$old.d.rej))) }
lines <- c(lines, "", "## Combination across cell types under the null (8 cell types sharing a sample-level factor)", "",
           "| shared variance | global max-T rej | any cell type FWER-sig | any BH-sig | any Bonferroni-sig | any raw p<.05 |", "|---|---|---|---|---|---|")
for (i in seq_len(nrow(msum))) { m <- msum[i, ]; lines <- c(lines, sprintf("| %.1f | %s | %s | %s | %s | %s |", m$shared, fmt(m$global.rej), fmt(m$any.fwer), fmt(m$any.bh), fmt(m$any.bonf), fmt(m$any.raw))) }
writeLines(lines, file.path(OUT, "REPORT.md"))
message("written: ", file.path(OUT, "REPORT.md"))
