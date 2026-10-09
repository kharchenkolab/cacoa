## Permutation inference for the pairwise-effects engine (Track A.4).
##
## All schemes permute SAMPLE labels (never pairs):
##  - "block": the existing block randomization, extended so that only the samples of the compared conditions
##    swap labels, within strata formed by the other discrete covariates (+ block.vars). Exact when every
##    adjusting covariate is discrete.
##  - "freedman-lane": fit the model without the tested term, permute residuals (within block.vars strata),
##    add the fitted part back; done on the Gower matrix (G* = K1 + K2[,p] + K3[p,] + K4[p,p]). Fallback for
##    continuous covariates; conservative when the effective dimension is large relative to n.
##  - "huh-jhun": permutation in the residual space of the reduced model; shift only; opt-in.
## One set of global permutations is shared by all cell types (induced on each cell type's sample subset),
## which allows max-statistic combination across cell types.

# ---- which samples may swap labels, within which strata -------------------------------------------

# Raw metadata variables of the design formula, split into discrete / numeric
designVariables <- function(design, meta) {
  vars <- all.vars(stats::as.formula(design$formula_used))
  vars <- intersect(vars, names(meta))
  isNum <- vapply(vars, function(v) is.numeric(meta[[v]]) && !is.matrix(meta[[v]]), logical(1))
  list(all = vars, numeric = vars[isNum], discrete = vars[!isNum])
}

# Describe the contrast in terms of raw variables: which variables are swapped ("diff"), which settings are
# held fixed, and the two labels being compared. type "labels": factor-level comparison; "free": numeric or
# coefficient contrast (every sample takes part, all labels distinct).
contrastSampleInfo <- function(design, meta, samples) {
  meta <- meta[samples, , drop = FALSE]; n <- nrow(meta)
  spec <- design$contrast_spec
  free <- list(type = "free", in.set = rep(TRUE, n), labels = as.character(seq_len(n)), diff.vars = character(0),
               fixed = list(), tested.vars = character(0))
  if (!is.null(spec) && identical(spec$type, "term")) {       # whole-factor test: every level takes part
    labels <- as.character(meta[[spec$term]])
    return(list(type = "labels", in.set = !is.na(labels), labels = labels, diff.vars = spec$term, fixed = list(),
                num = NA_character_, den = NA_character_, tested.vars = spec$term))
  }
  if (is.null(spec) || !spec$type %in% c("simple", "marginal")) {
    # coefficient-level contrasts: tested variables are those of the terms with non-zero weight
    cF <- design$contrast.F; assign <- attr(design$F, "assign"); tl <- attr(attr(design$F, "terms"), "term.labels")
    if (!is.null(assign) && length(tl)) {
      nz <- unique(assign[abs(cF[colnames(design$F)]) > 1e-12]); nz <- nz[nz > 0]
      free$tested.vars <- unique(unlist(lapply(tl[nz], function(t) all.vars(str2lang(t)))))
    }
    return(free)
  }
  if (grepl(":", spec$term, fixed = TRUE)) {
    vars <- parseTermVars(spec$term)
    numc <- parseCell(spec$num, vars); denc <- parseCell(spec$den, vars)
    diff.vars <- vars[unlist(numc) != unlist(denc)]
    fixed <- c(numc[setdiff(vars, diff.vars)], spec$at %||% list())
    numLab <- paste(unlist(numc[diff.vars]), collapse = ":"); denLab <- paste(unlist(denc[diff.vars]), collapse = ":")
  } else {
    diff.vars <- spec$term
    fixed <- if (spec$type == "simple") (spec$at %||% list()) else list()
    numLab <- spec$num; denLab <- spec$den
  }
  if (any(vapply(meta[diff.vars], is.numeric, logical(1)))) {
    # numeric tested variable (slope / step): every sample with a value takes part and all values are distinct labels
    ok <- stats::complete.cases(meta[, diff.vars, drop = FALSE])
    return(list(type = "labels", in.set = ok, labels = as.character(seq_len(n)), diff.vars = diff.vars, fixed = fixed,
                num = numLab, den = denLab, tested.vars = diff.vars))
  }
  labels <- apply(meta[, diff.vars, drop = FALSE], 1, function(r) paste(as.character(r), collapse = ":"))
  in.set <- labels %in% c(numLab, denLab)
  for (v in names(fixed)) {              # discrete fixed settings must match; numeric ones cannot be matched
    if (v %in% names(meta) && !is.numeric(meta[[v]])) in.set <- in.set & as.character(meta[[v]]) == as.character(fixed[[v]])
  }
  list(type = "labels", in.set = in.set, labels = labels, diff.vars = diff.vars, fixed = fixed,
       num = numLab, den = denLab, tested.vars = diff.vars)
}

# Strata (blocks) over `samples`: all discrete formula variables other than the swapped ones, plus block.vars
permutationStrata <- function(design, meta, samples, info, block.vars = NULL) {
  meta <- meta[samples, , drop = FALSE]
  dv <- designVariables(design, meta)
  svars <- unique(c(setdiff(dv$discrete, info$diff.vars), block.vars))
  svars <- intersect(svars, names(meta))
  if (!length(svars)) return(factor(rep("all", nrow(meta))))
  droplevels(interaction(lapply(meta[svars], function(x) factor(as.character(x))), drop = TRUE, lex.order = TRUE))
}

# log of the number of distinct relabelings: product over strata of multinomial coefficients of the label counts
logDistinctPermutations <- function(labels, in.set, strata) {
  tot <- 0
  for (s in levels(strata)) {
    idx <- which(in.set & strata == s)
    if (length(idx) < 2) next
    cnt <- table(labels[idx])
    tot <- tot + lgamma(length(idx) + 1) - sum(lgamma(cnt + 1))
  }
  tot
}

#' Permutation plan for a sample-level test
#'
#' Decides which samples may swap labels, within which strata, and which scheme to use.
#'
#' @param design output of `buildDesignMatrices()`
#' @param meta sample metadata (rows named by sample)
#' @param samples samples taking part (default: all rows of the design)
#' @param scheme `"auto"` (block when every adjusting covariate is discrete, else Freedman-Lane), `"block"`,
#'   `"freedman-lane"` or `"huh-jhun"`
#' @param block.vars additional metadata columns defining strata
#' @param n.permutations requested number of permutations; when fewer distinct relabelings exist they are
#'   enumerated exhaustively
#' @param max.enumerate largest number of relabelings that is enumerated
#' @return list with `scheme`, `in.set` (logical over `samples`), `strata` (factor), `labels`, `info`,
#'   `n.distinct` (may be `Inf`-like large; see `log.n.distinct`), `exhaustive`, `p.floor`, `notes`
#' @export
permutationPlan <- function(design, meta, samples = rownames(design$F),
                            scheme = c("auto", "block", "freedman-lane", "huh-jhun"),
                            block.vars = NULL, n.permutations = 999, max.enumerate = 20000) {
  scheme <- match.arg(scheme)
  info <- contrastSampleInfo(design, meta, samples)
  dv <- designVariables(design, meta[samples, , drop = FALSE])
  nuisance.numeric <- setdiff(dv$numeric, info$tested.vars)
  notes <- character(0)
  if (scheme == "auto") {
    scheme <- if (length(nuisance.numeric)) "freedman-lane" else "block"
    if (scheme == "freedman-lane")
      notes <- c(notes, sprintf("continuous covariate(s) %s: using Freedman-Lane residual permutation (approximate)",
                                paste(nuisance.numeric, collapse = ", ")))
  } else if (scheme == "block" && length(nuisance.numeric)) {
    notes <- c(notes, sprintf("block permutation with continuous covariate(s) %s is not exact",
                              paste(nuisance.numeric, collapse = ", ")))
  }
  if (scheme == "block") {
    strata <- permutationStrata(design, meta, samples, info, block.vars)
    in.set <- info$in.set; labels <- info$labels
  } else {
    # residual permutation: every sample takes part, strata from block.vars only
    strata <- if (length(block.vars)) permutationStrata(design, meta, samples, list(diff.vars = dv$all), block.vars)
              else factor(rep("all", length(samples)))
    in.set <- rep(TRUE, length(samples)); labels <- as.character(seq_along(samples))
  }
  log.nd <- logDistinctPermutations(labels, in.set, strata)
  n.distinct <- if (log.nd < log(.Machine$double.xmax)) exp(log.nd) else Inf
  exhaustive <- is.finite(n.distinct) && n.distinct <= n.permutations && n.distinct <= max.enumerate
  p.floor <- if (exhaustive) 1 / n.distinct else 1 / (n.permutations + 1)
  if (is.finite(n.distinct) && n.distinct < 200)
    notes <- c(notes, sprintf("only %d distinct permutations: the smallest attainable p-value is %.3g", round(n.distinct), p.floor))
  list(scheme = scheme, in.set = in.set, strata = strata, labels = labels, info = info, samples = samples,
       n.distinct = n.distinct, log.n.distinct = log.nd, exhaustive = exhaustive, p.floor = p.floor,
       n.permutations = if (exhaustive) round(n.distinct) else n.permutations, notes = notes)
}

# ---- generating permutations ---------------------------------------------------------------------

# all distinct permutations of a label vector; returns a list of index vectors p such that labels[p] runs over
# every distinct arrangement (one representative source index per label occurrence)
distinctLabelPermutations <- function(labels) {
  n <- length(labels)
  rec <- function(remaining) {               # remaining: named counts per label
    if (sum(remaining) == 0) return(list(character(0)))
    out <- list()
    for (l in names(remaining)[remaining > 0]) {
      r2 <- remaining; r2[l] <- r2[l] - 1
      for (tail in rec(r2)) out[[length(out) + 1]] <- c(l, tail)
    }
    out
  }
  arrangements <- rec(table(labels))
  pos <- split(seq_len(n), labels)          # source indices per label
  lapply(arrangements, function(arr) {
    used <- lapply(pos, function(x) 0L); p <- integer(n)
    for (i in seq_len(n)) { l <- arr[i]; used[[l]] <- used[[l]] + 1L; p[i] <- pos[[l]][used[[l]]] }
    p
  })
}

#' Permutation plan and draws for a model
#'
#' One set of permutations for every analysis on the same model: the plan of [permutationPlan()] for the model's
#' sample set, the drawn matrix `P` and the cells (stratum x permuted-set groups of row indices) the labels were
#' swapped in. Used by the per-column fitter (`performLMPermutations()`) so composition, density and cluster-free
#' DE share the shift tests' permutations.
#'
#' @param model a `cacoaModel` or a `buildDesignMatrices()` design with `F` and `meta`
#' @param scheme,n.permutations,block.vars see [permutationPlan()]; `block.vars` defaults to the model's
#' @param seed seed for the draws (`NULL`: current RNG state)
#' @param samples sample set (default: rows of the design)
#' @return list: `plan`, `P` (n x B integer matrix), `cells` (list of integer vectors, row indices), `samples`
#' @export
modelPermutations <- function(model, scheme = "auto", n.permutations = 999, block.vars = NULL, seed = NULL, samples = rownames(model$F)) {
  if (is.null(model$meta)) stop("the model carries no sample metadata; cannot build a permutation plan")
  plan <- permutationPlan(model, model$meta, samples, scheme = scheme, block.vars = block.vars %||% model$block.vars, n.permutations = n.permutations)
  P <- if (is.null(seed)) drawPermutations(plan) else withSeed(seed, drawPermutations(plan))
  storage.mode(P) <- "integer"
  cells <- split(which(plan$in.set), droplevels(plan$strata[plan$in.set]))
  cells <- unname(cells[lengths(cells) >= 2])
  list(plan = plan, P = P, cells = cells, samples = samples)
}

#' Draw permutations according to a plan
#'
#' @param plan output of [permutationPlan()]
#' @param B number of permutations (ignored when the plan is exhaustive)
#' @return integer matrix (n x B); column b is `p` such that sample i takes the labels of sample `p[i]`
#'   (i.e. the permuted design is `X[p, ]`); identity for samples outside the permuted set
#' @export
drawPermutations <- function(plan, B = plan$n.permutations) {
  n <- length(plan$in.set); cells <- split(which(plan$in.set), droplevels(plan$strata[plan$in.set]))
  cells <- cells[lengths(cells) >= 2]
  if (plan$exhaustive) {
    per.cell <- lapply(cells, function(idx) distinctLabelPermutations(plan$labels[idx]))
    grid <- expand.grid(lapply(per.cell, seq_along))
    P <- matrix(rep(seq_len(n), nrow(grid)), n, nrow(grid))
    for (b in seq_len(nrow(grid))) for (k in seq_along(cells)) P[cells[[k]], b] <- cells[[k]][per.cell[[k]][[grid[b, k]]]]
    return(P)
  }
  P <- matrix(rep(seq_len(n), B), n, B)
  for (b in seq_len(B)) for (idx in cells) P[idx, b] <- idx[sample.int(length(idx))]
  P
}

# Induce the permutation `p` of the full sample set onto the subset `sub` (positions in the full set) under the
# subset's plan: within each (stratum x permuted-set) cell of the subset, members are reassigned by the rank of
# their images under `p`. Uniform on the subset's relabelings when `p` is uniform, and coupled across subsets.
inducePermutation <- function(p, sub, sub.plan) {
  m <- length(sub); q <- seq_len(m)
  cells <- split(which(sub.plan$in.set), droplevels(sub.plan$strata[sub.plan$in.set]))
  for (idx in cells) {
    if (length(idx) < 2) next
    img <- p[sub[idx]]
    q[idx] <- idx[rank(img, ties.method = "first")]
  }
  q
}

# ---- statistics under a relabeling ---------------------------------------------------------------

# Precompute what does not change when rows of X and Z are permuted
inferencePrecompute <- function(G, X, Z, cvec) {
  hi <- hatInfo(X); n <- nrow(X)
  a <- drop(t(hi$A) %*% cvec)
  list(n = n, q = hi$rank, XtXi = hi$XtXi, A = hi$A, H = hi$H, a = a, cXc = drop(crossprod(cvec, hi$XtXi %*% cvec)),
       trG = sum(diag(G)), dG = diag(G))
}

# shift F, var, total for labels permuted by p (X[p, ], Z[p, ]) against a fixed Gower matrix G.
# O(n^2 q): uses that X'X is invariant under row permutation.
permutedStats <- function(G, X, Z, cvec, z.end, pre, p = NULL, bias.correct = TRUE, need.var = TRUE) {
  n <- pre$n
  if (is.null(p)) p <- seq_len(n)
  a <- pre$a[p]                                  # A_p' c = (A' c)[p]
  Hp <- pre$H[p, p]
  ss <- drop(crossprod(a, G %*% a)) / pre$cXc
  rss <- pre$trG - sum(Hp * G)
  Fs <- ss / (rss / (n - pre$q))
  if (!need.var) return(c(F = Fs, shift = NA, var = NA, total = NA))
  Xp <- X[p, , drop = FALSE]; Ap <- pre$A[, p, drop = FALSE]          # A_p = A[, p]
  AG <- Ap %*% G                                 # q x n
  M <- AG %*% t(Ap)                              # bias-uncorrected effect matrix
  hG <- rowSums(Xp * t(AG))                      # diag(H_p G)
  hGh <- rowSums((Xp %*% M) * Xp)                # diag(H_p G H_p)
  r <- pre$dG - 2 * hG + hGh                     # diag(R G R)
  Zp <- Z[p, , drop = FALSE]
  Rp <- -Hp; diag(Rp) <- diag(Rp) + 1            # residual projection under the relabeling
  gamma <- fitDispersion(Rp, Zp, r)
  s <- drop(Zp %*% gamma)
  if (bias.correct) M <- M - Ap %*% (s * t(Ap))
  shift <- drop(crossprod(cvec, M %*% cvec))
  s.alt <- sum(as.numeric(z.end$num) * gamma); s.ref <- sum(as.numeric(z.end$den) * gamma)
  var <- 2 * (s.alt - s.ref)
  c(F = Fs, shift = shift, var = var, total = shift + var / 2)
}

# Freedman-Lane parts for the reduced design Xr: G* = K1 + K2[, p] + K3[p, ] + K4[p, p]
flGowerParts <- function(G, Xr) {
  n <- nrow(G); Hr <- if (ncol(Xr)) hatInfo(Xr)$H else matrix(0, n, n); Rr <- diag(n) - Hr
  list(K1 = Hr %*% G %*% Hr, K2 = Hr %*% G %*% Rr, K3 = Rr %*% G %*% Hr, K4 = Rr %*% G %*% Rr)
}
flGowerPermute <- function(parts, p) parts$K1 + parts$K2[, p] + parts$K3[p, ] + parts$K4[p, p]

# Huh-Jhun: rotate into the residual space of the reduced model and permute there (shift F only)
hjPrecompute <- function(G, X, cvec, pre) {
  n <- nrow(G); Xr <- X %*% contrastNullBasis(cvec)
  qr_ <- qr(Xr)$rank
  Qr <- qr.Q(qr(Xr), complete = TRUE)[, seq_len(n - qr_) + qr_, drop = FALSE]
  aW <- drop(t(Qr) %*% pre$a)
  list(GW = t(Qr) %*% G %*% Qr, aW = aW, na2 = pre$cXc, k = ncol(Qr))
}
hjStat <- function(hj, p, n, q) {
  GW <- if (is.null(p)) hj$GW else hj$GW[p, p]
  num <- drop(crossprod(hj$aW, GW %*% hj$aW)) / hj$na2
  num / ((sum(diag(GW)) - num) / (n - q))
}

# ---- one cell type --------------------------------------------------------------------------------

# Permutation statistics for one cell type. `P` is the matrix of permutations of this cell type's samples
# (columns), already induced. Returns observed and permuted (F, shift, var, total).
permutationStatsForCellType <- function(eff, plan, P, bias.correct = TRUE) {
  G <- eff$G; X <- eff$X; Z <- eff$Z; cvec <- eff$contrast; z.end <- eff$z.end
  need.var <- all(is.finite(unlist(z.end)))
  pre <- inferencePrecompute(G, X, Z, cvec)
  obs <- permutedStats(G, X, Z, cvec, z.end, pre, NULL, bias.correct, need.var)
  B <- ncol(P)
  perm <- matrix(NA_real_, B, 4, dimnames = list(NULL, c("F", "shift", "var", "total")))
  if (B == 0) return(list(obs = obs, perm = perm))                       # nothing to permute (a single relabeling)
  storage.mode(P) <- "integer"
  znum <- if (need.var) as.numeric(z.end$num) else numeric(ncol(Z)); zden <- if (need.var) as.numeric(z.end$den) else numeric(ncol(Z))
  if (plan$scheme == "block") {                 # C++ kernel; R reference: permutedStats()
    perm[, ] <- permuted_contrast_stats(G, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, znum, zden, P, bias.correct, need.var)
  } else if (plan$scheme == "freedman-lane") {
    parts <- flGowerParts(G, X %*% contrastNullBasis(cvec))
    perm[, ] <- permuted_contrast_stats_fl(parts$K1, parts$K2, parts$K3, parts$K4, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, znum, zden, P, bias.correct, need.var)
  } else {  # huh-jhun: shift only
    hj <- hjPrecompute(G, X, cvec, pre)
    obs["F"] <- hjStat(hj, NULL, pre$n, pre$q)
    for (b in seq_len(B)) perm[b, "F"] <- hjStat(hj, sample.int(hj$k), pre$n, pre$q)
  }
  list(obs = obs, perm = perm)
}

permutationPValue <- function(obs, perm, alternative = c("greater", "two.sided"), exhaustive = FALSE) {
  alternative <- match.arg(alternative)
  perm <- perm[is.finite(perm)]
  if (!is.finite(obs) || !length(perm)) return(NA_real_)
  k <- if (alternative == "greater") sum(perm >= obs - 1e-12) else sum(abs(perm) >= abs(obs) - 1e-12)
  if (exhaustive) k / length(perm) else (k + 1) / (length(perm) + 1)
}

# ---- all cell types -----------------------------------------------------------------------------

#' Permutation tests of pairwise effects for a set of cell types
#'
#' Estimates shift / var / total per cell type with [pairwiseEffectsFromDesign()] and tests them by
#' permuting sample labels. One set of permutations is drawn for the full sample set and induced on each
#' cell type's samples, so the statistics are coupled across cell types and can be combined by their maximum
#' (global p-value and familywise-adjusted per-cell-type p-values).
#'
#' @param D.list named list of sample x sample distance matrices (one per cell type)
#' @param design sample-level design from `buildDesignMatrices()`
#' @param meta sample metadata (rows named by sample)
#' @param dispersion.formula,dist,bias.correct,influence,min.samp.per.level passed to [pairwiseEffectsFromDesign()]
#' @param permutation `"auto"`, `"block"`, `"freedman-lane"` or `"huh-jhun"` (see [permutationPlan()])
#' @param n.permutations number of permutations
#' @param block.vars metadata columns defining additional permutation strata
#' @param seed integer seed (default: drawn from R's RNG, so `set.seed()` applies)
#' @param alpha significance level used for flags
#' @param n.cores cores for the per-cell-type loop (results do not depend on it)
#' @param return.perm.stats keep the permutation statistics per cell type
#' @param verbose print notes
#' @return list: `results` (data.frame, one row per cell type), `global` (max-statistic p per effect),
#'   `fits` (per cell type), `skipped` (data.frame of skipped cell types and reasons), `plan` (global plan),
#'   `perm.stats` (optional), `call.info`
#' @export
testPairwiseEffects <- function(D.list, design, meta, dispersion.formula = NULL, dist = c("cor", "l2", "l1"),
                                permutation = c("auto", "block", "freedman-lane", "huh-jhun"), n.permutations = 999,
                                block.vars = NULL, bias.correct = TRUE, influence = FALSE, min.samp.per.level = 3,
                                seed = NULL, alpha = 0.05, n.cores = 1, return.perm.stats = FALSE, verbose = FALSE) {
  dist <- match.arg(dist); permutation <- match.arg(permutation)
  if (is.null(names(D.list))) names(D.list) <- paste0("CT", seq_along(D.list))
  all.samples <- intersect(rownames(design$F), rownames(meta))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)

  # global plan and permutations shared by all cell types (enumerated exactly when the design is small, so the
  # max-statistic p-values rest on the same relabelings as the per-cell-type p-values)
  gplan <- permutationPlan(design, meta, all.samples, scheme = permutation, block.vars = block.vars,
                           n.permutations = n.permutations)
  if (verbose && length(gplan$notes)) message(paste(gplan$notes, collapse = "\n"))
  P.global <- withSeed(seed, drawPermutations(gplan, n.permutations))

  # effects per cell type
  fits <- lapply(D.list, function(D) pairwiseEffectsFromDesign(D, design, meta, dispersion.formula = dispersion.formula,
                                                               dist = dist, bias.correct = bias.correct, influence = influence,
                                                               min.samp.per.level = min.samp.per.level))
  ok <- vapply(fits, function(f) isTRUE(f$ok), logical(1))
  skipped <- data.frame(celltype = names(fits)[!ok], reason = vapply(fits[!ok], function(f) f$reason, character(1)),
                        stringsAsFactors = FALSE, row.names = NULL)

  runOne <- function(ct) {
    eff <- fits[[ct]]
    plan <- permutationPlan(design, meta, eff$samples, scheme = gplan$scheme, block.vars = block.vars,
                            n.permutations = n.permutations)
    sub <- match(eff$samples, all.samples)
    P.ind <- apply(P.global, 2, inducePermutation, sub = sub, sub.plan = plan)
    if (!is.matrix(P.ind)) P.ind <- matrix(P.ind, ncol = ncol(P.global))          # a single permutation        # coupled permutations
    st.c <- permutationStatsForCellType(eff, plan, P.ind, bias.correct)
    # p-values: exhaustive enumeration when small, otherwise the coupled draws
    st.p <- if (plan$exhaustive) permutationStatsForCellType(eff, plan, drawPermutations(plan), bias.correct) else st.c
    list(plan = plan, obs = st.c$obs, perm.coupled = st.c$perm, perm = st.p$perm)
  }
  cts <- names(fits)[ok]
  runs <- if (n.cores > 1 && length(cts) > 1) sccore::plapply(cts, runOne, n.cores = n.cores, progress = verbose, fail.on.error = TRUE)
          else lapply(cts, runOne)
  names(runs) <- cts

  # per-cell-type p-values
  rows <- lapply(cts, function(ct) {
    r <- runs[[ct]]; e <- fits[[ct]]; ex <- r$plan$exhaustive
    data.frame(celltype = ct, shift = e$shift, var = e$var, total = e$total, ratio = e$ratio,
               shift.norm = e$shift.norm, var.norm = e$var.norm, total.norm = e$total.norm,
               s.ref = e$s.ref, s.alt = e$s.alt, F = unname(r$obs["F"]),
               p.shift = permutationPValue(r$obs["F"], r$perm[, "F"], "greater", ex),
               p.var = permutationPValue(r$obs["var"], r$perm[, "var"], "two.sided", ex),
               p.total = permutationPValue(r$obs["total"], r$perm[, "total"], "two.sided", ex),
               n = e$n, n.ref = e$n.ref, n.alt = e$n.alt,
               se.shift = if (!is.null(e$influence)) unname(e$influence$se["shift"]) else NA_real_,
               se.var = if (!is.null(e$influence)) unname(e$influence$se["var"]) else NA_real_,
               se.total = if (!is.null(e$influence)) unname(e$influence$se["total"]) else NA_real_,
               scheme = r$plan$scheme, n.perm = nrow(r$perm), n.perm.distinct = r$plan$n.distinct,
               exhaustive = ex, p.floor = r$plan$p.floor, stringsAsFactors = FALSE, row.names = NULL)
  })
  results <- do.call(rbind, rows)
  if (is.null(results)) results <- data.frame()
  global <- data.frame(effect = c("shift", "var", "total"), p = NA_real_, stringsAsFactors = FALSE)
  if (nrow(results)) {
    for (eff in c("shift", "var", "total")) results[[paste0("padj.", eff)]] <- stats::p.adjust(results[[paste0("p.", eff)]], "BH")
    # max-statistic combination over cell types from the coupled permutations
    for (eff in c("shift", "var", "total")) {
      col <- if (eff == "shift") "F" else eff
      obs <- vapply(runs, function(r) unname(r$obs[col]), numeric(1))
      perm <- sapply(runs, function(r) r$perm.coupled[, col])
      if (!is.matrix(perm)) perm <- matrix(perm, ncol = length(runs))
      mu <- colMeans(perm, na.rm = TRUE); sdv <- apply(perm, 2, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- NA
      z.obs <- (obs - mu) / sdv; z.perm <- sweep(sweep(perm, 2, mu), 2, sdv, "/")
      if (eff != "shift") { z.obs <- abs(z.obs); z.perm <- abs(z.perm) }
      usable <- is.finite(z.obs)
      if (any(usable)) {
        mx <- apply(z.perm[, usable, drop = FALSE], 1, max, na.rm = TRUE)
        add1 <- as.numeric(!gplan$exhaustive)                                  # exhaustive draws include the identity: no +1
        global$p[global$effect == eff] <- (sum(mx >= max(z.obs[usable]) - 1e-12) + add1) / (length(mx) + add1)
        results[[paste0("p.fwer.", eff)]] <- ifelse(usable, (vapply(z.obs, function(z) sum(mx >= z - 1e-12), numeric(1)) + add1) / (length(mx) + add1), NA_real_)
      } else results[[paste0("p.fwer.", eff)]] <- NA_real_
    }
    # flags
    results$flags <- vapply(seq_len(nrow(results)), function(i) {
      f <- character(0); r <- results[i, ]
      if (is.finite(r$p.var) && r$p.var < alpha && is.finite(r$n.ref) && is.finite(r$n.alt) &&
          max(r$n.ref, r$n.alt) > 2 * min(r$n.ref, r$n.alt)) f <- c(f, "unbalanced+dispersion")
      if (r$p.floor > alpha / 10) f <- c(f, "few-permutations")
      if (r$scheme == "freedman-lane") f <- c(f, "approximate-scheme")
      paste(f, collapse = ";")
    }, character(1))
    results <- results[order(results$shift, decreasing = TRUE), ]
    rownames(results) <- NULL
  }
  notes <- unique(unlist(c(gplan$notes, lapply(runs, function(r) r$plan$notes))))
  if (verbose && length(notes)) message(paste(notes, collapse = "\n"))
  list(results = results, global = global, fits = fits, skipped = skipped, plan = gplan, notes = notes,
       perm.stats = if (return.perm.stats) lapply(runs, `[[`, "perm") else NULL,
       call.info = list(dist = dist, permutation = permutation, n.permutations = n.permutations, seed = seed,
                        dispersion.formula = dispersion.formula, block.vars = block.vars, bias.correct = bias.correct))
}

# run expr with a given seed, restoring the RNG state afterwards
withSeed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  expr
}
