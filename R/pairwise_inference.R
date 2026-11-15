## Permutation inference for the pairwise-effects engine (Track A.4).
##
## All schemes permute SAMPLE labels (never pairs):
##  - "block": the existing block randomization, extended so that only the samples of the compared conditions
##    swap labels, within strata formed by the other discrete covariates (+ block.vars). Exact when every
##    adjusting covariate is discrete.
##  - "freedman-lane": fit the model without the tested term, permute residuals (within the levels of the discrete
##    covariates and block.vars), add the fitted part back; done on the Gower matrix (G* = K1 + K2[,p] + K3[p,] + K4[p,p]). Fallback for
##    continuous covariates; conservative when the effective dimension is large relative to n.
## One set of global permutations is shared by all cell types (induced on each cell type's sample subset),
## which allows max-statistic combination across cell types.
## Levels of a grouping variable that hold one value of the tested variable (a condition per donor, several samples per donor)
## are exchanged as wholes under "block": see wholeUnitCells().

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
logDistinctPermutations <- function(labels, in.set, strata, unit = NULL, ukey = NULL) {
  tot <- 0
  for (s in levels(strata)) {
    idx <- which(in.set & strata == s)
    if (length(idx) < 2) next
    # whole levels of the grouping variable: the arrangements of the levels' values over the levels
    cnt <- if (!is.null(unit) && unit[idx[1]] > 0) table(ukey[idx][!duplicated(unit[idx])]) else table(labels[idx])
    tot <- tot + lgamma(sum(cnt) + 1) - sum(lgamma(cnt + 1))
  }
  tot
}

# Levels of the grouping variable that hold one value of the tested variable (a condition per donor, several samples per donor) cannot
# change it by an exchange within the level. They are exchanged as wholes: among the levels with the same make-up (the same number of
# samples in every cell of the model's other discrete covariates) the j-th sample of a level takes the design row of the j-th sample of
# another level, so that a relabeling stays a permutation of the samples and the label, with every covariate, moves with the level.
# Levels with a make-up of their own stay as they are. Only for a single grouping variable. Returns NULL when there is nothing of the
# kind, otherwise `strata` (the samples of the levels exchanged with each other form one cell), `unit` (per sample: the index of its
# level when it is exchanged as a whole, 0 otherwise), `slot` (per sample: its cell of the other discrete covariates, by which the
# samples of two levels are matched) and `ukey` (the value of the tested variable).
wholeUnitCells <- function(design, meta, samples, info, block.vars, strata) {
  if (length(block.vars) != 1 || !block.vars %in% names(meta)) return(NULL)
  m <- meta[samples, , drop = FALSE]; n <- nrow(m); g <- as.character(m[[block.vars]])
  tv <- intersect(info$tested.vars, names(m))
  key <- if (length(tv)) do.call(paste, c(lapply(m[tv], as.character), sep = "\r")) else as.character(info$labels)
  slot <- as.integer(permutationStrata(design, meta, samples, info, NULL))
  ins <- info$in.set & !is.na(g); if (!any(ins)) return(NULL)
  lev <- split(which(ins), g[ins])
  single <- vapply(lev, function(ix) length(unique(key[ix])) == 1, logical(1))
  if (!any(single)) return(NULL)
  ord <- lapply(lev[single], function(ix) ix[order(slot[ix], ix)])
  sig <- vapply(ord, function(ix) paste(slot[ix], collapse = ","), character(1)); lkey <- vapply(ord, function(ix) key[ix[1]], character(1))
  unit <- integer(n); st <- as.character(strata); codes <- sort(unique(g)); k <- 0L; n.lev <- 0L
  for (sg in unique(sig)) {
    d <- names(ord)[sig == sg]
    if (length(d) < 2 || length(unique(lkey[d])) < 2) next          # nothing to exchange with, or the same value everywhere
    k <- k + 1L; n.lev <- n.lev + length(d)
    for (dd in d) { unit[ord[[dd]]] <- match(dd, codes); st[ord[[dd]]] <- paste0("whole levels ", k) }
  }
  if (!k) return(NULL)
  list(strata = factor(st), unit = unit, slot = slot, ukey = key, n.levels = n.lev, n.single = sum(single), n.classes = k)
}

#' Permutation plan for a sample-level test
#'
#' Decides which samples may swap labels, within which strata, and which scheme to use.
#'
#' @param design output of `buildDesignMatrices()`
#' @param meta sample metadata (rows named by sample)
#' @param samples samples taking part (default: all rows of the design)
#' @param scheme `"auto"`, `"block"` (the tested variable, i.e. the rows of the design, is exchanged within the strata) or
#'   `"freedman-lane"` (residuals of the model without the tested variable are exchanged). `"auto"`: block, except for a model
#'   with a continuous covariate and no grouping variable (`block.vars`, which include the grouping terms `(1 | g)` of the model),
#'   where residuals are exchanged. Residuals are exchanged within the levels of the model's discrete covariates (and `block.vars`),
#'   as the tested variable is under `"block"`: for a fitted factor with a level per pair of samples, the exchange among all samples
#'   is valid but found an added effect in 16 % of the tests where the exchange within the levels found it in 31 % (gsd) and
#'   15 % against 59 % (distance model). With a grouping variable the rows are exchanged within its levels also when a
#'   continuous covariate is present: the groups are not fitted, so the covariate's coefficient in the model without the tested
#'   variable is fitted partly to the group baselines, and residuals exchanged within groups are not comparable with the observed
#'   ones (9-11 % false positives at a nominal 5 % in calibration; 4-7 % with the exchange of rows)
#' @param block.vars additional metadata columns defining strata
#' @param residuals.within `TRUE` (default): under `"freedman-lane"` the residuals are exchanged within the levels of the discrete
#'   covariates of the model and of `block.vars`. `FALSE`: within the levels of `block.vars` only, i.e. among all samples when there
#'   are none (the behavior before 2026-11; of use when the levels are so small that too few relabelings exist). The default can be
#'   set for a session with `options(cacoa.residuals.within = FALSE)`
#' @param whole.levels `TRUE` (default): under `"block"` with one grouping variable, its levels that hold one value of the tested
#'   variable (a condition per donor, several samples per donor) are exchanged as wholes among the levels with the same make-up
#'   (the same number of samples in every cell of the other discrete covariates): the label, with every covariate, moves with the
#'   level. The test then rests on the levels as the independent units. `FALSE`: such levels keep their value in every relabeling
#' @param n.permutations requested number of permutations; when fewer distinct relabelings exist they are
#'   enumerated exhaustively
#' @param max.enumerate largest number of relabelings that is enumerated
#' @return list with `scheme`, `in.set` (logical over `samples`), `strata` (factor), `labels`, `unit` and `slot` (whole levels: per
#'   sample the index of its level, 0 when it is exchanged on its own, and the cell by which samples of two levels are matched), `info`,
#'   `n.distinct` (may be `Inf`-like large; see `log.n.distinct`), `exhaustive`, `p.floor`, `notes`
#' @export
permutationPlan <- function(design, meta, samples = rownames(design$F),
                            scheme = c("auto", "block", "freedman-lane"),
                            block.vars = NULL, n.permutations = 999, max.enumerate = 20000,
                            residuals.within = getOption("cacoa.residuals.within", TRUE), whole.levels = TRUE) {
  scheme <- match.arg(scheme)
  info <- contrastSampleInfo(design, meta, samples)
  dv <- designVariables(design, meta[samples, , drop = FALSE])
  nuisance.numeric <- setdiff(dv$numeric, info$tested.vars)
  notes <- character(0)
  if (scheme == "auto") {
    scheme <- if (length(nuisance.numeric) && !length(block.vars)) "freedman-lane" else "block"
    if (scheme == "freedman-lane")
      notes <- c(notes, sprintf("continuous covariate(s) %s: using Freedman-Lane residual permutation (approximate)",
                                paste(nuisance.numeric, collapse = ", ")))
    else if (length(nuisance.numeric))
      notes <- c(notes, sprintf("continuous covariate(s) %s with grouping variable(s) %s: the rows of the design are exchanged within the groups (the covariate moves with the tested variable; approximate)",
                                paste(nuisance.numeric, collapse = ", "), paste(block.vars, collapse = ", ")))
  } else if (scheme == "block" && length(nuisance.numeric)) {
    notes <- c(notes, sprintf("block permutation with continuous covariate(s) %s is not exact",
                              paste(nuisance.numeric, collapse = ", ")))
  } else if (scheme == "freedman-lane" && length(nuisance.numeric) && length(block.vars)) {
    notes <- c(notes, sprintf("exchange of residuals within the levels of %s with continuous covariate(s) %s is not valid when a covariate is related to the tested variable (the default exchanges the rows of the design within the groups)",
                              paste(block.vars, collapse = ", "), paste(nuisance.numeric, collapse = ", ")))
  }
  unit <- integer(length(samples)); slot <- unit; ukey <- NULL; crosses.levels <- FALSE
  if (scheme == "block") {
    strata <- permutationStrata(design, meta, samples, info, block.vars)
    in.set <- info$in.set; labels <- info$labels
    wu <- if (isTRUE(whole.levels)) wholeUnitCells(design, meta, samples, info, block.vars, strata) else NULL
    if (!is.null(wu)) {
      strata <- wu$strata; unit <- wu$unit; slot <- wu$slot; ukey <- wu$ukey
      notes <- c(notes, sprintf("%d of the %d levels of '%s' that hold one value of the tested variable are exchanged as wholes, among levels with the same make-up%s",
                                wu$n.levels, wu$n.single, block.vars, if (wu$n.levels < wu$n.single) "; the others keep their value" else ""))
    }
  } else {
    # residual permutation: every sample takes part, within the levels of the discrete covariates (not of the tested variable)
    # and of block.vars
    strata <- if (isTRUE(residuals.within)) permutationStrata(design, meta, samples, list(diff.vars = info$tested.vars), block.vars)
              else if (length(block.vars)) permutationStrata(design, meta, samples, list(diff.vars = dv$all), block.vars)
              else factor(rep("all", length(samples)))
    in.set <- rep(TRUE, length(samples)); labels <- as.character(seq_along(samples))
    nuisance.discrete <- setdiff(dv$discrete, info$tested.vars)
    if (isTRUE(residuals.within) && length(nuisance.discrete))
      notes <- c(notes, sprintf("residuals are exchanged within the levels of %s", paste(nuisance.discrete, collapse = ", ")))
    # residuals that leave the levels of a fitted factor (residuals.within = FALSE with discrete covariates): the relabeled data sets then
    # differ from the observed one in more than scale (see residualScale())
    crosses.levels <- !isTRUE(residuals.within) && length(nuisance.discrete) > 0
  }
  log.nd <- logDistinctPermutations(labels, in.set, strata, unit, ukey)
  n.distinct <- if (log.nd < log(.Machine$double.xmax)) exp(log.nd) else Inf
  exhaustive <- is.finite(n.distinct) && n.distinct <= n.permutations && n.distinct <= max.enumerate
  p.floor <- if (exhaustive) 1 / n.distinct else 1 / (n.permutations + 1)
  if (is.finite(n.distinct) && n.distinct < 200)
    notes <- c(notes, sprintf("only %d distinct permutations: the smallest attainable p-value is %.3g", round(n.distinct), p.floor))
  list(scheme = scheme, in.set = in.set, strata = strata, labels = labels, unit = unit, slot = slot, ukey = ukey, info = info, samples = samples,
       crosses.levels = crosses.levels, n.distinct = n.distinct, log.n.distinct = log.nd, exhaustive = exhaustive, p.floor = p.floor,
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
#' @return list: `plan`, `P` (n x B integer matrix), `cells` (list of integer vectors, row indices), `unit` and `slot` (whole levels of
#'   a grouping variable, see [permutationPlan()]), `samples`
#' @export
modelPermutations <- function(model, scheme = "auto", n.permutations = 999, block.vars = NULL, seed = NULL, samples = rownames(model$F)) {
  if (is.null(model$meta)) stop("the model carries no sample metadata; cannot build a permutation plan")
  plan <- permutationPlan(model, model$meta, samples, scheme = scheme, block.vars = block.vars %||% model$block.vars, n.permutations = n.permutations)
  P <- if (is.null(seed)) drawPermutations(plan) else withSeed(seed, drawPermutations(plan))
  storage.mode(P) <- "integer"
  # all distinct relabelings enumerated (the observed labeling among them): p-values are then counted over the set itself, without
  # adding the observed data set once more (see relabelingPValue())
  attr(P, "exhaustive") <- isTRUE(plan$exhaustive)
  cells <- split(which(plan$in.set), droplevels(plan$strata[plan$in.set]))
  cells <- unname(cells[lengths(cells) >= 2])
  list(plan = plan, P = P, cells = cells, unit = plan$unit, slot = plan$slot, samples = samples)
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
  # cells of whole levels: the members of every level in matching order (by slot, then by position); a permutation of the levels moves
  # the j-th sample of a level onto the j-th sample of its image
  unit <- plan$unit %||% integer(n); whole <- vapply(cells, function(idx) unit[idx[1]] > 0, logical(1))
  members <- lapply(seq_along(cells), function(k) if (whole[k]) wholeUnitMembers(cells[[k]], unit, plan$slot) else NULL)
  moveUnits <- function(mem, perm) { src <- unlist(mem, use.names = FALSE); img <- unlist(mem[perm], use.names = FALSE); list(src = src, img = img) }
  if (plan$exhaustive) {
    per.cell <- lapply(seq_along(cells), function(k) if (whole[k]) distinctLabelPermutations(plan$ukey[vapply(members[[k]], `[`, 1L, 1L)])
                                                      else distinctLabelPermutations(plan$labels[cells[[k]]]))
    grid <- expand.grid(lapply(per.cell, seq_along))
    P <- matrix(rep(seq_len(n), nrow(grid)), n, nrow(grid))
    for (b in seq_len(nrow(grid))) for (k in seq_along(cells)) {
      if (whole[k]) { mv <- moveUnits(members[[k]], per.cell[[k]][[grid[b, k]]]); P[mv$src, b] <- mv$img }
      else P[cells[[k]], b] <- cells[[k]][per.cell[[k]][[grid[b, k]]]]
    }
    return(P)
  }
  P <- matrix(rep(seq_len(n), B), n, B)
  for (b in seq_len(B)) for (k in seq_along(cells)) {
    if (whole[k]) { mv <- moveUnits(members[[k]], sample.int(length(members[[k]]))); P[mv$src, b] <- mv$img }
    else { idx <- cells[[k]]; P[idx, b] <- idx[sample.int(length(idx))] }
  }
  P
}

# the samples `idx` of a cell of whole levels, split by level (in the order of the level index) and ordered within the level by slot
# and position
wholeUnitMembers <- function(idx, unit, slot) {
  mem <- split(idx, unit[idx]); mem <- mem[order(as.integer(names(mem)))]
  lapply(mem, function(ix) ix[order(slot[ix], ix)])
}

# Induce the permutation `p` of the full sample set onto the subset `sub` (positions in the full set) under the
# subset's plan: within each (stratum x permuted-set) cell of the subset, members are reassigned by the rank of
# their images under `p`. Uniform on the subset's relabelings when `p` is uniform, and coupled across subsets.
#
# Whole levels (`plan` = the plan of the full sample set, whose `unit` marks them; see wholeUnitCells()): the levels present in the subset
# are grouped by the make-up of their present samples, and within a group they are reassigned by the rank of the levels they are mapped to
# under `p`, the j-th present sample of a level taking the j-th present sample of the other. A level never receives values from two levels.
inducePermutation <- function(p, sub, sub.plan, plan = NULL) {
  m <- length(sub); q <- seq_len(m)
  if (!is.null(plan)) {                  # the cells of the full plan among the present samples (those of the subset's plan, but for whole levels)
    ins <- plan$in.set[sub]; cells <- split(which(ins), droplevels(plan$strata[sub][ins])); unit <- plan$unit %||% integer(length(plan$in.set)); slot <- plan$slot
    for (idx in cells) {
      if (length(idx) < 2) next
      if (unit[sub[idx[1]]] == 0) { q[idx] <- idx[rank(p[sub[idx]], ties.method = "first")]; next }
      mem <- split(idx, unit[sub[idx]]); mem <- mem[order(as.integer(names(mem)))]
      mem <- lapply(mem, function(ix) ix[order(slot[sub[ix]], sub[ix])])
      sig <- vapply(mem, function(ix) paste(slot[sub[ix]], collapse = ","), character(1))
      for (sg in unique(sig)) {
        k <- which(sig == sg); if (length(k) < 2) next
        img <- vapply(mem[k], function(ix) unit[p[sub[ix[1]]]], integer(1))          # the level each level is mapped to
        q[unlist(mem[k], use.names = FALSE)] <- unlist(mem[k][rank(img, ties.method = "first")], use.names = FALSE)
      }
    }
    return(q)
  }
  if (any((sub.plan$unit %||% 0L) > 0)) stop("inducePermutation(): the plan of the full sample set is needed for whole levels")
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
permutedStats <- function(G, X, Z, cvec, z.end, pre, p = NULL, bias.correct = TRUE, need.var = TRUE, z.fixed = integer(0)) {
  n <- pre$n
  if (is.null(p)) p <- seq_len(n)
  a <- pre$a[p]                                  # A_p' c = (A' c)[p]
  Hp <- pre$H[p, p]
  ss <- drop(crossprod(a, G %*% a)) / pre$cXc
  rss <- pre$trG - sum(Hp * G)
  Fs <- ss / (rss / (n - pre$q))
  if (!need.var) return(c(F = Fs, shift = NA, var = NA, total = NA, var.t = NA))
  Xp <- X[p, , drop = FALSE]; Ap <- pre$A[, p, drop = FALSE]          # A_p = A[, p]
  AG <- Ap %*% G                                 # q x n
  M <- AG %*% t(Ap)                              # bias-uncorrected effect matrix
  hG <- rowSums(Xp * t(AG))                      # diag(H_p G)
  hGh <- rowSums((Xp %*% M) * Xp)                # diag(H_p G H_p)
  r <- pre$dG - 2 * hG + hGh                     # diag(R G R)
  Zp <- Z[p, , drop = FALSE]; if (length(z.fixed)) Zp[, z.fixed] <- Z[, z.fixed]   # sample properties stay with the sample
  Rp <- -Hp; diag(Rp) <- diag(Rp) + 1            # residual projection under the relabeling
  gamma <- fitDispersion(Rp, Zp, r)
  s <- drop(Zp %*% gamma)
  if (bias.correct) M <- M - Ap %*% (s * t(Ap))
  shift <- drop(crossprod(cvec, M %*% cvec))
  s.alt <- sum(as.numeric(z.end$num) * gamma); s.ref <- sum(as.numeric(z.end$den) * gamma)
  var <- 2 * (s.alt - s.ref)
  c(F = Fs, shift = shift, var = var, total = shift + var / 2, var.t = dispersionT((Rp * Rp) %*% Zp, r, gamma, as.numeric(z.end$num) - as.numeric(z.end$den)))
}

# studentized dispersion change d' gamma / se (R reference of disp_t() in src/gower_stats.h): the residual variance of the
# least-squares fit of the squared residual distances r on W = (R o R) Z, on n - ncol(Z) degrees of freedom
dispersionT <- function(W, r, gamma, d, w = rep(1, length(r))) {
  dfw <- sum(w) - ncol(W); if (dfw <= 0.5 || any(!is.finite(gamma))) return(NA_real_)
  res <- r - drop(W %*% gamma); s2 <- sum(w * res^2) / dfw
  Ws <- W * sqrt(w); v <- drop(t(d) %*% MASS::ginv(crossprod(Ws)) %*% d) * s2
  if (!is.finite(v) || v <= 0) return(NA_real_)
  sum(d * gamma) / sqrt(v)
}

# Freedman-Lane parts for the reduced design Xr: G* = K1 + K2[, p] + K3[p, ] + K4[p, p], the Gower matrix of the reduced fit plus the
# exchanged residuals. The residuals of the reduced model are smaller than the errors; with `kappa` (flKappa()) they are multiplied by
# sqrt(kappa), so that the statistics that are not ratios (the shift, the total, the null bands) are on the scale of the observed data.
flGowerParts <- function(G, Xr, kappa = 1) {
  n <- nrow(G); Hr <- if (ncol(Xr)) hatInfo(Xr)$H else matrix(0, n, n); Rr <- diag(n) - Hr
  list(K1 = Hr %*% G %*% Hr, K2 = sqrt(kappa) * (Hr %*% G %*% Rr), K3 = sqrt(kappa) * (Rr %*% G %*% Hr), K4 = kappa * (Rr %*% G %*% Rr))
}
# variance of the observed over the mean variance of the relabeled contrast estimate under the exchange of residuals of the reduced design
# Xr: |a|^2 / mean_b(|a|^2 - |Qr' a_b|^2), a = the contrast's weights over the samples, a_b = the weights moved with the residuals
# (sample i receives the residual of sample P[i, b]). 1 when the exchange stays within the levels of the reduced design's terms.
flKappa <- function(a, Xr, P) {
  if (!ncol(Xr) || !length(P) || !ncol(P)) return(1)
  q <- qr(Xr); if (q$rank < 1) return(1); Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]; aa <- sum(a^2); n <- length(a)
  g <- mean(apply(P, 2, function(p) { ab <- numeric(n); ab[p] <- a; sum(crossprod(Q, ab)^2) }))
  if (aa > 0 && aa - g > 1e-8 * aa) aa / (aa - g) else 1
}
flGowerPermute <- function(parts, p) parts$K1 + parts$K2[, p] + parts$K3[p, ] + parts$K4[p, p]

# Exchange of residuals: how the statistics of the relabeled data sets are put on the scale of the observed one (expression shifts: gsd and
# the distance model's F). The exchanged residuals are projections of the data, R0 x, so the variance of the relabeled contrast estimate and
# the expected residual sum of squares of the relabeled fit differ between relabelings, by factors common to all genes.
#   "average":     one factor for all relabelings (the behavior before this option existed);
#   "relabeling":  both quantities for every relabeling;
#   "independent" (default): the same, and the residuals of the relabeled fit are made independent of the relabeled contrast estimate
#                  (the two are correlated, which makes the relabeled statistics smaller and lighter-tailed than the observed ones);
#   "observed":    gsd only: as "relabeling", scaled to the observed arrangement (the distance model takes "relabeling").
# options(cacoa.residual.scale = ...) sets it for a session. When the residuals leave the levels of a fitted factor (the former exchange
# among all samples, options(cacoa.residuals.within = FALSE)) the default is "average": the exact scales are not valid there (gsd: 9-18 % of
# tests at p <= 0.05 without an effect, paired data).
residualScale <- function(crosses.levels = FALSE) {
  v <- getOption("cacoa.residual.scale", if (isTRUE(crosses.levels)) "average" else "independent")
  if (length(v) != 1 || !(v %in% c("average", "relabeling", "observed", "independent"))) stop("option cacoa.residual.scale must be 'independent', 'relabeling', 'observed' or 'average'")
  v
}

# Shift F of the distance model for every relabeling under the exchange of residuals, on the scale of the observed F. With K4 = R0 G R0 (the
# Gower matrix of the residuals of the reduced design Xr), a_b the contrast's weights moved with the residuals (sample i receives the
# residual of sample P[i, b]) and H_b the hat matrix of the relabeled arrangement:
#   ss_b = a_b' K4 a_b / cXc, expected in proportion to c_b = |R0 a_b|^2 / |a|^2;  rss_b = tr((I - H_b) K4), expected in proportion to
#   df_b = tr((I - H_b) R0)  (n - rank(X) for the observed arrangement);
#   "relabeling":  F_b = (ss_b / c_b) / (rss_b / df_b)  (the kernel's F_b, ss_b / (rss_b / df), times df_b / (df c_b));
#   "independent": the residuals of the relabeled fit without their part along the relabeled contrast estimate, k_b = (I - H_b) R0 a_b:
#                  rss_b - 2 k_b' K4 a_b / |R0 a_b|^2 + a_b' K4 a_b |k_b|^2 / |R0 a_b|^4 on df_b - |k_b|^2 / |R0 a_b|^2 degrees of freedom.
# Nothing changes when the exchange maps the reduced space onto itself (c_b = 1, df_b = df, k_b = 0).
flShiftF <- function(G, X, cvec, Xr, P, scale = c("independent", "relabeling")) {
  scale <- match.arg(scale); n <- nrow(X); hi <- hatInfo(X); H <- hi$H; a <- drop(t(hi$A) %*% cvec); aa <- sum(a^2); cXc <- drop(crossprod(cvec, hi$XtXi %*% cvec))
  Hr <- if (ncol(Xr)) hatInfo(Xr)$H else matrix(0, n, n); R0 <- diag(n) - Hr; K4 <- R0 %*% G %*% R0; trK <- sum(diag(K4)); tr0 <- sum(diag(R0))
  apply(P, 2, function(p) {
    ab <- numeric(n); ab[p] <- a; at <- drop(R0 %*% ab); att <- sum(at^2); if (!(att > 1e-12 * aa)) return(0)   # the contrast lies in the reduced space: no estimate
    K4a <- drop(K4 %*% ab); ss <- sum(ab * K4a); rss <- trK - sum(H * K4[p, p]); df <- tr0 - sum(H * R0[p, p])
    if (scale == "independent") { hb <- numeric(n); hb[p] <- drop(H %*% at[p]); k <- at - hb; kk <- sum(k^2); rss <- rss - 2 * sum(k * K4a) / att + ss * kk / att^2; df <- df - kk / att }
    if (!(rss > 0) || !(df > 1e-8)) return(NA_real_)
    (ss / cXc / (att / aa)) / (rss / df)
  })
}

# Can the relabeling plan test a change of dispersion? The dispersion contrast is estimated from the samples' residual distances; it
# is tested by moving the tested variable between samples within the exchange cells. A cell contributes only if the tested variable
# varies in it and its samples' residual distances are not forced equal: two samples that form a level of a fitted factor have mirror-image
# residuals (R[i, ] = -R[j, ]) for any data and covariates, hence equal distances, so a pair holding both values of the tested variable says
# nothing about their dispersions (paired designs with the pair fitted). Cells of whole levels (levels exchanged as wholes) contribute.
# When no cell contributes, the dispersion contrast rests on samples whose value never changes, and its permutation p-value is arbitrary.
# X, cvec: the location design of the cell type's samples and the contrast; plan: permutationPlan() of these samples
dispersionRandomizable <- function(X, cvec, plan) {
  hi <- hatInfo(X); n <- nrow(X); a <- drop(t(hi$A) %*% cvec); R <- diag(n) - hi$H; tol <- 1e-8 * max(1, max(abs(a)))
  unit <- plan$unit %||% integer(n)
  cells <- split(which(plan$in.set), droplevels(as.factor(plan$strata)[plan$in.set]))
  for (idx in cells) {
    if (length(idx) < 2) next
    if (any(unit[idx] > 0)) return(TRUE)                                    # whole levels are exchanged
    if (diff(range(a[idx])) <= tol) next                                     # the tested variable does not vary within the cell
    if (length(idx) >= 3 || max(abs(R[idx[1], ] + R[idx[2], ])) > 1e-8) return(TRUE)
  }
  FALSE
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
  perm <- matrix(NA_real_, B, 5, dimnames = list(NULL, c("F", "shift", "var", "total", "var.t")))
  if (B == 0) return(list(obs = obs, perm = perm))                       # nothing to permute (a single relabeling)
  storage.mode(P) <- "integer"
  znum <- if (need.var) as.numeric(z.end$num) else numeric(ncol(Z)); zden <- if (need.var) as.numeric(z.end$den) else numeric(ncol(Z))
  zfix <- as.integer(eff$z.fixed %||% integer(0)) - 1L   # dispersion-design columns that stay with their sample (0-based)
  if (isTRUE(eff$weighted)) {                   # weighted / robust path (C++; R reference: weightedContrastStats())
    code <- c(none = 0L, huber = 1L, winsor = 2L)[[eff$robust]]
    id <- permuted_contrast_stats_w(G, X, Z, eff$base.w, cvec, znum, zden, matrix(seq_len(nrow(X)), ncol = 1), FALSE, eff$w, zfix, bias.correct, need.var, code, eff$robust.k, 5L)
    obs <- c(F = eff$F.shift, shift = eff$shift, var = if (need.var) eff$var else NA, total = if (need.var) eff$total else NA, var.t = if (need.var) id[1, 5] else NA)
    perm[, ] <- permuted_contrast_stats_w(G, X, Z, eff$base.w, cvec, znum, zden, P, plan$scheme == "freedman-lane", eff$w, zfix, bias.correct, need.var, code, eff$robust.k, 5L)
    return(list(obs = obs, perm = perm))
  }
  if (plan$scheme == "block") {                 # C++ kernel; R reference: permutedStats()
    perm[, ] <- permuted_contrast_stats(G, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, znum, zden, P, zfix, bias.correct, need.var)
  } else {                                      # freedman-lane: Gower re-indexing of the reduced-model parts
    Xr <- X %*% contrastNullBasis(cvec); parts <- flGowerParts(G, Xr, flKappa(pre$a, Xr, P))
    perm[, ] <- permuted_contrast_stats_fl(parts$K1, parts$K2, parts$K3, parts$K4, X, Z, pre$A, pre$H, pre$a, cvec, pre$cXc, pre$q, znum, zden, P, bias.correct, need.var)
    # the shift F on the scale of the observed one for every relabeling (the kernel's F is right on average only; see flShiftF())
    sc <- residualScale(plan$crosses.levels); if (sc != "average") perm[, "F"] <- flShiftF(G, X, cvec, Xr, P, if (sc == "independent") "independent" else "relabeling")
  }
  list(obs = obs, perm = perm)
}

# 2.5-97.5 % of the normalized shift and of var.norm = log2(s.alt / s.ref) over relabelings (what noise alone gives), on the
# observed scale S = s.ref + s.alt: shift / S, and the arm dispersions (S +- var / 2) / 2
relabelingBands <- function(perm, S, var.ok = TRUE, probs = c(0.025, 0.975)) {
  q <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) c(NA_real_, NA_real_) else unname(stats::quantile(x, probs)) }
  if (!is.finite(S) || S <= 0) return(list(shift = c(NA_real_, NA_real_), var = c(NA_real_, NA_real_)))
  v <- perm[, "var"]; ok <- is.finite(v) & abs(v) / 2 < S
  list(shift = q(perm[, "shift"] / S), var = if (var.ok) q(log2((S + v[ok] / 2) / (S - v[ok] / 2))) else c(NA_real_, NA_real_))
}

# p-value from the number k of relabeled statistics at least as large as the observed one, out of B relabelings: (k + 1) / (B + 1)
# for random draws (the observed data set counts as one more); k / B when all distinct relabelings were enumerated, because the
# observed labeling is one of them (k >= 1 then; forced in case rounding lost the tie)
relabelingPValue <- function(k, B, exhaustive = FALSE) if (exhaustive) pmin(1, pmax(k, 1) / B) else (k + 1) / (B + 1)

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
#' @param permutation `"auto"`, `"block"` or `"freedman-lane"` (see [permutationPlan()]). Under `"freedman-lane"` (exchange of residuals) the shift F of
#'   every relabeled data set is put on the scale of the observed one (`options(cacoa.residual.scale = "average")` restores the single
#'   average factor)
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
#' @param robust,na.mode,robust.k robust fit (`"none"`, `"huber"`, `"winsor"`), treatment of samples absent from a cell type
#' @param cell.number.adjust per-sample 1 / n.cells term in the dispersion design, `var` at equal cell numbers (default TRUE; needs the
#'   `n.cells` attribute of the distance matrices)
#'   (`"drop"`, `"impute_weak"`) and robust tuning constant; recorded in `call.info`
testPairwiseEffects <- function(D.list, design, meta, dispersion.formula = NULL, dist = c("cor", "l2", "l1"),
                                permutation = c("auto", "block", "freedman-lane"), n.permutations = 999,
                                block.vars = NULL, bias.correct = TRUE, influence = FALSE, min.samp.per.level = 3,
                                seed = NULL, alpha = 0.05, n.cores = 1, return.perm.stats = FALSE, verbose = FALSE,
                                robust = c("none", "huber", "winsor"), na.mode = c("drop", "impute_weak"), robust.k = 1.345, cell.number.adjust = TRUE) {
  dist <- match.arg(dist); permutation <- match.arg(permutation); robust <- match.arg(robust); na.mode <- match.arg(na.mode)
  block.vars <- block.vars %||% design$block.vars   # the model's permutation strata unless overridden
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
                                                               min.samp.per.level = min.samp.per.level, robust = robust, na.mode = na.mode, robust.k = robust.k,
                                                               cell.number.adjust = cell.number.adjust))
  ok <- vapply(fits, function(f) isTRUE(f$ok), logical(1))
  skipped <- data.frame(celltype = names(fits)[!ok], reason = vapply(fits[!ok], function(f) f$reason, character(1)),
                        stringsAsFactors = FALSE, row.names = NULL)

  runOne <- function(ct) {
    eff <- fits[[ct]]
    plan <- permutationPlan(design, meta, eff$samples, scheme = gplan$scheme, block.vars = block.vars,
                            n.permutations = n.permutations)
    sub <- match(eff$samples, all.samples)
    P.ind <- apply(P.global, 2, inducePermutation, sub = sub, sub.plan = plan, plan = gplan)
    if (!is.matrix(P.ind)) P.ind <- matrix(P.ind, ncol = ncol(P.global))          # a single permutation        # coupled permutations
    st.c <- permutationStatsForCellType(eff, plan, P.ind, bias.correct)
    # p-values: exhaustive enumeration when small, otherwise the coupled draws
    st.p <- if (plan$exhaustive) permutationStatsForCellType(eff, plan, drawPermutations(plan), bias.correct) else st.c
    # relabelings for the noise band: the half least associated with the tested contrast (|cor| of the permuted and the original
    # per-sample contrast weights), so that relabelings close to the true labels (or to their mirror) do not carry the effect into it
    Pp <- if (plan$exhaustive) drawPermutations(plan) else P.ind
    a <- tryCatch(drop(t(hatInfo(eff$X)$A) %*% eff$contrast), error = function(e) NULL)
    band <- rep(TRUE, ncol(Pp))
    if (!is.null(a) && stats::sd(a) > 0 && ncol(Pp) >= 20) {
      assoc <- abs(apply(Pp, 2, function(q) suppressWarnings(stats::cor(a[q], a)))); assoc[!is.finite(assoc)] <- 1
      band <- assoc <= stats::median(assoc)
    }
    disp.rand <- dispersionRandomizable(eff$X, eff$contrast, plan)
    if (!disp.rand) plan$notes <- c(plan$notes, "the change of dispersion is not tested: every exchange cell in which the tested variable varies is a pair of samples forming a level of a fitted factor, whose residual distances are equal (fit the pairing variable as a grouping term to test it)")
    list(plan = plan, obs = st.c$obs, perm.coupled = st.c$perm, perm = st.p$perm, band = band, disp.rand = disp.rand)
  }
  cts <- names(fits)[ok]
  runs <- if (n.cores > 1 && length(cts) > 1) sccore::plapply(cts, runOne, n.cores = n.cores, progress = verbose, fail.on.error = TRUE)
          else lapply(cts, runOne)
  names(runs) <- cts

  # per-cell-type p-values
  rows <- lapply(cts, function(ct) {
    r <- runs[[ct]]; e <- fits[[ct]]; ex <- r$plan$exhaustive
    di <- (is.null(e$disp.identifiable) || isTRUE(e$disp.identifiable)) && !isFALSE(r$disp.rand)   # arm dispersions identifiable under the location model, and testable by the plan
    nb <- relabelingBands(r$perm[r$band, , drop = FALSE], e$s.ref + e$s.alt, di)
    data.frame(celltype = ct, shift = e$shift, var = if (di) e$var else NA_real_, total = if (di) e$total else NA_real_, ratio = e$ratio,
               shift.norm = e$shift.norm, var.norm = if (di) e$var.norm else NA_real_, total.norm = if (di) e$total.norm else NA_real_,
               s.ref = e$s.ref, s.alt = e$s.alt, F = unname(r$obs["F"]), disp.identifiable = di,
               p.shift = permutationPValue(r$obs["F"], r$perm[, "F"], "greater", ex),
               p.var = if (di) permutationPValue(r$obs["var.t"], r$perm[, "var.t"], "two.sided", ex) else NA_real_,
               p.total = if (di) permutationPValue(r$obs["total"], r$perm[, "total"], "two.sided", ex) else NA_real_,
               n = e$n, n.ref = e$n.ref, n.alt = e$n.alt,
               se.shift = if (!is.null(e$influence)) unname(e$influence$se["shift"]) else NA_real_,
               se.var = if (!is.null(e$influence)) unname(e$influence$se["var"]) else NA_real_,
               se.total = if (!is.null(e$influence)) unname(e$influence$se["total"]) else NA_real_,
               scheme = r$plan$scheme, n.perm = nrow(r$perm), n.perm.distinct = r$plan$n.distinct,
               exhaustive = ex, p.floor = r$plan$p.floor,
               null.shift.lo = nb$shift[1], null.shift.hi = nb$shift[2], null.var.lo = nb$var[1], null.var.hi = nb$var[2],
               cell.assoc = e$cell.assoc %||% NA_real_, stringsAsFactors = FALSE, row.names = NULL)
  })
  results <- do.call(rbind, rows)
  if (is.null(results)) results <- data.frame()
  global <- data.frame(effect = c("shift", "var", "total"), p = NA_real_, stringsAsFactors = FALSE)
  if (nrow(results)) {
    for (eff in c("shift", "var", "total")) results[[paste0("padj.", eff)]] <- stats::p.adjust(results[[paste0("p.", eff)]], "BH")
    # max-statistic combination over cell types from the coupled permutations
    for (eff in c("shift", "var", "total")) {
      col <- switch(eff, shift = "F", var = "var.t", eff)   # the studentized dispersion change, as its p-value
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
      if (isFALSE(r$disp.identifiable)) f <- c(f, "dispersion-not-identifiable")
      if (is.finite(r$cell.assoc) && r$cell.assoc > 0.5) f <- c(f, "cells-tied-to-contrast")   # var test anticonservative
      paste(f, collapse = ";")
    }, character(1))
    results <- results[order(results$shift, decreasing = TRUE), ]
    rownames(results) <- NULL
  }
  notes <- unique(unlist(c(gplan$notes, lapply(cts, function(ct) { nt <- setdiff(runs[[ct]]$plan$notes, gplan$notes); if (length(nt)) sprintf("%s: %s", ct, nt) }))))
  if (verbose && length(notes)) message(paste(notes, collapse = "\n"))
  list(results = results, global = global, fits = fits, skipped = skipped, plan = gplan, notes = notes,
       perm.stats = if (return.perm.stats) lapply(runs, `[[`, "perm") else NULL,
       call.info = list(dist = dist, permutation = permutation, n.permutations = n.permutations, seed = seed, robust = robust, na.mode = na.mode, robust.k = robust.k,
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
