## Expression shifts on the pairwise-effects engine (Track B.5): pseudobulk per cell type, distance
## matrices, per-test inference (contrast tests via testPairwiseEffects(), term tests via
## testTermEffects()), the long results table of the API note (§4.7), adjusted distances and influence.

# ---- pseudobulk ----------------------------------------------------------------------------------

#' Pseudobulk expression per cell type and sample
#'
#' Sums counts per (cell type, sample), drops (cell type, sample) cells with fewer than
#' `min.cells.per.sample` cells, keeps genes expressed (count > 1) in at least `min.gene.frac` of the cells
#' of at least 10% of the samples, and transforms to `log10(1e3 * count / total + 1)`.
#'
#' @param cms list (per sample) of cell x gene count matrices
#' @param cell.groups named factor of cell type per cell
#' @param sample.per.cell named factor of sample per cell
#' @param min.cells.per.sample minimum cells for a (cell type, sample) entry (default 10)
#' @param min.gene.frac gene filter (see above)
#' @param genes optional explicit gene set (skips the gene filter)
#' @return list: `cm.per.type` (per cell type: sample x gene matrix, absent samples as `NA` rows),
#'   `n.cells` (cell type x sample table of cell counts), `genes`
#' @export
pseudobulkPerCellType <- function(cms, cell.groups, sample.per.cell, min.cells.per.sample = 10, min.gene.frac = 0.01, genes = NULL) {
  stopifnot(is.list(cms), length(cms) > 0)
  cell.groups <- droplevels(as.factor(cell.groups))
  all.samples <- names(cms)
  cms <- lapply(cms, function(cm) cm[intersect(rownames(cm), names(cell.groups)), , drop = FALSE])
  if (is.null(genes)) {
    genes <- lapply(cms, function(cm) {
      if (!nrow(cm)) return(character(0))
      cm@x <- 1 * (cm@x > 1)
      names(which(Matrix::colMeans(cm) > min.gene.frac))
    }) %>% unlist() %>% table() %>% {. / length(cms) > 0.1} %>% which() %>% names()
  }
  if (!length(genes)) stop("no genes pass the expression filter")
  collapsed <- lapply(cms, function(cm) {
    if (!nrow(cm)) return(NULL)
    x <- sccore::collapseCellsByType(cm, groups = cell.groups, min.cell.count = min.cells.per.sample)
    x <- sccore::extendMatrix(x, genes)[, genes, drop = FALSE]
    as.matrix(x)
  })
  n.cells <- table(cell.groups, factor(as.character(sample.per.cell[names(cell.groups)]), levels = all.samples))
  types <- levels(cell.groups)
  cm.per.type <- lapply(sccore::sn(types), function(ct) {
    mat <- matrix(NA_real_, length(all.samples), length(genes), dimnames = list(all.samples, genes))
    for (s in all.samples) { x <- collapsed[[s]]; if (!is.null(x) && ct %in% rownames(x)) mat[s, ] <- x[ct, ] }
    tot <- pmax(1, rowSums(mat, na.rm = TRUE)); tot[rowSums(!is.na(mat)) == 0] <- NA
    log10(mat / tot * 1e3 + 1)
  })
  list(cm.per.type = cm.per.type, n.cells = n.cells, genes = genes)
}

# ---- term tests (K-level factors) ------------------------------------------------------------------

# Term-level effects for one cell type: location F / R2 for the term's columns (full vs reduced design),
# dispersion F for the term in the dispersion model, partial R^2 per term and the model-implied pairwise
# table over the term's levels.
termEffectsFromDesign <- function(D, design, meta, dispersion.formula = NULL, dist = c("cor", "l2", "l1"), min.samp.per.level = 3,
                                  robust = "none", na.mode = "drop", robust.k = 1.345) {
  dist <- match.arg(dist)
  skip <- function(reason) list(ok = FALSE, reason = reason)
  samples <- intersect(rownames(design$F), rownames(D))
  if (length(samples) < 3) return(skip("fewer than 3 samples with this cell type"))
  v <- design$term.variable
  present <- samples; w <- NULL
  D2 <- asSquaredDistance(D[samples, samples], dist)
  if (identical(na.mode, "impute_weak")) {
    all.s <- intersect(rownames(design$F), rownames(meta))
    if (length(samples) < length(all.s)) { imp <- imputeAbsentSamples(D2, samples, all.s); D2 <- imp$D2; w <- imp$w; samples <- all.s }
  }
  weighted <- !is.null(w) || !identical(robust, "none")
  G <- gowerCenter(D2)
  F <- design$F[samples, , drop = FALSE]
  keep <- colSums(abs(F)) > 1e-12
  Xf <- F[, keep, drop = FALSE]
  tcols <- intersect(design$term.cols, colnames(Xf))
  if (!length(tcols)) return(skip("no level of the tested term varies in this cell type"))
  Xr <- Xf[, setdiff(colnames(Xf), tcols), drop = FALSE]
  if (!ncol(Xr)) Xr <- matrix(1, nrow(Xf), 1, dimnames = list(rownames(Xf), "(Intercept)"))
  if (qr(Xf)$rank >= nrow(Xf)) return(skip("no residual degrees of freedom"))
  g <- as.character(meta[present, v]); tb <- table(g)
  if (sum(tb >= min.samp.per.level) < 2) return(skip(sprintf("fewer than 2 levels of %s with at least %d samples", v, min.samp.per.level)))
  if (is.null(dispersion.formula)) dispersion.formula <- stats::as.formula(paste("~", v))
  Zf.full <- buildFullDesign(dispersion.formula, meta[samples, , drop = FALSE])
  zk <- colSums(abs(Zf.full)) > 1e-12; Zf <- Zf.full[, zk, drop = FALSE]
  zassign <- attr(Zf.full, "assign")[zk]; ztl <- attr(attr(Zf.full, "terms"), "term.labels")
  zin <- vapply(ztl, function(t) v %in% all.vars(str2lang(t)), logical(1))
  zcols <- colnames(Zf)[zassign %in% which(zin)]
  Zr <- Zf[, setdiff(colnames(Zf), zcols), drop = FALSE]
  if (!ncol(Zr)) Zr <- matrix(1, nrow(Zf), 1, dimnames = list(rownames(Zf), "(Intercept)"))
  loc <- termTestGower(G, Xf, Xr)
  disp <- if (length(zcols)) dispersionTermTest(G, Xf, Zf, Zr) else c(F.disp = NA_real_, p.disp = NA_real_, df.disp = 0, nu.disp = NA_real_)
  base.w <- w %||% rep(1, nrow(Xf)); w.final <- base.w
  if (weighted) {                                            # weighted / robust statistics (R2 and the analytic p stay unweighted)
    wt <- weightedTermStats(G, Xf, Xr, Zf, Zr, w = base.w, robust = robust, k = robust.k)
    loc["F"] <- wt$F; loc["p.analytic"] <- NA_real_; if (length(zcols)) disp["F.disp"] <- wt$F.disp; w.final <- wt$w
  }
  names(base.w) <- names(w.final) <- rownames(Xf)
  disp.r2 <- if (length(zcols)) dispersionR2(G, Xf, Zf, Zr) else NA_real_
  # pairwise table over levels: shift_ab = (x_a - x_b)' M (x_a - x_b), with M bias-corrected
  hi <- hatInfo(Xf); n <- nrow(Xf); R <- diag(n) - hi$H
  r <- rowSums((R %*% G) * R)
  gamma <- fitDispersion(R, Zf, r); s <- drop(Zf %*% gamma)
  M <- hi$A %*% G %*% t(hi$A) - hi$A %*% (s * t(hi$A))
  levs <- names(tb)[tb >= 1]
  one <- function(Fm, at) .oneRowFromFormula(attr(Fm, "terms"), attr(Fm, "xlevels"), attr(Fm, "contrasts"), colnames(Fm),
                                             meta[samples, , drop = FALSE], at, design$numeric_ref_used)
  Xrows <- t(sapply(levs, function(l) one(design$F, setNames(list(l), v))[colnames(F)][keep]))
  Zrows <- t(sapply(levs, function(l) one(Zf.full, setNames(list(l), v))[zk]))
  if (length(levs) == 1) { Xrows <- matrix(Xrows, 1); Zrows <- matrix(Zrows, 1) }
  rownames(Xrows) <- rownames(Zrows) <- levs
  cell.table <- impliedCellTable(M, gamma, Xrows, Zrows)
  s.lev <- drop(Zrows %*% gamma); names(s.lev) <- levs
  pairs <- utils::combn(levs, 2)
  pair.table <- data.frame(a = pairs[1, ], b = pairs[2, ],
                           shift = apply(pairs, 2, function(ab) { d <- Xrows[ab[1], ] - Xrows[ab[2], ]; drop(crossprod(d, M %*% d)) }),
                           var = apply(pairs, 2, function(ab) 2 * (s.lev[ab[1]] - s.lev[ab[2]])), stringsAsFactors = FALSE)
  pair.table$total <- pair.table$shift + pair.table$var / 2
  list(ok = TRUE, kind = "term", variable = v, G = G, X = Xf, Xr = Xr, Z = Zf, Zr = Zr, samples = samples, n = n,
       location = loc, dispersion = disp, R2.disp.adj = disp.r2, cell.table = cell.table, pair.table = pair.table, s.level = s.lev,
       gamma = gamma, M = M, n.per.level = tb, dispersion.formula = dispersion.formula,
       weighted = weighted, base.w = base.w, w = w.final, robust = robust, robust.k = robust.k, na.mode = na.mode, present = present)
}

# chance-corrected R^2 of the dispersion regression sqrt(v) ~ Zf vs Zr
dispersionR2 <- function(G, Xloc, Zf, Zr) {
  n <- nrow(G); hi <- hatInfo(Xloc); R <- diag(n) - hi$H
  r <- rowSums((R %*% G) * R); ok <- hi$h < 0.99
  y <- sqrt(pmax(r[ok], 0) / (1 - hi$h[ok])); m <- sum(ok)
  Zf <- as.matrix(Zf)[ok, , drop = FALSE]; Zr <- as.matrix(Zr)[ok, , drop = FALSE]
  rss <- function(Z) { H <- hatInfo(Z)$H; sum(((diag(m) - H) %*% y)^2) }
  df <- qr(Zf)$rank - qr(Zr)$rank; nu <- m - qr(Zf)$rank
  if (df <= 0 || nu <= 0) return(NA_real_)
  rf <- rss(Zf); rr <- rss(Zr)
  (rr - rf - df * rf / nu) / (rr + rf / nu)
}

# hat matrices, ranks and degrees of freedom of a term test, as the C++ kernel needs them
termKernelInputs <- function(Xf, Xr, Zf, Zr, n) {
  hf <- hatInfo(Xf); hr <- hatInfo(Xr)
  Zf <- as.matrix(Zf); Zr <- as.matrix(Zr)
  list(Hf = hf$H, Hr = hr$H, Zf = Zf, Zr = Zr, df = hf$rank - hr$rank, nu = n - hf$rank, qZf = qr(Zf)$rank, qZr = qr(Zr)$rank)
}

# permutation statistics for a term test: F (location) and F.disp under relabelings
termPermutationStats <- function(eff, plan, P) {
  G <- eff$G; Xf <- eff$X; Xr <- eff$Xr; Zf <- eff$Z; Zr <- eff$Zr
  obs <- c(F = unname(eff$location["F"]), F.disp = unname(eff$dispersion["F.disp"]))
  B <- ncol(P); perm <- matrix(NA_real_, B, 2, dimnames = list(NULL, c("F", "F.disp")))
  if (B == 0) return(list(obs = obs, perm = perm))
  need.disp <- is.finite(obs["F.disp"])
  storage.mode(P) <- "integer"
  if (isTRUE(eff$weighted)) {                   # weighted / robust path (C++; R reference: weightedTermStats())
    code <- c(none = 0L, huber = 1L, winsor = 2L)[[eff$robust]]
    perm[, ] <- permuted_term_stats_w(G, Xf, Xr, Zf, Zr, eff$base.w, P, plan$scheme == "freedman-lane", eff$w, need.disp, code, eff$robust.k, 5L)
    return(list(obs = obs, perm = perm))
  }
  # C++ kernel (R references: termTestGower() / dispersionTermTest(), looped per relabeling)
  k <- termKernelInputs(Xf, Xr, Zf, Zr, nrow(G))
  perm[, ] <- if (plan$scheme == "block" || plan$scheme == "huh-jhun") {
    permuted_term_stats(G, k$Hf, k$Hr, k$Zf, k$Zr, k$df, k$nu, k$qZf, k$qZr, P, need.disp)
  } else {   # freedman-lane on the reduced location model
    parts <- flGowerParts(G, Xr)
    permuted_term_stats_fl(parts$K1, parts$K2, parts$K3, parts$K4, k$Hf, k$Hr, k$Zf, k$Zr, k$df, k$nu, k$qZf, k$qZr, P, need.disp)
  }
  list(obs = obs, perm = perm)
}

#' Permutation tests of term-level effects (K-level factor) for a set of cell types
#'
#' Companion of [testPairwiseEffects()] for tests of a whole factor: location (PERMANOVA-type F with
#' chance-corrected R^2) and dispersion (F of the per-sample dispersion on the factor), with permutations
#' shared across cell types for the max-statistic combination.
#'
#' @inheritParams testPairwiseEffects
#' @param dispersion.formula dispersion formula (default: the tested variable)
#' @param dist distance type of the matrices
#' @param min.samp.per.level minimum samples per level of the tested factor
#' @return list like [testPairwiseEffects()]: `results` (one row per cell type: `R2`, `R2.adj`, `R2.partial`,
#'   `F`, `df`, `r.eff`, `p.location`, `F.disp`, `R2.disp.adj`, `p.dispersion`, adjusted p-values, ...),
#'   `global`, `fits`, `skipped`, `plan`, `notes`, `call.info`
#' @export
testTermEffects <- function(D.list, design, meta, dispersion.formula = NULL, dist = c("cor", "l2", "l1"),
                            permutation = c("auto", "block", "freedman-lane", "huh-jhun"), n.permutations = 999,
                            block.vars = NULL, min.samp.per.level = 3, seed = NULL, alpha = 0.05, n.cores = 1,
                            return.perm.stats = FALSE, verbose = FALSE, robust = c("none", "huber", "winsor"), na.mode = c("drop", "impute_weak"), robust.k = 1.345) {
  dist <- match.arg(dist); permutation <- match.arg(permutation); robust <- match.arg(robust); na.mode <- match.arg(na.mode)
  if (permutation == "huh-jhun") { permutation <- "block"; warning("huh-jhun is a contrast-only scheme; using block permutations for the term test") }
  if (is.null(names(D.list))) names(D.list) <- paste0("CT", seq_along(D.list))
  all.samples <- intersect(rownames(design$F), rownames(meta))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  gplan <- permutationPlan(design, meta, all.samples, scheme = permutation, block.vars = block.vars, n.permutations = n.permutations)   # enumerated when small
  P.global <- withSeed(seed, drawPermutations(gplan, n.permutations))
  fits <- lapply(D.list, function(D) termEffectsFromDesign(D, design, meta, dispersion.formula = dispersion.formula, dist = dist,
                                                          min.samp.per.level = min.samp.per.level, robust = robust, na.mode = na.mode, robust.k = robust.k))
  ok <- vapply(fits, function(f) isTRUE(f$ok), logical(1))
  skipped <- data.frame(celltype = names(fits)[!ok], reason = vapply(fits[!ok], function(f) f$reason, character(1)), stringsAsFactors = FALSE, row.names = NULL)
  runOne <- function(ct) {
    eff <- fits[[ct]]
    plan <- permutationPlan(design, meta, eff$samples, scheme = gplan$scheme, block.vars = block.vars, n.permutations = n.permutations)
    sub <- match(eff$samples, all.samples)
    P.ind <- apply(P.global, 2, inducePermutation, sub = sub, sub.plan = plan)
    if (!is.matrix(P.ind)) P.ind <- matrix(P.ind, ncol = ncol(P.global))          # a single permutation
    st.c <- termPermutationStats(eff, plan, P.ind)
    st.p <- if (plan$exhaustive) termPermutationStats(eff, plan, drawPermutations(plan)) else st.c
    list(plan = plan, obs = st.c$obs, perm.coupled = st.c$perm, perm = st.p$perm)
  }
  cts <- names(fits)[ok]
  runs <- if (n.cores > 1 && length(cts) > 1) sccore::plapply(cts, runOne, n.cores = n.cores, progress = verbose, fail.on.error = TRUE) else lapply(cts, runOne)
  names(runs) <- cts
  rows <- lapply(cts, function(ct) {
    r <- runs[[ct]]; e <- fits[[ct]]; ex <- r$plan$exhaustive; L <- e$location
    data.frame(celltype = ct, R2 = unname(L["R2"]), R2.adj = unname(L["R2.adj"]), R2.partial = unname(L["R2.partial"]),
               F = unname(L["F"]), df = unname(L["df"]), nu = unname(L["nu"]), r.eff = unname(L["r.eff"]),
               p.location = permutationPValue(r$obs["F"], r$perm[, "F"], "greater", ex),
               F.disp = unname(e$dispersion["F.disp"]), R2.disp.adj = e$R2.disp.adj,
               p.dispersion = permutationPValue(r$obs["F.disp"], r$perm[, "F.disp"], "greater", ex),
               n = e$n, scheme = r$plan$scheme, n.perm = nrow(r$perm), n.perm.distinct = r$plan$n.distinct, exhaustive = ex,
               p.floor = r$plan$p.floor, stringsAsFactors = FALSE, row.names = NULL)
  })
  results <- do.call(rbind, rows); if (is.null(results)) results <- data.frame()
  global <- data.frame(effect = c("location", "dispersion"), p = NA_real_, stringsAsFactors = FALSE)
  if (nrow(results)) {
    results$padj.location <- stats::p.adjust(results$p.location, "BH"); results$padj.dispersion <- stats::p.adjust(results$p.dispersion, "BH")
    for (eff in c("location", "dispersion")) {
      col <- if (eff == "location") "F" else "F.disp"
      obs <- vapply(runs, function(r) unname(r$obs[col]), numeric(1))
      perm <- sapply(runs, function(r) r$perm.coupled[, col]); if (!is.matrix(perm)) perm <- matrix(perm, ncol = length(runs))
      mu <- colMeans(perm, na.rm = TRUE); sdv <- apply(perm, 2, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- NA
      z.obs <- (obs - mu) / sdv; z.perm <- sweep(sweep(perm, 2, mu), 2, sdv, "/")
      usable <- is.finite(z.obs)
      if (any(usable)) {
        mx <- apply(z.perm[, usable, drop = FALSE], 1, max, na.rm = TRUE)
        add1 <- as.numeric(!gplan$exhaustive)
        global$p[global$effect == eff] <- (sum(mx >= max(z.obs[usable]) - 1e-12) + add1) / (length(mx) + add1)
        results[[paste0("p.fwer.", eff)]] <- ifelse(usable, (vapply(z.obs, function(z) sum(mx >= z - 1e-12), numeric(1)) + add1) / (length(mx) + add1), NA_real_)
      } else results[[paste0("p.fwer.", eff)]] <- NA_real_
    }
    results$flags <- vapply(seq_len(nrow(results)), function(i) {
      f <- character(0); r <- results[i, ]
      if (r$p.floor > alpha / 10) f <- c(f, "few-permutations")
      if (r$scheme == "freedman-lane") f <- c(f, "approximate-scheme")
      if (is.finite(r$r.eff) && r$scheme == "freedman-lane" && r$r.eff > r$n) f <- c(f, "large-effective-dimension")
      paste(f, collapse = ";")
    }, character(1))
    results <- results[order(results$R2.adj, decreasing = TRUE), ]; rownames(results) <- NULL
  }
  notes <- unique(unlist(c(gplan$notes, lapply(runs, function(r) r$plan$notes))))
  list(results = results, global = global, fits = fits, skipped = skipped, plan = gplan, notes = notes,
       perm.stats = if (return.perm.stats) lapply(runs, `[[`, "perm") else NULL,
       call.info = list(dist = dist, permutation = permutation, n.permutations = n.permutations, seed = seed,
                        dispersion.formula = dispersion.formula, block.vars = block.vars))
}

# ---- adjusted distances and influence ------------------------------------------------------------

# Sample distances with the nuisance part of the design removed: G_adj = R_Z G R_Z (a proper Gram matrix),
# returned on the scale of the input distance (`cor`/`l1`: squared form as is, `l2`: square root).
adjustedDistanceMatrix <- function(G, Z, dist) {
  n <- nrow(G)
  Ga <- if (is.null(Z) || !ncol(Z)) G else { R <- diag(n) - hatInfo(Z)$H; R %*% G %*% R }
  D2 <- uncenterGower(Ga)
  if (dist == "l2") sqrt(D2) else D2
}

nuisanceColumns <- function(eff, design) {
  if (identical(eff$kind, "term")) return(eff$Xr)
  X <- eff$X; cvec <- eff$contrast
  N <- contrastNullBasis(cvec)
  if (is.null(N) || !ncol(N)) return(NULL)
  Zn <- X %*% N
  Zn[, colSums(abs(Zn)) > 1e-12, drop = FALSE]
}

# ---- results assembly ----------------------------------------------------------------------------

effectRowsContrast <- function(res, test.label, test.id) {
  w <- res$results
  if (!nrow(w)) return(NULL)
  z <- stats::qnorm(0.975)
  one <- function(eff) {
    stat <- if (eff == "shift") w$F else w[[eff]]
    se <- w[[paste0("se.", eff)]]
    data.frame(test = test.label, test.id = test.id, celltype = w$celltype, effect = eff,
               estimate = w[[eff]], estimate.norm = w[[paste0(eff, ".norm")]], se.jk = se,
               ci.low = w[[eff]] - z * se, ci.high = w[[eff]] + z * se,
               statistic = stat, p = w[[paste0("p.", eff)]], padj = w[[paste0("padj.", eff)]], p.fwer = w[[paste0("p.fwer.", eff)]],
               n = w$n, n.ref = w$n.ref, n.alt = w$n.alt, scheme = w$scheme, n.perm = w$n.perm, n.perm.distinct = w$n.perm.distinct,
               p.floor = w$p.floor, flags = w$flags, stringsAsFactors = FALSE)
  }
  do.call(rbind, lapply(c("shift", "var", "total"), one))
}

effectRowsTerm <- function(res, test.label, test.id) {
  w <- res$results
  if (!nrow(w)) return(NULL)
  mk <- function(eff, est, stat, p, padj, pf) data.frame(test = test.label, test.id = test.id, celltype = w$celltype, effect = eff,
    estimate = est, estimate.norm = est, se.jk = NA_real_, ci.low = NA_real_, ci.high = NA_real_, statistic = stat, p = p, padj = padj, p.fwer = pf,
    n = w$n, n.ref = NA_integer_, n.alt = NA_integer_, scheme = w$scheme, n.perm = w$n.perm, n.perm.distinct = w$n.perm.distinct,
    p.floor = w$p.floor, flags = w$flags, stringsAsFactors = FALSE)
  rbind(mk("location", w$R2.adj, w$F, w$p.location, w$padj.location, w$p.fwer.location),
        mk("dispersion", w$R2.disp.adj, w$F.disp, w$p.dispersion, w$padj.dispersion, w$p.fwer.dispersion))
}

#' Expression-shift tests for every test of a model
#'
#' Runs [testPairwiseEffects()] (contrast tests) or [testTermEffects()] (term tests) per test of a
#' `cacoaModel` on per-cell-type sample distance matrices and assembles the long results table.
#'
#' @param D.list named list of sample x sample distance matrices (one per cell type)
#' @param model a `cacoaModel` (see [buildCacoaModel()])
#' @param dist distance type of the matrices
#' @param permutation,n.permutations,block.vars,bias.correct,influence,min.samp.per.level,seed,alpha,n.cores,verbose
#'   passed to the engine
#' @param n.cells optional cell type x sample table of cell counts (attached to the distance matrices)
#' @return list: `results` (long table: one row per test x cell type x effect), `wide` (per test, the engine's
#'   table), `global` (per test x effect), `fits`, `skipped`, `adjusted.distances` (per test, per cell type),
#'   `influence` (per test: samples x cell types matrix of the change in shift when the sample is left out),
#'   `distances`, `model`, `notes`, `settings`
#' @export
expressionShiftsForModel <- function(D.list, model, dist = "cor", permutation = "auto", n.permutations = 999, block.vars = model$block.vars,
                                     bias.correct = TRUE, influence = TRUE, min.samp.per.level = 3, seed = NULL, alpha = 0.05,
                                     n.cores = 1, verbose = FALSE, n.cells = NULL, robust = "none", na.mode = "drop", robust.k = 1.345) {
  meta <- model$meta
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  per.test <- lapply(model$tests, function(t) {
    if (t$kind == "contrast") {
      r <- testPairwiseEffects(D.list, t$design, meta, dispersion.formula = model$dispersion.formula, dist = dist, permutation = permutation,
                               n.permutations = n.permutations, block.vars = block.vars, bias.correct = bias.correct, influence = influence,
                               min.samp.per.level = min.samp.per.level, seed = seed, alpha = alpha, n.cores = n.cores, verbose = verbose,
                               robust = robust, na.mode = na.mode, robust.k = robust.k)
      r$rows <- effectRowsContrast(r, t$label, t$id)
    } else {
      r <- testTermEffects(D.list, t$design, meta, dispersion.formula = model$dispersion.formula, dist = dist, permutation = permutation,
                           n.permutations = n.permutations, block.vars = block.vars, min.samp.per.level = min.samp.per.level, seed = seed,
                           alpha = alpha, n.cores = n.cores, verbose = verbose, robust = robust, na.mode = na.mode, robust.k = robust.k)
      r$rows <- effectRowsTerm(r, t$label, t$id)
    }
    r$global$test <- t$label; r$global$test.id <- t$id
    # adjusted distances and influence per cell type
    r$adjusted.distances <- Map(function(e, ct) if (isTRUE(e$ok)) {
      M <- adjustedDistanceMatrix(e$G, nuisanceColumns(e, t$design), dist)
      if (!is.null(n.cells) && ct %in% rownames(n.cells)) attr(M, "n.cells") <- setNames(as.numeric(n.cells[ct, rownames(M)]), rownames(M))
      M } else NULL, r$fits, names(r$fits))
    r$adjusted.distances <- Filter(Negate(is.null), r$adjusted.distances)
    if (t$kind == "contrast" && influence) {
      ok <- Filter(function(e) isTRUE(e$ok) && !is.null(e$influence), r$fits)
      if (length(ok)) {
        samples <- rownames(meta)
        infl <- lapply(c("shift", "var", "total"), function(eff) {
          m <- matrix(NA_real_, length(samples), length(ok), dimnames = list(samples, names(ok)))
          for (ct in names(ok)) { e <- ok[[ct]]; m[rownames(e$influence$effects), ct] <- e$influence$effects[, eff] - e[[eff]] }
          m
        })
        names(infl) <- c("shift", "var", "total"); r$influence <- infl
      }
    }
    r
  })
  names(per.test) <- vapply(model$tests, `[[`, character(1), "label")
  results <- do.call(rbind, Filter(Negate(is.null), lapply(per.test, `[[`, "rows")))
  if (is.null(results)) results <- data.frame()
  rownames(results) <- NULL
  global <- do.call(rbind, lapply(per.test, `[[`, "global")); rownames(global) <- NULL
  skipped <- do.call(rbind, Map(function(r, nm) if (nrow(r$skipped)) cbind(test = nm, r$skipped) else NULL, per.test, names(per.test)))
  if (is.null(skipped)) skipped <- data.frame(test = character(0), celltype = character(0), reason = character(0))
  rownames(skipped) <- NULL
  list(results = results, wide = lapply(per.test, `[[`, "results"), global = global, fits = lapply(per.test, `[[`, "fits"),
       skipped = skipped, adjusted.distances = lapply(per.test, `[[`, "adjusted.distances"), influence = lapply(per.test, `[[`, "influence"),
       distances = D.list, model = model, notes = unique(unlist(lapply(per.test, `[[`, "notes"))),
       settings = list(dist = dist, permutation = permutation, n.permutations = n.permutations, seed = seed, block.vars = block.vars,
                       bias.correct = bias.correct, influence = influence, min.samp.per.level = min.samp.per.level, alpha = alpha,
                       robust = robust, na.mode = na.mode, robust.k = robust.k))
}

# short settings key for caches
settingsKey <- function(...) {
  x <- list(...)
  if (requireNamespace("rlang", quietly = TRUE) && exists("hash", asNamespace("rlang"))) rlang::hash(x) else paste(deparse(x), collapse = "")
}
