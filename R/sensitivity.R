## Sensitivity of expression-shift results (Track D.1): alternative covariate sets (D31) and single-sample
## influence, summarized as a plain verdict per cell type.

# Candidate model formulas for a test: unadjusted, current, minus each adjusting covariate, plus each of the
# top-k screened covariates not in the model, all screened (if the df budget allows). The grouping terms `(1 | g)` of the model are kept
# in every formula (they set how the samples are relabeled); each one is also dropped on its own ("minus (1 | g)").
# key of a formula for spotting duplicates: its variables, the grouping variables marked as such
formulaKey <- function(f) { sp <- splitGroupingTerms(f); paste(c(sort(all.vars(sp$formula)), sort(paste0("|", sp$groups))), collapse = "+") }

sensitivityFormulas <- function(model, screen = NULL, meta, top.k = 3, min.resid.df = 5) {
  t <- model$tests[[1]]; tv <- t$variable
  adj <- setdiff(all.vars(model$formula), tv); grp <- model$random.vars %||% character(0)
  n <- length(model$samples$used)
  fml <- function(vars, g = grp) stats::as.formula(paste("~", paste(c(vars, if (length(g)) sprintf("(1 | %s)", g)), collapse = " + ")))
  out <- list()
  out[["unadjusted"]] <- fml(tv)
  out[["current"]] <- model$formula.full %||% model$formula
  for (v in adj) out[[paste("minus", v)]] <- fml(c(tv, setdiff(adj, v)))
  for (g in grp) out[[sprintf("minus (1 | %s)", g)]] <- fml(c(tv, adj), setdiff(grp, g))
  if (!is.null(screen)) {
    sm <- if ("partial" %in% screen$modes) "partial" else screen$modes[1]
    g <- screen$global[screen$global$mode == sm, ]
    g <- g[!g$covariate %in% c(tv, adj, grp) & is.finite(g$p.global), ]
    g <- g[order(-g$n.sig, g$p.global), ]
    cand <- g$covariate[g$n.sig > 0]
    assoc <- vapply(cand, function(v) { a <- tryCatch(covariateAssociation(meta[[tv]], meta[[all.vars(str2lang(v))[1]]]), error = function(e) NA_real_); if (is.na(a)) 0 else a }, numeric(1))
    cand <- utils::head(cand[assoc < 0.95], top.k)       # a covariate aliased with the test variable cannot be adjusted for
    for (v in cand) out[[paste("plus", v)]] <- fml(c(tv, adj, v))
    if (length(cand) > 1) {
      f.all <- fml(c(tv, adj, cand))
      q <- tryCatch(qr(buildFullDesign(stats::as.formula(paste("~", paste(c(tv, adj, cand), collapse = " + "))), meta))$rank, error = function(e) Inf)
      if (n - q >= min.resid.df) out[["all screened"]] <- f.all
    }
  }
  keys <- vapply(out, formulaKey, character(1))
  out <- out[setdiff(names(out), names(out)[names(out) != "current" & keys == keys[["current"]]])]   # "current" wins over duplicates
  keys <- keys[names(out)]
  out[!duplicated(keys)]
}

#' Sensitivity of expression-shift effects to the covariate set and to single samples
#'
#' @param D.list named list of sample distance matrices (per cell type)
#' @param model the current `cacoaModel` (its first test is examined)
#' @param meta sample metadata
#' @param formulas named list of alternative location formulas (default: the D31 sets of `sensitivityFormulas()`)
#' @param screen optional `cacoaCovariateScreen` used to propose added covariates
#' @param influence optional per-test influence matrices from the current result (`res$influence[[1]]`) and the
#'   current wide table (`res$wide[[1]]`), as `list(influence =, wide =)`; when `NULL` they are recomputed
#' @param dist,permutation,n.permutations,seed,alpha,min.samp.per.level,n.cores engine settings
#' @param top.k number of screened covariates to try adding
#' @param current.rows optional long results table of the stored model (`res$results`), used for the "current" model instead of a refit
#' @param rel.change relative change of the shift estimate across models that counts as "sensitive" (default 0.5)
#' @param influence.rel relative change of the shift estimate when one sample is left out that, together with being
#'   an outlier among the leave-one-out changes (> 3 MAD), marks the sample as driving the result (default 0.2)
#' @param significance optional function of a fitted model (`cacoaModel`, one test) returning a data.frame with `celltype`,
#'   `p` and `padj` of the shift: the test that decides significance when it is not the distance model's own (e.g. gsd, see
#'   `Cacoa$checkSensitivity()`); the distance model's p-values are then kept as `p.distance` / `padj.distance`.
#'   Rows of `current.rows` are taken as they are.
#' @param significance.label name of that test, for the printout and the plot (default "distance test")
#' @return object of class `cacoaSensitivity`: `table` (model x cell type x effect: estimate, ci, p, padj),
#'   `summary` (per cell type: same sign in all models, n significant of m, max relative change of the normalized shift, most
#'   influential sample, `verdict.magnitude` (the size of the shift: near zero (the current shift lies within the relabeling
#'   noise band of the main result, `influence$wide`), stable, changes by x%, sign depends on the model, driven by one sample), `verdict.significance` (significant in all / k of m / no models) and the combined `verdict`),
#'   `formulas`, `settings`
#' @export
#' @param robust,na.mode,robust.k robust fit (`"none"`, `"huber"`, `"winsor"`), treatment of absent samples (`"drop"`, `"impute_weak"`) and
#'   robust tuning constant, passed to [expressionShiftsForModel()] (default: the settings of the result being checked)
checkSensitivity <- function(D.list, model, meta, formulas = NULL, screen = NULL, influence = NULL, dist = "cor", permutation = "auto",
                             n.permutations = 199, seed = 1, alpha = 0.05, min.samp.per.level = 3, n.cores = 1, top.k = 3, rel.change = 0.5,
                             influence.rel = 0.2, robust = "none", na.mode = "drop", robust.k = 1.345, current.rows = NULL,
                             significance = NULL, significance.label = "distance test", jackknife = FALSE) {
  t <- model$tests[[1]]
  if (t$kind != "contrast") stop("sensitivity analysis is defined for contrast tests (two groups or a numeric step)")
  if (is.null(formulas)) formulas <- sensitivityFormulas(model, screen, meta, top.k = top.k)
  if (!"current" %in% names(formulas)) {            # the stored model is "current" whatever the caller named it
    cur <- model$formula.full %||% model$formula
    same <- which(vapply(formulas, formulaKey, character(1)) == formulaKey(cur))
    if (length(same)) names(formulas)[same[1]] <- "current" else formulas <- c(list(current = cur), formulas)
  }
  spec <- model$test.spec %||% list(test = NULL, contrast = t$contrast)
  runs <- lapply(names(formulas), function(nm) {
    m <- tryCatch(buildCacoaModel(meta, formula = formulas[[nm]], test = spec$test, contrast = spec$contrast, dispersion.formula = model$dispersion.formula,
                                  block.vars = setdiff(model$block.vars, model$random.vars), numeric.ref = model$numeric.ref, numeric.step = model$numeric.step,
                                  permutation = permutation, n.permutations = n.permutations), error = function(e) NULL)
    if (is.null(m) || any(m$issues$severity == "error")) return(NULL)
    m$tests <- m$tests[1]
    r <- expressionShiftsForModel(D.list, m, dist = dist, permutation = permutation, n.permutations = n.permutations, block.vars = m$block.vars,
                                  influence = jackknife || (nm == "current" && is.null(influence)), min.samp.per.level = min.samp.per.level, seed = seed, alpha = alpha,
                                  n.cores = n.cores, robust = robust, na.mode = na.mode, robust.k = robust.k)
    if (!nrow(r$results)) return(NULL)
    rows <- r$results; rows$model <- nm; rows$formula <- paste(deparse(formulas[[nm]]), collapse = "")
    rows <- addLeaveOneOutRange(rows, r$influence[[1]])   # the range of the estimate when one sample is left out (jackknife = TRUE: every model)
    list(rows = rows, res = r, model = m)
  })
  names(runs) <- names(formulas)
  runs <- Filter(Negate(is.null), runs)
  if (!is.null(current.rows) && "current" %in% names(runs) && setequal(current.rows$celltype, runs[["current"]]$rows$celltype)) {
    rows <- current.rows[, intersect(c(names(runs[["current"]]$rows), "p.distance", "padj.distance"), names(current.rows))]   # the main result's estimates and p-values
    rows$model <- "current"; rows$formula <- paste(deparse(formulas[["current"]]), collapse = "")
    rows <- addLeaveOneOutRange(rows, (influence %||% list(influence = runs[["current"]]$res$influence[[1]]))$influence)
    runs[["current"]]$rows <- rows; runs[["current"]]$stored <- TRUE
  }
  # significance of the shift by another test (the one of the stored result): refit it under every model not taken from the stored result
  if (!is.null(significance)) for (nm in names(runs)) {
    rows <- runs[[nm]]$rows
    if (!isTRUE(runs[[nm]]$stored)) {
      sg <- tryCatch(significance(runs[[nm]]$model), error = function(e) NULL)
      rows$p.distance <- rows$p; rows$padj.distance <- rows$padj
      sh <- rows$effect == "shift"; k <- if (is.null(sg)) rep(NA_integer_, nrow(rows)) else match(rows$celltype, sg$celltype)
      rows$p[sh] <- if (is.null(sg)) NA_real_ else sg$p[k[sh]]; rows$padj[sh] <- if (is.null(sg)) NA_real_ else sg$padj[k[sh]]
    }
    runs[[nm]]$rows <- rows
  }
  dropped <- setdiff(names(formulas), names(runs))
  if (!length(runs)) stop("no alternative model could be fitted")
  tab <- rbindFill(lapply(runs, `[[`, "rows")); rownames(tab) <- NULL
  tab$model <- factor(tab$model, levels = names(formulas)[names(formulas) %in% unique(tab$model)])
  # influence: from the current result when supplied, else from the "current" run
  cur <- runs[["current"]]
  infl <- if (!is.null(influence)) influence else if (!is.null(cur)) list(influence = cur$res$influence[[1]], wide = cur$res$wide[[1]]) else NULL
  summ <- sensitivitySummary(tab[tab$effect == "shift", ], infl, alpha = alpha, rel.change = rel.change, influence.rel = influence.rel)
  rownames(summ) <- NULL
  structure(list(table = tab, summary = summ, formulas = formulas, dropped = dropped, test = t$label, analysis = "expression.shifts", effect = "shift",
                 settings = list(dist = dist, permutation = permutation, n.permutations = n.permutations, alpha = alpha, rel.change = rel.change, influence.rel = influence.rel,
                                 significance = significance.label, jackknife = jackknife)),
            class = c("cacoaSensitivity", "list"))
}

# the range of the estimate over the leave-one-sample-out refits (`infl$shift`: samples x cell types, change of the estimate when the
# sample is left out), as `loo.min` / `loo.max` on the rows of the "shift" effect; NA when no influence was computed for the model
addLeaveOneOutRange <- function(rows, infl) {
  rows$loo.min <- NA_real_; rows$loo.max <- NA_real_
  if (is.null(infl) || is.null(infl$shift)) return(rows)
  sh <- which(rows$effect == "shift" & rows$celltype %in% colnames(infl$shift))
  for (i in sh) { ch <- infl$shift[, rows$celltype[i]]; ch <- ch[is.finite(ch)]; if (length(ch) > 1) { rows$loo.min[i] <- rows$estimate[i] + min(ch); rows$loo.max[i] <- rows$estimate[i] + max(ch) } }
  rows
}

# per cell type: the two verdicts (size across models and single samples; significance across models) from the rows of one effect
# (`sh`: the table restricted to that effect; `infl`: list(influence = list(shift = samples x cell types changes), wide = data.frame with
# `celltype`, `null.shift.lo`, `null.shift.hi`) of the current model, or NULL)
sensitivitySummary <- function(sh, infl, alpha = 0.05, rel.change = 0.5, influence.rel = 0.2) {
  summ <- do.call(rbind, lapply(split(sh, sh$celltype), function(d) {
    cur.row <- d[d$model == "current", ]
    est0 <- if (nrow(cur.row)) cur.row$estimate else d$estimate[1]
    same.sign <- length(unique(sign(d$estimate[is.finite(d$estimate)]))) == 1
    n.sig <- sum(d$padj < alpha, na.rm = TRUE); m <- sum(is.finite(d$p))
    # change across models: on the normalized shift (the magnitude of the main result and of the plot) when there is one
    v <- if (!is.null(d$estimate.norm) && any(is.finite(d$estimate.norm))) d$estimate.norm else d$estimate
    v0 <- if (nrow(cur.row)) v[d$model == "current"][1] else v[1]
    max.rel <- if (is.finite(v0) && v0 != 0) max(abs(v - v0) / abs(v0), na.rm = TRUE) else NA_real_
    worst <- if (nrow(d) > 1) { dd <- d[d$model != "current", ]; vv <- v[d$model != "current"]; dd$model[which.max(abs(vv - v0))] } else NA
    sig.change <- if (nrow(cur.row)) any((d$padj < alpha) != (cur.row$padj < alpha), na.rm = TRUE) else NA
    flipped.by <- if (!is.na(sig.change) && sig.change) paste(d$model[(d$padj < alpha) != (cur.row$padj < alpha) & !is.na(d$padj)], collapse = ", ") else ""
    # influence
    infl.sample <- NA_character_; infl.rel <- NA_real_
    if (!is.null(infl) && !is.null(infl$influence) && d$celltype[1] %in% colnames(infl$influence$shift)) {
      ch <- infl$influence$shift[, d$celltype[1]]; ch <- ch[is.finite(ch)]
      if (length(ch) > 2) {
        i <- which.max(abs(ch)); infl.rel <- abs(ch[i]) / abs(est0)
        others <- abs(ch[-i]); outlying <- abs(ch[i]) > stats::median(others) + 3 * stats::mad(others)
        if (isTRUE(outlying) && is.finite(infl.rel) && infl.rel > influence.rel) infl.sample <- names(ch)[i]
      }
    }
    cur.sig <- nrow(cur.row) && !is.na(cur.row$padj) && cur.row$padj < alpha
    # the two questions: does the size depend on the model (distance magnitudes, influence), does the significance
    m.est <- sum(is.finite(d$estimate))
    # near zero: the current normalized shift lies within the relabeling noise band of the main result (what shuffled labels give),
    # where relative changes across models are meaningless
    wb <- if (!is.null(infl$wide) && all(c("null.shift.lo", "null.shift.hi") %in% names(infl$wide))) infl$wide[match(d$celltype[1], infl$wide$celltype), ] else NULL
    near0 <- nrow(cur.row) == 1 && !is.null(wb) && nrow(wb) == 1 && is.finite(wb$null.shift.lo) && is.finite(wb$null.shift.hi) &&
      is.finite(cur.row$estimate.norm %||% NA_real_) && cur.row$estimate.norm >= wb$null.shift.lo && cur.row$estimate.norm <= wb$null.shift.hi
    v.mag <- if (m.est <= 1) "one model" else if (near0) "near zero (within the relabeling band)"
      else if (!is.na(infl.sample)) sprintf("driven by sample %s", infl.sample)
      else if (!same.sign) "sign depends on the model"
      else if (is.finite(max.rel) && max.rel > rel.change) sprintf("estimate changes by %.0f%% (%s)", 100 * max.rel, as.character(worst))
      else "stable"
    v.sig <- if (m <= 1) sprintf("one model (%s)", if (n.sig) "significant" else "not significant")
      else if (n.sig == 0) sprintf("not significant in any of %d models", m)
      else if (n.sig == m) sprintf("significant in all %d models", m)
      else if (nzchar(flipped.by)) sprintf("significant in %d of %d models (%s with: %s)", n.sig, m, if (cur.sig) "not" else "only", flipped.by)
      else sprintf("significant in %d of %d models", n.sig, m)
    verdict <- if (m <= 1) sprintf("only one model could be fitted (%s)", if (n.sig) "significant" else "not significant")
      else if (n.sig == 0) "no significant effect in any model"
      else if (!cur.sig && !is.na(infl.sample)) "no significant effect in the current model"
      else if (!is.na(infl.sample)) sprintf("driven by sample %s", infl.sample)
      else if (!same.sign) "sign sensitive to the adjustment set"
      else if (!is.na(sig.change) && sig.change) sprintf("significance sensitive to the adjustment set (%s)", flipped.by)
      else if (is.finite(max.rel) && max.rel > rel.change) sprintf("estimate changes by %.0f%% (%s)", 100 * max.rel, as.character(worst))
      else "robust"
    data.frame(celltype = d$celltype[1], estimate.current = est0, same.sign = same.sign, n.significant = n.sig, n.models = m,
               max.rel.change = max.rel, most.changing.model = as.character(worst), influential.sample = infl.sample,
               influence.rel = infl.rel, band.lo = if (!is.null(wb) && nrow(wb) == 1) wb$null.shift.lo else NA_real_,
               band.hi = if (!is.null(wb) && nrow(wb) == 1) wb$null.shift.hi else NA_real_,
               verdict.magnitude = v.mag, verdict.significance = v.sig, verdict = verdict, stringsAsFactors = FALSE)
  }))
  rownames(summ) <- NULL
  summ
}

#' Sensitivity of the compositional result to the covariate set and to single samples
#'
#' The contrast of [lmCoda()] refitted under alternative adjustment sets ([sensitivityFormulas()]: unadjusted; current; current minus
#' each covariate; current plus each of the top screened covariates; all screened), with the per-cell-type loadings (log2, the change of
#' a cell type's log proportion relative to the geometric mean of all cell types) and their permutation p-values; under every model the
#' loadings are also refitted with each sample left out (the ILR fit without permutations), which gives the jackknife standard error
#' and the leave-one-out range of each loading. The verdicts are those of [checkSensitivity()], on the loadings: the size (stable,
#' changes by x %, sign depends on the model, driven by one sample, near zero = within the current model's 95 % permutation band)
#' and the significance (significant in k of m models, by the primary test relative to the geometric mean).
#' @param cnts samples x cell types count matrix (as stored by `cao$estimateCellLoadings()`, `res$cnts`)
#' @param model the stored `cacoaModel` of the result
#' @param meta sample metadata
#' @param formulas named list of location formulas (default: [sensitivityFormulas()])
#' @param screen covariate screen result used to propose additions
#' @param current the stored [lmCoda()] result (its permutation band of the loadings gives the "near zero" verdict)
#' @param n.permutations permutations of each refit
#' @param seed seed of the relabelings (the same for every model)
#' @param alpha significance level on the adjusted p-values
#' @param top.k,rel.change,influence.rel as [checkSensitivity()]
#' @param jackknife leave-one-sample-out refits under every model (default TRUE; cheap: no permutations)
#' @param ... further arguments of [lmCoda()] (`ref.method`, `zero.pseudocount`, ...)
#' @return `cacoaSensitivity` with `analysis = "composition"` and `effect = "loading"`: `table` (one row per model x cell type:
#'   `estimate` = `estimate.norm` = loading in log2, `se.jk`, `loo.min`, `loo.max`, `p`, `padj`), `summary` (verdicts), `formulas`, `dropped`
#' @export
checkCompositionSensitivity <- function(cnts, model, meta, formulas = NULL, screen = NULL, current = NULL, n.permutations = 199, seed = 1, alpha = 0.05,
                                        top.k = 3, rel.change = 0.5, influence.rel = 0.2, jackknife = TRUE, ...) {
  t <- model$tests[[1]]
  if (t$kind != "contrast") stop("sensitivity analysis is defined for contrast tests (two groups or a numeric step)")
  if (is.null(formulas)) formulas <- sensitivityFormulas(model, screen, meta, top.k = top.k)
  if (!"current" %in% names(formulas)) {
    cur <- model$formula.full %||% model$formula
    same <- which(vapply(formulas, formulaKey, character(1)) == formulaKey(cur))
    if (length(same)) names(formulas)[same[1]] <- "current" else formulas <- c(list(current = cur), formulas)
  }
  spec <- model$test.spec %||% list(test = NULL, contrast = t$contrast)
  cnts <- as.matrix(cnts); dots <- list(...)
  runs <- lapply(names(formulas), function(nm) {
    m <- tryCatch(buildCacoaModel(meta, formula = formulas[[nm]], test = spec$test, contrast = spec$contrast, dispersion.formula = model$dispersion.formula,
                                  block.vars = setdiff(model$block.vars, model$random.vars), numeric.ref = model$numeric.ref, numeric.step = model$numeric.step,
                                  permutation = model$permutation %||% "auto", n.permutations = n.permutations), error = function(e) NULL)
    if (is.null(m) || any(m$issues$severity == "error")) return(NULL)
    m$tests <- m$tests[1]
    used <- rownames(m$F) %||% m$samples$used                       # the counts in the order of the model's design rows
    if (length(used) < 4 || !all(used %in% rownames(cnts))) return(NULL)
    x <- tryCatch(do.call(lmCoda, c(list(cnts = cnts[used, , drop = FALSE], model = m, n.permutations = n.permutations, seed = seed), dots)), error = function(e) NULL)
    if (is.null(x)) return(NULL)
    l <- x$contrast$loadings$obs / log(2); cts <- names(l)
    rows <- data.frame(test = t$label, celltype = cts, effect = "loading", estimate = unname(l), estimate.norm = unname(l), se.jk = NA_real_, loo.min = NA_real_, loo.max = NA_real_,
                       p = unname(x$contrast$per_cell$pval[cts]), padj = unname(x$contrast$per_cell$padj[cts]), n = length(used), model = nm,
                       formula = paste(deparse(formulas[[nm]]), collapse = ""), stringsAsFactors = FALSE)
    infl <- NULL
    if (jackknife && length(used) > 3) {   # leave-one-sample-out loadings from the ILR fit (the same contrast, no permutations)
      Fm <- m$F[used, , drop = FALSE]; cF <- m$contrast.F[colnames(Fm)]; ilr <- x$ilr[used, , drop = FALSE]; psi <- x$psi
      loo <- t(vapply(seq_along(used), function(i) { co <- qr.coef(qr(Fm[-i, , drop = FALSE]), ilr[-i, , drop = FALSE]); co[is.na(co)] <- 0
        as.numeric(psi %*% as.numeric(cF %*% co)) / log(2) }, numeric(length(cts))))
      colnames(loo) <- cts; rownames(loo) <- used; nj <- nrow(loo)
      rows$se.jk <- sqrt((nj - 1) / nj * colSums(sweep(loo, 2, colMeans(loo))^2))
      rows$loo.min <- apply(loo, 2, min); rows$loo.max <- apply(loo, 2, max)
      infl <- list(shift = sweep(loo, 2, unname(l)))   # change of the loading when the sample is left out
    }
    list(rows = rows, infl = infl, model = m)
  })
  names(runs) <- names(formulas); runs <- Filter(Negate(is.null), runs)
  dropped <- setdiff(names(formulas), names(runs))
  if (!length(runs)) stop("no alternative model could be fitted")
  tab <- rbindFill(lapply(runs, `[[`, "rows")); rownames(tab) <- NULL
  tab$model <- factor(tab$model, levels = names(formulas)[names(formulas) %in% unique(tab$model)])
  cur <- runs[["current"]]; wide <- NULL
  if (!is.null(cur) && !is.null(current$contrast$per_cell)) {   # the stored result stands in for the current model: its p-values (usually more permutations)
    k <- which(tab$model == "current"); pc <- current$contrast$per_cell
    tab$p[k] <- unname(pc$pval[tab$celltype[k]]); tab$padj[k] <- unname(pc$padj[tab$celltype[k]])
    runs[["current"]]$stored <- TRUE
  }
  # the current model's influence and the permutation band of its loadings (from the stored result when given)
  perm <- current$contrast$loadings$perm %||% NULL
  if (!is.null(perm)) { q <- apply(perm, 1, stats::quantile, c(0.025, 0.975), na.rm = TRUE) / log(2)
    wide <- data.frame(celltype = colnames(q), null.shift.lo = q[1, ], null.shift.hi = q[2, ], stringsAsFactors = FALSE) }
  infl <- if (!is.null(cur)) list(influence = cur$infl, wide = wide) else NULL
  sh <- tab; sh$effect <- "shift"   # the summary code is written for the "shift" effect
  summ <- sensitivitySummary(sh, infl, alpha = alpha, rel.change = rel.change, influence.rel = influence.rel)
  structure(list(table = tab, summary = summ, formulas = formulas, dropped = dropped, test = t$label, analysis = "composition", effect = "loading",
                 settings = list(n.permutations = n.permutations, alpha = alpha, rel.change = rel.change, influence.rel = influence.rel,
                                 significance = "permutation test (vs the geometric mean)", jackknife = jackknife)),
            class = c("cacoaSensitivity", "list"))
}

# short forms of the two verdicts for plot strips
shortMagnitudeVerdict <- function(v) {
  out <- v
  out[startsWith(v, "driven by sample")] <- "driven by one sample"
  out[startsWith(v, "estimate changes")] <- sub("^estimate changes by ([0-9]+%).*$", "changes by \\1", v[startsWith(v, "estimate changes")])
  out[startsWith(v, "near zero")] <- "near zero"
  out
}
shortSignificanceVerdict <- function(v) {
  out <- v
  k <- grepl("^significant in [0-9]+ of [0-9]+", v); out[k] <- sub("^significant in ([0-9]+) of ([0-9]+).*$", "significant in \\1/\\2", v[k])
  k <- grepl("^significant in all", v); out[k] <- sub("^significant in all ([0-9]+) models$", "significant in \\1/\\1", v[k])
  k <- grepl("^not significant in any", v); out[k] <- "not significant"
  out
}

# short form of a verdict for plot strips
shortVerdict <- function(v) {
  out <- v
  out[startsWith(v, "driven by sample")] <- "driven by one sample"
  out[startsWith(v, "sign sensitive")] <- "sign depends on the model"
  out[startsWith(v, "significance sensitive")] <- "significance depends on the model"
  out[startsWith(v, "estimate changes")] <- sub(" \\(.*\\)$", "", v[startsWith(v, "estimate changes")])
  out[startsWith(v, "only one model")] <- "only one model fitted"
  out[startsWith(v, "no significant effect")] <- "no significant effect"
  out
}

#' @exportS3Method base::print
print.cacoaSensitivity <- function(x, ...) {
  cat(sprintf("Sensitivity of '%s' (%s) across %d models: %s\n", x$test, x$effect %||% "shift", length(x$formulas),
              paste(sprintf("%s [%s]", names(x$formulas), vapply(x$formulas, function(f) paste(deparse(f), collapse = ""), character(1))), collapse = "; ")))
  if (length(x$dropped)) cat(sprintf("  not fitted: %s (no cell type could be tested with that covariate set)\n", paste(x$dropped, collapse = "; ")))
  s <- x$summary; s <- s[order(s$verdict != "robust", -abs(s$estimate.current)), ]
  if (is.null(s$verdict.magnitude)) {   # object of an earlier version
    for (i in seq_len(nrow(s))) cat(sprintf("  %-20s %s (significant in %d of %d models)\n", s$celltype[i], s$verdict[i], s$n.significant[i], s$n.models[i]))
    return(invisible(x))
  }
  cat(sprintf("  magnitude: %s; significance: %s (adjusted p < %.2g)\n", if (identical(x$analysis, "composition")) "loading (log2, relative to the geometric mean)" else "distance model", x$settings$significance %||% "distance test", x$settings$alpha))
  w <- max(nchar(s$celltype), 8)
  for (i in seq_len(nrow(s))) cat(sprintf("  %-*s %s\n", w, s$celltype[i], if (s$verdict.magnitude[i] == "one model") sprintf("only one model fitted: %s", sub("^one model \\((.*)\\)$", "\\1", s$verdict.significance[i]))
                                                                        else sprintf("magnitude %s; %s", s$verdict.magnitude[i], s$verdict.significance[i])))
  invisible(x)
}

#' Forest plot of a sensitivity analysis
#'
#' @param x `cacoaSensitivity`
#' @param effect effect to show (default: the result's, `"shift"` or `"loading"`)
#' @param cell.types optional subset
#' @param normalized plot normalized estimates (default TRUE)
#' @param interval `"jackknife"` (default): the jackknife 95 % interval as an error bar and the leave-one-sample-out range as a thin
#'   line, for every model that has them (`jackknife = TRUE` in `checkSensitivity()`; the current model always); `"none"`
#' @param plot.theme ggplot2 theme
#' @return ggplot2 object: cell types as panels (titled with the magnitude and significance verdicts), models on the y-axis,
#'   distance-model estimate with jackknife CI (when available) on the x-axis, filled points significant by the test named in
#'   the subtitle, the current model highlighted
#' @export
plotSensitivity <- function(x, effect = NULL, cell.types = NULL, normalized = TRUE, interval = c("jackknife", "none"), plot.theme = ggplot2::theme_bw()) {
  interval <- match.arg(interval); effect <- effect %||% x$effect %||% "shift"; coda <- identical(x$analysis, "composition")
  d <- x$table[x$table$effect == effect, ]
  if (!is.null(cell.types)) d <- d[d$celltype %in% cell.types, ]
  if (!nrow(d)) stop("nothing to plot")
  val <- if (normalized) "estimate.norm" else "estimate"
  d$value <- d[[val]]
  sc <- if (normalized) ifelse(is.finite(d$estimate) & d$estimate != 0, d$estimate.norm / d$estimate, NA) else 1
  d$lo <- d$value - 1.96 * d$se.jk * abs(sc); d$hi <- d$value + 1.96 * d$se.jk * abs(sc)
  has.loo <- all(c("loo.min", "loo.max") %in% names(d)) && any(is.finite(d$loo.min))
  if (has.loo) { d$loo.lo <- d$loo.min * abs(sc); d$loo.hi <- d$loo.max * abs(sc) }   # the leave-one-out range, on the plotted scale
  if (interval == "none") { d$lo <- d$hi <- NA_real_; has.loo <- FALSE }
  d$sig <- !is.na(d$padj) & d$padj < x$settings$alpha
  d$current <- d$model == "current"
  d$model <- factor(d$model, levels = rev(levels(x$table$model)))
  ord <- x$summary$celltype[order(x$summary$verdict != "robust", -abs(x$summary$estimate.current))]
  d$celltype <- factor(d$celltype, levels = intersect(ord, unique(d$celltype)))
  verd <- if (is.null(x$summary$verdict.magnitude)) shortVerdict(x$summary$verdict) else   # one line per verdict
    ifelse(x$summary$verdict.magnitude == "one model", shortSignificanceVerdict(x$summary$verdict.significance),
           paste(shortMagnitudeVerdict(x$summary$verdict.magnitude), shortSignificanceVerdict(x$summary$verdict.significance), sep = "\n"))
  verd <- setNames(verd, x$summary$celltype)
  d$panel <- factor(sprintf("%s\n%s", as.character(d$celltype), verd[as.character(d$celltype)]), levels = sprintf("%s\n%s", levels(d$celltype), verd[levels(d$celltype)]))
  bands <- NULL
  if (normalized && all(c("band.lo", "band.hi") %in% names(x$summary))) {   # the noise band of shuffled labels (current model), as in the shift plot
    bands <- x$summary[is.finite(x$summary$band.lo) & is.finite(x$summary$band.hi) & x$summary$celltype %in% levels(d$celltype), c("celltype", "band.lo", "band.hi")]
    bands$panel <- factor(sprintf("%s\n%s", bands$celltype, verd[bands$celltype]), levels = levels(d$panel))
  }
  ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$model)) +
    { if (!is.null(bands) && nrow(bands)) ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = .data$band.lo, xmax = .data$band.hi, ymin = -Inf, ymax = Inf),
                                                             inherit.aes = FALSE, fill = "grey90") } +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    { if (has.loo) ggplot2::geom_linerange(ggplot2::aes(xmin = .data$loo.lo, xmax = .data$loo.hi), colour = "grey78", linewidth = 2.4, na.rm = TRUE) } +   # the leave-one-out range, a light band inside the jackknife interval
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), width = 0.25, colour = "grey40", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$sig, colour = .data$current), size = 2.6, fill = "white", stroke = 0.9) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19), labels = c("not significant", "significant"), name = NULL) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "#d73027"), labels = c(`FALSE` = "alternative model", `TRUE` = "current model"), name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = 19)), shape = ggplot2::guide_legend(order = 1)) +
    ggplot2::facet_wrap(~ panel, scales = "free_x") + ggplot2::scale_x_continuous(n.breaks = 3) + plot.theme +
    ggplot2::theme(panel.spacing.x = grid::unit(1.2, "lines")) +
    ggplot2::labs(x = if (coda) "loading (log2, relative to the geometric mean)" else if (normalized) sprintf("normalized %s", effect) else effect, y = NULL,
                  title = sprintf("Sensitivity of the %s to the covariate set", if (coda) "composition loadings" else effect),
                  subtitle = paste(c(paste(c(if (coda) "position: loading" else "position: distance-model magnitude", if (!is.null(bands) && nrow(bands)) "gray band: shuffled labels (current model)",
                                             sprintf("filled: %s adjusted p < %.2g", x$settings$significance %||% "distance test", x$settings$alpha)), collapse = "; "),
                                     if (interval == "jackknife" && (any(is.finite(d$se.jk)) || has.loo))
                                       paste(c(if (any(is.finite(d$se.jk))) "error bars: jackknife 95% intervals", if (has.loo) "light bands: range over the leave-one-sample-out refits"), collapse = "; ")),
                                   collapse = "\n")) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"), strip.text = ggplot2::element_text(size = 7))
}
