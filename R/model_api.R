## Workflow API layer (Track B): options, test grammar, reference-level heuristics, the stored model
## object (`cacoaModel`) and the metadata audit (`describeMetadata`).
##
## Vocabulary used everywhere: a *test* is one question about the design. A test is either a
## "contrast" (two settings of the design are compared: alt vs ref, or a numeric step) and reports
## shift / var / total, or a "term" (a K-level factor tested as a whole) and reports location /
## dispersion / any difference.

# ---- options -------------------------------------------------------------------------------------

#' Default analysis options
#'
#' Options used by every analysis method of a `Cacoa` object unless the call sets the argument explicitly
#' (explicit argument > `setOptions()` > these defaults). See `Cacoa$setOptions()`.
#'
#' @return named list
#' @export
cacoaDefaultOptions <- function() {
  list(
    n.permutations = 999,          # permutations for every permutation test
    seed = 1,                      # integer seed; NULL draws from R's RNG
    alpha = 0.05,                  # significance level for flags and plot marks
    p.adjust.method = "BH",        # adjustment across cell types within a test
    n.cores = 1,
    verbose = TRUE,
    plot.provenance = TRUE,        # provenance subtitle on result plots
    min.samp.per.level = 3,        # minimum samples at each contrasted level per cell type
    dist = "cor",                  # expression distance: gene-centred cosine
    permutation = "auto",          # permutation scheme: auto / block / freedman-lane
    numeric.step = 1,              # step for numeric tests (per unit); per-SD is reported as well
    robust = "none",               # robust fit of the distance model: none / huber / winsor (samples with outlying residual distances)
    na.mode = "drop",              # samples absent from a cell type: drop (fit on the present ones) / impute_weak (kept with a near-zero weight)
    robust.k = 1.345               # robust tuning constant (robust standard deviations of the residual sizes)
  )
}

checkOptionValues <- function(opts) {
  chk <- function(cond, msg) if (!cond) stop(msg, call. = FALSE)
  if (!is.null(opts$n.permutations)) chk(is.numeric(opts$n.permutations) && opts$n.permutations >= 1, "n.permutations must be a positive number")
  if (!is.null(opts$seed)) chk(is.numeric(opts$seed) && length(opts$seed) == 1, "seed must be a single number or NULL")
  chk(is.numeric(opts$alpha) && opts$alpha > 0 && opts$alpha < 1, "alpha must be in (0, 1)")
  chk(opts$p.adjust.method %in% stats::p.adjust.methods, "unknown p.adjust.method")
  chk(is.numeric(opts$n.cores) && opts$n.cores >= 1, "n.cores must be >= 1")
  chk(is.logical(opts$verbose), "verbose must be logical")
  chk(is.logical(opts$plot.provenance), "plot.provenance must be logical")
  chk(is.numeric(opts$min.samp.per.level) && opts$min.samp.per.level >= 2, "min.samp.per.level must be >= 2")
  chk(opts$dist %in% c("cor", "l2", "l1"), "dist must be one of cor, l2, l1")
  chk(opts$permutation %in% c("auto", "block", "freedman-lane"), "unknown permutation scheme")
  chk(is.numeric(opts$numeric.step) && opts$numeric.step != 0, "numeric.step must be a non-zero number")
  chk(opts$robust %in% c("none", "huber", "winsor"), "robust must be one of none, huber, winsor")
  chk(opts$na.mode %in% c("drop", "impute_weak"), "na.mode must be drop or impute_weak")
  chk(is.numeric(opts$robust.k) && opts$robust.k > 0, "robust.k must be positive")
  invisible(opts)
}

# ---- metadata audit (B.2) -------------------------------------------------------------------------

idLikeName <- function(nm) grepl("(^|[._ ])(sample|donor|library|lib|patient|subject|individual|id|ids|barcode|run)([._ ]|$)|_id$|^id$|ids?$",
                                 tolower(nm))

#' Describe sample metadata columns
#'
#' Classifies every column of a sample metadata table by type and by its usefulness as a covariate
#' (`role`). Used by `Cacoa$new()` for the startup summary and by `screenCovariates()` /
#' `checkDesign()` to choose default covariate sets.
#'
#' @param meta data.frame of sample metadata (rows = samples)
#' @param max.missing.frac columns with more missing values than this fraction get role `mostly-missing`
#' @param max.levels.frac factors with more than `max.levels.frac * n` levels get role `high-cardinality`
#' @return data.frame with columns `column`, `type` (factor / numeric / logical / character / date / other),
#'   `n.levels`, `min.level.size`, `n.missing`, `role` (`usable` / `id-like` / `constant` /
#'   `high-cardinality` / `mostly-missing`) and `notes`
#' @export
describeMetadata <- function(meta, max.missing.frac = 0.3, max.levels.frac = 1/3) {
  stopifnot(is.data.frame(meta))
  n <- nrow(meta)
  rows <- lapply(names(meta), function(v) {
    x <- meta[[v]]
    type <- if (is.factor(x)) "factor" else if (is.logical(x)) "logical" else if (is.numeric(x)) "numeric"
            else if (is.character(x)) "character" else if (inherits(x, c("Date", "POSIXt"))) "date" else "other"
    miss <- sum(is.na(x)); xx <- x[!is.na(x)]
    isDiscrete <- type %in% c("factor", "logical", "character")
    nlev <- if (isDiscrete) length(unique(as.character(xx))) else length(unique(xx))
    minlev <- if (isDiscrete && nlev) min(table(as.character(xx))) else NA_integer_
    notes <- character(0)
    role <- "usable"
    if (nlev <= 1) role <- "constant"
    else if (miss > max.missing.frac * n) role <- "mostly-missing"
    else if (isDiscrete && nlev == length(xx) && length(xx) > 2) role <- "id-like"
    else if (isDiscrete && idLikeName(v) && nlev > 2 * length(xx) / 3) role <- "id-like"   # name says ID and nearly one level per sample
    else if (isDiscrete && nlev > max(2, n * max.levels.frac)) role <- "high-cardinality"
    else if (type == "numeric" && idLikeName(v) && nlev == length(xx)) role <- "id-like"
    if (role == "usable") {
      if (isDiscrete && !is.na(minlev) && minlev < 3) notes <- c(notes, sprintf("smallest level has %d sample%s", minlev, if (minlev == 1) "" else "s"))
      if (type == "numeric" && nlev > 5) {
        sk <- skewness(xx)
        if (is.finite(sk) && abs(sk) > 2 && all(xx > 0)) notes <- c(notes, sprintf("skewed (|skew| = %.1f): consider log", abs(sk)))
      }
      if (type == "numeric" && nlev <= 5 && nlev >= 2) notes <- c(notes, sprintf("numeric with %d distinct values: treat as factor?", nlev))
    }
    if (miss > 0 && role != "mostly-missing") notes <- c(notes, sprintf("%d missing", miss))
    data.frame(column = v, type = type, n.levels = nlev, min.level.size = minlev, n.missing = miss, role = role,
               notes = paste(notes, collapse = "; "), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame(column = character(0), type = character(0), n.levels = integer(0), min.level.size = integer(0),
                                      n.missing = integer(0), role = character(0), notes = character(0))
  rownames(out) <- NULL
  out
}

skewness <- function(x) { x <- x[is.finite(x)]; if (length(x) < 3) return(NA_real_); m <- mean(x); s <- stats::sd(x); if (s == 0) return(0); mean(((x - m) / s)^3) }

# compact one-screen summary used by Cacoa$new()
formatMetadataSummary <- function(desc, meta, max.show = 8) {
  n <- nrow(meta)
  usable <- desc[desc$role == "usable", , drop = FALSE]
  one <- function(i) {
    r <- usable[i, ]; x <- meta[[r$column]]
    if (r$type == "numeric") sprintf("%s [num, %s-%s]", r$column, format(min(x, na.rm = TRUE), digits = 3), format(max(x, na.rm = TRUE), digits = 3))
    else if (r$n.levels <= 3) { tb <- table(as.character(x)); sprintf("%s [%s]", r$column, paste(sprintf("%s %d", names(tb), as.integer(tb)), collapse = ", ")) }
    else sprintf("%s [%d levels]", r$column, r$n.levels)
  }
  lines <- sprintf("Sample metadata: %d samples, %d columns", n, ncol(meta))
  if (nrow(usable)) {
    shown <- vapply(seq_len(min(nrow(usable), max.show)), one, character(1))
    if (nrow(usable) > max.show) shown <- c(shown, "...")
    lines <- c(lines, sprintf("  usable (%d): %s", nrow(usable), paste(shown, collapse = ", ")))
  }
  excl <- desc[desc$role != "usable", , drop = FALSE]
  if (nrow(excl)) lines <- c(lines, sprintf("  excluded from screening (%d): %s", nrow(excl),
                                            paste(sprintf("%s (%s)", excl$column, excl$role), collapse = ", ")))
  miss <- desc[desc$n.missing > 0, , drop = FALSE]
  if (nrow(miss)) lines <- c(lines, sprintf("  missing values: %s", paste(sprintf("%s (%d)", miss$column, miss$n.missing), collapse = ", ")))
  paste(lines, collapse = "\n")
}

# ---- reference level (D24) ------------------------------------------------------------------------

controlLikeNames <- function() c("control", "controls", "ctrl", "ctl", "con", "healthy", "hc", "normal", "nl", "wt", "wildtype",
                                 "wild-type", "untreated", "vehicle", "veh", "mock", "baseline", "pre", "none", "no",
                                 "neg", "negative", "0", "false", "sham", "placebo", "ref", "reference")

#' Choose the reference level of a factor for a two-group test
#'
#' Rules (in order): 1) the factor's first level when the level order was set explicitly (not alphabetical);
#' 2) the unique level whose name looks like a control (`control`, `ctrl`, `healthy`, `wt`, `untreated`,
#' `vehicle`, `baseline`, `0`, `FALSE`, ...); 3) the most frequent level.
#'
#' @param x factor, character or logical vector
#' @return list with `level`, `reason` (one of `explicit-order`, `control-like`, `most-frequent`) and `levels`
#' @export
chooseReferenceLevel <- function(x) {
  xx <- x[!is.na(x)]
  levs <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(xx)))
  if (is.logical(x)) return(list(level = "FALSE", reason = "control-like", levels = c("FALSE", "TRUE")))
  if (length(levs) < 2) stop("the test variable needs at least two observed levels")
  if (is.factor(x) && !identical(levs, sort(levs))) return(list(level = levs[1], reason = "explicit-order", levels = levs))
  hit <- levs[tolower(levs) %in% controlLikeNames()]
  if (length(hit) == 1) return(list(level = hit, reason = "control-like", levels = levs))
  tb <- table(factor(as.character(xx), levels = levs))
  list(level = names(tb)[which.max(tb)], reason = "most-frequent", levels = levs)
}

referenceReason <- function(reason) switch(reason,
  "explicit-order" = "first level of the factor as set",
  "control-like"   = "matched a control-like name",
  "most-frequent"  = "most frequent level",
  "user"           = "as requested", reason)

# ---- test grammar (D23, §4.4) ---------------------------------------------------------------------

# "condition: IPF vs control" -> c("condition", "IPF", "control"); NULL when not of that form
parseTestString <- function(s) {
  m <- regmatches(s, regexec("^\\s*([^:]+?)\\s*:\\s*(.+?)\\s+vs\\.?\\s+(.+?)\\s*$", s))[[1]]
  if (length(m) == 4) return(c(m[2], m[3], m[4]))
  NULL
}

mainEffectTerms <- function(formula, meta) {
  tl <- attr(stats::terms(stats::as.formula(formula), data = meta), "term.labels")
  tl <- tl[!grepl(":", tl, fixed = TRUE)]
  vars <- unique(unlist(lapply(tl, function(t) all.vars(str2lang(t)))))
  intersect(vars, names(meta))
}

#' Resolve what to test from the `test` / `contrast` grammar
#'
#' @param test variable name(s), `"all"`, `"var: alt vs ref"`, a `c(var, alt, ref)` triple, or a structured
#'   contrast list (see [buildDesignMatrices()])
#' @param contrast expert synonym of `test` (a triple or structured contrast); either may be given
#' @param meta sample metadata
#' @param formula the location formula (needed for `"all"` and to check variables)
#' @param numeric.step step for numeric variables (per unit of the variable)
#' @return list of test specifications: each has `id`, `kind` (`contrast` / `term`), `variable`, `contrast`
#'   (a specification accepted by [buildDesignMatrices()], or `NULL` for terms), `levels` (`ref`, `alt`),
#'   `step`, `label`, `reference.reason`, `interpretation`
#' @export
resolveTests <- function(test = NULL, contrast = NULL, meta, formula = NULL, numeric.step = 1) {
  if (is.null(test) && is.null(contrast)) return(list())
  if (!is.null(test) && !is.null(contrast)) stop("give either `test` or `contrast`, not both")
  spec <- if (!is.null(test)) test else contrast
  out <- list()
  add <- function(t) { t$id <- length(out) + 1L; out[[length(out) + 1L]] <<- t }
  fromVariable <- function(v, ref = NULL, alt = NULL) {
    if (!v %in% names(meta)) stop(sprintf("test variable '%s' is not a column of the sample metadata", v))
    x <- meta[[v]]
    if (is.numeric(x) && !is.matrix(x)) {
      xx <- x[!is.na(x)]; m <- mean(xx); step <- numeric.step
      if (!is.null(ref) || !is.null(alt)) { ref <- as.numeric(ref); alt <- as.numeric(alt); step <- alt - ref; m <- ref }
      return(list(kind = "contrast", variable = v, step = step, sd = stats::sd(xx),
                  contrast = list(type = "simple", term = v, num = m + step, den = m),
                  levels = c(ref = format(m, digits = 4), alt = format(m + step, digits = 4)),
                  label = sprintf("%s: per %s unit%s", v, format(step, digits = 3), if (step == 1) "" else "s"),
                  reference.reason = NA_character_))
    }
    if (is.logical(x)) x <- factor(x, levels = c("FALSE", "TRUE"))
    levs <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x[!is.na(x)])))
    if (!is.null(ref) && !is.null(alt)) {
      bad <- setdiff(c(ref, alt), levs)
      if (length(bad)) stop(sprintf("level(s) %s not found in '%s' (levels: %s)", paste(sQuote(bad), collapse = ", "), v, paste(levs, collapse = ", ")))
      return(list(kind = "contrast", variable = v, contrast = c(v, alt, ref), levels = c(ref = ref, alt = alt),
                  label = sprintf("%s: %s vs %s", v, alt, ref), reference.reason = "user"))
    }
    if (length(levs) == 2) {
      r <- chooseReferenceLevel(x); alt <- setdiff(levs, r$level)
      return(list(kind = "contrast", variable = v, contrast = c(v, alt, r$level), levels = c(ref = r$level, alt = alt),
                  label = sprintf("%s: %s vs %s", v, alt, r$level), reference.reason = r$reason))
    }
    if (length(levs) < 2) stop(sprintf("'%s' has fewer than two observed levels", v))
    r <- chooseReferenceLevel(x)
    list(kind = "term", variable = v, contrast = NULL, levels = c(ref = r$level, alt = NA), n.levels = length(levs),
         all.levels = levs, label = sprintf("%s (%d levels)", v, length(levs)), reference.reason = r$reason)
  }
  if (is.list(spec) && !is.data.frame(spec)) {                       # structured contrast (expert)
    if (!is.null(spec$type)) {
      sp <- normalizeContrastSpec(spec)
      v <- switch(sp$type, simple = parseTermVars(sp$term)[1], marginal = sp$term, lincomb = parseTermVars(sp$term)[1], coef = NA_character_)
      add(list(kind = "contrast", variable = v, contrast = spec, levels = c(ref = as.character(sp$den %||% NA), alt = as.character(sp$num %||% NA)),
               label = contrastLabel(sp), reference.reason = "user"))
    } else {
      for (s in spec) out <- c(out, resolveTests(s, NULL, meta, formula, numeric.step))
      for (i in seq_along(out)) out[[i]]$id <- i
    }
  } else if (is.numeric(spec) && !is.null(names(spec))) {              # coefficient weights
    add(list(kind = "contrast", variable = NA_character_, contrast = spec, levels = c(ref = NA, alt = NA),
             label = paste(sprintf("%+g*%s", spec, names(spec)), collapse = " "), reference.reason = "user"))
  } else if (is.character(spec)) {
    if (length(spec) == 1 && identical(spec, "all")) {
      if (is.null(formula)) stop("test = \"all\" needs a formula")
      for (v in mainEffectTerms(formula, meta)) add(fromVariable(v))
    } else if (length(spec) == 3 && spec[1] %in% names(meta) && !all(spec %in% names(meta))) {
      add(fromVariable(spec[1], ref = spec[3], alt = spec[2]))         # c(var, alt, ref)
    } else {
      for (s in spec) {
        p <- parseTestString(s)
        if (!is.null(p)) add(fromVariable(p[1], ref = p[3], alt = p[2])) else add(fromVariable(trimws(s)))
      }
    }
  } else stop("unsupported `test` specification")
  for (i in seq_along(out)) { out[[i]]$id <- i; out[[i]]$interpretation <- testInterpretation(out[[i]]) }
  out
}

contrastLabel <- function(sp) {
  switch(sp$type,
         simple = sprintf("%s: %s vs %s%s", sp$term, sp$num, sp$den,
                          if (length(sp$at)) sprintf(" at %s", paste(sprintf("%s = %s", names(sp$at), unlist(sp$at)), collapse = ", ")) else ""),
         marginal = sprintf("%s: %s vs %s (averaged over %s)", sp$term, sp$num, sp$den, paste(sp$over, collapse = ", ")),
         lincomb = sprintf("%s: %s", sp$term, paste(sprintf("%+g*%s", sp$cells, names(sp$cells)), collapse = " ")),
         coef = paste(sprintf("%+g*%s", sp$coefs, names(sp$coefs)), collapse = " "))
}

# plain-language meaning of the three effects for a test
testInterpretation <- function(t) {
  if (t$kind == "term") {
    return(c(location = sprintf("location: the levels of %s differ in where their samples sit (a common direction per level), beyond within-level variability", t$variable),
             dispersion = sprintf("dispersion: the levels of %s differ in how spread out their samples are", t$variable),
             any = sprintf("any difference: samples from different levels of %s are farther apart than samples from the same level", t$variable)))
  }
  if (!is.null(t$step)) {
    return(c(shift = sprintf("shift > 0: samples move in a common direction as %s increases (per %s unit%s)", t$variable, format(t$step, digits = 3), if (t$step == 1) "" else "s"),
             var = sprintf("var > 0: samples become more heterogeneous as %s increases", t$variable),
             total = sprintf("total > 0: samples farther apart in %s are farther apart in expression", t$variable)))
  }
  a <- t$levels[["alt"]]; r <- t$levels[["ref"]]
  if (is.na(a) || is.na(r)) return(c(shift = "shift > 0: the compared settings differ in a common direction, beyond within-group variability",
                                     var = "var > 0: samples at the target setting are more heterogeneous than at the reference setting",
                                     total = "total > 0: samples at the two settings are farther apart than samples at the reference setting are from each other"))
  c(shift = sprintf("shift > 0: %s samples differ from %s samples in a common direction, beyond within-group variability", a, r),
    var = sprintf("var > 0: %s samples are more heterogeneous than %s samples", a, r),
    total = sprintf("total > 0: %s samples are farther from %s samples than %s samples are from each other", a, r, r))
}

# ---- the model object (D22) -----------------------------------------------------------------------

#' Build the stored model object
#'
#' Resolves the tests, builds the design for every test with [buildDesignMatrices()], determines the
#' permutation plan per test and collects design issues. The returned object also carries the primary
#' (first) test's design fields at the top level (`F`, `X`, `Z`, `contrast.F`, `contrast_spec`, ...), so
#' code that expects a `buildDesignMatrices()` result keeps working.
#'
#' @param meta sample metadata (rows named by sample)
#' @param formula location formula (default: `~ <test variables> + <block.vars>`)
#' @param test,contrast what to test (see [resolveTests()])
#' @param dispersion.formula dispersion formula (default: main effects of the test variables)
#' @param block.vars metadata columns defining permutation strata
#' @param numeric.ref `"auto"` or a named list of anchors for numeric covariates
#' @param numeric.step step for numeric tests
#' @param permutation permutation scheme (`"auto"`, `"block"`, `"freedman-lane"`)
#' @param n.permutations number of permutations the plan is made for
#' @param verbosity passed to [buildDesignMatrices()]
#' @return object of class `cacoaModel`
#' @export
buildCacoaModel <- function(meta, formula = NULL, test = NULL, contrast = NULL, dispersion.formula = NULL, block.vars = NULL,
                            numeric.ref = "auto", numeric.step = 1, permutation = "auto", n.permutations = 999,
                            verbosity = "none") {
  stopifnot(is.data.frame(meta))
  tests <- resolveTests(test, contrast, meta, formula, numeric.step)
  test.vars <- unique(stats::na.omit(vapply(tests, function(t) t$variable %||% NA_character_, character(1))))
  if (is.null(formula)) {
    if (!length(test.vars)) stop("cannot build a default formula: give `formula` (and `test`)")
    formula <- stats::as.formula(paste("~", paste(unique(c(test.vars, block.vars)), collapse = " + ")))
    formula.default <- TRUE
  } else {
    formula <- stats::as.formula(checkFormula(formula)); formula.default <- FALSE
  }
  if (!length(tests)) stop("no test given: use `test` (a variable name) or `contrast`")
  fvars <- all.vars(formula)
  miss <- setdiff(fvars, names(meta))
  if (length(miss)) stop("variable(s) not found in sample metadata: ", paste(miss, collapse = ", "))
  desc <- describeMetadata(meta)
  idlike <- intersect(fvars, desc$column[desc$role == "id-like"])
  if (length(idlike)) stop(sprintf("column '%s' has a unique value per sample and cannot be a covariate", idlike[1]))
  # listwise deletion (D29)
  used.vars <- unique(c(fvars, if (!is.null(dispersion.formula)) all.vars(stats::as.formula(dispersion.formula)), block.vars))
  used.vars <- intersect(used.vars, names(meta))
  complete <- stats::complete.cases(meta[, used.vars, drop = FALSE])
  dropped <- rownames(meta)[!complete]
  meta.used <- meta[complete, , drop = FALSE]
  if (is.null(dispersion.formula)) {
    dvars <- intersect(test.vars, names(meta))
    dispersion.formula <- if (length(dvars)) stats::as.formula(paste("~", paste(dvars, collapse = " + "))) else ~ 1
  } else dispersion.formula <- stats::as.formula(checkFormula(dispersion.formula))

  # per-test designs
  for (i in seq_along(tests)) {
    t <- tests[[i]]
    if (t$kind == "contrast") {
      d <- buildDesignMatrices(meta.used, contrast = t$contrast, formula = formula, numericRef = numeric.ref,
                               blockVars = block.vars, verbosity = verbosity)
    } else {   # term: design with an all-zero contrast on the term's columns marks the tested columns
      d <- buildTermDesign(meta.used, formula, t$variable, numeric.ref, block.vars)
    }
    t$design <- d
    plan <- tryCatch(permutationPlan(d, meta.used, scheme = permutation, block.vars = block.vars, n.permutations = n.permutations, max.enumerate = 0),
                     error = function(e) NULL)
    t$permutation <- if (!is.null(plan)) list(scheme = plan$scheme, n.distinct = plan$n.distinct, p.floor = plan$p.floor,
                                              n.strata = nlevels(plan$strata), notes = plan$notes) else NULL
    tests[[i]] <- t
  }
  issues <- designIssues(tests, meta.used, formula, dropped)
  m <- list(formula = formula, formula.default = formula.default, dispersion.formula = dispersion.formula, tests = tests,
            block.vars = block.vars, numeric.ref = numeric.ref, numeric.step = numeric.step, permutation = permutation,
            samples = list(used = rownames(meta.used), dropped = dropped), issues = issues, meta = meta.used)
  # primary design at top level (backward compatibility with buildDesignMatrices() consumers)
  prim <- tests[[1]]$design
  for (nm in names(prim)) m[[nm]] <- prim[[nm]]
  m$primary <- 1L
  class(m) <- c("cacoaModel", "list")
  m
}

# Design for a term test: the full location design plus `term.cols` (which columns belong to the tested
# variable) and a contrast carrying the term's columns (used only to locate them).
buildTermDesign <- function(meta, formula, variable, numeric.ref = "auto", block.vars = NULL) {
  F <- buildFullDesign(formula, meta)
  assign <- attr(F, "assign"); tl <- attr(attr(F, "terms"), "term.labels")
  in.term <- vapply(tl, function(t) variable %in% all.vars(str2lang(t)), logical(1))
  term.cols <- colnames(F)[assign %in% which(in.term)]
  if (!length(term.cols)) stop(sprintf("'%s' is not a term of the formula", variable))
  cF <- setNames(numeric(ncol(F)), colnames(F)); cF[term.cols] <- 1
  Z <- F[, setdiff(colnames(F), term.cols), drop = FALSE]
  list(F = F, X = F[, term.cols, drop = FALSE], Z = if (ncol(Z)) Z else NULL, contrast.F = cF, contrast.X = cF[term.cols],
       qrZ = if (ncol(Z)) qr(Z) else NULL,
       diagnostics = NULL, numeric_ref_used = list(), formula_used = formula, contrast_spec = list(type = "term", term = variable),
       baselines_used = NULL, contrast_endpoints_F = NULL, contrast_endpoints_X = NULL, contrast_label = sprintf("%s (term)", variable),
       contrast_endpoint_labels = NULL, contrast_endpoints_at = NULL, term.cols = term.cols, term.variable = variable)
}

# Issues of a model: structural (aliasing of the test term), warnings (small levels, df budget, permutation
# floor, missing samples) and notes. Returned as a data.frame(severity, message, suggestion).
designIssues <- function(tests, meta, formula, dropped = character(0)) {
  iss <- list()
  add <- function(sev, msg, sug = "") iss[[length(iss) + 1L]] <<- data.frame(severity = sev, message = msg, suggestion = sug, stringsAsFactors = FALSE)
  n <- nrow(meta)
  if (length(dropped)) add("warning", sprintf("%d sample(s) dropped for missing covariate values: %s", length(dropped), paste(dropped, collapse = ", ")),
                           "complete the metadata or remove the covariate")
  F <- tests[[1]]$design$F
  q <- qr(F)$rank
  if (q < ncol(F)) add("warning", sprintf("design is rank-deficient (rank %d < %d columns)", q, ncol(F)), "remove redundant covariates")
  if (n - q < 10) add(if (n - q < 1) "error" else "warning", sprintf("only %d residual degrees of freedom (n = %d, %d parameters)", n - q, n, q),
                      "fewer covariates, or restrict to the test variable")
  for (t in tests) {
    d <- t$design
    if (!is.null(d$Z) && ncol(d$Z) && !is.null(d$X)) {
      Xr <- qr.resid(qr(d$Z), as.matrix(d$X))
      rel <- sqrt(colSums(Xr^2)) / (sqrt(colSums(as.matrix(d$X)^2)) + 1e-15)
      if (any(rel < 1e-8)) {
        culprit <- aliasingCulprit(t, meta, formula)
        add("error", sprintf("'%s' is fully determined by %s: its effect cannot be separated from the adjustment", t$variable %||% t$label, culprit),
            "remove the aliased covariate or test it instead")
      } else if (any(rel < 0.3)) add("warning", sprintf("'%s' is strongly collinear with the other covariates (only %.0f%% of its variation is independent)",
                                                        t$variable %||% t$label, 100 * min(rel)), "check cao$checkDesign() for the association")
    }
    if (t$kind == "contrast" && !is.null(t$variable) && !is.na(t$variable) && !is.numeric(meta[[t$variable]])) {
      tb <- table(as.character(meta[[t$variable]]))[c(t$levels[["ref"]], t$levels[["alt"]])]
      tb <- tb[!is.na(tb)]
      if (length(tb) && min(tb) < 3) add("warning", sprintf("'%s' has only %d sample(s) at level %s", t$variable, min(tb), names(tb)[which.min(tb)]),
                                         "cell types with fewer than min.samp.per.level samples per level are skipped")
    }
    if (!is.null(t$permutation)) {
      if (is.finite(t$permutation$n.distinct) && t$permutation$n.distinct < 200)
        add("warning", sprintf("test '%s': only %d distinct permutations; the smallest attainable p-value is %.3g", t$label,
                               round(t$permutation$n.distinct), t$permutation$p.floor), "fewer strata (block.vars), or more samples")
      for (nt in t$permutation$notes) if (!grepl("distinct permutations", nt)) add("note", sprintf("test '%s': %s", t$label, nt))
    }
  }
  out <- if (length(iss)) do.call(rbind, iss) else data.frame(severity = character(0), message = character(0), suggestion = character(0), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

aliasingCulprit <- function(t, meta, formula) {
  v <- t$variable; if (is.null(v) || is.na(v)) return("the other covariates")
  others <- setdiff(all.vars(formula), v)
  hits <- Filter(function(o) { tb <- table(meta[[v]], meta[[o]]); all(rowSums(tb > 0) <= 1) || all(colSums(tb > 0) <= 1) }, others)
  if (length(hits)) sprintf("'%s' (each %s contains one %s)", hits[1], hits[1], v) else "the other covariates"
}

#' @exportS3Method base::print
print.cacoaModel <- function(x, ...) {
  cat(format(x, ...), sep = "\n"); invisible(x)
}

#' @exportS3Method base::format
format.cacoaModel <- function(x, ...) {
  fm <- function(f) paste(deparse(f), collapse = "")
  lines <- sprintf("Model: %s%s      dispersion: %s", fm(x$formula), if (isTRUE(x$formula.default)) " (default)" else "", fm(x$dispersion.formula))
  for (t in x$tests) {
    adj <- setdiff(all.vars(x$formula), t$variable)
    pl <- t$permutation
    perm.txt <- if (is.null(pl)) "" else sprintf("permutations: %s%s%s", pl$scheme,
                                                   if (pl$scheme == "block" && pl$n.strata > 1) sprintf(" within %d strata", pl$n.strata) else "",
                                                   if (is.finite(pl$n.distinct)) sprintf(" (%s distinct)", format(round(pl$n.distinct), big.mark = ",")) else "")
    ref.txt <- if (!is.null(t$reference.reason) && !is.na(t$reference.reason) && t$kind == "contrast" && is.null(t$step))
      sprintf("  (reference '%s': %s)", t$levels[["ref"]], referenceReason(t$reference.reason)) else ""
    lines <- c(lines, sprintf("Test%s: %s%s", if (length(x$tests) > 1) sprintf(" %d", t$id) else "", t$label, ref.txt))
    lines <- c(lines, sprintf("       %s%s", if (length(adj)) sprintf("adjusted for %s; ", paste(adj, collapse = ", ")) else "unadjusted; ", perm.txt))
    lines <- c(lines, paste0("       ", t$interpretation[1]))
  }
  dr <- x$samples$dropped
  iss <- x$issues
  cnt <- table(factor(iss$severity, levels = c("error", "warning", "note")))
  iss.txt <- if (!nrow(iss)) "none" else paste(sprintf("%d %s%s", cnt, names(cnt), ifelse(cnt == 1, "", "s"))[cnt > 0], collapse = ", ")
  lines <- c(lines, sprintf("Samples: %d used%s.   Issues: %s%s", length(x$samples$used),
                            if (length(dr)) sprintf(", %d dropped (%s)", length(dr), paste(dr, collapse = ", ")) else "",
                            iss.txt, if (nrow(iss)) " (see $issues)" else ""))
  lines
}

#' @exportS3Method base::summary
summary.cacoaModel <- function(object, ...) {
  print(object)
  if (nrow(object$issues)) { cat("\nIssues:\n"); for (i in seq_len(nrow(object$issues))) cat(sprintf("  [%s] %s%s\n", object$issues$severity[i], object$issues$message[i],
                                                                                              if (nzchar(object$issues$suggestion[i])) paste0(" -> ", object$issues$suggestion[i]) else "")) }
  invisible(object)
}

# index of a test within a model: by position, by label ("Group: Group2 vs Group1", "Group (4 levels)") or by variable name
matchModelTest <- function(model, test = NULL) {
  if (is.null(test)) return(1L)
  labels <- vapply(model$tests, `[[`, character(1), "label")
  if (is.numeric(test)) { if (any(test < 1 | test > length(labels))) stop("test index out of range (", length(labels), " tests)"); return(as.integer(test)) }
  vars <- vapply(model$tests, function(t) as.character(t$variable %||% NA_character_)[1], character(1))
  idx <- unique(c(which(labels %in% test), which(vars %in% test)))
  if (!length(idx)) stop("test not found; available: ", paste(labels, collapse = ", "))
  idx
}

# one-line provenance string for plots (D34)
modelProvenance <- function(model, test = NULL, extra = NULL) {
  if (is.null(model)) return(NULL)
  t <- model$tests[[matchModelTest(model, test)[1]]]
  adj <- setdiff(all.vars(model$formula), t$variable)
  parts <- c(t$label, if (length(adj)) sprintf("adjusted for %s", paste(adj, collapse = ", ")) else "unadjusted",
             if (!is.null(t$permutation)) sprintf("%s permutations%s", t$permutation$scheme,
                                                  if (is.finite(t$permutation$n.distinct) && t$permutation$n.distinct < 1e4) sprintf(" (floor %.3g)", t$permutation$p.floor) else ""),
             extra)
  paste(parts, collapse = "; ")
}

#' Regression weights of a contrast
#'
#' The sample weights `w = X (X'X)^- c` such that `w' y` is the least-squares estimate of the contrast `c'beta`
#' for any response `y`. For a balanced two-group comparison without covariates they are +1/n_alt and
#' -1/n_ref; with covariates they are the adjusted weights (which sum to zero within every stratum of a
#' discrete covariate).
#' @param design output of [buildDesignMatrices()] (or a `cacoaModel`)
#' @return named numeric vector over samples
#' @export
regressionWeights <- function(design) {
  X <- as.matrix(design$F); cvec <- design$contrast.F[colnames(X)]
  w <- drop(X %*% MASS::ginv(crossprod(X)) %*% cvec)
  names(w) <- rownames(X)
  w
}
