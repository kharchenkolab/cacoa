#' Sample-level design and contrast for one test
#'
#' @description
#' Builds the design matrix `F` of a sample-level model and the contrast vector `contrast.F` that every engine of
#' cacoa tests, from the sample metadata, a formula and a contrast specification. The workflow entry point is
#' [buildCacoaModel()] (which calls this for every test); this function is the expert constructor for one contrast.
#'
#' ## Coding
#' When the tested variable is a factor it is coded without an intercept and placed first, so its coefficients are
#' the level means at the reference setting of the other covariates (`~ 0 + group + batch` gives
#' `groupA, groupB, batchb2`). The reference level comes first. Other factors keep treatment coding; a numeric
#' tested variable keeps the intercept. Character and logical columns are treated as factors.
#'
#' ## Contrast
#' The contrast is the difference of two design rows (`row(alt) - row(ref)`), each evaluated with the other
#' factors at their reference level and the other numeric covariates at their anchors (`numericRef`: `"auto"` =
#' sample means, or a named list), so transformed covariates (`log(age)`, `poly(age, 2)`) and interactions are
#' handled uniformly. Forms:
#' - a triple `c("group", "B", "A")`: group B against group A; with interactions of `group` in the model, use
#'   one of the structured forms below;
#' - `list(type = "simple", term = "group", num = "B", den = "A", at = list(batch = "b2"))`: the comparison at a
#'   fixed setting of other factors; `term = "group:batch"` with `num = "B:b2"` compares interaction cells;
#' - `list(type = "marginal", term = "group", num = "B", den = "A", over = "batch", weights = "equal")`: the
#'   comparison averaged over the levels of other factors (`"equal"`, `"proportional"` or named weights);
#' - a named numeric vector over `colnames(F)`: coefficient weights (expert).
#'
#' ## Split for the per-column fitter
#' `X = F c / c'c` (one column whose coefficient is the contrast) and `Z = F N` with `N` an orthonormal basis of
#' the null space of `c`; `[X, Z]` spans the same space as `F`. Freedman-Lane permutations residualize on `Z`.
#'
#' @param data data.frame of sample-level covariates (rows = samples)
#' @param contrast contrast specification (see above)
#' @param formula right-hand-side formula; default `~ <contrast variables> + <blockVars>`
#' @param numericRef `"auto"` or a named list of numeric anchors
#' @param na.action NA handler for `model.frame`
#' @param verbosity `"none"`, `"warn"`, `"info"` or `"debug"`
#' @param blockVars metadata columns added to the default formula
#' @param ... deprecated arguments (`numericRefRows`, `tol`, `validate`, `computeQrZ`, `buildBlocks`), accepted and ignored
#' @return list: `F`, `contrast.F`, `X`, `Z`, `contrast.X`, `meta`, `formula_used` (the formula after pruning
#'   non-varying terms), `formula_coded` (with the coding applied), `contrast_spec`, `contrast_endpoints_F`,
#'   `contrast_endpoints_at`, `numeric_ref_used`, `contrast_label`, `contrast_endpoint_labels`, `tested.variable`,
#'   `reference`, `notes`
#' @export
buildDesignMatrices <- function(data, contrast, formula = NULL, numericRef = "auto", na.action = stats::na.pass,
                                verbosity = c("none", "warn", "info", "debug"), blockVars = NULL, ...) {
  verbosity <- match.arg(verbosity)
  dots <- list(...)
  dep <- intersect(names(dots), c("numericRefRows", "tol", "validate", "computeQrZ", "buildBlocks"))
  if (length(dep)) message("buildDesignMatrices(): argument(s) ", paste(dep, collapse = ", "), " are deprecated and ignored")
  unknown <- setdiff(names(dots), dep)
  if (length(unknown)) stop("unknown argument(s): ", paste(unknown, collapse = ", "))
  stopifnot(is.data.frame(data))

  if (is.null(formula)) {
    formula <- buildDefaultFormula(data, contrast = contrast, block.vars = blockVars)
    message("No formula supplied; using ", deparse(formula),
            " (contrast variables", if (length(blockVars)) " and block.vars" else "", "). Pass `formula` to adjust for other covariates.")
  } else formula <- checkFormula(formula)
  formula_used <- pruneFormulaByData(formula, data, na.action = na.action, verbosity = verbosity)
  notes <- character(0)
  dropped <- attr(formula_used, "dropped") %||% character(0)
  if (length(dropped)) notes <- c(notes, sprintf("term%s %s dropped from the formula: constant in the samples used", if (length(dropped) > 1) "s" else "", prettyJoin(dropped)))
  spec <- normalizeContrastSpec(contrast)
  # a plain test variable that interacts with other covariates: compare marginally over interacting factors and at
  # the anchors of interacting numeric covariates (structured contrasts with at = / over = are taken as given)
  idef <- interactionDefault(spec, formula_used, data)
  spec <- idef$spec; contrast <- spec

  # coding: the tested factor first and without an intercept, its reference level first
  tv <- testedVariable(spec, data)
  coding <- designCoding(formula_used, data, tested = tv$variable, ref = tv$ref)
  F <- buildFullDesign(coding$formula, data, na.action = na.action, baselines = coding$baselines)

  cF <- buildSyntheticContrast(F, data, contrast, numericRef = numericRef)
  if (all(abs(cF) < 1e-12)) stop("the contrast is zero on this design (the compared settings coincide)")
  sp <- splitByDirection(F, cF)
  nr <- attr(cF, "numeric_ref_used") %||% list()
  if (length(idef$numeric)) {
    anchors <- vapply(idef$numeric, function(v) format(nr[[v]] %||% mean(data[[v]], na.rm = TRUE), digits = 4), character(1))
    notes <- c(notes, sprintf("%s interacts with %s: compared at %s (%s); use a structured test with at = to change",
                              idef$term, paste(idef$numeric, collapse = ", "), paste(sprintf("%s = %s", idef$numeric, anchors), collapse = ", "),
                              if (is.list(numericRef)) "given anchors" else "the mean"))
  }
  if (length(idef$factors))
    notes <- c(notes, sprintf("%s interacts with %s: compared marginally (equal weights over %s levels); use a structured test with at = or over = to change",
                              idef$term, paste(idef$factors, collapse = ", "), paste(idef$factors, collapse = " x ")))

  contrast_label <- NULL
  contrast_endpoint_labels <- list(baseline = "baseline", target = "target")
  if (spec$type %in% c("simple", "marginal") && !grepl(":", spec$term, fixed = TRUE)) {
    contrast_label <- paste0(spec$term, ": ", spec$num, " vs ", spec$den)
    contrast_endpoint_labels <- list(baseline = paste0(spec$term, " = ", spec$den), target = paste0(spec$term, " = ", spec$num))
  }
  if (verbosity %in% c("info", "debug")) {
    message("Design: n=", nrow(F), ", p=", ncol(F), "; columns: ", prettyJoin(colnames(F)))
    nr <- attr(cF, "numeric_ref_used") %||% list()
    if (length(nr)) message("Numeric anchors: ", paste(sprintf("%s=%s", names(nr), format(unlist(nr), digits = 6)), collapse = ", "))
  }
  list(
    F = F, contrast.F = sp$contrast.F, X = sp$X, Z = sp$Z, contrast.X = sp$contrast.X,
    meta = data,                      # sample metadata (rows = samples of F), used by modelPermutations()
    formula_used = formula_used, formula_coded = coding$formula,
    contrast_spec = spec,
    contrast_endpoints_F = attr(cF, "endpoints_F") %||% NULL,
    contrast_endpoints_at = attr(cF, "endpoints_at") %||% NULL,
    numeric_ref_used = attr(cF, "numeric_ref_used") %||% list(),
    contrast_label = contrast_label, contrast_endpoint_labels = contrast_endpoint_labels,
    tested.variable = tv$variable, reference = tv$ref, notes = notes
  )
}

# Default comparison for a plain test variable (grammar: a name, a triple, "var: a vs b", or a numeric step) whose
# term interacts with other covariates: marginal with equal weights over the interacting factors; interacting numeric
# covariates sit at their anchors (which the endpoint rows do already). Returns the (possibly rewritten) spec and the
# interacting variables for the notes.
interactionDefault <- function(spec, formula, data) {
  out <- list(spec = spec, term = NULL, factors = character(0), numeric = character(0))
  if (!identical(spec$type, "simple") || !isTRUE(spec$plain_triple %||% spec$plain)) return(out)
  if (grepl(":", spec$term, fixed = TRUE) || length(spec$at %||% list())) return(out)
  tl <- attr(stats::terms(stats::as.formula(formula), data = data), "term.labels")
  ints <- tl[grepl(":", tl, fixed = TRUE) & vapply(tl, function(t) spec$term %in% all.vars(str2lang(t)), logical(1))]
  if (!length(ints)) return(out)
  partners <- setdiff(unique(unlist(lapply(ints, function(t) all.vars(str2lang(t))))), spec$term)
  partners <- intersect(partners, names(data))
  facs <- partners[vapply(partners, function(v) !is.numeric(data[[v]]), logical(1))]
  nums <- setdiff(partners, facs)
  sp <- spec; sp$plain_triple <- NULL; sp$plain <- NULL
  if (length(facs)) sp <- list(type = "marginal", term = spec$term, num = spec$num, den = spec$den, at = list(), over = facs, weights = "equal")
  list(spec = sp, term = spec$term, factors = facs, numeric = nums)
}

# The factor whose levels a contrast compares (NULL for numeric, coefficient-level or multi-variable cell contrasts),
# with its reference level (the denominator).
testedVariable <- function(spec, data) {
  none <- list(variable = NULL, ref = NULL)
  if (!spec$type %in% c("simple", "marginal")) return(none)
  if (grepl(":", spec$term, fixed = TRUE)) {
    vars <- parseTermVars(spec$term); numc <- parseCell(spec$num, vars); denc <- parseCell(spec$den, vars)
    diff <- vars[unlist(numc) != unlist(denc)]
    if (length(diff) != 1 || !diff %in% names(data) || is.numeric(data[[diff]])) return(none)
    return(list(variable = diff, ref = as.character(denc[[diff]])))
  }
  v <- spec$term
  if (!v %in% names(data) || is.numeric(data[[v]])) return(none)
  list(variable = v, ref = as.character(spec$den))
}

# Coding of the design for a tested factor: its main-effect term first, no intercept, reference level first.
# Other variables keep their coding. Numeric tested variables (or none) leave the formula unchanged.
designCoding <- function(formula, data, tested = NULL, ref = NULL) {
  f <- if (inherits(formula, "formula")) formula else stats::as.formula(formula)
  baselines <- list()
  if (is.null(tested) || !tested %in% names(data) || is.numeric(data[[tested]])) return(list(formula = f, baselines = baselines))
  if (!is.null(ref)) baselines[[tested]] <- ref
  trm <- stats::terms(f, data = data); tl <- attr(trm, "term.labels")
  main <- which(tl == tested)
  if (!length(main)) return(list(formula = f, baselines = baselines))      # tested variable only inside functions / interactions
  tl2 <- c(tl[main], tl[-main])
  f2 <- stats::as.formula(paste("~ 0 +", paste(tl2, collapse = " + ")), env = environment(f))
  list(formula = f2, baselines = baselines)
}

# Orthonormal basis of the null space of C' (C: p x k): the directions of the coefficient space that C does not test.
nullBasis <- function(C) {
  C <- as.matrix(C); p <- nrow(C)
  q <- qr(C); r <- q$rank
  if (r >= p) return(matrix(0, p, 0))
  Q <- qr.Q(q, complete = TRUE)
  Q[, seq.int(r + 1, p), drop = FALSE]
}

# Split of the design by the contrast direction: X = F c / c'c (coefficient = c'beta), Z = F N (the rest).
splitByDirection <- function(F, cF, name = "contrast") {
  cF <- cF[colnames(F)]
  X <- F %*% matrix(cF / sum(cF^2), ncol = 1); colnames(X) <- name; rownames(X) <- rownames(F)
  N <- nullBasis(matrix(cF, ncol = 1))
  Z <- if (ncol(N)) F %*% N else NULL
  if (!is.null(Z)) { colnames(Z) <- paste0("nuisance", seq_len(ncol(Z))); rownames(Z) <- rownames(F) }
  list(X = X, Z = Z, contrast.F = stats::setNames(as.numeric(cF), colnames(F)), contrast.X = stats::setNames(1, name))
}

# Design restricted to a set of samples and kept full rank: columns that vanish on those samples are dropped,
# dependent columns are removed by QR pivoting, and the contrast (vector or matrix) is re-expressed on the kept
# columns (w = F_k' (F_s^+)' c gives the same estimate as c on the full design for every response). NULL when a
# contrast is not estimable on the subset (a compared level absent, or confounded with a covariate there).
subsetDesign <- function(F, C, rows, tol = 1e-8) {
  Fs <- F[rows, , drop = FALSE]; C <- as.matrix(C)
  if (is.null(colnames(C))) colnames(C) <- paste0("c", seq_len(ncol(C)))
  if (!all(apply(C, 2, function(cv) isEstimable(Fs, cv, tol)))) return(NULL)
  keep <- colSums(abs(Fs)) > tol
  q <- qr(Fs[, keep, drop = FALSE]); r <- q$rank
  if (r < sum(keep)) { piv <- colnames(Fs)[keep][q$pivot[seq_len(r)]]; keep <- colnames(Fs) %in% piv }
  Fk <- Fs[, keep, drop = FALSE]
  Ck <- t(Fk) %*% t(MASS::ginv(Fs)) %*% C
  dimnames(Ck) <- list(colnames(Fk), colnames(C)); Ck[abs(Ck) < 1e-10] <- 0
  list(F = Fk, C = Ck, keep = keep)
}

# ---- Defaults & pruning ----

checkFormula <- function(formula) {
  # Accept formula or character; ensure "~" present
  if (inherits(formula, "formula")) {
    formula.str <- paste(deparse(formula), collapse = "")
  } else if (is.character(formula)) {
    formula.str <- formula
  } else {
    stop("Design formula must be a string or formula object.")
  }
  if (!grepl("~", formula.str)) {
    stop("Design formula must contain a '~' to separate response and predictors.")
  }
  
  # random-effect terms are not supported: say so instead of rewriting the formula
  if (grepl("|", formula.str, fixed = TRUE)) {
    rand.eff.vars <- trimws(unlist(regmatches(formula.str, gregexpr("(?<=\\|)[^\\)]+", formula.str, perl = TRUE))))
    stop(sprintf("random-effect terms (with '|') are not supported: use %s as a fixed effect in the formula or as block.vars (permutation strata)",
                 paste(sQuote(rand.eff.vars), collapse = ", ")), call. = FALSE)
  }
  formula.str
}

# Default sample-level formula: the variables named by the contrast plus block.vars (never every
# metadata column, which would pull in ID-like columns).
buildDefaultFormula <- function(data, contrast = NULL, block.vars = NULL) {
  stopifnot(is.data.frame(data))
  spec <- if (is.null(contrast)) NULL else normalizeContrastSpec(contrast)
  vars <- character(0)
  if (!is.null(spec)) {
    vars <- switch(spec$type,
                   simple   = c(parseTermVars(spec$term), names(spec$at)),
                   marginal = c(spec$term, spec$over, names(spec$at)),
                   coef     = stop("Cannot derive a default formula from a coefficient-level contrast; please supply `formula`."))
  }
  vars <- unique(c(vars, block.vars))
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("Variable(s) not found in sample metadata: ", paste(miss, collapse = ", "))
  if (!length(vars)) return(as.formula("~ 1"))
  as.formula(paste("~", paste(vars, collapse = " + ")))          # the tested factor's coding is applied by designCoding()
}

pruneFormulaByData <- function(formula, data, na.action = stats::na.pass, verbosity = c("none","warn","info","debug")) {
  verbosity <- match.arg(verbosity)
  f <- if (inherits(formula, "formula")) formula else as.formula(formula)
  
  trm <- stats::terms(f, data = data)
  fac <- attr(trm, "factors")
  tlabels <- colnames(fac)
  intercept <- attr(trm, "intercept")
  if (length(tlabels) == 0L) return(f)
  
  mf <- stats::model.frame(trm, data, na.action = na.action)
  rowvars <- rownames(fac)
  varEffective <- vapply(rowvars, function(v) {
    x <- mf[[v]]
    if (is.factor(x) || is.character(x)) {
      length(levels(droplevels(factor(x)))) >= 2L
    } else if (is.numeric(x) || is.logical(x)) {
      ux <- unique(x[!is.na(x)])
      length(ux) >= 2L
    } else TRUE
  }, logical(1))
  
  keep_col <- logical(ncol(fac))
  for (j in seq_len(ncol(fac))) {
    involved <- rowvars[ fac[, j] > 0 ]
    keep_col[j] <- all(varEffective[involved])
  }
  
  kept <- tlabels[keep_col]
  dropped <- tlabels[!keep_col]
  
  rhs <- if (length(kept)) paste(kept, collapse = " + ") else if (intercept == 0L) "0" else "1"
  if (length(kept) && intercept == 0L) rhs <- paste("0 +", rhs)
  f2 <- as.formula(paste("~", rhs), env = environment(f))
  
  if (length(dropped)) {
    msg <- paste0("Dropping non-varying term(s) from formula: ", prettyJoin(dropped))
    if (verbosity %in% c("warn")) warning(msg, call. = FALSE)
    if (verbosity %in% c("info","debug")) message(msg)
    attr(f2, "dropped") <- dropped                  # recorded as a model note by the builders
  }
  f2
}

# ---- Pretty helpers ----

`%||%` <- function(a, b) if (is.null(a)) b else a

prettyJoin <- function(x, maxShow = 12) {
  x <- as.character(x)
  if (!length(x)) return("(none)")
  if (length(x) <= maxShow) return(paste(x, collapse = ", "))
  paste0(paste(x[1:maxShow], collapse = ", "), sprintf(", ... (+%d more)", length(x) - maxShow))
}

parseTermVars <- function(termStr) strsplit(termStr, ":", fixed = TRUE)[[1]]

parseCell <- function(cell, vars) {
  parts <- strsplit(cell, ":", fixed = TRUE)[[1]]
  if (length(parts) != length(vars))
    stop(sprintf("Cell '%s' must have %d parts for {%s}.",
                 cell, length(vars), paste(vars, collapse = ", ")))
  as.list(stats::setNames(parts, vars))
}

# ---- Contrast spec + synthetic builder ----

normalizeContrastSpec <- function(contrast) {
  if (is.numeric(contrast) && !is.null(names(contrast)))
    return(list(type="coef", coefs=contrast))
  
  if (is.character(contrast) && length(contrast) == 3) {
    term <- contrast[[1]]; a <- contrast[[2]]; b <- contrast[[3]]
    if (grepl(":", term, fixed = TRUE))
      return(list(type="simple", term=term, num=a, den=b, at=list()))
    return(list(type="simple", term=term, num=a, den=b, at=list(), plain_triple=TRUE))
  }
  
  if (is.list(contrast) && length(contrast) > 0) {
    t <- contrast$type %||% stop("Structured contrast requires 'type'.")
    if (t == "coef") {
      co <- contrast$coefs %||% stop("type='coef' needs named numeric 'coefs'.")
      if (!is.numeric(co) || is.null(names(co))) stop("'coefs' must be named numeric.")
      return(list(type="coef", coefs=co))
    }
    if (t == "lincomb") stop("contrast type 'lincomb' is no longer supported: use type = 'marginal' with numeric weights, or named coefficient weights")
    if (t %in% c("simple","marginal")) {
      term <- contrast$term %||% stop("Field 'term' is required.")
      num  <- contrast$num  %||% stop("Field 'num' is required.")
      den  <- contrast$den  %||% stop("Field 'den' is required.")
      at   <- contrast$at   %||% list()
      if (t == "simple") return(list(type="simple", term=term, num=num, den=den, at=at, plain=isTRUE(contrast$plain)))
      over <- as.character(contrast$over %||% stop("type='marginal' needs 'over'."))
      w    <- contrast$weights %||% "equal"
      return(list(type="marginal", term=term, num=num, den=den, at=at, over=over, weights=w))
    }
  }
  stop("Unsupported contrast format.")
}

buildFullDesign <- function(formula, data, na.action = stats::na.pass,
                            contrasts.arg = NULL, baselines = NULL) {
  rhs <- if (inherits(formula, "formula")) formula else as.formula(formula)
  data <- prepareDesignData(data, baselines)

  # model.frame() on the raw data: its "terms" attribute carries `predvars` (e.g. poly() coefficients),
  # which is what lets .oneRowFromFormula() evaluate transformed covariates on new rows.
  mf  <- stats::model.frame(rhs, data, na.action = na.action)
  trm <- attr(mf, "terms")
  F   <- stats::model.matrix(trm, mf, contrasts.arg = contrasts.arg)

  # Attach metadata so downstream builders (.oneRowFromFormula, etc.) can reconstruct rows
  attr(F, "terms")     <- trm
  attr(F, "xlevels")   <- lapply(mf, function(x) if (is.factor(x)) levels(x) else NULL)
  attr(F, "contrasts") <- attr(F, "contrasts")

  F
}

# Character and logical columns become factors and requested baselines are applied, on the raw data so that
# transformed terms (log(age), poly(age, 2), factor(batch), ...) can still be evaluated from it.
prepareDesignData <- function(data, baselines = NULL) {
  for (v in names(data)) {
    if (is.character(data[[v]])) data[[v]] <- droplevels(factor(data[[v]]))
    else if (is.logical(data[[v]])) data[[v]] <- factor(data[[v]], levels = c("FALSE", "TRUE"))
  }
  for (v in names(baselines)) {
    x <- data[[v]]
    if (!is.null(x) && is.factor(x)) {
      lev <- as.character(baselines[[v]])
      x <- droplevels(x)
      if (length(lev) && lev %in% levels(x)) data[[v]] <- stats::relevel(x, ref = lev)
    }
  }
  data
}


resolveNumericRef <- function(F, data, contrast, numericRef = "auto", numericRefRows = NULL, tol = 1e-6) {
  if (is.list(numericRef)) return(numericRef)
  if (!identical(numericRef, "auto"))
    stop("numericRef must be a named list or 'auto'.")
  
  spec <- normalizeContrastSpec(contrast)
  trm  <- attr(F, "terms"); if (is.null(trm)) stop("F must carry 'terms'.")
  tlabels <- attr(trm, "term.labels")
  mf   <- stats::model.frame(trm, data, na.action = stats::na.pass)
  
  termVars <- switch(spec$type,
                     simple   = if (grepl(":", spec$term, fixed=TRUE)) parseTermVars(spec$term) else spec$term,
                     marginal = spec$term,
                     coef     = character(0))
  termVars <- unique(termVars)
  
  isNum <- vapply(mf, function(x) is.numeric(x) && !is.matrix(x), logical(1))
  numVars <- names(which(isNum))
  
  interacts <- function(a, b, terms) {
    any(grepl(":", terms) &
          grepl(paste0("(^|:)", a, "(:|$)"), terms) &
          grepl(paste0("(^|:)", b, "(:|$)"), terms))
  }
  
  need <- character(0)
  for (v in numVars) {
    if (spec$type == "simple" && length(termVars) == 1 && termVars[1] == v) next
    if (length(termVars) && any(sapply(termVars, function(tv) interacts(v, tv, tlabels))))
      need <- c(need, v)
  }
  need <- unique(need)
  
  rows <- if (is.null(numericRefRows)) seq_len(nrow(data)) else {
    if (is.logical(numericRefRows)) which(numericRefRows) else as.integer(numericRefRows)
  }
  
  out <- list()
  for (v in need) {
    x <- data[[v]]
    m <- mean(x[rows], na.rm = TRUE)
    out[[v]] <- if (is.finite(m) && abs(m) < tol) 0 else m
  }
  out
}

# One design row for a given covariate setting `at` (named by RAW metadata variables).
# Factors not in `at` sit at their first (reference) level; numerics not in `at` sit at `numericRef[[v]]`
# or, failing that, at their sample mean. New data is built from the raw variables so that transformed
# terms (log(age), poly(age, 2), factor(batch), ...) are evaluated the same way as on the original data.
.oneRowFromFormula <- function(trm, xlevels, contrastsArg, targetCols,
                               data, at = list(), numericRef = NULL) {
  trm  <- stats::delete.response(trm)
  vars <- all.vars(trm)
  data <- prepareDesignData(data)
  # xlevels are keyed by model-frame column names (e.g. "factor(batch)"); map them to raw variables
  xlevVar <- function(v) {
    if (!is.null(xlevels[[v]])) return(xlevels[[v]])
    for (nm in names(xlevels)) {
      if (!is.null(xlevels[[nm]]) && identical(all.vars(str2lang(nm)), v)) return(xlevels[[nm]])
    }
    NULL
  }
  newd <- vector("list", length(vars)); names(newd) <- vars
  for (v in vars) {
    x <- data[[v]]
    if (is.null(x)) stop(sprintf("Variable '%s' of the formula is not in the sample metadata.", v))
    levs <- xlevVar(v)
    if (is.factor(x)) {
      levs <- levs %||% levels(x)
      newd[[v]] <- factor(levs[1], levels = levs)
    } else if (!is.null(levs)) {                  # numeric/logical wrapped as factor(...) in the formula
      newd[[v]] <- if (is.logical(x)) as.logical(levs[1]) else as.numeric(levs[1])
    } else {
      ref <- numericRef[[v]]
      if (is.null(ref)) ref <- mean(x, na.rm = TRUE)
      newd[[v]] <- ref
    }
  }
  for (v in names(at)) {
    if (!v %in% vars) stop(sprintf("Variable '%s' not in model terms.", v))
    levs <- xlevVar(v)
    if (is.factor(newd[[v]])) {
      val <- as.character(at[[v]])
      if (!val %in% levs) stop(sprintf("Level '%s' not in levels(%s): {%s}", val, v, paste(levs, collapse=", ")))
      newd[[v]] <- factor(val, levels = levs)
    } else newd[[v]] <- at[[v]]
  }
  nd  <- as.data.frame(newd, stringsAsFactors = FALSE)
  mf1 <- stats::model.frame(trm, nd, xlev = xlevels, na.action = stats::na.pass)
  r1  <- stats::model.matrix(trm, mf1, contrasts.arg = contrastsArg)
  miss <- setdiff(targetCols, colnames(r1))
  if (length(miss)) r1 <- cbind(r1, matrix(0, 1, length(miss), dimnames = list(NULL, miss)))
  drop(r1[1, targetCols, drop = FALSE])
}

# Evaluate stored endpoint settings (attr "endpoints_at" of a synthetic contrast, see
# buildSyntheticContrast) on an arbitrary design `F` built by buildFullDesign() -- e.g. the dispersion
# design. Returns list(num =, den =) rows in the column space of `F`, each a weighted sum over settings.
endpointRowsFromDesign <- function(F, data, endpoints_at, numericRef = list()) {
  if (is.null(endpoints_at)) return(NULL)
  trm <- attr(F, "terms"); xlv <- attr(F, "xlevels"); ctr <- attr(F, "contrasts")
  if (is.null(trm) || is.null(xlv)) stop("F must carry 'terms' and 'xlevels' (use buildFullDesign()).")
  vars <- all.vars(stats::delete.response(trm))
  # settings for variables the design does not contain (e.g. the batch setting of an interaction-cell
  # contrast evaluated on a dispersion design without batch) do not affect the row and are dropped
  one <- function(at) .oneRowFromFormula(trm, xlv, ctr, colnames(F), data, at[intersect(names(at), vars)], numericRef)
  evalSide <- function(side) {
    r <- setNames(numeric(ncol(F)), colnames(F))
    for (e in side) r <- r + e$w * one(e$at)
    r
  }
  list(num = evalSide(endpoints_at$num), den = evalSide(endpoints_at$den))
}

buildSyntheticContrast <- function(F, data, contrast,
                                   numericRef = "auto",
                                   numericRefRows = NULL) {
  spec <- normalizeContrastSpec(contrast)
  trm  <- attr(F, "terms")
  xlv  <- attr(F, "xlevels")
  ctr  <- attr(F, "contrasts")
  if (is.null(trm) || is.null(xlv))
    stop("F must carry 'terms' and 'xlevels' (use buildFullDesign()).")
  
  numRef <- resolveNumericRef(
    F, data, contrast,
    numericRef     = numericRef,
    numericRefRows = numericRefRows
  )
  
  one <- function(at) .oneRowFromFormula(trm, xlv, ctr, colnames(F), data, at, numRef)
  
  cF <- setNames(numeric(ncol(F)), colnames(F))
  
  ## ------------------------------------------------------------
  ## 1) Coefficient-level contrast (no canonical endpoints)
  ## ------------------------------------------------------------
  if (spec$type == "coef") {
    miss <- setdiff(names(spec$coefs), colnames(F))
    if (length(miss))
      stop("Unknown coefficient(s) in contrast: ", paste(miss, collapse = ", "))
    
    cF[names(spec$coefs)] <- as.numeric(spec$coefs)
    attr(cF, "numeric_ref_used") <- numRef
    attr(cF, "endpoints_F")      <- NULL
    attr(cF, "endpoints_at")     <- NULL
    return(cF)
  }
  
  ## ------------------------------------------------------------
  ## 4) Simple contrasts
  ##    - Either single factor:   term = "Group", num="B", den="A"
  ##    - Or interaction cell:    term = "Group:Batch", num="B:Batch1", ...
  ##    We now also store endpoints in F-space.
  ## ------------------------------------------------------------
  if (spec$type == "simple") {
    if (grepl(":", spec$term, fixed = TRUE)) {
      ## interaction term: num/den are "A:B" style cells
      vars  <- parseTermVars(spec$term)
      atNum <- utils::modifyList(spec$at, parseCell(spec$num, vars))
      atDen <- utils::modifyList(spec$at, parseCell(spec$den, vars))
      rNum  <- one(atNum)
      rDen  <- one(atDen)
      cF    <- rNum - rDen
      atNum_ <- atNum; atDen_ <- atDen
    } else {
      ## simple single factor contrast
      atN  <- utils::modifyList(spec$at, setNames(list(spec$num), spec$term))
      atD  <- utils::modifyList(spec$at, setNames(list(spec$den), spec$term))
      rNum <- one(atN)
      rDen <- one(atD)
      cF   <- rNum - rDen
      atNum_ <- atN; atDen_ <- atD
    }
    attr(cF, "numeric_ref_used") <- numRef
    ## endpoints in F-space (num = alt, den = ref) and the covariate settings that produced them
    attr(cF, "endpoints_F")  <- list(num = rNum, den = rDen)
    attr(cF, "endpoints_at") <- list(num = list(list(at = atNum_, w = 1)), den = list(list(at = atDen_, w = 1)))
    return(cF)
  }
  
  ## ------------------------------------------------------------
  ## 5) Marginal contrasts
  ##    - Average (over=) across other factors, with weights
  ##    - We accumulate weighted endpoints (num/den) in F-space.
  ## ------------------------------------------------------------
  if (spec$type == "marginal") {
    ov  <- spec$over
    xlv <- attr(F, "xlevels")
    levs <- lapply(ov, function(v) {
      L <- xlv[[v]]
      if (is.null(L))
        stop(sprintf("'%s' in over= is not a factor.", v))
      L
    })
    names(levs) <- ov
    grid <- do.call(
      expand.grid,
      c(levs, stringsAsFactors = FALSE)
    )
    
    ## Weights over the grid
    ws <- spec$weights
    if (is.character(ws)) {
      if (!ws %in% c("equal", "proportional"))
        stop("weights must be 'equal', 'proportional', or named numeric.")
      w <- rep(1 / nrow(grid), nrow(grid))
      if (ws == "proportional") {
        tab <- data[, ov, drop = FALSE]
        for (v in ov) tab[[v]] <- factor(tab[[v]], levels = levs[[v]])
        idx <- do.call(interaction, c(tab, drop = TRUE, sep = ":"))
        key <- apply(grid, 1, function(r) paste(r, collapse = ":"))
        cnt <- tapply(rep(1, nrow(tab)), idx, sum)
        w   <- as.numeric(cnt[key])
        w[is.na(w)] <- 0
        if (sum(w) == 0) {
          w <- rep(1 / nrow(grid), nrow(grid))
        } else {
          w <- w / sum(w)
        }
      }
    } else if (is.numeric(ws)) {
      key <- apply(grid, 1, function(r) paste(r, collapse = ":"))
      if (is.null(names(ws)))
        stop("Numeric weights must be named by: ", paste(key, collapse = ", "))
      w <- as.numeric(ws[key])
      if (any(is.na(w)) || any(w < 0) || sum(w) == 0)
        stop("Invalid numeric weights.")
      w <- w / sum(w)
    } else stop("Unsupported weights spec.")
    
    ## Accumulate contrast and endpoints
    rNumTot <- cF  # same length / names as columns of F
    rDenTot <- cF
    atNumAll <- list(); atDenAll <- list()

    for (i in seq_len(nrow(grid))) {
      at_i <- as.list(grid[i, , drop = FALSE])
      atN  <- utils::modifyList(at_i, spec$at); atD <- atN
      atN[[spec$term]] <- spec$num
      atD[[spec$term]] <- spec$den
      
      rNum <- one(atN)
      rDen <- one(atD)
      
      cF      <- cF      + w[i] * (rNum - rDen)
      rNumTot <- rNumTot + w[i] * rNum
      rDenTot <- rDenTot + w[i] * rDen
      atNumAll[[i]] <- list(at = atN, w = w[i]); atDenAll[[i]] <- list(at = atD, w = w[i])
    }
    attr(cF, "numeric_ref_used") <- numRef
    ## weighted endpoints in F-space and the settings/weights that produced them
    attr(cF, "endpoints_F")  <- list(num = rNumTot, den = rDenTot)
    attr(cF, "endpoints_at") <- list(num = atNumAll, den = atDenAll)
    return(cF)
  }
  
  stop("Unhandled contrast type.")
}



# Level-difference contrasts of a factor on a design: K - 1 columns row(level_k) - row(ref), evaluated with the other
# factors at their reference level and the numeric covariates at their anchors (the main effect at the reference
# setting). `levels` orders the levels (reference first).
termContrastMatrix <- function(F, data, variable, levels, numericRef = list()) {
  trm <- attr(F, "terms"); xlv <- attr(F, "xlevels"); ctr <- attr(F, "contrasts")
  one <- function(l) .oneRowFromFormula(trm, xlv, ctr, colnames(F), data, stats::setNames(list(l), variable), numericRef)
  rows <- t(sapply(levels, one)); if (length(levels) == 1) rows <- matrix(rows, 1, dimnames = list(levels, colnames(F)))
  C <- t(rows[-1, , drop = FALSE] - matrix(rows[1, ], nrow(rows) - 1, ncol(rows), byrow = TRUE))
  dimnames(C) <- list(colnames(F), paste(levels[-1], "vs", levels[1]))
  C[abs(C) < 1e-12] <- 0
  C
}

# Extract sample-level variable names used by a formula
varsFromFormula <- function(formula, data, na.action = stats::na.pass) {
  f  <- if (inherits(formula, "formula")) formula else as.formula(formula)
  mf <- stats::model.frame(f, data, na.action = na.action)
  
  # RHS-only variables (no response in our use, but keep it robust)
  rhs_vars <- names(mf)
  
  # classify by *sample-level* types
  isNum <- rhs_vars[vapply(rhs_vars, function(v) is.numeric(data[[v]]) && !is.matrix(data[[v]]), logical(1))]
  isFac <- rhs_vars[vapply(rhs_vars, function(v) is.factor(data[[v]]) || is.character(data[[v]]), logical(1))]
  list(numeric = isNum, factors = isFac, all = rhs_vars)
}

