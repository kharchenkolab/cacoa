#' Prepare core/nuisance design matrices from metadata and a linear contrast
#'
#' @description
#' **`buildDesignMatrices()`** takes a sample-level metadata table, an optional
#' model formula, and a *synthetic* linear contrast, and returns:
#' - the **full model matrix** `F`,
#' - a **core** sub-matrix `X` containing the columns used by the contrast,
#' - a **nuisance** sub-matrix `Z` containing the remaining columns,
#' - the contrast vector in coefficient space (`contrast.F` and `contrast.X`),
#' - optional **FL plumbing** (`qrZ`) and **permutation blocks** diagnostics.
#'
#' ### What happens if `formula` is not supplied?
#' If `formula` is `NULL`, a safe default is built from the supplied `data`:
#'
#' - If **any factor/character** column is present: a **saturated** design
#'   is used (no intercept): `~ 0 + col1 + col2 + ...` so every factor level
#'   has its own dummy column and no level is dropped.
#' - If **all columns are numeric**: an **intercept-included** design is used:
#'   `~ col1 + col2 + ...` (i.e., includes `(Intercept)`).
#'
#' After that, the formula is **pruned** to drop non-varying variables and any
#' interaction terms that depend on them (safe to call on arbitrary metadata).
#'
#' ### How the split into `X` and `Z` works
#' The user-provided `contrast` is converted into a named numeric vector over
#' `colnames(F)` (synthetic; no data-weighting). Columns with
#' `|contrast.F| > tol` form **`X`** and the rest form **`Z`**.
#'
#' **No-nuisance promotion:** if `Z` is empty or only the intercept, we promote
#' to a full-core fit: `X <- F`, `Z <- NULL`, `qrZ <- NULL`. This guarantees the
#' intended design when there is nothing to adjust for.
#'
#' ### Numeric anchors (for interactions)
#' When your contrast involves *interactions with numerics*, those numeric
#' variables are evaluated at **means** (overall or within `numericRefRows`), or
#' you can provide fixed anchors via `numericRef = list(age=35, bmi=22)`.
#'
#' ### Supported contrast specifications
#' All contrasts must be **linear in coefficients**. Supported forms:
#'
#' 1) **DESeq2 triple** (single factor)
#' ```r
#' contrast <- c("group","B","A")  # mu(group=B) - mu(group=A)
#' ```
#'
#' 2) **Simple triple on an interaction**
#' ```r
#' contrast <- list(type="simple", term="group:batch", num="B:1", den="A:1")
#' ```
#'
#' 3) **Marginal contrast** (average over other factor(s))
#' ```r
#' contrast <- list(type="marginal", term="group", num="B", den="A",
#'                  over="batch", weights="equal")        # or "proportional" / named numeric
#' ```
#'
#' 4) **Linear combination over cells of a term**
#' ```r
#' contrast <- list(type="lincomb", term="group:batch",
#'                  cells=c("B:1"=1, "A:1"=-1, "A:2"=-0.5))
#' ```
#'
#' 5) **Direct coefficient-space weights** (named numeric over `colnames(F)`)
#' ```r
#' contrast <- c("(Intercept)"=0, "groupB"=1, "age"=-0.01, "groupB:age"=0.02)
#' ```
#'
#' **Ambiguous triples:** a plain triple like `c("group","B","A")` is ambiguous
#' if the model contains interactions with `group`. In that case, use
#' `type="simple"` with `at=` to fix other factor levels or `type="marginal"`
#' with `over=` to average across them.
#'
#' @param data A `data.frame` of sample-level covariates (factor/character/numeric).
#' @param contrast A linear contrast (see "Supported contrast specifications").
#' @param formula RHS formula for the design (`~ ...`). If `NULL`, a default
#'   is constructed as described above.
#' @param na.action NA handler for `model.frame`.
#' @param numericRef `"auto"` or named list of numeric anchors (e.g., `list(age=35)`).
#' @param numericRefRows Optional row indices / logical mask to compute mean anchors.
#' @param tol Numeric; threshold for selecting non-zero contrast columns into `X`.
#' @param tolRow Per-row activity threshold used to compute `core.rows`.
#' @param validate Logical; compute diagnostics (rank, aliasing, VIF, permutation checks).
#' @param verbosity `"none"|"warn"|"info"|"debug"`.
#' @param computeQrZ Logical; if `TRUE` return `qrZ` of `Z` for Freedman-Lane.
#' @param blockVars Optional character vector of factor names to define permutation blocks.
#' @param buildBlocks Logical; if `TRUE` compute `blocks`, permutation groups, and diagnostics.
#'
#' @return A list with:
#' \item{F}{Full model matrix (with attributes `terms`, `xlevels`, `contrasts`).}
#' \item{X}{Core submatrix (or `F` if no nuisance).}
#' \item{Z}{Nuisance submatrix (or `NULL` if none).}
#' \item{contrast.F, contrast.X}{Contrast vectors aligned to `F` and `X`.}
#' \item{core.rows}{Logical mask of rows with non-negligible activity in `X`.}
#' \item{qrZ}{QR decomposition of `Z` for Freedman-Lane (or `NULL`).}
#' \item{diagnostics}{If `buildBlocks=TRUE`, design diagnostics.}
#' \item{numeric_ref_used, formula_used, baselines_used, contrast_spec}{Metadata.}
#'
#' @examples
#' # 1) Simple factor triple (no interactions); default formula is saturated:
#' # out <- buildDesignMatrices(data=df, contrast=c("group","B","A"))
#'
#' # 2) Interaction triple (fix other factor):
#' # ctr <- list(type="simple", term="group:batch", num="B:1", den="A:1")
#' # out <- buildDesignMatrices(~ group*batch + age, data=df, contrast=ctr)
#'
#' # 3) Marginal across batch with proportional weights:
#' # ctr <- list(type="marginal", term="group", num="B", den="A",
#' #             over="batch", weights="proportional")
#' # out <- buildDesignMatrices(~ group*batch + age, data=df, contrast=ctr, buildBlocks=TRUE)
#'
#' @export
buildDesignMatrices <- function(data, contrast, 
                                formula = NULL,
                                numericRef = "auto",
                                numericRefRows = NULL,
                                tol = 1e-12,
                                na.action = stats::na.pass,
                                tolRow = sqrt(.Machine$double.eps),
                                validate = TRUE,
                                verbosity = c("none","warn","info","debug"),
                                computeQrZ = TRUE,
                                blockVars = NULL,
                                buildBlocks = TRUE) {
  verbosity <- match.arg(verbosity)
  
  # Default formula
  if (is.null(formula)) {
    formula <- buildDefaultFormula(data, contrast = contrast, block.vars = blockVars)
    message("No formula supplied; using ", deparse(formula),
            " (contrast variables", if (length(blockVars)) " and block.vars" else "", "). ",
            "Pass `formula` to adjust for other covariates.")
  } else { # check for mixed effect models
    formula <- checkFormula(formula)
  }
  
  # Prune non-varying terms
  formula_used <- pruneFormulaByData(formula, data, na.action = na.action, verbosity = verbosity)
  
  spec <- try(normalizeContrastSpec(contrast), silent = TRUE)
  
  # Pick baselines so contrasted levels are not dropped (when intercept is present)
  baselines <- chooseBaselinesForSpec(formula_used, data, spec)
  
  # Build F
  F <- buildFullDesign(formula_used, data, na.action = na.action, baselines = baselines)
  
  # Synthetic contrast over F
  cF <- buildSyntheticContrast(F, data, contrast,
                               numericRef = numericRef,
                               numericRefRows = numericRefRows)
  
  # Split by contrast (+ auto-promotion when no nuisance)
  sp <- splitByContrast(F, cF, tol = tol, tolRow = tolRow, promoteIfNoNuisance = TRUE)
  X <- sp$X; Z <- sp$Z
  
  ## extract endpoints in F-space, if present
  endpoints_F <- attr(cF, "endpoints_F") %||% NULL
  ## endpoints in X-space (restricted to core design columns)
  endpoints_X <- NULL
  if (!is.null(endpoints_F) && !is.null(sp$contrast.X)) {
    endpoints_X <- lapply(endpoints_F, function(v) v[colnames(X)])
  }
  
  ## human-readable labels for contrast and endpoints
  contrast_label <- NULL
  contrast_endpoint_labels <- list(
    baseline = "baseline",
    target   = "target"
  )
  if (!inherits(spec, "try-error") && !is.null(spec)) {
    if (is.list(spec) &&
        spec$type %in% c("simple", "marginal") &&
        !grepl(":", spec$term, fixed = TRUE)) {
      
      var <- spec$term
      num <- spec$num
      den <- spec$den
      
      contrast_label <- paste0(var, ": ", num, " vs ", den)
      contrast_endpoint_labels <- list(
        baseline = paste0(var, " = ", den),
        target   = paste0(var, " = ", num)
      )
    }
  }
  
  # qrZ for FL
  qrZ <- if (computeQrZ && !is.null(Z)) qr(as.matrix(Z)) else NULL
  
  # Optional blocks & permutation groups
 
  blocks <- NULL                      # permutation cells now come from permutationPlan() / modelPermutations()
  
  # Diagnostics
  diag <- NULL
  if (validate) {
    diag <- diagnoseDesign(F = F, X = X, Z = Z, qrZ = qrZ,
                           meta = data, blocks = blocks, core.rows = sp$core.rows,
                           contrastSpec = if (!inherits(spec, "try-error")) spec else NULL,
                           tol = tol)
    emitDiagnostics(diag, verbosity)
  }
  
  # Report
  reportContrastInfo(F, X, Z, sp$contrast.F,
                     numericRefUsed = attr(cF, "numeric_ref_used") %||% list(),
                     tol = tol, verbosity = verbosity)
  
  list(
    F = F, X = X, Z = Z,
    contrast.F = sp$contrast.F,
    contrast.X = sp$contrast.X,
    core.rows = sp$core.rows,
    qrZ = qrZ,
    diagnostics = diag,
    numeric_ref_used = attr(cF, "numeric_ref_used") %||% list(),
    formula_used = formula_used,
    contrast_spec = if (!inherits(spec, "try-error")) spec else NULL,
    baselines_used = baselines,
    contrast_endpoints_F = endpoints_F,
    contrast_endpoints_X = endpoints_X,
    contrast_endpoints_at = attr(cF, "endpoints_at") %||% NULL,
    contrast_label = contrast_label,
    contrast_endpoint_labels = contrast_endpoint_labels,
    meta = data                      # sample metadata (rows = samples of F), used by modelPermutations()
  )
}



# =====================================================================
# Helpers (shared by the two public constructors)
# =====================================================================

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
  
  # Demote random effects to fixed, preserving variables
  containsRandomEffects <- grepl("\\([^\\|]*\\|[^\\)]*\\)", formula.str)
  if (containsRandomEffects) {
    rand.eff.vars <- unlist(regmatches(formula.str, gregexpr("(?<=\\|)[^\\)]+", formula.str, perl = TRUE)))
    rand.eff.vars <- trimws(rand.eff.vars)
    warning(sprintf(
      "Random effects (terms with '|') are not supported in this workflow. The variable(s) '%s' will be treated as fixed effects.",
      paste(rand.eff.vars, collapse = ", ")
    ))
    formula.str <- gsub("\\([^\\|]*\\|[^\\)]*\\)", "", formula.str)
    rhs <- gsub("~", "", formula.str)
    rhs <- gsub("\\++", "+", rhs)
    rhs <- gsub("^\\s*\\+|\\+\\s*$", "", rhs)
    rhs <- trimws(rhs)
    rhs.terms <- trimws(unlist(strsplit(rhs, "\\+")))
    rhs.terms <- unique(c(rhs.terms, rand.eff.vars))
    rhs.terms <- rhs.terms[rhs.terms != ""]
    return( paste("~", paste(rhs.terms, collapse = " + ")) )
  } else {
    return(formula.str)
  }
}

# Default sample-level formula: the variables named by the contrast plus block.vars (never every
# metadata column, which would pull in ID-like columns). Uses `~ 0 + ...` coding when a factor is present,
# as the previous default did.
buildDefaultFormula <- function(data, contrast = NULL, block.vars = NULL) {
  stopifnot(is.data.frame(data))
  spec <- if (is.null(contrast)) NULL else normalizeContrastSpec(contrast)
  vars <- character(0)
  if (!is.null(spec)) {
    vars <- switch(spec$type,
                   simple   = c(parseTermVars(spec$term), names(spec$at)),
                   marginal = c(spec$term, spec$over, names(spec$at)),
                   lincomb  = c(parseTermVars(spec$term), names(spec$at)),
                   coef     = stop("Cannot derive a default formula from a coefficient-level contrast; please supply `formula`."))
  }
  vars <- unique(c(vars, block.vars))
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("Variable(s) not found in sample metadata: ", paste(miss, collapse = ", "))
  if (!length(vars)) return(as.formula("~ 1"))
  anyFac <- any(vapply(data[vars], function(x) is.factor(x) || is.character(x), logical(1)))
  if (anyFac) {
    as.formula(paste("~ 0 +", paste(vars, collapse = " + ")))
  } else {
    as.formula(paste("~", paste(vars, collapse = " + ")))
  }
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
  f2 <- as.formula(paste("~", rhs))
  
  if (length(dropped)) {
    msg <- paste0("Dropping non-varying term(s) from formula: ", prettyJoin(dropped))
    if (verbosity %in% c("warn")) warning(msg, call. = FALSE)
    if (verbosity %in% c("info","debug")) message(msg)
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
    if (t == "lincomb") {
      term  <- contrast$term  %||% stop("type='lincomb' needs 'term', e.g. 'group:batch'.")
      cells <- contrast$cells %||% stop("type='lincomb' needs named numeric 'cells'.")
      if (!is.numeric(cells) || is.null(names(cells))) stop("'cells' must be named numeric.")
      return(list(type="lincomb", term=term, cells=cells, at=contrast$at %||% list()))
    }
    if (t %in% c("simple","marginal")) {
      term <- contrast$term %||% stop("Field 'term' is required.")
      num  <- contrast$num  %||% stop("Field 'num' is required.")
      den  <- contrast$den  %||% stop("Field 'den' is required.")
      at   <- contrast$at   %||% list()
      if (t == "simple") return(list(type="simple", term=term, num=num, den=den, at=at))
      over <- as.character(contrast$over %||% stop("type='marginal' needs 'over'."))
      w    <- contrast$weights %||% "equal"
      return(list(type="marginal", term=term, num=num, den=den, at=at, over=over, weights=w))
    }
  }
  stop("Unsupported contrast format.")
}

chooseBaselinesForSpec <- function(formula, data, spec) {
  if (inherits(spec, "try-error")) return(list())
  
  # If the model is already saturated (~0 + ...), there is no intercept,
  # so we do NOT need to pick baselines at all.
  trm <- stats::terms(
    if (inherits(formula,"formula")) formula else as.formula(formula),
    data = data
  )
  if (attr(trm, "intercept", exact = TRUE) == 0L) {
    return(list())
  }
  
  # helper accessors
  vars_in_model <- all.vars(trm)
  levels_in_data <- lapply(
    intersect(vars_in_model, names(data)),
    function(v) {
      x <- data[[v]]
      if (is.factor(x) || is.character(x))
        levels(droplevels(factor(x)))
      else
        NULL
    }
  )
  names(levels_in_data) <- intersect(vars_in_model, names(data))
  
  isFacVar  <- function(v) !is.null(levels_in_data[[v]])
  varLevels <- function(v) levels_in_data[[v]] %||% character(0)
  
  # choose a baseline that:
  #  1. is NOT one of the contrasted levels if possible,
  #  2. otherwise falls back to old heuristic (den, or "not num", etc.).
  chooseBaselineAvoiding <- function(v, avoid_levels, fallback_first = NULL, fallback_not = NULL) {
    L <- varLevels(v)
    if (!length(L)) return(NULL)
    
    # First preference: any level not involved in the contrast at all
    cand <- setdiff(L, avoid_levels)
    if (length(cand)) return(cand[1])
    
    # Second preference: caller-supplied "fallback_first" (often denominator)
    if (!is.null(fallback_first) && fallback_first %in% L) {
      return(fallback_first)
    }
    
    # Third preference: "anything not fallback_not"
    if (!is.null(fallback_not)) {
      cand2 <- setdiff(L, fallback_not)
      if (length(cand2)) return(cand2[1])
    }
    
    # Last resort: just take the first level
    L[1]
  }
  
  bases <- list()
  
  # ---- Case 1: simple single-factor contrast (no ":")
  # contrast like list(type="simple", term="Group", num="Group2", den="Group1")
  # or DESeq2 style c("Group","Group2","Group1")
  if (is.list(spec) &&
      spec$type == "simple" &&
      !grepl(":", spec$term, fixed = TRUE)) {
    
    v <- spec$term
    if (isFacVar(v)) {
      # avoid BOTH contrasted levels so both show up explicitly
      avoid_levels <- unique(c(spec$num, spec$den))
      bases[[v]] <- chooseBaselineAvoiding(
        v,
        avoid_levels = avoid_levels,
        fallback_first = spec$den,   # old behavior fallback
        fallback_not   = spec$num
      )
    }
    return(bases)
  }
  
  # ---- Case 2: simple contrast on an interaction (term contains ":")
  # e.g. list(type="simple", term="Group:Batch",
  #           num="Group2:Batch1", den="Group1:Batch1")
  if (is.list(spec) &&
      spec$type == "simple" &&
      grepl(":", spec$term, fixed = TRUE)) {
    
    vs <- parseTermVars(spec$term)         # c("Group","Batch")
    numC <- parseCell(spec$num, vs)        # list(Group="Group2", Batch="Batch1")
    denC <- parseCell(spec$den, vs)        # list(Group="Group1", Batch="Batch1")
    
    for (v in vs) {
      if (isFacVar(v)) {
        # avoid BOTH numerator and denominator levels for this variable
        avoid_levels <- unique(c(numC[[v]], denC[[v]]))
        bases[[v]] <- chooseBaselineAvoiding(
          v,
          avoid_levels = avoid_levels,
          # fallbacks: pick denominator level first if needed
          fallback_first = denC[[v]],
          fallback_not   = numC[[v]]
        )
      }
    }
    return(bases)
  }
  
  # ---- Case 3: marginal contrast on a single factor
  # e.g. list(type="marginal", term="Group", num="Group2", den="Group1", over="Batch", ...)
  # For baseline purposes, 'marginal' is conceptually still "Group2 vs Group1".
  if (is.list(spec) &&
      spec$type == "marginal" &&
      !grepl(":", spec$term, fixed = TRUE)) {
    
    v <- spec$term
    if (isFacVar(v)) {
      avoid_levels <- unique(c(spec$num, spec$den))
      bases[[v]] <- chooseBaselineAvoiding(
        v,
        avoid_levels = avoid_levels,
        fallback_first = spec$den,
        fallback_not   = spec$num
      )
    }
    return(bases)
  }
  
  # ---- Case 4: lincomb
  # lincomb can hit multiple cells of something like "Group:Batch".
  # We'll keep your previous heuristic: pick the "largest-weight" cell,
  # then avoid that level for each factor. This is already decent.
  if (is.list(spec) && spec$type == "lincomb") {
    vs <- parseTermVars(spec$term)
    w  <- spec$cells
    pick <- if (any(w > 0)) {
      names(w)[which.max(w)]
    } else {
      names(w)[which.max(abs(w))]
    }
    if (length(pick)) {
      pickC <- parseCell(pick, vs)
      for (v in vs) {
        if (isFacVar(v)) {
          bases[[v]] <- chooseBaselineAvoiding(
            v,
            avoid_levels  = pickC[[v]],
            fallback_first = pickC[[v]],
            fallback_not   = NULL
          )
        }
      }
    }
    return(bases)
  }
  
  # ---- Case 5: marginal on interaction or other exotic structures
  # (No clear "two-level focus" to protect. Fall back to old behavior that
  #  nudges baselines away from the numerator if we can.)
  if (is.list(spec) && spec$type == "marginal") {
    v <- spec$term
    if (isFacVar(v)) {
      avoid_levels <- unique(c(spec$num, spec$den))
      bases[[v]] <- chooseBaselineAvoiding(
        v,
        avoid_levels = avoid_levels,
        fallback_first = spec$den,
        fallback_not   = spec$num
      )
    }
    return(bases)
  }
  
  # ---- Coef-style and anything else:
  # coef-style directly names columns; there is no concept of "baseline".
  # return empty -> let model.matrix pick defaults.
  list()
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

# Character columns become factors and requested baselines are applied, on the raw data so that
# transformed terms (log(age), poly(age, 2), factor(batch), ...) can still be evaluated from it.
prepareDesignData <- function(data, baselines = NULL) {
  for (v in names(data)) {
    if (is.character(data[[v]])) data[[v]] <- droplevels(factor(data[[v]]))
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
                     lincomb  = parseTermVars(spec$term),
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
                                   numericRefRows = NULL,
                                   stopOnAmbiguousTriple = TRUE) {
  spec <- normalizeContrastSpec(contrast)
  trm  <- attr(F, "terms")
  xlv  <- attr(F, "xlevels")
  ctr  <- attr(F, "contrasts")
  if (is.null(trm) || is.null(xlv))
    stop("F must carry 'terms' and 'xlevels' (use buildFullDesign()).")
  
  tlabels <- attr(trm, "term.labels")
  hasIntWith <- function(var) {
    any(
      grepl(":", tlabels) &
        grepl(paste0("(^|:)", var, "(:|$)"), tlabels)
    )
  }
  
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
  ## 2) Ambiguous triple guard (simple DESeq2-style triple with interactions)
  ## ------------------------------------------------------------
  if (identical(spec$plain_triple, TRUE) &&
      stopOnAmbiguousTriple &&
      hasIntWith(spec$term)) {
    stop(sprintf(
      "Ambiguous triple: interactions with '%s' present; ",
      spec$term
    ),
    "please specify type='simple' (add at=) or type='marginal' (add over=).")
  }
  
  ## ------------------------------------------------------------
  ## 3) General linear combination of cells (no canonical endpoints)
  ## ------------------------------------------------------------
  if (spec$type == "lincomb") {
    vars <- parseTermVars(spec$term)
    for (nm in names(spec$cells)) {
      at_i <- utils::modifyList(spec$at %||% list(), parseCell(nm, vars))
      cF   <- cF + as.numeric(spec$cells[[nm]]) * one(at_i)
    }
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


# ---- Split by contrast (+ promotion) ----

computeCoreRows <- function(X, tolRow = sqrt(.Machine$double.eps)) {
  if (!ncol(X)) return(rep(FALSE, nrow(X)))
  colScale <- pmax(1, apply(abs(X), 2, max, na.rm = TRUE))
  thr <- matrix(tolRow * colScale, nrow(X), ncol(X), byrow = TRUE)
  rowSums(abs(X) > thr) > 0
}

splitByContrast <- function(F, cF,
                            tol    = 1e-12,
                            tolRow = sqrt(.Machine$double.eps),
                            promoteIfNoNuisance = TRUE,
                            interceptName = "(Intercept)") {
  if (!all(colnames(F) %in% names(cF)))
    stop("F has columns absent from contrast vector: ",
         paste(setdiff(colnames(F), names(cF)), collapse = ", "))
  cF <- cF[colnames(F)]
  
  S  <- which(abs(cF) > tol)
  X  <- if (length(S)) F[, S, drop = FALSE] else F[, 0, drop = FALSE]
  Z  <- F[, setdiff(seq_len(ncol(F)), S), drop = FALSE]
  core.rows <- computeCoreRows(X, tolRow = tolRow)
  contrast.X <- cF[colnames(X)]
  
  if (promoteIfNoNuisance) {
    noNuisance <- (ncol(Z) == 0L) || (length(setdiff(colnames(Z), interceptName)) == 0L)
    if (noNuisance) {
      X <- F
      Z <- NULL
      contrast.X <- cF[colnames(X)]
      core.rows  <- computeCoreRows(X, tolRow = tolRow)
    }
  }
  
  list(X = X, Z = Z, contrast.F = cF, contrast.X = contrast.X, core.rows = core.rows)
}




# ---- Diagnostics ----

diagnoseDesign <- function(F = NULL, X = NULL, Z = NULL, qrZ = NULL,
                           meta = NULL, blocks = NULL, core.rows = NULL,
                           contrastSpec = NULL,
                           block.factors = NULL,
                           thresholds = list(
                             alias.tol     = 1e-8,
                             kappa.warn    = 1e3,
                             vif.warn      = 10,
                             min.core.size = 2L,
                             show.top      = 10,
                             min.eff.perm  = 100
                           ),
                           tol = 1e-12) {
  qrRankDiag <- function(M, name="F") {
    if (is.null(M) || !is.matrix(M) || !ncol(M)) {
      return(list(name=name, rank=0L, p=0L, dep=character(0), kappa=NA_real_))
    }
    q  <- qr(M, LAPACK = TRUE); p <- ncol(M); r <- q$rank
    dep <- if (r < p) colnames(M)[q$pivot[(r+1):p]] else character(0)
    s <- svd(M, nu=0, nv=0)$d
    k <- if (length(s)) (max(s) / max(1e-12, min(s))) else NA_real_
    list(name=name, rank=r, p=p, dep=dep, kappa=k)
  }
  aliasXbyZ <- function(X, Z, qrZ=NULL, tol=1e-8) {
    if (is.null(X) || !is.matrix(X) || !ncol(X) || is.null(Z) || !ncol(Z))
      return(list(aliased=character(0), X.r=X))
    if (is.null(qrZ)) qrZ <- qr(as.matrix(Z))
    Xr <- qr.resid(qrZ, as.matrix(X))
    rn <- sqrt(colSums(Xr^2)); xn <- sqrt(colSums(as.matrix(X)^2)) + 1e-15
    aliased <- colnames(X)[rn / xn < tol]
    list(aliased = aliased, X.r = Xr)
  }
  pinvSym <- function(A, eps=1e-8) {
    ev <- eigen(A, symmetric=TRUE); d <- ev$values; V <- ev$vectors
    di <- ifelse(d > eps * max(1, d[1]), 1/d, 0)
    V %*% (diag(di, nrow=length(di))) %*% t(V)
  }
  vif <- function(Xr, interceptName="(Intercept)", varTol=1e-12) {
    if (is.null(Xr) || !is.matrix(Xr)) return(numeric(0))
    if (ncol(Xr) <= 1) return(numeric(0))
    keep <- setdiff(colnames(Xr), interceptName)
    if (!length(keep)) return(numeric(0))
    Xk <- Xr[, keep, drop = FALSE]
    sds <- apply(Xk, 2, sd, na.rm = TRUE)
    keep <- keep[is.finite(sds) & sds > varTol]
    if (length(keep) < 2) return(numeric(0))
    Xc <- scale(Xr[, keep, drop = FALSE], center = TRUE, scale = TRUE)
    R  <- stats::cor(Xc)
    if (any(!is.finite(R))) return(numeric(0))
    VIF <- tryCatch(diag(solve(R)), error=function(e) diag(pinvSym(R)))
    setNames(as.numeric(VIF), keep)
  }
  
  msgs <- character(0); warns <- character(0)
  
  if (is.null(F) || !ncol(F)) stop("F has 0 columns.")
  if (!nrow(F)) stop("F has 0 rows.")
  if (!is.null(Z) && ncol(Z) == 0) warns <- c(warns, "Z has 0 columns (no nuisance predictors).")
  if (!is.null(X) && ncol(X) == 0) warns <- c(warns, "X has 0 columns (contrast is all zeros at the given tolerance).")
  
  if (!is.null(Z) && ncol(Z) > 0 && "(Intercept)" %in% colnames(F) && "(Intercept)" %in% colnames(X))
    warns <- c(warns, "Intercept is in X (usually you want it in Z).")
  
  dF  <- qrRankDiag(F, name="F")
  dX  <- qrRankDiag(X, name="X")
  dZ  <- qrRankDiag(Z, name="Z")
  if (dF$rank < dF$p) warns <- c(warns, sprintf("F is rank-deficient: rank %d < %d. Dependent: %s", dF$rank, dF$p, prettyJoin(dF$dep)))
  if (is.finite(dF$kappa) && dF$kappa >= thresholds$kappa.warn)
    warns <- c(warns, sprintf("F ill-conditioned (kappa=%.2e). Estimates may be unstable.", dF$kappa))
  
  aXZ <- aliasXbyZ(X, Z, qrZ = qrZ, tol = thresholds$alias.tol)
  if (length(aXZ$aliased)) warns <- c(warns, paste0("Core columns aliased by Z (not estimable after adjustment): ", prettyJoin(aXZ$aliased)))
  vifs <- vif(aXZ$X.r)
  if (length(vifs)) {
    bad <- vifs[vifs >= thresholds$vif.warn]
    if (length(bad)) warns <- c(warns, paste0("High VIFs in core after Z: ",
                                              paste(sprintf("%s=%.1f", names(bad), bad), collapse = ", "),
                                              ". Consider collapsing levels or using a single contrast regressor."))
  }
  
  perm <- NULL                      # permutation cells are described by permutationPlan() / modelPermutations()
  
  list(
    rankF = dF, rankX = dX, rankZ = dZ,
    aliased.by.Z = aXZ$aliased,
    vif = vifs,
    permutation = perm,
    warnings = warns,
    messages = msgs
  )
}

emitDiagnostics <- function(diag, verbosity) {
  if (verbosity %in% c("warn","info","debug")) {
    if (length(diag$warnings)) warning(paste(diag$warnings, collapse = "\n"), call. = FALSE)
  }
  invisible(NULL)
}

reportContrastInfo <- function(F, X, Z, cF, numericRefUsed, tol, verbosity) {
  if (verbosity %in% c("none","warn")) return(invisible(NULL))
  ncolZ <- if (is.null(Z)) 0L else ncol(Z)
  dims <- sprintf("n=%d, pF=%d, pX=%d, pZ=%d", nrow(F), ncol(F), ncol(X), ncolZ)
  if (length(numericRefUsed)) {
    message("Numeric anchors (auto): ",
            paste(sprintf("%s=%s", names(numericRefUsed),
                          format(unlist(numericRefUsed), digits=6)), collapse = ", "))
  }
  message("Design split: ", dims)
  message("X columns: ", prettyJoin(colnames(X)))
  message("Z columns: ", prettyJoin(if (is.null(Z)) character(0) else colnames(Z)))
  nz <- cF[abs(cF) > tol]; nz <- nz[order(abs(nz), decreasing = TRUE)]
  if (verbosity == "debug" && length(nz)) {
    top <- head(nz, 30)
    message("Non-zero contrast weights (top 30 by |weight|):")
    message(paste(sprintf("  %-30s % .6g", names(top), unclass(top)), collapse = "\n"))
    if (length(nz) > 30) message(sprintf("  ... (+%d more)", length(nz) - 30))
  }
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

