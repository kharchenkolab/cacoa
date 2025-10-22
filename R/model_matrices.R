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
#' contrast <- c("group","B","A")  # μ(group=B) − μ(group=A)
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
#' @param contrast A linear contrast (see “Supported contrast specifications”).
#' @param formula RHS formula for the design (`~ ...`). If `NULL`, a default
#'   is constructed as described above.
#' @param na.action NA handler for `model.frame`.
#' @param numericRef `"auto"` or named list of numeric anchors (e.g., `list(age=35)`).
#' @param numericRefRows Optional row indices / logical mask to compute mean anchors.
#' @param tol Numeric; threshold for selecting non-zero contrast columns into `X`.
#' @param tolRow Per-row activity threshold used to compute `core.rows`.
#' @param validate Logical; compute diagnostics (rank, aliasing, VIF, permutation checks).
#' @param verbosity `"none"|"warn"|"info"|"debug"`.
#' @param computeQrZ Logical; if `TRUE` return `qrZ` of `Z` for Freedman–Lane.
#' @param blockVars Optional character vector of factor names to define permutation blocks.
#' @param buildBlocks Logical; if `TRUE` compute `blocks`, permutation groups, and diagnostics.
#'
#' @return A list with:
#' \item{F}{Full model matrix (with attributes `terms`, `xlevels`, `contrasts`).}
#' \item{X}{Core submatrix (or `F` if no nuisance).}
#' \item{Z}{Nuisance submatrix (or `NULL` if none).}
#' \item{contrast.F, contrast.X}{Contrast vectors aligned to `F` and `X`.}
#' \item{core.rows}{Logical mask of rows with non-negligible activity in `X`.}
#' \item{qrZ}{QR decomposition of `Z` for Freedman–Lane (or `NULL`).}
#' \item{blocks, perm.groups, diagnostics}{If `buildBlocks=TRUE`, auxiliary info.}
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
buildDesignMatrices <- function(data, contrast, formula = NULL,
                                na.action = stats::na.pass,
                                numericRef = "auto",
                                numericRefRows = NULL,
                                tol = 1e-12,
                                tolRow = sqrt(.Machine$double.eps),
                                validate = TRUE,
                                verbosity = c("none","warn","info","debug"),
                                computeQrZ = TRUE,
                                blockVars = NULL,
                                buildBlocks = FALSE) {
  verbosity <- match.arg(verbosity)
  
  # Default formula
  if (is.null(formula)) {
    formula <- buildDefaultFormula(data)
    if (verbosity %in% c("info","debug")) {
      message("No formula supplied; using: ", deparse(formula))
    }
  }
  
  # Prune non-varying terms
  formula_used <- pruneFormulaByData(formula, data, na.action = na.action, verbosity = verbosity)
  
  # Pick baselines so contrasted levels are not dropped (when intercept is present)
  baselines <- chooseBaselinesForSpec(formula_used, data, contrast)
  
  # Build F
  F <- buildFullDesign(formula_used, data, na.action = na.action, baselines = baselines)
  
  # Synthetic contrast over F
  cF <- buildSyntheticContrast(F, data, contrast,
                               numericRef = numericRef,
                               numericRefRows = numericRefRows)
  
  # Split by contrast (+ auto-promotion when no nuisance)
  sp <- splitByContrast(F, cF, tol = tol, tolRow = tolRow, promoteIfNoNuisance = TRUE)
  X <- sp$X; Z <- sp$Z
  
  # qrZ for FL
  qrZ <- if (computeQrZ && !is.null(Z)) qr(as.matrix(Z)) else NULL
  
  # Optional blocks & permutation groups
  blocks <- NULL; perm.groups <- NULL
  if (buildBlocks) {
    spec <- try(normalizeContrastSpec(contrast), silent = TRUE)
    nuis.fac <- deriveNuisanceFactors(formula_used, data,
                                      contrastSpec = if (!inherits(spec, "try-error")) spec else NULL,
                                      explicit = blockVars)
    blocks <- makeBlocks(data, nuisance = nuis.fac,
                         block.vars = if (!is.null(blockVars)) blockVars else NULL)
    perm.groups <- permutationGroups(blocks = blocks,
                                     core.rows = if (is.null(sp$core.rows)) rep(TRUE, nrow(F)) else sp$core.rows)
  }
  
  # Diagnostics
  diag <- NULL
  if (validate) {
    spec <- try(normalizeContrastSpec(contrast), silent = TRUE)
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
  
  spec_out <- try(normalizeContrastSpec(contrast), silent = TRUE)
  
  list(
    F = F, X = X, Z = Z,
    contrast.F = sp$contrast.F,
    contrast.X = sp$contrast.X,
    core.rows = sp$core.rows,
    qrZ = qrZ,
    blocks = blocks,
    perm.groups = perm.groups,
    diagnostics = diag,
    numeric_ref_used = attr(cF, "numeric_ref_used") %||% list(),
    formula_used = formula_used,
    contrast_spec = if (!inherits(spec_out, "try-error")) spec_out else NULL,
    baselines_used = baselines
  )
}

#' Prepare pair-level design matrices from sample metadata and (inferred or explicit) paired contrast
#'
#' @description
#' **`buildPairDesignMatrices()`** lifts a sample-level metadata table to the
#' **pair level** by forming all unordered sample pairs \eqn{(i<j)}, creates a
#' pair-level metadata frame, and then builds a pair-level design using the same
#' machinery as `buildDesignMatrices()`.
#'
#' Concretely, it:
#' 1. Generates **all unordered pairs** (i.e., \eqn{n \choose 2} rows).
#' 2. Builds **pair-level metadata**:
#'    - If a *focus factor* is used (see below), a composition factor
#'      `"<focusVar>_pair"` with levels `cross` and `L:L` (one per level `L` of the focus factor).
#'      `cross` denotes pairs from different levels; `L:L` denotes same-level pairs.
#'    - For every numeric sample covariate, two **formula-safe** columns:
#'      `pair_<var>_mean` and `pair_<var>_diff` (with `diff = value[i] - value[j]`).
#' 3. Chooses a **paired contrast** (either **inferred** from the sample-level
#'    contrast or **explicitly provided**) and builds the design for the pair meta.
#'
#' ### What happens if `formula` is not supplied?
#' This call uses the **same defaulting logic as** `buildDesignMatrices()`:
#' if any factor appears in the pair meta (e.g., `group_pair`), a **saturated**
#' design is used (`~ 0 + ...`); otherwise numeric-only gets an intercept design.
#' Non-varying terms are pruned safely.
#'
#' ### How paired contrasts are obtained
#' - If `pairContrast` is **missing**, we try to **infer** the focus factor and
#'   ALT/REF levels from the **sample-level contrast**. Supported sample specs:
#'   *DESeq2 triple* `c("group","B","A")`, or `list(type="simple"| "marginal",
#'   term="group", num="B", den="A", ...)`. (A marginal sample contrast is treated
#'   like the corresponding simple `B vs A` at the pair level.)
#'
#'   With a focus factor `{ref, alt}`, we seed one of the standard **composition**
#'   contrasts over the `group_pair` factor:
#'   - **shift**: `cross − ½(REF:REF + ALT:ALT)`
#'   - **total**: `cross − REF:REF`
#'   - **var**:   `ALT:ALT − REF:REF`
#'
#' - If that inference is **not possible** (e.g., the sample contrast targets an
#'   interaction or a numeric-only slope), the function **errors** and prints the
#'   list of candidate pair columns so the user can specify `pairContrast` explicitly.
#'
#' ### Explicit paired contrasts
#' You can pass `pairContrast` explicitly in two ways:
#'
#' 1) **Named numeric vector** over pair-design coefficient names. Examples:
#' ```r
#' pairContrast <- c("group_paircross"=1, "group_pairA:A"=-0.5, "group_pairB:B"=-0.5)
#' pairContrast <- c("pair[age]_diff"=0.2, "pair[age]_mean"=0)  # friendly names OK
#' ```
#' Friendly numeric aliases like `pair[age]_diff` are normalized internally to
#' formula-safe names (e.g., `pair_age_diff`).
#'
#' 2) **List**: `list(var=, cells=, coefs=)` where:
#' - `var` is the focus factor (e.g., `"group"`),
#' - `cells` is a named numeric over `{'cross','L:L',...}`,
#' - `coefs` (optional) is a named numeric over other pair coefficients
#'   (e.g., `pair_age_diff`).
#'
#'
#' @param sample.meta A `data.frame` of **sample-level** covariates.
#' @param sampleContrast The **sample-level** contrast used to infer the paired
#'   focus factor and ALT/REF levels when `pairContrast` is not provided
#'   (supports DESeq2 triple and single-factor `type="simple"`/`"marginal"`).
#' @param dist.type One of `"shift"`, `"total"`, or `"var"`; used only when the
#'   paired contrast is inferred (i.e., `pairContrast` is `NULL`).
#' @param pairContrast Optional **explicit** paired contrast:
#'   - named numeric vector over pair-design coefficients, or
#'   - `list(var=, cells=, coefs=)` as described above.
#' @param na.action NA handling for `model.frame`.
#' @param verbosity `"none"|"warn"|"info"`.
#'
#' @examples
#' # 1) Infer paired contrast from sample triple (default SHIFT composition):
#' # p <- buildPairDesignMatrices(df, sampleContrast=c("group","B","A"), dist.type="shift")
#'
#' # 2) Explicit paired contrast via cells + numeric diff:
#' # pc <- list(var="group", cells=c("A:A"=-1, "C:C"=1), coefs=c("pair[age]_diff"=0.25))
#' # p  <- buildPairDesignMatrices(df, sampleContrast=c("group","B","A"), pairContrast=pc)
#'
#' # 3) Numeric-only paired contrast:
#' # p <- buildPairDesignMatrices(df, sampleContrast=c("group","B","A"),
#' #                              pairContrast=c("pair[age]_diff"=1))
#'
#' @return a list in the same format as the return of `buildDesignMatrices(pair.meta, contrast=...)` that also includes:
#' - `pairs`: the 2-column matrix of (i, j) indices (i < j),
#' - `pair.meta`: the constructed pair-level metadata,
#'
#' @export
buildPairDesignMatrices <- function(sample.meta, sampleContrast,
                                    dist.type = c("shift","total","var"),
                                    pairContrast = NULL,
                                    na.action = stats::na.pass,
                                    verbosity = c("none","warn","info")) {
  stopifnot(is.data.frame(sample.meta))
  dist.type <- match.arg(dist.type)
  verbosity <- match.arg(verbosity)
  
  # (1) All unordered pairs
  n <- nrow(sample.meta)
  pairs <- lowerTriIndices(n)
  if (!nrow(pairs)) stop("Not enough samples to form pairs (need at least 2).")
  
  # (2) Determine focus factor and/or explicit intent
  focusVar <- NULL
  explicitCoef  <- NULL  # final named-numeric vector
  explicitCells <- NULL  # list(var=,cells=,coefs=)
  
  if (!is.null(pairContrast)) {
    if (is.numeric(pairContrast) && !is.null(names(pairContrast))) {
      focusVar <- tryCatch(extractFocusVarFromNamedContrast(names(pairContrast)),
                           error = function(e) stop(conditionMessage(e), call. = FALSE))
    } else if (is.list(pairContrast)) {
      if (!is.null(pairContrast$cells)) {
        if (is.null(pairContrast$var))
          stop("Explicit list pairContrast requires 'var' when 'cells' are provided.", call. = FALSE)
        if (!is.numeric(pairContrast$cells) || is.null(names(pairContrast$cells)))
          stop("'cells' must be a named numeric over {'cross','<L>:<L>'}.", call. = FALSE)
        focusVar <- as.character(pairContrast$var)
        explicitCells <- pairContrast
      } else if (!is.null(pairContrast$coefs)) {
        if (!is.numeric(pairContrast$coefs) || is.null(names(pairContrast$coefs)))
          stop("'coefs' must be a named numeric vector.", call. = FALSE)
        explicitCoef <- pairContrast$coefs
      } else {
        stop("Unsupported 'pairContrast' list: provide either 'cells' (and 'var') or 'coefs'.", call. = FALSE)
      }
    } else {
      stop("Unsupported 'pairContrast' format. Use a named numeric vector, or list(var=,cells=,coefs=).", call. = FALSE)
    }
  }
  
  # If still no focus/explicit, infer from the sample-level contrast
  if (is.null(focusVar) && is.null(explicitCells) && is.null(explicitCoef)) {
    tgt <- tryCatch(inferPairTargetFromSampleContrast(sampleContrast, sample.meta),
                    error = function(e) e)
    if (inherits(tgt, "error")) {
      stop(paste0("Cannot infer paired contrast from the sample-level contrast: ",
                  conditionMessage(tgt), "\n",
                  "Provide 'pairContrast' or use a single-factor simple/marginal sample contrast.\n",
                  "Candidate pair columns by factor:\n  ",
                  candidatePairColumnsMsg(sample.meta)), call. = FALSE)
    }
    focusVar <- tgt$var
  }
  
  # (3) Pairify metadata
  pair.meta <- pairifyMeta(sample.meta, pairsIdx = pairs, focusVar = focusVar)
  
  # (4) Build explicit coefficient vector from cells/coefs if provided
  if (!is.null(explicitCells)) {
    coef_cells <- buildCoefFromCells(pair.meta, var = focusVar, cells = explicitCells$cells)
    explicitCoef <- coef_cells
    if (!is.null(explicitCells$coefs)) {
      explicitCoef <- c(explicitCoef, normalizeExplicitPairCoefNames(explicitCells$coefs, pair.meta, focusVar))
    }
  }
  
  # (5) Normalize names for direct named-numeric vector
  if (!is.null(pairContrast) && is.numeric(pairContrast) && !is.null(names(pairContrast))) {
    explicitCoef <- normalizeExplicitPairCoefNames(pairContrast, pair.meta, focusVar)
  }
  
  # (6) If no explicit vector, synthesize defaults from inferred sample contrast
  if (is.null(explicitCoef)) {
    tgt <- inferPairTargetFromSampleContrast(sampleContrast, sample.meta)  # validated
    if (!tgt$var %in% names(sample.meta))
      stop("Target factor '", tgt$var, "' not found in metadata.", call. = FALSE)
    f.levels <- levels(factor(sample.meta[[tgt$var]]))
    miss <- setdiff(c(tgt$ref, tgt$alt), f.levels)
    if (length(miss))
      stop("ALT/REF levels not found for factor '", tgt$var, "': ", paste(miss, collapse = ", "), call. = FALSE)
    
    pair.var    <- paste0(tgt$var, "_pair")
    if (!pair.var %in% names(pair.meta))
      stop("Internal error: '", pair.var, "' not found in pair.meta.", call. = FALSE)
    pair.levels <- levels(pair.meta[[pair.var]])              # e.g., cross, A:A, B:B
    w_by_level  <- pairContrastWeights(ref = tgt$ref, alt = tgt$alt,
                                       dist.type = dist.type, level.order = pair.levels)
    names(w_by_level) <- paste0(pair.var, pair.levels)        # e.g., group_paircross
    explicitCoef <- w_by_level
  } else {
    if (is.null(names(explicitCoef)) || !is.numeric(explicitCoef))
      stop("Explicit 'pairContrast' must be a named numeric vector.", call. = FALSE)
  }
  
  # (7) Build paired design using the sample-level constructor on pair meta
  out <- buildDesignMatrices(
    data        = pair.meta,
    contrast    = explicitCoef,   # named numeric over coefficient names
    formula     = NULL,           # default builder saturates if factor present
    na.action   = na.action,
    numericRef  = "auto",
    tol         = 1e-12,
    tolRow      = sqrt(.Machine$double.eps),
    validate    = TRUE,
    verbosity   = switch(verbosity, none="none", warn="warn", info="info"),
    computeQrZ  = TRUE,
    buildBlocks = FALSE
  )
  
  if (verbosity %in% c("info")) {
    if (is.null(pairContrast)) {
      message(sprintf("Paired design (auto): factor='%s', dist.type='%s'.",
                      focusVar %||% "(none)", dist.type))
    } else {
      message("Paired design: explicit pair-level contrast applied.")
    }
  }
  
  out$pairs <- pairs
  out$pair.meta <- pair.meta
  return(out)
}

# =====================================================================
# Helpers (shared by the two public constructors)
# =====================================================================

# ---- Defaults & pruning ----

buildDefaultFormula <- function(data) {
  stopifnot(is.data.frame(data))
  vars <- names(data)
  if (!length(vars)) return(as.formula("~ 1"))
  anyFac <- any(vapply(data, function(x) is.factor(x) || is.character(x), logical(1)))
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

prettyJoin <- function(x, maxShow = 12) {
  x <- as.character(x)
  if (!length(x)) return("(none)")
  if (length(x) <= maxShow) return(paste(x, collapse = ", "))
  paste0(paste(x[1:maxShow], collapse = ", "), sprintf(", … (+%d more)", length(x) - maxShow))
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

chooseBaselinesForSpec <- function(formula, data, contrast) {
  spec <- try(normalizeContrastSpec(contrast), silent = TRUE)
  if (inherits(spec, "try-error")) return(list())
  
  trm <- stats::terms(if (inherits(formula,"formula")) formula else as.formula(formula), data = data)
  if (attr(trm, "intercept", exact = TRUE) == 0L) return(list())
  
  vars_in_model <- all.vars(trm)
  levels_in_data <- lapply(intersect(vars_in_model, names(data)), function(v) {
    x <- data[[v]]
    if (is.factor(x) || is.character(x)) levels(droplevels(factor(x))) else NULL
  })
  names(levels_in_data) <- intersect(vars_in_model, names(data))
  
  isFacVar  <- function(v) !is.null(levels_in_data[[v]])
  varLevels <- function(v) levels_in_data[[v]] %||% character(0)
  
  chooseNot <- function(v, avoid) {
    L <- varLevels(v); if (!length(L)) return(NULL)
    if (!length(avoid) || !(avoid %in% L)) return(L[1])
    alt <- setdiff(L, avoid); if (length(alt)) alt[1] else L[1]
  }
  
  bases <- list()
  
  if (is.list(spec) && spec$type == "simple" && !grepl(":", spec$term, fixed = TRUE)) {
    v <- spec$term
    if (isFacVar(v)) {
      if (spec$den %in% varLevels(v)) bases[[v]] <- spec$den
      else { b <- chooseNot(v, spec$num); if (!is.null(b)) bases[[v]] <- b }
    }
    return(bases)
  }
  
  if (is.list(spec) && spec$type == "simple" && grepl(":", spec$term, fixed = TRUE)) {
    vs <- parseTermVars(spec$term)
    numC <- parseCell(spec$num, vs)
    for (v in vs) if (isFacVar(v)) { b <- chooseNot(v, numC[[v]]); if (!is.null(b)) bases[[v]] <- b }
    return(bases)
  }
  
  if (is.list(spec) && spec$type == "lincomb") {
    vs <- parseTermVars(spec$term)
    w  <- spec$cells
    pick <- if (any(w > 0)) names(w)[which.max(w)] else names(w)[which.max(w)]
    if (length(pick)) {
      numC <- parseCell(pick, vs)
      for (v in vs) if (isFacVar(v)) { b <- chooseNot(v, numC[[v]]); if (!is.null(b)) bases[[v]] <- b }
    }
    return(bases)
  }
  
  if (is.list(spec) && spec$type == "marginal") {
    v <- spec$term
    if (isFacVar(v)) {
      if (spec$den %in% varLevels(v)) bases[[v]] <- spec$den
      else { b <- chooseNot(v, spec$num); if (!is.null(b)) bases[[v]] <- b }
    }
    return(bases)
  }
  
  list()
}

buildFullDesign <- function(formula, data, na.action = stats::na.pass,
                            contrasts.arg = NULL, baselines = NULL) {
  rhs <- if (inherits(formula, "formula")) formula else as.formula(formula)
  mf  <- stats::model.frame(rhs, data, na.action = na.action)
  
  if (length(baselines)) {
    for (v in names(baselines)) {
      if (!is.null(mf[[v]]) && (is.factor(mf[[v]]) || is.character(mf[[v]]))) {
        lev <- as.character(baselines[[v]])
        if (length(lev) && lev %in% levels(droplevels(factor(mf[[v]])))) {
          mf[[v]] <- stats::relevel(droplevels(factor(mf[[v]])), ref = lev)
        }
      }
    }
  }
  
  trm <- stats::terms(rhs, data = data)
  F   <- stats::model.matrix(trm, mf, contrasts.arg = contrasts.arg)
  attr(F, "terms")   <- trm
  attr(F, "xlevels") <- lapply(mf, function(x) if (is.factor(x)) levels(x) else NULL)
  attr(F, "contrasts") <- attr(F, "contrasts")
  F
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

.oneRowFromFormula <- function(trm, xlevels, contrastsArg, targetCols,
                               data, at = list(), numericRef = NULL) {
  mf <- stats::model.frame(trm, data, na.action = stats::na.pass)
  vars <- names(mf)
  if (!is.null(attr(trm,"response")) && attr(trm,"response") > 0) {
    resp <- all.vars(stats::delete.response(trm))[1]
    vars <- setdiff(vars, resp)
  }
  newd <- vector("list", length(vars)); names(newd) <- vars
  for (v in vars) {
    if (!is.null(xlevels[[v]])) newd[[v]] <- factor(xlevels[[v]][1], levels = xlevels[[v]])
    else                        newd[[v]] <- (numericRef[[v]] %||% 0)
  }
  for (v in names(at)) {
    if (!v %in% vars) stop(sprintf("Variable '%s' not in model terms.", v))
    if (!is.null(xlevels[[v]])) {
      levs <- xlevels[[v]]; val <- as.character(at[[v]])
      if (!val %in% levs) stop(sprintf("Level '%s' not in levels(%s): {%s}", val, v, paste(levs, collapse=", ")))
      newd[[v]] <- factor(val, levels = levs)
    } else newd[[v]] <- at[[v]]
  }
  nd  <- as.data.frame(newd)
  mf1 <- stats::model.frame(trm, nd, xlev = xlevels, na.action = stats::na.pass)
  r1  <- stats::model.matrix(trm, mf1, contrasts.arg = contrastsArg)
  miss <- setdiff(targetCols, colnames(r1))
  if (length(miss)) r1 <- cbind(r1, matrix(0, 1, length(miss), dimnames = list(NULL, miss)))
  drop(r1[1, targetCols, drop = FALSE])
}

buildSyntheticContrast <- function(F, data, contrast,
                                   numericRef = "auto",
                                   numericRefRows = NULL,
                                   stopOnAmbiguousTriple = TRUE) {
  spec <- normalizeContrastSpec(contrast)
  trm  <- attr(F, "terms"); xlv <- attr(F, "xlevels"); ctr <- attr(F, "contrasts")
  if (is.null(trm) || is.null(xlv))
    stop("F must carry 'terms' and 'xlevels' (use buildFullDesign()).")
  
  tlabels <- attr(trm, "term.labels")
  hasIntWith <- function(var) any(grepl(":", tlabels) & grepl(paste0("(^|:)", var, "(:|$)"), tlabels))
  
  numRef <- resolveNumericRef(F, data, contrast,
                              numericRef = numericRef,
                              numericRefRows = numericRefRows)
  one <- function(at) .oneRowFromFormula(trm, xlv, ctr, colnames(F), data, at, numRef)
  
  cF <- setNames(numeric(ncol(F)), colnames(F))
  
  if (spec$type == "coef") {
    miss <- setdiff(names(spec$coefs), colnames(F))
    if (length(miss)) stop("Unknown coefficient(s) in contrast: ", paste(miss, collapse = ", "))
    cF[names(spec$coefs)] <- as.numeric(spec$coefs)
    attr(cF, "numeric_ref_used") <- numRef
    return(cF)
  }
  
  if (identical(spec$plain_triple, TRUE) && stopOnAmbiguousTriple && hasIntWith(spec$term)) {
    stop(sprintf("Ambiguous triple: interactions with '%s' present. Use type='simple' (add at=) or type='marginal' (add over=).", spec$term))
  }
  
  if (spec$type == "lincomb") {
    vars <- parseTermVars(spec$term)
    for (nm in names(spec$cells)) {
      at_i <- utils::modifyList(spec$at %||% list(), parseCell(nm, vars))
      cF <- cF + as.numeric(spec$cells[[nm]]) * one(at_i)
    }
    attr(cF, "numeric_ref_used") <- numRef
    return(cF)
  }
  
  if (spec$type == "simple") {
    if (grepl(":", spec$term, fixed = TRUE)) {
      vars  <- parseTermVars(spec$term)
      atNum <- utils::modifyList(spec$at, parseCell(spec$num, vars))
      atDen <- utils::modifyList(spec$at, parseCell(spec$den, vars))
      cF <- one(atNum) - one(atDen)
    } else {
      atN <- utils::modifyList(spec$at, setNames(list(spec$num), spec$term))
      atD <- utils::modifyList(spec$at, setNames(list(spec$den), spec$term))
      cF <- one(atN) - one(atD)
    }
    attr(cF, "numeric_ref_used") <- numRef
    return(cF)
  }
  
  if (spec$type == "marginal") {
    ov <- spec$over
    xlv <- attr(F, "xlevels")
    levs <- lapply(ov, function(v) { L <- xlv[[v]]; if (is.null(L)) stop(sprintf("'%s' in over= is not a factor.", v)); L })
    names(levs) <- ov
    grid <- do.call(expand.grid, c(levs, stringsAsFactors = FALSE))
    
    ws <- spec$weights
    if (is.character(ws)) {
      if (!ws %in% c("equal","proportional"))
        stop("weights must be 'equal', 'proportional', or named numeric.")
      w <- rep(1/nrow(grid), nrow(grid))
      if (ws == "proportional") {
        tab <- data[, ov, drop = FALSE]
        for (v in ov) tab[[v]] <- factor(tab[[v]], levels = levs[[v]])
        idx <- do.call(interaction, c(tab, drop=TRUE, sep=":"))
        key <- apply(grid, 1, function(r) paste(r, collapse=":"))
        cnt <- tapply(rep(1, nrow(tab)), idx, sum)
        w   <- as.numeric(cnt[key]); w[is.na(w)] <- 0
        if (sum(w) == 0) w <- rep(1/nrow(grid), nrow(grid)) else w <- w / sum(w)
      }
    } else if (is.numeric(ws)) {
      key <- apply(grid, 1, function(r) paste(r, collapse=":"))
      if (is.null(names(ws))) stop("Numeric weights must be named by: ", paste(key, collapse=", "))
      w <- as.numeric(ws[key]); if (any(is.na(w)) || any(w < 0) || sum(w) == 0) stop("Invalid numeric weights.")
      w <- w / sum(w)
    } else stop("Unsupported weights spec.")
    for (i in seq_len(nrow(grid))) {
      at_i  <- as.list(grid[i,,drop=FALSE])
      atN <- utils::modifyList(at_i, spec$at); atD <- atN
      atN[[spec$term]] <- spec$num; atD[[spec$term]] <- spec$den
      cF <- cF + w[i] * (one(atN) - one(atD))
    }
    attr(cF, "numeric_ref_used") <- numRef
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

# ---- Old-style helpers (blocks & FL plumbing) ----

makeBlocks <- function(meta, nuisance, block.vars = NULL) {
  if (!is.null(block.vars) && length(block.vars)) {
    return(interaction(lapply(meta[block.vars], as.factor), drop = TRUE, lex.order = TRUE))
  }
  fac.nuis <- if (length(nuisance)) {
    nuisance[vapply(meta[nuisance], function(x) is.factor(x) || is.character(x), logical(1))]
  } else character(0)
  if (length(fac.nuis)) interaction(lapply(meta[fac.nuis], as.factor), drop = TRUE, lex.order = TRUE)
  else factor(rep("all", nrow(meta)))
}

filterNuisance <- function(meta, nuisance.names) {
  if (!length(nuisance.names)) return(character(0))
  keep <- vapply(nuisance.names, function(v) {
    x <- meta[[v]]
    if (is.factor(x) || is.character(x)) length(levels(droplevels(factor(x)))) >= 2L
    else if (is.numeric(x)) length(unique(x[!is.na(x)])) >= 2L
    else FALSE
  }, logical(1))
  nuisance.names[keep]
}

permutationGroups <- function(blocks, core.rows = NULL) {
  stopifnot("blocks is not a factor"=is.factor(blocks), "no randomization blocks found"=length(blocks) >= 1L)
  n <- length(blocks)
  if (is.null(core.rows)) {
    core.rows <- rep(TRUE,n)
    groups.core <- groups.full <- split(seq_len(n), droplevels(blocks), drop = TRUE)
  } else {
    stopifnot("core.rows is not logi"=is.logical(core.rows), "core.rows length mismatch"=length(core.rows) == n)
    groups.full <- split(seq_len(n)[core.rows], droplevels(blocks[core.rows]), drop = TRUE)
    core.idx    <- which(core.rows)
    blocks.core <- droplevels(blocks[core.rows])
    pos.in.core <- integer(n); pos.in.core[core.idx] <- seq_along(core.idx)
    groups.core <- lapply(split(core.idx, blocks.core, drop = TRUE),
                          function(ids_full) pos.in.core[ids_full])
  }
  list(full = groups.full, core = groups.core)
}

# ---- FL residualization ----

#' @keywords internal
residualizeForFL <- function(y, qrZ, X) {
  X <- as.matrix(X)
  if (is.null(qrZ)) return(list(y.r = y, X.r = X))
  y.r <- qr.resid(qrZ, y)
  X.r <- qr.resid(qrZ, X)
  list(y.r = y.r, X.r = X.r)
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
  
  perm <- NULL
  if (!is.null(blocks)) {
    stopifnot(is.factor(blocks), nrow(F) == length(blocks))
    groups.core <- permutationGroups(blocks, if (is.null(core.rows)) rep(TRUE, nrow(F)) else core.rows)$core
    
    permCoreSummary <- function(meta, blocks, core.rows, contrastSpec, groups.core, block.factors, show.top) {
      labs   <- names(groups.core); if (is.null(labs)) labs <- as.character(seq_along(groups.core))
      n.core <- vapply(groups.core, length, integer(1))
      has.fac <- !is.null(contrastSpec) && contrastSpec$type == "simple" &&
        !grepl(":", contrastSpec$term, fixed = TRUE) &&
        !is.null(meta) && (contrastSpec$term %in% names(meta))
      n.levels <- integer(length(groups.core)); comp.str <- character(length(groups.core))
      if (has.fac) {
        keep <- c(contrastSpec$num, contrastSpec$den)
        for (i in seq_along(groups.core)) {
          ids <- which(core.rows)[groups.core[[i]]]
          v   <- as.character(meta[[contrastSpec$term]][ids]); v <- v[v %in% keep]
          tab <- sort(table(v), decreasing = TRUE)
          n.levels[i] <- length(tab)
          comp.str[i] <- if (length(tab)) paste(sprintf("%s:%d", names(tab), as.integer(tab)), collapse=",") else ""
        }
      } else {
        n.levels[] <- NA_integer_; comp.str[] <- ""
      }
      df <- data.frame(block = labs, n.core = n.core, n.levels = n.levels,
                       composition = comp.str, stringsAsFactors = FALSE)
      too.small <- df$n.core < thresholds$min.core.size
      no.var    <- !is.na(df$n.levels) & (df$n.levels < 2L)
      prob.idx  <- which(too.small | no.var)
      
      comb.df <- NULL
      if (length(prob.idx) && !is.null(block.factors) && length(block.factors) &&
          all(block.factors %in% names(meta))) {
        for (i in prob.idx) {
          rows.full <- which(blocks == levels(blocks)[match(df$block[i], levels(blocks))])
          if (length(rows.full)) {
            vals <- vapply(block.factors, function(v) as.character(meta[[v]][rows.full[1]]), character(1))
            comb.df <- rbind(comb.df,
                             data.frame(block = df$block[i],
                                        n.core = df$n.core[i],
                                        n.levels = df$n.levels[i],
                                        composition = df$composition[i],
                                        combo = paste(paste0(block.factors, "=", vals), collapse = ", "),
                                        stringsAsFactors = FALSE))
          }
        }
        if (!is.null(comb.df)) {
          ord <- order(comb.df$n.core, comb.df$n.levels)
          comb.df <- comb.df[ord, , drop = FALSE]
          comb.df <- head(comb.df, show.top)
        }
      }
      list(summary = df, prob.idx = prob.idx, comb.df = comb.df)
    }
    
    ps <- permCoreSummary(meta, blocks, core.rows, contrastSpec,
                          groups.core, block.factors = block.factors, show.top = thresholds$show.top)
    
    movable <- ps$summary$n.core[ps$summary$n.core >= thresholds$min.core.size]
    eff.perm.log <- if (length(movable)) sum(lfactorial(movable)) else 0
    
    perm <- list(groups.core = groups.core,
                 core.summary = ps$summary,
                 eff.perm.log = eff.perm.log,
                 problems = ps$prob.idx,
                 combos = ps$comb.df)
    
    too.small <- sum(ps$summary$n.core < thresholds$min.core.size)
    no.var    <- sum(!is.na(ps$summary$n.levels) & ps$summary$n.levels < 2L)
    if (too.small > 0)
      warns <- c(warns, sprintf("%d block(s) have < %d core rows; those rows will be frozen (no permutation).",
                                too.small, thresholds$min.core.size))
    if (no.var > 0)
      warns <- c(warns, sprintf("%d block(s) show no within-block variation in the contrasted variable; those rows will be frozen.",
                                no.var))
    if (exp(min(eff.perm.log, 50)) < thresholds$min.eff.perm)
      warns <- c(warns, sprintf("Effective permutations is small (~exp(%.1f)). Consider relaxing blocks or a wild bootstrap.",
                                eff.perm.log))
  }
  
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
    if (length(nz) > 30) message(sprintf("  … (+%d more)", length(nz) - 30))
  }
}

deriveNuisanceFactors <- function(formula, data, contrastSpec, explicit = NULL) {
  if (!is.null(explicit) && length(explicit)) {
    return(filterNuisance(data, explicit))
  }
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  allVars <- names(mf)
  isFac <- vapply(mf, function(x) is.factor(x) || is.character(x), logical(1))
  facVars <- allVars[isFac]
  coreVars <- character(0)
  if (!is.null(contrastSpec)) {
    coreVars
    