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
    formula <- buildDefaultFormula(data)
    if (verbosity %in% c("info","debug")) {
      message("No formula supplied; using: ", deparse(formula))
    }
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
  
  # qrZ for FL
  qrZ <- if (computeQrZ && !is.null(Z)) qr(as.matrix(Z)) else NULL
  
  # Optional blocks & permutation groups
 
  blocks <- NULL; perm.groups <- NULL
  if (buildBlocks) {
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
    blocks = blocks,
    perm.groups = perm.groups,
    diagnostics = diag,
    numeric_ref_used = attr(cF, "numeric_ref_used") %||% list(),
    formula_used = formula_used,
    contrast_spec = if (!inherits(spec, "try-error")) spec else NULL,
    baselines_used = baselines
  )
}

#' Build pair-level design matrices from sample metadata and an inferred or explicit paired contrast
#'
#' @description
#' **`buildPairDesignMatrices()`** lifts a sample-level metadata table to the
#' **pair level** by forming all unordered sample pairs \eqn{(i<j)}, constructs
#' a pair-level metadata frame, and then builds a pair-level design using the
#' same machinery as [`buildDesignMatrices()`].
#'
#' Concretely, it:
#' \enumerate{
#'   \item Generates **all unordered pairs** (i.e., \eqn{n \choose 2} rows);
#'   \item Builds **pair-level metadata** (see *Pair-level metadata* below);
#'   \item Chooses a **paired contrast** (either **inferred** from the sample-level
#'         design or **explicitly provided**) and builds the design for the pair meta.
#' }
#'
#' ### Pair-level metadata
#' For a chosen *focus factor* `var` (typically inferred from the sample-level
#' contrast), we create a **composition** factor `"<var>_pair"` with levels:
#' \itemize{
#'   \item `"L1:L2"` — pairs drawn from **different** focus levels of `var`;
#'   \item `"L:L"`   — one level per **same-level** class present in the sample data,
#'                     e.g. `"A:A"`, `"B:B"`, ... (only for levels seen in `data`).
#' }
#'
#' For **numeric** sample covariates we mirror (see *Numeric mirroring*):
#' \itemize{
#'   \item `pair_<safe(var)>_mean` — the average \eqn{(x_i + x_j)/2};
#'   \item `pair_<safe(var)>_diff` — the oriented difference \eqn{x_i - x_j} (with pairs `(i<j)`).
#' }
#' All names are **formula-safe**; a friendly alias like `pair[age]_diff` is
#' recognized (and mapped to `pair_age_diff`) when you pass an explicit contrast.
#'
#' ### How the paired contrast is obtained
#' The function prefers **inference** from the **sample-level design**:
#' \itemize{
#'   \item It inspects `sampleDesign$contrast_spec` (as produced by
#'         [`buildDesignMatrices()`]) and accepts **single-factor**
#'         contrasts of type `"simple"` or `"marginal"`. In both cases it extracts
#'         `{var, alt, ref}` and treats a marginal sample contrast as the corresponding
#'         *simple* \(alt vs ref\) at the **pair** level.
#'   \item With \{var, alt, ref\} in hand, it seeds one of the **composition defaults**:
#'         \describe{
#'           \item{`"shift"`}{`REF:ALT − ½(REF:REF + ALT:ALT)` (discordant vs mean of concordants)}
#'           \item{`"total"`}{`REF:ALT − REF:REF` (discordant vs reference-concordant)}
#'           \item{`"var"`}{`ALT:ALT − REF:REF` (between-class variance proxy)}
#'         }
#'   \item If inference is **not possible** (e.g., sample contrast targets an
#'         interaction `group:batch`, a numeric-only slope, a general `lincomb`,
#'         or `coef`-space weights), the function **errors** with a helpful list
#'         of candidate pair columns and asks you to pass `pairContrast` explicitly.
#' }
#'
#' You may **override** inference by supplying an explicit `pairContrast`:
#' \enumerate{
#'   \item **Named numeric vector** over pair-design coefficients (factor cells and/or
#'         numeric pair terms). Examples:
#'
#' \preformatted{
#' # composition default written explicitly
#' pairContrast <- c("group_paircross"=1, "group_pairA:A"=-0.5, "group_pairB:B"=-0.5)
#' # numeric-only contrast (friendly alias resolves to 'pair_age_diff')
#' pairContrast <- c("pair[age]_diff"=0.25)
#' }
#'
#'   \item **List form** `list(var=, cells=, coefs=)` where:
#'         \itemize{
#'           \item `var` is the focus factor name (e.g., `"group"`);
#'           \item `cells` is a named numeric over `{'cross','L:L',...}` specifying the
#'                 composition weights;
#'           \item `coefs` (optional) is a named numeric over any other pair coefficients,
#'                 e.g., `c("pair_age_diff"=0.2)`. Friendly aliases are accepted.
#'         }
#' }
#'
#' ### Formula handling (pair level)
#' If `pairFormula` is **not** supplied, the call uses the same defaulting logic
#' as [`buildDesignMatrices()`]: if any factor appears in the pair meta (e.g.,
#' `group_pair`), a **saturated** design is used (`~ 0 + ...`); otherwise a
#' numeric-only meta yields an **intercept** model. Non-varying terms are pruned safely.
#'
#' If you **do** supply `pairFormula`, it is used verbatim. Any numeric pair columns
#' referenced by `pairFormula` are auto-generated from the underlying sample numerics
#' (see *Numeric mirroring*).
#'
#' ### Numeric mirroring
#' Numeric pair terms are mirrored from the **sample-level formula** contained in
#' `sampleDesign$formula_used`. That is, we only create `pair_<var>_{mean|diff}` for
#' numeric variables that participated in the sample model (plus any additional ones
#' that `pairFormula` explicitly references).
#'
#' ### Relationship to the sample-level constructor
#' Internally, this function calls [`buildDesignMatrices()`] on the pair-level meta,
#' with your paired contrast (inferred or explicit). Therefore the returned object is
#' **identical in shape** to the sample-level output, with two extra fields:
#' \itemize{
#'   \item `pairs` — a two-column integer matrix of `(i, j)` indices (with `i < j`);
#'   \item `pair.meta` — the pair-level metadata used to build the design.
#' }
#'
#' The **split rule** for `X` vs `Z`, the **no-nuisance promotion** (promoting `X<-F`
#' when `Z` is empty or just the intercept), and the **diagnostics** are exactly the
#' same as at the sample level; see [`buildDesignMatrices()`] for details.
#'
#' @section Limitations and guidance:
#' \itemize{
#'   \item Inference requires a **single-factor** sample contrast (`"simple"` or `"marginal"`).
#'         If your sample contrast involves interactions, numeric-only slopes, or general
#'         linear combinations, specify `pairContrast` explicitly.
#'   \item Orientation of `pair_<var>_diff` is `i - j` for pairs `(i<j)`. This is consistent
#'         across the dataset and matters only if you include numeric differences in the contrast.
#'   \item Friendly aliases `pair[foo]_(mean|diff)` are accepted only in explicit contrasts;
#'         column names in the resulting design will always be **safe** (`pair_foo_mean`, etc.).
#' }
#'
#' @param sample.meta A `data.frame` of **sample-level** covariates.
#' @param sampleDesign **Required.** The result of [`buildDesignMatrices()`] at the
#'   sample level; we read `$contrast_spec` (for `{var, alt, ref}` inference) and
#'   `$formula_used` (to mirror numeric variables).
#' @param dist.type One of `"shift"`, `"total"`, or `"var"`; used only when
#'   the paired contrast is **inferred** (i.e., `pairContrast` is `NULL`).
#' @param pairContrast Optional **explicit** paired contrast:
#'   \itemize{
#'     \item A **named numeric** vector over pair-design coefficients; or
#'     \item `list(var=, cells=, coefs=)` (see *Explicit paired contrasts*).
#'   }
#' @param pairFormula Optional **RHS** formula over pair-level columns
#'   (e.g., `~ 0 + group_pair + pair_age_mean + pair_age_diff`). If missing,
#'   a default saturated or intercept model is used (see *Formula handling*).
#' @param na.action NA handling for `model.frame`/`model.matrix`.
#' @param verbosity `"none"` | `"warn"` | `"info"`; controls messages/warnings.
#'
#' @return A list identical to [`buildDesignMatrices()`] with two additional members:
#' \describe{
#'   \item{F}{full pair-level model matrix (with attributes `terms`, `xlevels`, `contrasts`).}
#'   \item{X}{core submatrix (columns where \eqn{|contrast.F| > tol}), or \code{F} when no nuisance.}
#'   \item{Z}{nuisance submatrix (complement of \code{X}), or \code{NULL} when no nuisance.}
#'   \item{contrast.F}{named numeric paired contrast over \code{colnames(F)}.}
#'   \item{contrast.X}{same, restricted to \code{colnames(X)}.}
#'   \item{core.rows}{logical mask of rows with non-negligible signal in \code{X}.}
#'   \item{qrZ}{QR decomposition of \code{Z} (or \code{NULL} if no nuisance).}
#'   \item{blocks, perm.groups, diagnostics, numeric_ref_used, formula_used,
#'         contrast_spec, baselines_used}{as in the sample-level constructor.}
#'   \item{pairs}{two-column integer matrix of `(i, j)` indices (with `i < j`).}
#'   \item{pair.meta}{the pair-level metadata used to build the design.}
#' }
#'
#' @examples
#' ## Sample-level: build a simple design and pass it in
#' # s_ctr <- c("group","B","A")                         # DESeq2 triple
#' # s     <- buildDesignMatrices(data=df, contrast=s_ctr, formula=~ group + age)
#'
#' ## [1] Default SHIFT paired contrast inferred from sample-level {group, B vs A}
#' # p1 <- buildPairDesignMatrices(df, sampleDesign = s, dist.type = "shift")
#' # names(p1$model$contrast.X)
#' # # -> c("group_paircross","group_pairA:A","group_pairB:B")
#'
#' ## [2] TOTAL and VAR variants
#' # p2 <- buildPairDesignMatrices(df, sampleDesign = s, dist.type = "total")
#' # p3 <- buildPairDesignMatrices(df, sampleDesign = s, dist.type = "var")
#'
#' ## [3] Explicit composition contrast (same as SHIFT above)
#' # pairCtr <- c("group_paircross"=1, "group_pairA:A"=-0.5, "group_pairB:B"=-0.5)
#' # p4 <- buildPairDesignMatrices(df, sampleDesign = s, pairContrast = pairCtr)
#'
#' ## [4] Numeric-only explicit contrast with a friendly alias
#' # p5 <- buildPairDesignMatrices(df, sampleDesign = s,
#' #        pairContrast = c("pair[age]_diff"=0.25))
#'
#' ## [5] List form: composition cells + extra numeric coefficient
#' # p6 <- buildPairDesignMatrices(df, sampleDesign = s,
#' #        pairContrast = list(
#' #          var   = "group",
#' #          cells = c(cross=1, "A:A"=-0.5, "B:B"=-0.5),
#' #          coefs = c("pair_age_diff"=0.1)
#' #        ))
#'
#' @seealso [buildDesignMatrices()] for the sample-level constructor and
#'   details on splitting `X` vs `Z`, numeric anchors, and diagnostics.
buildPairDesignMatrices <- function(sample.meta,
                                    sampleDesign,
                                    dist.type = c("shift","total","var"),
                                    pairContrast = NULL,
                                    pairFormula  = NULL,
                                    # NEW:
                                    pairBlockVars = NULL,
                                    buildBlocks   = TRUE,
                                    na.action = stats::na.pass,
                                    verbosity = c("none","warn","info")) {
  stopifnot(is.data.frame(sample.meta))
  if (missing(sampleDesign) || is.null(sampleDesign) || is.null(sampleDesign$formula_used)) {
    stop("sampleDesign with a valid $formula_used is required.", call. = FALSE)
  }
  dist.type <- match.arg(dist.type)
  verbosity <- match.arg(verbosity)
  
  ## 1) All unordered pairs
  n <- nrow(sample.meta)
  pairs <- {
    if (n < 2) stop("Not enough samples to form pairs (need at least 2).")
    idx <- which(lower.tri(matrix(NA, n, n)), arr.ind = TRUE)
    colnames(idx) <- c("i","j")
    idx
  }
  
  ## 2) Figure out the focus factor and any explicit contrast info
  focusVar      <- NULL
  explicitCoef  <- NULL
  explicitCells <- NULL
  
  # User-provided pairContrast
  if (!is.null(pairContrast)) {
    if (is.numeric(pairContrast) && !is.null(names(pairContrast))) {
      # try to extract focus from coefficient names like "Group_pair_Group1_and_Group2"
      focusVar <- tryCatch(
        extractFocusVarFromNamedContrast(names(pairContrast)),
        error = function(e) stop(conditionMessage(e), call. = FALSE)
      )
      # we'll normalize names below, after pair.meta exists
      explicitCoef <- pairContrast
    } else if (is.list(pairContrast)) {
      if (!is.null(pairContrast$cells)) {
        if (is.null(pairContrast$var))
          stop("Explicit list pairContrast requires 'var' when 'cells' are provided.", call. = FALSE)
        if (!is.numeric(pairContrast$cells) || is.null(names(pairContrast$cells)))
          stop("'cells' must be a named numeric vector over {'L1_and_L1','L1_and_L2',...}.", call. = FALSE)
        focusVar      <- as.character(pairContrast$var)
        explicitCells <- pairContrast  # we'll turn this into explicitCoef later
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
  
  # If still no focusVar, infer from sampleDesign$contrast_spec
  if (is.null(focusVar) && is.null(explicitCells) && is.null(explicitCoef)) {
    inferPairTargetFromDesign_local <- function(sampleDesign, data) {
      cs <- sampleDesign$contrast_spec
      if (is.null(cs)) stop("sampleDesign lacks $contrast_spec.")
      if (identical(cs$type, "simple") && !grepl(":", cs$term, fixed = TRUE))
        return(list(var = cs$term, alt = cs$num, ref = cs$den))
      if (identical(cs$type, "marginal"))
        return(list(var = cs$term, alt = cs$num, ref = cs$den))
      if (is.character(cs) && length(cs) == 3)
        return(list(var = cs[[1]], alt = cs[[2]], ref = cs[[3]]))
      stop("Cannot infer pair target from the provided sampleDesign$contrast_spec.")
    }
    tgt_try <- tryCatch(inferPairTargetFromDesign_local(sampleDesign, sample.meta),
                        error = function(e) e)
    if (inherits(tgt_try, "error")) {
      avail_vars <- names(sample.meta)
      msg <- paste0(
        "Cannot infer paired contrast from sampleDesign$contrast_spec: ",
        conditionMessage(tgt_try), "\n",
        "Either provide 'pairContrast', or use a single-factor simple/marginal sample contrast at the sample level.\n",
        "Available sample-level variables (columns) for contrasts: ",
        if (length(avail_vars)) paste(avail_vars, collapse = ", ") else "<none>"
      )
      stop(msg, call. = FALSE)
    }
    focusVar <- tgt_try$var
  }
  
  ## 3) Figure out which variables to mirror / pairify
  vf <- varsFromFormula(sampleDesign$formula_used, data = sample.meta, na.action = na.action)
  factorVarsFromSample  <- vf$factors
  numericFromSample     <- vf$numeric
  
  numericAlsoFromPairFormula <- extractPairNumericVarsFromFormula(pairFormula, sample.meta)
  numericFinal <- if (is.null(pairFormula)) {
    numericFromSample
  } else {
    unique(c(numericFromSample, numericAlsoFromPairFormula))
  }
  
  factorVarsAll <- unique(c(
    factorVarsFromSample,
    if (!is.null(focusVar)) focusVar else character(0)
  ))
  
  ## 4) Build pair.meta
  pair.meta <- pairifyMeta(sample.meta,
                           pairsIdx    = pairs,
                           factorVars  = factorVarsAll,
                           focusVar    = focusVar,
                           numericVars = numericFinal)
  
  ## 5) Turn explicitCells into explicitCoef, if needed
  if (!is.null(explicitCells)) {
    coef_cells <- buildCoefFromCells(pair.meta, var = focusVar, cells = explicitCells$cells)
    explicitCoef <- coef_cells
    if (!is.null(explicitCells$coefs)) {
      extra_norm <- normalizeExplicitPairCoefNames(explicitCells$coefs, pair.meta)
      explicitCoef <- c(explicitCoef, extra_norm)
    }
  }
  
  ## 6) If user gave a named numeric vector directly, normalize its names now
  if (!is.null(explicitCoef) && is.numeric(explicitCoef) && !is.null(names(explicitCoef))) {
    explicitCoef <- normalizeExplicitPairCoefNames(explicitCoef, pair.meta)
  }
  
  ## 7) If we STILL don't have explicitCoef, synthesize the default contrast
  if (is.null(explicitCoef)) {
    cs <- sampleDesign$contrast_spec
    if (is.null(cs)) stop("sampleDesign lacks $contrast_spec.")
    
    # extract {var, alt, ref} from simple/marginal/triple
    if (is.character(cs) && length(cs) == 3) {
      tgt <- list(var = cs[[1]], alt = cs[[2]], ref = cs[[3]])
    } else if (is.list(cs) &&
               identical(cs$type, "simple") &&
               !grepl(":", cs$term, fixed = TRUE)) {
      tgt <- list(var = cs$term, alt = cs$num, ref = cs$den)
    } else if (is.list(cs) && identical(cs$type, "marginal")) {
      tgt <- list(var = cs$term, alt = cs$num, ref = cs$den)
    } else {
      stop("Cannot synthesize default pair contrast from sampleDesign$contrast_spec.")
    }
    
    if (!tgt$var %in% names(sample.meta))
      stop("Target factor '", tgt$var, "' not found in metadata.", call. = FALSE)
    
    f.levels <- levels(factor(sample.meta[[tgt$var]]))
    miss <- setdiff(c(tgt$ref, tgt$alt), f.levels)
    if (length(miss))
      stop("ALT/REF levels not found for factor '", tgt$var, "': ",
           paste(miss, collapse = ", "), call. = FALSE)
    
    pair_fac <- paste0(tgt$var, "_pair")
    if (!pair_fac %in% names(pair.meta))
      stop("Internal error: '", pair_fac, "' not found in pair.meta.", call. = FALSE)
    
    levs <- levels(pair.meta[[pair_fac]])
    
    # Build canonical labels for ref-ref, alt-alt, and mixed(ref,alt)
    lab_refref <- paste0(tgt$ref, "_and_", tgt$ref)
    lab_altalt <- paste0(tgt$alt, "_and_", tgt$alt)
    lab_mixed  <- paste0(sort(c(tgt$ref, tgt$alt)), collapse = "_and_")
    
    # Set weights for each of these labels depending on dist.type
    w_refref <- 0
    w_altalt <- 0
    w_mixed  <- 0
    
    if (dist.type == "shift") {
      # mixed − 0.5(ref_ref + alt_alt)
      w_refref <- -0.5
      w_altalt <- -0.5
      w_mixed  <-  1
    } else if (dist.type == "total") {
      # mixed − ref_ref
      w_refref <- -1
      w_altalt <-  0
      w_mixed  <-  1
    } else if (dist.type == "var") {
      # alt_alt − ref_ref
      w_refref <- -1
      w_altalt <-  1
      w_mixed  <-  0
    }
    
    # Start with a fully named 0-vector over *all* levels (so Z can soak up "other")
    coef_names <- paste0(pair_fac, levs)
    explicitCoef <- setNames(numeric(length(coef_names)), coef_names)
    
    # Assign the weights we just computed
    # (only if those levels actually exist in this dataset)
    if (lab_refref %in% levs)
      explicitCoef[paste0(pair_fac, lab_refref)] <- w_refref
    if (lab_altalt %in% levs)
      explicitCoef[paste0(pair_fac, lab_altalt)] <- w_altalt
    if (lab_mixed %in% levs)
      explicitCoef[paste0(pair_fac, lab_mixed)]  <- w_mixed
  }
  
  ## 8) Pick the pair formula 
  pairFormulaUsed <- pairFormula %||% buildDefaultPairFormula(pair.meta)
  
  ## 8b) If the pair formula has an intercept and explicitCoef is expressed
  ##     as cell-level weights for the focus factor, rewrite it into
  ##     coefficient-space weights for (Intercept) + treatment-coded columns.
  if (!is.null(explicitCoef) &&
      is.numeric(explicitCoef) &&
      !is.null(names(explicitCoef)) &&
      !is.null(focusVar)) {
    
    pair_fac <- paste0(focusVar, "_pair")
    
    # only relevant if this factor exists in pair.meta
    if (pair_fac %in% colnames(pair.meta)) {
      
      trm_pair      <- stats::terms(pairFormulaUsed, data = pair.meta)
      has_intercept <- isTRUE(attr(trm_pair, "intercept") == 1L)
      
      if (has_intercept && !("(Intercept)" %in% names(explicitCoef))) {
        
        levs <- levels(pair.meta[[pair_fac]])
        if (!length(levs))
          stop("Focus pair factor '", pair_fac, "' has no levels.")
        
        # Expected full set of cell-level names for this factor
        fact_names_full <- paste0(pair_fac, levs)
        
        # Only transform if we truly have weights for *all* cells
        if (all(fact_names_full %in% names(explicitCoef))) {
          
          # Split factor (cell) part vs 'other' part (numeric mirrors, etc.)
          w_factor <- explicitCoef[fact_names_full]
          w_other  <- explicitCoef[setdiff(names(explicitCoef), fact_names_full)]
          
          # Sum of cell weights becomes the coefficient on the intercept
          w_sum <- sum(w_factor)
          
          coef_new <- w_other
          coef_new["(Intercept)"] <- w_sum
          
          if (length(levs) > 1L) {
            # Non-reference levels: keep their original cell weights
            nonref_levs  <- levs[-1L]                          # baseline = levs[1L]
            nonref_names <- paste0(pair_fac, nonref_levs)
            coef_new[nonref_names] <- w_factor[nonref_names]
          }
          
          explicitCoef <- coef_new
        }
      }
    }
  }
  
  
  ##  decide which pair-level block variables to use
  blockVarsPair <- character(0)
  if (isTRUE(buildBlocks)) {
    if (!is.null(pairBlockVars) && length(pairBlockVars)) {
      # allow both sample-level names ("Batch") and pair-level names ("Batch_pair")
      blockVarsPair <- normalizePairBlockVars(pairBlockVars, sample.meta, pair.meta)
    } else {
      # default: reuse the nuisance factor set that produced sample-level blocks
      blockVarsPair <- defaultPairBlockVarsFromSampleDesign(sampleDesign, sample.meta, pair.meta)
    }
  }
  
  ## 9) Build the actual model 
  out <- buildDesignMatrices(
    data        = pair.meta,
    contrast    = explicitCoef,
    formula     = pairFormulaUsed,
    na.action   = na.action,
    numericRef  = "auto",
    tol         = 1e-12,
    tolRow      = sqrt(.Machine$double.eps),
    validate    = TRUE,
    verbosity   = switch(verbosity, none="none", warn="warn", info="info"),
    computeQrZ  = TRUE,
    # CHANGED: pass blocks through to the pair level
    blockVars   = if (length(blockVarsPair)) blockVarsPair else NULL,
    buildBlocks = buildBlocks
  )
  
  ## 10) Attach extras
  out$pairs                 <- pairs
  out$pair.meta             <- pair.meta
  out$pair_formula_used     <- pairFormulaUsed
  out$focus_var             <- focusVar
  out$dist_type             <- dist.type
  out$pair_block_vars_used  <- blockVarsPair  # NEW: for transparency
  
  if (verbosity %in% c("info")) {
    fac_cols <- names(pair.meta)[grepl("_pair$", names(pair.meta))]
    num_cols <- names(pair.meta)[grepl("^pair_[A-Za-z0-9_]+_(mean|diff)$", names(pair.meta))]
    msg <- sprintf("Paired design: factor(s)=%s; numeric=%s; formula=%s",
                   if (length(fac_cols)) paste(sub("_pair$","",fac_cols), collapse=", ") else "(none)",
                   if (length(num_cols)) "(present)" else "(none)",
                   paste(deparse(pairFormulaUsed), collapse=""))
    if (length(blockVarsPair)) {
      msg <- paste0(msg, sprintf("; blocks=%s", paste(blockVarsPair, collapse = "+")))
    }
    message(msg)
  }
  
  out
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

`%||%` <- function(a, b) if (is.null(a)) b else a

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
  mf  <- stats::model.frame(rhs, data, na.action = na.action)
  
  # 1. Apply requested baselines (relevel)
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
  
  # 2. Force any remaining character columns to factors NOW.
  #    This makes the representation consistent with how model.matrix
  #    will treat them anyway, and guarantees that:
  #    - attr(F, "xlevels") has entries for them
  #    - .oneRowFromFormula() will treat them as factors
  #    - contrasts.arg will not be applied to a non-factor down the line
  for (v in names(mf)) {
    if (is.character(mf[[v]])) {
      mf[[v]] <- droplevels(factor(mf[[v]]))
    }
  }
  
  trm <- stats::terms(rhs, data = mf)   # NOTE: use mf, not the original data,
  # so terms() sees factors with baselines
  F   <- stats::model.matrix(trm, mf, contrasts.arg = contrasts.arg)
  
  # Attach metadata so downstream builders (.oneRowFromFormula, etc.) can reconstruct rows
  attr(F, "terms")     <- trm
  attr(F, "xlevels")   <- lapply(mf, function(x) if (is.factor(x)) levels(x) else NULL)
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
    coreVars <- switch(contrastSpec$type,
                       simple   = if (grepl(":", contrastSpec$term, fixed=TRUE)) parseTermVars(contrastSpec$term) else contrastSpec$term,
                       marginal = unique(c(contrastSpec$term, contrastSpec$over)),
                       lincomb  = parseTermVars(contrastSpec$term),
                       coef     = character(0))
  }
  filterNuisance(data, setdiff(facVars, unique(coreVars)))
}

# ---- Pair helpers ----

lowerTriIndices <- function(n) {
  if (n < 2) return(matrix(integer(0), ncol = 2, dimnames = list(NULL, c("i","j"))))
  idx <- which(lower.tri(matrix(NA, n, n)), arr.ind = TRUE)
  colnames(idx) <- c("i","j")
  idx
}

makeSafeVar <- function(s) {
  s <- as.character(s)
  s <- gsub("[^A-Za-z0-9_]+", "_", s)
  s <- gsub("_+", "_", s)
  s <- sub("^_", "", s)
  if (grepl("^[0-9]", s)) s <- paste0("v_", s)
  s
}

# Return variables that appear in the (sample-level) formula, split by type
varsFromFormula <- function(formula, data, na.action = stats::na.pass) {
  trm <- stats::terms(formula, data = data)
  # Build a model.frame once to evaluate factor/character vs numeric on raw columns
  mf <- stats::model.frame(trm, data, na.action = na.action)
  vars <- setdiff(names(mf), "(Intercept)")
  isFac <- vapply(mf, function(x) is.factor(x) || is.character(x), logical(1))
  isNum <- vapply(mf, function(x) is.numeric(x) && !is.matrix(x), logical(1))
  list(factors = vars[isFac], numeric = vars[isNum])
}


# Infer pair target (var, alt, ref) from a **sample-level design** object.
# Uses sampleDesign$contrast_spec, which is created by buildDesignMatrices().
inferPairTargetFromDesign <- function(sampleDesign, data) {
  spec <- sampleDesign$contrast_spec
  if (is.null(spec)) {
    stop("sampleDesign$contrast_spec is NULL; cannot infer paired contrast. ",
         "Provide 'pairContrast' explicitly.", call. = FALSE)
  }
  # Accept simple single-factor and marginal single-factor; reject interactions/lincomb/coef.
  if (is.list(spec)) {
    if (identical(spec$type, "simple") && !grepl(":", spec$term, fixed = TRUE)) {
      return(list(var = spec$term, alt = spec$num, ref = spec$den, source = "simple"))
    }
    if (identical(spec$type, "marginal")) {
      return(list(var = spec$term, alt = spec$num, ref = spec$den, source = "marginal"))
    }
    stop("Cannot infer paired contrast from sampleDesign$contrast_spec type='", spec$type,
         "'. Use a single-factor simple/marginal sample contrast or provide 'pairContrast'.",
         call. = FALSE)
  }
  stop("Unsupported contrast spec in sampleDesign; provide 'pairContrast' explicitly.", call. = FALSE)
}

inferPairTargetFromSampleContrast <- function(contrast, data) {
  if (is.character(contrast) && length(contrast) == 3) {
    return(list(var = contrast[[1]], alt = contrast[[2]], ref = contrast[[3]], source = "triple"))
  }
  if (is.list(contrast) && length(contrast)) {
    t <- contrast$type %||% stop("Structured contrast requires 'type'.")
    if (t == "simple" && is.character(contrast$term) && !grepl(":", contrast$term, fixed = TRUE)) {
      return(list(var = contrast$term, alt = contrast$num, ref = contrast$den, source = "simple"))
    }
    if (t == "marginal" && is.character(contrast$term)) {
      return(list(var = contrast$term, alt = contrast$num, ref = contrast$den, source = "marginal"))
    }
    stop("Cannot infer paired contrast from sample-level contrast type='", t,
         "'. Use a single-factor simple/marginal contrast or provide 'pairContrast'.")
  }
  stop("Unsupported sample-level contrast format. Use c('factor','ALT','REF'), ",
       "or list(type='simple'|'marginal', term='factor', num=, den=), ",
       "or provide 'pairContrast'.")
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

# Parse a pair-level RHS formula and discover requested pair_<safe(var)>_(mean|diff)
# Returns the *sample-level* variable names to include by matching safe names.
extractPairNumericVarsFromFormula <- function(pairFormula, data) {
  if (is.null(pairFormula)) return(character(0))
  rhs <- paste(deparse(pairFormula), collapse = "")
  # Matches pair_<anyword>_mean or pair_<anyword>_diff
  m <- gregexpr("\\bpair_([A-Za-z0-9_]+)_(mean|diff)\\b", rhs, perl = TRUE)
  hits <- regmatches(rhs, m)[[1]]
  if (length(hits) == 0) return(character(0))
  safe <- sub("^pair_([A-Za-z0-9_]+)_(mean|diff)$", "\\1", hits, perl = TRUE)
  # Try to unsafify (best-effort): convert underscores back to names if exact match exists
  # Prefer exact sample-column matches first; otherwise keep as-is (user may have already safe names)
  candidates <- names(data)
  out <- vapply(safe, function(s) {
    # If exact name exists in data, keep it; else try exact as-is; else return s as a fallback
    if (s %in% candidates) s else s
  }, character(1))
  unique(out)
}


# Build pair-level metadata for a set of factor vars (+ numeric mirrors)
# - factorVars: character vector of sample-level factor/character variables to pairify
# - focusVar: optional factor used later to seed default "shift/total/var" contrast
# - numericVars: character vector of sample-level numerics to mirror into pair_<var>_{mean,diff}
pairifyMeta <- function(meta,
                        pairsIdx,
                        factorVars = character(0),
                        focusVar = NULL,
                        numericVars = character(0)) {
  stopifnot(is.data.frame(meta))
  n <- nrow(meta); if (!n) stop("Empty metadata.")
  if (is.null(pairsIdx) || ncol(pairsIdx) != 2L)
    stop("pairsIdx must be 2-column (i,j).")
  i <- as.integer(pairsIdx[,1]); j <- as.integer(pairsIdx[,2])
  if (any(i < 1 | i > n | j < 1 | j > n))
    stop("pairsIdx out of bounds.")
  n_pairs <- length(i)
  
  cols <- list()
  friendly_map <- character(0)
  
  # helper: stable unordered label "A_and_B" with A,B ordered by factor level index
  join_levels <- function(levs, a, b) {
    ia <- match(a, levs); ib <- match(b, levs)
    if (is.na(ia) || is.na(ib))
      stop("Unknown level when forming pair labels.")
    if (ia <= ib) paste0(a, "_and_", b) else paste0(b, "_and_", a)
  }
  
  # ---- pairify categorical variables
  if (length(factorVars)) {
    for (v in factorVars) {
      if (!v %in% names(meta))
        stop("Factor '", v, "' not found in metadata.")
      raw <- meta[[v]]
      if (!(is.factor(raw) || is.character(raw)))
        stop("Variable '", v, "' must be factor/character to be pairified as factor.")
      
      f <- factor(raw)           # drop unused levels at sample level
      L <- levels(f)
      
      fi <- as.character(f[i])
      fj <- as.character(f[j])
      lev_lbl <- mapply(function(a,b) join_levels(L, a, b),
                        fi, fj,
                        SIMPLIFY = TRUE, USE.NAMES = FALSE)
      
      # build all possible unordered-with-replacement combos of levels
      all_pairs <- unique(c(outer(L, L, Vectorize(function(a,b) join_levels(L,a,b)))))
      
      fac_name <- paste0(v, "_pair")  # <-- no trailing underscore now
      cols[[fac_name]] <- factor(lev_lbl, levels = all_pairs)
    }
  }
  
  # ---- mirror numeric variables (mean + diff)
  if (length(numericVars)) {
    for (v in numericVars) {
      if (!v %in% names(meta))
        stop("Numeric variable '", v, "' not found in metadata.")
      x <- meta[[v]]
      if (!(is.numeric(x) && !is.matrix(x)))
        stop("Variable '", v, "' must be numeric for pair mirrors.")
      
      xi <- as.numeric(x[i])
      xj <- as.numeric(x[j])
      
      base_safe <- paste0("pair_", makeSafeVar(v))
      nm_mean   <- paste0(base_safe, "_mean")
      nm_diff   <- paste0(base_safe, "_diff")
      
      cols[[nm_mean]] <- 0.5 * (xi + xj)
      cols[[nm_diff]] <- (xi - xj)  # oriented i - j
      
      # record friendly aliases for explicit contrasts
      friendly_map[nm_mean] <- paste0("pair[", v, "]_mean")
      friendly_map[nm_diff] <- paste0("pair[", v, "]_diff")
    }
  }
  
  out <- if (length(cols)) {
    as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    data.frame(row.names = seq_len(n_pairs))[ , 0]
  }
  
  rownames(out) <- paste0("pair_", i, "_", j)
  attr(out, "friendly_map") <- friendly_map
  out
}


# Extract the focus variable from coefficient names in an explicit pairContrast.
# Works with colnames like "Group_pair_Group1_and_Group2".
extractFocusVarFromNamedContrast <- function(pc_names) {
  stopifnot(length(pc_names) > 0)
  m <- regexec("^([[:alnum:]_.]+)_pair_", pc_names)
  vars <- unique(vapply(
    regmatches(pc_names, m),
    function(mm) if (length(mm) > 1) mm[[2]] else NA_character_,
    character(1)
  ))
  vars <- vars[!is.na(vars)]
  if (!length(vars)) return(NULL)
  if (length(vars) > 1) {
    stop("Explicit pairContrast references multiple pair factors: ",
         paste(vars, collapse = ", "),
         ". Limit to a single factor or use list(var=,cells=,coefs=).")
  }
  vars[[1]]
}

# Map user-friendly coefficient names to real pair.meta column names.
# - Accepts:
#   * factor cell columns like "Group_pairGroup1_and_Group2"
#   * numeric columns like "pairdage_diff"
#   * friendly aliases like "pair[age]_diff"
normalizeExplicitPairCoefNames <- function(coefvec, pair_meta, focusVar = NULL) {
  stopifnot(is.numeric(coefvec), !is.null(names(coefvec)))
  nm_in <- names(coefvec)
  
  friendly_map <- attr(pair_meta, "friendly_map")
  if (is.null(friendly_map)) friendly_map <- character(0)
  
  # all columns that actually exist in pair.meta
  safe_numeric <- colnames(pair_meta)
  
  # expected factor-expanded coefficient names
  # For each pair factor column like 'Group_pair', model.matrix(~0+Group_pair)
  # will generate columns like 'Group_pairGroup1_and_Group2'.
  # We don't know the full expansion unless we look at levels().
  fact_expect <- character(0)
  fact_cols <- colnames(pair_meta)[grepl("_pair$", colnames(pair_meta))]
  if (length(fact_cols)) {
    fact_expect <- unlist(lapply(fact_cols, function(fc) {
      paste0(fc, levels(pair_meta[[fc]]))
    }))
  }
  
  out_names <- character(length(nm_in))
  for (k in seq_along(nm_in)) {
    nm <- nm_in[k]
    
    # case 1: already matches an expanded factor coef or a numeric pair term in pair.meta
    if (nm %in% fact_expect || nm %in% safe_numeric) {
      out_names[k] <- nm
      next
    }
    
    # case 2: friendly numeric alias like "pair[age]_diff"
    hit <- names(friendly_map)[friendly_map == nm]
    if (length(hit) == 1L) {
      out_names[k] <- hit
      next
    }
    
    # case 3: already safe numeric-style name
    if (grepl("^pair_[A-Za-z0-9_]+_(mean|diff)$", nm)) {
      out_names[k] <- nm
      next
    }
    
    stop("Unknown coefficient in explicit pairContrast: '", nm,
         "'. Available examples include factor columns like {",
         paste(head(fact_expect, 6), collapse=", "),
         "} or numeric mirrors like {",
         paste(head(names(friendly_map), 6), collapse=", "),
         "}.")
  }
  
  stats::setNames(unclass(coefvec), out_names)
}


# Build a coefficient vector from an explicit "cells" spec:
# cells is like c("Group1_and_Group1"= -0.5, "Group1_and_Group2"=1, ...)
# Returns a named numeric over the *design columns* for that pair factor.
buildCoefFromCells <- function(pair_meta, var, cells) {
  # var is e.g. "Group"
  pair_fac <- paste0(var, "_pair")
  if (!pair_fac %in% names(pair_meta))
    stop("Pair metadata lacks '", pair_fac, "'.")
  
  levs <- levels(pair_meta[[pair_fac]])
  
  # Validate provided cells
  bad <- setdiff(names(cells), levs)
  if (length(bad)) {
    stop("Unknown cells in explicit pair contrast: ",
         paste(bad, collapse = ", "),
         ". Available cells: ", paste(levs, collapse = ", "))
  }
  
  coef_names <- paste0(pair_fac, levs)   # e.g. "Group_pairGroup1_and_Group1"
  res <- setNames(numeric(length(coef_names)), coef_names)
  
  for (nm in names(cells)) {
    target_col <- paste0(pair_fac, nm)
    res[target_col] <- res[target_col] + as.numeric(cells[[nm]])
  }
  res
}


# Map user-supplied block vars to *pair-level* columns
# Accept both sample-level names ("Batch") and pair-level names ("Batch_pair").
normalizePairBlockVars <- function(vars, sample.meta, pair.meta) {
  if (!length(vars)) return(character(0))
  out <- character(0)
  for (v in vars) {
    if (v %in% names(pair.meta)) {
      out <- c(out, v)
    } else if (v %in% names(sample.meta) && paste0(v, "_pair") %in% names(pair.meta)) {
      out <- c(out, paste0(v, "_pair"))
    } else {
      stop("pairBlockVars entry '", v,
           "' is neither a pair-level column in pair.meta nor a sample-level factor present in metadata.")
    }
  }
  # must be factors on the pair meta
  bad <- out[!vapply(pair.meta[out], function(x) is.factor(x) || is.character(x), logical(1))]
  if (length(bad)) {
    stop("All pairBlockVars must resolve to factor pair columns. Offending: ",
         paste(bad, collapse=", "))
  }
  unique(out)
}

# Default: reuse sample-level nuisance factors that defined sampleDesign$blocks,
# then pairify them (e.g., 'Batch' -> 'Batch_pair') if present in pair.meta.
defaultPairBlockVarsFromSampleDesign <- function(sampleDesign, sample.meta, pair.meta) {
  # Derive the same nuisance set used at the sample level
  nuis <- deriveNuisanceFactors(sampleDesign$formula_used,
                                data = sample.meta,
                                contrastSpec = sampleDesign$contrast_spec,
                                explicit = NULL)
  if (!length(nuis)) return(character(0))
  cand <- paste0(nuis, "_pair")
  cand[cand %in% names(pair.meta)]
}


# Build a default pair-level formula if the user didn't supply one.
buildDefaultPairFormula <- function(pair.meta) {
  stopifnot(is.data.frame(pair.meta))
  cols <- names(pair.meta)
  
  # factor columns created by pairifyMeta()
  pair_factors <- cols[grepl("_pair$", cols)]
  # numeric mirror columns (pair_<num>_(mean|diff))
  pair_numeric <- cols[grepl("^pair_[A-Za-z0-9_]+_(mean|diff)$", cols)]
  
  rhs_terms <- c(pair_factors, pair_numeric)
  if (!length(rhs_terms)) return(~ 1)
  
  has_fac <- length(pair_factors) > 0
  rhs <- paste(rhs_terms, collapse = " + ")
  as.formula(paste("~", if (has_fac) paste("0 +", rhs) else rhs))
}

