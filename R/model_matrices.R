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
  } else { # check for mixed effect models
    formula <- checkFormula(formula)
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
#'   \item `"cross"` — pairs drawn from **different** levels of `var`;
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
#'           \item{`"shift"`}{`cross − ½(REF:REF + ALT:ALT)` (discordant vs mean of concordants)}
#'           \item{`"total"`}{`cross − REF:REF` (discordant vs reference-concordant)}
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
                                    na.action = stats::na.pass,
                                    verbosity = c("none","warn","info")) {
  stopifnot(is.data.frame(sample.meta))
  if (missing(sampleDesign) || is.null(sampleDesign) || is.null(sampleDesign$formula_used)) {
    stop("sampleDesign with a valid $formula_used is required.", call. = FALSE)
  }
  dist.type <- match.arg(dist.type)
  verbosity <- match.arg(verbosity)
  
  # 1) All unordered pairs
  n <- nrow(sample.meta)
  pairs <- lowerTriIndices(n)
  if (!nrow(pairs)) stop("Not enough samples to form pairs (need at least 2).")
  
  # 2) Determine focus factor / explicit contrast
  focusVar      <- NULL
  explicitCoef  <- NULL
  explicitCells <- NULL
  
  if (!is.null(pairContrast)) {
    if (is.numeric(pairContrast) && !is.null(names(pairContrast))) {
      focusVar <- tryCatch(extractFocusVarFromNamedContrast(names(pairContrast)),
                           error = function(e) stop(conditionMessage(e), call. = FALSE))
    } else if (is.list(pairContrast)) {
      if (!is.null(pairContrast$cells)) {
        if (is.null(pairContrast$var))
          stop("Explicit list pairContrast requires 'var' when 'cells' are provided.", call. = FALSE)
        if (!is.numeric(pairContrast$cells) || is.null(names(pairContrast$cells)))
          stop("'cells' must be a named numeric vector over {'cross','<L>:<L>'}.", call. = FALSE)
        focusVar      <- as.character(pairContrast$var)
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
  
  # If still no focus, try to infer from sampleDesign contrast_spec
  if (is.null(focusVar) && is.null(explicitCells) && is.null(explicitCoef)) {
    tgt <- tryCatch(inferPairTargetFromDesign(sampleDesign, sample.meta),
                    error = function(e) e)
    if (inherits(tgt, "error")) {
      stop(paste0("Cannot infer paired contrast from sampleDesign$contrast_spec: ",
                  conditionMessage(tgt), "\n",
                  "Either provide 'pairContrast', or use a single-factor simple/marginal sample contrast at the sample level."),
           call. = FALSE)
    }
    focusVar <- tgt$var
  }
  
  # 3) Mirror numeric variables from the sample-level formula (+ any required by pairFormula)
  vf <- varsFromFormula(sampleDesign$formula_used, data = sample.meta, na.action = na.action)
  useNumeric <- vf$numeric
  # Also auto-include sample numerics referenced by pairFormula via pair_<var>_(mean|diff)
  extraNums <- extractPairNumericVarsFromFormula(pairFormula, sample.meta)
  if (length(extraNums)) useNumeric <- unique(c(useNumeric, extraNums))
  
  # (optional) If pairFormula references a '<var>_pair' factor and we still lack focusVar, set it.
  if (is.null(focusVar) && !is.null(pairFormula)) {
    rhs <- paste(deparse(pairFormula), collapse = "")
    m <- regmatches(rhs, regexpr("\\b([A-Za-z0-9_.]+)_pair\\b", rhs, perl = TRUE))
    if (length(m) && nzchar(m)) {
      v <- sub("_pair$", "", m)
      if (v %in% names(sample.meta) && (is.factor(sample.meta[[v]]) || is.character(sample.meta[[v]]))) {
        focusVar <- v
      }
    }
  }
  
  # 4) Pairify metadata
  pair.meta <- pairifyMeta(sample.meta, pairsIdx = pairs, focusVar = focusVar, numericVars = useNumeric)
  
  # 5) Build explicit coefficient vector if provided as cells/coefs
  if (!is.null(explicitCells)) {
    coef_cells <- buildCoefFromCells(pair.meta, var = focusVar, cells = explicitCells$cells)
    explicitCoef <- coef_cells
    if (!is.null(explicitCells$coefs)) {
      explicitCoef <- c(explicitCoef, normalizeExplicitPairCoefNames(explicitCells$coefs, pair.meta, focusVar))
    }
  }
  
  # 6) If user passed a named-numeric vector directly, normalize names now
  if (!is.null(pairContrast) && is.numeric(pairContrast) && !is.null(names(pairContrast))) {
    explicitCoef <- normalizeExplicitPairCoefNames(pairContrast, pair.meta, focusVar)
  }
  
  # 7) If still no explicit vector, synthesize defaults from sampleDesign contrast_spec
  if (is.null(explicitCoef)) {
    tgt <- inferPairTargetFromDesign(sampleDesign, sample.meta)  # validated path
    if (!tgt$var %in% names(sample.meta))
      stop("Target factor '", tgt$var, "' not found in metadata.", call. = FALSE)
    f.levels <- levels(factor(sample.meta[[tgt$var]]))
    miss <- setdiff(c(tgt$ref, tgt$alt), f.levels)
    if (length(miss))
      stop("ALT/REF levels not found for factor '", tgt$var, "': ",
           paste(miss, collapse = ", "), call. = FALSE)
    
    pair.var    <- paste0(tgt$var, "_pair")
    if (!pair.var %in% names(pair.meta))
      stop("Internal error: '", pair.var, "' not found in pair.meta.", call. = FALSE)
    pair.levels <- levels(pair.meta[[pair.var]])
    w_by_level  <- pairContrastWeights(ref = tgt$ref, alt = tgt$alt,
                                       dist.type = dist.type, level.order = pair.levels)
    names(w_by_level) <- paste0(pair.var, pair.levels)  # match model.matrix(~0+pair.var) names
    explicitCoef <- w_by_level
  }
  
  # 8) Build paired design with your existing constructor
  out <- buildDesignMatrices(
    data        = pair.meta,
    contrast    = explicitCoef,             # named numeric over pair colnames
    formula     = pairFormula %||% NULL,    # if NULL, default saturated over pair.meta
    na.action   = na.action,
    numericRef  = "auto",
    tol         = 1e-12,
    tolRow      = sqrt(.Machine$double.eps),
    validate    = TRUE,
    verbosity   = switch(verbosity, none="none", warn="warn", info="info"),
    computeQrZ  = TRUE,
    buildBlocks = FALSE
  )
  
  # 9) Attach pair extras to mirror your sample-level return shape
  out$pairs     <- pairs
  out$pair.meta <- pair.meta
  
  if (verbosity %in% c("info")) {
    message(sprintf("Paired design (auto): factor='%s', dist.type='%s', mirrored numeric: %s.",
                    focusVar %||% "(none)", dist.type,
                    if (length(useNumeric)) paste(useNumeric, collapse=", ") else "(none)"))
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

# Parse a pair-level RHS formula and discover requested pair_<var>_(mean|diff)
# Returns the *sample-level* variable names to include by matching safe names.
extractPairNumericVarsFromFormula <- function(pairFormula, sample.meta) {
  if (is.null(pairFormula)) return(character(0))
  rhs <- paste(deparse(pairFormula), collapse = "")
  # capture tokens like pair_age_mean or pair_BMI_diff
  m <- gregexpr("\\bpair_([A-Za-z0-9_]+)_(mean|diff)\\b", rhs, perl = TRUE)
  hits <- regmatches(rhs, m)[[1]]
  if (!length(hits)) return(character(0))
  bases <- sub("^pair_([A-Za-z0-9_]+)_(mean|diff)$", "\\1", hits)
  # map safe bases back to sample columns by makeSafeVar() equality
  sampleNums <- names(which(vapply(sample.meta, function(x) is.numeric(x) && !is.matrix(x), logical(1))))
  want <- character(0)
  for (b in unique(bases)) {
    hit <- sampleNums[ vapply(sampleNums, function(v) identical(makeSafeVar(v), b), logical(1)) ]
    want <- c(want, hit)
  }
  unique(want)
}


pairContrastWeights <- function(ref, alt, dist.type = c("shift","total","var"),
                                level.order = NULL) {
  dist.type <- match.arg(dist.type)
  refref <- paste0(ref, ":", ref)
  altalt <- paste0(alt, ":", alt)
  key    <- c(refref, "cross", altalt)
  
  w <- setNames(numeric(3), key)
  if (dist.type == "shift") { w[refref] <- -0.5; w["cross"] <-  1;  w[altalt] <- -0.5 }
  if (dist.type == "total") { w[refref] <- -1;   w["cross"] <-  1;  w[altalt] <-  0   }
  if (dist.type == "var")   { w[refref] <- -1;   w["cross"] <-  0;  w[altalt] <-  1   }
  
  if (!is.null(level.order)) {
    out <- setNames(numeric(length(level.order)), level.order)
    m   <- intersect(names(out), names(w))
    out[m] <- w[m]
    return(out)
  }
  w
}

candidatePairColumnsMsg <- function(data) {
  facVars <- names(which(vapply(data, function(x) is.factor(x) || is.character(x), logical(1))))
  numVars <- names(which(vapply(data, function(x) is.numeric(x) && !is.matrix(x), logical(1))))
  pieces <- c()
  for (v in facVars) {
    L <- levels(droplevels(factor(data[[v]])))
    if (length(L) >= 2) {
      l1 <- L[1]; l2 <- L[2]
      pieces <- c(pieces,
                  paste0("If focus='", v, "': ",
                         paste(c(paste0(v, "_paircross"),
                                 paste0(v, "_pair", l1, ":", l1),
                                 paste0(v, "_pair", l2, ":", l2)), collapse = ", ")))
    }
  }
  if (length(numVars)) {
    pieces <- c(pieces, paste0("Numeric pair terms (friendly aliases accepted): ",
                               paste(c(paste0("pair[", numVars, "]_mean"),
                                       paste0("pair[", numVars, "]_diff")), collapse = ", ")))
  }
  if (!length(pieces)) "(no obvious pair columns)" else paste(pieces, collapse = " | ")
}

# Pairify sample metadata into one row per unordered pair (i<j).
# numericVars: character vector of sample-level numeric variables to include.
pairifyMeta <- function(meta, pairsIdx, focusVar = NULL, numericVars = NULL) {
  stopifnot(is.data.frame(meta))
  n <- nrow(meta); if (!n) stop("Empty metadata.")
  if (is.null(pairsIdx) || ncol(pairsIdx) != 2L) stop("pairsIdx must be 2-column (i,j).")
  i <- as.integer(pairsIdx[,1]); j <- as.integer(pairsIdx[,2])
  n_pairs <- length(i)
  
  cols <- list()
  friendly_map <- character(0)
  
  # Optional composition factor for the focus variable
  if (!is.null(focusVar)) {
    if (!focusVar %in% names(meta))
      stop("focusVar '", focusVar, "' not found in metadata.")
    fv <- meta[[focusVar]]
    if (!(is.factor(fv) || is.character(fv)))
      stop("focusVar '", focusVar, "' must be factor/character.")
    f  <- factor(fv)  # drop unused at sample level
    L  <- levels(f)
    fi <- as.character(f[i]); fj <- as.character(f[j])
    comp <- ifelse(fi == fj, paste0(fi, ":", fj), "cross")
    comp_levels <- unique(c("cross", paste0(L, ":", L)))
    cols[[paste0(focusVar, "_pair")]] <- factor(comp, levels = comp_levels)
  }
  
  # Numeric pair summaries restricted to numericVars
  if (length(numericVars)) {
    for (v in numericVars) {
      xi <- as.numeric(meta[[v]][i]); xj <- as.numeric(meta[[v]][j])
      base_safe <- paste0("pair_", makeSafeVar(v))
      nm_mean   <- paste0(base_safe, "_mean")
      nm_diff   <- paste0(base_safe, "_diff")
      cols[[nm_mean]] <- 0.5 * (xi + xj)
      cols[[nm_diff]] <- (xi - xj)           # order: i - j
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


extractFocusVarFromNamedContrast <- function(pc_names) {
  m <- regmatches(pc_names, regexec("^([[:alnum:]_.]+)_pair", pc_names))
  vars <- unique(vapply(m, function(mm) if (length(mm)) mm[[2]] else NA_character_, character(1)))
  vars <- vars[!is.na(vars)]
  if (!length(vars)) return(NULL)
  if (length(vars) > 1)
    stop("Explicit pairContrast references multiple pair factors: ",
         paste(vars, collapse = ", "),
         ". Limit to a single factor or use list(var=,cells=,coefs=).")
  vars[[1]]
}

normalizeExplicitPairCoefNames <- function(coefvec, pair_meta, focusVar = NULL) {
  stopifnot(is.numeric(coefvec), !is.null(names(coefvec)))
  names_in <- names(coefvec)
  
  friendly_map <- attr(pair_meta, "friendly_map"); if (is.null(friendly_map)) friendly_map <- character(0)
  safe_numeric <- colnames(pair_meta)
  
  fact_expect <- character(0)
  if (!is.null(focusVar) && paste0(focusVar, "_pair") %in% names(pair_meta)) {
    levs <- levels(pair_meta[[paste0(focusVar, "_pair")]])
    fact_expect <- paste0(paste0(focusVar, "_pair"), levs)
  }
  
  safe_out <- character(length(names_in))
  for (k in seq_along(names_in)) {
    nm <- names_in[k]
    if (nm %in% safe_numeric || nm %in% fact_expect) { safe_out[k] <- nm; next }
    hit <- names(friendly_map)[friendly_map == nm]
    if (length(hit) == 1L) { safe_out[k] <- hit; next }
    if (grepl("^pair_[A-Za-z0-9_]+_(mean|diff)$", nm)) { safe_out[k] <- nm; next }
    stop("Unknown coefficient in explicit pairContrast: '", nm,
         "'. Available examples: ",
         paste(c(head(fact_expect, 6), head(names(friendly_map), 6)), collapse = ", "),
         if (length(fact_expect) > 6 || length(friendly_map) > 6) ", …")
  }
  stats::setNames(unclass(coefvec), safe_out)
}

buildCoefFromCells <- function(pair_meta, var, cells) {
  pm <- as.data.frame(pair_meta)
  pair_var <- paste0(var, "_pair")
  if (!pair_var %in% names(pm)) stop("Pair metadata lacks '", pair_var, "'.")
  levs <- levels(pm[[pair_var]])
  valid_cells <- unique(c("cross", levs))
  if (!all(names(cells) %in% valid_cells)) {
    stop("Unknown cells in explicit pair contrast: ",
         paste(setdiff(names(cells), valid_cells), collapse = ", "),
         ". Available cells: ", paste(valid_cells, collapse = ", "))
  }
  coef_names <- paste0(pair_var, levs)
  res <- setNames(numeric(length(coef_names)), coef_names)
  for (nm in names(cells)) {
    target <- if (nm == "cross") paste0(pair_var, "cross") else paste0(pair_var, nm)
    res[target] <- res[target] + as.numeric(cells[[nm]])
  }
  res
}