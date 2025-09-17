#' @keywords internal
estimateGraphVarianceSignificance <- function(adj.mat, signal, n.permutations=5000) {
  if (!is.numeric(signal)) {
    signal <- as.factor(signal)
    signal[is.na(signal)] <- table(signal) %>% which.max() %>% names()
    comp.op <- "!="
  } else {
    signal[is.na(signal)] <- median(signal, na.rm=TRUE)
    comp.op <- "-"
  }

  obs.var <- signal %>% outer(., ., comp.op) %>% {. * . * adj.mat} %>% sum()
  perm.vars <- sapply(1:n.permutations, function(i) {
    sample(signal) %>% outer(., ., comp.op) %>% {. * . * adj.mat} %>% sum()
  })

  pvalue <- (sum(perm.vars <= obs.var) + 1) / (n.permutations + 1)
  pr2 <- 1 - obs.var / median(perm.vars)
  return(list(pvalue=pvalue, pr2=pr2))
}

#' @keywords internal
adjacencyMatrixFromPaiwiseDists <- function(p.dists, trim=0.05, k=NULL) {
  adj.mat <- p.dists %>% pmin(quantile(., 1 - trim)) %>%
    {pmax(0, . - quantile(., trim))} %>% {1 - . / max(.)} %>%
    matrix(ncol=ncol(p.dists))
  diag(adj.mat) <- 0

  if (!is.null(k) && (k < ncol(adj.mat))) { # remove edges for all but k nearest neighbors
    adj.mat %<>% apply(1, function(r) ifelse(r < sort(r, decreasing=TRUE)[k], 0, r)) %>% {(. + t(.)) / 2}
  }

  dimnames(adj.mat) <- dimnames(p.dists)

  return(adj.mat)
}

#' @keywords internal
estimateUMAPOnDistances <- function(p.dists, n.neighbors=15, verbose=FALSE, ...) {
  set.seed(42)
  idx <- (1:ncol(p.dists)) %>% lapply(function(i) head(order(p.dists[i,]), n.neighbors))
  dists <- (1:ncol(p.dists)) %>% lapply(function(i) p.dists[i, idx[[i]]])

  idx <- do.call(rbind, idx)
  dists <- do.call(rbind, dists)

  umap <- uwot::umap(data.frame(x=rep(0, nrow(dists))), nn_method=list(idx=idx, dist=dists), verbose=verbose, ...)

  return(umap)
}

#' @keywords internal
validateDesign <- function(formula, sample.meta = NULL, contrast = NULL, verbose = FALSE) {
  if (is.null(formula)) stop("Design formula must be provided.")

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
    formula.str <- paste("~", paste(rhs.terms, collapse = " + "))
  }
  parsedTerms <- terms(stats::as.formula(formula.str))
  termLabels  <- attr(parsedTerms, "term.labels")

  # default contrast if none specified (uses first term; last vs first level)
  if (is.null(contrast)) {
    if (length(termLabels) == 0) stop("No terms found in the design formula to use as contrast.")
    if (is.null(sample.meta)) stop("sample.meta must be given to infer a default contrast.")
    contrast.var <- termLabels[1]
    if (!contrast.var %in% colnames(sample.meta)) {
      stop(sprintf("Contrast variable '%s' not found in sample metadata.", contrast.var))
    }
    levels.contrast <- levels(factor(sample.meta[[contrast.var]]))
    if (length(levels.contrast) < 2) {
      stop(sprintf("Contrast variable '%s' must have at least two levels.", contrast.var))
    }
    contrast <- c(contrast.var, levels.contrast[length(levels.contrast)], levels.contrast[1]) # var,alt,ref
  }
  if (length(contrast) != 3) {
    stop("Contrast must be a vector of length 3: c('variable', 'level1', 'level2').")
  }
  if (verbose) {
    if (containsRandomEffects) message("Random effect terms provided will be treated as fixed effect terms.")
    message(sprintf("Final design formula: %s", formula.str))
  }
  list("formula" = parsedTerms, "contrast" = contrast)
}

#' @keywords internal
subsetMetadata <- function(sample.meta,
                           formula,                      # validateDesign() output OR formula OR character
                           drop.unused.levels = TRUE,
                           extra = NULL,
                           include.contrast = TRUE) {

  stopifnot(is.data.frame(sample.meta))

  # Extract a formula and (optionally) contrast from `design`
  if (is.list(formula) && !is.null(formula$formula)) {
    # validateDesign() returns a terms object in $formula
    terms.obj <- formula$formula
    fml <- formula(terms.obj)
    contrast  <- formula$contrast
  } else if (inherits(formula, "formula")) {
    fml <- formula
    contrast <- NULL
  } else if (is.character(formula)) {
    fml <- stats::as.formula(formula)
    contrast <- NULL
  } else {
    stop("`design` must be: the list returned by validateDesign(), or a formula/character formula.")
  }

  # Expand the formula using the data so '.' is resolved and term structure is known
  # The factors matrix rows = variables used; columns = term labels
  tt <- terms(fml, data = sample.meta)

  fac <- attr(tt, "factors")
  rhs.vars <- if (is.null(fac)) character(0) else rownames(fac)

  # Optionally ensure contrast variable is included
  if (include.contrast && !is.null(contrast) && length(contrast) >= 1) {
    rhs.vars <- union(rhs.vars, contrast[1])
  }

  # Add any extras
  if (!is.null(extra) && length(extra)) rhs.vars <- union(rhs.vars, extra)

  # Sanity checks and subsetting
  miss <- setdiff(rhs.vars, names(sample.meta))
  if (length(miss)) stop("Missing columns in 'meta': ", paste(miss, collapse = ", "))

  out <- sample.meta[, rhs.vars, drop = FALSE]
  rownames(out) <- rownames(sample.meta)

  if (drop.unused.levels) {
    for (nm in names(out)) {
      if (is.factor(out[[nm]])) out[[nm]] <- droplevels(out[[nm]])
    }
  }
  out
}

#' @keywords internal
getSampleGroups <- function(sample.meta, contrast, sample.id = NULL) {
  if (is.null(contrast)) return(NULL)
  var <- contrast[1]; alt <- contrast[2]; ref <- contrast[3]
  if (!var %in% colnames(sample.meta)) {
    stop(sprintf("Contrast variable '%s' not found in sample metadata.", var))
  }
  rn <- rownames(sample.meta)
  if ((is.null(rn) || anyNA(rn)) && !is.null(sample.id)) {
    if (!sample.id %in% colnames(sample.meta)) {
      stop(sprintf("`sample.id` column '%s' not found in sample metadata.", sample.id))
    }
    rn <- as.character(sample.meta[[sample.id]])
  }
  if (is.null(rn) || anyNA(rn)) {
    stop("Sample identifiers are unavailable: set rownames(sample_meta) OR provide a valid `sample.id` column.")
  }

  vals <- as.character(sample.meta[[var]])
  names(vals) <- rn
  vals <- vals[vals %in% c(ref, alt)]
  factor(vals, levels = c(ref, alt))
}


#' Prepare full/core/nuisance model matrices, a contrast vector, core row mask, and blocks
#' 
#' Defaults:
#' - Core X: only the factor levels that appear in the contrast (subset cell-means);
#'           numeric core variables get one column.
#' - Nuisance Z: all other columns in `meta`, with an intercept (good for FL).
#' - Blocks: interaction of nuisance *factor* variables (continuous ignored).
#' 
#' # Ensures the planned contrast sits entirely in X (zero on Z / intercept).
# - If the contrasted factor has ≥1 non-participating level present:
#     X: one column per participating level; Z: complement dummies (minus one baseline) + intercept.
# - If only the contrasted levels exist:
#     X: single contrast regressor g = sum_j w_j * 1_{level j}  (weights must sum to zero); Z as-is.
#' 
#' @param meta data.frame of per-sample covariates (factors or numeric)
#' @param contrast list(var, weights) OR named numeric like c("var=Level"=1)
#' @param core character vector of core variable names (default: contrasted var only)
#' @param nuisance character vector of nuisance variable names (default: all others)
#' @param block.vars optional character vector to override block definition
#' @param na.action NA handler passed to model.matrix (default stats::na.pass)
#' @return list(F, X, Z, contrast.coef, core.rows, blocks)
#'   F, X, Z - full, core, and nuisance model matrices
#'   qrZ - QR of Z (for residualization based on Z)
#'   contrast.F - linear combination implementing desired contrast on the
#'      current coefficient order relative to F. i.e. to get contrast value:
#'      crossprod(coef(fit), contrast.coef)
#'   contrast.X - same thing, but for X
#'   core.rows - logical vector specifying a subset of rows which are needed
#'      to fit the core coefficients (for optimizing FL)
#'   blocks - factor on rows specifying blocks in which randomization should be
#'      performed (if not using FL)
# helper used here; keep if not already in scope
# helper for "var=level" names
#' @keywords internal
buildDesignMatrices <- function(sample.meta, # data.frame of covariates (subsetted to the terms in the design)
               contrast,              # list(var,weights) or named "var=level"
               nuisance   = NULL,     # default: all other columns
               core.extra = NULL,     # optional other core variables
               block.vars = NULL,
               na.action  = stats::na.pass) {
  stopifnot(is.data.frame(sample.meta))

  ctr <- parseContrast(contrast, sample.meta)
  if (is.null(nuisance)) nuisance <- setdiff(names(sample.meta), c(ctr$var, core.extra))
  if (is.null(core.extra)) core.extra <- character(0)

  rn <- rownames(sample.meta); if (is.null(rn)) rn <- as.character(seq_len(nrow(sample.meta)))

  # keep only nuisance vars that actually vary (avoids 1-level factor error)
  nuisance.eff <- filterNuisance(sample.meta, nuisance)

  # ---- Z (with intercept iff any effective nuisance) ----
  Z <- NULL
  if (length(nuisance.eff)) {
  form.Z <- as.formula(paste("~ 1 +", paste(nuisance.eff, collapse = " + ")))
  Z <- model.matrix(form.Z, data = sample.meta, na.action = na.action)
  rownames(Z) <- rn
  }
  
  # ---- X for extra core vars (not the contrasted factor) ----
  X <- NULL
  if (length(core.extra)) {
  X.list <- lapply(core.extra, function(v) {
    x <- sample.meta[[v]]
    if (is.factor(x) || is.character(x)) {
    f  <- factor(x)
    mm <- model.matrix(~ 0 + f, na.action = na.action)
    colnames(mm) <- lvlKey(v, levels(f))
    mm
    } else {
    mm <- matrix(x, ncol = 1); colnames(mm) <- v; mm
    }
  })
  X.list <- Filter(function(m) ncol(m) > 0, X.list)
  if (length(X.list)) X <- do.call(cbind, X.list)
  if (!is.null(X)) rownames(X) <- rn
  }
  
  # ---- Handle contrasted variable ----
  if (ctr$type == "factor") {
  g <- ctr$var
  f <- factor(sample.meta[[g]])
  L <- levels(f)
  mm.G <- model.matrix(~ 0 + f, na.action = na.action)
  colnames(mm.G) <- lvlKey(g, L); rownames(mm.G) <- rn
  
  keep <- intersect(L, names(ctr$weights))     # levels in contrast
  comp <- setdiff(L, keep)                     # non-participating levels
  if (!length(keep)) stop("Contrast references no existing levels of '", g, "'.")
  
  if (length(comp) >= 1) {
    # Keep all contrasted dummies in X; complement (minus one baseline) in Z
    base <- comp[1]
    X.g   <- mm.G[, lvlKey(g, keep), drop = FALSE]
    Z.g   <- mm.G[, lvlKey(g, setdiff(comp, base)), drop = FALSE]

    # Ensure an n×1 intercept in Z (if Z is NULL or lacks one)
    
    if (is.null(Z)) {
    Z <- makeIntercept(nrow(sample.meta), rn)
    } else if (!("(Intercept)" %in% colnames(Z))) {
    Z <- cbind(makeIntercept(nrow(sample.meta), rn), Z)
    }
    if (!is.null(Z.g)) Z <- cbind(Z, Z.g)
    rownames(Z) <- rn
    
    X <- if (is.null(X)) X.g else cbind(X, X.g)
    
    # contrasts (zero on Z)
    contrast.X <- setNames(numeric(ncol(X)), colnames(X))
    contrast.X[lvlKey(g, names(ctr$weights))] <- as.numeric(ctr$weights)

    nm.F <- concatNamesXZ(X, Z)
    contrast.F <- setNames(numeric(length(nm.F)), nm.F)
    contrast.F[names(contrast.X)] <- contrast.X
    
  } else {
    # Only the contrasted levels exist → single contrast regressor (weights must sum to 0)
    w <- as.numeric(ctr$weights)
    if (abs(sum(w)) > 1e-12)
    stop("Only contrasted levels present: contrast weights must sum to zero.")
    g.col <- as.numeric(mm.G[, lvlKey(g, names(ctr$weights)), drop = FALSE] %*% matrix(w, ncol = 1))
    X.g <- matrix(g.col, ncol = 1); colnames(X.g) <- paste0(g, ".contrast"); rownames(X.g) <- rn
    X  <- if (is.null(X)) X.g else cbind(X, X.g)
    
    contrast.X <- setNames(numeric(ncol(X)), colnames(X)); contrast.X[colnames(X.g)] <- 1

    nm.F <- concatNamesXZ(X, Z)
    contrast.F <- setNames(numeric(length(nm.F)), nm.F)
    contrast.F[colnames(X.g)] <- 1
  }
  
  } else {  # numeric contrasted variable
  v <- ctr$var
  x.v <- matrix(meta[[v]], ncol = 1); colnames(x.v) <- v; rownames(x.v) <- rn
  X  <- if (is.null(X)) x.v else cbind(X, x.v)
  
  contrast.X <- setNames(numeric(ncol(X)), colnames(X)); contrast.X[v] <- as.numeric(ctr$weights[1])

  nm.F <- concatNamesXZ(X, Z)
  contrast.F <- setNames(numeric(length(nm.F)), nm.F)
  contrast.F[v] <- contrast.X[v]
  }
  
  # ---- Full design with X first ----
  F <- if (!is.null(Z)) cbind(X, Z) else X
  rownames(F) <- rn
  
  # Pre-compute qrZ once (NULL if Z is NULL) for FL residualization
  qr.Z <- if (!is.null(Z)) qr(as.matrix(Z)) else NULL
  
  # Masks / blocks
  core.rows <- if (!is.null(Z) && ctr$type == "factor") {
  z <- sample.meta[[ctr$var]] %in% names(ctr$weights); z[is.na(z)] <- FALSE; z
  } else NULL
  blocks <- makeBlocks(sample.meta, nuisance = nuisance.eff, block.vars = block.vars)

  list(F = F, X = X, Z = Z, qr.Z = qr.Z,
     contrast.F = contrast.F,   # over colnames(F)
     contrast.X = contrast.X,   # over colnames(X)
     core.rows = core.rows, blocks = blocks)
}

#' @keywords internal
#' Parse a contrast specification
#'
#' Accepts:
#' 1) list(var="<name>", weights = named numeric),
#' 2) named numeric like c("var=LevelA" = -1, "var=LevelB" = +1),
#' 3) character triplet for simple two-level factor contrasts:
#'    - unnamed: c(var, alt, ref)
#'    - or named: c(var="...", alt="...", ref="...")  (alt gets +1, ref gets -1)
#'
#' @param contrast list / named numeric / length-3 character vector
#' @param meta data.frame of covariates (to infer factor vs numeric)
#' @return list(var=<name>, weights=<named numeric by LEVEL>, type="factor"/"numeric")
#' @keywords internal
parseContrast <- function(contrast, meta) {
  ## --- character triplet c(var, alt, ref) or named triplet ---
  if (is.character(contrast) && length(contrast) == 3L) {
    nm <- names(contrast)
    
    # allow var name aliases when named
    if (!is.null(nm)) {
      var.key <- if ("var" %in% nm) "var" else if ("factor" %in% nm) "factor" else if ("f" %in% nm) "f" else NULL
      if (!is.null(var.key) && all(c("alt","ref") %in% nm)) {
        var <- contrast[[var.key]]
        alt <- contrast[["alt"]]
        ref <- contrast[["ref"]]
      } else {
        # fall back to positional if names are partial/missing
        var <- contrast[[1]]; alt <- contrast[[2]]; ref <- contrast[[3]]
      }
    } else {
      # purely positional: (var, alt, ref)
      var <- contrast[[1]]; alt <- contrast[[2]]; ref <- contrast[[3]]
    }
    
    if (!(var %in% names(meta)))
      stop("Variable '", var, "' not found in 'meta'.")
    
    # Triplet is for factor/character variables only
    is.fac <- is.factor(meta[[var]]) || is.character(meta[[var]])
    if (!is.fac)
      stop("Triplet contrasts are only supported for factor/character variables. ",
           "Use list(var=..., weights=...) for numeric contrasts.")
    
    if (identical(alt, ref))
      stop("Triplet contrast has identical alt and ref levels.")
    
    levs <- levels(droplevels(factor(meta[[var]])))
    miss <- setdiff(c(alt, ref), levs)
    if (length(miss))
      stop("Levels not found in meta[['", var, "']]: ", paste(miss, collapse = ", "))
    
    w <- setNames(c(+1, -1), c(alt, ref))  # alt = +1, ref = -1
    return(list(var = var, weights = w, type = "factor"))
  }
  
  ## --- list(var=..., weights=...) ---
  if (is.list(contrast) && length(contrast$var) == 1 && !is.null(contrast$weights)) {
    var <- contrast$var
    stopifnot(var %in% names(meta))
    type <- if (is.factor(meta[[var]]) || is.character(meta[[var]])) "factor" else "numeric"
    return(list(var = var, weights = contrast$weights, type = type))
  }
  
  ## --- named numeric like c("var=LevelA"=-1,"var=LevelB"=+1) ---
  if (is.numeric(contrast) && !is.null(names(contrast))) {
    parts <- strsplit(names(contrast), "=", fixed = TRUE)
    vars  <- vapply(parts, `[`, character(1), 1)
    lvls  <- vapply(parts, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
    if (length(unique(vars)) != 1) stop("Contrast must reference a single variable.")
    var <- unique(vars)
    if (!(var %in% names(meta))) stop("Variable '", var, "' not found in 'meta'.")
    type <- if (is.factor(meta[[var]]) || is.character(meta[[var]])) "factor" else "numeric"
    w <- tapply(contrast, lvls, sum)  # collapse duplicates if any
    return(list(var = var, weights = w, type = type))
  }
  
  stop("Unsupported 'contrast' format. Use:\n",
       " - list(var=..., weights=...), or\n",
       " - named numeric like c('var=LevelA'=-1,'var=LevelB'=+1), or\n",
       " - character triplet c(var, alt, ref) (or named: var=.., alt=.., ref=..).")
}

#' Add an intercept column
#' @keywords internal
makeIntercept <- function(n, rn) {
    M <- matrix(1, nrow = n, ncol = 1); colnames(M) <- "(Intercept)"; rownames(M) <- rn; M
    }


#' Build a canonical column name "var=level" for factor level columns
#' @param var character scalar, variable name
#' @param lvl character/atomic scalar, level
#' @return character scalar like "batch=A"
#' @keywords internal
lvlKey <- function(var, lvl) paste0(var, "=", as.character(lvl))

#' Check for non-varying nuisance terms
#' @keywords internal
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

#' Safe concat of X/Z column names (works when Z is NULL, and even if X were NULL)
#' @keywords internal
concatNamesXZ <- function(X, Z) {
  c(if (!is.null(X)) colnames(X) else character(0),
    if (!is.null(Z)) colnames(Z) else character(0))
} 

#' Build blocking factor for restricted permutations
#' 
#' Default: interact all *factor* nuisance variables (continuous nuisances ignored).
#' You can override by passing a character vector of block_vars (only those used).
#' 
#' @param meta data.frame
#' @param nuisance character vector of nuisance variable names
#' @param block_vars optional character vector to override default
#' @return factor with one level per block (or single 'all' level if no factor nuisances)
#' @keywords internal
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




# Build a p×K contrast matrix from:
#   - X: model matrix used for fitting (preferably with intercept kept)
#   - formula: the model formula used to build X
#   - contrasts: list of c(var, ref, alt) triplets, e.g. list(c("diagnosis","Control","ASD"))
#   - expand:
#       "marginal" -> only main contrasts for each requested var
#       "simple"   -> main + simple effects at each non-baseline level of partner factors found in interactions
#       "both"     -> "simple" plus, when a partner is 2-level, include the difference-of-simple-effects (the pure interaction)
#
# Notes:
# - Works with treatment coding (the default from model.matrix). If the intercept is dropped, it still works.
# - Supports factor×factor and continuous×factor interactions. For >2-level partners, it generates simple effects per level.
#' @keywords internal
buildContrastMatrix <- function(X, formula, contrasts, expand = c("marginal","simple","both")) {
  stopifnot(is.matrix(X))
  expand <- match.arg(expand)
  p  <- ncol(X); cn <- colnames(X)

  # helpers
  # Main-effect dummy columns for a factor var (exclude interactions)
  mainLevelToCol <- function(var) {
    idx <- grep(paste0("^", var), cn)
    idx <- idx[!grepl(":", cn[idx], fixed = TRUE)]
    lev <- sub(paste0("^", var), "", cn[idx])
    keep <- nzchar(lev)
    setNames(idx[keep], lev[keep])
  }
  # Check if a variable is continuous in X (has a single main column named exactly var)
  contCol <- function(var) {
    j <- which(cn == var)
    if (length(j) == 1L) j else integer(0)
  }
  # Build the regex for an interaction column (order-agnostic)
  toK  <- function(var, lev) if (missing(lev) || is.null(lev)) var else paste0(var, lev)
  find.inter.col <- function(t1, t2) {
    which(cn == paste0(t1, ":", t2) | cn == paste0(t2, ":", t1))
  }
  add <- function(v, j, w) { if (length(j)==1 && j>0) v[j] <- v[j] + w; v }

  # Parse formula to discover interactions containing each requested var
  tr <- terms(formula)
  term.labels <- attr(tr, "term.labels")
  inter.pairs <- strsplit(term.labels[grepl(":", term.labels, fixed = TRUE)], ":", fixed = TRUE)

  partnersFor <- function(var) {
    if (length(inter.pairs) == 0) return(character())
    unique(unlist(lapply(inter.pairs, function(ab) {
      if (var %in% ab) setdiff(ab, var) else character()
    })))
  }

  # Normalize contrasts input into list of lists
  specs <- lapply(contrasts, function(x) {
    stopifnot(is.character(x), length(x) == 3)
    list(var = x[1], ref = x[2], alt = x[3])
  })

  C <- NULL; cnames <- character()

  for (sp in specs) {
    A <- sp$var; ref <- sp$ref; alt <- sp$alt
    v.main <- numeric(p); names(v.main) <- cn

    # ----- main contrast for A (alt vs ref) -----
    mapA <- mainLevelToCol(A)
    jA.ref <- unname(mapA[ref]); if (length(jA.ref)==0) jA.ref <- NA_integer_
    jA.alt <- unname(mapA[alt]); if (length(jA.alt)==0) jA.alt <- NA_integer_
    jA.cont <- contCol(A)

    if (length(jA.cont) == 1L) {
      # A is continuous -> "main contrast" is just slope of A
      v.main <- add(v.main, jA.cont, +1)
      main.name <- paste0("Slope(", A, ")")
    } else {
      # A is factor
      if (!is.na(jA.alt) && !is.na(jA.ref)) {
        v.main <- add(v.main, jA.alt, +1); v.main <- add(v.main, jA.ref, -1)
      } else if (!is.na(jA.alt) && is.na(jA.ref)) {
        v.main <- add(v.main, jA.alt, +1)          # ref is baseline (intercept coding)
      } else if (is.na(jA.alt) && !is.na(jA.ref)) {
        v.main <- add(v.main, jA.ref, -1)          # alt is baseline
      } else {
        stop("For '", A, "': neither '", ref, "' nor '", alt, "' has a main column in X.")
      }
      main.name <- paste0(A, alt, "_vs_", ref)
    }

    # always include the marginal main contrast
    C <- cbind(C, v.main); cnames <- c(cnames, main.name)

    # ---------- expansions through interactions ----------
    if (expand != "marginal") {
      partners <- partnersFor(A)
      for (B in partners) {
        # continuous partner?
        jB.cont <- contCol(B)
        if (length(jB.cont) == 1L) {
          # A × continuous partner is uncommon for "simple effect"; skip by default
          next
        }

        # factor partner: simple effects at each non-baseline level of B
        mapB <- mainLevelToCol(B)
        if (length(mapB) == 0) next  # B not factor-coded in X

        # Determine B's non-baseline levels present as columns
        levB <- names(mapB)  # these are the dummy-coded levels (exclude baseline)
        for (b in levB) {
          v.simple <- v.main  # start from main effect of A
          # add interaction column for A_alt : B_b (or slope(A) : B_b if A is continuous)
          if (length(jA.cont) == 1L) {
            jInt <- find.inter.col(toK(A), toK(B, b))
            if (length(jInt) == 1L) v.simple <- add(v.simple, jInt, +1)
            simple.name <- paste0(main.name, " | ", B, "=", sub(paste0("^", B), "", b))
          } else {
            jInt <- find.inter.col(toK(A, alt), toK(B, b))
            if (length(jInt) == 1L) v.simple <- add(v.simple, jInt, +1)
            simple.name <- paste0("Simple(", A, " ", alt, "_vs_", ref, " | ",
                                  B, "=", sub(paste0("^", B), "", b), ")")
          }
          C <- cbind(C, v.simple); cnames <- c(cnames, simple.name)
        }

        # Differences of simple effects (pure interaction) when B has exactly 1 dummy (i.e., 2 levels)
        if (expand == "both" && length(levB) == 1L) {
          b <- levB[1]
          v.diff <- numeric(p); names(v.diff) <- cn
          if (length(jA.cont) == 1L) {
            # slope(A|B=b) - slope(A|B=baseline) = beta_{A:B_b}
            jInt <- find.inter.col(toK(A), toK(B, b))
            if (length(jInt) == 1L) v.diff <- add(v.diff, jInt, +1)
            diff.name <- paste0("DiffSlope(", A, " | ", B, "=", sub(paste0("^", B), "", b), " vs baseline)")
          } else {
            # [A_alt vs ref at B=b] - [A_alt vs ref at baseline] = beta_{A_alt:B_b}
            jInt <- find.inter.col(toK(A, alt), toK(B, b))
            if (length(jInt) == 1L) v.diff <- add(v.diff, jInt, +1)
            diff.name <- paste0("Interaction(", A, " ", alt, " × ", B, "=", sub(paste0("^", B), "", b), ")")
          }
          if (any(v.diff != 0)) { C <- cbind(C, v.diff); cnames <- c(cnames, diff.name) }
        }
      }
    }
  }

  if (is.null(C)) C <- matrix(numeric(0), nrow = p, ncol = 0, dimnames = list(colnames(X), character()))
  rownames(C) <- colnames(X); colnames(C) <- cnames
  C
}