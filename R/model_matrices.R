# code building model matrices

# -----------------------------------------------------------------------------
# Small design builder for permutation workflows (e.g., Freedman–Lane)
# - Core X defaults to only the factor levels participating in the contrast
# - Nuisance Z defaults to all other columns in `meta` (with intercept)
# - Blocks = interaction of nuisance *factor* variables (continuous nuisances ignored)
# - Returns: F (full), X, Z, contrast_coef (vector over coef(F)), core.rows, blocks
# -----------------------------------------------------------------------------


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
#' @param block_vars optional character vector to override block definition
#' @param na_action NA handler passed to model.matrix (default stats::na.pass)
#' @return list(F, X, Z, contrast_coef, core.rows, blocks)
#'   F, X, Z - full, core, and nuisance model matrices
#'   qrZ - QR of Z (for residualization based on Z)
#'   contrast_F - linear combination implementing desired contrast on the
#'      current coefficient order relative to F. i.e. to get contrast value:
#'      crossprod(coef(fit), contrast_coef)
#'   contrast_X - same thing, but for X
#'   core.rows - logical vector specifying a subset of rows which are needed
#'      to fit the core coefficients (for optimizing FL)
#'   blocks - factor on rows specifying blocks in which randomization should be
#'      performed (if not using FL)
# helper used here; keep if not already in scope
# helper for "var=level" names
prepare_design <- function(meta,
                           contrast,              # list(var,weights) or named "var=level"
                           nuisance   = NULL,     # default: all other columns
                           core_extra = NULL,     # optional other core variables
                           block_vars = NULL,
                           na_action  = stats::na.pass) {
  stopifnot(is.data.frame(meta))
  
  ctr <- .parse_contrast(contrast, meta)
  if (is.null(nuisance)) nuisance <- setdiff(names(meta), c(ctr$var, core_extra))
  if (is.null(core_extra)) core_extra <- character(0)
  
  rn <- rownames(meta); if (is.null(rn)) rn <- as.character(seq_len(nrow(meta)))
  
  # keep only nuisance vars that actually vary (avoids 1-level factor error)
  nuisance_eff <- .filter_effective_nuisance(meta, nuisance)
  
  # ---- Z (with intercept iff any effective nuisance) ----
  Z <- NULL
  if (length(nuisance_eff)) {
    formZ <- as.formula(paste("~ 1 +", paste(nuisance_eff, collapse = " + ")))
    Z <- model.matrix(formZ, data = meta, na.action = na_action)
    rownames(Z) <- rn
  }
  
  # ---- X for extra core vars (not the contrasted factor) ----
  X <- NULL
  if (length(core_extra)) {
    X_list <- lapply(core_extra, function(v) {
      x <- meta[[v]]
      if (is.factor(x) || is.character(x)) {
        f  <- factor(x)
        mm <- model.matrix(~ 0 + f, na.action = na_action)
        colnames(mm) <- .lvl_key(v, levels(f))
        mm
      } else {
        mm <- matrix(x, ncol = 1); colnames(mm) <- v; mm
      }
    })
    X_list <- Filter(function(m) ncol(m) > 0, X_list)
    if (length(X_list)) X <- do.call(cbind, X_list)
    if (!is.null(X)) rownames(X) <- rn
  }
  
  # ---- Handle contrasted variable ----
  if (ctr$type == "factor") {
    g <- ctr$var
    f <- factor(meta[[g]])
    L <- levels(f)
    mmG <- model.matrix(~ 0 + f, na.action = na_action)
    colnames(mmG) <- .lvl_key(g, L); rownames(mmG) <- rn
    
    keep <- intersect(L, names(ctr$weights))     # levels in contrast
    comp <- setdiff(L, keep)                     # non-participating levels
    if (!length(keep)) stop("Contrast references no existing levels of '", g, "'.")
    
    if (length(comp) >= 1) {
      # Keep all contrasted dummies in X; complement (minus one baseline) in Z
      base <- comp[1]
      Xg   <- mmG[, .lvl_key(g, keep), drop = FALSE]
      Zg   <- mmG[, .lvl_key(g, setdiff(comp, base)), drop = FALSE]
      
      # Ensure an n×1 intercept in Z (if Z is NULL or lacks one)
      make_intercept <- function(n, rn) {
        M <- matrix(1, nrow = n, ncol = 1); colnames(M) <- "(Intercept)"; rownames(M) <- rn; M
      }
      if (is.null(Z)) {
        Z <- make_intercept(nrow(meta), rn)
      } else if (!("(Intercept)" %in% colnames(Z))) {
        Z <- cbind(make_intercept(nrow(meta), rn), Z)
      }
      if (!is.null(Zg)) Z <- cbind(Z, Zg)
      rownames(Z) <- rn
      
      X <- if (is.null(X)) Xg else cbind(X, Xg)
      
      # contrasts (zero on Z)
      contrast_X <- setNames(numeric(ncol(X)), colnames(X))
      contrast_X[.lvl_key(g, names(ctr$weights))] <- as.numeric(ctr$weights)
      
      nmF <- .concat_names_XZ(X, Z)                                    
      contrast_F <- setNames(numeric(length(nmF)), nmF)                
      contrast_F[names(contrast_X)] <- contrast_X
      
    } else {
      # Only the contrasted levels exist → single contrast regressor (weights must sum to 0)
      w <- as.numeric(ctr$weights)
      if (abs(sum(w)) > 1e-12)
        stop("Only contrasted levels present: contrast weights must sum to zero.")
      gcol <- as.numeric(mmG[, .lvl_key(g, names(ctr$weights)), drop = FALSE] %*% matrix(w, ncol = 1))
      Xg <- matrix(gcol, ncol = 1); colnames(Xg) <- paste0(g, "_contrast"); rownames(Xg) <- rn
      X  <- if (is.null(X)) Xg else cbind(X, Xg)
      
      contrast_X <- setNames(numeric(ncol(X)), colnames(X)); contrast_X[colnames(Xg)] <- 1
      
      nmF <- .concat_names_XZ(X, Z)                                    
      contrast_F <- setNames(numeric(length(nmF)), nmF)                
      contrast_F[colnames(Xg)] <- 1
    }
    
  } else {  # numeric contrasted variable
    v <- ctr$var
    xv <- matrix(meta[[v]], ncol = 1); colnames(xv) <- v; rownames(xv) <- rn
    X  <- if (is.null(X)) xv else cbind(X, xv)
    
    contrast_X <- setNames(numeric(ncol(X)), colnames(X)); contrast_X[v] <- as.numeric(ctr$weights[1])
    
    nmF <- .concat_names_XZ(X, Z)                                       
    contrast_F <- setNames(numeric(length(nmF)), nmF)                   
    contrast_F[v] <- contrast_X[v]
  }
  
  # ---- Full design with X first ----
  F <- if (!is.null(Z)) cbind(X, Z) else X
  rownames(F) <- rn
  
  # Pre-compute qrZ once (NULL if Z is NULL) for FL residualization
  qrZ <- if (!is.null(Z)) qr(as.matrix(Z)) else NULL
  
  # Masks / blocks
  core.rows <- if (!is.null(Z) && ctr$type == "factor") {
    z <- meta[[ctr$var]] %in% names(ctr$weights); z[is.na(z)] <- FALSE; z
  } else NULL
  blocks <- .make_blocks(meta, nuisance = nuisance_eff, block_vars = block_vars)
  
  list(F = F, X = X, Z = Z, qrZ = qrZ,
       contrast_F = contrast_F,   # over colnames(F)
       contrast_X = contrast_X,   # over colnames(X)
       core.rows = core.rows, blocks = blocks)
}



# ---------- Residualize for FL (no add-back) ----------
# Uses qrZ computed once in prepare_design_dual(). If qrZ is NULL, returns inputs unchanged.
residualize_for_FL <- function(y, qrZ, X) {
  X <- as.matrix(X)
  if (is.null(qrZ)) return(list(y_r = y, X_r = X))
  y_r <- qr.resid(qrZ, y)
  X_r <- qr.resid(qrZ, X)
  list(y_r = y_r, X_r = X_r)
}


# Prepare pair-level design for modeling lower-tri distances
#
# Builds per-pair metadata (factor pairs "AA","AB",...; numeric diffs |xi-xj|),
# constructs a pair-level contrast from a sample-level triplet (var, alt, ref)
# according to dist.type, and calls `prepare_design()` to produce F, X, Z, qrZ,
# contrast vectors, core.rows, and blocks. Returns the (i,j) index table so you
# can vectorize ANY distance matrix later in the same order.
#
# @param meta data.frame of per-sample covariates; row order defines sample order
# @param triplet character triplet c(var, alt, ref) (or named var/alt/ref)
# @param dist.type one of "shift","total","var" (pair-level contrast form)
# @param block_vars optional character vector to override permutation blocks
# @param na_action NA handler for model.matrix
# @return list with:
#   - pairs: data.frame of (i,j) indices (1-based) in lower-tri order
#   - pair_meta: data.frame of per-pair covariates
#   - prep: design list from prepare_design() (F, X, Z, qrZ, contrast_F/_X, ...)
prepare_pair_design <- function(meta, triplet,
                                dist.type = c("shift","total","var"),
                                block_vars = NULL,
                                na_action = stats::na.pass) {
  stopifnot(is.data.frame(meta))
  dist.type <- match.arg(dist.type)
  
  # 1) Indices for lower triangle, derived from nrow(meta)
  n   <- nrow(meta)
  idx <- .lower_tri_indices(n)                 # expects helper already defined
  
  # 2) Pair-level metadata from sample-level meta
  pair_meta <- .pairify_meta(meta, idx)        # expects helper already defined
  
  # 3) Triplet (var, alt, ref) → pair-level contrast weights over AA/AB/BB
  tr <- .parse_triplet(triplet)                # var, alt, ref (your existing helper)
  if (!(tr$var %in% names(meta)))
    stop("Triplet 'var' not found in meta: ", tr$var)
  
  f_levels <- levels(factor(meta[[tr$var]]))
  if (!all(c(tr$ref, tr$alt) %in% f_levels))
    stop("alt/ref levels not found in meta[['", tr$var, "']]: ",
         tr$alt, ", ", tr$ref)
  
  weights <- .pair_contrast_weights(tr$ref, tr$alt, dist.type, level_order = f_levels)
  
  # 4) Build pair-level variable name for the contrasted factor (e.g., "group_pair")
  pair_var <- paste0(tr$var, "_pair")
  if (!(pair_var %in% names(pair_meta)))
    stop("Pair variable not found in pair_meta: ", pair_var,
         " (did '", tr$var, "' exist and was it factor-like?)")
  
  # 5) Call your prepare_design() on pair-level metadata
  prep <- prepare_design(
    meta       = pair_meta,
    contrast   = list(var = pair_var, weights = weights),
    #nuisance   = setdiff(names(pair_meta), pair_var),  # core is implicitly pair_var
    core_extra = NULL,
    block_vars = block_vars,
    na_action  = na_action
  )
  
  list(pairs = idx, pair_meta = pair_meta, model = prep)
}

# Vectorize the lower triangle of a square matrix using a precomputed index order
# If pairs is NULL, use standard lower.tri order for that matrix.
vectorize_lower_tri <- function(d, pairs = NULL) {
  stopifnot(is.matrix(d), nrow(d) == ncol(d))
  if (is.null(pairs)) {
    ij <- which(lower.tri(d), arr.ind = TRUE)
    return(d[cbind(ij[,1], ij[,2])])
  } else {
    stopifnot(all(c("i","j") %in% names(pairs)),
              nrow(d) >= max(pairs$i, pairs$j))
    return(d[cbind(pairs$i, pairs$j)])
  }
}



# Return minimal metadata for prepare_design(): only RHS vars from a formula
# - Expands '.' to all columns except LHS
# - Handles interactions/transforms via all.vars()
# - Optionally drop unused factor levels
select_minimal_meta <- function(meta, formula, drop_unused_levels = TRUE, extra = NULL) {
  stopifnot(is.data.frame(meta), inherits(formula, "formula"))
  lhs <- if (length(formula) >= 2) all.vars(formula[[2]]) else character(0)
  rhs_expr <- formula[[3]]
  rhs_vars <- all.vars(rhs_expr)
  
  if (any(rhs_vars == ".")) rhs_vars <- setdiff(names(meta), lhs)
  if (!is.null(extra) && length(extra)) rhs_vars <- union(rhs_vars, extra)
  
  miss <- setdiff(rhs_vars, names(meta))
  if (length(miss)) stop("Missing columns in 'meta': ", paste(miss, collapse = ", "))
  
  out <- meta[, rhs_vars, drop = FALSE]
  rownames(df) <- rownames(meta)
  if (drop_unused_levels) {
    for (nm in names(out)) if (is.factor(out[[nm]])) out[[nm]] <- droplevels(out[[nm]])
  }
  out
}


### Helper methods

#' Build a canonical column name "var=level" for factor level columns
#' @param var character scalar, variable name
#' @param lvl character/atomic scalar, level
#' @return character scalar like "batch=A"
.lvl_key <- function(var, lvl) paste0(var, "=", as.character(lvl))

# check for non-varying nuisance terms
.filter_effective_nuisance <- function(meta, nuisance_names) {
  if (!length(nuisance_names)) return(character(0))
  keep <- vapply(nuisance_names, function(v) {
    x <- meta[[v]]
    if (is.factor(x) || is.character(x)) length(levels(droplevels(factor(x)))) >= 2L
    else if (is.numeric(x)) length(unique(x[!is.na(x)])) >= 2L
    else FALSE
  }, logical(1))
  nuisance_names[keep]
}

# safe concat of X/Z column names (works when Z is NULL, and even if X were NULL)
.concat_names_XZ <- function(X, Z) {
  c(if (!is.null(X)) colnames(X) else character(0),
    if (!is.null(Z)) colnames(Z) else character(0))
}

# Choose which contrasted level to drop (the "reference")
# - If drop_level is supplied, use it (must be among referenced levels)
# - Else pick the first negative-weight level in factor-level order, or
#   if none are negative, the first referenced level in factor-level order
.choose_drop_level <- function(ctr_weights, factor_levels, drop_level = NULL) {
  keep <- intersect(factor_levels, names(ctr_weights))
  if (!length(keep)) stop("Contrast references no existing levels.")
  if (!is.null(drop_level)) {
    if (!drop_level %in% keep) stop("drop_level '", drop_level, "' is not among contrasted levels.")
    return(drop_level)
  }
  w <- ctr_weights[keep]
  neg <- keep[w < 0]
  if (length(neg)) {
    # choose first negative in factor-level order
    return(neg[order(match(neg, factor_levels))][1])
  } else {
    # no negatives: choose first referenced level in factor-level order
    return(keep[order(match(keep, factor_levels))][1])
  }
}

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
.parse_contrast <- function(contrast, meta) {
  ## --- character triplet c(var, alt, ref) or named triplet ---
  if (is.character(contrast) && length(contrast) == 3L) {
    nm <- names(contrast)
    
    # allow var name aliases when named
    if (!is.null(nm)) {
      var_key <- if ("var" %in% nm) "var" else if ("factor" %in% nm) "factor" else if ("f" %in% nm) "f" else NULL
      if (!is.null(var_key) && all(c("alt","ref") %in% nm)) {
        var <- contrast[[var_key]]
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
    is_fac <- is.factor(meta[[var]]) || is.character(meta[[var]])
    if (!is_fac)
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




#' Build blocking factor for restricted permutations
#' 
#' Default: interact all *factor* nuisance variables (continuous nuisances ignored).
#' You can override by passing a character vector of block_vars (only those used).
#' 
#' @param meta data.frame
#' @param nuisance character vector of nuisance variable names
#' @param block_vars optional character vector to override default
#' @return factor with one level per block (or single 'all' level if no factor nuisances)
.make_blocks <- function(meta, nuisance, block_vars = NULL) {
  if (!is.null(block_vars) && length(block_vars)) {
    return(interaction(lapply(meta[block_vars], as.factor), drop = TRUE, lex.order = TRUE))
  }
  fac_nuis <- if (length(nuisance)) {
    nuisance[vapply(meta[nuisance], function(x) is.factor(x) || is.character(x), logical(1))]
  } else character(0)
  if (length(fac_nuis)) interaction(lapply(meta[fac_nuis], as.factor), drop = TRUE, lex.order = TRUE)
  else factor(rep("all", nrow(meta)))
}


# -------------------------------------------------------------------------
# Pairwise design builder for lower-tri modeling of a distance matrix
# Depends on your existing `prepare_design()` defined earlier.
# -------------------------------------------------------------------------

#' Return lower-triangular (row, col) index pairs for an n x n matrix
#' @param n integer number of samples
#' @return data.frame with columns i, j (1-based), in row-major lower-tri order
.lower_tri_indices <- function(n) {
  ij <- which(lower.tri(matrix(NA_real_, n, n)), arr.ind = TRUE)
  data.frame(i = ij[,1], j = ij[,2])
}

#' Symmetric pair-code for a factor: "AA", "AB", "BB" (AB ≡ BA)
#' Uses the factor's level order to define A<B for concatenation.
#' @param f factor vector (length n)
#' @param i integer indices for first sample in each pair
#' @param j integer indices for second sample in each pair
#' @return character vector of codes ("AA","AB",...)
.pair_code_factor <- function(f, i, j) {
  lev <- levels(f); fi <- as.integer(f[i]); fj <- as.integer(f[j])
  any_na <- is.na(fi) | is.na(fj)
  lo <- pmin(fi, fj); hi <- pmax(fi, fj)
  code <- paste0(lev[lo], lev[hi])
  code[any_na] <- NA_character_
  code
}

#' Symmetric numeric pair difference |x_i - x_j|
#' @param x numeric vector (length n)
#' @param i, j integer index vectors
#' @return numeric vector of absolute differences
.pair_diff_numeric <- function(x, i, j) {
  xi <- x[i]; xj <- x[j]
  out <- abs(xi - xj)
  out
}

#' Parse a sample-level "triplet" contrast specification
#' @param triplet list(var="<factor>", ref="<A>", alt="<B>")
#' @return list(var, ref, alt)
.parse_triplet <- function(triplet) {
  stopifnot(length(triplet)==3)
  triplet <- as.list(triplet)
  names(triplet) <- c('var','alt','ref')
  if (triplet$ref == triplet$alt) stop("ref and alt must be different.")
  list(var = as.character(triplet$var),
       alt = as.character(triplet$alt),
       ref = as.character(triplet$ref))
}

#' Build pair-level metadata from sample-level metadata in lower-tri order
#' @param meta data.frame of per-sample covariates (rows correspond to d's order)
#' @param idx data.frame with columns i, j (lower-tri indices)
#' @return data.frame of per-pair covariates
.pairify_meta <- function(meta, idx) {
  stopifnot(is.data.frame(meta), all(c("i","j") %in% names(idx)))
  n <- nrow(meta)
  stopifnot(all(idx$i >= 1 & idx$i <= n & idx$j >= 1 & idx$j <= n & idx$i > idx$j))
  
  out <- list()
  for (v in names(meta)) {
    x <- meta[[v]]
    if (is.factor(x) || is.character(x) || is.logical(x)) {
      f <- factor(x)
      out[[paste0(v, "_pair")]] <- factor(.pair_code_factor(f, idx$i, idx$j))
    } else if (is.numeric(x)) {
      out[[paste0(v, "_diff")]] <- .pair_diff_numeric(x, idx$i, idx$j)
    } else {
      # unsupported type → coerce to factor then pair-code
      f <- factor(as.character(x))
      out[[paste0(v, "_pair")]] <- factor(.pair_code_factor(f, idx$i, idx$j))
    }
  }
  as.data.frame(out, stringsAsFactors = TRUE)
}

# Translate (ref, alt) + dist.type into pair-level contrast weights
# - level_order: optional character vector giving the factor's level order
#   (used to decide the label for the cross pair, e.g., "AB" vs "BA")
# - Returns a named numeric vector, e.g. c("AA"=-0.5, "AB"=1, "BB"=-0.5)
.pair_contrast_weights <- function(ref, alt, type = c("shift","total","var"), level_order = NULL) {
  type <- match.arg(type)
  
  # Cross label uses factor-level order if provided; otherwise alphabetical
  two <- c(ref, alt)
  if (!is.null(level_order)) {
    ord <- two[order(match(two, level_order))]
  } else {
    ord <- sort(two)  # fallback: alphabetical
  }
  AB <- paste0(ord[1], ord[2])
  
  AA <- paste0(ref, ref)
  BB <- paste0(alt, alt)
  
  if (type == "shift") {
    w <- c(-0.5, 1, -0.5)
    names(w) <- c(AA, AB, BB)
    return(w)
  } else if (type == "total") {
    w <- c(-1, 1)
    names(w) <- c(AA, AB)
    return(w)
  } else { # type == "var"
    w <- c(-1, 1)
    names(w) <- c(AA, BB)
    return(w)
  }
}



