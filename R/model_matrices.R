# code building model matrices

# -----------------------------------------------------------------------------
# Small design builder for permutation workflows (e.g., Freedman–Lane)
# - Core X defaults to only the factor levels participating in the contrast
# - Nuisance Z defaults to all other columns in `meta` (with intercept)
# - Blocks = interaction of nuisance *factor* variables (continuous nuisances ignored)
# - Returns: F (full), X, Z, contrast_coef (vector over coef(F)), core.rows, blocks, model diagnostics
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
  mmG <- model.matrix(~ 0 + f, na.action = na.action)
  colnames(mmG) <- lvlKey(g, L); rownames(mmG) <- rn
  
  keep <- intersect(L, names(ctr$weights))     # levels in contrast
  comp <- setdiff(L, keep)                     # non-participating levels
  if (!length(keep)) stop("Contrast references no existing levels of '", g, "'.")
  
  if (length(comp) >= 1) {
    # Keep all contrasted dummies in X; complement (minus one baseline) in Z
    base <- comp[1]
    X.g   <- mmG[, lvlKey(g, keep), drop = FALSE]
    Z.g   <- mmG[, lvlKey(g, setdiff(comp, base)), drop = FALSE]

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
    g.col <- as.numeric(mmG[, lvlKey(g, names(ctr$weights)), drop = FALSE] %*% matrix(w, ncol = 1))
    X.g <- matrix(g.col, ncol = 1); colnames(X.g) <- paste0(g, "_contrast"); rownames(X.g) <- rn
    X  <- if (is.null(X)) X.g else cbind(X, X.g)
    
    contrast.X <- setNames(numeric(ncol(X)), colnames(X)); contrast.X[colnames(X.g)] <- 1

    nm.F <- concatNamesXZ(X, Z)
    contrast.F <- setNames(numeric(length(nm.F)), nm.F)
    contrast.F[colnames(X.g)] <- 1
  }
  
  } else {  # numeric contrasted variable
  v <- ctr$var
  x.v <- matrix(sample.meta[[v]], ncol = 1); colnames(x.v) <- v; rownames(x.v) <- rn
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
  qrZ <- if (!is.null(Z)) qr(as.matrix(Z)) else NULL
  
  # Masks / blocks
  core.rows <- if (!is.null(Z) && ctr$type == "factor") {
  z <- sample.meta[[ctr$var]] %in% names(ctr$weights); z[is.na(z)] <- FALSE; z
  } else NULL
  blocks <- makeBlocks(sample.meta, nuisance = nuisance.eff, block.vars = block.vars)

  # Build permutation groups
  pg <- permutationGroups(blocks = blocks, core.rows = if (is.null(core.rows)) rep(TRUE, nrow(F)) else core.rows)
  
  # Run diagnostics & warnings (includes permutation diagnostics)
  diag <- diagnoseDesign(F = F, X = X, Z = Z, qrZ = qrZ,
                           meta = sample.meta, blocks = blocks, core.rows = core.rows, ctr = ctr,
                           verbose = TRUE)

  list(F = F, X = X, Z = Z, qrZ = qrZ,
     contrast.F = contrast.F,   # over colnames(F)
     contrast.X = contrast.X,   # over colnames(X)
     core.rows = core.rows, blocks = blocks, 
     perm.groups.full = pg$full,
     perm.groups.core = pg$core,
     model.diag = diag)
}

# Residualize for Freedman-Lane (no add-back) OLD
# Uses qrZ computed once. If qrZ is NULL, returns inputs unchanged.
#' @keywords internal
residualizeForFL <- function(y, qrZ, X) {
  X <- as.matrix(X)
  if (is.null(qrZ)) return(list(y.r = y, X.r = X))
  y.r <- qr.resid(qrZ, y)
  X.r <- qr.resid(qrZ, X)
  list(y.r = y.r, X.r = X.r)
}


# Split rows into permutation groups (freeze-only policy). No filtering here.
# - blocks: factor length n (interaction of nuisance factors)
# - core.rows: logical length n; if NULL, all rows are treated as core
# Returns indices both for the full data and for the core subset (relative to F[core.rows, ])
permutationGroups <- function(blocks, core.rows = NULL) {
  stopifnot("blocks is not a factor"=is.factor(blocks), "no randomization blocks found"=length(blocks) >= 1L)
  n <- length(blocks)
  
  # Full-data groups: indices are 1..n
  if(is.null(core.rows)) {
    core.rows <- rep(TRUE,n)
    groups.core <- groups_full <- split(seq_len(n), droplevels(blocks), drop = TRUE)
  } else {
    stopifnot("core.rows is not logi"=is.logical(core.rows), "core.rows length mismatch"=length(core.rows) == n)
    groups.full <- split(seq_len(n)[core.rows], droplevels(blocks[core.rows]), drop = TRUE)
    core.idx    <- which(core.rows)
    blocks.core <- droplevels(blocks[core.rows])
    # map full index -> position in core subset
    pos.in.core <- integer(n); pos.in.core[core.idx] <- seq_along(core.idx)
    groups.core <- lapply(split(core.idx, blocks.core, drop = TRUE),
                          function(ids_full) pos.in.core[ids_full])
  }

  list(full = groups.full, core = groups.core)
}



# Diagnose model matrices and permutation feasibility (freeze-only).
# Prints only warnings/hints; returns all metrics & messages.
diagnoseDesign <- function(F = NULL, X = NULL, Z = NULL, qrZ = NULL,
               meta = NULL, blocks = NULL, core.rows = NULL, ctr = NULL,
               block.factors = NULL,      # names of factor vars used to make 'blocks'
               verbose = TRUE,
               thresholds = list(
                 alias.tol     = 1e-8,
                 kappa.warn    = 1e3,    # warn if >= this (and stronger wording if >>)
                 vif.warn      = 10,     # warn if any VIF >= this
                 min.core.size = 2L,     # movable if >= 2
                 show.top      = 10,     # how many problematic blocks to list
                 min.eff.perm  = 100     # warn if effective permutations < this
               )) {
  # ---------- tiny helpers ----------
  emitInternal <- function(lines) { if (verbose && length(lines)) for (ln in lines) message(ln); invisible(NULL) }
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
  if (is.null(X) || !is.matrix(X) || !ncol(X)) return(list(aliased=character(0), X.r=X))
  if (is.null(Z)) return(list(aliased=character(0), X.r=X))
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
  vif <- function(Xr) {
  if (is.null(Xr) || !is.matrix(Xr) || ncol(Xr) <= 1) return(numeric(0))
  Xc <- scale(Xr, center=TRUE, scale=TRUE); R <- stats::cor(Xc)
  VIF <- tryCatch(diag(solve(R)), error=function(e) diag(pinvSym(R)))
  setNames(as.numeric(VIF), colnames(Xr))
  }

  permCoreSummary <- function(meta, blocks, core.rows, ctr, groups.core, block.factors, show.top) {
  labs   <- names(groups.core); if (is.null(labs)) labs <- as.character(seq_along(groups.core))
  n.core <- vapply(groups.core, length, integer(1))
  
  # Variation of contrast within each core-block
  has.ctr <- !is.null(ctr) && !is.null(meta) && (ctr$var %in% names(meta))
  if (has.ctr && ctr$type == "factor") {
    keep <- names(ctr$weights)
    n.levels <- integer(length(groups.core))
    comp.str <- character(length(groups.core))
    for (i in seq_along(groups.core)) {
    gidx <- groups.core[[i]]
    ids  <- if (is.null(core.rows) || all(core.rows)) gidx else which(core.rows)[gidx]
    v    <- as.character(meta[[ctr$var]][ids]); v <- v[v %in% keep]
    tab  <- sort(table(v), decreasing = TRUE)
    n.levels[i] <- length(tab)
    comp.str[i] <- if (length(tab)) paste(sprintf("%s:%d", names(tab), as.integer(tab)), collapse=",") else ""
    }
  } else if (has.ctr && ctr$type == "numeric") {
    n.levels <- vapply(groups.core, function(gidx) {
    ids <- if (is.null(core.rows) || all(core.rows)) gidx else which(core.rows)[gidx]
    length(unique(meta[[ctr$var]][ids]))
    }, integer(1))
    comp.str <- rep("", length(groups.core))
  } else {
    n.levels <- rep(NA_integer_, length(groups.core))
    comp.str <- rep("", length(groups.core))
  }
  
  df <- data.frame(block = labs, n.core = n.core, n.levels = n.levels,
           composition = comp.str, stringsAsFactors = FALSE)
  
  # Identify problematic blocks
  too.small <- df$n.core < 2L
  no.var    <- !is.na(df$n.levels) & (df$n.levels < 2L)
  prob.idx  <- which(too.small | no.var)
  
  # Decode block factor combinations if requested
  comb.df <- NULL
  if (length(prob.idx) && !is.null(block.factors) && length(block.factors) &&
    all(block.factors %in% names(meta))) {
    # pick any representative row in the *full* data for each block
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
    # keep at most show.top
    if (!is.null(comb.df)) {
    ord <- order(comb.df$n.core, comb.df$n.levels)
    comb.df <- comb.df[ord, , drop = FALSE]
    comb.df <- head(comb.df, show.top)
    }
  }
  
  list(summary = df,
     prob.idx = prob.idx,
     comb.df  = comb.df)
  }
  
  # ---------- normalize inputs ----------
  if (is.null(F)) {
  if (is.null(X) && is.null(Z)) stop("Provide F or (X and/or Z).")
  if (!is.null(X)) X <- as.matrix(X)
  if (!is.null(Z)) Z <- as.matrix(Z)
  F <- if (!is.null(X) && !is.null(Z)) cbind(X, Z) else if (!is.null(X)) X else Z
  } else F <- as.matrix(F)
  
  # ---------- model diagnostics ----------
  msgs <- character(0)
  dF  <- qrRankDiag(F, name = "F")
  aXZ <- aliasXbyZ(X, Z, qrZ = qrZ, tol = thresholds$alias.tol)
  dXr <- qrRankDiag(aXZ$X.r, name = "M_Z X")
  vif <- vif(aXZ$X.r)
  
  if (dF$rank < dF$p) {
  msgs <- c(msgs, sprintf("[WARN] Design is rank-deficient: rank(F)=%d < %d. Dependent columns: %s",
              dF$rank, dF$p, paste(dF$dep, collapse = ", ")))
  }
  if (is.finite(dF$kappa) && dF$kappa >= thresholds$kappa.warn) {
  msgs <- c(msgs, sprintf("[WARN] F ill-conditioned (kappa=%.2e). Estimates may be unstable.", dF$kappa))
  }
  if (!is.null(X)) {
  if (length(aXZ$aliased)) {
    msgs <- c(msgs, paste0("[WARN] Core columns aliased by Z (not estimable after adjustment): ",
               paste(aXZ$aliased, collapse = ", ")))
  }
  if (dXr$rank < dXr$p) {
    msgs <- c(msgs, sprintf("[WARN] Core after Z is rank-deficient: rank(M_Z X)=%d < %d. Dependent core columns: %s",
                dXr$rank, dXr$p, paste(dXr$dep, collapse = ", ")))
  }
  if (length(vif)) {
    bad <- vif[vif >= thresholds$vif.warn]
    if (length(bad)) {
    msgs <- c(msgs, paste0("[WARN] High VIFs in core after Z: ",
                 paste(sprintf("%s=%.1f", names(bad), bad), collapse = ", "),
                 ". Consider collapsing levels or using a single contrast regressor."))
    }
  }
  }
  
  # ---------- permutation diagnostics (freeze policy) ----------
  perm <- NULL
  if (!is.null(blocks)) {
  stopifnot(is.factor(blocks), nrow(F) == length(blocks))
  # groups must be built outside; but if not provided, recreate core groups here:
  groups.core <- permutationGroups(blocks, if (is.null(core.rows)) rep(TRUE, nrow(F)) else core.rows)$core

  ps <- permCoreSummary(meta, blocks, core.rows, ctr, groups.core,
               block.factors = block.factors,
               show.top = thresholds$show.top)
  # Effective permutations (consider only blocks with size >= 2)
  movable <- ps$summary$n.core[ps$summary$n.core >= thresholds$min.core.size]
  eff.perm.log <- if (length(movable)) sum(lfactorial(movable)) else 0
  
  # Warn only when there is a problem
  too.small <- sum(ps$summary$n.core < thresholds$min.core.size)
  no.var    <- sum(!is.na(ps$summary$n.levels) & ps$summary$n.levels < 2L)
  
  if (too.small > 0 || no.var > 0 || exp(min(eff.perm.log, 50)) < thresholds$min.eff.perm) {
    if (too.small > 0)
    msgs <- c(msgs, sprintf("[WARN] %d block(s) have < %d core rows; those rows will be frozen (no permutation).",
                too.small, thresholds$min.core.size))
    if (no.var > 0)
    msgs <- c(msgs, sprintf("[WARN] %d block(s) show no within-block variation in the contrasted variable; those rows will be frozen.",
                no.var))
    if (exp(min(eff.perm.log, 50)) < thresholds$min.eff.perm)
    msgs <- c(msgs, sprintf("[WARN] Effective number of permutations is very small (~exp(%.1f)). Consider relaxing blocks (drop/merge a blocking factor) or a wild bootstrap.",
                eff.perm.log))
    
    # List problematic blocks with combinations if available
    if (length(ps$prob.idx)) {
    msgs <- c(msgs, " Problematic blocks (examples):")
    if (!is.null(ps$comb.df)) {
      apply(ps$comb.df, 1, function(row) {
      msgs <<- c(msgs, sprintf("       %s: n_core=%s; levels=%s%s; combo: %s",
                   row[["block"]], row[["n.core"]], row[["n.levels"]],
                   ifelse(nchar(row[["composition"]])>0,
                      paste0(" [", row[["composition"]], "]"), ""),
                   row[["combo"]]))
      })
    } else {
      show <- ps$summary[ps$prob.idx, , drop = FALSE]
      ord  <- order(show$n.core, show$n.levels)
      show <- head(show[ord, , drop = FALSE], thresholds$show.top)
      apply(show, 1, function(row) {
      msgs <<- c(msgs, sprintf("       %s: n_core=%s; levels=%s%s",
                   row[["block"]], row[["n.core"]], row[["n.levels"]],
                   ifelse(nchar(row[["composition"]])>0,
                      paste0(" [", row[["composition"]], "]"), "")))
      })
    }
    }
    
    # Aggregate which factor levels dominate problematic blocks (to guide dropping)
    if (!is.null(block.factors) && length(block.factors) &&
      all(block.factors %in% names(meta)) && length(ps$prob.idx)) {
    msgs <- c(msgs, " Factor levels most often appearing in problematic blocks:")
    prob.blocks <- ps$summary$block[ps$prob.idx]
    # For each block, pick a representative row to read factor levels
    rep.row <- vapply(prob.blocks, function(b) which(blocks == b)[1], integer(1))
    for (v in block.factors) {
      levs <- as.character(meta[[v]][rep.row])
      tab  <- sort(table(levs), decreasing = TRUE)
      topK <- head(tab, 5)
      msgs <- c(msgs, paste0("       ", v, ": ",
                 paste(sprintf("%s (%d)", names(topK), as.integer(topK)), collapse = ", ")))
    }
    }
  }
  
  perm <- list(groups.core = groups.core,
         core.summary = ps$summary,
         eff.perm.log = eff.perm.log)
  } else {
  # no blocks → no permutation diagnostics
  }

  emitInternal(msgs)

  list(
  full.rank    = dF,
  core.after.Z = dXr,
  aliased.by.Z = aXZ$aliased,
  vif          = vif,
  perm         = perm,
  messages     = msgs
  )
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



###################### Paired setting  #########################

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
# @param block.vars optional character vector to override permutation blocks
# @param na.action NA handler for model.matrix
# @return list with:
#   - pairs: data.frame of (i,j) indices (1-based) in lower-tri order
#   - pair.meta: data.frame of per-pair covariates
#   - prep: design list from buildDesignMatrices() (F, X, Z, qrZ, contrast.F/.X, ...)
#' @keywords internal
buildPairDesignMatrices <- function(sample.meta, triplet,
                                dist.type = c("shift","total","var"),
                                block.vars = NULL,
                                na.action = stats::na.pass) {
  stopifnot(is.data.frame(sample.meta))
  dist.type <- match.arg(dist.type)
  
  # 1) Indices for lower triangle, derived from nrow(meta)
  n   <- nrow(sample.meta)
  idx <- lowerTriIndices(n)

  # 2) Pair-level metadata from sample-level meta
  pair.meta <- pairifyMeta(sample.meta, idx)

  # 3) Triplet (var, alt, ref) → pair-level contrast weights over AA/AB/BB
  tr <- parseTriplet(triplet)                # var, alt, ref
  if (!(tr$var %in% names(sample.meta)))
    stop("Triplet 'var' not found in meta: ", tr$var)

  f.levels <- levels(factor(sample.meta[[tr$var]]))
  if (!all(c(tr$ref, tr$alt) %in% f.levels))
    stop("alt/ref levels not found in meta[['", tr$var, "']]: ",
         tr$alt, ", ", tr$ref)

  weights <- pairContrastWeights(tr$ref, tr$alt, dist.type, level.order = f.levels)

  # 4) Build pair-level variable name for the contrasted factor (e.g., "group_pair")
  pair.var <- paste0(tr$var, "_pair")
  if (!(pair.var %in% names(pair.meta)))
    stop("Pair variable not found in pair_meta: ", pair.var,
         " (did '", tr$var, "' exist and was it factor-like?)")

  # 5) Call buildDesignMatrices() on pair-level metadata
  prep <- buildDesignMatrices(
    sample.meta       = pair.meta,
    contrast   = list(var = pair.var, weights = weights), 
    #nuisance   = setdiff(names(pair.meta), pair.var),  # core is implicitly pair.var
    core.extra = NULL,
    block.vars = block.vars,
    na.action  = na.action
  )

  list(pairs = idx, pair.meta = pair.meta, model = prep, model.diag = prep$model.diag)
}

#' Return lower-triangular (row, col) index pairs for an n x n matrix
#' @param n integer number of samples
#' @return data.frame with columns i, j (1-based), in row-major lower-tri order
#' @keywords internal
lowerTriIndices <- function(n) {
  ij <- which(lower.tri(matrix(NA_real_, n, n)), arr.ind = TRUE)
  data.frame(i = ij[,1], j = ij[,2])
}

#' Build pair-level metadata from sample-level metadata in lower-tri order
#' @param meta data.frame of per-sample covariates (rows correspond to d's order)
#' @param idx data.frame with columns i, j (lower-tri indices)
#' @return data.frame of per-pair covariates
#' @keywords internal
pairifyMeta <- function(meta, idx) {
  stopifnot(is.data.frame(meta), all(c("i","j") %in% names(idx)))
  n <- nrow(meta)
  stopifnot(all(idx$i >= 1 & idx$i <= n & idx$j >= 1 & idx$j <= n & idx$i > idx$j))
  
  out <- list()
  for (v in names(meta)) {
    x <- meta[[v]]
    if (is.factor(x) || is.character(x) || is.logical(x)) {
      f <- factor(x)
      out[[paste0(v, "_pair")]] <- factor(pairCodeFactor(f, idx$i, idx$j))
    } else if (is.numeric(x)) {
      out[[paste0(v, "_diff")]] <- pairDiffNumeric(x, idx$i, idx$j)
    } else {
      # unsupported type → coerce to factor then pair-code
      f <- factor(as.character(x))
      out[[paste0(v, "_pair")]] <- factor(pairCodeFactor(f, idx$i, idx$j))
    }
  }
  as.data.frame(out, stringsAsFactors = TRUE)
}

#' Symmetric pair-code for a factor: "AA", "AB", "BB" (AB ≡ BA)
#' Uses the factor's level order to define A<B for concatenation.
#' @param f factor vector (length n)
#' @param i integer indices for first sample in each pair
#' @param j integer indices for second sample in each pair
#' @return character vector of codes ("AA","AB",...)
#' @keywords internal
pairCodeFactor <- function(f, i, j) {
  lev <- levels(f); fi <- as.integer(f[i]); fj <- as.integer(f[j])
  check.na <- is.na(fi) | is.na(fj)
  lo <- pmin(fi, fj); hi <- pmax(fi, fj)
  code <- paste0(lev[lo], lev[hi])
  code[check.na] <- NA_character_
  code
}

#' Symmetric numeric pair difference |x_i - x_j|
#' @param x numeric vector (length n)
#' @param i, j integer index vectors
#' @return numeric vector of absolute differences
#' @keywords internal
pairDiffNumeric <- function(x, i, j) {
  xi <- x[i]; xj <- x[j]
  out <- abs(xi - xj)
  out
}

#' Parse a sample-level "triplet" contrast specification
#' @param triplet list(var="<factor>", ref="<A>", alt="<B>")
#' @return list(var, ref, alt)
#' @keywords internal
parseTriplet <- function(triplet) {
  stopifnot(length(triplet)==3)
  triplet <- as.list(triplet)
  names(triplet) <- c('var','alt','ref')
  if (triplet$ref == triplet$alt) stop("ref and alt must be different.")
  list(var = as.character(triplet$var),
       alt = as.character(triplet$alt),
       ref = as.character(triplet$ref))
}

# Translate (ref, alt) + dist.type into pair-level contrast weights
# - level.order: optional character vector giving the factor's level order
#   (used to decide the label for the cross pair, e.g., "AB" vs "BA")
# - Returns a named numeric vector, e.g. c("AA"=-0.5, "AB"=1, "BB"=-0.5)
#' @keywords internal
pairContrastWeights <- function(ref, alt, type = c("shift","total","var"), level.order = NULL) {
  type <- match.arg(type)
  
  # Cross label uses factor-level order if provided; otherwise alphabetical
  two <- c(ref, alt)
  if (!is.null(level.order)) {
    ord <- two[order(match(two, level.order))]
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

## Generic grouping for pair-level design matrices
## - collapses all "Var_pair*" columns into one group per Var
## - maps "Var_pair_contrast" to Var
## - keeps "(Intercept)" if present
## - anything else stays as-is (unless include.other=TRUE, then grouped under "Other")
#' @keywords internal
makeGroupsPair <- function(M, include.other = FALSE) {
  stopifnot(!is.null(colnames(M)))
  cn <- colnames(M)
  groups <- list()
  if ("(Intercept)" %in% cn) {
    groups$Intercept <- "(Intercept)"
  }
  # Handle all _pair terms
  pair.terms <- grep("_pair", cn, value = TRUE)
  if (length(pair.terms)) {
    # Extract the variable name before "_pair"
    varnames <- sub("_pair.*$", "", pair.terms)
    for (v in unique(varnames)) {
      cols <- pair.terms[varnames == v]
      # if there's a collapsed contrast, rename group nicely
      if (any(grepl("_pair_contrast$", cols))) {
        groups[[v]] <- cols[grepl("_pair_contrast$", cols)]
      } else {
        groups[[v]] <- cols
      }
    }
  }

  # Leftovers (e.g. numeric pair variables not matching _pair)
  used <- unlist(groups, use.names = FALSE)
  leftovers <- setdiff(cn, used)
  if (length(leftovers)) {
    if (include.other) {
      groups$Other <- leftovers
    } else {
      for (lf in leftovers) groups[[lf]] <- lf
    }
  }
  groups
}

vectorizeLowerTri <- function(d, pairs = NULL) {
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
