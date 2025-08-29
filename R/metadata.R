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
    formula_str <- paste(deparse(formula), collapse = "")
  } else if (is.character(formula)) {
    formula_str <- formula
  } else {
    stop("Design formula must be a string or formula object.")
  }
  if (!grepl("~", formula_str)) {
    stop("Design formula must contain a '~' to separate response and predictors.")
  }

  # Demote random effects to fixed, preserving variables
  containsRandomEffects <- grepl("\\([^\\|]*\\|[^\\)]*\\)", formula_str)
  if (containsRandomEffects) {
    rand_eff_vars <- unlist(regmatches(formula_str, gregexpr("(?<=\\|)[^\\)]+", formula_str, perl = TRUE)))
    rand_eff_vars <- trimws(rand_eff_vars)
    warning(sprintf(
      "Random effects (terms with '|') are not supported in this workflow. The variable(s) '%s' will be treated as fixed effects.",
      paste(rand_eff_vars, collapse = ", ")
    ))
    formula_str <- gsub("\\([^\\|]*\\|[^\\)]*\\)", "", formula_str)
    rhs <- gsub("~", "", formula_str)
    rhs <- gsub("\\++", "+", rhs)
    rhs <- gsub("^\\s*\\+|\\+\\s*$", "", rhs)
    rhs <- trimws(rhs)
    rhs_terms <- trimws(unlist(strsplit(rhs, "\\+")))
    rhs_terms <- unique(c(rhs_terms, rand_eff_vars))
    rhs_terms <- rhs_terms[rhs_terms != ""]
    formula_str <- paste("~", paste(rhs_terms, collapse = " + "))
  }
  parsedTerms <- terms(stats::as.formula(formula_str))
  termLabels  <- attr(parsedTerms, "term.labels")

  # default contrast if none specified (uses first term; last vs first level)
  if (is.null(contrast)) {
    if (length(termLabels) == 0) stop("No terms found in the design formula to use as contrast.")
    if (is.null(sample.meta)) stop("sample.meta must be given to infer a default contrast.")
    contrast_var <- termLabels[1]
    if (!contrast_var %in% colnames(sample.meta)) {
      stop(sprintf("Contrast variable '%s' not found in sample metadata.", contrast_var))
    }
    levels_contrast <- levels(factor(sample.meta[[contrast_var]]))
    if (length(levels_contrast) < 2) {
      stop(sprintf("Contrast variable '%s' must have at least two levels.", contrast_var))
    }
    contrast <- c(contrast_var, levels_contrast[1], levels_contrast[length(levels_contrast)])
  }
  if (length(contrast) != 3) {
    stop("Contrast must be a vector of length 3: c('variable', 'level1', 'level2').")
  }
  if (verbose) {
    if (containsRandomEffects) message("Random effect terms provided will be treated as fixed effect terms.")
    message(sprintf("Final design formula: %s", formula_str))
  }
  list("formula" = parsedTerms, "contrast" = contrast)
}

#' @keywords internal
buildModelMatrix <- function(sample.meta, formula, contrast = NULL, keep.intercept = FALSE, verbose = FALSE) {
  if (!is.null(formula) || !is.null(contrast)) {
    vd <- validateDesign(formula = formula, sample.meta = sample.meta, contrast = contrast, verbose = verbose)
    formula <- vd$formula
    contrast <- vd$contrast
  }
  predictorVars <- all.vars(formula)
  sample.meta   <- sample.meta[, predictorVars, drop = FALSE]

  valid.vars <- vapply(predictorVars, function(cov) { # remove single-level covariates
    x <- sample.meta[[cov]]
    if (is.factor(x) || is.character(x)) {
      length(unique(x)) > 1
    } else {
      vx <- stats::var(x, na.rm = TRUE)
      !is.na(vx) && vx > 0
    }
  }, logical(1))

  if (!all(valid.vars)) {
    removed.vars <- predictorVars[!valid.vars]
    warning(sprintf("Removing covariates with only one factor level or zero variance: %s", paste(removed.vars, collapse = ", ")))
    predictorVars   <- predictorVars[valid.vars]
    sample.meta     <- sample.meta[, predictorVars, drop = FALSE]
    intercept.flag  <- attr(terms(formula), "intercept") == 1
    formula         <- reformulate(predictorVars, intercept = intercept.flag)
  }

  if (length(predictorVars) > 1) { # flag identical columns, warn if correlated
    to.remove <- character(0)
    for (i in seq_along(predictorVars)) for (j in seq_along(predictorVars)) if (i < j) {
      col1 <- sample.meta[[predictorVars[i]]]
      col2 <- sample.meta[[predictorVars[j]]]
      if (is.factor(col1)) col1 <- as.character(col1)
      if (is.factor(col2)) col2 <- as.character(col2)
      if (all(col1 == col2, na.rm = TRUE)) {
        warning(sprintf("Covariates '%s' and '%s' are identical. Removing '%s'.", predictorVars[i], predictorVars[j], predictorVars[j]))
        to.remove <- c(to.remove, predictorVars[j])
      } else if (is.numeric(col1) && is.numeric(col2)) {
        cor.val <- suppressWarnings(stats::cor(col1, col2, use = "pairwise.complete.obs"))
        if (!is.na(cor.val) && abs(cor.val) > 0.95) {
          warning(sprintf("Covariates '%s' and '%s' are highly correlated (cor = %.2f).", predictorVars[i], predictorVars[j], cor.val))
        }
      }
    }
    if (length(to.remove)) {
      to.remove      <- unique(to.remove)
      sample.meta    <- sample.meta[, !(colnames(sample.meta) %in% to.remove), drop = FALSE]
      predictorVars  <- colnames(sample.meta)
      intercept.flag <- attr(terms(formula), "intercept") == 1
      formula        <- reformulate(predictorVars, intercept = intercept.flag)
    }
  }
  
  for (cov in predictorVars) { # coerce non-numeric/non-factor to factor
    if (!is.numeric(sample.meta[[cov]]) && !is.factor(sample.meta[[cov]])) {
      if (verbose) message(sprintf("Converting covariate '%s' to factor.", cov))
      lvls <- sort(unique(stats::na.omit(sample.meta[[cov]])))
      sample.meta[[cov]] <- factor(sample.meta[[cov]], levels = lvls)
    }
  }

  if (!is.null(contrast)) { # apply requested contrast’s reference/target level on its factor
    contrast.var <- contrast[1]; ref.level <- contrast[2]; target.level <- contrast[3]
    if (contrast.var %in% colnames(sample.meta)) {
      unique.levels <- unique(sample.meta[[contrast.var]])
      if (!ref.level %in% unique.levels || !target.level %in% unique.levels) {
        stop("ref.level or target.level not found in levels of ", contrast.var)
      }
      sample.meta[[contrast.var]] <- factor(sample.meta[[contrast.var]], levels = c(ref.level, setdiff(unique.levels, ref.level)))
    }
  }

  # Build the model matrix (keep intercept then optionally drop)
  mm <- stats::model.matrix(formula, data = sample.meta)
  assign.vec  <- attr(mm, "assign")
  contrasts.a <- attr(mm, "contrasts")

  if (!keep.intercept && "(Intercept)" %in% colnames(mm)) {
    keep <- colnames(mm) != "(Intercept)"
    mm   <- mm[, keep, drop = FALSE]
    if (!is.null(assign.vec))  attr(mm, "assign")    <- assign.vec[keep]
    if (!is.null(contrasts.a)) attr(mm, "contrasts") <- contrasts.a
  } else {
    if (!is.null(assign.vec))  attr(mm, "assign")    <- assign.vec
    if (!is.null(contrasts.a)) attr(mm, "contrasts") <- contrasts.a
  }

  qr.decomp <- qr(mm) # rank check
  if (qr.decomp$rank < ncol(mm)) {
    warning(sprintf(
      "Model matrix has linear dependencies: rank %d < number of columns %d. Possible confounding or redundant covariates.",
      qr.decomp$rank, ncol(mm)
    ))
  }
  dimnames(mm) <- list(rownames(mm), colnames(mm))
  mm
}

#' @keywords internal
getSampleGroups <- function(sample.meta, contrast, sample.id = NULL) {
  if (is.null(contrast)) return(NULL)
  var <- contrast[1]; ref <- contrast[2]; alt <- contrast[3]
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