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