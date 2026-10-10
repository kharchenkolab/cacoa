## Design check without expression data (Track B.4): associations among covariates, balance against the
## test variable, structural issues (aliasing, singleton cells, rank), collinearity (GVIF), degrees-of-freedom
## budget, permutation feasibility and over-adjustment candidates.

# ---- pairwise association measures ---------------------------------------------------------------

# bias-corrected Cramer's V (Bergsma 2013) between two discrete vectors
cramersV <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- factor(x[ok]); y <- factor(y[ok]); n <- length(x)
  if (n < 2 || nlevels(x) < 2 || nlevels(y) < 2) return(NA_real_)
  tb <- table(x, y); chi <- suppressWarnings(stats::chisq.test(tb, correct = FALSE)$statistic)
  phi2 <- unname(chi) / n; k <- ncol(tb); r <- nrow(tb)
  phi2c <- max(0, phi2 - (k - 1) * (r - 1) / (n - 1))
  kc <- k - (k - 1)^2 / (n - 1); rc <- r - (r - 1)^2 / (n - 1)
  d <- min(kc - 1, rc - 1)
  if (d <= 0) return(NA_real_)
  sqrt(phi2c / d)
}

# correlation ratio eta (factor x numeric): sqrt(between-group SS / total SS)
correlationRatio <- function(f, y) {
  ok <- !is.na(f) & !is.na(y); f <- factor(f[ok]); y <- y[ok]
  if (length(y) < 3 || nlevels(f) < 2 || stats::var(y) == 0) return(NA_real_)
  m <- tapply(y, f, mean); nk <- table(f)
  sqrt(sum(nk * (m - mean(y))^2) / sum((y - mean(y))^2))
}

isDiscreteVar <- function(x) is.factor(x) || is.character(x) || is.logical(x)

#' Association between two covariates on a common 0-1 scale
#'
#' Bias-corrected Cramer's V for two discrete variables, the correlation ratio eta for a discrete and a
#' numeric variable, and |Spearman rho| for two numeric variables.
#' @param x,y vectors of equal length
#' @return a number between 0 and 1 (NA when undefined)
#' @export
covariateAssociation <- function(x, y) {
  dx <- isDiscreteVar(x); dy <- isDiscreteVar(y)
  if (dx && dy) cramersV(x, y)
  else if (dx) correlationRatio(x, y)
  else if (dy) correlationRatio(y, x)
  else { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 3) NA_real_ else abs(suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))) }
}

associationMatrix <- function(meta, covariates) {
  k <- length(covariates); A <- matrix(NA_real_, k, k, dimnames = list(covariates, covariates))
  for (i in seq_len(k)) for (j in seq_len(k)) {
    if (i == j) A[i, j] <- 1 else if (j > i) A[i, j] <- A[j, i] <- covariateAssociation(meta[[covariates[i]]], meta[[covariates[j]]])
  }
  A
}

# generalized VIF per term of a formula (Fox & Monette 1992), reported as GVIF^(1/(2 df))
gvifPerTerm <- function(formula, meta) {
  F <- tryCatch(buildFullDesign(formula, meta), error = function(e) NULL)
  if (is.null(F)) return(numeric(0))
  assign <- attr(F, "assign"); tl <- attr(attr(F, "terms"), "term.labels")
  keep <- assign > 0 & apply(F, 2, stats::sd) > 0
  X <- F[, keep, drop = FALSE]; a <- assign[keep]
  if (ncol(X) < 2 || length(unique(a)) < 2) return(numeric(0))
  R <- stats::cor(X); detR <- det(R)
  out <- setNames(rep(NA_real_, length(unique(a))), tl[unique(a)])
  for (k in unique(a)) {
    A <- which(a == k); B <- which(a != k)
    g <- det(R[A, A, drop = FALSE]) * det(R[B, B, drop = FALSE]) / detR
    out[tl[k]] <- if (is.finite(g) && g > 0) g^(1 / (2 * length(A))) else NA_real_
  }
  out
}

# ---- the check -----------------------------------------------------------------------------------

#' Check a study design from the sample metadata alone
#'
#' @param meta sample metadata (rows = samples)
#' @param model optional `cacoaModel`; its tests, formula and permutation plans are checked
#' @param test.variable test variable name (default: from `model`)
#' @param covariates covariates to examine (default: all usable columns, see [describeMetadata()])
#' @param block.vars permutation strata variables
#' @param assoc.flag association level above which two covariates are reported as related (default 0.3)
#' @param overadjust.flag association with the test variable above which a covariate is an
#'   over-adjustment candidate (default 0.5)
#' @return object of class `cacoaDesignCheck`: `associations` (matrix), `balance` (per covariate), `issues`
#'   (data.frame: severity, message, suggestion), `df` (n, parameters, residual), `gvif`, `permutation` (per
#'   test), `test.variable`, `covariates`
#' @export
checkDesign <- function(meta, model = NULL, test.variable = NULL, covariates = NULL, block.vars = NULL,
                        assoc.flag = 0.3, overadjust.flag = 0.5) {
  stopifnot(is.data.frame(meta))
  desc <- describeMetadata(meta)
  if (is.null(test.variable) && !is.null(model)) test.variable <- model$tests[[1]]$variable
  if (!is.null(test.variable) && (is.na(test.variable) || !test.variable %in% names(meta))) test.variable <- NULL
  if (is.null(covariates)) covariates <- desc$column[desc$role == "usable"]
  covariates <- unique(c(test.variable, if (!is.null(model)) intersect(all.vars(model$formula), names(meta)), covariates, block.vars))
  covariates <- intersect(covariates, names(meta))
  iss <- list()
  add <- function(sev, msg, sug = "") iss[[length(iss) + 1L]] <<- data.frame(severity = sev, message = msg, suggestion = sug, stringsAsFactors = FALSE)
  n <- nrow(meta)

  A <- associationMatrix(meta, covariates)
  # related covariates
  if (length(covariates) > 1) {
    pairs <- which(upper.tri(A) & !is.na(A) & A > assoc.flag, arr.ind = TRUE)
    for (k in seq_len(nrow(pairs))) {
      i <- pairs[k, 1]; j <- pairs[k, 2]
      if (!is.null(test.variable) && test.variable %in% covariates[c(i, j)]) next   # reported below
      add("note", sprintf("'%s' and '%s' are associated (%.2f)", covariates[i], covariates[j], A[i, j]),
          "adjusting for both estimates each one's unique contribution only")
    }
  }

  balance <- list(); perm <- NULL
  if (!is.null(test.variable)) {
    tv <- meta[[test.variable]]
    for (v in setdiff(covariates, test.variable)) {
      x <- meta[[v]]
      if (isDiscreteVar(x)) {
        tb <- table(test = tv, covariate = x, useNA = "no")
        balance[[v]] <- as.data.frame.matrix(tb)
        if (isDiscreteVar(tv)) {
          if (all(colSums(tb > 0) <= 1) && ncol(tb) >= 2)
            add("error", sprintf("'%s' is fully determined by '%s' (each %s contains one %s): the %s effect cannot be separated from %s",
                                 test.variable, v, v, test.variable, test.variable, v), sprintf("do not adjust for '%s', or test it instead", v))
          else if (all(rowSums(tb > 0) <= 1) && nrow(tb) >= 2)
            add("warning", sprintf("'%s' is nested within '%s' (each %s has one %s)", v, test.variable, test.variable, v),
                sprintf("'%s' can only be a blocking factor (block.vars), not an adjustment", v))
          else {
            if (any(tb == 1)) add("note", sprintf("%d cell(s) of %s x %s hold a single sample", sum(tb == 1), test.variable, v), "estimates in those cells rest on one sample")
            if (any(tb == 0)) add("note", sprintf("%d empty cell(s) in %s x %s", sum(tb == 0), test.variable, v), "the design is unbalanced; a marginal contrast cannot average over empty cells")
          }
        }
      } else if (isDiscreteVar(tv)) {
        balance[[v]] <- do.call(rbind, lapply(split(x, tv), function(z) data.frame(n = sum(!is.na(z)), mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE))))
      }
      a <- A[test.variable, v]
      if (is.finite(a) && a > overadjust.flag)
        add("warning", sprintf("'%s' is strongly associated with '%s' (%.2f): adjusting for it estimates the %s effect not explained by %s",
                               v, test.variable, a, test.variable, v), "consider reporting both models (cao$checkSensitivity())")
      else if (is.finite(a) && a > assoc.flag)
        add("note", sprintf("'%s' differs between the levels of '%s' (association %.2f)", v, test.variable, a),
            sprintf("adjust for '%s' and compare with the unadjusted result (cao$checkSensitivity())", v))
    }
    if (isDiscreteVar(tv)) {
      tb <- table(tv); if (any(tb < 3)) add("warning", sprintf("'%s' has only %d sample(s) at level %s", test.variable, min(tb), names(tb)[which.min(tb)]), "cell types need min.samp.per.level samples per level")
    }
  }

  if (!is.null(test.variable) && isDiscreteVar(meta[[test.variable]])) {   # excluded factor columns nested in the test variable
    excl <- desc$column[desc$role %in% c("high-cardinality", "id-like") & !desc$column %in% covariates]
    for (v in excl) {
      x <- meta[[v]]; if (!isDiscreteVar(x)) next
      tb <- table(meta[[test.variable]], x, useNA = "no")
      if (ncol(tb) >= 2 && ncol(tb) < sum(!is.na(x)) && all(colSums(tb > 0) <= 1) && any(rowSums(tb > 0) > 1))
        add("warning", sprintf("'%s' (%s, not used as a covariate) is nested in '%s': each %s holds one %s, so a %s effect cannot be separated from the %s effect",
                               v, desc$role[desc$column == v], test.variable, v, test.variable, v, test.variable),
            sprintf("report this limitation; '%s' can serve as a permutation block only if it crosses '%s'", v, test.variable))
      else if (ncol(tb) >= 2 && all(tb <= 1) && mean(colSums(tb > 0) >= 2) >= 0.5)
        add("note", sprintf("'%s' pairs the samples across '%s' (%d of %d %ss have a sample at more than one level): a blocking factor",
                            v, test.variable, sum(colSums(tb > 0) >= 2), ncol(tb), v),
            sprintf("use block.vars = \"%s\" (permutations within %s), or add it to the model if the degrees of freedom allow", v, v))
    }
  }

  dfinfo <- NULL; gvif <- numeric(0)
  if (!is.null(model)) {
    F <- model$F; q <- qr(F)$rank
    dfinfo <- list(n = nrow(F), parameters = q, residual = nrow(F) - q)
    if (q < ncol(F)) add("warning", sprintf("the design is rank-deficient (rank %d < %d columns)", q, ncol(F)), "remove redundant covariates")
    if (dfinfo$residual < 10) add(if (dfinfo$residual < 1) "error" else "warning", sprintf("only %d residual degrees of freedom (n = %d, %d parameters)", dfinfo$residual, nrow(F), q),
                                  "fewer covariates")
    gvif <- gvifPerTerm(model$formula, model$meta)
    bad <- gvif[is.finite(gvif) & gvif > 2.5]
    for (nm in names(bad)) add("warning", sprintf("term '%s' is collinear with the rest of the model (GVIF^(1/2df) = %.1f)", nm, bad[[nm]]), "drop or merge related covariates")
    perm <- do.call(rbind, lapply(model$tests, function(t) if (is.null(t$permutation)) NULL else
      data.frame(test = t$label, scheme = t$permutation$scheme, n.strata = t$permutation$n.strata, n.distinct = t$permutation$n.distinct,
                 p.floor = t$permutation$p.floor, stringsAsFactors = FALSE)))
    for (t in model$tests) if (!is.null(t$permutation)) {
      if (is.finite(t$permutation$n.distinct) && t$permutation$n.distinct < 200)
        add("warning", sprintf("test '%s': only %d distinct permutations (smallest p-value %.3g)", t$label, round(t$permutation$n.distinct), t$permutation$p.floor), "fewer strata, or more samples")
      if (t$permutation$scheme == "freedman-lane") add("note", sprintf("test '%s': continuous covariate(s) -> Freedman-Lane residual permutation (approximate; conservative for shift)", t$label),
                                                      "exact block permutations need discrete covariates only")
    }
    dropped <- model$samples$dropped
    if (length(dropped)) add("warning", sprintf("%d sample(s) dropped for missing covariate values: %s", length(dropped), paste(dropped, collapse = ", ")), "complete the metadata or remove the covariate")
  } else {
    dfinfo <- list(n = n, parameters = NA_integer_, residual = NA_integer_)
  }
  excl <- desc[desc$role != "usable" & desc$column %in% covariates, ]
  for (i in seq_len(nrow(excl))) {
    v <- excl$column[i]; x <- meta[[v]]; pairing <- FALSE
    if (!is.null(test.variable) && isDiscreteVar(x) && isDiscreteVar(meta[[test.variable]])) {   # a pairing factor is a block, not a nuisance
      tb <- table(meta[[test.variable]], x, useNA = "no")
      pairing <- ncol(tb) >= 2 && all(tb <= 1) && mean(colSums(tb > 0) >= 2) >= 0.5
    }
    if (pairing) add("note", sprintf("'%s' pairs the samples across '%s' (%d of %d %ss have a sample at more than one level): a blocking factor",
                                     v, test.variable, sum(colSums(table(meta[[test.variable]], x) > 0) >= 2), nlevels(factor(x)), v),
                     sprintf("use block.vars = \"%s\" (permutations within %s), or add it to the model if the degrees of freedom allow", v, v))
    else add("note", sprintf("'%s' is %s", v, excl$role[i]), "not suitable as a covariate")
  }

  issues <- if (length(iss)) do.call(rbind, iss) else data.frame(severity = character(0), message = character(0), suggestion = character(0), stringsAsFactors = FALSE)
  issues$severity <- factor(issues$severity, levels = c("error", "warning", "note"))
  issues <- issues[order(issues$severity), ]; rownames(issues) <- NULL
  structure(list(associations = A, balance = balance, issues = issues, df = dfinfo, gvif = gvif, permutation = perm,
                 test.variable = test.variable, covariates = covariates, n = n, assoc.flag = assoc.flag),
            class = c("cacoaDesignCheck", "list"))
}

#' @exportS3Method base::print
print.cacoaDesignCheck <- function(x, ...) {
  cat(sprintf("Design check: %d samples, %d covariates%s\n", x$n, length(x$covariates),
              if (!is.null(x$test.variable)) sprintf(", test variable '%s'", x$test.variable) else ""))
  if (!is.null(x$df) && is.finite(x$df$parameters)) cat(sprintf("  model: %d parameters, %d residual degrees of freedom\n", x$df$parameters, x$df$residual))
  if (!is.null(x$permutation)) for (i in seq_len(nrow(x$permutation))) {
    p <- x$permutation[i, ]
    cat(sprintf("  permutations for '%s': %s%s%s\n", p$test, p$scheme, if (p$n.strata > 1) sprintf(" within %d strata", p$n.strata) else "",
                if (is.finite(p$n.distinct)) sprintf(", %s distinct (smallest p %.3g)", format(round(p$n.distinct), big.mark = ","), p$p.floor) else ""))
  }
  A <- x$associations
  if (!is.null(x$test.variable) && length(x$covariates) > 1) {
    cv <- setdiff(x$covariates, x$test.variable)
    a <- sort(stats::setNames(A[x$test.variable, cv], cv), decreasing = TRUE)
    a <- a[is.finite(a)]
    if (length(a)) cat(sprintf("  association with '%s': %s\n", x$test.variable, paste(sprintf("%s %.2f", names(a), a), collapse = ", ")))
  }
  if (!nrow(x$issues)) cat("  no issues found\n") else {
    for (sev in levels(x$issues$severity)) {
      rows <- x$issues[x$issues$severity == sev, ]
      if (!nrow(rows)) next
      cat(sprintf("  %s%s (%d):\n", toupper(substring(sev, 1, 1)), substring(sev, 2), nrow(rows)))
      for (i in seq_len(nrow(rows))) cat(sprintf("    - %s%s\n", rows$message[i], if (nzchar(rows$suggestion[i])) paste0(" -> ", rows$suggestion[i]) else ""))
    }
  }
  invisible(x)
}

# ---- plots ---------------------------------------------------------------------------------------

#' Plots of a design check
#'
#' @param x `cacoaDesignCheck` object
#' @param type `"associations"` (clustered heatmap of covariate associations), `"balance"` (covariates split by the
#'   test variable) or `"issues"` (table of issues)
#' @param meta sample metadata (needed for `"balance"`)
#' @param covariates optional subset
#' @param plot.theme ggplot2 theme
#' @param label.min label association cells above this value (default: the check's `assoc.flag`)
#' @return ggplot2 object (for `"balance"`: a cowplot grid)
#' @export
plotDesignCheck <- function(x, type = c("associations", "balance", "issues"), meta = NULL, covariates = NULL,
                            plot.theme = ggplot2::theme_bw(), label.min = NULL) {
  type <- match.arg(type)
  if (type == "associations") {
    A <- x$associations
    if (!is.null(covariates)) { cv <- intersect(covariates, rownames(A)); A <- A[cv, cv, drop = FALSE] }
    if (nrow(A) < 2) stop("at least two covariates are needed")
    ord <- if (nrow(A) > 2) { D <- 1 - A; D[is.na(D)] <- 1; rownames(A)[stats::hclust(stats::as.dist(D))$order] } else rownames(A)
    df <- data.frame(a = rep(rownames(A), ncol(A)), b = rep(colnames(A), each = nrow(A)), value = as.vector(A), stringsAsFactors = FALSE)
    df$a <- factor(df$a, levels = ord); df$b <- factor(df$b, levels = ord)
    df$label <- ifelse(!is.na(df$value) & df$value > (label.min %||% x$assoc.flag) & df$a != df$b, sprintf("%.2f", df$value), "")
    gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$a, y = .data$b, fill = .data$value)) + ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
      ggplot2::scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, 1), na.value = "grey90", name = "association") +
      plot.theme + ggplot2::labs(x = NULL, y = NULL, title = "Covariate associations",
                                 subtitle = "Cramer's V (two factors), correlation ratio (factor vs numeric), |Spearman| (two numerics)") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"))
    if (!is.null(x$test.variable) && x$test.variable %in% ord) {
      i <- match(x$test.variable, ord)
      gg <- gg + ggplot2::annotate("rect", xmin = i - 0.5, xmax = i + 0.5, ymin = 0.5, ymax = length(ord) + 0.5, fill = NA, colour = "black", linewidth = 0.8)
    }
    return(gg)
  }
  if (type == "issues") {
    iss <- x$issues
    if (!nrow(iss)) iss <- data.frame(severity = "note", message = "no issues found", suggestion = "")
    iss$y <- rev(seq_len(nrow(iss)))
    cols <- c(error = "#b2182b", warning = "#e08214", note = "#4393c3")
    gg <- ggplot2::ggplot(iss, ggplot2::aes(y = .data$y)) +
      ggplot2::geom_label(ggplot2::aes(x = 0, label = .data$severity, fill = .data$severity), colour = "white", hjust = 0, size = 3, label.size = 0) +
      ggplot2::geom_text(ggplot2::aes(x = 1, label = wrap_strings(paste0(.data$message, ifelse(nzchar(.data$suggestion), paste0("\n-> ", .data$suggestion), "")), 110)), hjust = 0, size = 3, lineheight = 0.9) +
      ggplot2::scale_fill_manual(values = cols, guide = "none") + ggplot2::xlim(0, 12) + ggplot2::theme_void() + ggplot2::labs(title = "Design issues")
    return(gg)
  }
  # balance
  if (is.null(x$test.variable)) stop("balance plots need a test variable")
  if (is.null(meta)) stop("`meta` is needed for balance plots")
  cv <- setdiff(covariates %||% x$covariates, x$test.variable)
  tv <- meta[[x$test.variable]]
  plots <- lapply(cv, function(v) {
    df <- data.frame(test = tv, value = meta[[v]])
    df <- df[!is.na(df$test) & !is.na(df$value), ]
    if (isDiscreteVar(df$value)) {
      df$value <- factor(df$value)
      ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$test), fill = .data$value)) + ggplot2::geom_bar(position = "fill", colour = "grey30") +
        plot.theme + ggplot2::labs(x = x$test.variable, y = "fraction of samples", fill = v, title = v)
    } else {
      ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$test), y = .data$value)) + ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey90") +
        ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) + plot.theme + ggplot2::labs(x = x$test.variable, y = v, title = v)
    }
  })
  if (!length(plots)) stop("no covariates to plot")
  cowplot::plot_grid(plotlist = plots, ncol = min(3, length(plots)))
}
