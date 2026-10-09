## Sensitivity of expression-shift results (Track D.1): alternative covariate sets (D31) and single-sample
## influence, summarized as a plain verdict per cell type.

# Candidate model formulas for a test: unadjusted, current, minus each adjusting covariate, plus each of the
# top-k screened covariates not in the model, all screened (if the df budget allows).
sensitivityFormulas <- function(model, screen = NULL, meta, top.k = 3, min.resid.df = 5) {
  t <- model$tests[[1]]; tv <- t$variable
  adj <- setdiff(all.vars(model$formula), tv)
  n <- length(model$samples$used)
  out <- list()
  out[["unadjusted"]] <- stats::as.formula(paste("~", tv))
  out[["current"]] <- model$formula
  for (v in adj) out[[paste("minus", v)]] <- stats::as.formula(paste("~", paste(c(tv, setdiff(adj, v)), collapse = " + ")))
  if (!is.null(screen)) {
    sm <- if ("partial" %in% screen$modes) "partial" else screen$modes[1]
    g <- screen$global[screen$global$mode == sm, ]
    g <- g[!g$covariate %in% c(tv, adj) & is.finite(g$p.global), ]
    g <- g[order(-g$n.sig, g$p.global), ]
    cand <- utils::head(g$covariate[g$n.sig > 0], top.k)
    for (v in cand) out[[paste("plus", v)]] <- stats::as.formula(paste("~", paste(c(tv, adj, v), collapse = " + ")))
    if (length(cand) > 1) {
      f.all <- stats::as.formula(paste("~", paste(c(tv, adj, cand), collapse = " + ")))
      q <- tryCatch(qr(buildFullDesign(f.all, meta))$rank, error = function(e) Inf)
      if (n - q >= min.resid.df) out[["all screened"]] <- f.all
    }
  }
  keep <- !duplicated(vapply(out, function(f) paste(sort(all.vars(f)), collapse = "+"), character(1)))
  out[keep]
}

#' Sensitivity of expression-shift effects to the covariate set and to single samples
#'
#' @param D.list named list of sample distance matrices (per cell type)
#' @param model the current `cacoaModel` (its first test is examined)
#' @param meta sample metadata
#' @param formulas named list of alternative location formulas (default: the D31 sets of `sensitivityFormulas()`)
#' @param screen optional `cacoaCovariateScreen` used to propose added covariates
#' @param influence optional per-test influence matrices from the current result (`res$influence[[1]]`) and the
#'   current wide table (`res$wide[[1]]`), as `list(influence =, wide =)`; when `NULL` they are recomputed
#' @param dist,permutation,n.permutations,seed,alpha,min.samp.per.level,n.cores engine settings
#' @param top.k number of screened covariates to try adding
#' @param rel.change relative change of the shift estimate across models that counts as "sensitive" (default 0.5)
#' @param influence.rel relative change of the shift estimate when one sample is left out that, together with being
#'   an outlier among the leave-one-out changes (> 3 MAD), marks the sample as driving the result (default 0.2)
#' @return object of class `cacoaSensitivity`: `table` (model x cell type x effect: estimate, ci, p, padj),
#'   `summary` (per cell type: same sign in all models, n significant of m, max relative change, most
#'   influential sample, verdict), `formulas`, `settings`
#' @export
#' @param robust,na.mode,robust.k robust fit (`"none"`, `"huber"`, `"winsor"`), treatment of absent samples (`"drop"`, `"impute_weak"`) and
#'   robust tuning constant, passed to [expressionShiftsForModel()] (default: the settings of the result being checked)
checkSensitivity <- function(D.list, model, meta, formulas = NULL, screen = NULL, influence = NULL, dist = "cor", permutation = "auto",
                             n.permutations = 199, seed = 1, alpha = 0.05, min.samp.per.level = 3, n.cores = 1, top.k = 3, rel.change = 0.5,
                             influence.rel = 0.2, robust = "none", na.mode = "drop", robust.k = 1.345) {
  t <- model$tests[[1]]
  if (t$kind != "contrast") stop("sensitivity analysis is defined for contrast tests (two groups or a numeric step)")
  if (is.null(formulas)) formulas <- sensitivityFormulas(model, screen, meta, top.k = top.k)
  spec <- model$test.spec %||% list(test = NULL, contrast = t$contrast)
  runs <- lapply(names(formulas), function(nm) {
    m <- tryCatch(buildCacoaModel(meta, formula = formulas[[nm]], test = spec$test, contrast = spec$contrast, dispersion.formula = model$dispersion.formula,
                                  block.vars = model$block.vars, numeric.ref = model$numeric.ref, numeric.step = model$numeric.step,
                                  permutation = permutation, n.permutations = n.permutations), error = function(e) NULL)
    if (is.null(m) || any(m$issues$severity == "error")) return(NULL)
    m$tests <- m$tests[1]
    r <- expressionShiftsForModel(D.list, m, dist = dist, permutation = permutation, n.permutations = n.permutations, block.vars = model$block.vars,
                                  influence = nm == "current" && is.null(influence), min.samp.per.level = min.samp.per.level, seed = seed, alpha = alpha,
                                  n.cores = n.cores, robust = robust, na.mode = na.mode, robust.k = robust.k)
    if (!nrow(r$results)) return(NULL)
    rows <- r$results; rows$model <- nm; rows$formula <- paste(deparse(formulas[[nm]]), collapse = "")
    list(rows = rows, res = r)
  })
  names(runs) <- names(formulas)
  runs <- Filter(Negate(is.null), runs)
  if (!length(runs)) stop("no alternative model could be fitted")
  tab <- do.call(rbind, lapply(runs, `[[`, "rows")); rownames(tab) <- NULL
  tab$model <- factor(tab$model, levels = names(formulas)[names(formulas) %in% unique(tab$model)])
  # influence: from the current result when supplied, else from the "current" run
  cur <- runs[["current"]]
  infl <- if (!is.null(influence)) influence else if (!is.null(cur)) list(influence = cur$res$influence[[1]], wide = cur$res$wide[[1]]) else NULL
  # summary per cell type on the shift effect
  sh <- tab[tab$effect == "shift", ]
  summ <- do.call(rbind, lapply(split(sh, sh$celltype), function(d) {
    cur.row <- d[d$model == "current", ]
    est0 <- if (nrow(cur.row)) cur.row$estimate else d$estimate[1]
    same.sign <- length(unique(sign(d$estimate[is.finite(d$estimate)]))) == 1
    n.sig <- sum(d$padj < alpha, na.rm = TRUE); m <- sum(is.finite(d$p))
    max.rel <- if (is.finite(est0) && est0 != 0) max(abs(d$estimate - est0) / abs(est0), na.rm = TRUE) else NA_real_
    worst <- if (nrow(d) > 1) { dd <- d[d$model != "current", ]; dd$model[which.max(abs(dd$estimate - est0))] } else NA
    sig.change <- if (nrow(cur.row)) any((d$padj < alpha) != (cur.row$padj < alpha), na.rm = TRUE) else NA
    flipped.by <- if (!is.na(sig.change) && sig.change) paste(d$model[(d$padj < alpha) != (cur.row$padj < alpha) & !is.na(d$padj)], collapse = ", ") else ""
    # influence
    infl.sample <- NA_character_; infl.rel <- NA_real_
    if (!is.null(infl) && !is.null(infl$influence) && d$celltype[1] %in% colnames(infl$influence$shift)) {
      ch <- infl$influence$shift[, d$celltype[1]]; ch <- ch[is.finite(ch)]
      if (length(ch) > 2) {
        i <- which.max(abs(ch)); infl.rel <- abs(ch[i]) / abs(est0)
        others <- abs(ch[-i]); outlying <- abs(ch[i]) > stats::median(others) + 3 * stats::mad(others)
        if (isTRUE(outlying) && is.finite(infl.rel) && infl.rel > influence.rel) infl.sample <- names(ch)[i]
      }
    }
    verdict <- if (!is.na(infl.sample)) sprintf("driven by sample %s", infl.sample)
      else if (!same.sign) "sign depends on adjustment"
      else if (!is.na(sig.change) && sig.change) sprintf("significance depends on adjustment (%s)", flipped.by)
      else if (is.finite(max.rel) && max.rel > rel.change) sprintf("estimate changes by %.0f%% (%s)", 100 * max.rel, as.character(worst))
      else "robust"
    data.frame(celltype = d$celltype[1], estimate.current = est0, same.sign = same.sign, n.significant = n.sig, n.models = m,
               max.rel.change = max.rel, most.changing.model = as.character(worst), influential.sample = infl.sample,
               influence.rel = infl.rel, verdict = verdict, stringsAsFactors = FALSE)
  }))
  rownames(summ) <- NULL
  structure(list(table = tab, summary = summ, formulas = formulas, test = t$label,
                 settings = list(dist = dist, permutation = permutation, n.permutations = n.permutations, alpha = alpha, rel.change = rel.change, influence.rel = influence.rel)),
            class = c("cacoaSensitivity", "list"))
}

#' @exportS3Method base::print
print.cacoaSensitivity <- function(x, ...) {
  cat(sprintf("Sensitivity of '%s' (shift) across %d models: %s\n", x$test, length(x$formulas),
              paste(sprintf("%s [%s]", names(x$formulas), vapply(x$formulas, function(f) paste(deparse(f), collapse = ""), character(1))), collapse = "; ")))
  s <- x$summary; s <- s[order(s$verdict != "robust", -abs(s$estimate.current)), ]
  for (i in seq_len(nrow(s))) cat(sprintf("  %-20s %s (significant in %d of %d models)\n", s$celltype[i], s$verdict[i], s$n.significant[i], s$n.models[i]))
  invisible(x)
}

#' Forest plot of a sensitivity analysis
#'
#' @param x `cacoaSensitivity`
#' @param effect effect to show (default `"shift"`)
#' @param cell.types optional subset
#' @param normalized plot normalized estimates (default TRUE)
#' @param plot.theme ggplot2 theme
#' @return ggplot2 object: cell types as panels, models on the y-axis, estimate with jackknife CI (when available)
#'   on the x-axis, filled points significant, the current model highlighted
#' @export
plotSensitivity <- function(x, effect = "shift", cell.types = NULL, normalized = TRUE, plot.theme = ggplot2::theme_bw()) {
  d <- x$table[x$table$effect == effect, ]
  if (!is.null(cell.types)) d <- d[d$celltype %in% cell.types, ]
  if (!nrow(d)) stop("nothing to plot")
  val <- if (normalized) "estimate.norm" else "estimate"
  d$value <- d[[val]]
  sc <- if (normalized) ifelse(is.finite(d$estimate) & d$estimate != 0, d$estimate.norm / d$estimate, NA) else 1
  d$lo <- d$value - 1.96 * d$se.jk * abs(sc); d$hi <- d$value + 1.96 * d$se.jk * abs(sc)
  d$sig <- !is.na(d$padj) & d$padj < x$settings$alpha
  d$current <- d$model == "current"
  d$model <- factor(d$model, levels = rev(levels(x$table$model)))
  ord <- x$summary$celltype[order(x$summary$verdict != "robust", -abs(x$summary$estimate.current))]
  d$celltype <- factor(d$celltype, levels = intersect(ord, unique(d$celltype)))
  verd <- setNames(x$summary$verdict, x$summary$celltype)
  d$panel <- factor(sprintf("%s\n%s", as.character(d$celltype), verd[as.character(d$celltype)]), levels = sprintf("%s\n%s", levels(d$celltype), verd[levels(d$celltype)]))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$model)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), width = 0.25, colour = "grey40", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$sig, colour = .data$current), size = 2.6, fill = "white", stroke = 0.9) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19), labels = c("not significant", "significant"), name = NULL) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "#d73027"), guide = "none") +
    ggplot2::facet_wrap(~ panel, scales = "free_x") + plot.theme +
    ggplot2::labs(x = if (normalized) sprintf("normalized %s", effect) else effect, y = NULL, title = sprintf("Sensitivity of %s to the covariate set", effect),
                  subtitle = "red: current model; intervals: jackknife 95%") +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"), strip.text = ggplot2::element_text(size = 7))
}
