## Plots for the long results table of expression-shift tests (Track B.6).

effectOrder <- c("shift", "var", "total", "location", "dispersion", "any")
effectTitle <- c(shift = "shift", var = "var (dispersion change)", total = "total", location = "location (adj. R2)",
                 dispersion = "dispersion (adj. R2)", any = "any difference")

#' Dot / bar plot of effects per cell type
#'
#' @param df long results table (`test`, `celltype`, `effect`, `estimate`, `estimate.norm`, `se.jk`, `ci.low`,
#'   `ci.high`, `p`, `padj`, `p.fwer`)
#' @param normalized use `estimate.norm` (default) or `estimate`
#' @param type `"dot"` or `"bar"`
#' @param order.by effect used to order cell types (default: first effect present)
#' @param show.ci show jackknife 95% intervals
#' @param significance p-value column marking significance (`padj`, `p.fwer`, `p`)
#' @param alpha significance level
#' @param palette named colours per cell type
#' @param plot.theme ggplot2 theme
#' @param subtitle provenance subtitle
#' @param cell.types optional subset / order of cell types
#' @return ggplot2 object
#' @export
plotEffectsPerCellType <- function(df, normalized = TRUE, type = c("dot", "bar"), order.by = NULL, show.ci = TRUE,
                                   significance = c("padj", "p.fwer", "p"), alpha = 0.05, palette = NULL, plot.theme = ggplot2::theme_bw(),
                                   subtitle = NULL, cell.types = NULL) {
  type <- match.arg(type); significance <- match.arg(significance)
  df <- as.data.frame(df)
  effects <- intersect(effectOrder, unique(df$effect))
  df$effect <- factor(df$effect, levels = effects, labels = effectTitle[effects])
  val <- if (normalized && "estimate.norm" %in% names(df)) "estimate.norm" else "estimate"
  df$value <- df[[val]]
  scale.se <- if (val == "estimate.norm") ifelse(is.finite(df$estimate) & df$estimate != 0, df$estimate.norm / df$estimate, NA) else 1
  se <- df$se.jk * abs(scale.se)
  if (val == "estimate.norm" && "se.jk.norm" %in% names(df)) se <- ifelse(is.finite(df$se.jk.norm), df$se.jk.norm, se)   # gsd: its own error of the standardized magnitude
  df$lo <- df$value - 1.96 * se; df$hi <- df$value + 1.96 * se
  df$pval <- df[[significance]]
  df$significant <- ifelse(is.na(df$pval), FALSE, df$pval < alpha)
  if (is.null(order.by)) order.by <- effects[1]
  ord.df <- df[df$effect == effectTitle[order.by] & df$test.id == df$test.id[1], ]
  ord <- ord.df$celltype[order(ord.df$value)]
  ord <- c(ord, setdiff(unique(df$celltype), ord))
  if (!is.null(cell.types)) ord <- intersect(cell.types, ord)
  df$celltype <- factor(df$celltype, levels = ord)
  df$test <- factor(df$test, levels = unique(df$test))
  n.tests <- length(unique(df$test))
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$celltype)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60")
  if (type == "bar") {
    gg <- gg + ggplot2::geom_col(ggplot2::aes(fill = .data$celltype, alpha = .data$significant), width = 0.7, colour = "grey30") +
      ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1), guide = "none")
    if (show.ci) gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), width = 0.3, colour = "grey30", na.rm = TRUE)
  } else {
    if (show.ci) gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi, colour = .data$celltype), width = 0.3, na.rm = TRUE)
    gg <- gg + ggplot2::geom_point(ggplot2::aes(colour = .data$celltype, fill = .data$celltype, shape = .data$significant), size = 2.6, stroke = 0.8) +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19), guide = "none", labels = c(`FALSE` = "not significant", `TRUE` = "significant"))
    gg <- gg + ggplot2::geom_point(data = df[!df$significant, ], ggplot2::aes(colour = .data$celltype), shape = 21, fill = "white", size = 2.6, stroke = 0.8)
  }
  if (!is.null(palette)) {
    pal <- palette[levels(df$celltype)]; pal[is.na(pal)] <- "grey50"
    gg <- gg + ggplot2::scale_colour_manual(values = pal, guide = "none") + ggplot2::scale_fill_manual(values = pal, guide = "none")
  } else gg <- gg + ggplot2::guides(colour = "none", fill = "none")
  gg <- gg + (if (n.tests > 1) ggplot2::facet_wrap(~ test + effect, nrow = n.tests, scales = "free_x") else ggplot2::facet_wrap(~ effect, nrow = 1, scales = "free_x")) +
    plot.theme + ggplot2::labs(x = if (val == "estimate.norm") "normalized effect" else "effect", y = NULL, subtitle = subtitle) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"))
  gg
}

# row order of the diagnostic heatmaps: by group, within a group the most deviating sample on top (ggplot draws the first level at the bottom)
sampleHeatmapOrder <- function(m, groups = NULL) {
  o <- rownames(m)[order(rowMeans(abs(m), na.rm = TRUE), decreasing = FALSE)]
  if (!is.null(groups)) o <- if (is.numeric(groups)) o[order(groups[o], na.last = TRUE)] else o[order(as.character(groups[o]), na.last = TRUE)]
  o
}

#' Heatmap of a per-sample diagnostic (samples x cell types)
#'
#' The common layout of the fit diagnostics of every analysis: samples in rows (grouped by the tested variable, the most
#' deviating first within a group), cell types in columns, a diverging colour scale centred on "as the model expects", and a
#' dot on the tiles beyond `mark` (the samples and cell types worth a look).
#'
#' @param m samples x cell types matrix
#' @param groups optional named vector of sample groups (annotation strip)
#' @param title,legend plot title and legend title
#' @param mark tiles with an absolute value beyond this are marked with a dot (NULL: none)
#' @param subtitle optional subtitle (how to read the values)
#' @param plot.theme ggplot2 theme
#' @param palette colours for the groups
#' @param sample.order optional order of the samples (bottom to top), to align several heatmaps
#' @return ggplot2 object
#' @export
plotSampleHeatmap <- function(m, groups = NULL, title = NULL, legend = "value", mark = NULL, subtitle = NULL, plot.theme = ggplot2::theme_bw(), palette = NULL, sample.order = NULL,
                              group.title = NULL, max.sample.labels = 60) {
  df <- data.frame(sample = rep(rownames(m), ncol(m)), celltype = rep(colnames(m), each = nrow(m)), value = as.vector(m), stringsAsFactors = FALSE)
  has.na <- anyNA(df$value)
  ord.s <- sampleHeatmapOrder(m, groups)
  if (!is.null(sample.order)) ord.s <- c(intersect(sample.order, rownames(m)), setdiff(ord.s, sample.order))
  df$sample <- factor(df$sample, levels = ord.s)
  strip <- " "
  df$celltype <- factor(df$celltype, levels = c(if (!is.null(groups)) strip, colnames(m)))
  lim <- max(abs(df$value), mark %||% 0, na.rm = TRUE)
  sub <- paste(c(subtitle, if (!is.null(mark)) sprintf("dots: beyond %g", mark), if (has.na) "grey: sample not in the cell type's estimate"), collapse = "; ")
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$celltype, y = .data$sample, fill = .data$value)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim), na.value = "grey80", name = legend) +
    ggplot2::scale_x_discrete(drop = FALSE) + plot.theme + ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = if (nzchar(sub)) sub else NULL) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(title.position = "top", order = 1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank(), legend.position = "bottom",
                   legend.title = ggplot2::element_text(size = 8), legend.text = ggplot2::element_text(size = 7), legend.key.width = ggplot2::unit(1.6, "lines"), legend.key.height = ggplot2::unit(0.6, "lines"),
                   plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"), plot.title = ggplot2::element_text(size = 11))
  if (!is.null(mark)) { mk <- df[is.finite(df$value) & abs(df$value) > mark, , drop = FALSE]
    if (nrow(mk)) gg <- gg + ggplot2::geom_point(data = mk, ggplot2::aes(x = .data$celltype, y = .data$sample), inherit.aes = FALSE, shape = 16, size = 1.1, colour = "black") }
  if (length(ord.s) > max.sample.labels) gg <- gg + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::labs(y = sprintf("%d samples", length(ord.s)))
  if (!is.null(groups)) {
    num <- is.numeric(groups)   # a continuous tested variable: a colour scale, samples ordered by it (sampleHeatmapOrder())
    gdf <- data.frame(sample = factor(ord.s, levels = ord.s), group = if (num) as.numeric(groups[ord.s]) else as.character(groups[ord.s]), x = strip, stringsAsFactors = FALSE)
    gdf <- gdf[!is.na(gdf$group), ]
    gg <- gg + ggplot2::geom_point(data = gdf, ggplot2::aes(x = .data$x, y = .data$sample, colour = .data$group), inherit.aes = FALSE, shape = 15, size = 3) +
      ggplot2::labs(colour = group.title)
    if (num) gg <- gg + ggplot2::scale_colour_viridis_c() + ggplot2::guides(colour = ggplot2::guide_colourbar(order = 2, title.position = "top"))
    else {
      gg <- gg + ggplot2::guides(colour = ggplot2::guide_legend(order = 2))
      if (!is.null(palette) && all(unique(gdf$group) %in% names(palette))) gg <- gg + ggplot2::scale_colour_manual(values = palette)
    }
  }
  gg
}

#' Heatmap of leave-one-sample-out influence
#'
#' @param m samples x cell types matrix: change in the effect when the sample is left out
#' @param se optional per-cell-type jackknife SE; when given, the change is shown in SE units
#' @param groups optional named vector of sample groups (annotation strip)
#' @param effect effect name (title)
#' @param plot.theme ggplot2 theme
#' @param palette colours for the groups
#' @return ggplot2 object
#' @export
plotInfluenceHeatmap <- function(m, se = NULL, groups = NULL, effect = "shift", plot.theme = ggplot2::theme_bw(), palette = NULL, group.title = "group",
                                 max.sample.labels = 60) {
  if (!is.null(se)) { se[!is.finite(se) | se == 0] <- NA; m <- sweep(m, 2, se, "/") }
  df <- data.frame(sample = rep(rownames(m), ncol(m)), celltype = rep(colnames(m), each = nrow(m)), value = as.vector(m), stringsAsFactors = FALSE)
  has.na <- anyNA(df$value)   # samples not in a cell type's estimate: grey tiles (white would read as "no change")
  ord.s <- rownames(m)[order(rowMeans(abs(m), na.rm = TRUE), decreasing = TRUE)]
  num <- !is.null(groups) && is.numeric(groups)   # a continuous tested variable: samples ordered by it, a colour scale instead of one colour per value
  if (!is.null(groups)) { ord.s <- if (num) ord.s[order(groups[ord.s], na.last = TRUE)] else ord.s[order(as.character(groups[ord.s]), na.last = TRUE)] }
  df$sample <- factor(df$sample, levels = ord.s)
  strip <- " "                                       # the group annotation gets its own column
  df$celltype <- factor(df$celltype, levels = c(if (!is.null(groups)) strip, colnames(m)))
  lim <- max(abs(df$value), na.rm = TRUE)
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$celltype, y = .data$sample, fill = .data$value)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim), na.value = "grey80",
                                  name = if (is.null(se)) sprintf("change in %s", effect) else sprintf("change in %s\n(SE units)", effect)) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    plot.theme + ggplot2::labs(x = NULL, y = NULL, title = sprintf("Influence of single samples on %s", effect),
                               subtitle = if (has.na) "grey: sample not in the cell type's estimate" else NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"))
  if (length(ord.s) > max.sample.labels) gg <- gg + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::labs(y = sprintf("%d samples%s", length(ord.s), if (num) sprintf(", ordered by %s", group.title) else ""))
  if (!is.null(groups)) {
    gdf <- data.frame(sample = factor(ord.s, levels = ord.s), group = if (num) as.numeric(groups[ord.s]) else as.character(groups[ord.s]), x = strip, stringsAsFactors = FALSE)
    gdf <- gdf[!is.na(gdf$group), ]
    gg <- gg + ggplot2::geom_point(data = gdf, ggplot2::aes(x = .data$x, y = .data$sample, colour = .data$group), inherit.aes = FALSE, shape = 15, size = 3) +
      ggplot2::labs(colour = group.title)
    if (num) gg <- gg + ggplot2::scale_colour_viridis_c()
    else if (!is.null(palette) && all(unique(gdf$group) %in% names(palette))) gg <- gg + ggplot2::scale_colour_manual(values = palette)
  }
  gg
}

#' Shift magnitudes per cell type for the gsd statistic
#'
#' One row per cell type: the magnitude (standardized by default) with its jackknife 95% interval, in front of two faint
#' backgrounds: the range of the magnitude under relabelings of the samples (grey band, 2.5-97.5%: the values expected
#' without a change) and the per-sample jackknife pseudo-values (small points coloured by group; their mean and spread give
#' the estimate and its interval, and a single outlying sample stands out). Significant cell types (by `significance`) come
#' first, ordered by the magnitude, and are drawn filled; the others follow, faded.
#' @param df long results rows of one test (`celltype`, `estimate`, `estimate.norm`, `se.jk`, `se.jk.norm`, p-value columns)
#' @param null relabelings x cell types matrix of the magnitude under relabelings (or NULL)
#' @param influence samples x cell types matrix of the leave-one-out change of the magnitude (or NULL)
#' @param groups named vector of sample groups (colours of the pseudo-values)
#' @param normalized standardized (default) or total shift
#' @param significance,alpha p-value column and level marking significance
#' @param order `"significance"` (default: significant types first, each part by magnitude) or `"magnitude"`
#' @param palette named colours per group; @param plot.theme ggplot2 theme; @param subtitle provenance subtitle
#' @param cell.types optional subset
#' @return ggplot2 object
#' @export
plotGsdShiftMagnitudes <- function(df, null = NULL, influence = NULL, groups = NULL, normalized = TRUE, significance = "padj", alpha = 0.05,
                                   order = c("significance", "magnitude"), palette = NULL, plot.theme = ggplot2::theme_bw(), subtitle = NULL, cell.types = NULL) {
  order <- match.arg(order)
  df <- as.data.frame(df); if (!is.null(cell.types)) df <- df[df$celltype %in% cell.types, ]
  val <- if (normalized) "estimate.norm" else "estimate"; se <- if (normalized) "se.jk.norm" else "se.jk"
  df$value <- df[[val]]; df$lo <- df$value - 1.96 * df[[se]]; df$hi <- df$value + 1.96 * df[[se]]
  df$significant <- is.finite(df[[significance]]) & df[[significance]] < alpha
  ord <- if (order == "significance") df$celltype[order(df$significant, df$value)] else df$celltype[order(df$value)]
  df$celltype <- factor(df$celltype, levels = ord)
  gg <- ggplot2::ggplot(df, ggplot2::aes(y = .data$celltype)) + ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60")
  if (!is.null(null)) {
    ct <- intersect(colnames(null), ord)
    nb <- data.frame(celltype = factor(ct, levels = ord), lo = apply(null[, ct, drop = FALSE], 2, stats::quantile, 0.025, na.rm = TRUE),
                     hi = apply(null[, ct, drop = FALSE], 2, stats::quantile, 0.975, na.rm = TRUE))
    gg <- gg + ggplot2::geom_linerange(data = nb, ggplot2::aes(xmin = .data$lo, xmax = .data$hi), linewidth = 4, colour = "grey88", na.rm = TRUE)
  }
  if (!is.null(influence)) {
    ct <- intersect(colnames(influence), ord)
    pv <- do.call(rbind, lapply(ct, function(c0) {
      d <- influence[, c0]; d <- d[is.finite(d)]; n <- length(d); if (n < 3) return(NULL)
      th <- df$value[df$celltype == c0]
      data.frame(celltype = c0, sample = names(d), value = th - (n - 1) * d)   # jackknife pseudo-values
    }))
    if (!is.null(pv) && nrow(pv)) {
      pv$celltype <- factor(pv$celltype, levels = ord)
      pv$group <- if (!is.null(groups)) as.character(groups[pv$sample]) else "sample"
      gg <- gg + ggplot2::geom_point(data = pv, ggplot2::aes(x = .data$value, colour = .data$group), size = 0.9, alpha = 0.35,
                                     position = ggplot2::position_jitter(width = 0, height = 0.15, seed = 1), na.rm = TRUE)
      if (!is.null(palette)) gg <- gg + ggplot2::scale_colour_manual(values = palette, name = NULL)
      else gg <- gg + ggplot2::scale_colour_brewer(palette = "Set1", name = NULL)
    }
  }
  gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi, alpha = .data$significant), width = 0.3, colour = "grey15", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(x = .data$value, shape = .data$significant, alpha = .data$significant), size = 2.6, stroke = 0.8, colour = "grey15", fill = "grey15") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none")
  # zoom on the estimates and their intervals (outlying pseudo-values are cut at the panel border)
  xr <- range(c(0, df$lo, df$hi, if (!is.null(null)) nb$hi), na.rm = TRUE); xr <- xr + c(-0.08, 0.08) * diff(xr)
  gg + ggplot2::coord_cartesian(xlim = xr) + plot.theme +
    ggplot2::labs(x = if (normalized) "standardized shift (change / donor-level variation)" else "total shift (squared log10 distance between the condition means)",
                  y = NULL, subtitle = subtitle) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"), legend.position = "bottom")
}

#' Top genes of the cell-type shifts (gsd)
#'
#' Dot plot of the genes with the largest per-gene z-scores in the selected cell types. A gene changed in a cell type at a
#' gene-level FDR of `gene.fdr` (BH on the normal p-values of its z-scores within the cell type) is a dot sized by |z| and
#' coloured by the change (log2 fold change of the profile, from the log10 change); otherwise a small grey dot (tested, not
#' changed) or nothing (not tested in that cell type). Cell types are ordered as given.
#' @param genes per cell type data.frames (`gene`, `beta`, `z`), as `res$genes[[test]]`
#' @param cell.types cell types to show (in this order)
#' @param n.genes top genes per cell type (default 10)
#' @param gene.fdr gene-level FDR of the coloured dots (default 0.1)
#' @param plot.theme ggplot2 theme
#' @return ggplot2 object
#' @export
plotGsdShiftGenes <- function(genes, cell.types, n.genes = 10, gene.fdr = 0.1, plot.theme = ggplot2::theme_bw()) {
  cell.types <- intersect(cell.types, names(genes)); if (!length(cell.types)) stop("no gene tables for these cell types")
  top <- unique(unlist(lapply(cell.types, function(ct) utils::head(genes[[ct]]$gene, n.genes))))
  d <- do.call(rbind, lapply(cell.types, function(ct) { g <- genes[[ct]]; q <- stats::p.adjust(2 * stats::pnorm(-abs(g$z)), "BH"); k <- g$gene %in% top
    if (!any(k)) return(NULL)
    data.frame(celltype = ct, gene = g$gene[k], lfc = g$beta[k] * log2(10), z = abs(g$z[k]), changed = q[k] < gene.fdr) }))
  rank.gene <- vapply(top, function(gn) match(TRUE, vapply(cell.types, function(ct) gn %in% utils::head(genes[[ct]]$gene, n.genes), logical(1))), integer(1))
  d$gene <- factor(d$gene, levels = rev(top[order(rank.gene)])); d$celltype <- factor(d$celltype, levels = cell.types)
  lim <- max(abs(d$lfc[d$changed]), 0.1, na.rm = TRUE); d$lfc <- pmax(-lim, pmin(lim, d$lfc))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$celltype, y = .data$gene)) +
    ggplot2::geom_point(data = d[!d$changed, , drop = FALSE], size = 0.7, colour = "grey70") +
    ggplot2::geom_point(data = d[d$changed, , drop = FALSE], ggplot2::aes(size = .data$z, colour = .data$lfc)) +
    ggplot2::scale_colour_gradient2(low = "#2166ac", mid = "grey80", high = "#b2182b", limits = c(-lim, lim), name = "log2 FC") +
    ggplot2::scale_size_area(max_size = 3.5, name = "|z|") + plot.theme +
    ggplot2::labs(x = NULL, y = NULL, caption = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), axis.text.y = ggplot2::element_text(size = 7),
                   panel.grid.major = ggplot2::element_line(colour = "grey94"))
}

#' Cell-type shift magnitudes of a contrast test (one effect)
#'
#' One row per cell type, ordered by the size of the effect. The point is the normalized effect of the distance model
#' (`shift`: the squared distance between the condition means over the summed within-group dispersions, shift / (s.ref + s.alt);
#' `var`: the dispersion change log2(s.alt / s.ref)), filled when significant (for the shift of a gsd run, by the gsd test).
#' Optional: the 95 % range of the effect over the test's relabelings (what noise alone gives), the jackknife 95 % interval (how
#' firmly the samples determine the size; not a test), the p-values, and a margin with the number of genes changed at a gene-level
#' FDR (gsd runs).
#' @param w wide results table of one contrast test (`res$wide[[test]]`)
#' @param effect `"shift"` or `"var"`
#' @param show.null,show.ci,show.genes see [Cacoa] `plotExpressionShiftMagnitudes()`
#' @param show.pvalues `"adjusted"`, `"raw"` or `"none"`
#' @param gene.counts named vector of changed genes per cell type (for `show.genes`)
#' @param gene.fdr FDR of the gene counts (label only)
#' @param significance column prefix marking significance: `"padj"`, `"p.fwer"` or `"p"`
#' @param alpha significance level
#' @param sig.label legend label of the significant points
#' @param palette cell-type colours
#' @param plot.theme ggplot2 theme
#' @param subtitle optional subtitle
#' @param cell.types optional subset of cell types
#' @param order `"magnitude"` (default) or `"significance"`
#' @param legend.position legend position (default "bottom"; "none" hides it)
#' @return ggplot2 object
#' @export
plotShiftMagnitudes <- function(w, effect = c("shift", "var"), show.null = TRUE, show.ci = FALSE, show.pvalues = c("adjusted", "raw", "none"),
                                show.genes = FALSE, gene.counts = NULL, gene.fdr = 0.1, significance = c("padj", "p.fwer", "p"), alpha = 0.05,
                                sig.label = NULL, palette = NULL, plot.theme = ggplot2::theme_bw(), subtitle = NULL, cell.types = NULL,
                                order = c("magnitude", "significance"), legend.position = "bottom") {
  effect <- match.arg(effect); show.pvalues <- match.arg(show.pvalues); significance <- match.arg(significance); order <- match.arg(order)
  w <- as.data.frame(w); if (!is.null(cell.types)) w <- w[w$celltype %in% cell.types, , drop = FALSE]
  S <- w$s.ref + w$s.alt
  if (effect == "shift") {
    val <- w$shift.norm; se <- w$se.shift / S; nlo <- w$null.shift.lo; nhi <- w$null.shift.hi
    xlab <- "normalized shift (distance between condition means / within-group spread)"
  } else {
    val <- w$var.norm; se <- w$se.var / (4 * log(2)) * (1 / w$s.alt + 1 / w$s.ref); nlo <- w$null.var.lo; nhi <- w$null.var.hi
    xlab <- "dispersion change, log2(s.alt / s.ref)"
    if (all(!is.finite(val))) stop("the dispersion change is not estimable for this design (e.g. paired samples: the dispersions are not identified)")
  }
  pcol <- paste0(significance, ".", effect); pv <- w[[pcol]] %||% rep(NA_real_, nrow(w))
  d <- data.frame(celltype = w$celltype, value = val, lo = val - 1.96 * se, hi = val + 1.96 * se, nlo = nlo %||% NA_real_, nhi = nhi %||% NA_real_,
                  p.show = if (show.pvalues == "raw") w[[paste0("p.", effect)]] %||% NA_real_ else w[[paste0("padj.", effect)]] %||% NA_real_,
                  sig = is.finite(pv) & pv < alpha, stringsAsFactors = FALSE)
  d$tied <- if (effect == "var" && "cell.assoc" %in% names(w)) is.finite(w$cell.assoc) & w$cell.assoc > 0.5 else FALSE
  if (!any(is.finite(d$value))) stop("no cell type with an estimate of the ", effect)
  # cell types tested (gsd) without a magnitude (too few samples with enough cells for the distance model): bottom rows, p-value only
  d$nomag <- !is.finite(d$value) & is.finite(d$p.show)
  d <- d[is.finite(d$value) | d$nomag, , drop = FALSE]
  if (any(d$tied)) subtitle <- paste(c(subtitle, "\u2020 cell numbers tied to the contrast: the dispersion test may be anticonservative"), collapse = "\n")
  key <- if (effect == "var") abs(d$value) else d$value
  ord <- if (order == "magnitude") d$celltype[order(!d$nomag, key)] else d$celltype[order(d$sig, -d$p.show, key)]
  d$celltype <- factor(d$celltype, levels = ord)
  sig.lab <- sig.label %||% sprintf("%s < %.2g", significance, alpha)
  d$sig.f <- factor(ifelse(d$sig, sig.lab, "not significant"), levels = c(sig.lab, "not significant"))
  pal <- palette[levels(d$celltype)] %||% stats::setNames(rep("grey25", nlevels(d$celltype)), levels(d$celltype)); pal[is.na(pal)] <- "grey40"
  gg <- ggplot2::ggplot(d, ggplot2::aes(y = .data$celltype)) + ggplot2::geom_vline(xintercept = 0, colour = "grey70")
  if (show.null && any(is.finite(d$nlo)))
    gg <- gg + ggplot2::geom_linerange(ggplot2::aes(xmin = .data$nlo, xmax = .data$nhi, alpha = "noise: 95 % of relabelings"), linewidth = 4.5, colour = "grey88", na.rm = TRUE) +
      ggplot2::scale_alpha_manual(values = 1, name = NULL, guide = ggplot2::guide_legend(override.aes = list(linewidth = 4.5, colour = "grey88"), order = 2))
  if (show.ci)
    gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi, linetype = "jackknife 95 % interval"), width = 0.25, colour = "grey45", na.rm = TRUE) +
      ggplot2::scale_linetype_manual(values = "solid", name = NULL, guide = ggplot2::guide_legend(order = 3))
  gg <- gg + ggplot2::geom_point(ggplot2::aes(x = .data$value, colour = .data$celltype, fill = .data$celltype, shape = .data$sig.f), size = 2.8, stroke = 0.9, show.legend = c(shape = TRUE, linetype = FALSE, alpha = FALSE, colour = FALSE, fill = FALSE), na.rm = TRUE) +
    ggplot2::geom_point(data = d[!d$sig, , drop = FALSE], ggplot2::aes(x = .data$value, colour = .data$celltype), shape = 21, fill = "white", size = 2.8, stroke = 0.9,
                        show.legend = FALSE, na.rm = TRUE) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(21, 21), c(sig.lab, "not significant")), drop = FALSE, name = NULL,
                                guide = ggplot2::guide_legend(override.aes = list(colour = "grey25", fill = c("grey25", "white")), order = 1)) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") + ggplot2::scale_fill_manual(values = pal, guide = "none")
  if (any(d$nomag))
    gg <- gg + ggplot2::geom_text(data = d[d$nomag, , drop = FALSE], ggplot2::aes(x = 0, label = "size n/a"), hjust = -0.1, size = 2.5, colour = "grey45", fontface = "italic")
  xr <- range(c(0, d$value, if (show.ci) c(d$lo, d$hi), if (show.null) c(d$nlo, d$nhi)), na.rm = TRUE)
  if (show.pvalues != "none") {
    xt <- xr[2] + 0.07 * diff(xr)   # clear of the largest point
    gg <- gg + ggplot2::geom_text(ggplot2::aes(x = xt, label = paste0(ifelse(is.finite(.data$p.show), format(signifBelow(.data$p.show, alpha), drop0trailing = TRUE), ""), ifelse(.data$tied, " \u2020", ""))),
                                  hjust = 0, size = 2.6, colour = "grey35") +
      ggplot2::annotate("text", x = xt, y = nlevels(d$celltype) + 0.7, label = if (show.pvalues == "raw") "p" else "padj", hjust = 0, size = 2.6, colour = "grey35", fontface = "italic")
    xr[2] <- xr[2] + 0.18 * diff(xr)
  }
  gg <- gg + ggplot2::coord_cartesian(xlim = xr, ylim = c(1, nlevels(d$celltype) + if (show.pvalues != "none") 0.7 else 0), clip = "off") +
    plot.theme + ggplot2::labs(x = xlab, y = NULL, subtitle = subtitle) +
    ggplot2::theme(legend.position = legend.position, plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
                   panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
  if (!show.genes) return(gg)
  if (is.null(gene.counts)) { warning("no gene tables (gene counts need a gsd run): show.genes ignored"); return(gg) }
  g <- data.frame(celltype = factor(levels(d$celltype), levels = levels(d$celltype)), n = as.numeric(gene.counts[levels(d$celltype)]))
  g$sig <- d$sig[match(g$celltype, d$celltype)]
  gm <- ggplot2::ggplot(g, ggplot2::aes(y = .data$celltype, x = .data$n)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$sig), width = 0.6, colour = "grey30", linewidth = 0.2, orientation = "y", na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = ifelse(is.finite(.data$n), round(.data$n), "")), hjust = -0.15, size = 2.3, colour = "grey15") +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "grey75"), guide = "none") +
    ggplot2::scale_x_continuous(trans = "sqrt", expand = ggplot2::expansion(mult = c(0, 0.05)), guide = ggplot2::guide_axis(check.overlap = TRUE),
                                breaks = function(lim) { b <- pretty(c(0, lim[2]), n = 2); b[b <= lim[2]] }) +
    ggplot2::coord_cartesian(ylim = c(1, nlevels(d$celltype) + if (show.pvalues != "none") 0.7 else 0)) +
    plot.theme + ggplot2::labs(x = sprintf("genes, FDR %g %%\n(grey: significant)", 100 * gene.fdr), y = NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 6), axis.title.x = ggplot2::element_text(size = 7))
  top <- cowplot::plot_grid(gg + ggplot2::theme(legend.position = "none"), gm, nrow = 1, rel_widths = c(1, 0.16), align = "h", axis = "tb")
  if (identical(legend.position, "none")) return(top)
  cowplot::plot_grid(top, cowplot::get_legend(gg), ncol = 1, rel_heights = c(1, 0.07))
}
