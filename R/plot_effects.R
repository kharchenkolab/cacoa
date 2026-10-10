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
  df$lo <- df$value - 1.96 * df$se.jk * abs(scale.se); df$hi <- df$value + 1.96 * df$se.jk * abs(scale.se)
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
      ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19), guide = "none", labels = c("not significant", "significant"))
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
plotInfluenceHeatmap <- function(m, se = NULL, groups = NULL, effect = "shift", plot.theme = ggplot2::theme_bw(), palette = NULL) {
  if (!is.null(se)) { se[!is.finite(se) | se == 0] <- NA; m <- sweep(m, 2, se, "/") }
  df <- data.frame(sample = rep(rownames(m), ncol(m)), celltype = rep(colnames(m), each = nrow(m)), value = as.vector(m), stringsAsFactors = FALSE)
  df <- df[!is.na(df$value), ]
  ord.s <- rownames(m)[order(rowMeans(abs(m), na.rm = TRUE), decreasing = TRUE)]
  if (!is.null(groups)) { g <- as.character(groups[ord.s]); ord.s <- ord.s[order(g, na.last = TRUE)] }
  df$sample <- factor(df$sample, levels = ord.s)
  strip <- " "                                       # the group annotation gets its own column
  df$celltype <- factor(df$celltype, levels = c(if (!is.null(groups)) strip, colnames(m)))
  lim <- max(abs(df$value), na.rm = TRUE)
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$celltype, y = .data$sample, fill = .data$value)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim),
                                  name = if (is.null(se)) sprintf("change in %s", effect) else sprintf("change in %s\n(SE units)", effect)) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    plot.theme + ggplot2::labs(x = NULL, y = NULL, title = sprintf("Influence of single samples on %s", effect)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
  if (!is.null(groups)) {
    gdf <- data.frame(sample = factor(ord.s, levels = ord.s), group = as.character(groups[ord.s]), x = strip, stringsAsFactors = FALSE)
    gdf <- gdf[!is.na(gdf$group), ]
    gg <- gg + ggplot2::geom_point(data = gdf, ggplot2::aes(x = .data$x, y = .data$sample, colour = .data$group), inherit.aes = FALSE, shape = 15, size = 3) +
      ggplot2::labs(colour = "group")
    if (!is.null(palette) && all(unique(gdf$group) %in% names(palette))) gg <- gg + ggplot2::scale_colour_manual(values = palette)
  }
  gg
}
