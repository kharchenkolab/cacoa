## Per-cell-type drill-down (api note §4.7, plotShiftDetail).

# squared pairwise distances of an adjusted distance matrix, labelled within-ref / within-alt / between
pairDistanceTable <- function(D2, groups) {
  s <- intersect(rownames(D2), names(groups)); g <- as.character(groups[s]); D2 <- D2[s, s]
  idx <- which(upper.tri(D2), arr.ind = TRUE)
  levs <- levels(groups)
  kind <- ifelse(g[idx[, 1]] == g[idx[, 2]], ifelse(g[idx[, 1]] == levs[1], sprintf("within %s", levs[1]), sprintf("within %s", levs[2])), "between")
  data.frame(a = s[idx[, 1]], b = s[idx[, 2]], d2 = D2[idx], kind = factor(kind, levels = c(sprintf("within %s", levs[1]), sprintf("within %s", levs[2]), "between")),
             stringsAsFactors = FALSE)
}

#' Detail panels for one cell type's fit
#'
#' @param eff a fit from `expressionShiftsForModel()$fits[[test]][[cell.type]]`
#' @param adjusted adjusted distance matrix for the cell type (optional; derived from the fit when `NULL`)
#' @param groups named factor / vector of sample groups (two-level contrast) or values (numeric test)
#' @param dist distance type (for the scale of the adjusted distances)
#' @param title title
#' @param palette colours for the groups
#' @param plot.theme ggplot2 theme
#' @param label.samples label the samples in the distance panel (default FALSE)
#' @return cowplot grid
#' @export
plotShiftDetailPanels <- function(eff, adjusted = NULL, groups = NULL, dist = "cor", title = NULL, palette = NULL, plot.theme = ggplot2::theme_bw(),
                                  label.samples = FALSE) {
  if (is.null(adjusted)) adjusted <- adjustedDistanceMatrix(eff$G, nuisanceColumns(eff, NULL), dist)
  D2 <- if (dist == "l2") adjusted^2 else adjusted
  samples <- rownames(D2)
  grp <- if (!is.null(groups)) groups[samples] else NULL
  is.two <- !is.null(grp) && is.factor(grp) && nlevels(grp) == 2
  # (a) MDS of adjusted distances
  mds <- stats::cmdscale(sqrt(pmax(D2, 0)), k = 2)
  df <- data.frame(x = mds[, 1], y = mds[, 2], sample = samples, group = if (is.null(grp)) "all" else grp, stringsAsFactors = FALSE)
  gg.a <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group)) + ggplot2::geom_point(size = 3) +
    plot.theme + ggplot2::labs(x = "MDS 1", y = "MDS 2", title = "adjusted sample distances", colour = NULL) +
    ggplot2::theme(legend.position = "bottom")
  if (label.samples) gg.a <- gg.a + ggrepel::geom_text_repel(ggplot2::aes(label = .data$sample), size = 2.5, show.legend = FALSE, max.overlaps = Inf)
  if (!is.null(palette) && is.two && all(levels(grp) %in% names(palette))) gg.a <- gg.a + ggplot2::scale_colour_manual(values = palette)
  panels <- list(gg.a)
  # (b) pair distances within / between, annotated with the three effects
  if (is.two) {
    pt <- pairDistanceTable(D2, grp)
    ann <- sprintf("shift = %.3g   var = %.3g   total = %.3g", eff$shift, eff$var, eff$total)
    kind.cols <- c("grey70", "grey70", "grey40"); names(kind.cols) <- levels(pt$kind)
    if (!is.null(palette) && all(levels(grp) %in% names(palette))) kind.cols[1:2] <- unname(palette[levels(grp)])
    gg.b <- ggplot2::ggplot(pt, ggplot2::aes(x = .data$kind, y = .data$d2, fill = .data$kind)) + ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      ggplot2::scale_fill_manual(values = kind.cols) +
      ggplot2::geom_jitter(width = 0.15, size = 1, alpha = 0.5) + plot.theme + ggplot2::guides(fill = "none") +
      ggplot2::labs(x = NULL, y = "adjusted squared distance", title = "pair distances", subtitle = ann) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8))
    panels <- c(panels, list(gg.b))
  }
  # (c) per-sample dispersion by group
  if (!is.null(eff$v)) {
    dv <- data.frame(sample = names(eff$v), v = eff$v, group = if (is.null(grp)) "all" else factor(as.character(grp[names(eff$v)]), levels = levels(grp)), stringsAsFactors = FALSE)
    dv <- dv[is.finite(dv$v), ]
    med <- stats::median(dv$v); madv <- stats::mad(dv$v)
    dv$label <- ifelse(abs(dv$v - med) > 3 * madv & madv > 0, dv$sample, "")
    gg.c <- ggplot2::ggplot(dv, ggplot2::aes(x = .data$group, y = .data$v, colour = .data$group)) + ggplot2::geom_boxplot(outlier.shape = NA, colour = "grey50") +
      ggplot2::geom_jitter(width = 0.15, size = 2) + ggrepel::geom_text_repel(ggplot2::aes(label = .data$label), size = 2.5, show.legend = FALSE) +
      plot.theme + ggplot2::guides(colour = "none") + ggplot2::labs(x = NULL, y = "per-sample dispersion", title = "dispersion (leverage-corrected)")
    if (!is.null(palette) && is.two && all(levels(grp) %in% names(palette))) gg.c <- gg.c + ggplot2::scale_colour_manual(values = palette)
    panels <- c(panels, list(gg.c))
  }
  # (d) level x level table for term tests
  if (identical(eff$kind, "term") && !is.null(eff$cell.table)) {
    ct <- eff$cell.table
    dd <- data.frame(a = rep(rownames(ct), ncol(ct)), b = rep(colnames(ct), each = nrow(ct)), value = as.vector(ct))
    gg.d <- ggplot2::ggplot(dd, ggplot2::aes(x = .data$a, y = .data$b, fill = .data$value)) + ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2g", .data$value)), size = 3) + ggplot2::scale_fill_gradient(low = "white", high = "#2c7fb8", name = "E[d2]") +
      plot.theme + ggplot2::labs(x = NULL, y = NULL, title = "model-implied pair distances between levels (diagonal: 2 x dispersion)")
    panels <- c(panels, list(gg.d))
  }
  gg <- cowplot::plot_grid(plotlist = panels, ncol = min(3, length(panels)))
  if (!is.null(title)) gg <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(title, fontface = "bold", x = 0, hjust = 0), gg, ncol = 1, rel_heights = c(0.07, 1))
  gg
}
