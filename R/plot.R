#' @import ggrepel
#' @import tibble
#' @import cowplot
#' @import dplyr
#' @import magrittr
#' @importFrom reshape2 melt
NULL

theme_legend_position <- function(position) {
  theme(legend.position=position, legend.justification=position)
}

plotNCellRegression <- function(n, n.total, y.lab="N", legend.position="right", label=T, size=5, palette=NULL) {
  p.df <- data.frame(N=n) %>% tibble::as_tibble(rownames="Type") %>%
    mutate(NCells=n.total[Type])

  gg <- ggplot(p.df, aes(x=NCells, y=N)) +
    geom_point(aes(color=Type)) +
    scale_x_log10() +
    ylim(0, max(p.df$N)) +
    labs(x="Number of cells", y=y.lab) +
    theme_bw()

  if(label) gg <- gg + geom_label_repel(aes(label=Type), size=size, min.segment.length=0.1, box.padding=0, label.size=0, max.iter=300, fill=alpha("white", 0.4))

  gg <- gg +
    theme(legend.background=element_rect(fill=alpha("white", 0.4))) +
    theme_legend_position(legend.position) +
    guides(color=guide_legend(title="Cell type"))

  if(!is.null(palette)) gg <- gg+ scale_color_manual(values=palette)

  return(gg)
}

#' @title Plot heatmap
#' @param df Data frame with plot data
#' @param color.per.group Colors per cell group (default=NULL)
#' @param row.order Forced row order (default=NULL)
#' @param col.order Forced column order (default=NULL)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param legend.key.width (default=unit(8, "pt))
#' @param legend.title Title on plot (default="-log10(p-value)")
#' @param x.axis.position Position of x axis (default="top")
#' @return A ggplot2 object
#' @export
plotHeatmap <- function(df, color.per.group=NULL, row.order=NULL, col.order=F, legend.position="right", legend.key.width=unit(8, "pt"), legend.title="-log10(p-value)", x.axis.position="top") {
  m <- max(df)

  if (is.null(row.order)) {
    row.order <- rownames(df)[dist(df) %>% hclust() %>% .$order]
  } else if (is.logical(row.order) && row.order) {
    row.order <- rownames(df)
  }

  if (is.null(col.order)) {
    col.order <- colnames(df)[dist(df) %>% hclust() %>% .$order]
  } else if (is.logical(col.order) && col.order) {
    col.order <- colnames(df)
  }

  df %<>% tibble::as_tibble(rownames="Pathway") %>%
    reshape2::melt(id.vars="Pathway", variable.name="Group", value.name="p.value")

  if (!is.logical(row.order)) {
    df %<>% dplyr::mutate(Pathway=factor(Pathway, levels=row.order))
  }

  if (!is.logical(col.order)) {
    df %<>% dplyr::mutate(Group=factor(Group, levels=col.order))
  }

  if (is.null(color.per.group)) {
    color.per.group <- "black"
  } else {
    color.per.group <- color.per.group[levels(df$Group)]
  }

  gg <- ggplot(df) + geom_tile(aes(x=Group, y=Pathway, fill=pmin(p.value, m)), colour = "grey50") +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, color=color.per.group),
          axis.text=element_text(size=8), axis.ticks=element_blank(), axis.title=element_blank()) +
    scale_fill_distiller(palette="RdYlBu", limits=c(0, m)) +
    guides(fill=guide_colorbar(title=legend.title, title.position="left", title.theme=element_text(angle=90, hjust=0.5))) +
    scale_y_discrete(position="right", expand=c(0, 0)) +
    scale_x_discrete(expand=c(0, 0), position=x.axis.position) +
    theme_legend_position(legend.position) +
    theme(legend.key.width=legend.key.width, legend.background=element_blank())

  return(gg)
}

#' @title Plot proportions
#' @description Plot the cell group proportions per sample
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param sample.per.cell Vector indicating sample name with cell names (default: stored vector)
#' @param sample.groups Vector indicating sample groups with sample names (default: stored vector)
#' @param alpha Transparency level on the data points (default: 0.1)
#' @return A ggplot2 object
plotProportions <- function(legend.position = "right", cell.groups, sample.per.cell, sample.groups, notch=FALSE, alpha=0.1, palette=NULL, show.significance = FALSE) {
  df.melt <- data.frame(anno=cell.groups, group=sample.per.cell[match(names(cell.groups), names(sample.per.cell))]) %>%
    table %>%
    rbind %>%
    t %>%
    as.data.frame %>%
    apply(2, as.numeric) %>%
    as.data.frame %>%
    magrittr::divide_by(rowSums(.)) %>%
    magrittr::multiply_by(1e2) %>%
    dplyr::mutate(group = sample.groups[match(levels(sample.per.cell), names(sample.groups))]) %>%
    reshape2::melt(., id.vars="group")

  gg <- ggplot(df.melt, aes(x=variable, y=value, by=group)) +
    geom_boxplot(position=position_dodge(), outlier.shape = NA, notch=notch) +
    ylab("% cells per sample") +
    xlab("") +
    theme_bw() +
    theme_legend_position(legend.position) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_blank()) +
    geom_point(position=position_jitterdodge(jitter.width=0.15), aes(col=group), alpha=alpha) +
    scale_y_continuous( expand=c(0, max(df.melt$value) * 0.1), limits=c(0, (max(df.melt$value) + max(df.melt$value) * 0.05 )))  #expand=c(0, 0),

  if(show.significance) gg <- gg + stat_compare_means(aes(group = group), label = "p.signif")  # willcox test


  if(!is.null(palette)) gg <- gg+ scale_color_manual(values=palette)
  gg
}


#' @title Plot proportions for a subset of cell types
#' @description Plot the cell group proportions per sample
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param sample.per.cell Vector indicating sample name with cell names (default: stored vector)
#' @param sample.groups Vector indicating sample groups with sample names (default: stored vector)
#' @param cells.to.remain Vector of cell types to remain in the composition
#' @param cells.to.remove Vector of cell types to remove from the composition
#' #' @param alpha Transparency level on the data points (default: 0.1)
#' @return A ggplot2 object
plotProportionsSubset <- function(legend.position = "right",
                                  cell.groups,
                                  sample.per.cell,
                                  sample.groups,
                                  cells.to.remove,
                                  cells.to.remain,
                                  notch = FALSE,
                                  alpha = 0.1, palette=NULL,
                                  show.significance = FALSE) {
  df.melt <- data.frame(anno=cell.groups, group=sample.per.cell[match(names(cell.groups), names(sample.per.cell))]) %>%
    table  %>%
    rbind

  if(!is.null(cells.to.remove)) df.melt = df.melt[!(rownames(df.melt) %in% cells.to.remove),]
  if(!is.null(cells.to.remain)) df.melt = df.melt[rownames(df.melt) %in% cells.to.remain,]

  df.melt <- df.melt %>%
    t %>%
    as.data.frame  %>%
    apply(2, as.numeric) %>%
    as.data.frame %>%
    magrittr::divide_by(rowSums(.)) %>%
    magrittr::multiply_by(1e2) %>%
    dplyr::mutate(group = sample.groups[match(levels(sample.per.cell), names(sample.groups))]) %>%
    reshape2::melt(., id.vars="group")

  gg <- ggplot(df.melt, aes(x=variable, y=value, by=group)) +
    geom_boxplot(position=position_dodge(), outlier.shape = NA, notch=notch) +
    ylab("% cells per sample") +
    xlab("") +
    theme_bw() +
    theme_legend_position(legend.position) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_blank()) +
    geom_point(position=position_jitterdodge(jitter.width=0.15), aes(col=group), alpha=alpha) +
    scale_y_continuous(limits=c(0, (max(df.melt$value) + 5))) +
    stat_compare_means(aes(group = group), label = "p.signif")

  if(show.significance) gg <- gg + stat_compare_means(aes(group = group), label = "p.signif")  # willcox test

  if(!is.null(palette)) gg <- gg+scale_color_manual(values=palette)
  return(gg)
}


#' @title Plot Expression Shift Magnitudes
#' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
#' @param name Test results to plot (default=expression.shifts)
#' @param size.norm Plot size normalized results. Requires cell.groups, and sample.per.cell (default=F)
#' @param notch Show notches in plot, see ggplot2::geom_boxplot for more info (default=T)
#' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
#' @param sample.per.cell Named sample factor with cell names (default: stored vector)
#' @return A ggplot2 object
plotExpressionShiftMagnitudes <- function(cluster.shifts, size.norm = F, notch = T, cell.groups = NULL, sample.per.cell = NULL, palette=NULL) {
  if (!size.norm) {
    m <- max(abs(cluster.shifts$value - 1))

    gg <- ggplot(na.omit(cluster.shifts), aes(x=as.factor(Type), y=value, fill=Type)) +
      geom_boxplot(notch=notch, outlier.shape=NA)  +
      geom_jitter(position=position_jitter(0.1), color='gray30', show.legend=FALSE,alpha=0.1,size=0.8) +
      theme_bw() +
      theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.5), legend.position = 'none') +
      labs(x="", y="normalized expression distance") +
      #ylim(c(1 - m, 1 + m)) +
      geom_hline(yintercept=1, linetype="dashed", color = "black")
    if(!is.null(palette)) { gg <- gg + scale_fill_manual(values=palette) }
  } else {
    if (length(setdiff(names(cell.groups), names(sample.per.cell)))>0) warning("Cell names in 'cell.groups' and 'sample.per.cell' are not identical, plotting intersect.")

    cct <- table(cell.groups, sample.per.cell[names(cell.groups)])
    x <- tapply(cluster.shifts$value, cluster.shifts$Type, median)
    odf <- data.frame(cell=names(x),size=rowSums(cct)[names(x)],md=x)

    #m <- max(abs(odf$md - 1))

    gg <- ggplot(odf, aes(size,md,color=cell,label=cell)) +
      ggrepel::geom_text_repel() +
      geom_point() +
      guides(color=F) +
      xlab("Cluster size") +
      theme_bw() +
      theme(legend.position = 'none')
      ylab("Median normalized distance") +
      #ylim(c(1 - m,1 + m)) +
      geom_hline(yintercept=1, linetype="dashed", color = "black")

    if(!is.null(palette)) { gg <- gg + scale_color_manual(values=palette) }
  }


  return(gg)
}

#' @title Plot Expression Shift Z Scores
#' @description  Plot results from estimateExpressionShiftZScores
#' @param plot.df Test results to plot
#' @param size.norm Plot size normalized results. Requires cell.groups, and sample.per.cell (default=F)
#' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
#' @param sample.per.cell Named sample factor with cell names (default: stored vector)
#' @param label Plot labels on size normalized plots (default=T)
#' @return A ggplot2 object
plotExpressionShiftZScores <- function(plot.df, size.norm = F, cell.groups = NULL, sample.per.cell = NULL) {
  if (!size.norm) {
    gg <- ggplot(plot.df, aes(x=Type, y=distance)) +
      geom_boxplot(outlier.alpha=0, show.legend=F) +
      geom_hline(aes(yintercept=split(distance, Type) %>% sapply(median) %>% median(), linetype="Median"), color="darkred", size=1) +
      labs(x="", y="Normalized distance") +
      scale_y_continuous(expand=c(0, 0)) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=9)) +
      scale_linetype_manual(name = "", values=1)
  } else {
    if(length(setdiff(names(cell.groups), names(sample.per.cell)))>0) warning("Cell names in 'cell.groups' and 'sample.per.cell' are not identical, plotting intersect.")

    cct <- table(cell.groups, sample.per.cell[names(cell.groups)])
    x <- tapply(plot.df$distance, plot.df$Type, median)
    odf <- data.frame(cell=names(x),size=rowSums(cct)[names(x)],md=x)

    gg <- ggplot(odf, aes(size,md,color=cell,label=cell)) +
      ggrepel::geom_text_repel() +
      geom_point() +
      guides(color=F) +
      xlab("Cluster size") +
      theme_bw() +
      ylab("Median normalized distance") +
      geom_hline(aes(yintercept=split(plot.df$distance, plot.df$Type) %>% sapply(median) %>% median(), linetype="Median"), color="darkred", size=1) +
      scale_linetype_manual(name = "", values=1)
  }
  return(gg)
}
