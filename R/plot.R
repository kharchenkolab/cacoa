#' @import ggrepel
#' @import tibble
#' @import cowplot
#' @importFrom reshape2 melt
NULL

plotNCellRegression <- function(n, n.total, y.lab="N", legend.position="right") {
  p.df <- data.frame(N=n) %>% as_tibble(rownames="Type") %>%
    mutate(NCells=n.total[Type])

  ggplot(p.df, aes(x=NCells, y=N)) +
    geom_point(aes(color=Type)) +
    geom_label_repel(aes(label=Type), size=2, min.segment.length=0.1, box.padding=0, label.size=0, max.iter=300, fill=alpha("white", 0.4)) +
    scale_x_log10() +
    ylim(0, max(p.df$N)) + labs(x="Number of cells", y=y.lab) +
    theme(legend.position=legend.position, legend.justification=legend.position, legend.background=element_rect(fill=alpha("white", 0.4))) +
    guides(color=guide_legend(title="Cell type"))
}

#' @title Plot onthology terms
#' @description Plot onthology terms as a function of both number of DE genes, and number of cells.
#' @param ont.res Onthology results in a data frame
#' @param type Onthology, must be either "GO" or "DO" (default=NULL)
#' @param de.genes.filtered
#' @param cell.groups Vector indicating cell group sizes with cell group names
#' @param show.legend Include legend in plot (default=T)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="bottom")
#' @param label.x.pos Plot label position on x axis (default=0.01)
#' @param label.y.pos Plot label position on y axis (default=1)
#' @param rel_heights Vector indicating relative heights for plots. Only relevant if show.legend=T. See cowplot::plot_grid for more info (default=c(2.5, 0.5))
#' @param scale Scaling of plots, adjust if e.g. label is misplaced. See cowplot::plot_grid for more info (default=0.93)
#' @return A ggplot2 object
#' @export
plotOnthologyTerms <- function(type=NULL, ont.res, de.genes.filtered, cell.groups, legend.position="bottom", label.x.pos=0.01, label.y.pos=1, rel_heights = c(2.5, 0.5), scale = 0.93, show.legend=T) {
  if(length(unique(ont.res$Type))==1) stop("The input only contains one cell type.")

  if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

  cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.genes.filtered)]
  leg <- cowplot::get_legend(plotNCellRegression(sapply(de.genes.filtered, length), cell.groups, legend.position=legend.position))

  pg <- cowplot::plot_grid(
    ont.res$Type %>% table %>% c %>%
      plotNCellRegression(sapply(de.genes.filtered, length), y.lab=NULL, legend.position="none") +
      geom_smooth(method=MASS::rlm, se=0, color="black", size=0.5) +
      scale_x_continuous(name="Number of highly-expressed DE genes"),
    ont.res$Type %>% table %>% c %>%
      plotNCellRegression(cell.groups, y.lab=NULL, legend.position="none") +
      geom_smooth(method=MASS::rlm, se=0, color="black", size=0.5),
    ncol=1, labels=c("a", "b"), label_x = label.x.pos, label_y = label.y.pos
  )

  if(type=="GO") {
    y_lab <- "Number of GO terms"
  } else if(type=="DO") {
    y_lab <- "Number of DO terms"
  }
  if(show.legend) {
    p <- cowplot::plot_grid(pg, leg, nrow=2, rel_heights = rel_heights, scale=scale) + draw_label(y_lab, x=0, y=0.6, vjust= 1.5, angle=90)
  } else {
    p <- cowplot::plot_grid(pg, nrow=1, scale=scale) + draw_label(y_lab, x=0, y=0.6, vjust= 1.5, angle=90)
  }
  return(p)
}

#' @title Plot DE genes
#' @description Plot number of DE genes as a function of number of cells
#' @param de.raw
#' @param de.genes.filtered
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param show.legend Include legend in plot (default=T)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="bottom")
#' @param p.adjust.cutoff Adjusted P cutoff (default=0.05)
#' @param label.x.pos Plot label position on x axis (default=0.01)
#' @param label.y.pos Plot label position on y axis (default=1)
#' @param rel_heights Relative heights for plots. Only relevant if show.legend=T. See cowplot::plot_grid for more info (default=c(2.5, 0.5))
#' @param scale Scaling of plots, adjust if e.g. label is misplaced. See cowplot::plot_grid for more info (defaul=0.93)
#' @return A ggplot2 object
#' @export
plotDEGenes <- function(de.raw, de.genes.filtered, cell.groups, legend.position="bottom", p.adjust.cutoff=0.05, label.x.pos=0.01, label.y.pos=1, rel_heights=c(2.5, 0.5), scale=0.93, show.legend=T) {
  cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.genes.filtered)]
  leg <- cowplot::get_legend(plotNCellRegression(1:length(cell.groups) %>% setNames(names(cell.groups)), cell.groups, legend.position=legend.position))

  pg <- cowplot::plot_grid(
    sapply(de.raw, function(n) n %>% dplyr::filter(padj <= p.adjust.cutoff) %>% nrow) %>%
      plotNCellRegression(cell.groups, y.lab="Significant DE genes", legend.position="none") +
      xlab("") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="black", size=0.5),
    sapply(de.genes.filtered, length) %>%
      plotNCellRegression(cell.groups, y.lab="Highly-expressed DE genes", legend.position="none") +
      geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="black", size=0.5),
    ncol=1, labels=c("a","b"), label_x = label.x.pos, label_y = label.y.pos
  )
  if(show.legend) {
    p <- cowplot::plot_grid(pg, leg, nrow=2, rel_heights = rel_heights, scale=scale)
  } else {
    p <- cowplot::plot_grid(pg, nrow=1, scale=scale)
  }
  return(p)
}

#' @title Plot heatmap
#' @param df Data frame with plot data
#' @param color.per.type Colors per cell type (default=NULL)
#' @param row.order Forced row order (default=NULL)
#' @param col.order Forced column order (default=NULL)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param legend.key.width (default=unit(8, "pt))
#' @param legend.title Title on plot (default="-log10(p-value)")
#' @param x.axis.position Position of x axis (default="top")
#' @return A ggplot2 object
#' @export
plotHeatmap <- function(df, color.per.type=NULL, row.order=NULL, col.order=F, legend.position="right", legend.key.width=unit(8, "pt"), legend.title="-log10(p-value)", x.axis.position="top") {
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

  df %<>% as_tibble(rownames="GO") %>%
    melt(id.vars="GO", variable.name="Type", value.name="p.value")

  if (!is.logical(row.order)) {
    df %<>% mutate(GO=factor(GO, levels=row.order))
  }

  if (!is.logical(col.order)) {
    df %<>% mutate(Type=factor(Type, levels=col.order))
  }

  if (is.null(color.per.type)) {
    color.per.type <- "black"
  } else {
    color.per.type <- color.per.type[levels(df$Type)]
  }

  gg <- ggplot(df) + geom_tile(aes(x=Type, y=GO, fill=pmin(p.value, m)), colour = "grey50") +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, color=color.per.type),
          axis.text=element_text(size=8), axis.ticks=element_blank(), axis.title=element_blank()) +
    scale_fill_distiller(palette="RdYlBu", limits=c(0, m)) +
    guides(fill=guide_colorbar(title=legend.title, title.position="left", title.theme=element_text(angle=90, hjust=0.5))) +
    scale_y_discrete(position="right", expand=c(0, 0)) +
    scale_x_discrete(expand=c(0, 0), position=x.axis.position) +
    theme(legend.position=legend.position, legend.key.width=legend.key.width, legend.background=element_blank())

  return(gg)
}

#' @title Plot pathway distribution
#' @description Bar plot of onthology pathways per cell type
#' @param ont.res Onthology resuls from estimateOnthology
#' @param type Onthology, must be either "GO" or "DO" (default=NULL)
#' @return A ggplot2 object
#' @export
plotPathwayDistribution <- function(type=NULL, ont.res, cell.groups) {
  if(type=="GO") {
    p_df <- table(ont.res$Type, ont.res$GO) %>%
      cbind %>%
      as_tibble(rownames="Type") %>%
      reshape2::melt(id.vars="Type", variable.name="GO", value.name="N") %>%
      mutate(Type=factor(Type, levels=unique(cell.groups) %>%
                           .[order(.)]))

    gg <- ggplot(p_df) +
      geom_bar(aes(x=Type, y=N, fill=GO), stat="identity")
  } else if(type=="DO") {
    p_df <- table(ont.res$Type) %>%
      cbind %>%
      as_tibble(rownames="Type") %>%
      reshape2::melt(id.vars="Type", variable.name="DO", value.name="N") %>%
      mutate(Type=factor(Type, levels=unique(cell.groups) %>%
                           .[order(.)]))

    gg <- ggplot(p_df) +
      geom_bar(aes(x=Type, y=N), stat="identity")
  } else {
    stop("'type' must be either 'GO' or 'DO'.")
  }

  m <- sapply(p_df$Type %>%
                unique, function(m) p_df %>%
                dplyr::filter(Type==m) %>%
                dplyr::select(N) %>%
                sum) %>%
    max %>%
    plyr::round_any(10, ceiling)

  gg +
    scale_y_continuous(expand=c(0, 0), limits=c(0, m)) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          panel.grid.major.x=element_blank(), legend.position=c(0.99, 0.99), legend.justification=c(1, 1)) +
    labs(x="", y="No. of pathways")
}

#' @title Plot onthology heatmap
#' @description Plot a heatmap of onthology P values per cell type
#' @param type Onthology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default=NULL)
#' @param ont.res Onthology resuls from estimateOnthology
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="left")
#' @param order Order of rows in heatmap. Can be 'unique' (only show pathways that are unique for any cell type); 'unique-max-row' (same as 'unique' but ordered by P value); 'all-max-rowsum' (all pathways ordered by cumulative P value for all cell types); 'all-max-row' (all pathways ordered by max P value) (default="all-max-row")
#' @param n Number of pathways to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
#' @export
plotOnthologyHeatmap <- function(type, ont.res, legend.position = "left", order = "all-max-row", n = 10) {
  if(type=="BP" | type=="CC" | type=="MF") {
    ont.sum <- getOnthologySummary(type, ont.res %>% filter(GO==type))
  } else if(type=="DO") {
    ont.sum <- getOnthologySummary(type, ont.res)
  }

  if(!order %in% c("unique","unique-max-row","all-max-rowsum","all-max-row")) stop("'order' must be one of the following: 'unique', 'unique-max-row', 'all-max-rowsum', 'all-max-row'.")

  if(order=="unique") {
    ont.sum %>% .[rowSums(. > 0.1) == 1,] %>%
      .[, colSums(.) > 0]  %>%
      plotHeatmap(legend.position=legend.position)
  } else if(order=="unique-max-row") {
    ont.sum %>% .[rowSums(. > 0.1) == 1,] %>%
      .[, colSums(.) > 0] %>%
      .[match(apply(., 1, max) %>%
                .[order(., decreasing = F)] %>%
                names, rownames(.)),] %>%
      plotHeatmap(legend.position=legend.position, row.order=T)
  } else if(order=="all-max-rowsum") {
    ont.sum %>% .[, colSums(.) > 0] %>%
      .[match(rowSums(.)[rowSums(.)>0] %>%
                .[order(., decreasing = F)] %>%
                names, rownames(.)),] %>%
      tail(n) %>%
      plotHeatmap(legend.position=legend.position, row.order=T)
  } else if(order=="all-max-row") {
    ont.sum %>% .[, colSums(.) > 0] %>%
      .[match(apply(., 1, max) %>%
                .[order(., decreasing = F)] %>%
                names, rownames(.)),] %>%
      tail(n) %>%
      plotHeatmap(legend.position=legend.position, row.order=T)
  }
}

#' @title Plot onthology correlations
#' @description Plot correlation matrix for onthologies between cell types
#' @param ont.res Onthology resuls from estimateOnthology
#' @param type Onthology, must be either "GO" or "DO" (default=NULL)
#' @return A ggplot2 object
#' @export
plotOnthologyCorrelations <- function(type=NULL, ont.res) {
  if(type=="GO") {
    pathway_df <- names(ont.res) %>%
      lapply(function(ng) lapply(ont.res[[ng]] %>% .$Type %>% levels, function(nt) {
        tibble(Pathway=ont.res[[ng]] %>%
                 dplyr::filter(Type==nt) %>%
                 .$Description, PathwayType=ng, Type=nt)
      }) %>% bind_rows) %>% bind_rows
  } else if(type=="DO") {
    pathway_df <- lapply(ont.res %>% names, function(nt) {
      tibble(Pathway=ont.res[[nt]] %>%
               .$Description, Type=nt)
    }) %>% bind_rows
  } else {
    stop("'type' must be either 'GO' or 'DO'.")
  }
  path_bin <- pathway_df %>%
    dplyr::select(Pathway, Type) %>%
    mutate(X=1) %>%
    tidyr::spread(Pathway, X) %>%
    as.data.frame() %>%
    set_rownames(.$Type) %>%
    .[, 2:ncol(.)] %>%
    as.matrix()
  path_bin[is.na(path_bin)] <- 0

  p_mat <- (1 - (path_bin %>% dist(method="binary") %>% as.matrix)) %>% pmin(0.5)

  t_tree <- dist(p_mat) %>% hclust()
  t_order <- t_tree %$% labels[order]
  t_cls <- cutree(t_tree, h=0.7) %>% .[t_order]
  t_cls[t_cls %in% names(which(table(t_cls) < 5))] <- max(t_cls) + 1

  t_cl_lengths <- rle(t_cls)$lengths %>% rev

  diag(p_mat) <- 1
  plotHeatmap(p_mat, color.per.type=NULL, row.order=t_order, col.order=rev(t_order), legend.title="Similarity") +
    scale_fill_distiller(palette="RdYlBu", limits=c(0, 0.5)) +
    geom_vline(aes(xintercept=x), data.frame(x=cumsum(t_cl_lengths)[t_cl_lengths > 1] + 0.5)) +
    geom_hline(aes(yintercept=x), data.frame(x=cumsum(t_cl_lengths)[t_cl_lengths > 1] + 0.5))
}
