#' @import tibble
#' @import cowplot
#' @import dplyr
#' @import magrittr
#' @importFrom reshape2 melt
#' @importFrom igraph get.adjacency graph_from_data_frame renderGraph
#' @importFrom Rgraphviz layoutGraph
#' @importFrom graph graphAM
#' @importFrom GOfuncR get_child_nodes
NULL

theme_legend_position <- function(position) {
  theme(legend.position=position, legend.justification=position)
}

plotNCellRegression <- function(n, n.total, x.lab="Number of cells", y.lab="N", legend.position="right", label=TRUE, size=5, palette=NULL) {
  p.df <- data.frame(N=n) %>% tibble::as_tibble(rownames="Type") %>%
    mutate(NCells=n.total[Type])

  gg <- ggplot(p.df, aes(x=NCells, y=N)) +
    geom_point(aes(color=Type)) +
    scale_x_log10() +
    ylim(0, max(p.df$N)) +
    labs(x=x.lab, y=y.lab) +
    theme_bw()

  if(label) {
    gg <- gg +
      ggrepel::geom_label_repel(aes(label=Type), size=size, min.segment.length=0.1, box.padding=0, label.size=0, max.iter=300, fill=alpha("white", 0.4))
  }

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
#' @param color.range Range for filling colors
#' @return A ggplot2 object
#' @export
plotHeatmap <- function(df, color.per.group=NULL, row.order=NULL, col.order=F, legend.position="right",
                        legend.key.width=unit(8, "pt"), legend.title="-log10(p-value)", x.axis.position="top",
                        color.range=NULL) {
  if (is.null(color.range)) {
    color.range <- c(min(0, min(df)), max(df))
  }

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

  df$p.value %<>% pmax(color.range[1]) %>% pmin(color.range[2])
  gg <- ggplot(df) + geom_tile(aes(x=Group, y=Pathway, fill=p.value), colour = "grey50") +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, color=color.per.group),
          axis.text=element_text(size=8), axis.ticks=element_blank(), axis.title=element_blank()) +
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

  if(show.significance) gg <- gg + ggpubr::stat_compare_means(aes(group = group), label = "p.signif")  # willcox test

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
    scale_y_continuous(limits=c(0, (max(df.melt$value) + 5)))

  if(show.significance) gg <- gg + ggpubr::stat_compare_means(aes(group = group), label = "p.signif")  # willcox test

  if(!is.null(palette)) gg <- gg + scale_color_manual(values=palette)

  return(gg)
}

#' Get Gene Scale
#' @param genes type of genes ("up", "down" or "all")
#' @param type type of scale ("fill" or "color")
#' @param high color for the highest value
#' @return ggplot2 fill or color scale
#' @export
getGeneScale <- function(genes=c("up", "down", "all"), type=c("fill", "color"), high="gray80", ...) {
  genes <- match.arg(genes)
  type <- match.arg(type)

  if (genes == "up") {
    low <- "red"
  } else if (genes == "down") {
    low <- "blue"
  } else {
    low <- "green"
  }

  if (type == "fill")
    return(scale_fill_gradient(low=low, high=high, ...))

  return(scale_color_gradient(low=low, high=high, ...))
}

prepareOntologyPlotDF <- function(ont.res, p.adj, n, log.colors) {
  ont.res$GeneRatio %<>% sapply(function(s) strsplit(s, "/")) %>%
    sapply(function(x) as.numeric(x[1])/as.numeric(x[2]) * 100)

  ont.res %<>% arrange(p.adjust) %>%
    filter(p.adjust <= p.adj) %>%
    {if(nrow(.) > n) .[1:n,] else .} %>%
    mutate(Description=as.factor(Description))

  if (log.colors) {
    ont.res$p.adjust %<>% log10()
  }

  return(ont.res)
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
    if (length(setdiff(names(cell.groups), names(sample.per.cell))) > 0)
      warning("Cell names in 'cell.groups' and 'sample.per.cell' are not identical, plotting intersect.")

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

plotOntologyFamily <- function(fam, data, plot.type = NULL, show.ids = F, string.length=18, legend.label.size = 1, legend.position = "topright", verbose = T, n.cores = 1) {
  parent.ids <- sapply(fam, function(x) data[[x]]$parent_go_id) %>% unlist() %>% unique()
  parent.names <- sapply(fam, function(x) data[[x]]$parent_name) %>% unlist() %>% unique()
  nodes <- data.frame(label = c(fam, parent.ids)) %>%
    mutate(., name = c(sapply(fam, function(x) data[[x]]$Description) %>% unlist(), parent.names))

  # Define edges
  edges <- lapply(nodes$label, function(y) {
    data.frame(from=y, to=nodes$label[nodes$label %in% (get_child_nodes(y)$child_go_id)])
  }) %>%
    bind_rows() %>%
    as.data.frame() %>%
    .[apply(., 1, function(x) x[1] != x[2]),] # Remove selves

  # Remove redundant inheritance
  edges %<>%
    pull(to) %>%
    unique() %>%
    plapply(function(x) {
      tmp.to <- edges[edges$to == x,]

      if(class(tmp.to) == "data.frame") {
        if(nrow(tmp.to) > 1) {
          tmp.children <- sapply(tmp.to$from, function(parent.id) {
            get_child_nodes(parent.id)$child_go_id
          })

          idx <- sapply(1:(tmp.children %>% length()), function(id) {
            any(tmp.children[-id] %>% names() %in% tmp.children[[id]])
          })

          res <- tmp.to[!idx,]
        } else if(nrow(tmp.to) == 1) {
          res <- tmp.to
        }
      } else {
        res <- NULL
      }
      return(res)
    }, progress = verbose, n.cores = n.cores) %>%
    .[!sapply(., is.null)] %>%
    bind_rows()

  # Minimize tree depending on plot.type
  if(plot.type == "dense") { # Plots all significanct terms and their 1st order relationships
    sig.parents <- c(fam, sapply(fam, function(x) data[[x]]$parents_in_IDs) %>% unlist())
    edges %<>% .[apply(., 1, function(r) any(unlist(r) %in% sig.parents)),]
  } else if(plot.type == "minimal") { # Plot significant terms only
    sig.parents <- c(fam, sapply(fam, function(x) data[[x]]$parents_in_IDs) %>% unlist())
    edges %<>% .[apply(., 1, function(r) all(unlist(r) %in% sig.parents)),]
  }

  # Convert IDs to names
  if(show.ids == F) {
    for(id in 1:nrow(nodes)) {
      edges[edges == nodes$label[id]] <- nodes$name[id]
    }
  }

  # Wrap strings for readability
  edges.wrapped <- edges %>% apply(2, function(x) wrap_strings(x, string.length))

  # Render graph
  p <- graphAM(get.adjacency(graph_from_data_frame(edges.wrapped)) %>% as.matrix(),
                      edgemode = 'directed') %>%
    layoutGraph()

  # Define layout
  ## Extract significance
  tmp.dat <- data.frame(id = c(fam, sapply(fam, function(x) data[[x]]$parents_in_IDs) %>% unlist()),
                        sig = c(sapply(fam, function(x) data[[x]]$Significance) %>% unlist(),
                                sapply(fam, function(x) data[[x]]$parents_in_IDs %>% names()) %>% unlist())) %>%
    .[match(unique(.$id), .$id),]

  ## Convert IDs to names
  for(id in tmp.dat$id) {
    tmp.dat[tmp.dat == id] <- nodes$name[nodes$label == id]
  }

  ## Wrap names
  tmp.dat$id %<>% wrap_strings(string.length)

  ## Color by significance
  node.color = tmp.dat$sig %>%
    as.numeric() %>%
    sapply(function(p) {
      if(p <= 0.05 & p > 0.01) {
        return("mistyrose1")
      } else if(p <= 0.01 & p > 0.001) {
        return("lightpink1")
      } else {
        return("indianred2")
      }
    }) %>%
    setNames(tmp.dat$id)

  ## Assign changes
  node.names <- p@renderInfo@nodes$fill %>% names()
  name.dif <- setdiff(node.names, names(node.color))
  node.color %<>% c(rep("transparent", length(name.dif)) %>% setNames(name.dif)) %>% .[match(node.names, names(.))]
  p@renderInfo@nodes$fill <- node.color %>% setNames(names(p@renderInfo@nodes$fill))
  p@renderInfo@nodes$shape <- rep("box", length(p@renderInfo@nodes$shape)) %>% setNames(names(p@renderInfo@nodes$shape))

  # Plot
  renderGraph(p)
  legend(legend.position,
         legend = c("P > 0.05","P < 0.05","P < 0.01","P < 0.001"),
         fill = c("white","mistyrose1","lightpink1","indianred2"),
         cex = legend.label.size)
}


