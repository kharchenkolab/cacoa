#' @import tibble
#' @import cowplot
#' @import dplyr
#' @import magrittr
#' @importFrom reshape2 melt
#' @importFrom sccore checkPackageInstalled
NULL

#' Helper function for creating color palettes
#' Syncs input to ggplot2::theme() params legend.position (i.e. the position of the legends)
#' and legend.justification (i.e. the anchor point for positioning legend inside the plot)
#'
#' @param position character string to pass into ggplot2::theme(legend.position=position, legend.justification=position)
#' @return returns ggplot2 theme() such that ggplot2::theme(legend.position=position, legend.justification=position)
#' @export
theme_legend_position <- function(position) {
  ggplot2::theme(legend.position=position, legend.justification=position)
}


#' Helper function for creating color palettes
#'
#' @param name character A palette name; please refer to RColorBrewer::brewer.pal()
#' @param n integer (default=NULL) Number of different colors in the palette, minimum 3, maximum depending on palette. Please refer to RColorBrewer::brewer.pal() for more details
#' @param rev boolean (default=TRUE) Whether to reverse the palette order
#' @return palette function
#'
#' @export
brewerPalette <- function(name, n=NULL, rev=TRUE) {
  checkPackageInstalled("RColorBrewer", cran=TRUE)
  if (is.null(n)) n <- RColorBrewer::brewer.pal.info[name,]$maxcolors
  pal <- RColorBrewer::brewer.pal(n, name)
  if (rev) pal <- rev(pal)
  return(grDevices::colorRampPalette(pal))
}

dark.red.palette <- colorRampPalette(c("gray95", "red", "#5A0000"), space="Lab")

#' Plot number of cells regression
#'
#' @param n the regressed variable
#' @param n.total integer Number of cells
#' @param x.lab character vector (default="Number of cells") x axis label
#' @param y.lab character vector (default="N") y axis label
#' @param legend.position character vector (default="right") legend position
#' @param label boolean (default=TRUE) Whether to show cell type labels
#' @param size integer (default=4) Label text size
#' @param palette color palette (default=NULL)
#' @param plot.line boolean (default=TRUE) Whether to plot the robust regression
#' @param line.width numeric (default=0.5) Regression line width
#' @return ggplot2 plot of regression
#'
#' @keywords internal
plotNCellRegression <- function(n, n.total, x.lab="Number of cells", y.lab="N", legend.position="right", label=TRUE, size=4,
                                palette=NULL, plot.line=TRUE, line.width=0.5, plot.theme=ggplot2::theme_get()) {

  p.df <- data.frame(N=n) %>% tibble::as_tibble(rownames="Type") %>%
    mutate(NCells=n.total[Type])

  gg <- ggplot(p.df, aes(x=NCells, y=N, color = Type)) +
    geom_point() +
    scale_x_log10() +
    ylim(0, max(p.df$N)) +
    labs(x=x.lab, y=y.lab)

  if(label) {
    gg <- gg +
      ggrepel::geom_label_repel(aes(label=Type), size=size, min.segment.length=0.1, box.padding=0, label.size=0, max.iter=300, fill=NA)
  }

  gg <- gg +
    plot.theme +
    theme(legend.background=element_rect(fill=alpha("white", 0.4))) +
    theme_legend_position(legend.position) +
    guides(color=guide_legend(title="Cell type"))

  if(!is.null(palette)) {
    gg <- gg+ scale_color_manual(values=palette)
  }

  if (plot.line) {
    gg <- gg +
      geom_smooth(method=MASS::rlm, formula = y~x, se=FALSE, color="gray", size=line.width, linetype=2)
  }

  return(gg)
}

#' Plot count boxplots per type
#'
#' @param count.df data.frame with columns `group`, `variable` and `value`
#' @param notch Whether to show notch in the boxplots
#' @param alpha Transparency level on the data points (default: 0.2)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param size marker size (default: 0.5)
#' @param jitter.width width of the point jitter (default: 0.15)
#' @keywords internal
plotCountBoxplotsPerType <- function(count.df, y.lab="count", x.lab="", y.expand=0.05, show.significance=FALSE,
                                     jitter.width=0.15, notch=FALSE, legend.position="right", alpha=0.2, size=0.5,
                                     palette=NULL, adjust.pvalues=TRUE, plot.theme=theme_get(), pvalue.y=NULL,
                                     ns.symbol="", p.adjust.method='BH') {
  gg <- ggplot(count.df, aes(x=variable, y=value, by=group, fill=group)) +
    geom_boxplot(position=position_dodge(), outlier.shape = NA, notch=notch) +
    labs(x=x.lab, y=y.lab) +
    plot.theme +
    theme_legend_position(legend.position) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_blank()) +
    geom_point(position=position_jitterdodge(jitter.width=jitter.width), color="black", size=size, alpha=alpha) +
    scale_y_continuous(expand=c(0, 0, y.expand, 0), limits=c(0, max(count.df$value)))

  if (show.significance) {
    suppressWarnings(
      pval.df <- count.df %>% group_by(variable) %>%
        summarise(pvalue=wilcox.test(value[group == group[1]], value[group != group[1]])$p.value)
    )

    if (adjust.pvalues) {
      pval.df$pvalue %<>% p.adjust(p.adjust.method)
    }

    if (is.null(pvalue.y)) pvalue.y <- max(count.df$value)

    pval.df$pvalue %<>% pvalueToCode(ns.symbol=ns.symbol)
    gg <- gg + geom_text(data=pval.df, mapping=aes(label=pvalue, by=NULL, fill=NULL), y=pvalue.y, color="black")
  }

  if (!is.null(palette)) gg <- gg + scale_fill_manual(values=palette)
  return(gg)
}

#' Plot heatmap
#'
#' @param df Data frame with plot data
#' @param color.per.group Colors per cell group (default=NULL)
#' @param row.order Forced row order (default=NULL)
#' @param col.order Forced column order (default=NULL)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param size.df (default=NULL)
#' @param size.range (default=c(1, 5))
#' @param size.legend.title (default="size")
#' @param legend.key.width (default=unit(8, "pt))
#' @param legend.title Title on plot (default="-log10(p-value)")
#' @param x.axis.position Position of x axis (default="top")
#' @param color.range Range for filling colors
#' @param plot.theme (default=ggplot2::theme_get())
#' @param symmetric boolean (default=FALSE)
#' @param palette (default=NULL)
#' @param font.size integer (default=8)
#' @param distance character string (default="manhattan")
#' @param clust.method character string (default="complete")
#' @param grid.color Color of the grid. Set to "transparent" to disable the grid. (default: "gray50")
#' @return A ggplot2 object
#'
#' @export
plotHeatmap <- function(df, color.per.group=NULL, row.order=TRUE, col.order=TRUE, legend.position="right",
                        size.df=NULL, size.range=c(1, 5), size.legend.title="size",
                        legend.key.width=unit(8, "pt"), legend.title="-log10(p-value)", x.axis.position="top",
                        color.range=NULL, plot.theme=ggplot2::theme_get(), symmetric=FALSE, palette=NULL, font.size=8,
                        distance="manhattan", clust.method="complete", grid.color="gray50") {
  if (is.null(color.range)) {
    if (prod(range(df, na.rm=TRUE)) < 0) {
      color.range <- c(-1, 1) * max(abs(df), na.rm=TRUE)
    } else {
      color.range <- c(min(0, min(df, na.rm=TRUE)), max(df, na.rm=TRUE))
    }
  }

  if (is.logical(row.order)) {
    if (row.order && (nrow(df) > 2)) {
      row.order <- rownames(df)[dist(df, method=distance) %>% hclust(method=clust.method) %>% .$order] %>% rev()
    } else {
      row.order <- rownames(df)
    }
  }

  if (is.logical(col.order)) {
    if (col.order && (ncol(df) > 2)) {
      col.order <- colnames(df)[dist(t(df), method=distance) %>% hclust(method=clust.method) %>% .$order] %>% rev()
    } else {
      col.order <- colnames(df)
    }
  }

  df %<>% tibble::as_tibble(rownames="G1") %>%
    reshape2::melt(id.vars="G1", variable.name="G2", value.name="value")

  df %<>% dplyr::mutate(G1=factor(G1, levels=row.order))
  df %<>% dplyr::mutate(G2=factor(G2, levels=col.order))

  if (is.null(color.per.group)) {
    color.per.group <- "black"
  } else {
    color.per.group <- color.per.group[levels(df$G2)]
  }

  df$value %<>% pmax(color.range[1]) %>% pmin(color.range[2])

  color.guide <- guide_colorbar(
    title=legend.title, title.position="left", title.theme=element_text(angle=90, hjust=0.5)
  )
  if (is.null(size.df)) {
    gg <- ggplot(df) + geom_tile(aes(x=G2, y=G1, fill=value), color=grid.color) + guides(fill=color.guide)
    return.fill <- TRUE
    expand <- expansion(0, 0)
  } else {
    df$size <- reshape2::melt(size.df, id.vars=NULL)$value
    size.guide <- guide_legend(
      title=size.legend.title, title.position="left", title.theme=element_text(angle=90, hjust=0.5)
    )

    gg <- ggplot(df) +
      geom_point(aes(x=G2, y=G1, color=value, size=size)) +
      scale_size_continuous(range=size.range) +
      guides(color=color.guide, size=size.guide)
    return.fill <- FALSE
    expand <- expansion(mult=0, add=0.5)
  }

  gg <- gg + plot.theme +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, color=color.per.group),
          axis.text=element_text(size=font.size), axis.ticks=element_blank(), axis.title=element_blank()) +
    scale_y_discrete(position="right", expand=expand) +
    scale_x_discrete(expand=expand, position=x.axis.position) +
    theme_legend_position(legend.position) +
    theme(legend.key.width=legend.key.width, legend.background=element_blank()) +
    val2ggcol(df$value, palette=palette, color.range=color.range, return.fill=return.fill)

  if (symmetric) {
    gg <- gg + theme(axis.text.y=element_text(color=color.per.group))
  }

  return(gg)
}

#' Get Gene Scale
#'
#' @param genes character vector Type of genes (default=c("up", "down" or "all"))
#' @param neutral.col character string (default="white")
#' @param bidirectional boolean (default=FALSE)
#' @return palette function
#'
#' @export
getGenePalette <- function(genes=c("up", "down", "all"), neutral.col="white", bidirectional=FALSE) {
  genes <- match.arg(genes)
  up.cols <- c("#d6604d", "#67001f")
  down.cols <- c("#4393c3", "#053061")
  all.cols <- c("#5aae61", "#00441b")

  if (bidirectional){
    return(colorRampPalette(c(rev(down.cols), neutral.col, up.cols)))
  }

  if (genes == "up") {
    sign.cols <- up.cols
  } else if (genes == "down") {
    sign.cols <- down.cols
  } else {
    sign.cols <- all.cols
  }

  return(colorRampPalette(c(neutral.col, sign.cols)))
}

#' @keywords internal
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

#' @keywords internal
getOntologyPlotTitle <- function(genes, cell.subgroup, type) {
  if(genes == "all"){
    return(ggtitle(paste(cell.subgroup, type, "terms, all DE genes")))
  }

  return(ggtitle(paste0(cell.subgroup, " ", type, " terms, ", genes,"-regulated DE genes")))
}

#' @keywords internal
estimateMeanCI <- function(arr, quant=0.05, n.samples=500, ...) {
  if (length(arr) == 1) return(c(NA, NA))
  s.means <- sapply(1:n.samples, function(i) mean(sample(arr, replace=TRUE), ...))
  return(quantile(s.means, c(quant, 1 - quant)))
}

#' Plot bar, point or boxplots showing mean/median values per cell type.
#' This is a generic function for plotting mean or median values per cell type
#' (used for expression shift distances and others).
#'
#' @param df dataframe containing the results, including $value and $Type slots which will be summarized
#' @param type character vector (default=c('box', 'point', 'bar')) Type of a plot "bar" (default), "point" (mean + sd), or "box" for boxplot
#' @param show.jitter boolean (default = FALSE) Whether to show individual data points
#' @param jitter.alpha transparency value for the data points (default: 0.05)
#' @param notch boolean Whether to show notches in the boxplot version (default=TRUE)
#' @param palette - cell type palette
#' @param coord.flip flip coordinates of the plot
#' @return A ggplot2 object
#' @keywords internal
plotMeanMedValuesPerCellType <- function(df, pvalues=NULL, type=c('box', 'point', 'bar'), show.jitter=TRUE,
                                         notch=TRUE, jitter.alpha=0.05, palette=NULL, ylab='expression distance',
                                         yline=1, plot.theme=theme_get(), jitter.size=1, line.size=0.75, trim=0,
                                         order.x=TRUE, pvalue.y=NULL, y.max=NULL, y.offset=NULL, ns.symbol="",
                                         coord.flip=FALSE) {
  type <- match.arg(type)
  df$Type %<>% as.factor()
  if (!is.null(y.offset)) {
    df$value <- df$value + y.offset
  }

  # calculate mean, se and median
  odf <- df <- na.omit(df); # full df is now in odf

  if (!is.null(y.max)) {
    odf$value %<>% pmin(y.max)
  }

  if (is.null(pvalue.y)) pvalue.y <- max(odf$value)

  conf.ints <- odf %$% split(value, Type) %>% lapply(estimateMeanCI, trim=trim)
  # calculate mean and CI
  df %<>% group_by(Type) %>%
    summarise(mean=mean(value, trim=trim), med=median(value)) %>%
    mutate(Type=as.character(Type)) %>%
    mutate(LI=as.numeric(sapply(Type, function(ct) conf.ints[[ct]][1])),
           UI=as.numeric(sapply(Type, function(ct) conf.ints[[ct]][2]))) %>%
    .[!is.na(.$mean),] %>%
    mutate(Type=factor(Type, levels=levels(odf$Type)))

  if (order.x) {
    if (type == 'box') {
      df %<>% arrange(med)
    } else {
      df %<>% arrange(mean)
    }

    # order cell types according to the mean
    odf$Type %<>% factor(levels=df$Type)
    df$Type %<>% factor(., levels=.)
  }

  if (type=='box') { # boxplot
    p <- ggplot(odf,aes(x=Type, y=value, fill=Type)) + geom_boxplot(notch=notch, outlier.shape=NA)
  } else if (type=='point') { # point + se
    p <- ggplot(df, aes(x=Type, y=mean, color=Type)) + geom_point(size=3) +
      geom_errorbar(aes(ymin=LI, ymax=UI), width=0.2, linewidth=line.size)
    if (!is.null(palette)) {p <- p + scale_color_manual(values=palette)}
  } else { # barplot
    p <- ggplot(df,aes(x=Type,y=mean,fill=Type)) + geom_bar(stat='identity') +
      geom_errorbar(aes(ymin=LI, ymax=UI), width=0.2, linewidth=line.size)
  }
  if (!is.na(yline) && !is.null(yline)) {p <- p + geom_hline(yintercept = yline, linetype=2, color='gray50')}
  p <- p +
    plot.theme +
    theme(
      panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
      legend.position="none", axis.text=element_text(size=12)
    ) +
    guides(fill="none") + labs(y=ylab)

  if (coord.flip) {
    p <- p + coord_flip() + theme(axis.title.y=element_blank())
  } else {
    p <- p + theme(
      axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
      axis.text.y=element_text(angle=90, hjust=0.5),
      axis.title.x=element_blank()
    )
  }

  if (show.jitter) {
    p <- p +
      geom_jitter(data=odf, aes(x=Type,y=value), color=1, position=position_jitter(0.1), show.legend=FALSE,
                  alpha=jitter.alpha, size=jitter.size)
  }

  if (!is.null(pvalues)) {
    pval.df <- pvalueToCode(pvalues, ns.symbol=ns.symbol) %>%
      tibble(Type=factor(names(.), levels=levels(df$Type)), pvalue=.) %>% na.omit()

    p <- p + geom_text(data=pval.df, mapping=aes(x=Type, label=pvalue), y=pvalue.y, color="black")
  }

  if(!is.null(palette)) {
    p <- p + scale_fill_manual(values=palette)
  }

  return(p)
}

# TODO: Improve speed of this function. No need to check BP/MF/CC all the time

#' @keywords internal
reduceEdges <- function(edges, verbose=TRUE, n.cores=1) {
  edges %>%
    pull(to) %>%
    unique() %>%
    plapply(function(x) {
      tmp.to <- edges[edges$to == x,]

      if(!("data.frame" %in% class(tmp.to)) || nrow(tmp.to) == 0) return(NULL)
      if(nrow(tmp.to) == 1) return(tmp.to)

      tmp.children <- sapply(tmp.to$from, function(parent.id) {
        GOfuncR::get_child_nodes(parent.id)$child_go_id
      })

      idx <- sapply(1:length(tmp.children), function(id) {
        any(names(tmp.children[-id]) %in% tmp.children[[id]])
      })

      return(tmp.to[!idx,])
    }, progress = verbose, n.cores = n.cores) %>%
    .[!sapply(., is.null)] %>%
    bind_rows()
}

#' Plot related ontologies in one hierarchical network plot
#'
#' @param fam List of ontology IDs for the chosen family
#' @param data Data frame if raw ontology data for the chosen cell type
#' @param plot.type character (default="complete") How much of the family network should be plotted. Can be "complete" (entire network), "dense" (show 1 parent for each significant term), or "minimal" (only show significant terms)
#' @param show.ids boolean (default=FALSE) Whether to show ontology IDs instead of names
#' @param string.length integer (default=18) Length of strings for wrapping in order to fit text within boxes
#' @param legend.label.size integer (default=1) Size of legend labels
#' @param legend.position character vector (default="topright") Position of legend
#' @param verbose boolean (default=TRUE) Print messages
#' @param n.cores integer (default=1) Number of cores to use
#' @param reduce.edges boolean (default=FALSE) Remove redundant edges in network
#' @param font.size integer (default=24) Size of the font
#' @return Rgraphviz object
#' @keywords internal
plotOntologyFamily <- function(fam, data, plot.type="complete", show.ids=FALSE, string.length=18, legend.label.size=1,
                               legend.position="topright", verbose=TRUE, n.cores=1, reduce.edges=FALSE, font.size=24) {
  checkPackageInstalled("Rgraphviz", bioc=TRUE)
  # Define nodes
  parent.ids <- sapply(fam, function(x) data[[x]]$parent_go_id) %>%
    unlist() %>% unique()
  parent.names <- sapply(fam, function(x) data[[x]]$parent_name) %>%
    unlist() %>% unique()
  nodes <- data.frame(label = c(fam, parent.ids)) %>%
    mutate(., name = c(sapply(fam, function(x) data[[x]]$Description) %>% unlist(), parent.names))

  # Define edges
  edges <- lapply(nodes$label, function(y) {
    child.nodes <- GOfuncR::get_child_nodes(y) %$% child_go_id[distance == 1]
    to.nodes <- nodes$label[nodes$label %in% child.nodes]
    if (length(to.nodes) == 0) return(data.frame())
    data.frame(from=y, to=to.nodes)
  }) %>%
    bind_rows() %>% as.data.frame() %$% .[from != to,]

  # Remove redundant inheritance
  if (reduce.edges) edges %<>% reduceEdges(verbose=verbose, n.cores=n.cores)

  # Minimize tree depending on plot.type
  if(plot.type == "dense") { # Plots all significanct terms and their 1st order relationships
    sig.parents <- c(fam, sapply(fam, function(x) data[[x]]$parents_in_IDs) %>% unlist())
    edges %<>%
      .[apply(., 1, function(r) any(unlist(r) %in% sig.parents)),]
  } else if(plot.type == "minimal") { # Plot significant terms only
    sig.parents <- c(fam, sapply(fam, function(x) data[[x]]$parents_in_IDs) %>% unlist())
    edges %<>%
      .[apply(., 1, function(r) all(unlist(r) %in% sig.parents)),]
  } else if (plot.type != "complete") stop("Unknown plot type: ", plot.type)

  # Check for remaining edges
  if (nrow(edges) == 0) stop("No significant terms left after filtering. Consider changing 'plot.type' to something less stringent, e.g., 'dense' or 'complete'.")
  
  # Convert IDs to names
  if(!show.ids) {
    for(id in 1:nrow(nodes)) {
      edges[edges == nodes$label[id]] <- nodes$name[id]
    }
  }

  # Wrap strings for readability
  edges.wrapped <- edges %>%
    apply(2, function(x) wrap_strings(x, string.length))
  
  # Check for correct format of edges.wrapped
  if (!inherits(edges.wrapped, "matrix")) edges.wrapped %<>% rbind()

  # Render graph
  p <- igraph::graph_from_data_frame(edges.wrapped) %>%
    igraph::get.adjacency() %>%
    as.matrix() %>%
    graph::graphAM(edgemode = 'directed') %>%
    Rgraphviz::layoutGraph()

  # Define layout
  ## Extract significance
  tmp.dat <- data.frame(id = c(fam, sapply(fam, function(x) data[[x]]$parents_in_IDs) %>% unlist()),
                        sig = c(sapply(fam, function(x) data[[x]]$Significance) %>% unlist(),
                                sapply(fam, function(x) data[[x]]$parents_in_IDs %>% names()) %>% unlist())) %>%
    .[match(unique(.$id), .$id),]

  ## Convert IDs to names
  # TODO: Dependent on show.ids?
  nodes %<>% .[match(unique(.$label), .$label),] # Must remove doublets

  for(id in tmp.dat$id) {
    tmp.dat[tmp.dat == id] <- nodes$name[nodes$label == id]
  }

  ## Wrap names
  tmp.dat$id %<>% wrap_strings(string.length)

  ## Color by significance
  node.color <- as.numeric(tmp.dat$sig) %>% sapply(function(p) {
    if(p <= 0.05 & p > 0.01) {
      return("mistyrose1")
    } else if (p <= 0.01 & p > 0.001) {
      return("lightpink1")
    } else {
      return("indianred2")
    }
  }) %>% setNames(tmp.dat$id)

  ## Assign changes
  node.names <- p@renderInfo@nodes$fill %>%
    names()
  name.dif <- setdiff(node.names, names(node.color))
  node.color %<>%
    c(rep("transparent", length(name.dif)) %>% setNames(name.dif)) %>%
    .[match(node.names, names(.))]
  p@renderInfo@nodes$fill <- node.color %>%
    setNames(names(p@renderInfo@nodes$fill))
  p@renderInfo@nodes$shape <- rep("box", length(p@renderInfo@nodes$shape)) %>%
    setNames(names(p@renderInfo@nodes$shape))
  p@renderInfo@nodes$fontsize <- font.size

  # Plot
  Rgraphviz::renderGraph(p)
  legend(legend.position,
         legend = c("P > 0.05","P < 0.05","P < 0.01","P < 0.001"),
         fill = c("white","mistyrose1","lightpink1","indianred2"),
         cex = legend.label.size)

  return(p)
}


#' @keywords internal
plotVolcano <- function(de.df, p.name='padj', color.var = 'CellFrac', legend.pos="none", palette=brewerPalette("RdYlBu"), lf.cutoff=1.5, p.cutoff=0.05,
                        cell.frac.cutoff=0.2, size=c(0.1, 1.0), lab.size=2, draw.connectors=TRUE, sel.labels=NULL, plot.theme=theme_get(), ...) {

  # fold changes far below the cutoff (a continuous tested variable: the change per unit, or per numeric step): the x axis follows the data,
  # no fold-change lines are drawn, and the labels are the most significant genes
  lfc <- de.df$log2FoldChange[is.finite(de.df$log2FoldChange)]; small <- length(lfc) > 0 && max(abs(lfc)) < lf.cutoff / 3
  if (small) { lf.cutoff <- 0; xl <- c(-1, 1) * max(abs(lfc)) * 1.1
    if (is.null(sel.labels)) { o <- order(de.df[[p.name]]); o <- o[de.df[[p.name]][o] <= p.cutoff]; sel.labels <- utils::head(de.df$Gene[o], 15) } }
  if (!requireNamespace("EnhancedVolcano", quietly = TRUE))
    return(plotVolcanoSimple(de.df, p.name = p.name, color.var = color.var, palette = palette, lf.cutoff = lf.cutoff, p.cutoff = p.cutoff,
                             cell.frac.cutoff = cell.frac.cutoff, size = size, lab.size = lab.size, sel.labels = sel.labels, plot.theme = plot.theme,
                             legend.pos = legend.pos))
  if (is.null(sel.labels) && (color.var == 'CellFrac')) {
    sel.labels <- de.df %$%
      Gene[(.[[p.name]] <= p.cutoff) & (abs(log2FoldChange) >= lf.cutoff) & (CellFrac >= cell.frac.cutoff)]
  } else if (is.null(sel.labels) && (color.var == 'Stability')) {
    sel.labels <- de.df %$%
      Gene[(.[[p.name]] <= p.cutoff) & (abs(log2FoldChange) >= lf.cutoff) & (Stability >= 0.5)]
  }
  gg <- EnhancedVolcano::EnhancedVolcano(
    de.df, lab=de.df$Gene, x='log2FoldChange', y=p.name, arrowheads=FALSE,
    pCutoff=p.cutoff, FCcutoff=lf.cutoff, title=NULL, subtitle=NULL, caption=NULL,
    selectLab=sel.labels, drawConnectors=draw.connectors, labSize=lab.size, ...
  )
  if (small) { gg$layers <- gg$layers[!vapply(gg$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))]; gg <- gg + ggplot2::coord_cartesian(xlim = xl) }

  point.id <- sapply(gg$layers, function(l) "GeomPoint" %in% class(l$geom)) %>%
    which() %>% .[1]
  if(color.var == 'CellFrac') {
    gg$layers[[point.id]] <- geom_point(aes(x=log2FoldChange, y=-log10(.data[[p.name]]), color=CellFrac, size=CellFrac))
    gg$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
    gg <- gg + val2ggcol(de.df$CellFrac, palette=palette, color.range=c(0, 1))
  } else if (color.var == 'Stability') {
    gg$layers[[point.id]] <- geom_point(aes(x=log2FoldChange, y=-log10(.data[[p.name]]), color=Stability, size=Stability))
    gg$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
    gg <- gg + val2ggcol(de.df$Stability, palette=palette, color.range=c(0, 1))
  }

  gg <- gg +
    scale_size_continuous(range=size, name="Expr. frac", limits=c(0, 1)) +
    guides(color=guide_colorbar(title="Expr. frac")) +
    plot.theme +
    theme_legend_position(legend.pos)
  return(gg)
}

#' @keywords internal
parseLimitRange <- function(lims, vals) {
  if (is.null(lims)) return(range(vals))
  if (!is.character(lims)) return(lims)
  lims[grep("%", lims)] %<>%
    sapply(function(q) {as.numeric(strsplit(q, "%")[[1]]) / 100}) %>%
    quantile(vals, ., na.rm = TRUE)
  return(as.numeric(lims))
}

#' Symmetric colour range for signed scores
#'
#' @param x the scores
#' @param range a quantile of `abs(x)` as a string (e.g. `"99.5%"`), or a numeric vector whose largest absolute value sets
#'   the range (`NULL`: the full range)
#' @param floor the smallest half-range (e.g. 1.96 for adjusted z-scores, so that a map without significant cells looks pale)
#' @param cap the largest half-range from a quantile or `NULL` range (a numeric `range` is taken as given); adjusted
#'   z-scores of q-values near 0 run to 37, which would wash out the rest of the map
#' @return `c(-lim, lim)`
#' @keywords internal
symmetricColorRange <- function(x, range = "99.5%", floor = 0, cap = Inf) {
  ax <- abs(x[is.finite(x)])
  lim <- if (is.null(range)) max(ax, 0) else if (is.character(range)) {
    q <- as.numeric(sub("%", "", range[grep("%", range)])) / 100
    if (!length(q) || any(!is.finite(q))) stop("range must be a quantile string such as '99.5%' or numeric")
    if (length(ax)) max(stats::quantile(ax, probs = max(q), names = FALSE)) else 0
  } else max(abs(range), na.rm = TRUE)
  lim <- max(min(lim, cap), floor, 1e-10)
  c(-lim, lim)
}

#' @keywords internal
prepareGeneExpressionComparisonPlotInfo <- function(de.info, genes, plots, smoothed, max.z, max.z.adj, max.lfc, z.palette, z.adj.palette, lfc.palette) {
  z.scores <- NULL
  z.adj <- NULL
  if (any(plots != "expression")) {
    z.scores <- de.info$z
    z.adj <- de.info$z.adj
    if (is.null(de.info)) {
      warning("Z-scores were not estimated. See estimateClusterFreeDE().")
      plots <- "expression"
    }

    if (!all(genes %in% colnames(de.info$z))) {
      missed.genes <- setdiff(genes, colnames(de.info$z))
      warning("Z-scores for genes ", paste(missed.genes, collapse=', '), " are not estimated. See estimateClusterFreeDE().")
      plots <- "expression"
    }

    if (is.null(max.z.adj)) {
      max.z.adj <- z.adj@x %>% c(1e-5) %>% range(z.adj@x, na.rm=TRUE) %>% abs() %>% max() %>% min(5)
    }

    z.scores@x %<>% pmin(max.z) %>% pmax(-max.z)
    z.adj@x %<>% pmin(max.z.adj) %>% pmax(-max.z.adj)
  }

  if (("z" %in% plots) && smoothed) {
    if ((is.null(de.info$z.smoothed) || !all(genes %in% colnames(de.info$z.smoothed)))){
      missed.genes <- setdiff(colnames(de.info$z.smoothed), genes)
      warning("Smoothed Z-scores for genes ", paste(missed.genes, collapse=', '), " are not estimated. See smoothClusterFreeZScores().")
    } else {
      z.scores <- de.info$z.smoothed
    }
  }

  if (("z.adj" %in% plots) && smoothed) {
    if ((is.null(de.info$z.adj.smoothed) || !all(genes %in% colnames(de.info$z.adj.smoothed)))){
      missed.genes <- setdiff(colnames(de.info$z.smoothed), genes)
      warning("Smoothed Z-scores for genes ", paste(missed.genes, collapse=', '), " are not estimated. See smoothClusterFreeZScores().")
    } else {
      z.adj <- de.info$z.adj.smoothed
    }
  }

  plot.parts <- list(
    z.adj=list(title="Z, adjusted", leg.title="Z,adj", max=max.z.adj, scores=z.adj, palette=z.adj.palette),
    z=list(title="Z-score", leg.title="Z", max=max.z, scores=z.scores, palette=z.palette),
    lfc=list(title="Log2(fold-change)", leg.title="LFC", max=max.lfc, scores=de.info$lfc, palette=lfc.palette)
  )[plots[plots != "expression"]]

  return(plot.parts)
}

#' @keywords internal
pvalueToCode <- function(pvals, ns.symbol="ns") {
  symnum(pvals, corr=FALSE, na=FALSE, legend=FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
         symbols = c("***", "**", "*", ns.symbol)) %>%
    as.character() %>% setNames(names(pvals))
}

#' @keywords internal
transferLabelLayer <- function(gg.target, gg.source, font.size) {
  ls <- gg.source$layers %>% .[sapply(., function(l) "GeomLabelRepel" %in% class(l$geom))]
  if (length(ls) != 1) {
    warning("Can't find annotation layer\n")
    return(gg.target)
  }

  gg.target <- gg.target + ls[[1]] +
    scale_size_continuous(range=font.size, trans='identity', guide='none')

  return(gg.target)
}

#' @keywords internal
getScaledZGradient <- function(min.z, palette, color.range) {
  if (all(abs(color.range) < 1e-10)) return(scale_color_gradientn(colors="grey80", limits=c(-1e-10, 1e-10), guide="none"))
  if (length(color.range) == 1) {
    if (min.z > (color.range - 1e-10))
      return(scale_color_gradientn(colors="grey80", limits=c(0, max(color.range, 1e-10)), guide="none"))

    col.vals <- c(0, seq(min.z, color.range, length.out=20))
    color.range <- c(0, color.range)
    pal.div <- palette; palette <- function(n) pal.div(2 * n - 1)[n:(2 * n - 1)]   # one-sided: the upper half of the diverging palette
  } else {
    col.vals <- c(seq(color.range[1], -min.z, length.out=10), 0, seq(min.z, color.range[2], length.out=10))
  }
  col.vals %<>% scales::rescale()
  scale <- scale_color_gradientn(colors=palette(21), values=col.vals, limits=color.range)
  return(scale)
}


#' @keywords internal
plotSampleDistanceMatrix <- function(p.dists, sample.labels=NULL, n.cells.per.samp, method='MDS', sample.colors=NULL,
                                     show.sample.size=TRUE, palette=NULL, font.size=NULL, show.ticks=FALSE, title=NULL,
                                     shape.labels=NULL, shape.title=NULL, show.labels=FALSE, size=5, color.title=NULL, 
                                     cont.palette = rev(RColorBrewer::brewer.pal(11, "Spectral")),
                                     perplexity=4, max.iter=1e3, plot.theme=theme_get(), n.neighbors=15, ...) {
      if (method == 'tSNE') {
        checkPackageInstalled('Rtsne', cran=TRUE, details='for `method="tSNE"`')
        emb <- Rtsne::Rtsne(p.dists, is_distance=TRUE, perplexity=perplexity, max_iter=max.iter)$Y
      } else if (method == "UMAP") {
        checkPackageInstalled('uwot', cran=TRUE, details='for `method="UMAP"`')
        emb <- estimateUMAPOnDistances(p.dists, n.neighbors=n.neighbors, n_epochs=max.iter)
      } else if (method == 'MDS') {
        # classical MDS of the geometry the models use: "cor" and "l1" distances are already squared distances there
        # (asSquaredDistance), so their square root is embedded; "l2" (and an untyped matrix) as is
        dt <- attr(p.dists, "dist.type")
        Dm <- if (!is.null(dt) && dt %in% c("cor", "l1")) sqrt(pmax(as.matrix(p.dists), 0)) else as.matrix(p.dists)
        emb <- suppressWarnings(cmdscale(Dm, eig=TRUE, k=2)$points)
        if (is.null(dim(emb)) || ncol(emb) < 2) {   # fewer than two positive eigenvalues: pad with zero coordinates
          if (!any(abs(Dm) > 1e-12)) warning("the samples do not differ in this space: all distances are zero")
          emb <- cbind(matrix(emb, nrow(Dm)), matrix(0, nrow(Dm), 2))[, 1:2, drop = FALSE]
        }
      } else if (method == 'heatmap') {
        color.per.group <- NULL
        if (!is.null(palette)) {
          color.per.group <- sample.groups %>% {setNames(palette[as.character(.)], names(.))}
        }
        gg <- plotHeatmap(p.dists, color.per.group=color.per.group, legend.title="Distance", symmetric=TRUE)
        return(gg)
      } else {
        stop("unknown embedding method")
      }

      df <- data.frame(emb) %>% set_rownames(rownames(p.dists)) %>% set_colnames(c("x", "y")) %>%
               mutate(sample = rownames(.), n.cells = as.vector(n.cells.per.samp[sample]))

      if (!is.null(sample.labels)) {
        # align by sample id
        df$label <- sample.labels[df$sample]
       } else {
        df$label <- factor("All")
       }

      if (!is.null(shape.labels)) {
        df$shape <- shape.labels[df$sample]
        } else {
          df$shape <- factor("All")
          }
          # shapes must be discrete; if numeric, coerce with a warning
          if (is.numeric(df$shape) || is.integer(df$shape)) {
            warning("shape.by is numeric; coercing to factor for discrete shapes.")
            df$shape <- as.factor(df$shape)
            } else {
              df$shape <- addNA(as.factor(df$shape))
              }
      
      if (is.null(color.title)) color.title <- "Covariate"
      if (is.null(shape.title)) shape.title <- "Covariate"

      if (!(is.numeric(df$label) || is.integer(df$label))) {
        df$label <- addNA(as.factor(df$label))
        }

      is.cont <- is.numeric(df$label) || is.integer(df$label)
      if (is.cont) {
        if (is.null(sample.colors)) {
          gg <- ggplot(df, aes(x, y, color=label, shape=shape))
          } else {
            df$color <- sample.colors[as.character(df$sample)]
            gg <- ggplot(df, aes(x, y, color=color, shape=shape))
            }
            } else {
              if (is.null(sample.colors)) {
                gg <- ggplot(df, aes(x, y, color=label, shape=shape))
                } else {
                  df$color <- sample.colors[as.character(df$sample)]
                  gg <- ggplot(df, aes(x, y, color=color, shape=shape))
                  }
                }
      # apply legend titles
      gg <- gg + labs(color = color.title, shape = shape.title)

      # Apply palette only when it is valid for the data type
      if (is.cont) {
        if (is.null(palette)) {
          # default continuous palette
          gg <- gg + scale_color_gradientn(colors = cont.palette)
        } else if (is.function(palette)) {
          gg <- gg + scale_color_gradientn(colors = palette(100))
        } else if (is.character(palette) && length(palette) > 1) {
          gg <- gg + scale_color_gradientn(colors = palette)
        } else {
          warning("Ignoring 'palette' for continuous covariate: provide a function like viridisLite::viridis or a vector of colors.")
        }
      } else {
        if (!is.null(palette)) {
          if (is.function(palette)) {
            gg <- gg + scale_color_manual(values = palette(length(levels(df$label))))
          } else if (is.character(palette) && length(palette) > 0) {
            nlev <- length(levels(df$label))
            if (is.null(names(palette)) && length(palette) < nlev) {
              warning("Ignoring 'palette': insufficient values (need ", nlev, ", got ", length(palette), ").")
            } else {
              gg <- gg + scale_color_manual(values = palette)
            }
          } else {
            warning("Ignoring 'palette': must be a non-empty character vector (discrete) or a function (continuous).")
          }
        }
      }

      gg <- gg + plot.theme
      # nothing mapped to colour / shape: one neutral colour, no legend entries
      if (is.null(sample.labels) && is.null(sample.colors)) gg <- gg + (if (is.null(palette)) scale_color_manual(values = "grey30", guide = "none") else guides(color = "none"))
      if (is.null(shape.labels)) gg <- gg + guides(shape = "none") else {
        n.sh <- nlevels(df$shape)                                      # solid shapes first, then the outlined ones
        if (n.sh <= 9) gg <- gg + scale_shape_manual(values = c(16, 17, 15, 18, 8, 4, 3, 1, 2)[seq_len(n.sh)], na.value = 1)
      }

      if (show.sample.size) {
        if (length(size) == 1) {
          size <- c(size * 0.5, size * 1.5)
        }
        gg <- gg + geom_point(aes(size=n.cells)) +
          scale_size_continuous(trans="log10", range=size, name="Num. cells")
        # readable legend keys (the point size of the keys would otherwise follow the smallest cell count)
        if (!is.null(sample.labels) && is.null(sample.colors) && !is.cont) gg <- gg + guides(color = guide_legend(override.aes = list(size = 4.5)))
        if (!is.null(shape.labels)) gg <- gg + guides(shape = guide_legend(override.aes = list(size = 4.5)))
      } else {
        gg <- gg + geom_point(size=size)
      }

      gg %<>% sccore::styleEmbeddingPlot(title=title, show.ticks=show.ticks, show.labels=show.labels, ...)
      gg <- gg + ggtitle(title) + theme(plot.title = element_text(face = "bold"))
      if (!is.null(font.size)) {
        gg <- gg + ggrepel::geom_text_repel(aes(label=sample), size=font.size, color="black")
      }

      return(gg)
    }

#' Performs ontology clustering by genes and then shows medoids of the clusters with
#'   enrichplot::dotplot
#'
#' @param ont.res something
#' @param p.adj numeric Adjusted p-value set (default=0.05)
#' @param min.genes integer Minimal number of genes (default=1)
#' @param cut.h numeric (default=0.66)
#' @param top.n integer Number of enriched terms to display (default=Inf). If Inf, no limit is set.
#' @param ... additional parameters passed to enrichplot::dotplot()
#' @return enrichplot::dotplot showing medoids of gene clustering
#' @export
clusteredOntologyDotplot <- function(ont.res, p.adj=0.05, min.genes=1, cut.h=0.66, top.n=Inf, ...) {
  checkPackageInstalled("enrichplot", bioc=TRUE)
  clusts <- ont.res@result %>% filter(p.adjust < p.adj) %$%
    setNames(strsplit(geneID, "/"), Description) %>% .[sapply(., length) >= min.genes] %>%
    estimateClusterPerGO(cut.h=cut.h) %>% {split(names(.), .)} %>%
    sapply(estimateOntologyClusterName, method='medoid')

  ont.res@result %<>% filter(Description %in% clusts)

  enrichplot::dotplot(ont.res, showCategory=top.n, ...)
}

#' Internal function for estimating CDA space
#' @keywords internal
estimateCdaSpace <- function(d.counts, d.groups, thresh.pc.var = 0.95, n.dim = 2){
  cell.loadings <- c()
  sample.pos <- c()
  
  bal <- getRndBalances(d.counts)
  d.used <- bal$norm
  
  if (ncol(d.used) != 2) {
    for (i in 1:n.dim) {
      res.remove <- removeGroupEffect(d.used, d.groups, thresh.pc.var = 0.9)
      cell.loadings <- cbind(cell.loadings, bal$psi %*% res.remove$rotation)
      sample.pos <- cbind(sample.pos, res.remove$scores)
      d.used <- d.used - res.remove$used.part
    }
  } else {
    cell.loadings <- bal$psi
    sample.pos <- d.used
  }
  
  cn <- paste('S', 1:n.dim, sep = '')
  df.cda <- as.data.frame(sample.pos) %>% set_colnames(cn)
  df.loadings <- as.data.frame(cell.loadings * 8) %>% set_colnames(cn)
  return(list(red=df.cda, loadings=df.loadings))
}

#' Internal plotting function to plot CoDA space
#' @keywords internal
plotCodaSpaceInner <- function(df.space, df.loadings, d.groups, ref.level, target.level, font.size=3, palette=NULL) {
  group.names <- c(ref.level, target.level)
  rda.plot <- ggplot(df.space, aes(x=S1, y=S2)) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(group.names[d.groups + 1] ))) +
    labs(colour="Condition")
  
  if (!is.null(palette)) rda.plot <- rda.plot + scale_color_manual(values=palette)
  
  rda.biplot <- rda.plot +
    geom_segment(data=df.loadings, aes(x=0, xend=S1, y=0, yend=S2),
                 color="grey", arrow=arrow(length=unit(0.01,"npc")))  +
    geom_text(data=df.loadings,
              aes(x=S1, y=S2, label=rownames(df.loadings)),
              color="black", size=font.size)
  
  dx <- max(diff(range(df.space$S1)), diff(range(df.loadings$S1)))
  dy <- max(diff(range(df.space$S2)), diff(range(df.loadings$S2)))
  
  rda.biplot <- rda.biplot + coord_fixed(ratio=dy/dx)
  
  return(rda.biplot)
}

#' helper function for creating dendograms
#' @keywords internal
createDendrogram <- function(dend.data, angle=90, plot.theme=theme_get(), font.size=3, hjust=1) {
  ggplot() +
    geom_segment(data = dend.data$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + plot.theme +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"),
          axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5, margin = margin()),
          panel.grid = element_blank(), panel.border=element_blank(),
          axis.line=element_blank()) +
    scale_x_continuous(breaks=dend.data$labels$x, labels=dend.data$labels$label) +
    scale_y_continuous(expand=c(0, 0, 0.1, 0))
}

#' helper function for creating dendograms
#' @keywords internal
distTreeOrder <- function(t, tree.order){
  t.labels <- t$tip.label[t$edge[t$edge[,2] < min(t$edge[,1]),2]]
  t.labels <- t.labels[t.labels %in% tree.order]
  s <- 0
  for(i in 1:length(t.labels)){
    j <- which(tree.order == t.labels[i])
    s <- s + abs(i - j) ^ 2
  }
  return(s)
}

#' Plot the contrast tree of a compositional result
#'
#' The tree of [contrastTree()] as a dendrogram: the cell types in the order of their loadings at the leaves (reference cell
#' types marked "(ref)"), the height of a node the number of cell types below it. Each node is labeled with the log2 ratio of
#' the geometric-mean proportions of its two sides between the contrast's endpoints (the right-hand side is the one that is
#' higher at the target), the significance code of its Benjamini-Hochberg adjusted permutation p-value and the share of
#' the compositional shift it explains; the point of a significant node is filled. Below the significant nodes the
#' samples' own log2 ratios of the two sides are shown (adjusted for the model's covariates, [codaSampleLogFreq()]), one
#' row per side of the contrast with the median and the quartiles, so that the separation behind a node is visible.
#' Tiles under the leaves show the reference-centered loadings (log2) of the cell types.
#'
#' @param res result of [lmCoda()] (as stored by `cao$estimateCellLoadings()`)
#' @param tree result of [contrastTree()] (default: computed from `res`, or `res$tree` when present)
#' @param alpha significance level on the adjusted node p-values (default 0.05)
#' @param show.samples show the samples' balances at the significant nodes (default TRUE)
#' @param show.text label the nodes (default TRUE)
#' @param leaf.values tiles under the leaves: `"reference-centered"` (default), `"clr"` or `"none"`
#' @param palette colors of the contrast's sides (named by the endpoint labels, by `"baseline"` / `"target"`, or two colors)
#' @param plot.theme ggplot2 theme
#' @param font.size size of the node labels and leaf labels
#' @param label.angle angle of the leaf labels
#' @param legend.position legend position (default "bottom")
#' @param min.share nodes explaining less than this share of the shift are not labeled unless significant (default 0.02)
#' @return ggplot2 object
#' @export
plotContrastTree <- function(res, tree = NULL, alpha = 0.05, show.samples = TRUE, show.text = TRUE, leaf.values = c("reference-centered", "clr", "none"),
                             palette = NULL, plot.theme = ggplot2::theme_bw(), font.size = 3, label.angle = 90, legend.position = "bottom", min.share = 0.02) {
  leaf.values <- match.arg(leaf.values)
  if (is.null(res$contrast)) stop("not a contrast result of lmCoda()")
  if (is.null(tree)) tree <- if (!is.null(res$tree) && identical(res$tree$order, names(res$contrast$loadings$obs)[order(res$contrast$loadings$obs)])) res$tree else contrastTree(res, alpha = alpha)
  nodes <- tree$nodes; nodes$significant <- is.finite(nodes$padj) & nodes$padj <= alpha
  cts <- tree$order; D <- length(cts); x.leaf <- stats::setNames(seq_len(D), cts); H <- D - 1
  members <- lapply(tree$sides, function(s) c(s$plus, s$minus))
  nodes$x <- vapply(members, function(m) mean(x.leaf[m]), numeric(1)); nodes$y <- vapply(members, length, integer(1)) - 1; nodes$range <- vapply(members, function(m) diff(range(x.leaf[m])), numeric(1))
  childOf <- function(side) { j <- which(vapply(members, function(m) setequal(m, side), logical(1))); if (length(j)) c(nodes$x[j[1]], nodes$y[j[1]]) else c(x.leaf[[side]], 0) }
  seg <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(k) { a <- childOf(tree$sides[[k]]$minus); b <- childOf(tree$sides[[k]]$plus)
    data.frame(x = c(a[1], b[1], a[1]), y = c(a[2], b[2], nodes$y[k]), xend = c(a[1], b[1], b[1]), yend = c(nodes$y[k], nodes$y[k], nodes$y[k])) }))
  ref <- res$reference$celltypes; lab <- ifelse(cts %in% ref, paste0(cts, " (ref)"), cts)
  gg <- createDendrogram(list(segments = seg, labels = data.frame(x = unname(x.leaf), label = lab)), plot.theme = plot.theme, font.size = font.size, angle = label.angle)
  # node labels: log2 ratio, significance code, share
  nodes$code <- pvalueToCode(nodes$padj, ns.symbol = ""); nodes$label <- sprintf("%+.2f%s  %.0f%%", nodes$log2.ratio, ifelse(nodes$code == "", "", paste0(" ", nodes$code)), 100 * nodes$share)
  gg <- gg + ggplot2::geom_point(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, shape = .data$significant), size = 2, colour = "black", show.legend = FALSE) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), drop = FALSE)
  nodes$labeled <- nodes$significant | (is.finite(nodes$share) & nodes$share >= min.share)
  if (show.text) gg <- gg + ggplot2::geom_text(data = nodes[nodes$labeled, , drop = FALSE], ggplot2::aes(x = .data$x, y = .data$y + 0.02 * H, label = .data$label), vjust = 0, hjust = 0.5, size = font.size, colour = "grey20")
  # the samples' balances at the significant nodes
  sig <- which(nodes$significant)
  if (show.samples && length(sig)) {
    sl <- codaSampleLogFreq(res); logf <- sl$logf; grp <- sl$group
    pts <- do.call(rbind, lapply(sig, function(k) { s <- tree$sides[[k]]; lr <- (rowMeans(logf[, s$plus, drop = FALSE]) - rowMeans(logf[, s$minus, drop = FALSE])) / log(2)
      z <- lr - mean(lr); w <- max(abs(z)); if (!is.finite(w) || w == 0) w <- 1; data.frame(node = k, sample = rownames(logf), x = nodes$x[k] + z / w * max(nodes$range[k] / 2, 0.6) * 0.9, lr = lr,
                 group = if (is.null(grp)) "sample" else if (sl$numeric) grp else as.character(grp), stringsAsFactors = FALSE) }))
    if (sl$numeric || is.null(grp)) {
      pts$y <- nodes$y[pts$node] - 0.045 * H
      gg <- if (sl$numeric) gg + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group), alpha = 0.6, size = 1.2) + ggplot2::labs(colour = res$model$contrast_spec$term %||% "sample")
            else gg + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$x, y = .data$y), colour = "grey40", alpha = 0.5, size = 1.2)
    } else {
      pts <- pts[pts$group != "other", , drop = FALSE]; pts$group <- factor(pts$group, levels = c("baseline", "target"))
      pts$y <- nodes$y[pts$node] - 0.03 * H - 0.035 * H * (as.integer(pts$group) - 1)
      summ <- do.call(rbind, lapply(split(pts, list(pts$node, pts$group), drop = TRUE), function(d) data.frame(node = d$node[1], group = d$group[1], y = d$y[1], med = stats::median(d$x), lo = stats::quantile(d$x, 0.25), hi = stats::quantile(d$x, 0.75))))
      labs <- c(baseline = sl$labels$baseline, target = sl$labels$target)
      pal <- if (is.null(palette)) c(baseline = "#5b7fa6", target = "#c0504d") else if (!is.null(names(palette)) && all(c("baseline", "target") %in% names(palette))) palette[c("baseline", "target")]
             else if (!is.null(names(palette)) && all(unlist(res$model$contrast_spec[c("den", "num")]) %in% names(palette))) stats::setNames(palette[unlist(res$model$contrast_spec[c("den", "num")])], c("baseline", "target"))
             else stats::setNames(rep_len(unname(palette), 2), c("baseline", "target"))
      gg <- gg + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$group), alpha = 0.35, size = 1) +
        ggplot2::geom_segment(data = summ, ggplot2::aes(x = .data$lo, xend = .data$hi, y = .data$y, yend = .data$y, colour = .data$group), linewidth = 0.7) +
        ggplot2::geom_point(data = summ, ggplot2::aes(x = .data$med, y = .data$y, colour = .data$group), size = 2.2, shape = 18) +
        ggplot2::scale_colour_manual(values = pal, labels = labs, name = NULL, drop = FALSE)
    }
  }
  if (leaf.values != "none") {
    v <- if (leaf.values == "clr") res$contrast$loadings$obs[cts] else res$contrast$loadings$ref_centered[cts]
    tiles <- data.frame(x = unname(x.leaf), y = -0.045 * H, value = unname(v) / log(2))
    gg <- gg + ggplot2::geom_tile(data = tiles, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value), width = 0.9, height = 0.05 * H) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B", midpoint = 0, name = if (leaf.values == "clr") "loading (log2)" else "reference-centered loading (log2)")
  }
  gg + ggplot2::coord_cartesian(ylim = c(-0.1 * H, H * 1.12), clip = "off") + ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(legend.position = legend.position, axis.text.x = ggplot2::element_text(angle = label.angle, hjust = 1, vjust = 0.5))
}

#' Internal function to plot global statistics from CoDA analysis (NEW)
#' @keywords internal
plotCodaGlobal <- function(x, style = c("bar", "violin"), fdr.threshold = 0.05, drop.intercept = TRUE,
                           use.contrast.title = FALSE) {
  style <- match.arg(style)
  
  ## -----------------------------
  ## 1. Collect global stats
  ## -----------------------------
  # Contrast
  contrast.label <- if (!is.null(x$contrast$label)) {
    x$contrast$label
  } else {
    "contrast"
  }
  
  terms   <- "contrast"
  labels  <- contrast.label
  stats   <- x$contrast$global$stat
  pvals   <- x$contrast$global$p
  
  # Coefficients
  coef.names <- names(x$coefficients)
  if (drop.intercept && "(Intercept)" %in% coef.names) {
    coef.names <- setdiff(coef.names, "(Intercept)")
  }
  
  if (length(coef.names) > 0) {
    terms  <- c(terms, coef.names)
    labels <- c(labels, coef.names)
    stats  <- c(stats,vapply(coef.names, function(nm) x$coefficients[[nm]]$global$stat, numeric(1)))
    pvals  <- c(pvals,vapply(coef.names, function(nm) x$coefficients[[nm]]$global$p, numeric(1)))
  }
  
  df.stats <- data.frame(term.id = terms, term.label = labels, stat = stats, pval = pvals, stringsAsFactors = FALSE)
  df.stats$padj   <- p.adjust(df.stats$pval, method = "fdr")
  df.stats$signif <- df.stats$padj < fdr.threshold
  
  ## 1–3 star scheme based on adjusted p-value
  df.stats$stars <- ifelse(
    df.stats$padj < 0.001, "***",
    ifelse(
      df.stats$padj < 0.01, "**",
      ifelse(df.stats$padj < 0.05, "*", "")
    )
  )
  
  # Order terms by decreasing variance explained
  df.stats$term.label <- factor(df.stats$term.label, levels = df.stats$term.label[order(-df.stats$stat)])
  
  # Optional title from contrast label
  plot.title <- NULL
  if (use.contrast.title && !is.null(x$contrast$label)) {
    plot.title <- x$contrast$label
  }
  
  ## ========================================================
  ## 2a. Bar-style global plot
  ## ========================================================
  if (style == "bar") {
    p <- ggplot(df.stats,
                aes(x = term.label,
                    y = stat,
                    fill = signif)) +
      geom_col() +
      geom_text(data = subset(df.stats, stars != ""), aes(label = stars), vjust = -0.0, size = 6, color = "black") +
      scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), name = paste0("FDR < ", fdr.threshold),
                        labels = c(`TRUE` = "yes", `FALSE` = "no")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = "Term (contrast / coefficient)", y = x$contrast$global$statistic %||% "Variance explained in ILR space", title = plot.title) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  }
  
  ## ========================================================
  ## 2b. Violin-style view of null distribution
  ## ========================================================
  # Build long data frame of permuted stats for each term
  
  perm.list <- list()
  
  # Contrast
  if (!is.null(x$contrast$global$stat_perm)) {
    perm.list[[length(perm.list) + 1]] <- data.frame(term.id = "contrast", term.label = contrast.label,
                                                     stat.perm  = as.numeric(x$contrast$global$stat_perm),
                                                     stringsAsFactors = FALSE)
  }
  
  # Coefficients
  if (length(coef.names) > 0) {
    for (nm in coef.names) {
      sp <- x$coefficients[[nm]]$global$stat_perm
      if (is.null(sp)) next  # if missing, skip that term
      perm.list[[length(perm.list) + 1]] <- data.frame(term.id = nm, term.label = nm, stat.perm  = as.numeric(sp),
                                                       stringsAsFactors = FALSE)
    }
  }
  
  if (length(perm.list) == 0) {
    stop("No permutation distributions found in x; ",
         "cannot produce violin-style global plot.")
  }
  
  df.perm <- do.call(rbind, perm.list)
  
  # Match factor levels to df_stats ordering
  df.perm$term.label <- factor(df.perm$term.label,
                               levels = levels(df.stats$term.label))
  
  p <- ggplot(df.perm, aes(x = term.label, y = stat.perm)) +
    geom_violin(fill = "grey85", color = NA, trim = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.01, size = 0.6) +  # lower alpha as you set
    # Observed stat overlaid
    geom_point(data = df.stats, aes(x = term.label, y = stat, color = signif), size = 2) +
    geom_text(data = subset(df.stats, stars != ""), aes(x = term.label, y = stat, label = stars),
              vjust = -0.0, size = 6, color = "black") +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black"), name = paste0("FDR < ", fdr.threshold),
                       labels = c(`TRUE` = "yes", `FALSE` = "no")) +
    labs(x = "Term (contrast or coefficient)", y = x$contrast$global$statistic %||% "Variance explained in ILR space", title = plot.title) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Internal function to plot heatmap of CoDA coefficients and contrast (NEW)
#' @keywords internal
plotCodaHeatmap <- function(x, coef.names = NULL, fdr.threshold = 0.05) {
  # -----------------------------
  # 1. Decide which coefficients
  # -----------------------------
  all.coef.names <- names(x$coefficients)
  
  if (is.null(coef.names)) {
    coef.names <- all.coef.names
  } else {
    missing <- setdiff(coef.names, all.coef.names)
    if (length(missing) > 0)
      stop("Some coef.names not found in x$coefficients: ",
           paste(missing, collapse = ", "))
  }
  
  # -----------------------------
  # 2. Build DF for coefficients
  # -----------------------------
  df.list <- lapply(coef.names, function(nm) {
    cf <- x$coefficients[[nm]]
    data.frame(celltype = names(cf$loadings), coef = nm, effect = as.numeric(cf$loadings),
               padj = as.numeric(cf$padj_cell), stringsAsFactors = FALSE)
  })
  df.coef <- do.call(rbind, df.list)
  
  # -----------------------------
  # 3. Build DF for contrast
  # -----------------------------
  if (is.null(x$contrast)) {
    stop("x$contrast not found in lmCoda result; cannot add contrast column.")
  }
  
  # Always use contrast label if available
  contrast.col.label <- if (!is.null(x$contrast$label)) {
    x$contrast$label
  } else {
    "contrast"
  }
  
  df.contrast <- data.frame(celltype = names(x$contrast$loadings$obs), coef = contrast.col.label,
                            effect = as.numeric(x$contrast$loadings$obs), padj = as.numeric(x$contrast$per_cell$padj),
                            stringsAsFactors = FALSE)
  
  # -----------------------------
  # 4. Combine contrast + coefs
  # -----------------------------
  df <- rbind(df.contrast, df.coef)
  
  # 1–3 star scheme based on adjusted p-value
  df$stars <- ifelse(
    df$padj < 0.001, "***",
    ifelse(
      df$padj < 0.01,  "**",
      ifelse(df$padj < 0.05, "*", "")
    )
  )
  df$signif <- df$padj < fdr.threshold
  
  # Order coef axis: contrast first, then chosen coef.names
  df$coef <- factor(df$coef, levels = c(contrast.col.label, coef.names))
  
  # Order cell types (alphabetical by default; change here if you want a custom order)
  df$celltype <- factor(df$celltype, levels = sort(unique(df$celltype)))
  
  # -----------------------------
  # 5. Plot
  # -----------------------------
  ggplot(df, aes(x = coef, y = celltype, fill = effect)) +
    geom_tile(color = "white") +
    # Stars for significance (larger text symbols)
    geom_text(data = subset(df, stars != ""), aes(label = stars), size = 4, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Term (contrast / coefficient)", y = "Cell type", fill = "Effect\n(log-ratio)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Internal function to plot cell-type effects from CoDA analysis (NEW)
#' @param x lmCoda result object
#' @param coef Coefficient to plot (default "contrast" for the main contrast)
#' @param show.ci `value = "predicted"`: error bars with the intervals of the fitted proportions (default TRUE; `lmCoda(ci.level)`,
#'   from the ILR fit: the endpoint's fitted ILR coordinates with their covariance, multivariate t draws mapped to proportions)
#' @param value Type of effect to plot: "predicted" (predicted baseline vs target compositions; contrast only),
#'   "reference-centered" (effects relative to reference clusters), or "volcano" (log2 effect vs -log10 p-value)    
#' @param style Style of plot for reference-centered effects only: "bar" or "violin"
#' @param fdr.threshold FDR threshold for significance
#' @param sortCellTypes Whether to sort cell types by adjusted p-value
#' @param labelCellTypes For volcano plots, which cell types to label: "significant", "all", or "none"
#' @param celltype.palette Optional named vector of colors for cell types
#' @param use.contrast.title Whether to use the contrast label as the plot title (if available)
#' @param significance p-values behind the stars: `"clr"` (default: the primary test, each cell type relative to the geometric mean of
#'   all cell types, which does not depend on the reference), or `"reference"` (relative to the reference cell types; the
#'   reference is chosen from the data, so this test is conditional on it)
#' @return ggplot2 object
#' @keywords internal
plotCodaCellTypes <- function(x, coef = c("contrast"), value = c("predicted", "reference-centered", "volcano"),
                              style = c("bar", "violin"), fdr.threshold = 0.05, sortCellTypes = FALSE,
                              labelCellTypes = c("significant", "all", "none"), celltype.palette = NULL,
                              use.contrast.title = FALSE, significance = c("clr", "reference"), show.ci = TRUE) {
  value <- match.arg(value); significance <- match.arg(significance); use.ref <- FALSE
  coef  <- match.arg(coef, choices = c("contrast", names(x$coefficients)))
  style <- match.arg(style)
  labelCellTypes <- match.arg(labelCellTypes)
  
  ## -----------------------------
  ## 1. Choose source of effects
  ## -----------------------------
  if (identical(coef, "contrast")) {
    beta.ilr  <- x$contrast$beta_ilr
    pval.cell <- x$contrast$per_cell$pval
    padj.cell <- x$contrast$per_cell$padj
    loadings  <- x$contrast$loadings$obs         # natural-log scale
    ref.loads <- x$contrast$loadings$ref_centered
    # the stars: the primary test (relative to the geometric mean of all cell types, reference-free) by default; the test relative to the
    # reference cell types only on request (the reference is chosen from the data and can change with the relabelings drawn)
    use.ref <- significance == "reference"
    if (use.ref) { if (is.null(x$contrast$per_cell$padj.ref)) stop("no reference-relative p-values in the result"); pval.cell <- x$contrast$per_cell$pval.ref; padj.cell <- x$contrast$per_cell$padj.ref }
  } else {
    if (!coef %in% names(x$coefficients)) {
      stop("Coefficient '", coef, "' not found in x$coefficients.")
    }
    beta.ilr  <- x$coefficients[[coef]]$beta_ilr
    pval.cell <- x$coefficients[[coef]]$pval_cell
    padj.cell <- x$coefficients[[coef]]$padj_cell
    loadings  <- x$coefficients[[coef]]$loadings  # natural-log scale
    
    # For coefficients, reference-centering uses the same reference cluster as contrast
    ref.idx   <- x$reference$idx
    ref.mean  <- mean(loadings[ref.idx])
    ref.loads <- loadings - ref.mean
  }
  
  ## -----------------------------
  ## 2. Common per-cell stats
  ## -----------------------------
  df.common <- data.frame(celltype = names(loadings), effect = as.numeric(loadings),    # natural-log scale
                          effect_ref = as.numeric(ref.loads),   # natural-log scale
    pval = as.numeric(pval.cell), padj = as.numeric(padj.cell), stringsAsFactors = FALSE)
  df.common$signif <- df.common$padj < fdr.threshold
  df.common$is_ref <- df.common$celltype %in% x$reference$celltypes
  
  # Convert to log2 scale (linear rescaling of natural-log effects)
  log2_factor <- 1 / log(2)
  df.common$effect_log2     <- df.common$effect     * log2_factor
  df.common$effect_ref_log2 <- df.common$effect_ref * log2_factor
  
  # 1–3 star scheme based on adjusted p-value
  df.common$stars <- ifelse(
    df.common$padj < 0.001, "***",
    ifelse(
      df.common$padj < 0.01,  "**",
      ifelse(df.common$padj < 0.05, "*", "")
    )
  )
  
  star.caption <- paste0("Stars: FDR relative to the ", if (use.ref) "reference cell types" else "geometric mean of all cell types", " (* < 0.05, ** < 0.01, *** < 0.001)")
  
  ## -----------------------------
  ## 3. Cell type ordering
  ## -----------------------------
  if (sortCellTypes) {
    ord <- order(df.common$padj, df.common$celltype, na.last = TRUE)
    df.common <- df.common[ord, , drop = FALSE]
    celltypes <- df.common$celltype
  } else {
    celltypes <- names(loadings)
  }
  celltypes <- unique(celltypes)
  
  ## Stable cell-type palette
  if (is.null(celltype.palette)) {
    cols <- scales::hue_pal()(length(celltypes))
    names(cols) <- celltypes
  } else {
    cols <- celltype.palette
    missing.cols <- setdiff(celltypes, names(cols))
    if (length(missing.cols) > 0) {
      stop("celltype.palette is missing colors for: ",
           paste(missing.cols, collapse = ", "))
    }
    cols <- cols[celltypes]
  }
  
  ## Optional plot title from contrast label
  plot.title <- NULL
  if (use.contrast.title && !is.null(x$contrast$label)) {
    plot.title <- x$contrast$label
  }
  
  ## ========================================================
  ## Mode 1: predicted baseline vs contrast (contrast only)
  ## ========================================================
  if (value == "predicted") {
    if (!identical(coef, "contrast")) {
      stop("value = 'predicted' is currently only implemented for coef = 'contrast'.")
    }
    if (is.null(x$contrast$predicted)) {
      stop("No predicted compositions found. Ensure model$contrast_endpoints_* is set and lmCoda was run with that model.")
    }
    
    baseline <- x$contrast$predicted$baseline
    target   <- x$contrast$predicted$target
    pr <- x$contrast$predicted
    show.ci <- show.ci && all(c("baseline.lo", "baseline.hi", "target.lo", "target.hi") %in% names(pr))
    if (show.ci) star.caption <- paste0(star.caption, sprintf("\nError bars: %.0f%% intervals of the fitted proportions (ILR fit)", 100 * (pr$ci.level %||% 0.95)))
    
    # human-readable labels for endpoints (from lmCoda/model)
    label.baseline <- "baseline"
    label.target   <- "target"
    if (!is.null(x$contrast$predicted$labels)) {
      labs <- x$contrast$predicted$labels
      if (!is.null(labs$baseline)) label.baseline <- labs$baseline
      if (!is.null(labs$target))   label.target   <- labs$target
    }
    
    # enforce consistent cell-type ordering
    celltypes.order <- celltypes
    df <- data.frame(celltype = celltypes.order, baseline = as.numeric(baseline[celltypes.order]),
                     target = as.numeric(target[celltypes.order]), stringsAsFactors = FALSE)
    
    df.long <- reshape2::melt(df, id.vars = "celltype", variable.name = "condition", value.name = "prop")
    df.long$lo <- df.long$hi <- NA_real_
    if (show.ci) {   # the intervals of the fitted proportions
      b <- df.long$condition == "baseline"
      df.long$lo[b] <- as.numeric(pr$baseline.lo[df.long$celltype[b]]); df.long$hi[b] <- as.numeric(pr$baseline.hi[df.long$celltype[b]])
      df.long$lo[!b] <- as.numeric(pr$target.lo[df.long$celltype[!b]]); df.long$hi[!b] <- as.numeric(pr$target.hi[df.long$celltype[!b]])
    }
    df.long$top <- ifelse(is.finite(df.long$hi), pmax(df.long$hi, df.long$prop), df.long$prop)   # where the stars go
    
    # merge padj/stars/signif from df_common
    df.long <- merge(df.long, df.common[, c("celltype", "padj", "signif", "stars")], by = "celltype", all.x = TRUE)
    
    df.long$condition <- factor(df.long$condition,
                                levels = c("baseline", "target"))
    
    p <- ggplot(df.long,aes(x = celltype, y = prop, fill = celltype, alpha = condition)) +
      geom_col(position = position_dodge(width = 0.8), color = NA) +
      { if (show.ci) geom_errorbar(aes(x = celltype, ymin = lo, ymax = hi, group = condition), position = position_dodge(width = 0.8), width = 0.3,
                                   colour = "grey30", linewidth = 0.4, inherit.aes = FALSE, na.rm = TRUE) } +
      # 1–3 stars only on target bars (explicit aes + no legend), above the bar or its interval
      geom_text(data = subset(df.long, condition == "target" & stars != ""),
        aes(x = celltype, y = top, label = stars), position = position_nudge(x = 0.2), vjust = 0, size = 6,
            color = "black", inherit.aes = FALSE, show.legend = FALSE) +
      scale_fill_manual(values = cols, guide = "none") +  # no cell-type legend
      scale_alpha_manual(values = c(baseline = 0.2, target = 1.0), name = "Condition",
                         labels = c(baseline = label.baseline, target = label.target)) +
      # Add headroom so stars are not clipped
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(x = "Cell type", y = "Predicted proportion", title = plot.title, caption = star.caption) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  }
  
  ## ========================================================
  ## Mode 2: reference-centered cell-type effect
  ## ========================================================
  if (value == "reference-centered") {
    # Enforce cell-type order
    df.common$celltype <- factor(df.common$celltype, levels = celltypes)
    
    if (style == "bar") {
      df <- df.common
      df$signif_flag <- ifelse(df$signif, "signif", "ns")
      df$is_ref_flag <- ifelse(df$is_ref, "ref", "nonref")
      
      p <- ggplot(df, aes(x = celltype, y = effect_ref_log2, fill = celltype, alpha = signif_flag,
                          color = is_ref_flag)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_col() +
        # 1–3 stars on significant bars (no legend participation)
        geom_text(data = subset(df, stars != ""), aes(x = celltype, y = effect_ref_log2, label = stars),
          vjust = 0, size = 6,color = "black", inherit.aes = FALSE, show.legend = FALSE) +
        # cell-type colors but no legend
        scale_fill_manual(values = cols, guide = "none") +
        # significance legend via alpha
        scale_alpha_manual(values = c(signif = 1.0, ns = 0.25),name = paste0("FDR < ", fdr.threshold)) +
        # Reference cell type legend: darker outline, very light fill
        scale_color_manual(values = c(ref = "black", nonref = "grey70"), name = "Reference cell types",
                           labels = c(ref = "Reference", nonref = "Other"),
                           guide  = guide_legend(override.aes = list(fill = c("grey95", "grey95")))) +
        # symmetric padding so both top and bottom have headroom
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        labs(x = "Cell type",y = "Reference-centered effect (log2-ratio)",title = plot.title, caption = star.caption) +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      return(p)
    }
    
    ## --------- style == "violin": null distribution per cell type ----------
    psi <- x$psi
    
    if (identical(coef, "contrast")) {
      if (is.null(x$fit$stats.perm)) {
        stop("No stats.perm found in x$fit; cannot build null for contrast.")
      }
      beta.perm <- as.matrix(x$fit$stats.perm)  # n_perm × K
      ell.perm  <- psi %*% t(beta.perm)         # D × n_perm
      ref.idx  <- x$reference$idx
      ref.mean <- mean(loadings[ref.idx])
      ell.perm.ref <- ell.perm - ref.mean       # natural-log scale
    } else {
      if (is.null(x$fit$sampled.fits)) {
        stop("No sampled.fits in x$fit; cannot build null for coefficients.")
      }
      coef.mat <- x$fit$coef
      j <- match(coef, rownames(coef.mat))
      if (is.na(j)) stop("Could not find coefficient index for ", coef)
      
      beta.j.perm <- sapply(x$fit$sampled.fits, function(mat) mat[, j])
      ell.perm <- psi %*% t(beta.j.perm)        # D × n_perm
      ref.idx  <- x$reference$idx
      ref.mean <- mean(loadings[ref.idx])
      ell.perm.ref <- ell.perm - ref.mean       # natural-log scale
    }
    
    # Convert permuted effects to log2 scale
    ell.perm.ref.log2 <- ell.perm.ref * log2_factor
    
    # Build long DF for violins
    celltypes_order <- rownames(ell.perm.ref.log2)
    df.perm <- reshape2::melt(
      ell.perm.ref.log2,
      varnames = c("celltype", "perm_id"),
      value.name = "effect_ref_perm_log2"
    )
    df.perm$celltype <- factor(df.perm$celltype, levels = celltypes)  # sorted / unsorted order
    
    # Observed summary DF (from df.common, log2 scale)
    df.obs <- df.common
    df.obs$celltype <- factor(df.obs$celltype, levels = celltypes)
    
    p <- ggplot(df.perm, aes(x = celltype, y = effect_ref_perm_log2)) +
      geom_violin(fill = "grey85", color = NA, trim = FALSE) +
      geom_jitter(width = 0.15, alpha = 0.01, size = 0.6) +
      geom_point(data = df.obs, aes(x = celltype, y = effect_ref_log2, color = signif), size = 2) +
      geom_text(data = subset(df.obs, stars != ""), aes(x = celltype, y = effect_ref_log2, label = stars),
        vjust = -0.4, size = 6, color = "black", inherit.aes = FALSE, show.legend = FALSE) +
      scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black"), name = paste0("FDR < ", fdr.threshold),
                         labels = c(`TRUE` = "yes", `FALSE` = "no")) +
      labs(x = "Cell type", y = "Reference-centered effect (log2-ratio)", title = plot.title, caption = star.caption) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
  }
  
  ## ========================================================
  ## Mode 3: volcano-style view
  ## ========================================================
  if (value == "volcano") {
    df <- df.common
    df$mlog10p     <- -log10(df$pval)
    df$signif_flag <- ifelse(df$signif, "signif", "ns")
    df$celltype    <- factor(df$celltype, levels = celltypes)
    
    # Decide which cell types to label
    if (labelCellTypes == "all") {
      df.label <- df
    } else if (labelCellTypes == "significant") {
      df.label <- subset(df, signif)
    } else { # "none"
      df.label <- df[0, , drop = FALSE]  # empty
    }
    
    p <- ggplot(df, aes(x = effect, y = mlog10p)) +
      geom_hline(yintercept = -log10(fdr.threshold), linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_point(aes(color = celltype, alpha = signif_flag), size = 2) +
      geom_text(data = df.label, aes(x = effect, y = mlog10p, label = celltype),
        vjust = -0.4, size = 3, inherit.aes = FALSE, show.legend = FALSE) +
      scale_color_manual(values = cols, name = "Cell type") +
      scale_alpha_manual(values = c(signif = 1.0, ns = 0.2), name = paste0("FDR < ", fdr.threshold)) +
      labs(x = "Effect (log-ratio / loading)", y = "-log10(p-value)", title = plot.title, caption = star.caption) +
      theme_bw()
    
    return(p)
  }
}

#' Internal plotting function for cell loadings (old two group version)
#' @keywords internal
plotCellLoadings <- function(loadings, pval, ref.level, target.level, signif.threshold=0.05, jitter.alpha=0.1,
                             palette=NULL, show.pvals=FALSE, plot.theme=theme_get(), jitter.size=1,
                             ordering=c("pvalue", "loading"), ref.load.level=0, annotation.x=NULL, annotation.y=1) {
  ordering <- match.arg(ordering)
  xintercept <- ref.load.level
  
  loading.order <- order(abs(rowMeans(loadings)))
  loadings <- loadings[loading.order, ]
  if (!is.null(pval)) {
    # if some p-values are the same - then order by mean, therefore a prior sorting is required
    pval <- pval[loading.order]
    
    if (ordering == "pvalue") {
      # additional ordering by p-value
      loadings <- loadings[order(-pval), ]
      pval <- pval[order(-pval)]
    }
    
    # Get significant cells
    n.significant.cells <- sum(pval <= signif.threshold)
  } else {
    n.significant.cells <- 0
  }
  
  # Normalization of loadings
  loadings <- loadings - xintercept
  
  res.ordered <- t(loadings) %>% as.data.frame()
  if (is.null(annotation.x)) {
    annotation.x <- max(loadings)
  }
  
  p <- ggplot(stack(res.ordered), aes(y=ind, x=values, fill=factor(ind))) +
    geom_boxplot(notch=TRUE, outlier.shape = NA)
  
  if ((jitter.alpha > 1e-5) && (jitter.size > 1e-5)) {
    p <- p + geom_jitter(alpha=jitter.alpha, size=jitter.size)
  }
  
  p <- p + geom_vline(xintercept=0, color = "grey37") +
    labs(y='', x='loadings') + plot.theme + theme(legend.position = "none") +
    scale_y_discrete(position = "right") + xlim(-1, 1)
  
  # Add text
  p <- p +
    annotate('text', x=-annotation.x, y=annotation.y, label=paste('\u2190', ref.level), hjust='left') +
    annotate('text', x=annotation.x, y=annotation.y, label=paste(target.level, '\u2192'), hjust='right')
  
  if (!is.null(palette)) p <- p + scale_fill_manual(values=palette)
  if ((n.significant.cells > 0) && (ordering == "pvalue")) {
    p <- p + geom_hline(yintercept=nrow(loadings) - n.significant.cells + 0.5, color='red')
  }
  
  if (show.pvals) {
    d <- data.frame(y=-log(pval, base=10), x=names(pval), row=1:length(pval))
    p.pval <- ggplot(d, aes(y=reorder(x, row), x=y, fill=factor(x))) +
      geom_bar(stat="identity") +
      geom_vline(xintercept=-log(signif.threshold, base=10)) +
      scale_x_continuous(expand=c(0, 0, 0.02, 0)) +
      labs(y='', x='-log10(adj. p-value)') +
      plot.theme + theme(legend.position="none") + theme(axis.text.y=element_blank())
    
    if (!is.null(palette)) p.pval <- p.pval + scale_fill_manual(values=palette)
    
    p.combo <- cowplot::plot_grid(p, p.pval, nrow=1, rel_widths=c(2, 1))
    return(p.combo)
  }
  
  return(p)
  
}

makeNamedPalette <- function(levels) {
    vals <- makeHCLPal(length(levels))
    stats::setNames(vals, levels)
}

makeHCLPal <- function(n) {
    base <- c("#0072B2","#E69F00","#009E73","#D55E00","#CC79A7","#56B4E9","#F0E442")
    if (n <= length(base)) base[seq_len(n)] else {
        extra_n <- n - length(base)
        extra <- grDevices::hcl(h = seq(15, 375, length.out = extra_n + 1)[-1], c = 100, l = 55)
        c(base, extra)
    }
}


pickColsForKey <- function(pm, key) {
  cn <- colnames(pm)
  # escape regex metacharacters
  k  <- gsub("([.|()\\^{}+$*?\\[\\]\\\\-])", "\\\\\\1", key)

  ## 1) factor covariate columns: "<key>_pair"
  exact_pair <- paste0(key, "_pair")
  if (exact_pair %in% cn)
    return(exact_pair)

  ## 2) numeric covariates: "pair_<key>_diff" or "pair_<key>_mean"
  # exact matches
  exact_diff <- paste0("pair_", key, "_diff")
  exact_mean <- paste0("pair_", key, "_mean")

  # return *_diff if exists, else *_mean
  if (exact_diff %in% cn)
    return(exact_diff)
  if (exact_mean %in% cn)
    return(exact_mean)

  ## 3) regex fallback: "pair_<something>_<key>_(mean|diff)"
  rx <- grep(paste0("^pair_.*", k, ".*_(mean|diff)$"), cn, value = TRUE)
  if (length(rx)) {
    diffs <- grep("_diff$", rx, value = TRUE)
    if (length(diffs)) return(diffs)
    means <- grep("_mean$", rx, value = TRUE)
    if (length(means)) return(means)
    return(rx)
  }

  ## none found
  character(0)
}

# ggplot2 volcano used when EnhancedVolcano is not installed
plotVolcanoSimple <- function(de.df, p.name = "padj", color.var = "CellFrac", palette = brewerPalette("RdYlBu"), lf.cutoff = 1.5, p.cutoff = 0.05,
                              cell.frac.cutoff = 0.2, size = c(0.1, 1.0), lab.size = 2, sel.labels = NULL, plot.theme = theme_get(), legend.pos = "none",
                              max.labels = 30) {
  df <- as.data.frame(de.df)
  if (!"Gene" %in% names(df)) df$Gene <- rownames(df)
  df <- df[is.finite(df$log2FoldChange) & is.finite(df[[p.name]]), ]
  df$nlp <- -log10(pmax(df[[p.name]], 1e-300))
  has.col <- color.var %in% names(df) && is.numeric(df[[color.var]])
  if (is.null(sel.labels)) {
    sig <- df[[p.name]] <= p.cutoff & abs(df$log2FoldChange) >= lf.cutoff
    if (has.col && color.var == "CellFrac") sig <- sig & df$CellFrac >= cell.frac.cutoff
    sel.labels <- unique(c(df$Gene[sig][order(df[[p.name]][sig])], df$Gene[order(df[[p.name]])][1:min(10, nrow(df))]))
  }
  sel.labels <- utils::head(sel.labels, max.labels)
  df$label <- ifelse(df$Gene %in% sel.labels, df$Gene, "")
  gg <- ggplot(df, aes(x = .data$log2FoldChange, y = .data$nlp)) +
    geom_vline(xintercept = c(-lf.cutoff, lf.cutoff), linetype = 2, colour = "grey60") +
    geom_hline(yintercept = -log10(p.cutoff), linetype = 2, colour = "grey60")
  if (has.col) {
    gg <- gg + geom_point(aes(colour = .data[[color.var]], size = .data[[color.var]]), alpha = 0.8) +
      scale_size_continuous(range = size, name = "Expr. frac", limits = c(0, 1)) + val2ggcol(df[[color.var]], palette = palette, color.range = c(0, 1)) +
      guides(color = guide_colorbar(title = "Expr. frac"))
  } else gg <- gg + geom_point(size = mean(size), colour = "grey40", alpha = 0.8)
  gg + ggrepel::geom_text_repel(aes(label = .data$label), size = lab.size, max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey50") +
    labs(x = "log2 fold change", y = sprintf("-log10 %s", p.name)) + plot.theme + theme_legend_position(if (legend.pos == "none") "right" else legend.pos)
}

#' Outline of the significant regions of a cell map
#'
#' Where on the embedding the significant cells are: the local fraction of significant cells, estimated by a two-dimensional
#' kernel density of the significant cells over that of all cells on a grid, is contoured at `level`. Grid points with too few
#' cells around them (fewer than `min.cells` expected in a grid cell) are left out, so that no contour runs through empty space.
#' @param embedding cells x 2 matrix with cell names as row names
#' @param significant logical vector named by cell (`NA` = not tested; counts as not significant)
#' @param level fraction of significant cells enclosed by the outline (default 0.5: regions where significant cells are the majority)
#' @param bandwidth kernel bandwidth on each axis (default: 1/30 of the axis range)
#' @param n grid size per axis (default 200)
#' @param min.cells expected cells per grid cell below which the grid point is masked (default 1)
#' @param min.area closed loops enclosing less than this fraction of the embedding's bounding box are dropped (default 0.001:
#'   specks around a handful of cells)
#' @return data.frame of contour paths (`x`, `y`, `group`), empty when there is nothing to outline
#' @export
significanceContours <- function(embedding, significant, level = 0.5, bandwidth = NULL, n = 200, min.cells = 1, min.area = 0.001) {
  checkPackageInstalled("MASS", cran = TRUE)
  emb <- as.matrix(embedding[, 1:2]); cells <- rownames(emb); ok <- is.finite(emb[, 1]) & is.finite(emb[, 2])
  sig <- as.logical(significant[cells]); sig[is.na(sig)] <- FALSE
  empty <- data.frame(x = numeric(0), y = numeric(0), group = character(0), stringsAsFactors = FALSE)
  if (sum(sig & ok) < 10) return(empty)
  rx <- range(emb[ok, 1]); ry <- range(emb[ok, 2]); pad <- c(diff(rx), diff(ry)) * 0.03
  lims <- c(rx[1] - pad[1], rx[2] + pad[1], ry[1] - pad[2], ry[2] + pad[2])
  h <- bandwidth %||% (c(diff(rx), diff(ry)) / 30)
  d.all <- MASS::kde2d(emb[ok, 1], emb[ok, 2], h = h, n = n, lims = lims)
  d.sig <- MASS::kde2d(emb[ok & sig, 1], emb[ok & sig, 2], h = h, n = n, lims = lims)
  cell.area <- diff(d.all$x[1:2]) * diff(d.all$y[1:2])
  expected <- d.all$z * sum(ok) * cell.area                       # cells expected in a grid cell
  frac <- (d.sig$z * sum(ok & sig)) / pmax(d.all$z * sum(ok), 1e-300)
  frac[expected < min.cells] <- 0                                  # masked: no contour through empty space
  cl <- grDevices::contourLines(d.all$x, d.all$y, frac, levels = level)
  if (!length(cl)) return(empty)
  box <- diff(range(d.all$x)) * diff(range(d.all$y))
  area <- vapply(cl, function(l) { x <- l$x; y <- l$y; abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2 }, numeric(1))   # shoelace; open paths count as closed
  cl <- cl[area >= min.area * box]
  if (!length(cl)) return(empty)
  do.call(rbind, lapply(seq_along(cl), function(i) data.frame(x = cl[[i]]$x, y = cl[[i]]$y, group = paste0("c", i), stringsAsFactors = FALSE)))
}

# the ggplot layer of significanceContours() for an embedding plot (NULL when there is nothing to outline)
significanceContourLayer <- function(embedding, significant, level = 0.5, color = "black", linewidth = 0.25, ...) {
  d <- significanceContours(embedding, significant, level = level, ...)
  if (!nrow(d)) return(NULL)
  ggplot2::geom_path(data = d, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group), colour = color, linewidth = linewidth, inherit.aes = FALSE)
}
