#' @import tibble
#' @import cowplot
#' @import dplyr
#' @import magrittr
#' @importFrom reshape2 melt
NULL

theme_legend_position <- function(position) {
  theme(legend.position=position, legend.justification=position)
}

#' @export
brewerPalette <- function(name, n=NULL, rev=TRUE) {
  checkPackageInstalled("RColorBrewer", cran=TRUE)
  if (is.null(n)) n <- RColorBrewer::brewer.pal.info[name,]$maxcolors
  pal <- RColorBrewer::brewer.pal(n, name)
  if (rev) pal <- rev(pal)
  return(grDevices::colorRampPalette(pal))
}

dark.red.palette <- colorRampPalette(c("gray95", "red", "#5A0000"), space="Lab")

#' @title Plot Number of Cells Regression
#' @param n the regressed variable
#' @param n.total number of cells
#' @param x.lab x axis label (default: "Number of cells").
#' @param y.lab y axis label (default: "N")
#' @param legend.position legend position (default: "right")
#' @param label show cell type labels (default: TRUE)
#' @param size label text size (default: 4)
#' @param palette color palette
#' @param plot.line plot the robust regression (default: TRUE)
#' @param line.width regression line width (default: 0.5)
plotNCellRegression <- function(n, n.total, x.lab="Number of cells", y.lab="N", legend.position="right", label=TRUE, size=4,
                                palette=NULL, plot.line=TRUE, line.width=0.5, plot.theme=theme_get()) {
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

#' @title Plot Count Boxplots Per Type
#' @param count.df data.frame with columns `group`, `variable` and `value`
#' @param notch Whether to show notch in the boxplots
#' @param alpha Transparency level on the data points (default: 0.2)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param size marker size (default: 0.5)
#' @param jitter.width width of the point jitter (default: 0.15)
plotCountBoxplotsPerType <- function(count.df, y.lab="count", x.lab="", y.expand=c(0, 0), show.significance=FALSE, palette=NULL,
                                     notch=FALSE, legend.position="right", alpha=0.2, plot.theme=theme_get(), size=0.5, jitter.width=0.15,
                                     adjust.pvalues=T) {
  gg <- ggplot(count.df, aes(x=variable, y=value, by=group, fill=group)) +
    geom_boxplot(position=position_dodge(), outlier.shape = NA, notch=notch) +
    labs(x=x.lab, y=y.lab) +
    plot.theme +
    theme_legend_position(legend.position) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_blank()) +
    geom_point(position=position_jitterdodge(jitter.width=jitter.width), color="black", size=size, alpha=alpha) +
    scale_y_continuous(expand=y.expand, limits=c(0, max(count.df$value) * 1.05))


  if (show.significance) {
    if (adjust.pvalues) {
      gg <- gg + ggpubr::stat_compare_means(aes(group = group), label = "p.adj")  # willcox test + adjustment
      # TODO
      # p.adjust.method = "fdr" ?
    } else {
      gg <- gg + ggpubr::stat_compare_means(aes(group = group), label = "p.signif")  # willcox test
    }
  }

  if (!is.null(palette)) gg <- gg + scale_fill_manual(values=palette)
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
plotHeatmap <- function(df, color.per.group=NULL, row.order=TRUE, col.order=TRUE, legend.position="right",
                        legend.key.width=unit(8, "pt"), legend.title="-log10(p-value)", x.axis.position="top",
                        color.range=NULL, plot.theme=theme_get(), symmetric=FALSE, palette=NULL, font.size=8,
                        distance="manhattan", clust.method="complete") {
  if (is.null(color.range)) {
    color.range <- c(min(0, min(df)), max(df))
  }

  if (is.logical(row.order) && row.order) {
    row.order <- rownames(df)[dist(df, method=distance) %>% hclust(method=clust.method) %>% .$order] %>% rev()
  } else if (is.logical(row.order) && !row.order) {
    row.order <- rownames(df)
  }

  if (is.logical(col.order) && col.order) {
    col.order <- colnames(df)[dist(t(df), method=distance) %>% hclust(method=clust.method) %>% .$order] %>% rev()
  } else if (is.logical(col.order) && !col.order) {
    col.order <- colnames(df)
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
  gg <- ggplot(df) + geom_tile(aes(x=G2, y=G1, fill=value), color="gray50") +
    plot.theme +
    theme(axis.text.x=element_text(angle=90, hjust=0, vjust=0.5, color=color.per.group),
          axis.text=element_text(size=font.size), axis.ticks=element_blank(), axis.title=element_blank()) +
    guides(fill=guide_colorbar(title=legend.title, title.position="left", title.theme=element_text(angle=90, hjust=0.5))) +
    scale_y_discrete(position="right", expand=c(0, 0)) +
    scale_x_discrete(expand=c(0, 0), position=x.axis.position) +
    theme_legend_position(legend.position) +
    theme(legend.key.width=legend.key.width, legend.background=element_blank()) +
    val2ggcol(df$value, palette=palette, color.range=color.range, return.fill=TRUE)

  if (symmetric) {
    gg <- gg + theme(axis.text.y=element_text(color=color.per.group))
  }

  return(gg)
}

#' Get Gene Scale
#' @param genes type of genes ("up", "down" or "all")
#' @param high color for the highest value (default: "gray80")
#' @return palette function
#' @export
getGenePalette <- function(genes=c("up", "down", "all"), high="gray80") {
  genes <- match.arg(genes)

  if (genes == "up") {
    low <- "red"
  } else if (genes == "down") {
    low <- "blue"
  } else {
    low <- "green"
  }

  return(colorRampPalette(c(high, low)))
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

getOntologyPlotTitle <- function(genes, cell.subgroup, type) {
  if(genes == "all")
    return(ggtitle(paste(cell.subgroup, type, "terms, all DE genes")))

  return(ggtitle(paste0(cell.subgroup, " ", type, " terms, ", genes,"-regulated DE genes")))
}

estimateMeanCI <- function(arr, quant=0.05, n.samples=500, ...) {
  if (length(arr) == 1) return(c(NA, NA))
  s.means <- sapply(1:n.samples, function(i) mean(sample(arr, replace=TRUE), ...))
  return(quantile(s.means, c(quant, 1 - quant)))
}

#' @title Plot bar, point or boxplots showing mean/median values per cell type
#' @description  Generic function for plotting mean or median values per cell type (used for expression shift distances and others)
#' @param df - data frame containing the results, including $value and $Type slots which will be summarized
#' @param type - type of a plot "bar" (default), "point" (mean + sd), or "box" for boxplot
#' @param show.jitter whether to show individual data points (default: FALSE)
#' @param jitter.alpha transparency value for the data points (default: 0.05)
#' @param notch - whether to show notches in the boxplot version (default=TRUE)
#' @param palette - cell type palette
#' @return A ggplot2 object
plotMeanMedValuesPerCellType <- function(df, pvalues=NULL, type=c('box', 'point', 'bar'), show.jitter=TRUE,
                                         notch=TRUE, jitter.alpha=0.05, palette=NULL, ylab='expression distance',
                                         yline=1, plot.theme=theme_get(), jitter.size=1, line.size=0.75, trim=0,
                                         order.x=TRUE, pvalue.y=NULL, y.max=NULL) {
  type <- match.arg(type)
  df$Type %<>% as.factor()

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

  if(type=='box') { # boxplot
    p <- ggplot(odf,aes(x=Type,y=value,fill=Type)) + geom_boxplot(notch=notch, outlier.shape=NA)
  } else if(type=='point') { # point + se
    p <- ggplot(df,aes(x=Type,y=mean,color=Type)) + geom_point(size=3) +
      geom_errorbar(aes(ymin=LI, ymax=UI), width=0.2, size=line.size)
    if(!is.null(palette)) {p <- p + scale_color_manual(values=palette)}
  } else { # barplot
    p <- ggplot(df,aes(x=Type,y=mean,fill=Type)) + geom_bar(stat='identity') +
      geom_errorbar(aes(ymin=LI, ymax=UI), width=0.2, size=line.size)
  }
  if(!is.na(yline) && !is.null(yline)) { p <- p + geom_hline(yintercept = yline, linetype=2, color='gray50') }
  p <- p +
    plot.theme +
    theme(
      axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=12),
      axis.text.y=element_text(angle=90, hjust=0.5, size=12),
      panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
      legend.position="none"
    ) +
    guides(fill="none") +
    labs(x="", y=ylab)

  if(show.jitter) {
    p <- p +
      geom_jitter(data=odf, aes(x=Type,y=value), color=1, position=position_jitter(0.1), show.legend=FALSE,
                  alpha=jitter.alpha, size=jitter.size)
  };

  if (!is.null(pvalues)) {
    pval.df <- pvalueToCode(pvalues) %>%
      tibble(Type=factor(names(.), levels=levels(df$Type)), pvalue=.) %>% na.omit()

    p <- p + geom_text(data=pval.df, mapping=aes(x=Type, label=pvalue), y=pvalue.y)
  }

  if(!is.null(palette)) {
    p <- p + scale_fill_manual(values=palette)
  }

  return(p)
}

##' show a scatter plot of cell-type values vs. number of cells per cell type
##'
##' @param df a data frame with $value and $Type columns, just like plotMeanValuesPerCellType
##' @param cell.groups a cell groups vector for calculating number of cells per cell type
##' @param show.whiskers whether se values should be plotted
##' @param palette cell type palette
##' @param ylab y axis label
##' @param yline value at which a horizontal reference value should be plotted
##' @return ggplot2 object
plotCellTypeSizeDep <- function(df, cell.groups, palette=NULL, font.size=4, ylab='expression distance', yline=1,
                                show.regression=TRUE, show.whiskers=TRUE, plot.theme=theme_get()) {
  cell.groups <- table(cell.groups) # %>% .[names(.) %in% names(de.raw)]

  # calculate mean, se and median
  odf <- na.omit(df); # full df is now in odf
  # calculate mean and se
  df$Type <- as.factor(df$Type)
  df <- data.frame(Type=levels(df$Type), mean=tapply(df$value,df$Type,mean), se=tapply(df$value, df$Type, function(x) sd(x)/sqrt(length(x))), stringsAsFactors=FALSE)
  df <- df[order(df$mean,decreasing=F),]
  df$Type <- factor(df$Type,levels=df$Type)
  df <- df[!is.na(df$mean),]
  df$size <- cell.groups[as.character(df$Type)]

  # order cell types according to the mean
  odf$Type <- factor(odf$Type,levels=df$Type)
  p <- ggplot(df,aes(x=size,y=mean,color=Type)) + geom_point(size=3)
  p <- p + ggrepel::geom_label_repel(aes(label=Type), size=font.size, min.segment.length=0.1, box.padding=0, label.size=0, max.iter=300, fill=NA)
  if(show.whiskers) p <- p+geom_errorbar(aes(ymin=mean-se*1.96, ymax=mean+se*1.96),width=0.2)
  if(!is.null(palette)) {p <- p+scale_color_manual(values=palette)}
  if(show.regression) {
    p <- p+geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="gray", size=0.5,linetype=2)
  }

  if(!is.na(yline)) { p <- p+ geom_hline(yintercept = 1,linetype=2,color='gray50') }
  p <- p +
    plot.theme +
    guides(fill="none")+
    theme(legend.position = "none")+
    labs(x="number of cells", y=ylab)
  if(!is.null(palette)) {
    p <- p + scale_fill_manual(values=palette)
  }

  return(p)
}

# TODO: Improve speed of this function. No need to check BP/MF/CC all the time
reduceEdges <- function(edges, verbose=TRUE, n.cores = 1) {
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

##' @title Plot ontology families
##' @description Plot related ontologies in one hierarchical network plot
##' @param fam List of ontology IDs for the chosen family
##' @param data Data frane if raw ontology data for the chosen cell type
##' @param plot.type How much of the family network should be plotted. Can be "complete" (entire network), "dense" (show 1 parent for each significant term), or "minimal" (only show significant terms) (default: complete)
##' @param show.ids Whether to show ontology IDs instead of names (default: F)
##' @param string.length Length of strings for wrapping in order to fit text within boxes (default: 18)
##' @param legend.label.size Size og legend labels (default: 1)
##' @param legend.position Position of legend (default: topright)
##' @param verbose Print messages (default: T)
##' @param n.cores Number of cores to use (default: 1)
##' @return Rgraphviz object
plotOntologyFamily <- function(fam, data, plot.type = "complete", show.ids=FALSE, string.length=18, legend.label.size = 1,
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

  # Convert IDs to names
  if(!show.ids) {
    for(id in 1:nrow(nodes)) {
      edges[edges == nodes$label[id]] <- nodes$name[id]
    }
  }

  # Wrap strings for readability
  edges.wrapped <- edges %>%
    apply(2, function(x) wrap_strings(x, string.length))

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


plotVolcano <- function(de.df, p.name='padj', color.var = 'CellFrac', legend.pos="none", palette=brewerPalette("RdYlBu"), lf.cutoff=1.5, p.cutoff=0.05,
                        cell.frac.cutoff=0.2, size=c(0.1, 1.0), lab.size=2, draw.connectors=TRUE, sel.labels=NULL, plot.theme=theme_get(), ...) {

  checkPackageInstalled("EnhancedVolcano", bioc=TRUE)
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

parseLimitRange <- function(lims, vals) {
  if (is.null(lims)) return(range(vals))
  if (!is.character(lims)) return(lims)
  lims[grep("%", lims)] %<>%
    sapply(function(q) {as.numeric(strsplit(q, "%")[[1]]) / 100}) %>%
    quantile(vals, .)
  return(as.numeric(lims))
}

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

pvalueToCode <- function(pvals) {
  symnum(pvals, corr=FALSE, na=FALSE, legend=FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
         symbols = c("***", "**", "*", "ns")) %>%
    as.character() %>% setNames(names(pvals))
}

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
