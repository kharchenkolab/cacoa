#' @import tibble
#' @import cowplot
#' @import dplyr
#' @import magrittr
#' @importFrom reshape2 melt
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

#' @title Plot raw DE genes
#' @description Plot number of DE genes as a function of number of cells
#' @param de.raw List of differentially expressed genes per cell type, results from getPerCellTypeDE (default: stored list)
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="none")
#' @param p.adjust.cutoff Adjusted P cutoff (default=0.05)
#' @return A ggplot2 object
#' @export
plotDEGenes <- function(de.raw, cell.groups, legend.position="none", p.adjust.cutoff = 0.05, label = T, palette=NULL, ...) {
  cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.raw)]

  sapply(de.raw, function(n) n %>% dplyr::filter(padj <= p.adjust.cutoff) %>% nrow) %>%
    plotNCellRegression(cell.groups, x.lab="Number of cells", y.lab="Significant DE genes", legend.position=legend.position, label=label, ...) +
    geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="black", size=0.5)
}

#' @title Plot ontology distribution
#' @description Bar plot of ontologies per cell type
#' @param type Ontology, must be either "GO" or "DO" (default=NULL)
#' @param ont.res Ontology resuls from estimateOntology
#' @return A ggplot2 object
#' @export
plotOntologyDistribution <- function(type = NULL, ont.res, cell.groups) {
  if(type=="GO") {
    p_df <- table(ont.res$Group, ont.res$Type, ont.res$direction) %>%
      as.data.frame() %>%
      setNames(c("Group","Type","direction","N")) %>%
      dplyr::arrange(Group)

    if(length(unique(p_df$direction)) > 1) {
      gg <- ggplot(p_df, aes(x=Group, y=N, fill=Type, group=Group)) +
        geom_bar(stat="identity") +
        facet_grid(~direction, switch="x")
    } else {
      gg <- ggplot(p_df) +
        geom_bar(aes(x=Group, y=N, fill=Type), stat="identity")
    }
  } else if(type=="DO") {
    p_df <- table(ont.res$Group, ont.res$direction) %>%
      as.data.frame() %>%
      setNames(c("Group","direction","N")) %>%
      dplyr::arrange(Group)

    if(length(unique(p_df$direction)) > 1) {
      gg <- ggplot(p_df) +
        geom_bar(aes(x=Group, y=N, fill=direction), stat="identity", position="dodge") +
        labs(fill="Gene set")
    } else {
      gg <- ggplot(p_df) +
        geom_bar(aes(x=Group, y=N), stat="identity")
    }
  } else {
    stop("'type' must be either 'GO' or 'DO'.")
  }

  gg +
    scale_y_continuous(expand=c(0, 0)) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          legend.position="right") +
    labs(x="", y=paste0("No. of ",type," terms"))
}

#' @title Plot ontology terms
#' @description Plot ontology terms as a function of both number of DE genes, and number of cells.
#' @param ont.res Ontology results in a data frame
#' @param type Ontology, must be either "GO" or "DO" (default=NULL)
#' @param de.filter List of filtered differentially expressed genes, results from prepareOntologyData (default: stored list)
#' @param cell.groups Vector indicating cell group sizes with cell group names
#' @param label.x.pos Plot label position on x axis (default=0.01)
#' @param label.y.pos Plot label position on y axis (default=1)
#' @param scale Scaling of plots, adjust if e.g. label is misplaced. See cowplot::plot_grid for more info (default=0.93)
#' @return A ggplot2 object
#' @export
plotOntologyTerms <- function(type=NULL, ont.res, de.filter, cell.groups, label.x.pos=0.01, label.y.pos=1, scale = 0.93) {
  if(length(unique(ont.res$Type))==1) stop("The input only contains one cell type.")

  if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

  cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.filter)]

  pg <- cowplot::plot_grid(
    ont.res$Group %>% table %>% c %>%
      plotNCellRegression(sapply(de.filter, length), x.lab="Number of highly-expressed DE genes",
                          y.lab=NULL, legend.position="none", label=T) +
      geom_smooth(method=MASS::rlm, formula = y~x, se=F, color="black", size=0.5),
    ont.res$Group %>% table %>% c %>%
      plotNCellRegression(cell.groups, y.lab=NULL, legend.position="none", label=T) +
      geom_smooth(method=MASS::rlm, formula = y~x, se=F, color="black", size=0.5),
    ncol=1, labels=c("a", "b"), label_x = label.x.pos, label_y = label.y.pos
  )

  if(type=="GO") {
    y_lab <- "Number of GO terms"
  } else if(type=="DO") {
    y_lab <- "Number of DO terms"
  }

  cowplot::plot_grid(pg, nrow=1, scale=scale) + draw_label(y_lab, x=0, y=0.6, vjust= 1.5, angle=90)
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

# TODO should depend on merged DF, not list
#' @title Plot ontology correlations
#' @description Plot correlation matrix for ontologies between cell types
#' @param ont.res Data frame with ontology resuls from estimateOntology
#' @param type Ontology, must be either "GO" or "DO" (default=NULL)
#' @param genes Specify which genes are plotted, can either be 'down', 'up' or 'all' (default=NULL)
#' @return A ggplot2 object
#' @export
plotOntologySimilarities <- function(type=NULL, ont.res, genes = NULL) {
  if(type=="GO") {
    pathway_df <- unique(ont.res$Group) %>%
      lapply(function(cell.group) {
        lapply(ont.res %>%
                 dplyr::filter(Group == cell.group) %>%
                 dplyr::pull(Type) %>%
                 as.factor() %>%
                 levels(), function(go) {
                   tibble::tibble(Pathway=ont.res %>%
                                    dplyr::filter(Group == cell.group) %>%
                                    dplyr::filter(Type==go) %>%
                                    dplyr::pull(Description),
                                  Group=cell.group,
                                  GO=go)
                   }) %>% dplyr::bind_rows()
        }) %>%
      dplyr::bind_rows()
    } else if(type=="DO") {
      pathway_df <- unique(ont.res$Group) %>%
        lapply(function(cell.group) {
          tibble::tibble(Pathway=ont.res %>%
                           dplyr::filter(Group == cell.group) %>%
                           dplyr::pull(Description),
                         Group=cell.group)
          }) %>%
        dplyr::bind_rows()
  } else {
    stop("'type' must be either 'GO' or 'DO'.")
  }
  path_bin <- pathway_df %>%
    dplyr::select(Pathway, Group) %>%
    dplyr::mutate(X=1) %>%
    tidyr::spread(Pathway, X) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$Group) %>%
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

  if(genes == "all") {
    l <- ggtitle(paste0(type," term similarities for all DE genes"))
  } else {
    l <- ggtitle(paste0(type," term similarities for ",genes,"-regulated DE genes"))
  }

  plotHeatmap(p_mat, color.per.group=NULL, row.order=t_order, col.order=rev(t_order), legend.title="Similarity", color.range=c(0, 0.5)) +
    scale_fill_distiller(palette="RdYlBu") +
    geom_vline(aes(xintercept=x), data.frame(x=cumsum(t_cl_lengths)[t_cl_lengths > 1] + 0.5)) +
    geom_hline(aes(yintercept=x), data.frame(x=cumsum(t_cl_lengths)[t_cl_lengths > 1] + 0.5)) + l
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


#' @title Plot proportions
#' @description Plot the cell group proportions per sample
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param sample.per.cell Vector indicating sample name with cell names (default: stored vector)
#' @param sample.groups Vector indicating sample groups with sample names (default: stored vector)
#' @return A ggplot2 object
plotCellNumbers <- function(legend.position = "right", cell.groups, sample.per.cell, sample.groups) {
  df.melt <- data.frame(anno=cell.groups, group=sample.per.cell[match(names(cell.groups), names(sample.per.cell))]) %>%
    table %>%
    rbind %>%
    t %>%
    as.data.frame %>%
    dplyr::mutate(group = sample.groups[match(levels(sample.per.cell), names(sample.groups))]) %>%
    reshape2::melt(., id.vars="group")

  ggplot(df.melt, aes(x=variable, y=value, by=group)) +
    geom_boxplot(position=position_dodge(), outlier.shape = NA) +
    ylab("Cells per sample") +
    xlab("") +
    theme_bw() +
    theme_legend_position(legend.position) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
          legend.title=element_blank()) +
    geom_point(position=position_jitterdodge(jitter.width=0.15), aes(col=group), alpha=0.4) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, (max(df.melt$value) + 50)))
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

getOntologyPlotTitle <- function(genes, cell.subgroup, type) {
  if(genes == "all")
    return(ggtitle(paste(cell.subgroup, type, "terms, all DE genes")))

  return(ggtitle(paste0(cell.subgroup, " ", type, " terms, ", genes,"-regulated DE genes")))
}


#' @title Plot bar, point or boxplots showing mean/median values per cell type
#' @description  Generic function for plotting mean or median values per cell type (used for expression shift distances and others)
#' @param df - data frame containing the results, including $val and $cell slots which will be summarized
#' @param type - type of a plot "bar" (default), "point" (mean + sd), or "box" for boxplot
#' @param show.jitter whether to show indiivudal data points (default: FALSE)
#' @param jitter.alpha transparency value for the data points (default: 0.05)
#' @param notch - whether to show notches in the boxplot version (default=TRUE)
#' @param palette - cell type palette
#' @return A ggplot2 object
plotMeanMedValuesPerCellType <- function(df, type='bar', show.jitter=TRUE, notch = T, jitter.alpha=0.05, palette=NULL, ylab='expression distance', yline=1) {

  # calculate mean, se and median
  odf <- na.omit(df); # full df is now in odf
  # calculate mean and se
  df$cell <- as.factor(df$cell)
  df <- data.frame(cell=levels(df$cell), mean=tapply(df$val,df$cell,mean), se=tapply(df$val,df$cell, function(x) sd(x)/sqrt(length(x))), stringsAsFactors=FALSE)
  df <- df[order(df$mean,decreasing=F),]
  df$cell <- factor(df$cell,levels=df$cell)
  df <- df[!is.na(df$mean),]

  # order cell types according to the mean
  odf$cell <- factor(odf$cell,levels=df$cell)
  
  if(type=='box') { # boxplot
    p <- ggplot(odf,aes(x=cell,y=val,fill=cell)) + geom_boxplot(notch=notch, outlier.shape=NA)
  } else if(type=='point') { # point + se
    p <- ggplot(df,aes(x=cell,y=mean,color=cell)) + geom_point(size=3)+ geom_errorbar(aes(ymin=mean-se*1.96, ymax=mean+se*1.96),width=0.2)
    if(!is.null(palette)) {p <- p+scale_color_manual(values=palette)}
  } else { # default to barplot
    p <- ggplot(df,aes(x=cell,y=mean,fill=cell)) + geom_bar(stat='identity')+ geom_errorbar(aes(ymin=mean-se*1.96, ymax=mean+se*1.96),width=0.2)
  }
  if(!is.na(yline)) { p <- p+ geom_hline(yintercept = 1,linetype=2,color='gray50') }
  p <- p+ theme_bw() +
    theme(axis.text.x=element_text(angle = 90, hjust=1, size=12), axis.text.y=element_text(angle=90, hjust=0.5, size=12))+ guides(fill=FALSE)+
    theme(legend.position = "none")+
    labs(x="", y=ylab)
  if(show.jitter) p <- p+geom_jitter(data=odf,aes(x=cell,y=val),color=1, position=position_jitter(0.1),show.legend=FALSE,alpha=jitter.alpha);
  if(!is.null(palette)) {
    p <- p+ scale_fill_manual(values=palette)
  }
  p
  
}

##' show a scatter plot of cell-type values vs. number of cells per cell type
##'
##' @param df a data frame with $val and $cell columns, just like plotMeanValuesPerCellType
##' @param cell.groups a cell groups vector for calculating number of cells per cell type
##' @param show.whiskers whether se values should be plotted
##' @param palette cell type palette
##' @param ylab y axis label
##' @param yline value at which a horizontal reference value should be plotted
##' @return ggplot2 object
plotCellTypeSizeDep <- function(df, cell.groups, palette=NULL, font.size=4, ylab='expression distance', yline=1, show.regression=TRUE, show.whiskers=TRUE) {
  cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.raw)]
  
  # calculate mean, se and median
  odf <- na.omit(df); # full df is now in odf
  # calculate mean and se
  df$cell <- as.factor(df$cell)
  df <- data.frame(cell=levels(df$cell), mean=tapply(df$val,df$cell,mean), se=tapply(df$val,df$cell, function(x) sd(x)/sqrt(length(x))), stringsAsFactors=FALSE)
  df <- df[order(df$mean,decreasing=F),]
  df$cell <- factor(df$cell,levels=df$cell)
  df <- df[!is.na(df$mean),]
  df$size <- cell.groups[as.character(df$cell)]

  # order cell types according to the mean
  odf$cell <- factor(odf$cell,levels=df$cell)
  p <- ggplot(df,aes(x=size,y=mean,color=cell)) + geom_point(size=3)
  p <- p+ggrepel::geom_label_repel(aes(label=cell), size=font.size, min.segment.length=0.1, box.padding=0, label.size=0, max.iter=300, fill=NA)
  if(show.whiskers) p <- p+geom_errorbar(aes(ymin=mean-se*1.96, ymax=mean+se*1.96),width=0.2)
  if(!is.null(palette)) {p <- p+scale_color_manual(values=palette)}
  if(show.regression) {
    p <- p+geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="gray", size=0.5,linetype=2)
  }

  if(!is.na(yline)) { p <- p+ geom_hline(yintercept = 1,linetype=2,color='gray50') }
  p <- p+ theme_bw() +
    guides(fill=FALSE)+
    theme(legend.position = "none")+
    labs(x="number of cells", y=ylab)
  if(!is.null(palette)) {
    p <- p+ scale_fill_manual(values=palette)
  }
  
  p
  
}
