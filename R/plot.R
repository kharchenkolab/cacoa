#' @import ggrepel
#' @import tibble
#' @import cowplot
#' @import dplyr
#' @import magrittr
#' @importFrom reshape2 melt
NULL

plotNCellRegression <- function(n, n.total, y.lab="N", legend.position="right", label=T) {
  p.df <- data.frame(N=n) %>% tibble::as_tibble(rownames="Type") %>%
    mutate(NCells=n.total[Type])

  gg <- ggplot(p.df, aes(x=NCells, y=N)) +
    geom_point(aes(color=Type)) +
    scale_x_log10() +
    ylim(0, max(p.df$N)) + labs(x="Number of cells", y=y.lab) +
    theme(legend.position=legend.position, legend.justification=legend.position, legend.background=element_rect(fill=alpha("white", 0.4))) +
    guides(color=guide_legend(title="Cell type"))

  if(label) gg <- gg + geom_label_repel(aes(label=Type), size=2, min.segment.length=0.1, box.padding=0, label.size=0, max.iter=300, fill=alpha("white", 0.4))

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
plotDEGenes <- function(de.raw, cell.groups, legend.position="none", p.adjust.cutoff = 0.05, label = T) {
  cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.raw)]

  sapply(de.raw, function(n) n %>% dplyr::filter(padj <= p.adjust.cutoff) %>% nrow) %>%
    plotNCellRegression(cell.groups, y.lab="Significant DE genes", legend.position=legend.position, label=label) +
    xlab("Number of cells") +
    geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="black", size=0.5) +
    theme_bw()
}

#' @title Plot filtered DE genes
#' @description Plot number of DE genes as a function of number of cells
#' @param de.filter List of filtered differentially expressed genes, results from prepareOntologyData (default: stored list)
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="bottom")
#' @param label Show labels on plot (default=T)
#' @export
plotFilteredDEGenes <- function(de.filter, cell.groups, legend.position="bottom", label = T) {
  cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.filter)]

  sapply(de.filter, length) %>%
    plotNCellRegression(cell.groups, y.lab="Highly-expressed DE genes", legend.position=legend.position, label=label) +
    xlab("Number of cells") +
    geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="black", size=0.5) +
    theme_bw()
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
#' @param rel_heights Vector indicating relative heights for plots. Only relevant if show.legend=T. See cowplot::plot_grid for more info (default=c(2.5, 0.5))
#' @param scale Scaling of plots, adjust if e.g. label is misplaced. See cowplot::plot_grid for more info (default=0.93)
#' @return A ggplot2 object
#' @export
plotOntologyTerms <- function(type=NULL, ont.res, de.filter, cell.groups, label.x.pos=0.01, label.y.pos=1, rel_heights = c(2.5, 0.5), scale = 0.93) {
  if(length(unique(ont.res$Type))==1) stop("The input only contains one cell type.")

  if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

  cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.filter)]

  if(legend.position != "none") leg <- cowplot::get_legend(plotNCellRegression(sapply(de.filter, length), cell.groups, legend.position=legend.position))

  pg <- cowplot::plot_grid(
    ont.res$Group %>% table %>% c %>%
      plotNCellRegression(sapply(de.filter, length), y.lab=NULL, legend.position="none", label=T) +
      geom_smooth(method=MASS::rlm, formula = y~x, se=F, color="black", size=0.5) +
      scale_x_continuous(name="Number of highly-expressed DE genes") +
      theme_bw(),
    ont.res$Group %>% table %>% c %>%
      plotNCellRegression(cell.groups, y.lab=NULL, legend.position="none", label=T) +
      geom_smooth(method=MASS::rlm, formula = y~x, se=F, color="black", size=0.5) +
      theme_bw(),
    ncol=1, labels=c("a", "b"), label_x = label.x.pos, label_y = label.y.pos
  )

  if(type=="GO") {
    y_lab <- "Number of GO terms"
  } else if(type=="DO") {
    y_lab <- "Number of DO terms"
  }
  if(legend.position != "none") {
    p <- cowplot::plot_grid(pg, leg, nrow=2, rel_heights = rel_heights, scale=scale) + cowplot::draw_label(y_lab, x=0, y=0.6, vjust= 1.5, angle=90)
  } else {
    p <- cowplot::plot_grid(pg, nrow=1, scale=scale) + draw_label(y_lab, x=0, y=0.6, vjust= 1.5, angle=90)
  }
  return(p)
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
    theme(legend.position=legend.position, legend.key.width=legend.key.width, legend.background=element_blank())

  return(gg)
}

#' @title Plot ontology heatmap
#' @description Plot a heatmap of ontology P values per cell type
#' @param type Ontology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default=NULL)
#' @param ont.res Ontology resuls from estimateOntology
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="left")
#' @param selection Order of rows in heatmap. Can be 'unique' (only show terms that are unique for any cell type); 'common' (only show terms that are not unique for any cell type); 'all' (all ontology terms) (default="all")
#' @param n Number of terms to show (default=10)
#' @export
plotOntologyHeatmap <- function(type = "GO", ont.res, legend.position = "left", selection = "all", n = 20, cell.subgroups = NULL) {
  if(is.null(selection) || (!selection %in% c("unique","common","all"))) stop("'selection' must be one of the following: 'unique', 'common', or 'all'.")

  if(!is.null(cell.subgroups) && (length(cell.subgroups) == 1)) stop("'cell.subgroups' must contain at least two groups. Please use plotOntologyBarplot or plotOntologyDotplot instead.")

  if(type=="GO") {
    ont.sum <- getOntologySummary(ont.res)
  } else if(type=="BP" | type=="CC" | type=="MF") {
    ont.sum <- getOntologySummary(ont.res %>% filter(Type==type))
  } else if(type=="DO") {
    ont.sum <- getOntologySummary(ont.res)
  }

  if(!is.null(cell.subgroups)) {
      ont.sum %<>% dplyr::select(all_of(cell.subgroups))
    }

  if(selection=="unique") {
      ont.sum %<>%
        .[rowSums(. > 0) == 1,]
  } else if(selection=="common") {
      ont.sum %<>%
        .[rowSums(. > 0) > 1,]
  }

  if(nrow(ont.sum) == 0) stop("Nothing to plot. Try another selection.")

  ont.sum %>%
    .[, colSums(.) > 0] %>%
    .[match(rowSums(.)[rowSums(.)>0] %>%
              .[order(., decreasing = F)] %>%
              names, rownames(.)),] %>%
    tail(n) %>%
    plotHeatmap(legend.position=legend.position, row.order=T)
}

# TODO should depend on merged DF, not list
#' @title Plot ontology correlations
#' @description Plot correlation matrix for ontologies between cell types
#' @param ont.res List with ontology resuls from estimateOntology
#' @param type Ontology, must be either "GO" or "DO" (default=NULL)
#' @return A ggplot2 object
#' @export
plotOntologySimilarities <- function(type=NULL, ont.res) {
  if(type=="GO") {
    pathway_df <- names(ont.res) %>%
      lapply(function(cell.group) lapply(ont.res[[cell.group]] %>% dplyr::pull(Type) %>% as.factor() %>% levels(), function(go) {
        tibble::tibble(Pathway=ont.res[[cell.group]] %>% dplyr::filter(Type==go) %>% dplyr::pull(Description),
                       Group=cell.group,
                       GO=go)
      }) %>% dplyr::bind_rows()) %>%
      dplyr::bind_rows()
    } else if(type=="DO") {
    pathway_df <- lapply(ont.res %>% names, function(cell.group) {
      tibble::tibble(Pathway=ont.res[[cell.group]] %>% dplyr::pull(Description),
                     Group=cell.group)
    }) %>% dplyr::bind_rows()
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
  plotHeatmap(p_mat, color.per.group=NULL, row.order=t_order, col.order=rev(t_order), legend.title="Similarity") +
    scale_fill_distiller(palette="RdYlBu", limits=c(0, 0.5)) +
    geom_vline(aes(xintercept=x), data.frame(x=cumsum(t_cl_lengths)[t_cl_lengths > 1] + 0.5)) +
    geom_hline(aes(yintercept=x), data.frame(x=cumsum(t_cl_lengths)[t_cl_lengths > 1] + 0.5))
}

#' @title Plot proportions
#' @description Plot the cell group proportions per sample
#' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param sample.per.cell Vector indicating sample name with cell names (default: stored vector)
#' @param sample.groups Vector indicating sample groups with sample names (default: stored vector)
#' @return A ggplot2 object
plotProportions <- function(legend.position = "right", cell.groups, sample.per.cell, sample.groups) {
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

  ggplot(df.melt, aes(x=variable, y=value, by=group)) +
    geom_boxplot(position=position_dodge(), outlier.shape = NA) +
    ylab("% cells per sample") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90),
          legend.position=legend.position,
          legend.title=element_blank()) +
    geom_point(position=position_jitterdodge(jitter.width=0.15), aes(col=group)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, (max(df.melt$value) + 1)))
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
    theme_classic() +
    theme(axis.text.x = element_text(angle=90),
          legend.position=legend.position,
          legend.title=element_blank()) +
    geom_point(position=position_jitterdodge(jitter.width=0.15), aes(col=group)) +
    scale_y_continuous(expand=c(0, 0), limits=c(0, (max(df.melt$value) + 50)))
}

#' @title Plot ontologies with barplot
#' @description Plot a barplot of ontologies with adj. P values for a specific cell subgroup
#' @param ont.res Data frame with ontology resuls from estimateOntology
#' @param genes Specify which genes are plotted, can either be 'down', 'up' or 'all' (default=NULL)
#' @param n Number of ontologies to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
#' @param p.adj Adjusted P cutoff (default=0.05)
#' @return A ggplot2 object
plotOntologyBarplot <- function(ont.res, genes = NULL, type = NULL, n = 20, p.adj = 0.05) {
  ont.res %<>% dplyr::mutate(., gratio = sapply(.$GeneRatio, function(s) strsplit(s, "/")) %>% sapply(function(x) as.numeric(x[1])/as.numeric(x[2]) * 100)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::filter(p.adjust <= p.adj) %>%
    {if(nrow(.) > n) .[1:n,] else .}
  ont.res$Description %<>% as.factor()

  if(is.null(type)) type <- "Ontology"

  gg <- ggplot(ont.res, aes(reorder(Description, -p.adjust), gratio, fill=p.adjust)) +
    geom_col() +
    coord_flip() +
    labs(title=paste0(type," terms"), y="% DE genes of total genes per pathway", x="", fill="Adj. P") +
    theme_bw() +
    scale_y_continuous(expand=c(0, 0), limits=c(0, (max(ont.res$gratio) + 1)))

  if(is.null(genes)) genes <- "all"

  if(genes == "up") {
    gg <- gg + scale_fill_gradient(low = "red", high = "gray80")
  } else if(genes == "down") {
    gg <- gg + scale_fill_gradient(low = "blue", high = "gray80")
  } else {
    gg <- gg + scale_fill_gradient(low = "green", high = "gray80")
  }

  gg
}

#' @title Plot ontologies with dotplot
#' @description Plot a dotplot of ontologies with adj. P values for a specific cell subgroup
#' @param ont.res Data frame with ontology resuls from estimateOntology
#' @param genes Specify which genes are plotted, can either be 'down', 'up' or 'all' (default=NULL)
#' @param n Number of ontologies to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
#' @param p.adj Adjusted P cutoff (default=0.05)
#' @return A ggplot2 object
plotOntologyDotplot <- function(ont.res, genes = NULL, type = NULL, n = 20, p.adj = 0.05) {
  ont.res %<>% dplyr::mutate(., gratio = sapply(.$GeneRatio, function(s) strsplit(s, "/")) %>% sapply(function(x) as.numeric(x[1])/as.numeric(x[2]) * 100)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::filter(p.adjust <= p.adj) %>%
    {if(nrow(.) > n) .[1:n,] else .}
  ont.res$Description %<>% as.factor()

  if(is.null(type)) type <- "Ontology"

  gg <- ggplot(ont.res, aes(reorder(Description, -p.adjust), gratio, col=p.adjust)) +
    geom_point(aes(size = Count)) +
    coord_flip() +
    labs(title=paste0(type," terms"), y="% DE genes of total genes per pathway", x="", col="Adj. P", size = "DE genes") +
    theme_bw() +
    scale_y_continuous(expand=c(0, 0), limits=c(0, (max(ont.res$gratio) + 1)))

  if(is.null(genes)) genes <- "all"

  if(genes == "up") {
    gg <- gg + scale_color_gradient(low = "red", high = "gray80")
  } else if(genes == "down") {
    gg <- gg + scale_color_gradient(low = "blue", high = "gray80")
  } else {
    gg <- gg + scale_color_gradient(low = "green", high = "gray80")
  }

  gg
}
