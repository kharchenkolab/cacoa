#' Estimate Local Z Scores
#'
#' @param graph Alignment graph adjacency matrix
#' @param genes Genes to be tested. Procedure is slow, so it's reasonable to restrict amount of genes.
#' @param is.ref Boolean vector per cell indicating whether it came from the reference condition
#' @param count.matrix.transposed Joint count matrix with cells by rows and genes by columns
#' @param ref.level Reference condition level, e.g., wt, ctrl, or healthy
estimateClusterFreeZScores <- function(graph, count.matrix.transposed, is.ref, genes=NULL, max.z=20, ...) {
  adj.mat <- igraph::as_adj(graph)

  cell.names <- intersect(rownames(count.matrix.transposed), rownames(adj.mat)) %>%
    intersect(names(is.ref))
  if (length(cell.names) == 0)
    stop("No cells in graph matching to count.matrix.transposed and is.ref")

  if (!is.null(genes)) {
    count.matrix.transposed <- count.matrix.transposed[,genes]
  }

  z.mat <- clusterFreeZScoreMat(adj.mat[cell.names, cell.names], count.matrix.transposed[cell.names,],
                                is.ref[cell.names], ...)

  z.mat@x %<>% pmin(max.z) %>% pmax(-max.z)

  return(z.mat)
}

### Selection

#' Get Top DE Genes
#'
#' @param z.scores Filtered z scores
#' @param min.z (default=0)
#' @param max.z (default=10)
#' @param cell.subset Cells to subset from z.scores (default=NULL)
#' @param top.quantile (default=NULL)
getTopDEGenes <- function(z.scores, min.z=0, max.z=10, cell.subset=NULL, top.quantile=NULL) {
  if (!is.null(cell.subset)) {
    z.scores %<>% .[intersect(rownames(.), cell.subset), ]
  }

  z.scores@x %<>% abs() %>% pmin(max.z)
  z.scores@x[z.scores@x < min.z] <- 0
  z.scores %<>% Matrix::drop0()

  if (!is.null(top.quantile)) {
    z.scores[Matrix::t(Matrix::t(z.scores) < apply(z.scores, 2, quantile, top.quantile, na.rm=TRUE))] <- 0
  }

  return(Matrix::colMeans(z.scores, na.rm=TRUE) %>% sort(decreasing=TRUE))
}

plotZScores <- function(gene, con, z.scores, cur.scores=NULL, color.range=c(-4, 4), show.legend=TRUE, legend.pos=c(1, 1), size=0.2, score.digits=3, ...) {
  cur.score <- if (is.null(cur.scores)) sum(z.scores[, gene], na.rm=TRUE) else cur.scores[gene]
  colors <- z.scores[,gene] %>% pmin(max(color.range)) %>% pmax(min(color.range))
  con$plotGraph(colors=colors, show.legend=show.legend, legend.pos=legend.pos, size=size,
                color.range=color.range, title=paste0(gene, ": ", round(cur.score, score.digits)), ...)
}

plotZScoreList <- function(con, z.scores, scores, n.genes=NULL, genes=NULL, ...) {
  if (is.null(n.genes) == is.null(genes))
    stop("Either n.genes or genes must be provided")

  if (is.null(genes)) {
    genes <- 1:n.genes
  }

  plots <- scores[genes] %>% names() %>%
    lapply(plotZScores, con, z.scores, cur.scores=scores, ...)

  return(plots)
}

#' Plot Gene Comparison Between Conditions
#'
#' @param genes Vector of genes to plot
#' @param con Conos object
#' @param condition.per.cell Named factor with cell names and condition
#' @param z.scores Z scores matrix. It is recommended to filter and adjust scores before plotting (default=NULL)
#' @param cur.scores Named numeric vector with gene names.
#' @param show.legend Plot legend (default=TRUE)
#' @param legend.pos Legend position, see ggplot2::theme
#' @param size Size of cells on plot
#' @param n.col Columns of plots. If NULL, will be number of conditions + 1 (default=NULL)
#' @param ... Plotting variables propagated to conos:::embeddingPlot
plotGeneComparisonBetweenCondition <- function(genes, con, condition.per.cell, z.scores=NULL, cur.scores=NULL, show.legend=TRUE, legend.pos=c(1, 1), size=0.2, n.col=NULL, z.lims=NULL, max.expr=NULL, adj.list=NULL, ...) {
  #TODO: Condition.per.cell should be derived from sample.groups and cell.groups
  if (!is.null(cur.scores) & !is.null(names(cur.scores))) {
    genes <- names(sort(cur.scores, decreasing=TRUE)[genes])
  }

  if (is.null(n.col)) n.col <- length(unique(condition.per.cell)) + 1

  lapply(genes, function(g) {
    lst <- lapply(unique(condition.per.cell), function(sg) {
      if (is.null(max.expr)) {
        m.expr <- sapply(con$samples, function(s) max(conos:::getGeneExpression(s, g), na.rm=TRUE)) %>% max
      } else {
        m.expr <- max.expr
      }
      con$plotGraph(gene=g, show.legend=show.legend, legend.pos=legend.pos, size=size,
                    title=paste(sg," ",g), groups=condition.per.cell, subgroups=sg,
                    color.range=c(0, m.expr), legend.title="Expression", ...)
    })

    if (!is.null(z.scores)) {
      lst <- plotZScores(g, con, z.scores, cur.scores=cur.scores, show.legend=show.legend, legend.title="Z-score",
                         legend.pos=legend.pos, size=size, color.range=z.lims, ...) %>% list() %>% c(lst)
    }

    if (!is.null(adj.list)) {
      lst %<>% lapply(function(p) Reduce(`+`, adj.list, p))
    }

    cowplot::plot_grid(plotlist=lst, ncol=n.col)
  })
}
