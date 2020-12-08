#' Estimate Gene Programmes
#' @param n.programmes maximal number of gene programmes to find (parameter `p` for fabia).
#' @param n.sampled.cells number of sub-sampled cells for estimating the gene programmes. If 0, all cells are used.
#' it is interpreted as a vector of cell
#' @inheritDotParams fabia::fabia -p -X -cyc -alpha -random
estimateGeneProgrammes <- function(z.scores, n.programmes, n.sampled.cells=15000, min.z=0.5, max.z=3,
                                   cyc=1500, alpha=0.2, random=-1, ...) {
  if (!requireNamespace("fabia", quietly=TRUE))
    stop("fabia package must be installed to run this function")

  z.scores@x[is.na(z.scores@x)] <- 0
  if (min.z > 1e-10) {
    z.scores.bin <- z.scores
    z.scores.bin@x %<>% {as.numeric(abs(.) >= min.z)}
    genes.filt <- which(colMeans(z.scores.bin) > 0.01)
    z.scores <- z.scores[, genes.filt]
  }

  if (n.sampled.cells > 0) {
    sample.ids <- unique(round(seq(1, nrow(z.scores), length.out=n.sampled.cells)))
    z.scores <- z.scores[sample.ids,]
  }

  z.scores@x %<>% pmax(-max.z) %<>% pmin(max.z)

  fabia.res <- z.scores %>% t() %>%
    fabia::fabia(p=n.programmes, alpha=alpha, cyc=cyc, random=random, ...)

  return(list(fabia=fabia.res, sample.ids=sample.ids))
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
