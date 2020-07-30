##' Local Z Scores
##'
##' @param graph Alignment graph (embedding)
##' @param condition.per.cell Named group factor with cell names. Must have exactly two levels.
##' @param genes Genes to be tested. Procedure is slow, so it's reasonable to restrict amount of genes.
##' @param count.matrix.transposed Joint count matrix with cells by rows and genes by columns
##' @param ref.level Reference condition level, e.g., wt, ctrl, or healthy
localZScores <- function(graph, condition.per.cell, count.matrix.transposed, ref.level, genes=NULL, max.z=20, ...) {
  #TODO: Condition.per.cell should be derived from sample.groups and cell.groups
  if (length(unique(condition.per.cell)) != 2)
    stop("Exactly two levels must be provided in 'condition.per.cell'")

  if (!(ref.level %in% unique(condition.per.cell)))
    stop("'ref.level' not present in 'condition.per.cell'")

  adj.mat <- igraph::as_adj(graph)

  cell.names <- intersect(rownames(count.matrix.transposed), names(condition.per.cell))
  if (length(cell.names) == 0)
    stop("No cells in condition.per.cell matching to count.matrix.transposed")

  if (!is.null(genes)) {
    count.matrix.transposed <- count.matrix.transposed[,genes]
  }

  z.scores <- localZScoreMat(adj.mat[cell.names, cell.names], count.matrix.transposed[cell.names,],
                             condition.per.cell[cell.names] == ref.level, ...)

  z.scores@x %<>% pmin(max.z) %>% pmax(-max.z)

  return(z.scores)
}

### Selection

##' Get Top DE Genes
##'
##' @param z.scores Filtered z scores
##' @param min.z (default=0)
##' @param max.z (default=10)
##' @param cell.subset Cells to subset from z.scores (default=NULL)
##' @param top.quantile (default=NULL)
getTopDEGenes <- function(z.scores, min.z=0, max.z=10, cell.subset=NULL, top.quantile=NULL) {
  if (!is.null(cell.subset)) {
    z.scores %<>% .[intersect(rownames(.), cell.subset), ]
  }

  z.scores %<>% abs()
  z.scores@x[z.scores@x < min.z] <- 0
  z.scores %<>% Matrix::drop0() %>% pmin(max.z)

  if (!is.null(top.quantile)) {
    z.scores[Matrix::t(Matrix::t(z.scores) < apply(z.scores, 2, quantile, top.quantile))] <- 0
  }

  return(Matrix::colMeans(z.scores) %>% sort(decreasing=T))
}

plotZScores <- function(gene, con, z.scores, cur.scores=NULL, color.range=c(-4, 4), show.legend=T, legend.pos=c(1, 1), size=0.2, score.digits=3, ...) {
  cur.score <- if (is.null(cur.scores)) sum(z.scores[, gene], na.rm=T) else cur.scores[gene]
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

##' Plot Gene Comparison Between Conditions
##'
##' @param genes Vector of genes to plot
##' @param con Conos object
##' @param condition.per.cell Named factor with cell names and condition
##' @param z.scores Z scores matrix. It is recommended to filter and adjust scores before plotting (default=NULL)
##' @param cur.scores Named numeric vector with gene names.
##' @param show.legend Plot legend (default=T)
##' @param legend.pos Legend position, see ggplot2::theme
##' @param size Size of cells on plot
##' @param n.col Columns of plots. If NULL, will be number of conditions + 1 (default=NULL)
##' @param ... Plotting variables propagated to conos:::embeddingPlot
plotGeneComparisonBetweenCondition <- function(genes, con, condition.per.cell, z.scores=NULL, cur.scores=NULL, show.legend=T, legend.pos=c(1, 1), size=0.2, n.col=NULL, ...) {
  #TODO: Condition.per.cell should be derived from sample.groups and cell.groups
  if (!is.null(cur.scores) & !is.null(names(cur.scores))) {
    genes <- names(sort(cur.scores, decreasing=T)[genes])
  }

  if (is.null(n.col)) n.col <- length(unique(condition.per.cell)) + 1

  lapply(genes, function(g) {
    lst <- lapply(unique(condition.per.cell), function(sg) {
      con$plotGraph(gene=g, show.legend=show.legend, legend.pos=legend.pos, size=size,
                      title=paste(sg," ",g), groups=condition.per.cell, subgroups=sg, ...)
    })

      if (!is.null(z.scores)) {
        lst <- plotZScores(g, con, z.scores, cur.scores=cur.scores, show.legend=show.legend,
                           legend.pos=legend.pos, size=size, ...) %>% list() %>% c(lst)
      }

      cowplot::plot_grid(plotlist=lst, ncol=n.col)
  })
}

#' @export
estimteWeightEntropyPerCell <- function(graph, factor.per.cell, annotation=NULL) {
  if (length(unique(factor.per.cell)) != 2)
    stop("factor.per.cell must have exactly two factors")

  adj.mat <- igraph::as_adjacency_matrix(graph, attr="weight") %>% as("dgTMatrix")
  factor.per.cell %<>% as.factor() %>% .[rownames(adj.mat)]
  if (!is.null(annotation)) {
    annotation %<>% as.factor() %>% .[rownames(adj.mat)]
  }

  weight.sum.per.fac.cell <- conos:::getSumWeightMatrix(adj.mat@x, adj.mat@i, adj.mat@j, as.integer(factor.per.cell)) %>%
    `colnames<-`(levels(adj.mat)) %>% `rownames<-`(rownames(adj.mat))

  if (is.null(annotation)) {
    xt <- table(factor.per.cell)
    max.ent <- (if (xt[1] > xt[2]) c(0, 1) else c(1, 0)) %>% entropy::KL.empirical(xt, unit='log2')
    entropy.per.cell <- apply(weight.sum.per.fac.cell, 1, entropy::KL.empirical, xt, unit='log2') / max.ent
  } else {
    xt.per.type <- factor.per.cell %>% split(annotation) %>% sapply(table) %>% t()
    max.ent.per.type <- apply(xt.per.type, 1, function(xt)
      (if (xt[1] > xt[2]) c(0, 1) else c(1, 0)) %>% entropy::KL.empirical(xt, unit='log2'))
    entropy.per.cell <- sapply(1:nrow(weight.sum.per.fac.cell), function(i)
      entropy::KL.empirical(weight.sum.per.fac.cell[i,], xt.per.type[annotation[i],], unit='log2') / max.ent.per.type[annotation[i]])
  }

  return(cbind(data.frame(weight.sum.per.fac.cell), entropy=entropy.per.cell))
}
