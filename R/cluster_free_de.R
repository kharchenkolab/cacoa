##' Local Z Scores
##'
##' @param graph Alignment graph (embedding)
##' @param condition.per.cell Named group factor with cell names. Must have exactly two levels.
##' @param genes Genes to be tested. Procedure is slow, so it's reasonable to restrict amount of genes.
##' @param count.matrix.transposed Joint count matrix with cells by rows and genes by columns
##' @param ref.level Reference condition level, e.g., wt, ctrl, or healthy
localZScores <- function(graph, condition.per.cell, genes, count.matrix.transposed, ref.level, min.expr=1e-10, min.std=1e-10, n.cores=1, verbose=T) {
  #TODO: Condition.per.cell should be derived from sample.groups and cell.groups
  if (length(unique(condition.per.cell)) != 2)
    stop("Exactly two levels must be provided in 'condition.per.cell'")

  if (!(ref.level %in% unique(condition.per.cell)))
    stop("'ref.level' not present in 'condition.per.cell'")

  # adj_mat <- igraph::as_adj(graph, attr="weight")
  if (verbose) cat("Estimating adjacency matrices... ")
  adj.mat <- igraph::as_adj(graph)
  adj.mat.by.cond <- colnames(adj.mat) %>% split(condition.per.cell[.]) %>%
    lapply(function(ns) adj.mat[,ns])
  if (verbose) cat("Done.\n")

  if (verbose) cat("Estimating local means... ")
  local.mean.mat <- sccore:::plapply(adj.mat.by.cond, function(mat)
    (mat %*% count.matrix.transposed[colnames(mat), genes]) / pmax(Matrix::rowSums(mat), min.expr),
    n.cores=n.cores, progress=verbose)
  if (verbose) cat("Done.\n")

  if (verbose) cat("Estimating stds... ")
  ref.adj.mat <- adj.mat.by.cond[[ref.level]]
  ref.cbs <- colnames(ref.adj.mat)
  stds <- sqrt(ref.adj.mat %*% (count.matrix.transposed[ref.cbs, genes] - local.mean.mat[[ref.level]][ref.cbs, genes]) ^ 2 /
                 pmax(Matrix::rowSums(ref.adj.mat), min.expr))
  if (verbose) cat("Done.\n")

  if (verbose) cat("Estimating z-scores... ")
  non.ref.level <- unique(condition.per.cell) %>% setdiff(ref.level)
  z.scores <- (local.mean.mat[[non.ref.level]] - local.mean.mat[[ref.level]]) / pmax(stds, min.std)
  if (verbose) cat("Done.\n")

  return(z.scores)
}

##' Filter Local Z Scores
##'
##' @param z.scores Output from localZScores function
##' @param min.row.z (default=1)
##' @param min.inidividual z (default=0.25)
##' @param n.cores Cores for parallel processing (default=1)
filterLocalZScores <- function(z.scores, min.row.z=1, min.individual.z=0.25, n.cores=1) {
  app.func <- if (requireNamespace("pbapply", quietly=T)) function(...) pbapply::pbapply(..., cl=n.cores) else apply
  z.scores.abs <- z.scores
  z.scores.abs@x %<>% abs()
  max.z <- app.func(z.scores.abs, 2, quantile, 0.99)
  z.scores <- z.scores[, max.z > min.row.z]
  z.scores@x[abs(z.scores@x) < min.individual.z] <- 0
  return(drop0(z.scores))
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
  cur.score <- if (is.null(cur.scores)) sum(z.scores[, gene]) else cur.scores[gene]
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
