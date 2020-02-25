##' Local Z Scores
##'
##' @param graph alignment graph
##' @param condition.per.cell conditions per cell. Must have two levels.
##' @param genes genes to be tested. Procedure is slow, so it's reasonable to restrict amount of genes.
##' @param count.matrix.transposed joint count matrix with cells by rows and genes by columns
##' @param reference.level reference condition level
localZScores <- function(graph, condition.per.cell, genes, count.matrix.transposed, reference.level, min.expr=1e-10, min.std=1e-10, n.cores=2, verbose=1) {
  if (length(unique(condition.per.cell)) != 2)
    stop("Exactly two condition levels must be provided")

  if (!(reference.level %in% unique(condition.per.cell)))
    stop("reference.level must be presented in condition.per.cell")

  # adj_mat <- igraph::as_adj(graph, attr="weight")
  if (verbose) cat("Estimating adjacency matrices... ")
  adj.mat <- igraph::as_adj(graph)
  adj.mat.by.cond <- colnames(adj.mat) %>% split(condition.per.cell[.]) %>%
    lapply(function(ns) adj.mat[,ns])
  if (verbose) cat("Done.\n")

  if (verbose) cat("Estimating local means... ")
  local.mean.mat <- sccore:::plapply(adj.mat.by.cond, function(mat)
    (mat %*% count.matrix.transposed[colnames(mat), genes]) / pmax(Matrix::rowSums(mat), min.expr),
    n.cores=n.cores, progress=(verbose > 1))
  if (verbose) cat("Done.\n")

  if (verbose) cat("Estimating stds... ")
  ref.adj.mat <- adj.mat.by.cond[[reference.level]]
  ref.cbs <- colnames(ref.adj.mat)
  stds <- sqrt(ref.adj.mat %*% (count.matrix.transposed[ref.cbs, genes] - local.mean.mat[[reference.level]][ref.cbs, genes]) ^ 2 /
                 pmax(Matrix::rowSums(ref.adj.mat), min.expr))
  if (verbose) cat("Done.\n")

  if (verbose) cat("Estimating z-scores... ")
  non.ref.level <- unique(condition.per.cell) %>% setdiff(reference.level)
  z.scores <- (local.mean.mat[[non.ref.level]] - local.mean.mat[[reference.level]]) / pmax(stds, min.std)
  if (verbose) cat("Done.\n")

  return(z.scores)
}

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

plotGeneComparisonBetweenCondition <- function(genes, con, condition.per.cell, z.scores=NULL, cur.scores=NULL, show.legend=T, legend.pos=c(1, 1), size=0.2, n.col=NULL, ...) {
  if (!is.null(cur.scores) & !is.null(names(cur.scores))) {
    genes <- names(sort(cur.scores, decreasing=T)[genes])
  }

  if (is.null(n.col)) n.col = length(unique(condition.per.cell)) + 1

  lapply(genes, function(g) {
    lst <- lapply(unique(condition.per.cell), function(sg) {
      con$plotGraph(gene=g, show.legend=show.legend, legend.pos=legend.pos, size=size,
                      title=paste(sg," ",g), groups=condition.per.cell, subgroups=sg, ...)
    })

      if (!is.null(z.scores)) {
        lst <- cacoa:::plotZScores(g, con, z.scores, cur.scores=cur.scores, show.legend=show.legend,
                                   legend.pos=legend.pos, size=size, ...) %>% list() %>% c(lst)
      }

      cowplot::plot_grid(plotlist=lst, ncol=n.col)
  })
}
