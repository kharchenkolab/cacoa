#' Estimate Gene Programs with FABIA
#' @param n.programs maximal number of gene programs to find (parameter `p` for fabia).
#' @param n.sampled.cells number of sub-sampled cells for estimating the gene programs. If 0, all cells are used.
#' it is interpreted as a vector of cell
#' @inheritDotParams fabia::fabia -p -X -cyc -alpha -random
estimateGeneProgramsFabia <- function(z.scores, n.programs, n.sampled.cells=15000, cyc=1500, alpha=0.2, random=-1, ...) {
  checkPackageInstalled("fabia", bioc=TRUE)

  sample.ids <- NULL
  if ((n.sampled.cells > 0) && (n.sampled.cells < nrow(z.scores))) {
    sample.ids <- unique(round(seq(1, nrow(z.scores), length.out=n.sampled.cells)))
    z.scores <- z.scores[sample.ids,]
  }

  fabia.res <- z.scores %>% t() %>%
    fabia::fabia(p=n.programs, alpha=alpha, cyc=cyc, random=random, ...)

  return(list(fabia=fabia.res, sample.ids=sample.ids))
}

estimateGeneClustersLeiden <- function(z.scores, resolution=1, n.pcs=100, k=30, n.cores=1, verbose=FALSE) {
  checkPackageInstalled(c("N2R", "leidenAlg"), cran=TRUE)
  z.scores %<>% t()
  pcs <- irlba::irlba(z.scores, nv=n.pcs, nu=0, right_only=FALSE, fastpath=TRUE, maxit=1000, reorth=TRUE, verbose=verbose)
  pcas <- as.matrix(z.scores %*% pcs$v)

  xn <- N2R::Knn(as.matrix(pcas), k, nThreads=n.cores, verbose=verbose, indexType='angular')
  xn@x <- 1 - xn@x
  xn <- (xn + t(xn)) / 2

  g <- igraph::graph_from_adjacency_matrix(xn, mode='undirected', weighted=TRUE)
  clusts <- leidenAlg::leiden.community(g, resolution=resolution, n.iterations=10) %$%
    setNames(membership, rownames(pcas))
  return(clusts)
}

estimateGeneClustersPam <- function(z.scores, n.programs) {
  checkPackageInstalled(c("cluster"), cran=TRUE)
  p.dists <- 1 - cor(as.matrix(z.scores))
  p.dists[is.na(p.dists)] <- 1
  pam.res <- cluster::pam(p.dists, k=n.programs, diss=TRUE)
  return(pam.res$clustering)
}

geneProgramSimilarityScores <- function(program.scores, gene.scores) {
  apply(gene.scores, 2, estimateCorrelationDistance, program.scores, centered=FALSE) %>%
      {1 - .} %>% sort(decreasing=TRUE)
}

geneProgramLoadingScores <- function(program.scores, gene.scores, min.score=0.05) {
  cell.subs <- which(abs(program.scores) > min.score) %>% names()
  scores <- gene.scores[cell.subs,, drop=FALSE] %>% abs() %>% colSums() %>% sort(decreasing=TRUE)
  return(scores)
}

geneProgramInfoByCluster <- function(clusters, z.scores, min.score=0.05, verbose=FALSE) {
  genes.per.clust <- clusters %>% {split(names(.), .)}
  program.scores <- genes.per.clust %>%
    plapply(function(ns) apply(as.matrix(z.scores[,ns,drop=FALSE]), 1, mean, trim=0.1), progress=verbose) %>%
    do.call(rbind, .) %>% set_colnames(rownames(z.scores))

  sim.scores <- lapply(1:nrow(program.scores), function(pid) {
    geneProgramSimilarityScores(program.scores[pid,], gene.scores=z.scores[, genes.per.clust[[pid]], drop=FALSE])
  })

  loading.scores <- lapply(1:nrow(program.scores), function(pid) {
    geneProgramLoadingScores(program.scores[pid,], z.scores[, genes.per.clust[[pid]], drop=FALSE], min.score=min.score)
  })

  return(list(program.scores=program.scores, genes.per.clust=genes.per.clust, clusters=clusters,
              sim.scores=sim.scores, loading.scores=loading.scores, n.progs=nrow(program.scores)))
}
