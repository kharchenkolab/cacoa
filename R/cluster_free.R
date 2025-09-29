#' @importFrom sccore checkPackageInstalled
NULL

#' Estimate Gene Programs with FABIA
#' 
#' @param n.programs maximal number of gene programs to find (parameter `p` for fabia).
#' @param n.sampled.cells number of sub-sampled cells for estimating the gene programs. If 0, all cells are used.
#' it is interpreted as a vector of cell
#' @inheritDotParams fabia::fabia -p -X -cyc -alpha -random
#' @keywords internal
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

#' @keywords internal
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

#' @keywords internal
estimateGeneClustersPam <- function(z.scores, n.programs) {
  checkPackageInstalled(c("cluster"), cran=TRUE)
  p.dists <- 1 - cor(as.matrix(z.scores))
  p.dists[is.na(p.dists)] <- 1
  pam.res <- cluster::pam(p.dists, k=n.programs, diss=TRUE)
  return(pam.res$clustering)
}

#' @keywords internal
geneProgramSimilarityScores <- function(program.scores, gene.scores) {
  apply(gene.scores, 2, estimateCorrelationDistance, program.scores, centered=FALSE) %>%
      {1 - .} %>% sort(decreasing=TRUE)
}

#' @keywords internal
geneProgramLoadingScores <- function(program.scores, gene.scores, min.score=0.05) {
  cell.subs <- which(abs(program.scores) > min.score) %>% names()
  scores <- gene.scores[cell.subs,, drop=FALSE] %>% abs() %>% colSums() %>% sort(decreasing=TRUE)
  return(scores)
}

#' @keywords internal
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


#' Estimate cluster-free shifts using linear model and permutations
#' @param cm gene-by-Cell count matrix (dgCMatrix or Matrix), genes x cells
#' @param sample.per.cell Vector of sample IDs for each cell
#' @param nns.per.cell List of nearest neighbors for each cell
#' @param sample.meta Data frame of sample-level metadata
#' @keywords internal
estimateClusterFreeShiftsLM <- function(cm, sample.per.cell, nns.per.cell,
                 sample.meta, contrast, dist = "cor", log.vecs = TRUE,
                 min.n.obs.per.samp = 2L, dist.type = "shift",
                 n.cores = 1,
                 x, n.permutations = 999,
                 perm.method = c("freedman-lane", "block"),
                 robust.method = c("none","huber", "winsor"),
                 na.mode = c("drop", "impute_weak"),
                 return.sampled.stats = FALSE,
                 return.residuals = FALSE,
                 adj.method = "BH",
                 verbose = FALSE) {
  perm.method <- match.arg(perm.method)
  robust.method <- match.arg(robust.method)
  na.mode <- match.arg(na.mode)

  Y.res <- getResponsePairMatrix(cm, sample.per.cell, nns.per.cell,
                                dist = dist, log.vecs = log.vecs,
                                min.n.obs.per.samp = min.n.obs.per.samp,
                                n.cores = n.cores)
  Y <- Y.res$Y
  samples <- Y.res$samples
  pairs.ij <- Y.res$pairs.ij

  stopifnot(all(rownames(sample.meta) %in% samples))
  sample.meta <- sample.meta[samples, , drop = FALSE]
  x <- buildPairDesignMatrices(sample.meta, triplet=contrast, dist.type = "shift")$model
  
  res <- performLMPermutations(
    x, Y, n.permutations = n.permutations,
    perm.method = perm.method,
    robust.method = robust.method,
    na.mode = na.mode,
    return.sampled.stats = return.sampled.stats,
    return.residuals = return.residuals,
    verbose = verbose
  )
  names(res$z.score) <- names(res$stat.obs) <- colnames(Y)
  res
}

#' Function for computing response pair matrix Y for pairwise sample distances across a set of neighborhoods
#' @param cm gene-by-Cell count matrix (dgCMatrix or Matrix), t(count.matrices)
#' @param sample.per.cell Vector of sample IDs for each cell
#' @param nns.per.cell List of nearest neighbors for each cell
#' @param dist Distance metric to use (default: "cor")
#' @param log.vecs Whether to log-transform expression vectors (default: TRUE)
#' @param min.n.obs.per.samp Minimum number of observations per sample (default: 1)
#' @examples 
#' \dontrun{
#'  b <- getResponsePairMatrix(cm, sample_per_cell, nns.per.cell,
#'  dist = "cor", log.vecs = TRUE, min.n.obs.per.samp = 3L, n.cores = 5)
#' }
#' @keywords internal
getResponsePairMatrix <- function(cm, sample.per.cell, nns.per.cell,
                 dist = "cor", log.vecs = TRUE,
                 min.n.obs.per.samp = 1L,
                 n.cores = 1) {
  stopifnot(inherits(cm, "dgCMatrix") || inherits(cm, "Matrix"))

  # samples & sample ids (0..n-1) once
  samples <- levels(factor(sample.per.cell))
  n <- length(samples)
  spc   <- factor(sample.per.cell, levels = samples)
  samp.0 <- as.integer(spc) - 1L
  
  # precompute lower-tri bookkeeping
  pairs.ij <- lowerTriIJ(n)
  n.pairs  <- nrow(pairs.ij)
  # prefix term T[i] = (i-1)*(i-2)/2 used in lt index; length n
  T.pref <- ((seq_len(n) - 1L) * (seq_len(n) - 2L)) %/% 2L
  
  # global -> local col map (use match; faster than named lookup)
  global.ids <- getGlobalCellIds(colnames(cm))
  if (anyNA(global.ids)) stop("Couldn't parse global IDs from cm colnames.")
  # to 1-based local
  map.to.local.1 <- function(ids.global.0) match(ids.global.0, global.ids)
  # to 0-based local
  map.to.local.0 <- function(ids.global.0) {
  pos.1 <- map.to.local.1(ids.global.0)
  pos.1 <- pos.1[!is.na(pos.1)]
  as.integer(pos.1 - 1L)
  }
  
  # worker per neighborhood
  per.nb <- function(ids.global.0) {
  ids.local.0 <- map.to.local.0(as.integer(ids.global.0))
  if (!length(ids.local.0)) return(rep(NA_real_, n.pairs))
  
  res <- estimateCellExpressionShift_export( 
    cm,
    as.integer(samp.0),                 # 0..n-1 samples
    as.integer(ids.local.0),            # 0-based local cells
    min.n.obs.per.samp = as.integer(min.n.obs.per.samp),
    dist = dist, log.vecs = log.vecs
  )
  
  # normalize possible field names from Rcpp
  d  <- res[["dists"]]; if (is.null(d))  d  <- res[["dist"]]
  s.1 <- res[["s1.ids"]]; if (is.null(s.1)) s.1 <- res[["s1"]]
  s.2 <- res[["s2.ids"]]; if (is.null(s.2)) s.2 <- res[["s2"]]
  if (is.null(d) || is.null(s.1) || is.null(s.2)) return(rep(NA_real_, n.pairs))
  
  keep <- is.finite(d)
  if (!any(keep)) return(rep(NA_real_, n.pairs))
  
  # back to 1-based sample indices and enforce lower-tri (i>j)
  i <- as.integer(s.1[keep]) + 1L
  j <- as.integer(s.2[keep]) + 1L
  i.L <- pmax(i, j); j.L <- pmin(i, j)
  
  # compute positions vectorized: pos = T.pref[i.L] + j.L
  pos <- T.pref[i.L] + j.L
  
  y.k <- rep(NA_real_, n.pairs)
  # write all distances at once (duplicates are rare; last wins)
  y.k[pos] <- as.numeric(d[keep])
  y.k
  }
  
  Y.cols <- sccore::plapply(nns.per.cell, per.nb, n.cores = n.cores, progress = TRUE, mc.preschedule = FALSE, 
                            mc.allow.recursive = TRUE, fail.on.error = FALSE)
  
  Y <- do.call(cbind, Y.cols)
  nbhd.names <- names(nns.per.cell)
  if (is.null(nbhd.names) || any(!nzchar(nbhd.names))) {
  nbhd.names <- paste0("cell.", seq_along(nns.per.cell))
  }
  colnames(Y) <- nbhd.names
  rownames(Y) <- paste0("(", pairs.ij$i, ",", pairs.ij$j, ")")
  list(Y = Y, samples = samples, pairs.ij = pairs.ij)
}

# helpers
# trailing digits in "cell.123" → 123
getGlobalCellIds <- function(cn) as.integer(sub("^.*?(\\d+)$", "\\1", cn))
  
# fixed lower-tri row order (i>j). We keep it for rownames only.
lowerTriIJ <- function(n) {
  ij <- which(lower.tri(matrix(NA_real_, n, n)), arr.ind = TRUE)
  data.frame(i = ij[,1], j = ij[,2])
  }
  
# vectorized lower-tri index: for i>j, pos = (i-1)*(i-2)/2 + j  (1-based)
.lt.index <- function(i, j) ((i - 1L) * (i - 2L)) %/% 2L + j
