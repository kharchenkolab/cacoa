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
#' @param cm            dgCMatrix (genes x cells)
#' @param sample.per.cell vector/factor of sample IDs for each cell (length = ncol(cm))
#' @param nns.per.cell  list of integer neighbor indices per cell (same length as ncol(cm))
#' @param x             output of buildDesignPairMatrices (must include $X, $Z, $contrast, $pairs)
#' @param dist,log.vecs, min.n.obs.per.samp as before
#' @param perm.method   "freedman-lane" or "block"
#' @param robust.method "none","huber","winsor"
#' @param na.mode       "drop" or "impute_weak"
#' @param alternative   "two-sided","greater","less"
#' @param smooth        logical; if TRUE, median-smooth effect sizes (and optionally z)
#' @param wins          winsorization fraction for permutation adjustment (e.g., 0.01)
#' @param return.sampled.stats  set TRUE if you want permutation-adjusted z's
#' @param return.residuals      pass-through to cpp
#' @param verbose              logical
#' @keywords internal
estimateClusterFreeExpressionShiftsLM <- function(cm, sample.per.cell, nns.per.cell, x, dist = "cor", log.vecs = TRUE,
                                                  min.n.obs.per.samp = 2L, n.cores = 1, perm.method = c("freedman-lane","block"), 
                                                  robust.method = c("none","huber","winsor"),n.permutations = 999,  adjust= TRUE,
                                                  na.mode = c("drop","impute_weak"), alternative = c("two-sided","greater","less"),
                                                  smooth = TRUE, wins = 0.025, return.sampled.stats = FALSE, 
                                                  return.residuals = FALSE, verbose = FALSE) {
  perm.method   <- match.arg(perm.method)
  robust.method <- match.arg(robust.method)
  na.mode       <- match.arg(na.mode)
  alternative   <- match.arg(alternative)

  #if (adjust) return.sampled.stats <- TRUE
  
  spc <- sample.per.cell
  if (is.factor(spc)) spc <- as.integer(spc) else spc <- as.integer(as.factor(spc))
  spc <- spc[match(colnames(cm), names(sample.per.cell))]
  stopifnot(length(spc) == ncol(cm))

  # neighbors as integers
  nn.list <- lapply(nns.per.cell, function(v) as.integer(v))

  # Pairs must correspond to the order of rows in x$X:
  stopifnot(!is.null(x$pairs))
  pairs.mat <- as.matrix(x$pairs)
  stopifnot(ncol(pairs.mat) == 2L)
  storage.mode(pairs.mat) <- "integer"

  ## ---- response matrix (pairwise distances per neighborhood) ----
  Y <- estimateExpressionShiftsPairsLM(cm = cm, sample_per_cell = spc, nn_ids = nn.list, pairs_mat = pairs.mat,
                                       min_n_obs_per_samp= min.n.obs.per.samp, dist = dist, log_vecs = log.vecs)
  if (!is.null(names(nn.list))) colnames(Y) <- names(nn.list)
  #na.cols <- which(colSums(is.na(Y)) > 0)

  ## ---- fit & permutations ----
  res <- performLMPermutations(x = x, y = Y, n.permutations = n.permutations, perm.method = perm.method,
                               robust.method = robust.method, na.mode = na.mode, alternative = alternative,
                               return.sampled.stats = return.sampled.stats, return.residuals = return.residuals,
                               n.cores = n.cores)

  valid <- if (!is.null(res$perm.valid)) {
    is.finite(res$stat.obs) & res$perm.valid
  } else {
    is.finite(res$stat.obs)
  }
  non.zero.ids <- which(valid)
  non.zero.ids.c <- as.integer(non.zero.ids - 1L)

  ## new cpp function; uses old functions internally; 
  z.adj <- if (adjust) adjustedZScoresMaxStat(z_obs = as.numeric(res$stat.obs), max_vals_in  = as.numeric(res$max.perm),   
                                  min_vals_in  = as.numeric(res$min.perm), alt = if (alternative == "two-sided") 0L 
                                  else if (alternative == "greater") 1L else 2L, wins = wins, smooth = smooth, 
                                  nn_ids = nn.list, non_zero_ids = non.zero.ids.c) else NULL

  ## ---- effect-size shifts (center by permutation mean) ----
  shifts <- res$stat.obs
  if (!is.null(res$stats.perm)) {
    shifts <- res$stat.obs - colMeans(res$stats.perm, na.rm = TRUE)
  }
  ## optional effect size smoothing
  shifts.smoothed <- if (smooth) applyMedianFilterES(res$stat.obs, nn_ids = nn.list,
                                           non_zero_ids = which(is.finite(res$stat.obs))) else NULL
  
  if (!is.null(colnames(Y))) names(shifts) <- names(res$z.score) <- colnames(Y)
  if (!is.null(z.adj))  names(z.adj) <- colnames(Y)
  if (!is.null(shifts.smoothed)) names(shifts.smoothed) <- colnames(Y)


  list(
    stat            = res$stat.obs,
    p.value         = res$pval,
    z.scores        = res$z.score,
    z.adj           = z.adj,
    shifts          = shifts,
    shifts.smoothed = shifts.smoothed,
    sampled.stats   = if (!is.null(res$sampled_stats)) res$sampled_stats else NULL,
    residuals       = if (!is.null(res$residuals)) res$residuals else NULL
  )
}