#' @keywords internal
estimateGraphVarianceSignificance <- function(adj.mat, signal) {
  if (!is.numeric(signal)) {
    signal <- as.factor(signal)
    signal[is.na(signal)] <- table(signal) %>% which.max() %>% names()
    comp.op <- "!="
  } else {
    signal[is.na(signal)] <- median(signal, na.rm=TRUE)
    comp.op <- "-"
  }

  obs.var <- signal %>% outer(., ., comp.op) %>% {. * . * adj.mat} %>% sum()
  perm.vars <- sapply(1:n.permutations, function(i) {
    sample(signal) %>% outer(., ., comp.op) %>% {. * . * adj.mat} %>% sum()
  })

  pvalue <- (sum(perm.vars <= obs.var) + 1) / (n.permutations + 1)
  pr2 <- 1 - obs.var / median(perm.vars)
  return(list(pvalue=pvalue, pr2=pr2))
}

#' @keywords internal
adjacencyMatrixFromPaiwiseDists <- function(p.dists, trim=0.05, k=NULL) {
  adj.mat <- p.dists %>% pmin(quantile(., 1 - trim)) %>%
    {pmax(0, . - quantile(., trim))} %>% {1 - . / max(.)} %>%
    matrix(ncol=ncol(p.dists))
  diag(adj.mat) <- 0

  if (!is.null(k) && (k < ncol(adj.mat))) { # remove edges for all but k nearest neighbors
    adj.mat %<>% apply(1, function(r) ifelse(r < sort(r, decreasing=TRUE)[k], 0, r)) %>% {(. + t(.)) / 2}
  }

  dimnames(adj.mat) <- dimnames(p.dists)

  return(adj.mat)
}

#' @keywords internal
estimateUMAPOnDistances <- function(p.dists, n.neighbors=15, verbose=FALSE, ...) {
  set.seed(42)
  idx <- (1:ncol(p.dists)) %>% lapply(function(i) head(order(p.dists[i,]), n.neighbors))
  dists <- (1:ncol(p.dists)) %>% lapply(function(i) p.dists[i, idx[[i]]])

  idx <- do.call(rbind, idx)
  dists <- do.call(rbind, dists)

  umap <- uwot::umap(data.frame(x=rep(0, nrow(dists))), nn_method=list(idx=idx, dist=dists), verbose=verbose, ...)

  return(umap)
}
