#' @import Matrix
#' @import ggplot2
#' @useDynLib cacoa
NULL


#' @keywords internal
.onUnload <- function (libpath) {
  library.dynam.unload("cacoa", libpath)
}

#' Iterate over tree (list of lists of lists, etc) for `length(ids)` level,
#' interpret each element as a data.frame and bind them appending all levels as columns with
#' colnames corresponding to `ids`
#'
#' @keywords internal
rblapply <- function(list, ids, func) {
  if (length(ids) == 1) return(bind_rows(lapply(list, func), .id=ids[1]))
  return(bind_rows(lapply(list, rblapply, tail(ids, -1), func), .id=ids[1]))
}

# block randomization 
#' @keywords internal
permuteWithinBlocks <- function(v, blocks) {
  if (is.null(blocks)) return(sample(v, length(v), replace = FALSE))
  ave(v, blocks, FUN = function(w) sample(w, length(w), replace = FALSE))
}

#' Extract DE result table from various possible formats
#' @keywords internal
getDeTable <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.list(x) && !is.null(x$res) && is.data.frame(x$res)) return(x$res)
  NULL
}

#' @keywords internal
pairTableToSquare <- function(pair.vec, pairs, sids, samp=NULL) {
    if (is.null(samp)) samp <- sids
    samp <- intersect(samp, sids)
    M <- matrix(NA_real_, length(samp), length(samp), dimnames=list(samp, samp))
    diag(M) <- 0
    for (k in seq_len(nrow(pairs))) {
      a <- sids[pairs[k,1]]
      b <- sids[pairs[k,2]]
      if (!(a %in% samp && b %in% samp)) next
      v <- pair.vec[k]
      if (!is.finite(v)) next
      M[a,b] <- v
      M[b,a] <- v
    }
    keep <- which(rowSums(is.finite(M), na.rm=TRUE) > 1)
    if (length(keep) >= 2) M <- M[keep, keep, drop=FALSE] else M <- M[0,0,drop=FALSE]
    M
  }

#' @keywords internal
attachNC <- function(M, ct, clust.info) {
  if (!is.matrix(M) || nrow(M) == 0) return(M)

  if (!is.null(clust.info$sample.table) && ct %in% rownames(clust.info$sample.table)) {
    nc <- clust.info$sample.table[ct, rownames(M)]
    nc <- as.numeric(nc)
    names(nc) <- rownames(M)
  } else {
    nc <- setNames(rep(0, nrow(M)), rownames(M))
  }

  attr(M, "n.cells") <- nc
  M
}

#' @keywords internal
samplesForCelltype <- function(clust.info, ct) {
  if (!is.null(clust.info$p.dist.info) && !is.null(clust.info$p.dist.info[[ct]])) {
    rownames(clust.info$p.dist.info[[ct]])
  } else {
    NULL
  }
}