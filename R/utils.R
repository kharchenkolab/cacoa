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

#' @keywords internal
summarizeLMResults <- function(res.list, adj.method="BH") {
    ct   <- names(res.list)
    eff  <- sapply(res.list, \(x) unname(x$res$stat.obs[1]))
    pvalues <- sapply(res.list, \(x) x$res$pval)
    out  <- data.frame(celltype=ct, obs.stat=eff, pvalue=pvalues, stringsAsFactors=FALSE)
    padjust <- p.adjust(pvalues, method=adj.method)
    out$p.adj <- padjust
    out <- out[order(out$obs.stat, decreasing=TRUE), ]
    p.dist.info    <- lapply(res.list, function(x) x$dist.mat)
    dists.per.type <- lapply(res.list, function(x) x$dists)
    r2 <- lapply(res.list, function(x) x$r2)
    list(results=out, p.dist.info=p.dist.info, dists.per.type=dists.per.type, r2=r2, pvalues=pvalues, padjust=padjust)
}
