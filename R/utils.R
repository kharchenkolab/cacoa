#' @import Matrix
#' @import ggplot2
NULL

#' @useDynLib cacoa
NULL

getColumnwiseCorrelations <- function(m1, m2) {
  .Deprecated() # TODO: remove lineup dependency
  common.cols <- intersect(colnames(m1), colnames(m2))
  return(lineup::corbetw2mat(m1[, common.cols], m2[, common.cols]))
}

.onUnload <- function (libpath) {
  library.dynam.unload("cacoa", libpath)
}
