#' @import Matrix
#' @import ggplot2
NULL

#' @useDynLib cacoa
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("cacoa", libpath)
}
