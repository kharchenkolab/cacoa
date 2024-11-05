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
