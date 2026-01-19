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