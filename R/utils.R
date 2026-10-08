#' @import Matrix
#' @import ggplot2
#' @useDynLib cacoa, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom rlang .data %||%
#' @importFrom stats ave lm.fit terms complete.cases p.adjust.methods
NULL

utils::globalVariables(c("celltype", "effect", "effect_ref_log2", "effect_ref_perm_log2", "go_name", "is_ref_flag", "label", "mlog10p",
                         "prop", "sample.groups", "shape", "signif_flag", "stars", "stat.perm", "term.label"))


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