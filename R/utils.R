#' @import Matrix
#' @import ggplot2
#' @importFrom magrittr %<>% %>% %$%
NULL

getColumnwiseCorrelations <- function(m1, m2) {
  common.cols <- intersect(colnames(m1), colnames(m2))
  return(lineup::corbetw2mat(m1[, common.cols], m2[, common.cols]))
}


