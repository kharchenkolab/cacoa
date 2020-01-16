# accessor methods for uniform access to data from Conos, Seurat and other objects

estimateExpressionShiftMagnitudes <- function(obj, ...) UseMethod("estimateExpressionShiftMagnitudes", obj)

##' @param con conos object
##' @inheritDotParams estimateExpressionShiftMagnitudes.default
##' @rdname estimateExpressionShiftMagnitudes
estimateExpressionShiftMagnitudes.Conos <- function(con, groups=NULL, ...) {
  if(is.null(groups)) {
    if(is.null(con$clusters)) stop('no groups specified and no clusterings found')
    groups <- as.factor(con$clusters[[1]]$groups)
  } else {
    groups <- as.factor(groups)
  }

  lapply(con$samples, conos:::getRawCountMatrix, transposed=T) %>%
    estimateExpressionShiftMagnitudes(..., groups=groups, transposed.matrices=T)
}
