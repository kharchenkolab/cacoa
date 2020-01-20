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


estimateExpressionShiftZScores <- function(obj, ...) UseMethod("expressionExpressionShiftZScores", obj)

##' @param con conos object
##' @inheritDotParams expressionExpressionShiftZScores.default
##' @rdname expressionExpressionShiftZScores
estimateExpressionShiftZScores.Conos <- function(con, n.od.genes=1000, n.pcs=100, pca.maxit=1000, ...) {
  cm.merged <- con$getJointCountMatrix(raw=F)
  od.genes <- conos:::getOdGenesUniformly(con$samples, n.od.genes)

  cm.merged[, od.genes] %>%
    irlba::irlba(nv=n.pcs, nu=0, center=Matrix::colMeans(.), right_only=F, fastpath=T, maxit=pca.maxit, reorth=T) %>%
    estimateExpressionShiftZScores.default(...)
}
