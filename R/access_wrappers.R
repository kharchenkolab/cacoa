# accessor methods for uniform access to data from Conos, Seurat and other objects

extractGroups <- function(obj) UseMethod("extractGroups", obj)
extractGroups.Conos <- function(con) {
  if(is.null(con$clusters)) stop('no groups specified and no clusterings found')
  return(as.factor(con$clusters[[1]]$groups))
}

extractRawCountMatrices <- function(obj, transposed=T) UseMethod("extractRawCountMatrices", obj)
extractRawCountMatrices.Conos <- function(con, transposed=T) {
  return(lapply(con$samples, conos:::getRawCountMatrix, transposed=transposed))
}

estimateExpressionShiftZScores <- function(obj, ...) UseMethod("estimateExpressionShiftZScores", obj)

##' @param con conos object
##' @inheritDotParams estimateExpressionShiftZScores.default
##' @rdname estimateExpressionShiftZScores
estimateExpressionShiftZScores.Conos <- function(con, n.od.genes=1000, n.pcs=100, pca.maxit=1000, ...) {
  cm.merged <- con$getJointCountMatrix(raw=F)
  od.genes <- conos:::getOdGenesUniformly(con$samples, n.od.genes)

  cm.merged[, od.genes] %>%
    irlba::irlba(nv=n.pcs, nu=0, center=Matrix::colMeans(.), right_only=F, fastpath=T, maxit=pca.maxit, reorth=T) %>%
    estimateExpressionShiftZScores.default(...)
}
