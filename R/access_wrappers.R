# accessor methods for uniform access to data from Conos, Seurat and other objects

extractCellGroups <- function(obj) UseMethod("extractCellGroups", obj)
extractCellGroups.Conos <- function(con) {
  if(is.null(con$clusters)) stop('No cell groups specified and no clusterings found')
  return(as.factor(con$clusters[[1]]$groups))
}

extractRawCountMatrices <- function(obj, transposed=T) UseMethod("extractRawCountMatrices", obj)
extractRawCountMatrices.Conos <- function(con, transposed=T) {
  return(lapply(con$samples, conos:::getRawCountMatrix, transposed=transposed))
}

extractJointCountMatrix <- function(obj, raw=T) UseMethod("extractJointCountMatrix", obj)
extractJointCountMatrix.Conos <- function(con, raw=T) {
  return(con$getJointCountMatrix(raw=raw))
}

extractOdGenes <- function(obj, n.genes=NULL) UseMethod("extractOdGenes", obj)
extractOdGenes.Conos <- function(con, n.genes=NULL) {
  return(conos:::getOdGenesUniformly(con$samples, n.genes))
}

extractSamplePerCell <- function(obj) UseMethod("extractSamplePerCell", obj)
extractSamplePerCell.Conos <- function(con) {
  return(con$getDatasetPerCell())
}

extractSampleGroups <- function(obj) UseMethod("extractSampleGroups", obj, ref.level, target.level)
extractSampleGroups.Conos <- function(con, ref.level=ref.level, target.level=target.level) {
  con.names <- names(con$samples)
  if(!any(grep(ref.level,con.names))) stop("'ref.level' not in Conos sample names.")
  if(!any(grep(target.level,con.names))) stop("'target.level' not in Conos sample names.")
  return(ifelse(grepl(ref.level,names(con$samples)),ref.level,target.level))
}
