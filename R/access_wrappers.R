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

extractSampleGroups <- function(obj, ref.level, target.level) UseMethod("extractSampleGroups", obj)
extractSampleGroups.Conos <- function(con, ref.level = ref.level, target.level = target.level) {
  if(!is.null(ref.level)) {
    con.names <- names(con$samples)

    if(!any(grep(ref.level,con.names))) stop("'ref.level' not in Conos sample names.")

    ifelse(grepl(ref.level,con.names),ref.level,target.level) %>% setNames(con.names)
  }
}
