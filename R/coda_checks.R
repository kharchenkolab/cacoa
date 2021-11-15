#' Check cell count data
#' 
#' @param d.counts Cell count table
#' @keywords internal
checkData <- function(d.counts){
  if(is.null(d.counts)) 
    stop('Cell count matrix is not provided')
  if(nrow(d.counts) == 0) 
    stop('Sample size is zero')
  if(ncol(d.counts) == 0) 
    stop('Cell count information is missed')
}


#' Check groups
#' 
#' @param d.groups Groups variable for samples
#' @keywords internal
checkGroups <- function(d.groups){
  if(is.null(d.groups)) 
    stop('Groups are not provided')
  if(length(d.groups) == 0) 
    stop('Group information is missed')
}


#' Check cell count data and groups
#' 
#' @param d.counts Cell count table
#' @param d.groups Groups variable for samples
#' @keywords internal
checkDataGroups <- function(d.counts, d.groups){
  checkData(d.counts)
  checkGroups(d.groups)
  if(nrow(d.counts) != length(d.groups)) 
    stop('Sample size in cell count matrix and in group variable is not the same')
}


#' Check sequential binary partitions(sbp)
#' 
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkSbp <- function(sbp){
  if(is.null(sbp)) 
    stop('Sequential binary partitions are not provided')
  if(nrow(sbp) == 0)
    stop('No partitions are provided')
}


#' Check sequential binary partitions(sbp) for the whole set of balances
#' 
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkSbpWhole <- function(sbp){
  checkSbp(sbp)
  if(nrow(sbp) < 2) 
    stop('Balances require at least 2 cell types')
  if((nrow(sbp) - ncol(sbp)) != 1) 
    stop('Wrong number of balances provided')
}


#' Check the agreement between sbp and cell count table
#' 
#' @param d.counts Cell count table
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkDataSbp <- function(d.counts, sbp){
  checkData(d.counts)
  checkSbp(sbp)
  if(ncol(sbp) != ncol(d.counts))
    stop('Number of cell types in data and sbp should be equal')
}


#' Check data and list of cells
#' 
#' @param d.counts Cell count table
#' @param cells List of cell types of interest
#' @keywords internal
checkDataAndCells <- function(d.counts, cells){
  if(is.null(cells))
    stop('Cell types are not provided')
  if(length(cells) == 0)
    stop('Cell types are not provided')
  if(sum(cells %in% colnames(d.counts)) != length(cells))
    stop('Some cell types do not exist in the data')
}


#' Check Tree
#' 
#' @param tree phylo tree
#' @keywords internal
checkTree <- function(tree){
  if(is.null(tree))
    stop('Tree is not provided')
  # typeof
}


#' Check agreement between Tree and the Cell count table
#' 
#' @param d.counts Cell count table
#' @param tree phylo tree
#' @keywords internal
checkDataTree <- function(d.counts, tree){
  checkTree(tree)
  checkData(d.counts)
  checkDataAndCells(d.counts, tree$tip.label)
}

