# accessor methods for uniform access to data from Conos, Seurat and other objects


#' Extract the cell groups from the object
#'
#' @param object object from which to extract the cell groups
#' @rdname extractCellGroups
#' @export
extractCellGroups <- function(object) UseMethod("extractCellGroups", object)

#' @rdname extractCellGroups
extractCellGroups.Conos <- function(object) {
  if(is.null(object$clusters)) {
    stop('No cell groups specified and no clusterings found')
  }
  return(as.factor(object$clusters[[1]]$groups))
}

#' @rdname extractCellGroups
extractCellGroups.Seurat <- function(object) {
  groups <- Seurat::Idents(object)
  if (length(unique(groups)) <= 1) {
    stop('No cell groups specified and no clusterings found')
  }
  return(as.factor(groups))
}




#' Extract the cell groups from the object
#'
#' @param object object from which to extract the cell groups
#' @rdname extractCellGraph
#' @export
extractCellGraph <- function(object) UseMethod("extractCellGraph", object)

#' @rdname extractCellGraph
extractCellGraph.Conos <- function(object) {
  if (is.null(object$graph)) {
    stop('No cell graph found in the object')
  }
  return(object$graph)
}

#' @rdname extractCellGraph
extractCellGraph.Seurat <- function(object) {
  graph.name <- object@misc$graph.name
  if (!is.null(graph.name)) {
    if (is.null(object@graphs[[graph.name]])){
      stop("Cannot find graph ", graph.name)
    }
    graph <- object@graphs[[graph.name]]
  } else {
    if (length(object@graphs) == 0) {
      stop('No cell graph found in the object')
    }
    graph <- object@graphs[[1]]
  }

  graph %<>% as("dgCMatrix") %>%
    igraph::graph_from_adjacency_matrix(mode="undirected", weighted=TRUE)
  return(graph)
}




#' Extract the cell groups from the object
#'
#' @param object object from which to extract the cell groups
#' @param transposed boolean Whether the raw count matrix should be transposed (default=TRUE)
#' @rdname extractRawCountMatrices
#' @export
extractRawCountMatrices <- function(object, transposed=TRUE) UseMethod("extractRawCountMatrices", object)

#' @rdname extractRawCountMatrices
extractRawCountMatrices.Conos <- function(object, transposed=TRUE) {
  return(lapply(object$samples, conos::getRawCountMatrix, transposed=transposed))
}

#' @rdname extractRawCountMatrices
extractRawCountMatrices.Seurat <- function(object, transposed=TRUE) {
  cms <- object$sample.per.cell %>% {split(names(.), .)} %>%
    lapply(function(cids) object@assays$RNA@counts[,cids])
  if (transposed) {
    cms %<>% lapply(Matrix::t)
  }
  return(cms)
}

#' @rdname extractRawCountMatrices
extractRawCountMatrices.dgCMatrix <- function(object, transposed=TRUE) {
  sample.per.cell <- attr(object, 'sample.per.cell')
  cms <- sample.per.cell %>% {split(names(.), .)} %>%
    lapply(function(cids) object[cids,])
  if (!transposed) {
    cms %<>% lapply(Matrix::t)
  }
  return(cms)
}




#' Extract the joint count matrix from the object
#'
#' @param object object from which to extract the cell groups
#' @param raw boolean If TRUE, return merged "raw" count matrices (default=TRUE)
#' @param ... additional parameters to be passed to extractJointCountMatrix
#' @rdname extractJointCountMatrix
#' @export
extractJointCountMatrix <- function(object, raw=TRUE, ...) UseMethod("extractJointCountMatrix", object)

#' @rdname extractJointCountMatrix
extractJointCountMatrix.Conos <- function(object, raw=TRUE) {
  return(object$getJointCountMatrix(raw=raw))
}

#' @param transposed boolean If TRUE, return merged transposed count matrices (default=TRUE)
#' @param raw boolean If TRUE, return merged the sparse dgCMatrix matrix (default=TRUE)
#' @rdname extractJointCountMatrix
extractJointCountMatrix.Seurat <- function(object, raw=TRUE, transposed=TRUE, sparse=TRUE) {
  if (raw) {
    dat <- object@assays$RNA@counts
    if (transposed){
      dat %<>% Matrix::t()
    }
    return(dat)
  }

  dat <- Seurat::GetAssayData(object, slot='scale.data', assay='RNA')
  dims <- dim(dat)
  dat.na <- all(dims == 1) && all(is.na(x = dat))
  if (all(dims == 0) || dat.na) {
    dat <- Seurat::GetAssayData(object, slot='data', assay='RNA')
  }

  if (transposed){
    dat %<>% Matrix::t()
  }
  if (is.matrix(dat) && sparse){
    dat %<>% as("dgCMatrix")
  }
  return(dat)
}

#' @rdname extractJointCountMatrix
extractJointCountMatrix.dgCMatrix <- function(object, raw=TRUE) {
  if (attr(object, 'raw') == raw){
    return(object)
  }

  if (!attr(object, 'raw') && raw){
    stop("Cannot extract raw matrix from a normalized dgCMatrix")
  }

  return(t(t(object) / pmax(colSums(object), 0.1)))
}



#' Extract the top overdispersed genes from the object
#'
#' @param object object from which to extract the top overdispersed genes
#' @param n.genes numeric Number of overdispersed genes to extract (default=NULL)
#' @rdname extractOdGenes
#' @export
extractOdGenes <- function(object, n.genes=NULL) UseMethod("extractOdGenes", object)

#' @rdname extractOdGenes
extractOdGenes.Conos <- function(object, n.genes=NULL) {
  getOdGenesUniformly <- utils::getFromNamespace("getOdGenesUniformly", "conos")
  return(getOdGenesUniformly(object$samples, n.genes))
}

#' @rdname extractOdGenes
extractOdGenes.Seurat <- function(object, n.genes=NULL) {
  if (is.null(object@assays$RNA@meta.features$vst.variance.standardized)){
    stop("The data object doesn't have gene variance info.",
         "Please, run FindVariableFeatures(assay='RNA, selection.method='vst') first")
  }
  genes <- object@assays$RNA@meta.features %>%
    {rownames(.)[order(.$vst.variance.standardized, decreasing=TRUE)]} %>%
    head(n.genes %||% length(.))

  return(genes)
}




#' Extract the sample/dataset per cell from the object
#'
#' @param object object from which to extract the sample/dataset per cell
#' @rdname extractSamplePerCell
#' @export
extractSamplePerCell <- function(object) UseMethod("extractSamplePerCell", object)

#' @rdname extractSamplePerCell
extractSamplePerCell.Conos <- function(object) {
  return(object$getDatasetPerCell())
}

#' @rdname extractSamplePerCell
extractSamplePerCell.Seurat <- function(object) {
  return(object$sample.per.cell)
}



extractSampleGroups <- function(object, ref.level, target.level) {
  samp.names <- unique(extractSamplePerCell(object))
  if(!any(grep(ref.level, con.names))) {
    stop("'ref.level' not in the data object sample names.")
  }
  sg <- grepl(ref.level, samp.names) %>% ifelse(ref.level, target.level) %>%
    setNames(samp.names)
  return(sg)
}




#' Extract embeddings from the object
#'
#' @param object object from which to extract the embeddings
#' @rdname extractEmbedding
#' @export
extractEmbedding <- function(object) UseMethod("extractEmbedding", object)

#' @rdname extractEmbedding
extractEmbedding.Conos <- function(object) {
  return(object$embedding)
}

#' @rdname extractEmbedding
extractEmbedding.Seurat <- function(object) {
  if (!is.null(object@reductions$umap)){
    return(object@reductions$umap@cell.embeddings)
  }
  if (length(object@reductions) == 0){
    return(NULL)
  }
  return(object@reductions[[1]]@cell.embeddings)
}




#' Extract the gene exrpession from the object
#'
#' @param object object from which to extract the cell groups
#' @param gene character vector of the specific gene names on which to subset
#' @rdname extractGeneExpression
#' @export
extractGeneExpression <- function(object, gene) UseMethod("extractGeneExpression", object)

#' @rdname extractGeneExpression
extractGeneExpression.Conos <- function(object, gene) {
  return(conos::getGeneExpression(object, gene))
}

#' @rdname extractGeneExpression
extractGeneExpression.Seurat <- function(object, gene) {
  return(extractJointCountMatrix(object, raw=FALSE, transposed=FALSE, sparse=FALSE)[gene,])
}

#' @rdname extractGeneExpression
extractGeneExpression.dgCMatrix <- function(object, gene) {
  return(object[,gene])
}
