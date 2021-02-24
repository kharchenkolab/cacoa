# accessor methods for uniform access to data from Conos, Seurat and other objects

extractCellGroups <- function(obj) UseMethod("extractCellGroups", obj)
extractCellGroups.Conos <- function(con) {
  if(is.null(con$clusters)) stop('No cell groups specified and no clusterings found')
  return(as.factor(con$clusters[[1]]$groups))
}

extractCellGroups.Seurat <- function(so) {
  groups <- Seurat::Idents(so)
  if(length(unique(groups)) <= 1) stop('No cell groups specified and no clusterings found')
  return(as.factor(groups))
}

extractCellGraph <- function(obj) UseMethod("extractCellGraph", obj)
extractCellGraph.Conos <- function(con) {
  if(is.null(con$graph)) stop('No cell graph found in the object')
  return(con$graph)
}

extractCellGraph.Seurat <- function(so) {
  graph.name <- so@misc$graph.name
  if (!is.null(graph.name)) {
    if (is.null(so@graphs[[graph.name]])) stop("Can't find graph ", graph.name)
    graph <- so@graphs[[graph.name]]
  } else {
    if(length(so@graphs) == 0) stop('No cell graph found in the object')
    graph <- so@graphs[[1]]
  }

  graph %<>% as("dgCMatrix") %>%
    igraph::graph_from_adjacency_matrix(mode="undirected", weighted=TRUE)
  return(graph)
}

extractRawCountMatrices <- function(obj, transposed=TRUE) UseMethod("extractRawCountMatrices", obj)
extractRawCountMatrices.Conos <- function(con, transposed=TRUE) {
  return(lapply(con$samples, conos:::getRawCountMatrix, transposed=transposed))
}

extractRawCountMatrices.Seurat <- function(so, transposed=TRUE) {
  cms <- so$sample.per.cell %>% {split(names(.), .)} %>%
    lapply(function(cids) so@assays$RNA@counts[,cids])
  if (transposed) {
    cms %<>% lapply(Matrix::t)
  }
  return(cms)
}

extractRawCountMatrices.dgCMatrix <- function(cm, transposed=TRUE) {
  sample.per.cell <- attr(cm, 'sample.per.cell')
  cms <- sample.per.cell %>% {split(names(.), .)} %>%
    lapply(function(cids) cm[cids,])
  if (!transposed) {
    cms %<>% lapply(Matrix::t)
  }
  return(cms)
}

extractJointCountMatrix <- function(obj, raw=TRUE, ...) UseMethod("extractJointCountMatrix", obj)
extractJointCountMatrix.Conos <- function(con, raw=TRUE) {
  return(con$getJointCountMatrix(raw=raw))
}

extractJointCountMatrix.Seurat <- function(so, raw=TRUE, transposed=TRUE, sparse=TRUE) {
  if (raw) {
    dat <- so@assays$RNA@counts
    if (transposed) dat %<>% Matrix::t()
    return(dat)
  }

  dat <- Seurat::GetAssayData(so, slot='scale.data', assay='RNA')
  dims <- dim(dat)
  dat.na <- all(dims == 1) && all(is.na(x = dat))
  if (all(dims == 0) || dat.na) {
    dat <- Seurat::GetAssayData(so, slot='data', assay='RNA')
  }

  if (transposed) dat %<>% Matrix::t()
  if (is.matrix(dat) && sparse) dat %<>% as("dgCMatrix")
  return(dat)
}

extractJointCountMatrix.dgCMatrix <- function(cm, raw=TRUE) {
  if (attr(cm, 'raw') == raw)
    return(cm)

  if (!attr(cm, 'raw') && raw)
    stop("Can't extract raw matrix from a normalized dgCMatrix")

  return(t(t(cm) / pmax(colSums(cm), 0.1)))
}

extractOdGenes <- function(obj, n.genes=NULL) UseMethod("extractOdGenes", obj)
extractOdGenes.Conos <- function(con, n.genes=NULL) {
  return(conos:::getOdGenesUniformly(con$samples, n.genes))
}

extractOdGenes.Seurat <- function(so, n.genes=NULL) {
  if (is.null(so@assays$RNA@meta.features$vst.variance.standardized))
    stop("The data object doesn't have gene variance info.",
         "Please, run FindVariableFeatures(assay='RNA, selection.method='vst') first")
  genes <- so@assays$RNA@meta.features %>%
    {rownames(.)[order(.$vst.variance.standardized, decreasing=TRUE)]} %>%
    head(n.genes %||% length(.))

  return(genes)
}

extractSamplePerCell <- function(obj) UseMethod("extractSamplePerCell", obj)
extractSamplePerCell.Conos <- function(con) {
  return(con$getDatasetPerCell())
}

extractSamplePerCell.Seurat <- function(so) {
  return(so$sample.per.cell)
}

extractSampleGroups <- function(obj, ref.level, target.level) {
  samp.names <- unique(extractSamplePerCell(obj))
  if(!any(grep(ref.level,con.names)))
    stop("'ref.level' not in the data object sample names.")
  sg <- grepl(ref.level, samp.names) %>% ifelse(ref.level, target.level) %>%
    setNames(samp.names)
  return(sg)
}

extractEmbedding <- function(obj) UseMethod("extractEmbedding", obj)
extractEmbedding.Conos <- function(con) {
  return(con$embedding)
}

extractEmbedding.Seurat <- function(so) {
  if (!is.null(so@reductions$umap)) return(so@reductions$umap@cell.embeddings)
  if (length(so@reductions) == 0) return(NULL)
  return(so@reductions[[1]]@cell.embeddings)
}

extractGeneExpression <- function(obj, gene) UseMethod("extractGeneExpression", obj)
extractGeneExpression.Conos <- function(obj, gene) {
  return(conos:::getGeneExpression(obj, gene))
}

extractGeneExpression.Seurat <- function(so, gene) {
  return(extractJointCountMatrix(so, raw=FALSE, transposed=FALSE, sparse=FALSE)[gene,])
}

extractGeneExpression.dgCMatrix <- function(cm, gene) {
  return(cm[,gene])
}

