# accessor methods for uniform access to data from Conos, Seurat and other objects

setGeneric("getCellNames", function(sample) standardGeneric("getCellNames"))
setMethod("getCellNames", signature("seurat"), function(sample) colnames(sample@data))
setMethod(f = 'getCellNames', signature = signature('Seurat'), definition = function(sample) return(colnames(x = sample)))
setMethod("getCellNames", signature("Conos"), function(sample) unlist(lapply(sample$samples,getCellNames)))

setGeneric("getGenes", function(sample) standardGeneric("getGenes"))
setMethod("getGenes", signature("seurat"), function(sample) rownames(sample@data))
setMethod(f = 'getGenes', signature = signature('Seurat'), definition = function(sample) return(rownames(x = sample)))
setMethod("getGenes", signature("Conos"), function(sample) unique(unlist(lapply(sample$samples,getGenes))))

setGeneric("getCountMatrix", function(sample, transposed=F) standardGeneric("getCountMatrix"))
setMethod("getCountMatrix", signature("Conos"), function(sample, transposed=F) {
  if (transposed)
    return(Matrix::t(sample$counts))

  return(sample$counts)
})
setMethod("getCountMatrix", signature("seurat"), function(sample, transposed=F) {
  cm <- if (is.null(sample@scale.data)) sample@data else sample@scale.data
  if (transposed)
    return(Matrix::t(cm))

  return(cm)
})
setMethod('getCountMatrix', signature('Seurat'), function(sample, transposed=F) {
    checkSeuratV3()
    dat <- Seurat::GetAssayData(object = sample, slot = 'scale.data')
    dims <- dim(x = dat)
    dat.na <- all(dims == 1) && all(is.na(x = dat))
    if (all(dims == 0) || dat.na) {
      dat <- Seurat::GetAssayData(object = sample, slot = 'data')
    }

    if (transposed)
      return(Matrix::t(dat))

    return(dat)
  }
)

setGeneric("getGeneExpression", function(sample, gene) standardGeneric("getGeneExpression"))
setMethod("getGeneExpression", signature("Conos"), function(sample, gene) {
  lapply(sample$samples, getGeneExpression, gene) %>% Reduce(c, .)
})

getGeneExpression.default <- function(sample, gene) {
  count.matrix <- getCountMatrix(sample)
  if(gene %in% rownames(count.matrix)) {
    return(count.matrix[gene,])
  }

  return(stats::setNames(rep(NA, ncol(count.matrix)), colnames(count.matrix)))
}

setGeneric("getRawCountMatrix", function(sample, transposed=F) standardGeneric("getRawCountMatrix"))
setMethod(
  f = "getRawCountMatrix",
  signature = signature("seurat"),
  definition = function(sample, transposed = F) {
    mi <- match(x = sample@cell.names, table = colnames(sample@raw.data))
    x <- sample@raw.data[, mi, drop = FALSE]
    if (transposed) {
      return(t(x = x))
    } else {
      return(x)
    }
  }
)
setMethod(
  f = 'getRawCountMatrix',
  signature = signature('Seurat'),
  definition = function(sample, transposed = FALSE) {
    checkSeuratV3()
    rd <- Seurat::GetAssayData(object = sample, slot = 'counts')
    # Raw data can be empty in Seurat v3
    # If it is, use data instead
    dims <- dim(x = rd)
    rd.na <- all(dims == 1) && all(is.na(x = rd))
    if (all(dims == 0) || rd.na) {
      rd <- Seurat::GetAssayData(object = sample, slot = 'data')
    }
    mi <- match(x = colnames(x = sample), table = colnames(x = rd))
    rd <- rd[, mi, drop = FALSE]
    if (transposed) {
      rd <- t(x = rd)
    }
    return(rd)
  }
)
setMethod("getRawCountMatrix", signature("Conos"), function(sample, transposed=F) { m <- sample$getJointCountMatrix(raw=TRUE); if (transposed) t(m) else m })

setGeneric("getEmbedding", function(sample, type) standardGeneric("getEmbedding"))
setMethod("getEmbedding", signature("seurat"), function(sample, type) if (is.null(sample@dr[[type]])) NULL else as.data.frame(sample@dr[[type]]@cell.embeddings))
setMethod(
  f = 'getEmbedding',
  signature = signature('Seurat'),
  definition = function(sample, type) {
    checkSeuratV3()
    emb <- tryCatch(
      expr = Seurat::Embeddings(object = sample, reduction = type),
      error = function(...) {
        return(NULL)
      }
    )
    return(emb)
  }
)
setMethod("getEmbedding", signature("Conos"), function(sample, type) sample$embedding)

setGeneric("getClustering", function(sample, type=NULL) standardGeneric("getClustering"))
setMethod("getClustering", signature("seurat"), function(sample, type) {if (!is.null(type)) warning("Seurat support only single type of clustering"); sample@ident})
setMethod(
  f = 'getClustering',
  signature = signature('Seurat'),
  definition = function(sample, type) {
    checkSeuratV3()
    if (missing(x = type)) {
      type <- NULL
    } else if (!is.null(x = type) && !type %in% colnames(x = sample[[]])) {
      warning(
        "Cannot find ",
        type,
        " in sample metadata, using normal identities",
        call. = FALSE,
        immediate. = TRUE
      )
      type <- NULL
    }
    idents <- if (is.null(x = type)) {
      Seurat::Idents(object = sample)
    } else {
      ids <- sample[[type]]
      if (!is.factor(x = ids)) {
        ids <- factor(x = ids)
      }
      ids
    }
    return(idents)
  }
)
setMethod("getClustering", signature("Conos"), function(sample, type=NULL) {if(is.null(type)) { type <- 1 }; cl <- sample$clusters[[type]]; if(is.null(cl)) NULL else cl$groups })
