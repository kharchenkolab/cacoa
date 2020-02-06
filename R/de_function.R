validatePerCellTypeParamsCacoa <- function(raw.mats, groups, sample.groups, ref.level, cluster.sep.chr) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("You have to install DESeq2 package to use differential expression")
  }
  
  if ( is.null(groups) ) stop('groups must be specified');
  if ( is.null(sample.groups) ) stop('sample.groups must be specified')
  if ( class(sample.groups) != 'list' ) stop('sample.groups must be a list');
  if ( length(sample.groups) != 2 ) stop('sample.groups must be of length 2');
  if ( ! all(unlist(lapply(sample.groups, function(x) class(x) == 'character'))) )
    stop('sample.groups must be a list of character vectors');
  if ( ! all(unlist(lapply(sample.groups, function(x) length(x) > 0))) )
    stop('sample.groups entries must be on length greater or equal to 1')
  if ( ! all(unlist(lapply(sample.groups, function(x) {all(x %in% names(raw.mats))}))) )
    stop('sample.groups entries must be names of samples in the raw.mats')
  if ( is.null(ref.level) ) stop('reference level is not defined')
  ## todo: check samplegrousp are named
  if(is.null(names(sample.groups))) stop('sample.groups must be named')
  if(class(groups) != 'factor') stop('groups must be a factor')
  if(any(grepl(cluster.sep.chr, names(raw.mats),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any sample name')
  if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any cluster name')
}

rawMatricesWithCommonGenesCacoa=function (raw.mats, sample.groups = NULL) 
{
  if (!is.null(raw.mats)) {
    raw.mats <- raw.mats[unlist(sample.groups)]
  }
  common.genes <- Reduce(intersect, lapply(raw.mats, colnames))
  return(lapply(raw.mats, function(x) {
    x[, common.genes]
  }))
}


#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' @param raw.mats  list of counts matrix; row for gene and column for cell 
#' @param groups factor specifying cell types
#' @param sample.groups a list of two character vector specifying the app groups to compare
#' @param cooks.cutoff cooksCutoff for DESeq2
#' @param independent.filtering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details return detals
#' @export getPerCellTypeDE
getPerCellTypeDEmat=function (raw.mats, groups = NULL, sample.groups = NULL, cooks.cutoff = FALSE, 
                              ref.level = NULL, min.cell.count = 10, independent.filtering = FALSE, 
                              n.cores = 1, cluster.sep.chr = "<!!>", return.details = TRUE) 
{
  validatePerCellTypeParamsCacoa(raw.mats, groups, sample.groups, 
                                    ref.level, cluster.sep.chr)
  dat <- rawMatricesWithCommonGenesCacoa(raw.mats, sample.groups)
  aggr2=lapply(dat,conos:::collapseCellsByType, groups = groups, min.cell.count = min.cell.count) %>% 
    conos:::rbindDEMatrices(cluster.sep.chr = cluster.sep.chr)
  
  gc()
  de.res <- conos:::papply(sn(levels(groups)), function(l) {
    tryCatch({
      cm <- aggr2[, conos:::strpart(colnames(aggr2), cluster.sep.chr, 
                                    2, fixed = TRUE) == l]
      meta <- data.frame(sample.id = colnames(cm), group = as.factor(unlist(lapply(colnames(cm), 
                                                                                   function(y) {
                                                                                     y <- conos:::strpart(y, cluster.sep.chr, 1, fixed = TRUE)
                                                                                     names(sample.groups)[unlist(lapply(sample.groups, 
                                                                                                                        function(x) any(x %in% y)))]
                                                                                   }))))
      if (!ref.level %in% levels(meta$group)) 
        stop("The reference level is absent in this comparison")
      meta$group <- relevel(meta$group, ref = ref.level)
      if (length(unique(as.character(meta$group))) < 2) 
        stop("The cluster is not present in both conditions")
      
      #cm=cm[rowMeans(as.matrix(cm))!=0,]
      dds1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta, 
                                             design = ~group)
      dds1 <- DESeq2::DESeq(dds1)
      res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, 
                              independentFiltering = independent.filtering)
      res1 <- as.data.frame(res1)
      res1 <- res1[order(res1$padj, decreasing = FALSE), 
                   ]
      
      # res1$ExpressionFraction=ExpressionFraction(dat,l,cm,groups)[rownames(res1)]
      
      if (return.details) {
        list(res = res1, cm = cm, sample.groups = sample.groups)
      }
      else {
        res1
      }
    }, error = function(err) NA)
  }, n.cores = n.cores)
  de.res
}




strpart <- function (x, split, n, fixed = FALSE) {
  sapply(strsplit(as.character(x), split, fixed = fixed), "[",n)
}



sn=function(x) { names(x) <- x; return(x); }

