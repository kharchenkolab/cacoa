##' @title Validate parameters per cell type
##' @param raw.mats List of raw count matrices
##' @param cell.groups Named clustering/annotation factor with cell names
##' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
##' @param ref.level Reference cluster level in 'sample.groups', e.g., ctrl, healthy, wt
##' @param cluster.sep.chr Character string of length 1 specifying a delimiter to separate cluster and app names
validatePerCellTypeParamsCacoa <- function(raw.mats, cell.groups, sample.groups, ref.level, cluster.sep.chr) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("You have to install DESeq2 package to use differential expression")
  }

  if ( is.null(cell.groups) ) stop('"cell.groups" must be specified');
  if ( is.null(sample.groups) ) stop('"sample.groups" must be specified')
  if ( class(sample.groups) != 'list' ) stop('"sample.groups" must be a list');
  if ( length(sample.groups) != 2 ) stop('"sample.groups" must be of length 2');
  if ( ! all(unlist(lapply(sample.groups, function(x) class(x) == 'character'))) )
    stop('"sample.groups" must be a list of character vectors');
  if ( ! all(unlist(lapply(sample.groups, function(x) length(x) > 0))) )
    stop('"sample.groups" entries must be on length greater or equal to 1')
  if ( ! all(unlist(lapply(sample.groups, function(x) {all(x %in% names(raw.mats))}))) )
    stop('"sample.groups" entries must be names of samples in the raw.mats')
  if ( is.null(ref.level) ) stop('"ref.level" is not defined')
  ## todo: check samplegrousp are named
  if(is.null(names(sample.groups))) stop('"sample.groups" must be named')
  if(class(cell.groups) != 'factor') stop('"cell.groups" must be a factor')
  if(any(grepl(cluster.sep.chr, names(raw.mats),fixed=TRUE)))
    stop('"cluster.sep.chr" must not be part of any sample name')
  if(any(grepl(cluster.sep.chr,levels(cell.groups),fixed=TRUE)))
    stop('"cluster.sep.chr" must not be part of any cluster name')
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
#' @param cell.groups factor specifying cell types (default=NULL)
#' @param sample.groups a list of two character vector specifying the app groups to compare (default=NULL)
#' @param cooks.cutoff cooksCutoff for DESeq2 (default=F)
#' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
#' @param min.cell.count (default=10)
#' @param independent.filtering independentFiltering for DESeq2 (default=F)
#' @param n.cores Number of cores (default=1)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details Return details
#' @export
getPerCellTypeDEmat=function (raw.mats, cell.groups = NULL, sample.groups = NULL, cooks.cutoff = FALSE,
                              ref.level = NULL, min.cell.count = 10, independent.filtering = FALSE,
                              n.cores = 1, cluster.sep.chr = "<!!>", return.details = F, verbose = T) {
  validatePerCellTypeParamsCacoa(raw.mats, cell.groups, sample.groups, ref.level, cluster.sep.chr)

  # TODO shouldn't depend on Conos
  aggr2<- rawMatricesWithCommonGenesCacoa(raw.mats, sample.groups) %>%
    lapply(conos:::collapseCellsByType, groups = cell.groups, min.cell.count = min.cell.count) %>%
    conos:::rbindDEMatrices(cluster.sep.chr = cluster.sep.chr)

  de.res <- sccore:::plapply(sn(levels(cell.groups)), function(l) {
    tryCatch({
      cm <- aggr2[, conos:::strpart(colnames(aggr2), cluster.sep.chr,
                                    2, fixed = TRUE) == l]
      meta <- data.frame(sample.id=colnames(cm), group=as.factor(unlist(lapply(colnames(cm), function(y) {
        y <- conos:::strpart(y, cluster.sep.chr, 1, fixed = TRUE)
        names(sample.groups)[unlist(lapply(sample.groups, function(x) any(x %in% y)))]}))))

      if (!ref.level %in% levels(meta$group))
        stop("The reference level is absent in this comparison")
      meta$group <- relevel(meta$group, ref = ref.level)
      if (length(unique(as.character(meta$group))) < 2)
        stop("The cluster is not present in both conditions")

      dds1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta,
                                             design = ~group)
      dds1 <- DESeq2::DESeq(dds1, quiet=T)
      res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff,
                              independentFiltering = independent.filtering)
      res1 <- as.data.frame(res1)
      res1 <- res1[order(res1$padj, decreasing = FALSE),]

      if (return.details) {
        list(res = res1, cm = cm)
      }
      else {
        res1
      }
    }, error = function(err) NA)
  }, n.cores = n.cores, progress=verbose) %>%
    .[!sapply(., `[[`, 1) %>% is.na]

  dif <- length(levels(cell.groups)) - length(de.res)
  if(dif > 0) warning(paste0("DEs not calculated for ",dif," cell group(s)."))

  return(de.res)
}

strpart <- function (x, split, n, fixed = FALSE) {
  sapply(strsplit(as.character(x), split, fixed = fixed), "[",n)
}
