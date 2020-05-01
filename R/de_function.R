##' @title Validate parameters per cell type
##' @param raw.mats List of raw count matrices
##' @param cell.groups Named clustering/annotation factor with cell names
##' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
##' @param ref.level Reference cluster level in 'sample.groups', e.g., ctrl, healthy, wt
##' @param cluster.sep.chr Character string of length 1 specifying a delimiter to separate cluster and app names
validatePerCellTypeParams <- function(raw.mats, cell.groups, sample.groups, ref.level, cluster.sep.chr) {
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

rawMatricesWithCommonGenes=function (raw.mats, sample.groups = NULL)
{
  if (!is.null(raw.mats)) {
    raw.mats <- raw.mats[unlist(sample.groups)]
  }
  common.genes <- Reduce(intersect, lapply(raw.mats, colnames))
  return(lapply(raw.mats, function(x) {
    x[, common.genes]
  }))
}

strpart <- function (x, split, n, fixed = FALSE) {
  sapply(strsplit(as.character(x), split, fixed = fixed), "[",n)
}

rbindDEMatrices <- function(mats, cluster.sep.chr) {
  mats <- lapply(names(mats), function(n) {
    rownames(mats[[n]]) <- paste0(n, cluster.sep.chr, rownames(mats[[n]]));
    return(mats[[n]])
  })

  return(t(do.call(rbind, mats)))
}

collapseCellsByType <- function(cm, groups, min.cell.count=10) {
  groups <- as.factor(groups);
  cl <- factor(groups[match(rownames(cm),names(groups))],levels=levels(groups));
  tc <- conos:::colSumByFactor(cm,cl);
  tc <- tc[-1,,drop=FALSE]  # omit NA cells
  tc[table(cl)>=min.cell.count,,drop=FALSE]
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
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
#' @param return.matrix Return merged matrix of results (default=F)
#' @export
getPerCellTypeDE=function (raw.mats, cell.groups = NULL, sample.groups = NULL, cooks.cutoff = FALSE,
                              ref.level = NULL, min.cell.count = 10, independent.filtering = FALSE,
                              n.cores = 1, cluster.sep.chr = "<!!>", return.matrix = F, verbose = T) {
  validatePerCellTypeParams(raw.mats, cell.groups, sample.groups, ref.level, cluster.sep.chr)

  # TODO genes should only be common per cell type - remember background per cell type
  aggr2<- rawMatricesWithCommonGenes(raw.mats, sample.groups) %>%
    lapply(collapseCellsByType, groups = cell.groups, min.cell.count = min.cell.count) %>%
    rbindDEMatrices(cluster.sep.chr = cluster.sep.chr)

  de.res <- sccore:::plapply(sccore:::sn(levels(cell.groups)), function(l) {
    tryCatch({
      cm <- aggr2[, nstrpart(colnames(aggr2), cluster.sep.chr,
                                    2, fixed = TRUE) == l]
      meta <- data.frame(sample.id=colnames(cm), group=as.factor(unlist(lapply(colnames(cm), function(y) {
        y <- strpart(y, cluster.sep.chr, 1, fixed = TRUE)
        names(sample.groups)[unlist(lapply(sample.groups, function(x) any(x %in% y)))]}))))

      if (!ref.level %in% levels(meta$group))
        stop("The reference level is absent in this comparison")
      meta$group <- relevel(meta$group, ref = ref.level)
      if (length(unique(as.character(meta$group))) < 2)
        stop("The cluster is not present in both conditions")

      res1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta, design = ~group) %>%
        DESeq2::DESeq(dds1, quiet=T) %>%
        DESeq2::results(cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering) %>%
        as.data.frame

      # add Z scores
      res1$Z <- -qnorm(res1$pval/2)
      res1$Z[is.na(res1$Z)] <- 0
      res1$Za <- -qnorm(res1$padj/2)
      res1$Za[is.na(res1$Za)] <- 0
      res1$Z <- res1$Z  * sign(res1$log2FoldChange)
      res1$Za <- res1$Za  * sign(res1$log2FoldChange)

      res1 <- list(down=rownames(res1)[order(res1$Z,decreasing=F)],
                  up=rownames(res1)[order(res1$Z,decreasing=T)],
                  all=rownames(res1)[order(abs(res1$Z),decreasing=T)])

      if (return.matrix) {
        res <- list(res = res, cm = cm)
      }
      else {
        res
      }
    }, error = function(err) NA)
  }, n.cores = n.cores, progress=verbose) %>%
    .[!sapply(., `[[`, 1) %>% is.na]

  dif <- length(levels(cell.groups)) - length(de.res)
  if(dif > 0) warning(paste0("DEs not calculated for ",dif," cell group(s)."))

  return(res)
}

