#' @import sccore
NULL

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

rawMatricesWithCommonGenes=function(raw.mats, sample.groups = NULL)
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
  # TODO remove dependency on Conos
  tc <- conos:::colSumByFactor(cm,cl);
  tc <- tc[-1,,drop=FALSE]  # omit NA cells
  tc[table(cl)>=min.cell.count,,drop=FALSE]
}

#' Add Z scores to DE results
#' @param df Data frame with the columns "pval", "padj" and "log2FoldChange"
#' @return Updated data frame with Z scores
#' @export
addZScores <- function(df) {
  df$Z <- -qnorm(df$pval/2)
  df$Z[is.na(df$Z)] <- 0
  df$Za <- -qnorm(df$padj/2)
  df$Za[is.na(df$Za)] <- 0
  df$Z <- df$Z * sign(df$log2FoldChange)
  df$Za <- df$Za * sign(df$log2FoldChange)

  return(df)
}

#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' @param raw.mats list of counts matrices; column for gene and row for cell
#' @param cell.groups factor specifying cell types (default=NULL)
#' @param sample.groups a list of two character vector specifying the app groups to compare (default=NULL)
#' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
#' @param common.genes Only investigate common genes across cell groups (default=F)
#' @param cooks.cutoff cooksCutoff for DESeq2 (default=F)
#' @param min.cell.count (default=10)
#' @param independent.filtering independentFiltering for DESeq2 (default=F)
#' @param n.cores Number of cores (default=1)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
#' @param return.matrix Return merged matrix of results (default=F)
#' @export
estimatePerCellTypeDE=function (raw.mats, cell.groups = NULL, sample.groups = NULL, ref.level = NULL,
                           common.genes = F, cooks.cutoff = FALSE, min.cell.count = 10, independent.filtering = T,
                           n.cores = 1, cluster.sep.chr = "<!!>", return.matrix = F, verbose = T) {
  validatePerCellTypeParams(raw.mats, cell.groups, sample.groups, ref.level, cluster.sep.chr)

  if(common.genes) {
    raw.mats <- rawMatricesWithCommonGenes(raw.mats, sample.groups)
  } else {
    gene.union <- lapply(raw.mats, colnames) %>% Reduce(union, .)
    raw.mats <- sccore:::plapply(raw.mats, sccore:::extendMatrix, gene.union, n.cores = n.cores)
  }

  aggr2 <- raw.mats %>%
    lapply(collapseCellsByType, groups = cell.groups, min.cell.count = min.cell.count) %>%
    rbindDEMatrices(cluster.sep.chr = cluster.sep.chr)

  de.res <- cell.groups %>%
    levels() %>%
    sccore:::sn() %>%
    sccore:::plapply(function(l) {
    tryCatch({
      cm <- aggr2[, strpart(colnames(aggr2), cluster.sep.chr, 2, fixed = TRUE) == l] %>%
        .[rowSums(.) > 0,] # Remove genes with no counts

      meta <- data.frame(sample.id=colnames(cm), group=as.factor(unlist(lapply(colnames(cm), function(y) {
        y <- strpart(y, cluster.sep.chr, 1, fixed = TRUE)
        names(sample.groups)[unlist(lapply(sample.groups, function(x) any(x %in% y)))]}))))

      if (!ref.level %in% levels(meta$group))
        stop("The reference level is absent in this comparison")
      meta$group <- relevel(meta$group, ref = ref.level)
      if (length(unique(as.character(meta$group))) < 2)
        stop("The cluster is not present in both conditions")

      res1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta, design = ~group) %>%
        DESeq2::DESeq(quiet=T) %>%
        DESeq2::results(cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering) %>%
        as.data.frame

      # add Z scores
      if(!is.na(res1[[1]][1])) {
        res1 <- addZScores(res1) %>%
          .[order(abs(.$pvalue),decreasing=F),]
      }

      if (return.matrix) {
        list(res = res1, cm = cm)
      }
      else {
        res1
      }
    }, error = function(err) NA)
  }, n.cores = n.cores, progress=verbose) %>%
    .[!sapply(., is.logical)]


  if(verbose) {
    dif <- setdiff(levels(cell.groups), names(de.res))
    if(length(dif) > 0) {
      message(paste0("\nDEs not calculated for ",length(dif)," cell group(s):"))
      print(dif)
    }
  }

  return(de.res)
}

#' @title Save DE results as JSON files
#' @param de.raw List of DE results
#' @param saveprefix Prefix for created files (default=NULL)
#' @param dir.name Name for directory with results. If it doesn't exist, it will be created. To disable, set as NULL (default="JSON")
#' @param gene.metadata (default=NULL)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
saveDEasJSON <- function(de.raw, saveprefix = NULL, dir.name = "JSON", gene.metadata = NULL, cluster.sep.chr = "<!!>")
{
  if(!is.null(dir.name)) {
    if(!dir.exists(dir.name)) dir.create(dir.name)
  } else {
    dir.name = "."
  }

  lapply(sccore:::sn(de.raw %>% names()), function(ncc) {
    res.table <- de.raw[[ncc]]
    res.table$gene <- rownames(res.table)
    res.table$significant <- res.table$padj < 0.05
    res.table$log2FoldChange[is.na(res.table$log2FoldChange)] <- 0
    res.table$rowid <- 1:nrow(res.table)
    if (!is.null(gene.metadata)) {
      mo <- match(as.character(gene.metadata$geneid), as.character(res.table$gene))
      keep.cols <- colnames(gene.metadata)[colnames(gene.metadata) != "geneid"]
      names(keep.cols) <- keep.cols
      res.table <- cbind(res.table, gene.metadata[mo, keep.cols, drop = FALSE])
    }

    tojson <- list(res = res.table, genes = rownames(res.table))
    y <- jsonlite::toJSON(tojson)
    file <- paste0(dir.name, "/", saveprefix, make.names(ncc), ".json")
    write(y, file)
    NULL
  })
  invisible(NULL)

  toc.file <- paste0(dir.name,"/toc.html")
  s <- paste(c(list('<html><head><style>
    table {
    font-family: arial, sans-serif;
    border-collapse: collapse;
    width: 100%;
    }

    td, th {
    border: 1px solid #dddddd;
    text-align: left;
    padding: 8px;
    }

    tr:nth-child(even) {
    background-color: #dddddd;
    }

    </style></head><body><table>'),lapply(names(de.raw),function(n) paste0('<tr><td><a href="deview.2.html?d=',saveprefix, make.names(n),'.json">',n,'</a></td></tr>')),list('</table></body></html>')),collapse='\n')

  write(s,file=toc.file)
}
