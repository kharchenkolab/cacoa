#' @import sccore
NULL

##' @title Validate parameters per cell type
##' @param raw.mats List of raw count matrices
##' @param cell.groups Named clustering/annotation factor with cell names
##' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
##' @param ref.level Reference cluster level in 'sample.groups', e.g., ctrl, healthy, wt
##' @param cluster.sep.chr Character string of length 1 specifying a delimiter to separate cluster and app names
validatePerCellTypeParams <- function(raw.mats, cell.groups, sample.groups, ref.level, cluster.sep.chr) {
  checkPackageInstalled("DESeq2", bioc=TRUE)

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
#' @param common.genes Only investigate common genes across cell groups (default=FALSE)
#' @param test which DESeq2 test to use (options: "LRT" (default), "Wald")
#' @param cooks.cutoff cooksCutoff for DESeq2 (default=FALSE)
#' @param min.cell.count (default=10)
#' @param max.cell.count maximal number of cells per cluster per sample to include in a comparison (useful for comparing the number of DE genes between cell types)
#' @param independent.filtering independentFiltering for DESeq2 (default=FALSE)
#' @param n.cores Number of cores (default=1)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
#' @param return.matrix Return merged matrix of results (default=TRUE)
#' @export
estimatePerCellTypeDE=function (raw.mats, cell.groups = NULL, sample.groups = NULL, ref.level = NULL,
                           common.genes = FALSE, test="LRT", cooks.cutoff = FALSE, min.cell.count = 10,max.cell.count=Inf, independent.filtering = TRUE,
                           n.cores = 1, cluster.sep.chr = "<!!>", return.matrix = TRUE, verbose = TRUE) {

  validatePerCellTypeParams(raw.mats, cell.groups, sample.groups, ref.level, cluster.sep.chr)

  if(common.genes) {
    raw.mats <- rawMatricesWithCommonGenes(raw.mats, sample.groups)
  } else {
    gene.union <- lapply(raw.mats, colnames) %>% Reduce(union, .)
    raw.mats <- sccore:::plapply(raw.mats, sccore:::extendMatrix, gene.union, n.cores = n.cores)
  }

  aggr2 <- raw.mats %>%
    .[sample.groups %>% unlist] %>% # Only consider samples in sample.groups
    lapply(collapseCellsByType, groups=cell.groups, min.cell.count=min.cell.count, max.cell.count=max.cell.count) %>%
    .[sapply(., nrow) > 0] %>% # Remove empty samples due to min.cell.count
    rbindDEMatrices(cluster.sep.chr = cluster.sep.chr)

  # Adjust sample.groups
  passed.samples <- strpart(colnames(aggr2), cluster.sep.chr, 1, fixed = TRUE) %>% unique()
  if(verbose) {
    if(length(passed.samples) != length(unlist(sample.groups))) {
      warning("Excluded ",length(unlist(sample.groups)) - length(passed.samples)," sample(s) due to 'min.cell.count'.")
    }
  }
  sample.groups %<>% lapply(function(n) n[n %in% passed.samples])

  ## For every cell type get differential expression results
  de.res <- sccore::plapply( sccore::sn( levels(cell.groups) ), function(l) {
    tryCatch({
      ## Get count matrix
      cm <- aggr2[, strpart(colnames(aggr2), cluster.sep.chr, 2, fixed = TRUE) == l] %>% .[rowSums(.) > 0,] # Remove genes with no counts
      ## Generate metadata
      meta <- data.frame(
        sample.id=colnames(cm),
        group=as.factor(unlist(lapply(colnames(cm), function(y) {
          y <- strpart(y, cluster.sep.chr, 1, fixed = TRUE)
          names(sample.groups)[unlist(lapply(sample.groups, function(x) any(x %in% y)))]})))
      )

      if (!ref.level %in% levels(meta$group))  stop("The reference level is absent in this comparison")
      meta$group <- relevel(meta$group, ref = ref.level)
      if (length(unique(as.character(meta$group))) < 2)  stop("The cluster is not present in both conditions")

      dds1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta, design=~group)
      if(test=="LRT") {
        dds1 <- DESeq2::DESeq(dds1,test="LRT", reduced = ~ 1,quiet=TRUE)
      } else { # defaults to Wald
        dds1 <- DESeq2::DESeq(dds1,quiet=TRUE)
      }
      res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering) %>% as.data.frame

      # add Z scores
      if(!is.na(res1[[1]][1])) {
        res1 <- addZScores(res1) %>%
          .[order(.$pvalue,decreasing=FALSE),]
      }

      if (return.matrix) {
        list(res = res1, cm = cm)
      }
      else {
        res1
      }
    }, error = function(err) {if (verbose) {warning(err)}; NA})
  }, n.cores = n.cores, progress=verbose) %>%  .[!sapply(., is.logical)]


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
#' @param sample.groups Sample groups as named list, each element containing a vector of samples
#' @param saveprefix Prefix for created files (default=NULL)
#' @param dir.name Name for directory with results. If it doesn't exist, it will be created. To disable, set as NULL (default="JSON")
#' @param gene.metadata (default=NULL)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
#' @param verbose Show progress (default=TRUE)
saveDEasJSON <- function(de.raw, saveprefix = NULL, dir.name = "JSON", gene.metadata = NULL, cluster.sep.chr = "<!!>", sample.groups = NULL, verbose=TRUE) {
  if(!is.null(dir.name)) {
    if(!dir.exists(dir.name)) dir.create(dir.name)
  } else {
    dir.name = "."
  }

  if(is.null(gene.metadata)) {
    gene.metadata <- data.frame()
    all.genes <- unique(unlist(lapply(de.raw, function(x) {
      if(!is.null(x)){
        rownames(as.data.frame(x$res))
      } else {
        NULL
      }
    })))
    gene.metadata <- data.frame(geneid=all.genes)
  } else {
    if(is.null(gene.metadata$gene.id)) stop("gene.metadata must contain $gene.id field")
  }

  lapply(sccore:::sn(de.raw %>% names()), function(ncc) {
    if(verbose) print(ncc)
    res.celltype <- de.raw[[ncc]]
    res.table <- res.celltype$res %>% as.data.frame()
    res.table$gene <- rownames(res.table)
    res.table$significant <- res.table$padj < 0.05
    res.table$log2FoldChange[is.na(res.table$log2FoldChange)] <- 0
    res.table$rowid <- 1:nrow(res.table)

    all.genes <- rownames(res.table)
    cm <- res.celltype$cm
    colnames(cm) <- strpart(colnames(cm),cluster.sep.chr,1,fixed=TRUE)

    ilev <- lapply(sample.groups, function(sg) {
      sg <- sg[sg %in% colnames(cm)]
      cm.tmp <- cm[,sg]
      cm.tmp <- as.matrix(cm.tmp)
      rownames(cm.tmp) <- rownames(cm)

      ## calculate cpm
      cpm <- sweep(cm.tmp, 2, apply(cm.tmp,2, sum), FUN='/')
      cpm <- log10(cpm * 1e6 + 1)
      snames1 <- colnames(cpm)

      ## Put genes in order
      cpm <- cpm[all.genes,]
      colnames(cpm) <- NULL;
      rownames(cpm) <- NULL;

      list(snames=snames1, val=as.matrix(cpm))
    })

    snames <- names(res.celltype$sample.groups)

    ## convert to json
    tojson <- list(
      res = res.table,
      genes = all.genes,
      ilev = ilev,
      snames = snames)
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


#' Differential expression using different methods (deseq2, edgeR, wilcoxon, ttest) with various covariates
#' @param raw.mats list of counts matrices; column for gene and row for cell
#' @param cell.groups factor specifying cell types (default=NULL)
#' @param sample.groups a list of two character vector specifying the app groups to compare (default=NULL)
#' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
#' @param common.genes Only investigate common genes across cell groups (default=FALSE)
#' @param cooks.cutoff cooksCutoff for DESeq2 (default=FALSE)
#' @param min.cell.count (default=10)
#' @param independent.filtering independentFiltering for DESeq2 (default=FALSE)
#' @param n.cores Number of cores (default=1)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
#' @param return.matrix Return merged matrix of results (default=TRUE)
#' @param covariates list of covariates to include; for example, cdr, sex or age
#' @param meta.info dataframe with possible covariates; for example, sex or age
#' @param test DE method: deseq2, edgeR, wilcoxon, ttest
#' @export
estimatePerCellTypeDEmethods=function (raw.mats,
                                       cell.groups=NULL,
                                       s.groups=NULL,
                                       ref.level=NULL,
                                       target.level=NULL,
                                       common.genes=FALSE,
                                       cooks.cutoff=FALSE,
                                       min.cell.count=10,
                                       max.cell.count=Inf,
                                       independent.filtering=TRUE,
                                       n.cores=4,
                                       cluster.sep.chr="<!!>",
                                       return.matrix=TRUE,
                                       verbose=TRUE,
                                       test='Wald',
                                       meta.info=NULL,
                                       gene.filter=NULL) {

  validatePerCellTypeParams(raw.mats, cell.groups, s.groups, ref.level, cluster.sep.chr)

  if(common.genes) {
    raw.mats <- rawMatricesWithCommonGenes(raw.mats, s.groups)
  } else {
    gene.union <- lapply(raw.mats, colnames) %>% Reduce(union, .)
    raw.mats <- sccore:::plapply(raw.mats, sccore:::extendMatrix, gene.union, n.cores = n.cores)
  }

  aggr2 <- raw.mats %>%
    .[s.groups %>% unlist] %>% # Only consider samples in s.groups
    lapply(collapseCellsByType, groups=cell.groups, min.cell.count=min.cell.count, max.cell.count=max.cell.count) %>%
    .[sapply(., nrow) > 0] %>% # Remove empty samples due to min.cell.count
    rbindDEMatrices(cluster.sep.chr = cluster.sep.chr)
  mode(aggr2) <- 'integer'

  # Adjust s.groups
  passed.samples <- strpart(colnames(aggr2), cluster.sep.chr, 1, fixed = TRUE) %>% unique()
  if(verbose && (length(passed.samples) != length(unlist(s.groups))))
    warning("Excluded ",length(unlist(s.groups)) - length(passed.samples)," sample(s) due to 'min.cell.count'.")

  s.groups %<>% lapply(function(n) n[n %in% passed.samples])
  
  ## For every cell type get differential expression results
  de.res <- sccore::plapply( sccore::sn( levels(cell.groups) ), function(l) {
    tryCatch({
      ## Get count matrix
      cm <- aggr2[, strpart(colnames(aggr2), cluster.sep.chr, 2, fixed = TRUE) == l] %>% .[rowSums(.) > 0,] # Remove genes with no counts
      if(!is.null(gene.filter)) {
        gene.to.remain <- gene.filter %>% {rownames(.)[.[,l]]} %>% intersect(rownames(cm))
        cm <- cm[gene.to.remain,]
      }
      ## Generate metadata
      meta <- data.frame(
        sample.id=colnames(cm),
        group=as.factor(unlist(lapply(colnames(cm), function(y) {
          y <- strpart(y, cluster.sep.chr, 1, fixed = TRUE)
          names(s.groups)[unlist(lapply(s.groups, function(x) any(x %in% y)))]})))
      )

      if (!ref.level %in% levels(meta$group))  stop("The reference level is absent in this comparison")
      meta$group <- relevel(meta$group, ref = ref.level)
      if (length(unique(as.character(meta$group))) < 2)  stop("The cluster is not present in both conditions")

      ## External covariates
      if(is.null(meta.info)) {
        design.formula = as.formula('~ group')
      } else {
        design.formula = as.formula(paste('~ ',
                                          paste(c(colnames(meta.info), 'group'),
                                                collapse=' + ')))
        meta.info.tmp = meta.info[gsub(paste("<!!>", l, sep=''),"",meta[,1]),, drop=F]
        
        meta = cbind(meta, meta.info.tmp)
      }
      

      if(grepl('wilcoxon', tolower(test)) || grepl('t-test', tolower(test))) {

        # ----- Wilcoxon and t-test -----

        # Normalization
        tmp = strsplit(test, split = '\\.')
        test = tolower(tmp[[1]][1])
        normalization = tolower(tmp[[1]][2])

        if(normalization == 'deseq2') {
          print('DESeq2 normalization')
          cnts.norm <- cm  %>%
            DESeq2::DESeqDataSetFromMatrix(colData = meta, design= ~ group)  %>%
            DESeq2::estimateSizeFactors()  %>% DESeq2::counts(normalized=TRUE)
        } else if(normalization == 'edger') {
          # EdgeR normalisation
          print('edgeR normalization')
          cnts.norm <- edgeR::DGEList(counts = cm) %>%
            edgeR::calcNormFactors() %>% edgeR::cpm()
        } else if((normalization == 'totcount') || (is.null(normalization))) {
          # the default should be normalization by the number of molecules!
          cnts.norm <- prop.table(cm, 2) # Should it be multiplied by median(colSums(cm)) ?
        }

        cm = cnts.norm #remove

        if(test == 'wilcoxon') {
          # Wilcoxon test
          res1 <- scran::pairwiseWilcox(cnts.norm, groups = meta$group)$statistics[[1]] %>%
            data.frame() %>%
            setNames(c("AUC","pvalue","padj"))
        } else if (test == 't-test') {
          res1 <- scran::pairwiseTTests(cnts.norm, groups = meta$group)$statistics[[1]] %>%
            data.frame() %>%
            setNames(c("AUC","pvalue","padj"))
        }

        res1$log2FoldChange <- apply(log(cnts.norm+1, base=2), 1, function(x) {
          mean(x[meta$group == target.level]) - mean(x[meta$group == ref.level])})
      } else if(grepl('deseq2', tolower(test))) {
        # ----- DESeq2 -----

        test.name = 'Wald'
        if(grepl('lrt', tolower(strsplit(test, split = '\\.')[[1]][2]))) test.name = 'LRT'

        if(test.name == 'Wald') {
          res1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta, design=design.formula) %>%
            DESeq2::DESeq(quiet=TRUE, test=test.name) %>%
            DESeq2::results(contrast=c('group', target.level, ref.level),
                            cooksCutoff = cooks.cutoff,
                            independentFiltering = independent.filtering) %>%
            as.data.frame
        } else {
          res1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta, design=design.formula) %>%
            DESeq2::DESeq(quiet=TRUE, test=test.name, reduced = ~ 1) %>%
            DESeq2::results(contrast=c('group', target.level, ref.level),
                            cooksCutoff = cooks.cutoff,
                            independentFiltering = independent.filtering) %>%
            as.data.frame
        }

        # Avoid NA padj values
        res1$padj[is.na(res1$padj)] <- 1

      } else if(tolower(test) == tolower('edgeR')) {

        # ----- EdgeR -----
        design <- model.matrix(design.formula, meta)

        qlf <- edgeR::DGEList(cm, group = meta$group) %>%
          edgeR::calcNormFactors() %>%
          edgeR::estimateDisp(design = design) %>%
          edgeR::glmQLFit(design = design) %>%
          edgeR::glmQLFTest(coef=ncol(design))

        res1 <- qlf$table %>% .[order(.$PValue),]
        colnames(res1) <- c("log2FoldChange","logCPM","stat","pvalue")
        res1$padj <- p.adjust(res1$pvalue, method = "BH")

      } else if(tolower(test) == tolower('limma-voom')) {

        # ----- limma-voom -----

        mm <- model.matrix(design.formula, meta)

        # cnts.norm <- DGEList(counts = cm) %>%
        #   edgeR::calcNormFactors() %>% cpm
        cnts.norm = cm

        y <- limma::voom(cnts.norm, mm, plot = FALSE)
        fit <- limma:lmFit(y, mm)

        contr <- makeContrasts(paste(c('group', target.level), collapse = ''),
                               levels = colnames(coef(fit)))

        tmp <- contrasts.fit(fit, contr)
        tmp <- eBayes(tmp)

        res1 <- topTable(tmp, sort.by = "P", n = Inf)
        colnames(res1) <- c('log2FoldChange', 'AveExpr', 'stat', 'pvalue', 'padj', 'B')

      }

      # add Z scores
      if(!is.na(res1[[1]][1])) {
        res1 <- addZScores(res1) %>% .[order(.$pvalue,decreasing=FALSE),]
      }

      if (return.matrix) {
        list(res = res1, cm = cm)
      }
      else {
        res1
      }
    }, error = function(err) NA)
  }, n.cores = n.cores, progress=verbose) %>%  .[!sapply(., is.logical)]


  if(verbose) {
    dif <- setdiff(levels(cell.groups), names(de.res))
    if(length(dif) > 0) {
      message(paste0("\nDEs not calculated for ",length(dif)," cell group(s):"))
      print(dif)
    }
  }

  return(de.res)
}

#' Summarize DE Resampling Results
#' @param var.to.sort Variable to calculate ranks
summarizeDEResamplingResults <- function(de.list, var.to.sort='pvalue') {
  de.res <- de.list[[1]]
  for(cell.type in names(de.res)) {
    genes.init <- genes.common <- rownames(de.res[[cell.type]]$res)
    mx.stat <- matrix(nrow = length(genes.common), ncol = 0, dimnames = list(genes.common,c()))
    for(i in 2:length(de.list)){
      if(!(cell.type %in% names(de.list[[i]]))) next
      genes.common <- intersect(genes.common, rownames(de.list[[i]][[cell.type]]))
      mx.stat <- cbind(mx.stat[genes.common,], de.list[[i]][[cell.type]][genes.common, var.to.sort])
    }

    mx.stat <- apply(mx.stat, 2, rank)
    stab.mean.rank <- rowMeans(mx.stat) # stab - for stability
    stab.median.rank <- apply(mx.stat, 1, median)
    stab.var.rank <- apply(mx.stat, 1, var)

    de.res[[cell.type]]$res$stab.median.rank <- stab.median.rank[genes.init]
    de.res[[cell.type]]$res$stab.mean.rank <- stab.mean.rank[genes.init]
    de.res[[cell.type]]$res$stab.var.rank <- stab.var.rank[genes.init]

    # Save subsamples
    de.res[[cell.type]]$subsamples <- lapply(de.list[2:length(de.list)], `[[`, cell.type)
  }

  return(de.res)
}

appendStatisticsToDE <- function(de.list, expr.frac.per.type) {
  for (n in names(de.list)) {
    de.list[[n]]$res %<>% mutate(Gene=rownames(.), CellFrac=expr.frac.per.type[Gene, n],
                                 SampleFrac=Matrix::rowMeans(de.list[[n]]$cm > 0)[Gene]) %>%
      as.data.frame(stringsAsFactors=FALSE) %>% set_rownames(.$Gene)
  }

  return(de.list)
}

getExpressionFractionPerGroup <- function(cm, cell.groups) {
  cm@x <- as.numeric(cm@x > 1e-10)
  fracs <- collapseCellsByType(cm, cell.groups, min.cell.count=0) %>%
    {. / as.vector(table(cell.groups)[rownames(.)])} %>% Matrix::t()
  return(fracs)
}
