#' @import sccore
NULL

##' @title Validate parameters per cell type
##' @param raw.mats List of raw count matrices
##' @param cell.groups Named clustering/annotation factor with cell names
##' @param sample.groups Named factor with cell names indicating condition/sample, e.g., ctrl/disease
##' @param ref.level Reference cluster level in 'sample.groups', e.g., ctrl, healthy, wt
validateDEPerCellTypeParams <- function(raw.mats, cell.groups, sample.groups, ref.level) {
  checkPackageInstalled("DESeq2", bioc=TRUE)

  if (is.null(cell.groups)) stop('"cell.groups" must be specified')
  if (is.null(sample.groups)) stop('"sample.groups" must be specified')
  if (class(sample.groups) != "list") stop('"sample.groups" must be a list')
  if (length(sample.groups) != 2) stop('"sample.groups" must be of length 2')
  if (!all(unlist(lapply(sample.groups, function(x) class(x) == "character"))))
    stop('"sample.groups" must be a list of character vectors')

  if (!all(unlist(lapply(sample.groups, function(x) length(x) > 0))))
    stop('"sample.groups" entries must be on length greater or equal to 1')

  if (!all(unlist(lapply(sample.groups, function(x) all(x %in% names(raw.mats))))))
    stop('"sample.groups" entries must be names of samples in the raw.mats')

  if (is.null(ref.level)) stop('"ref.level" is not defined')
  ## todo: check samplegrousp are named
  if (is.null(names(sample.groups))) stop('"sample.groups" must be named')
  if (class(cell.groups) != "factor") stop('"cell.groups" must be a factor')
}

subsetMatricesWithCommonGenes <- function(cms, sample.groups=NULL) {
  if (!is.null(sample.groups)) cms <- cms[unlist(sample.groups)]
  common.genes <- do.call(intersect, lapply(cms, colnames))
  cms %<>% lapply(function(m) m[, common.genes, drop=FALSE])
  return(cms)
}

strpart <- function(x, split, n, fixed = FALSE) {
  as.character(x) %>% strsplit(split, fixed=fixed) %>% sapply("[", n)
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

#' @title Save DE results as JSON files
#' @param de.raw List of DE results
#' @param sample.groups Sample groups as named list, each element containing a vector of samples
#' @param saveprefix Prefix for created files (default=NULL)
#' @param dir.name Name for directory with results. If it doesn't exist, it will be created. To disable, set as NULL (default="JSON")
#' @param gene.metadata (default=NULL)
#' @param verbose Show progress (default=TRUE)
saveDEasJSON <- function(de.raw, saveprefix=NULL, dir.name="JSON", gene.metadata=NULL,
                         sample.groups=NULL, verbose=TRUE) {
  if(!is.null(dir.name)) {
    if(!dir.exists(dir.name)) dir.create(dir.name)
  } else {
    dir.name = "."
  }

  if (is.null(gene.metadata)) {
    gene.metadata <- data.frame()
    all.genes <- unique(unlist(lapply(de.raw, function(x) {
      if (!is.null(x)) {
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

  toc.file <- paste0(dir.name, "/toc.html")
  s <- c(list('<html><head><style>
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

    </style></head><body><table>'),
    lapply(names(de.raw), function(n)
      paste0('<tr><td><a href="deview.2.html?d=', saveprefix, make.names(n),'.json">', n, '</a></td></tr>')),
    list('</table></body></html>')
    ) %>% paste(collapse='\n')

  write(s, file=toc.file)
}

prepareSamplesForDE <- function(sample.groups, resampling.method=c('loo', 'bootstrap', 'fix.cells', 'fix.samples'),
                                n.resamplings=30, n.biosamples=NULL) {
  resampling.method <- match.arg(resampling.method)

  if (resampling.method == 'loo') {
    samples <- unlist(sample.groups) %>% sn() %>% lapply(function(n) lapply(sample.groups, setdiff, n))
  } else if (resampling.method == 'bootstrap') {
    # TODO: Do we ever use bootstrap? It seems that including the same sample many times
    # reduces variation and skews the analysis
    samples <- (1:n.resamplings) %>% setNames(paste0('bootstrap.', .)) %>%
      lapply(function(i) lapply(sample.groups, function(x) sample(x, length(x), replace=TRUE)))
  } else { # 'fix.cells' or 'fix.samples'
    samples <- (1:n.resamplings) %>% setNames(., paste0('fix.', .)) %>% lapply(function(i) sample.groups)
  }

  return(samples)
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
#' @param return.matrix Return merged matrix of results (default=TRUE)
#' @param covariates list of covariates to include; for example, cdr, sex or age
#' @param meta.info dataframe with possible covariates; for example, sex or age
#' @param test DE method: deseq2, edgeR, wilcoxon, ttest
#' @export
estimateDEPerCellTypeInner=function(raw.mats, cell.groups=NULL, s.groups=NULL, ref.level=NULL, target.level=NULL,
                                    common.genes=FALSE, cooks.cutoff=FALSE, min.cell.count=10, max.cell.count=Inf,
                                    independent.filtering=TRUE, n.cores=4, return.matrix=TRUE, fix.n.samples=NULL,
                                    verbose=TRUE, test='Wald', meta.info=NULL, gene.filter=NULL) {
  # Validate input
  validateDEPerCellTypeParams(raw.mats, cell.groups, s.groups, ref.level)
  tmp <- tolower(strsplit(test, split='\\.')[[1]])
  test <- tmp[1]
  test.type <- ifelse(is.na(tmp[2]), '', tmp[2])

  # Filter data and convert to the right format
  if (verbose) message("Preparing matrices for DE")
  if (common.genes) {
    raw.mats %<>% subsetMatricesWithCommonGenes(s.groups)
  } else {
    gene.union <- lapply(raw.mats, colnames) %>% Reduce(union, .)
    raw.mats %<>% plapply(sccore:::extendMatrix, gene.union, n.cores=n.cores, progress=verbose, mc.preschedule=TRUE)
  }

  cm.bulk.per.samp <- raw.mats[unlist(s.groups)] %>% # Only consider samples in s.groups
    lapply(collapseCellsByType, groups=cell.groups, min.cell.count=min.cell.count, max.cell.count=max.cell.count) %>%
    .[sapply(., nrow) > 0] # Remove empty samples due to min.cell.count

  cm.bulk.per.type <- levels(cell.groups) %>% sn() %>% lapply(function(cg) {
    tcms <- cm.bulk.per.samp %>%
      lapply(function(cm) if (cg %in% rownames(cm)) cm[cg, , drop=FALSE] else NULL) %>%
      .[!sapply(., is.null)]
    if (length(tcms) == 0) return(NULL)

    tcms %>% {set_rownames(do.call(rbind, .), names(.))} %>% `mode<-`('integer') %>%
      .[,colSums(.) > 0,drop=FALSE]
  }) %>% .[sapply(., length) > 0] %>% lapply(t)

  ## Adjust s.groups
  passed.samples <- names(cm.bulk.per.samp)
  if (verbose && (length(passed.samples) != length(unlist(s.groups))))
    warning("Excluded ", length(unlist(s.groups)) - length(passed.samples), " sample(s) due to 'min.cell.count'.")

  s.groups %<>% lapply(intersect, passed.samples)

  # For every cell type get differential expression results
  if (verbose) message("Estimating DE per cell type")
  de.res <- names(cm.bulk.per.type) %>% sn()%>% plapply(function(l) {
    tryCatch({
      cm <- cm.bulk.per.type[[l]]
      if (!is.null(gene.filter)) {
        gene.to.remain <- gene.filter %>% {rownames(.)[.[,l]]} %>% intersect(rownames(cm))
        cm <- cm[gene.to.remain,,drop=FALSE]
      }

      cur.s.groups <- lapply(s.groups, intersect, colnames(cm))
      if (!is.null(fix.n.samples)) {
        if (min(sapply(s.groups, length)) < fix.n.samples) stop("The cluster does not have enough samples")
        cur.s.groups %<>% lapply(sample, fix.n.samples)
        cm <- cm[, unlist(cur.s.groups), drop=FALSE]
      }

      ## Generate metadata
      meta.groups <- colnames(cm) %>% lapply(function(y) {
        names(cur.s.groups)[sapply(cur.s.groups, function(x) any(x %in% y))]
      }) %>% unlist() %>% as.factor()

      if (length(levels(meta.groups)) < 2) stop("The cluster is not present in both conditions")
      if (!ref.level %in% levels(meta.groups)) stop("The reference level is absent in this comparison")

      meta <- data.frame(sample.id=colnames(cm), group=relevel(meta.groups, ref=ref.level))

      ## External covariates
      if (is.null(meta.info)) {
        design.formula <- as.formula('~ group')
      } else {
        design.formula <- c(colnames(meta.info), 'group') %>%
          paste(collapse=' + ') %>% {paste('~', .)} %>% as.formula()
        meta %<>% cbind(meta.info[meta$sample.id, , drop=FALSE])
      }

      if (test %in% c('wilcoxon', 't-test')) {
        cm <- normalizePseudoBulkMatrix(cm, meta=meta, design.formula=design.formula, type=test.type)
        res <- estimateDEForTypePairwiseStat(cm, meta=meta, target.level=target.level, test=test)
      } else if (test == 'deseq2') {
        res <- estimateDEForTypeDESeq(
          cm, meta=meta, design.formula=design.formula, ref.level=ref.level, target.level=target.level,
          test.type=test.type, cooksCutoff=cooks.cutoff, independentFiltering=independent.filtering
        )
      } else if (test == 'edger') {
        res <- estimateDEForTypeEdgeR(cm, meta=meta, design.formula=design.formula)
      } else if (test == 'limma-voom') {
        res <- estimateDEForTypeLimma(cm, meta=meta, design.formula=design.formula, target.level=target.level)
      }

      if (!is.na(res[[1]][1])) {
        res <- addZScores(res) %>% .[order(.$pvalue, decreasing=FALSE),]
      }

      if (return.matrix)
        return(list(res = res, cm = cm, meta=meta))

      return(res)
    }, error = function(err) NA)
  }, n.cores=n.cores, progress=verbose) %>%  .[!sapply(., is.logical)]


  if (verbose) {
    dif <- setdiff(levels(cell.groups), names(de.res))
    if (length(dif) > 0) {
      message("DEs not calculated for ", length(dif), " cell group(s): ", paste(dif, collapse=', '))
    }
  }

  return(de.res)
}

normalizePseudoBulkMatrix <- function(cm, meta=NULL, design.formula=NULL, type='totcount') {
  if (type == 'deseq2') {
    cnts.norm <- DESeq2::DESeqDataSetFromMatrix(cm, meta, design=design.formula)  %>%
      DESeq2::estimateSizeFactors()  %>% DESeq2::counts(normalized=TRUE)
  } else if (type == 'edger') {
    cnts.norm <- edgeR::DGEList(counts=cm) %>% edgeR::calcNormFactors() %>% edgeR::cpm()
  } else if (type == 'totcount') {
    # the default should be normalization by the number of molecules!
    cnts.norm <- prop.table(cm, 2) # Should it be multiplied by median(colSums(cm)) ?
  }

  return(cnts.norm)
}

estimateDEForTypePairwiseStat <- function(cm.norm, meta, target.level, test) {
  if (test == 'wilcoxon') {
    res <- scran::pairwiseWilcox(cm.norm, groups = meta$group)$statistics[[1]] %>%
      data.frame() %>% setNames(c("AUC", "pvalue", "padj"))
  } else if (test == 't-test') {
    res <- scran::pairwiseTTests(cm.norm, groups = meta$group)$statistics[[1]] %>%
      data.frame() %>% setNames(c("AUC", "pvalue", "padj"))
  }

  # TODO: log2(x + 1) does not work for total-count normalization
  res$log2FoldChange <- log2(cm.norm + 1) %>% apply(1, function(x) {
    mean(x[meta$group == target.level]) - mean(x[meta$group != target.level])})

  return(res)
}

estimateDEForTypeDESeq <- function(cm, meta, design.formula, ref.level, target.level, test.type, ...) {
  res <- DESeq2::DESeqDataSetFromMatrix(cm, meta, design=design.formula)
  if (test.type == 'wald') {
      res %<>% DESeq2::DESeq(quiet=TRUE, test='Wald')
  } else {
    res %<>% DESeq2::DESeq(quiet=TRUE, test='LRT', reduced = ~ 1)
  }

  res %<>% DESeq2::results(contrast=c('group', target.level, ref.level), ...) %>% as.data.frame()

  res$padj[is.na(res$padj)] <- 1

  return(res)
}

estimateDEForTypeEdgeR <- function(cm, meta, design.formula) {
  design <- model.matrix(design.formula, meta)

  qlf <- edgeR::DGEList(cm, group = meta$group) %>%
    edgeR::calcNormFactors() %>%
    edgeR::estimateDisp(design = design) %>%
    edgeR::glmQLFit(design = design) %>%
    edgeR::glmQLFTest(coef=ncol(design))

  res <- qlf$table %>% .[order(.$PValue),] %>% set_colnames(c("log2FoldChange", "logCPM", "stat", "pvalue"))
  res$padj <- p.adjust(res$pvalue, method="BH")

  return(res)
}

estimateDEForTypeLimma <- function(cm, meta, design.formula, target.level) {
  mm <- model.matrix(design.formula, meta)
  fit <- limma::voom(cm, mm, plot = FALSE) %>% limma::lmFit(mm)

  contr <- paste0('group', target.level) %>% limma::makeContrasts(levels=colnames(coef(fit)))
  res <- limma::contrasts.fit(fit, contr) %>% limma::eBayes() %>% limma::topTable(sort.by="P", n=Inf) %>%
    set_colnames(c('log2FoldChange', 'AveExpr', 'stat', 'pvalue', 'padj', 'B'))

  return(res)
}

#' Summarize DE Resampling Results
#' @param var.to.sort Variable to calculate ranks
summarizeDEResamplingResults <- function(de.list, var.to.sort='pvalue') {
  de.res <- de.list[[1]]
  for (cell.type in names(de.res)) {
    genes.init <- genes.common <- rownames(de.res[[cell.type]]$res)
    mx.stat <- matrix(nrow = length(genes.common), ncol = 0, dimnames = list(genes.common,c()))
    for (i in 2:length(de.list)) {
      if (!(cell.type %in% names(de.list[[i]]))) next
      genes.common <- intersect(genes.common, rownames(de.list[[i]][[cell.type]]))
      mx.stat <- cbind(mx.stat[genes.common,,drop=FALSE],
                       de.list[[i]][[cell.type]][genes.common, var.to.sort,drop=FALSE])
    }

    if (ncol(mx.stat) == 0) {
      warning("Cell type ", cell.type, " was not present in any subsamples")
      next
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
