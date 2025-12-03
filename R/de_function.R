#' @import sccore
NULL

#' Validate parameters per cell type
#'
#' @param raw.mats List of raw count matrices
#' @param cell.groups Named clustering/annotation factor with cell names
#' @param model Model object returned by buildDesignMatrices
#' @keywords internal
validateDEPerCellTypeParams <- function(raw.mats, cell.groups, model) {
  checkPackageInstalled("DESeq2", bioc=TRUE)

  if (is.null(cell.groups)) stop('"cell.groups" must be specified')
  if (class(cell.groups) != "factor") stop('"cell.groups" must be a factor')

  ## raw.mats checks 
  if (is.null(raw.mats)) stop('"raw.mats" must be provided')
  if (!is.list(raw.mats))stop('"raw.mats" must be a list (e.g. list of count matrices)')
  
  if (is.null(names(raw.mats)) || any(names(raw.mats) == "")) {
    stop('"raw.mats" must be a *named* list; names should be sample IDs')
  }
  # model checks
  if (is.null(model)) stop('model must be provided')
  if (is.null(model$F) || is.null(model$contrast.F)) {
      stop("Model object must contain model matrix F and contrast.F")
  }
  if (is.null(rownames(model$F))) stop("model$F has no rownames (sample IDs)")
  missing.in.raw <- setdiff(rownames(model$F), names(raw.mats))
  if (length(missing.in.raw) > 0) {
    stop(
      "The following samples are in model$F but not in raw.mats: ",
      paste(missing.in.raw, collapse = ", ")
    )
  }
  # contrast checks
  if (!identical(names(model$contrast.F), colnames(model$F))) {
    stop('names(model$contrast.F) must exactly match colnames(model$F)')
  }
  if (all(abs(model$contrast.F) < 1e-12)) {
    stop("model$contrast.F is (numerically) all zeros; invalid contrast")
  }
  invisible(TRUE)
}

#' Subset matrices with common genes
#' @param cms List with count matrices
#' @param sample.groups (default=NULL)
#' @return list with count matrices with common genes
#' @keywords internal
subsetMatricesWithCommonGenes <- function(cms, sample.groups=NULL) {
  if (!is.null(sample.groups)) cms <- cms[unlist(sample.groups)]
  common.genes <- do.call(intersect, lapply(cms, colnames))
  cms %<>% lapply(function(m) m[, common.genes, drop=FALSE])
  return(cms)
}

#' Split strings and extract nth element
#' @description Function building on base::strsplit to extract nth element of a character string after splitting
#' 
#' @param x input character vector
#' @param split character vector containing a regular expression for splitting
#' @param n element to extract
#' @param fixed passed to strsplit
#' @keywords internal
strpart <- function(x, split, n, fixed = FALSE) {
  as.character(x) %>% strsplit(split, fixed=fixed) %>% sapply("[", n)
}

#' Add Z scores to DE results
#'
#' @param df Data.frame with the columns "pval", "padj" and "log2FoldChange"
#' @return Updated data.frame with Z scores
#' @examples 
#' \dontrun{
#' df_adj <- addZScores(df)
#' }
#' 
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

#' Prepare samples for DE analysis
#' @param sample.groups named list containing sample names
#' @param resampling.method one of "loo" (leave-one-out, remove one sample per iteration), "bootstrap", "fix.cells" (fixed number of cells per subsample), or "fix.samples" (fixed number of samples per iteration)
#' @param n.resamplings number of iterations (default=30)
#' @keywords internal
prepareSamplesForDE <- function(sample.groups, resampling.method=c('loo', 'bootstrap', 'fix.cells', 'fix.samples'),
                                n.resamplings=30) {
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

#' Differential expression using different methods (DESeq2, edgeR, wilcoxon, ttest) with various covariates
#'
#' @param raw.mats list of counts matrices; column for gene and row for cell
#' @param cell.groups factor specifying cell types (default=NULL)
#' @param model Model object returned by buildDesignMatrices (default=NULL)
#' @param common.genes boolean Only investigate common genes across cell groups (default=FALSE)
#' @param cooks.cutoff boolean cooksCutoff for DESeq2 (default=FALSE)
#' @param min.cell.count numeric Minimum cell count (default=10)
#' @param max.cell.count numeric Maximum cell count (default=Inf). If Inf, there is no limit set.
#' @param fix.n.samples Number of samples to fix (default=NULL). If greater the the length of the s.groups, an error is thrown.
#' @param verbose boolean Whether to output verbose messages (default=TRUE)
#' @param independent.filtering boolean independentFiltering for DESeq2 (default=FALSE)
#' @param n.cores numeric Number of cores (default=1)
#' @param return.matrix Return merged matrix of results (default=TRUE)
#' @param sample.meta dataframe with sample metadata (default=NULL)
#' @param formula design formula (default=NULL)
#' @param contrast character vector of length 3 specifying the comparison. e.g c("group", "control", "treatment")
#' @param test DE method: DESeq2, edgeR, wilcoxon, ttest
#' @param gene.filter matrix/boolean Genes to omit (rows) per cluster (cols) (default=NULL)
#' @return differential expression for each cell type
#'
#' @export
estimateDEPerCellTypeInner <- function(raw.mats, cell.groups=NULL, s.groups=NULL, sample.meta=NULL, model=NULL, 
                                       common.genes=FALSE, cooks.cutoff=FALSE, min.cell.count=10, max.cell.count=Inf,
                                       independent.filtering=TRUE, n.cores=4, return.matrix=TRUE, fix.n.samples=NULL,
                                       verbose=TRUE, test='Wald', gene.filter=NULL) {
  # Validate input
  validateDEPerCellTypeParams(raw.mats, cell.groups, model)
  tmp <- tolower(strsplit(test, split='\\.')[[1]])
  test <- tmp[1]
  test.type <- ifelse(is.na(tmp[2]), '', tmp[2])

  # Filter data and convert to the right format
  if (verbose) message("Preparing matrices for DE")
  if (common.genes) {
    raw.mats %<>% subsetMatricesWithCommonGenes(s.groups)
  } else {
    gene.union <- lapply(raw.mats, colnames) %>% Reduce(union, .)
    raw.mats %<>% lapply(sccore::extendMatrix, gene.union)
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
    cm <- cm.bulk.per.type[[l]]
    if (!is.null(gene.filter)) {
      gene.to.remain <- gene.filter %>% {rownames(.)[.[,l]]} %>% intersect(rownames(cm))
      cm <- cm[gene.to.remain,,drop=FALSE]
    }

    cur.s.groups <- lapply(s.groups, intersect, colnames(cm))
    if (!is.null(fix.n.samples)) {
      if (min(sapply(s.groups, length)) < fix.n.samples) {
        warning("The cluster does not have enough samples")
        return(NULL)
      }
      cur.s.groups %<>% lapply(sample, fix.n.samples)
      cm <- cm[, unlist(cur.s.groups), drop=FALSE]
    }

    # subset/reorder design rows to match cm columns
  model.ct <- model
  model.ct$F <- model.ct$F[match(colnames(cm), rownames(model.ct$F)), , drop = FALSE]
  if (!identical(rownames(model.ct$F), colnames(cm))) {
    warning("Failed to align model$F rows to cm columns for cell type ", l)
    return(NULL)
  }

  # DE Testing 
  if (verbose) message("Running DE for cell type: ", l)
  res <- tryCatch({
    if (test %in% c('wilcoxon', 't-test')) {
    ## Derive simple triplet from model$contrast_spec
    triplet <- tryCatch(
      extractSimpleTripletFromSpec(model.ct),
      error = function(e) {
      warning("Cannot use Wilcoxon/t-test for cell type ", l,
          ": ", conditionMessage(e))
      return(NULL)
      }
    )
    if (is.null(triplet)) return(NULL)

    group.var <- triplet[1]
    num.level <- triplet[2]
    den.level <- triplet[3]

    ## Check that the model formula is "simple enough" (no covariates)
    covars <- getCovariatesForTriplet(model, group.var)
    if (length(covars) > 0L) {
      stop(
      "Wilcoxon/t-test are unadjusted two-group tests, ",
      "but your model formula includes additional terms: ",
      paste(covars, collapse = ", "), ".\n",
      "Either:\n",
      "  - use a model-based method (DESeq2 / edgeR / limma) with this formula, or\n",
      "  - provide a simpler design formula (e.g. ~ ", group.var,
      ") if you explicitly want unadjusted Wilcoxon/t-tests."
      )
    }

    ## Proceed if the design is effectively ~ group.var only
    if (!group.var %in% colnames(sample.meta)) {
      warning("Grouping variable '", group.var,
          "' from contrast_spec not found in metadata for cell type ", l)
      return(NULL)
    }
    
    sample.meta <- sample.meta[colnames(cm), , drop = FALSE]
    sample.meta[[group.var]] <- factor(sample.meta[[group.var]])

    # keep only the two contrast levels
    keep <- sample.meta[[group.var]] %in% c(num.level, den.level)
    if (sum(keep) < 2L) {
      warning("Fewer than two samples in contrast levels for cell type ", l)
      return(NULL)
    }

    cm2   <- cm[, keep, drop = FALSE]
    meta2 <- sample.meta[keep, , drop = FALSE]

    # drop unused and relevel reference
    meta2[[group.var]] <- droplevels(meta2[[group.var]])
    if (!all(c(num.level, den.level) %in% levels(meta2[[group.var]]))) {
      warning("Contrast levels not found after subsetting for cell type ", l)
      return(NULL)
    }
    meta2[[group.var]] <- stats::relevel(meta2[[group.var]], ref = den.level)

    # require at least 2 samples per group
    tab <- table(meta2[[group.var]])
    if (any(tab < 2L)) {
      warning("Each group must be present in at least two samples (Wilcoxon/t-test) — skipping cell type: ", l)
      return(NULL)
    }

    # simple design for normalization: ~ group.var
    design.formula <- stats::reformulate(group.var)

    cm.norm <- normalizePseudoBulkMatrix(cm2, meta = meta2, design.formula = design.formula, type = test.type)
    estimateDEForTypePairwiseStat(cm.norm, meta = meta2, group.var = group.var, target.level = num.level, test = test)
    } else if (test == 'deseq2') {
    estimateDEForTypeDESeq(cm, sample.meta, model.ct, test.type=test.type)
    } else if (test == 'edger') {
    estimateDEForTypeEdgeR(cm, model.ct)
    } else if (test == 'limma-voom') {
    estimateDEForTypeLimma(cm, model.ct)
    } else {
    stop("Unknown test: ", test)
    }
  }, error = function(e) {
    warning("DE failed for cell type ", l, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(res)) return(NULL)
  res$Gene <- rownames(res)

  if (!is.null(res) && !is.na(res[[1]][1])) {
    res <- addZScores(res) %>% .[order(.$pvalue, decreasing = FALSE), ]
  }

  return(res)
  }, n.cores = n.cores, progress = verbose, mc.preschedule = TRUE, mc.allow.recursive = TRUE) %>%
  .[!sapply(., is.null)]

  if (verbose) {
    dif <- setdiff(levels(cell.groups), names(de.res))
    if (length(dif) > 0) {
      message("DEs not calculated for ", length(dif), " cell group(s): ", paste(dif, collapse=', '))
    }
  }

  return(de.res)
}

#' Filter DE metadata
#' @param meta data frame containing metadata in columns
#' @return cleaned data frame with metadata
#' @keywords internal
filterDEMetadata <- function(meta) {
  # Remove unique columns
  i.m.remain <- c()
  for(i.m in 1:ncol(meta)){
    if(length(unique(meta[, i.m])) != 1) i.m.remain <- c(i.m.remain, i.m)
  }
  meta <- meta[, i.m.remain, drop=F]
  if(ncol(meta) == 1) return(meta)
  # The same columns
  for(i in 2:ncol(meta)){
    for(j in i:ncol(meta)){
      if(i == j) next
      if((nrow(unique(meta[,c(i, j)])) == length(unique(meta[,i]))))
        return(NULL)
    }
  }
  return(meta)
}

#' Normalize pseudo-bulk matrix
#' @param cm count matrix
#' @param meta (default=NULL)
#' @param design.formula (default=NULL)
#' @param type (default="totcount")
#' @return normalized count matrix
#' @keywords internal
# Requires availability test for DESeq2 or edgeR
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

#' Estimate pair-wise DEGs
#' @param cm.norm normalized count matrix
#' @param meta data frame with meta data
#' @param group.var grouping variable
#' @param target.level target level, e.g., disease group
#' @param test type of test, either "wilcoxon" or "t-test"
#' @return data frame containing DEGs using a pair-wise test
#' @keywords internal
# Requires availability test for "scran"
estimateDEForTypePairwiseStat <- function(cm.norm, meta, group.var, target.level, test) {
  if (!group.var %in% colnames(meta)) {
    stop("Grouping variable '", group.var, "' not found in meta.")
  }
  groups <- droplevels(factor(meta[[group.var]]))
  if (test == 'wilcoxon') {
    res <- scran::pairwiseWilcox(cm.norm, groups = groups)$statistics[[1]] %>%
      data.frame() %>% setNames(c("AUC", "pvalue", "padj"))
    res$log2FoldChange <- log2(cm.norm + 1) %>% apply(1, function(x) {
    mean(x[groups == target.level]) - mean(x[groups != target.level])})
  } else if (test == 't-test') {
    cm.norm <- log2(cm.norm * 1e6 + 1)
    res <- scran::pairwiseTTests(cm.norm, groups = groups)$statistics[[1]] %>%
      data.frame() %>% setNames(c("AUC", "pvalue", "padj"))
    res$log2FoldChange <- cm.norm %>% apply(1, function(x) {
    mean(x[groups == target.level]) - mean(x[groups != target.level])})
  }
  return(res)
}

#' Estimate DE using DESeq2
#' 
#' @param cm count matrix
#' @param sample.meta sample metadata data frame
#' @param model list returned by buildDesignMatrices (must contain F and contrast.F)
#' @param test.type test type incorporated in DESeq2, either "wald" or "LRT"
#' @param ... additional parameters forwarded to DESeq2
#' 
#' @keywords internal
#' Requires "ashr" package for lfcShrink with ashr
estimateDEForTypeDESeq <- function(cm, sample.meta, model, test.type='wald', cooks.cutoff=FALSE, independent.filtering=TRUE, ...) {
  dds <- DESeq2::DESeqDataSetFromMatrix(cm, sample.meta[colnames(cm), ], design = ~ 1)
  DESeq2::design(dds) <- model$F
  if (test.type == 'wald') {
      dds <- DESeq2::DESeq(dds, quiet = TRUE, test = 'Wald')
      cF <- model$contrast.F

      # Make sure order matches what DESeq2 expects
      n <- DESeq2::resultsNames(dds)
      n[1] <- "(Intercept)"
      cF <- cF[n]

      res0 <- DESeq2::results(dds, contrast = cF)
      res <- DESeq2::lfcShrink(dds, contrast=cF, res = res0, type = "ashr") # TODO: add 'ashr' package to dependencies
    } else { # lrt
        reduced <- model$Z
        if (is.null(reduced)) {
             n <- nrow(model$F)
             reduced <- matrix(1, nrow = n, ncol = 1)
             colnames(reduced) <- "Intercept"
             rownames(reduced) <- rownames(model$F)
         } else {
           reduced <- reduced[rownames(model$F), , drop = FALSE]
         }
        dds <- DESeq2::DESeq(dds, quiet = TRUE, test = 'LRT', full = model$F, reduced = reduced)
         n <- DESeq2::resultsNames(dds)
         n[1] <- "(Intercept)"
         cF <- cF[n]

        res <- DESeq2::results(dds, contrast = cF, cooksCutoff = cooks.cutoff,
                               independentFiltering = independent.filtering)  # No lfcShrink for LRT
    }
  res <- as.data.frame(res) 
  res$padj[is.na(res$padj)] <- 1

  return(res)
}

#' Estimate DE using edgeR (design + contrast from buildDesignMatrices)
#'
#' @param cm count matrix (genes x samples)
#' @param model list returned by buildDesignMatrices (must contain F and contrast.F)
#'
#' @keywords internal
estimateDEForTypeEdgeR <- function(cm, model) {
    # Extract design and contrast
    if (is.null(model$F) || is.null(model$contrast.F)) {
        stop("model must contain components 'F' (design) and 'contrast.F' (contrast vector).")
    }
    design   <- model$F
    contrast <- model$contrast.F
    # Align design rows to columns of cm (samples)
    if (is.null(rownames(design))) {
        stop("design matrix (model$F) must have rownames corresponding to sample IDs.")
    }
    if (is.null(colnames(cm))) {
        stop("count matrix 'cm' must have column names corresponding to sample IDs.")
    }
    # Reorder / subset design to match cm columns
    design <- design[match(colnames(cm), rownames(design)), , drop = FALSE]
    if (!identical(rownames(design), colnames(cm))) {
        stop("Row names of model$F (samples) must match column names of cm (samples).")
    }
    # Sanity check for contrast
    if (!identical(names(contrast), colnames(design))) {
        stop("names(model$contrast.F) must exactly match colnames(model$F).")
    }
    # edgeR pipeline
    dge <- edgeR::DGEList(counts = cm)
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge, design = design)
    fit <- edgeR::glmQLFit(dge, design = design)
    
    # Test the supplied contrast
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
    
    # --- Compute per-gene stats (already in qlf$table)
    tab <- qlf$table
    tab$FDR <- p.adjust(tab$PValue, method = "BH")
    
    # --- Build an edgeR-style object
    #     topTags() only requires: table, comparison, genes (optional)
    result <- list(
        table       = tab,
        comparison  = contrast     # printed by topTags()
       # genes       = data.frame(Gene = rownames(tab), row.names = rownames(tab))
    )
    class(result) <- c("DGELRT", "DGEExact")  # important for topTags()
    return(result)
}

#' Estimate DE using limma-voom (design + contrast from buildDesignMatrices)
#'
#' @param cm count matrix (genes x samples)
#' @param model list returned by buildDesignMatrices (must contain F and contrast.F)
#'
#' @keywords internal
estimateDEForTypeLimma <- function(cm, model) {
    # Extract design and contrast
    if (is.null(model$F) || is.null(model$contrast.F)) {
        stop("model must contain components 'F' (design) and 'contrast.F' (contrast vector).")
    }
    design   <- model$F
    contrast <- model$contrast.F
    
    # Align design rows to columns of cm (samples)
    if (is.null(rownames(design))) {
        stop("design matrix (model$F) must have rownames corresponding to sample IDs.")
    }
    if (is.null(colnames(cm))) {
        stop("count matrix 'cm' must have column names corresponding to sample IDs.")
    }
    
    design <- design[match(colnames(cm), rownames(design)), , drop = FALSE]
    
    if (!identical(rownames(design), colnames(cm))) {
        stop("Row names of model$F (samples) must match column names of cm (samples).")
    }
    # Sanity check for contrast
    if (!identical(names(contrast), colnames(design))) {
        stop("names(model$contrast.F) must exactly match colnames(model$F).")
    }
    
    # limma-voom pipeline
    dge <- edgeR::DGEList(counts = cm)
    dge <- edgeR::calcNormFactors(dge)
    v   <- limma::voom(dge, design = design, plot = FALSE)
    
    fit <- limma::lmFit(v, design = design)
  
    # Build 1-column contrast matrix from model$contrast.F
    Cmat <- matrix(contrast,
                   ncol = 1,
                   dimnames = list(colnames(design), "contrast1"))
    fit2 <- limma::contrasts.fit(fit, Cmat) %>%
        limma::eBayes()
    res <- limma::topTable(fit2,
                           coef    = "contrast1",
                           sort.by = "P",
                           n       = Inf) %>%
        setNames(c("log2FoldChange", "AveExpr", "stat", "pvalue", "padj", "B"))
    
    return(res)
}

#' Summarize DE Resampling Results
#'
#' @param de.list list with DE results. Data frame with DE results are found in the first element
#' @param var.to.sort Variable to calculate ranks (default="pvalue")
#' 
#' @keywords internal
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

#' Append statistics to DE results
#' 
#' @param de.list list with DE results
#' @param expr.frac.per.type 
#' 
#' @keywords internal
appendStatisticsToDE <- function(de.list, expr.frac.per.type) {
  for (n in names(de.list)) {
    de.list[[n]]$res %<>%
      mutate(CellFrac=expr.frac.per.type[Gene, n], SampleFrac=Matrix::rowMeans(de.list[[n]]$cm > 0)[Gene]) %>%
      as.data.frame(stringsAsFactors=FALSE) %>% set_rownames(.$Gene)
  }

  return(de.list)
}

#' get expression fraction per cell group
#' 
#' @param cm sparse count matrix 
#' @param cell.groups factor containing cell groups with cell names as names
#' 
#' @keywords internal
getExpressionFractionPerGroup <- function(cm, cell.groups) {
  cm@x <- as.numeric(cm@x > 1e-10)
  fracs <- collapseCellsByType(cm, cell.groups, min.cell.count=0) %>%
    {. / as.vector(table(cell.groups)[rownames(.)])} %>% Matrix::t()
  return(fracs)
}

#' filter genes by expression per celltype not globally
#' 
#' @param counts sparse count matrix 
#' @param cell.groups factor containing cell groups with cell names as names
#' @param threshold minimum fraction of expressing cells per cell type (default=0.05)
#
#' @keywords internal
getPerCellTypeGeneFilter <- function(counts, cell.groups, threshold = 0.05) {
  counts@x <- as.numeric(counts@x > 1)
  cell.types <- unique(cell.groups)
  cell.type.indicator <- Matrix(0, nrow = length(cell.groups), ncol = length(cell.types),
                               dimnames = list(names(cell.groups), cell.types))
  row.inds <- match(names(cell.groups), rownames(cell.type.indicator))
  col.inds <- match(cell.groups, colnames(cell.type.indicator))
  cell.type.indicator[cbind(row.inds, col.inds)] <- 1
  
  counts.per.type <- counts %*% cell.type.indicator
  n.cells <- colSums(cell.type.indicator)
  expr.frac <- sweep(counts.per.type, 2, n.cells, FUN = "/")
  filters <- expr.frac > threshold
  
  filters.list <- lapply(seq_len(ncol(filters)), function(i) {
    setNames(as.vector(filters[, i]), rownames(counts))
  })
  names(filters.list) <- colnames(filters)
  return(filters.list)
}

#' Extract simple triplet contrast from model contrast specification
#' @param model model object returned by buildDesignMatrices
#' @return character vector of length 3 with c("group", "num", "den")
#' @keywords internal
extractSimpleTripletFromSpec <- function(model) {
    spec <- model$contrast_spec
    if (is.null(spec)) {
        stop("model$contrast_spec is NULL; cannot derive a simple triplet contrast.")
    }
    if (is.character(spec) && length(spec) == 3L) {
        # DESeq2-style triple: c("group", "B", "A")
        return(spec)
    }
    if (is.list(spec) &&
        identical(spec$type, "simple") &&
        !grepl(":", spec$term, fixed = TRUE)) {
        # simple main-effect contrast:
        # list(type="simple", term="group", num="B", den="A")
        return(c(spec$term, spec$num, spec$den))
    }
    
    stop("model$contrast_spec is not representable as a simple triplet (main-effect) contrast.")
}

#' Get covariates from model formula for triplet tests
#' @param model model object returned by buildDesignMatrices
#' @param group.var grouping variable in
#' @return character vector with covariate names
#' @keywords internal
getCovariatesForTriplet <- function(model, group.var) {
    # prefer formula_used if present, else fall back to whatever you pass in
    if (!is.null(model$formula_used)) {
        f <- model$formula_used
    } else {
        stop("model$formula_used is NULL; cannot inspect covariates for triplet tests.")
    }
    
    tt <- terms(f)
    # RHS term labels, e.g. c("group", "batch", "age", "group:batch")
    term.labels <- attr(tt, "term.labels")
    if (is.null(term.labels)) term.labels <- character(0)
    
    # treat any term that involves ':' as an interaction / covariate in this context
    # we only allow a pure main effect on group.var
    # So: anything not equal to group.var is a "covariate" for Wilcoxon/t
    covars <- setdiff(term.labels, group.var)
    covars
}

