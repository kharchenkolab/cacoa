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

#' Subset matrices with common genes (by sample vector)
#' @param cms named list of count matrices (cells x genes)
#' @param samples character vector of sample IDs to keep (optional)
#' @keywords internal
subsetMatricesWithCommonGenesSamples <- function(cms, samples = NULL) {
  if (!is.null(samples)) cms <- cms[samples]
  common.genes <- Reduce(intersect, lapply(cms, colnames))
  lapply(cms, function(m) m[, common.genes, drop = FALSE])
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

#' Prepare sample sets for DE resampling (model-based)
#' - loo: drop ONE contrast-side sample (core) at a time, keep all other samples (incl non-core levels/nuisance)
#' - fix.*: return all samples; variation happens inside DE inner (fix.samples / fix.cells)
#'
#' @param model list returned by buildDesignMatrices (must contain F and contrast_spec)
#' @param sample.meta data.frame with rownames = sample IDs
#' @param resampling.method "loo", "fix.cells", "fix.samples"
#' @param n.resamplings number of iterations (used for fix.*)
#' @param min.core.per.block minimum CORE (contrast-side) samples per block to allow dropping from that block
#' @param min.core.per.side minimum CORE samples per contrast side AFTER dropping (safety; default 2)
#'
#' @keywords internal
prepareSamplesForDE_model <- function(model, sample.meta, resampling.method = c("loo", "fix.cells", "fix.samples"),
                                      n.resamplings = 30, min.core.per.block = 3, min.core.per.side = 2) {
  resampling.method <- match.arg(resampling.method)

  all.samples <- rownames(model$F)
  if (is.null(all.samples) || length(all.samples) == 0L) {
    stop("model$F must have non-empty rownames (sample IDs)")
  }
  if (is.null(sample.meta) || is.null(rownames(sample.meta))) {
    stop("sample.meta must be a data.frame with rownames = sample IDs")
  }

  # fix.cells / fix.samples: sample set unchanged; the variation is induced later
  if (resampling.method != "loo") {
    out <- lapply(seq_len(n.resamplings), function(i) all.samples)
    names(out) <- paste0("fix.", seq_len(n.resamplings))
    return(out)
  }

  # LOO - core only
  triplet <- tryCatch(extractSimpleTripletFromSpec(model), error = function(e) NULL)
  if (is.null(triplet)) {
    stop("loo requires a simple triplet contrast_spec (e.g. c('group','B','A')).")
  }
  group.var <- triplet[1]; num.level <- triplet[2]; den.level <- triplet[3]

  sm <- sample.meta[all.samples, , drop = FALSE]
  if (!group.var %in% colnames(sm)) {
    stop("Grouping var '", group.var, "' not found in sample.meta")
  }

  sm[[group.var]] <- factor(sm[[group.var]])
  core.mask <- sm[[group.var]] %in% c(num.level, den.level)

  core.samples <- rownames(sm)[core.mask]
  if (length(core.samples) == 0L) {
    stop("No core samples found for levels {", den.level, ", ", num.level, "}")
  }
  # ensure there are enough samples per side to do LOO at all
  n.den <- sum(sm[[group.var]] == den.level, na.rm = TRUE)
  n.num <- sum(sm[[group.var]] == num.level, na.rm = TRUE)
  if (n.den <= min.core.per.side || n.num <= min.core.per.side) {
    stop("Not enough core samples per contrast side for loo (need > ", min.core.per.side,
         " per side; have ", n.den, " and ", n.num, ").")
  }

  drop.candidates <- core.samples

  # Optional block-aware: only drop from blocks that keep enough CORE samples after drop
  if (!is.null(model$blocks)) {
    b <- model$blocks
    names(b) <- rownames(model$F)

    # blocks that have enough total core to be "drop-eligible"
    core.counts <- tapply(core.mask, b, sum)
    ok.blocks <- names(core.counts)[core.counts >= min.core.per.block]

    drop.candidates <- names(b)[(b %in% ok.blocks) & core.mask]
    if (length(drop.candidates) == 0L) {
      warning("No blocks satisfy min.core.per.block; falling back to dropping from all core samples.")
      drop.candidates <- core.samples
    }
  }

  # Filter drops that would violate min.core.per.side (per-side after drop)
  is.valid.drop <- vapply(drop.candidates, function(s) {
    g <- as.character(sm[s, group.var])
    if (is.na(g)) return(FALSE)
    if (g == den.level) return((n.den - 1) >= min.core.per.side)
    if (g == num.level) return((n.num - 1) >= min.core.per.side)
    FALSE
  }, logical(1))

  drop.candidates <- drop.candidates[is.valid.drop]
  if (length(drop.candidates) == 0L) {
    warning("No core samples can be dropped while keeping >= ", min.core.per.side,
            " per contrast side. Falling back to dropping from all core samples.")
    drop.candidates <- core.samples
  }

  out <- lapply(drop.candidates, function(s) setdiff(all.samples, s))
  names(out) <- paste0("loo.core.", drop.candidates)
  out
}


#' Prepare samples for DE analysis (old; two-group)
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

#' Differential expression per cell type (model-matrix aware; NO two-group lists)
#'
#' NOTE: This replaces the existing estimateDEPerCellTypeInner() that uses s.groups.
#' The new argument is `samples` (character vector of sample IDs).
#'
#' @export
estimateDEPerCellTypeInner_model <- function(raw.mats, cell.groups = NULL, samples = NULL, sample.meta = NULL, model = NULL,
                                             sample.per.cell = NULL, n.cells.subsample.core = NULL, seed = NULL,
                                             common.genes = FALSE, cooks.cutoff = FALSE, min.cell.count = 10,
                                             max.cell.count = Inf, independent.filtering = TRUE, n.cores = 1,
                                             return.matrix = TRUE, fix.n.samples = NULL, max.resample.tries = 50,
                                             verbose = TRUE, test = "Wald", gene.filter = NULL) {
    validateDEPerCellTypeParams(raw.mats, cell.groups, model)
    
    tmp <- tolower(strsplit(test, split = "\\.")[[1]])
    test <- tmp[1]
    test.type <- ifelse(is.na(tmp[2]), "", tmp[2])
    
    if (is.null(samples)) samples <- rownames(model$F)
    samples <- intersect(samples, names(raw.mats))
    if (length(samples) == 0L) stop("No overlap between requested samples and raw.mats names")
    
    if (verbose) message("Preparing matrices for DE")
    raw.mats.sub <- raw.mats[samples]
    if (!is.null(n.cells.subsample.core)) {
     if (is.null(sample.per.cell)) stop("fix.cells requires sample.per.cell")
     core.samples <- getCoreSamples(model, samples)

    raw.mats.sub <- subsampleCellsPerTypePerSample(raw.mats.sub = raw.mats.sub, cell.groups = cell.groups,
                                                   sample.per.cell = sample.per.cell, core.samples = core.samples,
                                                   n.cells.subsample = n.cells.subsample.core,
                                                   min.cell.count = min.cell.count, seed = seed)
  }

  if (common.genes) {
    raw.mats.sub <- subsetMatricesWithCommonGenesSamples(raw.mats.sub)
  } else {
    gene.union <- Reduce(union, lapply(raw.mats.sub, colnames))
    raw.mats.sub <- lapply(raw.mats.sub, sccore::extendMatrix, gene.union)
  }

  cm.bulk.per.samp <- raw.mats.sub %>%
    lapply(collapseCellsByType, groups = cell.groups,
           min.cell.count = min.cell.count,
           max.cell.count = max.cell.count) %>%
    .[sapply(., nrow) > 0]

  passed.samples <- names(cm.bulk.per.samp)
  if (length(passed.samples) == 0L) return(list())

  cm.bulk.per.type <- levels(cell.groups) %>% sn() %>%
    lapply(function(cg) {
      tcms <- cm.bulk.per.samp %>%
        lapply(function(cm) {
          if (cg %in% rownames(cm)) cm[cg, , drop = FALSE] else NULL
        }) %>%
        .[!sapply(., is.null)]
      if (length(tcms) == 0) return(NULL)

      tcms %>%
        { set_rownames(do.call(rbind, .), names(.)) } %>%
        `mode<-`("integer") %>%
        .[, colSums(.) > 0, drop = FALSE]
    }) %>%
        .[sapply(., length) > 0] %>%
        lapply(t)
    
    if (verbose) message("Estimating DE per cell type")
    
    # DE celltype loop
    de.res <- names(cm.bulk.per.type) %>% sn() %>% plapply(function(l) {
        
        cm <- cm.bulk.per.type[[l]]
        
        if (!is.null(gene.filter)) {
            gene.to.remain <- rownames(gene.filter)[as.logical(gene.filter[, l])]
            gene.to.remain <- intersect(gene.to.remain, rownames(cm))
            cm <- cm[gene.to.remain, , drop = FALSE]
        }
        
        samp.ct <- colnames(cm)
        if (length(samp.ct) == 0L) return(NULL)
        
        if (is.null(sample.meta) || is.null(rownames(sample.meta))) {
            warning("sample.meta with rownames(sample IDs) required; skipping cell type ", l)
            return(NULL)
        }
        
        # eligible = samples where this cell type exists
        meta.ct <- sample.meta[samp.ct, , drop = FALSE]
        
        model.ct <- model

      if (!is.null(fix.n.samples)) {
        triplet <- tryCatch(
          extractSimpleTripletFromSpec(model),
          error = function(e) NULL)
        if (is.null(triplet)) {
          warning("fix.samples requested but contrast_spec is not a ",
                  "simple triplet; skipping cell type ", l)
          return(NULL)
        }

        group.var <- triplet[1]
        num.level <- triplet[2]
        den.level <- triplet[3]

        if (is.null(sample.meta) || is.null(rownames(sample.meta))) {
          warning("sample.meta with rownames(sample IDs) required; ",
                  "skipping cell type ", l)
          return(NULL)
        }

        meta.ct <- sample.meta[samp.ct, , drop = FALSE]
        if (!group.var %in% colnames(meta.ct)) {
          warning("Grouping var '", group.var,
                  "' not found in sample.meta; ",
                  "skipping cell type ", l)
          return(NULL)
        }

        meta.ct[[group.var]] <- factor(meta.ct[[group.var]])

        is.core <- meta.ct[[group.var]] %in% c(den.level, num.level)
        core.ct <- rownames(meta.ct)[is.core]
        nuis.ct <- rownames(meta.ct)[!is.core]

        pool.den <- rownames(meta.ct)[
          meta.ct[[group.var]] == den.level]
        pool.num <- rownames(meta.ct)[
          meta.ct[[group.var]] == num.level]

        if (length(pool.den) < fix.n.samples ||
            length(pool.num) < fix.n.samples) {
          warning("Not enough CORE samples per contrast side for ",
                  "fix.samples in cell type ", l, " (need ",
                  fix.n.samples, " per side; have ", length(pool.den),
                  " and ", length(pool.num), ").")
          return(NULL)
        }

        found <- FALSE
        keep.core <- NULL
        model.keep <- NULL

        for (try.i in seq_len(max.resample.tries)) {
          keep.den <- sample(pool.den, fix.n.samples)
          keep.num <- sample(pool.num, fix.n.samples)

          keep.core.try <- unique(c(keep.den, keep.num))
          keep.all <- unique(c(keep.core.try, nuis.ct))

          cm.try <- cm[, keep.all, drop = FALSE]

        # subset design rows to match kept samples
        F.sub <- model$F[match(keep.all, rownames(model$F)), , drop = FALSE]
        if (!identical(rownames(F.sub), keep.all)) next

          rep <- tryCatch(
            repairDesignAfterRowSubset(F.sub, model$contrast.F),
            error = function(e) NULL)
          if (is.null(rep)) next

          if (!isContrastNonzero(rep$contrast)) next

          model.try <- model
          model.try$F <- rep$F
          model.try$contrast.F <- rep$contrast

          found <- TRUE
          keep.core <- keep.core.try
          model.keep <- model.try

          cm <- cm.try
          samp.ct <- colnames(cm)
          model.ct <- model.keep
          meta.ct <- sample.meta[samp.ct, , drop = FALSE]
          break
        }

        if (!found) {
          warning("fix.samples: could not find a usable ",
                  "(contrast-estimable) repaired design for cell type ",
                  l, " after ", max.resample.tries, " tries; skipping.")
          return(NULL)
        }
      } else {
        F.sub <- model$F[match(samp.ct, rownames(model$F)), ,
                         drop = FALSE]
        if (!identical(rownames(F.sub), samp.ct)) {
          warning("Failed to align model$F rows to cm columns for ",
                  "cell type ", l)
          return(NULL)
        }

        rep <- tryCatch(
          repairDesignAfterRowSubset(F.sub, model$contrast.F),
          error = function(e) NULL)
        if (is.null(rep) || !isContrastNonzero(rep$contrast)) {
          warning("Contrast not estimable after design repair for ",
                  "cell type ", l)
          return(NULL)
        }

        model.ct$F <- rep$F
        model.ct$contrast.F <- rep$contrast
        meta.ct <- if (!is.null(sample.meta)) {
          sample.meta[samp.ct, , drop = FALSE]
        } else {
          NULL
        }
      }

      sample.meta.ct <- sample.meta[rownames(model.ct$F), ,
                                    drop = FALSE]

      if (verbose) message("Running DE for cell type: ", l)

      res <- tryCatch({
        if (test %in% c("wilcoxon", "t-test")) {
          triplet <- extractSimpleTripletFromSpec(model.ct)
          group.var <- triplet[1]
          num.level <- triplet[2]
          den.level <- triplet[3]

          covars <- getCovariatesForTriplet(model.ct, group.var)
          if (length(covars) > 0L) {
            stop("Wilcoxon/t-test require no covariates; found: ",
                 paste(covars, collapse = ", "))
          }
          if (!group.var %in% colnames(sample.meta.ct)) {
            stop("sample.meta must contain grouping variable '",
                 group.var, "' for Wilcoxon/t-test")
          }

          meta2 <- sample.meta.ct
          meta2[[group.var]] <- factor(meta2[[group.var]])
          keep <- meta2[[group.var]] %in% c(num.level, den.level)
          cm2 <- cm[, keep, drop = FALSE]
          meta2 <- meta2[keep, , drop = FALSE]

          meta2[[group.var]] <- droplevels(meta2[[group.var]])
          meta2[[group.var]] <- stats::relevel(meta2[[group.var]],
                                               ref = den.level)

          tab <- table(meta2[[group.var]])
          if (any(tab < 2L)) {
            stop("Need >=2 samples per group for Wilcoxon/t-test")
          }

          design.formula <- stats::reformulate(group.var)
          cm.norm <- normalizePseudoBulkMatrix(
            cm2, meta = meta2, design.formula = design.formula,
            type = test.type)
          estimateDEForTypePairwiseStat(
            cm.norm, meta = meta2, group.var = group.var,
            target.level = num.level, test = test)

        } else if (test == "deseq2") {
          estimateDEForTypeDESeq(
            cm, sample.meta.ct, model.ct, test.type = test.type,
            cooksCutoff = cooks.cutoff,
            independentFiltering = independent.filtering)

        } else if (test == "edger") {
          estimateDEForTypeEdgeR(cm, model.ct)

        } else if (test == "limma-voom") {
          estimateDEForTypeLimma(cm, model.ct)

        } else {
          stop("Unknown test: ", test)
        }
      }, error = function(e) {
        warning("DE failed for cell type ", l, ": ",
                conditionMessage(e))
        NULL
      })

      res <- coerceDEResultToDF(res)
      if (is.null(res) || nrow(res) == 0L) return(NULL)

      if (is.null(rownames(res)) || anyNA(rownames(res)) ||
          any(rownames(res) == "")) {
        gene.col <- intersect(c("Gene", "gene", "genes"),
                              colnames(res))
        if (length(gene.col) > 0) {
          rownames(res) <- as.character(res[[gene.col[1]]])
        }
      }

      res$Gene <- rownames(res)
      res <- standardizeDEColumns(res)

      if ("pvalue" %in% names(res) &&
          "log2FoldChange" %in% names(res) &&
          is.numeric(res$log2FoldChange)) {
        ok <- !is.na(res$pvalue)
        if (any(ok)) {
          res <- addZScores(res)
          res <- res[order(res$pvalue, decreasing = FALSE,
                           na.last = TRUE), , drop = FALSE]
        }
      } else {
        if ("pvalue" %in% names(res)) {
          res <- res[order(res$pvalue, decreasing = FALSE,
                           na.last = TRUE), , drop = FALSE]
        }
      }

      attr(res, "core.samples.used") <-
        attr(sample.meta.ct, "fix.samples.core")

      if (return.matrix) {
        return(list(res = res, cm = cm, meta = sample.meta.ct))
      }
      res
    }, n.cores = n.cores, progress = verbose, mc.preschedule = TRUE,
    mc.allow.recursive = TRUE) %>%
    .[!sapply(., is.null)]

  if (verbose) {
    dif <- setdiff(levels(cell.groups), names(de.res))
    if (length(dif) > 0) {
      message("DEs not calculated for ", length(dif),
              " cell group(s): ", paste(dif, collapse = ", "))
    }
  }

  de.res
}

# keep-or-drop columns that break rank after row subsetting
# returns list(F = ..., contrast = ...)
#' @keywords internal
repairDesignAfterRowSubset <- function(F, contrast, tol = 1e-12) {
    
    if (is.null(F) || nrow(F) == 0L) stop("F is empty")
    if (is.null(contrast)) stop("contrast is NULL")
    
    # align contrast to design columns
    contrast <- contrast[colnames(F)]
    contrast[is.na(contrast)] <- 0
    
    # drop columns that are all-zero or constant
    is.zero <- apply(F, 2, function(x) all(abs(x) < tol))
    is.const <- apply(F, 2, function(x) max(x) - min(x) < tol)
    
    keep <- !(is.zero | is.const)
    # keep intercept even if constant? Usually intercept is constant by definition.
    # If you have (Intercept), do NOT drop it just because it's constant.
    if ("(Intercept)" %in% colnames(F)) keep[colnames(F) == "(Intercept)"] <- TRUE
    
    F <- F[, keep, drop = FALSE]
    contrast <- contrast[colnames(F)]
    
    if (ncol(F) == 0L) stop("All columns dropped from design after repair")
    
    # ensure full rank by dropping dependent columns using QR pivoting
    qrF <- qr(F)
    rnk <- qrF$rank
    
    if (rnk < ncol(F)) {
        piv <- qrF$pivot[seq_len(rnk)]
        F <- F[, piv, drop = FALSE]
        contrast <- contrast[colnames(F)]
    }
    
    list(F = F, contrast = contrast)
}

# check whether contrast is usable in this repaired design
# if contrast is all zeros after repair -> nothing to test
#' @keywords internal
isContrastNonzero <- function(contrast, tol = 1e-12) {
    any(abs(contrast) > tol)
}


#' @keywords internal
coerceDEResultToDF <- function(res) {
  if (is.null(res)) return(NULL)
  if (inherits(res, "TopTags")) {
    return(as.data.frame(res))
  }
  if (is.list(res) && !is.data.frame(res) && "table" %in% names(res) && is.data.frame(res$table)) {
    return(res$table)
  }
  if (is.matrix(res)) return(as.data.frame(res))
  if (!is.data.frame(res)) {
    res <- tryCatch(as.data.frame(res), error = function(e) NULL)
  }
  res
}

#' Standardize DE result column names across DE methods
#' @keywords internal
standardizeDEColumns <- function(res) {
    if (is.null(res) || !is.data.frame(res)) return(res)
    if (!"pvalue" %in% names(res)) {
        if ("PValue" %in% names(res)) res$pvalue <- res$PValue
        if ("p.value" %in% names(res)) res$pvalue <- res[["p.value"]]
    }
    if (!"padj" %in% names(res)) {
        if ("FDR" %in% names(res)) res$padj <- res$FDR
        if ("adj.P.Val" %in% names(res)) res$padj <- res[["adj.P.Val"]]
    }
    if (!"log2FoldChange" %in% names(res)) {
        if ("logFC" %in% names(res)) res$log2FoldChange <- res$logFC
        if ("log_fold_change" %in% names(res)) res$log2FoldChange <- res$log_fold_change
    }
    # Force numeric 
    for (nm in intersect(c("log2FoldChange", "pvalue", "padj"), names(res))) {
        if (!is.numeric(res[[nm]])) {
            res[[nm]] <- suppressWarnings(as.numeric(as.character(res[[nm]])))
        }
    }
    
    res
}

#' Subsample cells per cell type per sample
#' @param raw.mats.sub named list of raw count matrices (cells x genes)
#' @keywords internal
subsampleCellsPerTypePerSample <- function(raw.mats.sub, cell.groups, sample.per.cell, core.samples,
                                           n.cells.subsample, min.cell.count = 1L, seed = NULL) {
    if (is.null(n.cells.subsample) || is.na(n.cells.subsample)) return(raw.mats.sub)
    if (length(core.samples) == 0L) return(raw.mats.sub)
    
    if (!is.null(seed)) set.seed(seed)
    
    # cell.groups must be named by cell IDs
    if (is.null(names(cell.groups))) stop("cell.groups must be a named factor/vector (names = cell IDs)")
    if (is.null(names(sample.per.cell))) stop("sample.per.cell must be named (names = cell IDs)")
    
    # Only cells that exist in both mappings
    common.cells <- intersect(names(cell.groups), names(sample.per.cell))
    if (length(common.cells) == 0L) stop("No overlap between names(cell.groups) and names(sample.per.cell)")
    
    # For each core sample, pick cells within each cell type
    for (s in intersect(names(raw.mats.sub), core.samples)) {
        m <- raw.mats.sub[[s]]
        if (is.null(m)) next
        
        # infer cell IDs present in this sample's matrix
        cell.ids <- rownames(m) # assumes cells x genes format
        if (is.null(cell.ids)) stop("raw.mats matrices must have rownames = cell IDs")
        
        # only keep cells mapped to this sample
        cell.ids <- intersect(cell.ids, common.cells)
        cell.ids <- cell.ids[sample.per.cell[cell.ids] == s]
        if (length(cell.ids) == 0L) next
        
        # subset groups for these cells
        cg <- cell.groups[cell.ids]
        cg <- droplevels(as.factor(cg))
        
        keep.cells <- character(0)
        
        for (ct in levels(cg)) {
            cells.ct <- cell.ids[cg == ct]
            n.avail <- length(cells.ct)
            
            # keep all if small; else subsample
            if (n.avail >= n.cells.subsample) {
              if (!is.null(seed)) {set.seed(hash.seed(seed, ct, s, "fix.cells"))}
              keep.cells <- c(keep.cells, sample(cells.ct, n.cells.subsample))
            } else if (n.avail >= min.cell.count) {
                keep.cells <- c(keep.cells, cells.ct)
            } else {
                # too few cells for this cell type; drop them (forces collapseCellsByType to drop this type in this sample)
            }
        }
        
        # apply
        if (length(keep.cells) == 0L) {
            # drop sample entirely 
            raw.mats.sub[[s]] <- m[0, , drop = FALSE]
        } else {
            raw.mats.sub[[s]] <- m[keep.cells, , drop = FALSE]
        }
    }
    
    raw.mats.sub
}

#' Get core samples from model
#' @param model list returned by buildDesignMatrices (must contain F and core.rows)
#' @param samples character vector of sample IDs
#' @return character vector of core sample IDs
#' @keywords internal
getCoreSamples <- function(model, samples) {
    all.samples <- rownames(model$F)
    if (is.null(all.samples)) return(samples)
    
    core.mask <- rep(TRUE, length(all.samples)); names(core.mask) <- all.samples
    if (!is.null(model$core.rows)) {
        core.mask <- as.logical(model$core.rows)
        names(core.mask) <- all.samples
    }
    core.samples <- names(core.mask)[core.mask]
    intersect(samples, core.samples)
}

#' @keywords internal
hash.seed <- function(..., mod = .Machine$integer.max) {
  x <- paste0(list(...), collapse = "||")
  # simple deterministic hash
  s <- sum(utf8ToInt(x) * seq_along(utf8ToInt(x)))
  as.integer(abs(s) %% mod)
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
estimateDEForTypeDESeq <- function(cm, sample.meta, model, test.type = "wald",
                                   cooksCutoff = FALSE, independentFiltering = TRUE, ...) {

  if (is.null(model$F)) stop("model must contain 'F' (design matrix).")
  if (is.null(rownames(model$F))) stop("model$F must have rownames (sample IDs).")
  if (is.null(colnames(cm))) stop("cm must have colnames (sample IDs).")
  if (is.null(rownames(sample.meta))) stop("sample.meta must have rownames (sample IDs).")

  # Align samples
  samp <- colnames(cm)
  if (!all(samp %in% rownames(model$F))) {
    stop("Some cm columns are missing from model$F rownames.")
  }
  if (!all(samp %in% rownames(sample.meta))) {
    stop("Some cm columns are missing from sample.meta rownames.")
  }

  X <- model$F[match(samp, rownames(model$F)), , drop = FALSE]
  if (!identical(rownames(X), samp)) stop("Failed to align model$F rows to cm columns.")

  # Build DESeq2 object (design is dummy; we pass X explicitly)
  coldata <- sample.meta[samp, , drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cm, colData = coldata, design = ~ 1)

  # Fit size factors + dispersions 
  dds <- DESeq2::estimateSizeFactors(dds, quiet = TRUE)
  dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)

  if (tolower(test.type) %in% c("wald", "")) {
    # Wald fit using custom model matrix
    dds <- DESeq2::nbinomWaldTest(dds, modelMatrix = X, quiet = TRUE)
    if (is.null(model$contrast.F)) {
      stop("model must contain 'contrast.F' for Wald test.")
    }
    cF <- model$contrast.F
    # resultsNames(dds) for a custom modelMatrix are the coefficient names DESeq2 stored
    rn <- DESeq2::resultsNames(dds)
    # numeric contrast vector aligned to resultsNames
    contrast.vec <- setNames(rep(0, length(rn)), rn)
    # assumes names(model$contrast.F) correspond to coef names in X
    common <- intersect(names(cF), rn)
    if (length(common) == 0L) {
      stop("names(model$contrast.F) do not match DESeq2 resultsNames(dds). ",
           "Have: ", paste(rn, collapse = ", "))
    }
    contrast.vec[common] <- cF[common]

    res <- DESeq2::results(dds, contrast = contrast.vec, cooksCutoff = cooksCutoff,
                           independentFiltering = independentFiltering)

    return(as.data.frame(res))

  } else {    # Reduced model: prefer model$Z if available; else intercept-only
    if (!is.null(model$Z) && is.matrix(model$Z)) {
      if (is.null(rownames(model$Z))) stop("model$Z must have rownames (sample IDs).")
      Z <- model$Z[match(samp, rownames(model$Z)), , drop = FALSE]
      if (!identical(rownames(Z), samp)) stop("Failed to align model$Z rows to cm columns.")
    } else {
      Z <- matrix(1, nrow = nrow(X), ncol = 1,
                  dimnames = list(rownames(X), "(Intercept)"))
    }
    dds <- DESeq2::nbinomLRT(dds, full = X, reduced = Z, quiet = TRUE)
    res <- DESeq2::results(dds, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)

    return(as.data.frame(res))
  }
  stop("Unknown DESeq2 test.type: ", test.type)
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
    
    # Ensure rownames are genes
    if (is.null(rownames(tab))) {
      stop("edgeR result table has no rownames (gene IDs).")
    }
    res <- as.data.frame(tab, stringsAsFactors = FALSE)
    if ("logFC" %in% names(res)) {
      res$log2FoldChange <- as.numeric(res$logFC)
     } else if (!"log2FoldChange" %in% names(res)) {
      res$log2FoldChange <- NA_real_
    }
    # Standard p-values + adjusted p-values
    res$pvalue <- as.numeric(res$PValue)
    res$padj   <- as.numeric(res$FDR)
    res$Gene <- rownames(res)

    return(res)
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

#' helper: get the DE table (data.frame) from either a df or a list(res=...)
#' @keywords internal
getResDf <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.list(x) && !is.null(x$res) && is.data.frame(x$res)) return(x$res)
  NULL
}

#' helper: choose an available statistic column (with fallbacks)
#' @keywords internal
pickStatCol <- function(df, var.to.sort) {
  if (is.null(df) || !is.data.frame(df)) return(NULL)
  if (var.to.sort %in% colnames(df)) return(var.to.sort)

  # p-value fallbacks (edgeR / misc)
  if (identical(var.to.sort, "pvalue")) {
    for (alt in c("PValue", "p.value", "p_val", "pval", "p")) {
      if (alt %in% colnames(df)) return(alt)
    }
  }

  # adjusted p-value fallbacks
  if (identical(var.to.sort, "padj")) {
    for (alt in c("FDR", "adj.P.Val", "qvalue", "q_value")) {
      if (alt %in% colnames(df)) return(alt)
    }
  }

  # last resort: common p-value-ish columns
  for (alt in c("PValue", "p.value", "p_val", "pval", "p", "FDR", "adj.P.Val")) {
    if (alt %in% colnames(df)) return(alt)
  }

  NULL
}

#' Summarize DE Resampling Results
#'
#' @param de.list list with DE results. Data frame with DE results are found in the first element
#' @param var.to.sort Variable to calculate ranks (default="pvalue")
#' 
#' @keywords internal
summarizeDEResamplingResults <- function(de.list, var.to.sort = "pvalue") {

  # initial run
  de.res <- de.list[[1]]

  for (cell.type in names(de.res)) {

    # initial result table must exist
    res0 <- getResDf(de.res[[cell.type]])
    if (is.null(res0)) {
      warning("Initial DE result for cell type ", cell.type, " is missing/invalid; skipping.")
      next
    }

    genes.init <- rownames(res0)
    if (is.null(genes.init) || length(genes.init) == 0L) next

    # Start with all genes from initial; we'll intersect with each resample
    genes.common <- genes.init

    # We'll collect one column per resample; start with 0 columns
    mx.stat <- matrix(
      nrow = length(genes.common),
      ncol = 0,
      dimnames = list(genes.common, character(0))
    )

    used.subsamples <- list()

    # Loop over resamples (exclude the initial)
    for (i in 2:length(de.list)) {

      # get this cell type's result, could be missing
      if (!(cell.type %in% names(de.list[[i]]))) next

      ri <- getResDf(de.list[[i]][[cell.type]])
      if (is.null(ri)) next

      # must have the statistic column (or skip)
      stat.col <- pickStatCol(ri, var.to.sort)
      if (is.null(stat.col)) next

      # intersect genes
      genes.common <- intersect(genes.common, rownames(ri))
      if (length(genes.common) == 0L) {
        # if no overlap anymore, can't proceed for this cell type
        mx.stat <- mx.stat[character(0), , drop = FALSE]
        break
      }

      # align mx.stat to updated genes.common and bind this resample stat
      mx.stat <- mx.stat[genes.common, , drop = FALSE]

      # extract as 2D column matrix
      v <- ri[genes.common, stat.col, drop = FALSE]

      # keep names for traceability
      col.name <- names(de.list)[i] %||% paste0("resample.", i)
      colnames(v) <- col.name

      mx.stat <- cbind(mx.stat, v)

      # store the whole resample object for this cell type (optional)
      used.subsamples[[col.name]] <- de.list[[i]][[cell.type]]
    }

    # If no subsample contributed, skip
    if (ncol(mx.stat) == 0L || nrow(mx.stat) == 0L) {
      warning("Cell type ", cell.type, " was not present (or had no usable '", var.to.sort,
              "') in any subsamples.")

      # still record that there were 0 usable subsamples
      if (is.list(de.res[[cell.type]]) && !is.null(de.res[[cell.type]]$res)) {
        de.res[[cell.type]]$res$n.subsamples.used <- 0L
        de.res[[cell.type]]$subsamples <- list()
      } else if (is.data.frame(de.res[[cell.type]])) {
        de.res[[cell.type]]$n.subsamples.used <- 0L
        attr(de.res[[cell.type]], "subsamples") <- list()
      } else if (is.list(de.res[[cell.type]])) {
        de.res[[cell.type]]$n.subsamples.used <- 0L
        de.res[[cell.type]]$subsamples <- list()
      }
      next
    }

    # Ensure 2D numeric matrix (protect against drop-to-vector issues)
    mx.stat <- as.matrix(mx.stat)
    storage.mode(mx.stat) <- "numeric"

    # Rank within each subsample column
    mx.rank <- apply(mx.stat, 2, rank)
    mx.rank <- as.matrix(mx.rank)

    stab.mean.rank <- rowMeans(mx.rank)
    stab.median.rank <- apply(mx.rank, 1, median)
    stab.var.rank <- if (ncol(mx.rank) >= 2) apply(mx.rank, 1, var) else rep(0, nrow(mx.rank))

    # Attach stability stats back to the *initial* result table order
    out.median <- setNames(rep(NA_real_, length(genes.init)), genes.init)
    out.mean   <- setNames(rep(NA_real_, length(genes.init)), genes.init)
    out.var    <- setNames(rep(NA_real_, length(genes.init)), genes.init)

    out.median[genes.common] <- stab.median.rank
    out.mean[genes.common]   <- stab.mean.rank
    out.var[genes.common]    <- stab.var.rank

    # Write back depending on storage shape of initial
    if (is.list(de.res[[cell.type]]) && !is.null(de.res[[cell.type]]$res)) {
      de.res[[cell.type]]$res$stab.median.rank <- out.median
      de.res[[cell.type]]$res$stab.mean.rank   <- out.mean
      de.res[[cell.type]]$res$stab.var.rank    <- out.var
      de.res[[cell.type]]$res$n.subsamples.used <- ncol(mx.stat)

      # Save only subsamples that were used
      de.res[[cell.type]]$subsamples <- used.subsamples

    } else if (is.data.frame(de.res[[cell.type]])) {
      de.res[[cell.type]]$stab.median.rank <- out.median
      de.res[[cell.type]]$stab.mean.rank   <- out.mean
      de.res[[cell.type]]$stab.var.rank    <- out.var
      de.res[[cell.type]]$n.subsamples.used <- ncol(mx.stat)

      # keep parity with old output:
      attr(de.res[[cell.type]], "subsamples") <- used.subsamples

    } else {
      # unknown shape; do best effort
      de.res[[cell.type]]$stab.median.rank <- out.median
      de.res[[cell.type]]$stab.mean.rank   <- out.mean
      de.res[[cell.type]]$stab.var.rank    <- out.var
      de.res[[cell.type]]$n.subsamples.used <- ncol(mx.stat)
      de.res[[cell.type]]$subsamples <- used.subsamples
    }
  }

  de.res
}


#' Append statistics to DE results
#' 
#' @param de.list list with DE results
#' @param expr.frac.per.type 
#' 
#' @keywords internal
appendStatisticsToDE <- function(de.list, expr.frac.per.type) {

  for (n in names(de.list)) {

    # Skip cell types with no DE table
    if (is.null(de.list[[n]]$res) || nrow(de.list[[n]]$res) == 0) {
      next
    }

    res <- de.list[[n]]$res

    # Ensure there is a Gene column
    if (!("Gene" %in% colnames(res))) {
      res <- tibble::rownames_to_column(as.data.frame(res), var = "Gene")
    }

    # Compute CellFrac safely (match genes)
    cf <- as.numeric(expr.frac.per.type[res$Gene, n])

    # Compute SampleFrac if cm exists
    if (!is.null(de.list[[n]]$cm)) {
      sf_all <- Matrix::rowMeans(de.list[[n]]$cm > 0)
      sf <- as.numeric(sf_all[res$Gene])
    } else {
      sf <- rep(NA_real_, nrow(res))
    }

    # Add columns and restore rownames
    res <- dplyr::mutate(res, CellFrac = cf, SampleFrac = sf) |>
      as.data.frame(stringsAsFactors = FALSE)

    rownames(res) <- res$Gene
    de.list[[n]]$res <- res
  }

  de.list
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

  # --- checks + alignment  ---
  if (is.null(colnames(counts))) {
    stop("getPerCellTypeGeneFilter(): counts must have colnames (cell IDs) for alignment.")
  }
  if (is.null(names(cell.groups))) {
    stop("getPerCellTypeGeneFilter(): cell.groups must be a named vector (names are cell IDs).")
  }

  common.cells <- intersect(colnames(counts), names(cell.groups))
  if (length(common.cells) == 0) {
    stop("getPerCellTypeGeneFilter(): no overlapping cell IDs between counts and cell.groups.")
  }

  counts <- counts[, common.cells, drop = FALSE]
  cell.groups <- cell.groups[common.cells]
  stopifnot(identical(colnames(counts), names(cell.groups)))

  # binarize counts (>1)
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  counts@x <- as.numeric(counts@x > 1)

  # build indicator 
  cell.types <- sort(unique(cell.groups))
  cell.type.indicator <- Matrix::sparseMatrix(
    i = seq_along(cell.groups),
    j = match(cell.groups, cell.types),
    x = 1,
    dims = c(length(cell.groups), length(cell.types)),
    dimnames = list(names(cell.groups), cell.types)
  )

  # compute fraction of expressing cells per gene per celltype
  counts.per.type <- counts %*% cell.type.indicator   # (genes x cells) %*% (cells x types) = genes x types
  n.cells <- Matrix::colSums(cell.type.indicator)
  expr.frac <- sweep(counts.per.type, 2, n.cells, FUN = "/")
  filters <- expr.frac > threshold

  filters.list <- lapply(seq_len(ncol(filters)), function(i) {
    setNames(as.vector(filters[, i]), rownames(counts))
  })
  names(filters.list) <- colnames(filters)
  filters.list
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

