#' @title Cacoa R6 class
#'
#' @import methods enrichplot dplyr
#' @export Cacoa
#' @exportClass Cacoa
#' @param sample.groups a two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param n.cores number of cores for parallelization
#' @param verbose show progress (default: stored value)
#' @param name field name where the test results are stored
#' @param n.top.genes number of top genes for estimation
#' @param gene.selection a method to select top genes, "z" selects genes by cluster-free Z-score change, "lfc" uses log2(fold-change) instead,
#' "expression" picks the most expressed genes and "od" picks overdispersed genes.  Default: "z".
#' @param excluded.genes list of genes to exclude during estimation. For example, a list of mitochondrial genes.
Cacoa <- R6::R6Class("Cacoa", lock_objects=FALSE,
  public = list(
    n.cores = 1,

    #' @field verbose print diagnostic messages
    verbose = FALSE,

    #' @field test.results results of the estimations, ready to use
    test.results = list(),

    #' @field cache intermediate results of the estimations, which can be used during some other computations
    cache = list(),

    #' @field data.object the main object storing data (Conos or Seurat)
    data.object = list(),

    #' @field sample.groups 2-factor vector with annotation of groups/condition per sample
    sample.groups = NULL,

    #' @field cell.groups named factor with cell names with cluster per cell
    cell.groups = NULL,

    #' @field embedding a 2D embedding to visualize the cells in
    embedding = NULL,

    #' @field sample.per.cell named factor with cell names
    sample.per.cell = NULL,

    #' @field ref.level reference level for sample.group vector
    ref.level = NULL,

    #' @field target.level target/disease level for sample.group vector
    target.level = NULL,

    #' @field sample.groups.palette a color palette for the sample.groups
    sample.groups.palette = NULL,

    #' @field cell.groups.palette a color palette for the cell.groups
    cell.groups.palette = NULL,

    #' @field plot.theme ggplot2 theme for all plots
    plot.theme=NULL,

    #' @field plot.params parameters, forwarded to all `plotEmbedding` calls
    plot.params=NULL,

    initialize=function(data.object, sample.groups=NULL, cell.groups=NULL, sample.per.cell=NULL, ref.level=NULL,
                        target.level=NULL, sample.groups.palette=NULL, cell.groups.palette=NULL,
                        embedding=extractEmbedding(data.object), graph.name=NULL, n.cores=1, verbose=TRUE,
                        plot.theme=theme_bw(), plot.params=NULL) {
      if ('Cacoa' %in% class(data.object)) { # copy constructor
        for (n in ls(data.object)) {
          if (!is.function(get(n, data.object))) assign(n, get(n, data.object), self)
        }

        return()
      }

      if (!is.null(sample.groups)) {
        if (length(unique(sample.groups)) != 2) stop("sample.groups must have exactly two levels")
        if (is.null(ref.level) && !is.null(target.level)) ref.level <- setdiff(unique(sample.groups), target.level)[1]
        if (!is.null(ref.level) && is.null(target.level)) target.level <- setdiff(unique(sample.groups), ref.level)[1]
        if (length(setdiff(sample.groups, c(ref.level, target.level)) > 0))
          stop("sample.groups must contain only ref.level '", ref.level, "' and target.level '", target.level, "'")
      }

      if (is.null(ref.level) || is.null(target.level))
        stop("Both ref.level and target.level must be provided")

      self$n.cores <- n.cores
      self$verbose <- verbose
      self$ref.level <- ref.level
      self$target.level <- target.level

      if ("Seurat" %in% class(data.object)) {
        if (is.null(sample.groups) || is.null(sample.per.cell))
          stop("Both sample.groups and sample.per.cell must be provided for Seurat objects")
        data.object$sample.per.cell <- sample.per.cell
        if (is.null(graph.name)) warning("No graph.name provided. The algorithm will use the first available graph.")
        data.object@misc$graph.name <- graph.name
      } else if (('Conos' %in% class(data.object))) {
        if (!is.null(graph.name)) warning("graph.name is not supported for Conos objects")
      } else {
        warning("Many function may be not supported for an object of class ", class(data.object));
        if (is.null(sample.groups) || is.null(sample.per.cell) || is.null(cell.groups))
          stop("All sample.groups, sample.per.cell and cell.groups must be provided")
      }

      if (any(c("dgCMatrix", "dgTMatrix", "dgEMatrix", "matrix") %in% class(data.object))) {
        data.object %<>% as("dgCMatrix") %>% Matrix::t()
        if (max(abs(round(data.object@x) - data.object@x)) < 1e-10) {
          message("Interpreting data.object as a raw count matrix")
          attr(data.object, "raw") <- TRUE
        } else {
          message("Interpreting data.object as a normalized count matrix")
          attr(data.object, "raw") <- FALSE
        }

        if (length(setdiff(names(sample.per.cell), rownames(data.object))) > 0)
          stop("All cells in the count matrix columns must be present in sample.per.cell")
        attr(data.object, 'sample.per.cell') <- as.factor(sample.per.cell)
      }

      self$data.object <- data.object

      if(is.null(sample.groups)) {
        self$sample.groups <- extractSampleGroups(data.object, ref.level, self$target.level)
      } else {
        self$sample.groups <- sample.groups <- as.factor(sample.groups)
      }

      if(is.null(cell.groups)) {
        self$cell.groups <- extractCellGroups(data.object)
      } else {
        self$cell.groups <- cell.groups <- as.factor(cell.groups)
      }

      if(is.null(sample.per.cell)) {
        self$sample.per.cell <- extractSamplePerCell(data.object) %>% as.factor()
      } else {
        self$sample.per.cell <- as.factor(sample.per.cell)
      }

      if(is.null(sample.groups.palette)) {
        self$sample.groups.palette <- c("#d73027", "#4575b4") %>%
          setNames(c(self$target.level, self$ref.level))
      } else {
        self$sample.groups.palette <- sample.groups.palette
      }

      if(is.null(cell.groups.palette)) {
        self$cell.groups.palette <- levels(cell.groups) %>% length() %>%
          rainbow(s=0.9,v=0.9) %>% setNames(levels(cell.groups))
      } else {
        self$cell.groups.palette <- cell.groups.palette
      }

      self$plot.theme <- plot.theme
      self$embedding <- embedding;
      self$plot.params <- plot.params
    },

    ### Expression shifts

    #' @description  Calculate expression shift magnitudes of different clusters between conditions
    #' @param cell.groups Named cell group factor with cell names (default: stored vector)
    #' @param dist 'cor' - correlation distance, 'l1' - manhattan distance or 'l2' - euclidean (default depends on dimensionality)
    #' @param within.group.normalization Normalize the shift magnitude by the mean magnitude of within-group variation (default=`TRUE`)
    #' @param n.cells Number of cells to subsmaple across all samples (if not specified, defaults to the total size of the smallest cell cluster)
    #' @param n.top.genes Number of top highest-expressed genes to consider (default: all genes)
    #' @param n.subsamples Number of samples to draw (default=100)
    #' @param min.cells Minimum number of cells per cluster/per sample to be included in the analysis (default=10)
    #' @param n.cores Number of cores (default: stored integer)
    #' @param name Test name (default="expression.shifts")
    #' @return a list include
    #'   - `dist.df`: a table with cluster distances (normalized if within.gorup.normalization=TRUE), cell type and the number of cells # TODO: update
    #'   - `p.dist.info`: list of distance matrices per cell type
    #'   - `sample.groups`: filtered sample groups
    #'   - `cell.groups`: filtered cell groups
    estimateExpressionShiftMagnitudes=function(cell.groups=self$cell.groups, dist=NULL, dist.type="cross.both",
                                               min.cells.per.sample=10, min.samp.per.type=2, min.gene.frac=0.01,
                                               ref.level=self$ref.level, sample.groups=self$sample.groups,
                                               verbose=self$verbose, n.cores=self$n.cores, name="expression.shifts",
                                               n.permutations=1000, p.adjust.method='BH', genes=NULL, n.pcs=NULL,
                                               top.n.genes=NULL, ...) {
      count.matrices <- extractRawCountMatrices(self$data.object, transposed=TRUE)

      if (verbose) cat("Filtering data... ")
      shift.inp <- filterExpressionDistanceInput(
        count.matrices, cell.groups=cell.groups,
        sample.per.cell=self$sample.per.cell, sample.groups=self$sample.groups,
        min.cells.per.sample=min.cells.per.sample, min.samp.per.type=min.samp.per.type,
        min.gene.frac=min.gene.frac, genes=genes, verbose=verbose
      )
      if (verbose) cat("done!\n")

      if (!is.null(n.pcs)) {
        if (!is.null(top.n.genes) && n.pcs > top.n.genes) {
          n.pcs <- top.n.genes - 1
          warning("n.pcs can't be larger than top.n.genes - 1, setting it to ", n.pcs)
        }

        n.samps.per.type <- shift.inp$cm.per.type %>% sapply(nrow)
        affected.types <- which(n.samps.per.type <= n.pcs)
        if (length(affected.types) > 0) {
          affected.types %<>% names() %>% paste(collapse=", ")
          n.pcs <- min(n.samps.per.type) - 1
          warning("Cell types '", affected.types, "' don't have enough samples present. Setting n.pcs to ", n.pcs,
                  ". Consider increasing min.samp.per.type.")
        }
      }

      self$test.results[[name]] <- shift.inp %$%
        estimateExpressionShiftMagnitudes(
          cm.per.type, sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=self$sample.per.cell,
          dist=dist, dist.type=dist.type, verbose=verbose, ref.level=ref.level,
          n.permutations=n.permutations, top.n.genes=top.n.genes, n.pcs=n.pcs, n.cores=n.cores, ...
        )

      return(invisible(self$test.results[[name]]))
    },

    estimateCommonExpressionShiftMagnitudes=function(cell.groups=self$cell.groups, name='common.expression.shifts',
                                                     min.cells.per.sample=10, min.samp.per.type=2, min.gene.frac=0.01,
                                                     genes=NULL, n.permutations=1000, trim=0.2, p.adjust.method="BH",
                                                     verbose=self$verbose, n.cores=self$n.cores, ...) {
      if (verbose) cat("Filtering data... ")
      shift.inp <- extractRawCountMatrices(self$data.object, transposed=TRUE) %>%
        filterExpressionDistanceInput(
          cell.groups=cell.groups, sample.per.cell=self$sample.per.cell, sample.groups=self$sample.groups,
          min.cells.per.sample=min.cells.per.sample, min.samp.per.type=min.samp.per.type, min.gene.frac=min.gene.frac,
          genes=genes
        )
      if (verbose) cat("done!\n")

      sample.groups <- shift.inp$sample.groups

      if (verbose) cat('Calculating distances ... ')

      dists.norm <- list()
      res.per.type <- levels(shift.inp$cell.groups) %>% sccore:::sn() %>% plapply(function(ct) {
        cm.norm <- t(shift.inp$cm.per.type[[ct]])

        dists <- consensusShiftDistances(cm.norm, sample.groups, ...)
        obs.diff <- mean(dists, trim=trim)
        randomized.dists <- lapply(1:n.permutations, function(i) {
          sg.shuff <- sample.groups[colnames(cm.norm)] %>% as.factor() %>%
            droplevels() %>% {setNames(sample(.), names(.))}
          consensusShiftDistances(cm.norm, sg.shuff, ...) %>% mean(trim=trim)
        }) %>% unlist()

        pvalue <- (sum(randomized.dists >= obs.diff, na.rm=TRUE) + 1) / (sum(!is.na(randomized.dists)) + 1)
        dists <- dists - median(randomized.dists, na.rm=TRUE)
        list(dists=dists, pvalue=pvalue)
      }, progress=verbose, n.cores=n.cores, mc.preschedule=TRUE, fail.on.error=TRUE)

      if (verbose) cat("done!\n")

      pvalues <- sapply(res.per.type, `[[`, "pvalue")
      dists.per.type <- lapply(res.per.type, `[[`, "dists")
      padjust <- p.adjust(pvalues, method=p.adjust.method)

      self$test.results[[name]] <- list(dists.per.type=dists.per.type, pvalues=pvalues, padjust=padjust)
      return(invisible(self$test.results[[name]]))
    },

    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes() (shift.type="normal") or
    #'   cao$estimateCommonExpressionShiftMagnitudes() (shift.type="common")
    #' @param name - results slot name (default: 'expression.shifts')
    #' @param show.jitter whether to show indivudal data points (default: FALSE)
    #' @param jitter.alpha transparency value for the data points (default: 0.05)
    #' @param type - type of a plot "bar" (default) or "box"
    #' @param notch - whether to show notches in the boxplot version (default=TRUE)
    #' @return A ggplot2 object
    plotExpressionShiftMagnitudes=function(name=NULL, type='box', notch=TRUE, show.jitter=TRUE, jitter.alpha=0.05,
                                           show.pvalues=c("adjusted", "raw", "none"), shift.type=c("normal", "common"),
                                           ylab='normalized expression distance', ...) {
      show.pvalues <- match.arg(show.pvalues)
      shift.type <- match.arg(shift.type)
      if (is.null(name)) {
        name <- if (shift.type == "normal") "expression.shifts" else "common.expression.shifts"
      }

      func.name <- if (shift.type == "normal") {"estimateExpressionShiftMagnitudes()"}
        else {"estimateCommonExpressionShiftMagnitudes()"}
      res <- private$getResults(name, func.name)
      df <- names(res$dists.per.type) %>%
        lapply(function(n) data.frame(value=res$dists.per.type[[n]], Type=n)) %>%
        do.call(rbind, .) %>% na.omit()

      if (show.pvalues == "adjusted") {
        pvalues <- res$padjust
      } else if (show.pvalues == "raw") {
        pvalues <- res$pvalues
      } else {
        pvalues <- NULL
      }

      plotMeanMedValuesPerCellType(
        df, pvalues=pvalues, show.jitter=show.jitter,jitter.alpha=jitter.alpha, notch=notch, type=type,
        palette=self$cell.groups.palette, ylab=ylab, plot.theme=self$plot.theme, yline=0.0, ...
      )
    },

    #' @description Alias for estimateDEPerCellType
    estimatePerCellTypeDE=function(...) {
      .Deprecated("cao$estimateDEPerCellType")
      return(self$estimateDEPerCellType(...))
    },

    #' @description Estimate differential gene expression per cell type between conditions
    #' @param cell.groups factor specifying cell types (default from self)
    #' @param sample.groups 2-factor vector with annotation of groups/condition per sample (default from self)
    #' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy (default from self)
    #' @param target.level Reference level in 'sample.groups', e.g., case, diseased (default from self)
    #' @param common.genes Only investigate common genes across cell groups (default=FALSE)
    #' @param cooks.cutoff cooksCutoff for DESeq2 (default=FALSE)
    #' @param test which DESeq2 test to use (options: "LRT" (default), "Wald")
    #' @param cooks.cutoff cooksCutoff for DESeq2 (default=FALSE)
    #' @param min.cell.count minimum number of cells that need to be present in a given cell type in a given sample in order to be taken into account (default=10)
    #' @param independent.filtering independentFiltering parameter for DESeq2 (default=FALSE)
    #' @param resampling.method which resampling method should be used "loo" for leave-one-out or "bootstrap", (default:NULL no resampling)
    #' @param name slot in which to save the results (default: 'de')
    #' @return A list of DE genes
    estimateDEPerCellType=function(cell.groups=self$cell.groups, sample.groups=self$sample.groups,
                                   ref.level=self$ref.level, target.level=self$target.level, name='de',
                                   test='DESeq2.Wald', resampling.method=NULL, n.resamplings=30, seed.resampling=239,
                                   min.cell.frac=0.05, covariates=NULL, common.genes=FALSE, n.cores=self$n.cores,
                                   cooks.cutoff=FALSE, independent.filtering=FALSE, min.cell.count=10,
                                   n.cells.subsample=NULL, verbose=self$verbose, fix.n.samples=NULL, ...) {
      set.seed(seed.resampling)
      if (!is.list(sample.groups)) {
        sample.groups %<>% {split(names(.), . == ref.level)} %>% setNames(c(target.level, ref.level))
      }

      possible.tests <- c('DESeq2.Wald', 'DESeq2.LRT', 'edgeR',
                          'Wilcoxon.edgeR', 'Wilcoxon.DESeq2', 'Wilcoxon.totcount',
                          't-test.edgeR', 't-test.DESeq2', 't-test.totcount',
                          'limma-voom')

      if (tolower(test) == tolower('DESeq2')) test <- paste(test, 'Wald', sep='.')
      if (tolower(test) %in% tolower(c('Wilcoxon', 't-test')))  test <- paste(test, 'edgeR', sep='.')

      if (!(tolower(test) %in% tolower(possible.tests)))
        stop('Test ', test, ' is not supported. Available tests: ', paste(possible.tests, collapse=', '))

      # s.groups.new contains list of case/control groups of samples to run DE on.
      # First element in s.groups.new corresponds to the initial grouping.

      if (!is.null(n.cells.subsample) && is.null(resampling.method)) resampling.method <- 'fix.cells'
      s.groups.new <- list(initial=sample.groups)
      max.cell.count <- Inf
      fix.samples <- NULL
      # If resampling is defined, new contrasts will append to s.groups.new
      if (!is.null(resampling.method) && (n.resamplings != 0)) {
        s.groups.new %<>% c(
          prepareSamplesForDE(sample.groups, resampling.method=resampling.method, n.resamplings=n.resamplings)
        )

        if (resampling.method == 'fix.samples') {
          if (is.null(fix.n.samples)) stop("fix.n.samples must be provided for resampling.method='fix.samples'")
          fix.samples <- fix.n.samples
        }
      }

      if (!is.null(n.cells.subsample)) {
        if (verbose) message('Number of cell counts is fixed to ', n.cells.subsample)
        max.cell.count <- min.cell.count <- n.cells.subsample
      }

      raw.mats <- extractRawCountMatrices(self$data.object, transposed=TRUE)

      expr.fracs <- self$getJointCountMatrix() %>% getExpressionFractionPerGroup(cell.groups)
      gene.filter <- (expr.fracs > min.cell.frac)

      # parallelize the outer loop if subsampling is on
      outer.multicore <- (length(s.groups.new) >= n.cores) && (n.cores > 1)
      inner.verbose <- (length(s.groups.new) == 1) || (!outer.multicore && verbose > 1)
      de.res <- names(s.groups.new) %>% sn() %>% plapply(function(resampling.name) {
        estimateDEPerCellTypeInner(
          raw.mats=raw.mats, cell.groups=cell.groups, s.groups=s.groups.new[[resampling.name]],
          ref.level=ref.level, target.level=target.level, common.genes=common.genes,
          cooks.cutoff=cooks.cutoff, min.cell.count=min.cell.count, max.cell.count=max.cell.count,
          independent.filtering=independent.filtering, test=test, meta.info=covariates, gene.filter=gene.filter,
          fix.n.samples=(if (resampling.name == 'initial') NULL else fix.samples),
          n.cores=ifelse(outer.multicore, 1, n.cores),
          return.matrix=(resampling.name == 'initial'),
          verbose=(inner.verbose & verbose), ...
        )
      }, n.cores=ifelse(outer.multicore, n.cores, 1), progress=(!inner.verbose & verbose))

      # if resampling: calculate median and variance on ranks after resampling
      de.res <- if (length(de.res) > 1) summarizeDEResamplingResults(de.res) else de.res[[1]]
      de.res %<>% appendStatisticsToDE(expr.fracs)
      self$test.results[[name]] <- de.res

      # TODO: add overall p-adjustment

      return(invisible(self$test.results[[name]]))
    },

    #' @description  Plot DE stability per cell type
    #' @param name - results slot name (default: 'de')
    #' @param show.pairs transparency value for the data points (default: 0.05)
    #' @param notch - whether to show notches in the boxplot version (default=TRUE)
    #' @return A ggplot2 object
    estimateDEStabilityPerCellType=function(de.name='de', name='de.jaccards', top.n.genes=NULL, p.val.cutoff=NULL) {
      if(!is.null(p.val.cutoff) & !is.null(top.n.genes))
        stop('Only one threshold (top.n.genes or p.val.cutoff) should be provided')
      if(is.null(p.val.cutoff) & is.null(top.n.genes))
        stop('At least one threshold (top.n.genes or p.val.cutoff) should be provided')

      de.res <- private$getResults(de.name, 'estimateDEPerCellType()')

      if(!all(sapply(names(de.res), function(x) 'subsamples' %in% names(de.res[[x]]))))
        stop('Resampling was not performed')

      jaccards <- estimateStabilityPerCellType(de.res, top.n.genes, p.val.cutoff)
      self$test.results[[name]] <- jaccards
    },

    estimateDEStabilityPerTest=function(de.names, name='jacc.per.test', top.n.genes=NULL, p.val.cutoff=NULL) {
      if(!is.null(p.val.cutoff) & !is.null(top.n.genes))
        stop('Only one threshold (top.n.genes or p.val.cutoff) should be provided')

      if(is.null(p.val.cutoff) & is.null(top.n.genes))
        stop('At least one threshold (top.n.genes or p.val.cutoff) should be provided')

      jaccards.all <- data.frame()
      for(de.name in de.names){
        if(!(de.name %in% names(self$test.results)) || (length(self$test.results[[de.name]]) == 0) ){
          message(paste0('DE by ', de.name, ' was not estimated', collapse = ' '))
          next
        }
        de.res <- private$getResults(de.name)

        if(!all(sapply(names(de.res), function(x) 'subsamples' %in% names(de.res[[x]]))))
          stop('Resampling was not performed')
        jaccards <- estimateStabilityPerCellType(de.res=de.res, top.n.genes=top.n.genes, p.val.cutoff=p.val.cutoff)
        # print(jaccards)

        jacc.medians <- sapply(unique(jaccards$group), function(x) median(jaccards$value[jaccards$group == x]))
        jaccards.tmp <- data.frame(group = de.name, value = jacc.medians,
                                   cmp = names(jacc.medians))
        jaccards.all <- rbind(jaccards.all, jaccards.tmp)
      }
      self$test.results[[name]] <- jaccards.all
    },


    estimateDEStabilityPerGene=function(de.name,
                                        top.n.genes = 500,
                                        p.adj.cutoff = NULL,
                                        visualize=FALSE) {
      de.res <- self$test.results[[de.name]]
      for(cell.type in names(de.res)) {
        if(!is.null(p.adj.cutoff)) {
          top.n.genes <- sum(de.res[[cell.type]]$res$padj <= p.adj.cutoff)
        }
        genes.tmp <- de.res[[cell.type]]$res$Gene
        tmp <- rep(0, length(genes.tmp))
        n.tmp <- 0
        for(i in 1:length(de.res[[cell.type]]$subsamples)){
          if(is.null(de.res[[cell.type]]$subsamples[[i]])) next
          tmp <- tmp + 1*(rank(de.res[[cell.type]]$subsamples[[i]][genes.tmp,'padj']) < top.n.genes)
          n.tmp <- n.tmp + 1
        }
        de.res[[cell.type]]$res$Stability <- tmp / n.tmp
      }
      self$test.results[[de.name]] <- de.res
      if(visualize){
        stability.all <- c()
        for(cell.type in names(de.res)) {
          tmp <- de.res[[cell.type]]$res
          stability.all <- rbind(stability.all,
                                 data.frame(rank = 1:nrow(tmp),
                                            stability = tmp$Stability, cell.type = cell.type))
        }

        p <- ggplot(stability.all, aes(rank, stability, colour = cell.type)) +
          xlim(0, top.n.genes) + geom_smooth(method = "loess") + self$plot.theme +
          xlab('Gene rank by p-value') + ylab('Fraction of LOOs') +
          scale_color_manual(values=self$cell.groups.palette) +
          guides(color=guide_legend(override.aes=list(fill=NA), ncol=1)) +
          theme(legend.text = element_text(size=6),legend.key.height= unit(0.1, 'cm'),
                legend.key.width = unit(0.2, 'cm'))
        return(p)
      }

    },


    estimateDEStabilityBetweenTests=function(de.names,
                                             name='jacc.bw.tests',
                                             top.n.genes = NULL,
                                             p.val.cutoff = NULL){

      if( (!is.null(p.val.cutoff)) & (!is.null(top.n.genes)) ) stop('Only one threshold (top.n.genes or p.val.cutoff) should be provided')
      if(is.null(p.val.cutoff) & is.null(top.n.genes)) stop('At least one threshold (top.n.genes or p.val.cutoff) should be provided')

      data.all = data_frame()
      cell.types = levels(self$cell.groups)
      for(cell.type in cell.types) {
        data.cell.type = data_frame()
        for(i in 1:length(de.names)) {
          for(j in 1:length(de.names)) {
            if(j <= i) next

            common.resamplings = intersect(names(self$test.results[[de.names[i]]][[cell.type]]$subsamples),
                                           names(self$test.results[[de.names[j]]][[cell.type]]$subsamples))
            if (length(common.resamplings) == 0){
              # message(paste('There is no corresponding resamplings for',
              #               de.names[i], 'and', de.names[j], 'in', cell.type))
              next
            }
            for(sample.name in common.resamplings) {

              s1 = self$test.results[[de.names[i]]][[cell.type]]$subsamples[[sample.name]]
              s2 = self$test.results[[de.names[j]]][[cell.type]]$subsamples[[sample.name]]
              if(is.null(p.val.cutoff)) {
                jac = jaccardPwTop(list(s1, s2), top.n.genes)
              } else {
                jac = jaccardPwPval(list(s1, s2), p.val.cutoff)
              }

              data.tmp = data_frame(s1 = de.names[i], s2 = de.names[j], val = jac$jac, cell.type = cell.type)
              data.cell.type = rbind(data.cell.type, data.tmp)
              data.tmp = data_frame(s1 = de.names[j], s2 = de.names[i], val = jac$jac, cell.type = cell.type)
              data.cell.type = rbind(data.cell.type, data.tmp)

            }
          }
          data.tmp = data_frame(s1 = de.names[i], s2 = de.names[i], val = 1, cell.type = cell.type)
          data.cell.type = rbind(data.cell.type, data.tmp)
        }
        data.cell.type <- data.cell.type[!duplicated(data.cell.type), ]
        data.all <- rbind(data.all, data.cell.type)
      }

      self$test.results[[name]] <- data.all
    },

    estimateDEStabilityTrend=function(de.name='de',
                                      name='de.trend',
                                      top.n.genes = c(100,200,300),
                                      p.val.cutoffs = NULL) {

      de.res <- private$getResults(de.name, 'estimateDEPerCellType()')

      data.all = data.frame()

      if(!is.null(p.val.cutoffs)) {
        cutoffs = p.val.cutoffs
      } else {
        cutoffs = top.n.genes
      }
      cutoffs = sort(cutoffs)

      for(i in 1:length(cutoffs)) {

        if(!is.null(p.val.cutoffs)) {
          jaccards <- estimateStabilityPerCellType(de.res = de.res, top.n.genes = NULL,
                                                   p.val.cutoff = cutoffs[i])
        } else {
          jaccards <- estimateStabilityPerCellType(de.res = de.res, top.n.genes = cutoffs[i],
                                                   p.val.cutoff = NULL)
        }

        jacc.medians <- sapply(unique(jaccards$group), function(x) mean(jaccards$value[jaccards$group == x]))

        if(length(jacc.medians) == 0) next

        data.tmp <- data.frame(group = cutoffs[i],
                              value = jacc.medians,
                              cmp = names(jacc.medians))

        data.all <- rbind(data.all, data.tmp)

      }
      self$test.results[[name]] <- data.all

    },

    estimateGOStabilityPerCellTypePermDir = function(de.test = 'de', top.n.de = 500, padj.go = 0.05,
                                                  name = 'go.jaccards') {

      test.name = de.test
      pref = paste('go', test.name, top.n.de, sep = '.')

      pref.init <- paste(pref, 'init', sep = '.')
      pref.resampling <- paste(pref.init, 'res', sep = '.')
      cell.types <- names(self$test.results[[pref.init]]$res)

      names.test.results <- names(self$test.results)

      jc.all <- c()
      # for(de.type in c('up', 'down', 'all')) {
      for(de.type in c('up', 'down')) {
        go.stability = list()
        for(cell.type in cell.types) {
          for(go.type in c('BP', 'MF', 'CC')) {
            if(!('res' %in% names(self$test.results[[pref.init]]))) next
            result <- self$test.results[[pref.init]]$res[[cell.type]][[go.type]][[de.type]]@result
            result <- result[,c('ID', 'pvalue', 'p.adjust')]
            colnames(result) <- c('ID', 'pvalue', 'padj')  # rename to be compatible with DE analysis
            if(go.type %in% c('MF', 'CC')){
              go.stability[[cell.type]]$res <- rbind(go.stability[[cell.type]]$res, result)
            } else {
              go.stability[[cell.type]]$res <- result
            }

            go.names.res <- unique(names.test.results[grepl(pref.resampling, names.test.results, fixed = TRUE)])

            for(go.name in go.names.res) {
              if(!('res' %in% names(self$test.results[[go.name]]))) next
              id.res <- gsub(paste(pref.init, '.', sep = ''), "", go.name)
              result <- self$test.results[[go.name]]$res[[cell.type]][[go.type]][[de.type]]@result
              result <- result[,c('ID', 'pvalue', 'p.adjust')]
              colnames(result) <- c('ID', 'pvalue', 'padj')  # rename to be compatible with DE analysis
              if(go.type %in% c('MF', 'CC')){
                go.stability[[cell.type]]$subsamples[[id.res]] <- rbind(go.stability[[cell.type]]$subsamples[[id.res]],
                                                                        result)
              } else {
                go.stability[[cell.type]]$subsamples[[id.res]] <- result
              }
            }
          }
        }
        # self$test.results[['tmp']] <- go.stability
        jc <- estimateStabilityPerCellType(go.stability, top.n.genes = NULL, p.val.cutoff = padj.go)
        if(nrow(jc) == 0) next
        jc$de.type <- de.type

        jc.all <- rbind(jc.all, jc)

      }
      self$test.results[[name]] <- jc.all
    },

    estimateGOStabilityPerCellTypePerm = function(de.test = 'de', top.n.de = 500, padj.go = 0.05,
                                              name = 'go.jaccards') {

      test.name = de.test
      pref = paste('go', test.name, top.n.de, sep = '.')

      pref.init <- paste(pref, 'init', sep = '.')
      pref.resampling <- paste(pref.init, 'res', sep = '.')
      cell.types <- names(self$test.results[[pref.init]]$res)

      names.test.results <- names(self$test.results)

      jc.all <- c()
      # for(de.type in c('up', 'down', 'all')) {

      go.stability = list()
      for(cell.type in cell.types) {
        for(de.type in c('up', 'down')) {
          for(go.type in c('BP', 'MF', 'CC')) {
            if(!('res' %in% names(self$test.results[[pref.init]]))) next
            result <- self$test.results[[pref.init]]$res[[cell.type]][[go.type]][[de.type]]@result
            result <- result[,c('ID', 'pvalue', 'p.adjust')]
            colnames(result) <- c('ID', 'pvalue', 'padj')  # rename to be compatible with DE analysis
            if(go.type %in% c('MF', 'CC')){
              go.stability[[cell.type]]$res <- rbind(go.stability[[cell.type]]$res, result)
            } else {
              go.stability[[cell.type]]$res <- result
            }

            go.names.res <- unique(names.test.results[grepl(pref.resampling, names.test.results, fixed = TRUE)])

            for(go.name in go.names.res) {
              if(!('res' %in% names(self$test.results[[go.name]]))) next
              id.res <- gsub(paste(pref.init, '.', sep = ''), "", go.name)
              result <- self$test.results[[go.name]]$res[[cell.type]][[go.type]][[de.type]]@result
              result <- result[,c('ID', 'pvalue', 'p.adjust')]
              colnames(result) <- c('ID', 'pvalue', 'padj')  # rename to be compatible with DE analysis
              if(go.type %in% c('MF', 'CC')){
                go.stability[[cell.type]]$subsamples[[id.res]] <- rbind(go.stability[[cell.type]]$subsamples[[id.res]],
                                                                        result)
              } else {
                go.stability[[cell.type]]$subsamples[[id.res]] <- result
              }
            }
          }
        }
        # self$test.results[['tmp']] <- go.stability
        jc <- estimateStabilityPerCellType(go.stability, top.n.genes = NULL, p.val.cutoff = padj.go)
        if(nrow(jc) == 0) next
        jc.all <- rbind(jc.all, jc)

      }
      self$test.results[[name]] <- jc.all
    },

    estimateGOStabilityPerCellTypePermGseaDir = function(de.test = 'de',padj.go = 0.05,
                                                  name = 'go.jaccards') {

      test.name = de.test
      pref = paste('gsea', test.name, sep = '.')

      pref.init <- paste(pref, 'init', sep = '.')
      pref.resampling <- paste(pref.init, 'res', sep = '.')
      cell.types <- names(self$test.results[[pref.init]]$res)

      names.test.results <- names(self$test.results)

      jc.all <- c()
      # for(de.type in c('up', 'down', 'all')) {
      for(de.type in c('up', 'down')) {
        go.stability = list()
        for(cell.type in cell.types) {
          for(go.type in c('BP', 'MF', 'CC')) {
            if(!('res' %in% names(self$test.results[[pref.init]]))) next
            result <- self$test.results[[pref.init]]$res[[cell.type]][[go.type]][[de.type]]
            result <- result[,c('pathway', 'pval', 'padj')]
            colnames(result) <- c('ID', 'pvalue', 'padj')  # rename to be compatible with DE analysis
            if(go.type %in% c('MF', 'CC')){
              go.stability[[cell.type]]$res <- rbind(go.stability[[cell.type]]$res, result)
            } else {
              go.stability[[cell.type]]$res <- result
            }

            go.names.res <- unique(names.test.results[grepl(pref.resampling, names.test.results, fixed = TRUE)])

            for(go.name in go.names.res) {
              if(!('res' %in% names(self$test.results[[go.name]]))) next
              id.res <- gsub(paste(pref.init, '.', sep = ''), "", go.name)
              result <- self$test.results[[go.name]]$res[[cell.type]][[go.type]][[de.type]]
              result <- result[,c('pathway', 'pval', 'padj')]
              colnames(result) <- c('ID', 'pvalue', 'padj')  # rename to be compatible with DE analysis
              if(go.type %in% c('MF', 'CC')){
                go.stability[[cell.type]]$subsamples[[id.res]] <- rbind(go.stability[[cell.type]]$subsamples[[id.res]],
                                                                        result)
              } else {
                go.stability[[cell.type]]$subsamples[[id.res]] <- result
              }
            }
          }
        }
        # self$test.results[['tmp']] <- go.stability
        jc <- estimateStabilityPerCellType(go.stability, top.n.genes = NULL, p.val.cutoff = padj.go)
        if(nrow(jc) == 0) next
        jc$de.type <- de.type

        jc.all <- rbind(jc.all, jc)

      }
      self$test.results[[name]] <- jc.all
    },

    estimateGOStabilityPerCellTypePermGsea = function(de.test = 'de',padj.go = 0.05,
                                                      name = 'go.jaccards') {

      test.name = de.test
      pref = paste('gsea', test.name, sep = '.')

      pref.init <- paste(pref, 'init', sep = '.')
      pref.resampling <- paste(pref.init, 'res', sep = '.')
      cell.types <- names(self$test.results[[pref.init]]$res)

      names.test.results <- names(self$test.results)

      jc.all <- c()


      go.stability = list()
      for(cell.type in cell.types) {
        for(de.type in c('up', 'down')) {
          for(go.type in c('BP', 'MF', 'CC')) {
            if(!('res' %in% names(self$test.results[[pref.init]]))) next
            result <- self$test.results[[pref.init]]$res[[cell.type]][[go.type]][[de.type]]
            result <- result[,c('pathway', 'pval', 'padj')]
            colnames(result) <- c('ID', 'pvalue', 'padj')  # rename to be compatible with DE analysis
            if(go.type %in% c('MF', 'CC')){
              go.stability[[cell.type]]$res <- rbind(go.stability[[cell.type]]$res, result)
            } else {
              go.stability[[cell.type]]$res <- result
            }

            go.names.res <- unique(names.test.results[grepl(pref.resampling, names.test.results, fixed = TRUE)])

            for(go.name in go.names.res) {
              if(!('res' %in% names(self$test.results[[go.name]]))) next
              id.res <- gsub(paste(pref.init, '.', sep = ''), "", go.name)
              result <- self$test.results[[go.name]]$res[[cell.type]][[go.type]][[de.type]]
              result <- result[,c('pathway', 'pval', 'padj')]
              colnames(result) <- c('ID', 'pvalue', 'padj')  # rename to be compatible with DE analysis
              if(go.type %in% c('MF', 'CC')){
                go.stability[[cell.type]]$subsamples[[id.res]] <- rbind(go.stability[[cell.type]]$subsamples[[id.res]],
                                                                        result)
              } else {
                go.stability[[cell.type]]$subsamples[[id.res]] <- result
              }
            }
          }
        }
        # self$test.results[['tmp']] <- go.stability
        jc <- estimateStabilityPerCellType(go.stability, top.n.genes = NULL, p.val.cutoff = padj.go)
        if(nrow(jc) == 0) next
        jc$de.type <- de.type

        jc.all <- rbind(jc.all, jc)

      }
      self$test.results[[name]] <- jc.all
    },

    estimateGOStabilityTrend=function(de.name='go',
                                      name='go.trend',
                                      top.n.genes = c(100,200,300),
                                      p.val.cutoffs = NULL) {

      de.res <- private$getResults(de.name, 'estimateOntology()')

      data.all = data.frame()

      if(!is.null(p.val.cutoffs)) {
        cutoffs = p.val.cutoffs
      } else {
        cutoffs = top.n.genes
      }
      cutoffs = sort(cutoffs)

      for(i in 1:length(cutoffs)) {

        if(!is.null(p.val.cutoffs)) {
          jaccards <- estimateStabilityPerCellType(de.res = de.res, top.n.genes = NULL,
                                                   p.val.cutoff = cutoffs[i])
        } else {
          jaccards <- estimateStabilityPerCellType(de.res = de.res, top.n.genes = cutoffs[i],
                                                   p.val.cutoff = NULL)
        }

        jacc.medians <- sapply(unique(jaccards$group), function(x) median(jaccards$value[jaccards$group == x]))

        if(length(jacc.medians) == 0) next

        data.tmp <- data.frame(group = cutoffs[i],
                               value = jacc.medians,
                               cmp = names(jacc.medians))

        data.all <- rbind(data.all, data.tmp)

      }
      self$test.results[[name]] <- data.all

    },

    plotGOStabilityPerCellType=function(name='go.jaccards',
                                        notch = FALSE,
                                        show.jitter = TRUE,
                                        jitter.alpha = 0.05,
                                        show.pairs = FALSE,
                                        sort.order = FALSE) {

      jaccards <- private$getResults(name, 'estimateGOStabilityPerCellType()')

      p <- plotStability(jaccards = jaccards,
                         notch = notch,
                         show.jitter = show.jitter,
                         jitter.alpha = jitter.alpha,
                         show.pairs = show.pairs,
                         sort.order = sort.order,
                         xlabel = 'Cell Type',
                         ylabel = 'Jaccard Index',
                         palette=self$cell.groups.palette,
                         plot.theme=self$plot.theme)
      # p <- p + facet_grid(rows = vars(go.type), cols = vars(de.type))
      p <- p + facet_grid(cols = vars(de.type)) + ylim(0, 1)
      return(p)
    },

    plotDEStabilityTrend=function(name='de.trend',
                                        notch = FALSE,
                                        show.jitter = TRUE,
                                        jitter.alpha = 0.05,
                                        show.pairs = FALSE,
                                        sort.order = FALSE) {

      jaccards <- private$getResults(name, 'estimateDEStabilityTrend()')
      p <- plotStability(jaccards = jaccards,
                         notch = notch,
                         show.jitter = show.jitter,
                         jitter.alpha = jitter.alpha,
                         show.pairs = show.pairs,
                         sort.order = sort.order,
                         xlabel = 'Cutoffs',
                         ylabel = 'Jaccard Index',
                         palette=self$cell.groups.palette,
                         plot.theme=self$plot.theme) +
        theme(legend.position = "right") +
        geom_jitter(aes(color = cmp)) + labs(color = "Cell types")
      return(p)
    },

    plotDEStabilityFDR=function(de.name='de',
                                p.adj.cutoffs = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2 ),
                                cell.types = NULL, type = 'relative'){
      de.res <- private$getResults(de.name, 'estimateDEPerCellType()')

      df.n.genes <- estimateDEStabilityFDR(de.res, p.adj.cutoffs)

      if(!is.null(cell.types)) df.n.genes <- df.n.genes[df.n.genes$Var2 %in% cell.types,]

      # df.n.genes <- df.n.genes[df.n.genes$type %in% c('common'),]

      df.n.common.genes <- df.n.genes[df.n.genes$type %in% c('common'),]
      df.n.all.genes <- df.n.genes[df.n.genes$type %in% c('all'),]
      df.n.init.genes <- df.n.genes[df.n.genes$type %in% c('init'),]

      if(type == 'relative') {
        de.plot <- df.n.common.genes
        de.plot$value <- de.plot$value / df.n.all.genes$value
      } else if (type == 'relative.to.init') {
        de.plot <- df.n.common.genes
        de.plot$value <- de.plot$value / df.n.init.genes$value
      } else if (type == 'common') {
        de.plot <- df.n.common.genes
      } else if (type == 'all') {
        de.plot <- df.n.all.genes
      } else if (type == 'init') {
        de.plot <- df.n.init.genes
      } else {
        stop('Incorrect')
      }

      ggplot(de.plot, aes(Var1, value, colour = Var2,
                             group = interaction(type, Var2))) +
        geom_line() + scale_y_continuous(trans='log10') +
        ggtitle(type)


    },

    plotDEStabilityFDR1loo=function(cell.type, de.name='de',
                                p.adj.cutoffs = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2 ),
                                cell.types = NULL, type = 'relative'){
      de.res <- private$getResults(de.name, 'estimateDEPerCellType()')

      df.n.genes <- estimateDEStabilityFDR1loo(de.res, p.adj.cutoffs)
      self$test.results[['fdr.stability']] <- df.n.genes
      df.n.genes <- df.n.genes[df.n.genes$cell.type == cell.type,]

      ggplot(df.n.genes, aes(x=fdr, y=frac,
                           group=fdr)) + geom_boxplot() + geom_jitter()

    },

    plotDEStabilityBetweenTests=function(name='jacc.bw.tests',
                                         cell.types = NULL) {
      data.all <- private$getResults(name, 'estimateDEStabilityBetweenTests()')

      if(is.null(cell.types)) cell.types = levels(self$cell.groups)
      data.all.median <- data.frame()

      unique.combinations <- data.all[!duplicated(data.all[,1:2]),1:2]
      for(i in 1:nrow(unique.combinations)) {
        idx = (data.all$s1 == unique.combinations$s1[i]) & (data.all$s2 == unique.combinations$s2[i]) &
          (data.all$cell.type %in% cell.types)
        data.tmp = data.frame(s1 = unique.combinations$s1[i],
                              s2 = unique.combinations$s2[i],
                              val = median(data.all$val[idx]))
        data.all.median = rbind(data.all.median, data.tmp)
      }

      p <- ggplot(data.all.median, aes(x=s1, y=s2, fill= val)) +
        geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylab(' ') + xlab(' ') + coord_equal() + scale_fill_gradient(low="white", high="seagreen4") +
        labs(fill='Jaccard index')

      return(p)
    },

    #' @description  Plot DE stability per cell type
    #' @param name - results slot name (default: 'de')
    #' @param show.pairs transparency value for the data points (default: 0.05)
    #' @param notch - whether to show notches in the boxplot version (default=TRUE)
    #' @return A ggplot2 object
    plotDEStabilityPerCellType=function(name='de.jaccards',
                                       notch = FALSE,
                                       show.jitter = TRUE,
                                       jitter.alpha = 0.05,
                                       show.pairs = FALSE,
                                       sort.order = TRUE) {

      jaccards <- private$getResults(name, 'estimateDEStability()')
      p <- plotStability(jaccards = jaccards,
                         notch = notch,
                         show.jitter = show.jitter,
                         jitter.alpha = jitter.alpha,
                         show.pairs = show.pairs,
                         sort.order = sort.order,
                         xlabel = 'Cell Type',
                         ylabel = 'Jaccard Index',
                         palette=self$cell.groups.palette,
                         plot.theme=self$plot.theme) + theme(legend.position = "none")
      p <- p + theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_color_manual(values=self$cell.groups.palette)
      return(p)
    },



    plotDEStabilityPerTest=function(name='jacc.per.test',
                                        notch = FALSE,
                                        show.jitter = TRUE,
                                        jitter.alpha = 0.5,
                                        show.pairs = FALSE,
                                        sort.order = FALSE) {

      jaccards <- private$getResults(name, 'estimateDEStability()')
      p <- plotStability(jaccards = jaccards,
                         notch = notch,
                         show.jitter = show.jitter,
                         jitter.alpha = jitter.alpha,
                         show.pairs = show.pairs,
                         sort.order = sort.order,
                         xlabel = 'Test name',
                         ylabel = 'Jaccard Index',
                         palette=self$cell.groups.palette,
                         plot.theme=self$plot.theme) + ylim(0, 1) +
        labs(color = "Cell types")
      return(p)
    },

    # plotNumberOfDEGenesStability=function(name = 'de',
    #                                       p.adjust = TRUE,
    #                                       pvalue.cutoff=0.05,
    #                                       notch = TRUE,
    #                                       show.jitter = TRUE,
    #                                       jitter.alpha = 0.1,
    #                                       show.pairs = FALSE,
    #                                       sort.order = FALSE,
    #                                       log.y.axis = FALSE){
    #   de.res <- private$getResults(name, 'estimateDEPerCellType()')
    #   if(!all(sapply(names(de.res), function(x) 'subsamples' %in% names(de.res[[x]])))) stop('Resampling was not performed')
    #
    #   de.numbers <- estimateNumberOfTermsStability(de.res,
    #                                                p.adjust = p.adjust,
    #                                                pvalue.cutoff = pvalue.cutoff)
    #
    #   p <- plotStability(jaccards = de.numbers,
    #                      notch = notch,
    #                      show.jitter = show.jitter,
    #                      jitter.alpha = jitter.alpha,
    #                      show.pairs = show.pairs,
    #                      sort.order = sort.order,
    #                      xlabel = 'Cell Type',
    #                      ylabel = 'Number of Significant DE genes',
    #                      log.y.axis = log.y.axis,
    #                      palette=self$cell.groups.palette,
    #                      plot.theme=self$plot.theme)
    #   return(p)
    # },


    plotDEStabilityPerGene = function(name = 'de',
                                      cell.type = NULL,
                                      stability.score = 'stab.median.rank'){
      de.res <- private$getResults(name, 'estimateDEPerCellType()')
      possible.scores = c('stab.median.rank', 'stab.mean.rank', 'stab.var.rank')
      if ( !(stability.score %in% possible.scores) ) stop('Please provide correct name of the stability core')
      if ( !(cell.type %in% names(de.res)) ) stop('Please provide correct cell type to visualise')
      if ( !(stability.score %in% names(de.res[[cell.type]]$res)) ) stop('Stability score was not calculated')

      idx = rank(de.res[[cell.type]]$res$pvalue) < 500
      p <- smoothScatter(x = rank(de.res[[cell.type]]$res$pvalue)[idx],
                         y = rank(de.res[[cell.type]]$res[[stability.score]])[idx],
                         xlab = 'rank of DE gene', ylab = 'stability score',
                         colramp = colorRampPalette(c("white", 'seagreen4')))
    },

    #' @description Plot number of significant DE genes
    #' @param name results slot in which the DE results should be stored (default: 'de')
    #' @param pvalue.cutoff P value cutoff (default=0.05)
    #' @param p.adjust whether the cutoff should be based on the adjusted P value (default: TRUE)
    #' @param show.resampling.results whether to show uncertainty based on resampling results (default: TRUE)
    #' @return A ggplot2 object
    plotNumberOfDEGenes=function(name='de', p.adjust=TRUE, pvalue.cutoff=0.05, show.resampling.results=TRUE,
                                 show.jitter=FALSE, jitter.alpha=0.05, type='bar', notch=TRUE, ...) {
      de.raw <- private$getResults(name, 'estimateDEPerCellType()')

      if (show.resampling.results) {
        subsamples <- lapply(de.raw, `[[`, 'subsamples')
        miss.subsamples <- names(subsamples)[sapply(subsamples, is.null)]
        if (length(miss.subsamples) == length(subsamples)) {
          warning("resampling results are missing for all cell types, falling back to point estimates.",
                  "Please rerun estimateDEPerCellType() with resampling != NULL")
          rl <- lapply(de.raw, `[[`, 'res')
        } else {
          warning("Subtypes ", paste(miss.subsamples, collapse=","), " are missed form sampling, ignoring those")
          rl <- unlist(subsamples, recursive=FALSE) %>%
            setNames(rep(names(de.raw), sapply(subsamples, length)))
        }
      } else {
        rl <- lapply(de.raw, `[[`, 'res')
      }

      # convert to dataframe for plotting
      p.col <- if (p.adjust) "padj" else "pvalue"
      df <- lapply(rl, function(d) data.frame(value=sum(d[[p.col]] <= pvalue.cutoff, na.rm=TRUE))) %>%
        bind_rows(.id="Type")

      plotMeanMedValuesPerCellType(
        df, show.jitter=show.jitter, jitter.alpha=jitter.alpha, notch=notch, type=type,
        palette=self$cell.groups.palette, ylab='number of DE genes', yline=NA, plot.theme=self$plot.theme, ...
      )
    },

    plotVolcano=function(name='de', cell.types=NULL, palette=NULL, build.panel=TRUE, n.col=3,
                         color.var = 'CellFrac', ...) {
      de <- private$getResults(name, 'estimateDEPerCellType()') %>% lapply(`[[`, 'res')
      if (is.null(palette)) {
        palette <- c("#5e4fa2", "#3288bd", "#abdda4", "#fdae61", "#f46d43", "#9e0142") %>%
          colorRampPalette()
      }

      if(!(color.var %in% names(de[[1]])))
        stop(paste(color.var, 'is not calculated'))

      if (is.null(cell.types))
        cell.types <- names(de)
      de <- de[intersect(cell.types, names(de))]

      if (length(de) == 0)
        stop("No cell types left after the filtering")

      if (length(cell.types) == 1)
        return(plotVolcano(de[[cell.types]], color.var=color.var, palette=palette, plot.theme=self$plot.theme, ...))

      if (!build.panel)
        return(lapply(de, plotVolcano, color.var=color.var, palette=palette, plot.theme=self$plot.theme, ...))

      gg <- lapply(de, plotVolcano, color.var=color.var, palette=palette,
                   xlab=NULL, ylab=NULL, plot.theme=self$plot.theme, ...) %>%
        cowplot::plot_grid(plotlist=., ncol=n.col, labels=paste0(names(.)), label_x=0.14,
                           label_y=0.99, label_size=10, align="hv", axis="lrtb", hjust=0)

      gg <- gg +
        theme(plot.margin=margin(b=12, l=12)) +
        draw_label("Log2(Fold Change)", size=12, y=-0.01, angle = 0) +
        draw_label("-Log10(P)", size=12, x=-0.01, angle = 90)

      return(gg)
    },

    #' @description Save DE results as JSON files
    #' @param saveprefix Prefix for created files (default=NULL)
    #' @param dir.name Name for directory with results (default="JSON")
    #' @param de.raw List of DE results
    #' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
    #' @param gene.metadata (default=NULL)
    saveDEasJSON=function(saveprefix=NULL, dir.name="JSON", de.raw=NULL, sample.groups=self$sample.groups, de.name='de',
                          ref.level=self$ref.level, gene.metadata=NULL, verbose=TRUE) {
      if (is.null(de.raw)) {
        de.raw <- private$getResults(de.name, "estimateDEPerCellType")
      }

      if (!is.list(sample.groups)) {
        sample.groups <- list(names(sample.groups[sample.groups == ref.level]),
                              names(sample.groups[sample.groups != ref.level])) %>%
          setNames(c(ref.level, self$target.level))
      }

      if (class(de.raw[[1]]) != "list") stop("Please rerun 'estimateDEPerCellType' with return.matrix=T")

      if (is.null(saveprefix)) saveprefix <- ""

      saveDEasJSON(de.raw=de.raw, saveprefix=saveprefix, dir.name=dir.name, gene.metadata=gene.metadata,
                   sample.groups=sample.groups, verbose=verbose)
    },

    ### Ontology analysis

    #' @description  Plot embedding
    #' @param embedding A cell embedding to use (two-column data frame with rownames corresponding to cells) (default: stored embedding object)
    #' @param plot.theme plot theme to use (default: `self$plot.theme`)
    #' @param color.by color cells by 'cell.groups', 'condition' or 'sample'. Overrides `groups` and `palette`. (default: NULL)
    #' @param ... other parameters passed to \link[sccore:embeddingPlot]{embeddingPlot}
    plotEmbedding=function(embedding=self$embedding, color.by=NULL, plot.theme=self$plot.theme, ...) {
      new.params <- list(...)
      params <- if (is.null(self$plot.params)) list() else self$plot.params
      params[names(new.params)] <- new.params

      if(is.null(embedding)) stop("embedding must be provided to Cacoa constructor or to this method.")
      private$checkCellEmbedding(embedding)
      if (!is.null(color.by)) {
        if (color.by == 'cell.groups') {
          params$groups <- self$cell.groups
          params$palette <- self$cell.groups.palette
        } else if (color.by == 'condition') {
          params$groups <- self$getConditionPerCell()
          params$palette <- self$sample.groups.palette
        } else if (color.by == 'sample') {
          params$groups <- self$sample.per.cell
        } else stop("Unknown color.by option: ", color.by)
      }

      params$embedding <- embedding
      params$plot.theme <- plot.theme
      if (is.null(params$show.legend)) {
        params$show.legend <- !is.null(params$colors)
      }

      rlang::invoke(sccore::embeddingPlot, params)
    },

    #' @description Estimate ontology terms based on DEs
    #' @param type Ontology type, either GO (gene ontology) or DO (disease ontology). Please see DOSE package for more information (default="GO")
    #' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
    #' @param p.adjust.method Method for calculating adj. P. Please see DOSE package for more information (default="BH")
    #' @param readable Mapping gene ID to gene name (default=TRUE)
    #' @param min.genes Minimum number of input genes overlapping with ontologies (default=0)
    #' @param qvalue.cutoff Q value cutoff, please see clusterProfiler package for more information (default=0.2)
    #' @param min.gs.size Minimal geneset size, please see clusterProfiler package for more information (default=5)
    #' @param max.gs.size Minimal geneset size, please see clusterProfiler package for more information (default=5e2)
    #' @return A list containing a list of terms per ontology, and a data frame with merged results
    estimateOntology=function(type="GO", name=NULL, de.name='de', org.db, n.top.genes=500, p.adj=1,
                              p.adjust.method="BH", readable=TRUE, verbose=TRUE, qvalue.cutoff=0.2, min.gs.size=10,
                              max.gs.size=5e2, keep.gene.sets=FALSE, ignore.cache=NULL, de.raw=NULL, ...) {
      if(!is.null(type) & !type %in% c("GO", "DO", "GSEA"))
        stop("'type' must be 'GO', 'DO', or 'GSEA'.")

      if(is.null(name)) {
        name <- type
      }
      if (is.null(de.raw)) {
        de.raw <- private$getResults(de.name, "estimateDEPerCellType()")
      }

      # If estimateDEPerCellType was run with return.matrix = TRUE, remove matrix before calculating
      if(class(de.raw[[1]]) == "list") de.raw %<>% lapply(`[[`, "res")

      de.gene.ids <- getDEEntrezIdsSplitted(de.raw, org.db=org.db, p.adj=p.adj)

      go.environment <- self$getGOEnvironment(org.db, verbose=verbose, ignore.cache=ignore.cache)
      res <- estimateOntologyFromIds(
        de.gene.ids, type=type, org.db=org.db, n.top.genes=n.top.genes, go.environment=go.environment,
        verbose=verbose, qvalue.cutoff=qvalue.cutoff, pAdjustMethod=p.adjust.method, readable=readable,
        minGSSize=min.gs.size, maxGSSize=max.gs.size, keep.gene.sets=keep.gene.sets, ...
      )

      self$test.results[[name]] <- list(res=res, de.gene.ids=de.gene.ids) # redundancy needed
      return(invisible(self$test.results[[name]]))
    },

    estimateOntologyGsea=function(name = NULL, de.name='de', org.db, p.adj=1,
                              p.adjust.method="BH",
                              readable=TRUE, verbose=TRUE,
                              qvalue.cutoff=0.2, min.gs.size=10,
                              max.gs.size=5e2, keep.gene.sets=FALSE,
                              ignore.cache=NULL, de.raw=NULL, ...) {

      if(is.null(name)) {
        name <- type
      }
      if (is.null(de.raw)) {
        de.raw <- private$getResults(de.name, "estimateDEPerCellType()")
      }

      # If estimateDEPerCellType was run with return.matrix = TRUE, remove matrix before calculating
      if(class(de.raw[[1]]) == "list") de.raw %<>% lapply(`[[`, "res")

      de.gene.ids <- getDEEntrezIdsSplitted(de.raw, org.db=org.db, p.adj=p.adj)

      go.environment <- self$getGOEnvironment(org.db, verbose=verbose, ignore.cache=ignore.cache)

      options(warn=-1)
      res <- list()
      for(cell.type in names(de.gene.ids)){
        print(cell.type)
        # stats - Named numeric vector with gene-level statistics sorted in decreasing order
        stats <- sort(de.gene.ids[[cell.type]]$all)
        for(go.type in c('BP', 'MF', 'CC')){
          # print(go.type)
          pathways <- go.environment[[go.type]]$PATHID2EXTID

          # All genes together
          # fgsea.res <- fgsea(pathways = pathways, stats = stats)
          suppressMessages(fgsea.res <- fgsea::fgsea(pathways = pathways, stats = stats, nperm = 10000))
          fgsea.res$pathw.name <- go.environment[[go.type]]$PATHID2NAME[fgsea.res$pathway]

          res.up <- fgsea.res[fgsea.res$ES > 0,]
          res.up$padj.new <- p.adjust(res.up$pval, method = 'fdr')
          res.down <- fgsea.res[fgsea.res$ES < 0,]
          res.down$padj.new <- p.adjust(res.down$pval, method = 'fdr')
          res[[cell.type]][[go.type]]$up <- res.up[order(res.up$pval),]
          res[[cell.type]][[go.type]]$down <- res.down[order(res.down$pval),]

          # # Up and down regulated genes separately
          # suppressMessages(fgsea.res.up <- fgsea::fgsea(pathways = pathways, stats = stats[stats>0], nperm = 10000))
          # suppressMessages(fgsea.res.down <- fgsea::fgsea(pathways = pathways, stats = stats[stats<0], nperm = 10000))
          #
          # res.up <- fgsea.res.up[fgsea.res.up$ES > 0,]
          # res.up$padj.new <- p.adjust(res.up$pval, method = 'fdr')
          # res.down <- fgsea.res.down[fgsea.res.down$ES < 0,]
          # res.down$padj.new <- p.adjust(res.down$pval, method = 'fdr')
          # res[[cell.type]][[go.type]]$up.sep <- res.up[order(res.up$pval),]
          # res[[cell.type]][[go.type]]$down.sep <- res.down[order(res.down$pval),]
        }
      }
      options(warn=0)

      self$test.results[[name]] <- list(res=res, de.gene.ids=de.gene.ids)
      return(invisible(self$test.results[[name]]))
    },

    #' @title Estimate ontology families
    #' @description Estimate ontology families based on ontology results
    #' @param type Type of ontology result, i.e., GO, GSEA, or DO (default: GO)
    #' @param p.adj Cut-off for adj. P values (default: 0.05)
    #' @return List of families and ontology data per cell type
    estimateOntologyFamilies=function(type = "GO", p.adj = 0.05, name = NULL) {
      checkPackageInstalled("GOfuncR", bioc=TRUE)
      if(!type %in% c("GO", "DO", "GSEA")) stop("'type' must be 'GO', 'DO', or 'GSEA'.")

      if(is.null(name)) {
        name <- type
      }

      # TODO: Test DO
      if (type == "GO") {
        ont.list <- self$test.results[[type]]$res %>%
          lapply(lapply, lapply, function(x) {
            tmp <- x@result %>% filter(p.adjust <= p.adj)
            if (nrow(tmp) > 0) return(tmp)
          }) %>%
          lapply(lapply, plyr::compact) %>%
          lapply(plyr::compact)
      } else {
        ont.list <- self$test.results[[type]]$res %>%
          lapply(lapply, function(x) {
            tmp <- x@result %>% filter(p.adjust <= p.adj)
            if (nrow(tmp) > 0) return(tmp)
          }) %>%
          lapply(plyr::compact)
      }
      self$test.results[[name]]$families <- estimateOntologyFamilies(ont.list=ont.list, type=type)
      return(invisible(self$test.results[[name]]))
    },

    #' @description Identify families containing a specific ontology term or ID
    #' @param go.term Character vector with term description(s)
    #' @param go.id Character vector with ID(s)
    #' @param common Boolean, only identify families with all the supplied terms or IDs (default = FALSE)
    #' @return Data frame
    getFamiliesPerGO=function(go.term=NULL, go.id=NULL, common = FALSE) {
      type <- "GO" # Maybe it would be nice to support other types
      if (is.null(go.term) && is.null(go.id)) stop("Please specify either 'go.term' or 'go.id'.")
      if (!is.null(go.term) && !is.null(go.id))
        warning("Both 'go.term' and 'go.id' specified, will only use 'go.term'.")

      # Extract data
      test.res <- self$test.results[[type]]
      ont.fam <- test.res$families
      if (is.null(ont.fam)) stop("No family data found.")

      # Get mapping of IDs to Descriptions
      desc.per.id <- rapply(test.res$res, function(x) x@result %$% setNames(Description, ID), how="list") %>%
        do.call(c, .) %>% do.call(c, .) %>% unname() %>% do.call(c, .)

      # Make data.frame
      res <- ont.fam %>% rblapply(c("CellType", "Ontology", "Genes"), function (fs) {
          fpg <- aggregate(f ~ value, rblapply(fs$families, "f", as_tibble), paste, collapse=",")
          mapply(function(go, f) {
              tibble(ID=c(go, fs$data[[go]]$parents_in_IDs), Families=strsplit(gsub("Family", "", f), ",", TRUE))
            }, fpg$value, fpg$f, SIMPLIFY=FALSE
          ) %>% bind_rows()
        }) %>%
        mutate(Description=desc.per.id[ID]) %>% filter((ID %in% go.id) | (Description %in% go.term))

      # Check and get result
      if (nrow(res) == 0) {
        warning("No families matching the query")
        return(res)
      }

      # Find common families
      if (common) {
        n.terms <- max(length(go.term), length(go.id))
        filt.fams <- res %$% split(Families, ID) %>%
          lapply(function(x) unique(unlist(x))) %>%
          unlist() %>% table() %>% {names(.)[. == n.terms]}
        res$Families %<>% lapply(intersect, filt.fams)
        res <- res[sapply(res$Families, length) > 0,]

        if (nrow(res) == 0) warning("No families having all terms present")
      }
      return(mutate(res, Families=sapply(Families, paste, collapse=", ")))
    },

    #' @description Estimate Gene Ontology clusters
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all'. Default: "up".
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types), "GO" or "DO". Default: "GO".
    #' @param name Name of the field to store the results. Default: `cacoa:::getOntClustField(subtype, genes)`.
    #' @param ind.h Cut height for hierarchical clustering of terms per cell type.
    #' Approximately equal to the fraction of genes, shared between the GOs. Default: 0.66.
    #' @param total.h Cut height for hierarchical clustering of GOs across all subtypes.
    #' Approximately equal to the fraction of subtypes for which two GOs should belong to the same cluster.
    #'   Default: 0.5.
    #' @return List containing:
    #'   - `df`: data.frame with information about individual gene ontolodies and columns `Cluster` and `ClusterName`
    #'     for the clustering info
    #'   - `hclust`: the object of class \link[stats:hclust]{hclust} with hierarchical clustering of GOs across all
    #'     subtypes
    estimateOntologyClusters=function(type="GO", subtype=NULL, genes="all", ind.h=0.66, total.h=0.5, p.adj=0.05,
                                      min.genes=1, name=getOntClustField(subtype, genes), clust.naming="medoid",
                                      verbose=self$verbose) {
      ont.df <- private$getOntologyPvalueResults(genes, type=type, subtype=subtype, p.adj=p.adj, min.genes=min.genes)
      if (nrow(ont.df) == 0) {
        res <- list(df=ont.df)
        self$test.results[[type]][[name]] <- res
        return(invisible(res))
      }

      genes.per.go.per.type <- ont.df$geneID %>% strsplit("/") %>%
        setNames(ont.df$Description) %>% split(ont.df$Group)

      clust.df <- clusterIndividualGOs(genes.per.go.per.type, cut.h=ind.h)
      clusts <- clusterGOsPerType(clust.df, cut.h=total.h, verbose=verbose)

      ont.df$Cluster <- clusts$clusts[as.character(ont.df$Description)]

      name.per.clust <- estimateOntologyClusterNames(ont.df, clust.naming=clust.naming)

      ont.df$ClusterName <- name.per.clust[ont.df$Cluster]
      res <- list(df=ont.df, hclust=clusts$hclust)
      self$test.results[[type]][[name]] <- res

      return(invisible(res))
    },

    #' @description Bar plot of ontology terms per cell type
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "GO" or "DO" (default="GO")
    #' @return A ggplot2 object
    plotNumOntologyTermsPerType=function(genes="up", type="GO", p.adj=0.05, min.genes=1,
                                         cell.groups=self$cell.groups, name=NULL) {
      if (length(genes) > 0) {
        ont.res <- genes %>% setNames(., .) %>%
          lapply(private$getOntologyPvalueResults, type=type, p.adj=p.adj, min.genes=min.genes)

        if (type != "GSEA") {
          classes <- sapply(ont.res[genes], class)
          if(any(classes == "character")) {
            message(paste0("No significant results found for genes = '", names(classes[classes == "character"]), "'."))
            genes <- names(classes[classes == "data.frame"])
            if(length(genes) == 0) stop("No results to plot.")
          }

          ont.res %<>% .[genes]
          ont.res %<>% names() %>%
            lapply(function(d) ont.res[[d]] %>% dplyr::mutate(direction = d))
        }

        ont.res %<>% Reduce(rbind, .)
      } else {
        ont.res <- private$getOntologyPvalueResults(genes=genes, type=type, p.adj=p.adj, min.genes=min.genes)
        if (type != "GSEA") {
          ont.res %<>% dplyr::mutate(direction = genes)
        }
      }

      # Prepare data further
      if(type=="GO") {
        p.df <- table(ont.res$Group, ont.res$Type, ont.res$direction) %>%
          as.data.frame() %>%
          setNames(c("Group","Type","direction","N")) %>%
          dplyr::arrange(Group)

        if(length(unique(p.df$direction)) > 1) {
          gg <- ggplot(p.df, aes(x=Group, y=N, fill=Type, group=Group)) +
            geom_bar(stat="identity") +
            facet_grid(~direction, switch="x")
        } else {
          gg <- ggplot(p.df) +
            geom_bar(aes(x=Group, y=N, fill=Type), stat="identity")
        }
      } else if(type=="DO") {
        p.df <- table(ont.res$Group, ont.res$direction) %>%
          as.data.frame() %>%
          setNames(c("Group","direction","N")) %>%
          dplyr::arrange(Group)

        if(length(unique(p.df$direction)) > 1) {
          gg <- ggplot(p.df) +
            geom_bar(aes(x=Group, y=N, fill=direction), stat="identity", position="dodge") +
            labs(fill="Gene set")
        } else {
          gg <- ggplot(p.df) +
            geom_bar(aes(x=Group, y=N), stat="identity")
        }
      } else if(type == "GSEA") {
        p.df <- table(ont.res$Group, ont.res$Type) %>%
          as.data.frame() %>%
          setNames(c("Group","Type","N")) %>%
          dplyr::arrange(Group)

          gg <- ggplot(p.df) +
            geom_bar(aes(x=Group, y=N, fill=Type), stat="identity")
      }

      gg <- gg +
        scale_y_continuous(expand=c(0, 0)) +
        self$plot.theme +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
              legend.position="right") +
        labs(x="", y=paste0("No. of ", type, " terms"))
      return(gg)
    },

    #' @description Plot a dotplot of ontology terms with adj. P values for a specific cell subgroup
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types), "GO" or "DO" (default="GO")
    #' @param cell.subgroup Cell group to plot
    #' @param n Number of ontology terms to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
    #' @param p.adj Adjusted P cutoff (default=0.05)
    #' @param log.colors Use log10 p-values for coloring (default=FALSE)
    #' @return A ggplot2 object
    plotOntology = function(cell.subgroups, plot="dot", genes="up", type="GO", subtype="BP", n=20, p.adj=0.05,
                            min.genes=1, ...) {
      # Checks
      checkPackageInstalled("enrichplot", bioc=TRUE)
      if (is.null(type) || (!type %in% c("GO","DO","GSEA"))) stop("'type' must be 'GO', 'DO', or 'GSEA'.")
      if (is.null(subtype) || (!subtype %in% c("BP","CC","MF"))) stop("'subtype' must be 'BP', 'CC', or 'MF'.")
      if (is.null(genes) || (!genes %in% c("down","up","all"))) stop("'genes' must be 'down', 'up', or 'all'.")
      if (is.null(cell.subgroups)) stop("Please define which cells to plot using 'cell.subgroups'.")
      if (type == "GSEA" & plot == "bar") stop("No 'enrichplot' method exists for making barplots of GSEA results.")

      # Extract results
      ont.res <- self$test.results[[type]]$res
      if (is.null(ont.res)) stop("No results found for ", type)
      if (!cell.subgroups %in% names(ont.res)) stop("'cell.subgroups' not found in results.")
      ont.res %<>% .[[cell.subgroups]]
      if (type != "DO") ont.res %<>% .[[subtype]]
      if (is.null(ont.res))
        stop("No results found for ", type, ", ", subtype, " for ", cell.subgroups)

      if (type %in% c("GO","DO")) ont.res %<>% .[[genes]]
      if (is.null(ont.res))
        stop("No results found for ", genes, " genes for ", type, ", ", subtype, " for ", cell.subgroups)

      # Prepare data
      df <- ont.res@result %>% filter(p.adjust <= p.adj)
      if (nrow(df) == 0)
        stop("Nothing to plot. Try relaxing 'p.adj'. The lowest adj. P value is ",
             formatC(min(ont.res@result$p.adjust), digits=3))

      # Allow plotting of terms with p.adj > 0.05
      if (p.adj != 0.05) {
        ont.res@pvalueCutoff <- 1
        ont.res@qvalueCutoff <- 1
      }

      if (min.genes > 1) {
        idx <- df$GeneRatio %>%
          strsplit("/", fixed=TRUE) %>%
          sapply(`[[`, 1)

        df <- df[idx > min.genes,]
      }
      ont.res@result <- df

      # Plot
      if (plot == "dot")
        return(dotplot(ont.res, showCategory=n, orderBy="x", ...))

      if (plot == "bar")
        return(barplot(ont.res, showCategory=n, ...))

      stop("Unknown plot type: ", plot)
    },

    #' @description Plot a heatmap of ontology P values per cell type
    #' @param genes Specify which genes to plot, can either be 'down' for downregulated genes, 'up' or 'all'
    #'   (default="up")
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default="GO")
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="left")
    #' @param selection Order of rows in heatmap. Can be 'unique' (only show terms that are unique for any cell type);
    #'   'common' (only show terms that are present in at least two cell types); 'all' (all ontology terms)
    #'   (default="all")
    #' @param top.n Number of terms to show (default=Inf)
    #' @param clusters Whether to show GO clusters or raw GOs (default=TRUE)
    #' @param cluster.name Field with the results for GO clustering. Ignored if `clusters == FALSE`.
    #' @param cell.subgroups Cell groups to plot (default=NULL). This affects only visualization, but not clustering.
    #' @param color.range vector with two values for min/max values of p-values
    #' @param ... parameters forwarded to \link{plotHeatmap}
    #' @return A ggplot2 object
    plotOntologyHeatmap=function(genes="up", type="GO", subtype="BP", min.genes=1, p.adj=0.05, legend.position="left",
                                 selection="all", top.n=Inf, clusters=TRUE, cluster.name=NULL,
                                 cell.subgroups=NULL, palette=NULL, row.order=TRUE, col.order=TRUE, max.log.p=10,
                                 only.family.children=FALSE, description.regex=NULL, clust.naming="medoid", ...) {
      ont.sum <- private$getOntologyHeatmapInfo(
        genes=genes, type=type, subtype=subtype, min.genes=min.genes, p.adj=p.adj, selection=selection,
        clusters=clusters, cluster.name=cluster.name, cell.subgroups=cell.subgroups,
        only.family.children=only.family.children, description.regex=description.regex, clust.naming=clust.naming
      )

      if (is.null(ont.sum)) return(ggplot())
      if (is.null(palette)) palette <- getGenePalette(genes, high="white")

      ont.sum[ont.sum > max.log.p] <- max.log.p
      ont.sum %<>% .[order(rowSums(.), decreasing=TRUE),,drop=FALSE] %>% head(top.n)

      plt <- plotHeatmap(
        ont.sum, legend.position=legend.position, row.order=row.order, col.order=col.order,
        plot.theme=self$plot.theme, palette=palette, ...
      )
      # plt <- as.matrix(ont.sum) %>%
      #   ComplexHeatmap::Heatmap(border=TRUE, show_row_dend=FALSE, show_column_dend=FALSE,
      #                           row_names_max_width=unit(8, "cm"), row_names_gp=grid::gpar(fontsize=10),
      #                           cluster_rows=row.order, cluster_columns=col.order,
      #                           col=palette(100), ...)
      return(plt)
    },

    plotOntologyHeatmapCollapsed=function(genes="up", type="GO", subtype="BP", min.genes=1, p.adj=0.05,
                                          legend.position="left", selection="all", n=20, clusters=TRUE,
                                          cluster.name=NULL, cell.subgroups=NULL, palette=NULL, row.order=TRUE,
                                          col.order=TRUE, max.log.p=10, only.family.children=FALSE,
                                          description.regex=NULL, distance="manhattan", clust.method="complete",
                                          clust.naming="consensus", n.words=5, exclude.words=NULL, ...) {
      ont.sum <- private$getOntologyHeatmapInfo(
        genes=genes, type=type, subtype=subtype, min.genes=min.genes, p.adj=p.adj, selection=selection,
        clusters=clusters, cluster.name=cluster.name, cell.subgroups=cell.subgroups,
        only.family.children=only.family.children, description.regex=description.regex, clust.naming=clust.naming,
        return.descriptions=TRUE
      )

      desc.per.clust <- ont.sum$desc.per.clust
      ont.sum <- ont.sum$ont.sum

      if (is.null(ont.sum)) return(ggplot())
      if (is.null(palette)) palette <- getGenePalette(genes, high="white")

      ont.sum.raw <- ont.sum
      ont.sum[ont.sum > max.log.p] <- max.log.p

      gos.per.clust <- dist(ont.sum, method=distance) %>%
        hclust(method=clust.method) %>% cutree(k=n) %>% {split(names(.), .)}

      clust.names <- sapply(gos.per.clust, function(gos) {
        n.gos <- length(gos)
        if (n.gos == 1) return(paste("1: ", gos))
        if (clusters) gos <- unlist(desc.per.clust[gos], use.names=FALSE)

        estimateOntologyClusterName(gos, method=clust.naming, n.words=n.words, exclude.words=exclude.words) %>%
          {paste0(n.gos, ": ", .)}
      })

      ont.sum <- lapply(gos.per.clust, function(ns) colMeans(ont.sum.raw[ns,,drop=FALSE])) %>%
        do.call(rbind, .) %>% as.data.frame() %>% set_rownames(clust.names)
      ont.sum[ont.sum > max.log.p] <- max.log.p

      ont.freqs <- gos.per.clust %>%
        lapply(function(gos) colMeans(ont.sum.raw[gos,] > 1e-5) * 100) %>%
        do.call(rbind, .) %>% as.data.frame()

      plt <- plotHeatmap(
        ont.sum, size.df=ont.freqs, legend.position=legend.position, row.order=row.order, col.order=col.order,
        plot.theme=self$plot.theme, palette=palette, distance=distance, clust.method=clust.method,
        size.legend.title="Cluster GOs %", ...
      )
      return(plt)
    },

    #' @description Plot correlation matrix for ontology terms between cell types
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "GO" or "DO" (default="GO")
    #' @return A ggplot2 object
    plotOntologySimilarities=function(genes = "up", type = "GO", p.adj = 0.05, min.genes = 1) {
      ont.res <- private$getOntologyPvalueResults(genes=genes, type=type, p.adj=p.adj, min.genes=min.genes)

      if((ont.res$Group %>% unique() %>% length()) == 1)
        stop("Only one group present, correlation cannot be performed.")

      if(nrow(ont.res) == 0)
        stop("No significant ontology terms identified. Try relaxing p.adj.")

      if(type %in% c("GO", "GSEA")) {
        pathway.df <- ont.res[c("Description", "Group", "Type")] %>% rename(Pathway=Description, GO=Type)
      } else if(type=="DO") {
        pathway.df <- unique(ont.res$Group) %>%
          lapply(function(cell.group) {
            tibble::tibble(Pathway=ont.res %>%
                             dplyr::filter(Group == cell.group) %>%
                             dplyr::pull(Description),
                           Group=cell.group)
          }) %>%
          dplyr::bind_rows()
      }

      path.bin <- pathway.df %>%
        dplyr::select(Pathway, Group) %>%
        dplyr::mutate(X=1) %>%
        tidyr::spread(Pathway, X) %>%
        as.data.frame() %>%
        magrittr::set_rownames(.$Group) %>%
        .[, 2:ncol(.)] %>%
        as.matrix()
      path.bin[is.na(path.bin)] <- 0

      # TODO: currently we use binary distance. Probably, checking z-scores would give bitter results.
      p.mat <- (1 - (path.bin %>% dist(method="binary") %>% as.matrix)) %>% pmin(0.5)
      cl.tree <- dist(p.mat) %>% hclust()
      clust.order <- cl.tree %$% labels[order]
      clusts <- cutree(cl.tree, h=0.7)[clust.order]
      clusts[clusts %in% names(which(table(clusts) < 5))] <- max(clusts) + 1
      clust.lengths <- rle(clusts)$lengths %>% rev
      diag(p.mat) <- 1

      # Plot
      plotHeatmap(p.mat, color.per.group=NULL, row.order=clust.order, col.order=rev(clust.order),
                  legend.title="Similarity", plot.theme=self$plot.theme) +
        scale_fill_distiller(palette="RdYlBu", limits=c(0, 0.5)) +
        geom_vline(aes(xintercept=x), data.frame(x=cumsum(clust.lengths)[clust.lengths > 1] + 0.5)) +
        geom_hline(aes(yintercept=x), data.frame(x=cumsum(clust.lengths)[clust.lengths > 1] + 0.5))
    },

    #' @title Plot ontology families
    #' @description Plot related ontologies in one hierarchical network plot
    #' @param type Type of ontology result, i.e., GO, GSEA, or DO (default: GO)
    #' @param cell.subgroups Cell subtype to plot
    #' @param family Family within cell subtype to plot (numeric value)
    #' @param genes Only for GO results: Direction of genes, must be "up", "down", or "all" (default: up)
    #' @param subtype Only for GO results: Type of result, must be "BP", "MF", or "CC" (default: BP)
    #' @param plot.type How much of the family network should be plotted. Can be "complete" (entire network), "dense" (show 1 parent for each significant term), or "minimal" (only show significant terms) (default: complete)
    #' @param show.ids Whether to show ontology IDs instead of names (default: FALSE)
    #' @param string.length Length of strings for wrapping in order to fit text within boxes (default: 14)
    #' @param legend.label.size Size og legend labels (default: 1)
    #' @param legend.position Position of legend (default: topright)
    #' @param verbose Print messages (default: stored value)
    #' @param n.cores Number of cores to use (default: stored value)
    #' @return Rgraphviz object
    plotOntologyFamily=function(type = "GO", cell.subgroups, family=NULL, genes="up", subtype="BP",
                                plot.type="complete", show.ids=FALSE, string.length=14, legend.label.size=1,
                                legend.position="topright", verbose=self$verbose, n.cores=self$n.cores, ...) {
      #Checks
      checkPackageInstalled(c("GOfuncR", "graph", "Rgraphviz"), bioc=TRUE)

      ont.fam.res <- self$test.results[[type]]$families
      if(is.null(ont.fam.res))
        stop("No results found for type '", type, "'. Please run 'estimateOntologyFamilies' first.")

      ont.fam.res %<>% .[[cell.subgroups]]
      if (is.null(ont.fam.res)) stop("No results found for cell.subgroups '", cell.subgroups, "'.")
      # TODO: Test for GSEA/GO. Update description!
      ont.fam.res %<>% .[[subtype]]
      if (is.null(ont.fam.res)) stop("No results found for subtype '", subtype, "'.")

      if (type != "GSEA") {
        ont.fam.res %<>% .[[genes]]
        if (is.null(ont.fam.res)) stop("No results found for genes '", genes, "'.")
      }

      if (is.null(family)) {
        fam.names <- names(ont.fam.res$families)
      } else {
        if (!is.numeric(family)) {
          fam.names <- family
        } else {
          fam.names <- paste0("Family",family)
        }
      }

      if(!all(fam.names %in% names(ont.fam.res$families))) stop("Not all families are found in 'ont.fam.res'.")
      families <- lapply(fam.names, function(n) ont.fam.res$families[[n]]) %>% unlist() %>% unique()

      plotOntologyFamily(fam=families, data=ont.fam.res$data, plot.type=plot.type, show.ids=show.ids,
                         string.length=string.length, legend.label.size=legend.label.size, legend.position=legend.position,
                         verbose=verbose, n.cores=n.cores, ...)
    },

    #' @title Save ontology results
    #' @description Save ontology results as a table
    #' @param file File name. Set to NULL to return the table instead of saving
    #' @param type Type of ontology result, i.e., GO, GSEA, or DO (default: GO)
    #' @param subtype Only for GO results: Type of result to filter by, must be "BP", "MF", or "CC" (default: NULL)
    #' @param genes Only for GO results: Direction of genes to filter by, must be "up", "down", or "all" (default: NULL)
    #' @param p.adj Adjusted P to filter by (default: 0.05)
    #' @param sep Separator (default: tab)
    #' @return Table for import into text editor
    saveOntologyAsTable=function(file, type="GO", subtype=NULL, genes=NULL, p.adj=0.05, sep="\t", ...) {
      res <- self$test.results[[type]]$res %>%
        rblapply(c("CellType", "Subtype", "Genes"), function (r) r@result) %>%
        filter(p.adjust <= p.adj) %>% as_tibble()

      if (!is.null(subtype)) res %<>% filter(Subtype %in% subtype)
      if (!is.null(genes)) res %<>% filter(Genes %in% genes)

      if (is.null(file)) return(res)
      write.table(res, file=file, sep=sep, row.names=FALSE, ...)
    },

    #' @title Save family results
    #' @description Save family results as a table
    #' @param file File name (default: Family_results.tsv)
    #' @param type Type of ontology result, i.e., GO, GSEA, or DO (default: GO)
    #' @param subtype Only for GO results: Type of result to filter by, must be "BP", "MF", or "CC" (default: NULL)
    #' @param genes Only for GO results: Direction of genes to filter by, must be "up", "down", or "all" (default: NULL)
    #' @param p.adj Adjusted P to filter by (default: 0.05)
    #' @param sep Separator (default: tab)
    #' @return Table for import into text editor
    saveFamiliesAsTable=function(file, type="GO", subtype=NULL, genes=NULL, p.adj=0.05, sep="\t", ...) {
      # Extract results
      tmp <- lapply(self$test.results[[type]]$families, lapply, lapply, function(x) {
        lapply(x$families, function(y) {
          if (length(y) <= 1) return(NULL) # What about == 1?
          tmp.data <- x$data[y] %>%
            lapply(lapply, function(z) if (length(z) > 1) paste0(z, collapse="/") else z) %>%
            lapply(function(z) z[c("Description", "Significance", "parents_in_IDs", "parent_go_id")]) %>%
            bind_rows() %>%
            data.frame() %>%
            mutate(No_parents=sapply(.$parents_in_IDs, function(z) length(strsplit(z, split="/", fixed=TRUE)[[1]])),
                    Child_terms = nrow(.))

          data.frame(y, tmp.data) %>%
            setNames(c("Lonely_child_IDs", "Description", "P.adj", "Significant_parents", "All_parents",
                        "#_sig_parents", "#_child_terms"))
        }) %>% .[!sapply(., is.null)]
      })

      # Convert to data frame
      res <- rblapply(tmp, c("CellType", "Subtype", "Genes", "Family"), function (r) r)

      # Filtering
      if (!is.null(subtype)) res %<>% filter(subtype %in% subtype)
      if (!is.null(genes)) res %<>% filter(genes %in% genes)

      # Write table
      if (is.null(file)) return(res)
      write.table(res, file=file, sep=sep, row.names=FALSE, ...)
    },

    #' @description Plot the cell group sizes or proportions per sample
    #' @param palette color palette to use for conditions (default: stored $sample.groups.palette)
    #' @param show.significance whether to show statistical significance betwwen sample groups. wilcox.test was used; (`*` < 0.05; `**` < 0.01; `***` < 0.001)
    #' @param proportions whether to plot proportions or absolute numbers (default: TRUE)
    #' @param ... additional plot parameters, forwarded to \link{plotCountBoxplotsPerType}
    #' @return A ggplot2 object
    plotCellGroupSizes=function(cell.groups=self$cell.groups, show.significance=FALSE, filter.empty.cell.types=TRUE,
                                proportions=TRUE, palette=self$sample.groups.palette, ...) {
      df.melt <- private$extractCodaData(cell.groups=cell.groups, ret.groups=FALSE)

      if (proportions) {
        df.melt %<>% {100 * . / rowSums(.)}
        y.lab <- "% cells per sample"
      } else {
        y.lab <- "Num. cells per sample"
      }

      df.melt %<>% as.data.frame() %>%
        dplyr::mutate(group=self$sample.groups[levels(self$sample.per.cell)]) %>%
        reshape2::melt(id.vars="group")

      # Filtration
      if (filter.empty.cell.types) {
        cell.types.counts <- table(df.melt$variable[df.melt$value>0], df.melt$group[df.melt$value>0])
        cell.types.to.remain <- cell.types.counts %>% {rownames(.)[rowSums(. == 0) == 0]}
        df.melt <- df.melt[df.melt$variable %in% cell.types.to.remain,]

        for (tmp.level in colnames(cell.types.counts)) {
          cell.types.tmp <- cell.types.counts %>% {rownames(.)[(rowSums(. == 0) != 0) & (.[, tmp.level] > 0)]}
          if (length(cell.types.tmp) > 0)
            message('Cell types {', paste(cell.types.tmp, collapse=", "), '} are present only in ', tmp.level, ' samples')
        }
      }

      gg <- plotCountBoxplotsPerType(df.melt, y.lab=y.lab, palette=palette, show.significance=show.significance,
                                     plot.theme=self$plot.theme, ...)

      return(gg)
    },

    plotCellGroupAbundanceVariation=function(cell.groups=self$cell.groups, type='mad', rotate.xticks=TRUE, min.rel.abundance=0.05) {
      n.cells.per.samp <- table(self$sample.per.cell)
      vars.per.group <- cell.groups %>%
        {split(names(.), .)} %>% sapply(function(ids) {
          sample.fracs <- table(self$sample.per.cell[ids]) / n.cells.per.samp / length(ids)
          if (type == 'mad') return(mad(sample.fracs))
          if (type == 'sd') return(sd(sample.fracs))

          is.missed <- sample.fracs %>% {. / mean(.) < min.rel.abundance}
          if (type == 'sample.num') return(sum(is.missed))
          if (type == 'sample.frac') return(mean(is.missed))
          stop("Unknown type: ", type)
        }) %>% {data.frame(var=., type=factor(names(.), levels=names(.)[order(.)]))}

      gg <- ggplot(vars.per.group) +
        geom_bar(aes(x=type, y=var), stat="identity") +
        self$plot.theme +
        scale_y_continuous(expand=c(0, 0), limits=c(0, max(vars.per.group$var) * 1.05), name=type) +
        theme(axis.title.x=element_blank())

      if (rotate.xticks) {
        gg <- gg + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
      }

      return(gg)
    },

    #' @description Plot compositions in CoDA space (PCA or CDA)
    #' @param space either 'PCA' or 'CDA'
    #' @return A ggplot2 object
    plotCodaSpace=function(space='CDA', cell.groups=self$cell.groups, font.size=3, cells.to.remain=NULL, cells.to.remove=NULL,
                           samples.to.remove=NULL, palette=self$sample.groups.palette) {
      tmp <- private$extractCodaData(cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain, samples.to.remove=samples.to.remove, cell.groups=cell.groups)

      if (space == 'PCA') {
        bal <- getRndBalances(tmp$d.counts)
        pca.res <- prcomp(bal$norm)
        pca.loadings <- bal$psi %*% pca.res$rotation

        dfs <- list(
          red=as.data.frame(pca.res$x) %>% set_colnames(c("S1", "S2")),
          loadings=(10 * as.data.frame(pca.loadings[,1:2])) %>% set_colnames(c("S1", "S2"))
        )
        gg.labs <- labs(x="PC1", y="PC2")
      } else if (space == 'CDA') {
        dfs <- estimateCdaSpace(tmp$d.counts, tmp$d.groups)
        gg.labs <- labs(x="Score 1", y="Score 2")
      }

      gg <- plotCodaSpaceInner(dfs$red, dfs$loadings, d.groups=tmp$d.groups, ref.level=self$ref.level, target.level=self$target.level,
                               font.size=font.size, palette=palette)
      return(gg + gg.labs + self$plot.theme)
    },

    #' @description Plot contrast tree
    #' @return A ggplot2 object
    plotContrastTree=function(cell.groups=self$cell.groups, palette=self$sample.groups.palette,
                              cells.to.remain = NULL, cells.to.remove = NULL, filter.empty.cell.types = TRUE,
                              adjust.pvalues = TRUE, h.method=c('both', 'up', 'down'),
                              reorder.tree = TRUE, ...) {
      h.method <- match.arg(h.method)
      tmp <- private$extractCodaData(cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain,
                                     cell.groups=cell.groups)
      if(filter.empty.cell.types) {
        cell.type.to.remain <- (colSums(tmp$d.counts[tmp$d.groups,]) > 0) &
          (colSums(tmp$d.counts[!tmp$d.groups,]) > 0)
        tmp$d.counts <- tmp$d.counts[,cell.type.to.remain]
      }

      tree.order = NULL
      pval.cell.types = NULL
      if(reorder.tree){
        if ("coda" %in% names(self$test.results)){
          loadings.mean <- rowMeans(self$test.results$coda$loadings) - self$test.results$coda$ref.load.level
          pval <- self$test.results$coda$pval
          loadings.mean <- loadings.mean[order(pval)]
          tree.order <- c(names(loadings.mean)[loadings.mean < 0],
                         rev(names(loadings.mean)[loadings.mean > 0]))
          pval.cell.types <-self$test.results$coda$padj
        }
      }

      gg <- plotContrastTree(tmp$d.counts, tmp$d.groups, self$ref.level, self$target.level,
                             plot.theme=self$plot.theme, adjust.pvalues=adjust.pvalues,
                             h.method=h.method, tree.order=tree.order,
                             pval.cell.types=pval.cell.types,...)

      if (!is.null(palette)) {
        gg <- gg + scale_color_manual(values=palette)
      }
      return(gg)
    },


    #' @description Plot composition similarity
    #' @return A ggplot2 object
    plotCompositionSimilarity=function(cell.groups=self$cell.groups, cells.to.remain=NULL, cells.to.remove=NULL) {
      tmp <- private$extractCodaData(cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain, cell.groups=cell.groups)
      res <- referenceSet(tmp$d.counts, tmp$d.groups)
      heatmap(res$mx.first, scale = 'none')
    },



    estimateCellLoadings=function(n.perm=1000, n.boot=1000, coda.test='significance', # TODO: n.perm, coda.test, criteria and define.ref.cell.type are never used
                                  ref.cell.type=NULL, criteria='cda.std', name='coda',
                                  n.seed=239, cells.to.remove=NULL, cells.to.remain=NULL,
                                  samples.to.remove=NULL, filter.empty.cell.types=TRUE,
                                  define.ref.cell.type=FALSE, n.cores=self$n.cores, verbose=self$verbose) {

      if (!(coda.test %in% c('significance', 'confidence'))) stop('Test is not supported')
      checkPackageInstalled(c("coda.base", "psych"), cran=TRUE)

      if ((!is.null(ref.cell.type)) && (!(ref.cell.type %in% levels(self$cell.groups))))
        stop('Incorrect reference cell type')
      if (define.ref.cell.type & (!is.null(ref.cell.type))) define.ref.cell.type <- FALSE

      # Get cell counts and groups
      tmp <- private$extractCodaData(cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain, samples.to.remove=samples.to.remove)
      # self$test.results[['tmp']] <- tmp

      if(filter.empty.cell.types) {
        cell.type.to.remain <- (colSums(tmp$d.counts[tmp$d.groups,]) > 0) &
          (colSums(tmp$d.counts[!tmp$d.groups,]) > 0)
        tmp$d.counts <- tmp$d.counts[,cell.type.to.remain]
      }
      cnts <- tmp$d.counts
      groups <- tmp$d.groups

      # if(coda.test == 'confidence'){
      #   n.perm <- 1
      # }

      res <- runCoda(cnts, groups, n.boot=n.boot, n.seed=n.seed, ref.cell.type=ref.cell.type)
      loadings.init <- res$loadings.init
      padj <- res$padj
      pval <- res$pval
      ref.load.level <- res$ref.load.level

      ## Calculate normalized counts
      ref.cell.type <- res$ref.cell.type

      ref.cnts <- cnts[, ref.cell.type, drop=FALSE]
      ref.cnts[ref.cnts == 0] <- 0.5
      norm.val <- 1 / nrow(ref.cnts) * rowSums(log(ref.cnts))
      cnts.nonzero <- cnts
      cnts.nonzero[cnts.nonzero == 0] <- 0.5
      norm.cnts <- log(cnts.nonzero) - norm.val

      self$test.results[[name]] <- list(
        loadings=loadings.init,
        pval=pval,
        padj = padj,
        ref.load.level = ref.load.level,
        cell.list=res$cell.list,
        cnts=cnts,
        groups=groups,
        ref.cell.types=ref.cell.type,
        norm.cnts=norm.cnts
      )

      return(invisible(self$test.results[[name]]))
    },

    estimateGaPartition=function(cells.to.remain=NULL, cells.to.remove=NULL, samples.to.remove=NULL, ...){ # TODO: do we ever use this?
      tmp <- private$extractCodaData(cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain, samples.to.remove=samples.to.remove)
      ga.res <- gaPartition(tmp$d.counts, tmp$d.groups, ...)

      self$test.results[['ga.partition']] <- rownames(t(ga.res[1,ga.res[1,] != 0,drop=FALSE]))

      return(invisible(self$test.results[['ga.partition']]))
    },

    #' @description Plot Loadings
    #' @param palette palette specification for cell types (default: stored $cell.groups.palette)
    #' @return A ggplot2 object
    plotCellLoadings=function(alpha=0.01, palette=self$cell.groups.palette, font.size=NULL, name='coda',
                              ordering='pvalue', show.pvals=TRUE, ...) {

      loadings <- private$getResults(name, 'estimateCellLoadings()')
      p <- loadings %$% plotCellLoadings(
        loadings, pval=padj, jitter.alpha=alpha, palette=palette, show.pvals=show.pvals,
        ref.level=self$ref.level, target.level=self$target.level, plot.theme=self$plot.theme,
        ref.load.level=ref.load.level, ordering=ordering, ...
      )

      return(p)
    },

    estimateWilcoxonTest = function(cell.groups=self$cell.groups, cells.to.remove = NULL, cells.to.remain = NULL) { # TODO: do we ever use this?
      tmp <- private$extractCodaData(cell.groups=cell.groups, cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain)
      p.vals <- calcWilcoxonTest(tmp$d.counts, tmp$d.groups)
      self$test.results[['p.vals.balances']] <- p.vals
      return(invisible(p.vals))
    },


    ### Segmentation-free cell density

    #' @description Estimate cell density in giving embedding
    #' @param bins number of bins for density estimation, default 400
    #' @param method density estimation method, graph: graph smooth based density estimation. kde: embedding grid based density  estimation. (default: 'kde')
    #' @param beta smoothing strength parameter of the \link[sccore:heatFilter]{heatFilter} for graph based cell density (default: 30)
    #' @param bandwidth KDE bandwidth multiplier. The full bandwidth is estimated by multiplying this value on the difference between 90% and 10%
    #' of the corresponding embedding dimension. Set it to NULL to use \link[MASS:bandwidth.nrd]{bandwidth.nrd} estimator. (default: 0.05)
    #' @param name slot in which to save the results (default: 'cell.density')
    estimateCellDensity = function(bins=400, method='kde', name='cell.density', beta=30, estimate.variation=TRUE,
                                   sample.groups=self$sample.groups, verbose=self$verbose, n.cores=self$n.cores,
                                   bandwidth=0.05, ...){
      sample.per.cell <- self$sample.per.cell

      if (method == 'kde') {
        private$checkCellEmbedding()
        res <- self$embedding %>% estimateCellDensityKde(
          sample.per.cell=sample.per.cell, sample.groups=sample.groups, bins=bins, bandwidth=bandwidth, ...
        )
      } else if (method == 'graph') {
        res <-  extractCellGraph(self$data.object) %>%
          estimateCellDensityGraph(sample.per.cell=sample.per.cell, sample.groups=sample.groups,
                                   n.cores=n.cores, beta=beta, verbose=verbose, ...)
      } else stop("Unknown method: ", method)

      sg.ids <- sample.groups %>% {split(names(.), .)}

      if (estimate.variation) {
        res$density.mad <- res %$% lapply(sg.ids, function(ids) apply(density.mat[,ids,drop=FALSE], 1, mad)) %>%
          Reduce(`+`, .)
        res$density.sd <-res %$% lapply(sg.ids, function(ids) apply(density.mat[,ids,drop=FALSE], 1, sd)) %>%
          Reduce(`+`, .)
        res$missed.sample.frac <- res$density.mat %>% {(. / rowMeans(.)) < 0.05} %>% rowMeans()
      }

      self$test.results[[name]] <- res

      return(invisible(res))
    },

    #' @description Plot cell density depending on the method that was used for estimating `cao$test.resulst[[name]]`
    #' @param add.points default is TRUE, add points to cell density figure
    #' @param contours specify cell types for contour, multiple cell types are also supported
    #' @param contour.color color for contour line
    #' @param contour.conf confidence interval of contour
    #' @param name slot in which to saved results from estimateCellDensity (default: 'cell.density')
    #' @param ... plot style parameters forwarded to \link[sccore:styleEmbeddingPlot]{sccore::styleEmbeddingPlot}.
    #' @return A ggplot2 object
    plotCellDensity = function(show.grid=TRUE, add.points=TRUE, size=0.1, show.legend=FALSE, palette=NULL,
                               point.col='#313695', contours=NULL, contour.color='black', contour.conf='10%',
                               name='cell.density', show.cell.groups=TRUE, cell.groups=self$cell.groups,
                               font.size=c(2, 4), color.range=c(0, "99%"), ...) {
      dens.res <- private$getResults(name, 'estimateCellDensity()')
      private$checkCellEmbedding()

      if (is.null(palette)) {
        palette <- c("#FFFFFF", brewerPalette("YlOrRd", rev=FALSE)(9)) %>% colorRampPalette()
      }

      ps <- names(dens.res$cond.densities) %>% sn() %>% lapply(function(l) {
        color.range %<>% parseLimitRange(unlist(dens.res$cond.densities))
        if (dens.res$method =='graph') {
          p <- self$plotEmbedding(colors=dens.res$cond.densities[[l]], size=size, title=l, show.legend=show.legend,
                                  legend.title='Density', color.range=color.range, palette=palette, ...)
        } else {# dens.res$method =='kde'
          condition.per.cell <- self$getConditionPerCell()
          p <- dens.res %$% data.frame(density.emb, z=cond.densities[[l]]) %>%
            plotDensityKde(bins=dens.res$bins, lims=color.range, title=l, show.legend=show.legend,
                           show.grid=show.grid, plot.theme=self$plot.theme, palette=palette, ...)

          if (add.points) {
            emb <- as.data.frame(self$embedding) %>% set_colnames(c('x', 'y')) %>% cbind(z=1)
            nnames <- condition.per.cell %>% {names(.)[. == l]} %>% sample(min(2000, length(.)))
            p <- p + geom_point(data=emb[nnames, ], aes(x=x, y=y), col=point.col, size=0.00001, alpha=0.2)
          }
        }

        if (show.cell.groups) {
          p %<>% transferLabelLayer(self$plotEmbedding(groups=cell.groups), font.size=font.size)
        }
        p
      })

      if (!is.null(contours)) {
        cn.geoms <- private$getDensityContours(groups=contours, conf=contour.conf, color=contour.color)
        ps %<>% lapply(`+`, cn.geoms)
      }

      ps %<>% lapply(`+`, theme(legend.background=element_blank()))
      return(ps)
    },

    plotCellDensityVariation = function(type='mad', plot.type='embedding', name='cell.density', cutoff=NULL,
                                        condition=c('both', 'ref', 'target'), ...) {
      dens.res <- private$getResults(name, 'estimateCellDensity(estimate.variation=TRUE)')
      condition <- match.arg(condition)
      if (type == 'mad') {
        name <- 'MAD'
        scores <- dens.res$density.mad
      } else if (type == 'sd') {
        name <- 'SD'
        scores <- dens.res$density.sd
      } else if (type == 'sample.frac') {
        name <- 'Missed sample frac.'
        scores <- dens.res$missed.sample.frac
      } else stop("Unknown type: ", type)

      if (is.null(scores))
        stop("To use this function, please re-run estimateCellDensity() with estimate.variation=TRUE")

      if ((dens.res$method == 'graph') && (condition != 'both')) {
        subgr <- if (condition == 'ref') self$ref.level else self$target.level
        scores %<>% .[self$getConditionPerCell()[names(.)] == subgr]
      }

      if (plot.type == 'hist') {
        gg <- ggplot(data.frame(var=scores), aes(x=scores)) + geom_histogram() +
          xlab(name) + self$plot.theme
      } else if (plot.type == 'embedding') {
        if (!is.null(cutoff)) {
          scores <- scores[abs(scores) > cutoff]
        }

        if (dens.res$method == 'graph') {
          gg <- self$plotEmbedding(colors=scores, plot.na=FALSE, legend.title=name, ...)
        } else {
          gg <- dens.res %$% data.frame(density.emb, z=scores) %>%
            plotDensityKde(bins=dens.res$bins, plot.theme=self$plot.theme, ...)
        }

      } else stop("Unknown plot.type: ", plot.type)

      return(gg)
    },

    #' @description estimate differential cell density
    #' @param type method to calculate differential cell density; permutation, t.test, wilcox or subtract (target subtract ref density);
    #' @param adjust.pvalues whether to adjust Z-scores for multiple comparison using BH method (default: FALSE for type='sutract', TRUE for everything else)
    #' @param name slot with results from estimateCellDensity. New results will be appended there. (Default: 'cell.density')
    estimateDiffCellDensity=function(type='permutation', adjust.pvalues=(type != 'subtract'), name='cell.density',
                                     n.permutations=400, smooth=TRUE, verbose=self$verbose, n.cores=self$n.cores, ...){
      dens.res <- private$getResults(name, 'estimateCellDensity')
      density.mat <- dens.res$density.mat
      if (dens.res$method == 'kde'){
        density.mat <- density.mat[dens.res$density.emb$counts > 0,]
        bins <- dens.res$bins

        sig.ids <- rownames(density.mat) %>% as.integer()
        graph <- lapply(sig.ids, function (i) c(i+1, i-1, i+bins, i-bins)) %>%
          mapIds(sig.ids) %>% igraph::graph_from_adj_list()
        igraph::V(graph)$name <- sig.ids
      } else if (adjust.pvalues && smooth) {
        graph <- extractCellGraph(self$data.object)
      }

      if (!adjust.pvalues) {
        if (type %in% c('permutation', 'permutation.mean')) {
          score <- density.mat %>%
            diffCellDensityPermutations(sample.groups=self$sample.groups, ref.level=self$ref.level,
                                        target.level=self$target.level, type=type, verbose=verbose,
                                        n.permutations=n.permutations, n.cores=n.cores) %>% .$score
        } else {
          score <- density.mat %>%
            diffCellDensity(self$sample.groups, ref.level=self$ref.level, target.level=self$target.level, type=type)
        }
        res <- list(raw=score)
      } else {
        perm.res <- density.mat %>%
          diffCellDensityPermutations(sample.groups=self$sample.groups, ref.level=self$ref.level,
                                      target.level=self$target.level, type=type, verbose=verbose,
                                      n.permutations=n.permutations, n.cores=n.cores)

        res <- list(
          raw=perm.res$score,
          adj=perm.res %$% adjustZScoresByPermutations(
            score, permut.scores, smooth=smooth, graph=graph, n.cores=n.cores, verbose=verbose, ...
          )
        )
      }

      self$test.results[[name]]$diff[[type]] <- res

      return(invisible(self$test.results[[name]]))
    },

    #' @description estimate differential cell density
    #' @param palette color palette, default is c('blue','white','red')
    #' @param type method to calculate differential cell density; t.test, wilcox or subtract (target subtract ref density);
    #' @param contours specify cell types for contour, multiple cell types are also supported (default: NULL)
    #' @param contour.color color for contour line (default: 'black')
    #' @param contour.conf confidence interval of contour (default: '10%')
    #' @param name slot with results from estimateCellDensity. New results will be appended there. (Default: 'cell.density')
    plotDiffCellDensity=function(type='permutation', name='cell.density', size=0.2, palette=NULL,
                                 adjust.pvalues=NULL, contours=NULL, contour.color='black', contour.conf='10%',
                                 plot.na=FALSE, color.range=NULL, mid.color='gray95',
                                 scale.z.palette=adjust.pvalues, min.z=qnorm(0.9), ...) {
      if (is.null(palette)) {
        if (is.null(self$sample.groups.palette)) {
          palette <- c('blue', mid.color,'red')
        } else {
          palette <- self$sample.groups.palette %>% {c(.[self$ref.level], mid.color, .[self$target.level])}
        }
        palette %<>% grDevices::colorRampPalette(space="Lab")
      }
      private$checkCellEmbedding()
      dens.res <- private$getResults(name, 'estimateCellDensity')
      scores <- dens.res$diff[[type]]
      if (is.null(adjust.pvalues)) {
        scores <- if (!is.null(scores$adj)) scores$adj else scores$raw
        adjust.pvalues <- TRUE
      } else if (adjust.pvalues) {
        if (is.null(scores$adj)) {
          warning("Adjusted scores are not estimated. Using raw scores. ",
                  "Please, run estimateCellDensity with adjust.pvalues=TRUE")
          scores <- scores$raw
        } else {
          scores <- scores$adj
        }
      } else {
        scores <- scores$raw
      }

      if (is.null(scores)) {
        warning("Can't find results for name, '", name, "' and type '", type,
                "'. Running estimateDiffCellDensity with default parameters.")
        self$estimateDiffCellDensity(type=type, name=name, adjust.pvalues=adjust.pvalues)
        dens.res <- self$test.results[[name]]
        scores <- dens.res$diff[[type]]
      }

      if (dens.res$method == 'graph'){
        density.emb <- self$embedding
        scores %<>% .[intersect(names(.), rownames(density.emb))]
      } else if (dens.res$method == 'kde') {
        density.emb <- dens.res$density.emb[,1:2]
      } else stop("Unknown method: ", dens.res$method)

      density.mat <- dens.res$density.mat[names(scores),]
      density.emb <- density.emb[names(scores),]

      if (is.null(color.range)) {
        color.range <- c(-1, 1) * max(abs(scores))
      } else {
        color.range %<>% parseLimitRange(scores)
        color.range <- max(abs(color.range)) %>% {c(-., .)}
        scores %<>% pmin(color.range[2]) %>% pmax(color.range[1])
      }

      leg.title <- if (type == 'subtract') 'Prop. change' else {if (adjust.pvalues) 'Z adj.' else 'Z-score'}
      gg <- self$plotEmbedding(density.emb, colors=scores, size=size, legend.title=leg.title, palette=palette,
                               midpoint=0, plot.na=plot.na, color.range=color.range, ...)

      if (scale.z.palette && (type != 'subtract')) {
        gg$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
        gg <- gg + getScaledZGradient(min.z=min.z, palette=palette, color.range=color.range)
      }

      if (!is.null(contours)) {
        gg <- gg + private$getDensityContours(groups=contours, conf=contour.conf, color=contour.color)
      }
      return(gg)
    },

    #' @title Plot inter-sample expression distance
    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
    #' @param name Test results to plot (default=expression.shifts)
    #' @param joint whether to show joint boxplot with the expression distance weighed by the sizes of cell types (default: TRUE), or show distances for each individual cell type
    #' @param show.significance whether to show statistical significance between sample groups. wilcox.test was used; (`*` < 0.05; `**` < 0.01; `***` < 0.001)
    #' @param ... other plot parameters, forwarded to \link{plotCountBoxplotsPerType}
    #' @return A ggplot2 object
    plotExpressionDistance = function(name='expression.shifts', joint=FALSE, palette=self$sample.groups.palette,
                                      show.significance=FALSE, ...) {
      cluster.shifts <- private$getResults(name, 'estimateExpressionShiftMagnitudes()')
      if (!joint) {
        df <- cluster.shifts %$%
          lapply(p.dist.info, subsetDistanceMatrix, sample.groups, cross.factor=FALSE, build.df=TRUE) %>%
          joinExpressionShiftDfs(sample.groups=cluster.shifts$sample.groups) %>%
          rename(group=Condition, variable=Type)
        plot.theme <- self$plot.theme
      } else {
        df <- cluster.shifts %$%
          prepareJointExpressionDistance(p.dist.info, sample.groups=sample.groups, return.dists=FALSE) %>%
          group_by(Var1, Var2, type1) %>%
          summarize(value=median(value)) %>%
          mutate(group=type1, variable="")
        plot.theme <- self$plot.theme +
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())
      }

      gg <- plotCountBoxplotsPerType(df, y.lab="expression distance", show.significance=show.significance,
                                     plot.theme=plot.theme, palette=palette, ...)

      return(gg)
    },

    getSampleDistanceMatrix=function(space=c('expression.shifts', 'coda', 'pseudo.bulk'), cell.type=NULL,
                                     dist=NULL, name=NULL, verbose=self$verbose, ...) {
      space <- match.arg(space)
      if ((space != 'pseudo.bulk') && (length(list(...)) > 0)) stop("Unexpected arguments: ", names(list(...)))
      sample.groups <- self$sample.groups
      if (space == 'expression.shifts') {
        if (is.null(name)) name <- 'expression.shifts'
        clust.info <- private$getResults(name, 'estimateExpressionShiftMagnitudes()')
        if (!is.null(cell.type)) { # use distances based on the specified cell type
          title <- cell.type
          p.dists <- clust.info$p.dist.info[[cell.type]]
          if (is.null(p.dists)) {
            warning("Distances were not estimated for cell type ", cell.type)
            return(NULL)
          }
        } else { # weighted expression distance across all cell types
          p.dists <- prepareJointExpressionDistance(clust.info$p.dist.info)
        }
        if (any(is.na(p.dists))) { # NA imputation
          p.dists %<>% ape::additive() %>% `dimnames<-`(dimnames(p.dists))
        }

        return(p.dists)
      }

      if (space=='coda') {
        n.cells.per.samp <- table(self$sample.per.cell)
        mat <- private$extractCodaData() %$% getRndBalances(d.counts) %$% prcomp(norm) %$% as.data.frame(x)
      } else { # space == 'pseudo.bulk'
        stop("Not implemented!")
      }

      dist %<>% parseDistance(ncol(mat), NULL, verbose=FALSE)
      if (dist == 'cor') {
        p.dists <- 1 - cor(t(mat))
      } else if (dist == 'l2') {
        p.dists <- dist(mat, method="euclidean") %>% as.matrix()
      } else if (dist == 'l1') {
        p.dists <- dist(mat, method="manhattan") %>% as.matrix()
      } else {
        stop("Unknown distance: ", dist)
      }

      return(p.dists)
    },

    #' @title Project samples to 2D space with MDS
    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes() or cao$estimateCellLoadings()
    #' @param space expression.shifts- expression shifts results from cao$estimateExpressionShiftMagnitudes(); CDA- cell composition shifts result from cao$estimateCellLoadings(); sudo.bulk- expression distance of sudo bulk
    #' @param cell.type If a name of a cell type is specified, the sample distances will be assessed based on this cell type alone. Otherwise (cell.type=NULL, default), sample distances will be estimated as an average distance across all cell types (weighted by the minimum number of cells of that cell type between any two samples being compared)
    #' @param dist 'cor' - correlation distance, 'l1' - manhattan distance or 'l2' - euclidean (default correlation distance)
    #' @param palette a set of colors to use for conditions (default: stored $sample.groups.palette)
    #' @param font.size font size of the sample labels. If NULL, the labels are not shown. (default: NULL)
    #' @param show.sample.size make point size proportional to the log10 of the number of cells per sample (default: FALSE)
    #' @param show.ticks show tick labels on axes (default: FALSE)
    #' @param show.labels show axis labels (default: FALSE)
    #' @param size point size. For `show.sample.size==TRUE`, it can be vector of length 2.  (default: 5)
    #' @return A ggplot2 object
    plotSampleDistances=function(space='expression.shifts', method='MDS', dist=NULL, name=NULL, cell.type=NULL,
                                 palette=NULL, show.sample.size=FALSE, sample.colors=NULL, color.title=NULL,
                                 title=NULL, n.permutations=2000, show.pvalues=FALSE, ...) {
      if (is.null(cell.type)) {
        n.cells.per.samp <- table(self$sample.per.cell)
      } else {
        if (is.null(title)) title <- cell.type
        n.cells.per.samp <- self$sample.per.cell %>% .[self$cell.groups[names(.)] == cell.type] %>% table()
      }

      p.dists <- self$getSampleDistanceMatrix(space=space, cell.type=cell.type, dist=dist, name=name)
      if (is.null(p.dists)) return(NULL)

      if (is.null(sample.colors) && is.null(palette)) palette <- self$sample.groups.palette
      gg <- plotSampleDistanceMatrix(
        p.dists=p.dists, sample.groups=self$sample.groups, n.cells.per.samp=n.cells.per.samp, method=method,
        sample.colors=sample.colors, show.sample.size=show.sample.size, palette=palette, color.title=color.title,
        title=title, plot.theme=self$plot.theme, ...
      )

      return(gg)
    },

    estimateMetadataSeparation=function(sample.meta, space='expression.shifts', dist=NULL, space.name=NULL,
                                        name='metadata.separation', n.permutations=5000, trim=0.05, show.warning=TRUE,
                                        verbose=self$verbose, n.cores=self$n.cores, adjust.pvalues=TRUE,
                                        p.adjust.method="BH", pvalue.cutoff=0.05) {
      p.dists <- self$getSampleDistanceMatrix(space=space, cell.type=NULL, dist=dist, name=space.name)
      if (is.null(p.dists)) return(NULL)

      if (is.data.frame(sample.meta)) {
        sample.meta %<>% lapply(setNames, rownames(.))
      } else if (!is.list(sample.meta)) {
        sample.meta %<>% list()
      }

      adj.mat <- p.dists %>% pmin(quantile(., 1 - trim)) %>%
        {pmax(0, . - quantile(., trim))} %>% {1 - . / max(.)} %>%
        matrix(ncol=ncol(p.dists))
      diag(adj.mat) <- 0

      pvalues <- plapply(sample.meta, function(mg) {
        mg <- mg[colnames(p.dists)]
        if (!is.numeric(mg)) {
          mg <- as.factor(mg)
          mg[is.na(mg)] <- table(mg) %>% which.max() %>% names()
          comp.op <- "!="
        } else {
          mg[is.na(mg)] <- median(mg, na.rm=TRUE)
          comp.op <- "-"
        }

        obs.var <- mg %>% outer(., ., comp.op) %>% {. * . * adj.mat} %>% sum()
        perm.vars <- sapply(1:n.permutations, function(i) {
          sample(mg) %>% outer(., ., comp.op) %>% {. * . * adj.mat} %>% sum()
        })

        (sum(perm.vars <= obs.var) + 1) / (n.permutations + 1)
      }, progress=(verbose && (length(sample.meta) > 1)), n.cores=n.cores, mc.preschedule=TRUE) %>% unlist()

      res <- list(metadata=sample.meta, pvalues=pvalues)
      if (adjust.pvalues) {
        pvalues %<>% p.adjust(method=p.adjust.method)
        res$padjust <- pvalues
      }

      if (show.warning && any(pvalues < pvalue.cutoff))
        warning("Significant separation by: ", paste(names(pvalues)[pvalues < pvalue.cutoff], collapse=', '))

      self$test.results[[name]] <- res

      # gg <- (-log10(pvalues)) %>% {tibble(Type=names(.), value=.)} %>%
      #   plotMeanMedValuesPerCellType(type="bar", yline=-log10(0.05), palette=palette, ylab="-log10(separation P-value)")

      return(invisible(res))
    },

    ### Cluster-free differential expression

    #' @description Estimate differential expression Z-scores between two conditions per individual cell
    #' @param max.z z-score value to winsorize the estimates for reducing impact of outliers. Default: 20.
    #' @param min.expr.frac minimal fraction of cell expressing a gene for estimating z-scores for it. Default: 0.001.
    #' @param min.n.samp.per.cond minimul number of samples per condition for estimating z-scores (default: 2)
    #' @param min.n.obs.per.samp minimul number of cells per samples for estimating z-scores (default: 2)
    #' @param robust whether to use median estimates instead of mean. Using median is more robust,
    #' but greatly increase the number of zeros in the data, leading to bias towards highly-express genes. (Default: FALSE)
    #' @param lfc.pseudocount pseudocount value for estimation of log2(fold-change)
    #' @return list with sparce matrices containing various DE metrics with genes as columns and cells as rows:
    #'   - `z`: DE Z-scores
    #'   - `reference` mean or median expression in reference samples
    #'   - `target` mean or median expression in target samples
    #'   - `lfc`: log2(fold-change) of expression
    #' Cells that have only one condition in their expression neighborhood have NA Z-scores for all genes.
    #' Results are also stored in the `cluster.free.de` field.
    estimateClusterFreeDE=function(n.top.genes=Inf, genes=NULL, max.z=20, min.expr.frac=0.01, min.n.samp.per.cond=2,
                                   min.n.obs.per.samp=2, robust=FALSE, norm.both=TRUE, adjust.pvalues=FALSE,
                                   smooth=TRUE, wins=0.01, n.permutations=200, lfc.pseudocount=1e-5,
                                   min.edge.weight=0.6, verbose=self$verbose, n.cores=self$n.cores,
                                   name="cluster.free.de") {
      if (is.null(genes)) {
        genes <- private$getTopGenes(n.top.genes, gene.selection="expression", min.expr.frac=min.expr.frac)
      }

      if (verbose)
        message("Estimating cluster-free Z-scores for ", length(genes), " most expressed genes")

      de.inp <- private$getClusterFreeDEInput(genes, min.edge.weight=min.edge.weight)
      mats <- de.inp %$% clusterFreeZScoreMat(
        cm, sample_per_cell=self$sample.per.cell[rownames(cm)], nn_ids=nns.per.cell, is_ref=is.ref,
        min_n_samp_per_cond=min.n.samp.per.cond, min_n_obs_per_samp=min.n.obs.per.samp, robust=robust,
        norm_both=norm.both, adjust_pvalues=adjust.pvalues, smooth=smooth, wins=wins, n_permutations=n.permutations,
        verbose=verbose, n_cores=n.cores
      )

      mats$z@x %<>% pmin(max.z) %>% pmax(-max.z)
      if (length(mats$z.adj@x) > 0) {
        mats$z.adj@x %<>% pmin(max.z) %>% pmax(-max.z)
      }

      lf.mat <- mats$reference
      lf.mat@x <- log2(mats$target@x + lfc.pseudocount) - log2(lf.mat@x + lfc.pseudocount)
      mats$lfc <- lf.mat

      self$test.results[[name]] <- mats
      return(invisible(self$test.results[[name]]))
    },

    getMostChangedGenes=function(n, method=c("z", "z.adj", "lfc"), min.z=0.5, min.lfc=1, max.score=20,
                                 cell.subset=NULL, excluded.genes=NULL, included.genes=NULL, name="cluster.free.de") {
      method <- match.arg(method)
      de.info <- private$getResults(name, "estimateClusterFreeDE")
      score.mat <- de.info[[method]]
      z.scores <- if (method == "lfc") de.info$z else de.info[[method]]

      score.mat@x[(abs(z.scores@x) < min.z) | abs(de.info$lfc@x) < min.lfc] <- 0

      if (!is.null(cell.subset)) {
        score.mat <- score.mat[cell.subset,]
      }

      if (is.null(included.genes)) {
        included.genes <- colnames(score.mat)
      }

      score.mat@x %<>% abs() %>% pmin(max.score)
      scores <- colSums(score.mat, na.rm=TRUE) %>% sort(decreasing=TRUE) %>%
        .[setdiff(names(.), excluded.genes)] %>% .[intersect(names(.), included.genes)] %>%
        .[1:min(n, length(.))]
      return(scores)
    },

    #' @description Estimate Cluster-free Expression Shift
    #' @param n.top.genes number of top genes for the distance estimation (default: 3000)
    #' @param min.n.between minimal number of pairs between condition for distance estimation (default: 2)
    #' @param min.n.within minimal number of pairs within one condition for distance estimation (default: `min.n.between`)
    #' @param min.n.obs.per.samp minimal number of cells per sample for using it in distance estimation (default: 3)
    #' @param normalize.both whether to normalize results relative to distances within both conditions (TRUE) or only to the control (FALSE)
    #' @param dist distance measure. Options: "cor" (correlation), "cosine" or "js" (Jensen–Shannon)
    #' @param log.vectors whether to use log10 on the normalized expression before estimating the distance.
    #' In most cases, must be TRUE for "cosine" and "cor" distances and always must be FALSE for "js". (default: `dist != 'js'`)
    #' @return Vector of cluster-free expression shifts per cell. Values above 1 correspond to difference between conditions.
    #' Results are also stored in the `cluster.free.expr.shifts` field.
    estimateClusterFreeExpressionShifts=function(n.top.genes=3000, gene.selection="z", name="cluster.free.expr.shifts",
                                                 min.n.between=2, min.n.within=max(min.n.between, 1),
                                                 min.expr.frac=0.0, min.n.obs.per.samp=3, normalize.both=FALSE,
                                                 dist="cor", log.vectors=(dist != "js"), wins=0.025,
                                                 n.permutations=500, verbose=self$verbose, n.cores=self$n.cores, ...) {
      cm <- self$getJointCountMatrix(raw=FALSE)
      genes <- private$getTopGenes(n.top.genes, gene.selection=gene.selection, cm.joint=cm, min.expr.frac=min.expr.frac)
      cm <- Matrix::t(cm[, genes])

      is.ref <- (self$sample.groups[levels(self$sample.per.cell)] == self$ref.level)

      nns.per.cell <- extractCellGraph(self$data.object) %>%
        igraph::as_adjacency_matrix() %>% as("dgTMatrix") %>%
        {setNames(split(.@j, .@i + 1), rownames(.))}

      shifts <- estimateClusterFreeExpressionShiftsC(
        cm, self$sample.per.cell[names(nns.per.cell)], nn_ids=nns.per.cell, is_ref=is.ref, min_n_between=min.n.between,
        min_n_within=min.n.within, min_n_obs_per_samp=min.n.obs.per.samp, norm_all=normalize.both, verbose=verbose,
        n_cores=n.cores, dist=dist, log_vecs=log.vectors, wins=wins, n_permutations=n.permutations, ...
      )
      self$test.results[[name]] <- shifts

      return(invisible(shifts))
    },

    #' @description Performs graph smoothing of the cluster-free DE Z-scores
    #' @param smoothing `beta` parameter of the \link[sccore:heatFilter]{heatFilter}. Default: 20.
    #' @param filter graph filter function. Default: \link[sccore:heatFilter]{heatFilter}.
    #' @param ... parameters forwarded to \link[sccore:smoothSignalOnGraph]{smoothSignalOnGraph}
    #' @return Sparse matrix of smoothed Z-scores. Results are also stored in the `cluster.free.de$z.smoothed` field.
    smoothClusterFreeZScores = function(n.top.genes=1000, smoothing=20, filter=NULL, z.adj=FALSE, gene.selection=ifelse(z.adj, "z.adj", "z"),
                                        excluded.genes=NULL, n.cores=self$n.cores, verbose=self$verbose, name="cluster.free.de", ...) {
      z.scores <- private$getResults(name, "estimateClusterFreeDE")
      z.scores <- if (z.adj) z.scores$z.adj else z.scores$z

      genes <- private$getTopGenes(n.top.genes, gene.selection=gene.selection,
                                   excluded.genes=excluded.genes, included.genes=colnames(z.scores))
      z.scores <- z.scores[,genes]

      if (verbose) message("Smoothing Z-scores for ", ncol(z.scores), " genes passed filtration")
      if (is.null(filter)) {
        filter <- function(...) heatFilter(..., beta=smoothing)
      }

      z.smoothed <- z.scores
      z.smoothed@x[is.na(z.smoothed@x)] <- 0
      z.smoothed %<>% smoothSignalOnGraph(filter=filter, graph=extractCellGraph(self$data.object),
                                          n.cores=n.cores, progress=verbose, ...)

      z.smoothed[is.na(z.scores)] <- NA
      if (z.adj) {
        self$test.results[[name]]$z.adj.smoothed <- z.smoothed
      } else {
        self$test.results[[name]]$z.smoothed <- z.smoothed
      }

      return(invisible(z.smoothed))
    },

    #' @description Estimate Gene Programs based on cluster-free Z-scores on a subsample of
    #' cells using \link[fabia:fabia]{fabia}. # TODO: update it
    #' @param n.programs maximal number of gene programs to find (parameter `p` for fabia). Default: 15.
    #' @param ... keyword arguments forwarded to \link{estimateGenePrograms}
    #' @return a list includes:
    #'   - `fabia`: \link[fabia:Factorization]{fabia::Factorization} object, result of the
    #'       \link[fabia:fabia]{fabia::fabia} call
    #'   - `sample.ids`: ids of the subsampled cells used for fabia estimates
    #'   - `scores.exact`: vector of fabia estimates of gene program scores per cell. Estimated only for the
    #'     subsampled cells.
    #'   - `scores.approx`: vector of approximate gene program scores, estimated for all cells in the dataset
    #'   - `loadings`: matrix with fabia gene loadings per program
    #'   - `gene.scores`: list of vectors of gene scores per program. Contains only genes, selected for
    #'     the program usin fabia biclustering.
    #'   - `bi.clusts` fabia biclustering information, result of the \link[fabia:extractBic]{fabia::extractBic} call
    estimateGenePrograms = function(method=c("pam", "leiden", "fabia"), n.top.genes=Inf, genes=NULL, n.programs=15,
                                    z.adj=FALSE, gene.selection=ifelse(z.adj, "z.adj", "z"), smooth=TRUE,
                                    abs.scores=FALSE, name="gene.programs", cell.subset=NULL, n.cores=self$n.cores,
                                    verbose=self$verbose, max.z=5, min.z=0.5, min.change.frac=0.01, de.name="cluster.free.de", ...) {
      z.scores <- private$getResults(de.name, "estimateClusterFreeDE")
      method <- match.arg(method)
      if (smooth) {
        z.scores <- if (z.adj) z.scores$z.adj.smoothed else z.scores$z.smoothed
        if (is.null(z.scores)) stop("Smoothed z.scores were not estimated. Please, run smoothClusterFreeZScores first.")
      } else {
        warning("Gene programs without smoothing often produce intractable results, especially with z.adj=FALSE\n")
        z.scores <- if (z.adj) z.scores$z.adj else z.scores$z
      }

      if (!is.null(genes)) {
        z.scores <- z.scores[,genes]
      } else if (!is.infinite(n.top.genes)) {
        genes <- private$getTopGenes(n.top.genes, gene.selection=gene.selection, included.genes=colnames(z.scores),
                                     cell.subset=cell.subset)
        z.scores <- z.scores[,genes]
      }

      if (!is.null(cell.subset)) {
        z.scores <- z.scores[cell.subset,]
      }

      z.scores@x %<>% pmax(-max.z) %<>% pmin(max.z)
      z.scores@x[abs(z.scores@x) < min.z] <- 0

      if (min.change.frac > 0) {
        z.scores.bin <- drop0(z.scores)
        z.scores.bin@x[] <- 1
        z.scores %<>% .[, colMeans(z.scores.bin, na.rm=TRUE) >= min.change.frac]
      }

      if (abs.scores) z.scores@x %<>% abs()

      z.scores %<>% as.matrix() %>% na.omit()

      if (method == "fabia") {
        warning("fabia method is deprecated")
        checkPackageInstalled('fabia', bioc=TRUE)
        res <- estimateGeneProgramsFabia(z.scores, n.programs, ...)
        fr <- res$fabia

        bi.clusts <- fabia::extractBic(fr)
        mask <- bi.clusts$bic[,"bixv"] %>% {sapply(., length) > 0}
        res$scores.exact <- t(fr@Z[mask,, drop=FALSE])
        res$scores.approx <- t(t(fr@L[,mask, drop=FALSE]) %*% (Matrix::t((z.scores[,rownames(fr@L)]) - fr@center) / fr@scaleData))
        res$loadings <- fr@L[,mask, drop=FALSE]
        res$gene.scores <- apply(bi.clusts$bic, 1, `[[`, "bixv")[mask, drop=FALSE]
        res$bi.clusts <- bi.clusts

        colnames(res$scores.exact) <- colnames(res$scores.approx) <-
          colnames(res$loadings) <- names(res$gene.scores) <- paste0("P", 1:ncol(res$scores.exact))
      } else {
        if (method == "pam") {
          clusters <- estimateGeneClustersPam(z.scores, n.programs=n.programs)
        } else { # method == "leiden"
          clusters <- estimateGeneClustersLeiden(z.scores, n.cores=n.cores, verbose=verbose, ...)
        }

        res <- geneProgramInfoByCluster(clusters, z.scores, verbose=verbose)
      }

      res$method <- method
      self$test.results[[name]] <- res

      return(invisible(self$test.results[[name]]))
    },

    plotGeneProgramScores = function(name="gene.programs", prog.ids=NULL, build.panel=TRUE, nrow=NULL, legend.title="Score", adj.list=NULL,
                                     palette=NULL, min.genes.per.prog=10, ...) {
      gene.progs <- private$getResults(name, "estimateGenePrograms")
      if (gene.progs$method == "fabia")
        stop("fabia is deprecated and plotting is not supported anymore")

      if (is.null(prog.ids)) {
        prog.ids <- (sapply(gene.progs$genes.per.clust, length) > min.genes.per.prog) %>% which()
      }

      if (all(gene.progs$program.scores >= 0)) palette <- dark.red.palette

      ggs <- lapply(prog.ids, function(i) {
        title <- paste0("Program ", i, ". ", length(gene.progs$genes.per.clust[[i]]), " genes.")
        gg <- self$plotEmbedding(colors=gene.progs$program.scores[i,], title=title, legend.title=legend.title, palette=palette, ...) +
          theme(legend.background = element_blank())
        if (!is.null(adj.list)) gg <- gg + adj.list
        gg
      })

      if (length(ggs) == 1)
        return(ggs[[1]])

      if (build.panel)
        return(cowplot::plot_grid(plotlist=ggs, nrow=nrow))

      return(ggs)
    },

    plotGeneProgramGenes = function(program.id, name="gene.programs", ordering=c("similarity", "loading"), max.genes=9, plots="z.adj", ...) {
      ordering <- match.arg(ordering)
      gene.progs <- private$getResults(name, "estimateGenePrograms")
      if (gene.progs$method == "fabia")
        stop("fabia is deprecated and plotting is not supported anymore")

      scores <- if (ordering == "similarity") gene.progs$sim.scores else gene.progs$loading.scores
      scores %<>% .[[program.id]]
      if (is.null(scores))
        stop("Can't find program", program.id)

      scores %<>% head(max.genes)
      return(self$plotGeneExpressionComparison(scores=scores, plots=plots, ...))
    },

    #' @description Plot cluster-free expression shift z-scores
    #' @param cell.groups cell type labels. Set to NULL if it shouldn't be shown
    #' @param plot.both.conditions show both case and control cells. Normally, showing control cells doesn't
    #' make sense, as control cells always have small distance from control.
    #' @param font.size size range for cell type labels
    #' @param ... parameters forwarded to \link[sccore:embeddingPlot]{embeddingPlot}
    plotClusterFreeExpressionShifts = function(cell.groups=self$cell.groups, smooth=TRUE, plot.na=FALSE,
                                               name="cluster.free.expr.shifts",
                                               color.range=c("0", "97.5%"), alpha=0.2, font.size=c(3,5), adj.list=NULL,
                                               palette=brewerPalette("YlOrRd", rev=FALSE), build.panel=TRUE, ...) {
      shifts <- private$getResults(name, "estimateClusterFreeExpressionShifts")
      private$checkCellEmbedding()
      z.scores <- shifts$z_scores
      shifts <- if (smooth) shifts$shifts_smoothed else shifts$shifts

      shifts %<>% na.omit()
      color.range %<>% parseLimitRange(shifts)
      shifts %<>% pmax(color.range[1]) %>% pmin(color.range[2])

      ggs <- mapply(function(cls, lt) {
        self$plotEmbedding(colors=cls, plot.na=plot.na, alpha=alpha, palette=palette, legend.title=lt, ...) +
          theme(legend.background = element_blank())
      }, list(shifts, z.scores), c("Distance", "Z-score"), SIMPLIFY=FALSE)

      if (!is.null(cell.groups)) {
        ggs %<>% lapply(transferLabelLayer, self$plotEmbedding(groups=cell.groups), font.size=font.size)
      }

      if (!is.null(adj.list)) ggs %<>% lapply(`+`, adj.list)
      if (build.panel) ggs %<>% cowplot::plot_grid(plotlist=., ncol=2, labels=c("Shifts", "Adj. z-scores"))

      return(ggs)
    },

    plotMostChangedGenes = function(n.top.genes, method="z", min.z=0.5, min.lfc=1, max.score=20, cell.subset=NULL, excluded.genes=NULL, ...) {
      scores <- self$getMostChangedGenes(n.top.genes, method=method, min.z=min.z, min.lfc=min.lfc, max.score=max.score,
                                         cell.subset=cell.subset, excluded.genes=excluded.genes)
      self$plotGeneExpressionComparison(scores=scores, cell.subset=cell.subset, ...)
    },

    plotGeneExpressionComparison=function(genes=NULL, scores=NULL, max.expr="97.5%", plots=c("z.adj", "z", "expression"),
                                          min.z=qnorm(0.9), max.z=4, max.z.adj=NULL, max.lfc=3, smoothed=FALSE,
                                          gene.palette=dark.red.palette, z.palette=NULL, z.adj.palette=z.palette,
                                          lfc.palette=NULL, scale.z.palette=TRUE, plot.na=-1, adj.list=NULL,
                                          build.panel=TRUE, nrow=1, cell.subset=NULL, groups=NULL, subgroups=NULL,
                                          keep.limits=NULL, name="cluster.free.de", ...) {
      unexpected.plots <- setdiff(plots, c("z.adj", "z", "lfc", "expression"))
      if (length(unexpected.plots) > 0) stop("Unexpected values in `plots`: ", unexpected.plots)

      if (is.null(genes)) {
        if (is.null(scores)) stop("Either 'genes' or 'scores' must be provided")
        genes <- names(scores)
      }

      if (length(plots) == 0) return(NULL)

      de.info <- private$getResults(name, "estimateClusterFreeDE")

      if (is.null(z.palette)) {
        if (is.null(self$sample.groups.palette)) {
          z.palette <- c('blue', 'grey90', 'red')
        } else {
          z.palette <- c(self$sample.groups.palette[self$ref.level], 'grey90',
                         self$sample.groups.palette[self$target.level])
        }
        z.palette %<>% grDevices::colorRampPalette(space="Lab")
      }

      if (is.null(z.adj.palette)) z.adj.palette <- z.palette

      plot.parts <- prepareGeneExpressionComparisonPlotInfo(
        de.info, genes=genes, plots=plots, smoothed=smoothed, max.z=max.z, max.z.adj=max.z.adj, max.lfc=max.lfc,
        z.palette=z.palette, z.adj.palette=z.adj.palette, lfc.palette=lfc.palette
      )

      condition.per.cell <- self$getConditionPerCell()

      if (is.null(groups) && !is.null(cell.subset)) {
        groups <- rownames(self$embedding) %>% {setNames(. %in% cell.subset, .)}
        subgroups <- TRUE
      }

      if (is.null(keep.limits)) {keep.limits <- FALSE}

      ggs <- lapply(genes, function(g) {
        lst <- list()

        title <- if (is.null(scores)) paste0(g, ". ") else paste0(g, ": ", signif(scores[g], 3), ". ")
        for (n in names(plot.parts)) {
          pp <- plot.parts[[n]]
          gg <- self$plotEmbedding(colors=pp$scores[,g], title=paste0(title, pp$title),
                                   color.range=c(-pp$max, pp$max), plot.na=plot.na, legend.title=pp$leg.title,
                                   palette=pp$palette, groups=groups, subgroups=subgroups, keep.limits=keep.limits, ...)

          if ((n == "z.adj") && scale.z.palette) {
            gg$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
            gg <- gg + getScaledZGradient(min.z=min.z, palette=pp$palette, color.range=c(-pp$max, pp$max))
          }

          lst <- c(lst, list(gg))
          title <- ""
        }

        if ("expression" %in% plots) {
          expr <- extractGeneExpression(self$data.object, g)
          m.expr <- parseLimitRange(c(0, max.expr), expr)[2]
          lst <- lapply(unique(condition.per.cell), function(sg) {
            self$plotEmbedding(colors=expr[condition.per.cell[names(expr)] == sg], title=paste(title, sg),
                               color.range=c(0, m.expr), legend.title="Expression", plot.na=FALSE, palette=gene.palette,
                               groups=groups, subgroups=subgroups, keep.limits=keep.limits, ...)
          }) %>% {c(lst, .)}
          title <- ""
        }

        lst <- lapply(lst, function(x) x + theme(legend.background = element_blank()) + adj.list)
        if (build.panel) {
          lst <- if ((length(lst) > 1)) cowplot::plot_grid(plotlist=lst, nrow=nrow) else lst[[1]]
        }

        lst
      })

      if (length(genes) == 1) return(ggs[[1]])
      return(ggs)
    },

    getConditionPerCell=function() {
      self$sample.per.cell %>%
        {setNames(as.character(self$sample.groups[as.character(.)]), names(.))} %>%
        as.factor()
    },

    getJointCountMatrix=function(force=FALSE, raw=TRUE) {
      cache.name <- if (raw) "joint.count.matrix" else "joint.count.matrix.norm"
      if (force || is.null(self$cache[[cache.name]])) {
        self$cache[[cache.name]] <- extractJointCountMatrix(self$data.object, raw=raw)
      }

      return(self$cache[[cache.name]])
    },

    getGOEnvironment=function(org.db, verbose=FALSE, ignore.cache=NULL) {
      checkPackageInstalled("clusterProfiler", bioc=TRUE)
      if (!is.null(self$cache$go.environment)) {
        if (is.null(ignore.cache)) {
          message("Using stored GO environment. Use `ignore.cache=TRUE` if you want to re-estimate it. ",
                  "Set `ignore.cache=FALSE` to suppress this message.")
          return(self$cache$go.environment)
        }

        if (!ignore.cache) return(self$cache$go.environment)
      }

      if(class(org.db) != "OrgDb")
        stop("'org.db' must be of class 'OrgDb'. Please input an organism database.")

      self$cache$go.environment <- c("BP", "CC", "MF") %>% sn() %>%
        plapply(function(n) clusterProfiler:::get_GO_data(org.db, n, "ENTREZID") %>%
                  as.list() %>% as.environment(), n.cores=1, progress=verbose)
      return(self$cache$go.environment)
    }
  ),

  private = list(
    getResults=function(name, suggested.function=NULL) {
      if (!is.null(self$test.results[[name]]))
        return(self$test.results[[name]])

      msg <- paste0("A result named \"", name, "\" cannot be found.");
      if(!is.null(suggested.function)) {
        msg <- paste(msg, "Please first run", suggested.function)
      }
      stop(msg)
    },

    getTopGenes = function(n, gene.selection=c("z", "z.adj", "lfc", "expression", "od"), cm.joint=NULL,
                           min.expr.frac=0.0, excluded.genes=NULL, included.genes=NULL, name="cluster.free.de", ...) {
      gene.selection <- match.arg(gene.selection)
      if ((gene.selection %in% c("z", "lfc")) && is.null(self$test.results[[name]])) {
        warning("Please run estimateClusterFreeDE() first to use gene.selection='z' or 'lfc'. Fall back to gene.selection='expression'.")
        gene.selection <- "expression"
      }

      if (min.expr.frac > 0) {
        if (is.null(cm.joint)) {
          cm.joint <- self$getJointCountMatrix()
        }

        cm.joint@x <- 1 * (cm.joint@x > 0)
        excluded.genes %<>% union(colnames(cm.joint)[colMeans(cm.joint, na.rm=TRUE) < min.expr.frac])
      }

      if (gene.selection %in% c("z", "z.adj", "lfc")) {
        genes <- names(self$getMostChangedGenes(Inf, method=gene.selection, ...))
      } else if (gene.selection == "od") {
        genes <- extractOdGenes(self$data.object)
      } else { # expression
        if (is.null(cm.joint)) {
          cm.joint <- self$getJointCountMatrix()
        }

        cm.joint@x <- 1 * (cm.joint@x > 0)
        genes <- colMeans(cm.joint, na.rm=TRUE) %>% sort(decreasing=TRUE) %>% names()
      }

      if (is.null(included.genes)) {
        included.genes <- genes
      }
      genes %<>% setdiff(excluded.genes) %>% intersect(included.genes) %>% .[1:min(length(.), n)]
      if (is.finite(n) && (length(genes) < n))
        warning("Not enough genes pass the criteria. Returning ", length(genes), " genes.")

      return(genes)
    },

    getOntologyPvalueResults=function(genes, type, p.adj=0.05, min.genes=1,
                                      subtype=NULL, cell.subgroups=NULL, name = NULL) {
      if (!type %in% c("GO", "DO", "GSEA"))
        stop("'type' must be 'GO', 'DO', or 'GSEA'.")
      if (is.null(name)) {
        name <- type
      }

      if (!is.null(subtype) && !all(subtype %in% c("BP", "CC", "MF")))
        stop("'subtype' must be 'BP', 'CC', or 'MF'.")

      if ((length(genes) != 1) || (!genes %in% c("down","up","all")))
        stop("'genes' must be 'down', 'up', or 'all'.")

      ont.res <- self$test.results[[name]]$res
      if (is.null(ont.res)) stop(paste0("No results found for '", type, "'. Please run 'estimateOntology' first."))

      ont.res %<>% prepareOntologyPlotData(type, p.adj, min.genes)

      # Extract genes and subgroups
      if (type == "GSEA") {
        ont.res %<>% addGseaGroup() %>% rename(geneID=core_enrichment)
        if (genes == "up") {
          ont.res %<>% filter(enrichmentScore > 0)
        } else if (genes == "down") {
          ont.res %<>% filter(enrichmentScore < 0)
        }
      } else {
        ont.res %<>% .[[genes]]
      }

      if (is.null(ont.res))
        stop(paste0("No results found for ", genes, " genes for ", type, ".")) # TODO: Include GSEA

      if (!is.null(cell.subgroups)) {
        if (!any(cell.subgroups %in% unique(ont.res$Group)))
          stop("None of 'cell.subgroups' was found in results.")

        ont.res %<>% filter(Group %in% cell.subgroups)
      }

      if (!is.null(subtype)) {
        ont.res %<>% filter(Type %in% subtype)
      }

      return(ont.res)
    },

    getOntologyHeatmapInfo=function(genes="up", type="GO", subtype="BP", min.genes=1, p.adj=0.05,
                                    selection=c("all", "common", "unique"), clusters=TRUE, cluster.name=NULL,
                                    cell.subgroups=NULL, only.family.children=FALSE, description.regex=NULL,
                                    clust.naming="medoid", return.descriptions=FALSE) {
      # Checks
      selection <- match.arg(selection)
      if (!is.null(cell.subgroups) && (length(cell.subgroups) == 1))
        stop("'cell.subgroups' must contain at least two groups. Please use plotOntology instead.")

      if (only.family.children) {
        fams <- self$test.results[[type]]$families
        if (is.null(fams))
          stop("No ontology family results found, please run 'estimateOntologyFamilies' first",
               " or set only.family.children=FALSE")
      }

      # Extract results
      desc.per.clust <- NULL
      if (!clusters) {
        ont.sum <- private$getOntologyPvalueResults(
          genes=genes, type=type, subtype=subtype, cell.subgroups=cell.subgroups, p.adj=p.adj, min.genes=min.genes
        )
        group.field <- "Description"
      } else {
        name <- if (is.null(cluster.name)) getOntClustField(subtype, genes) else cluster.name
        if (is.null(self$test.results[[type]][[name]])) {
          if (!is.null(cluster.name))
            stop("Can't find the results for ", cluster.name) # stop if user specified a wrong cluster.name

          warning("Can't find the results for ", name, ". Running estimateOntologyClusters()...\n")
          self$estimateOntologyClusters(type=type, subtype=subtype, genes=genes, name=name, p.adj=p.adj,
                                        min.genes=min.genes, clust.naming=clust.naming)
        }

        ont.sum <- self$test.results[[type]][[name]]$df

        if (!is.null(cell.subgroups)) ont.sum %<>% filter(Group %in% cell.subgroups)
        group.field <- "ClusterName"
        desc.per.clust <- ont.sum %$% split(Description, ClusterName) %>% lapply(unique)
      }

      if (only.family.children) {
        ont.sum %<>% getOntologyFamilyChildren(fams=fams, subtype=subtype, genes=genes)
      }

      if (!is.null(description.regex)) ont.sum %<>% .[grep(description.regex, .$Description),]

      ont.sum <- -groupOntologiesByCluster(ont.sum, field=group.field)

      if (nrow(ont.sum) == 0) {
        warning("No ontologies pass the filtration for type=", type, ", subtype=", subtype, " and genes=", genes)
        return(NULL)
      }

      if (selection=="unique") {
        ont.sum %<>% .[rowSums(abs(.) > 0) == 1,,drop=FALSE]
      } else if(selection=="common") {
        ont.sum %<>% .[rowSums(abs(.) > 0) > 1,,drop=FALSE]
      }
      if (nrow(ont.sum) == 0) {
        warning("Nothing to plot. Try another selection.")
        return(NULL)
      }

      ont.sum %<>% .[, colSums(abs(.)) > 0, drop=FALSE]
      if (return.descriptions)
        return(list(ont.sum=ont.sum, desc.per.clust=desc.per.clust))

      return(ont.sum)
    },

    extractCodaData = function(ret.groups=TRUE, cell.groups=self$cell.groups, cells.to.remove=NULL, cells.to.remain=NULL, samples.to.remove=NULL) {
      d.counts <- cell.groups %>% data.frame(anno=., group=self$sample.per.cell[names(.)]) %>%
        table() %>% rbind() %>% t()

      if (!is.null(cells.to.remove)) d.counts %<>% .[,!(colnames(.) %in% cells.to.remove)]
      if (!is.null(cells.to.remain)) d.counts %<>% .[,colnames(.) %in% cells.to.remain]
      if (!is.null(samples.to.remove)) d.counts %<>% .[!(rownames(.) %in% samples.to.remove),]

      if (!ret.groups)
        return(d.counts)

      d.groups <- (self$sample.groups[rownames(d.counts)] == self$target.level) %>%
        setNames(rownames(d.counts))

      return(list(d.counts = d.counts,
                  d.groups = d.groups))
    },

    #' @description Extract contours from embedding
    #' @param groups specify cell groups for contour, multiple cell groups are also supported
    #' @param conf confidence interval of contour
    getDensityContours = function(groups, color='black', linetype=2, conf="10%", n.cores=1, verbose=FALSE, ...) {
      cnl <- sn(groups) %>% plapply(function(g) {
        getDensityContour(
          self$embedding, cell.groups=self$cell.groups, linetype=linetype, group=g, conf=conf, color=color, ...
        )
      }, n.cores=n.cores, progress=verbose, mc.preschedule=TRUE) %>% do.call(c, .)
      return(cnl)
    },

    checkCellEmbedding = function(embedding=self$embedding) {
      if(is.null(embedding) || ncol(embedding) != 2)
        stop("self$embedding must contain 2D cell embedding")

      if(is.null(rownames(embedding)))
        stop("self$embedding must have rownames, equal to cell ids")
    },

    getClusterFreeDEInput = function(genes, min.edge.weight=0.0) {
      cm <- self$getJointCountMatrix(raw=FALSE)
      is.ref <- (self$sample.groups[levels(self$sample.per.cell)] == self$ref.level)

      adj.mat <- extractCellGraph(self$data.object) %>% igraph::as_adj()
      diag(adj.mat) <- 1
      cell.names <- intersect(rownames(cm), rownames(adj.mat))

      adj.mat %<>% .[cell.names, cell.names, drop=FALSE] %>% as("dgTMatrix")

      if (min.edge.weight > 1e-10) {
        samp.per.cell <- self$sample.per.cell[cell.names]
        adj.mat@x[(samp.per.cell[adj.mat@i + 1] != samp.per.cell[adj.mat@j + 1]) & (adj.mat@x < min.edge.weight)] <- 0.0
        adj.mat %<>% drop0() %>% as("dgTMatrix")
      }

      nns.per.cell <- split(adj.mat@j, adj.mat@i) %>% setNames(cell.names)
      cm %<>% .[cell.names, genes, drop=FALSE]

      return(list(cm=cm, adj.mat=adj.mat, is.ref=is.ref, nns.per.cell=nns.per.cell))
    }
  )
)
