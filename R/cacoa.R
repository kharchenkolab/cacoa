#' @import methods dplyr R6
#' @importFrom grDevices colorRampPalette contourLines
#' @importFrom graphics legend
#' @importFrom stats as.dendrogram as.dist as.formula as.hclust cmdscale coef
#' @importFrom stats cor cutree dist hclust lm loadings median model.matrix na.omit p.adjust pf prcomp
#' @importFrom stats qnorm quantile relevel reorder rmultinom runif sd setNames symnum var wilcox.test
#' @importFrom utils stack
NULL


## for magrittr and dplyr functions in the package
if (getRversion() >= "2.15.1"){
  utils::globalVariables(c("%||%", ".", "CellFrac", "Cluster", "Condition", "Description", "ENTREZID",
  "Freq", "G1", "G2", "Gene", "Group", "ID", "LI", "N", "NCells",
  "S1", "S2", "SYMBOL", "Sample", "Sign", "Stability", "Type", "UI",
  "Var1", "Z", "child_distance", "child_go_id", "cmp", "color", "con.names", "condition", "cont",
  "distance", "estimate", "eval.points", "fill", "from", "geneID", "group", "ind",
  "log2FoldChange", "med", "membership", "n.cells", "node", "padj", "parent_distance",
  "pvalue", "qvalue", "se", "size", "to", "value", "values", "variable", "x", "xend", "y", "yend", "z"))
}


#' @title Cacoa R6 class
#'
#' @description The class encompasses etc etc
#' @param sample.groups a two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param cell.groups vector Indicates cell groups with cell IDs as names (default: stored vector)
#' @param n.cores numeric Number of cores for parallelization
#' @param verbose boolean Whether to show progress
#' @param ref.level reference sample group level (default=self$ref.level)
#' @param name string Field name where the test results are stored
#' @param n.top.genes numeric Number of top genes for estimation
#' @param p.adj numeric Cut-off for adjusted p-values (default=0.05)
#' @param p.adjust.method character string Method for calculating adjusted p-values (default="BH")
#' @param gene.selection character string Method to select top genes, "z" selects genes by cluster-free Z-score change, "lfc" uses log2(fold-change) instead,
#' "expression" picks the most expressed genes and "od" picks overdispersed genes.  Default: "z".
#' @param excluded.genes List of genes to exclude during estimation. For example, a list of mitochondrial genes.
#' @param sample.subset subset data for analysis only to the given samples
#' @export Cacoa
Cacoa <- R6::R6Class("Cacoa", lock_objects=FALSE,
  public = list(
    #' @field n.cores Number of cores (default=1)
    n.cores = 1,

    #' @field verbose boolean Whether to provide verbose output with diagnostic messages (default=FALSE)
    verbose = FALSE,

    #' @field test.results list Results of the estimations, ready to use (default=list())
    test.results = list(),

    #' @field cache list Intermediate results of the estimations, which can be used during some other computations (default=list())
    cache = list(),

    #' @field data.object list The main object storing data (Conos or Seurat) (default=list())
    data.object = list(),

    #' @field sample.groups 2-factor vector with annotation of groups/condition per sample (default=NULL)
    sample.groups = NULL,

    #' @field cell.groups Named factor with cell names with cluster per cell (default=NULL)
    cell.groups = NULL,

    #' @field embedding 2D embedding to visualize the cells in (default=NULL)
    embedding = NULL,

    #' @field sample.per.cell Named factor with cell names (default=NULL)
    sample.per.cell = NULL,

    #' @field ref.level Reference level for sample.group vector (default=NULL)
    ref.level = NULL,

    #' @field target.level Target/disease level for sample.group vector
    target.level = NULL,

    #' @field sample.groups.palette Color palette for the sample.groups (default=NULL)
    sample.groups.palette = NULL,

    #' @field cell.groups.palette Color palette for the cell.groups (default=NULL)
    cell.groups.palette = NULL,

    #' @field plot.theme ggplot2 theme for all plots (default=NULL)
    plot.theme = NULL,

    #' @field plot.params list with parameters, forwarded to all `plotEmbedding` calls (default=NULL)
    plot.params = NULL,


    #' @description Initialize Cacoa class
    #'
    #' @param data.object Object used to initialize the Cacoa class. Either a raw or normalized count matrix, Conos object, or Seurat object.
    #' @param sample.groups a two-level factor on the sample names describing the conditions being compared (default: extracted from `data.object`)
    #' @param cell.groups vector Indicates cell groups with cell names (default: extracted from `data.object`)
    #' @param sample.per.cell vector Sample name per cell (default: extracted from `data.object`)
    #' @param ref.level reference sample group level
    #' @param target.level target sample group level
    #' @param sample.groups.palette Color palette for the sample.groups (default=NULL)
    #' @param cell.groups.palette Color palette for the cell.groups (default=NULL)
    #' @param embedding embedding 2D embedding to visualize the cells in (default: extracted from `data.object`)
    #' @param graph.name graph name for Seurat object, ignored otherwise (default=NULL)
    #' @param assay.name assay name for Seurat object, ignored otherwise (default="RNA")
    #' @param n.cores Number of cores for parallelization (default=1)
    #' @param verbose boolean Whether to provide verbose output with diagnostic messages (default=TRUE)
    #' @param plot.theme ggplot2 plot theme (default=ggplot2::theme_bw())
    #' @param plot.params list with parameters, replacing defaults from \link[sccore:embeddingPlot]{embeddingPlot} (default=NULL)
    #' @return a new 'Cacoa' object
    #' @examples 
    #' # Is it highly recommended that sample.groups and cell.groups are assigned in the initialization call. 
    #' # Here, "con" is a Conos object.
    #' \dontrun{
    #' sample.groups <- c("control","control","disease","disease")
    #' names(sample.groups) <- names(con$samples)
    #' }
    #' # cell.groups should be a named factor where names are cell names corresponding to cell names in the data object. 
    #' # For Conos objects, they should overlap with rownames(con$embedding)
    #' \dontrun{
    #' cell.groups <- my.named.annotation.factor
    #' cao <- Cacoa$new(data.object = con, sample.groups = sample.groups, cell.groups =  ref.level = "control", target.level = "disease")
    #' }
    initialize=function(data.object, sample.groups=NULL, cell.groups=NULL, sample.per.cell=NULL,
                        ref.level=NULL, target.level=NULL, sample.groups.palette=NULL,
                        cell.groups.palette=NULL, embedding=NULL,
                        graph.name=NULL, assay.name="RNA", n.cores=1, verbose=TRUE,
                        plot.theme=ggplot2::theme_bw(), plot.params=NULL) {

      if ('Cacoa' %in% class(data.object)) { # copy constructor
        for (n in ls(data.object)) {
          if (!is.function(get(n, data.object))) assign(n, get(n, data.object), self)
        }

        return(NULL)
      }

      if (!is.null(sample.groups)) {
        if (length(unique(sample.groups)) != 2) {
          stop("sample.groups must have exactly two levels")
        }
        if (is.null(ref.level) && !is.null(target.level)){
          ref.level <- setdiff(unique(sample.groups), target.level)[1]
        }
        if (!is.null(ref.level) && is.null(target.level)){
          target.level <- setdiff(unique(sample.groups), ref.level)[1]
        }
        if (length(setdiff(sample.groups, c(ref.level, target.level)) > 0)){
          stop("sample.groups must contain only ref.level '", ref.level, "' and target.level '", target.level, "'")
        }
      }

      if (is.null(ref.level) || is.null(target.level)){
        stop("Both ref.level and target.level must be provided")
      }

      self$n.cores <- n.cores
      self$verbose <- verbose
      self$ref.level <- ref.level
      self$target.level <- target.level

      if ("Seurat" %in% class(data.object)) {
        if (is.null(sample.groups) || is.null(sample.per.cell)){
          stop("Both sample.groups and sample.per.cell must be provided for Seurat objects")
        }
        data.object$sample.per.cell <- sample.per.cell
        if (is.null(graph.name)){
          warning("No graph.name provided. The algorithm will use the first available graph.")
        }
        data.object@misc$graph.name <- graph.name
        data.object@misc$assay.name <- assay.name
      } else if (('Conos' %in% class(data.object))) {
        if (!is.null(graph.name)) {
          warning("graph.name is not supported for Conos objects")
        }
      } else {
        warning("Many function may be not supported for an object of class ", class(data.object));
        if (is.null(sample.groups) || is.null(sample.per.cell) || is.null(cell.groups)){
          stop("All sample.groups, sample.per.cell and cell.groups must be provided")
        }
      }

      if (any(c("dgCMatrix", "dgTMatrix", "dgEMatrix", "matrix") %in% class(data.object))) {
        data.object %<>% as("CsparseMatrix") %>% Matrix::t()
        if (max(abs(round(data.object@x) - data.object@x)) < 1e-10) {
          message("Interpreting data.object as a raw count matrix")
          attr(data.object, "raw") <- TRUE
        } else {
          message("Interpreting data.object as a normalized count matrix")
          attr(data.object, "raw") <- FALSE
        }

        if (length(setdiff(names(sample.per.cell), rownames(data.object))) > 0){
          stop("All cells in the count matrix columns must be present in sample.per.cell")
        }
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

      if (is.null(embedding)) {
        try({ # In case extractEmbedding doesn't exist for this type of objects
          embedding <- extractEmbedding(data.object)
        }, silent=TRUE)
      }

      self$plot.theme <- plot.theme
      self$embedding <- embedding
      self$plot.params <- plot.params
    },


    #' @description Calculate expression shift magnitudes of different clusters between conditions
    #'
    #' @param top.n.genes character vector Vector of top genes to show (default=NULL)
    #' @param dist.type type of expression distance: 'shift' (linear shift) 'var' (variance change) or 'total' (both) (default="shift")
    #' @param sample.per.cell Sample per cell (default=self$sample.per.cell)
    #' @param min.cells.per.sample numeric minimum cells per sample (default=10)
    #' @param ref.level character Reference level, e.g. "control" (default=self$ref.level)
    #' @param sample.groups named vector indicating sample groups with sample IDs as names (default: stored sample.groups)
    #' @param n.cores integer Number of cores for parallelization (default: stored integer)
    #' @param name character Test name (default="expression.shifts")
    #' @param dist distance metric: 'cor' - correlation distance, 'l1' - manhattan distance or 'l2' - euclidean (default=NULL, depends on dimensionality)
    #' @param min.samp.per.type numeric minimal number of samples per cell type for it to be included (default=2)
    #' @param min.gene.frac numeric minimal number of cells per cell type expressing a gene for it to be included (default=0.01)
    #' @param n.permutations numeric number of permutations for estimating normalization coefficient (default=1000)
    #' @param genes character if provided, the expression distance is estimated only based on these genes (default=NULL)
    #' @param n.pcs numeric Number of principal components for estimating expression distance (default=NULL, no PCA)
    #' @param ... extra parameters to \link{estimateExpressionChange}
    #' @return List including:
    #'   `dist.df`: a table with cluster distances (normalized if within.group.normalization=TRUE), cell type and the number of cells # TODO: update
    #'   `p.dist.info`: list of distance matrices per cell type
    #'   `sample.groups`: filtered sample groups
    #'   `cell.groups`: filtered cell groups
    #' @examples
    #' \dontrun{
    #' cao$estimateExpressionShiftMagnitudes()
    #' }
    estimateExpressionShiftMagnitudes=function(cell.groups=self$cell.groups,
      sample.per.cell=self$sample.per.cell, dist=NULL, dist.type="shift",
      min.cells.per.sample=10, min.samp.per.type=2, min.gene.frac=0.01,
      ref.level=self$ref.level, sample.groups=self$sample.groups,
      verbose=self$verbose, n.cores=self$n.cores, name="expression.shifts",
      n.permutations=1000, genes=NULL, n.pcs=NULL, top.n.genes=NULL, ...) {

      count.matrices <- extractRawCountMatrices(self$data.object, transposed=TRUE)

      if (verbose) message("Filtering data... ")
      shift.inp <- filterExpressionDistanceInput(
        count.matrices, cell.groups=cell.groups,
        sample.per.cell=self$sample.per.cell, sample.groups=self$sample.groups,
        min.cells.per.sample=min.cells.per.sample, min.samp.per.type=min.samp.per.type,
        min.gene.frac=min.gene.frac, genes=genes, verbose=verbose
      )
      if (verbose) message("done!\n")

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
        estimateExpressionChange(
          cm.per.type, sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=self$sample.per.cell,
          dist=dist, dist.type=dist.type, verbose=verbose, ref.level=ref.level,
          n.permutations=n.permutations, top.n.genes=top.n.genes, n.pcs=n.pcs, n.cores=n.cores, ...
        )

      return(invisible(self$test.results[[name]]))
    },

    #' @description Plot results from cao$estimateExpressionShiftMagnitudes() (shift.type="normal") or
    #'   cao$estimateCommonExpressionShiftMagnitudes() (shift.type="common")
    #'
    #' @param name character Results slot name (default="expression.shifts")
    #' @param type character type of a plot "bar" or "box" (default="bar")
    #' @param notch boolean Whether to show notches in the boxplot version (default=TRUE)
    #' @param show.jitter boolean Whether to show individual data points (default=FALSE)
    #' @param jitter.alpha numeric Transparency value for the data points (default=0.05)
    #' @param show.pvalues character string Which p-values to plot. Accepted values are "none", "raw", or "adjusted". (default=c("adjusted", "raw", "none"))
    #' @param ylab character string Label of the y-axis (default="normalized expression distance")
    #' @param ... additional arguments
    #' @return A ggplot2 object
    #' @examples
    #' \dontrun{
    #' cao$estimateExpressionShiftMagnitudes()
    #' cao$plotExpressionShiftMagnitudes()
    #' }
    plotExpressionShiftMagnitudes=function(name="expression.shifts", type='box', notch=TRUE, show.jitter=TRUE,
                                           jitter.alpha=0.05, show.pvalues=c("adjusted", "raw", "none"),
                                           ylab='normalized expression distance', ...) {
      show.pvalues <- match.arg(show.pvalues)

      res <- private$getResults(name, "estimateExpressionShiftMagnitudes()")
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

      plotMeanMedValuesPerCellType(df, pvalues=pvalues, show.jitter=show.jitter,jitter.alpha=jitter.alpha, notch=notch, type=type,
        palette=self$cell.groups.palette, ylab=ylab, plot.theme=self$plot.theme, yline=0.0, ...
      )
    },

    #' @description Alias for estimateDEPerCellType
    #' @param ... parameters fed to estimateDEPerCellType
    #' @return A list of DE genes
    #' @examples 
    #' \dontrun{
    #' cao$estimatePerCellTypeDE() # Deprecated
    #' }
    estimatePerCellTypeDE=function(...) {
      .Deprecated("cao$estimateDEPerCellType")
      return(self$estimateDEPerCellType(...))
    },

    #' @description Estimate differential gene expression per cell type between conditions
    #' @param cell.groups factor specifying cell types (default=self$cell.groups)
    #' @param sample.groups 2-factor vector with annotation of groups/condition per sample (default=self$sample.groups)
    #' @param ref.level character Reference level in 'sample.groups', e.g., ctrl, healthy (default=self$ref.level)
    #' @param target.level character Target level in 'sample.groups', e.g., case, diseased (default=self$target.level)
    #' @param name character string Slot in which to save the results (default='de')
    #' @param test character string Which DESeq2 test to use. The available options are "LRT", "Wald". (default="DESeq2.Wald")
    #' @param resampling.method character which resampling method should be used "loo" for leave-one-out or "bootstrap", (default=NULL, i.e. no resampling)
    #' @param n.resamplings numeric Number of resamplings to perform (default=30)
    #' @param seed.resampling numeric Seed to use for resamplings, input to set.seed() (default=239)
    #' @param min.cell.frac numeric Minimum fraction of cells to use to perform DE (default=0.05)
    #' @param covariates (default=NULL)
    #' @param common.genes boolean Whether to investigate common genes across cell groups (default=FALSE)
    #' @param cooks.cutoff boolean cooksCutoff for DESeq2 (default=FALSE)
    #' @param independent.filtering boolean independentFiltering parameter for DESeq2 (default=FALSE)
    #' @param min.cell.count numeric minimum number of cells that need to be present in a given cell type in a given sample in order to be taken into account (default=10)
    #' @param n.cells.subsample integer Number of cells to subsample (default=NULL)
    #' @param fix.n.samples Samples to be provided if resampling.method='fix.samples'.
    #' @param ... additional parameters
    #' @return A list of DE genes
    #' @examples
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' }
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

      if (!(tolower(test) %in% tolower(possible.tests))) {
        stop('Test ', test, ' is not supported. Available tests: ', paste(possible.tests, collapse=', '))
      }

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
          if (is.null(fix.n.samples)) {
            stop("fix.n.samples must be provided for resampling.method='fix.samples'")
          }
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
      n.cores.outer <- min(length(s.groups.new), n.cores)
      n.cores.inner <- max(n.cores %/% length(s.groups.new), 1)
      verbose.inner <- (verbose & (n.cores.outer == 1))

      de.res <- names(s.groups.new) %>% sn() %>% plapply(function(resampling.name) {
        estimateDEPerCellTypeInner(
          raw.mats=raw.mats, cell.groups=cell.groups, s.groups=s.groups.new[[resampling.name]],
          ref.level=ref.level, target.level=target.level, common.genes=common.genes,
          cooks.cutoff=cooks.cutoff, min.cell.count=min.cell.count, max.cell.count=max.cell.count,
          independent.filtering=independent.filtering, test=test, meta.info=covariates, gene.filter=gene.filter,
          fix.n.samples=(if (resampling.name == 'initial') NULL else fix.samples),
          n.cores=n.cores.inner, verbose=verbose.inner, return.matrix=(resampling.name == 'initial'), ...
        )
      }, n.cores=n.cores.outer, progress=(!verbose.inner & verbose), mc.preschedule=TRUE, mc.allow.recursive=TRUE)



      # if resampling: calculate median and variance on ranks after resampling
      de.res <- if (length(de.res) > 1) summarizeDEResamplingResults(de.res) else de.res[[1]]
      de.res %<>% appendStatisticsToDE(expr.fracs)
      self$test.results[[name]] <- de.res

      # TODO: add overall p-adjustment

      return(invisible(self$test.results[[name]]))
    },

    #' @description Estimate DE stability per cell type
    #' @param de.name character string DE results slot name (default='de')
    #' @param name character string Name for storing results (default='de.jaccards')
    #' @param top.n.genes numeric Number of top DE genes to return (default=NULL)
    #' @param p.val.cutoff numeric The p-value cutoff to apply for returned DE values (default=NULL)
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateDEStabilityPerCellType()
    #' }
    estimateDEStabilityPerCellType=function(de.name='de', name='de.jaccards', top.n.genes=NULL, p.val.cutoff=NULL) {
      if(!is.null(p.val.cutoff) & !is.null(top.n.genes)){
        stop('Only one threshold (top.n.genes or p.val.cutoff) should be provided')
      }
      if(is.null(p.val.cutoff) & is.null(top.n.genes)){
        stop('At least one threshold (top.n.genes or p.val.cutoff) should be provided')
      }

      de.res <- private$getResults(de.name, 'estimateDEPerCellType()')

      if(length(de.res) < 3){
        stop('Resampling was not performed')
      }

      jaccards <- estimateStabilityPerCellType(de.res, top.n.genes, p.val.cutoff)
      self$test.results[[name]] <- jaccards
    },

    #' @description Estimate DE stability per gene
    #' @param de.name character string DE results slot name (default='de')
    #' @param top.n.genes numeric Number of top DE genes to return (default=500)
    #' @param p.adj.cutoff numeric The adjusted p-value cutoff to apply for returned DE values (default=NULL)
    #' @param visualize boolean Whether to visualize results (default=FALSE)
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateDEStabilityPerGene()
    #' }
    estimateDEStabilityPerGene=function(de.name="de", top.n.genes=500, p.adj.cutoff=NULL, visualize=FALSE) {
      de.res <- self$test.results[[de.name]]
      for(cell.type in names(de.res$initial)) {
        if(!is.null(p.adj)) {
          top.n.genes <- sum(de.res$initial[[cell.type]]$res$padj <= p.adj)
        }
        genes.tmp <- rownames(de.res$initial[[cell.type]]$res)
        tmp <- rep(0, length(genes.tmp))
        n.tmp <- 0
        for(resampling.name in setdiff(names(de.res), 'initial')){
          if(is.null(de.res[[resampling.name]][[cell.type]])) next
          tmp <- tmp + 1*(rank(de.res[[resampling.name]][[cell.type]][genes.tmp,'padj']) < top.n.genes)
          n.tmp <- n.tmp + 1
        }
        de.res$initial[[cell.type]]$res$Stability <- tmp / n.tmp
      }
      self$test.results[[de.name]] <- de.res
      if(visualize){
        stability.all <- c()
        for(cell.type in names(de.res$initial)) {
          tmp <- de.res$initial[[cell.type]]$res
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
                legend.key.width = unit(0.2, 'cm')) + ylim(0, 1)
        return(p)
      }
    },

    #' @description Plot DE stability per cell type
    #' @param name character string DE stability results slot name (default='de.jaccards')
    #' @param notch boolean Whether to show notches on plot (default=FALSE)
    #' @param show.jitter boolean Whether to show jitter on plots (default=TRUE)
    #' @param jitter.alpha numeric Parameter for jitter (default=0.05)
    #' @param show.pairs boolean Whether to show pairs (default=FALSE)
    #' @param sort.order boolean Whether to show notches in the boxplot version (default=TRUE)
    #' @param pallete plot palette (default=self$cell.groups.palette)
    #' @param set.fill (default=TRUE)
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateDEStability()
    #' cao$plotDEStabilityPerCellType)
    #' }
    plotDEStabilityPerCellType=function(name='de.jaccards', notch=FALSE, show.jitter=TRUE, jitter.alpha = 0.05,
                                        show.pairs = FALSE, sort.order = TRUE, pallete=self$cell.groups.palette,
                                        set.fill=TRUE) {
      jaccards <- private$getResults(name, 'estimateDEStability()')
      p <- plotStability(jaccards = jaccards,
                         notch = notch,
                         show.jitter = show.jitter,
                         jitter.alpha = jitter.alpha,
                         show.pairs = show.pairs,
                         sort.order = sort.order,
                         xlabel = 'Cell Type',
                         ylabel = 'Jaccard Index',
                         palette=NULL,
                         plot.theme=self$plot.theme, set.fill=set.fill) + theme(legend.position = "none")

      p <- p + theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylim(0, 1)  + xlab('')
      if(!is.null(palette) && !(show.pairs) ){
        p <- p + scale_color_manual(values=self$cell.groups.palette)
      }

      if(!set.fill){
        p <- p + scale_fill_manual(values=c('white'))
      }

      return(p)
    },

    #' @description Plot DE stability per gene
    #' @param name character string DE results slot name (default='de')
    #' @param cell.type character If set only show stability for a specific cell type in DE results (default=NULL)
    #' @param stability.score character string Any of "stab.median.rank", "stab.mean.rank", or "stab.var.rank" (default='stab.median.rank')
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateDEStabilityPerGene()
    #' cao$plotDEStabilityPerGene()
    #' }
    plotDEStabilityPerGene=function(name='de', cell.type=NULL, stability.score='stab.median.rank') {
      de.res <- private$getResults(name, 'estimateDEPerCellType()')
      possible.scores <- c('stab.median.rank', 'stab.mean.rank', 'stab.var.rank')
      if ( !(stability.score %in% possible.scores) ) stop('Please provide correct name of the stability core')
      if ( !(cell.type %in% names(de.res$initial)) ) stop('Please provide correct cell type to visualise')
      if ( !(stability.score %in% names(de.res$initial[[cell.type]]$res)) ) stop('Stability score was not estimated')

      idx <- (rank(de.res$initial[[cell.type]]$res$pvalue) < 500)
      p <- smoothScatter(x = rank(de.res$initial[[cell.type]]$res$pvalue)[idx],
                         y = rank(de.res$initial[[cell.type]]$res[[stability.score]])[idx],
                         xlab = 'rank of DE gene', ylab = 'stability score',
                         colramp = grDevices::colorRampPalette(c("white", 'seagreen4')))
      return(p)
    },

    #' @description Plot number of significant DE genes
    #' @param p.adjust boolean Whether the cutoff should be based on the adjusted P value (default=TRUE)
    #' @param pvalue.cutoff numeric P-value cutoff (default=0.05)
    #' @param show.resampling.results boolean Whether to show uncertainty based on resampling results (default=TRUE)
    #' @param show.jitter boolean Whether to apply jitter to the ggplot (default=FALSE)
    #' @param jitter.alpha numeric Opacity setting (default=0.05)
    #' @param type character string Any of 'box', 'point', or 'bar' (default='bar')
    #' @param notch boolean Whether to show notches (default=TRUE)
    #' @param ... additional parameters passed to plotMeanMedValuesPerCellType()
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$plotNumberOfDEGenes()
    #' }
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
          if (length(miss.subsamples) > 0) {
            warning("Subtypes ", paste(miss.subsamples, collapse=","), " are missed form sampling, ignoring those")
          }

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

    #' @description Make volcano plots
    #' @param name character string DE results slot name (default='de')
    #' @param cell.types character If set will plot only for selected cell types in DE results (default=NULL)
    #' @param palette plot palette If NULL will use standard palette (default=NULL)
    #' @param build.panel boolean (default=TRUE)
    #' @param n.col numeric Number of columns (default=3)
    #' @param color.var character string (default='CellFrac')
    #' @param ... additional parameters fed to plotVolcano
    #' @return A ggplot2 object
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$plotVolcano()
    #' }
    plotVolcano=function(name='de', cell.types=NULL, palette=NULL, build.panel=TRUE, n.col=3,
                         color.var = 'CellFrac', ...) {
      de <- private$getResults(name, 'estimateDEPerCellType()') %>% lapply(`[[`, 'res')
      if (is.null(palette)) {
        palette <- c("#5e4fa2", "#3288bd", "#abdda4", "#fdae61", "#f46d43", "#9e0142") %>%
          grDevices::colorRampPalette()
      }

      if(!(color.var %in% names(de[[1]]))) {
        stop(paste(color.var, 'is not calculated'))
      }

      if (is.null(cell.types)) {
        cell.types <- names(de)
      }
      de <- de[intersect(cell.types, names(de))]

      if (length(de) == 0) {
        stop("No cell types left after the filtering")
      }

      if (length(cell.types) == 1) {
        return(plotVolcano(de[[cell.types]], color.var=color.var, palette=palette, plot.theme=self$plot.theme, ...))
      }

      if (!build.panel) {
        return(lapply(de, plotVolcano, color.var=color.var, palette=palette, plot.theme=self$plot.theme, ...))
      }

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
    #' @param saveprefix character Prefix for created files (default=NULL)
    #' @param dir.name character Name for directory with results (default="JSON")
    #' @param de.raw List of DE results. If NULL will use stored DE results defined by "de.name" (default=NULL)
    #' @param de.name character string DE results slot name (default='de')
    #' @param ref.level character Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
    #' @param gene.metadata (default=NULL)
    #' @param verbose boolean Whether to output verbose messages (default=self$verbose)
    #' @return saved JSON objects
    #' @examples
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$saveDEasJSON()
    #' }
    saveDEasJSON=function(saveprefix=NULL, dir.name="JSON", de.raw=NULL, sample.groups=self$sample.groups, de.name='de',
                          ref.level=self$ref.level, gene.metadata=NULL, verbose=self$verbose) {
      if (is.null(de.raw)) {
        de.raw <- private$getResults(de.name, "estimateDEPerCellType")
      }

      if (!is.list(sample.groups)) {
        sample.groups <- list(names(sample.groups[sample.groups == ref.level]),
                              names(sample.groups[sample.groups != ref.level])) %>%
          setNames(c(ref.level, self$target.level))
      }

      if (!("list" %in% class(de.raw[[1]]))) stop("Please rerun 'estimateDEPerCellType' with return.matrix=T")

      if (is.null(saveprefix)) saveprefix <- ""

      saveDEasJSON(de.raw=de.raw, saveprefix=saveprefix, dir.name=dir.name, gene.metadata=gene.metadata,
                   sample.groups=sample.groups, verbose=verbose)
    },

    #' @description Plot embedding
    #' @param embedding A cell embedding to use (two-column data frame with rownames corresponding to cells) (default: stored embedding object)
    #' @param plot.theme plot theme to use (default: `self$plot.theme`)
    #' @param color.by color cells by 'cell.groups', 'condition' or 'sample'. Overrides `groups` and `palette`. (default: NULL)
    #' @param ... other parameters passed to \link[sccore:embeddingPlot]{embeddingPlot}
    #' @return embedding plot as output by sccore::embeddingPlot
    #' @examples 
    #' \dontrun{
    #' cao$plotEmbedding()
    #' }
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

      rlang::exec(sccore::embeddingPlot, !!!params)
    },

    #' @description Estimate ontology terms based on DEs
    #' @param type character Ontology type, either GO (gene ontology) or DO (disease ontology). Please see DOSE package for more information (default="GO")
    #' @param name character If NULL will use `type` to look for ontology results (default=NULL)
    #' @param de.name character string DE results slot name (default='de')
    #' @param org.db Organism database, e.g., org.Hs.eg.db::org.Hs.eg.db for human or org.Ms.eg.db::org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
    #' @param n.top.genes numeric Number of top highest-expressed genes to consider (default=500)
    #' @param p.adj numeric adjust-pvalues cutoff fed to getDEEntrezIdsSplitted() (default: 1)
    #' @param readable boolean Mapping gene ID to gene name (default=TRUE)
    #' @param keep.gene.sets boolean (default=FALSE)
    #' @param ignore.cache (default=NULL)
    #' @param de.raw (default=NULL)
    #' @param min.genes numeric Minimum number of input genes overlapping with ontologies (default=0)
    #' @param qvalue.cutoff numeric Q value cutoff, please see clusterProfiler package for more information (default=0.2)
    #' @param min.gs.size numeric Minimal geneset size, please see clusterProfiler package for more information (default=5)
    #' @param max.gs.size numeric Minimal geneset size, please see clusterProfiler package for more information (default=5e2)
    #' @param ... further argument for ontology estimation. Pass `nPerm` with `type='GSEA'` to use fgseaSimple method
    #' @return A list containing a list of terms per ontology, and a data frame with merged results
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' library(org.Hs.eg.db)
    #' cao$estimateOntology(type = "GSEA", org.db = org.Hs.eg.db)
    #' }
    estimateOntology=function(type=c("GO", "DO", "GSEA"), name=NULL, de.name='de', org.db, n.top.genes=500, p.adj=1,
                              p.adjust.method="BH", readable=TRUE, min.gs.size=10, max.gs.size=500,
                              keep.gene.sets=FALSE, ignore.cache=NULL, de.raw=NULL, verbose=self$verbose,
                              n.cores=self$n.cores, ...) {
      type <- match.arg(type)
      if (is.null(name)) {
        name <- type
      }

      if (is.null(de.raw)) {
        de.raw <- private$getResults(de.name, "estimateDEPerCellType()")
      }
      # If estimateDEPerCellType was run with return.matrix = TRUE, remove matrix before calculating
      if ("list" %in% class(de.raw[[1]])) de.raw %<>% lapply(`[[`, "res")

      de.gene.ids <- getDEEntrezIdsSplitted(de.raw, org.db=org.db, p.adj=p.adj)

      go.environment <- self$getGOEnvironment(org.db, verbose=verbose, ignore.cache=ignore.cache)
      res <- estimateOntologyFromIds(
        de.gene.ids, type=type, org.db=org.db, n.top.genes=n.top.genes, go.environment=go.environment,
        pAdjustMethod=p.adjust.method, readable=readable, minGSSize=min.gs.size, maxGSSize=max.gs.size,
        keep.gene.sets=keep.gene.sets, verbose=verbose, n.cores=n.cores, ...
      )

      self$test.results[[name]] <- list(res=res, de.gene.ids=de.gene.ids, type=type) # redundancy needed
      return(invisible(self$test.results[[name]]))
    },

    #' @description Estimate ontology families based on ontology results
    #' @param name character Type of ontology result: "GO", "GSEA", or "DO" (default="GO")
    #' @param p.adj double Cutoff for adjusted p (default=0.05)
    #' @return List of families and ontology data per cell type
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$estimateOntologyFamilies(name = "GSEA")
    #' }
    estimateOntologyFamilies=function(name="GO", p.adj=0.05) {
      checkPackageInstalled("GOfuncR", bioc=TRUE)
      ont.res <- private$getResults(name, 'estimateOntology()')

      if (ont.res$type == "GO") {
        ont.list <- lapply(ont.res$res, lapply, lapply, function(x) {
          tmp <- x@result %>% filter(p.adjust <= p.adj)
          if (nrow(tmp) > 0) return(tmp)
        }) %>% lapply(lapply, plyr::compact) %>% lapply(plyr::compact)
      } else {
        ont.list <- lapply(ont.res$res, lapply, function(x) {
          tmp <- x@result %>% filter(p.adjust <= p.adj)
          if (nrow(tmp) > 0) return(tmp)
        }) %>% lapply(plyr::compact)
      }
      self$test.results[[name]]$families <- estimateOntologyFamilies(ont.list=ont.list, type=ont.res$type)
      return(invisible(self$test.results[[name]]))
    },

    #' @description Identify families containing a specific ontology term or ID
    #' @param name character string Type of ontology result: "GO", "GSEA", or "DO" (default="GO")
    #' @param go.term character vector with term description(s) (default=NULL)
    #' @param go.id character vector with ID(s) (default=NULL)
    #' @param common boolean Only identify families with all the supplied terms or IDs (default = FALSE)
    #' @return Data frame
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$estimateOntologyFamilies(name = "GSEA")
    #' cao$getFamiliesPerGO(name = "GSEA", go.term = "antigen presentation") # Either go.term og go.id has to be specified 
    #' }
    getFamiliesPerGO=function(name="GO", go.term=NULL, go.id=NULL, common=FALSE) {

      ont.res <- private$getResults(name, 'estimateOntology()')

      if (is.null(go.term) && is.null(go.id)) stop("Please specify either 'go.term' or 'go.id'.")
      if (!is.null(go.term) && !is.null(go.id))
        warning("Both 'go.term' and 'go.id' specified, will only use 'go.term'.")

      # Extract data
      ont.fam <- ont.res$families
      if (is.null(ont.fam)) stop("No family data found.")

      # Get mapping of IDs to Descriptions
      desc.per.id <- rapply(ont.res$res, function(x) x@result %$% setNames(Description, ID), how="list") %>%
        do.call(c, .)

      if (ont.res$type == "GO") desc.per.id %<>% do.call(c, .)
      desc.per.id %<>% unname() %>% do.call(c, .)

      # Make data.frame
      list.levels <- getOntologyListLevels(ont.res$type)
      res <- ont.fam %>% rblapply(list.levels, function(fs) {
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

    #' @description Bar plot of ontology terms per cell type
    #' @param name character string Type of ontology result: "GO", "GSEA", or "DO" (default="GO")
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param p.adj numeric adjusted p-value cutoff (default=0.05)
    #' @param q.value numeric Q value used for filtering (default=0.2)
    #' @param min.genes integer Minimum number of overlapping genes in terms (default=1)
    #' @param families boolean Plot family terms (default=FALSE)
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$plotNumOntologyTermsPerType()
    #' }
    plotNumOntologyTermsPerType=function(name="GO", genes="up", p.adj=0.05, q.value=0.2, min.genes=1, families=FALSE) {
      type <- private$getResults(name, 'estimateOntology()')$type
      
      if (families) {
        tmp <- private$getResults(name, 'estimateOntology()')$families
        
        if (type == "GO") {
          p.df <- tmp %>% 
            lapply(lapply, lapply, '[[', "families") %>% 
            lapply(lapply, lapply, unlist) %>% 
            lapply(lapply, lapply, unique) %>% 
            lapply(lapply, sapply, length) %>%
            lapply(lapply, \(type) data.frame(Direction = names(type), N = unname(type))) %>%
            lapply(bind_rows, .id = "Type") %>% 
            bind_rows(.id = "Group") %>% 
            filter(Direction %in% genes)
        } else {
          p.df <- tmp %>% 
            lapply(lapply, '[[', "families") %>% 
            lapply(lapply, unlist) %>% 
            lapply(lapply, unique)
          
          # Filter based on 'genes'
          if (type == "GSEA" && any(genes != "all")) {
            p.df <- genes %>% 
              lapply(\(gene) {
                ont.df <- private$getOntologyPvalueResults(
                  name=name, genes=gene, p.adj=p.adj, q.value=q.value,
                  min.genes=min.genes
                ) %>% 
                  select(Group, Type, ID) %>% 
                  split(., .$Group) %>% 
                  lapply(\(x) split(x, x$Type) %>% 
                           lapply(pull, ID))
                
                # This is really ugly, sorry
                tmp <- p.df %>%
                  names() %>%
                  lapply(\(ct) {
                    p.df[[ct]] %>%
                      names() %>%
                      lapply(\(type) {
                        p.df[[ct]][[type]] %>%
                          .[. %in% ont.df[[ct]][[type]]]
                      }) %>% 
                      `names<-`(p.df[[ct]] %>% names()) %>%
                      .[sapply(., length) > 0]
                  }) %>% 
                  `names<-`(p.df %>% names()) %>%
                  .[sapply(., length) > 0]
                
                tmp %<>% lapply(sapply, length) %>% 
                  {lapply(names(.), \(x) data.frame(Group = x, Type = names(.[[x]]), N = unname(.[[x]])))} %>% 
                  bind_rows() %>% 
                  mutate(Direction = gene)
                
                return(tmp)
              }) %>%
              bind_rows()
          } else {
            p.df %<>% 
              lapply(sapply, length) %>% 
              {lapply(names(.), \(x) data.frame(Group = x, Type = names(.[[x]]), N = unname(.[[x]])))} %>% 
              bind_rows()
          }
          
          if (type == "DO") {
            p.df %<>%
              rename(Direction = Type) %>% 
              filter(Direction %in% genes)
          }
        }
      } else {
        if (length(genes) > 1) {
          ont.res <- genes %>% setNames(., .) %>% lapply(function(g) {
            private$getOntologyPvalueResults(name=name, gene=g, p.adj=p.adj, q.value=q.value, min.genes=min.genes)
          })
          
          classes <- sapply(ont.res[genes], class)
          if (any(classes == "character")) {
            message("No significant results found for genes = '",
                    paste(names(classes[classes == "character"]), collapse=','), "'.")
            genes <- names(classes[classes == "data.frame"])
            if (length(genes) == 0) stop("No results to plot.")
            
            ont.res %<>% .[genes]
          }
          
          ont.res %<>% names() %>% 
            lapply(function(d) dplyr::mutate(ont.res[[d]], Direction=d)) %>%
            Reduce(rbind, .)
        } else {
          ont.res <- private$getOntologyPvalueResults(
            name=name, genes=genes, p.adj=p.adj, q.value=q.value, min.genes=min.genes
          ) %>% dplyr::mutate(Direction=genes)
        }
        
        # Prepare data further
        if (type %in% c("GO", "GSEA")) {
          p.df <- ont.res %$% table(Group=Group, Type=Type, Direction=Direction) %>% as.data.frame(responseName='N')
        } else if (type == "DO") {
          p.df <- ont.res %$% table(Group=Group, Direction=Direction) %>% as.data.frame(responseName='N')
        } else stop("Unexpected type ", type)
      }
      
      # Create plot
      if (type %in% c("GO", "GSEA")) {
        gg <- ggplot(p.df, aes(x=Group, y=N, fill=Type, group=Group)) +
          geom_bar(stat="identity")
        if (length(unique(p.df$Direction)) > 1) {
          gg <- gg + facet_grid(~Direction, switch="x")
        }
      } else if (type == "DO") {
        if (length(unique(p.df$Direction)) > 1) {
          gg <- ggplot(p.df) +
            geom_bar(aes(x=Group, y=N, fill=Direction), stat="identity", position="dodge") +
            labs(fill="Gene set")
        } else {
          gg <- ggplot(p.df) +
            geom_bar(aes(x=Group, y=N), stat="identity")
        }
      }
      
      if (families) type %<>% paste0(.," family")
      
      gg <- gg +
        scale_y_continuous(expand=c(0, 0, 0.05, 0)) +
        self$plot.theme +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position="right",
              axis.title.x=element_blank(), panel.grid.major.x=element_blank()) +
        labs(x="", y=paste0("No. of ", type, " terms"))
      return(gg)
    },

    #' @description Plot a dotplot of ontology terms with adj. P values for a specific cell subgroup
    #' @param cell.type character string Cell type to plot
    #' @param name character string Type of ontology result: "GO", "GSEA", or "DO" (default="GO")
    #' @param plot character string Type of plot to return (default="dot"). Either "dot" or "bar".
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param subtype character string Ontology type, must be either "BP", "CC", or "MF" (GO types), "GO" or "DO" (default="GO")
    #' @param cell.subgroup character Specific cell group to plot
    #' @param n integer Number of ontology terms to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
    #' @param p.adj numeric Adjusted P cutoff (default=0.05)
    #' @param min.genes integer Minimum overlapping genes per term (default=1)
    #' @param ... additional parameters passed to enrichplot::dotplot() or enrichplot::barplot()
    #' @return A ggplot2 object
    #' @examples
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$plotOntology(name = "GSEA", cell.type = "Neurons") # "cell.type" is a cell type in self$cell.groups used for calculating ontologies
    #' }
    plotOntology = function(cell.type, name="GO", plot="dot", genes=c("up", "down", "all"),
                            subtype=c("BP", "CC", "MF"), n=20, p.adj=0.05, min.genes=1, ...) {
      # Checks
      checkPackageInstalled("enrichplot", bioc=TRUE)
      subtype <- match.arg(subtype)
      genes <- match.arg(genes)

      ont.res <- private$getResults(name, 'estimateOntology()')
      type <- ont.res$type
      ont.res <- ont.res$res
      if ((type == "GSEA") && (plot == "bar")) stop("No 'enrichplot' method exists for making barplots of GSEA results.")

      if (plot != "dot" && plot != "bar") {
        stop(paste("Unknown plot type: ", plot, ". The plot parameter must be specified as either 'dot' or 'bar'."))
      }
      
      # Extract results
      if (!cell.type %in% names(ont.res)) stop("'cell.type' not found in results.")
      ont.res %<>% .[[cell.type]]
      if (type != "DO") ont.res %<>% .[[subtype]]
      if (is.null(ont.res)){
        stop("No results found for ", name, ", ", subtype, " for ", cell.type)
      }

      if (type %in% c("GO", "DO")) {
        ont.res %<>% .[[genes]]
      } else {
        if (genes == "up") ont.res %<>% filter(NES >= 0) else if (genes == "down") ont.res %<>% filter(NES <= 0)
      } 
      if (is.null(ont.res)){
        stop("No results found for ", genes, " genes for ", name, ", ", subtype, " for ", cell.type)
      }

      # Prepare data
      df <- ont.res@result %>% filter(p.adjust <= p.adj)
      if (nrow(df) == 0){
        stop("Nothing to plot. Try relaxing 'p.adj'. The lowest adj. P value is ",
             formatC(min(ont.res@result$p.adjust), digits=3))
      }

      # Allow plotting of terms with p.adj > 0.05
      if (p.adj > 0.05) {
        if (type == "GSEA") {
          ont.res@params$pvalueCutoff <- 1
        } else {
          ont.res@pvalueCutoff <- 1
          ont.res@qvalueCutoff <- 1
        }
      }

      if (min.genes > 1 && type != "GSEA") { # GeneRatio is not provided for GSEA
        idx <- df$GeneRatio %>%
          strsplit("/", fixed=TRUE) %>%
          sapply(`[[`, 1)

        df <- df[idx > min.genes,]
      }
      ont.res@result <- df

      # Plot
      if (plot == "dot"){
        return(enrichplot::dotplot(ont.res, showCategory=n, orderBy="x", ...))
      } else if (plot == "bar"){
        return(enrichplot::barplot(ont.res, showCategory=n, ...))
      } else {
        stop("Unknown plot type: ", plot)
      }

    },

    #' @description Plot a heatmap of ontology P values per cell type
    #' @param genes Specify which genes to plot, can either be 'down' for downregulated genes, 'up' or 'all'
    #'   (default="up")
    #' @param subtype character string (default="BP")
    #' @param q.value numeric (default=0.2)
    #' @param min.genes integer Minimum genes (default=1)
    #' @param top.n Number of terms to show (default=Inf)
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default="GO")
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="left")
    #' @param selection Order of rows in heatmap. Can be 'unique' (only show terms that are unique for any cell type);
    #'   'common' (only show terms that are present in at least two cell types); 'all' (all ontology terms)
    #'   (default="all")
    #' @param cluster boolean Whether to show GO clusters or raw GOs (default=TRUE)
    #' @param clust.naming Field with the results for GO clustering. Ignored if `clusters == FALSE`.
    #' @param cell.subgroups Cell groups to plot (default=NULL). This affects only visualization, but not clustering.
    #' @param color.range vector with two values for min/max values of p-values
    #' @param row.order boolean Whether to order rows (default=TRUE)
    #' @param col.order boolean Whether to order columns (default=TRUE)
    #' @param max.log.p numeric Maximum log P value, used for coloring (default=10)
    #' @param only.family.children boolean Whether to only include family children (default=FALSE)
    #' @param description.regex (default=NULL)
    #' @param description.exclude.regex (default=NULL)
    #' @param readjust.p boolean Whether to re-adjust p-values (default=TRUE)
    #' @param p.adjust.method character string Method used to adjust p-values (default="BH")
    #' @param palette plot palette. If NULL default will be used (default=NULL)
    #' @param return.info boolean (default=FALSE)
    #' @param ... parameters forwarded to \link{plotHeatmap}
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$plotOntologyHeatmap()
    #' 
    #' cao$estimateOntologyFamilies(name = "GSEA")
    #' cao$plotOntologyHeatmap(name = "GSEA", only.family.children = TRUE)
    #' }
    plotOntologyHeatmap=function(name="GO", genes="up", subtype="BP", p.adj=0.05, q.value=0.2, min.genes=1, top.n=Inf,
                                 legend.position="left", selection="all", cluster=TRUE, cell.subgroups=NULL,
                                 row.order=TRUE, col.order=TRUE, max.log.p=10, only.family.children=FALSE,
                                 description.regex=NULL, description.exclude.regex=NULL, clust.naming="medoid",
                                 readjust.p=TRUE, p.adjust.method="BH",
                                 palette=NULL, color.range=NULL, return.info=FALSE, ...) {
      ont.info <- private$getOntologyHeatmapInfo(
        name=name, genes=genes, subtype=subtype, p.adj=p.adj, q.value=q.value, min.genes=min.genes,
        selection=selection, cluster=cluster, cell.subgroups=cell.subgroups, only.family.children=only.family.children,
        clust.naming=clust.naming, description.regex=description.regex,
        description.exclude.regex=description.exclude.regex, readjust.p=readjust.p, p.adjust.method=p.adjust.method
      )

      ont.sum <- ont.info$ont.sum
      if (is.null(ont.sum)) return(NULL)
      if (is.null(palette)) {
        palette <- getGenePalette(genes, bidirectional=(prod(range(ont.sum, na.rm=TRUE)) < 0))
      }

      ont.sum[abs(ont.sum) > max.log.p] %<>% {max.log.p * sign(.)}
      ont.sum %<>% .[order(rowSums(abs(.)), decreasing=TRUE),,drop=FALSE] %>% head(top.n)

      gg <- plotHeatmap(
        ont.sum, legend.position=legend.position, row.order=row.order, col.order=col.order,
        plot.theme=self$plot.theme, palette=palette, color.range=color.range, ...
      )

      if (return.info) {
        ont.info$ont.sum <- ont.sum
        ont.info$gg <- gg
        return(ont.info)
      }

      return(gg)
    },

    #' @description Plot a heatmap of ontology p-values per cell type heatmap, collapsed
    #' @param name character string Type of ontology result: "GO", "GSEA", or "DO" (default="GO")
    #' @param genes character Specify which genes to plot, can either be 'down' for downregulated genes, 'up' or 'all'
    #'   (default="up")
    #' @param subtype character Ontology type, must be either "BP", "CC", or "MF" (GO types) or "DO" (default="GO")
    #' @param p.adj numeric Adj. P value cutoff (default=0.05)
    #' @param q.value numeric Q value cutoff (default=0.2)
    #' @param min.genes integer Minimum number of overlapping genes per term (default=1)
    #' @param n integer Number of terms to plot (default=20)
    #' @param legend.position character Position of legend in plot. See ggplot2::theme (default="left")
    #' @param selection character Order of rows in heatmap. Can be 'unique' (only show terms that are unique for any cell type);
    #'   'common' (only show terms that are present in at least two cell types); 'all' (all ontology terms)
    #'   (default="all")
    #' @param max.log.p numeric Maximum log P value, used for coloring (default=10)
    #' @param top.n Number of terms to show (default=Inf)
    #' @param cluster Whether to show GO clusters or raw GOs (default=TRUE)
    #' @param clust.naming Field with the results for GO clustering. Ignored if `clusters == FALSE`.
    #' @param cell.subgroups Cell groups to plot (default=NULL). This affects only visualization, but not clustering.
    #' @param color.range vector with two values for min/max values of p-values
    #' @param row.order boolean Whether to order rows (default=TRUE)
    #' @param col.order boolean Whether to order columns (default=TRUE)
    #' @param n.words integer (default=5)
    #' @param exclude.words (default=NULL)
    #' @param return.info boolean Whether to return the info (default=FALSE)
    #' @param palette (default=NULL)
    #' @param only.family.children boolean (default=FALSE)
    #' @param distance character string (default="manhattan")
    #' @param readjust.p boolean Whether to re-adjust p-values (default=TRUE)
    #' @param clust.method character string (default="complete")
    #' @param description.regex (default=NULL)
    #' @param description.exclude.regex (default=NULL)
    #' @param ... parameters forwarded to \link{plotHeatmap}
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$plotOntologyHeatmapCollapsed(name = "GSEA")
    #' 
    #' cao$estimateOntologyFamilies(name = "GSEA")
    #' cao$plotOntologyHeatmapCollapsed(name = "GSEA", only.family.children = TRUE)
    #' }
    plotOntologyHeatmapCollapsed=function(name="GO", genes="up", subtype="BP", p.adj=0.05, q.value=0.2, min.genes=1,
                                          n=20, legend.position="left", selection="all", max.log.p=10, cluster=TRUE,
                                          cell.subgroups=NULL, palette=NULL, row.order=TRUE, col.order=TRUE,
                                          only.family.children=FALSE, readjust.p=TRUE, p.adjust.method='BH',
                                          description.regex=NULL, description.exclude.regex=NULL,
                                          distance="manhattan", clust.method="complete", clust.naming="consensus",
                                          n.words=5, exclude.words=NULL, return.info=FALSE, ...) {
      ont.info <- private$getOntologyHeatmapInfo(
        name=name, genes=genes, subtype=subtype, p.adj=p.adj, q.value=q.value, min.genes=min.genes,
        selection=selection, cluster=cluster, cell.subgroups=cell.subgroups, only.family.children=only.family.children,
        clust.naming=clust.naming, description.regex=description.regex,
        description.exclude.regex=description.exclude.regex, readjust.p=readjust.p, p.adjust.method=p.adjust.method
      )

      desc.per.clust <- ont.info$desc.per.clust
      ont.sum <- ont.info$ont.sum

      if (is.null(ont.sum)) return(NULL)
      if (is.null(palette)) {
        palette <- getGenePalette(genes, bidirectional=(prod(range(ont.sum, na.rm=TRUE)) < 0))
      }

      ont.sum.raw <- ont.sum
      ont.sum[abs(ont.sum) > max.log.p] %<>% 
        {max.log.p * sign(.)}

      gos.per.clust <- dist(ont.sum, method=distance) %>%
        hclust(method=clust.method) 
      
      # Adding check for n to avoid errors
      max.gpc <- max(gos.per.clust$order)
      if (max.gpc < n) {
        warning(paste0("Reducing 'n' to ",max.gpc))
        n <- max.gpc
      }
      
      gos.per.clust %<>% 
        cutree(k=n) %>% 
        {split(names(.), .)}
      
      clust.names <- sapply(gos.per.clust, function(gos) {
        n.gos <- length(gos)
        if (n.gos == 1) return(paste("1: ", gos))
        if (cluster) gos <- unlist(desc.per.clust[gos], use.names=FALSE)

        estimateOntologyClusterName(gos, method=clust.naming, n.words=n.words, exclude.words=exclude.words) %>%
          {paste0(n.gos, ": ", .)}
      })

      ont.sum <- lapply(gos.per.clust, function(ns) colMeans(ont.sum.raw[ns,,drop=FALSE])) %>%
        do.call(rbind, .) %>% 
        as.data.frame() %>% 
        set_rownames(clust.names)
      ont.sum[abs(ont.sum) > max.log.p] %<>% 
        {max.log.p * sign(.)}

      ont.freqs <- gos.per.clust %>%
        lapply(function(gos) colMeans(abs(ont.sum.raw[gos,]) > 1e-5) * 100) %>%
        do.call(rbind, .) %>% 
        as.data.frame()

      gg <- plotHeatmap(
        ont.sum, size.df=ont.freqs, legend.position=legend.position, row.order=row.order, col.order=col.order,
        plot.theme=self$plot.theme, palette=palette, distance=distance, clust.method=clust.method,
        size.legend.title="Cluster GOs %", ...
      )

      if (return.info) {
        ont.info$ont.sum <- ont.sum
        ont.info$ont.freqs <- ont.freqs
        ont.info$gg <- gg
        return(ont.info)
      }

      return(gg)
    },

    #' @description Plot correlation matrix for ontology terms between cell types
    #' @param name character string Type of ontology result: "GO", "GSEA", or "DO" (default="GO")
    #' @param subtype character Type of ontology result, must be "BP", "MF", or "CC" (default="BP")
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param p.adj numeric Cut-off for adjusted p-values (default=0.05)
    #' @param only.family.children boolean Plot similarities for ontology family lonely children (default=FALSE)
    #' @param q.value numeric Q value for filtering (default=0.2)
    #' @param min.genes numeric Minimum number of overlapping genes per term (default=1)
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology()
    #' cao$plotOntologySimilarities()
    #' }
    plotOntologySimilarities=function(name="GO", subtype=c("BP","MF","CC"), genes="up", p.adj=0.05, only.family.children=FALSE, q.value=0.2, min.genes=1) {
      subtype <- match.arg(subtype)
      
      if (only.family.children) {
        fams <- private$getResults(name, 'estimateOntology()')$families
        if (is.null(fams))
          stop("No ontology family results found, please run 'estimateOntologyFamilies' first, or set only.family.children=FALSE")
      }
      
      ont.res <- private$getOntologyPvalueResults(
        name=name, genes=genes, p.adj=p.adj, q.value=q.value, min.genes=min.genes, subtype=subtype
      )
      
      if (nrow(ont.res) == 0) {
        stop("No ontologies found for name=", name, ", subtype=", subtype, " and genes=", genes,". You could also consider relaxing p.adj.")
      }
      
      type <- private$getResults(name, 'estimateOntology()')$type

      if ((ont.res$Group %>% unique() %>% length()) == 1)
        stop("Only one group present, correlation cannot be performed.")

      if (only.family.children) {
        ont.res %<>% getOntologyFamilyChildren(fams=fams, subtype=subtype, genes=genes, type=private$getResults(name, 'estimateOntology()')$type)
        if (nrow(ont.res) == 0) {
          stop("No ontology family children found.")
        }
      }
      
      if (type %in% c("GO", "GSEA")) {
        pathway.df <- ont.res[c("Description", "Group", "Type")] %>% rename(Pathway=Description, GO=Type)
      } else if (type=="DO") {
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

      # TODO: currently we use binary distance. Probably, checking z-scores would give better results.
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

    #' @description Plot related ontologies in one hierarchical network plot
    #' @param name character string Type of ontology result: "GO", "GSEA", or "DO" (default="GO")
    #' @param cell.type character Cell subtype to plot
    #' @param family numeric Family within cell subtype to plot (default=NULL)
    #' @param genes character string Only for GO results: Direction of genes, must be "up", "down", or "all" (default="up")
    #' @param subtype character Type of ontology result, must be "BP", "MF", or "CC" (default="BP")
    #' @param plot.type character Extend of family network Can be "complete" (entire network), "dense" (show 1 parent for each significant term), or "minimal" (only show significant terms) (default="complete")
    #' @param show.ids boolean Whether to show ontology IDs instead of names (default=FALSE)
    #' @param string.length integer Length of strings for wrapping in order to fit text within boxes (default: 14)
    #' @param legend.label.size numeric Size og legend labels (default: 1)
    #' @param legend.position numeric Position of legend (default: topright)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param ... additional parameters passed to plotOntologyFamily
    #' @return Rgraphviz object
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$estimateOntologyFamilies(name = "GSEA")
    #' cao$plotOntologyFamily(name = "GSEA", cell.type = "Neurons") # "cell.type" is a cell type in self$cell.groups used for calculating ontologies
    #' }
    plotOntologyFamily=function(name="GO", cell.type, family=NULL, genes="up", subtype="BP",
                                plot.type="complete", show.ids=FALSE, string.length=14, legend.label.size=1,
                                legend.position="topright", verbose=self$verbose, n.cores=self$n.cores, ...) {
      #Checks
      checkPackageInstalled(c("GOfuncR", "graph", "Rgraphviz"), bioc=TRUE)

      plot.type <- match.arg(plot.type)
      
      ont.res <- private$getResults(name, 'estimateOntology()')
      ont.fam.res <- ont.res$families
      if (is.null(ont.fam.res))
        stop("No results found for '", name, "'. Please run 'estimateOntologyFamilies' first.")

      ont.fam.res %<>% .[[cell.type]]
      if (is.null(ont.fam.res)) stop("No results found for cell.type '", cell.type, "'.")
      
      if (ont.res$type != "DO") {
        ont.fam.res %<>% .[[subtype]]
        if (is.null(ont.fam.res)) stop("No results found for subtype '", subtype, "'.")
      }

      if (ont.res$type != "GSEA") {
        ont.fam.res %<>% .[[genes]]
        if (is.null(ont.fam.res)) stop("No results found for genes '", genes, "'.")
      }

      if (is.null(family)) {
        fam.names <- names(ont.fam.res$families)
      } else {
        if (!is.numeric(family)) {
          fam.names <- family
        } else {
          fam.names <- paste0("Family", family)
        }
      }

      if (!all(fam.names %in% names(ont.fam.res$families))) stop("Not all families are found in 'ont.fam.res'.")
      families <- lapply(fam.names, function(n) ont.fam.res$families[[n]]) %>% unlist() %>% unique()

      plotOntologyFamily(fam=families, data=ont.fam.res$data, plot.type=plot.type, show.ids=show.ids,
                         string.length=string.length, legend.label.size=legend.label.size,
                         legend.position=legend.position, verbose=verbose, n.cores=n.cores, ...)
    },

    #' @description Save ontology results as a table
    #' @param file character string File name passed to write.table(). Set to NULL to return the table instead of saving.
    #' @param subtype character string Only for GO results: Type of result to filter by, must be "BP", "MF", or "CC" (default: NULL)
    #' @param genes character Direction of genes to filter by, must be "up", "down", or "all" (default: NULL)
    #' @param sep character Separator (default: "\t", tab)
    #' @param ... additional arguments passed to write.table()
    #' @return table for import into text editor
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$saveOntologyAsTable(name = "GSEA", file = "Ontologies.tsv")
    #' }
    saveOntologyAsTable=function(file, name="GO", subtype=NULL, genes=NULL, p.adj=0.05, sep="\t", ...) {
      ont.res <- private$getResults(name, 'estimateOntology()')

      list.levels <- getOntologyListLevels(ont.res$type)
      res <- rblapply(ont.res$res, list.levels, function(r) r@result) %>%
        filter(p.adjust <= p.adj) %>% as_tibble()

      if (!is.null(subtype)) res %<>% filter(Subtype %in% subtype)
      if (!is.null(genes)) res %<>% filter(Genes %in% genes)

      if (is.null(file)) return(res)
      write.table(res, file=file, sep=sep, row.names=FALSE, ...)
    },

    #' @description Save family results as a table
    #' @param file character string File name passed to write.table(). Set to NULL to return the table instead of saving.
    #' @param type character string Type of ontology result, i.e., GO, GSEA, or DO (default='GO')
    #' @param subtype character Type of result to filter by, must be "BP", "MF", or "CC" (default=NULL)
    #' @param genes character Direction of genes to filter by, must be "up", "down", or "all" (default=NULL)
    #' @param p.adj numeric Adjusted P to filter by (default=0.05)
    #' @param sep character Separator (default = "\t", tab)
    #' @param ... additional arguments passed to write.table()
    #' @return table for import into text editor
    #' @examples 
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$estimateOntologyFamilies(name = "GSEA")
    #' cao$saveFamiliesAsTable(name = "GSEA", file = "Families.tsv")
    #' }
    saveFamiliesAsTable=function(file, name="GO", subtype=NULL, genes=NULL, p.adj=0.05, sep="\t", ...) {
      # Extract results
      ont.res <- private$getResults(name, 'estimateOntology()')
      ont.fam.res <- ont.res$families
      if (is.null(ont.fam.res))
        stop("No results found for '", name, "'. Please run 'estimateOntologyFamilies' first.")

      list.levels <- getOntologyListLevels(ont.res$type)
      res <- rblapply(ont.fam.res, list.levels, function(x) {
        lapply(x$families, function(y) {
          if (length(y) <= 1) return(NULL) # TODO: What about == 1?
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
        }) %>% .[!sapply(., is.null)] %>% bind_rows(.id='Family')
      })

      # Filtering
      if (!is.null(subtype)) res %<>% filter(subtype %in% subtype)
      if (!is.null(genes)) res %<>% filter(genes %in% genes)

      # Write table
      if (is.null(file)) return(res)
      write.table(res, file=file, sep=sep, row.names=FALSE, ...)
    },

    #' @description Plot the cell group sizes or proportions per sample
    #' @param cell.groups factor Cell annotations with cell IDs as names (default=self$cell.groups)
    #' @param palette color palette to use for conditions (default: stored $sample.groups.palette)
    #' @param show.significance boolean show statistical significance between sample groups. wilcox.test was used; (`*` < 0.05; `**` < 0.01; `***` < 0.001) (default=FALSE)
    #' @param filter.empty.cell.types boolean Remove cell types without cells (default=TRUE)
    #' @param proportions boolean Plot proportions or absolute numbers (default=TRUE)
    #' @param ... additional plot parameters, forwarded to \link{plotCountBoxplotsPerType}
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$plotCellGroupSizes()
    #' }
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

    #' @description Plot the cell group sizes or proportions per sample
    #' @param cell.groups character Cell annotations with cell IDs as names(default=self$cell.groups)
    #' @param type character string Must be "mad", "sd", "sample.num", or "sample.frac" (default='mad')
    #' @param rotate.xticks boolean Turn x labels 90 degrees (default=TRUE)
    #' @param min.rel.abundance numeric Minimum relative abundance to plot (default=0.05)
    #' @return ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$plotCellGroupAbundanceVariation()
    #' }
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
    #' @param space character either 'PCA' or 'CDA' (default="CDA")
    #' @param font.size numeric Font size (default=3)
    #' @param cells.to.remain character Specific cell types to keep (default=NULL)
    #' @param cells.to.remove character Specific cell types to remove (default=NULL)
    #' @param samples.to.remove character Specific samples to remove (default=NULL)
    #' @param palette plot palette (default=self$sample.groups.palette)
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellLoadings()
    #' cao$plotCodaSpace()
    #' }
    plotCodaSpace=function(space='CDA', cell.groups=self$cell.groups, font.size=3,
                            cells.to.remain=NULL, cells.to.remove=NULL,
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
    #' @param cell.groups character Cell annotations with cell IDs as name (default=self$cell.groups)
    #' @param palette plot palette (default=self$sample.groups.palette)
    #' @param name character Results name slot (default='coda')
    #' @param cells.to.remain character Specific cell types to keep (default=NULL)
    #' @param cells.to.remove character Specific cell types to remove (default=NULL)
    #' @param filter.empty.cell.types boolean Remove cell types without cells (default=TRUE)
    #' @param adjust.pvalues boolean Adjust P values or not (default=TRUE)
    #' @param h.method character Must be one of 'both', 'up', 'down' (default='both')
    #' @param reorder.tree boolean Reorder tree or not (default=TRUE)
    #' @param ... additional parameters
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellLoadings()
    #' cao$plotContratsTree()
    #' }
    plotContrastTree=function(cell.groups=self$cell.groups, palette=self$sample.groups.palette, name='coda',
                              cells.to.remain=NULL, cells.to.remove=NULL, filter.empty.cell.types=TRUE,
                              adjust.pvalues=TRUE, h.method=c('both', 'up', 'down'), reorder.tree=TRUE, ...) {
      h.method <- match.arg(h.method)
      tmp <- private$extractCodaData(cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain,
                                     cell.groups=cell.groups)
      if (filter.empty.cell.types) {
        cell.type.to.remain <- (colSums(tmp$d.counts[tmp$d.groups,]) > 0) &
          (colSums(tmp$d.counts[!tmp$d.groups,]) > 0)
        tmp$d.counts <- tmp$d.counts[,cell.type.to.remain]
      }

      loadings.mean <- NULL
      if (reorder.tree) {
        if (name %in% names(self$test.results)) {
          loadings.mean <- rowMeans(self$test.results[[name]]$loadings)
        } else {
          message('To show significance, please run estimateCellLoadings()')
        }
      }

      gg <- plotContrastTree(tmp$d.counts, tmp$d.groups, self$ref.level, self$target.level,
                             plot.theme=self$plot.theme, adjust.pvalues=adjust.pvalues,
                             h.method=h.method, loadings.mean=loadings.mean, palette=palette, ...)

      return(gg)
    },


    #' @description Plot composition similarity
    #' @param cell.groups character Cell annotations with cell IDs as name (default=self$cell.groups)
    #' @param cells.to.remain character Specific cell types to keep (default=NULL)
    #' @param cells.to.remove character Specific cell types to remove (default=NULL)
    #' @param palette plot palette (default=brewerPalette("YlOrRd", rev=FALSE))
    #' @param ... parameters passed to plotHeatmap()
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellLoadings()
    #' cao$plotCompositionSimilarity()
    #' }
    plotCompositionSimilarity=function(cell.groups=self$cell.groups, cells.to.remain=NULL, cells.to.remove=NULL,
                                       palette=brewerPalette("YlOrRd", rev=FALSE), ...) {
      tmp <- private$extractCodaData(
        cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain, cell.groups=cell.groups
      )
      mat <- referenceSet(tmp$d.counts, tmp$d.groups)$mx.first
      diag(mat) <- 1
      gg <- plotHeatmap(mat, legend.title="Similarity", palette=palette, ...)
      return(gg)
    },

    #' @description Estimate cell loadings
    #' @param n.boot numeric Number of boot straps (default=1000)
    #' @param ref.cell.type character Reference cell type (default=NULL)
    #' @param name character Results name slot (default='coda')
    #' @param n.seed numeric Seed number for reproducibility (default=239)
    #' @param cells.to.remove character Specific cell types to keep (default=NULL)
    #' @param cells.to.remain character Specific cell types to remove (default=NULL)
    #' @param samples.to.remove character Specific samples to remove (default=NULL)
    #' @param filter.empty.cell.types boolean Remove cell types without cells (default=TRUE)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @return resulting cell loadings
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellLoadings()
    #' }
    estimateCellLoadings=function(n.boot=1000, ref.cell.type=NULL, name='coda', n.seed=239,
                                  cells.to.remove=NULL, cells.to.remain=NULL, samples.to.remove=NULL,
                                  filter.empty.cell.types=TRUE, n.cores=self$n.cores, verbose=self$verbose) {
      checkPackageInstalled(c("coda.base", "psych"), cran=TRUE)

      if ((!is.null(ref.cell.type)) && (!(ref.cell.type %in% levels(self$cell.groups))))
        stop('Incorrect reference cell type')

      # Get cell counts and groups
      tmp <- private$extractCodaData(cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain,
                                     samples.to.remove=samples.to.remove)

      if (filter.empty.cell.types) {
        cell.type.to.remain <- (colSums(tmp$d.counts[tmp$d.groups,]) > 0) &
          (colSums(tmp$d.counts[!tmp$d.groups,]) > 0)
        tmp$d.counts <- tmp$d.counts[,cell.type.to.remain]
      }
      cnts <- tmp$d.counts
      groups <- tmp$d.groups

      res <- runCoda(cnts, groups, n.boot=n.boot, n.seed=n.seed, ref.cell.type=ref.cell.type)
      res$cnts <- cnts
      res$groups <- groups

      ## Calculate normalized counts
      ref.cell.type <- res$ref.cell.type

      ref.cnts <- cnts[, ref.cell.type, drop=FALSE]
      ref.cnts[ref.cnts == 0] <- 0.5
      norm.val <- 1 / nrow(ref.cnts) * rowSums(log(ref.cnts))
      cnts.nonzero <- cnts
      cnts.nonzero[cnts.nonzero == 0] <- 0.5
      res$norm.cnts <- log(cnts.nonzero) - norm.val

      self$test.results[[name]] <- res

      return(invisible(res))
    },

    #' @description Plot Loadings
    #' @param alpha numeric Transparency (default=0.01)
    #' @param palette plot palette specification for cell types (default: stored $cell.groups.palette)
    #' @param font.size numeric Font size (default=NULL)
    #' @param name character Results slot name (default='coda')
    #' @param ordering character Must be one of "pvalue", "loadings"  (default='pvalue')
    #' @param show.pvals boolean Show P values (default=TRUE)
    #' @param ... additional parameters plotCellLoadings()
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellLoadings()
    #' cao$plotCellLoadings()
    #' }
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

    ### Cluster-free cell density

    #' @description Estimate cell density in giving embedding
    #' @param bins numeric Number of bins for density estimation (default=400)
    #' @param method character string Density estimation method, graph: graph smooth based density estimation. kde: embedding grid based density  estimation. (default: 'kde')
    #' @param beta numeric Smoothing strength parameter of the \link[sccore:heatFilter]{heatFilter} for graph based cell density (default=30)
    #' @param estimate.variation boolean Estimate variation (default=TRUE)
    #' @param sample.groups 2-factor vector with annotation of groups/condition per sample (default=self$sample.groups)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param bandwidth numeric KDE bandwidth multiplier (default=0.5). The full bandwidth is estimated by multiplying this value on the difference between 90% and 10%
    #' of the corresponding embedding dimension. Set it to NULL to use \link[MASS:bandwidth.nrd]{bandwidth.nrd} estimator. (default=0.05)
    #' @param ... additional arguments
    #' @return estimated cell densities
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellDensity()
    #' }
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
    #' @param show.grid boolean Whether to show grid (default=TRUE)
    #' @param add.points boolean Add points to cell density figure (default=TRUE)
    #' @param size numeric Point size (default=0.1)
    #' @param show.legend boolean Show legend (default=FALSE)
    #' @param palette plot palette (default=NULL)
    #' @param point.col character Point color (default='#313695')
    #' @param contours character Specify cell types for contour, multiple cell types are also supported (default=NULL)
    #' @param contour.color character Color for contour line (default='black')
    #' @param contour.conf character Confidence interval of contour  (default='10%')
    #' @param name character Slot in which to saved results from estimateCellDensity (default='cell.density')
    #' @param show.cell.groups boolean Plot cell group names (default=TRUE)
    #' @param cell.groups character Cell annotations with cell IDs as name (default=self$cell.groups)
    #' @param font.size numeric Font size (default=c(2, 4))
    #' @param color.range character Color range (default=c(0, "99%"))
    #' @param ... plot style parameters forwarded to \link[sccore:styleEmbeddingPlot]{sccore::styleEmbeddingPlot}.
    #' @return A ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellDensity()
    #' cao$plotCellDensity()
    #' }
    plotCellDensity = function(show.grid=TRUE, add.points=TRUE, size=0.1, show.legend=FALSE, palette=NULL,
                               point.col='#313695', contours=NULL, contour.color='black', contour.conf='10%',
                               name='cell.density', show.cell.groups=TRUE, cell.groups=self$cell.groups,
                               font.size=c(2, 4), color.range=c(0, "99%"), ...) {
      dens.res <- private$getResults(name, 'estimateCellDensity()')
      private$checkCellEmbedding()

      if (is.null(palette)) {
        palette <- c("#FFFFFF", brewerPalette("YlOrRd", rev=FALSE)(9)) %>% grDevices::colorRampPalette()
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

    #' @description Plot cell density variation
    #' @param type character Must be one of "mad", "sd", "sample.frac" (default='mad')
    #' @param plot.type character Must be one of "hist", "embedding" (default='embedding')
    #' @param name character Results slot name (default='cell.density')
    #' @param cutoff numeric Score cutoff (default=NULL)
    #' @param condition character Must be one of 'both', 'ref', 'target' (default="both")
    #' @param ... additional arguments
    #' @return ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellDensity(estimate.variation=TRUE)
    #' cao$plotCellDensityVariation()
    #' }
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

      if (is.null(scores)) {
        stop("To use this function, please re-run estimateCellDensity() with estimate.variation=TRUE")
      }

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

    #' @description Estimate differential cell density
    #' @param type character method to calculate differential cell density; permutation, t.test, wilcox or subtract (target subtract ref density);
    #' @param adjust.pvalues boolean Whether to adjust Z-scores for multiple comparison using BH method (default: FALSE for type='subtract', TRUE for everything else)
    #' @param name character Slot with results from estimateCellDensity. New results will be appended there. (Default: 'cell.density')
    #' @param n.permutations numeric Number of permutations (default=400)
    #' @param smooth boolean Smooth results (default=TRUE)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param ... additional arguments to the function
    #' @return estimated differential cell densities
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellDensity()
    #' cao$estimateDiffCellDensity()
    #' }
    estimateDiffCellDensity=function(type='permutation', adjust.pvalues=NULL, name='cell.density',
                                     n.permutations=400, smooth=TRUE, verbose=self$verbose, n.cores=self$n.cores, ...){
      dens.res <- private$getResults(name, 'estimateCellDensity')
      if (is.null(adjust.pvalues)) adjust.pvalues <- (type != 'subtract') # NULL can be forwarded here
      density.mat <- dens.res$density.mat
      if (dens.res$method == 'kde'){
        if (adjust.pvalues) {
          # For p-value adjustment matrix filtration can result in disconnected grid and break estimating l.max
          l.max <- rownames(density.mat) %>% as.integer() %>% graphFromGrid(n.bins=dens.res$bins) %>%
            igraph::laplacian_matrix(sparse=TRUE) %>% irlba::partial_eigen(n=1) %>% .$values
        }

        density.mat <- density.mat[dens.res$density.emb$counts > 0,]
        graph <- rownames(density.mat) %>% as.integer() %>% graphFromGrid(n.bins=dens.res$bins)
      } else if (adjust.pvalues && smooth) {
        graph <- extractCellGraph(self$data.object)
        l.max <- NULL
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
            score, permut.scores, smooth=smooth, graph=graph, n.cores=n.cores, verbose=verbose,
            l.max=l.max, ...
          )
        )
      }

      self$test.results[[name]]$diff[[type]] <- res

      return(invisible(self$test.results[[name]]))
    },

    #' @description Estimate differential cell density
    #' @param type character method to calculate differential cell density; t.test, wilcox or subtract (target subtract ref density);
    #' @param name character Slot with results from estimateCellDensity. New results will be appended there. (Default: 'cell.density')
    #' @param size numeric (default=0.2)
    #' @param palette color palette, default is c('blue','white','red')
    #' @param adjust.pvalues boolean Adjust P values (default=NULL)
    #' @param contours character Specify cell types for contour, multiple cell types are also supported (default: NULL)
    #' @param contour.color character color for contour line (default: 'black')
    #' @param contour.conf character confidence interval of contour (default: '10%')
    #' @param plot.na boolean Plot NAs (default=FALSE)
    #' @param color.range numeric, e.g. c(0,90) (default=NULL)
    #' @param mid.color character Color code for medium value in color range (default='gray95')
    #' @param scale.z.palette boolean Scale plot palette for Z scores (default=adjust.pvalues)
    #' @param min.z numeric Minimum Z score to plot (default=qnorm(0.9))
    #' @param ... additional parameters
    #' @return ggplot2 object
    #' @examples 
    #' \dontrun{
    #' cao$estimateCellDensity()
    #' cao$estimateDiffCellDensity()
    #' cao$plotDiffCellDensity()
    #' }
    plotDiffCellDensity=function(type=NULL, name='cell.density', size=0.2, palette=NULL,
                                 adjust.pvalues=NULL, contours=NULL, contour.color='black', contour.conf='10%',
                                 plot.na=FALSE, color.range=NULL, mid.color='gray95',
                                 scale.z.palette=adjust.pvalues, min.z=qnorm(0.9), ...) {
      if (is.null(palette)) {
        if (is.null(self$sample.groups.palette)) {
          palette <- c('blue', mid.color, 'red')
        } else {
          palette <- self$sample.groups.palette %>% {c(.[self$ref.level], mid.color, .[self$target.level])}
        }
        palette %<>% grDevices::colorRampPalette(space="Lab")
      }
      private$checkCellEmbedding()
      dens.res <- private$getResults(name, 'estimateCellDensity')

      if (is.null(type)) {
        type <- if (length(dens.res$diff) == 0) 'permutation' else names(dens.res$diff)[1]
      }
      scores <- dens.res$diff[[type]]
      if (is.null(scores)) {
        warning("Can't find results for name, '", name, "' and type '", type,
                "'. Running estimateDiffCellDensity with default parameters.")
        self$estimateDiffCellDensity(type=type, name=name, adjust.pvalues=adjust.pvalues)
        dens.res <- self$test.results[[name]]
        scores <- dens.res$diff[[type]]
      }

      if (is.null(adjust.pvalues)) {
        if (!is.null(scores$adj)) {
          scores <- scores$adj
          adjust.pvalues <- TRUE
        } else {
          scores <- scores$raw
          adjust.pvalues <- FALSE
        }
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

      if (dens.res$method == 'graph') {
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


    #' @description Plot inter-sample expression distance. The inputs to this function are the results from cao$estimateExpressionShiftMagnitudes()
    #' @param name character Test results to plot (default=expression.shifts)
    #' @param joint boolean Whether to show joint boxplot with the expression distance weighed by the sizes of cell types (default: TRUE), or show distances for each individual cell type
    #' @param palette plot palette (default=self$sample.groups.palette)
    #' @param show.significance boolean Whether to show statistical significance between sample groups. wilcox.test was used; (`*` < 0.05; `**` < 0.01; `***` < 0.001)
    #' @param ... other plot parameters, forwarded to \link{plotCountBoxplotsPerType}
    #' @return A ggplot2 object
    #' @examples
    #' \dontrun{
    #' cao$estimateExpressionShiftMagnitudes()
    #' cao$plotExpressionDistance()
    #' }
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

    #' @description Plot inter-sample expression distance. The inputs to this function are the results from cao$estimateExpressionShiftMagnitudes()
    #' @param space character One of 'expression.shifts', 'coda', 'pseudo.bulk' (default="expression.shifts")
    #' @param cell.type character Cell type reference for distancing (default=NULL)
    #' @param dist character Must be one of "cor", "l1" (manhattan), "l2" (euclidian) (default=NULL)
    #' @param name character Results slot name (default=NULL)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param ... additional arguments
    #' @return sample distance matrix
    #' @examples
    #' \dontrun{
    #' cao$getSampleDistanceMatrix()
    #' }
    getSampleDistanceMatrix=function(space=c('expression.shifts', 'coda', 'pseudo.bulk'), cell.type=NULL,
                                     dist=NULL, name=NULL, verbose=self$verbose, sample.subset=NULL, ...) {
      space <- match.arg(space)
      if ((space != 'pseudo.bulk') && (length(list(...)) > 0)) stop("Unexpected arguments: ", names(list(...)))

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
      } else if (space=='coda') {
        n.cells.per.samp <- table(self$sample.per.cell)
        mat <- private$extractCodaData() %$% getRndBalances(d.counts) %$% prcomp(norm) %$% as.data.frame(x)

        dist %<>% parseDistance(top.n.genes=ncol(mat), n.pcs=NULL)
        if (dist == 'cor') {
          p.dists <- 1 - cor(t(mat))
        } else if (dist == 'l2') {
          p.dists <- dist(mat, method="euclidean") %>% as.matrix()
        } else if (dist == 'l1') {
          p.dists <- dist(mat, method="manhattan") %>% as.matrix()
        } else {
          stop("Unknown distance: ", dist)
        }
      } else {
        stop("Not implemented space: ", space, "!")
      }

      if (!is.null(sample.subset)) {
        p.dists <- p.dists[sample.subset, sample.subset]
      }

      return(p.dists)
    },

    #' @description Project samples to 2D space with MDS. Plots results from cao$estimateExpressionShiftMagnitudes() or cao$estimateCellLoadings()
    #' @param space character string "expression.shifts" Results from cao$estimateExpressionShiftMagnitudes(); CDA- cell composition shifts result from cao$estimateCellLoadings(); sudo.bulk- expression distance of sudo bulk
    #' @param method character string "MDS"
    #' @param dist 'cor' - correlation distance, 'l1' - manhattan distance or 'l2' - euclidean (default correlation distance)
    #' @param cell.type If a name of a cell type is specified, the sample distances will be assessed based on this cell type alone. Otherwise (cell.type=NULL, default), sample distances will be estimated as an average distance across all cell types (weighted by the minimum number of cells of that cell type between any two samples being compared)
    #' @param palette a set of colors to use for conditions (default: stored $sample.groups.palette)
    #' @param show.sample.size make point size proportional to the log10 of the number of cells per sample (default: FALSE)
    #' @param sample.colors (default=NULL)
    #' @param color.title (default=NULL)
    #' @param title (default=NULL)
    #' @param n.permutations numeric (default=2000)
    #' @param show.pvalues boolean (default=FALSE)
    #' @param ... additional parameters passed to plotSampleDistanceMatrix()
    #' @return A ggplot2 object
    #' @examples
    #' \dontrun{
    #' cao$estimateExpressionShiftMagnitudes()
    #' cao$plotSampleDistances()
    #' }
    plotSampleDistances=function(space='expression.shifts', method='MDS', dist=NULL, name=NULL, cell.type=NULL,
                                 palette=NULL, show.sample.size=FALSE, sample.colors=NULL, color.title=NULL,
                                 title=NULL, n.permutations=2000, show.pvalues=FALSE, sample.subset=NULL,
                                 n.cores=self$n.cores, ...) {
      if (is.null(cell.type)) {
        n.cells.per.samp <- table(self$sample.per.cell)
      } else {
        if (is.null(title)) title <- cell.type
        n.cells.per.samp <- self$sample.per.cell %>% .[self$cell.groups[names(.)] == cell.type] %>% table()
      }

      p.dists <- self$getSampleDistanceMatrix(
        space=space, cell.type=cell.type, dist=dist, name=name, sample.subset=sample.subset
      )
      if (is.null(p.dists)) return(NULL)

      if (is.null(sample.colors) && is.null(palette)) {
        # Has to be in the same order, or ggplot separates shape and color legends into two
        palette <- self$sample.groups.palette[levels(as.factor(self$sample.groups))]
      }
      gg <- plotSampleDistanceMatrix(
        p.dists=p.dists, sample.groups=self$sample.groups, n.cells.per.samp=n.cells.per.samp, method=method,
        sample.colors=sample.colors, show.sample.size=show.sample.size, palette=palette, color.title=color.title,
        title=title, plot.theme=self$plot.theme, ...
      )

      return(gg)
    },

    #' @description Estimate metadata separation using variance on the sample distance graph
    #' @param sample.meta sample metadata is a list or data.frame with metadata per sample
    #' @param space (default="expression shifts")
    #' @param dist (default=NULL)
    #' @param space.name (default=NULL)
    #' @param n.permutations number permutations for the test (default=5000)
    #' @param trim trim distance matrix above the given quantile (default=0.05)
    #' @param k if this parameter is supplied, k-NN graph is used for variance estimation, otherwise
    #' the function uses a fully-connected graph (default=20)
    #' @param show.warning boolean (default=TRUE)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param adjust.pvalues boolean (default=TRUE)
    #' @param pvalue.cutoff numeric (default=0.05)
    #' @return results
    #' @examples 
    #' \dontrun{
    #' cao$estimateExpressionShiftMagnitudes() # or estimateCellLoadings()
    #' cao$estimateMetadataSeparation(sample.meta = meta.data) # meta.data is a list or data.frame with metadata per sample
    #' }
    estimateMetadataSeparation=function(sample.meta, space='expression.shifts', dist=NULL, space.name=NULL,
                                        sample.subset=NULL,
                                        name='metadata.separation', n.permutations=5000, trim=0.05, k=20,
                                        show.warning=TRUE, verbose=self$verbose, n.cores=self$n.cores,
                                        adjust.pvalues=TRUE, p.adjust.method="BH", pvalue.cutoff=0.05) {
      p.dists <- self$getSampleDistanceMatrix(
        space=space, cell.type=NULL, dist=dist, name=space.name, sample.subset=sample.subset
      )
      # Check whether results are empty
      if (is.null(p.dists)) {
        warning("An empty sample distance matrix was returned. Consider changing 'space', 'cell.type', 'name',  or 'sample.subset'.")
        return(NULL)
      }

      # Check whether any sample names are present in sample.meta and p.dists
      if (!any(rownames(sample.meta) %in% rownames(p.dists))) {
        cat("Printing the first three rownames of sample.meta:\n",head(rownames(sample.meta), 3),"\nPrinting the first three sample names:\n",head(rownames(p.dists), 3),"\n"); stop("The rownames of the sample.meta object doesn't match any sample names.")
      } 
      
      if (is.data.frame(sample.meta)) {
        sample.meta %<>% lapply(setNames, rownames(.))
      } else if (!is.list(sample.meta)) {
        sample.meta %<>% list()
      }

      adj.mat <- adjacencyMatrixFromPaiwiseDists(p.dists, trim=trim, k=k)
      sep.info <- sample.meta %>% plapply(
        function(mg) estimateGraphVarianceSignificance(adj.mat, signal=mg[colnames(adj.mat)]),
        progress=(verbose && (length(sample.meta) > 1)), n.cores=n.cores, mc.preschedule=TRUE, fail.on.error=TRUE
      )

      pvalues <- sapply(sep.info, `[[`, 'pvalue')
      pseudo.r2 <- sapply(sep.info, `[[`, 'pr2')

      res <- list(metadata=sample.meta, pvalues=pvalues, pseudo.r2=pseudo.r2)
      if (adjust.pvalues) {
        pvalues %<>% p.adjust(method=p.adjust.method)
        res$padjust <- pvalues
      }

      if (any(is.na(pvalues))) warning(paste0(paste(names(pvalues)[is.na(pvalues)], sep = "\t")," resulted in NAs when calculating p values. Is the metadata defined for all samples?"))
      
      if (show.warning && any(pvalues < pvalue.cutoff, na.rm = TRUE)){
        warning("Significant separation by: ", paste(names(pvalues)[pvalues < pvalue.cutoff], collapse=', '))
      }

      self$test.results[[name]] <- res
      return(invisible(res))
    },

    #' @description Plot metadata separation
    #' @param name character Name for storage in test.results (default="metadata.separation")
    #' @param pvalue.y numeric (default=0.93)
    #' @param ... additional parameters forwarded to \link[plotMeanMedValuesPerCellType]{plotMeanMedValuesPerCellType}
    #' @examples 
    #' \dontrun{
    #' cao$estimateExpressionShiftMagnitudes() # or estimateCellLoadings()
    #' cao$estimateMetadataSeparation(sample.meta = meta.data)
    #' cao$plotMetadataSeparation()
    #' }
    plotMetadataSeparation=function(name='metadata.separation', pvalue.y=0.93, ...) {
      res <- private$getResults(name, "estimateMetadataSeparation()")
      pvals <- if (is.null(res$padjust)) res$pvalues else res$padjust
      gg <- res$pseudo.r2 %>% {tibble(Type=names(.), value=ifelse(is.na(.), 0, .))} %>%
        plotMeanMedValuesPerCellType(type="bar", ylab=expression(Pseudo-R^2), jitter.alpha=0, pvalues=pvals,
                                     pvalue.y=pvalue.y, ...) +
        scale_y_continuous(expand=c(0, 0)) +
        scale_fill_manual(values=rep("#2b8cbe", length(pvals)))
      return(gg)
    },

    #' @description Estimate differential expression Z-scores between two conditions per individual cell
    #' @param n.top.genes (default=Inf)
    #' @param genes (default=NULL)
    #' @param max.z z-score value to winsorize the estimates for reducing impact of outliers. Default: 20.
    #' @param min.expr.frac minimal fraction of cell expressing a gene for estimating z-scores for it. Default: 0.001.
    #' @param min.n.samp.per.cond minimal number of samples per condition for estimating z-scores (default: 2)
    #' @param min.n.obs.per.samp minimal number of cells per samples for estimating z-scores (default: 2)
    #' @param robust whether to use median estimates instead of mean. Using median is more robust,
    #' but greatly increase the number of zeros in the data, leading to bias towards highly-express genes. (Default: FALSE)
    #' @param norm.both boolean (default=TRUE)
    #' @param adjust.pvalues boolean (default=FALSE)
    #' @param smoooth boolean (default=TRUE)
    #' @param wins numeric (default=0.01)
    #' @param n.permutations numeric (default=200)
    #' @param lfc.pseudocount pseudocount value for estimation of log2(fold-change)
    #' @param min.edge.weight numeric (default=0.6)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param name character string (default='cluster.free.de')
    #' @param smooth boolean Whether to apply smoothing (default=TRUE)
    #' @return list with sparce matrices containing various DE metrics with genes as columns and cells as rows:
    #'   - `z`: DE Z-scores
    #'   - `reference` mean or median expression in reference samples
    #'   - `target` mean or median expression in target samples
    #'   - `lfc`: log2(fold-change) of expression
    #' Cells that have only one condition in their expression neighborhood have NA Z-scores for all genes.
    #' Results are also stored in the `cluster.free.de` field.
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' }
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

    #' @description Get most changed genes
    #' @param n numeric Number of genes to retrieve
    #' @param method character Must be one of "z", "z.adj", "lfc" (default="z")
    #' @param min.z numeric Minimum Z score (default=0.5)
    #' @param min.lfc numeric Minimum log fold change (default=1)
    #' @param max.score numeric Maximum Z score (default=20)
    #' @param cell.subset character Cells to subset (default=NULL)
    #' @param excluded.genes character Genes to exclude (default=NULL)
    #' @param included.genes character Genes to include (default=NULL)
    #' @param name character Results slot name (default="cluster.free.de")
    #' @return named numeric with scores and gene symbols as names
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$getMostChangedGenes(n = 10) # n can be any number of genes to extract
    #' }
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
    #' @param min.expr.frac numeric (default=0.0)
    #' @param min.n.obs.per.samp minimal number of cells per sample for using it in distance estimation (default: 3)
    #' @param normalize.both whether to normalize results relative to distances within both conditions (TRUE) or only to the control (FALSE)
    #' @param dist distance measure. Options: "cor" (correlation), "cosine" or "js" (Jensen–Shannon)
    #' @param log.vectors whether to use log10 on the normalized expression before estimating the distance.
    #' In most cases, must be TRUE for "cosine" and "cor" distances and always must be FALSE for "js". (default: `dist != 'js'`)
    #' @param wins numeric (default=0.025)
    #' @param genes character vector Genes to include (default=NULL)
    #' @param n.permutations numeric Number of permutations (default=500)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param min.edge.weight numeric Minimum edge weight (default=0.0)
    #' @param ... additional parameters passed to estimateClusterFreeExpressionShiftsC()
    #' @return Vector of cluster-free expression shifts per cell. Values above 1 correspond to difference between conditions.
    #' Results are also stored in the `cluster.free.expr.shifts` field.
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$estimateClusterFreeExpressionShifts()
    #' }
    estimateClusterFreeExpressionShifts=function(n.top.genes=3000, gene.selection="z", name="cluster.free.expr.shifts",
                                                 min.n.between=2, min.n.within=max(min.n.between, 1),
                                                 min.expr.frac=0.0, min.n.obs.per.samp=3, normalize.both=FALSE,
                                                 dist="cor", log.vectors=(dist != "js"), wins=0.025, genes=NULL,
                                                 n.permutations=500, verbose=self$verbose, n.cores=self$n.cores,
                                                 min.edge.weight=0.0, ...) {
      if (normalize.both)
        warning("Setting normalize.both=TRUE likely leads to wrong results for cluster-free shifts")

      if (is.null(genes)) {
        genes <- private$getTopGenes(n.top.genes, gene.selection=gene.selection, min.expr.frac=min.expr.frac)
      }

      inp <- private$getClusterFreeDEInput(genes, min.edge.weight=min.edge.weight)
      shifts <- inp %$% estimateClusterFreeExpressionShiftsC(
        t(cm), self$sample.per.cell[rownames(cm)], nn_ids=nns.per.cell, is_ref=is.ref, min_n_between=min.n.between,
        min_n_within=min.n.within, min_n_obs_per_samp=min.n.obs.per.samp, norm_all=normalize.both, verbose=verbose,
        n_cores=n.cores, dist=dist, log_vecs=log.vectors, wins=wins, n_permutations=n.permutations, ...
      )
      self$test.results[[name]] <- shifts

      return(invisible(shifts))
    },

    #' @description Performs graph smoothing of the cluster-free DE Z-scores
    #' @param n.top.genes numeric Number of top ranked genes to include (default=1000)
    #' @param smoothing `beta` parameter of the \link[sccore:heatFilter]{heatFilter}. (default=20)
    #' @param filter graph filter function. (default=NULL)
    #' @param z.adj boolean Adjust Z scores (default=FALSE)
    #' @param gene.selection character Must be one of "z.adj" or "z", default is based on the "z.adj" parameter (default=ifelse(z.adj, "z.adj", "z"))
    #' @param exluded.genes character Genes to exclude (default=NULL)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param name character Results slot name (default='cluster.free.de')
    #' @param ... parameters forwarded to \link[sccore:smoothSignalOnGraph]{smoothSignalOnGraph}
    #' @return Sparse matrix of smoothed Z-scores. Results are also stored in the `cluster.free.de$z.smoothed` field.
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$smoothClusterFreeZScores()
    #' }
    smoothClusterFreeZScores = function(n.top.genes=1000, smoothing=20, filter=NULL, z.adj=FALSE, gene.selection=ifelse(z.adj, "z.adj", "z"),
                                        excluded.genes=NULL, n.cores=self$n.cores, verbose=self$verbose, name="cluster.free.de", ...) {
      z.scores <- private$getResults(name, "estimateClusterFreeDE")
      z.scores <- if (z.adj) z.scores$z.adj else z.scores$z

      genes <- private$getTopGenes(n.top.genes, gene.selection=gene.selection,
                                   excluded.genes=excluded.genes, included.genes=colnames(z.scores))
      z.scores <- z.scores[,genes]

      if (verbose) message("Smoothing Z-scores for ", ncol(z.scores), " genes passed filtration")
      if (is.null(filter)) {
        filter <- function(...) sccore::heatFilter(..., beta=smoothing)
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

    #' @description Estimate Gene Programs based on cluster-free Z-scores
    #' @param method character String Method to use (default=c("pam", "leiden", "fabia"))
    #' @param n.top.genes (default=Inf)
    #' @param genes (default=NULL)
    #' @param n.programs maximal number of gene programs to find (parameter `p` for fabia). (default=15)
    #' @param z.adj boolean (default=FALSE)
    #' @param smooth boolean (default=TRUE)
    #' @param abs.scores boolean (default=FALSE)
    #' @param cell.subset (default=NULL)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param max.z numeric (default=5)
    #' @param min.z numeric (default=0.5)
    #' @param min.change.frac numeric (default=0.01)
    #' @param de.name character string (default="cluster.free.de")
    #' @param ... keyword arguments forwarded to estimateGenePrograms
    #' @return a list includes:
    #'   - `fabia`: \link[fabia:Factorization]{fabia::Factorization} object, result of the
    #'       \link[fabia:fabia]{fabia::fabia} call
    #'   - `sample.ids`: ids of the subsampled cells used for fabia estimates
    #'   - `scores.exact`: vector of fabia estimates of gene program scores per cell. Estimated only for the
    #'     subsampled cells.
    #'   - `scores.approx`: vector of approximate gene program scores, estimated for all cells in the dataset
    #'   - `loadings`: matrix with fabia gene loadings per program
    #'   - `gene.scores`: list of vectors of gene scores per program. Contains only genes, selected for
    #'     the program using fabia biclustering.
    #'   - `bi.clusts` fabia biclustering information, result of the \link[fabia:extractBic]{fabia::extractBic} call
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$estimateGenePrograms()
    #' }
    estimateGenePrograms = function(method=c("pam", "leiden", "fabia"), n.top.genes=Inf, genes=NULL, n.programs=15,
                                    z.adj=FALSE, gene.selection=ifelse(z.adj, "z.adj", "z"), smooth=TRUE,
                                    abs.scores=FALSE, name="gene.programs", cell.subset=NULL, n.cores=self$n.cores,
                                    verbose=self$verbose, max.z=5, min.z=0.5, min.change.frac=0.01,
                                    de.name="cluster.free.de", ...) {
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

    #' @description Plot gene program scores
    #' @param prog.ids (default=NULL)
    #' @param build.panel boolean (default=TRUE)
    #' @param nrow (default=NULL)
    #' @param adj.list (default=NULL)
    #' @param legend.title character string (default="Score")
    #' @param palette (default=NULL)
    #' @param min.genes.per.prog numeric (default=10)
    #' @param color.range (default=c("0.5%", "99.5%"))
    #' @param ... additional parameters
    #' @return gene program scores
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$estimateGenePrograms()
    #' cao$plotGeneProgramScores()
    #' }
    plotGeneProgramScores=function(name="gene.programs", prog.ids=NULL, build.panel=TRUE, nrow=NULL,
                                    adj.list=NULL, legend.title="Score", palette=NULL, min.genes.per.prog=10,
                                    color.range=c("0.5%", "99.5%"), ...) {
      gene.progs <- private$getResults(name, "estimateGenePrograms")
      if (gene.progs$method == "fabia")
        stop("fabia is deprecated and plotting is not supported anymore")

      if (is.null(prog.ids)) {
        prog.ids <- (sapply(gene.progs$genes.per.clust, length) > min.genes.per.prog) %>% which()
      }

      if (is.null(palette) && all(gene.progs$program.scores >= 0)) palette <- dark.red.palette

      ggs <- lapply(prog.ids, function(i) {
        title <- paste0("Program ", i, ". ", length(gene.progs$genes.per.clust[[i]]), " genes.")
        cur.scores <- gene.progs$program.scores[i,]
        cur.range <- parseLimitRange(color.range, cur.scores)
        cur.scores %<>% pmax(cur.range[1]) %>% pmin(cur.range[2])
        gg <- self$plotEmbedding(colors=cur.scores, title=title, legend.title=legend.title, palette=palette, ...) +
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



    #' @description Plot gene program genes
    #' @param program.id program id
    #' @param ordering character vector (default=c("similarity", "loading"))
    #' @param max.genes integer (default=9)
    #' @param plots character string (default="z.adj")
    #' @param ... additional parameters passed to plotGeneExpressionComparison()
    #' @return plotGeneExpressionComparison
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$estimateGenePrograms()
    #' cao$plotGeneProgramGenes(program.id = 1) # program.id is any gene program ID in 1:cao$test.results$gene.programs$n.progs
    #' }
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
    #' @param cell.groups Indicates cell groups with cell names. Set to NULL if it shouldn't be shown. (default: stored vector)
    #' @param smooth boolean (default=TRUE)
    #' @param plot.na boolean (default=FALSE)
    #' @param scale.z.palette boolean (default=TRUE)
    #' @param min.z (default=qnorm(0.9))
    #' @param color.range (default=c("0", "97.5%"))
    #' @param alpha numeric (default=0.2)
    #' @param font.size size range for cell type labels
    #' @param adj.list (default=NULL)
    #' @param palette (default=brewerPalette("YlOrRd", rev=FALSE))
    #' @param build.panel boolean (default=TRUE)
    #' @param ... parameters forwarded to \link[sccore:embeddingPlot]{embeddingPlot}
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$estimateClusterFreeExpressionShifts()
    #' cao$plotClusterFreeExpressionShifts()
    #' }
    plotClusterFreeExpressionShifts = function(cell.groups=self$cell.groups, smooth=TRUE, plot.na=FALSE,
                                               name="cluster.free.expr.shifts", scale.z.palette=TRUE, min.z=qnorm(0.9),
                                               color.range=c("0", "97.5%"), alpha=0.2, font.size=c(3, 5), adj.list=NULL,
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

      if (scale.z.palette) {
        ggs[[2]]$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
        max.score <- max(z.scores, na.rm=TRUE)
        ggs[[2]] <- ggs[[2]] + getScaledZGradient(min.z=min.z, palette=palette, color.range=max.score)
      }

      if (!is.null(cell.groups)) {
        ggs %<>% lapply(transferLabelLayer, self$plotEmbedding(groups=cell.groups), font.size=font.size)
      }

      if (!is.null(adj.list)) ggs %<>% lapply(`+`, adj.list)
      if (build.panel) ggs %<>% cowplot::plot_grid(plotlist=., ncol=2, labels=c("Shifts", "Adj. z-scores"))

      return(ggs)
    },

    #' @description Plot most changed genes
    #' @param n.top.genes numeric
    #' @param method character string (default='z')
    #' @param min.z numeric (default=0.5)
    #' @param min.lfc numeric (default=1)
    #' @param max.score numeric (default=20)
    #' @param cell.subset (default=NULL)
    #' @param ... additional parameters input to self$plotGeneExpressionComparison()
    #' @return plot of the most changed genes via plotGeneExpressionComparison()
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$plotMostChangedGenes(n.top.genes = 10) # n.top.genes is any number of genes to plot
    #' }
    plotMostChangedGenes = function(n.top.genes, method="z", min.z=0.5, min.lfc=1, max.score=20, cell.subset=NULL, excluded.genes=NULL, ...) {
      scores <- self$getMostChangedGenes(n.top.genes, method=method, min.z=min.z, min.lfc=min.lfc, max.score=max.score,
                                         cell.subset=cell.subset, excluded.genes=excluded.genes)
      self$plotGeneExpressionComparison(scores=scores, cell.subset=cell.subset, ...)
    },

    #' @description Plot gene expression comparison
    #' @param genes (default=NULL)
    #' @param scores (default=NULL)
    #' @param max.expr character string (default="97.5%")
    #' @param plots (default=c("z.adj", "z", "expression"))
    #' @param min.z (default=qname(0.9))
    #' @param max.z numeric (default=4)
    #' @param max.z.adj (default=NULL)
    #' @param max.lfc numeric (default=3)
    #' @param smoothed boolean (default=FALSE)
    #' @param gene.palette (default=dark.red.palette)
    #' @param z.palette (default=NULL)
    #' @param z.adj.palette (default=z.palette)
    #' @param lfc.palette (default=NULL)
    #' @param scale.z.palette boolean (default=TRUE)
    #' @param plot.na (default=-1)
    #' @param adj.list (default=NULL)
    #' @param build.panel boolean (default=TRUE)
    #' @param nrow (default=1)
    #' @param cell.subset (default=NULL)
    #' @param groups (default=NULL)
    #' @param subgroups (default=NULL)
    #' @param keep.limits (default=NULL)
    #' @param ... additional parameters
    #' @return list
    #' @examples 
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$plotGeneExpressionComparison()
    #' }
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

    #' @description Get condition per cell
    #' @return conditions per cell
    #' @examples 
    #' \dontrun{
    #' cao$getConditionPerCell()
    #' }
    getConditionPerCell=function() {
      self$sample.per.cell %>%
        {setNames(as.character(self$sample.groups[as.character(.)]), names(.))} %>%
        as.factor()
    },

    #' @description Get joint count matrix
    #' @param force boolean, if TRUE the joint count matrix will be recalculated even though it already exists in self$cache (default=FALSE)
    #' @param raw boolean, return raw counts (default=TRUE)
    #' @return joint count matrix
    #' @examples 
    #' \dontrun{
    #' cao$getJointCountMatrix()
    #' }
    getJointCountMatrix=function(force=FALSE, raw=TRUE) {
      cache.name <- if (raw) "joint.count.matrix" else "joint.count.matrix.norm"
      if (force || is.null(self$cache[[cache.name]])) {
        self$cache[[cache.name]] <- extractJointCountMatrix(self$data.object, raw=raw)
      }

      return(self$cache[[cache.name]])
    },

    #' @description Get GO environment
    #' @param org.db object of class OrgDB from Bioconductor (e.g. org.Hs.eg.db::org.Hs.eg.db)
    #' @param verbose boolean, print progress (default=FALSE)
    #' @param ignore.cache ignore GO environments already in self$cache (default=NULL)
    #' @return GO environment
    #' @examples
    #' \dontrun{
    #' cao$getGOEnvironment(org.db = org.Hs.eg.db::org.Hs.eg.db)
    #' }
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

      if (!("OrgDb" %in% class(org.db))){
        stop("'org.db' must be of class 'OrgDb'. Please input an organism database.")
      }

      get_GO_data <- utils::getFromNamespace("get_GO_data", "clusterProfiler")
      self$cache$go.environment <- c("BP", "CC", "MF") %>% sn() %>%
        plapply(function(n) get_GO_data(org.db, n, "ENTREZID") %>%
                  as.list() %>% as.environment(), n.cores=1, progress=verbose)
      return(self$cache$go.environment)
    }
  ),

  private = list(
    getResults=function(name, suggested.function=NULL) {
      if (!is.null(self$test.results[[name]])){
        return(self$test.results[[name]])
      }

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

    getOntologyPvalueResults=function(name, genes=c('up', 'down', 'all'), readjust.p=FALSE, p.adjust.method='BH',
                                      p.adj=0.05, q.value=0.2, min.genes=1, subtype=NULL, cell.subgroups=NULL) {
      genes <- match.arg(genes)
      ont.res <- private$getResults(name, 'estimateOntology()')
      type <- ont.res$type
      ont.res %<>% .$res

      if (type == "DO") {
        ont.res <- lapply(ont.res, `[[`, genes) %>%
          lapply(function(x) if (length(x)) x@result else x) %>%
          plyr::compact() %>%
          lapply(mutate, Type='DO')
      } else if (type == "GO") {
        ont.res <- lapply(ont.res, lapply, `[[`, genes)
      }

      if (type %in% c("GO", "GSEA")) {
        ont.res %<>%
          lapply(lapply, function(x) if (length(x)) x@result else x) %>%
          lapply(plyr::compact) %>%
          lapply(bind_rows, .id='Type')
      }

      ont.res %<>% bind_rows(.id='Group')

      if (type == "GSEA") {
        if ((nrow(ont.res) == 0) && !('core_enrichment' %in% colnames(ont.res))) {
          # Sometimes, empty results from GSEA don't have core_enrichment (don't know when, see issue #21)
          ont.res$core_enrichment <- character()
        }
        ont.res %<>% rename(geneID=core_enrichment)
      }

      if (readjust.p) {
        ont.res$p.adjust <- p.adjust(ont.res$pvalue, method=p.adjust.method)
      }

      if ((type == "GSEA") && (genes != "all")) {
        ont.res %<>% {if (genes == "up") filter(., enrichmentScore > 0) else filter(., enrichmentScore < 0)}
      }

      ont.res %<>% filterOntologyDf(
        p.adj=p.adj, q.value=q.value, min.genes=min.genes, subtype=subtype, cell.subgroups=cell.subgroups
      )

      return(ont.res)
    },

    getOntologyHeatmapInfo=function(name="GO", genes="up", subtype="BP", p.adj=0.05, q.value=0.2, min.genes=1,
                                    selection=c("all", "common", "unique"), cell.subgroups=NULL,
                                    only.family.children=FALSE, description.regex=NULL, description.exclude.regex=NULL,
                                    cluster=TRUE, clust.naming="medoid", readjust.p=TRUE, p.adjust.method='BH') {
      # Checks
      selection <- match.arg(selection)
      if (!is.null(cell.subgroups) && (length(cell.subgroups) == 1))
        stop("'cell.subgroups' must contain at least two groups. Please use plotOntology instead.")

      if (only.family.children) {
        fams <- private$getResults(name, 'estimateOntology()')$families
        if (is.null(fams))
          stop("No ontology family results found, please run 'estimateOntologyFamilies' first, or set only.family.children=FALSE")
      }

      # Extract results
      ont.df <- private$getOntologyPvalueResults(
        name=name, genes=genes, subtype=subtype, cell.subgroups=cell.subgroups, p.adj=p.adj, q.value=q.value,
        min.genes=min.genes, readjust.p=readjust.p, p.adjust.method=p.adjust.method
      )

      if (nrow(ont.df) == 0) {
        warning("No ontologies found for name=", name, ", subtype=", subtype, " and genes=", genes)
        return(NULL)
      }

      if (only.family.children) {
        ont.df %<>% getOntologyFamilyChildren(fams=fams, subtype=subtype, genes=genes, type=self$test.results[[name]]$type)
      }

      if (!is.null(description.regex)) ont.df %<>% .[grep(description.regex, .$Description),]
      if (!is.null(description.exclude.regex)) {
        ont.df %<>% .[grep(description.exclude.regex, .$Description, invert=TRUE),]
      }

      if (nrow(ont.df) == 0) {
        warning("No ontologies pass the regex filtration for name=", name, ", subtype=", subtype, " and genes=", genes)
        return(NULL)
      }

      desc.per.clust <- NULL
      group.field <- "Description"
      if (cluster) {
        ont.df %<>% clusterOntologyDF(clust.naming=clust.naming) %>% .$df
        desc.per.clust <- ont.df %$% split(Description, ClusterName) %>% lapply(unique)
        group.field <- "ClusterName"
      }

      sign.field <- if ((genes == "all") && !is.null(ont.df$enrichmentScore)) 'enrichmentScore' else 'p.adjust'
      ont.sum <- -groupOntologiesByCluster(ont.df, field=group.field, sign.field=sign.field)

      if (nrow(ont.sum) == 0) {
        warning("No ontologies pass the filtration for name=", name, ", subtype=", subtype, " and genes=", genes)
        return(NULL)
      }

      if (selection=="unique") {
        ont.sum %<>% .[rowSums(abs(.) > 0) == 1,,drop=FALSE]
      } else if (selection=="common") {
        ont.sum %<>% .[rowSums(abs(.) > 0) > 1,,drop=FALSE]
      }

      if (nrow(ont.sum) == 0) {
        warning("Nothing to plot. Try another selection.")
        return(NULL)
      }

      ont.sum %<>% .[, colSums(abs(.)) > 0, drop=FALSE]
      return(list(ont.sum=ont.sum, ont.df=ont.df, desc.per.clust=desc.per.clust, group.field=group.field))
    },

    ## Extract coda data function
    ##
    ## ret.groups boolean Whether to return groups (default=TRUE). If FASLE, returns table of cell counts.
    ## cells.to.remove character vector Cells to remove (default=NULL)
    ## cells.to.remain character vector Cells which should remain (default=NULL)
    ## samples.to.remove character vector Samples which should remain (default=NULL)
    ## list of cell counts and cell groups
    extractCodaData=function(ret.groups=TRUE, cell.groups=self$cell.groups, cells.to.remove=NULL, cells.to.remain=NULL,
                             samples.to.remove=NULL) {
      d.counts <- cell.groups %>% data.frame(anno=., group=self$sample.per.cell[names(.)]) %>%
        table() %>% rbind() %>% t()

      if (!is.null(cells.to.remove)) d.counts %<>% .[,!(colnames(.) %in% cells.to.remove)]
      if (!is.null(cells.to.remain)) d.counts %<>% .[,colnames(.) %in% cells.to.remain]
      if (!is.null(samples.to.remove)) d.counts %<>% .[!(rownames(.) %in% samples.to.remove),]

      if (!ret.groups) {
        return(d.counts)
      }

      d.groups <- (self$sample.groups[rownames(d.counts)] == self$target.level) %>%
        setNames(rownames(d.counts))

      return(list(d.counts = d.counts,
                  d.groups = d.groups))
    },

    ## Extract contours from embedding
    ##
    ## groups specify cell groups for contour, multiple cell groups are also supported
    ## color character string (default='black')
    ## linetype numeric (default=2)
    ## conf character string confidence interval of contour (default='10%')
    ## ... additional parameters
    ## density contours
    getDensityContours = function(groups, color='black', linetype=2, conf="10%", ...) {
      cnl <- sn(groups) %>% plapply(function(g) {
        getDensityContour(
          self$embedding, cell.groups=self$cell.groups, linetype=linetype, group=g, conf=conf, color=color, ...
        )
      }, n.cores=self$n.cores, progress=self$verbose, mc.preschedule=TRUE) %>% do.call(c, .)
      return(cnl)
    },

    ## Check cell embedding. The function throws errors in the embedding is empty or incorrectly formatted.
    ##
    ## embedding (default=self$embedding)
    ## the function only checks if there are errors; nothing is returned
    checkCellEmbedding = function(embedding=self$embedding) {
      if (is.null(embedding) || ncol(embedding) != 2) {
        stop("self$embedding must contain 2D cell embedding")
      }

      if (is.null(rownames(embedding))){
        stop("self$embedding must have rownames, equal to cell ids")
      }
    },

    ## Get cluster-free DE input
    ##
    ## genes
    ## min.edge.weight numeric (default=0.6)
    ## list with fields 'cm', 'adj.mat', 'is.ref', 'nns.per.cell'
    getClusterFreeDEInput = function(genes, min.edge.weight=0.0) {
      cm <- self$getJointCountMatrix(raw=FALSE)
      is.ref <- (self$sample.groups[levels(self$sample.per.cell)] == self$ref.level)

      adj.mat <- extractCellGraph(self$data.object) %>% igraph::as_adj()
      diag(adj.mat) <- 1
      cell.names <- intersect(rownames(cm), rownames(adj.mat))

      adj.mat %<>% .[cell.names, cell.names, drop=FALSE] %>% as("TsparseMatrix")

      if (min.edge.weight > 1e-10) {
        samp.per.cell <- self$sample.per.cell[cell.names]
        adj.mat@x[(samp.per.cell[adj.mat@i + 1] != samp.per.cell[adj.mat@j + 1]) & (adj.mat@x < min.edge.weight)] <- 0.0
        adj.mat %<>% drop0() %>% as("TsparseMatrix")
      }

      nns.per.cell <- adj.mat %>% {split(.@j, .@i + 1)} %>%
        .[paste(1:nrow(adj.mat))] %>% setNames(cell.names)
      cm %<>% .[cell.names, genes, drop=FALSE]

      for (id in which(sapply(nns.per.cell, is.null))) {
         # If there are cells without neighbors. May be a consequence of graph modification
        nns.per.cell[[id]] <- numeric()
      }

      return(list(cm=cm, adj.mat=adj.mat, is.ref=is.ref, nns.per.cell=nns.per.cell))
    }
  )
)
