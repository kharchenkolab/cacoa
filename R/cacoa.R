#' @title Cacoa R6 class
#'
#' @import methods
#' @export Cacoa
#' @exportClass Cacoa
Cacoa <- R6::R6Class("Cacoa", lock_objects=F,
  public = list(
    #' @field n.cores number of cores
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

    #' @field sample.per.cell named factor with cell names
    sample.per.cell = NULL,

    #' @field ref.level reference level for sample.group vector
    ref.level = NULL,

    #' @field target.level target/disease level for sample.group vector
    target.level = NULL,

    initialize=function(data.object, sample.groups=NULL, cell.groups=NULL, sample.per.cell=NULL, ref.level=NULL, target.level=NULL, n.cores=1, verbose=TRUE) {
      self$n.cores <- n.cores
      self$verbose <- verbose
      self$ref.level <- ref.level

      if(is.null(target.level)) {
        self$target.level <- "target"
      } else {
        self$target.level <- target.level
      }

      if('Cacoa' %in% class(data.object)) { # copy constructor
        for(n in ls(data.object)) {
          if (!is.function(get(n, data.object))) assign(n, get(n, data.object), self)
        }
      } else {
        # TODO: would be nice to support a list of count matrices as input
        if (!('Conos' %in% class(data.object)) && !('Seurat' %in% class(data.object)))
          stop("only Conos or Seurat v3 data objects are currently supported");

        self$data.object <- data.object
      }

      if(is.null(sample.groups) && !is.null(ref.level)) {
        self$sample.groups <- extractSampleGroups(data.object, ref.level, self$target.level)
      } else {
        self$sample.groups <- sample.groups
      }

      if(is.null(cell.groups)) {
        self$cell.groups <- extractCellGroups(data.object)
      } else {
        self$cell.groups <- cell.groups
      }

      if(is.null(sample.per.cell)) {
        self$sample.per.cell <- extractSamplePerCell(data.object)
      } else {
        self$sample.per.cell <- sample.per.cell
      }
    },

    #' @description  Calculate expression shift magnitudes of different clusters between conditions
    #' @param cell.groups Named cell group factor with cell names (default: stored vector)
    #' @param dist 'JS' - Jensen Shannon divergence, or 'cor' - correlation distance (default="JS")
    #' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
    #' @param within.group.normalization Normalize the shift magnitude by the mean magnitude of within-group variation (default=T)
    #' @param valid.comparisons A logical matrix (rows and columns are samples) specifying valid between-sample comparisons. Note that if within.group.normalization=T, the method will automatically include all within-group comparisons of the samples for which at least one valid pair is included in the valid.comparisons (default=NULL)
    #' @param n.cells Number of cells to subsmaple across all samples (if not specified, defaults to the total size of the smallest cell cluster)
    #' @param n.top.genes Number of top highest-expressed genes to consider (default: all genes)
    #' @param n.subsamples Number of samples to draw (default=100)
    #' @param min.cells Minimum number of cells per cluster/per sample to be included in the analysis (default=10)
    #' @param verbose Print progress
    #' @param n.cores Number of cores (default: stored integer)
    #' @param name Test name (default="expression.shifts")
    #' @return a list include \itemize{
    #'   \item{df: a table with cluster distances (normalized if within.gorup.normalization=T), cell type, number of cells}
    #'   \item{ctdml: raw list of `n.subsamples` sampled distance matrices (cells were subsampled)}
    #'   \item{sample.groups: same as the provided variable}
    #'   \item{valid.comparisons: a matrix of valid comparisons (in this case all should be valid, since we're not restricting samples that should be compared)}
    #' }
    estimateExpressionShiftMagnitudes=function(cell.groups=self$cell.groups, dist='JS', within.group.normalization=TRUE, valid.comparisons=NULL,
                                               n.cells=NULL, n.top.genes=Inf, n.subsamples=100, min.cells=10,
                                               sample.groups=self$sample.groups, n.cores=self$n.cores, verbose=self$verbose,
                                               name="expression.shifts") {
      if (is.null(sample.groups))
        stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if (is.null(cell.groups)) {
        stop("'cell.groups' must be provided either during the object initialization or during this function call")
      }

      count.matrices <- extractRawCountMatrices(self$data.object, transposed=T)

      self$test.results[[name]] <- count.matrices %>%
        estimateExpressionShiftMagnitudes(sample.groups, cell.groups, dist=dist, within.group.normalization=within.group.normalization,
                                          valid.comparisons=valid.comparisons, n.cells=n.cells, n.top.genes=n.top.genes, n.subsamples=n.subsamples,
                                          min.cells=min.cells, n.cores=n.cores, verbose=verbose, transposed.matrices=T)

      return(invisible(self$test.results[[name]]))
    },

    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
    #' @param name Test results to plot (default=expression.shifts)
    #' @param size.norm Plot size normalized results. Requires cell.groups, and sample.per.cell (default=F)
    #' @param notch Show notches in plot, see ggplot2::geom_boxplot for more info (default=T)
    #' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
    #' @param sample.per.cell Named sample factor with cell names (default: stored vector)
    #' @return A ggplot2 object
    plotExpressionShiftMagnitudes=function(name="expression.shifts", size.norm=F, notch = T, cell.groups=self$cell.groups, sample.per.cell=self$sample.per.cell) {
      private$checkTestResults(name)

      if (is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if (is.null(sample.per.cell)) stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      plotExpressionShiftMagnitudes(cluster.shifts = self$test.results[[name]]$df, size.norm = size.norm, notch = notch, cell.groups = cell.groups, sample.per.cell = sample.per.cell)
    },

    #' @description  Calculate expression shift Z scores of different clusters between conditions
    #' @param cell.groups Named cell group factor with cell names (default: stored vector)
    #' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
    #' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
    #' @param sample.per.cell (default: stored vector)
    #' @param n.od.genes (default=1000)
    #' @param n.pcs (default=100)
    #' @param pca.maxit (default=1000)
    #' @param ignore.cache (default=F)
    #' @param name (default="expression.z.scores")
    estimateExpressionShiftZScores=function(cell.groups=self$cell.groups, ref.level=self$ref.level, sample.groups=self$sample.groups,
                                            sample.per.cell=self$sample.per.cell, n.od.genes=1000, n.pcs=100, pca.maxit=1000, ignore.cache=F,
                                            name="expression.z.scores", verbose = self$verbose) {
      if (is.null(cell.groups))
        stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if (is.null(ref.level))
        stop("'ref.level' must be provided either during the object initialization or during this function call")

      if (is.null(sample.groups))
        stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if (is.null(sample.per.cell))
        stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      if(!ref.level %in% sample.groups) stop(paste0("Reference group '",ref.level,"' not in 'sample.groups': ",paste(unique(sample.groups), collapse=" ")))

      if(verbose) cat("Merging count matrices ... ")
      mtx <- extractJointCountMatrix(self$data.object, raw=F)

      if (n.od.genes > 0) {
        if(verbose) cat("done!\nExtracting OD genes ... ")
        mtx %<>% .[, extractOdGenes(self$data.object, n.od.genes)]
      }

      # TODO: use scaled count matrix here
      if (n.pcs < ncol(mtx)) {
        if(verbose) cat("done!\nCalculating PCs ... ")
        pca.name <- paste("joint.pca", n.od.genes, n.pcs, sep="_")
        if (ignore.cache || !is.null(self$cache[[pca.name]])) {
          mtx <- self$cache[[pca.name]]
        } else {
          centers <- Matrix::colMeans(mtx)
          pcs <- mtx %>% irlba::irlba(nv=n.pcs, nu=0, center=centers, right_only=F,
                                      fastpath=T, maxit=pca.maxit, reorth=T)
          mtx <- as.matrix(t(as(t(mtx %*% pcs$v), "dgeMatrix") - t(centers %*% pcs$v)))
          self$cache[[pca.name]] <- mtx
        }
      }

      if(verbose) cat("done!\n")
      self$test.results[[name]] <- estimateExpressionShiftZScores(mtx, sample.per.cell, sample.groups, cell.groups, ref.level, verbose = verbose)
      return(invisible(self$test.results[[name]]))
    },

    #' @description  Plot results from cao$estimateExpressionShiftZScores()
    #' @param type.order (default=NULL)
    #' @param name Test results to plot (default=expression.z.shifts)
    #' @param size.norm Show cluster size-normalized results
    #' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
    #' @param sample.per.cell Named sample factor with cell names (default: stored vector)
    #' @return A ggplot2 object
    plotExpressionShiftZScores=function(type.order=NULL, name="expression.z.scores", size.norm = F, cell.groups = self$cell.groups, sample.per.cell = self$sample.per.cell) {
      private$checkTestResults(name)

      plot.df <- self$test.results[[name]] %>% dplyr::filter(complete.cases(.))
      if (!is.null(type.order)) {
        plot.df %<>% dplyr::filter(Type %in% type.order) %>%
          dplyr::mutate(Type=factor(Type, levels=type.order))
      }

      plotExpressionShiftZScores(plot.df = plot.df, size.norm = size.norm, cell.groups = cell.groups, sample.per.cell = sample.per.cell)
    },

    #' @description Estimate differential gene expression per cell type between conditions
    #' @param cell.groups factor specifying cell types (default=NULL)
    #' @param sample.groups a list of two character vector specifying the app groups to compare (default=NULL)
    #' @param cooks.cutoff cooksCutoff for DESeq2 (default=F)
    #' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
    #' @param common.genes Only investigate common genes across cell groups (default=F)
    #' @param min.cell.count (default=10)
    #' @param independent.filtering independentFiltering for DESeq2 (default=F)
    #' @param n.cores Number of cores (default=1)
    #' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
    #' @param return.matrix Return merged matrix of results (default=T)
    #' @param verbose Show progress (default=T)
    #' @return A list of DE genes
    estimatePerCellTypeDE=function(cell.groups = self$cell.groups, sample.groups = self$sample.groups, ref.level = self$ref.level,
                              common.genes = F, n.cores = self$n.cores, cooks.cutoff = FALSE, min.cell.count = 10, independent.filtering = FALSE,
                              cluster.sep.chr = "<!!>", return.matrix = T, verbose=T) {
      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if(is.null(ref.level)) stop("'ref.level' must be provided either during the object initialization or during this function call")

      if(is.null(sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if(!is.list(sample.groups)) {
        sample.groups <- list(names(sample.groups[sample.groups == ref.level]),
                              names(sample.groups[sample.groups != ref.level])) %>%
          setNames(c(ref.level, self$target.level))
      }

      self$test.results[["de"]] <- extractRawCountMatrices(self$data.object, transposed=T) %>%
        estimatePerCellTypeDE(cell.groups = cell.groups, sample.groups = sample.groups, ref.level = ref.level, n.cores = n.cores,
                            cooks.cutoff = cooks.cutoff, min.cell.count = min.cell.count, independent.filtering = independent.filtering,
                            cluster.sep.chr = cluster.sep.chr, return.matrix = return.matrix)
      return(invisible(self$test.results[["de"]]))
    },

    #' @description Plot number of significant DE genes as a function of number of cells
    #' @param de.raw List with differentially expressed genes per cell group (default: stored list, results from estimatePerCellTypeDE)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="none")
    #' @param label Show labels on plot (default=T)
    #' @param p.adjust.cutoff Adjusted P cutoff (default=0.05)
    #' @return A ggplot2 object
    plotDEGenes=function(de.raw = self$test.results$de, cell.groups = self$cell.groups, legend.position = "none", label = T, p.adjust.cutoff = 0.05) {
      if(is.null(de.raw)) stop("Please run 'estimatePerCellTypeDE' first.")

      # If estimatePerCellTypeDE was run with return.matrix = T, remove matrix before plotting
      if(class(de.raw[[1]]) == "list") de.raw %<>% lapply(`[[`, 1)

      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotDEGenes(de.raw = de.raw, cell.groups = cell.groups, legend.position = legend.position, p.adjust.cutoff = p.adjust.cutoff)
    },

    #' @description Save DE results as JSON files
    #' @param saveprefix Prefix for created files (default=NULL)
    #' @param dir.name Name for directory with results (default="JSON")
    #' @param de.raw List of DE results
    #' @param sample.groups Sample groups (default: stored vector)
    #' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
    #' @param gene.metadata (default=NULL)
    #' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
    #' @param verbose Show progress (default=T)
    saveDEasJSON=function(saveprefix = NULL, dir.name = "JSON", de.raw = self$test.results$de, sample.groups = self$sample.groups, ref.level = self$ref.level, gene.metadata = NULL, cluster.sep.chr = "<!!>", verbose = T) {
      if(is.null(sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if(is.null(ref.level) & !is.list(sample.groups)) stop("'ref.level' must be provided either during the object initialization or during this function call")

      if(!is.list(sample.groups)) {
        sample.groups <- list(names(sample.groups[sample.groups == ref.level]),
                              names(sample.groups[sample.groups != ref.level])) %>%
          setNames(c(ref.level, self$target.level))
      }

      if(is.null(de.raw)) stop("Please run 'estimatePerCellTypeDE' first.")

      if(class(de.raw[[1]]) != "list") stop("Please rerun 'estimatePerCellTypeDE' with return.matrix=T")

      if(is.null(saveprefix)) saveprefix <- ""

      saveDEasJSON(de.raw = de.raw, saveprefix = saveprefix, dir.name = dir.name, gene.metadata = gene.metadata, cluster.sep.chr = cluster.sep.chr, sample.groups = sample.groups, verbose = verbose)
    },

    #' @description  Filter and prepare DE genes for ontology calculations
    #' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
    #' @param n.top.genes Number of most different genes to take as input. If less are left after filtering for p.adj.cutoff, additional genes are included. To disable, set n.top.genes=0 (default=1e2)
    #' @param p.adj.cutoff Cutoff for filtering highly-expressed DE genes (default=0.05)
    #' @param expr.cutoff Cutoff for cells per group expressing a DE gene, i.e., cutoff for highly-expressed genes (default=0.05)
    #' @param de.raw Differentially expressed genes per cell group, results from estimatePerCellTypeDE (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param universe Only set this if a common background gene set is desired for all cell groups (default: NULL)
    #' @param transposed Whether count matrices should be transposed (default=T)
    #' @param verbose Print progress (default=T)
    #' @param n.cores Number of cores to use (default: stored integer)
    #' @return A list containing DE gene IDs, filtered DE genes, and input DE genes
    prepareOntologyData=function(org.db, n.top.genes = 1e2, p.adj.cutoff = 0.05, expr.cutoff = 0.05, de.raw = self$test.results$de, cell.groups = self$cell.groups, universe = NULL, transposed = T, verbose = T, n.cores = self$n.cores) {
      if(is.null(de.raw)) stop("Please run 'estimatePerCellTypeDE' first.")

      # If estimatePerCellTypeDE was run with return.matrix = T, remove matrix before calculating
      if(class(de.raw[[1]]) == "list") de.raw %<>% lapply(`[[`, 1)

      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      self$test.results[["ontology"]] <- extractRawCountMatrices(self$data.object, transposed = transposed) %>%
        prepareOntologyData(org.db = org.db, n.top.genes = n.top.genes, p.adj.cutoff = p.adj.cutoff, expr.cutoff = expr.cutoff, de.raw = de.raw, cell.groups = cell.groups, universe = universe, transposed = transposed, verbose = verbose, n.cores = n.cores)
      return(invisible(self$test.results[["ontology"]]))
    },

    #' @description Plot number of highly-expressed DE genes as a function of number of cells
    #' @param de.filter Filtered DE genes, results from prepareOntologyData (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="none")
    #' @param label Show labels on plot (default=T)
    #' @return A ggplot2 object
    plotFilteredDEGenes=function(de.filter = self$test.results$ontology$de.filter, cell.groups = self$cell.groups, legend.position = "none", label = T) {
      if(is.null(de.filter)) stop("Please run 'estimatePerCellTypeDE' first.")

      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotFilteredDEGenes(de.filter = de.filter, cell.groups = cell.groups, legend.position = legend.position, label = label)
    },

    #' @description  Estimate ontology terms based on DEs
    #' @param type Ontology type, either GO (gene ontology) or DO (disease ontology). Please see DOSE package for more information (default="GO")
    #' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
    #' @param de.gene.ids List containing DE gene IDs, and filtered DE genes (default: stored list, results from prepareOntologyData)
    #' @param universe All measured genes (default: stored vector)
    #' @param go.environment Extracted GO environment. If set to NULL, the environment will be re-extracted (default: stored environment)
    #' @param p.adj Adjusted P cutoff (default=0.05)
    #' @param p.adjust.method Method for calculating adj. P. Please see DOSE package for more information (default="BH")
    #' @param readable Mapping gene ID to gene name (default=T)
    #' @param verbose Print progress (default=T)
    #' @param min.genes Minimum number of input genes overlapping with ontologies (default=0)
    #' @param qvalueCutoff Q value cutoff, please see clusterProfiler package for more information (default=0.2)
    #' @param minGSSize Minimal geneset size, please see clusterProfiler package for more information (default=5)
    #' @param maxGSSize Minimal geneset size, please see clusterProfiler package for more information (default=5e2)
    #' @param ... Additional parameters for sccore:::plapply function
    #' @return A list containing a list of terms per ontology, and a data frame with merged results
    estimateOntology=function(type = "GO", org.db, de.gene.ids = self$test.results$ontology$de.gene.ids, go.environment = self$test.results$GO$go.environment, p.adj = 0.05, p.adjust.method = "BH", readable = T, verbose = T, min.genes = 0, qvalueCutoff = 0.2, minGSSize = 5, maxGSSize = 5e2, ...) {
      if(!is.null(type) & !type %in% c("GO", "DO")) stop("'type' must be 'GO' or 'DO'.")

      if(is.null(de.gene.ids)) stop("Please run 'prepareOntologyData' first.")

      self$test.results[[type]] <- estimateOntology(type = type, org.db = org.db, de.gene.ids = de.gene.ids, go.environment = go.environment, p.adj = p.adj, p.adjust.method = p.adjust.method, readable = readable, verbose = verbose, min.genes = min.genes, qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize, ...)
      return(invisible(self$test.results[[type]]))
    },

    #' @description Bar plot of ontology terms per cell type
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "GO" or "DO" (default="GO")
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @return A ggplot2 object
    plotOntologyDistribution=function(genes = c("up","down"), type = "GO", cell.groups = self$cell.groups) {
      if(is.null(type) || (!type %in% c("GO", "DO"))) stop("'type' must be 'GO' or 'DO'.")

      if(is.null(genes) || (!all(genes %in% c("down","up","all")))) stop("'genes' must be 'down', 'up', 'all', or a combination of these.")

      ont.res <- self$test.results[[type]][["df"]]

      if(is.null(ont.res)) stop(paste0("No results found for '",type,"'. Please run 'estimateOntology' first and specify type='",type,"'."))

      classes <- sapply(ont.res[genes], class)
      if(any(classes == "character")) {
        message(paste0("No results found for genes = '",names(classes[classes == "character"]),"'."))
        genes <- names(classes[classes == "data.frame"])
        if(length(genes) == 0) stop("No results to plot.")
      }

      ont.res %<>% .[names(.) %in% genes]

      if(length(genes) > 1) {
        ont.res %<>% names() %>%
          lapply(function(d) ont.res[[d]] %>% dplyr::mutate(direction = d)) %>%
          Reduce(rbind, .)
      } else {
        ont.res %<>% .[[1]] %>% dplyr::mutate(direction = genes)
      }

      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotOntologyDistribution(type = type, ont.res = ont.res, cell.groups = cell.groups)
    },

    #' @description Plot ontology terms as a function of both number of DE genes, and number of cells.
    #' @param genes Specify which genes to plot, can either be 'down', 'up', 'all' or a combination of these (default=c("up","down"))
    #' @param type Ontology, must be either "GO" or "DO" (default="GO")
    #' @param de.filter Filtered DE genes, results from prepareOntologyData (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param label.x.pos Plot label position on x axis (default=0.01)
    #' @param label.y.pos Plot label position on y axis (default=1)
    #' @param scale Scaling of plots, adjust if e.g. label is misplaced. See cowplot::plot_grid for more info (default=0.93)
    #' @return A ggplot2 object
    plotOntologyTerms=function(genes = c("up","down"), type = "GO", de.filter = self$test.results$ontology$de.filter, cell.groups = self$cell.groups, label.x.pos = 0.01, label.y.pos = 1, scale = 0.93) {
      if(is.null(type) || (!type %in% c("GO", "DO"))) stop("'type' must be 'GO' or 'DO'.")

      if(is.null(genes) || (!all(genes %in% c("down","up","all")))) stop("'genes' must be 'down', 'up', 'all', or a combination of these.")

      if(is.null(de.filter)) stop("Please run 'prepareOntologyData' first.")

      ont.res <- self$test.results[[type]][["df"]]

      if(is.null(ont.res)) stop(paste0("No results found for ",type,". Please run estimateOntology first and specify type='",type,"'."))

      classes <- sapply(ont.res[genes], class)
      if(any(classes == "character")) {
        message(paste0("No results found for genes = '",names(classes[classes == "character"]),"'."))
        genes <- names(classes[classes == "data.frame"])
        if(length(genes) == 0) stop("No results to plot.")
      }

      ont.res %<>% .[names(.) %in% genes]

      if(length(genes) > 1) ont.res %<>% Reduce(rbind, .) else ont.res %<>% .[[1]]

      if (is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotOntologyTerms(type = type, ont.res = ont.res, de.filter = de.filter, cell.groups = cell.groups, label.x.pos = label.x.pos, label.y.pos = label.y.pos, scale = scale)
    },

    #' @description Plot a barplot of ontology terms with adj. P values for a specific cell subgroup
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default="GO")
    #' @param cell.subgroups Cell group to plot (default=NULL)
    #' @param n Number of ontology terms to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
    #' @param p.adj Adjusted P cutoff (default=0.05)
    #' @return A ggplot2 object
    plotOntologyBarplot = function(genes = "up", type = "BP", cell.subgroups = NULL, n = 20, p.adj = 0.05) {
      if(is.null(type) || (!type %in% c("GO","BP","CC","MF","DO"))) stop("'type' must be 'GO', BP', 'CC', 'MF', or 'DO'.")

      if(is.null(genes) || (!genes %in% c("down","up","all"))) stop("'genes' must be 'down', 'up', or 'all'.")

      if(is.null(cell.subgroups)) stop("Please define 'cell.subgroups'.")

      if(type=="DO") {
        ont.res <- self$test.results[["DO"]][["df"]]
      } else {
        ont.res <- self$test.results[["GO"]][["df"]]
      }

      if(is.null(ont.res)) {
        if(type == "DO") t <- "DO" else t <- "GO"
        stop(paste0("No results found for ",type,". Please run estimateOntology first and specify type='",t,"'."))
      }

      ont.res %<>% .[[genes]]

      if(is.null(ont.res)) stop("No results found for genes = '",genes,"'.")

      if(!cell.subgroups %in% unique(ont.res$Group)) stop("'cell.subgroups' not found in results.")

      ont.res %<>% dplyr::filter(Group == cell.subgroups)

      if(type %in% c("BP","CC","MF")) ont.res %<>% dplyr::filter(Type == type)

      plotOntologyBarplot(ont.res = ont.res, genes = genes, type = type, cell.subgroups = cell.subgroups, n = n, p.adj = p.adj)
    },

    #' @description Plot a dotplot of ontology terms with adj. P values for a specific cell subgroup
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default="GO")
    #' @param cell.subgroups Cell group to plot (default=NULL)
    #' @param n Number of ontology terms to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
    #' @param p.adj Adjusted P cutoff (default=0.05)
    #' @return A ggplot2 object
    plotOntologyDotplot = function(genes = "up", type = "BP", cell.subgroups = NULL, n = 20, p.adj = 0.05) {
      if(is.null(type) || (!type %in% c("GO","BP","CC","MF","DO"))) stop("'type' must be 'GO', BP', 'CC', 'MF', or 'DO'.")

      if(is.null(genes) || (!genes %in% c("down","up","all"))) stop("'genes' must be 'down', 'up', or 'all'.")

      if(is.null(cell.subgroups)) stop("Please define 'cell.subgroups'.")

      if(type=="DO") {
        ont.res <- self$test.results[["DO"]][["df"]]
      } else {
        ont.res <- self$test.results[["GO"]][["df"]]
      }

      if(is.null(ont.res)) {
        if(type == "DO") t <- "DO" else t <- "GO"
        stop(paste0("No results found for ",type,". Please run estimateOntology first and specify type='",t,"'."))
      }

      ont.res %<>% .[[genes]]

      if(is.null(ont.res)) stop("No results found for genes = '",genes,"'.")

      if(!cell.subgroups %in% unique(ont.res$Group)) stop("'cell.subgroups' not found in results.")

      ont.res %<>% dplyr::filter(Group == cell.subgroups)

      if(type %in% c("BP","CC","MF")) ont.res %<>% dplyr::filter(Type == type)

      plotOntologyDotplot(ont.res = ont.res, genes = genes, type = type, cell.subgroups = cell.subgroups, n = n, p.adj = p.adj)
    },

    #' @description Plot a heatmap of ontology P values per cell type
    #' @param genes Specify which genes to plot, can either be 'down' for downregulated genes, 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default="GO")
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="left")
    #' @param selection Order of rows in heatmap. Can be 'unique' (only show terms that are unique for any cell type); 'common' (only show terms that are present in at least two cell types); 'all' (all ontology terms) (default="all")
    #' @param n Number of terms to show (default=10)
    #' @param cell.subgroups Cell groups to plot (default=NULL)
    #' @return A ggplot2 object
    plotOntologyHeatmap=function(genes = "up", type = "GO", legend.position = "left", selection = "all", n = 10, cell.subgroups = NULL) {
      if(is.null(type) || (!type %in% c("GO","BP","CC","MF","DO"))) stop("'type' must be 'BP', 'CC', 'MF', or 'DO'.")

      if(is.null(genes) || (!genes %in% c("down","up","all"))) stop("'genes' must be 'down', 'up', or 'all'.")

      if(!is.null(cell.subgroups)) {
        if(length(cell.subgroups) == 1) stop("'cell.subgroups' must contain at least two groups. Please use plotOntologyBarplot or plotOntologyDotplot instead.")
      }

      if(is.null(selection) || (!selection %in% c("unique","common","all"))) stop("'selection' must be one of the following: 'unique', 'common', or 'all'.")

      if(type=="DO") {
        ont.res <- self$test.results[["DO"]][["df"]]
      } else {
        ont.res <- self$test.results[["GO"]][["df"]]
      }

      if(is.null(ont.res)) stop(paste0("No results found for '",type,"'. Please run 'estimateOntology' first and specify type='",type,"'."))

      ont.res %<>% .[[genes]]

      if(!is.null(cell.subgroups) && !cell.subgroups %in% unique(ont.res$Group)) stop("'cell.subgroups' not found in results.")

      plotOntologyHeatmap(type = type, ont.res = ont.res, legend.position = legend.position, selection = selection, n = n, cell.subgroups = cell.subgroups, genes = genes)
    },

    #' @description Plot correlation matrix for ontology terms between cell types
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "GO" or "DO" (default="GO")
    #' @return A ggplot2 object
    plotOntologySimilarities=function(genes = "up", type = "GO") {
      if(is.null(type) || (!type %in% c("GO", "DO"))) stop("'type' must be 'GO' or 'DO'.")

      if(is.null(genes) || (!genes %in% c("down","up","all"))) stop("'genes' must be 'down', 'up', or 'all'.")

      ont.res <- self$test.results[[type]][["df"]]

      if(is.null(ont.res)) stop(paste0("No results found for '",type,"'. Please run 'estimateOntology' first and specify type='",type,"'."))

      ont.res %<>% .[[genes]]

      if(length(ont.res) == 1) stop("Only one group present, correlation cannot be performed.")
      if(length(ont.res) == 0) stop("No significant ontology terms identified. Try relaxing p.adj.")

      plotOntologySimilarities(type = type, ont.res = ont.res, genes = genes)
    },

    #' @description Plot the cell group proportions per sample
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param sample.per.cell Vector indicating sample name with cell names (default: stored vector)
    #' @param sample.groups Vector indicating sample groups with sample names (default: stored vector)
    #' @param cells.to.remove Vector of cell types to remove from the composition
    #' @param cells.to.remain Vector of cell types to remain in the composition
    #' @return A ggplot2 object
    plotProportions=function(legend.position = "right",
                             cell.groups = self$cell.groups,
                             sample.per.cell = self$sample.per.cell,
                             sample.groups = self$sample.groups,
                             cells.to.remove = NULL,
                             cells.to.remain = NULL) {
      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if(is.null(sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if(is.null(sample.per.cell)) stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      if(is.null(cells.to.remove) && is.null(cells.to.remain)){
        plotProportions(legend.position = legend.position, cell.groups = cell.groups, sample.per.cell = sample.per.cell, sample.groups = sample.groups)
      }else{  # Anna modified
        plotProportionsSubset(legend.position = legend.position,
                        cell.groups = cell.groups,
                        sample.per.cell = sample.per.cell,
                        sample.groups = sample.groups,
                        cells.to.remove = cells.to.remove,
                        cells.to.remain = cells.to.remain)
      }

    },

    #' @description Plot the cell numbers per sample
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param sample.per.cell Vector indicating sample name with cell names (default: stored vector)
    #' @param sample.groups Vector indicating sample groups with sample names (default: stored vector)
    #' @return A ggplot2 object
    plotCellNumbers=function(legend.position = "right", cell.groups = self$cell.groups, sample.per.cell = self$sample.per.cell, sample.groups = self$sample.groups) {
      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if(is.null(sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if(is.null(sample.per.cell)) stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      plotCellNumbers(legend.position = legend.position, cell.groups = cell.groups, sample.per.cell = sample.per.cell, sample.groups = sample.groups)
    },

    #' @description Plot compositions in CoDA-PCA space
    #' @return A ggplot2 object
    plotPcaSpace=function(cells.to.remove = NULL) {
      # Cope with levels
      if(is.null(self$ref.level) && is.null(self$target.level)) stop('Target or Reference levels must be provided')
      if((is.null(self$ref.level) || is.null(self$target.level)) && (length(levels(self$sample.groups)) != 2)) stop('Only two levels should be provided')
      if(is.null(self$ref.level)) self$ref.level = setdiff(levels(self$sample.groups), self$target.level)
      if(is.null(self$target.level)) self$target.level = setdiff(levels(self$sample.groups), self$ref.level)

      # Construct sample groups and count data
      # ---- The following can be significantly reduced
      d.counts <- data.frame(anno=self$cell.groups,
                             group=self$sample.per.cell[match(names(self$cell.groups), names(self$sample.per.cell))]) %>%
        table  %>% rbind %>% t
      if(!is.null(cells.to.remove)) d.counts = d.counts[,!(colnames(d.counts) %in% cells.to.remove)]

      d.groups = self$sample.groups[rownames(d.counts)] == self$target.level
      names(d.groups) <- rownames(d.counts)
      # ----

      plotPcaSpace(d.counts, d.groups)
    },

    #' @description Plot compositions in CoDA-CDA space
    #' @return A ggplot2 object
    plotCdaSpace=function(cells.to.remain = NULL,
                          cells.to.remove = NULL,
                          samples.to.remove = NULL) {
      # Cope with levels
      if(is.null(self$ref.level) && is.null(self$target.level)) stop('Target or Reference levels must be provided')
      if((is.null(self$ref.level) || is.null(self$target.level)) && (length(levels(self$sample.groups)) != 2)) stop('Only two levels should be provided')
      if(is.null(self$ref.level)) self$ref.level = setdiff(levels(self$sample.groups), self$target.level)
      if(is.null(self$target.level)) self$target.level = setdiff(levels(self$sample.groups), self$ref.level)

      # Construct sample groups and count data
      # ---- The following can be significantly reduced
      d.counts <- data.frame(anno=self$cell.groups,
                             group=self$sample.per.cell[match(names(self$cell.groups), names(self$sample.per.cell))]) %>%
        table  %>% rbind %>% t

      if(!is.null(cells.to.remain)) d.counts = d.counts[,colnames(d.counts) %in% cells.to.remain]
      if(!is.null(cells.to.remove)) d.counts = d.counts[,!(colnames(d.counts) %in% cells.to.remove)]
      if(!is.null(samples.to.remove)) d.counts = d.counts[!(rownames(d.counts) %in% samples.to.remove),]

      d.groups = self$sample.groups[rownames(d.counts)] == self$target.level
      names(d.groups) <- rownames(d.counts)
      # ----

      plotCdaSpace(d.counts, d.groups)

    },

    #' @description Plot contrast tree
    #' @return A ggplot2 object
    plotContrastTree=function(cells.to.remain = NULL,
                              cells.to.remove = NULL) {
      # Cope with levels
      if(is.null(self$ref.level) && is.null(self$target.level)) stop('Target or Reference levels must be provided')
      if((is.null(self$ref.level) || is.null(self$target.level)) && (length(levels(self$sample.groups)) != 2)) stop('Only two levels should be provided')
      if(is.null(self$ref.level)) self$ref.level = setdiff(levels(self$sample.groups), self$target.level)
      if(is.null(self$target.level)) self$target.level = setdiff(levels(self$sample.groups), self$ref.level)

      # Construct sample groups and count data
      # ---- The following can be significantly reduced
      d.counts <- data.frame(anno=self$cell.groups,
                             group=self$sample.per.cell[match(names(self$cell.groups), names(self$sample.per.cell))]) %>%
        table  %>% rbind %>% t

      if(!is.null(cells.to.remain)) d.counts = d.counts[,colnames(d.counts) %in% cells.to.remain]
      if(!is.null(cells.to.remove)) d.counts = d.counts[,!(colnames(d.counts) %in% cells.to.remove)]

      d.groups = self$sample.groups[rownames(d.counts)] == self$target.level
      names(d.groups) <- rownames(d.counts)
      # ----

      plotContrastTree(d.counts, d.groups)
    },

    #' @description Plot Loadings
    #' @return A ggplot2 object
    estimateCellLoadings=function(n.cell.counts = 1000,
                              n.seed = 239,
                              aplha = 0.01,
                              cells.to.remove = NULL,
                              cells.to.remain = NULL,
                              samples.to.remove = NULL,
                              signif.threshold = 0.2){
      # Cope with levels
      if(is.null(self$ref.level) && is.null(self$target.level)) stop('Target or Reference levels must be provided')
      if((is.null(self$ref.level) || is.null(self$target.level)) && (length(levels(self$sample.groups)) != 2)) stop('Only two levels should be provided')
      if(is.null(self$ref.level)) self$ref.level = setdiff(levels(self$sample.groups), self$target.level)
      if(is.null(self$target.level)) self$target.level = setdiff(levels(self$sample.groups), self$ref.level)

      # Construct sample groups and count data
      # ---- The following can be significantly reduced
      d.counts <- data.frame(anno=self$cell.groups,
                             group=self$sample.per.cell[match(names(self$cell.groups), names(self$sample.per.cell))]) %>%
        table  %>% rbind %>% t
      if(!is.null(cells.to.remove)) d.counts = d.counts[,!(colnames(d.counts) %in% cells.to.remove)]
      if(!is.null(cells.to.remain)) d.counts = d.counts[,colnames(d.counts) %in% cells.to.remain]
      if(!is.null(samples.to.remove)) d.counts = d.counts[!(rownames(d.counts) %in% samples.to.remove),]

      d.groups = self$sample.groups[rownames(d.counts)] == self$target.level
      names(d.groups) <- rownames(d.counts)
      # ----

      self$test.results[['cda']] <- resampleContrast(d.counts, d.groups,
                                                     n.cell.counts = n.cell.counts,
                                                     n.seed = n.seed)


      balances = self$test.results$cda$balances
      n.bal = ncol(balances)
      n.plus = rowSums(balances[,2:n.bal] > 0)
      n.minus = rowSums(balances[,2:n.bal] < 0)
      perm.frac = mapply(function(num1, num2) min(num1, num2), n.plus, n.minus) /
        mapply(function(num1, num2) max(num1, num2), n.plus, n.minus)

      cda.top.cells = names(perm.frac)[perm.frac < signif.threshold]

      self$test.results[['cda.top.cells']] = cda.top.cells

      return(invisible(self$test.results[['cda']]))
    },

    estimateGaPartiotion=function(n.cell.counts = 1000,
                                  n.seed = 239,
                                  aplha = 0.01,
                                  cells.to.remain = NULL,
                                  cells.to.remove = NULL,
                                  samples.to.remove = NULL){
      # Cope with levels
      if(is.null(self$ref.level) && is.null(self$target.level)) stop('Target or Reference levels must be provided')
      if((is.null(self$ref.level) || is.null(self$target.level)) && (length(levels(self$sample.groups)) != 2)) stop('Only two levels should be provided')
      if(is.null(self$ref.level)) self$ref.level = setdiff(levels(self$sample.groups), self$target.level)
      if(is.null(self$target.level)) self$target.level = setdiff(levels(self$sample.groups), self$ref.level)

      # Construct sample groups and count data
      # ---- The following can be significantly reduced
      d.counts <- data.frame(anno=self$cell.groups,
                             group=self$sample.per.cell[match(names(self$cell.groups), names(self$sample.per.cell))]) %>%
        table  %>% rbind %>% t

      if(!is.null(cells.to.remain)) d.counts = d.counts[,colnames(d.counts) %in% cells.to.remain]
      if(!is.null(cells.to.remove)) d.counts = d.counts[,!(colnames(d.counts) %in% cells.to.remove)]
      if(!is.null(samples.to.remove)) d.counts = d.counts[!(rownames(d.counts) %in% samples.to.remove),]

      d.groups = self$sample.groups[rownames(d.counts)] == self$target.level
      names(d.groups) <- rownames(d.counts)
      # ----

      ga.res = gaPartition(d.counts, d.groups)

      self$test.results[['ga.partition']] <- rownames(t(ga.res[1,ga.res[1,] != 0,drop=FALSE]))

      return(invisible(self$test.results[['ga.partition']]))
    },

    #' @description Plot Loadings
    #' @return A ggplot2 object
    plotCellLoadings = function(n.cell.counts = 1000,
                              n.seed = 239,
                              aplha = 0.01,
                              font.size = 16,
                              cells.to.remove = NULL,
                              samples.to.remove = NULL){
      # Cope with levels
      if(is.null(self$ref.level) && is.null(self$target.level)) stop('Target or Reference levels must be provided')
      if((is.null(self$ref.level) || is.null(self$target.level)) && (length(levels(self$sample.groups)) != 2)) stop('Only two levels should be provided')
      if(is.null(self$ref.level)) self$ref.level = setdiff(levels(self$sample.groups), self$target.level)
      if(is.null(self$target.level)) self$target.level = setdiff(levels(self$sample.groups), self$ref.level)

      if(is.null(self$test.results$cda)) self$estimateCellLoadings(n.cell.counts = n.cell.counts,
                                                     n.seed = n.seed,
                                                     aplha = aplha,
                                                     cells.to.remove = cells.to.remove,
                                                     samples.to.remove = samples.to.remove)

      plotCellLoadings(self$test.results[['cda']],
                       aplha = aplha,
                       n.significant.cells = length(self$test.results$cda.top.cells),
                       font.size = font.size)
    },

    extimateWilcoxonTest = function(cell.groups = self$cell.groups,
                                     sample.per.cell = self$sample.per.cell,
                                     sample.groups = self$sample.groups,
                                     cells.to.remove = NULL,
                                     cells.to.remain = NULL){

      p.vals <- calcWilcoxonTest(cell.groups = cell.groups,
                                 sample.per.cell = sample.per.cell,
                                 sample.groups = sample.groups,
                                 cells.to.remove = cells.to.remove,
                                 cells.to.remain = cells.to.remain)

      self$test.results[['p.vals.balances']] <- p.vals
      return(self$test.results[['p.vals.balances']])

      cda = resampleContrast(d.counts, d.groups,
                             n.cell.counts = n.cell.counts,
                             n.seed = n.seed)
      plotCellLoadings(cda$balances, aplha = aplha)
    },

    #' @description Estimate cell density in giving embedding
    #' @param emb cell embedding matrix
    #' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
    #' @param sample.groups @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
    #' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
    #' @param target.level target/disease level for sample.group vector
    #' @param bins number of bins for density esitmation, default 400
    #' @param by.sample  if TRUE, density will esitmated by sample and quantiles normlization will applied to indivisual sample. If FALSE, cell fraction need to be provided and density will simply esitmated by fraction.
    #' @add.ponits add.ponits  show cells in density plot
    estimateCellDensity=function(emb,fraction=NULL,anoSample=NULL,ref.level,target.level,bins=400,by.sample=TRUE){

      if(is.null(emb)) stop("'emb' must be provided either during the object initialization or during this function call")

      if(is.null(self$sample.per.cell)) stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      if(is.null(self$sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if(is.null(self$ref.level)) stop("'ref.level' must be provided either during the object initialization or during this function call")

      if(is.null(self$target.level)) stop("'target.level' must be provided either during the object initialization or during this function call")

      anoSample=self$sample.per.cell
      ref.level=self$ref.level
      target.level=self$target.level
      sample.groups=self$sample.groups

      self$test.results[['bins']]=bins
      res=estimateCellDensity(emb,anoSample=anoSample,sample.groups=sample.groups,bins=bins,ref.level=ref.level,target.level=target.level,fraction=fraction,by.sample=by.sample)
      self$test.results[['DensityMatrix']]=res[['denMatrix.nor']]
      self$test.results[['targetDensity']]=res[['density.fraction']][[target.level]]
      self$test.results[['refDensity']]=res[['density.fraction']][[ref.level]]
      self$test.results[['DiffDensity']]=res[['density.fraction']][[target.level]]-res[['density.fraction']][[ref.level]]

      # adjust embedding to the same density space
      x=emb[,1]
      y=emb[,2]

      x=(x-range(x)[1])
      x=(x/max(x))*bins

      y=(y-range(y)[1])
      y=(y/max(y))*bins
      emb2=data.frame(x=x,y=y)
      self$test.results[['density.emb']]=emb2
      return(invisible(self$density[['bins']]))

    },

    ##' @description Plot cell density
    ##' @param add.ponits default is TRUE, add points to cell density figure
    ##' @param fraction fraction must be provided when add.ponits is TRUE
    plotCellDensity=function(legend=NULL,title=NULL,grid=NULL,add.ponits=TRUE,fraction=NULL){
      bins=self$test.results$bins
      ref=self$ref.level
      target=self$target.level
      p1=plotDensity(cao$test.results$refDensity,bins=bins,col='B',legend=NULL,title=self$ref.level,grid=TRUE)
      p2=plotDensity(cao$test.results$targetDensity,bins=bins,col='B',legend=NULL,title=self$target.level,grid=TRUE)

      if (add.ponits){
        if(is.null(fraction)) stop("'fraction' must be provided when add points")

        emb=self$test.results$density.emb
        emb$Z=1
        nname1=names(fraction[fraction==ref])
        nname1=sample(nname1,min(2000,nrow(emb[nname1,])))

        nname2=names(fraction[fraction==target])
        nname2=sample(nname2,min(2000,nrow(emb[nname2,])))

        p1=p1+geom_point(data=emb[nname1,],aes(x=x,y=y), col='#FCFDBFFF',size=0.00001,alpha=0.2)
        p2=p2+geom_point(data=emb[nname2,],aes(x=x,y=y), col='#FCFDBFFF',size=0.00001,alpha=0.2)
        return(list('ref'=p1,'target'=p2))
      }
    },

    ##' @description esitmate differential cell density
    ##' @param col color palettes, 4 different color palettes are supported; default is yellow-black-magenta; BWR: blue-white-red;  WR: white-read; B: magma in viridi;
    ##' @param fraction A two-level factor on the cell names describing the conditions being compared (default: stored vector)
    ##' @method method to cacuated differential cell density of each bin; substract: target density minus ref density; entropy: estimated kl divergence entropy betwwen sample grapups ; t.test: zscore of t-test,global variacen is setting for t.test;
    DiffCellDensity=function(fraction=NULL,method='substract',legend=NULL,grid=TRUE,col='YBM',title=NULL){
      ref.level=self$ref.level
      target.level=self$target.level
      sample.groups=self$sample.groups
      bins=self$test.results$bins

      DensityMatrix=self$test.results$DensityMatrix
      if (method=='entropy'){
        if(is.null(fraction)) stop("'fraction' must be provided when entropy was used")
      }
      p=DiffCellDensity(DensityMatrix,fraction,sample.groups,bins=bins,target.level=target.level,ref.level=ref.level,method=method,
                        title=title,grid=grid,legend=legend)
      return(p)
    },


    #' @description esitmate expression distance between samples of each cell type
    ##' @param dist what distance measure to use: 'JS' - Jensen-Shannon divergence, 'cor' - Pearson's linear correlation on log transformed values
    ##' @param min.cluster.size minimum number of cells in a cluster (in a sample) for the distance to be estimated. default: 10
    esitmate.expression.dsiatnce=function(dist='JS',min.cells =10){
      cell.type <- self$cell.groups
      samplef=self$sample.per.cell
      count.matrices <- extractRawCountMatrices(self$data.object, transposed=T)
      self$test.results[['exp_disance']]=esitmate.expression.dsiatnce(count.matrices,samplef,cell.type,dist=dist,min.cells =min.cells)
      return(invisible(self$test.results[['exp_disance']]))
    },


    #' @description Plot expression distance
    #' @param disData esitmated cell expression distance with esitmate.expression.dsiatnce
    #' @param type output figure format to measure inter sample distance or intra sample distance;
    #' @param cell.Type default is null for weigeted expression distance across mutiple cell types, if setting, draw plot for specific cell type.
    #' @return A ggplot2 object
    plotExpressionDistance=function(disData=NULL,cell.Type=NULL,type='inter'){

      sample.groups=self$sample.groups
      if (is.null(disData)){
        disData=self$test.results[['exp_disance']]
      }
      plotExpressionDistance(disData,sample.groups,cell.Type=cell.Type,type=type)
    }

  ),
  private = list(
    checkTestResults=function(name) {
      if (is.null(self$test.results[[name]]))
        stop("Test result for ", name, " wasn't found")
    }
  )
)
