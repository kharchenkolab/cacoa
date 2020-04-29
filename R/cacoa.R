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
      self$target.level <- target.level

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
        self$sample.groups <- extractSampleGroups(data.object, ref.level, target.level)
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
    #' @param normalized.distance Plot the absolute median distance (default=F)
    #' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
    #' @param sample.per.cell Named sample factor with cell names (default: stored vector)
    #' @param label Plot labels on size normalized plots (default=T)
    #' @return A ggplot2 object
    plotExpressionShiftMagnitudes=function(name="expression.shifts", size.norm=F, normalized.distance=F, cell.groups=self$cell.groups, sample.per.cell=self$sample.per.cell, label=T) {
      private$checkTestResults(name)
      if (is.null(cell.groups)) {
        stop("'cell.groups' must be provided either during the object initialization or during this function call")
      }

      if (is.null(sample.per.cell)) {
        stop("'sample.per.cell' must be provided either during the object initialization or during this function call")
      }

      if (!size.norm) {
        if(normalized.distance) message("'normalized.distance' has no effect when 'size.norm' = F.")
        if(label) message("'label' has no effect when 'size.norm' = F.")

        m <- max(abs(df$value - 1))

        gg <- ggplot(na.omit(self$test.results[[name]]$df), aes(x=as.factor(Type), y=value)) +
          geom_boxplot(notch=T, outlier.shape=NA) +
          geom_jitter(position=position_jitter(0.1), aes(color=patient), show.legend=FALSE,alpha=0.1) +
          theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.5)) +
          labs(x="", y="Normalized distance") +
          ylim(c(1 - m, 1 + m)) +
          geom_hline(yintercept=1, linetype="dashed", color = "black")
      } else {
        if (length(setdiff(names(cell.groups), names(sample.per.cell)))>0) warning("Cell names in 'cell.groups' and 'sample.per.cell' are not identical, plotting intersect.")

        cct <- table(cell.groups, sample.per.cell[names(cell.groups)])
        cluster.shifts <- cao$test.results[[name]]$df
        x <- tapply(cluster.shifts$value, cluster.shifts$Type, median)

        if(normalized.distance) {
          odf <- data.frame(cell=names(x),size=rowSums(cct)[names(x)],md=abs(1-x))
        } else {
          odf <- data.frame(cell=names(x),size=rowSums(cct)[names(x)],md=x)
        }

        if (label) {
          gg <- ggplot(odf, aes(size,md,color=cell,label=cell)) +
            ggrepel::geom_text_repel()
        } else {
          gg <- ggplot(odf, aes(size,md,color=cell))
        }

        gg <- gg +
          geom_point() +
          guides(color=F) +
          xlab("Cluster size")

        if(normalized.distance) {
          gg <- gg +
            ylab("Absolute median distance")
        } else {
          m <- max(abs(odf$md - 1))

          gg <- gg +
            ylab("Median distance") +
            ylim(c(1 - m,1 + m)) +
            geom_hline(yintercept=1, linetype="dashed", color = "black")
        }
        return(gg)
      }
      return(gg)
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
                                            sample.per.cell=self$sample.per.cell,
                                            n.od.genes=1000, n.pcs=100, pca.maxit=1000, ignore.cache=F,
                                            name="expression.z.scores") {
      if (is.null(cell.groups))
        stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if (is.null(ref.level))
        stop("'ref.level' must be provided either during the object initialization or during this function call")

      if (is.null(sample.groups))
        stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if (is.null(sample.per.cell))
        stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      if(!ref.level %in% sample.groups) stop(paste0("Reference group '",ref.level,"' not in 'sample.groups': ",paste(unique(sample.groups), collapse=" ")))
      mtx <- extractJointCountMatrix(self$data.object, raw=F)

      if (n.od.genes > 0) {
        mtx %<>% .[, extractOdGenes(self$data.object, n.od.genes)]
      }

      # TODO: use scaled count matrix here
      if (n.pcs < ncol(mtx)) {
        pca.name <- paste("joint.pca", n.od.genes, n.pcs, sep="_")
        if (ignore.cache || !is.null(self$cache[[pca.name]])) {
          mtx <- self$cache[[pca.name]]
        } else {
          centers <- Matrix::colMeans(mtx)
          pcs <- mtx %>% irlba::irlba(nv=n.pcs, nu=0, center=Matrix::colMeans(.), right_only=F,
                                      fastpath=T, maxit=pca.maxit, reorth=T)
          mtx <- as.matrix(t(as(t(mtx %*% pcs$v), "dgeMatrix") - t(centers %*% pcs$v)))
          self$cache[[pca.name]] <- mtx
        }
      }

      self$test.results[[name]] <- estimateExpressionShiftZScores(mtx, sample.per.cell, sample.groups, cell.groups, ref.level)
      return(invisible(self$test.results[[name]]))
    },

    #' @description  Plot results from cao$estimateExpressionShiftZScores()
    #' @param type.order (default=NULL)
    #' @param name Test results to plot (default=expression.z.shifts)
    #' @return A ggplot2 object
    plotExpressionShiftZScores=function(type.order=NULL, name="expression.z.scores") {
      private$checkTestResults(name)

      plot.df <- self$test.results[[name]]
      if (!is.null(type.order)) {
        plot.df %<>% dplyr::filter(Type %in% type.order) %>%
          dplyr::mutate(Type=factor(Type, levels=type.order))
      }

      ggplot(plot.df, aes(x=Type, y=distance)) +
        geom_boxplot(outlier.alpha=0, show.legend=F) +
        geom_hline(aes(yintercept=split(distance, Type) %>% sapply(median) %>% median()), color="darkred", size=1) +
        geom_hline(aes(yintercept=0), color="darkgreen", size=1) +
        labs(x="", y="normalized distance") +
        scale_y_continuous(expand=c(0, 0)) +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=9))
    },

    #' @description Estimate differential gene expression per cell type between conditions
    #' @param cell.groups factor specifying cell types (default=NULL)
    #' @param sample.groups a list of two character vector specifying the app groups to compare (default=NULL)
    #' @param cooks.cutoff cooksCutoff for DESeq2 (default=F)
    #' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
    #' @param min.cell.count (default=10)
    #' @param independent.filtering independentFiltering for DESeq2 (default=F)
    #' @param n.cores Number of cores (default=1)
    #' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
    #' @param return.details Return details (default=T)
    #' @param verbose Show progress (default=T)
    #' @return A list of DE genes
    getPerCellTypeDE=function(cell.groups = self$cell.groups, sample.groups = self$sample.groups, ref.level = self$ref.level,
                              n.cores = self$n.cores, cooks.cutoff = FALSE, min.cell.count = 10, independent.filtering = FALSE,
                              cluster.sep.chr = "<!!>", return.details = TRUE, verbose=T) {
      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if(is.null(ref.level)) stop("'ref.level' must be provided either during the object initialization or during this function call")

      if(is.null(sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if(!is.list(sample.groups)) {
        sample.groups <- list(names(sample.groups[sample.groups == ref.level]),
                              names(sample.groups[sample.groups != ref.level])) %>%
          setNames(c(ref.level, "target"))
      }

      self$test.results[["de"]] <- extractRawCountMatrices(self$data.object, transposed=T) %>%
        getPerCellTypeDEmat(cell.groups = cell.groups, sample.groups = sample.groups, ref.level = ref.level, n.cores = n.cores,
                            cooks.cutoff = cooks.cutoff, min.cell.count = min.cell.count, independent.filtering = independent.filtering,
                            cluster.sep.chr = cluster.sep.chr, return.details = return.details)
      return(invisible(self$test.results[["de"]]))
    },

    #' @description Plot number of significant DE genes as a function of number of cells
    #' @param de.raw List with differentially expressed genes per cell group (default: stored list, results from getPerCellTypeDE)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="bottom")
    #' @param p.adjust.cutoff Adjusted P cutoff (default=0.05)
    #' @return A ggplot2 object
    plotDEGenes=function(de.raw=self$test.results$de, cell.groups=self$cell.groups, legend.position="bottom", p.adjust.cutoff=0.05) {
      if(is.null(de.raw)) stop("Please run 'getPerCellTypeDE' first.")

      if (is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotDEGenes(de.raw=de.raw, cell.groups=cell.groups, legend.position=legend.position, p.adjust.cutoff=p.adjust.cutoff)
    },

    #' @description  Filter and prepare DE genes for onthology calculations
    #' @param de.raw Differentially expressed genes per cell group, results from getPerCellTypeDE (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param org Organism, can be "human", "mouse", "zebrafish", "worm", or "fly" (default="human")
    #' @param stat.cutoff Cutoff for filtering highly-expressed DE genes (default=3)
    #' @param verbose Print progress (default=T)
    #' @return A list containing DE gene IDs, filtered DE genes, and input DE genes
    prepareOnthologyData=function(de.raw=self$test.results$de, cell.groups=self$cell.groups, org="human", stat.cutoff=3, verbose=T) {
      if (is.null(de)) stop("Please run 'getPerCellTypeDE' first.")

      if (is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      self$test.results[["onthology"]] <- extractRawCountMatrices(self$data.object, transposed=T) %>%
        prepareOnthologyData(de.raw=de.raw, transpose=T, cell.groups=cell.groups, org=org, verbose=verbose, stat.cutoff=stat.cutoff)
      return(invisible(self$test.results[["onthology"]]))
    },

    #' @description Plot number of highly-expressed DE genes as a function of number of cells
    #' @param de.filter Filtered DE genes, results from prepareOnthologyData (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="bottom")
    #' @return A ggplot2 object
    plotFilteredDEGenes=function(de.filter=self$test.results$onthology$de.filter, cell.groups=self$cell.groups, legend.position="bottom") {
      if(is.null(de.filter)) stop("Please run 'getPerCellTypeDE' first.")

      if (is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotFilteredDEGenes(de.filter=de.filter, cell.groups=cell.groups, legend.position=legend.position)
    },

    #' @description  Calculate onthologies based on DEs
    #' @param type Onthology type, either GO (gene onthology) or DO (disease onthology). Please see DOSE package for more information.
    #' @param de.gene.ids List containing DE gene IDs, and filtered DE genes (default: stored list, results from prepareOnthologyData)
    #' @param org Organism, can be "human", "mouse", "zebrafish", "worm", or "fly" (default="human")
    #' @param p.adj Adjusted P cutoff (default=0.05)
    #' @param p.adjust.method Method for calculating adj. P. Please see DOSE package for more information (default="BH")
    #' @param readable Mapping gene ID to gene name (default=T)
    #' @param n.cores Number of cores used (default: stored vector)
    #' @param verbose Print progress (default=T)
    #' @param ... Additional parameters for sccore:::plapply function
    #' @return A list containing a list of onthologies per type of onthology, and a data frame with merged results
    estimateOnthology=function(type=NULL, de.gene.ids=self$test.resuls$onthology$de.gene.ids, background=self$test.results$onthology$background, org="human", p.adj=0.05, p.adjust.method="BH", readable=T, verbose=T, ...) {
      if(is.null(type)) stop("'type' must be 'GO' or 'DO'.")

      if(is.null(de.gene.ids)) stop("Please run 'prepareOnthologyData' first.")

      self$test.results[[type]] <- estimateOnthology(type=type, de.gene.ids=de.gene.ids, background = background, org=org, p.adj=p.adj, p.adjust.method=p.adjust.method, readable=readable, verbose=verbose, ...)
      return(invisible(self$test.results[[type]]))
    },

    #' @description Bar plot of onthologies per cell type
    #' @param type Onthology, must be either "GO" or "DO" (default=NULL)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @return A ggplot2 object
    plotOnthologyDistribution=function(type=NULL, cell.groups=self$cell.groups) {
      if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

      ont.res <- self$test.results[[type]][["df"]]

      if(is.null(ont.res)) stop(paste0("No results found for '",type,"'. Please run 'estimateOnthology' first and specify type='",type,"'."))

      if((ont.res$Type %>% unique %>% length) == 1) stop("The input only contains one cell type.")

      plotOnthologyDistribution(type=type, ont.res=ont.res, cell.groups=cell.groups)
    },

    #' @description Plot onthology terms as a function of both number of DE genes, and number of cells.
    #' @param type Onthology, must be either "GO" or "DO" (default=NULL)
    #' @param de.filter Filtered DE genes, results from prepareOnthologyData (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param show.legend Include legend in plot (default=T)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="bottom")
    #' @param label.x.pos Plot label position on x axis (default=0.01)
    #' @param label.y.pos Plot label position on y axis (default=1)
    #' @param rel_heights Relative heights for plots. Only relevant if show.legend=T. See cowplot::plot_grid for more info (default=c(2.5, 0.5))
    #' @param scale Scaling of plots, adjust if e.g. label is misplaced. See cowplot::plot_grid for more info (default=0.93)
    #' @return A ggplot2 object
    plotOnthologyTerms=function(type=NULL, de.filter = self$test.results$onthology$de.filter, cell.groups=self$cell.groups, show.legend=T, legend.position="bottom", label.x.pos=0.01, label.y.pos=1, rel_heights = c(2.5, 0.5), scale = 0.93) {
      if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

      if(is.null(de.filter)) stop("Please run 'prepareOnthologyData' first.")

      ont.res <- self$test.results[[type]][["df"]]

      if(is.null(ont.res)) stop(paste0("No results found for ",type,". Please run estimateOnthology first and specify type='",type,"'."))

      if (is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotOnthologyTerms(type=type, ont.res=ont.res, de.filter=de.filter, cell.groups=cell.groups, show.legend=show.legend, legend.position=legend.position, label.x.pos=label.x.pos, label.y.pos=label.y.pos, rel_heights=rel_heights, scale=scale)
    },

    #' @description Plot a heatmap of onthology P values per cell type
    #' @param type Onthology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default=NULL)
    #' @param ont.data List containing a list of results from estimateOnthology
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="left")
    #' @param order Order of rows in heatmap. Can be 'unique' (only show onthologies that are unique for any cell type); 'unique-max-row' (same as 'unique' but ordered by P value); 'all-max-rowsum' (all onthologies ordered by cumulative P value for all cell types); 'all-max-row' (all onthologies ordered by max P value) (default="all-max-row")
    #' @param n Number of onthologies to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
    #' @return A ggplot2 object
    plotOnthologyHeatmap=function(type=NULL, legend.position = "left", order = "all-max-row", n = 10) {
      if(type=="BP" | type=="CC" | type=="MF") {
        ont.res <- self$test.results[["GO"]][["df"]]
      } else if(type=="DO") {
        ont.res <- self$test.results[["DO"]][["df"]]
      } else {
        stop("'type' must be 'BP', 'CC', 'MF', or 'DO'.")
      }

      if(is.null(ont.res)) stop(paste0("No results found for '",type,"'. Please run 'estimateOnthology' first and specify type='",type,"'."))

      plotOnthologyHeatmap(type=type, ont.res=ont.res, legend.position=legend.position, order=order, n=n)
    },

    #' @description Plot correlation matrix for onthologies between cell types
    #' @param type Onthology, must be either "GO" or "DO" (default=NULL)
    #' @param ont.data List containing a list of results from estimateOnthology
    #' @return A ggplot2 object
    plotOnthologyCorrelations=function(type=NULL) {
      if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

      ont.res <- self$test.results[[type]][["list"]]

      if(is.null(ont.res)) stop(paste0("No results found for '",type,"'. Please run 'estimateOnthology' first and specify type='",type,"'."))

      plotOnthologyCorrelations(type=type, ont.res=ont.res)
    }
  ),
  private = list(
    checkTestResults=function(name) {
      if (is.null(self$test.results[[name]]))
        stop("Test result for ", name, " wasn't found")
    }
  )
)
