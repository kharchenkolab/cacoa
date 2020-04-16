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

    #' @field reference level for sample.group vector
    ref.level = NULL,

    #' @field target/disease level for sample.group vector
    target.level = NULL,

    initialize=function(data.object, sample.groups=NULL, cell.groups=NULL, sample.per.cell=NULL, ref.level=NULL, target.level=NULL, n.cores=parallel::detectCores(logical=F), verbose=TRUE) {
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

      if(is.null(sample.groups) && (!is.null(ref.level) && !is.null(target.level))) {
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
    #'
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

      if (is.null(sample.per.cel))
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

    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
    #' @param name Test results to plot (default=expression.shifts)
    #' @param size.norm Plot size normalized results. Requires cell.groups, and sample.per.cell (default=F)
    #' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
    #' @param sample.per.cell Named sample factor with cell names (default: stored vector)
    #' @param label Plot labels on size normalized plots (default=T)
    #' @return A ggplot2 object
    plotExpressionShiftMagnitudes=function(name="expression.shifts", size.norm=F, cell.groups=self$cell.groups, sample.per.cell=self$sample.per.cell, label=T) {
      private$checkTestResults(name)
      if (is.null(cell.groups)) {
        stop("'cell.groups' must be provided either during the object initialization or during this function call")
      }

      if (is.null(sample.per.cell)) {
        stop("'sample.per.cell' must be provided either during the object initialization or during this function call")
      }

      if (!size.norm) {
        gg <- ggplot(na.omit(self$test.results[[name]]$df),aes(x=as.factor(Type), y=value)) +
          geom_boxplot(notch=T, outlier.shape=NA) +
          geom_jitter(position=position_jitter(0.1), aes(color=patient), show.legend=FALSE,alpha=0.1) +
          theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.5)) +
          labs(x="", y="Normalized distance") +
          geom_hline(yintercept=1, linetype="dashed", color = "black")
      } else {
        if (length(setdiff(names(cell.groups), names(sample.per.cell)))>0) warning("Cell names in 'cell.groups' and 'sample.per.cell' are not identical, plotting intersect.")

        cct <- table(cell.groups, sample.per.cell[names(cell.groups)])
        cluster.shifts <- cao$test.results[[name]]$df
        x <- tapply(cluster.shifts$value, cluster.shifts$Type, median)
        odf <- data.frame(cell=names(x),size=rowSums(cct)[names(x)],md=x)

        if (label) {
          gg <- ggplot(odf, aes(size,md,color=cell,label=cell)) +
            ggrepel::geom_text_repel()
        } else {
          gg <- ggplot(odf, aes(size,md,color=cell))
        }

        gg <- gg +
          geom_point() +
          guides(color=F) + geom_hline(yintercept=1, linetype="dashed", color = "black") +
          ylab("Median distance")

        return(gg)
      }
      return(gg)
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

    #' @title Prepare pathway data
    #' @description  Filter and prepare DE genes for onthology calculations
    #' @param de List with differentially expressed genes per cell group
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param OrgDB Genome-wide annotation (default=org.Hs.eg.db)
    #' @param stat.cutoff Cutoff for filtering highly-expressed DE genes (default=3)
    #' @param verbose Print progress (default=T)
    #' @return A list containing DE gene IDs, filtered DE genes, and input DE genes
    preparePathwayData <- function(de, cell.groups=self$cell.groups, OrgDB=org.Hs.eg.db, stat.cutoff=3, verbose=T) {
      if (is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if(verbose) cat("Extracting raw count matrices ... ")
      count.matrices <- extractRawCountMatrices(self$data.object, transposed=T)
      if(verbose) cat("done!\n")

      self$pathway.data <- preparePathwayData(cms=count.matrices, de=de, transpose=T, cell.groups=cell.groups, OrgDB=OrgDB, verbose=verbose, stat.cutoff=stat.cutoff)
      return(invisible(self$pathway.data))
    },

    #' @title Estimate onthology
    #' @description  Calculate onthologies based on DEs
    #' @param type Onthology type, either GO (gene onthology) or DO (disease onthology). Please see DOSE package for more information.
    #' @param pathway.data List containing DE gene IDs, and filtered and unfiltered DE genes
    #' @param OrgDB Genome-wide annotation (default=org.Hs.eg.db)
    #' @param p.adj Adjusted P cutoff (default=0.05)
    #' @param p.adjust.method Method for calculating adj. P. Please see DOSE package for more information (default="BH")
    #' @param readable Mapping gene ID to gene name (default=T)
    #' @param n.cores Number of cores used (default: stored vector)
    #' @param verbose Print progress (default=T)
    #' @param ... Additional parameters for sccore:::plapply function
    #' @return A list containing a list of onthologies per type of onthology, and a data frame with merged results
    estimateOnthology <- function(type=NULL, pathway.data=self$pathway.data, OrgDB=org.Hs.eg.db, p.adj=0.05, p.adjust.method="BH", readable=T, n.cores=self$n.cores, verbose=T, ...) {
      if(is.null(type)) stop("'type' must be 'GO' or 'DO'.")

      if(is.null(pathway.data)) stop("Please run 'preparePathwayData' first.")

      self$pathway.data[[type]] <- estimateOnthology(type=type, pathway.data=pathway.data, OrgDB=OrgDB, p.adj=p.adj, p.adjust.method=p.adjust.method, readable=readable, n.cores=n.cores, verbose=verbose, ...)
      return(invisible(self$pathway.data[[type]]))
    },

    #' @title Plot onthology terms
    #' @description Plot onthology terms as a function of both number of DE genes, and number of cells.
    #' @param type Onthology, must be either "GO" or "DO" (default=NULL)
    #' @param pathway.data Results from preparePathwayData (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param show.legend Include legend in plot (default=T)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="bottom")
    #' @param label.x.pos Plot label position on x axis (default=0.01)
    #' @param label.y.pos Plot label position on y axis (default=1)
    #' @param rel_heights Relative heights for plots. Only relevant if show.legend=T. See cowplot::plot_grid for more info (default=c(2.5, 0.5))
    #' @param scale Scaling of plots, adjust if e.g. label is misplaced. See cowplot::plot_grid for more info (default=0.93)
    #' @return A ggplot2 object
    plotOnthologyTerms <- function(type=NULL, pathway.data=self$pathway.data, cell.groups=self$cell.groups, show.legend=T, legend.position="bottom", label.x.pos=0.01, label.y.pos=1, rel_heights = c(2.5, 0.5), scale = 0.93) {
      if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

      if(is.null(pathway.data)) stop("Please run 'preparePathwayData' first.")

      ont.res <- pathway.data[[type]][["df"]]

      if(is.null(ont.res)) stop(paste0("No results found for ",type,". Please run estimateOnthology first and specify type='",type,"'."))

      if (is.null(cell.groups))
        stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotOnthologyTerms(ont.res=ont.res, type=type, de.genes.filtered=pathway.data$de.genes.filtered, cell.groups=cell.groups, show.legend=show.legend, legend.position=legend.position, label.x.pos=label.x.pos, label.y.pos=label.y.pos, rel_heights=rel_heights, scale=scale)
    },

    #' @title Plot DE genes
    #' @description Plot number of DE genes as a function of number of cells
    #' @param pathway.data Results from preparePathwayData (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param show.legend Include legend in plot (default=T)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="bottom")
    #' @param p.adj Adjusted P cutoff (default=0.05)
    #' @param label.x.pos Plot label position on x axis (default=0.01)
    #' @param label.y.pos Plot label position on y axis (default=1)
    #' @param rel_heights Relative heights for plots. Only relevant if show.legend=T. See cowplot::plot_grid for more info (default=c(2.5, 0.5))
    #' @param scale Scaling of plots, adjust if e.g. label is misplaced. See cowplot::plot_grid for more info (defaul=0.93)
    #' @return A ggplot2 object
    plotDEGenes <- function(pathway.data=self$pathway.data, cell.groups=self$cell.groups, show.legend=T, legend.position="bottom", p.adj=0.05, label.x.pos=0.01, label.y.pos=1, rel_heights=c(2.5, 0.5), scale=0.93) {
      if(is.null(pathway.data)) stop("Please run 'preparePathwayData' first.")

      if (is.null(cell.groups))
        stop("'cell.groups' must be provided either during the object initialization or during this function call")

      plotDEgenes(de.raw=pathway.data$de.raw, de.genes.filtered=pathway.data$de.genes.filtered, cell.groups=cell.groups, show.legend=show.legend, legend.position=legend.position, p.adj=p.adj, label.x.pos=label.x.pos, label.y.pos=label.y.pos, rel_heights=rel_heights, scale=scale)
    },

    #' @title Plot pathway distribution
    #' @description Bar plot of onthology pathways per cell type
    #' @param type Onthology, must be either "GO" or "DO" (default=NULL)
    #' @param pathway.data List containing a list of results from estimateOnthology
    #' @return A ggplot2 object
    plotPathwayDistribution <- function(type=NULL, pathway.data=self$pathway.data) {
      if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

      ont.res <- pathway.data[[type]]

      if(is.null(ont.res)) stop(paste0("No results found for ",type,". Please run estimateOnthology first and specify type='",type,"'."))

      plotPathwayDistribution(ont.res=ont.res, type=type)
    },

    #' @title Plot onthology heatmap
    #' @description Plot a heatmap of onthology P values per cell type
    #' @param type Onthology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default=NULL)
    #' @param pathway.data List containing a list of results from estimateOnthology
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="left")
    #' @param order Order of rows in heatmap. Can be 'unique' (only show pathways that are unique for any cell type); 'unique-max-row' (same as 'unique' but ordered by P value); 'all-max-rowsum' (all pathways ordered by cumulative P value for all cell types); 'all-max-row' (all pathways ordered by max P value) (default="all-max-row")
    #' @param n Number of pathways to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
    #' @return A ggplot2 object
    plotOnthologyHeatmap <- function(type=NULL, pathway.data=self$pathway.data, legend.position = "left", order = "all-max-row", n = 10) {
      if(type=="BP" | type=="CC" | type=="MF") {
        ont.res <- pathway.data[["GO"]]
      } else if(type=="DO") {
        ont.res <- pathway.data[["DO"]]
      } else {
        stop("'type' must be 'BP', 'CC', 'MF', or 'DO'.")
      }

      plotOnthologyHeatmap(type=type, ont.res=ont.res, legend.position=legend.position, order=order, n=n)
    },

    #' @title Plot onthology correlations
    #' @description Plot correlation matrix for onthologies between cell types
    #' @param type Onthology, must be either "GO" or "DO" (default=NULL)
    #' @param pathway.data List containing a list of results from estimateOnthology
    #' @return A ggplot2 object
    plotOnthologyCorrelations <- function(type=NULL, pathway.data=self$pathway.data) {
      if(is.null(type) & type!="GO" & type!="DO") stop("'type' must be 'GO' or 'DO'.")

      ont.res <- pathway.data[[type]][["list"]]

      if(is.null(ont.res)) stop(paste0("No results found for ",type,". Please run estimateOnthology first and specify type='",type,"'."))

      plotOnthologyCorrelations(ont.res=ont.res, type=type)
    }
  ),
  private = list(
    checkTestResults=function(name) {
      if (is.null(self$test.results[[name]]))
        stop("Test result for ", name, " wasn't found")
    }
  )
)
