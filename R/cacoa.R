#' @title Cacoa R6 class
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

    #' @field sample.groups 2-factor vector with annotation of groups per sample
    sample.groups = NULL,

    initialize=function(data.object, sample.groups=NULL, n.cores=parallel::detectCores(logical=F), verbose=TRUE) {
      self$n.cores <- n.cores
      self$verbose <- verbose
      self$sample.groups <- sample.groups

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
    },

    #' @description  Calculate expression shift magnitudes of different clusters between conditions
    #' @param sample.groups a two-level factor on the sample names describing the conditions being compared
    #' @param groups cell cluster factor
    #' @param dist 'JS' - Jensen Shannon divergence, or 'cor' - correlation distance
    #' @param within.group.normalization normalize the shift magnitude by the mean magnitude of within-group variation
    #' @param valid.comparisons a logical matrix (rows and columns are samples) specifying valid between-sample comparisons. Note that if within.group.normalization=T, the method will automatically include all within-group comparisons of the samples for which at least one valid pair is included in the valid.comparisons
    #' @param n.cells number of cells to subsmaple across all samples (if not specified, defaults to the total size of the smallest cell cluster)
    #' @param n.top.genes number of top highest-expressed genes to consider (default: all genes)
    #' @param n.subsamples number of samples to draw (default:100)
    #' @param min.cells minimum number of cells per cluster/per sample to be included in the analysis
    #' @param ref.level reference sample group, e.g., ctrl, healthy, or untreated.
    #' @param n.cores number of cores to use
    #' @param verbose
    #'
    #' @return a list include \itemize{
    #'   \item{df: a table with cluster distances (normalized if within.gorup.normalization=T), cell type, number of cells}
    #'   \item{ctdml: raw list of `n.subsamples` sampled distance matrices (cells were subsampled)}
    #'   \item{sample.groups: same as the provided variable}
    #'   \item{valid.comparisons: a matrix of valid comparisons (in this case all should be valid, since we're not restricting samples that should be compared)}
    #' }
    estimateExpressionShiftMagnitudes=function(groups=NULL, dist='JS', within.group.normalization=TRUE, valid.comparisons=NULL,
                                               n.cells=NULL, n.top.genes=Inf, n.subsamples=100, min.cells=10,
                                               sample.groups=self$sample.groups, n.cores=self$n.cores, verbose=self$verbose,
                                               name="expression.shifts") {
      if (is.null(sample.groups))
        stop("sample.groups must be provided either during the object initialization or during this function call")

      if (is.null(groups)) {
        groups <- extractGroups(self$data.object)
      }

      count.matrices <- extractRawCountMatrices(self$data.object, transposed=T)

      self$test.results[[name]] <- count.matrices %>%
        estimateExpressionShiftMagnitudes(sample.groups, groups, dist=dist, within.group.normalization=within.group.normalization,
                                          valid.comparisons=valid.comparisons, n.cells=n.cells, n.top.genes=n.top.genes, n.subsamples=n.subsamples,
                                          min.cells=min.cells, n.cores=n.cores, verbose=verbose, transposed.matrices=T)

      return(invisible(self$test.results[[name]]))
    },

    estimateExpressionShiftZScores=function(groups, ref.level, sample.groups=self$sample.groups,
                                            n.od.genes=1000, n.pcs=100, pca.maxit=1000, ignore.cache=F,
                                            name="expression.z.scores") {
      if(!ref.level %in% sample.groups) stop(paste0("Reference group '",ref.level,"' not in sample groups: ",paste(unique(sample.groups), collapse=" ")))
      sample.per.cell <- extractSamplePerCell(self$data.object)
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

      self$test.results[[name]] <- estimateExpressionShiftZScores(mtx, sample.per.cell, sample.groups, groups, ref.level)
      return(invisible(self$test.results[[name]]))
    },

    plotExpressionShiftMagnitudes=function(name="expression.shifts") {
      private$checkTestResults(name)

      ggplot(na.omit(self$test.results[[name]]$df),aes(x=as.factor(Type), y=value)) +
        geom_boxplot(notch=T, outlier.shape=NA) +
        geom_jitter(position=position_jitter(0.1), aes(color=patient), show.legend=FALSE,alpha=0.1) +
        theme(axis.text.x=element_text(angle = 90, hjust=1), axis.text.y=element_text(angle=90, hjust=0.5)) +
        labs(x="", y="normalized distance") +
        geom_hline(yintercept=1, linetype="dashed", color = "black")
    },

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
    }
  ),
  private = list(
    checkTestResults=function(name) {
      if (is.null(self$test.results[[name]]))
        stop("Test result for ", name, " wasn't found")
    }
  )
)
