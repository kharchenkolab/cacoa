#' @title Cacoa R6 class
#'
#' @import methods enrichplot dplyr
#' @export Cacoa
#' @exportClass Cacoa
#' @param sample.groups a two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param n.cores number of cores for parallelisation
#' @param verbose show progress (default: stored value)
#' @param name field name where the test results are stored
#' @param n.top.genes number of top genes for estimation
#' @param gene.selection a method to select top genes, "change" selects genes by cluster-free Z-score change,
#' "expression" picks the most expressed genes and "od" picks overdispersed genes.  Default: "change".
#' @param excluded.genes list of genes to exclude during estimation. For example, a list of mitochondrial genes.
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

    initialize=function(data.object, sample.groups=NULL, cell.groups=NULL, sample.per.cell=NULL, ref.level=NULL, target.level=NULL, sample.groups.palette=NULL, cell.groups.palette=NULL,
                        embedding=extractEmbedding(data.object), n.cores=1, verbose=TRUE) {
      if ('Cacoa' %in% class(data.object)) { # copy constructor
        for (n in ls(data.object)) {
          if (!is.function(get(n, data.object))) assign(n, get(n, data.object), self)
        }

        return()
      }

      self$n.cores <- n.cores
      self$verbose <- verbose
      self$ref.level <- ref.level

      if (is.null(ref.level))
        stop("ref.level must be provided")

      if(is.null(target.level)) {
        self$target.level <- "target"
      } else {
        self$target.level <- target.level
      }

      # TODO: would be nice to support a list of count matrices as input
      if (!('Conos' %in% class(data.object)))
        stop("only Conos data objects are currently supported");

      self$data.object <- data.object

      if(is.null(sample.groups)) {
        self$sample.groups <- extractSampleGroups(data.object, ref.level, self$target.level)
      } else {
        self$sample.groups <- sample.groups <- as.factor(sample.groups)
      }

      if(is.null(cell.groups)) {
        self$cell.groups <- extractCellGroups(data.object)
      } else {
        self$cell.groups <- as.factor(cell.groups)
      }

      if(is.null(sample.per.cell)) {
        self$sample.per.cell <- extractSamplePerCell(data.object)
      } else {
        self$sample.per.cell <- sample.per.cell
      }

      if(is.null(sample.groups.palette)) {
        self$sample.groups.palette <- setNames(rev(scales::hue_pal()(length(levels(sample.groups)))), levels(sample.groups))
      } else {
        self$sample.groups.palette <- sample.groups.palette
      }

      if(is.null(cell.groups.palette)) {
        self$cell.groups.palette <- setNames(rainbow(length(levels(cell.groups)),s=0.9,v=0.9), levels(cell.groups))
      } else {
        self$cell.groups.palette <- cell.groups.palette
      }

      self$embedding <- embedding;
    },

    ### Expression shifts

    #' @description  Calculate expression shift magnitudes of different clusters between conditions
    #' @param cell.groups Named cell group factor with cell names (default: stored vector)
    #' @param dist 'JS' - Jensen Shannon divergence, or 'cor' - correlation distance (default="JS")
    #' @param within.group.normalization Normalize the shift magnitude by the mean magnitude of within-group variation (default=T)
    #' @param valid.comparisons A logical matrix (rows and columns are samples) specifying valid between-sample comparisons. Note that if within.group.normalization=T, the method will automatically include all within-group comparisons of the samples for which at least one valid pair is included in the valid.comparisons (default=NULL)
    #' @param n.cells Number of cells to subsmaple across all samples (if not specified, defaults to the total size of the smallest cell cluster)
    #' @param n.top.genes Number of top highest-expressed genes to consider (default: all genes)
    #' @param n.subsamples Number of samples to draw (default=100)
    #' @param min.cells Minimum number of cells per cluster/per sample to be included in the analysis (default=10)
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
                                               name="expression.shifts", ...) {
      count.matrices <- extractRawCountMatrices(self$data.object, transposed=T)


      self$test.results[[name]] <- count.matrices %>%
        estimateExpressionShiftMagnitudes(sample.groups, cell.groups, dist=dist, within.group.normalization=within.group.normalization,
                                          valid.comparisons=valid.comparisons, n.cells=n.cells, n.top.genes=n.top.genes, n.subsamples=n.subsamples,
                                          min.cells=min.cells, n.cores=n.cores, verbose=verbose, transposed.matrices=T, ...)

      return(invisible(self$test.results[[name]]))
    },

    estimateCommonExpressionShiftMagnitudes=function(sample.groups=self$sample.groups, cell.groups=self$cell.groups, n.cells=NULL, n.randomizations=50, n.subsamples=30, min.cells=10, n.cores=self$n.cores, verbose=self$verbose,  mean.trim=0.1, name='common.expression.shifts') {
      if(length(levels(sample.groups))!=2) stop("'sample.groups' must be a 2-level factor describing which samples are being contrasted")

      count.matrices <- extractRawCountMatrices(self$data.object, transposed=T)

      common.genes <- Reduce(intersect, lapply(count.matrices, colnames))
      count.matrices %<>% lapply(`[`, , common.genes)

      comp.matrix <- outer(sample.groups,sample.groups,'!='); diag(comp.matrix) <- FALSE


      # get a cell sample factor, restricted to the samples being contrasted
      cl <- lapply(count.matrices[names(sample.groups)], rownames)
      cl <- rep(names(cl), sapply(cl, length)) %>% setNames(unlist(cl)) %>%  as.factor()

      # cell factor
      cf <- cell.groups
      cf <- cf[names(cf) %in% names(cl)]

      # if(is.null(n.cells)) {
      #   n.cells <- min(table(cf)) # use the size of the smallest group
      #   if(verbose) cat('setting group size of ',n.cells,' cells for comparisons\n')
      # }

      if(!is.null(n.cells) && n.cells>=length(cf)) {
        if(n.subsamples>1) {
          warning('turning off subsampling, as n.cells exceeds the total number of cells')
          n.subsamples <- 1;
        }
      }

      if(verbose) cat('running',n.subsamples,'subsamples,',n.randomizations,'randomizations each ... \n')
      ctdll <- sccore:::plapply(1:n.subsamples,function(i) {
        # subsample cells

        # draw cells without sample stratification - this can drop certain samples, particularly those with lower total cell numbers
        if(!is.null(n.cells)) {
          cf <- tapply(names(cf),cf,function(x) {
            if(length(x)<=n.cells) { return(cf[x]) } else { setNames(rep(cf[x[1]],n.cells), sample(x,n.cells)) }
          })
        }

        # calculate expected mean number of cells per sample and aim to sample that
        n.cells.scaled <- max(min.cells,ceiling(n.cells/length(sample.groups)));
        cf <- tapply(names(cf),list(cf,cl[names(cf)]),function(x) {
          if(length(x)<=n.cells.scaled) { return(cf[x]) } else { setNames(rep(cf[x[1]],n.cells.scaled), sample(x,n.cells.scaled)) }
        })

        cf <- as.factor(setNames(unlist(lapply(cf,as.character)),unlist(lapply(cf,names))))

        # table of sample types and cells
        cct <- table(cf,cl[names(cf)])
        caggr <- lapply(count.matrices, conos:::collapseCellsByType, groups=as.factor(cf), min.cell.count=1)[names(sample.groups)]

        ctdl <- sccore::plapply(sccore:::sn(levels(cf)),function(ct) { # for each cell type
          tcm <- na.omit(do.call(rbind,lapply(caggr,function(x) x[match(ct,rownames(x)),])))
          tcm <- log(tcm/rowSums(tcm)*1e3+1) # log transform
          tcm <- t(tcm[cct[ct,rownames(tcm)]>=min.cells,,drop=F])

          # an internal function to calculate consensus change direction and distances between samples along this axis
          t.consensus.shift.distances <- function(tcm,sample.groups, useCpp=TRUE) {
            if(min(table(sample.groups[colnames(tcm)]))<1) return(NA); # not enough samples
            if(useCpp) {
              g1 <- which(sample.groups[colnames(tcm)]==levels(sample.groups)[1])-1
              g2 <- which(sample.groups[colnames(tcm)]==levels(sample.groups)[2])-1
              as.numeric(projdiff(tcm,g1,g2))
            } else {
              # calculate consensus expression shift
              s1 <- colnames(tcm)[sample.groups[colnames(tcm)]==levels(sample.groups)[1]]
              s2 <- colnames(tcm)[sample.groups[colnames(tcm)]==levels(sample.groups)[2]]
              dm <- do.call(rbind,lapply(s1,function(n1) {
                do.call(rbind,lapply(s2,function(n2) {
                  tcm[,n1]-tcm[,n2]
                }))
              }))
              dmm <- apply(dm,2,mean,trim=mean.trim)
              dmm <- dmm/sqrt(sum(dmm^2)) # normalize

              # project samples and calculate distances
              as.numeric(dm %*% dmm)
            }
          }

          # true distances
          tdist <- t.consensus.shift.distances(tcm,sample.groups)
          # randomized distances

          rdist <- lapply(1:n.randomizations,function(i) {
            t.consensus.shift.distances(tcm,as.factor(setNames(sample(as.character(sample.groups)), names(sample.groups))))
          })

          # normalize true distances by the mean of the randomized ones
          return(abs(tdist)/mean(abs(unlist(rdist))))
        },n.cores = 1,mc.preschedule = FALSE, progress=(n.subsamples<=1))

      },n.cores=ifelse(n.subsamples>1,n.cores,1), mc.preschedule=FALSE, progress=(verbose && n.subsamples>1))

      return(invisible(self$test.results[[name]] <- ctdll))

    },

    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
    #' @param name - results slot name (default: 'expression.shifts')
    #' @param show.jitter whether to show indiivudal data points (default: FALSE)
    #' @param jitter.alpha transparency value for the data points (default: 0.05)
    #' @param type - type of a plot "bar" (default) or "box"
    #' @param notch - whether to show notches in the boxplot version (default=TRUE)
    #' @param show.size.depenency whether to show mean vs. number of cells in a cell type instead
    #' @param font.size font size for the cell type labels in the size dependency plot
    #' @param show.regression whether to show a slope line in the size dependency plot
    #' @param show.regression whether to show a whiskers in the size dependency plot
    #' @return A ggplot2 object
    plotExpressionShiftMagnitudes=function(name="expression.shifts", type='box', notch = TRUE, show.jitter=TRUE, jitter.alpha=0.05, show.size.dependency=FALSE, show.whiskers=TRUE, show.regression=TRUE, font.size=5) {
      df <- private$getResults(name)$df
      df <- data.frame(cell=df$Type, val=df$value)

      if(show.size.dependency) {
        plotCellTypeSizeDep(df, self$cell.groups, palette=self$cell.groups.palette,ylab='normalized expression distance', yline=NA, show.whiskers=show.whiskers, show.regression=show.regression)
      } else {
        plotMeanMedValuesPerCellType(df,show.jitter=show.jitter,jitter.alpha=jitter.alpha, notch=notch, type=type, palette=self$cell.groups.palette, ylab='normalized expression distance')
      }
    },

    ##' Plot common expression shift estimates across cell types
    ##'
    ##' @param name result slot name (default: common.expression.shifts
    ##' @param show.subsampling.variability  - whether the spread should illustrate subsampling variability instead of the inter-sample variability (default: FALSE)
    ##' @param show.jitter whether to show indiivudal data points (default: FALSE)
    ##' @param jitter.alpha transparency value for the data points (default: 0.05)
    ##' @param type - type of a plot "bar" (default) or "box"
    ##' @param notch - whether to show notches in the boxplot version (default=TRUE)
    ##' @param show.size.depenency whether to show mean vs. number of cells in a cell type instead
    ##' @param font.size font size for the cell type labels in the size dependency plot
    ##' @param show.regression whether to show a slope line in the size dependency plot
    ##' @param show.regression whether to show a whiskers in the size dependency plot
    ##' @return A ggplot2 object
    plotCommonExpressionShiftMagnitudes=function(name='common.expression.shifts', show.subsampling.variability=FALSE, show.jitter=FALSE, jitter.alpha=0.05, type='bar', notch=TRUE, show.size.dependency=FALSE, show.whiskers=TRUE, show.regression=TRUE, font.size=5) {
      res <- private$getResults(name)
      cn <- setNames(names(res[[1]]),names(res[[1]]))
      if(show.subsampling.variability) { # average across patient pairs
        if(length(res)<2) stop('the result has only one subsample; please set show.sampling.variability=FALSE')
        df <- do.call(rbind,lapply(res,function(d) data.frame(val=unlist(lapply(d,mean)),cell=names(d))))
      } else { # average across subsampling rounds
        df <- do.call(rbind,lapply(cn,function(n) data.frame(val=colMeans(do.call(rbind,lapply(res,function(x) x[[n]]))),cell=n)))
      }

      if(show.size.dependency) {
        plotCellTypeSizeDep(df, self$cell.groups, palette=self$cell.groups.palette,ylab='common expression distance', yline=NA, show.whiskers=show.whiskers, show.regression=show.regression)
      } else {
        plotMeanMedValuesPerCellType(df,show.jitter=show.jitter,jitter.alpha=jitter.alpha, notch=notch, type=type, palette=self$cell.groups.palette, ylab='common expression distance')
      }
    },

    ### Differential expression

    #' @description Estimate differential gene expression per cell type between conditions
    #' @param cell.groups factor specifying cell types (default=NULL)
    #' @param cooks.cutoff cooksCutoff for DESeq2 (default=F)
    #' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
    #' @param common.genes Only investigate common genes across cell groups (default=F)
    #' @param test which DESeq2 test to use (options: "LRT" (default), "Wald")
    #' @param cooks.cutoff cooksCutoff for DESeq2 (default=FALSE)
    #' @param min.cell.count minimum number of cells that need to be present in a given cell type in a given sample in order to be taken into account (default=10)
    #' @param max.cell.count maximal number of cells per cluster per sample to include in a comparison (useful for comparing the number of DE genes between cell types) (default: Inf)
    #' @param independent.filtering independentFiltering parameter for DESeq2 (default=FALSE)
    #' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
    #' @param return.matrix Return merged matrix of results (default=TRUE)
    #' @param name slot in which to save the results (default: 'de')
    #' @return A list of DE genes
    estimatePerCellTypeDE=function(cell.groups = self$cell.groups, sample.groups = self$sample.groups, ref.level = self$ref.level, common.genes = FALSE, n.cores = self$n.cores, cooks.cutoff = FALSE, min.cell.count = 10, max.cell.count= Inf, test='Wald', independent.filtering = FALSE, cluster.sep.chr = "<!!>", return.matrix = T, verbose=self$verbose, name ='de') {
      if(!is.list(sample.groups)) {
        sample.groups <- list(names(sample.groups[sample.groups == ref.level]),
                              names(sample.groups[sample.groups != ref.level])) %>%
          setNames(c(ref.level, self$target.level))
      }

      self$test.results[[name]] <- extractRawCountMatrices(self$data.object, transposed=T) %>%
        estimatePerCellTypeDE(cell.groups = cell.groups, sample.groups = sample.groups, ref.level = ref.level, n.cores = n.cores,
                              cooks.cutoff = cooks.cutoff, min.cell.count = min.cell.count, max.cell.count=max.cell.count, test=test, independent.filtering = independent.filtering,
                              cluster.sep.chr = cluster.sep.chr, return.matrix = return.matrix)
      return(invisible(self$test.results[[name]]))
    },

    #' @description Estimate differential gene expression per cell type between conditions
    #' @param cell.groups factor specifying cell types (default=NULL)
    #' @param cooks.cutoff cooksCutoff for DESeq2 (default=F)
    #' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
    #' @param common.genes Only investigate common genes across cell groups (default=F)
    #' @param test which DESeq2 test to use (options: "LRT" (default), "Wald")
    #' @param cooks.cutoff cooksCutoff for DESeq2 (default=FALSE)
    #' @param min.cell.count minimum number of cells that need to be present in a given cell type in a given sample in order to be taken into account (default=10)
    #' @param max.cell.count maximal number of cells per cluster per sample to include in a comparison (useful for comparing the number of DE genes between cell types) (default: Inf)
    #' @param independent.filtering independentFiltering parameter for DESeq2 (default=FALSE)
    #' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
    #' @param return.matrix Return merged matrix of results (default=TRUE)
    #' @param resampling.method which resampling method should be used "loo" for leave-one-out or "bootstrap", (default:NULL no resampling)
    #' @param name slot in which to save the results (default: 'de')
    #' @return A list of DE genes
    estimatePerCellTypeDEnew=function(cell.groups = self$cell.groups,
                                      sample.groups = self$sample.groups,
                                      ref.level = self$ref.level,
                                      target.level = self$target.level,
                                      common.genes = FALSE,
                                      n.cores = self$n.cores,
                                      cooks.cutoff = FALSE,
                                      min.cell.count = 10,
                                      max.cell.count= Inf,
                                      independent.filtering = FALSE,
                                      cluster.sep.chr = "<!!>",
                                      verbose=self$verbose,
                                      name ='de',
                                      test='DESeq2.Wald',
                                      resampling.method=NULL, # default - without resampling
                                      max.resamplings=30,
                                      seed.resampling=239, # shouldn't this be external?
                                      covariates = NULL) {

      if(!is.list(sample.groups)) {
        s.groups <- list(names(sample.groups[sample.groups == ref.level]),
                         names(sample.groups[sample.groups != ref.level])) %>%
          setNames(c(ref.level, target.level))
      } else {
        s.groups = sample.groups
      }

      possible.tests <- c('DESeq2.Wald', 'DESeq2.LRT', 'edgeR',
                          'Wilcoxon.edgeR', 'Wilcoxon.DESeq2', 'Wilcoxon.totcount',
                          't-test.edgeR', 't-test.DESeq2', 't-test.totcount',
                          'limma-voom')

      # Default test for DESeq2 is Wald
      if(tolower(test) == tolower('DESeq2')) test = paste(test, 'Wald', sep='.')
      # Default normalization for Wilcoxon and t-test is edgeR
      if(tolower(test) %in% tolower(c('Wilcoxon', 't-test')) )  test = paste(test, 'edgeR', sep='.')

      if(!(tolower(test) %in% tolower(possible.tests))) stop(paste('Test', test, 'is not supported. Available tests:',paste(possible.tests,collapse=', ')))

      test <- possible.tests[tolower(test) == tolower(possible.tests)]
      message(paste0(c('DE method ', test, ' is used'), collapse = ''))

      # s.groups.new contains list of case/control groups of samples to run DE on.
      # First element in s.groups.new corresponds to the initial grouping.

      s.groups.new = list(initial = s.groups)
      # If resampling is defined, new contrasts will append to s.groups.new
      if (is.null(resampling.method)){
      } else if (resampling.method == 'loo') {
        s.groups.new = c(s.groups.new, lapply(unlist(s.groups), function(name)
          lapply(s.groups, function(group) setdiff(group, name)))  %>%
            setNames(unlist(s.groups)))
      } else if (resampling.method == 'bootstrap') {
        if(max.resamplings < 2) {
          warning('Bootstrap was not applied, because the number of resamplings was less than 2')
        } else {
          n.bootstrap <- max.resamplings
          set.seed(seed.resampling)
          s.groups.new <- c(s.groups.new, lapply(setNames(1:n.bootstrap,paste0('bootstrap.',1:n.bootstrap)),function(i) lapply(s.groups,function(x) sample(x,length(x),replace=T))))
        }
      } else stop(paste('Resampling method', resampling.method, 'is not supported'))

      raw.mats <- extractRawCountMatrices(self$data.object, transposed=T)

      de.res = sccore::plapply(names(s.groups.new), function(resampling.name) {

        estimatePerCellTypeDEmethods(raw.mats=raw.mats,
                                     cell.groups = cell.groups,
                                     s.groups = s.groups.new[[resampling.name]],
                                     ref.level = ref.level,
                                     target.level = target.level,
                                     common.genes = common.genes,
                                     cooks.cutoff = cooks.cutoff,
                                     min.cell.count = min.cell.count,
                                     max.cell.count = max.cell.count,
                                     independent.filtering = independent.filtering,
                                     n.cores = ifelse(length(s.groups.new)>=n.cores,1,n.cores),
                                     cluster.sep.chr = cluster.sep.chr,
                                     return.matrix = ifelse(resampling.name == 'initial',T,F),
                                     verbose = length(s.groups.new)<n.cores,
                                     useT = useT,
                                     minmu = minmu,
                                     test = test,
                                     meta.info = covariates)
      },n.cores=ifelse(length(s.groups.new)>=n.cores,n.cores,1),progress=length(s.groups.new)>=n.cores) %>% setNames(names(s.groups.new)) # parallelize the outer loop if subsampling is on


      # if resampling: calculate median and variance on ranks after resampling
      if(length(de.res) > 1) {
        var.to.sort = 'pvalue' # Variable to calculate ranks
        for(cell.type in names(de.res[[1]])) {
          genes.init <- genes.common <- rownames(de.res[[1]][[cell.type]]$res)
          mx.stat <- matrix(nrow = length(genes.common), ncol = 0, dimnames = list(genes.common,c()))
          for(i in 2:length(de.res)){
            if(!(cell.type %in% names(de.res[[i]]))) next
            genes.common = intersect(genes.common, rownames(de.res[[i]][[cell.type]]))
            mx.stat = cbind(mx.stat[genes.common,], de.res[[i]][[cell.type]][genes.common, var.to.sort])
          }

          mx.stat = apply(mx.stat, 2, rank)
          stab.mean.rank = rowMeans(mx.stat) # stab - for stability
          stab.median.rank = apply(mx.stat, 1, median)
          stab.var.rank = apply(mx.stat, 1, var)

          de.res[[1]][[cell.type]]$res$stab.median.rank = stab.median.rank[genes.init]
          de.res[[1]][[cell.type]]$res$stab.mean.rank = stab.mean.rank[genes.init]
          de.res[[1]][[cell.type]]$res$stab.var.rank = stab.var.rank[genes.init]

          # Save subsamples
          de.res[[1]][[cell.type]]$subsamples <- lapply(de.res[2:length(de.res)], function(de) de[[cell.type]])
        }
      }

      # Overwrite 'de' slot with the last result
      return(invisible(self$test.results[[name]] <- de.res[[1]]))
    },

    #' @description Plot number of significant DE genes as a function of number of cells
    #' @param name results slot in which the DE results should be stored (default: 'de')
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param palette cell group palette (default: stored $cell.groups.palette)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="none")
    #' @param label Show labels on plot (default=T)
    #' @param p.adjust.cutoff Adjusted P cutoff (default=0.05)
    #' @return A ggplot2 object
    plotDEGenes=function(name='de', cell.groups = self$cell.groups, legend.position = "none", label = TRUE, p.adj = 0.05, size=4, palette=self$cell.groups.palette) {
      de.raw <- private$getResults(name, 'estimatePerCellTypeDE()')

      if(class(de.raw[[1]]) == "list") de.raw %<>% lapply(`[[`, 1) # If estimatePerCellTypeDE was run with return.matrix = T
      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.raw)]
      gg <- sapply(de.raw, function(d) sum(d$padj <= p.adj)) %>%
        plotNCellRegression(cell.groups, x.lab="Number of cells", y.lab="Significant DE genes",
                            legend.position=legend.position, label=label, size=size, palette=palette) +
        geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="black", size=0.5)

      return(gg)
    },

    #' @description Plot number of significant DE genes
    #' @param name results slot in which the DE results should be stored (default: 'de')
    #' @param palette cell group palette (default: stored $cell.groups.palette)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="none")
    #' @param label Show labels on plot (default=T)
    #' @param pvalue.cutoff P value cutoff (default=0.05)
    #' @param p.adjust whether the cutoff should be based on the adjusted P value (default: TRUE)
    #' @param show.resampling.results whether to show uncertainty based on resampling results (default: TRUE)
    #' @param show.size.depenency whether to show mean vs. number of cells in a cell type instead (default: FALSE)
    #' @param font.size font size for the cell type labels in the size dependency plot
    #' @param show.regression whether to show a slope line in the size dependency plot (default:TRUE)
    #' @param show.whiskers whether to show a whiskers in the size dependency plot (default:TRUE)
    #' @return A ggplot2 object
    plotNumberOfDEGenes=function(name='de', legend.position="none", label=TRUE, p.adjust=TRUE, pvalue.cutoff=0.05, show.resampling.results=TRUE, show.jitter=FALSE, jitter.alpha=0.05, type='bar', notch=TRUE, show.size.dependency=FALSE, show.whiskers=TRUE, show.regression=TRUE, font.size=5) {

      de.raw <- private$getResults(name, 'estimatePerCellTypeDE()')

      if(show.resampling.results) {
        if(!all(unlist(lapply(de.raw,function(x) !is.null(x$subsamples))))) {
          warning("resampling results are missing for at least some cell types, falling back to point estimates. Please rerun estimatePerCellTypeDE() with resampling='bootstrap' or resampling='loo'")
          rl <- lapply(de.raw,function(x) x$res)
        } else {
          rl <- setNames( unlist(lapply(de.raw,function(x) x$subsamples),recursive=F), rep(names(de.raw),unlist(lapply(de.raw,function(x) length(x$subsamples)))))
        }
      } else {
        rl <- lapply(de.raw,function(x) x$res)
      }
      # convert to dataframe for plotting
      df <- do.call(rbind,lapply(1:length(rl),function(i) {
        if(p.adjust) {
          ndiff <- sum(na.omit(rl[[i]]$padj<=pvalue.cutoff))
        } else {
          ndiff <- sum(na.omit(rl[[i]]$pvalue<=pvalue.cutoff))
        }
        data.frame(cell=names(rl)[i],val=ndiff,stringsAsFactors=FALSE)
      }))

      if(show.size.dependency) {
        plotCellTypeSizeDep(df, self$cell.groups, palette=self$cell.groups.palette,ylab='number of DE genes', yline=NA, show.whiskers=show.whiskers, show.regression=show.regression)
      } else {
        plotMeanMedValuesPerCellType(df,show.jitter=show.jitter,jitter.alpha=jitter.alpha, notch=notch, type=type, palette=self$cell.groups.palette, ylab='number of DE genes',yline=NA)
      }
    },

    #' @description Save DE results as JSON files
    #' @param saveprefix Prefix for created files (default=NULL)
    #' @param dir.name Name for directory with results (default="JSON")
    #' @param de.raw List of DE results
    #' @param ref.level Reference level in 'sample.groups', e.g., ctrl, healthy, wt (default=NULL)
    #' @param gene.metadata (default=NULL)
    #' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default="<!!>")
    saveDEasJSON=function(saveprefix = NULL, dir.name = "JSON", de.raw = NULL, sample.groups = self$sample.groups, ref.level = self$ref.level, gene.metadata = NULL, cluster.sep.chr = "<!!>", verbose = T) {
      if (is.null(de.raw)) {
        de.raw <- private$getResults("de", "estimatePerCellTypeDE")
      }

      if(!is.list(sample.groups)) {
        sample.groups <- list(names(sample.groups[sample.groups == ref.level]),
                              names(sample.groups[sample.groups != ref.level])) %>%
          setNames(c(ref.level, self$target.level))
      }

      if(class(de.raw[[1]]) != "list") stop("Please rerun 'estimatePerCellTypeDE' with return.matrix=T")

      if(is.null(saveprefix)) saveprefix <- ""

      saveDEasJSON(de.raw = de.raw, saveprefix = saveprefix, dir.name = dir.name, gene.metadata = gene.metadata, cluster.sep.chr = cluster.sep.chr, sample.groups = sample.groups, verbose = verbose)
    },

    #' @description Plot number of highly-expressed DE genes as a function of number of cells
    #' @param de.filter Filtered DE genes, results from prepareOntologyData (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="none")
    #' @param label Show labels on plot (default=T)
    #' @return A ggplot2 object
    plotFilteredDEGenes=function(de.filter=self$test.results$ontology$de.filter, cell.groups = self$cell.groups, legend.position = "none", label = T) {
      if(is.null(de.filter)) stop("Please run 'estimatePerCellTypeDE' first.")

      cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.filter)]

      gg <- sapply(de.filter, length) %>%
        plotNCellRegression(cell.groups, x.lab="Number of cells", y.lab="Highly-expressed DE genes", legend.position=legend.position, label=label) +
        geom_smooth(method=MASS::rlm, formula=y~x, se=0, color="black", size=0.5)

      return(gg)
    },

    ### Onthology analysis

    #' @description  Filter and prepare DE genes for ontology calculations
    #' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
    #' @param n.top.genes Number of most different genes to take as input. If less are left after filtering for p.adj.cutoff, additional genes are included. To disable, set n.top.genes=0 (default=1e2)
    #' @param p.adj Cutoff for filtering highly-expressed DE genes (default=0.05)
    #' @param expr.cutoff Cutoff for cells per group expressing a DE gene, i.e., cutoff for highly-expressed genes (default=0.05)
    #' @param de.raw Differentially expressed genes per cell group, results from estimatePerCellTypeDE (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param universe Only set this if a common background gene set is desired for all cell groups (default: NULL)
    #' @param transposed Whether count matrices should be transposed (default=T)
    #' @return A list containing DE gene IDs, filtered DE genes, and input DE genes
    prepareOntologyData=function(org.db, p.adj=1, expr.cutoff=0.05, de.raw=NULL, cell.groups=self$cell.groups, universe=NULL, transposed=TRUE, verbose=self$verbose, n.cores=self$n.cores) {
      if (is.null(de.raw)) {
        de.raw <- private$getResults("de", "estimatePerCellTypeDE")
      }

      # If estimatePerCellTypeDE was run with return.matrix = T, remove matrix before calculating
      if(class(de.raw[[1]]) == "list") de.raw %<>% lapply(`[[`, 1)

      if(p.adj < 1) warning("You are filtering based on adj. P value through the 'p.adj' parameter. We do not recommend this. Proceed with caution.")

      if(is.null(cell.groups))
        stop("'cell.groups' must be provided either to Cacoa constructor or to this method.")

      self$test.results[["gene.ids"]] <- extractRawCountMatrices(self$data.object, transposed = transposed) %>%
        prepareOntologyData(org.db = org.db, p.adj = p.adj, expr.cutoff = expr.cutoff, de.raw = de.raw, cell.groups = cell.groups, universe = universe, transposed = transposed, verbose = verbose, n.cores = n.cores)
      return(invisible(self$test.results[["gene.ids"]]))
    },

    #' @description  Plot embedding
    #' @param embedding A cell embedding to use (two-column data frame with rownames corresponding to cells) (default: stored embedding object)
    #' @param plot.theme plot theme to use (default: ggplot2::theme_bw())
    #' @param ... other parameters are passed to \link[sccore:embeddingPlot]{embeddingPlot}
    plotEmbedding=function(embedding=self$embedding, plot.theme=ggplot2::theme_bw(), show.legend=TRUE, ...) {
      if(is.null(embedding)) stop("embedding must be provided to Cacoa constructor or to this method.")
      sccore::embeddingPlot(embedding, plot.theme=plot.theme, show.legend=show.legend, ...)
    },

    #' @description Estimate ontology terms based on DEs
    #' @param type Ontology type, either GO (gene ontology) or DO (disease ontology). Please see DOSE package for more information (default="GO")
    #' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
    #' @param de.gene.ids List containing DE gene IDs, and filtered DE genes (default: stored list, results from prepareOntologyData)
    #' @param go.environment Extracted GO environment. If set to NULL, the environment will be re-extracted (default: stored environment)
    #' @param p.adjust.method Method for calculating adj. P. Please see DOSE package for more information (default="BH")
    #' @param readable Mapping gene ID to gene name (default=T)
    #' @param min.genes Minimum number of input genes overlapping with ontologies (default=0)
    #' @param qvalue.cutoff Q value cutoff, please see clusterProfiler package for more information (default=0.2)
    #' @param min.gs.size Minimal geneset size, please see clusterProfiler package for more information (default=5)
    #' @param max.gs.size Minimal geneset size, please see clusterProfiler package for more information (default=5e2)
    #' @return A list containing a list of terms per ontology, and a data frame with merged results
    estimateOntology=function(type = "GO", org.db, n.top.genes = 500, de.gene.ids=NULL, go.environment=self$test.results$GO$go.environment, p.adjust.method="BH",
                              readable=TRUE, verbose=TRUE, qvalue.cutoff=0.2, min.gs.size=10, max.gs.size=5e2, ...) {
      if(!is.null(type) & !type %in% c("GO", "DO", "GSEA"))
        stop("'type' must be 'GO', 'DO', or 'GSEA'.")

      if(is.null(de.gene.ids)) {
        de.gene.ids <- private$getResults("gene.ids", "prepareOntologyData")
      }

      self$test.results[[type]] <- estimateOntology(type=type, org.db=org.db, n.top.genes=n.top.genes, de.gene.ids=de.gene.ids, go.environment=go.environment,
                                                    verbose=verbose, qvalue.cutoff=qvalue.cutoff, pAdjustMethod=p.adjust.method, readable=readable,
                                                    minGSSize=min.gs.size, maxGSSize=max.gs.size, ...)
      return(invisible(self$test.results[[type]]))
    },

    #' @description Estimate ontology families
    #' @return A list of results
    estimateOntologyFamilies=function(type = "GO", p.adj = 0.05) {
      # TODO: Checks
      if (!requireNamespace("GOfuncR", quietly = TRUE)) stop("You need 'GOfuncR' to perform the ontology family analysis.")
      if(!type %in% c("GO","DO","GSEA")) stop("'type' must be 'GO', or 'GSEA'.")

      # TODO: Test DO
      if(type == "GO") {
        ont.list <- self$test.results[[type]]$res %>%
          lapply(lapply, lapply, function(x) {
            tmp <- x@result %>% filter(p.adjust <= p.adj)
            if(nrow(tmp) > 0) return(tmp)
          }) %>%
          lapply(lapply, plyr::compact) %>%
          lapply(plyr::compact)
      } else {
        ont.list <- self$test.results[[type]]$res %>%
          lapply(lapply, function(x) {
            tmp <- x@result %>% filter(p.adjust <= p.adj)
            if(nrow(tmp) > 0) return(tmp)
          }) %>%
          lapply(plyr::compact)
      }
      self$test.results[[type]]$families <- estimateOntologyFamilies(ont.list = ont.list, type = type, p.adj = p.adj)
    },

    #' @description Estimate Gene Ontology clusters
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all'. Default: "up".
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types), "GO" or "DO". Default: "GO".
    #' @param name Name of the field to store the results. Default: `cacoa:::getOntClustField(type, subtype, genes)`.
    #' @param ind.h Cut height for hierarchical clustering of terms per cell type.
    #' Approximately equal to the fraction of genes, shared between the GOs. Default: 0.66.
    #' @param total.h Cut height for hierarchical clustering of GOs across all subtypes.
    #' Approximately equal to the fraction of subtypes for which two GOs should belong to the same cluster. Default: 0.5.
    #' @return List containing:
    #'   - `df`: data.frame with information about individual gene ontolodies and columns `Cluster` and `ClusterName` for the clustering info
    #'   - `hclust`: the object of class \link[stats:hclust]{hclust} with hierarchical clustering of GOs across all subtypes
    estimateOntologyClusters=function(type="GO", subtype=NULL, genes="all", ind.h=0.66, total.h=0.5, verbose=self$verbose,
                                      p.adj=p.adj, min.genes=min.genes, name=getOntClustField(type, subtype, genes)) {
      ont.df <- private$getOntologyPvalueResults(genes=genes, type=type, p.adj=p.adj, min.genes=min.genes)
      clust.mat <- ont.df %>% split(.$Group) %>% clusterIndividualGOs(cut.h=ind.h) %>%
        as.matrix() %>% t()

      apply.fun <- if (verbose) pbapply::pbapply else apply
      cl.dists <- apply.fun(clust.mat, 2, function(ct1) apply(clust.mat, 2, function(ct2) {
        mask <- !is.na(ct1) & !is.na(ct2)
        if (sum(mask) == 0) 1 else (1 - mean(ct1[mask] == ct2[mask]))
      }))

      cl.clusts <- as.dist(cl.dists) %>% hclust(method="average")
      clusts <- cutree(cl.clusts, h=total.h)
      ont.df$Cluster <- clusts[as.character(ont.df$Description)]

      name.per.clust <- ont.df %>% group_by(Cluster, Description) %>% summarise(pvalue=exp(mean(log(pvalue)))) %>%
        split(.$Cluster) %>% sapply(function(df) df$Description[which.min(df$pvalue)])

      ont.df$ClusterName <- name.per.clust[ont.df$Cluster]
      self$test.results[[name]] <- list(df=ont.df, hclust=cl.clusts)

      return(invisible(self$test.results[[name]]))
    },

    #' @description Bar plot of ontology terms per cell type
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="all")
    #' @param type Ontology, must be either "GO" or "DO" (default="GO")
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @return A ggplot2 object
    plotOntologyDistribution=function(genes="all", type="GO", p.adj=0.05, min.genes=1, cell.groups=self$cell.groups) {
      if (length(genes) > 0) {
        ont.res <- genes %>% setNames(., .) %>% lapply(private$getOntologyPvalueResults, type=type, p.adj=p.adj, min.genes=min.genes)

        if (type != "GSEA") {
          classes <- sapply(ont.res[genes], class)
          if(any(classes == "character")) {
            message(paste0("No significant results found for genes = '",names(classes[classes == "character"]),"'."))
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
        p_df <- table(ont.res$Group, ont.res$Type, ont.res$direction) %>%
          as.data.frame() %>%
          setNames(c("Group","Type","direction","N")) %>%
          dplyr::arrange(Group)

        if(length(unique(p_df$direction)) > 1) {
          gg <- ggplot(p_df, aes(x=Group, y=N, fill=Type, group=Group)) +
            geom_bar(stat="identity") +
            facet_grid(~direction, switch="x")
        } else {
          gg <- ggplot(p_df) +
            geom_bar(aes(x=Group, y=N, fill=Type), stat="identity")
        }
      } else if(type=="DO") {
        p_df <- table(ont.res$Group, ont.res$direction) %>%
          as.data.frame() %>%
          setNames(c("Group","direction","N")) %>%
          dplyr::arrange(Group)

        if(length(unique(p_df$direction)) > 1) {
          gg <- ggplot(p_df) +
            geom_bar(aes(x=Group, y=N, fill=direction), stat="identity", position="dodge") +
            labs(fill="Gene set")
        } else {
          gg <- ggplot(p_df) +
            geom_bar(aes(x=Group, y=N), stat="identity")
        }
      } else if(type == "GSEA") {
        p_df <- table(ont.res$Group, ont.res$Type) %>%
          as.data.frame() %>%
          setNames(c("Group","Type","N")) %>%
          dplyr::arrange(Group)

          gg <- ggplot(p_df) +
            geom_bar(aes(x=Group, y=N, fill=Type), stat="identity")
      }

      gg <- gg +
        scale_y_continuous(expand=c(0, 0)) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
              legend.position="right") +
        labs(x="", y=paste0("No. of ",type," terms"))
      return(gg)
    },

    #' @description Plot ontology terms as a function of both number of DE genes, and number of cells.
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default='all')
    #' @param type Ontology, must be either "GO" or "DO" (default="GO")
    #' @param de.filter Filtered DE genes, results from prepareOntologyData (default: stored list)
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param label.x.pos Plot label position on x axis (default=0.01)
    #' @param label.y.pos Plot label position on y axis (default=1)
    #' @param scale Scaling of plots, adjust if e.g. label is misplaced. See \link[cowplot:plot_grid]{cowplot::plot_grid} for more info (default=1.0)
    #' @param ... parameters forwarded to \link{plotNCellRegression}
    #' @return A ggplot2 object
    plotOntologyTerms=function(genes='all', type="GO", p.adj=0.05, min.genes=1, de.filter=self$test.results$gene.ids,
                               cell.groups=self$cell.groups, label.x.pos=0.01, label.y.pos=1, scale=1.0, ...) {
      ont.res <- private$getOntologyPvalueResults(genes=genes, type=type, p.adj=p.adj, min.genes=min.genes)

      if(length(unique(ont.res$Group))==1) stop("The input only contains one cell type.")
      cell.groups <- table(cell.groups) %>% .[names(.) %in% names(de.filter)]

      n.go.per.type <- table(ont.res$Group) %>% c()
      n.de.per.type <- sapply(de.filter, function(x) length(x$all))

      y.lab <- paste("Number of", type, "terms")
      pg <- cowplot::plot_grid(
        plotNCellRegression(n.go.per.type, n.de.per.type, x.lab="Number of highly-expressed DE genes",
                            y.lab=y.lab, legend.position="none", label=TRUE, ...),
        plotNCellRegression(n.go.per.type, cell.groups, y.lab=y.lab, legend.position="none", label=TRUE, ...),
        ncol=1, labels=c("a", "b"), label_x=label.x.pos, label_y=label.y.pos, scale=scale
      )

      return(pg)
    },

    #' @description Plot a dotplot of ontology terms with adj. P values for a specific cell subgroup
    #' @param genes Specify which genes to plot, can either be 'down', 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types), "GO" or "DO" (default="GO")
    #' @param cell.subgroup Cell group to plot (default=NULL)
    #' @param n Number of ontology terms to show. Not applicable when order is 'unique' or 'unique-max-row' (default=10)
    #' @param p.adj Adjusted P cutoff (default=0.05)
    #' @param log.colors Use log10 p-values for coloring (default=FALSE)
    #' @return A ggplot2 object
    plotOntology = function(cell.subgroups, plot = "dot", genes = "up", type = "GO", subtype = "BP", n = 20, p.adj = 0.05, min.genes = 1, ...) {
      # Checks
      if (!requireNamespace("enrichplot", quietly = TRUE)) stop("You need 'enrichplot' to plot.")
      if(is.null(type) || (!type %in% c("GO","DO","GSEA"))) stop("'type' must be 'GO', 'DO', or 'GSEA'.")
      if(is.null(subtype) || (!subtype %in% c("BP","CC","MF"))) stop("'subtype' must be 'BP', 'CC', or 'MF'.")
      if(is.null(genes) || (!genes %in% c("down","up","all"))) stop("'genes' must be 'down', 'up', or 'all'.")
      if(is.null(cell.subgroups)) stop("Please define which cells to plot using 'cell.subgroups'.")
      if(type == "GSEA" & plot == "bar") stop("No 'enrichplot' method exists for making barplots of GSEA results.")

      # Extract results
      ont.res <- self$test.results[[type]]$res
      if(is.null(ont.res)) stop(paste0("No results found for ",type,"."))
      if(!cell.subgroups %in% names(ont.res)) stop("'cell.subgroups' not found in results.")
      ont.res %<>% .[[cell.subgroups]]
      if(type != "DO") ont.res %<>% .[[subtype]]
      if(is.null(ont.res)) stop(paste0("No results found for ",type,", ",subtype," for ",cell.subgroups,"."))
      if(type %in% c("GO","DO")) ont.res %<>% .[[genes]]
      if(is.null(ont.res)) stop(paste0("No results found for ",genes," genes for ",type,", ",subtype," for ",cell.subgroups,"."))

      # Prepare data
      df <- ont.res@result %>% filter(p.adjust <= p.adj)
      if(nrow(df) == 0) stop(paste0("Nothing to plot. Try relaxing 'p.adj'. The lowest adj. P value is ",min(ont.res@result$p.adjust) %>% formatC(digits = 3),"."))

      # Allow plotting of terms with p.adj > 0.05
      if(p.adj != 0.05) {
        ont.res@pvalueCutoff = 1
        ont.res@qvalueCutoff = 1
      }

      if(min.genes > 1) {
        idx <- df$GeneRatio %>%
          strsplit("/", fixed=T) %>%
          sapply(`[[`, 1)

        df <- df[idx > min.genes,]
      }
      ont.res@result <- df

      # Plot
      if(plot == "dot")
        return(dotplot(ont.res, showCategory=n, orderBy="x", ...))

      if(plot == "bar")
        return(barplot(ont.res, showCategory=n, ...))

      stop("Unknown plot type: ", plot)
    },

    #' @description Plot a heatmap of ontology P values per cell type
    #' @param genes Specify which genes to plot, can either be 'down' for downregulated genes, 'up' or 'all' (default="up")
    #' @param type Ontology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default="GO")
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="left")
    #' @param selection Order of rows in heatmap. Can be 'unique' (only show terms that are unique for any cell type); 'common' (only show terms that are present in at least two cell types); 'all' (all ontology terms) (default="all")
    #' @param n Number of terms to show (default=10)
    #' @param clusters Whether to show GO clusters or raw GOs (default=TRUE)
    #' @param cluster.name Field with the results for GO clustering. Ignored if `clusters == FALSE`.
    #' @param cell.subgroups Cell groups to plot (default=NULL)
    #' @param color.range vector with two values for min/max values of p-values
    #' @param ... parameters forwarded to \link{plotHeatmap}
    #' @return A ggplot2 object
    plotOntologyHeatmap=function(genes="up", type="GO", subtype="BP", min.genes=1, p.adj=0.05, legend.position="left", selection="all", n=10,
                                 clusters=TRUE, cluster.name=NULL, cell.subgroups=NULL, color.range=NULL, ...) {
      # Checks
      if(!is.null(cell.subgroups) && (length(cell.subgroups) == 1))
        stop("'cell.subgroups' must contain at least two groups. Please use plotOntology instead.")

      if(is.null(selection) || (!selection %in% c("unique","common","all")))
        stop("'selection' must be one of the following: 'unique', 'common', or 'all'.")

      # Extract results
      if (!clusters) {
        ont.sum <- private$getOntologyPvalueResults(genes=genes, type=type, subtype=subtype, cell.subgroups=cell.subgroups,
                                                    p.adj=p.adj, min.genes=min.genes) %>%
          groupOntologiesByCluster(field="Description")
      } else {
        name <- if (is.null(cluster.name)) getOntClustField(type, subtype, genes) else cluster.name
        if (is.null(self$test.results[[name]])) {
          if (!is.null(cluster.name))
            stop("Can't find the results for ", cluster.name) # stop if user specified a wrong cluster.name

          warning("Can't find the results for ", name, ". Running estimateOntologyClusters()...\n")
          ont.sum <- self$estimateOntologyClusters(type=type, genes=genes, name=name, p.adj=p.adj, min.genes=min.genes)$df
        } else {
          ont.sum <- self$test.results[[name]]$df
        }

        ont.sum %<>% groupOntologiesByCluster(field="ClusterName")
      }

      if(selection=="unique") {
        ont.sum %<>% .[rowSums(abs(.) > 0) == 1,]
      } else if(selection=="common") {
        ont.sum %<>% .[rowSums(abs(.) > 0) > 1,]
      }
      if(nrow(ont.sum) == 0) stop("Nothing to plot. Try another selection.")

      # Plot
      gg <- ont.sum %>%
        .[, colSums(abs(.)) > 0] %>%
        .[match(rowSums(.)[rowSums(abs(.)) > 0] %>% .[order(., decreasing=TRUE)] %>% names, rownames(.)),] %>%
        tail(n) %>%
        plotHeatmap(legend.position=legend.position, row.order=TRUE, color.range=color.range, ...) +
        getGeneScale(genes=genes, type="fill", high="white", limits=color.range)

      return(gg)
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
        pathway_df <- unique(ont.res$Group) %>%
          lapply(function(cell.group) {
            lapply(ont.res %>%
                     dplyr::filter(Group == cell.group) %>%
                     dplyr::pull(Type) %>%
                     as.factor() %>%
                     levels(), function(go) {
                       tibble::tibble(Pathway=ont.res %>%
                                        dplyr::filter(Group == cell.group) %>%
                                        dplyr::filter(Type==go) %>%
                                        dplyr::pull(Description),
                                      Group=cell.group,
                                      GO=go)
                     }) %>% dplyr::bind_rows()
          }) %>%
          dplyr::bind_rows()
      } else if(type=="DO") {
        pathway_df <- unique(ont.res$Group) %>%
          lapply(function(cell.group) {
            tibble::tibble(Pathway=ont.res %>%
                             dplyr::filter(Group == cell.group) %>%
                             dplyr::pull(Description),
                           Group=cell.group)
          }) %>%
          dplyr::bind_rows()
      }

      path_bin <- pathway_df %>%
        dplyr::select(Pathway, Group) %>%
        dplyr::mutate(X=1) %>%
        tidyr::spread(Pathway, X) %>%
        as.data.frame() %>%
        magrittr::set_rownames(.$Group) %>%
        .[, 2:ncol(.)] %>%
        as.matrix()
      path_bin[is.na(path_bin)] <- 0

      p_mat <- (1 - (path_bin %>% dist(method="binary") %>% as.matrix)) %>% pmin(0.5)
      t_tree <- dist(p_mat) %>% hclust()
      t_order <- t_tree %$% labels[order]
      t_cls <- cutree(t_tree, h=0.7) %>% .[t_order]
      t_cls[t_cls %in% names(which(table(t_cls) < 5))] <- max(t_cls) + 1
      t_cl_lengths <- rle(t_cls)$lengths %>% rev
      diag(p_mat) <- 1

      # Plot
      plotHeatmap(p_mat, color.per.group=NULL, row.order=t_order, col.order=rev(t_order), legend.title="Similarity") +
        scale_fill_distiller(palette="RdYlBu", limits=c(0, 0.5)) +
        geom_vline(aes(xintercept=x), data.frame(x=cumsum(t_cl_lengths)[t_cl_lengths > 1] + 0.5)) +
        geom_hline(aes(yintercept=x), data.frame(x=cumsum(t_cl_lengths)[t_cl_lengths > 1] + 0.5))
    },

    #' #' @description Plot ontology families heatmap
    #' #' @return A ggplot2 object
    #' plotOntologyFamilyHeatmap=function(genes = "up", type = "GO", subtype = "BP", min.genes = 1, selection = "unique", n = 20, legend.position = "right") {
    #'   # Checks
    #'   if(!is.null(cell.subgroups) && (length(cell.subgroups) == 1))
    #'     stop("'cell.subgroups' must contain at least two groups. Please use plotOntology instead.")
    #'
    #'   if(is.null(selection) || (!selection %in% c("unique","common","all")))
    #'     stop("'selection' must be one of the following: 'unique', 'common', or 'all'.")
    #'
    #'   # Extract results
    #'   if (!clusters) {
    #'     ont.sum <- private$getOntologyPvalueResults(genes=genes, type=type, subtype=subtype, cell.subgroups=cell.subgroups,
    #'                                                 p.adj=p.adj, min.genes=min.genes) %>%
    #'       groupOntologiesByCluster(field="Description")
    #'   } else {
    #'     name <- if (is.null(cluster.name)) getOntClustField(type, subtype, genes) else cluster.name
    #'     if (is.null(self$test.results[[name]])) {
    #'       if (!is.null(cluster.name))
    #'         stop("Can't find the results for ", cluster.name) # stop if user specified a wrong cluster.name
    #'
    #'       warning("Can't find the results for ", name, ". Running estimateOntologyClusters()...\n")
    #'       ont.sum <- self$estimateOntologyClusters(type=type, genes=genes, name=name, cell.subgroups=cell.subgroups,
    #'                                                p.adj=p.adj, min.genes=min.genes)$df
    #'     } else {
    #'       ont.sum <- self$test.results[[name]]$df
    #'     }
    #'
    #'     ont.sum %<>% groupOntologiesByCluster(field="ClusterName")
    #'   }
    #'
    #'   if(selection=="unique") {
    #'     ont.sum %<>% .[rowSums(. > 0) == 1,]
    #'   } else if(selection=="common") {
    #'     ont.sum %<>% .[rowSums(. > 0) > 1,]
    #'   }
    #'   if(nrow(ont.sum) == 0) stop("Nothing to plot. Try another selection.")
    #'
    #'   # Plot
    #'   gg <- ont.sum %>%
    #'     .[, colSums(abs(.)) > 0] %>%
    #'     .[match(rowSums(.)[rowSums(abs(.)) > 0] %>% .[order(., decreasing=TRUE)] %>% names, rownames(.)),] %>%
    #'     tail(n) %>%
    #'     plotHeatmap(legend.position=legend.position, row.order=TRUE, color.range=color.range, ...) +
    #'     getGeneScale(genes=genes, type="fill", high="white", limits=color.range)
    #'
    #'   return(gg)
    #' },

    #' @description Plot ontology family tree
    #' @return An Rgraphviz object
    plotOntologyFamily=function(type = "GO", cell.subgroups, family, genes = "up", subtype = "BP", plot.type = "complete", show.ids = FALSE,
                                string.length=18, legend.label.size = 1, legend.position = "topright", verbose = self$verbose, n.cores = self$n.cores) {
      #Checks
      if (!requireNamespace("GOfuncR", quietly = TRUE))
        stop("You need 'GOfuncR' to plot ontology families. Use `BiocManager::install('GOfuncR')`.")
      if (!requireNamespace("graph", quietly = TRUE))
        stop("You need 'graph' to plot ontology families. Use `BiocManager::install('graph')`.")
      if (!requireNamespace("Rgraphviz", quietly = TRUE))
        stop("You need 'Rgraphviz' to plot ontology families. Use `BiocManager::install('Rgraphviz')`.")

      if(!is.numeric(family)) stop("'family' must be numeric.")
      if(!is.null(plot.type) && !plot.type %in% c("complete","dense","minimal")) stop("'plot.type' must be 'complete', 'dense', or 'minimal'.")

      fam.name <- paste0("Family",family)
      ont.fam.res <- self$test.results[[type]]$families
      if(is.null(ont.fam.res)) stop(paste0("No results found for type '",type,"'."))
      ont.fam.res %<>% .[[cell.subgroups]]
      if(is.null(ont.fam.res)) stop(paste0("No results found for cell.subgroups '",cell.subgroups,"'."))
      ont.fam.res %<>% .[[subtype]]
      if(is.null(ont.fam.res)) stop(paste0("No results found for subtype '",subtype,"'."))
      ont.fam.res %<>% .[[genes]]
      if(is.null(ont.fam.res)) stop(paste0("No results found for genes '",genes,"'."))
      if(!fam.name %in% names(ont.fam.res$families)) stop("'family' not in 'ont.fam.res'.")

      plotOntologyFamily(fam = ont.fam.res$families[[fam.name]], data = ont.fam.res$data, plot.type = plot.type, show.ids = show.ids, string.length = string.length, legend.label.size = legend.label.size, legend.position = legend.position, verbose = verbose, n.cores = n.cores)
    },

    #' @description Plot the cell group proportions per sample
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param sample.per.cell Vector indicating sample name with cell names (default: stored vector)
    #' @param cells.to.remove Vector of cell types to remove from the composition
    #' @param cells.to.remain Vector of cell types to remain in the composition
    #' @param notch Whether to show notch in the boxplots
    #' @param alpha Transparency level on the data points (default: 0.2)
    #' @param palette color palette to use for conditions (default: stored $sample.groups.palette)
    #' @param show.significance whether to show statistical significance betwwen sample groups. wilcox.test was used; (* < 0.05; ** < 0.01; *** < 0.001)
    #' @return A ggplot2 object
    plotProportions=function(legend.position = "right", cell.groups = self$cell.groups, sample.per.cell = self$sample.per.cell,
                             sample.groups = self$sample.groups, cells.to.remove = NULL, cells.to.remain = NULL, notch = FALSE,
                             alpha=0.2, palette=self$sample.groups.palette, show.significance = FALSE) {
      if(is.null(cells.to.remove) && is.null(cells.to.remain)){
        plotProportions(legend.position = legend.position, cell.groups = cell.groups, sample.per.cell = sample.per.cell, sample.groups = sample.groups,
                        notch=notch, alpha = alpha, palette=palette, show.significance=show.significance)
      } else {  # Anna modified
        plotProportionsSubset(legend.position = legend.position,
                        cell.groups = cell.groups,
                        sample.per.cell = sample.per.cell,
                        sample.groups = sample.groups,
                        cells.to.remove = cells.to.remove,
                        cells.to.remain = cells.to.remain,
                        notch=notch,
                        alpha = alpha, palette=palette,
                        show.significance=show.significance)
      }
    },

    #' @description Plot the cell numbers per sample
    #' @param legend.position Position of legend in plot. See ggplot2::theme (default="right")
    #' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
    #' @param sample.per.cell Vector indicating sample name with cell names (default: stored vector)
    #' @return A ggplot2 object
    plotCellNumbers=function(legend.position = "right", cell.groups = self$cell.groups, sample.per.cell = self$sample.per.cell, sample.groups = self$sample.groups) {
      df.melt <- data.frame(anno=cell.groups, group=sample.per.cell[match(names(cell.groups), names(sample.per.cell))]) %>%
        table %>%
        rbind %>%
        t %>%
        as.data.frame %>%
        dplyr::mutate(group = sample.groups[match(levels(sample.per.cell), names(sample.groups))]) %>%
        reshape2::melt(., id.vars="group")

      ggplot(df.melt, aes(x=variable, y=value, by=group)) +
        geom_boxplot(position=position_dodge(), outlier.shape = NA) +
        ylab("Cells per sample") +
        xlab("") +
        theme_bw() +
        theme_legend_position(legend.position) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
              legend.title=element_blank()) +
        geom_point(position=position_jitterdodge(jitter.width=0.15), aes(col=group), alpha=0.4) +
        scale_y_continuous(expand=c(0, 0), limits=c(0, (max(df.melt$value) + 50)))
    },

    #' @description Plot compositions in CoDA-PCA space
    #' @return A ggplot2 object
    plotPcaSpace=function(cells.to.remove = NULL, font.size=NULL, palette=self$sample.groups.palette) {
      # Construct sample groups and count data
      # ---- The following can be significantly reduced
      d.counts <- data.frame(anno=self$cell.groups,
                             group=self$sample.per.cell[match(names(self$cell.groups), names(self$sample.per.cell))]) %>%
        table  %>% rbind %>% t
      if(!is.null(cells.to.remove)) d.counts = d.counts[,!(colnames(d.counts) %in% cells.to.remove)]

      d.groups = self$sample.groups[rownames(d.counts)] == self$target.level
      names(d.groups) <- rownames(d.counts)
      # ----

      plotPcaSpace(d.counts, d.groups, self$ref.level, self$target.level, font.size, palette=palette)
    },

    #' @description Plot compositions in CoDA-CDA space
    #' @return A ggplot2 object
    plotCdaSpace=function(cells.to.remain = NULL, cells.to.remove = NULL, samples.to.remove = NULL, font.size=NULL) {
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

      plotCdaSpace(d.counts, d.groups, self$ref.level, self$target.level, font.size)

    },

    #' @description Plot contrast tree
    #' @return A ggplot2 object
    plotContrastTree=function(cells.to.remain = NULL, cells.to.remove = NULL) {
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

      plotContrastTree(d.counts, d.groups, self$ref.level, self$target.level)
    },

    #' @description Plot Loadings
    #' @return A ggplot2 object
    estimateCellLoadings=function(n.cell.counts = 1000, n.seed = 239, cells.to.remove = NULL,
                                  cells.to.remain = NULL, samples.to.remove = NULL, n.iter=1000){
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
                                                     n.seed = n.seed, n.iter = n.iter)

      self$test.results$cda$pvals = getCellSignificance(self$test.results$cda$balances)

      return(invisible(self$test.results[['cda']]))
    },

    estimateGaPartiotion=function(cells.to.remain = NULL, cells.to.remove = NULL, samples.to.remove = NULL, ...){
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

      ga.res <- gaPartition(d.counts, d.groups, ...)

      self$test.results[['ga.partition']] <- rownames(t(ga.res[1,ga.res[1,] != 0,drop=FALSE]))

      return(invisible(self$test.results[['ga.partition']]))
    },

    #' @description Plot Loadings
    #' @param palette palette specification for cell types (default: stored $cell.groups.palette)
    #' @return A ggplot2 object
    plotCellLoadings = function(alpha = 0.01, palette=self$cell.groups.palette, font.size=NULL,
                                ordering='by.pvalue', signif.threshold=0.05, show.pvals=F) {
      possible.ordering = c('by.pvalue', 'by.mean', 'by.median')
      if(!(ordering %in% possible.ordering)){
        warning('Defaulf ordiring \'by.pvalue\' is applied ' )
        ordering = 'by.pvalue'
      }

      cda <- private$getResults('cda', 'estimateCellLoadings()')
      p <- plotCellLoadings(cda, ordering, signif.threshold, font.size, alpha, palette, show.pvals,
                            self$ref.level, self$target.level)

      return(p)
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

      cda <- resampleContrast(d.counts, d.groups,
                             n.cell.counts = n.cell.counts,
                             n.seed = n.seed)
      plotCellLoadings(cda$balances, aplha = aplha)
    },


    ### Segmentation-free cell density

    #' @description Estimate cell density in giving embedding
    #' @param emb cell embedding matrix
    #' @param bins number of bins for density estimation, default 400
    #' @param method density estimation method, graph: graph smooth based density estimation. kde: embedding grid based density  estimation. (default: embGrid)
    #' @param m numeric Maximum order of Chebyshev coeff to compute (default=50) for graph based cell density
    #' @param name slot in which to save the results (default: 'cell.density')
    estimateCellDensity = function(embedding=self$embedding, bins = 400, method = 'kde',
                                   verbose=self$verbose, m=50, n.cores=self$n.cores,name='cell.density'){
      if(is.null(embedding))
        stop("'embedding' must be provided either during the object initialization or during this function call")
      
      sample.per.cell <- self$sample.per.cell
      sample.groups <- self$sample.groups
      
      if (method == 'kde'){
        self$test.results[[name]] <- embedding %>%
          estimateCellDensityKde(sample.per.cell=sample.per.cell, sample.groups=sample.groups, bins=bins)
        return(invisible(self$test.results[[name]]))
      }
      if (method == 'graph'){
        self$test.results[[name]] <-  extractCellGraph(self$data.object) %>%
          estimateCellDensityGraph(sample.per.cell=sample.per.cell, sample.groups=sample.groups,
                                   n.cores=n.cores, m=m, verbose=verbose)
        return(invisible(self$test.results[[name]]))
      }
      
      stop("Unknown method: ", method)
    },

    
    
    #' @description Plot cell density
    #' @param method density estimation method (graph, ked)
    #' @param add.points default is TRUE, add points to cell density figure
    #' @param contours specify cell types for contour, multiple cell types are also supported
    #' @param contour.color color for contour line
    #' @param contour.conf confidence interval of contour
    #' @param name slot in which to saved results from estimateCellDensity (default: 'cell.density')
    #' @return A ggplot2 object
    plotCellDensity = function(method ='kde', show.legend=FALSE, legend.position=NULL, show.grid=TRUE, add.points=TRUE,size=0.1,
                               point.col='#FCFDBFFF', contours=NULL, contour.color='white', contour.conf='10%', name='cell.density') {
      dens.res <- private$getResults(name)
      
      if (method == 'kde'){
        if (dens.res$method!='kde') stop('please estimate cell density with estimateCellDensity(method="kde")')
        # calculate sample.per.cell
        condition.per.cell <- as.factor(setNames( as.character(self$sample.groups[ as.character(self$sample.per.cell)]), names(self$sample.per.cell) ))
        
        target.density <- dens.res %$% data.frame(density.emb, z=density.fraction[[self$target.level]])
        ref.density <- dens.res %$% data.frame(density.emb, z=density.fraction[[self$ref.level]])
        
        mi <- min(min(ref.density$z), min(target.density$z))
        ma <- max(max(ref.density$z), max(target.density$z))
        
        p1 <- plotDensity(target.density, bins=dens.res$bins, legend.position=legend.position, show.legend=show.legend, title=self$ref.level, show.grid=show.grid, mi=mi, ma=ma)
        p2 <- plotDensity(ref.density, bins=dens.res$bins, legend.position=legend.position, show.legend=show.legend, title=self$target.level, show.grid=show.grid, mi=mi, ma=ma)
        
        if (add.points){
          emb <- self$embedding %>% as.data.frame()
          colnames(emb) <- c('x','y')
          emb$z <- 1
          nname1 <- names(condition.per.cell)[condition.per.cell == self$ref.level] %>%
            sample(min(2000, length(.)))
          
          nname2 <- names(condition.per.cell)[condition.per.cell == self$target.level] %>%
            sample(min(2000, length(.)))
          
          p1 <- p1 + geom_point(data=emb[nname1, ], aes(x=x, y=y), col=point.col, size=0.00001, alpha=0.2)
          p2 <- p2 + geom_point(data=emb[nname2, ], aes(x=x, y=y), col=point.col, size=0.00001, alpha=0.2)
        }
      }  
      else if (method =='graph'){
        if (dens.res$method != 'graph') stop('please estimate cell density with estimateCellDensity(method="graph")')
        emb <- self$embedding
        target.density <- dens.res %$% density.fraction[[self$target.level]]
        ref.density <- dens.res %$% density.fraction[[self$ref.level]]        
        p1 <- sccore::embeddingPlot(emb, plot.theme=ggplot2::theme_bw(), colors = target.density, size = size,title = self$target.level, legend.position = legend.position, show.legend = show.legend) + 
          #scale_fill_gradient2(low = col[1], high = col[3], mid = col[2], midpoint = 0, limits = c(mi, ma)) +
          theme(legend.background = element_blank())
        p2 <- sccore::embeddingPlot(emb, plot.theme=ggplot2::theme_bw(), colors = ref.density, size = size,, title=self$ref.level, legend.position = legend.position, show.legend = show.legend) + 
          #scale_fill_gradient2(low = col[1], high = col[3], mid = col[2], midpoint = 0, limits = c(mi, ma)) +
          theme(legend.background = element_blank())
      }else stop("Unknown method: ", method)
      if(!is.null(contours)){
        cnl <- do.call(c, lapply(sn(contours), function(x) getContour(self$embedding, cell.type=self$cell.groups , cell=x ,conf=contour.conf, color=contour.color)))
        p1 <- p1 + cnl
        p2 <- p2 + cnl
      }
      return(list(ref=p1, target=p2))
    },


    #' @description estimate differential cell density
    #' @param method density estimation method (graph or ked)
    #' @param col color palettes,  default is c('blue','white','red')
    #' @param type method to calculate differential cell density; t.test, wilcox or subtract (target subtract ref density);
    #' @param contours specify cell types for contour, multiple cell types are also supported
    #' @param contour.color color for contour line
    #' @param z.cutoff absolute z score cutoff
    #' @param contour.conf confidence interval of contour
    #' @param name slot in which to saved results from estimateCellDensity (default: 'cell.density')
    diffCellDensity=function(method = 'kde', type='subtract', col=c('blue','white','red'), show.legend=FALSE, legend.position=NULL, title=NULL, show.grid=NULL, plot=TRUE, contours=NULL, contour.color='white', contour.conf='10%', z.cutoff=NULL, size =0.2, adjust.pvalues=TRUE, name = 'cell.density', ...){
      # TODO: rename it to start with estimate*
      dens.res <- private$getResults(name)
      if (method == 'graph'){
        if (dens.res$method != 'graph') stop('please estimate cell density with estimateCellDensity(method="graph")')
        density.emb <-  self$embedding
        density.mat <-  dens.res$density.mat
        cname <- intersect(rownames(density.emb),rownames(density.mat))
        density.emb <- density.emb[cname,]
        density.mat <- density.mat[cname,]
      }
      else if (method == 'kde'){
        if (dens.res$method != 'kde') stop('please estimate cell density with estimateCellDensity(method="kde")')
        density.emb <- dens.res$density.emb
        density.mat <- dens.res$density.mat
        # remove empty bins
        index = density.emb$counts > 0
        density.emb <- density.emb[index,1:2]
        density.mat <- density.mat[index,]
      }else stop("Unknown method: ", method)
      
      mat <- diffCellDensity(density.emb, density.mat, self$sample.groups, bins=bins, target.level=self$target.level,ref.level=self$ref.level, type=type, z.cutoff=z.cutoff, adjust.pvalues=adjust.pvalues)
      emb <- mat[,1:2]
      score <- mat$z
      names(score) <- rownames(mat)
      
      if (plot){
        fig <- sccore::embeddingPlot(emb, plot.theme=ggplot2::theme_bw(), colors = score, size=size,title = title, legend.position = legend.position, show.legend = show.legend, ...) + scale_color_gradient2(low = col[1], high = col[3], mid = col[2],, midpoint = 0) +theme(legend.background = element_blank()) +  labs(color='Zscore') 
        
        if(!is.null(contours)){
          cnl <- do.call(c, lapply(sn(contours), function(x)
            getContour(self$embedding, cell.type=self$cell.groups , cell=x, conf=contour.conf, color=contour.color)))
          fig <- fig + cnl
        }
        return(fig)
      }
      return(mat)
    },

    #' @title Plot inter-sample expression distance
    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
    #' @param name Test results to plot (default=expression.shifts)
    #' @param notch Show notches in plot, see ggplot2::geom_boxplot for more info (default=T)
    #' @param cell.groups Named factor with cell names defining groups/clusters (default: stored $cell.groups vector)
    #' @param weighted.distance whether to weigh the expression distance by the sizes of cell types (default: TRUE), or show distances for each individual cell type
    #' @param show.significance whether to show statistical significance between sample groups. wilcox.test was used; (\* < 0.05; \*\* < 0.01; \*\*\* < 0.001)
    #' @param alpha dot transparency
    #' @return A ggplot2 object
    plotExpressionDistance = function(name='expression.shifts', notch=TRUE, sample.groups=self$sample.groups, weighted.distance=TRUE,
                                      min.cells=10, palette=self$sample.groups.palette, show.significance=FALSE, alpha=0.2) {
      cluster.shifts <- private$getResults(name, 'estimateExpressionShiftMagnitudes()')
      ctdml <- cluster.shifts$ctdml
      valid.comparisons <- cluster.shifts$valid.comparisons
      if (!weighted.distance) {
        gg <- plotExpressionDistanceIndividual(ctdml, valid.comparisons, sample.groups=sample.groups, notch=notch, alpha=alpha, min.cells=min.cells, show.significance=show.significance)
      } else {
        gg <- plotExpressionDistanceJoint(ctdml, valid.comparisons, sample.groups=sample.groups, notch=notch, alpha=alpha, show.significance=show.significance)
      }

      if(!is.null(palette)) {
        gg <- gg + scale_fill_manual(values=palette)
      }

      if(show.significance) {
        gg <- gg + ggpubr::stat_compare_means(aes(group = group), label = "p.signif", label.x.npc="centre")  # willcox test
      }

      return(gg)
    },

    #' @title Plot sample-sample expression distance as a 2D embedding
    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
    #' @param cell.type If a name of a cell type is specified, the sample distances will be assessed based on this cell type alone. Otherwise (cell.type=NULL, default), sample distances will be estimated as an average distance across all cell types (weighted by the minimum number of cells of that cell type between any two samples being compared)
    #' @param method dimension reduction methods (MDS or tSNE ) , default is MDS
    #' @param perplexity tSNE perpexity (default: 4)
    #' @param max_iter tSNE max_iter (default: 1e3)
    #' @param palette a set of colors to use for conditions (default: stored $sample.groups.palette)
    #' @return A ggplot2 object
    plotExpressionDistanceEmbedding = function(name='expression.shifts', sample.groups = self$sample.groups, cell.type = NULL, method = 'tSNE', perplexity=4, max_iter=1e3, palette=self$sample.groups.palette) {
      cluster.shifts <- private$getResults(name, 'estimateExpressionShiftMagnitudes()')
      plotExpressionDistancetSNE(cluster.shifts, sample.groups = sample.groups, cell.type = cell.type, method = method, perplexity=perplexity, max_iter=max_iter, palette=palette)
    },

    #' @description Extract contour from embedding
    #' @param cell specify cell types for contour, mutiple cell types are also suported
    #' @param conf confidence interval of contour
    getContour = function(cells,  color = 'white', linetype = 2, conf = "10%") {
      cnl <- do.call(c, lapply(sn(cells), function(x) getContour(self$embedding, cell.type=self$cell.groups, linetype = linetype,
                                                                 cell=x ,conf = conf, color = color)))
      return(cnl)
    },


    ### Cluster-free differential expression

    #' @description Estimate differential expression Z-scores between two conditions per individual cell
    #' @param max.z z-score value to winsorize the estimates for reducing impact of outliers. Default: 20.
    #' @param min.expr.frac minimal fraction of cell expressing a gene for estimating z-scores for it. Default: 0.001.
    #' @param normalize whether to normalize z-scores over std for reference ("ref") or both ("both"). Default: "ref".
    #' @return Sparse matrix of z-scores with genes as columns and cells as rows.
    #' Cells that have only one condition in their expression neighborhood have NA Z-scores for all genes.
    #' Results are also stored in the `cluster.free.z` field.
    estimateClusterFreeZScores = function(n.top.genes=NULL, max.z=20, min.expr.frac=0.001, normalize=c("ref", "both"),
                                          verbose=self$verbose, n.cores=self$n.cores) {
      normalize <- match.arg(normalize)
      cm <- extractJointCountMatrix(self$data.object)
      genes <- private$getTopGenes(n.top.genes, gene.selection="expression", cm.joint=cm, min.expr.frac=min.expr.frac)

      if (verbose)
        message("Estimating cluster-free Z-scores for ", length(genes), " top genes")

      is.ref <- self$sample.per.cell %>%
        {setNames(self$sample.groups[as.character(.)] == self$ref.level, names(.))}

      adj.mat <- extractCellGraph(self$data.object) %>% igraph::as_adj()
      cell.names <- intersect(rownames(cm), rownames(adj.mat)) %>%
        intersect(names(is.ref))

      z.mat <- clusterFreeZScoreMat(adj.mat[cell.names, cell.names], cm[cell.names, genes, drop=FALSE],
                                    is.ref[cell.names], normalize_both=(normalize == "both"), verbose=verbose, n_cores=n.cores)
      z.mat@x %<>% pmin(max.z) %>% pmax(-max.z)

      self$test.results[["cluster.free.z"]] <- z.mat
      return(invisible(self$test.results[["cluster.free.z"]]))
    },

    getMostChangedGenes = function(n, min.z=0.5, max.z=20, cell.subset=NULL, excluded.genes=NULL, included.genes=NULL) {
      z.scores <- private$getResults("cluster.free.z", "estimateClusterFreeZScores")
      if (!is.null(cell.subset)) {
        z.scores <- z.scores[cell.subset,]
      }

      if (is.null(included.genes)) {
        included.genes <- colnames(z.scores)
      }

      z.scores@x %<>% abs() %>% pmin(max.z)
      z.scores@x[z.scores@x < min.z] <- 0
      scores <- colMeans(z.scores, na.rm=TRUE) %>% sort(decreasing=TRUE) %>%
        .[setdiff(names(.), excluded.genes)] %>% .[intersect(names(.), included.genes)] %>%
        .[1:min(n, length(.))]
      return(scores)
    },

    #' @description Estimate Cluster-free Expression Shift
    #' @return Vector of cluster-free expression shifts per cell. Values above 1 correspond to difference between conditions.
    #' Results are also stored in the `cluster.free.expr.shifts` field.
    estimateClusterFreeExpressionShifts = function(n.top.genes=3000, gene.selection="change", excluded.genes=NULL,
                                                   verbose=self$verbose, n.cores=self$n.cores) {
      cm <- extractJointCountMatrix(self$data.object)
      genes <- private$getTopGenes(n.top.genes, gene.selection=gene.selection, cm.joint=cm, excluded.genes=excluded.genes)
      cm <- Matrix::t(cm[, genes])

      is.ref <- (self$sample.groups[levels(self$sample.per.cell)] == self$ref.level)

      nns.per.cell <- extractCellGraph(self$data.object) %>%
        igraph::as_adjacency_matrix() %>% as("dgTMatrix") %>%
        {setNames(split(.@j, .@i + 1), rownames(.))}

      shifts <- estimateClusterFreeExpressionShifts(cm, self$sample.per.cell[names(nns.per.cell)], nns.per.cell,
                                                    is.ref, verbose=verbose, n_cores=n.cores)
      self$test.results[["cluster.free.expr.shifts"]] <- shifts

      return(invisible(shifts))
    },

    #' @description Performs graph smoothing of the cluster-free DE Z-scores
    #' @param smoothing `beta` parameter of the \link[sccore:heatFilter]{heatFilter}. Default: 20.
    #' @param filter graph filter function. Default: \link[sccore:heatFilter]{heatFilter}.
    #' @param ... parameters forwarded to \link[sccore:smoothSignalOnGraph]{smoothSignalOnGraph}
    #' @return Sparse matrix of smoothed Z-scores. Results are also stored in the `cluster.free.z.smoothed` field.
    smoothClusterFreeZScores = function(n.top.genes=1000, smoothing=20, filter=NULL, gene.selection="change", excluded.genes=NULL,
                                        n.cores=self$n.cores, verbose=self$verbose, ...) {
      z.scores <- private$getResults("cluster.free.z", "estimateClusterFreeZScores")
      genes <- private$getTopGenes(n.top.genes, gene.selection=gene.selection,
                                   excluded.genes=excluded.genes, included.genes=colnames(z.scores))
      z.scores <- z.scores[,genes]

      if (verbose) message("Smoothing Z-scores for ", ncol(z.scores), " genes passed filtration")
      if (is.null(filter)) {
        filter <- function(...) heatFilter(..., beta=smoothing)
      }

      z.scores@x[is.na(z.scores@x)] <- 0
      z.smoothed <- z.scores %>%
        smoothSignalOnGraph(extractCellGraph(self$data.object), filter, n.cores=n.cores,
                            progress=verbose, ...)

      z.smoothed[is.na(z.scores)] <- NA
      self$test.results[["cluster.free.z.smoothed"]] <- z.smoothed
      return(invisible(z.smoothed))
    },

    #' @description Estimate Gene Programmes based on cluster-free Z-scores on a subsample of
    #' cells using \link[fabia:fabia]{fabia}.
    #' @param n.programmes maximal number of gene programmes to find (parameter `p` for fabia). Default: 15.
    #' @param ... keyword arguments forwarded to \link{estimateGeneProgrammes}
    #' @return a list includes:
    #'   - `fabia`: \link[fabia:Factorization]{fabia::Factorization} object, result of the
    #'       \link[fabia:fabia]{fabia::fabia} call
    #'   - `sample.ids`: ids of the subsampled cells used for fabia estimates
    #'   - `scores.exact`: vector of fabia estimates of gene programme scores per cell. Estimated only for the
    #'     subsampled cells.
    #'   - `scores.approx`: vector of approximate gene programme scores, estimated for all cells in the dataset
    #'   - `loadings`: matrix with fabia gene loadings per programme
    #'   - `gene.scores`: list of vectors of gene scores per programme. Contains only genes, selected for
    #'     the programme usin fabia biclustering.
    #'   - `bi.clusts` fabia biclustering information, result of the \link[fabia:extractBic]{fabia::extractBic} call
    estimateGeneProgrammes = function(n.top.genes=NULL, n.programmes=15, gene.selection="change", name="gene.programmes", ...) {
      if (!requireNamespace("fabia", quietly=TRUE))
        stop("fabia package must be installed to run this function")

      z.scores <- private$getResults("cluster.free.z.smoothed", "smoothClusterFreeZScores")
      if (!is.null(n.top.genes)) {
        genes <- private$getTopGenes(n.top.genes, gene.selection=gene.selection, included.genes=colnames(z.scores))
        z.scores <- z.scores[,genes]
      }

      res <- estimateGeneProgrammesFabia(z.scores, n.programmes, ...)
      fr <- res$fabia

      bi.clusts <- fabia::extractBic(fr)
      mask <- bi.clusts$bic[,"bixv"] %>% {sapply(., length) > 0}
      res$scores.exact <- t(fr@Z[mask,])
      res$scores.approx <- t(t(fr@L[,mask]) %*% (Matrix::t((z.scores[,rownames(fr@L)]) - fr@center) / fr@scaleData))
      res$loadings <- fr@L[,mask]
      res$gene.scores <- apply(bi.clusts$bic, 1, `[[`, "bixv")[mask]
      res$bi.clusts <- bi.clusts

      colnames(res$scores.exact) <- colnames(res$scores.approx) <-
        colnames(res$loadings) <- names(res$gene.scores) <- paste0("P", 1:ncol(res$scores.exact))
      self$test.results[[name]] <- res

      return(invisible(self$test.results[[name]]))
    },

    plotGeneProgrammeScores = function(name="gene.programmes", approximate=FALSE, build.panel=TRUE, nrow=NULL,
                                       gradient.range.quantile=0.975, plot.na=approximate, legend.title="Score", ...) {
      fr <- private$getResults(name, "estimateGeneProgrammes")
      scores <- if (approximate) fr$scores.approx else fr$scores.exact

      ggs <- lapply(1:ncol(scores), function(i) {
        self$plotEmbedding(colors=scores[,i], gradient.range.quantile=gradient.range.quantile, plot.na=plot.na,
                           legend.title=legend.title, title=colnames(scores)[i], ...)
      })
      ggs <- lapply(ggs,function(x) x+theme(legend.background = element_blank()))
      if (build.panel)
        return(cowplot::plot_grid(plotlist=ggs, nrow=nrow))

      return(ggs)
    },

    plotGeneProgrammeGenes = function(programme.id, name="gene.programmes", max.genes=9, plot.expression=FALSE, ...) {
      fr <- private$getResults(name, "estimateGeneProgrammes")
      scores <- fr$gene.scores[[programme.id]]
      if (is.null(scores))
        stop("Can't find programme", programme.id)

      scores %<>% .[1:min(length(.), max.genes)]
      return(self$plotGeneExpressionComparison(scores=scores, plot.expression=plot.expression, ...))
    },

    #' @description Plot cluster-free expression shift z-scores
    #' @param cell.groups cell type labels. Set to NULL if it shouldn't be shown
    #' @param plot.both.conditions show both case and control cells. Normally, showing control cells doesn't
    #' make sense, as control cells always have small distance from control.
    #' @param max.shift all shift values above `max.shift` are set to this value when plotting. Default: 95% of the shifts.
    #' @param font.size size range for cell type labels
    #' @param ... parameters forwarded to \link[sccore:embeddingPlot]{embeddingPlot}
    plotClusterFreeExpressionShifts = function(cell.groups=self$cell.groups, plot.both.conditions=FALSE, plot.na=FALSE, max.shift=NULL,
                                               alpha=0.2, font.size=c(3,5), ...) {
      shifts <- private$getResults("cluster.free.expr.shifts", "estimateClusterFreeExpressionShifts")
      if (is.null(self$embedding))
        stop("embedding must not be NULL. Please, set the 'embedding' field.")

      if (!plot.both.conditions) {
        shifts %<>%  .[self$sample.groups[self$sample.per.cell[names(.)]] != self$ref.level]
      }

      if (is.null(max.shift)) {
        max.shift <- quantile(shifts, 0.95, na.rm=TRUE)
      }

      if (max.shift > 0) {
        shifts %<>% pmin(max.shift)
      }

      gg <- self$plotEmbedding(colors=shifts, plot.na=plot.na, alpha=alpha, ...)
      gg$scales$scales %<>% .[sapply(., function(s) s$aesthetics != "colour")]

      if (!is.null(cell.groups)) {
        ann.ls <- self$plotEmbedding(groups=cell.groups)$layers
        gg <- gg + ann.ls[[which(sapply(ann.ls, function(l) "GeomLabelRepel" %in% class(l$geom)))]]
      }

      gg <- gg +
        scale_size_continuous(range=font.size, trans='identity', guide='none') +
        scale_color_gradientn(colours=c('#2c7bb6', '#abd9e9', '#ffffbf', '#fdae61', '#d7191c'),
                              values=c(0, 0.5 / max.shift, 1 / max.shift, 0.5 - 0.5 * max.shift, 1.0), # Ensure that 1.0 has yellow color
                              name="Ratio") +theme(legend.background = element_blank())
      return(gg)
    },

    plotMostChangedGenes = function(n.top.genes, min.z=0.5, max.z=20, max.z.plot=max.z, cell.subset=NULL, excluded.genes=NULL, ...) {
      scores <- self$getMostChangedGenes(n.top.genes, min.z=min.z, max.z=max.z,
                                         cell.subset=cell.subset, excluded.genes=excluded.genes)
      self$plotGeneExpressionComparison(scores=scores, max.z=max.z.plot, ...)
    },

    plotGeneExpressionComparison = function(genes=NULL, scores=NULL, max.expr=NULL, plot.z=TRUE, plot.expression=TRUE, max.z=5, smoothed=FALSE, plot.na=-1, ...) {
      if (is.null(genes)) {
        if (is.null(scores)) stop("Either 'genes' or 'scores' must be provided")
        genes <- names(scores)
      }

      if (!plot.z && !plot.expression) return()

      if (plot.z) {
        if (smoothed) {
          z.scores <- self$test.results$cluster.free.z.smoothed
          if ((is.null(z.scores) || !all(genes %in% colnames(z.scores)))){
            missed.genes <- setdiff(colnames(z.scores), genes)
            warning("Smoothed Z-scores for genes ", paste(missed.genes, collapse=', '), " are not estimated. See smoothClusterFreeZScores().")
            smoothed <- FALSE
          }
        }

        if (!smoothed) {
          z.scores <- self$test.results$cluster.free.z
          if ((is.null(z.scores) || !all(genes %in% colnames(z.scores)))) {
            missed.genes <- setdiff(colnames(z.scores), genes)
            warning("Z-scores for genes ", paste(missed.genes, collapse=', '), " are not estimated. See estimateClusterFreeZScores().")
            plot.z <- FALSE
          }
        }
      }

      condition.per.cell <- self$sample.per.cell %>%
        {setNames(as.character(self$sample.groups[as.character(.)]), names(.))} %>%
        as.factor()

      ggs <- lapply(genes, function(g) {
        lst <- list()
        if (plot.expression) {
          expr <- extractGeneExpression(self$data.object, g)
          m.expr <- if (is.null(max.expr)) max(expr) else max.expr
          lst <- lapply(unique(condition.per.cell), function(sg) {
            self$plotEmbedding(colors=expr, title=paste(sg, " ",g), groups=condition.per.cell, subgroups=sg,
                               color.range=c(0, m.expr), legend.title="Expression", plot.na=FALSE, ...)
          }) %>% c(lst)
        }

        if (plot.z) {
          title <- if (is.null(scores)) g else paste0(g, ": ", signif(scores[g], 3))
          lst <- self$plotEmbedding(colors=z.scores[,g], title=title, color.range=c(-max.z, max.z), plot.na=plot.na, ...) %>%
            list() %>% c(lst)
        }
        lst <- lapply(lst,function(x) x+theme(legend.background = element_blank()))
        if (length(lst) > 1) cowplot::plot_grid(plotlist=lst, ncol=3) else lst[[1]]
      })

      if (length(genes) == 1) return(ggs[[1]])
      return(ggs)
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

    getTopGenes = function(n, gene.selection=c("change", "expression", "od"), cm.joint=NULL,
                           min.expr.frac=0.0, excluded.genes=NULL, included.genes=NULL, ...) {
      gene.selection <- match.arg(gene.selection)
      if ((gene.selection == "change") && is.null(self$test.results$cluster.free.z)) {
        warning("Please run estimateClusterFreeZScores() first to use gene.selection='change'. Fall back to gene.selection='expression'.")
        gene.selection <- "expression"
      }

      if (gene.selection == "change") {
        genes <- names(self$getMostChangedGenes(Inf, ...))
      } else if (gene.selection == "od") {
        genes <- extractOdGenes(self$data.object, n + length(excluded.genes))
      } else {
        if (is.null(cm.joint)) {
          cm.joint <- extractJointCountMatrix(self$data.object)
        }

        cm.joint@x <- 1 * (cm.joint@x > 0)
        col.means <- colMeans(cm.joint, na.rm=TRUE)
        gene.mask <- (col.means >= min.expr.frac)
        genes <- col.means %>% order(decreasing=TRUE) %>%
          .[gene.mask[.]] %>% colnames(cm.joint)[.]
      }

      if (is.null(included.genes)) {
        included.genes <- genes
      }
      genes %<>% setdiff(excluded.genes) %>% intersect(included.genes) %>% .[1:min(length(.), n)]
      return(genes)
    },

    getOntologyPvalueResults=function(genes, type, p.adj=0.05, min.genes=1, subtype=NULL, cell.subgroups=NULL) {
      if(!type %in% c("GO", "DO", "GSEA"))
        stop("'type' must be 'GO', 'DO', or 'GSEA'.")

      if(!is.null(subtype) && !all(subtype %in% c("BP", "CC", "MF")))
        stop("'subtype' must be 'BP', 'CC', or 'MF'.")

      if((length(genes) != 1) || (!genes %in% c("down","up","all")))
        stop("'genes' must be 'down', 'up', or 'all'.")

      ont.res <- self$test.results[[type]][["res"]]
      if(is.null(ont.res)) stop(paste0("No results found for '", type, "'. Please run 'estimateOntology' first."))

      ont.res %<>% preparePlotData(type, p.adj, min.genes)

      # Extract genes and subgroups
      if(type == "GSEA") {
        ont.res %<>% addGseaGroup() %>% rename(geneID=core_enrichment)
      } else {
        ont.res %<>% .[[genes]]
      }

      if(is.null(ont.res))
        stop(paste0("No results found for ", genes, " genes for ", type, ".")) # TODO: Include GSEA

      if(!is.null(cell.subgroups)) {
        if(!cell.subgroups %in% unique(ont.res$Group))
          stop("'cell.subgroups' not found in results.")

        ont.res %<>% .[cell.subgroups]
      }

      if(!is.null(subtype)) {
        ont.res %<>% filter(Type %in% subtype)
      }

      return(ont.res)
    }
  )
)
