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

    initialize=function(data.object, sample.groups=NULL, cell.groups=NULL, sample.per.cell=NULL, ref.level=NULL, target.level=NULL, sample.groups.palette=NULL, cell.groups.palette=NULL, embedding=NULL, n.cores=1, verbose=TRUE) {
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

        return()
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
      
      if(is.null(sample.groups.palette)) {
        self$sample.groups.palette <- setNames(rev(scales::hue_pal()(length(levels(sample.groups)))), levels(sample.groups))
      } else {
        self$sample.groups.palette <- sample.groups.palette
      }
      
      if(is.null(cell.groups.palette)) {
        self$cell.groups.palette <- setNames( rainbow(length(levels(cell.groups)),s=0.9,v=0.9), levels(cell.groups))
      } else {
        self$cell.groups.palette <- cell.groups.palette
      }
      
      if(is.null(embedding)) {
        # TODO: extract from the object
      } else {
        self$embedding <- embedding;
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
    
    estimateCommonExpressionShiftMagnitudes=function(sample.groups=self$sample.groups, cell.groups=self$cell.groups, n.cells=NULL, n.randomizations=50, n.subsamples=30, min.cells=10, n.cores=1, verbose=TRUE,  mean.trim=0.1, name='common.expression.shifts') {
      if (is.null(sample.groups))
        stop("'sample.groups' must be provided either during the object initialization or during this function call")
      
      if (is.null(cell.groups)) {
        stop("'cell.groups' must be provided either during the object initialization or during this function call")
      }
      
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
    #' @param name Test results to plot (default=expression.shifts)
    #' @param size.norm Plot size normalized results. Requires cell.groups, and sample.per.cell (default=F)
    #' @param notch Show notches in plot, see ggplot2::geom_boxplot for more info (default=T)
    #' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
    #' @param sample.per.cell Named sample factor with cell names (default: stored vector)
    #' @return A ggplot2 object
    plotExpressionShiftMagnitudes=function(name="expression.shifts", size.norm=F, notch = T, cell.groups=self$cell.groups, sample.per.cell=self$sample.per.cell, palette=self$cell.groups.palette) {
      cluster.shifts <- private$getResults(name)$df

      if (is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if (is.null(sample.per.cell)) stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      plotExpressionShiftMagnitudes(cluster.shifts = cluster.shifts, size.norm = size.norm, notch = notch, cell.groups = cell.groups, sample.per.cell = sample.per.cell, palette=palette)
    },
    
    
    plotCommonExpressionShiftMagnitudes=function(name='common.expression.shifts', show.subsampling.variability=FALSE, show.jitter=FALSE, jitter.alpha=0.05, palette=self$cell.groups.palette) { 
      res <- private$getResults(name)
      cn <- setNames(names(res[[1]]),names(res[[1]]))
      if(show.subsampling.variability) { # average across patient pairs
        if(length(res)<2) stop('the result has only one subsample; please set show.sampling.variability=FALSE')
        df <- do.call(rbind,lapply(res,function(d) data.frame(val=unlist(lapply(d,mean)),cell=names(d))))
      } else { # average across subsampling rounds
        df <- do.call(rbind,lapply(cn,function(n) data.frame(val=colMeans(do.call(rbind,lapply(res,function(x) x[[n]]))),cell=n)))
      }
      odf <- na.omit(df);
      df <- data.frame(cell=levels(df$cell),mean=tapply(df$val,df$cell,mean),se=tapply(df$val,df$cell,function(x) sd(x)/sqrt(length(x))),stringsAsFactors=F)
      df <- df[order(df$mean,decreasing=F),]
      df$cell <- factor(df$cell,levels=df$cell)
      df <- na.omit(df);
      
      p <- ggplot(df,aes(x=cell,y=mean,fill=cell))+
        geom_bar(stat='identity') + 
        geom_errorbar(aes(ymin=mean-se*1.96, ymax=mean+se*1.96),width=0.2)+
        geom_hline(yintercept = 1,linetype=2,color='gray50')+
        theme_bw() +
        theme(axis.text.x=element_text(angle = 90, hjust=1, size=12), axis.text.y=element_text(angle=90, hjust=0.5, size=12))+ guides(fill=FALSE)+
        theme(legend.position = "none")+
        labs(x="", y="normalized distance (common)")
      if(show.jitter) p <- p+geom_jitter(data=odf,aes(x=cell,y=val),position=position_jitter(0.1),show.legend=FALSE,alpha=jitter.alpha);
      if(!is.null(palette)) {
        p <- p+ scale_fill_manual(values=palette)
      }
      p
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
    
    #' @description  Plot embedding
    #' @param embedding A cell embedding to use (two-column data frame with rownames corresponding to cells) (default: stored embedding object)
    #' @param plot.theme plot theme to use (default: ggplot2::theme_bw())
    #' @param ... other parameters are passed to sccore::embeddingPlot()
    plotEmbedding=function( embedding=self$embedding, plot.theme=ggplot2::theme_bw(), ... ) {
      if(is.null(embedding)) stop("embedding must be provided to cacoa constructor or to this method")
      sccore::embeddingPlot(embedding, plot.theme=plot.theme, ...)
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
    #' @param notch Whether to show notch in the boxplots
    #' @param alpha Transparency level on the data points (default: 0.2)
    #' @param palette color palette to use for conditions (default: stored $sample.groups.palette)
    #' @return A ggplot2 object
    plotProportions=function(legend.position = "right",
                             cell.groups = self$cell.groups,
                             sample.per.cell = self$sample.per.cell,
                             sample.groups = self$sample.groups,
                             cells.to.remove = NULL,
                             cells.to.remain = NULL,
                             notch = FALSE,
                             alpha=0.2, palette=self$sample.groups.palette) {
      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")

      if(is.null(sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")

      if(is.null(sample.per.cell)) stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      if(is.null(cells.to.remove) && is.null(cells.to.remain)){
        plotProportions(legend.position = legend.position, cell.groups = cell.groups, sample.per.cell = sample.per.cell, sample.groups = sample.groups, notch=notch, alpha = alpha, palette=palette)
      }else{  # Anna modified
        plotProportionsSubset(legend.position = legend.position,
                        cell.groups = cell.groups,
                        sample.per.cell = sample.per.cell,
                        sample.groups = sample.groups,
                        cells.to.remove = cells.to.remove,
                        cells.to.remain = cells.to.remain,
                        notch=notch,
                        alpha = alpha, palette=palette)
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
    plotPcaSpace=function(cells.to.remove = NULL, palette=self$sample.groups.palette) {
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

      plotPcaSpace(d.counts, d.groups, palette=palette)
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
    #' @param palette palette specification for cell types (default: stored $cell.groups.palette)
    #' @return A ggplot2 object
    plotCellLoadings = function(n.cell.counts = 1000,
                              n.seed = 239,
                              aplha = 0.01,
                              font.size = NULL,
                              cells.to.remove = NULL,
                              samples.to.remove = NULL,
                              palette=self$cell.groups.palette){
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
                       font.size = font.size, palette=palette)
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
    #' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
    #' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
    #' @param target.level target/disease level for sample.group vector
    #' @param bins number of bins for density esitmation, default 400
    #' @param condition.per.cell Named group factor with cell names. Must have exactly two levels.
    #' @param by.sample  if TRUE, density will esitmated by sample and quantiles normlization will applied to indivisual sample. If FALSE, cell condition.per.cell need to be provided and density will simply esitmated by condition.per.cell. 
    #' @add.ponits add.ponits  show cells in density plot     
    estimateCellDensity = function(embedding=self$embedding, cell.groups=self$cell.groups, sample.groups=self$sample.groups, sample.per.cell = self$sample.per.cell, ref.level=self$ref.level, target.level=self$target.level, bins = 400, by.sample = TRUE){

      if(is.null(embedding)) stop("'embedding' must be provided either during the object initialization or during this function call")

      if(is.null(cell.groups)) stop("'cell.groups' must be provided either during the object initialization or during this function call")
      
      if(is.null(sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")
      
      if(is.null(sample.per.cell)) stop("'sample.per.cell' must be provided either during the object initialization or during this function call")

      if(is.null(ref.level)) stop("'ref.level' must be provided either during the object initialization or during this function call")

      if(is.null(target.level)) stop("'target.level' must be provided either during the object initialization or during this function call")
      
      # calculate sample.per.cell
      condition.per.cell <- as.factor(setNames( as.character(sample.groups[ as.character(sample.per.cell)]), names(sample.per.cell) ))

      ref.level <- self$ref.level
      target.level <- self$target.level
      sample.groups <- self$sample.groups
      
      self$test.results[['bins']] <- bins
      res <- estimateCellDensity(embedding, sample.per.cell = sample.per.cell, sample.groups = sample.groups, bins = bins, ref.level =
                                   ref.level, target.level = target.level, condition.per.cell = condition.per.cell, by.sample =
                                   by.sample)
      self$test.results[['density.mat']] <- res[['density.mat']]
      self$test.results[['target.density']] <- res[['density.fraction']][[target.level]]
      self$test.results[['ref.density']] <- res[['density.fraction']][[ref.level]]
      #self$test.results[['diff.density']] <- res[['density.fraction']][[target.level]] - res[['density.fraction']][[ref.level]]

      # adjust embedding to the same density space 
      x <- embedding[, 1]
      y <- embedding[, 2]
      x <- (x - range(x)[1])
      x <- (x / max(x)) * bins
      
      y <- (y - range(y)[1])
      y <- (y / max(y)) * bins
      emb2 <- data.frame(x = x, y = y)
      self$test.results[['density.emb']] <- emb2
      
      
      #count cell number in each bin
      x=emb2[,1]
      y=emb2[,2]
      s1 = seq(from = min(x),
               to = max(x),
               length.out = bins + 1)
      s2 = seq(from = min(y),
               to = max(y),
               length.out = bins + 1)
      dcounts = table(cut(x, breaks = s1), cut(y, breaks = s2)) %>% as.matrix.data.frame
      self$test.results[['density.counts']] <- dcounts
      
      return(invisible(self$density[['bins']]))
    },

    ##' @description Plot cell density
    ##' @param add.ponits default is TRUE, add points to cell density figure
    ##' @param condition.per.cell Named group factor with cell names. Must have exactly two levels. condition.per.cell must be provided when add.ponits is TRUE
    plotCellDensity = function(legend = NULL, title = NULL, grid = NULL, add.ponits = TRUE, condition.per.cell = NULL, color='B', point.col='#FCFDBFFF') {
      bins <- private$getResults('bins', 'estimateCellDensity()')
      ref <- self$ref.level
      target <- self$target.level
      
      # calculate sample.per.cell
      condition.per.cell <- as.factor(setNames( as.character(self$sample.groups[ as.character(self$sample.per.cell)]), names(self$sample.per.cell) ))
      
      
      target.density <- private$getResults('target.density', 'estimateCellDensity()')
      ref.density <- private$getResults('ref.density', 'estimateCellDensity()')
      
      mi <- min(c(min(ref.density), min(target.density)))
      ma <- max(c(max(ref.density), max(target.density)))
      
      p1 <- plotDensity(target.density, bins = bins, col = color, legend = legend, title =self$ref.level, grid = grid, mi = mi, ma = ma)
      p2 <- plotDensity(ref.density, bins = bins, col = color, legend = legend, title =self$target.level, grid = grid, mi = mi, ma = ma)
      
      if (add.ponits){
        if(is.null(condition.per.cell)) stop("'condition.per.cell' must be provided when add points")
    
        emb <- private$getResults('density.emb', 'estimateCellDensity()')
        emb$Z <- 1
        nname1 <- names(condition.per.cell)[condition.per.cell == ref]
        nname1 <- sample(nname1, min(2000, nrow(emb[nname1, ])))
        
        nname2 <- names(condition.per.cell)[condition.per.cell == target]
        nname2 <- sample(nname2, min(2000, nrow(emb[nname2, ])))
        
        p1 <- p1 + geom_point(data = emb[nname1, ], aes(x = x, y = y), col = point.col, size = 0.00001, alpha = 0.2)  
        p2 <- p2 + geom_point(data = emb[nname2, ], aes(x = x, y = y), col = point.col, size = 0.00001, alpha = 0.2) 
      }
      return(list('ref' = p1, 'target' = p2))
    },
    
    

    ##' @description esitmate differential cell density
    ##' @param col color palettes, 4 different color palettes are supported; default is blue-white-red; BWR: blue-white-red;  WR: white-read; B: magma in viridi;
    ##' @param condition.per.cell A two-level factor on the cell names describing the conditions being compared (default: stored vector)
    ##' @method method to cacuated differential cell density of each bin; substract: target density minus ref density; entropy: estimated kl divergence entropy betwwen sample grapups ; t.test: zscore of t-test,global variacen is setting for t.test;     
    diffCellDensity = function(condition.per.cell = NULL, method = 'substract', legend = NULL, grid = TRUE, col = 'BWR', title = NULL, plot = TRUE){
      ref.level <- self$ref.level
      target.level <- self$target.level
      sample.groups <- self$sample.groups
      bins <- private$getResults('bins', 'estimateCellDensity()')
      density.matrix <- private$getResults('density.mat', 'estimateCellDensity()')
      dcounts <- private$getResults('density.counts', 'estimateCellDensity()')

      if (method == 'entropy'){
        if(is.null(condition.per.cell)) stop("'condition.per.cell' must be provided when entropy was used")
      }
      p <- diffCellDensity(density.matrix, dcounts = dcounts, condition.per.cell = condition.per.cell, sample.groups, bins = bins, col = col, target.level = target.level, ref.level =
                             ref.level, method = method, title = title, grid = grid, legend = legend)
      if (plot)
        return(p$fig)
      
      return(p$score)
    },

    #' @title Plot inter-sample expression distance 
    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
    #' @param name Test results to plot (default=expression.shifts)
    #' @param notch Show notches in plot, see ggplot2::geom_boxplot for more info (default=T)
    #' @param cell.groups Named factor with cell names defining groups/clusters (default: stored $cell.groups vector)
    #' @param sample.groups Named sample factor with cell names (default: stored $sample.groups vector)
    #' @param weighted.distance whether to weigh the expression distance by the sizes of cell types (default: TRUE), or show distances for each individual cell type
    #' @return A ggplot2 object
    plotExpressionDistance = function(name='expression.shifts', notch = TRUE, cell.groups = self$cell.groups, sample.groups = self$sample.groups, weighted.distance = TRUE,  min.cells = 10, palette=self$sample.groups.palette) {
      plotExpressionDistance(private$getResults(name, 'estimateExpressionShiftMagnitudes()'), notch = notch, cell.groups = cell.groups, sample.groups = sample.groups, weighted.distance = weighted.distance,  min.cells = min.cells, palette=palette)
    },

    #' @title Plot sample-sample expression distance as a 2D embedding
    #' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
    #' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
    #' @param cell.type If a name of a cell type is specified, the sample distances will be assessed based on this cell type alone. Otherwise (cell.type=NULL, default), sample distances will be estimated as an average distance across all cell types (weighted by the minimum number of cells of that cell type between any two samples being compared)
    #' @param sample.groups Named sample factor with cell names (default: stored $sample.groups )
    #' @method dimension reduction methods (MDS or tSNE ) , default is MDS
    #' @param perplexity tSNE perpexity (default: 4)
    #' @param max_iter tSNE max_iter (default: 1e3)
    #' @param palette a set of colors to use for conditions (default: stored $sample.groups.palette)
    #' @return A ggplot2 object
    plotExpressionDistanceEmbedding = function(name='expression.shifts', sample.groups = self$sample.groups, cell.type = NULL, method = 'tSNE', perplexity=4, max_iter=1e3, palette=self$sample.groups.palette) {
      cluster.shifts <- private$getResults(name, 'estimateExpressionShiftMagnitudes()')
      plotExpressionDistancetSNE(cluster.shifts, sample.groups = sample.groups, cell.type = cell.type, method = method, perplexity=perplexity, max_iter=max_iter, palette=palette)
    },
    
    #' @title get a cell group contour for 
    getContour = function(cell, embedding=self$embedding, cell.groups=self$cell.groups, ...) {
      getContour(cell=cell, emb=embedding,cell.type=cell.groups, ...)
    }
  ),
  private = list(
    checkTestResults=function(name) {
      if (is.null(self$test.results[[name]])) {
        stop("Test result for ", name, " wasn't found")
      }
    },
    
    getResults=function(name,suggestedFunction=NULL) {
      if (is.null(self$test.results[[name]])) {
        msg <- paste0("A result named \"", name, "\" cannot be found.");
        if(!is.null(suggestedFunction)) {
          msg <- paste(msg,"Please first run",suggestedFunction)
        }
        stop(msg)
      } else {
        return(self$test.results[[name]])
      }
    }
  )
)
