#' @import methods dplyr R6
#' @importFrom grDevices colorRampPalette contourLines
#' @importFrom graphics legend
#' @importFrom stats as.dendrogram as.dist as.formula as.hclust cmdscale coef
#' @importFrom stats cor cutree dist hclust lm loadings median model.matrix na.omit p.adjust pf prcomp
#' @importFrom stats qnorm quantile relevel reorder rmultinom runif sd setNames symnum var wilcox.test
#' @importFrom utils stack
#' @importFrom sccore checkPackageInstalled
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

    #' @field sample.meta Data frame with annotation of covariates per sample (default=NULL)
    sample.meta = NULL,

    #' @field full.meta Full sample metadata (default=NULL)
    full.meta = NULL,

    #' @field excluded.samples samples left out of all analyses by `excludeSamples()` (default=NULL: none)
    excluded.samples = NULL,

    #' @field cell.groups Named factor with cell names with cluster per cell (default=NULL)
    cell.groups = NULL,

    #' @field embedding 2D embedding to visualize the cells in (default=NULL)
    embedding = NULL,

    #' @field sample.per.cell Named factor with cell names (default=NULL)
    sample.per.cell = NULL,



    #' @field model.matrices Design model matrices for the analysis (default=NULL)
    model.matrices = NULL,

    #' @field ref.level Reference level for sample.group vector (default=NULL)
    ref.level = NULL,

    #' @field target.level Target/disease level for sample.group vector
    target.level = NULL,

    #' @field sample.ids Character vector of sample identifiers (row names of `sample.meta`)
    sample.ids = NULL,

    #' @field block.vars Character of column name containing variable name to restrict permutations to
    block.vars = NULL,

    #' @field numeric.ref Reference points for numeric covariates in contrasts: "auto" or a named list
    numeric.ref = "auto",

    #' @field perm.method Legacy field (unused): permutation method of the previous implementation
    perm.method = NULL,

    #' @field model The stored model (class `cacoaModel`, see `setModel()`): formula, tests, designs, issues (default=NULL)
    model = NULL,

    #' @field options Analysis options shared by all methods (see `setOptions()`)
    options = NULL,

    #' @field sample.groups Named factor of the primary test's groups per sample (reference, target); derived from the model
    sample.groups = NULL,

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
    #' @param sample.metadata data.frame; rows = samples, columns = covariates. Row names must be sample IDs.
    #' @param formula character or formula (e.g. "~ condition + batch"). Optional; the default is `~ <test variables> + <block.vars>`.
    #' @param test what to test: a metadata column name (`"condition"`), several (`c("condition", "sex")`), `"all"`,
    #'   an explicit comparison (`"condition: disease vs control"` or `c("condition", "disease", "control")`), or a
    #'   structured contrast list (see `buildDesignMatrices()`). Optional: without a model the object can still be
    #'   explored (`checkDesign()`, `screenCovariates()`); set the model later with `setModel()`.
    #' @param contrast expert synonym of `test` (a `c(var, alt, ref)` triple or a structured contrast).
    #' @param dispersion.formula formula for the dispersion (heterogeneity) model (default: the test variables).
    #' @param numeric.ref reference points for numeric covariates in the contrast specification. `"auto"` or named list of numeric anchors (e.g., `list(age=35)`).
    #' @param sample.groups deprecated: named vector of group per sample (old two-group API); becomes a `condition` metadata column and `test = "condition"`.
    #' @param ref.level deprecated: the reference level of `sample.groups` (old two-group API).
    #' @param target.level deprecated: the target level of `sample.groups` (old two-group API).
    #' @param block.vars optional list of covariates in sample.metadata on which to form randomizaton blocks
    #' @param sample.ids character scalar naming the column in `sample.metadata` that contains sample IDs.
    #' @param cell.groups vector Indicates cell groups with cell names (default: extracted from `data.object`)
    #' @param sample.per.cell vector Sample name per cell (default: extracted from `data.object`)
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
    #' # Is it highly recommended that sample.metadata and cell.groups are assigned in the initialization call.
    #' # Here, "con" is a Conos object.
    #' \dontrun{
    #' sample.metadata <- data.frame(condition=c("control","control","disease","disease"), batch=c("batch1","batch1","batch2","batch2"), sample=names(con$samples))
    #' sample.id <- "sample"
    #' formula <- "~ condition + batch"
    #' contrast <- c("condition","disease", "control")
    #' }
    #' # cell.groups should be a named factor where names are cell names corresponding to cell names in the data object.
    #' # For Conos objects, they should overlap with rownames(con$embedding)
    #' \dontrun{
    #' cell.groups <- my.named.annotation.factor
    #' cao <- Cacoa$new(data.object = con, sample.metadata = sample.metadata, sample.id=sample.id, formula = formula, contrast = contrast, cell.groups = cell.groups)
    #' }
    #' @param data.layer Seurat data layer to use for normalized expression (default='scale.data'; 'data' is usually better)
    initialize=function(
      data.object, sample.metadata=NULL, sample.ids=NULL, formula=NULL, test=NULL, contrast=NULL, dispersion.formula=NULL,
      numeric.ref = "auto", block.vars=NULL, cell.groups=NULL, sample.per.cell=NULL, sample.groups.palette=NULL,
      cell.groups.palette=NULL, embedding=NULL, n.cores=1, verbose=TRUE,
      graph.name=NULL, assay.name="RNA", data.layer='scale.data',
      plot.theme=ggplot2::theme_bw(), plot.params=NULL,
      sample.groups=NULL, ref.level=NULL, target.level=NULL
    ) {

      if ('Cacoa' %in% class(data.object)) {
        for (n in ls(data.object)) {
          if (!is.function(get(n, data.object))) assign(n, get(n, data.object), self)
        }
        if (is.null(self$options)) self$options <- cacoaDefaultOptions()
        return(NULL)
      }

      # deprecated two-group API (D35): sample.groups + ref.level / target.level -> a `condition` column and test
      if (!is.null(sample.groups)) {
        .Deprecated(msg = "`sample.groups` / `ref.level` / `target.level` are deprecated: give `sample.metadata` with a condition column and `test = \"condition\"`.")
        if (is.null(names(sample.groups))) stop("`sample.groups` must be named by sample")
        sg <- data.frame(condition = as.character(sample.groups), row.names = names(sample.groups), stringsAsFactors = FALSE)
        if (is.null(sample.metadata)) sample.metadata <- sg
        else {
          if (is.null(rownames(sample.metadata)) && !is.null(sample.ids)) rownames(sample.metadata) <- as.character(sample.metadata[[sample.ids]])
          if (!"condition" %in% names(sample.metadata)) sample.metadata$condition <- sg[rownames(sample.metadata), "condition"]
        }
        if (is.null(test) && is.null(contrast)) {
          levs <- unique(sg$condition)
          if (is.null(ref.level)) ref.level <- chooseReferenceLevel(sg$condition)$level
          if (is.null(target.level)) target.level <- setdiff(levs, ref.level)[1]
          test <- c("condition", target.level, ref.level)
        }
      }

      if (is.null(sample.metadata)) stop("sample.metadata must be provided")
      if (!is.data.frame(sample.metadata)) stop("sample.metadata must be a data.frame")

      # check sample IDs

      samp.names <- rownames(sample.metadata)
      if ((is.null(samp.names) || anyNA(samp.names)) && !is.null(sample.ids)) {
          if (!sample.ids %in% colnames(sample.metadata)) {
            stop(sprintf("`sample.ids` column '%s' not found in sample metadata.", sample.ids))
          }
          samp.names <- as.character(sample.metadata[[sample.ids]])
          if (any(duplicated(samp.names))) {
            stop("duplicate sample identifiers found in sample.metadata. Please ensure all sample IDs are unique.")
          }
          rownames(sample.metadata) <- samp.names
        }
        if (is.null(samp.names) || anyNA(samp.names)) {
          stop("sample identifiers are unavailable: set rownames(sample.metadata) OR provide a valid `sample.id` column.")
        }


      self$options <- cacoaDefaultOptions()
      self$options$n.cores <- n.cores; self$options$verbose <- verbose
      self$numeric.ref <- numeric.ref
      self$block.vars <- block.vars
      self$sample.meta <- self$full.meta <- sample.metadata
      self$sample.ids <- samp.names
      self$n.cores <- n.cores
      self$verbose <- verbose
      self$sample.groups.palette <- sample.groups.palette

      # the model is optional (D25): set it now when asked, otherwise print the metadata audit
      if (!is.null(formula) || !is.null(test) || !is.null(contrast)) {
        self$setModel(formula = formula, test = test, contrast = contrast, dispersion.formula = dispersion.formula,
                      block.vars = block.vars, numeric.ref = numeric.ref, verbose = verbose)
      } else if (verbose) {
        message(formatMetadataSummary(describeMetadata(sample.metadata), sample.metadata))
        message("Next: cao$checkDesign() and cao$screenCovariates(), or cao$setModel(test = \"<variable>\").")
      }

      # interpret different types of data objects

      if ("Seurat" %in% class(data.object)) {
        if (is.null(sample.metadata) || is.null(sample.per.cell)){
          stop("Both sample.metadata and sample.per.cell must be provided for Seurat objects")
        }
        data.object$sample.per.cell <- sample.per.cell
        
        if (is.null(graph.name)){
          warning("No graph.name provided. The algorithm will use the first available graph.")
        }

        data.object@misc$graph.name <- graph.name
        is.connected <- extractCellGraph(data.object) %>% igraph::is_connected()
        if (!is.connected) {
          warning("The provided graph is not connected. Cluster-free methods will not work. Using NN graph instead of SNN may help.")
          data.object@misc$graph.name <- NULL
        }

        data.object@misc$assay.name <- assay.name
        data.object@misc$data.layer <- data.layer
        if (data.layer == 'scale.data') {
          warning("Using data.layer='scale.data' is not recommended, as it usually contains only highly-variable genes. Consider using data.layer='data' instead.")
        }
      } else if (('Conos' %in% class(data.object))) {
        if (!is.null(graph.name)) {
          warning("graph.name is not supported for Conos objects")
        }
      } else {
        warning("Many function may be not supported for an object of class ", class(data.object));
        if (is.null(sample.metadata) || is.null(sample.per.cell) || is.null(cell.groups)){
          stop("All sample.metadata, sample.per.cell and cell.groups must be provided")
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
      if(length(setdiff(unique(as.character(self$sample.per.cell)), samp.names)) > 0) {
          stop("sample identifiers found in sample.metadata do not match sample.per.cell")
      }

      if (is.null(self$sample.groups.palette)) {
        self$sample.groups.palette <- c("#d73027", "#4575b4")
        if (!is.null(self$target.level)) names(self$sample.groups.palette) <- c(self$target.level, self$ref.level)
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


    #' @description Set analysis options shared by all methods (explicit call arguments still win)
    #' @param ... options to change: `n.permutations` (999), `seed` (1; `NULL` to draw from R's RNG), `alpha` (0.05),
    #'   `p.adjust.method` ("BH"), `n.cores`, `verbose`, `plot.provenance` (TRUE), `min.samp.per.level` (3),
    #'   `dist` ("cor"), `permutation` ("auto"), `numeric.step` (1)
    #' @return the full option list, invisibly
    #' @examples
    #' \dontrun{
    #' cao$setOptions(n.permutations = 1999, seed = 7, n.cores = 8)
    #' }
    setOptions = function(...) {
      new <- list(...)
      if (is.null(self$options)) self$options <- cacoaDefaultOptions()
      def <- cacoaDefaultOptions(); miss <- setdiff(names(def), names(self$options))   # options added after the object was created
      if (length(miss)) self$options[miss] <- def[miss]
      unknown <- setdiff(names(new), names(self$options))
      if (length(unknown) || (length(new) && is.null(names(new)))) stop("unknown option(s): ", paste(unknown, collapse = ", "),
                                                                      ". Known: ", paste(names(self$options), collapse = ", "))
      opts <- self$options
      for (nm in names(new)) opts[nm] <- list(new[[nm]])   # list() so that NULL (seed) is stored
      checkOptionValues(opts)
      self$options <- opts
      if ("n.cores" %in% names(new)) self$n.cores <- new$n.cores
      if ("verbose" %in% names(new)) self$verbose <- new$verbose
      invisible(self$options)
    },

    #' @description Get analysis options
    #' @param name optional option name; all options when `NULL`
    #' @return the option value, or the list of all options
    getOptions = function(name = NULL) {
      if (is.null(self$options)) self$options <- cacoaDefaultOptions()
      if (is.null(name)) self$options else self$options[[name]]
    },

    #' @description Set (or change) the stored model used by every analysis
    #'
    #' Resolves what to test (`test` grammar), builds the design, checks it, determines the permutation plan and
    #' stores the result as `cao$model`. The constructor calls this when it receives `formula` / `test` / `contrast`;
    #' call it directly to change the model after exploration (`checkDesign()`, `screenCovariates()`).
    #' @param formula location formula (default: `~ <test variables> + <block.vars>`)
    #' @param test what to test (see `Cacoa$new()`)
    #' @param contrast expert synonym of `test`
    #' @param dispersion.formula dispersion formula (default: the test variables)
    #' @param block.vars metadata columns defining permutation strata (default: stored)
    #' @param numeric.ref anchors for numeric covariates (default: stored)
    #' @param numeric.step step for numeric tests (default: option)
    #' @param permutation permutation scheme (default: option `"auto"`)
    #' @param verbose print the model summary (default: option)
    #' @return the model (class `cacoaModel`), invisibly
    #' @examples
    #' \dontrun{
    #' cao$setModel(~ condition + batch + sex, test = "condition")
    #' cao$setModel(test = "condition: disease vs control", block.vars = "batch")
    #' }
    setModel = function(formula = NULL, test = NULL, contrast = NULL, dispersion.formula = NULL, block.vars = self$block.vars,
                        numeric.ref = self$numeric.ref, numeric.step = NULL, permutation = NULL, verbose = NULL) {
      verbose <- private$opt("verbose", verbose)
      if (is.null(test) && is.null(contrast) && !is.null(self$model)) { test <- self$model$test.spec$test; contrast <- self$model$test.spec$contrast }
      m <- private$buildModel(formula = formula, test = test, contrast = contrast, dispersion.formula = dispersion.formula,
                              block.vars = block.vars, permutation = permutation, numeric.ref = numeric.ref, numeric.step = numeric.step)
      self$model <- m
      self$block.vars <- block.vars; self$numeric.ref <- numeric.ref
      private$syncLegacyFields()
      # a tested variable with one value in every level of a grouping term: no analysis can run (reported here, where the formula is given)
      btw <- m$issues[m$issues$severity == "warning" & grepl("takes one value within every level", m$issues$message, fixed = TRUE), , drop = FALSE]
      if (nrow(btw)) warning(paste0(btw$message[1], ". ", sub("^c", "C", btw$suggestion[1]), "."), call. = FALSE)
      if (verbose) print(m)
      invisible(m)
    },

    #' @description Get the stored model
    #' @return the model (class `cacoaModel`) or `NULL`
    getModel = function() self$model,

    #' @description Check the study design from the sample metadata (no expression data needed)
    #'
    #' Reports associations among covariates, balance of each covariate against the test variable, structural
    #' problems (a covariate that fully determines the test variable, nesting, single-sample or empty cells, rank
    #' deficiency), collinearity (GVIF), the degrees-of-freedom budget, permutation feasibility per test and
    #' over-adjustment candidates (covariates strongly tied to the test variable).
    #' @param formula optional model to check instead of the stored one (a temporary model is built)
    #' @param test optional model to check instead of the stored one (a temporary model is built)
    #' @param contrast optional model to check instead of the stored one (a temporary model is built)
    #' @param covariates covariates to examine (default: all usable metadata columns)
    #' @param block.vars permutation strata (default: stored)
    #' @param verbose print the report (default: option)
    #' @return object of class `cacoaDesignCheck` (see [checkDesign()]), also stored in `cao$cache$design.check`
    #' @examples
    #' \dontrun{
    #' cao$checkDesign()
    #' cao$plotDesign("associations")
    #' }
    checkDesign = function(formula = NULL, test = NULL, contrast = NULL, covariates = NULL, block.vars = NULL, verbose = NULL) {
      verbose <- private$opt("verbose", verbose)
      bare <- is.null(formula) && is.null(contrast) && is.character(test) && length(test) == 1 && test %in% names(self$sample.meta)
      model <- if (bare) NULL else if (!is.null(formula) || !is.null(test) || !is.null(contrast)) {
        tryCatch(private$resolveModel(formula = formula, test = test, contrast = contrast, block.vars = block.vars, verbose = FALSE, what = "checkDesign()"),
                 error = function(e) { if (verbose) message("Model could not be built: ", conditionMessage(e)); NULL })
      } else self$model
      test.variable <- if (is.null(model) && !is.null(test) && is.character(test) && length(test) == 1 && test %in% names(self$sample.meta)) test else NULL
      chk <- checkDesign(self$sample.meta, model = model, test.variable = test.variable, covariates = covariates,
                         block.vars = block.vars %||% self$block.vars)
      self$cache$design.check <- chk
      if (verbose) print(chk)
      invisible(chk)
    },

    #' @description Plot the design check
    #' @param type `"associations"` (default), `"balance"` or `"issues"`
    #' @param covariates optional subset of covariates
    #' @param ... passed to [plotDesignCheck()]
    #' @return ggplot2 object
    plotDesign = function(type = c("associations", "balance", "issues"), covariates = NULL, ...) {
      type <- match.arg(type)
      chk <- self$cache$design.check
      if (is.null(chk)) chk <- self$checkDesign(verbose = FALSE)
      plotDesignCheck(chk, type = type, meta = self$sample.meta, covariates = covariates, plot.theme = self$plot.theme, ...)
    },

    #' @description Describe the sample metadata columns (type, levels, missing values, role as a covariate)
    #' @return data.frame, see [describeMetadata()]
    describeMetadata = function() describeMetadata(self$sample.meta),

    #' @description Leave samples out of every analysis run afterwards (expression shifts, composition, differential expression,
    #'   the cluster-free maps): their cells are treated as unassigned and their rows are removed from the working sample
    #'   metadata; the model is rebuilt on the remaining samples. Nothing is excluded automatically: the sensitivity and fit
    #'   diagnostics name samples, the decision is the user's. Called without samples, it restores the full set. Stored results
    #'   are not recomputed.
    #' @param samples character vector of sample ids to exclude (the complete list: a second call replaces the first); `NULL`
    #'   or empty restores all samples
    #' @param verbose print what was done
    #' @return the object, invisibly; `cao$excluded.samples` holds the excluded ids
    excludeSamples = function(samples = NULL, verbose = NULL) {
      verbose <- private$opt("verbose", verbose)
      if (is.null(private$all.samples)) private$all.samples <- list(sample.per.cell = self$sample.per.cell, sample.meta = self$sample.meta)
      full <- private$all.samples; samples <- unique(as.character(samples))
      unknown <- setdiff(samples, rownames(full$sample.meta))
      if (length(unknown)) stop("not in the sample metadata: ", paste(unknown, collapse = ", "))
      spc <- full$sample.per.cell; nm <- names(spc); spc <- as.character(spc); spc[spc %in% samples] <- NA
      keep.lv <- setdiff(if (is.factor(full$sample.per.cell)) levels(full$sample.per.cell) else unique(stats::na.omit(as.character(full$sample.per.cell))), samples)
      self$sample.per.cell <- setNames(factor(spc, levels = keep.lv), nm)
      self$sample.meta <- full$sample.meta[!rownames(full$sample.meta) %in% samples, , drop = FALSE]
      self$sample.ids <- rownames(self$sample.meta)
      self$excluded.samples <- if (length(samples)) samples else NULL
      if (!is.null(self$model)) { m <- self$model
        self$setModel(formula = if (isTRUE(m$formula.default)) NULL else (m$formula.full %||% m$formula), test = m$test.spec$test, contrast = m$test.spec$contrast, dispersion.formula = m$dispersion.formula, verbose = FALSE) }
      if (verbose) message(if (length(samples)) sprintf("Excluded from all analyses run from now on: %s (%d samples remain).", paste(samples, collapse = ", "), nrow(self$sample.meta)) else "All samples are included again.",
                           " Stored results were computed with the previous sample set: rerun the estimates you need.")
      invisible(self)
    },

    #' @description Sample groups of a test: a factor (reference, target) for two-level contrasts, the variable
    #'   itself for numeric tests or term tests
    #' @param test test index or label within the stored model (default: the primary test)
    #' @param model model to take the test from (default: stored)
    #' @return named vector over samples (`NULL` when the test has no variable, e.g. coefficient contrasts)
    getSampleGroups = function(test = NULL, model = self$model) {
      if (is.null(model)) return(NULL)
      t <- model$tests[[matchModelTest(model, test)[1]]]
      if (is.null(t) || is.null(t$variable) || is.na(t$variable)) return(NULL)
      x <- self$sample.meta[[t$variable]]
      if (t$kind == "contrast" && is.null(t$step)) {
        x <- setNames(as.character(x), rownames(self$sample.meta)); keep <- x %in% t$levels[c("ref", "alt")]
        return(factor(x[keep], levels = unname(t$levels[c("ref", "alt")])))
      }
      names(x) <- rownames(self$sample.meta)
      if (is.factor(x) || is.character(x)) droplevels(factor(x)) else x
    },

    #' @description Robustness of a result: alternative covariate sets and single-sample influence
    #'
    #' Reruns the first test of the stored model under alternative adjustment sets (unadjusted; current; current
    #' minus each covariate; current plus each of the top screened covariates not yet in the model; all screened,
    #' if the degrees of freedom allow) and refits every model with each sample left out. The two questions are kept apart:
    #' the size of the effect (stable, changes by x%, sign depends on the model, or driven by sample S) and its significance
    #' (significant in all / k of m / no models). One line per cell type gives both.
    #'
    #' `analysis = "expression.shifts"`: the expression shifts on the cached distances; the size is the distance-model
    #' magnitude, as in the main result, the significance by the test of the main result (with `statistic = "gsd"` the gsd test
    #' is rerun under every alternative model, with the main result's gsd settings).
    #' `analysis = "composition"`: the compositional contrast ([lmCoda()]) refitted on the stored cell-type counts; the size
    #' is the loading of each cell type (log2, relative to the geometric mean of all cell types), the significance its
    #' permutation test (see [checkCompositionSensitivity()]).
    #' @param analysis `"expression.shifts"` (default) or `"composition"`
    #' @param covariate.sets `"auto"` (default, see above) or a named list of location formulas
    #' @param screen.name results slot of a covariate screen used to propose additions (default "covariate.screen")
    #' @param top.k number of screened covariates to try adding (default 3)
    #' @param name results slot (default "sensitivity" for the expression shifts, "composition.sensitivity" for the composition)
    #' @param shifts.name results slot of the expression shifts the sensitivity refers to (default "expression.shifts")
    #' @param coda.name results slot of the compositional result the sensitivity refers to (default "coda")
    #' @param jackknife leave-one-sample-out refits under every alternative model (default TRUE), which give each estimate its
    #'   jackknife standard error and leave-one-out range (shown by `plotSensitivity()`); the current model always has them
    #' @param n.permutations number of permutations of the refits (default: the object option `n.permutations`, 999); the
    #'   gsd refits use the main result's count, so that their p-values have the resolution of the current model's
    #' @param significance expression shifts only: `"main"` (default: the test of the stored result, gsd or distance) or `"distance"`
    #'   (the distance model's own test, faster: no gsd refits)
    #' @param cell.groups cell annotation of the gsd refits (default: the object's; the one the stored result was computed with)
    #' @param verbose defaults from options (permutations capped at 499)
    #' @param n.cores defaults from options (permutations capped at 499)
    #' @param seed defaults from options (permutations capped at 499)
    #' @return `cacoaSensitivity` (see [checkSensitivity()]), also stored in `cao$test.results[[name]]`
    checkSensitivity = function(analysis = c("expression.shifts", "composition"), covariate.sets = "auto", screen.name = "covariate.screen", top.k = 3, name = NULL,
                                shifts.name = "expression.shifts", coda.name = "coda", jackknife = TRUE, n.permutations = NULL, significance = c("main", "distance"),
                                cell.groups = self$cell.groups, verbose = NULL, n.cores = NULL, seed = NULL) {
      analysis <- match.arg(analysis); significance <- match.arg(significance)
      name <- name %||% if (analysis == "composition") "composition.sensitivity" else "sensitivity"
      verbose <- private$opt("verbose", verbose); n.cores <- private$opt("n.cores", n.cores)
      n.permutations <- n.permutations %||% private$opt("n.permutations")
      seed <- if (missing(seed)) private$opt("seed") else seed
      screen <- self$test.results[[screen.name]]
      formulas <- if (identical(covariate.sets, "auto")) NULL else covariate.sets
      if (analysis == "composition") {
        res <- self$test.results[[coda.name]]
        if (is.null(res)) stop("run estimateCellLoadings() first")
        if (identical(res$kind, "term") || is.null(res$contrast)) stop("the compositional sensitivity is defined for contrast tests")
        st <- res$settings %||% list()
        if (verbose) message(sprintf("checkSensitivity(): composition loadings refitted under alternative covariate sets%s", if (jackknife) " and with each sample left out" else ""))
        out <- checkCompositionSensitivity(res$cnts, res$model, self$sample.meta, formulas = formulas, screen = screen, current = res,
                                           n.permutations = n.permutations, seed = seed %||% sample.int(.Machine$integer.max, 1), alpha = private$opt("alpha"),
                                           top.k = top.k, jackknife = jackknife, perm.method = res$perm.method %||% "block",
                                           zero.pseudocount = st$zero.pseudocount %||% 0.1, scale.pseudocount = st$scale.pseudocount %||% FALSE, basis.type = st$basis.type %||% "default",
                                           ref.method = st$ref.method %||% "block", ref.block.thresh = st$ref.block.thresh %||% 0.05, ref.size.frac = st$ref.size.frac %||% 0.5)
        self$test.results[[name]] <- out
        if (verbose) print(out)
        return(invisible(out))
      }
      res <- self$test.results[[shifts.name]]
      if (is.null(res)) stop("run estimateExpressionShiftMagnitudes() first")
      model <- res$model
      D.list <- res$distances
      infl <- if (!is.null(res$influence[[1]])) list(influence = res$influence[[1]], wide = res$wide[[1]]) else NULL
      cur <- res$results[res$results$test.id == res$results$test.id[1], ]        # the stored result stands in for the "current" model
      gsd.run <- identical(res$settings$statistic, "gsd") && !is.null(res$settings$gsd) && "p.distance" %in% names(cur)
      sig.fun <- NULL; sig.label <- "distance test"
      if (gsd.run && significance == "distance") {   # compare the distance model's own tests across the models
        i <- which(is.finite(cur$p.distance)); cur$p[i] <- cur$p.distance[i]; cur$padj[i] <- cur$padj.distance[i]
        if (verbose) message("checkSensitivity(): magnitudes and significance from the distance model (the stored result was tested by gsd)")
      } else if (gsd.run) {                          # significance by gsd, as in the stored result, refitted under every alternative model
        gs <- res$settings$gsd; sig.label <- "gsd"
        inp <- private$getGsdCellTypeInput(res$genes.tested)
        types <- unique(cur$celltype[cur$effect == "shift"])
        cg <- cell.groups[!is.na(cell.groups) & as.character(cell.groups) %in% types]
        n.gsd <- gs$n.permutations
        if (verbose) message(sprintf("checkSensitivity(): magnitudes from the distance model, significance from gsd (as the stored result; %d permutations per model)", n.gsd))
        sig.fun <- function(m) {
          g <- gsdShiftsPerCellType(inp$cm, self$sample.per.cell, cg, m, inp$lib, min.cells.per.sample = res$settings$min.cells.per.sample %||% 5,
                                    min.samp.per.level = gs$min.samp.per.level, neighbourhood.gene.frac = gs$neighbourhood.gene.frac,
                                    permutation = gs$permutation, n.permutations = n.gsd, block.vars = gs$block.vars, seed = gs$seed,
                                    robust = gs$robust %||% "none", winsor.k = gs$winsor.k %||% 3, gsd.d0 = gs$gsd.d0 %||% 3, influence = FALSE, n.null = 0,
                                    n.cores = n.cores, verbose = FALSE)
          g$results[, c("celltype", "p", "padj")]
        }
      }
      out <- checkSensitivity(D.list, model, self$sample.meta, formulas = formulas, screen = screen, influence = infl, current.rows = cur, dist = res$settings$dist,
                              permutation = res$settings$permutation, n.permutations = n.permutations, seed = seed %||% sample.int(.Machine$integer.max, 1),
                              alpha = private$opt("alpha"), min.samp.per.level = res$settings$min.samp.per.level, n.cores = n.cores, top.k = top.k,
                              robust = res$settings$robust %||% "none", na.mode = res$settings$na.mode %||% "drop", robust.k = res$settings$robust.k %||% 1.345,
                              significance = sig.fun, significance.label = sig.label, jackknife = jackknife)
      self$test.results[[name]] <- out
      if (verbose) print(out)
      invisible(out)
    },

    #' @description Forest plot of a sensitivity analysis (see [plotSensitivity()])
    #' @param analysis `"expression.shifts"` (default) or `"composition"`: which stored sensitivity to plot when `name` is not given
    #' @param name results slot (default "sensitivity" or "composition.sensitivity", by `analysis`)
    #' @param interval `"jackknife"` (default): jackknife 95 % intervals and leave-one-sample-out ranges of the estimates, where
    #'   computed; `"none"`
    #' @param ... passed to [plotSensitivity()] (`effect`, `cell.types`, `normalized`)
    #' @return ggplot2 object
    plotSensitivity = function(analysis = c("expression.shifts", "composition"), name = NULL, interval = c("jackknife", "none"), ...) {
      analysis <- match.arg(analysis); interval <- match.arg(interval)
      name <- name %||% if (analysis == "composition") "composition.sensitivity" else "sensitivity"
      plotSensitivity(private$getResults(name, "checkSensitivity()"), interval = interval, plot.theme = self$plot.theme, ...)
    },

    #' @description Screen covariates: which metadata variables are associated with sample-level variation, per cell type
    #'
    #' Exploratory, term-level and contrast-free. For every covariate and cell type the between-sample variation it
    #' explains is reported alone (`marginal`) and after adjusting for the other screened covariates (`partial`),
    #' as chance-corrected R2 with permutation p-values (BH over the whole grid) and a per-covariate global p-value
    #' across cell types; a dispersion association is reported alongside. Ends with a suggested model (not applied).
    #' @param covariates covariates to screen (default: all usable metadata columns); transform strings allowed
    #' @param space `"expression"` (per cell type; default) or `"composition"` (ILR coordinates of the cell-type counts, one
    #'   row); `"expression.shifts"` is accepted as the earlier name
    #' @param summary.rows with `space = "expression"`, rows screened next to the cell types on the combined sample spaces
    #'   (see `getSampleDistanceMatrix()`): any of `"expression"` (all cell types, weighted by cells), `"expression.equal"`
    #'   (all cell types, equal weights), `"pseudobulk"` (all cells pooled) and `"composition"`; default
    #'   `c("expression", "composition")`. They are adjusted among themselves and do not enter the global p-values or the
    #'   suggested model.
    #' @param mode `"both"` (default), `"partial"` or `"marginal"`
    #' @param adjust.for optional formula of covariates always adjusted for
    #' @param dist as in `estimateExpressionShiftMagnitudes()`
    #' @param n.pcs as in `estimateExpressionShiftMagnitudes()`
    #' @param cell.groups as in `estimateExpressionShiftMagnitudes()`
    #' @param min.cells.per.sample as in `estimateExpressionShiftMagnitudes()`
    #' @param min.gene.frac as in `estimateExpressionShiftMagnitudes()`
    #' @param genes as in `estimateExpressionShiftMagnitudes()`
    #' @param min.samples.per.type cell types with fewer samples are skipped (default 6)
    #' @param p.values `"permutation"` (default) or `"analytic"` (instant preview, approximate)
    #' @param n.permutations number of permutations (default: the object option `n.permutations`, 999)
    #' @param name results slot (default "covariate.screen")
    #' @param verbose defaults from options
    #' @param n.cores defaults from options
    #' @param seed defaults from options
    #' @return `cacoaCovariateScreen` (see [screenCovariates()]), also stored in `cao$test.results[[name]]`
    #' @examples
    #' \dontrun{
    #' cao$screenCovariates(); cao$plotCovariateScreen(); cao$plotCovariateSummary()
    #' }
    #' @param robust `"none"` (default, option `robust`), `"huber"` or `"winsor"`: robust location / dispersion statistics
    screenCovariates = function(covariates = NULL, space = c("expression", "composition"), mode = c("both", "partial", "marginal"),
                                adjust.for = NULL, dist = NULL, n.pcs = NULL, cell.groups = self$cell.groups, min.cells.per.sample = 10,
                                min.gene.frac = 0.01, genes = NULL, min.samples.per.type = 6, summary.rows = c("expression", "composition"),
                                n.permutations = NULL, name = "covariate.screen", verbose = NULL, n.cores = NULL, seed = NULL, robust = NULL) {
      space <- normalizeSampleSpace(space[1]); mode <- match.arg(mode); robust <- private$opt("robust", robust)
      if (space == "pseudobulk") stop("space must be 'expression' or 'composition' (the pooled pseudobulk is a summary row: summary.rows = 'pseudobulk')")
      bad <- setdiff(summary.rows, c("expression", "expression.equal", "pseudobulk", "composition"))
      if (length(bad)) stop("unknown summary.rows: ", paste(bad, collapse = ", "))
      verbose <- private$opt("verbose", verbose); n.cores <- private$opt("n.cores", n.cores); dist <- private$opt("dist", dist)
      n.permutations <- n.permutations %||% private$opt("n.permutations")
      seed <- if (missing(seed)) private$opt("seed") else seed
      D.list <- private$distancesForSpace(space, dist = dist, n.pcs = n.pcs, cell.groups = cell.groups, min.cells.per.sample = min.cells.per.sample,
                                          min.gene.frac = min.gene.frac, genes = genes, verbose = verbose)
      if (space == "composition") attr(D.list[[1]], "dist.type") <- "l2"
      if (space == "expression" && length(summary.rows)) {   # combined sample spaces, screened as summary rows
        sums <- list()
        for (sr in summary.rows) {
          D <- tryCatch(switch(sr,
            expression = private$sampleSpaceFromList(D.list, NULL, NULL, "cells"),
            expression.equal = private$sampleSpaceFromList(D.list, NULL, NULL, "equal"),
            pseudobulk = private$distancesForSpace("pseudobulk", dist = dist, n.pcs = n.pcs, min.cells.per.sample = min.cells.per.sample,
                                                   min.gene.frac = min.gene.frac, genes = genes, verbose = verbose)[[1]],
            composition = private$distancesForSpace("composition", dist = "l2")[[1]]), error = function(e) {
              if (verbose) message("summary row '", sr, "' skipped: ", conditionMessage(e)); NULL })
          if (is.null(D) || !any(abs(D) > 1e-12, na.rm = TRUE)) next   # absent, or no variation between the samples
          attr(D, "dist.type") <- if (sr == "composition") "l2" else dist; attr(D, "summary") <- TRUE
          sums[[SAMPLE.SPACE.LABELS[[sr]]]] <- D
        }
        D.list <- c(sums, D.list)
      }
      test.variable <- if (!is.null(self$model)) self$model$tests[[1]]$variable else self$cache$design.check$test.variable
      if (!is.null(test.variable) && is.na(test.variable)) test.variable <- NULL
      res <- screenCovariates(D.list, self$sample.meta, covariates = covariates, mode = mode, adjust.for = adjust.for, dist = if (space == "composition") "l2" else dist,
                              test.variable = test.variable, n.permutations = n.permutations,
                              min.samples.per.type = min.samples.per.type, alpha = private$opt("alpha"), seed = seed %||% sample.int(.Machine$integer.max, 1),
                              n.cores = n.cores, verbose = verbose, robust = robust, robust.k = private$opt("robust.k"))
      res$settings$space <- space
      self$test.results[[name]] <- res
      invisible(res)
    },

    #' @description Screen covariates against the residuals of a stored cell-density or cluster-free DE result
    #'
    #' The estimator must have been run with `return.residuals = TRUE`. Residuals are turned into sample distances
    #' and screened marginally for every usable covariate that is not in the model (see
    #' [screenResidualCovariates()]); `plotCovariateScreen(name = <result.name>)` plots the result.
    #' @param name name of the stored result (`"cell.density"` or `"cluster.free.de"`)
    #' @param covariates covariates to screen (default: all usable columns not in the model)
    #' @param n.permutations number of permutations (default: option, at most 499)
    #' @param seed integer seed (default: option)
    #' @param n.cores cores (default: option)
    #' @param verbose print the summary (default: option)
    #' @param robust `"none"`, `"huber"` or `"winsor"` (default: option)
    #' @param result.name where to store the screen (default: `<name>.residual.screen`)
    #' @return object of class `cacoaCovariateScreen`
    #' @examples
    #' \dontrun{
    #' cao$estimateDiffCellDensity(return.residuals = TRUE); cao$screenResiduals("cell.density")
    #' cao$plotCovariateScreen(name = "cell.density.residual.screen")
    #' }
    screenResiduals = function(name = "cell.density", covariates = NULL, n.permutations = NULL, seed = NULL, n.cores = NULL, verbose = NULL,
                               robust = NULL, result.name = paste0(name, ".residual.screen")) {
      verbose <- private$opt("verbose", verbose); n.cores <- private$opt("n.cores", n.cores); robust <- private$opt("robust", robust)
      res <- private$getResults(name, "the estimator with return.residuals = TRUE")
      if (is.null(res$residuals) && !is.null(res$diff$permutation)) res <- res$diff$permutation     # density results keep the test under $diff
      resid <- res$residuals
      if (is.null(resid)) stop("result '", name, "' carries no residuals: re-run the estimator with return.residuals = TRUE")
      model <- res$model %||% self$model
      if (is.null(model)) stop("no model is set")
      tv <- model$tests[[1]]$variable; if (!is.null(tv) && is.na(tv)) tv <- NULL
      out <- screenResidualCovariates(resid, self$sample.meta, covariates = covariates, exclude = all.vars(model$formula),
                                      n.permutations = n.permutations %||% private$opt("n.permutations"),
                                      seed = (if (missing(seed)) private$opt("seed") else seed) %||% sample.int(.Machine$integer.max, 1),
                                      n.cores = n.cores, test.variable = tv, robust = robust, robust.k = private$opt("robust.k"), verbose = verbose)
      out$settings$source <- name
      self$test.results[[result.name]] <- out
      invisible(out)
    },

    #' @description Heatmap of the covariate screen (see [plotCovariateScreen()])
    #' @param name results slot (default "covariate.screen")
    #' @param ... passed to [plotCovariateScreen()] (`effect`, `mode`, `value`, `cell.types`, `covariates`, `cluster`, `alpha`)
    #' @return ggplot2 object
    plotCovariateScreen = function(name = "covariate.screen", ...) {
      plotCovariateScreen(private$getResults(name, "screenCovariates()"), plot.theme = self$plot.theme, ...)
    },

    #' @description Summary bars of the covariate screen (see [plotCovariateSummary()])
    #' @param name results slot (default "covariate.screen")
    #' @param ... passed to [plotCovariateSummary()]
    #' @return ggplot2 object
    plotCovariateSummary = function(name = "covariate.screen", ...) {
      plotCovariateSummary(private$getResults(name, "screenCovariates()"), plot.theme = self$plot.theme, ...)
    },

    #' @description Variance partition of 2-4 covariates per cell type: unique, shared and residual fractions
    #'   of the between-sample variation (chance-corrected). For questions like "is this disease or medication?"
    #' @param covariates 2-4 covariates
    #' @param space `"expression.shifts"` (default) or `"composition"`
    #' @param cell.types optional subset
    #' @param dist as in `estimateExpressionShiftMagnitudes()`
    #' @param n.pcs as in `estimateExpressionShiftMagnitudes()`
    #' @param return.table return the partition table instead of the plot
    #' @return ggplot2 object (or a matrix: cell types x components)
    plotVariancePartition = function(covariates, space = c("expression.shifts", "composition"), cell.types = NULL, dist = NULL, n.pcs = NULL,
                                     return.table = FALSE) {
      space <- match.arg(space)
      if (length(covariates) < 2 || length(covariates) > 4) stop("give 2 to 4 covariates")
      dist <- private$opt("dist", dist)
      D.list <- private$distancesForSpace(space, dist = dist, n.pcs = n.pcs, cell.groups = self$cell.groups, verbose = FALSE)
      if (!is.null(cell.types)) D.list <- D.list[intersect(cell.types, names(D.list))]
      parts <- lapply(D.list, function(D) variancePartition(D, self$sample.meta, covariates, dist = if (space == "composition") "l2" else dist))
      skipped <- names(parts)[vapply(parts, is.null, logical(1))]; parts <- Filter(Negate(is.null), parts)
      if (length(skipped)) message("skipped (fewer than 6 samples or a constant covariate): ", paste(skipped, collapse = ", "))
      if (!length(parts)) stop("no cell type could be partitioned")
      if (return.table) return(do.call(rbind, parts))
      plotVariancePartition(parts, plot.theme = self$plot.theme)
    },

    #' @description Drill-down for one cell type: MDS of the adjusted distances, the pair-distance boxplots that
    #'   define shift / var / total, per-sample dispersion by group, and the level x level table for term tests
    #' @param cell.type cell type
    #' @param test test label or index (default: first test)
    #' @param label.samples label the samples in the distance panel (default FALSE)
    #' @param name results slot (default "expression.shifts")
    #' @return cowplot grid of ggplot2 panels
    plotShiftDetail = function(cell.type, test = NULL, name = "expression.shifts", label.samples = FALSE) {
      res <- private$getResults(name, "estimateExpressionShiftMagnitudes()")
      if (is.null(test)) test <- 1
      tn <- if (is.numeric(test)) names(res$wide)[test] else test
      eff <- if (tn %in% names(res$fits)) res$fits[[tn]][[cell.type]] else NULL
      if (is.null(eff) && identical(res$settings$statistic, "gsd")) stop("plotShiftDetail() shows the distance test: rerun with statistic = 'distance' (gsd results: per-gene changes in $genes)")
      if (is.null(eff) || !isTRUE(eff$ok)) stop("no fit for cell type '", cell.type, "' in test '", tn, "'")
      groups <- self$getSampleGroups(tn, model = res$model)
      w <- res$wide[[tn]]; w <- if (!is.null(w)) w[match(cell.type, w$celltype), , drop = FALSE] else NULL
      effects <- if (!is.null(w) && nrow(w) == 1 && all(c("shift.norm", "var.norm", "total.norm") %in% names(w)))
        c(shift = w$shift.norm, var = w$var.norm, total = w$total.norm) else NULL
      plotShiftDetailPanels(eff, adjusted = res$adjusted.distances[[tn]][[cell.type]], groups = groups, dist = res$settings$dist,
                            title = sprintf("%s: %s", cell.type, tn), palette = self$sample.groups.palette, plot.theme = self$plot.theme,
                            label.samples = label.samples, effects = effects)
    },

    #' @description Expression shift magnitudes per cell type
    #'
    #' With `statistic = "gsd"` (default), every gene of a cell type is tested by a per-sample linear model of its
    #' expression profile, standardized by a moderated standard error (sampling-noise model and residual variance), and
    #' the gene z-scores are combined adaptively into one statistic per cell type (see [gsdShiftsPerCellType()]), tested
    #' by relabeling the samples. The magnitude of the shift is reported separately: `estimate` is the bias-corrected mean
    #' squared change in log10 expression over the type's genes (independent of the numbers of cells), `estimate.norm`
    #' the same change relative to the donor-level variance; per-gene changes and z-scores are in `genes`. Whole-factor
    #' (term) tests use the distance test below; `plotShiftDetail()` and `checkSensitivity()` need the distance test.
    #'
    #' With `statistic = "distance"`, pseudobulk profiles per sample are compared through sample-sample distances
    #' (gene-centred cosine by default). An individual-level model (location: the design formula; dispersion:
    #' `dispersion.formula`) is fitted to the distances and, for each test of the model, three effects are
    #' estimated and tested by permuting sample labels:
    #' \itemize{
    #'   \item \code{shift}: a common change in direction between the compared settings (squared distance between
    #'     the group means, corrected for within-group variability);
    #'   \item \code{var}: the change in within-group heterogeneity (dispersion);
    #'   \item \code{total}: how much farther samples of the target setting are from the reference than reference
    #'     samples are from each other (\code{shift + var/2}).
    #' }
    #' Term tests (a factor with more than two levels) report \code{location} and \code{dispersion} instead, plus
    #' a pairwise table over the levels. Unset arguments fall back to the stored options (`cao$getOptions()`).
    #'
    #' @param test what to test (default: the stored model's tests). Giving `test`, `formula` or `contrast` builds
    #'   a temporary model for this call; `cao$setModel()` changes the stored one.
    #' @param formula location formula for a temporary model
    #' @param contrast expert synonym of `test`
    #' @param dispersion.formula dispersion formula for a temporary model
    #' @param dist `"cor"` (gene-centred cosine; default from options), `"l2"` or `"l1"` (`l1` does not separate shift
    #'   from dispersion and is discouraged)
    #' @param permutation `"auto"`, `"block"` or `"freedman-lane"` (default from options)
    #' @param n.permutations number of permutations (default from options)
    #' @param block.vars metadata columns defining additional permutation strata (default: stored)
    #' @param cell.groups cell annotations (default: stored)
    #' @param sample.per.cell cell annotations (default: stored)
    #' @param n.pcs optional number of principal components to reduce the pseudobulk profiles to before the distance
    #' @param min.cells.per.sample minimum cells per (cell type, sample) of the test (default 5 with gsd, 10 with the distance test)
    #' @param min.cells.magnitude with `statistic = "gsd"`: minimum cells per (cell type, sample) of the distance model that gives
    #'   the magnitudes, the `var` / `total` tests, the noise bands and the jackknife (default 10). Profiles of fewer cells are
    #'   dominated by sampling noise, which shrinks the distance-based shift; the gsd test models that noise gene by gene and
    #'   keeps the thinner samples. Cell types tested by gsd with too few such samples get a p-value without a magnitude
    #' @param min.samp.per.level minimum samples at each compared level per cell type (default 2 with gsd, from options with the distance test)
    #' @param min.gene.frac gene filter: fraction of cells expressing the gene (default 0.01)
    #' @param genes optional gene subset
    #' @param bias.correct correct the shift for within-group dispersion (default TRUE)
    #' @param cell.number.adjust model the per-sample sampling noise in the dispersions with a 1 / (number of cells) term, so that the
    #'   dispersion change (`var`) is compared at equal cell numbers (default TRUE; otherwise a condition with fewer cells per sample
    #'   appears more dispersed)
    #' @param influence compute leave-one-sample-out effects and jackknife standard errors (default TRUE)
    #' @param seed integer seed (default from options; `NULL` draws from R's RNG)
    #' @param name results slot (default "expression.shifts")
    #' @param verbose defaults from options
    #' @param n.cores defaults from options
    #' @param ... deprecated arguments of the previous implementation (`dist.type`, `perm.method`, `pairFormula`,
    #'   `pairContrast`, `alternative`, `robust.method`, `na.mode`, `top.n.genes`, `gene.selection`, ...) are
    #'   accepted with a message and otherwise ignored
    #' @return list, also stored in `cao$test.results[[name]]`: \code{results} (one row per test x cell type x effect:
    #'   `estimate`, `estimate.norm`, `se.jk`, `ci.low`, `ci.high`, `statistic`, `p`, `padj` (BH across cell types),
    #'   `p.fwer` (max-statistic, familywise), `n`, `scheme`, `n.perm.distinct`, `p.floor`, `flags`), \code{global}
    #'   (per test x effect: p-value across all cell types), \code{wide} (per test, one row per cell type),
    #'   \code{fits}, \code{skipped} (cell types not tested, with the reason), \code{distances},
    #'   \code{adjusted.distances}, \code{influence}, \code{model}, \code{settings}, \code{n.cells}
    #' @examples
    #' \dontrun{
    #' cao$setModel(~ condition + batch, test = "condition")
    #' cao$estimateExpressionShiftMagnitudes()
    #' cao$plotExpressionShiftMagnitudes()
    #' }
    #' @param robust `"none"` (default, option `robust`), `"huber"` or `"winsor"`: down-weight samples whose residual distance is
    #'   outlying, re-estimated under every relabeling
    #' @param na.mode samples absent from a cell type: `"drop"` (default, option `na.mode`) fits the present samples; `"impute_weak"` keeps
    #'   them with a near-zero weight (same sample set in every cell type, levels stay represented)
    #' @param robust.k robust tuning constant (default: option `robust.k`, 1.345)
    #' @param statistic `"gsd"` (default: per-gene test combined over genes, see above) or `"distance"` (the distance test)
    #' @param type.gene.frac per cell type, the genes (of those passing `min.gene.frac`) detected in at least this fraction of the
    #'   type's cells, for both the gsd test and the distance model (default: option `type.gene.frac`, 0.05; see [typeGeneSets()])
    #' @param neighbourhood.gene.frac earlier name of `type.gene.frac`
    #' @param winsor.k gsd with `robust = "winsor"`: winsorization threshold in robust sd units (default 3)
    #' @param gsd.d0 gsd: prior degrees of freedom of the moderated per-gene variance (default 3); see [gsdShiftsPerCellType()]
    estimateExpressionShiftMagnitudes = function(test = NULL, formula = NULL, contrast = NULL, dispersion.formula = NULL,
                                                 statistic = c("gsd", "distance"),
                                                 dist = NULL, permutation = NULL, n.permutations = NULL, block.vars = NULL,
                                                 cell.groups = self$cell.groups, sample.per.cell = self$sample.per.cell, n.pcs = NULL,
                                                 min.cells.per.sample = NULL, min.cells.magnitude = 10, min.samp.per.level = NULL, min.gene.frac = 0.01, genes = NULL,
                                                 type.gene.frac = NULL, neighbourhood.gene.frac = NULL, winsor.k = 3, gsd.d0 = 3,
                                                 bias.correct = TRUE, cell.number.adjust = TRUE, influence = TRUE, seed = NULL, name = "expression.shifts",
                                                 verbose = NULL, n.cores = NULL, robust = NULL, na.mode = NULL, robust.k = NULL, ...) {
      statistic <- match.arg(statistic)
      min.cells.per.sample <- min.cells.per.sample %||% (if (statistic == "gsd") 5 else 10)
      # the distance model (magnitudes, var / total, noise bands, jackknife) needs well-estimated profiles: with statistic = "gsd" it
      # takes the samples with at least min.cells.magnitude cells (the gsd test models the sampling noise of thinner samples per gene)
      min.cells.dist <- if (statistic == "gsd") max(min.cells.magnitude, min.cells.per.sample) else min.cells.per.sample
      type.gene.frac <- type.gene.frac %||% neighbourhood.gene.frac %||% private$opt("type.gene.frac")   # one gene set per type: gsd and distances
      verbose <- private$opt("verbose", verbose); n.cores <- private$opt("n.cores", n.cores)
      robust <- private$opt("robust", robust); na.mode <- private$opt("na.mode", na.mode); robust.k <- private$opt("robust.k", robust.k)
      legacy <- private$legacyShiftArgs(list(...), verbose)
      permutation <- legacy$permutation %||% permutation
      dist <- private$opt("dist", dist); permutation <- private$opt("permutation", permutation)
      n.permutations <- private$opt("n.permutations", n.permutations)
      msl.gsd <- min.samp.per.level %||% 2; min.samp.per.level <- private$opt("min.samp.per.level", min.samp.per.level)   # gsd: 2 (as calibrated), distance: option
      seed <- if (missing(seed)) private$opt("seed") else seed
      alpha <- private$opt("alpha")
      if (dist == "l1") warning("dist = 'l1' does not separate location from dispersion: shift estimates absorb dispersion changes")

      model <- private$resolveModel(formula = formula, test = test, contrast = contrast, dispersion.formula = dispersion.formula,
                                    block.vars = block.vars, permutation = permutation, verbose = verbose,
                                    what = "estimateExpressionShiftMagnitudes()")
      block.vars <- block.vars %||% model$block.vars

      pb <- private$getPseudobulk(cell.groups = cell.groups, sample.per.cell = sample.per.cell, min.cells.per.sample = min.cells.per.sample,
                                  min.gene.frac = min.gene.frac, genes = genes, verbose = verbose)
      # expression shifts are estimated for contrasts (two levels, a combination of levels, or a numeric variable); whole-factor
      # tests of the model (used by the DE, composition and screen analyses) are skipped
      is.term <- vapply(model$tests, function(t) t$kind != "contrast", logical(1))
      if (all(is.term)) stop("expression shifts are estimated for contrasts; whole-factor tests (",
                             paste(vapply(model$tests, `[[`, character(1), "label"), collapse = ", "), ") are not supported: name the comparison, ",
                             "e.g. test = \"stage: II vs I\", or a combination of levels", call. = FALSE)
      if (any(is.term)) {
        if (verbose) message("Skipping whole-factor test(s) ", paste(vapply(model$tests[is.term], `[[`, character(1), "label"), collapse = ", "),
                             ": expression shifts are estimated for contrasts")
        model$tests <- model$tests[!is.term]; is.term <- is.term[!is.term]
      }
      res.gsd <- NULL
      if (statistic == "gsd") {
        if (identical(robust, "huber")) stop("robust = 'huber' is available for statistic = 'distance'; statistic = 'gsd' takes robust = 'winsor'")
        if (verbose && !identical(na.mode, "drop")) message("statistic = 'gsd' fits the samples present in each cell type (na.mode is not used)")
        if (any(!is.term)) {
          inp <- private$getGsdCellTypeInput(pb$genes)
          mc <- model; mc$tests <- model$tests[!is.term]
          if (verbose) message(sprintf("Testing %d cell type(s), %d test(s) with gsd (%d genes, per type those in >= %g %% of its cells), %d permutations...",
                                       nlevels(droplevels(as.factor(cell.groups))), length(mc$tests), length(pb$genes), 100 * type.gene.frac, n.permutations))
          res.gsd <- tryCatch(gsdShiftsPerCellType(inp$cm, sample.per.cell, cell.groups, mc, inp$lib, min.cells.per.sample = min.cells.per.sample,
                                          min.samp.per.level = msl.gsd, neighbourhood.gene.frac = type.gene.frac,
                                          permutation = permutation, n.permutations = n.permutations, block.vars = block.vars, seed = seed,
                                          robust = if (identical(robust, "winsor")) "winsor" else "none", winsor.k = winsor.k, gsd.d0 = gsd.d0, influence = FALSE, n.null = 0,
                                          n.cores = n.cores, verbose = FALSE),
                              error = function(e) {
                                if (!grepl("too few expressed genes", conditionMessage(e))) stop(e)
                                warning("statistic = 'gsd' needs more expressed genes to fit its noise model (", conditionMessage(e),
                                        "); using statistic = 'distance'", call. = FALSE)
                                NULL
                              })
          if (is.null(res.gsd)) statistic <- "distance"
        }
      }
      model.dist <- model   # the distance model gives the effect sizes of every test (and the significance of var, total and term tests)
      if (min.cells.dist != min.cells.per.sample)
        pb <- private$getPseudobulk(cell.groups = cell.groups, sample.per.cell = sample.per.cell, min.cells.per.sample = min.cells.dist,
                                    min.gene.frac = min.gene.frac, genes = genes, verbose = verbose)
      res <- NULL
      if (length(model.dist$tests)) {
        D.list <- private$getSampleDistances(pb, dist = dist, n.pcs = n.pcs, type.gene.frac = type.gene.frac, verbose = verbose)
        if (verbose) message(sprintf("Testing %d cell type(s), %d test(s), %d permutations (%s)...", length(D.list), length(model.dist$tests),
                                     n.permutations, permutation))
        res <- expressionShiftsForModel(D.list, model.dist, dist = dist, permutation = permutation, n.permutations = n.permutations,
                                        block.vars = block.vars, bias.correct = bias.correct, influence = influence,
                                        min.samp.per.level = if (statistic == "gsd") msl.gsd else min.samp.per.level, seed = seed, alpha = alpha, n.cores = n.cores,
                                        verbose = FALSE, n.cells = pb$n.cells, robust = robust, na.mode = na.mode, robust.k = robust.k,
                                        cell.number.adjust = cell.number.adjust)
      }
      res <- mergeShiftResults(res.gsd, res, model, statistic)
      if (is.null(res$distances)) res$distances <- private$getSampleDistances(pb, dist = dist, n.pcs = n.pcs, type.gene.frac = type.gene.frac, verbose = verbose)   # for the sample-space plots
      res$settings$statistic <- statistic; res$settings$dist <- dist; res$settings$alpha <- alpha; res$settings$type.gene.frac <- type.gene.frac
      res$n.cells <- pb$n.cells; res$genes.tested <- pb$genes
      res$settings$n.pcs <- n.pcs; res$settings$min.cells.per.sample <- min.cells.per.sample; res$settings$min.cells.magnitude <- min.cells.dist
      res$settings$temporary.model <- !is.null(formula) || !is.null(test) || !is.null(contrast)
      res$sample.ids <- rownames(model$meta)
      if (verbose) {
        if (nrow(res$skipped)) message(sprintf("Skipped %d cell type/test combination(s): %s", nrow(res$skipped),
                                               paste(unique(res$skipped$reason), collapse = "; ")))
        for (nt in res$notes) message("Note: ", nt)
        if (nrow(res$results)) {
          sig <- res$results[res$results$effect %in% c("shift", "location") & !is.na(res$results$padj) & res$results$padj < alpha, ]
          message(sprintf("%d of %d cell types with a significant %s at FDR %.0f%% (global p = %s).", length(unique(sig$celltype)),
                          length(unique(res$results$celltype)), if (any(res$results$effect == "shift")) "shift" else "location effect", 100 * alpha,
                          paste(format(res$global$p[res$global$effect %in% c("shift", "location")], digits = 2), collapse = ", ")))
        }
      }
      class(res) <- c("cacoaExpressionShifts", "list")
      self$test.results[[name]] <- res
      invisible(res)
    },

    #' @description Top genes of the cell-type expression shifts (gsd runs)
    #'
    #' Dot plot of the genes with the largest per-gene z-scores of the shift test in the selected cell types: genes changed at
    #' gene-level FDR `gene.fdr` as dots sized by |z| and coloured by the log2 fold change, the others as small grey dots.
    #' @param name results slot (default "expression.shifts")
    #' @param test test label or index (default: the first test)
    #' @param cell.types cell types to show, in this order (default: the significant ones by the size of the shift, at most 6)
    #' @param n.genes top genes per cell type (default 8)
    #' @param gene.fdr gene-level FDR of the coloured dots (default 0.1)
    #' @param alpha significance level for the default cell types (default from options)
    #' @return ggplot2 object
    plotExpressionShiftGenes = function(name = "expression.shifts", test = NULL, cell.types = NULL, n.genes = 8, gene.fdr = 0.1, alpha = NULL) {
      alpha <- private$opt("alpha", alpha)
      res <- private$getResults(name, "estimateExpressionShiftMagnitudes()")
      if (!length(res$genes)) stop("no per-gene results: run estimateExpressionShiftMagnitudes() with statistic = 'gsd'")
      tid <- if (is.null(test)) res$results$test.id[1] else matchModelTest(res$model, test)[1]
      tl <- res$results$test[match(tid, res$results$test.id)]
      if (is.na(tl) || is.null(res$genes[[tl]])) stop("no per-gene results for this test")
      if (is.null(cell.types)) {
        w <- res$wide[[tl]]; w <- w[is.finite(w$padj.shift) & w$padj.shift < alpha, , drop = FALSE]
        if (!nrow(w)) stop("no cell type with a significant shift: pass cell.types")
        cell.types <- utils::head(w$celltype[order(-w$shift.norm)], 6)
      }
      plotGsdShiftGenes(res$genes[[tl]], cell.types = cell.types, n.genes = n.genes, gene.fdr = gene.fdr, plot.theme = self$plot.theme)
    },

    #' @description Plot the expression-shift magnitudes per cell type
    #'
    #' For a contrast test: one row per cell type, ordered by the size of the effect. The point is the normalized effect of
    #' the distance model (`effect = "shift"`: the squared distance between the condition means relative to the summed
    #' within-group dispersions; `effect = "var"`: the dispersion change log2(s.alt / s.ref)), in the cell-type colour, filled
    #' when significant. With `statistic = "gsd"` (the default of [Cacoa] `estimateExpressionShiftMagnitudes()`) the
    #' significance of the shift is the gsd test's; the dispersion change is tested by the distance model. The grey band is the
    #' 95 % range of the effect over the test's relabelings (the size noise alone gives); the optional bar is the jackknife 95 %
    #' interval (how firmly the samples determine the size, not a test).
    #' @param name results slot (default "expression.shifts")
    #' @param test test label or index (default: the first test)
    #' @param effect `"shift"` (default) or `"var"`
    #' @param show.null show the relabeling band (default TRUE)
    #' @param show.ci show the jackknife 95 % interval (default FALSE)
    #' @param show.pvalues print the p-values: `"adjusted"` (default), `"raw"` or `"none"`
    #' @param show.genes add a margin with the number of genes changed at gene-level FDR `gene.fdr` (gsd runs; default FALSE).
    #'   The count depends on power (samples, cells): a cross-check, not a ranking
    #' @param gene.fdr FDR for `show.genes` (default 0.1)
    #' @param significance which p-value marks significance: `"padj"` (BH across cell types; default), `"p.fwer"` or `"p"`
    #' @param alpha significance level (default from options)
    #' @param order `"magnitude"` (default) or `"significance"`
    #' @param show.provenance add the provenance subtitle (default from options)
    #' @param legend.position legend position (default "bottom"; "none" hides it)
    #' @param cell.types optional subset of cell types
    #' @param ... ignored (accepted for backward compatibility: `effects`, `normalized`, `type`, `order.by`, ...)
    #' @return ggplot2 object
    plotExpressionShiftMagnitudes = function(name = "expression.shifts", test = NULL, effect = c("shift", "var"), show.null = TRUE, show.ci = FALSE,
                                             show.pvalues = c("adjusted", "raw", "none"), show.genes = FALSE, gene.fdr = 0.1,
                                             significance = c("padj", "p.fwer", "p"), alpha = NULL, order = c("magnitude", "significance"),
                                             show.provenance = NULL, legend.position = "bottom", cell.types = NULL, ...) {
      effect <- match.arg(effect); significance <- match.arg(significance); show.pvalues <- match.arg(show.pvalues); order <- match.arg(order)
      alpha <- private$opt("alpha", alpha); show.provenance <- private$opt("plot.provenance", show.provenance)
      res <- private$getResults(name, "estimateExpressionShiftMagnitudes()")
      if (is.null(res$results) || !nrow(res$results)) stop("no cell type could be tested (see cao$test.results[['", name, "']]$skipped)")
      tid <- if (is.null(test)) res$results$test.id[1] else matchModelTest(res$model, test)[1]
      if (is.na(tid) || !(tid %in% res$results$test.id)) stop("test not found; available: ", paste(unique(res$results$test), collapse = ", "))
      df <- res$results[res$results$test.id == tid, , drop = FALSE]; tl <- df$test[1]
      gsd <- identical(res$settings$statistic, "gsd") && "p.distance" %in% names(df) && any(is.finite(df$p.distance))
      subtitle <- if (show.provenance) {
        prov <- permutationPhrase(modelProvenance(res$model, tid), df)
        rob <- c(if (!identical(res$settings$robust %||% "none", "none")) sprintf("robust: %s", res$settings$robust),
                 if (identical(res$settings$na.mode, "impute_weak")) "absent samples weakly imputed")
        paste(c(prov, sprintf("%s; distance: %s%s", if (gsd && effect == "shift") "shift tested by gsd" else "tested by the distance model", res$settings$dist,
                              if (length(rob)) paste0("; ", paste(rob, collapse = "; ")) else "")), collapse = "\n")
      } else NULL
      if (!any(df$effect %in% c("shift", "var"))) stop("no shift / var effects in these results (whole-factor results of an older version?): rerun estimateExpressionShiftMagnitudes()")
      w <- res$wide[[tl]]
      counts <- NULL
      if (show.genes && length(res$genes[[tl]])) counts <- vapply(res$genes[[tl]], function(g) sum(stats::p.adjust(2 * stats::pnorm(-abs(g$z)), "BH") < gene.fdr), numeric(1))
      sig.label <- sprintf("%s%s < %.2g", if (gsd && effect == "shift") "gsd " else "", significance, alpha)
      plotShiftMagnitudes(w, effect = effect, show.null = show.null, show.ci = show.ci, show.pvalues = show.pvalues, show.genes = show.genes,
                          gene.counts = counts, gene.fdr = gene.fdr, significance = significance, alpha = alpha, sig.label = sig.label,
                          palette = self$cell.groups.palette, plot.theme = self$plot.theme, subtitle = subtitle, cell.types = cell.types,
                          order = order, legend.position = legend.position)
    },

    #' @description Leave-one-sample-out influence on the expression-shift effects (replaces `plotExpressionShiftResiduals`)
    #' @param name results slot (default "expression.shifts")
    #' @param effect `"shift"` (default), `"var"` or `"total"`
    #' @param test test label or index (default: the first contrast test)
    #' @param cell.types optional subset of cell types
    #' @param normalized express the change relative to the effect's jackknife SE (default TRUE)
    #' @return ggplot2 heatmap: samples x cell types, the change in the effect when the sample is left out
    plotSampleInfluence = function(name = "expression.shifts", effect = c("shift", "var", "total"), test = NULL, cell.types = NULL, normalized = TRUE) {
      effect <- match.arg(effect)
      res <- private$getResults(name, "estimateExpressionShiftMagnitudes()")
      if (is.null(test)) test <- which(vapply(res$model$tests, function(t) t$kind == "contrast", logical(1)))[1]
      if (is.numeric(test)) test <- names(res$influence)[test]
      infl <- res$influence[[test]]
      if (is.null(infl)) stop("no influence results for test '", test, "' (run estimateExpressionShiftMagnitudes(influence = TRUE))")
      m <- infl[[effect]]
      if (is.null(m)) stop("no '", effect, "' influence for test '", test, "' (gsd results have the shift only)")
      if (!is.null(cell.types)) m <- m[, intersect(cell.types, colnames(m)), drop = FALSE]
      w <- res$wide[[test]]
      se <- setNames(w[[paste0("se.", effect)]], w$celltype)
      groups <- self$getSampleGroups(test, model = res$model)
      plotInfluenceHeatmap(m, se = if (normalized) se[colnames(m)] else NULL, groups = groups, effect = effect, plot.theme = self$plot.theme,
                           palette = self$sample.groups.palette)
    },

    #' @description Fit diagnostics in one layout for every analysis: two heatmaps of samples by cell types. Left, how far each
    #'   sample lies from the model's fit (composition: the studentized residual of the cell type's log ratio; expression shifts:
    #'   the distance of the sample's profile from its fitted position relative to a typical sample, log2; cell density: the mean
    #'   standardized residual of the cell type's cells). Right, how much the cell type's effect changes when the sample is left
    #'   out, in standard errors of the effect. Marked tiles are beyond `resid.thresh` and `influence.se`.
    #' @param analysis `"composition"`, `"expression.shifts"` or `"cell.density"`
    #' @param name result name (default: the analysis' default slot)
    #' @param what `"both"` (default), `"residual"` or `"influence"`
    #' @param cell.types optional subset of cell types
    #' @param resid.thresh,influence.se marking thresholds (defaults: residual 3, or a sample twice as far from the fit as typical for
    #'   the expression shifts; 1 standard error)
    #' @return ggplot2 object (a two-panel grid for `what = "both"`)
    plotFitDiagnostics = function(analysis = c("composition", "expression.shifts", "cell.density"), name = NULL, what = c("both", "residual", "influence"),
                                  cell.types = NULL, resid.thresh = NULL, influence.se = 1) {
      analysis <- match.arg(analysis); what <- match.arg(what)
      name <- name %||% c(composition = "coda", expression.shifts = "expression.shifts", cell.density = "cell.density")[[analysis]]
      sub <- function(m) if (is.null(cell.types) || is.null(m)) m else m[, intersect(cell.types, colnames(m)), drop = FALSE]
      rm <- im <- NULL; r.title <- r.leg <- r.sub <- NULL
      if (analysis == "composition") {
        res <- private$getResults(name, "estimateCellLoadings()"); if (identical(res$kind, "term")) stop("the whole-factor composition test has no per-cell-type fit")
        d <- codaSampleDiagnostics(res); sm <- rownames(res$cnts); cts <- unique(d$celltype)
        mk <- function(col) matrix(d[[col]], length(sm), length(cts), dimnames = list(sm, cts))
        rm <- mk("residual"); im <- mk("influence"); groups <- self$getSampleGroups(model = res$model)
        r.title <- "Composition: residual of the log ratio"; r.leg <- "studentized residual"; resid.thresh <- resid.thresh %||% 3; i.title <- "Composition: influence on the loading"
      } else if (analysis == "expression.shifts") {
        res <- private$getResults(name, "estimateExpressionShiftMagnitudes()"); test <- which(vapply(res$model$tests, function(t) t$kind == "contrast", logical(1)))[1]
        rm <- shiftSampleResiduals(res$distances, res$model, dist = res$settings$dist %||% "cor", test = test)
        infl <- res$influence[[test]]; w <- res$wide[[test]]
        if (!is.null(infl$shift)) { se <- setNames(w$se.shift %||% w$se.jk %||% rep(NA_real_, nrow(w)), w$celltype)
          if (all(!is.finite(se))) { rr <- res$results[res$results$effect == "shift" & res$results$test.id == res$model$tests[[test]]$id, ]; se <- setNames(rr$se.jk, rr$celltype) }
          im <- sweep(infl$shift, 2, se[colnames(infl$shift)], "/") }
        rm <- rm[, colSums(is.finite(rm)) > 0, drop = FALSE]; groups <- self$getSampleGroups(test, model = res$model)
        r.title <- "Expression: distance from the fitted position"; r.leg <- "log2 distance relative to a typical sample"; resid.thresh <- resid.thresh %||% 1; i.title <- "Expression: influence on the shift"
      } else {
        res <- private$getResults(name, "estimateDiffCellDensity()"); dr <- NULL
        for (tp in names(res$diff)) if (!is.null(res$diff[[tp]]$residuals)) { dr <- res$diff[[tp]]$residuals.pearson %||% res$diff[[tp]]$residuals; break }
        if (is.null(dr)) stop("the density residuals are not stored; run estimateDiffCellDensity(..., return.residuals = TRUE)")
        cg <- self$cell.groups[colnames(dr)]; rm <- t(apply(dr, 1, function(v) tapply(v, cg, mean, na.rm = TRUE))); groups <- self$getSampleGroups()
        r.title <- "Cell density: mean residual of the cell type's cells"; r.leg <- "mean standardized residual"; resid.thresh <- resid.thresh %||% 1; i.title <- NULL
      }
      rm <- sub(rm); im <- sub(im)
      if (!is.null(rm) && !is.null(im)) { cc <- intersect(colnames(rm), colnames(im)); ss <- intersect(rownames(rm), rownames(im)); rm <- rm[ss, cc, drop = FALSE]; im <- im[ss, cc, drop = FALSE] }   # the same samples and cell types, in the same order
      so <- sampleHeatmapOrder(rm %||% im, groups)
      g.r <- if (what != "influence" && !is.null(rm)) plotSampleHeatmap(rm, groups = groups, title = r.title, legend = r.leg, mark = resid.thresh, plot.theme = self$plot.theme, palette = self$sample.groups.palette, sample.order = so)
      g.i <- if (what != "residual" && !is.null(im)) plotSampleHeatmap(im, groups = groups, title = i.title, legend = "change without the sample (standard errors)", mark = influence.se, plot.theme = self$plot.theme, palette = self$sample.groups.palette, sample.order = so)
      if (is.null(g.r) || is.null(g.i)) return(g.r %||% g.i)
      checkPackageInstalled("cowplot", cran = TRUE)
      cowplot::plot_grid(g.r, g.i + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()), ncol = 2, rel_widths = c(1.12, 1), align = "h", axis = "tb")   # the same rows in both panels: sample names once
    },

    #' @description Deprecated: use `plotSampleInfluence()`
    #' @param ... passed to `plotSampleInfluence()`
    plotExpressionShiftResiduals = function(...) {
      .Deprecated("plotSampleInfluence", msg = "plotExpressionShiftResiduals() is deprecated: residual plots of the pair model no longer exist; showing plotSampleInfluence() instead.")
      self$plotSampleInfluence(...)
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
    #' @param sample.meta data.frame Sample-level metadata (rows = samples, columns = covariates) (default=self$sample.meta)
    #' @param formula formula|character Design formula specifying covariates to model (default=NULL)
    #' @param contrast character length-3 Contrast triplet c(var, ref, target) indicating the grouping variable and the two levels to compare (default=NULL)
    #' @param name character string Slot in which to save the results (default='de')
    #' @param test character string Which DESeq2 test to use. The available options are "LRT", "Wald". (default="DESeq2.Wald")
    #' @param resampling.method character which resampling method should be used "loo" for leave-one-out or "bootstrap", (default=NULL, i.e. no resampling)
    #' @param n.resamplings numeric Number of resamplings to perform (default=30)
    #' @param seed.resampling numeric Seed to use for resamplings, input to set.seed() (default=239)
    #' @param min.cell.frac numeric Minimum fraction of cells to use to perform DE (default=0.05)
    #' @param common.genes boolean Whether to investigate common genes across cell groups (default=FALSE)
    #' @param n.cores integer Number of cores for parallelization (default=self$n.cores)
    #' @param cooks.cutoff boolean cooksCutoff for DESeq2 (default=FALSE)
    #' @param independent.filtering boolean independentFiltering parameter for DESeq2 (default=FALSE)
    #' @param min.cell.count numeric minimum number of cells that need to be present in a given cell type in a given sample in order to be taken into account (default=10)
    #' @param n.cells.subsample integer Number of cells to subsample (default=NULL)
    #' @param verbose boolean Whether to provide verbose output with diagnostic messages (default=self$verbose)
    #' @param fix.n.samples Samples to be provided if resampling.method='fix.samples'.
    #' @param genes.to.omit character Genes to omit from calculations (default = NULL)
    #' @param ... additional parameters
    #' @return A list of DE genes
    #' @examples
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' }
    #' @param sample.per.cell named factor of sample per cell (default: stored)
    #' @param test.spec what to test (see `Cacoa$new()`); builds a temporary model for this call
    estimateDEPerCellType=function(cell.groups=self$cell.groups, sample.meta=self$sample.meta, sample.per.cell=self$sample.per.cell,
                                   formula=NULL, contrast=NULL, test.spec=NULL, name='de', test='DESeq2.Wald', resampling.method=NULL, 
                                   n.resamplings=30, seed.resampling=239, min.cell.frac=0.05, common.genes=FALSE, 
                                   n.cores=self$n.cores, cooks.cutoff=FALSE, independent.filtering=FALSE, min.cell.count=10,
                                   n.cells.subsample=NULL, verbose=self$verbose, fix.n.samples=NULL, genes.to.omit = NULL, ...) {
      set.seed(seed.resampling)
      sample.model <- private$resolveModel(formula = formula, test = test.spec, contrast = contrast, verbose = verbose, what = "estimateDEPerCellType()")
      if (length(sample.model$random.vars)) {   # tests with model-based p-values cannot use a grouping term: it is fitted as a fixed term
        gf <- groupingAsFixedFormula(sample.model)
        warning(sprintf("estimateDEPerCellType(): '%s' takes its p-values from the model and cannot use a grouping term; %s %s fitted as %s with one parameter per level (model %s).%s",
                        test, paste(sprintf("'(1 | %s)'", gf$groups), collapse = ", "), if (length(gf$groups) > 1) "are" else "is", if (length(gf$groups) > 1) "fixed terms" else "a fixed term", paste(deparse(gf$formula), collapse = ""),
                        if (length(gf$dropped)) sprintf(" Covariates that are constant within its levels are absorbed by it and were dropped: %s.", paste(gf$dropped, collapse = ", ")) else ""), call. = FALSE)
        sample.model <- tryCatch(private$buildModel(formula = gf$formula, test = sample.model$test.spec$test, contrast = sample.model$test.spec$contrast),
                                 error = function(e) stop("the model with the grouping term as a fixed term cannot be fitted: ", conditionMessage(e), call. = FALSE))
      }
      if (identical(sample.model$contrast_spec$type, "term")) {   # whole-factor test: F-type tests only
        if (tolower(test) %in% c("deseq2", "deseq2.wald")) { test <- "DESeq2.LRT"; if (verbose) message("Whole-factor test: using DESeq2 LRT (full vs reduced model)") }
        if (grepl("^wilcoxon|^t-test", tolower(test))) stop("Wilcoxon / t-test need a two-group contrast")
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

      # If n.cells.subsample is set and no resampling specified -> fix.cells
      if (!is.null(n.cells.subsample) && is.null(resampling.method)) resampling.method <- 'fix.cells'

      raw.mats <- extractRawCountMatrices(self$data.object, transposed=TRUE)

      # Determine initial sample universe from model and raw.mats intersection
      all.samples <- intersect(rownames(sample.model$F), names(raw.mats))
      if (length(all.samples) == 0L) stop("No overlap between model samples and raw.mats")

      expr.fracs <- self$getJointCountMatrix() %>% getExpressionFractionPerGroup(cell.groups)
      
      gene.filter <- getPerCellTypeGeneFilter(t(self$getJointCountMatrix(raw=TRUE)), cell.groups, threshold=min.cell.frac)
      gene.filter.df <- do.call(cbind, gene.filter) 
      if (is.null(rownames(gene.filter.df))) { rownames(gene.filter.df) <- names(gene.filter[[1]])}
      gene.filter <- as.data.frame(gene.filter.df)

      if (!is.null(genes.to.omit)) {
        gene.filter[genes.to.omit, ] <- FALSE
        gene.filter %<>%
          apply(2, as.logical) %>% as.matrix() %>% `rownames<-`(gene.filter %>% rownames())
      }

      max.cell.count <- Inf
      if (!is.null(n.cells.subsample)) {
       if (verbose) message("Subsampling up to ", n.cells.subsample,
                             " cells per (cell type, core samples) in resamples.")
      }
      n.cells.subsample.core <- n.cells.subsample

      # Resampling plan = list of sample vectors (NO two-group lists)
      sample.sets <- list(initial = all.samples)

      fix.samples <- NULL
      if (!is.null(resampling.method) && (n.resamplings != 0)) {
        sample.sets <- c(sample.sets,
           prepareSamplesForDE_model(sample.model, sample.meta=sample.meta, resampling.method=resampling.method, n.resamplings=n.resamplings)
         )

        if (resampling.method == 'fix.samples') {
          if (is.null(fix.n.samples)) {
            stop("fix.n.samples must be provided for resampling.method='fix.samples'")
          }
          fix.samples <- fix.n.samples
        }
      }

      # parallelize the outer loop if subsampling is on
      n.cores.outer <- min(length(sample.sets), n.cores)
      n.cores.inner <- max(n.cores %/% length(sample.sets), 1)
      verbose.inner <- (verbose & (n.cores.outer == 1))

      de.res <- names(sample.sets) %>% sn() %>% plapply(function(resampling.name) {
        estimateDEPerCellTypeInner_model(
          raw.mats=raw.mats, cell.groups=cell.groups, samples=sample.sets[[resampling.name]],
          common.genes=common.genes, sample.meta=sample.meta, model=sample.model,  sample.per.cell = sample.per.cell,
          cooks.cutoff=cooks.cutoff, min.cell.count=min.cell.count, max.cell.count=max.cell.count,
          n.cells.subsample.core = (if (resampling.name == "initial") NULL else n.cells.subsample.core),
          independent.filtering=independent.filtering, test=test,  gene.filter=gene.filter,
          fix.n.samples=(if (resampling.name == 'initial') NULL else fix.samples), 
          seed = seed.resampling + match(resampling.name, names(sample.sets)), 
          n.cores=n.cores.inner, verbose=verbose.inner, return.matrix=is.null(resampling.method), ...)
      }, n.cores=n.cores.outer, progress=(!verbose.inner & verbose), mc.preschedule=TRUE, mc.allow.recursive=TRUE)

      # if resampling: calculate median and variance on ranks after resampling
      de.res <- if(length(de.res) > 1) summarizeDEResamplingResults(de.res) else de.res[[1]]
      de.res %<>% appendStatisticsToDE(expr.fracs[, names(de.res)])
      attr(de.res, "model") <- sample.model          # provenance (an attribute, so the per-cell-type list keeps its shape)
      attr(de.res, "test") <- test
      
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
    estimateDEStabilityPerGene = function(de.name = "de", top.n.genes = 500, p.adj.cutoff = NULL, visualize = FALSE) {
      de.res <- self$test.results[[de.name]]

      for (cell.type in names(de.res)) {
        df0 <- de.res[[cell.type]]
        subs <- attr(df0, "subsamples")

        if (is.null(df0) || !is.data.frame(df0)) next
        if (is.null(subs) || length(subs) == 0L) {
           df0$Stability <- NA_real_
           df0$n.subsamples.used <- 0L
           de.res[[cell.type]] <- df0
           next
        }

        # Decide threshold for "top" genes
        top.n <- top.n.genes
        if (!is.null(p.adj.cutoff)) {
          if (!("padj" %in% names(df0))) stop("padj column missing in initial DE table for ", cell.type)
          top.n <- sum(df0$padj <= p.adj.cutoff, na.rm = TRUE)
          top.n <- max(top.n, 1L)
        }

        genes <- rownames(df0)
        if (is.null(genes) || length(genes) == 0L) next

        hits <- setNames(numeric(length(genes)), genes)
        n.used <- 0L

        for (s in subs) {
          if (is.null(s) || !is.data.frame(s)) next
          if (!("padj" %in% names(s))) next

          common <- intersect(genes, rownames(s))
          if (length(common) == 0L) next

          r <- rank(s[common, "padj"], ties.method = "average", na.last = "keep")
          hits[common] <- hits[common] + as.numeric(r <= top.n)
          n.used <- n.used + 1L
        }

        df0$Stability <- if (n.used > 0L) hits / n.used else NA_real_
        df0$n.subsamples.used <- n.used
        de.res[[cell.type]] <- df0
      }

      self$test.results[[de.name]] <- de.res

      if (visualize) {
        stability.all <- do.call(rbind, lapply(names(de.res), function(ct) {
          d <- de.res[[ct]]
          if (is.null(d) || !is.data.frame(d) || is.null(d$Stability)) return(NULL)
          if (!("pvalue" %in% names(d))) return(NULL)

          data.frame(
            rank = rank(d$pvalue, na.last = "keep"),
            stability = d$Stability,
            cell.type = ct
          )
        }))

        p <- ggplot(stability.all, aes(rank, stability, colour = cell.type)) +
          xlim(0, top.n.genes) +
          geom_smooth(method = "loess") + self$plot.theme +
          xlab("Gene rank by p-value") + ylab("Fraction of resamples") +
          scale_color_manual(values = self$cell.groups.palette) +
          guides(color = guide_legend(override.aes = list(fill = NA), ncol = 1)) +
          theme(legend.text = element_text(size = 6),
               legend.key.height = unit(0.1, "cm"),
                legend.key.width = unit(0.2, "cm")) +
          ylim(0, 1)

        return(p)
      }
      invisible(self$test.results[[de.name]])
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
      p <- plotStability(jaccards = jaccards,notch = notch, show.jitter = show.jitter,
                         jitter.alpha = jitter.alpha, show.pairs = show.pairs, sort.order = sort.order,
                         xlabel = 'Cell Type', ylabel = 'Jaccard Index', palette=NULL, 
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
    #' @param top.n number of top genes to show
    plotDEStabilityPerGene=function(name='de', cell.type=NULL, stability.score='stab.median.rank', top.n=500) {
      de.res <- private$getResults(name, 'estimateDEPerCellType()')
      possible.scores <- c('stab.median.rank', 'stab.mean.rank', 'stab.var.rank')
      if ( !(stability.score %in% possible.scores) ) stop('Please provide correct name of the stability core')
      if (!(cell.type %in% names(de.res))) stop('Please provide correct cell type to visualise')
      df <- de.res[[cell.type]]
      if (!is.data.frame(df) || !(stability.score %in% colnames(df))) stop("Stability score not found in DE table.")
       
      if (!("pvalue" %in% colnames(df))) stop("Need 'pvalue' column in DE result.")
      idx <- rank(df$pvalue, na.last = NA) <= top.n
      if (!any(idx)) stop("No genes found in top.n filter.")
      p <- smoothScatter(x = rank(df$pvalue)[idx],
                         y = rank(df[[stability.score]])[idx],
                         xlab = 'rank of DE gene', ylab = paste0("rank of ", stability.score),
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
    #' @param name results slot name
    plotNumberOfDEGenes=function(name='de', p.adjust=TRUE, pvalue.cutoff=0.05, show.resampling.results=TRUE,
                                 show.jitter=FALSE, jitter.alpha=0.05, type='bar', notch=TRUE, ...) {
      de.raw <- private$getResults(name, 'estimateDEPerCellType()')

      if (show.resampling.results) {
        subsamples <- lapply(de.raw, function(df) attr(df, "subsamples"))
        miss.subsamples <- names(subsamples)[sapply(subsamples, function(x) is.null(x) || length(x) == 0L)]
        if (length(miss.subsamples) == length(subsamples)) {
          message("resampling results are missing for all cell types: showing point estimates ",
                  "(rerun estimateDEPerCellType() with resampling != NULL for uncertainty)")
          rl <- de.raw
        } else {
          if (length(miss.subsamples) > 0) {
            warning("Subtypes ", paste(miss.subsamples, collapse=","), " are missed form sampling, ignoring those")
          }

          # combine: resampling dfs for cell types that have them
          rl <- list()

          for (ct in names(de.raw)) {
            subs <- subsamples[[ct]]
            if (!is.null(subs) && length(subs) > 0L) {
              # name resamples for plotting
              nm <- names(subs)
              if (is.null(nm) || any(nm == "")) nm <- paste0("resample.", seq_along(subs))
              names(subs) <- paste0(ct, ".", nm)
              rl <- c(rl, subs)
            } else {
              # fallback to point estimate df
              rl[[paste0(ct, ".initial")]] <- de.raw[[ct]]
             }
          }
        }

      } else {
        rl <- de.raw
        names(rl) <- paste0(names(rl), ".initial")
      }

      # convert to dataframe for plotting
      p.col <- if (p.adjust) "padj" else "pvalue"
      df <- lapply(rl, function(d) {
        if (is.list(d) && !is.data.frame(d) && !is.null(d$res)) d <- d$res
        if (is.null(d) || !(p.col %in% colnames(d))) return(data.frame(value = NA_real_))
        data.frame(value = sum(d[[p.col]] <= pvalue.cutoff, na.rm = TRUE))
      }) %>% dplyr::bind_rows(.id = "Type")
      df$Type <- sub("\\.initial$", "", df$Type)                     # cell type labels without the resampling suffix
      df$Type <- sub("\\.resample\\.[0-9]+$", "", df$Type)

      plotMeanMedValuesPerCellType(
        df, show.jitter=show.jitter, jitter.alpha=jitter.alpha, notch=notch, type=type,
        palette=self$cell.groups.palette, ylab='number of DE genes', yline=NA, plot.theme=self$plot.theme, ...)
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
    #' @examples
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$plotVolcano()
    #' }
    plotVolcano=function(name='de', cell.types=NULL, palette=NULL, build.panel=TRUE, n.col=3,
                         color.var = 'CellFrac', ...) {
      de <- private$getResults(name, 'estimateDEPerCellType()') 
      if (is.null(palette)) {
        palette <- c("#5e4fa2", "#3288bd", "#abdda4", "#fdae61", "#f46d43", "#9e0142") %>%
          grDevices::colorRampPalette()
      }

      de <- lapply(de, getResDf)
      de <- de[!vapply(de, is.null, logical(1))]
      if (!length(de)) stop("no DE tables in cao$test.results[['", name, "']]")
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
        return(plotVolcano(de[[cell.types]], color.var=color.var, palette=palette, plot.theme=self$plot.theme, ...) + ggtitle(cell.types))
      } else if (length(cell.types) == 2) {
        n.col <- 2
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
    #' cao$saveDeAsJson()
    #' }
    #' @param sample.groups named vector of group per sample (default: derived from the stored model)
    saveDeAsJson=function(saveprefix=NULL, dir.name="JSON", de.raw=NULL, sample.groups=self$sample.groups, de.name='de',
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

      sccore::saveDeAsJson(de.raw=de.raw, saveprefix=saveprefix, dir.name=dir.name, gene.metadata=gene.metadata,
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
        } else if (color.by == 'sample') {
          params$groups <- self$sample.per.cell
        } else if(color.by %in% colnames(self$sample.meta)) {
          # map from sample values to cells
          params$groups <- setNames(self$sample.meta[as.character(self$sample.per.cell), color.by], names(self$sample.per.cell))
          if (is.null(params$show.legend)) params$show.legend <- TRUE
          if (is.null(params$legend.title)) params$legend.title <- color.by
          legend.keys <- TRUE                                             # readable keys for small, transparent points
          if (is.null(params$palette) && !is.null(names(self$sample.groups.palette)) && all(unique(stats::na.omit(as.character(params$groups))) %in% names(self$sample.groups.palette)))
            params$palette <- self$sample.groups.palette
        } else {
          stop("Unknown color.by option: '", color.by, "'. It must be 'cell.groups', 'sample', or a column in self$sample.meta",
               if (is.null(self$model)) " (no model is set: name the metadata column directly)" else "", ".")
        }
      }

      params$object <- embedding
      params$plot.theme <- plot.theme
      if (is.null(params$show.legend)) {
        params$show.legend <- !is.null(params$colors)
      }

      gg <- rlang::exec(sccore::embeddingPlot, !!!params)
      if (exists("legend.keys", inherits = FALSE) && isTRUE(params$show.legend)) gg <- gg + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
      gg
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
    #' @param p.adjust.method p-value adjustment method (default='BH')
    #' @param verbose print progress messages (default: stored option)
    #' @param n.cores number of cores (default: stored option)
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
      de.raw <- lapply(de.raw, getDeTable)
      de.raw <- de.raw[!vapply(de.raw, is.null, logical(1))]
      if (length(de.raw) == 0L) {
        stop("No usable DE tables found in de.raw.")
      }

      de.gene.ids <- getDEEntrezIdsSplitted(de.raw, org.db=org.db, p.adj=p.adj)

      go.environment <- self$getGOEnvironment(org.db, verbose=verbose, ignore.cache=ignore.cache)
      res <- estimateOntologyFromIds(
        de.gene.ids, type=type, org.db=org.db, n.top.genes=n.top.genes, go.environment=go.environment,
        pAdjustMethod=p.adjust.method, readable=readable, minGSSize=min.gs.size, maxGSSize=max.gs.size,
        keep.gene.sets=keep.gene.sets, verbose=verbose, n.cores=n.cores, ...)

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
        if (genes == "up") ont.res@result %<>% filter(NES >= 0) else if (genes == "down") ont.res@result %<>% filter(NES <= 0)
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
    #' @param name results slot name
    #' @param p.adj adjusted p-value cutoff
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
    #' @param cluster Whether to show GO clusters or raw GOs (default=TRUE)
    #' @param clust.naming Field with the results for GO clustering. Ignored if `clusters == FALSE`.
    #' @param cell.subgroups Cell groups to plot (default=NULL). This affects only visualization, but not clustering.
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
    #' @param p.adjust.method p-value adjustment method (default='BH')
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
                                plot.type=c("complete","dense","minimal"), show.ids=FALSE, string.length=14, legend.label.size=1,
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
    #' @param sep character Separator (default: tab)
    #' @param ... additional arguments passed to write.table()
    #' @return table for import into text editor
    #' @examples
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$saveOntologyAsTable(name = "GSEA", file = "Ontologies.tsv")
    #' }
    #' @param name results slot name
    #' @param p.adj adjusted p-value cutoff
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
    #' @param subtype character Type of result to filter by, must be "BP", "MF", or "CC" (default=NULL)
    #' @param genes character Direction of genes to filter by, must be "up", "down", or "all" (default=NULL)
    #' @param p.adj numeric Adjusted P to filter by (default=0.05)
    #' @param sep character Separator (default: tab)
    #' @param ... additional arguments passed to write.table()
    #' @return table for import into text editor
    #' @examples
    #' \dontrun{
    #' cao$estimateDEPerCellType()
    #' cao$estimateOntology(name = "GSEA")
    #' cao$estimateOntologyFamilies(name = "GSEA")
    #' cao$saveFamiliesAsTable(name = "GSEA", file = "Families.tsv")
    #' }
    #' @param name results slot name
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
    #' @param condition character Metadata column to group samples by (default: the variable of the stored contrast)
    #' @param palette color palette to use for conditions (default: stored $sample.groups.palette)
    #' @param show.significance boolean show statistical significance between sample groups. wilcox.test was used; (`*` < 0.05; `**` < 0.01; `***` < 0.001) (default=FALSE)
    #' @param filter.empty.cell.types boolean Remove cell types without cells (default=TRUE)
    #' @param proportions boolean Plot proportions or absolute numbers (default=TRUE)
    #' @param ... additional plot parameters, forwarded to \link{plotCountBoxplotsPerType}
    #' @return A ggplot2 object
    plotCellGroupSizes=function(cell.groups=self$cell.groups, condition=NULL, show.significance=FALSE, filter.empty.cell.types=TRUE,
                                proportions=TRUE, palette=self$sample.groups.palette, ...) {
      df.melt <- private$extractCodaData(cell.groups=cell.groups, ret.groups=FALSE)
      
      if (proportions) {
        df.melt %<>% {100 * . / rowSums(.)}
        y.lab <- "% cells per sample"
      } else {
        y.lab <- "Num. cells per sample"
      }
      
      # grouping variable: the primary test's variable
      if (is.null(condition)) condition <- private$testVariable()
      
      df.melt <- as.data.frame(df.melt)
      df.melt$group <- self$sample.meta[rownames(df.melt), condition]
      df.melt <- reshape2::melt(df.melt, id.vars="group")
      
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
    #' @param cell.groups named factor of cell type per cell (default: stored)
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

    #' @description Plot the contrast tree of the compositional result (see [contrastTree()] and [plotContrastTree()]): the
    #'   sequential binary partition of the cell types that explains the compositional shift of the contrast with the fewest
    #'   balances, each node labeled with the log2 ratio of its two sides between the contrast's endpoints, the significance
    #'   code of its permutation p-value and the share of the shift it explains; the samples' own ratios at the significant
    #'   nodes; the reference-centered loadings under the leaves.
    #' @param name name of the `estimateCellLoadings()` result (default: `"coda"`)
    #' @param alpha significance level on the adjusted node p-values (default: option `alpha`)
    #' @param split how the tree is built: `"effect"` (default) or `"evidence"`, see [contrastTree()]
    #' @param show.samples show the samples' balances at the significant nodes (default TRUE)
    #' @param leaf.values tiles under the leaves: `"reference-centered"` (default), `"clr"` or `"none"`
    #' @param palette colors of the two sides of the contrast (default: `sample.groups.palette`)
    #' @param ... passed to [plotContrastTree()] (`show.text`, `min.share`, `font.size`, `label.angle`, `legend.position`)
    #' @return ggplot2 object
    #' @examples
    #' \dontrun{
    #' cao$estimateCellLoadings()
    #' cao$plotContrastTree()
    #' }
    plotContrastTree=function(name='coda', alpha=NULL, split=c("effect", "evidence"), show.samples=TRUE, leaf.values=c("reference-centered", "clr", "none"),
                              palette=self$sample.groups.palette, ...) {
      split <- match.arg(split); leaf.values <- match.arg(leaf.values); alpha <- private$opt("alpha", alpha)
      res <- private$getResults(name, 'estimateCellLoadings()')
      if (identical(res$kind, "term")) stop("the whole-factor composition test has no contrast; the contrast tree needs a contrast result")
      tree <- if (split == "effect" && !is.null(res$tree)) res$tree else contrastTree(res, split = split, alpha = alpha)
      plotContrastTree(res, tree = tree, alpha = alpha, show.samples = show.samples, leaf.values = leaf.values, palette = palette, plot.theme = self$plot.theme, ...)
    },

    #' @description The contrast tree (the former name of `plotContrastTree()`)
    #' @param name name of the `estimateCellLoadings()` result (default: `"coda"`)
    #' @param ... passed to `plotContrastTree()`
    #' @return ggplot2 object
    plotCodaContrastTree=function(name = "coda", ...) {
      self$plotContrastTree(name = name, ...)
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
    #' @param name character Results name slot (default='coda')
    #' @param n.seed numeric Seed number for reproducibility (default=239)
    #' @param cells.to.remove character Specific cell types to keep (default=NULL)
    #' @param cells.to.remain character Specific cell types to remove (default=NULL)
    #' @param filter.empty.cell.types boolean Remove cell types without cells (default=TRUE)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param method character One of 'lda', svm', 'cda', or 'cda.std' (default=lda)
    #' @return resulting cell loadings
    #' @examples
    #' \dontrun{
    #' cao$estimateCellLoadings()
    #' }
    #' @param n.permutations number of permutations
    #' @param formula optional location formula for a temporary model (see `estimateExpressionShiftMagnitudes()`)
    #' @param contrast optional contrast for a temporary model (expert synonym of `test`)
    #' @param test optional test for a temporary model (see `Cacoa$new()`)
    #' @param perm.method relabeling scheme: `"auto"` (default: the scheme of the model's test, as for the other analyses: the tested
    #'   variable is exchanged within the levels of discrete covariates and grouping terms; residuals are exchanged, within the levels of
    #'   the discrete covariates, when the model has a continuous covariate), `"block"` or `"freedman-lane"`
    #' @param random.baseline fit a grouping term `(1 | g)` of the model as a random baseline per level of `g` (default `FALSE`: the
    #'   grouping term only restricts the relabeling to within the levels of `g`, as in every other analysis). With `TRUE` the levels
    #'   that contain one value of the tested variable contribute to the estimate, and a design in which every level contains one
    #'   value can be tested; see [lmCoda()]
    #' @param statistic test statistic of the per-column fitter: `"t"` (studentized contrast, default) or `"coef"` (raw contrast estimate)
    #' @param zero.pseudocount pseudocount added to the counts before the log-ratio transform (default 0.5 cells in a sample of median
    #'   size, scaled with the sample's cell number; see [lmCoda()])
    #' @param basis.type ILR basis type (default 'default')
    #' @param ref.method reference cell types for the presentation of the loadings: `"block"` (default) the largest set of cell
    #'   types whose proportions do not change relative to each other ([referenceBlock()]), or `"pvalue"` (the cell types with the
    #'   weakest evidence of a change relative to the geometric mean; `ref.p.thresh`, `ref.min.size`, `ref.max.size`)
    #' @param ref.block.thresh p-value below which two groups of cell types are kept apart in the block search (default 0.05)
    #' @param ref.size.frac blocks at least this fraction of the largest block's size compete for the reference, the one closest to
    #'   the center of the composition winning (default 0.5)
    #' @param ref.p.thresh p-value threshold for picking the reference cell types (`ref.method = "pvalue"`; default 0.3)
    #' @param ref.min.size minimum number of reference cell types (default 1)
    #' @param ref.max.size maximum number of reference cell types (default 3)
    #' @param block.vars metadata columns defining permutation strata (default: stored)
    #' @param ... further arguments passed to the underlying function
    estimateCellLoadings=function(n.permutations=NULL, name='coda', n.seed=239,
                                  cells.to.remove=NULL, cells.to.remain=NULL, 
                                  filter.empty.cell.types=TRUE, n.cores=self$n.cores, verbose=self$verbose, method="lda",
                                  formula=NULL, contrast=NULL, test=NULL, perm.method="auto", random.baseline=FALSE, statistic=c("t", "coef"), zero.pseudocount=0.5,
                                  basis.type = c("default"), ref.method = c("block", "pvalue"), ref.block.thresh = 0.05, ref.size.frac = 0.5, ref.p.thresh = 0.3, ref.min.size = 1, ref.max.size = 3,
                                  block.vars = NULL, ...) {
      n.permutations <- private$opt("n.permutations", n.permutations)
      # Checks
      checkPackageInstalled(c("coda.base", "psych"), cran=TRUE)
      #method <- match.arg(method)
      #if (method == "svm") checkPackageInstalled("e1071", cran=TRUE)
      #if (method == "cda") checkPackageInstalled("candisc", cran=TRUE)
      #if (method == "lda") checkPackageInstalled("quadprog", cran=TRUE)

      statistic   <- match.arg(statistic)
      basis.type  <- match.arg(basis.type)

      # Get cell counts and groups
      cnts <- private$extractCodaData(cells.to.remove=cells.to.remove, cells.to.remain=cells.to.remain,
                                     samples.to.remove=NULL)

      if (ncol(cnts) < 3) stop("Cell loadings cannot be estimated for less than 3 cell types.")

      if (filter.empty.cell.types) {
        keep.ct <- colSums(cnts, na.rm = TRUE) > 0
        cnts <- cnts[, keep.ct, drop = FALSE]
      }

      # Build model
      sample.model <- private$resolveModel(formula = formula, test = test, contrast = contrast, block.vars = block.vars, verbose = verbose,
                                           what = "estimateCellLoadings()")
      # relabeling scheme of the model's test, as for the other analyses
      if (identical(perm.method, "auto")) perm.method <- sample.model$tests[[sample.model$primary %||% 1L]]$permutation$scheme %||% "block"
      perm.method <- match.arg(perm.method, c("freedman-lane", "block"))
      if (!isTRUE(random.baseline) && !is.null(msg <- groupingBetweenMessage(sample.model, "estimateCellLoadings()"))) stop(msg, call. = FALSE)
      if (identical(sample.model$contrast_spec$type, "term")) {
        # whole-factor test on the ILR coordinates (Euclidean): location / dispersion through the distance engine
        if (verbose) message("Whole-factor test on composition: location and dispersion of the ILR coordinates")
        ilr <- computeILRMatrix(cnts, zero.pseudocount = zero.pseudocount, basis.type = basis.type)
        D <- as.matrix(stats::dist(ilr$ilr))
        tr <- testTermEffects(list(composition = D), sample.model$tests[[1]]$design, sample.model$meta, dispersion.formula = sample.model$dispersion.formula,
                              dist = "l2", permutation = private$opt("permutation"), n.permutations = n.permutations, block.vars = block.vars %||% self$block.vars,
                              seed = private$opt("seed") %||% sample.int(.Machine$integer.max, 1), alpha = private$opt("alpha"), n.cores = n.cores,
                              robust = private$opt("robust"), robust.k = private$opt("robust.k"))
        res <- list(kind = "term", results = tr$results, fits = tr$fits, ilr = ilr$ilr, psi = ilr$psi, cnts = cnts, model = sample.model, notes = tr$notes)
        self$test.results[[name]] <- res
        return(invisible(res))
      }

      if (verbose) message("Running lmCoda with design='", deparse(sample.model$formula_used), "' and ", perm.method, " permutations")

      #res <- runCoda(tmp$d.counts, tmp$d.groups, n.boot=n.boot, n.seed=n.seed, ref.cell.type=ref.cell.type, method=method, n.cores=n.cores, verbose=verbose)
      res <- lmCoda(cnts, sample.model, seed = private$opt("seed"), perm.method=perm.method, random.baseline=random.baseline, statistic=statistic, n.permutations=n.permutations,
                    zero.pseudocount=zero.pseudocount, basis.type = basis.type, ref.p.thresh = ref.p.thresh,
                    ref.min.size = ref.min.size, ref.max.size = ref.max.size, ref.method = ref.method, ref.block.thresh = ref.block.thresh, ref.size.frac = ref.size.frac, ...)
      res$cnts <- cnts
      res$model <- sample.model
      #res$groups <- tmp$d.groups

      ## Calculate normalized counts
      #ref.cell.type <- res$ref.cell.type
      #ref.cnts <- tmp$d.counts[, ref.cell.type, drop=FALSE]
      #ref.cnts[ref.cnts == 0] <- 0.5
      #norm.val <- 1 / nrow(ref.cnts) * rowSums(log(ref.cnts))
      #cnts.nonzero <- tmp$d.counts
      #cnts.nonzero[cnts.nonzero == 0] <- 0.5
      #res$norm.cnts <- log(cnts.nonzero) - norm.val

      self$test.results[[name]] <- res

      return(invisible(res))
    },

    #' @description Plot global (composition-level) effects from lmCoda
    #' @param name character Results slot name (default='coda')
    #' @param style character One of "bar" or "violin" (default="bar")
    #' @param fdr.threshold numeric FDR threshold for marking significance (default=0.05)
    #' @param drop.intercept logical Drop "(Intercept)" from coefficient list (default=TRUE)
    #' @param use.contrast.title logical Use contrast label as plot title (default=FALSE)
    #' @param ... passed to plotCodaGlobal()
    #' @return A ggplot2 object
    plotCodaGlobal = function(name = "coda", style = c("bar", "violin"), fdr.threshold = 0.05,
                              drop.intercept = TRUE, use.contrast.title = FALSE, ...) {
      style <- match.arg(style)

      x <- private$getResults(name, "estimateCellLoadings()")

      # Basic structure check
      if (is.null(x$contrast)) {
        stop("Result '", name, "' does not look like an lmCoda output. ",
             "Run estimateCellLoadings() first.")
      }

      p <- plotCodaGlobal(x, style = style, fdr.threshold = fdr.threshold, drop.intercept = drop.intercept,
                          use.contrast.title = use.contrast.title, ...)
      return(p)
    },

    #' @description Plot a heatmap of cell-type effects for contrast and coefficients (lmCoda)
    #' @param name character Results slot name (default='coda')
    #' @param coef.names character vector of coefficient names to include (default=NULL = all)
    #' @param fdr.threshold numeric FDR threshold for marking significance (default=0.05)
    #' @param ... passed to plotCodaHeatmap()
    #' @return A ggplot2 object
    plotCodaHeatmap = function(name = "coda", coef.names = NULL, fdr.threshold = 0.05, ...) {

      x <- private$getResults(name, "estimateCellLoadings()")

      if (is.null(x$contrast)) {
         stop("Result '", name, "' does not look like an lmCoda output. ",
              "Run estimateCellLoadings() first.")
      }

      plotCodaHeatmap(x, coef.names = coef.names, fdr.threshold = fdr.threshold, ...)
    },

    #' @description Plot cell-type–level effects from lmCoda (contrast or coefficient)
    #' @param name character Results slot name (default='coda')
    #' @param coef character "contrast" or a coefficient name in x$coefficients
    #' @param value character One of "predicted", "reference-centered", "volcano"
    #' @param style character "bar" or "violin" (used when value="reference-centered")
    #' @param fdr.threshold numeric FDR threshold (default=0.05)
    #' @param sortCellTypes logical Sort cell types (default=FALSE)
    #' @param labelCellTypes character Label mode for volcano: "significant","all","none"
    #' @param palette named vector of colors for cell types (default: self$cell.groups.palette)
    #' @param use.contrast.title logical Use contrast label as plot title (default=FALSE)
    #' @param significance the test behind the stars: `"clr"` (default: the primary test, relative to the geometric mean of all cell
    #'   types, reference-free) or `"reference"` (relative to the reference cell types, conditional on that choice)
    #' @param ... passed to plotCodaCellTypes() (`show.ci`: error bars of the predicted proportions, default TRUE)
    #' @return A ggplot2 object
    plotCodaCellTypes = function(name = "coda", coef = "contrast", style = c("bar", "violin"),
                                 value = c("predicted", "reference-centered", "volcano"), 
                                 fdr.threshold = 0.05, sortCellTypes = FALSE,
                                 labelCellTypes = c("significant", "all", "none"),
                                 palette = self$cell.groups.palette, use.contrast.title = FALSE, significance = c("clr", "reference"), ...) {

      value <- match.arg(value); significance <- match.arg(significance)
      style <- match.arg(style)
      labelCellTypes <- match.arg(labelCellTypes)

      x <- private$getResults(name, "estimateCellLoadings()")

      if (is.null(x$contrast) || is.null(x$psi)) {
        stop("Result '", name, "' does not look like an lmCoda output. ",
            "Run estimateCellLoadings() first.")
      }

      # Ensure palette is either NULL or named with celltypes
      celltype.palette <- palette
      if (!is.null(celltype.palette) && is.null(names(celltype.palette))) {
       # unnamed palettes are ambiguous; let plotCodaCellTypes generate stable colors
        celltype.palette <- NULL
      }

      p <- plotCodaCellTypes(x, coef = coef, value = value, style = style, fdr.threshold = fdr.threshold,
                        sortCellTypes = sortCellTypes, labelCellTypes = labelCellTypes, 
                        celltype.palette = celltype.palette, use.contrast.title = use.contrast.title, significance = significance, ...)
      return(p)
    },




    #' @description Plot the compositional loadings of the contrast per cell type (see [plotCodaLoadings()])
    #' @param alpha significance threshold on the adjusted per-cell-type p-values (default: option `alpha`)
    #' @param palette named colours per cell type (default: `cell.groups.palette`)
    #' @param name name of the `estimateCellLoadings()` result (default: `"coda"`)
    #' @param ordering `"loading"` (default) or `"pvalue"`
    #' @param show.pvals print the adjusted p-values next to the points (default TRUE)
    #' @param show.null show the 95 % permutation range of each loading as a grey band (default TRUE)
    #' @param significance the test behind the filled points and p-values: `"clr"` (default: the primary test, relative to the
    #'   geometric mean of all cell types, reference-free) or `"reference"` (relative to the reference cell types, conditional on
    #'   that data-chosen set); see [plotCodaLoadings()]
    #' @param center `"reference"` (default): loadings centred on the reference cell types; `"none"`: as they are (zero = the
    #'   geometric mean of all cell types); `"both"`: the two side by side (cowplot)
    #' @param ... passed to [plotCodaLoadings()] (e.g. `legend.position`)
    #' @return ggplot2 object (`center = "both"`: a cowplot grid)
    #' @examples
    #' \dontrun{
    #' cao$estimateCellLoadings()
    #' cao$plotCellLoadings()
    #' cao$plotCellLoadings(center = "both")
    #' }
    plotCellLoadings=function(alpha=NULL, palette=self$cell.groups.palette, name='coda', ordering=c('loading', 'pvalue'),
                              show.pvals=TRUE, show.null=TRUE, significance=c("clr", "reference"), center=c("reference", "none", "both"), ...) {
      ordering <- match.arg(ordering); significance <- match.arg(significance); center <- match.arg(center); alpha <- private$opt("alpha", alpha)
      res <- private$getResults(name, 'estimateCellLoadings()')
      if (identical(res$kind, "term")) {
        stop("the whole-factor composition test has no per-cell-type loadings; see cao$test.results[['", name, "']]$results ",
             "(location / dispersion R2 and p-values)")
      }
      lab <- res$contrast$predicted$labels
      one <- function(ctr, title = NULL) plotCodaLoadings(res, ref.level = self$ref.level %||% lab$baseline %||% "reference", target.level = self$target.level %||% lab$target %||% "target",
                       alpha = alpha, palette = palette, ordering = ordering, show.pvals = show.pvals, show.null = show.null, significance = significance, center = ctr, title = title, plot.theme = self$plot.theme, ...)
      if (center != "both") return(one(center))
      checkPackageInstalled("cowplot", cran = TRUE)
      cowplot::plot_grid(one("none", "relative to the geometric mean"), one("reference", "relative to the reference cell types"), ncol = 2)
    },

    #' @description The samples behind the composition result: per cell type the samples' proportions by level of the tested variable,
    #'   the fitted proportions with their intervals, and the samples that carry the result named (see [plotCodaSamples()])
    #' @param name result name (default 'coda')
    #' @param cell.types cell types to show (default: the `n.top` with the smallest p-values)
    #' @param n.top number of cell types shown by default
    #' @param pair.by optional column of the sample metadata joining samples of the same unit (e.g. patient)
    #' @param palette colours of the levels of the tested variable
    #' @param ... passed to [plotCodaSamples()] (`influence.se`, `resid.thresh`, `label.samples`, `ncol`)
    #' @return A ggplot2 object
    plotCodaSamples = function(name = "coda", cell.types = NULL, n.top = 12, pair.by = NULL, palette = self$sample.groups.palette, ...) {
      res <- private$getResults(name, "estimateCellLoadings()")
      if (identical(res$kind, "term")) stop("the whole-factor composition test has no per-cell-type fit to show")
      if (!is.null(pair.by) && !pair.by %in% names(res$model$meta) && pair.by %in% names(self$sample.meta)) res$model$meta[[pair.by]] <- self$sample.meta[rownames(res$model$meta), pair.by]
      plotCodaSamples(res, cell.types = cell.types, n.top = n.top, pair.by = pair.by, palette = palette, plot.theme = self$plot.theme, ...)
    },

    ### Cluster-free cell density

    #' @description Diagnostics of the kNN neighbourhood units used by the cluster-free maps: cells per sample per
    #'   unit (power), unit sizes, and purity against the cell annotation (resolution), for a size `k`; no test is run.
    #' @param k neighbourhood size (default: [defaultNeighbourhoodSize()] for this object)
    #' @param min.cells cells per sample below which a sample is dropped from a unit's test (default 3)
    #' @param verbose print the one-sentence summary
    #' @return invisibly, the list of [neighbourhoodDiagnostics()]
    #' @examples
    #' \dontrun{
    #' cao$checkNeighbourhoods()
    #' }
    checkNeighbourhoods = function(k = NULL, min.cells = 3, verbose = self$verbose) {
      k <- k %||% defaultNeighbourhoodSize(length(self$sample.per.cell), length(unique(self$sample.per.cell)))
      d <- neighbourhoodDiagnostics(private$knnUnits(k, self$n.cores), self$sample.per.cell, self$cell.groups, min.cells = min.cells, k = k)
      if (verbose) message(d$text)
      invisible(d)
    },

    #' @description Estimate cell density in giving embedding
    #' @param bins numeric Number of bins for density estimation (default=400)
    #' @param method character string Density estimation method, graph: graph smooth based density estimation. kde: embedding grid based density  estimation. (default: 'kde')
    #' @param beta numeric Smoothing strength parameter of the \link[sccore:heatFilter]{heatFilter} for graph based cell density (default=30)
    #' @param estimate.variation boolean Estimate variation (default=TRUE)
    #' @param contrast Linear contrast specification (required).
    #' @param formula RHS formula for the design (~ ...). If NULL, uses default.
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
    #' @param name results slot name
    #' @param test optional test for a temporary model (see `Cacoa$new()`)
    #' @param sample.metadata sample metadata (default: stored)
    #' @param block.vars metadata columns defining permutation strata (default: stored)
    #' @param smoothing graph method only: `"knn"` (default): per-sample counts in each cell's kNN neighbourhood of
    #'   `k` cells (the unit with the best power at matched error in the benchmark, `misc/benchmark/README.md`);
    #'   `"heat"`: the heat-kernel smooth, width set by `beta` or by `target.size`; `"neighbourhood"`: counts within
    #'   `n.hops` graph hops. The count units are tested on the log scale by default (`estimateDiffCellDensity(transform)`)
    #' @param k neighbourhood size for `smoothing = "knn"` (default: the larger of 100 and 10 cells per sample, capped
    #'   at 1 % of the cells; [defaultNeighbourhoodSize()])
    #' @param target.size heat kernel: median effective neighbourhood size in cells, from which `beta` is solved
    #'   ([solveKernelBeta()]); the same `beta` is 20x wider on one graph than another, so a width in cells is the
    #'   comparable setting. The effective size per cell is stored as `unit.size`
    #' @param n.hops neighbourhood radius for `smoothing = "neighbourhood"` (default 1)
    estimateCellDensity = function(bins=400, method='kde', name='cell.density', beta=30, estimate.variation=TRUE, contrast=NULL,
                                   formula=NULL, test=NULL, verbose=self$verbose, n.cores=self$n.cores, sample.metadata=self$sample.meta,
                                   block.vars=NULL, bandwidth=0.05, smoothing=c("knn", "heat", "neighbourhood"), k=NULL, target.size=NULL, n.hops=1, ...){
      smoothing <- match.arg(smoothing)
      sample.per.cell <- self$sample.per.cell
      if (anyNA(sample.per.cell)) { sample.per.cell <- sample.per.cell[!is.na(sample.per.cell)]; if (is.factor(sample.per.cell)) sample.per.cell <- droplevels(sample.per.cell) }   # cells of excluded samples

      sample.model <- private$resolveModel(formula = formula, test = test, contrast = contrast, block.vars = block.vars, verbose = verbose,
                                           what = "estimateCellDensity()")
      if (identical(sample.model$contrast_spec$type, "term")) stop("cell density needs a contrast test (two groups or a numeric step), not a whole-factor test")

      # Get sample weights
      # regression weights w = X (X'X)^- c: +1/n_alt and -1/n_ref for a plain two-group comparison, adjusted weights with covariates
      sample.weights <- regressionWeights(sample.model)

      # Ensure weights are named by sample and restrict to samples present in data
      sample.weights <- sample.weights[intersect(names(sample.weights), unique(sample.per.cell))]

      if (length(sample.weights) == 0L) {
       stop("No sample weights overlap with sample.per.cell; check sample IDs.")
      }

      if (method == 'kde') {
        private$checkCellEmbedding()
        res <- self$embedding %>% estimateCellDensityKde(
          sample.per.cell=sample.per.cell, sample.weights=sample.weights, bins=bins, bandwidth=bandwidth, ...
        )
      } else if (method == 'graph') {
        res <-  private$cellGraph() %>%
          estimateCellDensityGraph(sample.per.cell=sample.per.cell, sample.weights=sample.weights,
                                   n.cores=n.cores, beta=beta, verbose=verbose, smoothing=smoothing, k=k, target.size=target.size, n.hops=n.hops,
                                   neighbourhoods=if (smoothing == "knn") private$knnUnits(k %||% defaultNeighbourhoodSize(length(sample.per.cell), length(unique(sample.per.cell))), n.cores) else NULL, ...)
      } else stop("Unknown method: ", method)
      

      if (estimate.variation) {
        # groups defined by sign of weights
        w <- sample.weights[colnames(res$density.mat)]
        sg.ids <- list( num = names(w)[w >  0], den = names(w)[w <  0])

        res$density.mad <- res %$% lapply(sg.ids, function(ids) apply(density.mat[,ids,drop=FALSE], 1, mad)) %>%
          Reduce(`+`, .)
        res$density.sd <-res %$% lapply(sg.ids, function(ids) apply(density.mat[,ids,drop=FALSE], 1, sd)) %>%
          Reduce(`+`, .)
        res$missed.sample.frac <- res$density.mat %>% {(. / rowMeans(.)) < 0.05} %>% rowMeans()
      }

      res$sample.weights  <- sample.weights
      res$contrast_spec   <- sample.model$contrast_spec
      res$model           <- sample.model

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

      color.range %<>% parseLimitRange(unlist(dens.res$cond.densities))

      ps <- names(dens.res$cond.densities) %>% sn() %>% lapply(function(l) {
        #color.range %<>% parseLimitRange(unlist(dens.res$cond.densities))
        if (dens.res$method =='graph') {
          col.vec <- dens.res$cond.densities[[l]]

          # Align to embedding cells if names are present
          if (!is.null(names(col.vec))) {
           col.vec <- col.vec[rownames(self$embedding)]
          }
          p <- self$plotEmbedding(colors=col.vec, size=size, title=l, show.legend=show.legend,
                                  legend.title='Density', color.range=color.range, palette=palette, ...)
        } else {# dens.res$method =='kde'
          p <- data.frame(density.emb, z=dens.res$cond.densities[[l]]) %>%
            plotDensityKde(bins=dens.res$bins, lims=color.range, title=l, show.legend=show.legend,
                           show.grid=show.grid, plot.theme=self$plot.theme, palette=palette, ...)

          # ------------------- overlay points by weights -------------------
        if (add.points && !is.null(dens.res$sample.weights)) {
          w <- dens.res$sample.weights

        # Samples belonging to numerator / denominator side
          samps <- switch(l,num = names(w)[w >  0],den = names(w)[w <  0], character(0))
          
          if (length(samps) > 0L) {
            sample.per.cell <- self$sample.per.cell
            # cells whose sample is in the chosen side
            cells <- names(sample.per.cell)[sample.per.cell %in% samps]

            if (length(cells) > 0L) {
              nnames <- sample(cells, min(2000, length(cells)))
              emb <- as.data.frame(self$embedding) %>%
              set_colnames(c('x', 'y')) %>%
              cbind(z = 1)

              p <- p + ggplot2::geom_point(data = emb[nnames, ],
                ggplot2::aes(x = x, y = y), col = point.col, size  = 0.00001, alpha = 0.2)
            }
           }
          }
          # -----------------------------------------------------------------
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
    #' @param adjust boolean Whether to adjust Z-scores for multiple comparison using BH method (default: FALSE for type='subtract', TRUE for everything else)
    #' @param name character Slot with results from estimateCellDensity. New results will be appended there. (Default: 'cell.density')
    #' @param sample.metadata data.frame Sample metadata (default=self$sample.meta)
    #' @param formula formula Model formula (default=NULL)
    #' @param contrast character Vector with two elements specifying the reference and target levels (default=NULL)
    #' @param n.permutations numeric Number of permutations (default=400)
    #' @param smooth boolean Smooth results (default=TRUE)
    #' @param verbose boolean Print messages (default=self$verbose)
    #' @param n.cores integer Number of cores to use for parallelization (default=self$n.cores)
    #' @param ... additional arguments to the function
    #' @return estimated differential cell densities: `pooled` (unadjusted signed z-scores), `adj` (adjusted), `p.adj`,
    #'   `effect` (the contrast estimate per cell: the density difference, or the log2 fold change with
    #'   `transform = "log"`), `p.perm`, `p.pooled`, `p.param`, `t`, `group.test` (when cell groups are stored);
    #'   with `scale`: `scale`, `scale.mode` and `scales` (the per-scale maps, see `scale`)
    #' @examples
    #' \dontrun{
    #' cao$estimateCellDensity()
    #' cao$estimateDiffCellDensity()
    #' }
    #' @param test optional test for a temporary model (see `Cacoa$new()`)
    #' @param block.vars metadata columns defining permutation strata (default: stored)
    #' @param perm.method relabeling scheme: `"auto"` (default: the scheme of the model's test, i.e. sample labels relabeled within the
    #'   levels of the discrete covariates and grouping variables (`"block"`) when every adjusting covariate is discrete, residual
    #'   permutation (`"freedman-lane"`) when a continuous covariate is present), `"block"` or `"freedman-lane"`. With a covariate of many
    #'   levels (patients) residual permutation across all samples is not valid for the density map: in data without a label effect it
    #'   called cells in 67-78 % of the data sets, relabeling within the levels in 0-5 % (residuals are now exchanged within the levels
    #'   of the discrete covariates, see [permutationPlan()])
    #' @param statistic test statistic of the fitter: `"t"` (studentized contrast, default) or `"coef"` (raw contrast estimate)
    #' @param robust.method robust fitting: 'none', 'huber' or 'winsor'
    #' @param na.mode handling of missing responses: 'drop' or 'impute_weak'
    #' @param alternative alternative hypothesis for the z-scores: 'two-sided', 'greater' or 'less'. The unadjusted z-scores
    #'   (`pooled`) are the signed tail z of the per-cell p-values of the pooled permutation null (`p.pooled`, see
    #'   [pooledPermutationP()]; calibrated and not floored), two-sided: `sign(effect) * qnorm(1 - p / 2)`; `adj` holds
    #'   the signed tail z of the adjusted p-values
    #' @param return.residuals also return residuals
    #' @param return.sampled.stats also return the permutations x cells matrix of the permuted statistics (default FALSE:
    #'   the pooled quantities are accumulated in C++ without it; max-t keeps it regardless)
    #' @param adjust.method multiplicity adjustment across cells, see [adjustCellScores()]: `"perm-fdr"` (default: the
    #'   pooled permutation FDR of the studentized contrast over cells, see [permutationFDR()]; no per-cell floor, uses
    #'   the permutation blocks of the design), `"bh"` (Benjamini-Hochberg over the tested cells' p-values; the
    #'   parametric p-values ignore permutation blocks and are conservative where few samples carry density),
    #'   `"spatial-fdr"` (Milo's spatially weighted BH), `"neighbourhood"` (a covering set of index neighbourhoods is
    #'   tested; liberal) or `"max-t"` (the maximum statistic under permutations; the former default). Benchmarks:
    #'   `misc/benchmark/README.md`
    #' @param p.source per-cell p-values for the BH-type alternatives: `"parametric"` (t statistic of the contrast, default),
    #'   `"permutation"` (floored at 1 / (n.permutations + 1)) or `"pooled"` (the cell's standardized statistic against
    #'   the pooled permutation null of all cells, see [pooledPermutationP()]: no floor, calibrated with enough
    #'   permutations); `p.pooled` is returned whenever the permuted statistics are kept
    #' @param moderate shrink the per-cell residual variances towards a prior across cells (a scaled-F prior with a trend on the
    #'   cell's mean density, [fitFPrior()]; a moderated t with extra degrees of freedom), in the permutation statistic of the
    #'   diffusion neighborhoods and before the parametric p-values of the stored units (default TRUE)
    #' @param moderate.d0.max cap on the prior degrees of freedom of the moderation (default 4; see [fitFPrior()]): the cell's own
    #'   variance keeps a share of at least `df / (df + moderate.d0.max)` of the moderated variance. Without the cap a homogeneous set of
    #'   variances gives an infinite prior df, which replaces every cell's variance by the trend and costs the cells with a real
    #'   difference (whose within-group variance is below the trend) part of their statistic
    #' @param transform `"auto"` (default): `"log"` for the graph units (`smoothing = "knn"`, `"neighbourhood"` or
    #'   `"heat"`), `"none"` for the kde; `"none"`: test the densities as they are; `"log"`: for the count units `log2`
    #'   counts per 10,000 cells of the sample (`log2((density * n + 0.5) / (n + 1) * 1e4)`, `n` the sample's cell
    #'   number), the log-count model of a neighbourhood abundance test; for the heat kernel `log2(d + kernel.alpha /
    #'   n.cells)` of the per-sample densities `d` per cell of the sample, the cell itself counted as one cell of its
    #'   neighborhood (a count pseudo-count would swamp a kernel density; [heatLogDensity()], [ownCellDensity()]); the
    #'   effect is then a log2 fold change of the density
    #' @param pooled.statistic the per-cell statistic pooled over cells for the pooled p-values and FDR: `"raw"`
    #'   (default; the studentized contrast t, as benchmarked) or `"standardized"` (t minus its permutation mean over
    #'   its permutation sd, the shift map's convention)
    #' @param post post-processing of the per-cell test on the cell graph (graph density only): `"averaged"`
    #'   (default): the statistic is averaged over each cell and its graph neighbours, under every permutation alike,
    #'   before the pooled FDR, an exact test of the one-hop-larger unit ([averageOverNeighbours()]; +10-40 % recall
    #'   at matched error in the benchmark, q-values stay calibrated); `"inherit"`: every cell takes the smallest
    #'   q among itself and its neighbours after the test ([inheritNeighbourQ()]; Milo's rule, more reach but
    #'   the realized false-call rate at a nominal threshold is about the nominal rate); `"none"`. The unadjusted
    #'   z-scores, effects and q-values of the result are those of the averaged statistic under `"averaged"`
    #' @param scale the neighborhoods of the test. `NULL`: the units of the stored density (`estimateCellDensity()`, the
    #'   kNN count units by default). Otherwise diffusion neighborhoods on the cell graph, whose width is the share of
    #'   the cells a typical neighborhood covers (the kernel strength is solved per graph so that the median cell's
    #'   effective neighborhood, [kernelEffectiveSize()], holds that share; [solveKernelBeta()]; the diffusion follows
    #'   the graph, so a neighborhood stays mostly within its own population). `"auto"` (default): the scale scan over
    #'   `auto.scales` chooses the scale whose own map calls the most cells at q < 0.1 (the finest of the scales calling at
    #'   least 90 % of the cells when there are such; the scale of its strongest 1 % of the cells when none calls any), and
    #'   the result is the single-scale map at that scale, so that its q-values and
    #'   effects belong to one scale. One value (`scale = 0.1`) tests a single scale; several values (`scale = c(0.01,
    #'   0.03, 0.1)`) run the scale scan itself ([scanDensityScales()]): every scale is tested under the same relabelings,
    #'   a cell's statistic is its smallest p over the scales and its q-value the pooled permutation false discovery rate
    #'   of that statistic (against the same quantity under every relabeling, no Bonferroni factor), so that each cell
    #'   may be called at its own scale. Each map is the same test on `log2(d + kernel.alpha / n.cells)` of the per-sample
    #'   kernel densities ([heatLogDensity()]; at each cell the heat arriving from every sample's cells, the cell itself
    #'   counted as one cell of its neighborhood, so that a weakly connected cell, whose own heat stays on it, is tested
    #'   through the neighborhood its arriving heat comes from) with the trended moderation of the statistic
    #'   (`moderate`; needs limma).
    #'   The scan needs `adjust.method = "perm-fdr"`. The result records `scale` (the scale(s) used), `scale.mode`
    #'   (`"single"`, `"scan"` or `"auto"`) and `scales` (per scale: strength, median size, and per cell the pooled p,
    #'   the effect, the single-scale q and the best scale; under `"auto"` also the scan's q and p, `scan.q`, `scan.p`)
    #' @param auto.scales the scales `scale = "auto"` chooses from (default 1, 3 and 5 % of the cells, the expression-shift map's set;
    #'   10 % as the coarsest scale gives this map 0.007 more ranking recall on the benchmark at a higher false share on PF; ignored otherwise)
    #' @param effect.scale under the scan, the scale whose log2 fold change is reported as `effect`: `"finest"` (default;
    #'   the smallest scale, the least attenuated estimate, noisy only where the scan is not significant) or `"best"`
    #'   (each cell's scale of the smallest p)
    #' @param kernel.alpha pseudo-count of a heat-kernel density on the log scale, as a fraction of the mean density
    #'   (default 0.3; the diffusion neighborhoods, and `smoothing = "heat"` with `transform = "log"`)
    estimateDiffCellDensity=function(type='permutation', adjust=NULL, name='cell.density', sample.metadata=self$sample.meta, 
                                     formula=NULL, contrast=NULL, test=NULL, block.vars=NULL, perm.method="auto", statistic=c("t", "coef"),
                                     robust.method="none", na.mode="drop", alternative="two-sided", return.residuals=FALSE, return.sampled.stats=FALSE,
                                     n.permutations=NULL, smooth=TRUE, adjust.method=c("perm-fdr", "bh", "spatial-fdr", "neighbourhood", "max-t"),
                                     p.source=c("parametric", "permutation", "pooled"), moderate=TRUE, moderate.d0.max=4, transform=c("auto", "none", "log"),
                                     post=c("averaged", "none", "inherit"), pooled.statistic=c("raw", "standardized"), scale="auto", auto.scales=c(0.01, 0.03, 0.05), effect.scale=c("finest", "best"), kernel.alpha=0.3,
                                     verbose=self$verbose, n.cores=self$n.cores, ...){
      n.permutations <- private$opt("n.permutations", n.permutations)
      dens.res <- private$getResults(name, 'estimateCellDensity')
      effect.scale <- match.arg(effect.scale); scale.given <- !missing(scale)
      statistic <- match.arg(statistic); adjust.method <- match.arg(adjust.method); p.source <- match.arg(p.source); transform <- match.arg(transform); post <- match.arg(post); pooled.statistic <- match.arg(pooled.statistic)
      if (is.null(adjust)) adjust <- (type != 'subtract') # NULL can be forwarded here
      if (transform == "auto") transform <- if (isTRUE(dens.res$smoothing %in% c("knn", "neighbourhood", "heat"))) "log" else "none"
      if (post != "none" && dens.res$method != "graph") { if (verbose) message("post-processing on the cell graph applies to the graph density only: post = 'none'"); post <- "none" }
      if (post != "none" && !adjust) post <- "none"
      density.mat <- dens.res$density.mat
      if (transform == "log" && isTRUE(dens.res$smoothing == "heat")) {   # a kernel density, per cell of the sample (estimateCellDensityGraph()): relative pseudo-count (heatLogDensity())
        density.mat[!is.finite(density.mat) | density.mat < 0] <- 0
        density.mat <- log2(density.mat + kernel.alpha / nrow(density.mat))
      } else if (transform == "log") {
        sizes <- dens.res$sample.sizes %||% table(self$sample.per.cell)
        sizes <- as.numeric(sizes[colnames(density.mat)])
        if (any(is.na(sizes))) stop("sample sizes are missing for some density columns")
        density.mat <- log2(t((t(density.mat) * sizes + 0.5) / (sizes + 1)) * 1e4)
      }
      if (dens.res$method == 'kde'){
        if (adjust) {
          # For p-value adjustment matrix filtration can result in disconnected grid and break estimating l.max
          l.max <- rownames(density.mat) %>% as.integer() %>% graphFromGrid(n.bins=dens.res$bins) %>%
            igraph::laplacian_matrix(sparse=TRUE) %>% irlba::partial_eigen(n=1) %>% .$values
        }

        density.mat <- density.mat[dens.res$density.emb$counts > 0,]
        graph <- rownames(density.mat) %>% as.integer() %>% graphFromGrid(n.bins=dens.res$bins)
      } else if (adjust && smooth) {
        graph <- private$cellGraph()
        l.max <- NULL
      }

      sample.model <- if (!is.null(formula) || !is.null(contrast) || !is.null(test)) {
        private$resolveModel(formula = formula, test = test, contrast = contrast, block.vars = block.vars, verbose = verbose, what = "estimateDiffCellDensity()")
      } else if (!is.null(self$model)) {   # the current model, as in every other analysis (the density estimate keeps the model it was made under)
        cur <- private$resolveModel(verbose = verbose, what = "estimateDiffCellDensity()"); old <- dens.res$model
        sig <- function(m) paste(paste(deparse(m$formula.full %||% m$formula), collapse = ""), paste(vapply(m$tests, `[[`, character(1), "label"), collapse = ";"), paste(m$block.vars, collapse = ","))
        if (verbose && !is.null(old) && !identical(sig(cur), sig(old)))
          message("The model has changed since estimateCellDensity() (", paste(deparse(old$formula.full %||% old$formula), collapse = ""), "): the test uses the current model, ",
                  paste(deparse(cur$formula.full %||% cur$formula), collapse = ""), ". Re-run estimateCellDensity() to refresh the condition summaries of the density estimate.")
        cur
      } else dens.res$model %||% private$resolveModel(verbose = verbose, what = "estimateDiffCellDensity()")
      # relabeling scheme of the model's test: within the levels of discrete covariates / grouping variables when all adjusting covariates are discrete
      if (identical(perm.method, "auto")) perm.method <- sample.model$tests[[sample.model$primary %||% 1L]]$permutation$scheme %||% "freedman-lane"
      perm.method <- match.arg(perm.method, c("freedman-lane", "block"))

      if (!is.null(scale)) {   # diffusion neighborhoods at the given scale(s): the scale scan (scanDensityScales()); the stored units when the request does not fit the scan and scale was not given
        why <- if (!(adjust && adjust.method == "perm-fdr")) "the diffusion neighborhoods are adjusted by the pooled permutation FDR (adjust = TRUE, adjust.method = 'perm-fdr')" else if (dens.res$method != "graph") "the diffusion neighborhoods need the graph density (estimateCellDensity(method = 'graph'))" else NULL
        if (!is.null(why) && scale.given) stop(why, "; scale = NULL for the stored units")
        if (!is.null(why)) { if (verbose) message(why, ": the stored units"); scale <- NULL }
      }
      if (!is.null(scale)) {
        auto <- identical(scale, "auto"); if (auto) scale <- checkAutoScales(auto.scales)
        scale <- resolveScales(scale, length(self$sample.per.cell), floor.cells = 10, verbose = verbose); K <- length(scale)
        g.s <- extractCellGraph(self$data.object); spc.s <- self$sample.per.cell
        g.s <- igraph::induced_subgraph(g.s, intersect(igraph::V(g.s)$name, names(spc.s)[!is.na(spc.s)]))
        kb <- lapply(scale, function(sc) private$kernelBetaForScale(g.s, sc, n.cores))
        mp <- if (moderate) list(moderate.perm = "trend", moderate.d0.max = moderate.d0.max) else list()
        if (verbose) message(sprintf("diffusion neighborhoods at %s of the cells%s", paste(sprintf("%g %%", 100 * scale), collapse = ", "), if (auto) " (auto: one scale is chosen)" else if (K > 1) " (scale scan)" else ""))
        sc <- do.call(scanDensityScales, c(list(g.s, spc.s, sample.model, scale = scale, beta = sapply(kb, `[[`, "beta"), self = lapply(kb, `[[`, "self"), off.size = lapply(kb, `[[`, "off.size"), alpha = kernel.alpha, post = if (post == "none") "none" else "averaged",
                                                alternative = alternative, n.permutations = n.permutations, seed = private$opt("seed"), groups = self$cell.groups, n.cores = n.cores, verbose = verbose,
                                                perm.method = perm.method, robust.method = robust.method, na.mode = na.mode), mp))
        cells <- sc$cells; msize <- sapply(kb, `[[`, "median.size")
        if (auto || K == 1) {   # one scale: its own map (auto: the scale calling the most cells at 0.1, the finest when a scale calls nearly all; see chooseScanScale)
          k <- if (K == 1) 1L else chooseScanScale(sc)
          q <- sc$Q[, k]; p <- sc$P[, k]; effect <- sc$E[, k]; scale.used <- scale[k]; gt <- sc$group.test[[k]]
        } else { q <- sc$q; p <- sc$p; effect <- if (effect.scale == "finest") sc$E[, 1L] else sc$E[cbind(seq_along(cells), sc$best)]; scale.used <- scale; gt <- sc$group.test[[1L]] }
        names(q) <- names(p) <- names(effect) <- cells
        tz <- function(pv, sg) { pa <- pmin(pmax(pv, 1e-300), 1 - 1e-16)
          z <- switch(alternative, "two-sided" = sg * stats::qnorm(pa / 2, lower.tail = FALSE), greater = pmax(stats::qnorm(pa, lower.tail = FALSE), 0), less = pmin(stats::qnorm(pa), 0))
          z[!is.finite(pv)] <- NA_real_; stats::setNames(z, cells) }
        sg <- sign(effect); sg[!is.finite(sg) | sg == 0] <- 1
        res <- list(pooled = tz(p, sg), adj = tz(q, sg), p.adj = q, p.pooled = p, effect = effect, model = sample.model, perm.method = perm.method, statistic = statistic, robust = robust.method,
                    adjust.method = "perm-fdr", p.source = "permutation", post = post, transform = "log", unit = "diffusion", scale = scale.used,
                    scale.mode = if (auto) "auto" else if (K == 1) "single" else "scan", effect.scale = if (!auto && K > 1) effect.scale else NULL,
                    scales = list(scale = scale, beta = sc$scales$beta, median.size = msize, best = stats::setNames(scale[sc$best], cells), p = sc$P, q = sc$Q, effect = sc$E, scan.q = sc$q, scan.p = sc$p),
                    unit.size = stats::setNames(msize, sprintf("%g", scale)), kernel.alpha = kernel.alpha, moderate = if (length(mp)) "trend" else "none", group.test = gt, n.permutations = sc$B)
        if (post == "inherit") {
          ih <- inheritNeighbourQ(res$p.adj, res$effect, g.s); res$p.adj.cell <- res$p.adj; res$p.adj <- ih$q; res$effect.inherited <- ih$effect; res$adj <- tz(ih$q, sign(ih$effect))
        }
        res$z.pooled <- res$pooled; res$z.adj <- res$adj
        if (verbose) { n.adj <- sum(res$p.adj < 0.1, na.rm = TRUE)   # the two-sided adjusted z of q = 0.1 is 1.64; count the q-values directly
          message(sprintf("%d of %d cells at FDR 10 %%%s", n.adj, sum(is.finite(res$p.adj)), if (auto) sprintf(" (scale %g %% chosen)", 100 * scale.used) else "")) }
        self$test.results[[name]]$diff[[type]] <- res
        return(invisible(self$test.results[[name]]))
      }

      # the C++ accumulator (exact p, standardized z, pooled histograms, group sums, optional neighbour averaging) replaces
      # the permutations x cells matrix for every adjustment but max-t, which smooths the permuted maps on the graph
      use.acc <- adjust.method != "max-t"
      cells.d <- rownames(density.mat); grp.levels <- NULL; cg.codes <- NULL
      if (!is.null(self$cell.groups) && dens.res$method == "graph") { gf <- as.factor(self$cell.groups[cells.d]); grp.levels <- levels(gf); cg.codes <- as.integer(gf) - 1L; cg.codes[is.na(cg.codes)] <- -1L
        if (!any(cg.codes >= 0)) { grp.levels <- NULL; cg.codes <- NULL } }
      nn.graph <- if (post != "none" && dens.res$method == "graph") graphNeighbourLists0(private$cellGraph(), cells.d) else NULL
      perm.res <- density.mat %>%
          diffCellDensityPermutations(sample.model=sample.model, seed = private$opt("seed"), perm.method=perm.method, statistic=statistic, robust.method = robust.method,
                                      na.mode = na.mode, alternative = alternative, n.permutations=n.permutations,
                                      return.residuals=return.residuals, return.sampled.stats=return.sampled.stats || !use.acc,
                                      n.cores=n.cores, verbose=verbose, moderate=moderate, moderate.d0.max = moderate.d0.max,
                                      pooled = use.acc, group = if (use.acc) cg.codes else NULL, average = use.acc && post == "averaged", nn.graph = nn.graph, pooled.statistic = pooled.statistic)
      if (post == "averaged" && !use.acc) { if (verbose) message("post = 'averaged' is not available with max-t: post = 'none'"); post <- "none" }
      p.pooled <- perm.res$p.pooled   # the unadjusted z-scores (`pooled`) are the signed tail z of this p
      if (is.null(p.pooled) && verbose) message("the permuted statistics were not kept (return.sampled.stats = FALSE): the unadjusted z-scores come from the floored per-cell permutation p-values")
      if (p.source == "pooled" && is.null(p.pooled)) stop("p.source = 'pooled' needs the permutation statistics (return.sampled.stats = TRUE)")
      if(!adjust){
        score <- perm.res %>% .$z.score
        res <- list(pooled=score, model=sample.model, perm.method=perm.method, statistic=statistic, robust = robust.method, effect = perm.res$effect, se = perm.res$se,
                    t = perm.res$t, p.param = perm.res$p.param)
      } else if (adjust.method == "max-t") {
        res <- list(pooled=perm.res$z.score, adj=perm.res %$% adjustZScoresByPermutations(
                    score, permut.scores, smooth=smooth, graph=graph, n.cores=n.cores, verbose=verbose,
                    l.max=l.max, ...), model = sample.model, perm.method=perm.method, statistic=statistic, robust = robust.method, effect = perm.res$effect, se = perm.res$se,
                    t = perm.res$t, p.param = perm.res$p.param, adjust.method = adjust.method)
      } else if (adjust.method == "perm-fdr") {
        if (!is.null(perm.res$pooled.hist)) {   # pooled FDR of the pooled statistic (|t| by default) from the accumulated histograms
          pf <- permutationFDR(perm.res$pooled.hist$s.obs, perm.res$pooled.hist, alternative = "greater")
          pa <- pmin(pmax(pf$q, 1e-300), 1 - 1e-16)
          z.adj <- switch(alternative, "two-sided" = sign(perm.res$effect) * stats::qnorm(pa / 2, lower.tail = FALSE), greater = pmax(stats::qnorm(pa, lower.tail = FALSE), 0), less = pmin(stats::qnorm(pa), 0))
          z.adj[!is.finite(pf$q)] <- NA_real_; names(z.adj) <- names(pf$q)
          ad <- list(z.adj = z.adj, p.adj = pf$q, fdr.curve = pf$curve)
        } else {
          if (is.null(perm.res$stats.perm)) stop("perm-fdr needs the permutation statistics (return.sampled.stats = TRUE)")
          ad <- adjustCellScores(perm.res$stat.obs, NULL, method = "perm-fdr", perm = perm.res$stats.perm,
                                 alternative = if (alternative == "two-sided") "two.sided" else "greater")
        }
        res <- list(pooled=perm.res$z.score, adj=ad$z.adj, p.adj = ad$p.adj, fdr.curve = ad$fdr.curve,
                    model = sample.model, perm.method=perm.method, statistic=statistic, robust = robust.method, effect = perm.res$effect, se = perm.res$se,
                    t = perm.res$t, p.param = perm.res$p.param, adjust.method = adjust.method, p.source = "permutation")
      } else {
        # false-discovery-rate alternatives on per-cell p-values (parametric t by default: no permutation floor)
        z.in <- if (p.source == "parametric") perm.res$t else perm.res$z.score
        p.in <- switch(p.source, parametric = perm.res$p.param, permutation = perm.res$p.perm, pooled = p.pooled)
        nn <- NULL
        if (adjust.method != "bh") {
          cells <- names(z.in)
          g <- if (dens.res$method == 'graph') private$cellGraph() else graph
          nn <- graphNeighbourList(g, cells)
          if (length(nn) != length(cells)) stop("the graph does not cover all tested cells")
        }
        ad <- adjustCellScores(z.in, p.in, method = adjust.method, nn.list = nn, seed = private$opt("seed") %||% 1)
        res <- list(pooled=perm.res$z.score, adj=ad$z.adj, p.adj = ad$p.adj, index.cells = if (!is.null(ad$index)) names(z.in)[ad$index] else NULL,
                    model = sample.model, perm.method=perm.method, statistic=statistic, robust = robust.method, effect = perm.res$effect, se = perm.res$se,
                    t = perm.res$t, p.param = perm.res$p.param, adjust.method = adjust.method, p.source = p.source)
      }

      if (post == "inherit" && !is.null(res$p.adj)) {
        ih <- inheritNeighbourQ(res$p.adj, res$effect, private$cellGraph())
        res$p.adj.cell <- res$p.adj; res$p.adj <- ih$q; res$effect.inherited <- ih$effect
        pa <- pmin(pmax(ih$q, 1e-300), 1 - 1e-16)
        res$adj <- if (alternative == "two-sided") sign(ih$effect) * stats::qnorm(pa / 2, lower.tail = FALSE) else pmax(stats::qnorm(pa, lower.tail = FALSE), 0)
        res$adj[!is.finite(ih$q)] <- NA_real_; names(res$adj) <- names(ih$q)
      }
      res$post <- post; res$unit.size <- dens.res$unit.size
      res$transform <- transform; res$p.perm <- perm.res$p.perm; res$p.pooled <- p.pooled; res$z.std <- perm.res$z.std; res$closed.form <- perm.res$closed.form; res$pooled.statistic <- perm.res$pooled.statistic
      res$z.pooled <- res$pooled; res$z.adj <- res$adj   # the names of the shift map
      if (!is.null(perm.res$group.acc) && !is.null(grp.levels) && length(perm.res$group.acc$group_n) == length(grp.levels)) {   # exact group-level permutation test from the accumulated group sums
        res$group.test <- groupTestFromAccum(perm.res$group.acc, grp.levels, alternative = if (alternative == "two-sided") "two.sided" else "greater", negate = alternative == "less")
      } else if (!is.null(self$cell.groups) && !is.null(perm.res$stats.perm)) {
        res$group.test <- groupPermutationTest(perm.res$stat.obs, perm.res$stats.perm, self$cell.groups, alternative = if (alternative == "two-sided") "two.sided" else "greater")
      }
      if (return.residuals) {
        res$residuals <- perm.res$residuals
        res$residuals.pearson <- perm.res$residuals.pearson
      }
      if (verbose && adjust) {
        n.adj <- sum(res$p.adj < 0.1, na.rm = TRUE)   # the two-sided adjusted z of q = 0.1 is 1.64; count the q-values directly
        message(sprintf("%d of %d cells at FDR 10 %%%s", n.adj, sum(is.finite(res$p.adj)),
                        if (n.adj == 0) " (none: more permutations, or plotDiffCellDensity(adjust = FALSE) for the raw scores)" else ""))
      }

      self$test.results[[name]]$diff[[type]] <- res

      return(invisible(self$test.results[[name]]))
    },

    #' @description Estimate differential cell density
    #' @param name character Slot with results from estimateCellDensity. New results will be appended there. (Default: 'cell.density')
    #' @param size numeric (default=0.2)
    #' @param palette color palette, default is c('blue','white','red')
    #' @param adjust plot the adjusted z-scores (TRUE), the unadjusted pooled-null z-scores (FALSE) or the adjusted ones when
    #'   they were estimated (NULL, default)
    #' @param what `"z"` (default): the signed z-scores; `"effect"`: the contrast estimate per cell (the density
    #'   difference, or the log2 fold change of the density with `transform = "log"`)
    #' @param contours character Specify cell types for contour, multiple cell types are also supported (default: NULL)
    #' @param contour.color character color for contour line (default: 'black')
    #' @param contour.conf character confidence interval of contour (default: '10%')
    #' @param outline draw a thin outline around the regions where the significant cells are (default FALSE; see
    #'   [significanceContours()]): the cells with adjusted p below `outline.alpha`, outlined where they make up at least
    #'   `outline.level` of the cells around a point, in `outline.color` with line width `outline.linewidth`
    #' @param outline.alpha,outline.level,outline.color,outline.linewidth settings of the outline
    #' @param plot.na boolean Plot NAs (default=FALSE)
    #' @param color.range symmetric colour range: a quantile of the absolute values as a string (default `"99.5%"`, so that
    #'   a few extreme cells do not wash out the map; adjusted z-scores keep at least +-1.96) or a numeric vector whose
    #'   largest absolute value sets the range; values beyond it are clipped to the range
    #' @param mid.color character Color code for medium value in color range (default='gray95')
    #' @param scale.z.palette boolean Scale plot palette for Z scores (default: TRUE for z-scores, FALSE for the effect)
    #' @param min.z numeric Minimum Z score to plot (default=qnorm(0.9))
    #' @param max.z cap of the colour range of a z map taken from a quantile (default 8: z-scores of smaller p-values are
    #'   not distinguishable and would wash out the map); a numeric `color.range` is taken as given
    #' @param ... additional parameters
    #' @return ggplot2 object
    #' @examples
    #' \dontrun{
    #' cao$estimateCellDensity()
    #' cao$estimateDiffCellDensity()
    #' cao$plotDiffCellDensity()
    #' }
    plotDiffCellDensity=function(name='cell.density', size=0.2, palette=NULL,
                                 adjust=NULL, what=c("z", "effect"), contours=NULL, contour.color='black', contour.conf='10%',
                                 plot.na=FALSE, color.range="99.5%", max.z=8, mid.color='gray95',
                                 scale.z.palette=NULL, min.z=qnorm(0.9), outline=FALSE, outline.alpha=0.05, outline.level=0.5, outline.color="black", outline.linewidth=0.25, ...) {
      what <- match.arg(what)
      if (is.null(palette)) {
        palette <- grDevices::colorRampPalette(c("blue", mid.color, "red"), space = "Lab")
      }
      private$checkCellEmbedding()
      dens.res <- private$getResults(name, 'estimateCellDensity')
           
      scores <- dens.res$diff[[1]]
      if (is.null(scores)) {                                             
             warning("Can't find differential density results for name '", name,
                 "'. Running estimateDiffCellDensity() with default parameters.")
         self$estimateDiffCellDensity(name = name, adjust = adjust)       
         dens.res <- self$test.results[[name]]
         scores  <- dens.res$diff[[1]]                                    
      }

      res <- scores
      unadj <- res$pooled %||% res$raw   # `raw` in results of earlier versions
      if (what == "effect") {
        if (is.null(res$effect)) stop("no per-cell effect in the results of '", name, "' (type = 'permutation' is needed)")
        scores <- res$effect; adjust <- FALSE
        leg.title <- if (identical(res$transform, "log")) "log2 FC" else "density diff."
      } else if (is.null(adjust)) {
        if (!is.null(res$adj)) {
          scores <- res$adj
          adjust <- TRUE
        } else {
          scores <- unadj
          adjust <- FALSE
        }
      } else if (adjust) {
        if (is.null(res$adj)) {
          warning("Adjusted scores are not estimated. Using the unadjusted scores. ",
                  "Please, run estimateDiffCellDensity with adjust=TRUE")
          scores <- unadj
        } else {
          scores <- res$adj
        }
      } else {
        scores <- unadj
      }
      if (what == "z") leg.title <- if (adjust) 'Z adj.' else 'Z (pooled)'
      if (is.null(scale.z.palette)) scale.z.palette <- (what == "z")

      if (dens.res$method == 'graph') {
        density.emb <- self$embedding
        scores %<>% .[intersect(names(.), rownames(density.emb))]
      } else if (dens.res$method == 'kde') {
        density.emb <- dens.res$density.emb[,1:2]
      } else stop("Unknown method: ", dens.res$method)

      density.emb <- density.emb[names(scores),]

      # adjusted scores: keep at least +-1.96 in the range so that a map without significant cells looks pale
      color.range <- symmetricColorRange(scores, color.range, floor = if (adjust) stats::qnorm(0.975) else 0, cap = if (what == "z") max.z else Inf)
      scores %<>% pmin(color.range[2]) %>% pmax(color.range[1])

      gg <- self$plotEmbedding(density.emb, colors=scores, size=size, legend.title=leg.title, palette=palette,
                               midpoint=0, plot.na=plot.na, color.range=color.range, ...)

      if (scale.z.palette) {                               
        gg$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
        gg <- gg + getScaledZGradient(min.z=min.z, palette=palette, color.range=color.range)
     }

      if (outline) {   # the regions where the cells significant at `outline.alpha` are
        q <- res$p.adj %||% (if (!is.null(res$adj)) 2 * stats::pnorm(-abs(res$adj)) else NULL)
        if (is.null(q)) warning("no adjusted p-values in the results of '", name, "': no outline") else {
          ol <- significanceContourLayer(self$embedding, stats::setNames(q < outline.alpha, names(q)), level = outline.level, color = outline.color, linewidth = outline.linewidth)
          if (!is.null(ol)) gg <- gg + ol }
      }
      if (!is.null(contours)) {
        gg <- gg + private$getDensityContours(groups=contours, conf=contour.conf, color=contour.color)
      }
      return(gg)
    },
                                       
    #' @description Plot residual diagnostics for differential cell density
    #' @param name character Slot with results from estimateCellDensity (default: 'cell.density')
    #' @param type character key in diff list, e.g. 'permutation' (default)
    #' @param size numeric point size for embedding (default=0.2)
    #' @param palette color palette for embedding plots (default: brewerPalette("BuGn"))
    #' @param alpha numeric transparency for embedding points (default=0.2)
    #' @param sample.metadata data.frame with sample-level covariates (default = self$sample.meta).
    #'        Row names should match residual row names (sample IDs).
    #' @param cov.plot.keys character vector of column names in sample.metadata to plot against
    #'        sample-level residual metrics (default = NULL, meaning no sample-covariate plots).
    #' @param metric character, which sample-level residual metric to use on y-axis:
    #'        "mean_abs" (mean |residual| per sample) or "sd" (SD of residuals per sample).
    #' @param build.panel logical; if TRUE and length(cov.plot.keys) == 1, returns a combined panel
    #'        (2 embedding plots + 1 sample-covariate plot). Otherwise returns a list of ggplots.
    #' @param cont.palette palette for continuous covariates (numeric)
    #' @param ... passed to plotEmbedding
    #' @return ggplot or list of ggplots
    #' @param jitter.size point size
    #' @param jitter.alpha point transparency
    #' @param residual.type 'pearson' or 'raw' residuals
    plotDiffCellDensityResiduals = function(name = "cell.density", type = "permutation", size = 0.2, palette = brewerPalette("GnBu"),
                                            alpha = 0.2, sample.metadata = self$sample.meta, cov.plot.keys = NULL, 
                                            metric = c("mean_abs", "sd"), build.panel = TRUE, jitter.size = 1, jitter.alpha = 0.8,
                                            cont.palette  = rev(RColorBrewer::brewer.pal(11, "Spectral")), residual.type = c("pearson","raw"),
                                             ...) {
      metric <- match.arg(metric)
      residual.type <- match.arg(residual.type)
      private$checkCellEmbedding()

      dens.res <- private$getResults(name, "estimateCellDensity")

      diff.res <- dens.res$diff[[type]]
      if (is.null(diff.res)) {
          stop("Can't find differential density results for name '", name,
              "' and type '", type, "'. Run estimateDiffCellDensity(..., type='", type,
              "', return.residuals=TRUE).")
      }

      if (is.null(diff.res$residuals)) {
          stop("Residuals are not stored. Re-run estimateDiffCellDensity(..., return.residuals=TRUE).")
      }

      res <- if (residual.type == "pearson") {
               if (!is.null(diff.res$residuals.pearson)) diff.res$residuals.pearson else diff.res$residuals
             } else {
               diff.res$residuals # matrix: samples 
             }
  
      ## ---- per-cell metrics (across samples) ----
      mean.abs.resid.cell <- colMeans(abs(res), na.rm = TRUE)
      sd.resid.cell <- apply(res, 2, sd, na.rm = TRUE)

      ## align cells with embedding / density embedding
      if (dens.res$method == "graph") {
        density.emb <- self$embedding
        common.cells <- intersect(colnames(res), rownames(density.emb))
        mean.abs.resid.cell <- mean.abs.resid.cell[common.cells]
        sd.resid.cell <- sd.resid.cell[common.cells]
        density.emb  <- density.emb[common.cells, , drop = FALSE]
      } else if (dens.res$method == "kde") {
        density.emb <- dens.res$density.emb[, 1:2]
        common.cells <- intersect(colnames(res), rownames(density.emb))
        mean.abs.resid.cell <- mean.abs.resid.cell[common.cells]
        sd.resid.cell <- sd.resid.cell[common.cells]
        density.emb <- density.emb[common.cells, , drop = FALSE]
      } else {
        stop("Unknown density method: ", dens.res$method)
      }

      ## ---- embedding plots: where the density model fits badly ----
      g.mean <- self$plotEmbedding(density.emb, colors = mean.abs.resid.cell, size = size,
                               alpha = alpha, legend.title = "Density: mean |residual|",
                               palette = palette, ...)
      g.sd <- self$plotEmbedding(density.emb, colors = sd.resid.cell, size = size, alpha = alpha, 
                             legend.title = "Density: residual SD", palette = palette, ...)
      ## ---- sample-level metrics ----
      sample.ids <- rownames(res)
      if (is.null(sample.ids)) {
         stop("Residual matrix must have rownames corresponding to samples.")
      }
      # subset / reorder sample.metadata to match residuals
      if (!all(sample.ids %in% rownames(sample.metadata))) {
        stop("Not all residual sample IDs are present as rownames in sample.metadata.")
      }
      sm <- sample.metadata[sample.ids, , drop = FALSE]

      mean.abs.resid.sample <- rowMeans(abs(res), na.rm = TRUE)
      sd.resid.sample <- apply(res, 1, sd, na.rm = TRUE)

      metric.vec <- if (metric == "mean_abs") mean.abs.resid.sample else sd.resid.sample
      ylab <- if (metric == "mean_abs") {
                if (residual.type == "pearson") "Mean |Pearson residual| per sample" else "Mean |residual| per sample"
              } else {
                if (residual.type == "pearson") "Pearson residual SD per sample" else "Residual SD per sample"
              }

      # avoid duplicate column names (sample is commonly present in metadata)
      sm2 <- sm
      if ("sample" %in% names(sm2)) sm2 <- sm2[, names(sm2) != "sample", drop = FALSE]

      sample.df <- data.frame(sample = sample.ids, resid_metric = metric.vec, sm2, check.names = FALSE)
      names(sample.df) <- make.unique(names(sample.df))

      ## ---- sample-covariate plots ----
      sample.plots <- list()
      if (!is.null(cov.plot.keys)) {
         for (covkey in cov.plot.keys) {
           if (!covkey %in% colnames(sample.df)) {
            warning("Covariate '", covkey, "' not found in sample.metadata; skipping.")
           next
           }

          cov.vec <- sample.df[[covkey]]
          is.cont <- is.numeric(cov.vec)

          if (is.cont) {
            p <- ggplot(sample.df, aes(x = .data[[covkey]], y = resid_metric, color = .data[[covkey]])) +
              geom_point(size = jitter.size, alpha = jitter.alpha) + scale_colour_gradientn(colors = cont.palette) +
              theme_bw() + labs(x = covkey,y = ylab,color = covkey, title = paste("Residual metric vs", covkey))
          } else {
            cov.fac <- factor(cov.vec)
            p <- ggplot(sample.df, aes(x = cov.fac, y = resid_metric, color = cov.fac)) +
              geom_boxplot(notch = TRUE, outlier.shape = NA,
                           fill = "grey90", colour = "grey30", alpha = 0.9) +
              geom_jitter(width = 0.15, height = 0,
                          size = jitter.size, alpha = jitter.alpha) + theme_bw() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
              labs(x = covkey,y = ylab, color = covkey, title = paste("Residual metric by", covkey))
          }

          sample.plots[[covkey]] <- p
         }
      }

      ## ---- return / panel assembly ----
      if (build.panel && length(sample.plots) == 1) {
       # if exactly one covariate, build a 3-panel figure
        g.sample <- sample.plots[[1]]
        return(cowplot::plot_grid(g.mean, g.sd, g.sample, labels = c("Mean |Residual|", "Residual SD", "Sample-level"),
                                  ncol = 3))
      }

      # otherwise return as a list
      out <- list(mean_cell = g.mean, sd_cell   = g.sd, sample = sample.plots)
      return(out)
    },

    #' @description Plot inter-sample expression distance. The inputs to this function are the results from cao$estimateExpressionShiftMagnitudes()
    #' @param name character Test results to plot (default="expression.shifts")
    #' @param joint boolean Whether to show joint boxplot with the expression distance weighed by the sizes of cell types (default: TRUE), or show distances for each individual cell type
    #' @param condition character Metadata column to group samples by (default: inferred from contrast)
    #' @param values character One of "unadjusted" (observed raw distances), "adjusted" (after regressing out nuissance covariates) (default="unadjusted")
    #' @param palette plot palette (default=self$sample.groups.palette)
    #' @param show.significance boolean Whether to show statistical significance between sample groups. wilcox.test was used; (`*` < 0.05; `**` < 0.01; `***` < 0.001)
    #' @param ... other plot parameters, forwarded to \link{plotCountBoxplotsPerType}
    #' @return A ggplot2 object
    plotExpressionDistance = function(name='expression.shifts', joint=FALSE, condition=NULL,
                                      values=c("unadjusted", "adjusted"), palette=NULL,
                                      show.significance=FALSE, ...) {
      values <- match.arg(values)
      clust.info <- private$getResults(name, 'estimateExpressionShiftMagnitudes()')
      
      # grouping variable: the primary test's variable
      if (is.null(condition)) condition <- private$testVariable()
      if (is.null(palette)) palette <- self$sample.groups.palette
      
      # Helper to melt square matrix and filter for within-group distances
      get_within_group_dists <- function(mat, var_name) {
        df <- reshape2::melt(as.matrix(mat))
        
        # Keep unique pairs only (avoid duplicates & self-distances)
        df <- df[as.character(df$Var1) < as.character(df$Var2), ]
        
        # Annotate with metadata
        df$cond1 <- self$sample.meta[as.character(df$Var1), condition]
        df$cond2 <- self$sample.meta[as.character(df$Var2), condition]
        
        # Filter for within-condition comparisons only
        df <- df[df$cond1 == df$cond2, ]
        df$group <- df$cond1
        df$variable <- var_name
        
        # Drop NAs
        df <- df[!is.na(df$value), ]
        return(df[, c("Var1", "Var2", "value", "group", "variable")])
      }
      
      if (joint) {
        dist_mat <- self$getSampleDistanceMatrix(space="expression.shifts", name=name, values=values)
        if (is.null(dist_mat)) stop("Could not compute joint distance matrix.")
        
        df <- get_within_group_dists(dist_mat, var_name="")
        
        plot.theme <- self$plot.theme +
          theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                axis.ticks.x=element_blank(), panel.grid.major.x=element_blank())
      } else {
        cell_types <- colnames(clust.info$res$sample.distances)
        
        df_list <- lapply(cell_types, function(ct) {
          dist_mat <- self$getSampleDistanceMatrix(space="expression.shifts", name=name, cell.type=ct, values=values)
          if (is.null(dist_mat)) return(NULL)
          get_within_group_dists(dist_mat, var_name=ct)
        })
        
        df <- do.call(rbind, df_list)
        plot.theme <- self$plot.theme
      }
      
      if (nrow(df) == 0) {
        stop("No within-group distances could be calculated. Please check your condition and sample metadata.")
      }
      
      # Safety net for wilcox.test
      if (show.significance) {
        # Check how many distinct groups have data for each cell type
        counts <- table(df$variable, df$group)
        valid_vars <- rownames(counts)[rowSums(counts > 0) >= 2]
        invalid_vars <- setdiff(unique(df$variable), valid_vars)
        
        if (length(invalid_vars) > 0) {
          warning("The following cell types do not have enough pairs in both conditions for Wilcoxon testing and will be excluded from the plot: ", 
                  paste(invalid_vars, collapse = ", "))
          df <- df[df$variable %in% valid_vars, ]
        }
      }
      
      gg <- plotCountBoxplotsPerType(df, y.lab="expression distance", show.significance=show.significance,
                                     plot.theme=plot.theme, palette=palette, ...)
      
      return(gg)
    },
    

    #' @description Sample x sample distance matrix of a sample space
    #'
    #' The spaces: `"expression"` (pseudobulk expression of one cell type, `cell.type = "X"`, or, with `cell.type = NULL`,
    #' all cell types combined: the average of the per-type distances, each type weighted by its cells, `joint.weighting =
    #' "cells"`, or equally, `"equal"`); `"pseudobulk"` (one pseudobulk of all cells per sample, so differences in
    #' composition enter the distance); `"composition"` (Euclidean distance of the ILR coordinates of the cell-type
    #' counts). Expression distances are taken from the stored results of `estimateExpressionShiftMagnitudes()` when there
    #' are any (`name`), so that they are the distances behind the shift estimates; otherwise they are computed with the
    #' same defaults. Nothing requires a model, except `values = "adjusted"`.
    #' @param space `"expression"` (default), `"pseudobulk"` or `"composition"` (`"expression.shifts"` and `"coda"` are
    #'   accepted as the earlier names)
    #' @param values `"unadjusted"` (default) or `"adjusted"`: the model's nuisance part removed (the covariates of the
    #'   model other than the tested variable; for stored expression-shift results, the adjusted distances of the fit)
    #' @param cell.type cell type for `space = "expression"` (default NULL: all cell types combined)
    #' @param joint.weighting `"cells"` (default) or `"equal"`: weights of the cell types in the combined expression distance
    #' @param dist distance of the expression profiles (`"cor"`, `"l1"`, `"l2"`; default: the object option)
    #' @param name results of `estimateExpressionShiftMagnitudes()` to take the expression distances from (default
    #'   "expression.shifts" when present)
    #' @param verbose print messages (default: the object option)
    #' @param sample.subset optional subset of samples
    #' @param adjust.for optional formula of covariates regressed out of the distances (G_adj = R G R)
    #' @param ... `test` (which test's adjusted distances of stored shift results)
    #' @return sample x sample distance matrix with attributes `dist.type`, `space` and `n.cells`
    #' @examples
    #' \dontrun{
    #' cao$getSampleDistanceMatrix(space = "composition")
    #' }
    getSampleDistanceMatrix=function(space=c("expression", "pseudobulk", "composition"), values=c("unadjusted", "adjusted"),
                                     cell.type=NULL, joint.weighting=c("cells", "equal"), dist=NULL, name=NULL, verbose=NULL,
                                     sample.subset=NULL, adjust.for=NULL, ...) {
      space <- normalizeSampleSpace(space); values <- match.arg(values); joint.weighting <- match.arg(joint.weighting)
      verbose <- private$opt("verbose", verbose)
      extra <- setdiff(names(list(...)), "test")
      if (length(extra)) stop("Unexpected arguments: ", paste(extra, collapse = ", "))
      if (!is.null(cell.type) && space != "expression") stop("cell.type applies to space = 'expression' only")
      res <- if (space == "expression") self$test.results[[name %||% "expression.shifts"]] else NULL
      if (!is.null(name) && is.null(res) && space == "expression") stop("no results named '", name, "': run estimateExpressionShiftMagnitudes() first")
      if (!is.null(res) && is.null(res$distances)) res <- NULL
      p.dists <- NULL
      if (values == "adjusted" && !is.null(res)) {   # the fit's own adjusted distances
        tst <- list(...)$test
        mats <- if (is.null(tst)) res$adjusted.distances[[1]] else res$adjusted.distances[[tst]]
        if (is.null(mats)) stop("no adjusted distances for the requested test")
        p.dists <- private$sampleSpaceFromList(mats, res, cell.type, joint.weighting)
        if (is.null(p.dists)) return(NULL)
        attr(p.dists, "dist.type") <- res$settings$dist
      } else {
        p.dists <- if (!is.null(res)) private$sampleSpaceFromList(res$distances, res, cell.type, joint.weighting) else
          private$sampleSpaceDistance(space, cell.type = cell.type, joint.weighting = joint.weighting, dist = dist, verbose = verbose)
        if (is.null(p.dists)) return(NULL)
        if (!is.null(res)) attr(p.dists, "dist.type") <- res$settings$dist
        if (values == "adjusted") {
          nuis <- private$modelNuisanceFormula()
          if (is.null(nuis)) { if (verbose) message("the model has no covariates to adjust for: unadjusted distances"); }
          else adjust.for <- if (is.null(adjust.for)) nuis else stats::as.formula(paste("~", paste(unique(c(all.vars(nuis), all.vars(stats::as.formula(adjust.for)))), collapse = " + ")))
        }
      }
      dist.type <- attr(p.dists, "dist.type"); nc <- attr(p.dists, "n.cells")
      if (!is.null(sample.subset)) p.dists <- p.dists[intersect(rownames(p.dists), sample.subset), intersect(rownames(p.dists), sample.subset), drop = FALSE]
      if (!is.null(adjust.for)) {   # regress the covariates of `adjust.for` out of the Gower matrix: G_adj = R G R
        meta <- self$sample.meta[rownames(p.dists), , drop = FALSE]
        Z <- do.call(cbind, lapply(intersect(all.vars(stats::as.formula(adjust.for)), names(meta)), function(v) covariateColumns(meta, v)))
        if (!is.null(Z)) {
          ok <- stats::complete.cases(Z)
          G <- gowerCenter(asSquaredDistance(p.dists[ok, ok], dist.type))
          p.dists <- adjustedDistanceMatrix(G, cbind(1, Z[ok, , drop = FALSE]), dist.type)
        }
      }
      attr(p.dists, "dist.type") <- dist.type; attr(p.dists, "space") <- space
      if (!is.null(nc)) attr(p.dists, "n.cells") <- nc[rownames(p.dists)]
      p.dists
    },

    #' @description MDS map of the samples in a sample space, colored by a metadata column
    #'
    #' Draws the samples in two dimensions (classical MDS of the distances of `getSampleDistanceMatrix()`), colored and
    #' shaped by metadata columns. The subtitle reports how much of the between-sample variation each of them explains on
    #' exactly these distances: the chance-corrected R2 and permutation p-value of the covariate screen (marginal, as in
    #' `screenCovariates()`). Works before a model is set; `values = "adjusted"` or `"both"` need one.
    #' @param space `"expression"` (default), `"pseudobulk"` or `"composition"` (see `getSampleDistanceMatrix()`)
    #' @param method `"MDS"` (default), `"tSNE"`, `"UMAP"` or `"heatmap"`
    #' @param values `"unadjusted"` (default), `"adjusted"` or `"both"` (side by side)
    #' @param dist distance of the expression profiles (default: the object option)
    #' @param name results of `estimateExpressionShiftMagnitudes()` to take the expression distances from
    #' @param cell.type cell type (`space = "expression"`; default NULL: all cell types combined)
    #' @param joint.weighting `"cells"` (default) or `"equal"` (see `getSampleDistanceMatrix()`)
    #' @param sample.meta sample metadata (default: the object's)
    #' @param color.by,shape.by metadata columns for colour and shape; `color.by = "test"` colours by the tested variable
    #' @param palette colours (default: the stored palette for the tested variable)
    #' @param show.sample.size point size proportional to log10 of the cells per sample (default FALSE)
    #' @param sample.colors per-sample colors (overrides the `color.by` colors)
    #' @param color.title,shape.title,title legend and plot titles
    #' @param both.ncol,both.align layout of `values = "both"`
    #' @param show.r2 report the chance-corrected R2 of `color.by` / `shape.by` in the subtitle (default TRUE)
    #' @param n.permutations permutations of the R2 p-values (default: the object option)
    #' @param sample.subset optional subset of samples
    #' @param n.cores cores (default: the object option)
    #' @param adjust.for optional formula of covariates regressed out of the distances
    #' @param ... passed to the embedding plot (`font.size`, `size`, `show.labels`, ...)
    #' @return ggplot2 object
    #' @examples
    #' \dontrun{
    #' cao$plotSampleDistances(color.by = "batch")
    #' cao$plotSampleDistances(space = "composition", color.by = "condition")
    #' }
    plotSampleDistances=function(space="expression", method="MDS", values="unadjusted", dist=NULL, name=NULL, cell.type=NULL,
                                 joint.weighting=c("cells", "equal"), sample.meta=NULL, color.by=NULL, shape.by=NULL, palette=NULL,
                                 show.sample.size=FALSE, sample.colors=NULL, color.title=NULL, shape.title=NULL, title=NULL,
                                 both.ncol=2, both.align="hv", show.r2=TRUE, n.permutations=NULL, sample.subset=NULL,
                                 n.cores=NULL, adjust.for=NULL, ...) {
      space <- normalizeSampleSpace(space); joint.weighting <- match.arg(joint.weighting)
      values <- match.arg(values, c("unadjusted", "adjusted", "both"))
      n.permutations <- private$opt("n.permutations", n.permutations); n.cores <- private$opt("n.cores", n.cores)
      dots <- list(...)
      if (!is.null(dots$show.pvalues)) { show.r2 <- isTRUE(dots$show.pvalues) || show.r2; dots$show.pvalues <- NULL }   # earlier name
      sample.meta <- sample.meta %||% self$sample.meta
      if (identical(color.by, "test") || identical(shape.by, "test")) {   # the active test's groups / values
        if (is.null(self$model)) stop("color.by = 'test' needs a model (setModel())")
        tv <- private$testVariable(); sg <- self$getSampleGroups()
        sample.meta$test <- if (!is.null(sg)) as.vector(sg[rownames(sample.meta)]) else sample.meta[[tv]]
        if (identical(color.by, "test")) color.title <- color.title %||% tv
        if (identical(shape.by, "test")) shape.title <- shape.title %||% tv
      }
      for (v in c(color.by, shape.by)) if (!v %in% colnames(sample.meta)) stop("column not found in sample.meta: ", v)
      color.title <- color.title %||% color.by; shape.title <- shape.title %||% shape.by
      if (is.null(palette) && !is.null(color.by) && !is.null(names(self$sample.groups.palette))) {   # the condition colors of the other plots
        is.test <- !is.null(self$model) && (identical(color.by, "test") || identical(color.by, private$testVariable()))
        vals <- if (identical(color.by, "test")) NULL else (sample.meta %||% self$sample.meta)[[color.by]]
        if (is.test || (!is.null(vals) && !is.numeric(vals) && length(stats::na.omit(unique(as.character(vals)))) > 0 &&
                        all(stats::na.omit(unique(as.character(vals))) %in% names(self$sample.groups.palette))))
          palette <- self$sample.groups.palette
      }
      if (!is.null(adjust.for) && values != "unadjusted") stop("adjust.for replaces the model-adjusted view: use values = 'unadjusted' with it")
      if (values != "unadjusted" && is.null(self$model) && is.null(self$test.results[[name %||% "expression.shifts"]]))
        stop("values = '", values, "' needs a model (setModel()); use values = 'unadjusted'")
      if (values == "both" && is.null(self$test.results[[name %||% "expression.shifts"]]) && is.null(private$modelNuisanceFormula())) {
        message("the model has no covariates to adjust for: showing the unadjusted distances")
        values <- "unadjusted"
      }
      n.cells.per.samp <- if (is.null(cell.type)) table(self$sample.per.cell) else table(self$sample.per.cell[self$cell.groups[names(self$sample.per.cell)] %in% cell.type])
      space.title <- switch(space, expression = if (is.null(cell.type)) sprintf("expression, all cell types%s", if (joint.weighting == "equal") " (equal weights)" else "") else cell.type,
                            pseudobulk = "pooled pseudobulk", composition = "composition")
      onePanel <- function(vals, panel.title) {
        D <- self$getSampleDistanceMatrix(space = space, cell.type = cell.type, joint.weighting = joint.weighting, dist = dist, name = name,
                                          sample.subset = sample.subset, values = vals, adjust.for = adjust.for)
        if (is.null(D)) return(NULL)
        lab <- function(v) if (is.null(v)) NULL else setNames(sample.meta[rownames(D), v], rownames(D))
        adj.vars <- c(if (!is.null(adjust.for)) all.vars(stats::as.formula(adjust.for)),
                      if (vals == "adjusted") { nf <- private$modelNuisanceFormula(); if (!is.null(nf)) all.vars(nf) })
        sub <- if (show.r2) private$sampleSpaceR2Label(D, sample.meta, unique(c(color.by, shape.by)), c(color.title, shape.title)[!duplicated(c(color.by, shape.by))],
                                                       n.permutations = n.permutations, n.cores = n.cores, adjusted = adj.vars) else NULL
        args <- c(list(p.dists = D, sample.labels = lab(color.by), n.cells.per.samp = n.cells.per.samp, method = method, sample.colors = sample.colors,
                       show.sample.size = show.sample.size, palette = palette, color.title = color.title, shape.labels = lab(shape.by), shape.title = shape.title,
                       title = panel.title, plot.theme = self$plot.theme), dots)
        gg <- do.call(plotSampleDistanceMatrix, args)
        if (!is.null(sub) && method != "heatmap") gg <- gg + ggplot2::labs(subtitle = sub) + ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"))
        gg
      }
      if (values != "both") {
        ttl <- title %||% if (!is.null(adjust.for)) sprintf("%s | adjusted for %s", space.title, paste(all.vars(stats::as.formula(adjust.for)), collapse = ", ")) else
          if (values == "adjusted") sprintf("%s | model-adjusted", space.title) else space.title
        return(onePanel(values, ttl))
      }
      checkPackageInstalled("cowplot", cran = TRUE, details = "for values = 'both'")
      g1 <- onePanel("unadjusted", "Unadjusted"); g2 <- onePanel("adjusted", "Model-adjusted")
      if (is.null(g1) || is.null(g2)) return(NULL)
      gg <- cowplot::plot_grid(g1, g2, ncol = both.ncol, align = both.align)
      cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(title %||% space.title, fontface = "bold", x = 0, hjust = 0), gg, ncol = 1, rel_heights = c(0.08, 1))
    },

    #' @description Deprecated: use `screenCovariates()`, whose summary row "all cell types" gives the same R2
    #'
    #' Chance-corrected R2 and permutation p-value of each metadata column on the combined sample distances (the marginal
    #' covariate screen of one matrix).
    #' @param sample.meta sample metadata (default: the object's)
    #' @param space sample space (see `getSampleDistanceMatrix()`; default "expression")
    #' @param dist distance (default: the object option)
    #' @param space.name results to take the expression distances from
    #' @param sample.subset optional subset of samples
    #' @param name results slot (default "metadata.separation")
    #' @param n.permutations permutations (default: the object option)
    #' @param trim,k unused (the graph-based statistic was replaced by the covariate screen)
    #' @param show.warning warn about significant covariates (default TRUE)
    #' @param verbose,n.cores defaults from options
    #' @param adjust.pvalues BH-adjust the p-values (default TRUE)
    #' @param p.adjust.method adjustment method (default "BH")
    #' @param pvalue.cutoff threshold of the warning (default 0.05)
    #' @param mode `"marginal"` (default), `"partial"` or `"both"`
    #' @return list: `metadata`, `pvalues`, `pseudo.r2`, `padjust`, `screen`
    estimateMetadataSeparation=function(sample.meta=self$sample.meta, space="expression", dist=NULL, space.name=NULL,
                                        sample.subset=NULL, name="metadata.separation", n.permutations=NULL, trim=NULL, k=NULL,
                                        show.warning=TRUE, verbose=NULL, n.cores=NULL,
                                        adjust.pvalues=TRUE, p.adjust.method="BH", pvalue.cutoff=0.05, mode=c("marginal", "partial", "both")) {
      mode <- match.arg(mode); verbose <- private$opt("verbose", verbose); n.cores <- private$opt("n.cores", n.cores)
      n.permutations <- n.permutations %||% private$opt("n.permutations")
      if (verbose) message("estimateMetadataSeparation() is deprecated: screenCovariates() reports the same R2 in its 'all cell types' row")
      if (!is.null(trim) || !is.null(k)) message("`trim` and `k` are deprecated (the graph-based pseudo-R2 was replaced by the covariate screen)")
      p.dists <- self$getSampleDistanceMatrix(space=space, cell.type=NULL, dist=dist, name=space.name, sample.subset=sample.subset)
      if (is.null(p.dists)) {
        warning("An empty sample distance matrix was returned. Consider changing 'space', 'cell.type', 'name',  or 'sample.subset'.")
        return(NULL)
      }
      if (!is.data.frame(sample.meta)) sample.meta <- as.data.frame(sample.meta)
      if (!any(rownames(sample.meta) %in% rownames(p.dists))) stop("The rownames of the sample.meta object don't match any sample names.")
      meta <- sample.meta[intersect(rownames(p.dists), rownames(sample.meta)), , drop = FALSE]
      sc <- screenCovariates(list(joint = p.dists[rownames(meta), rownames(meta)]), meta, covariates = names(meta), mode = mode,
                             dist = attr(p.dists, "dist.type") %||% "l2", n.permutations = n.permutations,
                             min.samples.per.type = 4, alpha = pvalue.cutoff, seed = private$opt("seed") %||% sample.int(.Machine$integer.max, 1),
                             n.cores = n.cores, verbose = FALSE)
      tb <- sc$table[sc$table$mode == (if (mode == "both") "marginal" else mode), ]
      pvalues <- setNames(tb$p, tb$covariate); pseudo.r2 <- setNames(tb$R2.adj, tb$covariate)
      res <- list(metadata=meta, pvalues=pvalues, pseudo.r2=pseudo.r2, screen=sc)
      if (adjust.pvalues) res$padjust <- p.adjust(pvalues, method=p.adjust.method)
      pv <- if (adjust.pvalues) res$padjust else pvalues
      if (any(is.na(pv))) warning(paste(names(pv)[is.na(pv)], collapse = ", "), " could not be tested (constant or too many missing values?)")
      if (show.warning && any(pv < pvalue.cutoff, na.rm = TRUE)) warning("Significant separation by: ", paste(names(pv)[pv < pvalue.cutoff], collapse=', '))
      self$test.results[[name]] <- res
      return(invisible(res))
    },

    #' @description Deprecated: bars of `estimateMetadataSeparation()` (see `plotCovariateScreen()`)
    #' @param name character Name for storage in test.results (default="metadata.separation")
    #' @param pvalue.y numeric (default=0.93)
    #' @param ... additional parameters forwarded to \code{plotMeanMedValuesPerCellType}
    #' @return ggplot2 object
    plotMetadataSeparation=function(name='metadata.separation', pvalue.y=0.93, ...) {
      res <- private$getResults(name, "estimateMetadataSeparation()")
      pvals <- if (is.null(res$padjust)) res$pvalues else res$padjust
      gg <- res$pseudo.r2 %>% {tibble(Type=names(.), value=ifelse(is.na(.), 0, pmax(., 0)))} %>%
        plotMeanMedValuesPerCellType(type="bar", ylab=expression(adjusted~R^2), jitter.alpha=0, pvalues=pvals,
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
    #' @param adjust.pvalues boolean (default=FALSE)
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
    #' @param sample.per.cell named factor of sample per cell (default: stored)
    #' @param formula optional location formula for a temporary model (see `estimateExpressionShiftMagnitudes()`)
    #' @param test optional test for a temporary model (see `Cacoa$new()`); the test must be a contrast (two groups or a numeric step)
    #' @param contrast optional contrast for a temporary model (expert synonym of `test`)
    #' @param perm.method relabeling scheme: `"auto"` (default: the scheme of the model's test, as for the other analyses), `"block"` or
    #'   `"freedman-lane"`
    #' @param statistic test statistic of the fitter: `"t"` (studentized contrast, default) or `"coef"` (raw contrast estimate)
    #' @param robust.method robust fitting: 'none', 'huber' or 'winsor'
    #' @param na.mode handling of missing responses: 'drop' or 'impute_weak'
    #' @param alternative alternative hypothesis for the z-scores: 'two-sided', 'greater' or 'less'
    #' @param keep.means also return per-condition mean expression
    #' @param block.vars metadata columns defining permutation strata (default: stored)
    #' @param ... further arguments passed to the underlying function
    estimateClusterFreeDE=function(n.top.genes=Inf, genes=NULL, max.z=20, min.expr.frac=0.01, min.n.samp.per.cond=2,
                                   sample.per.cell=self$sample.per.cell, formula=NULL, test=NULL, contrast=NULL,
                                   perm.method = "auto", statistic = c("t", "coef"), robust.method = "none", na.mode = "drop", alternative = "two-sided",
                                   min.n.obs.per.samp=2, adjust.pvalues=FALSE, keep.means=FALSE, robust=FALSE,
                                   smooth=TRUE, wins=0.01, n.permutations=NULL, lfc.pseudocount=1e-5, block.vars=NULL,
                                   min.edge.weight=0.6, verbose=self$verbose, n.cores=self$n.cores, 
                                   name="cluster.free.de", ...){
      statistic <- match.arg(statistic); n.permutations <- private$opt("n.permutations", n.permutations)
      if (is.null(genes)) {
        genes <- private$getTopGenes(n.top.genes, gene.selection="expression", min.expr.frac=min.expr.frac)
      }

      sample.model <- private$resolveModel(formula = formula, test = test, contrast = contrast, block.vars = block.vars, verbose = verbose, what = "estimateClusterFreeDE()")
      if (identical(sample.model$contrast_spec$type, "term")) stop("cluster-free DE needs a contrast test (two groups or a numeric step), not a whole-factor test")
      if (identical(perm.method, "auto")) perm.method <- sample.model$tests[[sample.model$primary %||% 1L]]$permutation$scheme %||% "block"
      perm.method <- match.arg(perm.method, c("freedman-lane", "block"))

      if (verbose)
        message("Estimating cluster-free Z-scores for ", length(genes), " most expressed genes")

      de.inp <- private$getClusterFreeDEInput(genes, min.edge.weight=min.edge.weight)
      #mats <- de.inp %$% clusterFreeZScoreMat(
      #  cm, sample_per_cell=self$sample.per.cell[rownames(cm)], nn_ids=nns.per.cell, is_ref=is.ref,
      #  min_n_samp_per_cond=min.n.samp.per.cond, min_n_obs_per_samp=min.n.obs.per.samp, robust=robust,
      #  norm_both=norm.both, adjust_pvalues=adjust.pvalues, smooth=smooth, wins=wins, n_permutations=n.permutations,
      #  verbose=verbose, n_cores=n.cores
      #)
      mats <- estimateClusterFreeDE_LM(seed = private$opt("seed"), genes, de.inp, sample.per.cell = sample.per.cell, design = sample.model, 
                                                  perm.method = perm.method, statistic = statistic, robust.method = robust.method, 
                                                  na.mode = na.mode, alternative = alternative, n.cores = n.cores,
                                                  min.n.samp.per.cond = min.n.samp.per.cond, min.n.obs.per.samp = min.n.obs.per.samp,
                                                  keep.means = keep.means, lfc.pseudocount = lfc.pseudocount, 
                                                  n.permutations = n.permutations, max.z= max.z, ...) 
      # TODO: add z-score adjustment

      #mats$z@x %<>% pmin(max.z) %>% pmax(-max.z)
      #if (length(mats$z.adj@x) > 0) {
      #  mats$z.adj@x %<>% pmin(max.z) %>% pmax(-max.z)
      #}
      #lf.mat <- mats$reference
      #lf.mat@x <- log2(mats$target@x + lfc.pseudocount) - log2(lf.mat@x + lfc.pseudocount)
      #mats$lfc <- lf.mat

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

    #' @description Cluster-free expression shifts: for every cell, the contrast's shift between samples' mean
    #'   expression over the cell's neighbourhood, with permutation z-scores adjusted across cells
    #'
    #' Uses the stored model (or a temporary one) and the same engine as `estimateExpressionShiftMagnitudes()`:
    #' one set of sample permutations shared by all cells, the max-statistic adjustment of the z-scores, and
    #' median smoothing over the graph.
    #' @param statistic `"gsd"` (default: per-gene test combined over the unit's genes, see [clusterFreeExpressionShifts()]) or
    #'   `"distance"` (the sample-distance pseudo-F, the previous method)
    #' @param n.top.genes optional: restrict the gene universe to this many top genes by `gene.selection` (default NULL: no
    #'   restriction; the distance map used the 3000 most expressed genes, with no per-neighbourhood selection, before 2026-10)
    #' @param gene.selection ranking used with `n.top.genes` (`"expression"` or `"od"`)
    #' @param genes explicit gene universe (default: see `min.gene.frac` and `n.top.genes`)
    #' @param min.expr.frac minimum fraction of cells expressing a gene (default 0)
    #' @param min.gene.frac the gene universe is every gene detected in at least this fraction of all cells (default 0.01)
    #' @param neighbourhood.gene.frac both statistics: in each neighbourhood, the genes of the universe detected in at least this
    #'   fraction of its cells (default 0.05; 0: every gene of the universe)
    #' @param gsd.d0 gsd: prior degrees of freedom of the moderated per-gene variance (the noise model's variance shrunk towards the
    #'   residual variance; default 3); see [clusterFreeExpressionShifts()]
    #' @param centres.per.cell gsd: test the neighbourhoods of a subset of cells, about this many centres among each cell's graph
    #'   neighbours (default 2; `Inf` tests every cell's neighbourhood), see
    #'   [clusterFreeExpressionShifts()]
    #' @param name results slot (default "cluster.free.expr.shifts")
    #' @param test optional temporary model (see `estimateExpressionShiftMagnitudes()`)
    #' @param formula optional temporary model (see `estimateExpressionShiftMagnitudes()`)
    #' @param contrast optional temporary model (see `estimateExpressionShiftMagnitudes()`)
    #' @param block.vars optional temporary model (see `estimateExpressionShiftMagnitudes()`)
    #' @param min.n.obs.per.samp minimum cells of a sample in a neighbourhood (default 3)
    #' @param min.samp.per.level minimum samples per compared level in a neighbourhood (default 2)
    #' @param permutation defaults from options (permutations capped at 499)
    #' @param n.permutations number of permutations (default: the object option `n.permutations`, 999)
    #' @param seed defaults from options (permutations capped at 499)
    #' @param verbose defaults from options (permutations capped at 499)
    #' @param n.cores defaults from options (permutations capped at 499)
    #' @param dist `"cor"` (default) or `"euclidean"` with the normalized profiles; `"cor"`, `"cosine"` or `"js"` with
    #'   `profile = "raw"`; see [clusterFreeExpressionShifts()]
    #' @param profile `"normalized"` (default: depth-normalized, noise-weighted neighbourhood profiles) or `"raw"` (the
    #'   previous mean raw counts, unweighted); see [clusterFreeExpressionShifts()]
    #' @param profile.scale transform scale of the normalized profiles (default `"median"`: the median library size)
    #' @param weighting `"noise"` (default) or `"none"`; see [clusterFreeExpressionShifts()]
    #' @param sample.floor between-sample variance floor in the gene weights (default TRUE)
    #' @param noise.correction known-noise correction (default FALSE); see [clusterFreeExpressionShifts()]
    #' @param adjust see [clusterFreeExpressionShifts()]
    #' @param smooth see [clusterFreeExpressionShifts()]
    #' @param wins see [clusterFreeExpressionShifts()]
    #' @param min.edge.weight minimum graph edge weight between cells of different samples (default 0)
    #' @param k neighbourhood unit: NULL (default) uses kNN neighbourhoods sized by [defaultNeighbourhoodSize()]
    #'   (the larger of 100 and 10 cells per sample, capped at 1 % of the cells; the unit with the best power at
    #'   matched error in the benchmark, `misc/benchmark/README.md`); an integer sets the size
    #'   ([knnNeighbourhoodMatrix()]); `"1-hop"` tests each cell with its graph neighbours only (the former unit);
    #'   with diffusion neighbourhoods (`scale` given) it is the smallest neighbourhood: a weakly connected centre, whose kernel stays on
    #'   itself, gets its `k` cells of largest kernel weight, so that its neighbourhood can be tested
    #' @param scale the neighbourhoods of the test. `NULL`: the kNN neighbourhoods (`k`). Otherwise diffusion
    #'   neighbourhoods on the cell graph ([diffusionNeighbourhoodLists()]; the scale is the share of the cells a typical
    #'   neighbourhood covers, the kernel strength solved per graph, [solveKernelBeta()]; each tested neighbourhood is
    #'   the scale's number of cells its centre's heat reaches most, [diffusionNeighbourhoodLists()]): `"auto"` (default)
    #'   runs the scale scan over `auto.scales`, chooses the scale whose own map calls the most cells at q < 0.1 (the finest
    #'   of the scales calling at least 90 % of the cells when there are such; its strongest 1 % of the cells when none
    #'   calls any) and reports that single-scale map; one value tests a single
    #'   scale; several values (`scale = c(0.01, 0.03, 0.05)`) report the scale scan itself ([scanShiftScales()]: every
    #'   scale under the same relabelings, a cell's statistic its smallest p over the scales, its q-value the pooled
    #'   permutation false discovery rate of that statistic; no Bonferroni factor). `statistic = "gsd"` only. A scale below the floor of
    #'   [defaultNeighbourhoodSize()] cells is raised to it (enough cells per sample). The result records `scale`, `scale.mode`
    #'   (`"single"`, `"scan"`, `"auto"`) and `scales` (per scale: strength, median size, centres, and per cell the pooled
    #'   p, the shift, the single-scale q, the best scale); as for the density map, `estimateDiffCellDensity()`
    #' @param auto.scales the scales `scale = "auto"` chooses from (default 1, 3 and 5 % of the cells, as for the density map: at 10 % the
    #'   scale-sized neighbourhoods of this map take in the halo around a changed population and call it; ignored otherwise)
    #' @param effect.scale under the scan, the scale whose shift estimate is reported as `shifts`: `"finest"` (default) or
    #'   `"best"` (each cell's scale of the smallest p)
    #' @param post post-processing on the cell graph: `"averaged"` (default; the shift statistic averaged over each
    #'   cell and its graph neighbours under every permutation, exact and calibrated; needs `adjust.method =
    #'   "perm-fdr"`), `"inherit"` (smallest q among the cell and its neighbours; more reach, realized error at the
    #'   nominal threshold above the nominal rate) or `"none"`; see [clusterFreeExpressionShifts()]
    #' @param ... deprecated arguments of the previous implementation are accepted and ignored with a message
    #' @return list (see [clusterFreeExpressionShifts()]), also stored in `cao$test.results[[name]]`
    #' @param adjust.method multiplicity adjustment across cells: `"bh"` (default with gsd: BH on the pooled p), `"perm-fdr"`
    #'   (default with the distance map; pooled permutation FDR, see [permutationFDR()]), `"spatial-fdr"`, `"neighbourhood"` or
    #'   `"max-t"`, see [adjustCellScores()]
    #' @param p.source per-cell p-values for the BH-type adjustments (default `"pooled"` with gsd, `"permutation"` with the
    #'   distance map): `"permutation"` (calibrated, floored at
    #'   1 / (n.permutations + 1)), `"pooled"` (pooled permutation null of all cells, see [pooledPermutationP()]: no
    #'   floor, calibrated with enough permutations), `"z"` (normal tail of the permutation-standardized z: no floor,
    #'   anti-conservative) or `"parametric"` (the F of the neighbourhood shift; not calibrated)
    #' @param robust `"none"` (default, option `robust`), `"huber"` or `"winsor"`: robust per-cell shift test
    estimateClusterFreeExpressionShifts=function(statistic=c("gsd", "distance"), n.top.genes=NULL, gene.selection="expression", name="cluster.free.expr.shifts",
                                                 test=NULL, formula=NULL, contrast=NULL, block.vars=NULL,
                                                 min.n.obs.per.samp=3, min.samp.per.level=2, permutation=NULL, n.permutations=NULL, seed=NULL,
                                                 min.gene.frac=0.01, neighbourhood.gene.frac=0.05, centres.per.cell=2, gsd.d0=3,
                                                 dist=c("cor", "euclidean", "cosine", "js"), profile=c("normalized", "raw"), profile.scale="median", weighting=c("noise", "none"),
                                                 sample.floor=TRUE, noise.correction=FALSE, adjust=TRUE, smooth=TRUE, wins=0.025, genes=NULL, min.expr.frac=0.0,
                                                 min.edge.weight=0.0, adjust.method=NULL, p.source=NULL,
                                                 k=NULL, scale="auto", auto.scales=c(0.01, 0.03, 0.05), effect.scale=c("finest", "best"), post=c("averaged", "inherit", "none"), verbose=NULL, n.cores=NULL, robust=NULL, ...) {
      statistic <- match.arg(statistic); gsd <- statistic == "gsd"
      dist <- match.arg(dist); robust <- private$opt("robust", robust); post <- match.arg(post)
      adjust.method <- match.arg(adjust.method %||% (if (gsd) "bh" else "perm-fdr"), c("perm-fdr", "bh", "spatial-fdr", "neighbourhood", "max-t"))
      p.source <- match.arg(p.source %||% (if (gsd) "pooled" else "permutation"), c("permutation", "pooled", "z", "parametric"))
      profile <- match.arg(profile); weighting <- match.arg(weighting)
      verbose <- private$opt("verbose", verbose); n.cores <- private$opt("n.cores", n.cores)
      if (post == "averaged" && !(adjust.method == "perm-fdr" || (adjust.method == "bh" && p.source == "pooled"))) {
        if (verbose) message("post = 'averaged' needs adjust.method = 'perm-fdr' or 'bh' with p.source = 'pooled': post = 'none'"); post <- "none" }
      legacy <- private$legacyShiftArgs(list(...), verbose)
      permutation <- private$opt("permutation", legacy$permutation %||% permutation)
      n.permutations <- n.permutations %||% private$opt("n.permutations")
      seed <- if (missing(seed)) private$opt("seed") else seed
      model <- private$resolveModel(formula = formula, test = test, contrast = contrast, block.vars = block.vars, verbose = verbose,
                                    what = "estimateClusterFreeExpressionShifts()")
      # both statistics: every gene detected in >= min.gene.frac of the cells (optionally the top n.top.genes); per neighbourhood,
      # those detected in >= neighbourhood.gene.frac of its cells (the rule of the per-cell-type analyses)
      if (is.null(genes)) genes <- if (is.null(n.top.genes)) private$getTopGenes(Inf, gene.selection="expression", min.expr.frac=max(min.gene.frac, min.expr.frac))
        else private$getTopGenes(n.top.genes, gene.selection=gene.selection, min.expr.frac=max(min.gene.frac, min.expr.frac))
      if (verbose) message(sprintf("%d genes%s", length(genes), if (neighbourhood.gene.frac > 0) sprintf(" (per neighborhood: those detected in >= %g %% of its cells)", 100 * neighbourhood.gene.frac) else ""))
      if (identical(k, "1-hop")) k <- NULL else if (is.null(k)) k <- defaultNeighbourhoodSize(length(self$sample.per.cell), nlevels(droplevels(as.factor(self$sample.per.cell))))
      inp <- private$getClusterFreeDEInput(genes, raw=TRUE, min.edge.weight=min.edge.weight, k=k)
      if (!is.null(scale) && !gsd) { if (!missing(scale)) stop("the diffusion neighborhoods are tested with statistic = 'gsd'; scale = NULL for the distance statistic")
        if (verbose) message("the distance statistic runs on the kNN neighborhoods (scale = NULL)"); scale <- NULL }
      if (!is.null(scale)) {   # diffusion neighbourhoods at the given scale(s): the scale scan (scanShiftScales())
        auto <- identical(scale, "auto"); if (auto) scale <- checkAutoScales(auto.scales)
        n.s <- nlevels(droplevels(as.factor(self$sample.per.cell[rownames(inp$cm)])))
        scale <- resolveScales(scale, nrow(inp$cm), floor.cells = defaultNeighbourhoodSize(nrow(inp$cm), n.s), verbose = verbose); K <- length(scale); effect.scale <- match.arg(effect.scale)
        k.floor <- k %||% defaultNeighbourhoodSize(nrow(inp$cm), n.s)   # the smallest neighbourhood (a weakly connected centre's): the kNN unit size, so that it can be tested
        cells <- rownames(inp$cm); g.s <- extractCellGraph(self$data.object); g.s <- igraph::induced_subgraph(g.s, cells); g.s <- igraph::permute(g.s, match(igraph::V(g.s)$name, cells))
        kb <- lapply(scale, function(sc) private$kernelBetaForScale(g.s, sc, n.cores))
        if (verbose) message(sprintf("diffusion neighborhoods at %s of the cells%s, at least %d cells each (%.1f per sample)", paste(sprintf("%g %%", 100 * scale), collapse = ", "), if (auto) " (auto: one scale is chosen)" else if (K > 1) " (scale scan)" else "", k.floor, k.floor / n.s))
        sc <- scanShiftScales(cm = Matrix::t(inp$cm), sample.per.cell = self$sample.per.cell[cells], graph = g.s, design = model, meta = model$meta, scale = scale, beta = sapply(kb, `[[`, "beta"), min.size = k.floor,
                              n.permutations = n.permutations, seed = seed %||% sample.int(.Machine$integer.max, 1), n.cores = n.cores, verbose = verbose,
                              neighbourhood.gene.frac = neighbourhood.gene.frac, gsd.d0 = gsd.d0, profile = profile, lib = inp$lib, profile.scale = profile.scale, weighting = weighting,
                              sample.floor = sample.floor, noise.correction = noise.correction, min.n.obs.per.samp = min.n.obs.per.samp, min.samp.per.level = min.samp.per.level,
                              permutation = permutation, block.vars = block.vars %||% model$block.vars, robust = robust, robust.k = private$opt("robust.k"), groups = self$cell.groups)
        msize <- sapply(kb, `[[`, "median.size")
        if (auto || K == 1) { j <- if (K == 1) 1L else chooseScanScale(sc); q <- sc$Q[, j]; p <- sc$P[, j]; eff <- sc$E[, j]; st <- sc$Z[, j]; scale.used <- scale[j]; gt <- sc$group.test[[j]] }
        else { q <- sc$q; p <- sc$p; ks <- if (effect.scale == "finest") rep(1L, length(cells)) else sc$best; eff <- sc$E[cbind(seq_along(cells), ks)]; st <- sc$Z[cbind(seq_along(cells), ks)]; scale.used <- scale; gt <- sc$group.test[[1L]] }
        names(q) <- names(p) <- names(eff) <- names(st) <- cells
        tz <- function(pv) { z <- stats::qnorm(pmin(pmax(pv, 1e-300), 1 - 1e-16), lower.tail = FALSE); z[!is.finite(pv)] <- NA_real_; stats::setNames(pmax(z, 0), cells) }
        shifts <- list(stat = st, z.std = st, p.value = p, p.adj = q, p.pooled = p, z.pooled = tz(p), z.adj = tz(q), pooled = tz(p), adj = tz(q), shifts = eff, post = "averaged", statistic = "gsd", unit = "diffusion",
                       scale = scale.used, scale.mode = if (auto) "auto" else if (K == 1) "single" else "scan", effect.scale = if (!auto && K > 1) effect.scale else NULL,
                       scales = list(scale = scale, beta = sc$scales$beta, median.size = msize, n.centres = sc$scales$n.centres, best = stats::setNames(scale[sc$best], cells), p = sc$P, q = sc$Q, shifts = sc$E, scan.q = sc$q, scan.p = sc$p),
                       group.test = gt, settings = sc$settings, adjust.method = "perm-fdr", p.source = "permutation", n.permutations = sc$B, model = model)
        if (post == "inherit") { ih <- inheritNeighbourQ(shifts$p.adj, shifts$shifts, g.s); shifts$p.adj.cell <- shifts$p.adj; shifts$p.adj <- ih$q; shifts$z.adj <- tz(ih$q); shifts$adj <- shifts$z.adj }
        if (verbose) { n.adj <- sum(shifts$z.adj >= stats::qnorm(0.9), na.rm = TRUE)
          message(sprintf("%d of %d cells keep adjusted z >= %.2f%s", n.adj, sum(is.finite(shifts$z.adj)), stats::qnorm(0.9), if (auto) sprintf(" (scale %g %% chosen)", 100 * scale.used) else "")) }
        self$test.results[[name]] <- shifts
        return(invisible(shifts))
      }
      if (verbose && !is.null(k)) message(sprintf("kNN neighborhoods of at least %d cells (%.1f per sample)", k, k / nlevels(droplevels(as.factor(self$sample.per.cell)))))
      shifts <- clusterFreeExpressionShifts(cm = Matrix::t(inp$cm), sample.per.cell = self$sample.per.cell[rownames(inp$cm)], nns.per.cell = inp$nns.per.cell, nns.graph = inp$nns.graph, post = post,
                                            statistic = statistic, neighbourhood.gene.frac = neighbourhood.gene.frac, centres.per.cell = centres.per.cell, gsd.d0 = gsd.d0,
                                            design = model, meta = model$meta, dist = dist, profile = profile, lib = inp$lib, profile.scale = profile.scale, weighting = weighting,
                                            sample.floor = sample.floor, noise.correction = noise.correction, min.n.obs.per.samp = min.n.obs.per.samp,
                                            min.samp.per.level = min.samp.per.level, permutation = permutation, n.permutations = n.permutations,
                                            block.vars = block.vars %||% model$block.vars, seed = seed %||% sample.int(.Machine$integer.max, 1),
                                            adjust = adjust, smooth = smooth, wins = wins, n.cores = n.cores, verbose = verbose,
                                            robust = robust, robust.k = private$opt("robust.k"), adjust.method = adjust.method, p.source = p.source,
                                            groups = self$cell.groups)
      shifts$model <- model
      self$test.results[[name]] <- shifts
      return(invisible(shifts))
    },

    #' @description Group-level summary of a cluster-free test: one Cauchy-combined p-value per cell group
    #'   (see [summarizeCellScoresByGroup()]) from the per-cell p-values of the differential density or of the
    #'   cluster-free expression shifts, BH-adjusted across groups, next to the median effect and the fraction of
    #'   cells significant after the per-cell adjustment. The combination holds under the dependence between
    #'   neighbouring cells and needs no per-cell adjustment, so it gives a group verdict where the per-cell map is empty.
    #' @param name results slot: `"cell.density"` (after `estimateDiffCellDensity()`) or `"cluster.free.expr.shifts"`
    #' @param cell.groups cell groups (default: stored annotation)
    #' @param alpha significance level (default 0.05)
    #' @param min.cells groups with fewer tested cells are not combined (default 10)
    #' @param type for the density: the `diff` result to summarize (default: the first stored)
    #' @param method `"permutation"` (default: the exact permutation test of the group's mean statistic, computed by the
    #'   estimators when cell groups are stored) or `"cauchy"` (Cauchy combination of the per-cell p-values)
    #' @return data.frame (see [summarizeCellScoresByGroup()]), also stored in `cao$test.results[[name]]$group.summary`
    #' @examples
    #' \dontrun{
    #' cao$estimateDiffCellDensity(); cao$estimateClusterFreeSummary("cell.density")
    #' cao$estimateClusterFreeExpressionShifts(); cao$estimateClusterFreeSummary("cluster.free.expr.shifts")
    #' }
    estimateClusterFreeSummary = function(name=c("cell.density", "cluster.free.expr.shifts"), cell.groups=self$cell.groups, alpha=0.05, min.cells=10, type=NULL,
                                          method=c("permutation", "cauchy")) {
      name <- match.arg(name); method <- match.arg(method)
      if (is.null(cell.groups)) stop("cell.groups are needed")
      if (name == "cell.density") {
        dens.res <- private$getResults(name, 'estimateDiffCellDensity')
        r <- if (is.null(type)) dens.res$diff[[1]] else dens.res$diff[[type]]
        if (is.null(r) || (is.null(r$p.perm) && is.null(r$p.pooled) && is.null(r$p.param))) stop("run estimateDiffCellDensity(type = 'permutation') first")
        p <- r$p.perm %||% r$p.pooled %||% r$p.param; eff <- r$effect; z <- r$pooled %||% r$raw; padj <- r$p.adj; gt <- r$group.test   # diffusion maps store the pooled permutation p only
      } else {
        r <- private$getResults(name, 'estimateClusterFreeExpressionShifts')
        p <- r$p.value; eff <- r$shifts; z <- r$z.pooled; padj <- r$p.adj; gt <- r$group.test
      }
      groups <- as.factor(cell.groups[intersect(names(cell.groups), names(p))])
      if (!is.null(gt) && !identical(sort(unique(as.character(groups))), sort(gt$group[gt$n.tested > 0]))) gt <- NULL   # other groups than the stored test: Cauchy only
      if (method == "permutation" && is.null(gt)) message("no stored group permutation test for these groups: using the Cauchy combination")
      res <- summarizeCellScoresByGroup(p, groups, effect = eff, z = z, p.adj = padj, alpha = alpha, min.cells = min.cells, group.test = gt, method = method)
      res$test <- name
      self$test.results[[name]]$group.summary <- res
      return(invisible(res))
    },

    #' @description Plot the group-level summary of a cluster-free test (see `estimateClusterFreeSummary()`): bars of
    #'   -log10 of the group p-value per cell group, signed by the direction of the median effect for the
    #'   density, filled by significance after BH across groups, labelled with the fraction of cells significant in the
    #'   per-cell map
    #' @param name results slot (see `estimateClusterFreeSummary()`)
    #' @param alpha significance level for the fill (default 0.05)
    #' @param show.frac label the bars with the fraction of significant cells (default TRUE)
    #' @param palette two colours: not significant, significant (default grey and blue)
    #' @return ggplot2 object
    plotClusterFreeSummary = function(name=c("cell.density", "cluster.free.expr.shifts"), alpha=0.05, show.frac=TRUE, palette=c("grey70", "#2a78d6")) {
      name <- match.arg(name)
      res <- self$test.results[[name]]$group.summary
      if (is.null(res)) res <- self$estimateClusterFreeSummary(name, alpha = alpha)
      df <- res[is.finite(res$p), ]
      if (!nrow(df)) stop("no group with enough tested cells")
      sgn <- if (name == "cell.density") ifelse(df$direction == "down", -1, 1) else rep(1, nrow(df))
      df$score <- sgn * pmin(-log10(pmax(df$p, 1e-300)), 300)
      df$group <- factor(df$group, levels = df$group[order(df$score)])
      df$sig <- factor(ifelse(is.finite(df$padj) & df$padj < alpha, "significant", "n.s."), levels = c("n.s.", "significant"))
      df$label <- if (show.frac && any(is.finite(df$frac.sig))) sprintf("%.0f%%", 100 * df$frac.sig) else ""
      lab <- if (identical(res$method[1], "permutation")) "group permutation test" else "Cauchy combination"
      xl <- if (name == "cell.density") sprintf("signed -log10 p (%s) of the density change", lab) else sprintf("-log10 p (%s) of the expression shift", lab)
      gg <- ggplot(df, aes(x = score, y = group, fill = sig)) + geom_col(width = 0.7) +
        geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
        scale_fill_manual(values = stats::setNames(palette, c("n.s.", "significant")), drop = FALSE, name = sprintf("BH across groups, %d%%", round(100 * alpha))) +
        labs(x = xl, y = NULL) + theme_bw() + theme(legend.position = "bottom", panel.grid.minor = element_blank())
      if (show.frac) gg <- gg + geom_text(aes(label = label, hjust = ifelse(score >= 0, -0.1, 1.1)), size = 3, colour = "grey30") +
        scale_x_continuous(expand = expansion(mult = 0.15))
      gg
    },

    #' @description Performs graph smoothing of the cluster-free DE Z-scores
    #' @param n.top.genes numeric Number of top ranked genes to include (default=1000)
    #' @param smoothing `beta` parameter of the \link[sccore:heatFilter]{heatFilter}. (default=20)
    #' @param filter graph filter function. (default=NULL)
    #' @param z.adj boolean Adjust Z scores (default=FALSE)
    #' @param gene.selection character Must be one of "z.adj" or "z", default is based on the "z.adj" parameter (default=ifelse(z.adj, "z.adj", "z"))
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
    #' @param excluded.genes genes to exclude
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
    #' @param gene.selection gene selection method: 'z' or 'z.adj'
    #' @param name results slot name
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
    #' @param name results slot name
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
    #' @param build.panel boolean Plot in a grid (default=TRUE)
    #' @param ncol numeric Number of columns for build.panel (default=3)
    #' @param plots character string (default="z.adj")
    #' @param ... additional parameters passed to plotGeneExpressionComparison()
    #' @return plotGeneExpressionComparison
    #' @examples
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$estimateGenePrograms()
    #' cao$plotGeneProgramGenes(program.id = 1) # program.id is any gene program ID in 1:cao$test.results$gene.programs$n.progs
    #' }
    #' @param name results slot name
    plotGeneProgramGenes = function(program.id, name="gene.programs", ordering=c("similarity", "loading"), max.genes=9, build.panel=TRUE, ncol=3, plots="z.adj", ...) {
      ordering <- match.arg(ordering)
      gene.progs <- private$getResults(name, "estimateGenePrograms")
      if (gene.progs$method == "fabia")
        stop("fabia is deprecated and plotting is not supported anymore")

      scores <- if (ordering == "similarity") gene.progs$sim.scores else gene.progs$loading.scores
      scores %<>% .[[program.id]]
      if (is.null(scores))
        stop("Can't find program", program.id)

      scores %<>% head(max.genes)

      if (build.panel && max.genes == 2) ncol <- 2
      return(self$plotGeneExpressionComparison(scores=scores, plots=plots, ncol=ncol, ...))
    },

    #' @description Plot the cluster-free expression shifts: the shift magnitude next to the signed z-scores
    #' @param cell.groups Indicates cell groups with cell names. Set to NULL if it shouldn't be shown. (default: stored vector)
    #' @param smooth boolean (default=TRUE)
    #' @param plot.na boolean (default=FALSE)
    #' @param scale.z.palette boolean (default=TRUE)
    #' @param min.z (default=qnorm(0.9))
    #' @param color.range range of the shift panel (default=c("0", "97.5%"))
    #' @param z.range symmetric colour range of the z panels: a quantile of the absolute z-scores as a string (default
    #'   `"99.5%"`, so that a few extreme cells do not wash out the map; adjusted z-scores keep at least +-1.96) or a
    #'   numeric vector whose largest absolute value sets the range; values beyond it are clipped
    #' @param max.z cap of a z range taken from a quantile (default 8; a numeric `z.range` is taken as given)
    #' @param z.floor smallest half-range of the z panels (default `qnorm(0.99)` = 2.33): the colours saturate only beyond it, so
    #'   a map whose z-scores barely reach significance stays pale (was 1.96 for the adjusted and 0 for the unadjusted panel)
    #' @param alpha numeric (default=0.2)
    #' @param font.size size range for cell type labels
    #' @param adj.list (default=NULL)
    #' @param build.panel boolean (default=TRUE)
    #' @param ncol columns of the panel (default: all panels in one row); each map keeps the aspect ratio of the embedding
    #' @param ... parameters forwarded to \link[sccore:embeddingPlot]{embeddingPlot}
    #' @examples
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$estimateClusterFreeExpressionShifts()
    #' cao$plotClusterFreeExpressionShifts()
    #' }
    #' @param adjusted `"both"` (default): the unadjusted (pooled-null) and the adjusted z-scores side by side; `TRUE`:
    #'   the adjusted z-scores only; `FALSE`: the unadjusted only
    #' @param name results slot name
    #' @param pal.seq sequential palette for the shift magnitudes
    #' @param pal.div diverging palette for the z-scores
    #' @param outline draw a thin outline around the regions where the cells significant at `outline.alpha` (adjusted p) are,
    #'   on every panel (default FALSE; see [significanceContours()])
    #' @param outline.alpha,outline.level,outline.color,outline.linewidth settings of the outline
    #' @return a panel (`build.panel = TRUE`) or a list of ggplot objects: the shift, then the z-score panel(s)
    plotClusterFreeExpressionShifts = function(cell.groups=self$cell.groups, smooth=TRUE, plot.na=FALSE, adjusted="both",
                                               name="cluster.free.expr.shifts", scale.z.palette=TRUE, min.z=qnorm(0.9),
                                               color.range=c("0", "97.5%"), z.range="99.5%", max.z=8, z.floor=stats::qnorm(0.99), alpha=0.2, font.size=c(3, 5), adj.list=NULL,
                                               pal.seq=brewerPalette("YlOrRd",rev=FALSE), pal.div=brewerPalette("RdBu",rev=TRUE),
                                               build.panel=TRUE, ncol=NULL, outline=FALSE, outline.alpha=0.05, outline.level=0.5, outline.color="black", outline.linewidth=0.25, ...) {
      res <- private$getResults(name, "estimateClusterFreeExpressionShifts")
      private$checkCellEmbedding()
      z.unadj <- res$z.pooled %||% res$z.score   # `z.score` in results of earlier versions
      if (identical(adjusted, "both") && is.null(res$z.adj)) {
        message("adjusted z-scores were not estimated (adjust = FALSE): showing the unadjusted z-scores only")
        adjusted <- FALSE
      }
      panels <- if (identical(adjusted, "both")) list(z.unadj, res$z.adj) else if (isTRUE(adjusted)) list(res$z.adj) else list(z.unadj)
      is.adj <- if (identical(adjusted, "both")) c(FALSE, TRUE) else isTRUE(adjusted)
      leg.titles <- ifelse(is.adj, "z (adj.)", "z (pooled)")
      titles <- ifelse(is.adj, "Adjusted z-score", "z-score (pooled null)")

      shifts <- if (smooth) res$shifts.smoothed %||% res$shifts else res$shifts   # diffusion neighbourhoods are averages already (no smoothed copy)
      shifts %<>% na.omit()
      color.range %<>% parseLimitRange(shifts)
      shifts %<>% pmax(color.range[1]) %>% pmin(color.range[2])
      gg.shift <- self$plotEmbedding(colors=shifts, plot.na=plot.na, alpha=alpha, palette=pal.seq, legend.title="shift", ...) +
        theme(legend.background = element_blank())
      for (aes.name in c("colour", "fill")) { sc <- gg.shift$scales$get_scales(aes.name); if (!is.null(sc) && inherits(sc, "ScaleContinuous")) sc$n.breaks <- 3 }   # narrow ranges crowd the labels

      gg.z <- mapply(function(z, adj, lt) {
        degenerate <- !any(abs(z) >= min.z, na.rm = TRUE)
        if (degenerate) message(sprintf("no cell reaches |z| >= %.2f (%s z-scores): the z-score panel is uniform", min.z, if (adj) "adjusted" else "unadjusted"))
        cr <- if (degenerate) c(-1, 1) * max(min.z, z.floor) else symmetricColorRange(z, z.range, floor = z.floor, cap = max.z)   # the range reaches at least z.floor, so that colours saturate only beyond it   # a uniform panel has no range to take
        z %<>% pmin(cr[2]) %>% pmax(cr[1])
        gg <- self$plotEmbedding(colors=z, plot.na=plot.na, alpha=alpha, palette=pal.div, legend.title=lt, ...) +
          theme(legend.background = element_blank())
        if (scale.z.palette) {
          gg$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
          # one-sided adjusted z-scores are non-negative: use the upper half of the range only
          gg <- gg + getScaledZGradient(min.z = min.z, palette = pal.div, color.range = if (any(z < 0, na.rm = TRUE)) cr else cr[2])
        }
        if (degenerate) gg <- gg + theme(legend.position = "none")
        gg
      }, panels, is.adj, leg.titles, SIMPLIFY=FALSE)
      ggs <- c(list(gg.shift), gg.z)
      if (outline && !is.null(res$p.adj)) {   # the regions where the cells significant at `outline.alpha` are
        ol <- significanceContourLayer(self$embedding, stats::setNames(res$p.adj < outline.alpha, names(res$p.adj)), level = outline.level, color = outline.color, linewidth = outline.linewidth)
        if (!is.null(ol)) ggs %<>% lapply(`+`, ol)
      }

      if (!is.null(cell.groups)) {
        ggs %<>% lapply(transferLabelLayer, self$plotEmbedding(groups=cell.groups), font.size=font.size)
      }
      if (!is.null(adj.list)) ggs %<>% lapply(`+`, adj.list)
      if (build.panel) {
        ggs[[1]] <- ggs[[1]] + ggtitle("Expression shift")
        for (i in seq_along(gg.z)) ggs[[i + 1]] <- ggs[[i + 1]] + ggtitle(titles[i])
        # panels keep the aspect ratio of the embedding; legends below, so that the maps get the width
        emb <- self$embedding; ar <- if (!is.null(emb) && ncol(emb) >= 2) diff(range(emb[, 2], na.rm = TRUE)) / diff(range(emb[, 1], na.rm = TRUE)) else NULL
        ggs %<>% lapply(function(g) g + theme(legend.position = "bottom", legend.key.height = grid::unit(0.3, "cm"), legend.title = element_text(vjust = 0.8),
                                             aspect.ratio = if (is.finite(ar %||% NA)) ar else NULL))
        ggs %<>% cowplot::plot_grid(plotlist=., ncol=ncol %||% length(ggs), align = "hv")
      }

      return(ggs)
    },
                                            
    #' @description Plot most changed genes
    #' @param n.top.genes numeric
    #' @param method character string (default='z')
    #' @param min.z numeric (default=0.5)
    #' @param min.lfc numeric (default=1)
    #' @param max.score numeric (default=20)
    #' @param cell.subset (default=NULL)
    #' @param excluded.genes (default=NULL)
    #' @param build.panel boolean Plot in grid (default=TRUE)
    #' @param ncol numeric Number of columns for build.panel (default=1)
    #' @param ... additional parameters input to self$plotGeneExpressionComparison()
    #' @return plot of the most changed genes via plotGeneExpressionComparison()
    #' @examples
    #' \dontrun{
    #' cao$estimateClusterFreeDE()
    #' cao$plotMostChangedGenes(n.top.genes = 10) # n.top.genes is any number of genes to plot
    #' }
    plotMostChangedGenes = function(n.top.genes, method="z", min.z=0.5, min.lfc=1, max.score=20, cell.subset=NULL, excluded.genes=NULL, build.panel=TRUE, ncol = 1, ...) {
      scores <- self$getMostChangedGenes(n.top.genes, method=method, min.z=min.z, min.lfc=min.lfc, max.score=max.score,
                                         cell.subset=cell.subset, excluded.genes=excluded.genes)
      self$plotGeneExpressionComparison(scores=scores, cell.subset=cell.subset, build.panel=build.panel, ncol=ncol, ...)
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
    #' @param nrow numeric Number of rows for build.panel (default=1)
    #' @param ncol numeric Number of columns for build.panel (default = 1)
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
    #' @param name results slot name
    plotGeneExpressionComparison=function(genes=NULL, scores=NULL, max.expr="97.5%", plots=c("z.adj", "z", "expression"),
                                          min.z=qnorm(0.9), max.z=4, max.z.adj=NULL, max.lfc=3, smoothed=FALSE,
                                          gene.palette=dark.red.palette, z.palette=NULL, z.adj.palette=z.palette,
                                          lfc.palette=NULL, scale.z.palette=TRUE, plot.na=-1, adj.list=NULL,
                                          build.panel=TRUE, nrow=1, ncol=1, cell.subset=NULL, groups=NULL, subgroups=NULL,
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

      if (build.panel && length(ggs) > 1) ggs <- cowplot::plot_grid(plotlist = ggs, ncol = ncol)

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

    # per-"cell type" distance matrices for a space: expression shifts (cached pseudobulk + distances) or composition
    distancesForSpace = function(space, dist, n.pcs = NULL, cell.groups = self$cell.groups, min.cells.per.sample = 10, min.gene.frac = 0.01,
                                 genes = NULL, verbose = FALSE) {
      space <- normalizeSampleSpace(space)
      if (space == "pseudobulk") {   # all cells of a sample pooled into one profile
        key <- settingsKey(names(self$sample.per.cell), as.character(self$sample.per.cell), min.cells.per.sample, min.gene.frac, genes, dist, n.pcs)
        hit <- self$cache$pooled.pseudobulk
        if (!is.null(hit) && identical(hit$key, key)) return(hit$value)
        if (verbose) message("Building pooled pseudobulk profiles...")
        cms <- extractRawCountMatrices(self$data.object, transposed = TRUE)
        cms <- cms[intersect(names(cms), rownames(self$sample.meta))]
        cg1 <- setNames(factor(rep("pseudobulk", length(self$sample.per.cell))), names(self$sample.per.cell))
        pb <- pseudobulkPerCellType(cms, cell.groups = cg1, sample.per.cell = self$sample.per.cell, min.cells.per.sample = min.cells.per.sample,
                                    min.gene.frac = min.gene.frac, genes = genes)
        D.list <- Filter(function(D) !is.null(D) && nrow(D) >= 3, sampleDistanceMatrices(pb$cm.per.type, dist = dist, n.pcs = n.pcs))
        for (ct in names(D.list)) attr(D.list[[ct]], "n.cells") <- setNames(as.numeric(pb$n.cells[ct, rownames(D.list[[ct]])]), rownames(D.list[[ct]]))
        self$cache$pooled.pseudobulk <- list(key = key, value = D.list)
        return(D.list)
      }
      if (space == "composition") {
        checkPackageInstalled("coda.base", cran = TRUE)
        cnts <- private$extractCodaData(cells.to.remove = NULL, cells.to.remain = NULL, samples.to.remove = NULL)
        cnts <- cnts[, colSums(cnts, na.rm = TRUE) > 0, drop = FALSE]
        ilr <- computeILRMatrix(cnts)$ilr
        D <- as.matrix(stats::dist(ilr)); attr(D, "n.cells") <- setNames(rowSums(cnts)[rownames(D)], rownames(D))
        return(list(composition = D))
      }
      pb <- private$getPseudobulk(cell.groups = cell.groups, sample.per.cell = self$sample.per.cell, min.cells.per.sample = min.cells.per.sample,
                                  min.gene.frac = min.gene.frac, genes = genes, verbose = verbose)
      private$getSampleDistances(pb, dist = dist, n.pcs = n.pcs, verbose = verbose)
    },

    # gene x cell counts of the given genes and per-cell library sizes over all genes (cell-type gsd), cached on the genes
    getGsdCellTypeInput = function(genes) {
      key <- settingsKey(genes)
      hit <- self$cache$gsd.celltype.input
      if (!is.null(hit) && identical(hit$key, key)) return(hit$value)
      cms <- extractRawCountMatrices(self$data.object, transposed = TRUE)
      cms <- cms[intersect(names(cms), rownames(self$sample.meta))]
      lib <- unlist(lapply(unname(cms), Matrix::rowSums))
      cm <- Matrix::t(do.call(rbind, lapply(unname(cms), function(x) sccore::extendMatrix(x, genes)[, genes, drop = FALSE])))
      cm <- methods::as(cm, "CsparseMatrix")
      val <- list(cm = cm, lib = lib[colnames(cm)])
      self$cache$gsd.celltype.input <- list(key = key, value = val)
      val
    },

    # pseudobulk per cell type, cached on the settings that produce it
    getPseudobulk = function(cell.groups, sample.per.cell, min.cells.per.sample, min.gene.frac, genes, verbose = FALSE) {
      key <- settingsKey(as.character(cell.groups), names(cell.groups), as.character(sample.per.cell), min.cells.per.sample, min.gene.frac, genes)
      hit <- self$cache$pseudobulk
      if (!is.null(hit) && identical(hit$key, key)) return(hit$value)
      if (verbose) message("Building pseudobulk profiles per cell type...")
      cms <- extractRawCountMatrices(self$data.object, transposed = TRUE)
      cms <- cms[intersect(names(cms), rownames(self$sample.meta))]
      pb <- pseudobulkPerCellType(cms, cell.groups = cell.groups, sample.per.cell = sample.per.cell,
                                  min.cells.per.sample = min.cells.per.sample, min.gene.frac = min.gene.frac, genes = genes)
      self$cache$pseudobulk <- list(key = key, value = pb)
      pb
    },

    # per-cell-type sample distances, each on the type's genes (typeGeneSets(): detected in >= type.gene.frac of its cells)
    getSampleDistances = function(pb, dist, n.pcs = NULL, type.gene.frac = NULL, verbose = FALSE) {
      type.gene.frac <- private$opt("type.gene.frac", type.gene.frac) %||% 0
      key <- settingsKey(self$cache$pseudobulk$key, dist, n.pcs, type.gene.frac)
      hit <- self$cache$sample.distances
      if (!is.null(hit) && identical(hit$key, key)) return(hit$value)
      if (verbose) message(sprintf("Computing sample distances (%s%s%s)...", dist, if (!is.null(n.pcs)) sprintf(", %d PCs", n.pcs) else "",
                                   if (type.gene.frac > 0) sprintf(", per cell type the genes detected in >= %g %% of its cells", 100 * type.gene.frac) else ""))
      gs <- typeGeneSets(pb, type.gene.frac)
      D.list <- sampleDistanceMatrices(lapply(sccore::sn(names(pb$cm.per.type)), function(ct) pb$cm.per.type[[ct]][, gs[[ct]], drop = FALSE]),
                                       dist = dist, n.pcs = n.pcs)
      for (ct in names(D.list)) if (!is.null(D.list[[ct]])) attr(D.list[[ct]], "n.genes") <- length(gs[[ct]])
      D.list <- Filter(function(D) !is.null(D) && nrow(D) >= 3, D.list)
      for (ct in names(D.list)) attr(D.list[[ct]], "n.cells") <- setNames(as.numeric(pb$n.cells[ct, rownames(D.list[[ct]])]), rownames(D.list[[ct]]))
      self$cache$sample.distances <- list(key = key, value = D.list)
      D.list
    },

    # arguments of the previous expression-shift implementation: translate what has an equivalent, announce the rest
    legacyShiftArgs = function(dots, verbose = TRUE) {
      out <- list()
      if (!length(dots)) return(out)
      nms <- names(dots)
      if (is.null(nms) || any(!nzchar(nms))) stop("unnamed extra arguments are not accepted")
      if (any(c("pairFormula", "pairContrast") %in% nms))
        stop("`pairFormula` / `pairContrast` no longer exist: the model is specified at the sample level (formula + test). ",
             "Pair-level terms are derived from it; see ?Cacoa for `setModel()` and the `test` grammar.", call. = FALSE)
      if ("perm.method" %in% nms) { out$permutation <- dots$perm.method; if (verbose) message("`perm.method` is deprecated: use `permutation`") }
      if ("dist.type" %in% nms && verbose) message("`dist.type` is deprecated: shift, var and total are all reported; choose the effect when plotting")
      if ("min.samp.per.type" %in% nms && verbose) message("`min.samp.per.type` is deprecated: use `min.samp.per.level` (samples per compared level)")
      known <- c("perm.method", "dist.type", "pairFormula", "pairContrast", "alternative", "robust.method", "na.mode", "return.residuals",
                 "return.sampled.stats", "return.sampled.fits", "top.n.genes", "gene.selection", "sample.metadata", "sample.ids",
                 "min.samp.per.type", "p.adjust.method", "trim", "keep.all")
      ignored <- setdiff(intersect(nms, known), c("perm.method", "dist.type", "min.samp.per.type"))
      if (length(ignored) && verbose) message("Ignoring deprecated argument(s): ", paste(ignored, collapse = ", "))
      unknown <- setdiff(nms, known)
      if (length(unknown)) stop("unknown argument(s): ", paste(unknown, collapse = ", "))
      out
    },

    # option lookup: an explicit (non-NULL) argument wins over the stored options
    opt = function(name, value = NULL) {
      if (!is.null(value)) return(value)
      if (is.null(self$options)) self$options <- cacoaDefaultOptions()
      if (!name %in% names(self$options)) return(cacoaDefaultOptions()[[name]])   # an option added after the object was created
      self$options[[name]]
    },

    # The model a method should use: the stored one, or a temporary one when the call gives formula / test /
    # contrast (announced with a message; the stored model is left unchanged).
    resolveModel = function(formula = NULL, test = NULL, contrast = NULL, dispersion.formula = NULL, block.vars = NULL,
                            permutation = NULL, verbose = TRUE, what = "this call") {
      if (is.null(formula) && is.null(test) && is.null(contrast) && is.null(dispersion.formula) && is.null(block.vars)) {
        if (is.null(self$model)) stop("no model is set: call cao$setModel(test = \"<variable>\") or pass `test` / `formula` to ", what)
        if (!grepl("^estimateCellLoadings", what) && !is.null(msg <- groupingBetweenMessage(self$model, what))) stop(msg, call. = FALSE)
        return(self$model)
      }
      m <- self$model
      if (is.null(test) && is.null(contrast)) {
        if (is.null(m)) stop("`formula` given without `test`: say what to test (e.g. test = \"condition\")")
        test <- m$test.spec$test; contrast <- m$test.spec$contrast
      }
      if (is.null(formula) && !is.null(m)) {
        tv <- unique(stats::na.omit(vapply(resolveTests(test, contrast, self$sample.meta, m$formula, private$opt("numeric.step")),
                                           function(t) t$variable %||% NA_character_, character(1))))
        if (all(tv %in% all.vars(m$formula))) formula <- m$formula.full %||% m$formula
      }
      tmp <- private$buildModel(formula = formula, test = test, contrast = contrast,
                                dispersion.formula = dispersion.formula %||% (if (!is.null(m) && is.null(test) && is.null(contrast)) m$dispersion.formula else NULL),
                                block.vars = block.vars %||% self$block.vars, permutation = permutation)
      if (verbose) message("Using a temporary model for ", what, ": ", paste(deparse(tmp$formula), collapse = ""), ", test ",
                           paste(vapply(tmp$tests, `[[`, character(1), "label"), collapse = "; "),
                           ". The stored model is unchanged (cao$setModel() changes it).")
      if (!grepl("^estimateCellLoadings", what) && !is.null(msg <- groupingBetweenMessage(tmp, what))) stop(msg, call. = FALSE)
      tmp
    },

    buildModel = function(formula = NULL, test = NULL, contrast = NULL, dispersion.formula = NULL, block.vars = NULL,
                          permutation = NULL, numeric.ref = NULL, numeric.step = NULL) {
      m <- buildCacoaModel(self$sample.meta, formula = formula, test = test, contrast = contrast, dispersion.formula = dispersion.formula,
                           block.vars = block.vars, numeric.ref = numeric.ref %||% self$numeric.ref,
                           numeric.step = private$opt("numeric.step", numeric.step), permutation = private$opt("permutation", permutation),
                           n.permutations = private$opt("n.permutations"))
      m$test.spec <- list(test = test, contrast = contrast)
      err <- m$issues[m$issues$severity == "error", , drop = FALSE]
      if (nrow(err)) stop(paste(c(err$message, if (nzchar(err$suggestion[1])) paste("Suggestion:", err$suggestion[1])), collapse = "\n"), call. = FALSE)
      m
    },

    # the variable of the primary test (first metadata column when no model is set)
    # one sample-space distance matrix computed from the data (no stored results needed): attributes dist.type, n.cells
    sampleSpaceDistance = function(space, cell.type = NULL, joint.weighting = "cells", dist = NULL, verbose = FALSE) {
      dist <- private$opt("dist", dist)
      if (space == "composition") {
        D <- private$distancesForSpace("composition", dist = "l2")[[1]]; attr(D, "dist.type") <- "l2"; return(D)
      }
      if (space == "pseudobulk") {
        D <- private$distancesForSpace("pseudobulk", dist = dist, verbose = verbose)[[1]]
        if (is.null(D)) stop("too few samples with at least 10 cells for the pooled pseudobulk")
        attr(D, "dist.type") <- dist; return(D)
      }
      D.list <- private$distancesForSpace("expression", dist = dist, verbose = verbose)
      D <- private$sampleSpaceFromList(D.list, NULL, cell.type, joint.weighting)
      if (!is.null(D)) attr(D, "dist.type") <- dist
      D
    },

    # one cell type's matrix, or the combination over cell types, from a per-type list (n.cells from `res` when missing)
    sampleSpaceFromList = function(mats, res, cell.type, joint.weighting) {
      if (!is.null(cell.type)) {
        if (!cell.type %in% names(mats)) { warning("Distances were not estimated for cell type ", cell.type); return(NULL) }
        D <- mats[[cell.type]]
        if (is.null(attr(D, "n.cells")) && !is.null(res)) attr(D, "n.cells") <- setNames(as.numeric(res$n.cells[cell.type, rownames(D)]), rownames(D))
        return(D)
      }
      mats <- mats[vapply(mats, function(M) is.matrix(M) && nrow(M) >= 2, logical(1))]
      if (!length(mats)) return(NULL)
      for (ct in names(mats)) if (is.null(attr(mats[[ct]], "n.cells")))
        attr(mats[[ct]], "n.cells") <- setNames(as.numeric(res$n.cells[ct, rownames(mats[[ct]])]), rownames(mats[[ct]]))
      D <- prepareJointExpressionDistance(mats, weighting = joint.weighting)
      if (any(is.na(D))) D <- ape::additive(D) %>% `dimnames<-`(dimnames(D))   # pairs of samples that share no cell type
      nc <- table(self$sample.per.cell); attr(D, "n.cells") <- setNames(as.numeric(nc[rownames(D)]), rownames(D))
      D
    },

    # the model's covariates other than the tested variable, as a formula (NULL when there are none)
    modelNuisanceFormula = function() {
      if (is.null(self$model)) return(NULL)
      vs <- setdiff(all.vars(self$model$formula), private$testVariable())
      vs <- intersect(vs, names(self$sample.meta))
      if (!length(vs)) return(NULL)
      stats::as.formula(paste("~", paste(vs, collapse = " + ")))
    },

    # subtitle of a sample map: chance-corrected R2 and permutation p of each covariate on the shown distances
    sampleSpaceR2Label = function(D, meta, covariates, titles = covariates, n.permutations = 999, n.cores = 1, adjusted = character(0)) {
      if (!length(covariates)) return(NULL)
      meta <- meta[rownames(D), , drop = FALSE]
      out <- character(0)
      for (i in seq_along(covariates)) {
        v <- covariates[i]
        if (v %in% adjusted) { out[i] <- sprintf("%s: removed by the adjustment", titles[i]); next }
        sc <- tryCatch(suppressWarnings(screenCovariates(list(map = D), meta, covariates = v, mode = "marginal", dist = attr(D, "dist.type") %||% "l2",
                                                         n.permutations = n.permutations, min.samples.per.type = 4, seed = private$opt("seed") %||% 1,
                                                         n.cores = 1, verbose = FALSE)), error = function(e) NULL)
        r <- if (!is.null(sc)) sc$table[sc$table$covariate == v, ] else NULL
        out[i] <- if (is.null(r) || !nrow(r) || !is.finite(r$R2.adj[1])) sprintf("%s: not testable", titles[i]) else
          sprintf("%s: adj. R\u00b2 %.2f, p %.2g", titles[i], r$R2.adj[1], r$p[1])
      }
      paste(out, collapse = "; ")
    },

    testVariable = function() {
      v <- if (!is.null(self$model)) self$model$tests[[1]]$variable else NULL
      if (is.null(v) || is.na(v)) colnames(self$sample.meta)[1] else v
    },

    # legacy two-group fields (ref.level, target.level, sample.groups, palette names) from the primary test
    syncLegacyFields = function() {
      m <- self$model
      self$ref.level <- self$target.level <- self$sample.groups <- NULL
      t <- m$tests[[1]]
      if (t$kind == "contrast" && is.null(t$step) && !is.null(t$variable) && !is.na(t$variable) && all(!is.na(t$levels))) {
        self$ref.level <- unname(t$levels[["ref"]]); self$target.level <- unname(t$levels[["alt"]])
        self$sample.groups <- self$getSampleGroups()
      }
      pal <- self$sample.groups.palette
      if (!is.null(self$target.level) && (is.null(pal) || is.null(names(pal)) || !all(c(self$target.level, self$ref.level) %in% names(pal)))) {
        if (is.null(pal) || length(pal) < 2) pal <- c("#d73027", "#4575b4")
        self$sample.groups.palette <- setNames(unname(pal[1:2]), c(self$target.level, self$ref.level))
      }
      invisible(NULL)
    },

    getTopGenes = function(n, gene.selection=c("expression", "od"), cm.joint=NULL,
                           min.expr.frac=0.0, excluded.genes=NULL, included.genes=NULL, name="cluster.free.de", ...) {
      gene.selection <- match.arg(gene.selection)

      if (min.expr.frac > 0) {
        if (is.null(cm.joint)) {
          cm.joint <- self$getJointCountMatrix()
        }

        cm.joint@x <- 1 * (cm.joint@x > 0)
        excluded.genes %<>% union(colnames(cm.joint)[colMeans(cm.joint, na.rm=TRUE) < min.expr.frac])
      }

      if (gene.selection == "od") {
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
    extractCodaData=function(ret.groups=FALSE, cell.groups=self$cell.groups, cells.to.remove=NULL, cells.to.remain=NULL,
                             samples.to.remove=NULL) {
      d.counts <- cell.groups %>% data.frame(anno=., group=self$sample.per.cell[names(.)]) %>%
        table() %>% rbind() %>% t()
      
      if (!is.null(cells.to.remove)) d.counts %<>% .[,!(colnames(.) %in% cells.to.remove), drop=FALSE]
      if (!is.null(cells.to.remain)) d.counts %<>% .[,colnames(.) %in% cells.to.remain, drop=FALSE]
      if (!is.null(samples.to.remove)) d.counts %<>% .[!(rownames(.) %in% samples.to.remove), , drop=FALSE]
      
      s.order <- intersect(rownames(self$sample.meta), rownames(d.counts))
      d.counts <- d.counts[s.order, , drop=FALSE]
      
      if (!ret.groups) {
        return(d.counts)
      }
      
      # grouping variable and target level of the primary test
      condition <- private$testVariable()
      target <- self$target.level
      
      sample.groups <- setNames(self$sample.meta[[condition]], rownames(self$sample.meta))
      
      d.groups <- (sample.groups[rownames(d.counts)] == target) %>%
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
    ## kNN neighbourhood lists of the cell graph for a size k, computed once per object and size (both cluster-free maps use them)
    cellGraph = function() {   # the cell graph of the analyses: without the cells of excluded samples (excludeSamples())
      g <- extractCellGraph(self$data.object); spc <- self$sample.per.cell
      if (!anyNA(spc)) return(g)
      igraph::induced_subgraph(g, intersect(igraph::V(g)$name, names(spc)[!is.na(spc)]))
    },
    knnUnits = function(k, n.cores = 1) {
      g <- private$cellGraph(); cells <- igraph::V(g)$name; key <- sprintf("k%d_n%d", as.integer(k), length(cells))
      if (is.null(private$knn.cache[[key]])) private$knn.cache[[key]] <- knnNeighbourhoodLists(g, k, cells = cells, n.cores = n.cores)
      private$knn.cache[[key]]
    },
    knn.cache = list(),
    all.samples = NULL,   # the full sample assignment and metadata, kept by excludeSamples()
    kernelBetaForScale = function(graph, scale, n.cores = 1) {   # kernel strength of a diffusion neighborhood at a scale (label-free, cached per graph and scale)
      key <- sprintf("s%g_n%d_e%d", scale, igraph::vcount(graph), igraph::ecount(graph))
      if (is.null(private$knn.cache[[key]])) { sb <- solveKernelBeta(graph, target.size = round(scale * igraph::vcount(graph)), n.vec = 16, n.cores = n.cores)
        private$knn.cache[[key]] <- list(beta = sb$beta, median.size = stats::median(sb$size), self = sb$self, off.size = sb$off.size) }   # the probes' kernel diagonal and arriving-heat size: the cell as one cell of its neighborhood (heatLogDensity())
      private$knn.cache[[key]]
    },
    getClusterFreeDEInput = function(genes, raw=FALSE, min.edge.weight=0.0, k=NULL) {
      cm <- self$getJointCountMatrix(raw=raw)
      genes <- intersect(genes, colnames(cm))
      is.ref <- if (!is.null(self$sample.groups)) (as.character(self$sample.groups[levels(self$sample.per.cell)]) == self$ref.level) else rep(NA, nlevels(self$sample.per.cell))

      adj.mat <- private$cellGraph() %>% igraph::as_adj()
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
      lib <- Matrix::rowSums(cm[cell.names, , drop=FALSE])   # library size over all genes
      cm %<>% .[cell.names, genes, drop=FALSE]

      for (id in which(sapply(nns.per.cell, is.null))) {
         # If there are cells without neighbors. May be a consequence of graph modification
        nns.per.cell[[id]] <- numeric()
      }
      nns.graph <- nns.per.cell
      if (!is.null(k)) {   # kNN unit of at least k cells per cell (self, every graph neighbour, then two-hop by shared neighbours); the 1-hop graph is kept for post-processing
        nn <- private$knnUnits(k)
        if (!identical(names(nn), cell.names)) { idx <- match(cell.names, names(nn)); if (anyNA(idx)) stop("the kNN units do not cover the tested cells")
          remap <- match(names(nn), cell.names); nns.per.cell <- lapply(nn[idx], function(v) { w <- remap[v + 1L]; as.integer(w[!is.na(w)] - 1L) }) } else nns.per.cell <- nn
        names(nns.per.cell) <- cell.names
      }

      return(list(cm=cm, lib=lib, adj.mat=adj.mat, is.ref=is.ref, nns.per.cell=nns.per.cell, nns.graph=nns.graph))
    }
  )
)
