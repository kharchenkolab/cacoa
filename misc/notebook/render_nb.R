## Render the cacoa walkthrough as a Jupyter notebook (.ipynb) with executed outputs, plus a Markdown copy (nb_to_md.R).
## usage: R_LIBS=<lib> Rscript render_nb.R vignettes/walkthrough_short.ipynb <path to shifts_sim_objects.rds>
## Code cells are evaluated with `evaluate`; text goes to stream outputs, plots to PNG display_data.
suppressPackageStartupMessages({library(cacoa); library(conos); library(ggplot2); library(cowplot)})
args <- commandArgs(TRUE); out.file <- args[1]; sim.file <- args[2]
cells <- list()
md <- function(...) cells[[length(cells) + 1]] <<- list(type = "markdown", text = paste0(..., collapse = ""))
code <- function(text, fig = c(8, 4.5)) cells[[length(cells) + 1]] <<- list(type = "code", text = text, fig = fig)

md("# Getting started with Cacoa\n\n",
   "Cacoa compares samples (patients, donors, replicates) across conditions. One **model** is set once and used by every ",
   "analysis: expression shifts per cell type, composition, differential expression, cell density and the cluster-free analyses.\n\n",
   "The recommended order is\n\n",
   "1. `Cacoa$new()` with the sample metadata. No model is needed yet; a summary of the metadata is printed.\n",
   "2. `cao$checkDesign()`, `cao$plotDesign()`: confounding, balance, degrees of freedom and permutation feasibility, from the metadata alone.\n",
   "3. `cao$screenCovariates()`, `cao$plotCovariateScreen()`: which covariates are associated with sample-level variation, and in which cell types (exploratory).\n",
   "4. `cao$setModel(~ condition + batch, test = \"condition\")`: the model used from here on.\n",
   "5. `cao$estimateExpressionShiftMagnitudes()`, `cao$plotExpressionShiftMagnitudes()` and the other analyses.\n",
   "6. `cao$checkSensitivity()`, `cao$plotSensitivity()`: is the result robust to the covariate set and to single samples?\n\n",
   "This notebook runs the workflow on a **simulated** dataset: a Conos object with 40 samples (4 groups x 2 batches, 5 samples each) and ",
   "9 cell types of 500 cells. Expression differences between `Group2` and `Group1` were planted with a known strength per cell type: ",
   "strong in `IN-SST`, `IN-PV` and `L2/3`; moderate in `OPC`, `L5/6-CC` and `AST-PP`; weak in `AST-FB`; none in `Microglia` and `Neu-NRGN`. ",
   "The object is not shipped with the package (720 MB); the notebook is provided for reading, with its outputs. ",
   "The same code runs on any Conos or Seurat object with a sample metadata table.")
md("> GitHub's notebook preview runs its math renderer over code cells, so R code containing `$` is displayed incorrectly there. For reading on GitHub use [walkthrough_short.md](walkthrough_short.md), the same notebook rendered as Markdown. The `.ipynb` is for Jupyter, VS Code and nbviewer.")
code("library(conos)\nlibrary(cacoa)\nlibrary(ggplot2)\nlibrary(cowplot)")

md("## Load data\n\n",
   "`sample.metadata` is a data frame with one row per sample (row names are the sample names) and one column per covariate. ",
   "Here it has the group, the batch and two technical columns from the simulation. A random `age` column is added to show how an unrelated covariate is treated.")
code(sprintf('sim <- readRDS("%s")$with_batch                      # a Cacoa object from an earlier version
con <- sim$data.object                                 # the Conos object inside it
sample.meta <- con$misc$sample_meta
set.seed(1); sample.meta$age <- round(rnorm(nrow(sample.meta), 55, 8))   # unrelated to the groups, for illustration
head(sample.meta)', sim.file))
code(sprintf('cao <- Cacoa$new(con, sample.metadata = sample.meta, cell.groups = sim$cell.groups, sample.per.cell = sim$sample.per.cell, n.cores = %d)
cao$plot.params <- list(size = 0.3, alpha = 0.3, font.size = c(2, 3))', min(16L, parallel::detectCores())))
md("The constructor prints which metadata columns can be used as covariates and which cannot (ID-like, constant, mostly missing, too many levels). ",
   "The previous two-group constructor (`sample.groups`, `ref.level`, `target.level`) still works and is translated into a `condition` column and `test = \"condition\"`.\n\n",
   "Options shared by all methods are set once; explicit arguments of a method call still win.")
code("cao$setOptions(n.permutations = 499, seed = 1)")

md("## Check the design\n\n",
   "Before fitting anything, the metadata alone tells whether the planned comparison can be adjusted for the other covariates: ",
   "associations between covariates (bias-corrected Cramer's V, correlation ratio, Spearman correlation), balance of the groups over the covariates, aliasing, ",
   "degrees of freedom and how many distinct permutations the design allows.")
code('chk <- cao$checkDesign(test = "Group")\nchk')
code('cao$plotDesign("associations")', fig = c(5, 4))
code('cao$plotDesign("balance")', fig = c(8, 3.5))
md("The report lists errors (a covariate that fully determines the condition), warnings (strongly associated covariates, small groups, few distinct permutations) and notes, each with a suggestion. Here the design is balanced and nothing is flagged.")

md("## Screen covariates\n\n",
   "The screen tests every usable covariate against the sample-to-sample expression distances of each cell type, ",
   "**marginally** (alone) and **partially** (adjusted for all other covariates). ",
   "A covariate that is significant only marginally is explained by the others (shown as a ring in the plot). ",
   "The last column is a global test over all cell types. The screen ends with a suggested model; it is not applied automatically.")
code('scr <- cao$screenCovariates(n.permutations = 199)\nscr')
code('cao$plotCovariateScreen()', fig = c(8, 4))
code('cao$plotCovariateSummary()', fig = c(7, 3.5))
md("Variance partition: how much of the sample-level variation is explained uniquely by each covariate, shared between them, or left unexplained.")
code('cao$plotVariancePartition(c("Group", "Batch"))', fig = c(8, 3.5))

md("## Set the model\n\n",
   "`test` accepts a variable name (two levels: a contrast with an automatically chosen reference; numeric: a slope per unit; a factor with more levels needs the comparison, e.g. `\"Group: Group2 vs Group1\"`, or `\"Group: all\"` for a whole-factor test), ",
   "an explicit comparison such as `\"Group: Group2 vs Group1\"`, several variables, `\"all\"`, or a structured contrast. ",
   "The printout states the reference level and how it was chosen, the adjustment set, the permutation scheme and the number of distinct permutations.\n\n",
   "Here two tests are set: the planted contrast `Group2 vs Group1`, and the whole-factor test of `Group` over all four groups. Both are adjusted for `Batch`.")
code('cao$setModel(~ Group + Batch, test = c("Group: Group2 vs Group1", "Group: all"))')

md("## Expression shifts per cell type\n\n",
   "The expression shift is estimated with an individual-level model fitted to the sample-sample distances of each cell type. Three effects are reported for a contrast:\n\n",
   "| effect | meaning |\n|---|---|\n",
   "| `shift` | the groups differ in a common direction, beyond within-group variability |\n",
   "| `var`   | the target group is more (or less) heterogeneous than the reference group |\n",
   "| `total` | target samples are farther from reference samples than reference samples are from each other (`shift + var/2`) |\n\n",
   "For a whole-factor test the effects are `location` (any group differs) and `dispersion` (the groups differ in heterogeneity). ",
   "P-values come from permutations (block randomization within batches here); they are adjusted across cell types by Benjamini-Hochberg (`padj`) and by the max-statistic (`p.fwer`), and a global p-value per test is reported.")
code('res <- cao$estimateExpressionShiftMagnitudes()\nsubset(res$results, effect == "shift")[, c("celltype", "estimate", "estimate.norm", "se.jk", "p", "padj", "p.fwer", "n")]')
code('res$global')
code('cao$plotExpressionShiftMagnitudes()', fig = c(9, 4.5))
md("Points are normalized effects with jackknife intervals; filled points are significant after Benjamini-Hochberg adjustment across cell types. ",
   "The subtitle records the test, the adjustment set and the permutation scheme. The planted strong cell types come out on top and the two unchanged ones at the bottom.\n\n",
   "The second test (all four groups) is plotted by name.")
code('cao$plotExpressionShiftMagnitudes(test = "Group")', fig = c(9, 4.5))
md("### Behind one cell type\n\n`plotShiftDetail()` shows the sample-level picture: the distance distribution within and between groups, and the samples in the adjusted distance space.")
code('cao$plotShiftDetail("IN-SST")', fig = c(10, 4))
md("`plotSampleInfluence()` shows how much each sample changes the shift estimate when it is left out.")
code('cao$plotSampleInfluence()', fig = c(9, 4.5))

md("## Robustness\n\n",
   "The expression shifts are re-estimated without adjustment, with the current model, without each covariate, and with the top covariates proposed by the screen. ",
   "Per cell type, one verdict summarizes the result: robust, sign or significance depends on the adjustment, the estimate changes by more than half, or the result is driven by one sample.")
code('sens <- cao$checkSensitivity(n.permutations = 199)\nsens')
code('cao$plotSensitivity()', fig = c(9, 5))

md("## Sample structure\n\n",
   "The sample distances of one cell type (or of all cell types jointly) can be shown as is, or after removing the part explained by a covariate (`adjust.for`).")
code('plot_grid(
  cao$plotSampleDistances(space = "expression.shifts", cell.type = "IN-SST", color.by = "test", values = "unadjusted"),
  cao$plotSampleDistances(space = "expression.shifts", cell.type = "IN-SST", color.by = "Batch", adjust.for = ~ Group),
  ncol = 2)', fig = c(10, 4))

md("## Differential expression and composition on the same model\n\n",
   "The stored model (formula, test, reference level) is used by the per-cell-type differential expression, by the compositional analysis and by the cell density analysis, ",
   "and is recorded with each result. Whole-factor tests are routed to the corresponding multi-group tests (DESeq2 LRT, edgeR/limma F).")
code('de <- cao$estimateDEPerCellType(test = "limma-voom", min.cell.count = 5)
sapply(de, function(d) sum(d$res$padj < 0.05, na.rm = TRUE))             # DE genes per cell type')
md("The top genes of one cell type (`cao$plotVolcano()` draws them when the EnhancedVolcano package is installed):")
code('r <- de[["IN-SST"]]$res; r <- r[order(r$padj), intersect(c("Gene", "log2FoldChange", "stat", "pvalue", "padj", "CellFrac"), names(r))]
head(r, 8)')
code('cao$estimateCellLoadings()\ncao$plotCellLoadings(show.pvals = FALSE)', fig = c(5, 5))

md("## Cluster-free analyses\n\n",
   "Cell density is compared between the two groups using the regression weights of the contrast, so adjusted covariates are accounted for. ",
   "Cluster-free expression shifts test the same contrast in every cell's neighbourhood, sharing the permutations across cells and adjusting the z-scores by the maximum statistic.")
code('cao$estimateCellDensity(method = "graph")
cao$estimateDiffCellDensity(type = "permutation")
plot_grid(cao$plotEmbedding(color.by = "cell.groups"), cao$plotDiffCellDensity(), ncol = 2)', fig = c(10, 4.5))
code('cao$estimateClusterFreeExpressionShifts(n.top.genes = 500, n.permutations = 99)
cao$plotClusterFreeExpressionShifts(font.size = 2)', fig = c(5.5, 4.5))

md("## Session info")
code("sessionInfo()")

## ---- evaluation -------------------------------------------------------------------------------------------
env <- new.env(parent = globalenv())
tmp <- tempfile("fig"); dir.create(tmp)
b64 <- function(f) jsonlite::base64_enc(readBin(f, "raw", file.size(f)))
lines <- function(s) { x <- strsplit(s, "\n", fixed = TRUE)[[1]]; if (!length(x)) x <- ""; as.list(paste0(x, c(rep("\n", length(x) - 1), ""))) }
empty <- setNames(list(), character(0))
nb.cells <- list(); k <- 0L
for (cell in cells) {
  if (cell$type == "markdown") { nb.cells[[length(nb.cells) + 1]] <- list(cell_type = "markdown", metadata = empty, source = lines(cell$text)); next }
  k <- k + 1L
  cat(sprintf("[%s] cell %d: %s\n", format(Sys.time(), "%H:%M:%S"), k, substr(gsub("\n", " | ", cell$text), 1, 90)))
  outs <- list(); txt <- character(); err.txt <- character()
  flush <- function() { if (length(txt)) outs[[length(outs) + 1]] <<- list(output_type = "stream", name = "stdout", text = lines(paste(txt, collapse = ""))); txt <<- character()
                        if (length(err.txt)) outs[[length(outs) + 1]] <<- list(output_type = "stream", name = "stderr", text = lines(paste(err.txt, collapse = ""))); err.txt <<- character() }
  w <- cell$fig[1]; h <- cell$fig[2]
  ev <- evaluate::evaluate(cell$text, envir = env, new_device = TRUE, stop_on_error = 0L, log_echo = FALSE)
  for (it in ev) {
    if (inherits(it, "source")) next
    if (is.character(it)) { txt <- c(txt, it); next }
    if (inherits(it, "message")) { err.txt <- c(err.txt, conditionMessage(it)); next }
    if (inherits(it, "warning")) { err.txt <- c(err.txt, paste0("Warning: ", conditionMessage(it), "\n")); next }
    if (inherits(it, "error")) { flush(); outs[[length(outs) + 1]] <- list(output_type = "error", ename = "Error", evalue = conditionMessage(it), traceback = list(paste0("Error: ", conditionMessage(it))))
      cat("   ERROR:", conditionMessage(it), "\n"); next }
    if (inherits(it, "recordedplot")) {
      flush(); f <- file.path(tmp, sprintf("c%03d_%02d.png", k, length(outs)))
      ragg::agg_png(f, width = w, height = h, units = "in", res = 110); grDevices::replayPlot(it); grDevices::dev.off()
      outs[[length(outs) + 1]] <- list(output_type = "display_data", metadata = list(`image/png` = list(width = round(w * 110), height = round(h * 110))),
                                       data = list(`image/png` = b64(f), `text/plain` = list("plot without title")))
    }
  }
  flush()
  nb.cells[[length(nb.cells) + 1]] <- list(cell_type = "code", execution_count = k, metadata = empty, outputs = outs, source = lines(cell$text))
}
nb <- list(cells = nb.cells,
           metadata = list(kernelspec = list(display_name = "R", language = "R", name = "ir"),
                           language_info = list(name = "R", codemirror_mode = "r", file_extension = ".r", mimetype = "text/x-r-source",
                                                pygments_lexer = "r", version = paste(R.version$major, R.version$minor, sep = "."))),
           nbformat = 4L, nbformat_minor = 4L)
writeLines(jsonlite::toJSON(nb, auto_unbox = TRUE, pretty = TRUE, null = "null"), out.file)
cat("written", out.file, file.size(out.file), "bytes\n")
system2("Rscript", c(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "nb_to_md.R"), out.file, sub("\\.ipynb$", ".md", out.file), sub("\\.ipynb$", "_files", out.file)))
