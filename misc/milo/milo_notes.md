# Milo (miloR) cluster-free differential abundance on the cacoa benchmark datasets

Working directory: `/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad/milo/`

## Software

- R 4.2.2, Bioconductor 3.16, miloR 1.6.0 (installed into `rlib_milo` with `BiocManager::install("miloR")`; no compilation problems, the gnu++14 Makevars was not needed), edgeR 3.40.2, igraph 2.x, Seurat 4.3.0.1, ggplot2 4.0.2 (Milo's `plotDAbeeswarm` was not used; the beeswarm is own ggplot2 + ggbeeswarm code).

## Scripts

| file | purpose |
|---|---|
| `install_milo.R` | installs miloR into the private library |
| `extract_inputs.R <ds>` | loads the large objects once and saves light-weight inputs `inputs_<ds>.rds` (Conos joint graph, 2-D embedding, raw joint count matrix, cell.groups, sample.per.cell, sample.meta) |
| `run_milo.R <ds> [n.cores]` | the Milo pipeline (both routes, all designs, both FDR weightings, per-cell mapping, per-cell-type summaries, beeswarm PNGs, timing); writes `milo_<ds>.rds` |
| `sim_empirical_truth.R` | per-sample composition of the simulation, treated vs control (see "Simulation truth" below) |
| `diag_purity.R <ds>` | condition purity and cell-type purity of neighbourhoods per route |
| `make_report.R` | builds `summary_tables.md` (markdown tables) from the `milo_<ds>.rds` files |

## What was done (identical for the three datasets)

1. **Inputs.** Cells = vertices of the Conos joint graph (`con$graph`), in graph order; `cell.groups`, `sample.per.cell`, 2-D embedding (`con$embedding`) restricted to those cells. MS: 48,918 cells (one cell of the 48,919 in the joint count matrix is not in the graph), SCC: 35,640, simulation: 2,557. Sample metadata: simulation `treatment` (control = reference, 10 vs 10), MS `diagnosis` (Control = reference; 9 Control vs 12 MS samples), SCC built from sample names (`P<i>_Tumor` / `P<i>_Normal`: `condition` with Normal = reference, `patient`; 9 patients x 2 = 18 samples).
2. **Milo object.** `Milo(SingleCellExperiment(counts = joint raw counts))`, `colData` = celltype, sample; `reducedDim(milo, "UMAP") <- con$embedding` (required by `testNhoods` even for graph-overlap weighting, and used by the k-distance weighting).
3. **Route `conos` (primary, same cell graph as cacoa).** The Conos joint graph is used as Milo's kNN graph: edge weights dropped, graph made undirected and simplified, then `graph(milo) <- g`. `buildFromAdjacency()` was not used because for a non-binary square matrix it densifies the matrix (49k x 49k) before binarising; assigning the igraph directly is equivalent to `buildFromAdjacency(adjacency > 0, is.binary = TRUE)`. `milo@.k <- 30` (Milo only uses `.k` for the k-distance FDR weighting: the distance to the 30th nearest neighbour of the index cell in the 2-D embedding).
   Neighbourhoods: `makeNhoods(prop, refined = TRUE, refinement_scheme = "graph")`, seed 1. Decision: the graph refinement scheme (available since miloR 1.3) was preferred over refinement in the 2-D embedding because it uses the same graph as the neighbourhoods and because refined sampling in a 2-D UMAP is not what Milo's reduced-dimension scheme is meant for. A Milo neighbourhood is the full order-1 graph neighbourhood of the index cell, so its size is set by the Conos graph degree, not by k. prop = 0.1 (MS, SCC), 0.2 (simulation, 2.5k cells).
4. **Route `pca` (standard Milo, for reference).** Seurat 4 on the joint raw counts: `NormalizeData`, 2000 variable features, `ScaleData`, `RunPCA(npcs = 30)`; `buildGraph(k = 30, d = 30)`; `makeNhoods(prop, k = 30, d = 30, refined = TRUE, reduced_dims = "PCA")`, seed 1. No sample/batch integration (plain PCA), which is how Milo is often run but which matters for the simulation (see below).
5. `countCells` by sample; `testNhoods(design, design.df = sample.meta, fdr.weighting = "graph-overlap")` (primary) and `"k-distance"` (secondary; uses the index cell's 30-NN distance in the 2-D embedding for the conos route and in the 30 PCs for the pca route); TMM normalisation, QL F-test on the last design coefficient (treated vs control / MS vs Control / Tumor vs Normal), `robust = TRUE`. Designs: `~ condition` for all; SCC additionally `~ patient + condition`. `calcNhoodDistance` was not run (not needed for these two weightings).
6. `annotateNhoods(milo, res, "celltype")`: dominant cell type and its fraction per neighbourhood.
7. **Per-cell mapping** (`cells` data frame in the rds: cell, celltype, sample, milo.logFC, milo.fdr, milo.pval, n.nhoods): for each cell the minimum SpatialFDR over the neighbourhoods containing it and the logFC of that most significant neighbourhood (ties by PValue); cells in no neighbourhood get NA. Done for the graph-overlap SpatialFDR (`cells`) and for the k-distance SpatialFDR (`cells.kdist`).
8. **Summaries** per cell type (`summary`, `summary.kdist`): number of neighbourhoods (by dominant type) and of cells, counts with SpatialFDR < 0.1 and < 0.05 split by logFC sign, fractions, median logFC, cells in no neighbourhood.
9. Beeswarm: one point per neighbourhood, logFC on x, dominant cell type on y, blue = SpatialFDR < 0.1 (graph-overlap), grey otherwise; label = significant / total neighbourhoods of that type. `milo_<ds>_beeswarm.png` is the primary run (conos route, `~ condition`); the other runs are `milo_<ds>_beeswarm_<route>_<design>.png`.

## Output `milo_<ds>.rds` (list)

`dataset`, `params`, `sample.meta` (with the `condition` factor used), `primary` = `"conos|cond"`, `runs[[<route>|<design>]]` = list(route, design, n.nhoods, `nhoods` (per-neighbourhood table: logFC, logCPM, F, PValue, FDR, Nhood, SpatialFDR (= graph-overlap), SpatialFDR.graph_overlap, SpatialFDR.k_distance, celltype, celltype_fraction, nhood.size, index.cell), `cells`, `summary`, `cells.kdist`, `summary.kdist`), `routes[[route]]` = list(`nhood.matrix` (cells x neighbourhoods sparse 0/1), nhood.size, k, prop, graph.degree), `timing` (seconds per step).

