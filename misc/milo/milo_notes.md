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

## Results overview

Primary run for the comparison with cacoa = **conos route, `~ condition`, graph-overlap SpatialFDR** (`milo_<ds>.rds$runs[["conos|cond"]]`, `milo_<ds>_beeswarm.png`).

| dataset | cells | route | nhoods (median size) | cells in >= 1 nhood | SpatialFDR < 0.1 (up / down) | < 0.05 | cells with milo.fdr < 0.1 |
|---|---|---|---|---|---|---|---|
| sim | 2,557 | conos (degree median 341) | 293 (405) | 2557 | 117 (52 / 65) | 82 | 2354 |
| sim | 2,557 | pca k=30 | 140 (88) | 2542 | 81 (42 / 39) | 71 | 2296 |
| MS | 48,918 | conos (degree median 58) | 3339 (150) | 46,545 | 357 (302 / 55) | 0 | 16,123 |
| MS | 48,918 | pca k=30 | 3424 (64) | 48,343 | 1467 (1051 / 416) | 1233 | 31,029 |
| SCC `~condition` | 35,640 | conos (degree median 61) | 2945 (106) | 34,618 | 816 (297 / 519) | 655 | 20,224 |
| SCC `~patient+condition` | 35,640 | conos | 2945 (106) | 34,618 | 609 (228 / 381) | 406 | 18,279 |
| SCC `~condition` | 35,640 | pca k=30 | 2799 (58) | 35,047 | 2390 (1045 / 1345) | 2258 | 33,223 |
| SCC `~patient+condition` | 35,640 | pca k=30 | 2799 (58) | 35,047 | 2396 (1007 / 1389) | 2312 | 33,632 |

### Simulation (sim_density)

**The planted truth does not describe the per-sample composition of the dataset.** `misc/notebook/prepare_density_sim.R` removes 60 % of IN-SST cells in treated samples and 60 % of OPC cells in control samples from `test/sim_treatment.rds`, but the base simulation already has arm-dependent cell numbers (per-sample counts, t-test on log counts, 10 vs 10 samples):

| cell type | base sim_treatment: control / treated mean cells | sim_density: control / treated mean cells | log2FC | p (log counts) | stated truth | empirical |
|---|---|---|---|---|---|---|
| AST-FB | 11.7 / 18.8 | 11.7 / 18.8 | +0.68 | 2e-4 | none | up |
| AST-PP | 14.4 / 15.9 | 14.4 / 15.9 | +0.14 | 0.44 | none | none |
| IN-PV | 10.0 / 21.3 | 10.0 / 21.3 | +1.09 | 3e-6 | none | up |
| IN-SST | 10.0 / 21.3 | 10.0 / 8.5 | -0.23 | 0.19 | down | none (n.s.) |
| L2/3 | 19.6 / 10.9 | 19.6 / 10.9 | -0.85 | 2e-3 | none | down |
| L5/6-CC | 19.6 / 10.9 | 19.6 / 10.9 | -0.85 | 2e-3 | none | down |
| Microglia | 14.4 / 15.9 | 14.4 / 15.9 | +0.14 | 0.42 | none | none |
| Neu-NRGN | 14.4 / 15.9 | 14.4 / 15.9 | +0.14 | 0.45 | none | none |
| OPC | 11.7 / 18.8 | 4.7 / 18.8 | +2.00 | 5e-9 | up | up |

So in terms of sample-level abundance, OPC is a 4-fold increase (planted on top of an existing 1.6x), the planted IN-SST decrease roughly cancels a pre-existing 2x increase (net -0.23 log2, p = 0.19 on counts, p = 0.02 on logit proportions, the latter because treated samples have more cells: 137 vs 119 on average), and IN-PV (2x up), L2/3 and L5/6-CC (2x down) and AST-FB (1.6x up) differ between arms although the truth table says "none". Any sample-level abundance test (Milo, or a cluster-level cacoa test) must call these; whether cacoa's cluster-free density test "should" is a question for the lead (if its density is relative within a sample, the same composition differences apply). The `truth` table should be regenerated from the actual per-sample counts, or the base simulation should be re-drawn with equal expected cell-type proportions across arms. Script: `sim_empirical_truth.R` (`sim_empirical_truth.rds`).

Milo, conos route (`~ treatment`): IN-SST 18/18 neighbourhoods down (median logFC -0.88; 17 at FDR < 0.05), OPC 17/17 up (logFC +1.6). Against the stated truth there are 82/258 "false positive" neighbourhoods in the other 7 types (AST-FB 19/42 up, IN-PV 7/30 up, L2/3 22/36 down, L5/6-CC 23/46 down, Microglia 7/42 up, Neu-NRGN 2 up/2 down, AST-PP 0/26); all except Microglia/Neu-NRGN go in the direction of the empirical composition difference. The IN-SST logFC of -0.9 (larger than the -0.23 per-sample count change) comes from the library-size (TMM) normalisation across samples of different size and from impure neighbourhoods: the Conos graph in this simulation is very dense (median degree 341 in 2557 cells; a Milo neighbourhood = full graph neighbourhood, median 405 cells), so neighbourhoods are only 72 % pure on median (min 40 %) and IN-SST neighbourhoods contain L2/3 and L5/6-CC cells. Per-cell mapping is correspondingly coarse: every cell is in >= 1 neighbourhood and 92 % of cells have milo.fdr < 0.1 (in AST-PP 33 %).

Milo, pca route (plain PCA, no integration, k = 30): the planted expression shifts (truth `expression` = strong for IN-PV, IN-SST, L2/3; moderate for L5/6-CC) split those cell types by condition in the kNN graph: 100 % of their neighbourhoods are condition-pure (all treated or all control cells; `diag_purity.R`), giving |logFC| ~ 6 in both directions and "significant" calls for every neighbourhood of those types in both directions. This is the well-known Milo failure mode without batch/condition integration and makes the pca route uninformative for this benchmark (IN-SST 6/8 down but also 2/8 up; L2/3 8 up / 9 down). OPC 13/13 up is correct; AST-FB 4/22 up.

### MS (Schirmer 2019; 9 Control vs 12 MS, 19 cell types)

Conos route, `~ diagnosis`: 357/3339 neighbourhoods at SpatialFDR < 0.1, none at < 0.05. The signal is marginal overall: min PValue 2.8e-5 over 3339 tests, plain BH FDR >= 0.094 everywhere, SpatialFDR minimum 0.066 (the 357 calls all lie between 0.066 and 0.1). The calls are: Microglia 76/109 up (median logFC +1.5), Phagocytes 21/22 up, OL 106/811 up, Astrocytes 51 up / 8 down, OPC 19 up / 4 down, Glia-MIX 8 up, EN-L2-3 27 down (a tail of neighbourhoods with logFC < -1.5, 5 % of EN-L2-3 neighbourhoods), EN-L4 12 down, EN-L5-6 8 up, Stromal 4 up, B cells 3/4 up; inhibitory neurons essentially nothing (IN-SST, IN-SV2C, IN-VIP 0; IN-PV 2). Per cell: 16,123/48,918 cells have milo.fdr < 0.1 (Microglia 94 %, Phagocytes 97 %, OL 42 %, Astrocytes 51 %, EN-L2-3 25 % down, EN-L4 38 % down). The k-distance weighting (30-NN distance in the 2-D embedding) calls 1 neighbourhood; see Problems.

pca route (no integration): 1467/3424 at FDR < 0.1, 1233 at < 0.05; 48-60 % of OL, Astrocyte, EN-L2-3, EN-PYR, Microglia neighbourhoods are condition-extreme (> 90 % or < 10 % MS cells), EN-L2-3 has 100 up and 228 down. Not a fair comparison with cacoa, kept for reference only.

### SCC (Ji 2020; 9 patients, Tumor vs Normal)

Conos route, `~ condition`: 816/2945 neighbourhoods at SpatialFDR < 0.1 (297 up in tumour, 519 down), 655 at < 0.05; `~ patient + condition`: 609 (228 up, 381 down), 406 at < 0.05. Both designs agree qualitatively: down in tumour = KC_Diff (337/881, resp. 247), KC_Basal (124/385, resp. 89), KC_Cyc (21-23/239), Pilosebaceous (9/167); up in tumour = Mac (81-91/228), PDC (17/17, resp. 12), TSK (18/26, resp. 11), ASDC (20/31, resp. 16), B cells (12/13, resp. 8), CD1C (49-50/360), Fibroblast (28/83 with `~condition`, 2/83 with the patient term, i.e. the fibroblast increase is between-patient), Melanocyte (13 up with `~condition`, 5 down with the patient term). Per cell: 20,224 (resp. 18,279)/35,640 cells have milo.fdr < 0.1 (KC_Diff 67 %, KC_Basal 64 %, Mac 66 %, PDC 95 %, TSK 88 %). Conos-route neighbourhoods are 91 % pure (median); 24-28 % of KC_Basal/KC_Diff neighbourhoods are condition-extreme, which here reflects genuinely tumour-specific keratinocyte states (TSK/tumour KC) rather than a failed alignment.

pca route (no integration): 85 % of neighbourhoods significant with either design (keratinocyte neighbourhoods 95-98 % condition-pure, |logFC| 5-6 in both directions within KC_Basal/KC_Diff); kept for reference only.

## Problems and decisions

1. **Simulation truth vs composition** (above): the main benchmark cannot be scored against the stated truth as is; against the empirical composition Milo (conos route) has 100 % sensitivity for OPC and IN-SST (the latter helped by neighbourhood impurity), 45-61 % of neighbourhoods called in the empirically changed types AST-FB, L2/3, L5/6-CC, 23 % in IN-PV, and 0-17 % in the three empirically unchanged types (AST-PP 0/26, Microglia 7/42, Neu-NRGN 4/36 at FDR < 0.1; 0/26, 2/42, 1/36 at < 0.05).
2. **Dense Conos graph in the simulation**: median degree 341 (20 samples of ~128 cells, all pairs aligned) gives Milo neighbourhoods of ~400 cells (a whole cell type) with 72 % purity; the sim's "cluster-free" comparison is therefore nearly a cluster-level one on the Milo side. MS/SCC graphs are normal (median degree ~60, neighbourhoods 100-150 cells, purity 0.91-0.95).
3. **k-distance SpatialFDR on the 2-D embedding is unreliable** (MS conos: 1 significant neighbourhood vs 357 with graph-overlap, because the 1/distance weights in a 2-D UMAP span orders of magnitude and the weighted BH is dominated by neighbourhoods in dense embedding regions). It is also slow in miloR 1.6.0 (one `findKNN` call over all cells per neighbourhood: 15-25 min per test for the conos route, 70-130 min for the pca route in 30 PCs). graph-overlap is the weighting used for all reported numbers; k-distance results are kept in the rds (`SpatialFDR.k_distance`, `cells.kdist`, `summary.kdist`) and agree with graph-overlap for sim/SCC.
4. **pca route without integration** is dominated by condition-pure neighbourhoods in all three datasets (expression differences between conditions, not abundance). It is included because it is the "standard" Milo workflow, but the fair head-to-head with cacoa is the conos route (same graph).
5. `annotateNhoods` dominant type with fraction: for the summaries a neighbourhood counts for its dominant type only; the per-cell table (`cells`) does not depend on this.
6. `SpatialFDR < 0.05` in MS conos is empty (see above); the lead may prefer to compare at matched thresholds or on the p-value/FDR rank scale (`PValue`, `FDR`, `SpatialFDR` are all in the per-neighbourhood table).
7. Minor: Seurat renamed gene names with underscores (warning only); one MS cell present in the joint count matrix but not in the graph was dropped; cells in no neighbourhood (MS 2373, SCC 1022, sim 0 for conos) have NA in the per-cell table.
8. ggplot2 4.0.2 caused no problems with own beeswarm code (`ggbeeswarm::geom_quasirandom`); Milo's `plotDAbeeswarm` was not tried.

## Runtime (seconds; single core, Milo does not parallelise; the n.cores argument is unused)

| step | sim | MS | SCC |
|---|---|---|---|
| load inputs (from `inputs_<ds>.rds`; the initial extraction from the raw objects took ~1-2 min per dataset) | 4 | 21 | 21 |
| conos: set graph | 0.4 | 3 | 2 |
| conos: makeNhoods (graph-refined) | 42 | 23 | 13 |
| countCells | 0.2 | 0.2 | 0.2 |
| testNhoods graph-overlap (per design) | 0.4 | 1.3 | 1.1 / 4.8 |
| testNhoods k-distance (per design) | 5 | 1471 | 886 / 876 |
| pca: Seurat PCA | 8 | 53 | 45 |
| pca: buildGraph k=30 | 0.3 | 9 | 6 |
| pca: makeNhoods | 0.4 | 19 | 12 |
| pca: testNhoods graph-overlap / k-distance | 0.2 / 4 | 1.2 / 7845 | 0.9 / 4024 (and 5.6 / 4409) |

Total wall time: sim ~1.5 min; MS ~2.7 h; SCC ~2.9 h (almost all of it the k-distance weighting; without it each dataset runs in ~3 min).

## Per-cell-type summary tables (all runs; generated by make_report.R)

Columns: n.nhoods = neighbourhoods whose dominant type is the cell type; nh.fdr10.up/down = of these, SpatialFDR (graph-overlap) < 0.1 with logFC > 0 / < 0; nh.fdr05.* = same at < 0.05; frac.nh.fdr10 = fraction of neighbourhoods at < 0.1; n.cells = cells of the type; n.cells.no.nhood = cells in no neighbourhood; cells.fdr10.up/down = cells whose most significant neighbourhood has SpatialFDR < 0.1 and logFC > 0 / < 0; frac.cells.fdr10 = fraction of cells at < 0.1.

### Dataset: sim

Parameters: prop = 0.2; k (conos route, .k for k-distance weighting) = 30; k (pca route) = 30, d = 30; condition column = treatment (reference = control); miloR 1.6.0

- route conos: 293 neighbourhoods; size median 405 (IQR 345-442, range 185-567); cells in >= 1 neighbourhood: 2557 / 2557; graph degree median 341
- route pca: 140 neighbourhoods; size median 88 (IQR 73-109.25, range 48-182); cells in >= 1 neighbourhood: 2542 / 2557; graph degree median 39

#### Run `conos|cond`: route conos, design `~condition`, SpatialFDR weighting graph-overlap (primary)

Neighbourhoods: 293; SpatialFDR < 0.1: 117 (up 52, down 65); SpatialFDR < 0.05: 82; with k-distance weighting: FDR < 0.1: 112, < 0.05: 79; cells with milo.fdr < 0.1: 2354 / 2557; dominant-type fraction median 0.72

| celltype | n.nhoods | nh.fdr10.up | nh.fdr10.down | nh.fdr05.up | nh.fdr05.down | frac.nh.fdr10 | median.logFC | n.cells | n.cells.no.nhood | cells.fdr10.up | cells.fdr10.down | cells.fdr05.up | cells.fdr05.down | frac.cells.fdr10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AST-FB | 42 | 19 | 0 | 15 | 0 | 0.452 | 0.433 | 305 | 0 | 295 | 10 | 295 | 10 | 1.000000 |
| AST-PP | 26 | 0 | 0 | 0 | 0 | 0.0000000 | 0.109 | 303 | 0 | 38 | 62 | 33 | 47 | 0.33 |
| IN-PV | 30 | 7 | 0 | 4 | 0 | 0.233 | 0.284 | 313 | 0 | 296 | 17 | 296 | 17 | 1.000000 |
| IN-SST | 18 | 0 | 18 | 0 | 17 | 1.0000000 | -0.883 | 185 | 0 | 2 | 183 | 2 | 183 | 1.000000 |
| L2/3 | 36 | 0 | 22 | 0 | 14 | 0.611 | -0.491 | 305 | 0 | 107 | 198 | 107 | 198 | 1.000000 |
| L5/6-CC | 46 | 0 | 23 | 0 | 12 | 0.500 | -0.462 | 305 | 0 | 88 | 217 | 88 | 217 | 1.000000 |
| Microglia | 42 | 7 | 0 | 2 | 0 | 0.167 | 0.162 | 303 | 0 | 205 | 98 | 199 | 98 | 1.000000 |
| Neu-NRGN | 36 | 2 | 2 | 1 | 0 | 0.111 | 0.166 | 303 | 0 | 229 | 74 | 229 | 74 | 1.000000 |
| OPC | 17 | 17 | 0 | 17 | 0 | 1.0000000 | 1.590 | 235 | 0 | 235 | 0 | 235 | 0 | 1.000000 |

k-distance weighting, neighbourhoods with SpatialFDR < 0.1 (up/down) per cell type: AST-FB 17/0; AST-PP 0/0; IN-PV 7/0; IN-SST 0/17; L2/3 0/22; L5/6-CC 0/21; Microglia 7/0; Neu-NRGN 2/2; OPC 17/0

Truth comparison (stated truth: IN-SST down, OPC up, others none):

- sensitivity IN-SST: 18 / 18 neighbourhoods (100%) significant in the planted direction at FDR < 0.1; 17 at FDR < 0.05
- sensitivity OPC: 17 / 17 neighbourhoods (100%) significant in the planted direction at FDR < 0.1; 17 at FDR < 0.05
- stated-truth false positives in the other 7 types: 82 / 258 neighbourhoods (32%) at FDR < 0.1; 48 at FDR < 0.05; by type: AST-FB 19/42, AST-PP 0/26, IN-PV 7/30, L2/3 22/36, L5/6-CC 23/46, Microglia 7/42, Neu-NRGN 4/36
- empirical composition truth (per-sample count t-test, see notes): AST-FB up, AST-PP none, IN-PV up, IN-SST none, L2/3 down, L5/6-CC down, Microglia none, Neu-NRGN none, OPC up
- versus empirical truth: AST-FB 19/42 sig in empirical direction; AST-PP 0/26 sig; IN-PV 7/30 sig in empirical direction; IN-SST 18/18 sig; L2/3 22/36 sig in empirical direction; L5/6-CC 23/46 sig in empirical direction; Microglia 7/42 sig; Neu-NRGN 4/36 sig; OPC 17/17 sig in empirical direction

#### Run `pca|cond`: route pca, design `~condition`, SpatialFDR weighting graph-overlap (primary)

Neighbourhoods: 140; SpatialFDR < 0.1: 81 (up 42, down 39); SpatialFDR < 0.05: 71; with k-distance weighting: FDR < 0.1: 77, < 0.05: 69; cells with milo.fdr < 0.1: 2296 / 2557; dominant-type fraction median 1

| celltype | n.nhoods | nh.fdr10.up | nh.fdr10.down | nh.fdr05.up | nh.fdr05.down | frac.nh.fdr10 | median.logFC | n.cells | n.cells.no.nhood | cells.fdr10.up | cells.fdr10.down | cells.fdr05.up | cells.fdr05.down | frac.cells.fdr10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AST-FB | 22 | 4 | 0 | 2 | 0 | 0.182 | 0.144 | 305 | 0 | 212 | 0 | 165 | 0 | 0.695 |
| AST-PP | 12 | 0 | 3 | 0 | 2 | 0.250 | -0.256 | 303 | 4 | 0 | 241 | 0 | 121 | 0.795 |
| IN-PV | 14 | 11 | 3 | 11 | 3 | 1.0000000 | 6.176 | 313 | 4 | 210 | 99 | 210 | 99 | 0.987 |
| IN-SST | 8 | 2 | 6 | 2 | 6 | 1.0000000 | -6.249 | 185 | 1 | 84 | 100 | 84 | 100 | 0.995 |
| L2/3 | 17 | 8 | 9 | 8 | 9 | 1.0000000 | -5.895 | 305 | 1 | 109 | 195 | 109 | 195 | 0.997 |
| L5/6-CC | 13 | 4 | 9 | 4 | 9 | 1.0000000 | -6.304 | 305 | 4 | 108 | 193 | 108 | 193 | 0.987 |
| Microglia | 23 | 0 | 3 | 0 | 0 | 0.130 | -0.308 | 303 | 0 | 0 | 239 | 0 | 0 | 0.789 |
| Neu-NRGN | 18 | 0 | 6 | 0 | 2 | 0.333 | -0.437 | 303 | 1 | 0 | 271 | 0 | 216 | 0.894 |
| OPC | 13 | 13 | 0 | 13 | 0 | 1.0000000 | 1.383 | 235 | 0 | 235 | 0 | 235 | 0 | 1.0000000 |

k-distance weighting, neighbourhoods with SpatialFDR < 0.1 (up/down) per cell type: AST-FB 2/0; AST-PP 0/3; IN-PV 11/3; IN-SST 2/6; L2/3 8/9; L5/6-CC 4/9; Microglia 0/1; Neu-NRGN 0/6; OPC 13/0

Truth comparison (stated truth: IN-SST down, OPC up, others none):

- sensitivity IN-SST: 6 / 8 neighbourhoods (75%) significant in the planted direction at FDR < 0.1; 6 at FDR < 0.05
- sensitivity OPC: 13 / 13 neighbourhoods (100%) significant in the planted direction at FDR < 0.1; 13 at FDR < 0.05
- stated-truth false positives in the other 7 types: 60 / 119 neighbourhoods (50%) at FDR < 0.1; 50 at FDR < 0.05; by type: AST-FB 4/22, AST-PP 3/12, IN-PV 14/14, L2/3 17/17, L5/6-CC 13/13, Microglia 3/23, Neu-NRGN 6/18
- empirical composition truth (per-sample count t-test, see notes): AST-FB up, AST-PP none, IN-PV up, IN-SST none, L2/3 down, L5/6-CC down, Microglia none, Neu-NRGN none, OPC up
- versus empirical truth: AST-FB 4/22 sig in empirical direction; AST-PP 3/12 sig; IN-PV 11/14 sig in empirical direction; IN-SST 8/8 sig; L2/3 9/17 sig in empirical direction; L5/6-CC 9/13 sig in empirical direction; Microglia 3/23 sig; Neu-NRGN 6/18 sig; OPC 13/13 sig in empirical direction

Runtime (s): load_inputs 4.3; conos_build_graph 0.4; conos_makeNhoods 41.5; conos_countCells 0.2; conos_cond_testNhoods_graph-overlap 0.4; conos_cond_testNhoods_k-distance 5.2; pca_seurat_pca 7.5; pca_buildGraph 0.3; pca_makeNhoods 0.4; pca_countCells 0.1; pca_cond_testNhoods_graph-overlap 0.2; pca_cond_testNhoods_k-distance 3.7

### Dataset: ms

Parameters: prop = 0.1; k (conos route, .k for k-distance weighting) = 30; k (pca route) = 30, d = 30; condition column = diagnosis (reference = Control); miloR 1.6.0

- route conos: 3339 neighbourhoods; size median 150 (IQR 88-212, range 8-438); cells in >= 1 neighbourhood: 46545 / 48918; graph degree median 58
- route pca: 3424 neighbourhoods; size median 64 (IQR 50-86.25, range 31-396); cells in >= 1 neighbourhood: 48343 / 48918; graph degree median 39

#### Run `conos|cond`: route conos, design `~condition`, SpatialFDR weighting graph-overlap (primary)

Neighbourhoods: 3339; SpatialFDR < 0.1: 357 (up 302, down 55); SpatialFDR < 0.05: 0; with k-distance weighting: FDR < 0.1: 1, < 0.05: 0; cells with milo.fdr < 0.1: 16123 / 48918; dominant-type fraction median 0.95

| celltype | n.nhoods | nh.fdr10.up | nh.fdr10.down | nh.fdr05.up | nh.fdr05.down | frac.nh.fdr10 | median.logFC | n.cells | n.cells.no.nhood | cells.fdr10.up | cells.fdr10.down | cells.fdr05.up | cells.fdr05.down | frac.cells.fdr10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Astrocytes | 306 | 51 | 8 | 0 | 0 | 0.193 | 0.712 | 5381 | 246 | 1903 | 841 | 0 | 0 | 0.510 |
| B cells | 4 | 3 | 0 | 0 | 0 | 0.750 | 2.169 | 53 | 4 | 41 | 0 | 0 | 0 | 0.774 |
| EN-L2-3 | 518 | 0 | 27 | 0 | 0 | 0.052 | -0.359 | 6120 | 210 | 7 | 1555 | 0 | 0 | 0.255 |
| EN-L4 | 285 | 1 | 12 | 0 | 0 | 0.046 | -0.312 | 3125 | 54 | 51 | 1126 | 0 | 0 | 0.377 |
| EN-L5-6 | 284 | 8 | 1 | 0 | 0 | 0.032 | 0.027 | 3058 | 84 | 100 | 295 | 0 | 0 | 0.129 |
| EN-MIX | 18 | 0 | 0 | 0 | 0 | 0.00000000 | 0.016 | 2402 | 243 | 63 | 323 | 0 | 0 | 0.161 |
| EN-PYR | 119 | 5 | 0 | 0 | 0 | 0.042 | 0.267 | 2082 | 190 | 163 | 3 | 0 | 0 | 0.080 |
| Endo cells | 25 | 0 | 0 | 0 | 0 | 0.00000000 | 0.579 | 375 | 12 | 115 | 0 | 0 | 0 | 0.307 |
| Glia-MIX | 78 | 8 | 0 | 0 | 0 | 0.103 | 0.967 | 1226 | 128 | 316 | 43 | 0 | 0 | 0.293 |
| IN-PV | 139 | 0 | 2 | 0 | 0 | 0.014 | -0.103 | 1897 | 46 | 2 | 449 | 0 | 0 | 0.238 |
| IN-SST | 100 | 0 | 0 | 0 | 0 | 0.00000000 | -0.015 | 1486 | 40 | 9 | 10 | 0 | 0 | 0.013 |
| IN-SV2C | 62 | 0 | 0 | 0 | 0 | 0.00000000 | -0.177 | 779 | 12 | 2 | 0 | 0 | 0 | 0.003 |
| IN-VIP | 226 | 0 | 0 | 0 | 0 | 0.00000000 | -0.070 | 2755 | 41 | 14 | 1 | 0 | 0 | 0.005 |
| Microglia | 109 | 76 | 0 | 0 | 0 | 0.697 | 1.493 | 1425 | 38 | 1333 | 0 | 0 | 0 | 0.935 |
| OL | 811 | 106 | 1 | 0 | 0 | 0.132 | 0.826 | 12394 | 839 | 5141 | 51 | 0 | 0 | 0.419 |
| OPC | 219 | 19 | 4 | 0 | 0 | 0.105 | 0.401 | 3844 | 167 | 1021 | 707 | 0 | 0 | 0.450 |
| Phagocytes | 22 | 21 | 0 | 0 | 0 | 0.955 | 2.326 | 258 | 7 | 251 | 0 | 0 | 0 | 0.973 |
| Stromal cells | 12 | 4 | 0 | 0 | 0 | 0.333 | 1.228 | 153 | 0 | 141 | 0 | 0 | 0 | 0.922 |
| T cells | 2 | 0 | 0 | 0 | 0 | 0.00000000 | 2.255 | 105 | 12 | 46 | 0 | 0 | 0 | 0.438 |

k-distance weighting, neighbourhoods with SpatialFDR < 0.1 (up/down) per cell type: Astrocytes 0/0; B cells 0/0; EN-L2-3 0/0; EN-L4 0/0; EN-L5-6 0/0; EN-MIX 0/0; EN-PYR 0/0; Endo cells 0/0; Glia-MIX 0/0; IN-PV 0/0; IN-SST 0/0; IN-SV2C 0/0; IN-VIP 0/0; Microglia 0/0; OL 1/0; OPC 0/0; Phagocytes 0/0; Stromal cells 0/0; T cells 0/0

#### Run `pca|cond`: route pca, design `~condition`, SpatialFDR weighting graph-overlap (primary)

Neighbourhoods: 3424; SpatialFDR < 0.1: 1467 (up 1051, down 416); SpatialFDR < 0.05: 1233; with k-distance weighting: FDR < 0.1: 1468, < 0.05: 1227; cells with milo.fdr < 0.1: 31029 / 48918; dominant-type fraction median 1

| celltype | n.nhoods | nh.fdr10.up | nh.fdr10.down | nh.fdr05.up | nh.fdr05.down | frac.nh.fdr10 | median.logFC | n.cells | n.cells.no.nhood | cells.fdr10.up | cells.fdr10.down | cells.fdr05.up | cells.fdr05.down | frac.cells.fdr10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Astrocytes | 367 | 148 | 39 | 129 | 33 | 0.510 | 1.851 | 5381 | 79 | 3076 | 1038 | 2816 | 897 | 0.765 |
| B cells | 2 | 2 | 0 | 2 | 0 | 1.00000000 | 5.236 | 53 | 4 | 49 | 0 | 49 | 0 | 0.925 |
| EN-L2-3 | 515 | 100 | 228 | 85 | 192 | 0.637 | -1.918 | 6120 | 34 | 1591 | 3598 | 1436 | 3312 | 0.848 |
| EN-L4 | 221 | 33 | 41 | 28 | 27 | 0.335 | -0.766 | 3125 | 49 | 737 | 1204 | 646 | 923 | 0.621 |
| EN-L5-6 | 235 | 72 | 19 | 68 | 10 | 0.387 | 0.681 | 3058 | 22 | 1174 | 746 | 1107 | 455 | 0.628 |
| EN-MIX | 150 | 15 | 22 | 13 | 19 | 0.247 | -0.256 | 2402 | 33 | 440 | 635 | 391 | 531 | 0.448 |
| EN-PYR | 145 | 61 | 18 | 52 | 9 | 0.545 | 1.918 | 2082 | 30 | 1177 | 496 | 1073 | 302 | 0.804 |
| Endo cells | 22 | 3 | 0 | 2 | 0 | 0.136 | 0.587 | 375 | 10 | 101 | 0 | 73 | 0 | 0.269 |
| Glia-MIX | 72 | 18 | 0 | 15 | 0 | 0.250 | 0.805 | 1226 | 36 | 621 | 24 | 512 | 17 | 0.526 |
| IN-PV | 124 | 9 | 5 | 8 | 2 | 0.113 | -0.216 | 1897 | 21 | 208 | 220 | 184 | 71 | 0.226 |
| IN-SST | 124 | 9 | 2 | 8 | 1 | 0.089 | 0.077 | 1486 | 9 | 210 | 85 | 192 | 46 | 0.199 |
| IN-SV2C | 56 | 0 | 1 | 0 | 0 | 0.018 | -0.026 | 779 | 11 | 0 | 23 | 0 | 10 | 0.030 |
| IN-VIP | 243 | 9 | 11 | 8 | 10 | 0.082 | -0.283 | 2755 | 19 | 258 | 406 | 235 | 388 | 0.241 |
| Microglia | 117 | 78 | 0 | 65 | 0 | 0.667 | 3.155 | 1425 | 15 | 1212 | 0 | 1123 | 0 | 0.851 |
| OL | 767 | 390 | 25 | 351 | 9 | 0.541 | 2.902 | 12394 | 132 | 7645 | 1422 | 7234 | 625 | 0.732 |
| OPC | 223 | 77 | 5 | 64 | 4 | 0.368 | 1.339 | 3844 | 53 | 1922 | 237 | 1731 | 209 | 0.562 |
| Phagocytes | 24 | 16 | 0 | 12 | 0 | 0.667 | 2.835 | 258 | 3 | 253 | 0 | 243 | 0 | 0.981 |
| Stromal cells | 11 | 5 | 0 | 1 | 0 | 0.455 | 2.450 | 153 | 6 | 130 | 0 | 38 | 0 | 0.850 |
| T cells | 6 | 6 | 0 | 6 | 0 | 1.00000000 | 4.456 | 105 | 9 | 91 | 0 | 78 | 0 | 0.867 |

k-distance weighting, neighbourhoods with SpatialFDR < 0.1 (up/down) per cell type: Astrocytes 148/39; B cells 2/0; EN-L2-3 100/228; EN-L4 33/41; EN-L5-6 72/19; EN-MIX 15/22; EN-PYR 61/18; Endo cells 3/0; Glia-MIX 18/0; IN-PV 9/5; IN-SST 9/2; IN-SV2C 0/1; IN-VIP 9/11; Microglia 78/0; OL 391/25; OPC 77/5; Phagocytes 16/0; Stromal cells 5/0; T cells 6/0

Runtime (s): load_inputs 21.1; conos_build_graph 3.1; conos_makeNhoods 22.8; conos_countCells 0.2; conos_cond_testNhoods_graph-overlap 1.3; conos_cond_testNhoods_k-distance 1470.7; pca_seurat_pca 52.8; pca_buildGraph 9.4; pca_makeNhoods 18.5; pca_countCells 0.2; pca_cond_testNhoods_graph-overlap 1.2; pca_cond_testNhoods_k-distance 7845.3

### Dataset: scc

Parameters: prop = 0.1; k (conos route, .k for k-distance weighting) = 30; k (pca route) = 30, d = 30; condition column = condition (reference = Normal); miloR 1.6.0

- route conos: 2945 neighbourhoods; size median 106 (IQR 67-154, range 11-273); cells in >= 1 neighbourhood: 34618 / 35640; graph degree median 61
- route pca: 2799 neighbourhoods; size median 58 (IQR 48-72, range 31-176); cells in >= 1 neighbourhood: 35047 / 35640; graph degree median 39

#### Run `conos|cond`: route conos, design `~condition`, SpatialFDR weighting graph-overlap (primary)

Neighbourhoods: 2945; SpatialFDR < 0.1: 816 (up 297, down 519); SpatialFDR < 0.05: 655; with k-distance weighting: FDR < 0.1: 606, < 0.05: 463; cells with milo.fdr < 0.1: 20224 / 35640; dominant-type fraction median 0.91

| celltype | n.nhoods | nh.fdr10.up | nh.fdr10.down | nh.fdr05.up | nh.fdr05.down | frac.nh.fdr10 | median.logFC | n.cells | n.cells.no.nhood | cells.fdr10.up | cells.fdr10.down | cells.fdr05.up | cells.fdr05.down | frac.cells.fdr10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ASDC | 31 | 20 | 0 | 16 | 0 | 0.645 | 2.571 | 283 | 7 | 255 | 3 | 240 | 3 | 0.912 |
| B Cell | 13 | 12 | 0 | 12 | 0 | 0.923 | 4.166 | 205 | 13 | 181 | 3 | 180 | 3 | 0.898 |
| CD1C | 360 | 49 | 5 | 36 | 4 | 0.150 | 0.349 | 3985 | 120 | 1765 | 203 | 1395 | 146 | 0.494 |
| CLEC9A | 38 | 6 | 0 | 3 | 0 | 0.158 | 0.638 | 556 | 22 | 336 | 1 | 260 | 1 | 0.606 |
| Endothelial Cell | 33 | 2 | 1 | 1 | 1 | 0.091 | -0.335 | 421 | 5 | 144 | 74 | 57 | 74 | 0.518 |
| Eccrine | 4 | 0 | 0 | 0 | 0 | 0.00000000 | -1.239 | 59 | 0 | 3 | 1 | 3 | 0 | 0.068 |
| Fibroblast | 83 | 28 | 0 | 25 | 0 | 0.337 | 0.253 | 840 | 28 | 490 | 29 | 409 | 29 | 0.618 |
| KC_Basal | 385 | 5 | 124 | 2 | 109 | 0.335 | -0.893 | 4808 | 156 | 272 | 2810 | 158 | 2619 | 0.641 |
| KC_Cyc | 239 | 3 | 21 | 3 | 15 | 0.100 | -0.221 | 2459 | 39 | 175 | 778 | 111 | 541 | 0.388 |
| KC_Diff | 881 | 20 | 337 | 11 | 280 | 0.405 | -1.307 | 10496 | 328 | 700 | 6303 | 455 | 5595 | 0.667 |
| LC | 221 | 12 | 9 | 9 | 7 | 0.095 | 0.032 | 3072 | 113 | 551 | 464 | 386 | 377 | 0.330 |
| Mac | 228 | 81 | 1 | 62 | 0 | 0.360 | 1.075 | 3020 | 72 | 1889 | 115 | 1735 | 44 | 0.664 |
| MDSC | 62 | 2 | 4 | 2 | 3 | 0.097 | -0.171 | 711 | 23 | 101 | 145 | 84 | 129 | 0.346 |
| Melanocyte | 55 | 13 | 3 | 11 | 2 | 0.291 | -0.301 | 732 | 22 | 332 | 137 | 288 | 57 | 0.641 |
| NK | 3 | 0 | 0 | 0 | 0 | 0.00000000 | 0.824 | 111 | 4 | 5 | 2 | 5 | 2 | 0.063 |
| PDC | 17 | 17 | 0 | 16 | 0 | 1.00000000 | 3.416 | 245 | 10 | 232 | 0 | 232 | 0 | 0.947 |
| Pilosebaceous | 167 | 4 | 9 | 0 | 5 | 0.078 | -0.076 | 1988 | 23 | 219 | 557 | 53 | 352 | 0.390 |
| Tcell | 99 | 5 | 5 | 4 | 4 | 0.101 | 0.081 | 1327 | 32 | 366 | 300 | 270 | 277 | 0.502 |
| TSK | 26 | 18 | 0 | 12 | 0 | 0.692 | 2.261 | 322 | 5 | 274 | 9 | 241 | 3 | 0.879 |

k-distance weighting, neighbourhoods with SpatialFDR < 0.1 (up/down) per cell type: ASDC 14/0; B Cell 12/0; CD1C 31/4; CLEC9A 2/0; Endothelial Cell 1/1; Eccrine 0/0; Fibroblast 24/0; KC_Basal 2/101; KC_Cyc 3/11; KC_Diff 9/263; LC 8/7; Mac 56/0; MDSC 2/3; Melanocyte 11/2; NK 0/0; PDC 16/0; Pilosebaceous 0/4; Tcell 4/4; TSK 11/0

#### Run `conos|patient_cond`: route conos, design `~patient + condition`, SpatialFDR weighting graph-overlap (primary)

Neighbourhoods: 2945; SpatialFDR < 0.1: 609 (up 228, down 381); SpatialFDR < 0.05: 406; with k-distance weighting: FDR < 0.1: 634, < 0.05: 428; cells with milo.fdr < 0.1: 18279 / 35640; dominant-type fraction median 0.91

| celltype | n.nhoods | nh.fdr10.up | nh.fdr10.down | nh.fdr05.up | nh.fdr05.down | frac.nh.fdr10 | median.logFC | n.cells | n.cells.no.nhood | cells.fdr10.up | cells.fdr10.down | cells.fdr05.up | cells.fdr05.down | frac.cells.fdr10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ASDC | 31 | 16 | 0 | 8 | 0 | 0.516 | 3.072 | 283 | 7 | 236 | 0 | 186 | 0 | 0.834 |
| B Cell | 13 | 8 | 0 | 7 | 0 | 0.615 | 4.014 | 205 | 13 | 167 | 3 | 166 | 3 | 0.829 |
| CD1C | 360 | 50 | 2 | 28 | 1 | 0.144 | 0.538 | 3985 | 120 | 2143 | 116 | 1516 | 82 | 0.567 |
| CLEC9A | 38 | 3 | 0 | 1 | 0 | 0.079 | 0.956 | 556 | 22 | 289 | 0 | 199 | 0 | 0.520 |
| Endothelial Cell | 33 | 0 | 2 | 0 | 1 | 0.061 | -0.607 | 421 | 5 | 46 | 169 | 39 | 152 | 0.511 |
| Eccrine | 4 | 0 | 0 | 0 | 0 | 0.00000000 | -1.115 | 59 | 0 | 2 | 6 | 2 | 1 | 0.136 |
| Fibroblast | 83 | 2 | 0 | 0 | 0 | 0.024 | 0.202 | 840 | 28 | 80 | 77 | 0 | 61 | 0.187 |
| KC_Basal | 385 | 1 | 89 | 0 | 65 | 0.234 | -1.054 | 4808 | 156 | 155 | 2493 | 102 | 2011 | 0.551 |
| KC_Cyc | 239 | 3 | 23 | 2 | 6 | 0.109 | -0.249 | 2459 | 39 | 127 | 899 | 83 | 401 | 0.417 |
| KC_Diff | 881 | 9 | 247 | 3 | 177 | 0.291 | -1.565 | 10496 | 328 | 398 | 5591 | 171 | 4704 | 0.571 |
| LC | 221 | 12 | 2 | 9 | 0 | 0.063 | 0.058 | 3072 | 113 | 462 | 257 | 316 | 12 | 0.234 |
| Mac | 228 | 91 | 0 | 70 | 0 | 0.399 | 1.546 | 3020 | 72 | 2033 | 21 | 1813 | 13 | 0.680 |
| MDSC | 62 | 0 | 1 | 0 | 0 | 0.016 | 0.031 | 711 | 23 | 58 | 73 | 42 | 0 | 0.184 |
| Melanocyte | 55 | 0 | 5 | 0 | 1 | 0.091 | -0.755 | 732 | 22 | 122 | 242 | 114 | 137 | 0.497 |
| NK | 3 | 0 | 0 | 0 | 0 | 0.00000000 | 0.636 | 111 | 4 | 16 | 1 | 4 | 0 | 0.153 |
| PDC | 17 | 12 | 0 | 10 | 0 | 0.706 | 3.569 | 245 | 10 | 219 | 0 | 202 | 0 | 0.894 |
| Pilosebaceous | 167 | 0 | 9 | 0 | 6 | 0.054 | -0.376 | 1988 | 23 | 46 | 719 | 15 | 504 | 0.385 |
| Tcell | 99 | 10 | 1 | 3 | 1 | 0.111 | 0.319 | 1327 | 32 | 609 | 164 | 342 | 147 | 0.583 |
| TSK | 26 | 11 | 0 | 7 | 0 | 0.423 | 3.131 | 322 | 5 | 237 | 3 | 185 | 1 | 0.745 |

k-distance weighting, neighbourhoods with SpatialFDR < 0.1 (up/down) per cell type: ASDC 18/0; B Cell 9/0; CD1C 51/2; CLEC9A 3/0; Endothelial Cell 0/2; Eccrine 0/0; Fibroblast 2/0; KC_Basal 2/96; KC_Cyc 3/24; KC_Diff 9/252; LC 14/3; Mac 94/0; MDSC 0/1; Melanocyte 0/5; NK 0/0; PDC 12/0; Pilosebaceous 0/10; Tcell 10/1; TSK 11/0

#### Run `pca|cond`: route pca, design `~condition`, SpatialFDR weighting graph-overlap (primary)

Neighbourhoods: 2799; SpatialFDR < 0.1: 2390 (up 1045, down 1345); SpatialFDR < 0.05: 2258; with k-distance weighting: FDR < 0.1: 2397, < 0.05: 2267; cells with milo.fdr < 0.1: 33223 / 35640; dominant-type fraction median 0.99

| celltype | n.nhoods | nh.fdr10.up | nh.fdr10.down | nh.fdr05.up | nh.fdr05.down | frac.nh.fdr10 | median.logFC | n.cells | n.cells.no.nhood | cells.fdr10.up | cells.fdr10.down | cells.fdr05.up | cells.fdr05.down | frac.cells.fdr10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ASDC | 21 | 20 | 0 | 20 | 0 | 0.952 | 4.717 | 283 | 5 | 254 | 11 | 254 | 10 | 0.936 |
| B Cell | 3 | 3 | 0 | 3 | 0 | 1.0000000 | 4.334 | 205 | 14 | 184 | 3 | 184 | 3 | 0.912 |
| CD1C | 347 | 219 | 62 | 201 | 57 | 0.810 | 3.509 | 3985 | 55 | 2646 | 1106 | 2545 | 1046 | 0.942 |
| CLEC9A | 33 | 16 | 3 | 13 | 3 | 0.576 | 2.375 | 556 | 8 | 419 | 77 | 360 | 77 | 0.892 |
| Endothelial Cell | 30 | 0 | 7 | 0 | 4 | 0.233 | -1.303 | 421 | 1 | 0 | 189 | 0 | 133 | 0.449 |
| Eccrine | 4 | 0 | 3 | 0 | 2 | 0.750 | -2.785 | 59 | 1 | 1 | 43 | 0 | 43 | 0.746 |
| Fibroblast | 69 | 56 | 13 | 56 | 12 | 1.0000000 | 5.280 | 840 | 7 | 615 | 218 | 615 | 213 | 0.992 |
| KC_Basal | 363 | 97 | 261 | 93 | 259 | 0.986 | -5.722 | 4808 | 86 | 1300 | 3397 | 1286 | 3391 | 0.977 |
| KC_Cyc | 211 | 51 | 133 | 46 | 124 | 0.872 | -3.593 | 2459 | 67 | 700 | 1600 | 689 | 1534 | 0.935 |
| KC_Diff | 852 | 188 | 634 | 181 | 630 | 0.965 | -5.816 | 10496 | 166 | 2240 | 8042 | 2226 | 8027 | 0.980 |
| LC | 223 | 53 | 63 | 42 | 49 | 0.520 | -0.124 | 3072 | 54 | 1279 | 1203 | 1173 | 1088 | 0.808 |
| Mac | 253 | 156 | 64 | 144 | 62 | 0.870 | 3.641 | 3020 | 45 | 2109 | 748 | 2043 | 742 | 0.946 |
| MDSC | 56 | 10 | 26 | 8 | 25 | 0.643 | -1.874 | 711 | 6 | 246 | 355 | 208 | 350 | 0.845 |
| Melanocyte | 53 | 31 | 12 | 31 | 10 | 0.811 | 3.744 | 732 | 3 | 456 | 235 | 456 | 214 | 0.944 |
| NK | 6 | 3 | 0 | 2 | 0 | 0.500 | 2.547 | 111 | 4 | 79 | 4 | 59 | 4 | 0.748 |
| PDC | 16 | 16 | 0 | 15 | 0 | 1.0000000 | 3.732 | 245 | 3 | 242 | 0 | 242 | 0 | 0.988 |
| Pilosebaceous | 145 | 68 | 45 | 60 | 41 | 0.779 | 1.271 | 1988 | 35 | 1056 | 784 | 1000 | 770 | 0.926 |
| Tcell | 92 | 37 | 19 | 28 | 16 | 0.609 | 1.623 | 1327 | 26 | 690 | 378 | 613 | 341 | 0.805 |
| TSK | 22 | 21 | 0 | 21 | 0 | 0.955 | 5.371 | 322 | 7 | 314 | 0 | 313 | 0 | 0.975 |

k-distance weighting, neighbourhoods with SpatialFDR < 0.1 (up/down) per cell type: ASDC 20/0; B Cell 3/0; CD1C 219/63; CLEC9A 16/3; Endothelial Cell 0/7; Eccrine 0/3; Fibroblast 56/13; KC_Basal 98/261; KC_Cyc 51/134; KC_Diff 188/634; LC 53/63; Mac 157/64; MDSC 10/26; Melanocyte 31/12; NK 3/0; PDC 16/0; Pilosebaceous 68/46; Tcell 39/19; TSK 21/0

#### Run `pca|patient_cond`: route pca, design `~patient + condition`, SpatialFDR weighting graph-overlap (primary)

Neighbourhoods: 2799; SpatialFDR < 0.1: 2396 (up 1007, down 1389); SpatialFDR < 0.05: 2312; with k-distance weighting: FDR < 0.1: 2398, < 0.05: 2316; cells with milo.fdr < 0.1: 33632 / 35640; dominant-type fraction median 0.99

| celltype | n.nhoods | nh.fdr10.up | nh.fdr10.down | nh.fdr05.up | nh.fdr05.down | frac.nh.fdr10 | median.logFC | n.cells | n.cells.no.nhood | cells.fdr10.up | cells.fdr10.down | cells.fdr05.up | cells.fdr05.down | frac.cells.fdr10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ASDC | 21 | 20 | 0 | 20 | 0 | 0.952 | 4.573 | 283 | 5 | 253 | 12 | 253 | 12 | 0.936 |
| B Cell | 3 | 3 | 0 | 3 | 0 | 1.0000000 | 3.679 | 205 | 14 | 184 | 4 | 184 | 3 | 0.917 |
| CD1C | 347 | 201 | 68 | 187 | 64 | 0.775 | 2.643 | 3985 | 55 | 2539 | 1204 | 2505 | 1199 | 0.939 |
| CLEC9A | 33 | 14 | 6 | 12 | 3 | 0.606 | 1.345 | 556 | 8 | 379 | 126 | 363 | 81 | 0.908 |
| Endothelial Cell | 30 | 0 | 9 | 0 | 5 | 0.300 | -1.323 | 421 | 1 | 1 | 268 | 0 | 181 | 0.639 |
| Eccrine | 4 | 0 | 4 | 0 | 4 | 1.0000000 | -2.963 | 59 | 1 | 0 | 50 | 0 | 49 | 0.847 |
| Fibroblast | 69 | 56 | 11 | 56 | 11 | 0.971 | 4.618 | 840 | 7 | 615 | 181 | 615 | 181 | 0.948 |
| KC_Basal | 363 | 97 | 260 | 94 | 260 | 0.983 | -5.124 | 4808 | 86 | 1298 | 3404 | 1291 | 3401 | 0.978 |
| KC_Cyc | 211 | 55 | 136 | 53 | 132 | 0.905 | -3.672 | 2459 | 67 | 725 | 1628 | 712 | 1620 | 0.957 |
| KC_Diff | 852 | 190 | 627 | 187 | 620 | 0.959 | -5.038 | 10496 | 166 | 2244 | 8022 | 2239 | 8006 | 0.978 |
| LC | 223 | 50 | 91 | 46 | 86 | 0.632 | -0.528 | 3072 | 54 | 1243 | 1511 | 1206 | 1466 | 0.896 |
| Mac | 253 | 158 | 57 | 153 | 55 | 0.850 | 3.480 | 3020 | 45 | 2161 | 711 | 2136 | 705 | 0.951 |
| MDSC | 56 | 12 | 28 | 10 | 27 | 0.714 | -1.621 | 711 | 6 | 291 | 371 | 272 | 361 | 0.931 |
| Melanocyte | 53 | 10 | 17 | 6 | 17 | 0.509 | 0.724 | 732 | 3 | 267 | 322 | 215 | 322 | 0.805 |
| NK | 6 | 2 | 0 | 2 | 0 | 0.333 | 1.505 | 111 | 4 | 68 | 4 | 68 | 4 | 0.649 |
| PDC | 16 | 16 | 0 | 15 | 0 | 1.0000000 | 3.685 | 245 | 3 | 242 | 0 | 242 | 0 | 0.988 |
| Pilosebaceous | 145 | 60 | 53 | 59 | 49 | 0.779 | 0.486 | 1988 | 35 | 959 | 898 | 955 | 864 | 0.934 |
| Tcell | 92 | 42 | 22 | 35 | 20 | 0.696 | 1.412 | 1327 | 26 | 736 | 398 | 719 | 391 | 0.855 |
| TSK | 22 | 21 | 0 | 21 | 0 | 0.955 | 5.264 | 322 | 7 | 313 | 0 | 313 | 0 | 0.972 |

k-distance weighting, neighbourhoods with SpatialFDR < 0.1 (up/down) per cell type: ASDC 20/0; B Cell 3/0; CD1C 202/68; CLEC9A 14/6; Endothelial Cell 0/9; Eccrine 0/4; Fibroblast 56/11; KC_Basal 97/260; KC_Cyc 55/136; KC_Diff 190/627; LC 50/91; Mac 158/57; MDSC 12/28; Melanocyte 10/17; NK 2/0; PDC 16/0; Pilosebaceous 60/53; Tcell 42/23; TSK 21/0

Runtime (s): load_inputs 21.4; conos_build_graph 1.6; conos_makeNhoods 12.9; conos_countCells 0.2; conos_cond_testNhoods_graph-overlap 1.1; conos_cond_testNhoods_k-distance 885.5; conos_patient_cond_testNhoods_graph-overlap 4.8; conos_patient_cond_testNhoods_k-distance 876.4; pca_seurat_pca 45.0; pca_buildGraph 5.8; pca_makeNhoods 11.6; pca_countCells 0.2; pca_cond_testNhoods_graph-overlap 0.9; pca_cond_testNhoods_k-distance 4023.6; pca_patient_cond_testNhoods_graph-overlap 5.6; pca_patient_cond_testNhoods_k-distance 4408.7

