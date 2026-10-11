## Dataset: sim

Parameters: prop = 0.2; k (conos route, .k for k-distance weighting) = 30; k (pca route) = 30, d = 30; condition column = treatment (reference = control); miloR 1.6.0

- route conos: 293 neighbourhoods; size median 405 (IQR 345-442, range 185-567); cells in >= 1 neighbourhood: 2557 / 2557; graph degree median 341
- route pca: 140 neighbourhoods; size median 88 (IQR 73-109.25, range 48-182); cells in >= 1 neighbourhood: 2542 / 2557; graph degree median 39

### Run `conos|cond`: route conos, design `~condition`, SpatialFDR weighting graph-overlap (primary)

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

### Run `pca|cond`: route pca, design `~condition`, SpatialFDR weighting graph-overlap (primary)

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

## ms: no results file

## scc: no results file

