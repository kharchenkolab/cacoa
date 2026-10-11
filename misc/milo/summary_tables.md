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

## Dataset: ms

Parameters: prop = 0.1; k (conos route, .k for k-distance weighting) = 30; k (pca route) = 30, d = 30; condition column = diagnosis (reference = Control); miloR 1.6.0

- route conos: 3339 neighbourhoods; size median 150 (IQR 88-212, range 8-438); cells in >= 1 neighbourhood: 46545 / 48918; graph degree median 58
- route pca: 3424 neighbourhoods; size median 64 (IQR 50-86.25, range 31-396); cells in >= 1 neighbourhood: 48343 / 48918; graph degree median 39

### Run `conos|cond`: route conos, design `~condition`, SpatialFDR weighting graph-overlap (primary)

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

### Run `pca|cond`: route pca, design `~condition`, SpatialFDR weighting graph-overlap (primary)

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

## Dataset: scc

Parameters: prop = 0.1; k (conos route, .k for k-distance weighting) = 30; k (pca route) = 30, d = 30; condition column = condition (reference = Normal); miloR 1.6.0

- route conos: 2945 neighbourhoods; size median 106 (IQR 67-154, range 11-273); cells in >= 1 neighbourhood: 34618 / 35640; graph degree median 61
- route pca: 2799 neighbourhoods; size median 58 (IQR 48-72, range 31-176); cells in >= 1 neighbourhood: 35047 / 35640; graph degree median 39

### Run `conos|cond`: route conos, design `~condition`, SpatialFDR weighting graph-overlap (primary)

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

### Run `conos|patient_cond`: route conos, design `~patient + condition`, SpatialFDR weighting graph-overlap (primary)

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

### Run `pca|cond`: route pca, design `~condition`, SpatialFDR weighting graph-overlap (primary)

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

### Run `pca|patient_cond`: route pca, design `~patient + condition`, SpatialFDR weighting graph-overlap (primary)

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

