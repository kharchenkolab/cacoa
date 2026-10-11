# cacoa cluster-free density test versus Milo (miloR 1.6.0), 2026-10-11

Setup (scripts in this directory): Milo run by a subagent on the same Conos joint graph that cacoa uses (the graph
assigned as Milo's kNN graph, neighbourhoods by graph-refined sampling, `prop` 0.2 for the simulation and 0.1 for the
real data; design `~ condition`; TMM + quasi-likelihood F test; spatial FDR with graph-overlap weighting). A second
Milo route with its own PCA + kNN graph (no sample integration) is reported in the agent's notes
(`scratchpad/milo/milo_notes.md`, `summary_tables.md`). cacoa: graph density per sample (`estimateCellDensity(method =
"graph")`), per-cell contrast of the treated-versus-control densities with the stored model, 999 permutations, per-cell
p from the t statistic, adjusted by BH (new default) or by the max statistic (former default). Per-cell Milo results:
the smallest spatial FDR over the neighbourhoods containing the cell and the logFC of that neighbourhood.

## Planted density simulation (`test/sim_density.rds`; 2,557 cells, 10 vs 10 samples)

Empirical truth (per-sample proportions, treated vs control; now stored as `truth$composition` in `test/sim_density.rds`, Wilcoxon p < 0.01): IN-SST, L2/3, L5/6-CC down; IN-PV, OPC, AST-FB up;
AST-PP, Microglia, Neu-NRGN unchanged. Fractions of cells called:

| cell type | change | cacoa BH (|z| > 1.96) | cacoa max-t | Milo cells (FDR < 0.1) | Milo neighbourhoods (FDR < 0.1) |
|---|---|---|---|---|---|
| IN-SST | down | 1.00 (down) | 0.00 | 1.00 | 18/18 |
| L2/3 | down | 1.00 (down) | 0.00 | 1.00 | 22/36 |
| L5/6-CC | down | 1.00 (down) | 0.00 | 1.00 | 23/46 |
| IN-PV | up | 1.00 (up) | 0.00 | 1.00 | 7/30 |
| OPC | up | 1.00 (up) | 1.00 | 1.00 | 17/17 |
| AST-FB | up (1.35x) | 0.00 | 0.00 | 1.00 | 19/42 |
| AST-PP | none | 0.00 | 0.00 | 0.33 | 0/26 |
| Microglia | none | 0.04 | 0.00 | 1.00 | 7/42 |
| Neu-NRGN | none | 0.00 | 0.00 | 1.00 | 4/36 |

Cell-level Spearman between cacoa's raw z and Milo's logFC: 0.56. Milo's per-cell calls cover 92 % of all cells,
including the unchanged types, because the Conos-graph neighbourhoods are large (median 405 cells) and overlapping and
every cell inherits its most significant neighbourhood; at the neighbourhood level Milo still calls 17 % of the
Microglia and 11 % of the Neu-NRGN neighbourhoods. Milo's own PCA route (no integration) splits the shifted cell types
by condition and then calls neighbourhoods in both directions within IN-SST, L2/3 and L5/6-CC (see the agent's notes):
an unintegrated embedding confounds abundance with expression shifts. cacoa with BH recovers every change with
p < 0.01 in the per-sample test, misses the mild AST-FB change, and stays at or below 4 % in the unchanged types;
max-t finds one of six. Maps: `compare_sim_maps.png`.

## MS (Schirmer 2019; 48,918 nuclei, 12 MS vs 9 control)

Reference: the cluster-based loadings (neurons of all layers and interneurons down, Microglia / Phagocytes / T cells /
Stromal cells up), reproduced by cacoa and the preprint.

| cell type | loadings | cacoa BH (fraction of cells) | Milo cells FDR < 0.1 | Milo neighbourhoods FDR < 0.1 |
|---|---|---|---|---|
| EN-L2-3 | down | 0.71 (down) | 0.26 | 0.05 |
| EN-L4 | down | 0.70 (down) | 0.38 | 0.05 |
| EN-L5-6 | down | 0.41 (down) | 0.13 | 0.03 |
| EN-MIX | down | 0.35 (down) | 0.18 | 0.00 |
| IN-PV / IN-SST / IN-SV2C / IN-VIP | down | 0.00-0.004 | 0.00-0.24 | 0.00-0.01 |
| Microglia | up | 0.76 (up) | 0.96 | 0.70 |
| Phagocytes | up | 0.72 (up) | 1.00 | 0.96 |
| Astrocytes | n.s. | 0.17 (up) | 0.53 | 0.19 |
| OL | n.s. | 0.01 | 0.45 | 0.13 |
| OPC | n.s. | 0.003 | 0.47 | 0.11 |
| Stromal cells | up | 0.00 | 0.92 | 0.33 |

Overall: cacoa BH calls 23 % of the nuclei, max-t none, Milo 35 % at FDR < 0.1 and none at FDR < 0.05 (every Milo call
lies between 0.05 and 0.1). Cell-level Spearman between cacoa's z and Milo's logFC: 0.69. Milo recovers the glial /
immune gains but not the neuronal depletion (5 % of the excitatory-neuron neighbourhoods), which the cluster-based
analysis, the original study and cacoa's BH map all show; cacoa misses the interneuron and stromal changes (small
populations: 105-1,900 nuclei spread over 21 samples; the per-cell density contrast has little power there). Maps:
`compare_ms_maps.png`.

## SCC (Ji 2020; 35,640 cells, 9 tumour vs 9 normal, 2,945 Milo neighbourhoods)

Reference: tumour-only populations (TSK, PDC, B cells, ASDC) in the loadings. cacoa's per-cell density test calls
nothing at 5 % FDR with any adjustment (BH on the parametric t, on the moderated t, on the permutation p; max-t):
the smallest BH-adjusted p is 0.07. Milo calls 58 % of the cells at FDR < 0.1 and 49 % at 0.05, with the tumour-only
populations on top (PDC 99 %, B cells 96 %, ASDC 94 %, TSK 89 % of cells; 100 %, 92 %, 65 %, 69 % of their
neighbourhoods) and the keratinocyte states down (KC_Diff 41 % of neighbourhoods). Cell-level Spearman between
cacoa's z and Milo's logFC: 0.82 (the two statistics agree on direction and order; the difference is power).

Why cacoa loses here: (1) the graph density is a heat-kernel smooth (beta = 30) of each sample's cell indicator, so a
small population embedded in a large one (TSK inside KC_Basal) has its tumour-only signal diluted by the normal
neighbours (TSK median raw z 1.1, PDC 3.3); (2) the per-cell contrast rests on 18 sample densities, so even a
perfectly separated population reaches only the permutation floor (0.001 at 999 permutations; 310 cells at the floor
in SCC), and BH over 35,640 cells needs about 700 such cells before any passes; the parametric t is no better (zero
densities in one arm and heterogeneous tumours give t about 3). Shrinking the per-cell variances across cells (a
moderated t, limma's empirical Bayes, now `moderate = TRUE` by default) changes nothing here: the limitation is the
statistic, not its variance estimate. Milo pools counts per neighbourhood and borrows dispersion across
neighbourhoods, and tests about 3,000 neighbourhoods instead of 35,000 cells. Maps: `compare_scc_maps.png`.

## Reading

- The former max-statistic adjustment is unusable on real data; BH over the per-cell parametric p-values is calibrated on
  the simulation and recovers the cluster-based picture on MS.
- Milo's advantage is in small populations where pooling a neighbourhood's counts gives power; its weakness on the same
  graph is the coverage of its calls (large overlapping neighbourhoods, minimum-FDR propagation), which paints whole
  unchanged populations at the 10 % threshold, and its dependence on an integrated embedding when built in its usual way.
- cacoa's per-cell test remains a per-cell contrast of smoothed densities: it wins on the simulation and on MS (large
  populations, calibrated calls) and loses on SCC (small condition-specific populations among many cells). A
  neighbourhood-level count model on index cells (Milo's unit of testing) with cacoa's design and permutation
  machinery, and a smaller smoothing bandwidth for the density, are the natural next steps; the present
  `adjust.method = "neighbourhood"` only reduces the number of tests (it reuses the per-cell smoothed statistics) and is
  liberal on the simulation.

Milo agent's full write-up: `milo_notes.md` (methods, per-cell-type tables for every route and design, runtimes) and
`summary_tables.md`. Its two caveats stand: the simulation's Conos graph is dense (median degree 341 of 2,557 cells), so
Milo's neighbourhoods there are nearly cluster-sized (median 405 cells, 72 % type purity); and Milo's PCA route
without integration is uninformative on all three datasets (the expression shifts split the planted types into
condition-pure neighbourhoods), so the Conos-graph route is the fair comparison. The k-distance spatial FDR on the
2-D embedding is unusable (weights span orders of magnitude; 1 neighbourhood called on MS); graph-overlap is used.
