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

Empirical truth (per-sample proportions, treated vs control): IN-SST, L2/3, L5/6-CC down; IN-PV, OPC, AST-FB up;
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

## SCC (Ji 2020)

Pending (Milo run in progress).

## Reading

- The former max-statistic adjustment is unusable on real data; BH over the per-cell parametric p-values is calibrated on
  the simulation and recovers the cluster-based picture on MS.
- Milo's advantage is in small populations where pooling a neighbourhood's counts gives power; its weakness on the same
  graph is the coverage of its calls (large overlapping neighbourhoods, minimum-FDR propagation), which paints whole
  unchanged populations at the 10 % threshold, and its dependence on an integrated embedding when built in its usual way.
- cacoa's per-cell test remains a per-cell contrast of smoothed densities; a neighbourhood-level count model (Milo's) with
  cacoa's design and permutation machinery would combine both, and is the natural next step if small populations matter.
