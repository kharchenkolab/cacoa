# Cluster-free multiplicity adjustment: benchmark on the planted density simulation (2026-10-11)

Motivation: the max-statistic (familywise) adjustment of the cluster-free density and expression-shift z-scores left no
significant cell in any real dataset (AD, MS, SCC) and Milo's neighbourhood testing with spatial FDR looked more
sensitive. Alternatives implemented in `adjustCellScores()` (`R/cell_adjust.R`) and offered through
`estimateDiffCellDensity(adjust.method = )` and `estimateClusterFreeExpressionShifts(adjust.method = )`:
`"bh"` (BH over tested cells), `"spatial-fdr"` (Milo's weighted BH, weights 1/(1 + degree)), `"neighbourhood"` (a
covering set of index cells is tested, weighted BH over the index neighbourhoods, every cell inherits the smallest
adjusted p of the neighbourhoods containing it; Milo-like without its NB model), `"max-t"` (the former default).
Per-cell p-values: density uses the parametric t of the contrast (well calibrated); the cluster-free shift uses the
permutation p (its parametric F is not calibrated: with the parametric p nothing is significant at all).

Simulation: `test/sim_density.rds` from `misc/notebook/prepare_density_sim.R` (20 samples; 60 % of the treated IN-SST
and of the control OPC cells removed). Note that the base simulation already differs in composition between the arms
(per-sample proportions, Wilcoxon): L2/3 and L5/6-CC at 0.47x in treated (p 2e-4, 6e-4), IN-PV 1.8x (p 0.001), IN-SST
1.8x before the removal (p 2e-4; 0.72x after), AST-FB 1.35x (p 0.026), OPC 1.34x before (p 0.017; up after),
AST-PP / Microglia / Neu-NRGN unchanged. Expression truth: L2/3, IN-PV, IN-SST strong; OPC, L5/6-CC, AST-PP
moderate; AST-FB weak; Microglia, Neu-NRGN none.

Density (999 permutations; fraction of cells with |z.adj| > 1.96):

| cell type | change | max-t | bh | spatial-fdr | neighbourhood |
|---|---|---|---|---|---|
| IN-SST | down | 0.00 | 1.00 (down) | 1.00 | 1.00 |
| OPC | up | 1.00 | 1.00 (up) | 1.00 | 1.00 |
| L2/3 | down (0.47x) | 0.00 | 1.00 (down) | 1.00 | 1.00 |
| L5/6-CC | down (0.47x) | 0.00 | 1.00 (down) | 1.00 | 1.00 |
| IN-PV | up (1.8x) | 0.00 | 1.00 (up) | 1.00 | 1.00 |
| AST-FB | up (1.35x) | 0.00 | 0.00 | 0.00 | 0.68 |
| AST-PP | none | 0.00 | 0.00 | 0.00 | 0.05 |
| Microglia | none | 0.00 | 0.04 | 0.04 | 0.60 |
| Neu-NRGN | none | 0.00 | 0.00 | 0.00 | 0.32 |

BH and spatial FDR recover every change with p < 0.01 in the per-sample test and stay quiet in the unchanged types;
max-t finds one of six; the neighbourhood variant (45 index neighbourhoods) paints whole unchanged populations.
Runtime: max-t 7 s, the others < 1 s.

Cluster-free expression shifts (500 genes, 499 permutations, permutation p; fraction with one-sided z.adj > 1.64):

| cell type | planted | max-t | bh | spatial-fdr | neighbourhood |
|---|---|---|---|---|---|
| IN-PV | strong | 0.58 | 0.91 | 0.92 | 1.00 |
| L2/3 | strong | 0.44 | 0.80 | 0.82 | 1.00 |
| IN-SST | strong | 0.42 | 0.78 | 0.80 | 1.00 |
| AST-PP | moderate | 0.00 | 0.77 | 0.77 | 1.00 |
| OPC | moderate | 0.00 | 0.80 | 0.80 | 1.00 |
| L5/6-CC | moderate | 0.07 | 0.51 | 0.55 | 0.25 |
| AST-FB | weak | 0.00 | 0.19 | 0.21 | 0.98 |
| Microglia | none | 0.00 | 0.24 | 0.26 | 0.98 |
| Neu-NRGN | none | 0.00 | 0.09 | 0.09 | 0.32 |

The FDR adjustments find the moderate shifts that max-t misses entirely. The Microglia false positives (24 %) sit at the
AST-FB boundary of the graph: neighbourhoods that mix cell types whose composition differs between the arms carry a
"shift" (the per-cell permutation p there has median 0.11), a known confound of neighbourhood-level shifts that the
cluster-based test does not have. Defaults set 2026-10-11: density `"bh"` with parametric p; shifts `"bh"` with
permutation p; `"max-t"` kept as an option. Head-to-head with Milo on the same simulation and on MS / SCC: see
`misc/milo_comparison.md` (pending).

Addendum (same day): (1) `p.source = "z"` for the cluster-free shifts (normal tail of the permutation-standardized z,
no floor) raises sensitivity (strong 79-91 %, moderate 72-84 %) but is anti-conservative on the null types (Microglia
37 %, Neu-NRGN 15 % at a nominal 5 %), so the default stays `"permutation"`; on MS (30,000 cells, sparse signal) no
cell passes BH with the permutation p at 199-999 permutations because too few cells sit at the floor (108 at 999),
the notebook shows the raw z map instead. (2) Density: a moderated t (limma `squeezeVar` across cells, `moderate =
TRUE`) is the default parametric p; it leaves the simulation and MS results unchanged and does not rescue SCC. See
`misc/milo/milo_comparison.md` for the head-to-head with Milo.

## Calibration and power benchmark (2026-10-11, second pass): see `misc/benchmark/README.md`

Added to the package after the benchmark: `adjust.method = "perm-fdr"` (pooled permutation FDR, `permutationFDR()`;
now the default of both cluster-free tests), `estimateCellDensity(smoothing = "neighbourhood")` with
`estimateDiffCellDensity(transform = "log")` (1-hop neighbourhood counts on the log scale, the unit of a Milo-type
test), and the per-cell-type summaries `estimateClusterFreeSummary()` / `plotClusterFreeSummary()` (exact group
permutation test from the shared permutations, `groupPermutationTest()`, or the Cauchy combination of the cells'
p-values, `cauchyCombination()`). On the planted simulation the pooled permutation FDR calls the same cells as BH for the
density (every change with p < 0.01 in the per-sample test, nothing in the unchanged types) and, for the shifts,
fractions within a few percent of BH on the permutation p (strong 82-94 %, moderate 57-80 %, weak 18 %, null types
Microglia 24 % at the AST-FB boundary and Neu-NRGN 9 %). The group permutation test marks the six changed types of the
density simulation and none of the three unchanged ones (p 0.001-0.002 versus 0.40-0.70).

Kernel fix (2026-10-11, from the null study): the shift kernel induced each shared permutation within a cell's present
samples, which keeps the per-cell test exact but makes the relabelings incoherent across cells and the group
permutation test anti-conservative (9-22 % of the types at the 5 % level under relabeling; 31 % on iid synthetic data
with 6 of 16 samples absent). Permutations of the labels are now global relabelings restricted to each cell's present
samples (cells are skipped under relabelings that leave too few samples per level; `n.perm` per cell), Freedman-Lane
keeps the within-cell induction. The low-level `clusterFreeExpressionShifts()`, `testTermEffects()` and
`testPairwiseEffects()` now default `block.vars` to the model's, as the class methods already did. Details and the
synthetic checks in `misc/benchmark/README.md`.
