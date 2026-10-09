# Fit and randomization: where each test does it, and how to converge on two kernels

State as of 2026-10-09 (dev_lm, after the Track B-E session). The original strategy was one Rcpp fit/randomize core
used by every test, for speed, uniformity and simplicity. The redesign did not land there. This note lists where
fitting and permutation happen today, why it diverged, what it costs, and a concrete path to two C++ kernels fed by
one permutation source.

## 1. Inventory

| # | analysis (entry point) | statistic | observed fit | permutations drawn by | permuted statistic computed in | schemes | file |
|---|---|---|---|---|---|---|---|
| 1 | expression shifts, contrast test (`estimateExpressionShiftMagnitudes`) | shift F, shift, var, total from the Gower matrix G | R `estimatePairwiseEffects` (+ `fitDispersion`) | R `permutationPlan` / `drawPermutations` (restricted: contrasted levels only, within strata; exhaustive when small) | **R** `permutedStats` per permutation (block: relabel X; FL: reindex K1..K4). C++ `permuted_contrast_F(_fl)` is used **only when var/total are not needed**, which is not the default path | block, FL, HJ (HJ in R `hjStat`) | `pairwise_inference.R`, `perm_stats.cpp` |
| 2 | expression shifts, whole-factor test (same entry) | location F (X_full vs X_reduced), dispersion F | R `termTestGower`, `dispersionTermTest` | same as 1 | **R** `termPermutationStats` (per permutation: relabel or reindex, then the two term tests) | block, FL (HJ falls back to relabeling) | `expression_shifts_api.R` |
| 3 | covariate screen (`screenCovariates`, `estimateMetadataSeparation`) | term F and dispersion F per covariate x cell type, marginal and partial | R, same as 2 | global P from R (`drawPermutations`) induced per complete-case subset (`inducePermutationSimple`), or `sample.int` | **R** loop inside `screenOneMatrix` (FL only) | FL only | `covariate_screen.R` |
| 4 | sensitivity (`checkSensitivity`) | rows 1-2 re-run per covariate set; leave-one-out refits are analytic (no permutation) | R | as 1 | as 1-2 | as 1-2 | `sensitivity.R` |
| 5 | cluster-free expression shifts (`estimateClusterFreeExpressionShifts`) | shift F per cell neighbourhood (+ shift estimate), max-statistic z adjustment | distances: C++ `estimateExpressionShiftsPairsLM`; shift estimate: R `estimatePairwiseEffects` per cell | global P from R, induced per neighbourhood **in R** (`apply(P, 2, inducePermutation)` per cell) | **C++** `permuted_contrast_F(_fl)` per cell, called from an R loop over cells (`oneCell`); max-stat: C++ `adjustedZScoresMaxStat`; smoothing: C++ `applyMedianFilterES` | block, FL | `cluster_free_shifts.R` |
| 6 | composition, contrast (`estimateCellLoadings`) | per-ILR-coordinate linear model, contrast statistic, loadings | **C++** `fit_and_randomize` (block) or `fl_fwl_cpp` (FL) via `performLMPermutations` | **C++ fitter's own RNG** over `perm_groups` (all samples swapped within blocks; built in R by `makeBlocks`/`permutationGroups`) | C++ | block, FL (feature-level FL) | `coda.R`, `model_fits.R`, `lm_fit.cpp` |
| 7 | composition, whole-factor | row 2 on Euclidean distances of the ILR coordinates | R | as 1 | R | block, FL | `cacoa.R` (`estimateCellLoadings`) |
| 8 | cell density (`estimateDiffCellDensity`) | per-bin / per-cell linear model, contrast z | C++ as 6 (default FL) | C++ fitter's own RNG | C++ | block, FL | `cell_density.R` |
| 9 | cluster-free DE (`estimateClusterFreeDE`) | per-gene per-neighbourhood linear model, contrast z | C++ as 6 | C++ fitter's own RNG | C++; z adjustment `adjustZScoresByPermutations` | block, FL | `cluster_free.R` |
| 10 | DE per cell type (`estimateDEPerCellType`) | DESeq2 / edgeR / limma | parametric | none | none | - | `de_function.R` |
| - | unused C++ exports | `fit_density_lm`, `perm_FL_contrast_mat`, `perm_full_contrast_mat` (`projdiff.cpp`); `pca_project`, `estimateCorrelationDistance` (`expression_shifts.cpp`) | | | | | |

Two permutation generators exist and they do **not** draw from the same distribution: the R plan swaps only the
samples of the compared levels within strata (exact for a two-level contrast with discrete nuisance, and the one
whose p-value floors and enumeration counts are reported in the model), while `fit_and_randomize` shuffles all
samples within each `perm_groups` block with its own Mersenne stream. With more than two levels (the simulated
object has four groups) rows 6, 8 and 9 therefore permute differently from rows 1-5, and their permutations are not
coupled with the shift tests.

## 2. Why it diverged

- The redesigned statistics are quadratic forms of an n x n Gower matrix (a'Ga, tr(HG), diag(RGR) for the
  dispersion fit), not per-column linear-model fits of a feature matrix. `fit_and_randomize` fits Y (n x m) by
  columns and returns a contrast t/z per column; it cannot compute the distance statistics, so the only reusable
  part was the design. The handoff's R reference (`proposed.R`) was ported first, in R, and validated; only the
  hottest scalar (the shift F) was then moved to C++.
- Freedman-Lane exists twice with different mathematics: Gower-form reindexing of the residual kernel (rows 1-5)
  and residual permutation at the feature level (rows 6, 8, 9).
- The permutation plan (restricted swaps, strata, enumeration, induced sub-permutations) was written in R because
  the model object needs it before any fit (distinct-permutation counts, p-value floors); the C++ fitter never
  learned to consume it. That is the "fit_and_randomize hook" left undone.

## 3. What it costs

Per cell type, n = 40 samples, q = 3 columns, 999 permutations (this machine, single core):

| path | time |
|---|---|
| R `permutedStats`, shift + var + total (default cell-type test) | 0.33 s |
| R `permutedStats`, F only | 0.04 s |
| C++ `permuted_contrast_F`, F only | 0.007 s |
| R term F (`termTestGower` loop) | 0.72 s |
| R dispersion term F (`dispersionTermTest` loop) | 0.95 s |

Per-permutation work is O(n^2 q) everywhere, so for n <= 100 the kernels themselves are microseconds; what
costs is R call overhead multiplied by the number of (test x cell type x permutation) units:

- Cell-type tests: 9 cell types x 999 permutations stay under a few seconds. Not a problem.
- Screen: covariates x cell types x 2 modes x (term + dispersion) R loops; 5 covariates x 9 cell types at 999
  permutations is a few minutes. Acceptable, but it is the slowest exploratory step.
- Cluster-free shifts: cells x (R induction of P over B columns + C++ kernel + R shift estimate). Measured in the
  notebook render: 4,500 cells x 99 permutations x 500 genes took 8.7 minutes on 16 cores (the C++ kernel itself
  accounts for seconds of that); the R per-cell overhead dominates. **This is the real performance problem.**
- Rows 6, 8, 9 are fast (C++ end to end) but on a different permutation distribution.

## 4. Target: two kernels, one permutation source

The two statistic families are different mathematical objects, so one routine is not achievable without
pretending; two is the honest minimum. Everything else in R is design, plan, estimability, result assembly,
multiplicity and plots.

**P (R, exists).** `permutationPlan()` + `drawPermutations()` produce an integer n x B matrix for a sample set,
plus `inducePermutation()` for subsets. This becomes the single source of randomness for every test.

**Kernel A: `gower_perm_stats()` (C++, extend `perm_stats.cpp`).**
Inputs: G, X_full, X_reduced (for term tests), Z_full, Z_reduced, contrast c, dispersion endpoints, P, scheme
(block / FL / HJ), which statistics. Output: observed row + B x k matrix (shift F, shift, var, total, location F,
dispersion F). Pieces: the a'Ga / tr(HG) part exists; port `fitDispersion` (least squares of diag(RGR) on
(R o R) Z, q_z x q_z solve per permutation), `termTestGower`, `dispersionTermTest`, and `hjStat`. FL via the
K1..K4 reindexing already in C++. Add a batched entry `gower_perm_stats_batch()` taking a list of (G_k, sample
subset_k, strata_k) and the global P, inducing the sub-permutations in C++ and looping over units with the sccore
thread pool: this serves the cluster-free shifts (units = cells) and the screen (units = covariate x cell type).

**Kernel B: `fit_and_randomize()` (C++, exists).** Add a `perm_matrix` argument (n x B, 1-based); when given,
ignore `perm_groups` and the internal RNG, and let `fl_fwl_cpp` pass it through to its inner call. Delete
`projdiff.cpp`. `makeBlocks` / `permutationGroups` are then only needed for the legacy path and can go.

Consumers after convergence:

| rows | kernel | permutations |
|---|---|---|
| 1, 2, 4, 7 (cell-type and composition term tests, sensitivity) | A | P from the model's plan, shared across cell types |
| 3 (screen) | A, batched | one global P induced per subset |
| 5 (cluster-free shifts) | A, batched | one global P induced per neighbourhood, in C++ |
| 6, 8, 9 (composition contrast, density, cluster-free DE) | B with `perm_matrix` | the same P as rows 1-5 on the same object |
| 10 (DE) | none (parametric) | - |

The R implementations (`permutedStats`, the `termPermutationStats` and `screenOneMatrix` loops, `hjStat`) move
to `tests/testthat/helper-reference.R` and exist only to validate the kernels, as `test-track-d.R` already does
for the shift F (agreement to 1e-10).

## 5. Steps, each independently testable

1. **P hook for kernel B** (`perm_matrix` in `fit_and_randomize` / `fl_fwl_cpp`; `performLMPermutations(P =)`);
   route CoDA, density and cluster-free DE through `drawPermutations()` of the stored model. Tests: identical
   statistics with a given P versus the R block generator restricted to two levels; CoDA and shift tests on one
   object share `P`; exhaustive enumeration respected. Removes the "not coupled" leftover.
2. **Kernel A, contrast**: extend `permuted_contrast_F` to shift / var / total (port `fitDispersion`); switch the
   default cell-type path (`permutationStatsForCellType`) to it. Test versus `permutedStats` to 1e-10 in block and FL.
3. **Kernel A, terms**: location F and dispersion F for X_full vs X_reduced; switch `termPermutationStats` and
   `screenOneMatrix`. Test versus the R loops; E1 calibration test unchanged.
4. **Batched kernel** with C++ induction and the thread pool; replace `oneCell` in `cluster_free_shifts.R`
   (shift estimate per cell moves into the kernel too: it is M = A G A' with the dispersion correction, already
   computed for the permuted statistics). Test: identical to the current per-cell results on the toy; timing on
   the simulated object (target: 4,500 cells x 199 permutations in well under a minute).
5. **Retire**: `projdiff.cpp`, `pca_project`, the unused `estimateCorrelationDistance` export, `makeBlocks` /
   `permutationGroups`, R reference loops out of the namespace, `huh-jhun` either in the kernel or dropped.
6. Optional: fold `fl_fwl_cpp` into `fit_and_randomize` behind a `scheme` argument so each family has one entry point.

Steps 1-3 are each about half a day with tests; step 4 is the largest (a day) and the one with the visible payoff.
