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

Two permutation generators exist. Both restrict the swaps to the samples of the compared levels within nuisance
strata (see §2a), but they are not the same draws: the R plan draws one matrix P per test, shared by all cell
types and coordinates, enumerates when few relabelings exist and reports the p-value floor in the model; the C++
fitter draws its own Mersenne stream per response column and never enumerates. Rows 6, 8 and 9 are therefore not
coupled with the shift tests, nor with each other's columns.

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


## 2a. Why the randomization logic differs (in detail)

The two generators answer different questions and were written at different times for different statistics.

**C++ `fit_and_randomize` (rows 6, 8, 9).** The design matrix F and the response matrix Y (n x m) are fixed; for
each column j a job draws `n_randomizations` permutations and refits. Three properties matter:

1. *Which samples move.* `perm_groups` are built by `permutationGroups(blocks, core.rows)`: the blocks are the
   interaction of the discrete nuisance factors (`makeBlocks()`), and only the **core rows** are placed in them,
   where `core.rows` marks the samples with non-zero weight in the contrast (Group1 and Group2 for a
   Group2-vs-Group1 test; verified for a 4 x 2 design: two blocks of 5 + 5 Group1/Group2 samples, one per batch).
   `generate_permutation()` starts from the identity and shuffles inside each block, so Group3 and Group4 rows keep
   their labels. The fitter's block scheme is therefore contrast-restricted, like the R plan, and an earlier
   version of this note claiming otherwise was wrong. For a whole-factor (term) test the fitter has no equivalent
   (it tests one contrast at a time), and for interaction-cell contrasts the restriction follows the contrast's
   weights rather than the cell labels, which can differ from the R plan's `in.set`; otherwise the two agree on
   which samples move.
2. *One RNG stream per column.* `make_rng(seed, j)` seeds the Mersenne generator with the column index, so
   permutation b of column 1 is not permutation b of column 2. Anything that combines columns at a fixed
   permutation index assumes they were relabelled together, and that assumption is false here:
   - `lmCoda()` forms the permuted loadings as `psi %*% t(beta.contrast.perm)`, i.e. row b of the permuted ILR
     coefficients across all K coordinates is back-transformed as one vector. Those K values come from K different
     relabelings, so the null of each cell-type loading ignores the correlation between ILR coordinates, and the
     global statistic `T.contrast.perm` (variance explained across coordinates) has the same problem.
   - `estimateDiffCellDensity()` adjusts z-scores with `adjustZScoresByPermutations(score, permut.scores)`, a
     max-statistic over bins per permutation row; neighbouring bins are strongly correlated, but the row mixes
     independent draws, so the max-null is the one for independent bins (conservative when statistics are
     positively correlated, and in any case not the permutation distribution of the maximum).
   This is the most consequential difference and the strongest reason for step 1 below.
3. *No enumeration, no plan.* The fitter always samples with replacement from the relabelings, even when there are
   only a handful (2 vs 2 in a block gives 6). It cannot report the number of distinct permutations or the smallest
   attainable p-value, which the model object prints, and the Freedman-Lane variant (`fl_fwl_cpp`) permutes
   residuals of Y under the reduced model within `core_perm_groups` with the same per-column streams.

**R plan (rows 1-5, 7).** `permutationPlan()` is computed from the metadata before any fit, because the model
object has to describe the test (scheme, distinct permutations, p-value floor) when it is set. It draws one
integer matrix P (n x B) for the sample set and everything consumes it:

1. *Which samples move.* `contrastSampleInfo()` marks `in.set`: only the samples whose label is one of the two
   compared levels (for a cell contrast, the samples matching the fixed settings as well); for a whole-factor test
   every level of the factor. Strata are all discrete formula variables other than the swapped one, plus
   `block.vars`. Labels are swapped within (stratum x in.set) cells, so the other levels keep their labels and the
   null is specific to the contrast.
2. *Coupling.* P is drawn once per test and shared by all cell types (and by both tests of a two-test model). A
   cell type missing some samples gets the induced permutation (`inducePermutation()`: within each of its
   stratum cells the members are reassigned by the ranks of their images under P), which is uniform on the subset's
   relabelings and coupled with the full set. That is what makes the across-cell-type max-statistic and the global
   p-value valid. The cluster-free shifts use the same mechanism per neighbourhood.
3. *Exactness.* When the number of distinct relabelings is at most `n.permutations` (and 20,000) they are
   enumerated and the p-values are exact (the SCC paired design, 2^8 = 256). Otherwise Monte Carlo with the floor
   1/(B+1). Freedman-Lane in Gower form permutes all samples (strata from `block.vars` only) and re-indexes the
   residual kernel; Huh-Jhun is opt-in for the shift.

So the two sides agree on the restricted set of samples and on the strata for a simple contrast; they differ in
coupling (one shared P versus one stream per column), in enumeration and reporting (plan before the fit, exact
p-values when few relabelings exist), in the Freedman-Lane mechanics, and in the seed. The R side was written in R
because the plan has to exist before any fit and because its statistics are not what the fitter computes; the C++
side kept its own generator because nobody taught the fitter to accept a permutation matrix. Once it does (step 1),
rows 6, 8 and 9 inherit the coupling across columns (which fixes the two defects in point 2 above), the enumeration
and the shared seed, with no change to their statistics.

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
   route CoDA, density and cluster-free DE through `drawPermutations()` of the stored model. Robust fits and both
   NA modes are preserved: per NA pattern the global P is induced on the observed rows, then used by the OLS or
   the iterative robust refit exactly as the internally generated permutation is used today. Tests: identical
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
4b. **Weighted Gower form** (serves both §8 and §10): weighted centring, H_w = X (X'WX)^- X'W, weighted dispersion
   fit, in kernel A and its R reference; then `robust = c("none", "huber", "winsor")` and
   `na.mode = c("drop", "impute_weak")` for every distance-based test, recomputed per relabeling. Tests: weights of
   one reproduce the unweighted results; a planted outlying sample is down-weighted and the null stays calibrated;
   weak imputation of an absent sample reproduces the dropped-sample estimates to numerical tolerance while keeping
   the full sample set. About a day plus the calibration runs.
5. **Retire**: `projdiff.cpp`, `pca_project`, the unused `estimateCorrelationDistance` export, `makeBlocks` /
   `permutationGroups`, R reference loops out of the namespace, `huh-jhun` either in the kernel or dropped.
6. Optional: fold `fl_fwl_cpp` into `fit_and_randomize` behind a `scheme` argument so each family has one entry point.

Steps 1-3 are each about half a day with tests; step 4 is the largest (a day) and the one with the visible payoff.
## 3a. Profile of the cluster-free loop

Toy case, 1,000 cells x 10 samples x 500 genes, neighbourhoods of 200 cells, 99 permutations (`Rprof`):

| part | time |
|---|---|
| C++ `estimateExpressionShiftsPairsLM` (all neighbourhood distances) | 0.35 s total, 0.35 ms per cell |
| whole `clusterFreeExpressionShifts()` | 31.7 s, 31.7 ms per cell |
| of which `apply(P.global, 2, inducePermutation)` per cell (split / droplevels / rank / factor / order) | ~90 % |
| `match.arg`, `pairVectorToMatrix`, `permutationPlan` for unseen subsets, `estimatePairwiseEffects` | most of the rest |
| C++ `permuted_contrast_F` per cell | not visible in the profile |

The statistic itself is negligible; inducing B sub-permutations in R per cell (B R-level calls each doing a split
and a rank) is the cost, multiplied by the number of cells, and on the simulated object (4,500 cells, 40 samples,
B = 99) it came to 8.7 minutes on 16 cores through `sccore::plapply`, with the forking overhead on top. The work
is O(cells x B x n log n) integer operations, which is a fraction of a second in C++; the whole cluster-free step
should take seconds, dominated by the neighbourhood distances.

Measured on the simulated object (2026-10-09, `cao$estimateClusterFreeExpressionShifts(n.top.genes = 500,
n.permutations = 99)`, graph neighbourhoods of median 638 cells, 40 samples): 349 s on 16 cores for 3,371 tested
cells (the rest had too few cells per sample), i.e. about 1.7 core-seconds per cell for 99 permutations of a
40 x 3 design, against 0.35 ms per cell for the C++ distances.

The batched kernel of step 4 therefore takes: the pair-layout distance matrix Y (pairs x cells, already C++), the
sample strata and `in.set` flags, the global P, the design X and contrast, and does per cell in C++: square the
pair column into a distance matrix over the samples present, Gower-centre, induce the sub-permutation by ranks
within stratum cells, compute the observed and permuted shift F (existing kernel), the shift estimate (M = A G A'
with the dispersion correction, the same quantities), and return stat / p / shift / n per cell plus the per-
permutation maxima for the max-statistic adjustment; the loop over cells runs on the sccore thread pool.

## 6. Scaling the cluster-free tests to ~10^6 cells

The batched kernel of step 4 handles this if two things are built in from the start; the current R loop cannot
(10^6 cells x 1.7 core-seconds is about 20 core-days, and it keeps every cell's B permuted statistics in memory).

**Stream per cell; never materialize a cells-wide matrix.** The pair-layout distance matrix (pairs x cells) is
780 x 10^6 doubles = 6 GB for 40 samples, and a B x cells matrix of permuted statistics is 8 GB at B = 999. The
kernel must fuse the steps per cell: collapse the neighbourhood to sample profiles, distances, Gower matrix,
induced sub-permutation, observed and permuted statistics, then keep only stat / p / z / shift / n for the cell
and update running maxima. Memory is then O(cells) for the outputs plus O(n^2 + nB) per thread.

**Accumulate the max-statistic null inside the kernel.** The adjustment needs, per permutation b, the maximum
(and minimum) z over cells. Each thread keeps its own length-B vectors of running extremes and they are reduced
at the end; results do not depend on the number of threads because P is shared and fixed.

Cost per cell (n samples, g genes, k neighbourhood cells, B permutations):

| part | flops | n = 40, g = 1000, k = 300, B = 199 |
|---|---|---|
| collapse neighbourhood to n profiles | k x nnz per cell | ~ 2 x 10^5 |
| pairwise distances (n(n-1)/2 pairs x g) | ~ n^2 g / 2 | ~ 8 x 10^5 |
| Gower centring, precompute a = A'c | n^2 | ~ 2 x 10^3 |
| permuted F for all b at once (G A_P, BLAS-3) | n^2 B | ~ 3 x 10^5 |
| induced sub-permutation | B n log n | ~ 4 x 10^4 |

About 1.5 x 10^6 flops per cell, i.e. roughly 1 ms single-threaded; 10^6 cells is 15-20 minutes on one core and
under a minute on 32 threads, dominated by the distances rather than the test. With n = 200 samples and B = 999 the
test term grows to n^2 B = 4 x 10^7 per cell (4 x 10^13 in total), still minutes on a many-core machine because it
is a dense matrix product. Freedman-Lane adds four n x n kernel parts per cell and costs about 4x the block scheme.

The same streaming applies to cell density through kernel B: `fit_and_randomize` already works per column, but it
returns `sampled_stats` (B x columns) for the max-statistic step; for 10^6 columns the kernel should instead return
the per-permutation extremes (optionally winsorized) and the per-column z, computed inside the loop.

Statistical side: with 10^6 cells the max-statistic adjustment is what controls the family-wise error, the raw
per-cell p-values floor at 1/(B+1), and the winsorized extremes (`wins`) keep a few extreme cells from dominating
the null. Cells whose neighbourhood lacks enough cells per sample are skipped per cell (induced permutation on the
present samples), as now.

## 7. What becomes unnecessary in `fit_and_randomize` once P comes from R

Dead already (no R caller passes it): graph mode. `x$pairs` is never set by `buildDesignMatrices()` any more, so
`pair_indices` / `core_pair_indices`, `PairLookup`, `subset_blocks()` and the node-shuffling branch of
`generate_permutation()` can go (`lm_common.h` 23-53, 75-127 and the `is_graph` branches in `lm_fit.cpp`).

Replaced by P: `perm_groups` parsing (`lm_fit.cpp` 280-290), `generate_permutation()`, `make_rng()`, the `seed`
argument, and on the R side `makeBlocks()`, `deriveNuisanceFactors()`, `filterNuisance()`, `permutationGroups()`,
`computeCoreRows()` and the `blocks` / `perm.groups` / `core.rows` fields of the design (`model_matrices.R`
190-202, 950-1010). `core_rows` in `fl_fwl_cpp` and `parse_core_rows()` become the identity part of P.

Duplicated: `z_from_p()` and the `z_score` output; `performLMPermutations()` recomputes z from p in R.

With robust fits dropped, the permutation loop itself collapses. For OLS the permuted contrast is linear in the
permuted response: stat_b = alpha' y[p_b] with alpha = (X'X)^-1 X' projected on c, so for all columns and all
permutations at once S = Y' A_P, where A_P (n x B) holds alpha re-indexed by each p_b. Freedman-Lane is the same
product on the residuals of the reduced model, R_red Y. Per NA pattern that is one BLAS-3 product instead of a
per-column loop with per-permutation refits; p-values and winsorized extremes are read off S in one pass (and for
10^6 columns S is streamed in column chunks). Kernel A is the quadratic analogue, F_b = a[p_b]' G a[p_b], so the two
kernels share the induced-permutation step and the A_P construction.

**Kept, by user decision (2026-10-09): robust fitting and missing-data handling are requirements, not options
to trim.** Kernel B therefore keeps two paths that both consume P:
- *OLS path* (`robust = "none"`): the matrix-product form above, per NA pattern.
- *Robust path* (`huber`, `winsor`): per column, per permutation iterative refit (IRLS / winsorized OLS), as now;
  only the source of the permutation changes.
Missing data: with `na_mode = "drop"` each NA pattern defines its own observed-row set, and the global P is
**induced** on that set with the rank trick already used for cell types (this replaces `subset_blocks()`, it
does not simply delete it); with `impute_weak` all rows take part and P applies directly; Freedman-Lane
residualizes on the reduced model per pattern before the induced P is applied. Residual outputs (`residuals`,
`residuals.pearson`, `partial_core`) are consumed by plots and are independent of the permutations; they stay but
belong in a plain fit function, not in the randomization kernel.
Open question for kernel A: the distance engine has no robust variant (robustness there comes from the
leave-one-out influence and the sensitivity report); whether a robust dispersion fit or down-weighting of
outlying samples is wanted is for the user to decide.

Rough size: `lm_fit.cpp` 749 + `lm_common.h` 207 + `projdiff.cpp` 187 lines today. What goes is graph mode, the
internal generator and the duplicated z computation (roughly 250 lines across the three files); the robust and
NA-handling code stays, the OLS path gets the matrix-product form, and the shared induced-permutation / A_P
helper (~60 lines) is added. `performLMPermutations()` (329 lines) loses the block / core-row plumbing but keeps
the robust and NA arguments.

## 8. Robust fitting in every test (user decision 2026-10-09)

The user wants `robust` options in all tests, not only in the per-column fitter. For kernel A the unit that can be
outlying is a sample (its whole distance profile), so robustness acts on samples, and the natural pair of
options mirrors kernel B:

- `robust = "huber"`: sample weights by IRLS. Each sample's residual size is its leverage-corrected residual
  distance r_i = diag(RGR)_i (the quantity the dispersion model already fits); Huber weights w_i from r_i scaled
  by a robust scale (MAD), then the model is refitted with the weighted Gower form (weighted centring, weighted hat
  matrix H_w = X (X'WX)^- X'W, weighted dispersion fit), iterated a few times. Shift, var, total, the term F and
  the dispersion F all come out of the weighted fit; the leave-one-out influence is reported alongside as now.
- `robust = "winsor"`: clip the residual distance contributions (or the off-diagonal entries of the residual
  kernel RGR) at a quantile before the statistics are formed; a one-step, non-iterative alternative.

Inference: the weights (or clipping) are **recomputed under every relabeling** inside the kernel, exactly as the
robust path of kernel B refits per permutation; weights frozen from the observed data would make the test
anti-conservative. Cost: a few IRLS iterations x O(n^2 q) per permutation, so roughly 3-5x the OLS path; acceptable
for cell-type tests, the screen and sensitivity, and still feasible per cell in the streaming kernel (robust
cluster-free tests simply take a few times longer). Calibration of both options is added to the slow simulation
suite (null rejection rates with and without planted outlying samples; power loss under clean data).

Robust options therefore appear uniformly: `estimateExpressionShiftMagnitudes(robust =)`, the whole-factor tests,
`screenCovariates(robust =)`, `checkSensitivity()` (inherits), composition term tests, cluster-free shifts, and the
existing `robust.method` of density / cluster-free DE / CoDA contrast through kernel B. Default stays `"none"`.

## 9. Residual diagnostics: what the old residuals did and where that lives now

The pair model exposed Pearson residuals and `plotExpressionShiftResiduals(cov.plot.keys = ...)` plotted them
against covariates that were not in the model, to spot relevant omitted covariates; residual magnitudes were also
used to visualize how large an effect was at the sample level. In the Gower engine the residual object is the
residual kernel G_adj = R G R (stored per cell type as `res$adjusted.distances`), and the same two uses are:

- omitted covariates: `screenCovariates(adjust.for = <model formula>)` tests every remaining covariate against the
  residual structure (location and dispersion), with `plotCovariateScreen()`; `plotSampleDistances(adjust.for =)`
  shows the residual sample configuration coloured by a candidate covariate;
- effect size at the sample level: `plotShiftDetail()` (within / between-group distance distributions and the
  adjusted sample map) and `plotSampleInfluence()` (leave-one-out change of the estimate).

For density and cluster-free DE the per-sample residuals of kernel B's observed fit stay available
(`return.residuals`, `plotDiffCellDensityResiduals()`), and a covariate diagnostic on them (residual against
candidate covariate per bin / gene, summarized) is a small addition worth making for parity with the shift engine.
Requirement for the convergence work: keep residual outputs first-class in the observed-fit function of kernel B
and keep `adjusted.distances` in the shift results; the dead cluster-free-shift residual branch in
`plotClusterFreeExpressionShifts` is removed.

## 10. Missing data in the Gower engine

Today there is one mode, sample-level complete cases per unit:

- Metadata: `buildCacoaModel()` drops samples with a missing value in any model variable (listwise deletion),
  reported in the model's issues.
- Pseudobulk: a sample with fewer than `min.cells.per.sample` cells in a cell type is absent from that cell
  type's distance matrix; the unit is fitted on the samples present (design subset, estimability check,
  `min.samp.per.level`), skipped with a reason otherwise (`res$skipped`), and its permutations are induced from the
  shared P on the present samples. The distance matrix itself never holds NA: a sample is either present with all
  its distances or absent.
- Screen: per covariate, the samples with complete covariate values are used (`n.used` reported) and the Gower
  matrix is re-centred on that subset; the shared P is induced likewise.
- Cluster-free: per cell, samples with fewer than `min.n.obs.per.samp` cells in the neighbourhood are absent,
  same mechanism.

This is the analogue of kernel B's `na_mode = "drop"`, applied at the sample level because the response of this
model is a sample's whole distance profile. **Decision (2026-10-09): `impute_weak` is added to the Gower fit.** It falls out of
the weighted Gower form planned for the robust fit (§8): keep the absent sample in the design, fill its distances with the
within-level mean distance, and give it a near-zero weight. That keeps the sample set identical across cell types
(no induced permutations, one P for everything), keeps design levels represented that a dropped sample would have
removed (estimability), and matches kernel B's option set; `na.mode = c("drop", "impute_weak")` is offered
uniformly (default `"drop"`), with calibration of the weak mode in the slow suite. Partial missingness of single distances does not arise (pseudobulk
profiles are either present or not), so there is no pairwise-deletion mode to support.

### 10a. Where missing-data awareness lives in the permutation machinery

The source itself is deliberately NA-agnostic: P is drawn once for the test's sample set (after the model's
listwise deletion of samples with missing covariates). Awareness of what is missing per unit lives in the
**induction step**, so the kernel interface is not P alone but (P, stratum id per sample, in-set flag per
sample): for a unit (cell type, NA pattern of columns in kernel B, neighbourhood, covariate subset in the screen)
the present samples are selected, and within each (stratum x in-set) cell of the present samples the labels are
reassigned by the ranks of their images under P. This is uniform on the subset's relabelings and coupled with the
full set, and it is the same code for both kernels (today `inducePermutation()` in R; in the batched kernels it
runs in C++ with the same inputs). Per unit, the number of distinct induced relabelings can be smaller than the
global one, so the per-unit `n.perm.distinct` and `p.floor` are reported in the results table rather than assumed
from the global plan.

Under `impute_weak` no induction is needed: every sample is present, P applies directly, and the near-zero weight
stays with the sample's response row while the labels move, which is exactly the exchangeability the weak mode is
meant to preserve. Under Freedman-Lane with `drop`, residualization on the reduced model is done per NA pattern on
its observed rows before the induced P is applied.
