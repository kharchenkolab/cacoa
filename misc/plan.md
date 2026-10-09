# cacoa: expression-shift redesign and workflow API — implementation plan

Status: working plan, 2026-10-08. Untracked (`misc/` is not committed; consider `^misc$` in `.Rbuildignore`).
Sources: `handoff/pairwise_shift_redesign.md` (engine; bugs B1–B20; decisions D1–D20),
`handoff/api_workflow_design.md` (workflow/API; decisions D21–D35), prototypes in `handoff/scripts/`.
Code base: branch `dev_lm` @ 99bc90a, with one uncommitted change in `R/cacoa.R`
(`plotExpressionDistance(values = "unadjusted"/"adjusted")`).

---

## 0. Facts verified on this machine (not in the handoff)

- The package does not link: `src/expression_shifts.cpp` redefines helpers that already live in
  `src/cluster_free.cpp` (`median`, `mad`, `var`, `range`, `count_values`, `collapseMatrixNorm`,
  `estimateKLDivergence`, `estimateJSDivergence`, `estimateVectorDistance`, `estimateCellExpressionShift`,
  `estimateCorrelationDistance`, `mapIds`, `applyMedianFilter` ×2, `adjustZScoresWithPermutations` ×2).
- `R/RcppExports.R` / `src/RcppExports.cpp` are stale: 14 `[[Rcpp::export]]` functions are unregistered
  (`fit_and_randomize`, `fl_fwl_cpp`, `fit_with_focusing`, `estimateExpressionShiftsPairsLM`,
  `clusterFreeGeneMat`, `mapIds`, `pca_project`, `applyMedianFilter`, `applyMedianFilterES`,
  `adjustedZScoresMaxStat`, `estimateCorrelationDistance`, `fit_density_lm`, `perm_FL_contrast_mat`,
  `perm_full_contrast_mat`); `adjustZScoresWithPermutations` is registered but no longer exported.
- C++ compiles under `CXX_STD = CXX11` (RcppArmadillo falls back to bundled Armadillo 14.6.3 with a
  pragma message). B18 is a warning here; the fix is to drop the standard pin, not raise it.
- Threading: `lm_fit.cpp`/`lm_focusing.cpp` use OpenMP (`<omp.h>` unguarded in `lm_common.h`);
  `cluster_free.cpp` uses `sccore::runTaskParallelFor` (std::thread). R side uses forked workers.
- `fit_and_randomize` seeds its RNG from `time(0)`; per-job streams `make_rng(seed, j)` already exist.
- Installed `cacoa` 0.5.0 in the user library is an older build without any lm/pair code.
- Companion design note `pairwise_model_design.md` is not on this machine.
- Toolchain: R 4.2.2 (`g++ -std=gnu++14` default), Rcpp 1.1.0, RcppArmadillo 15.2, testthat 3.3.2,
  sccore 1.0.6, conos 1.5.4. Handoff script 04 reproduces exactly.
- Fixtures: `test/shifts_sim_objects.rds` (720 MB, simulated Cacoa objects with known per-cell-type
  effect strength) — slow integration tests only. `examples/{ms,scc,alz}` — real-data notebooks.

## 1. Decisions taken (beyond the two notes)

- Gene focusing is REMOVED (user decision): `top.n.genes`, `gene.selection` (shift methods),
  `fitWithFocusingWrapper`, `fit_with_focusing`, `src/lm_focusing.cpp`. Closes B6, B10, D14.
  Gene centering (per-gene centering across the samples present for a cell type, before cosine)
  stays as an explicit, default-on, tested step.
- OpenMP is removed; remaining parallel C++ loops use `sccore::runTaskParallelFor` (already a dependency,
  already used in `cluster_free.cpp`). No new thread pool.
- C++ standard: delete `CXX_STD` from both Makevars and `C++11` from `SystemRequirements`; do not pin.
- Engine in R (D19); only the cluster-free permutation loop needs C++.
- Pair API (`pairFormula`, `pairContrast`) removed with an informative error (D11); `dist.type` aliased
  to `effect` for one release.
- Default screens do NOT switch to PCs; report `r.eff` and warn (API open question 5).
- Over-adjustment candidates: allowed with a warning plus an automatic sensitivity row (open question 1).
- D24 "user-set levels" is detected by the proxy "levels not in alphabetical order".
- Jackknife uncertainty reported as estimate ± jackknife SE (not called a CI).
- Reviewed 2026-10-08 and KEPT as is: additive plotted scale (as in `main`; ratio scale is an option), gene-centred
  cosine (dev_lm), pseudobulk transform `log10(1e3*c/S+1)` (S = row total over kept genes).
- Test grammar: vector/list forms canonical; `"var: alt vs ref"` string is sugar with documented limits.

Consolidated list of remaining deficiencies, autonomous decisions and cleanup candidates: `misc/issues.md`.

## 2. Open items needing a spec line

- "Any difference" statistic for K-level TERM tests (contrast case is shift + ½ var). Proposal: max-T over
  the location F and dispersion F under shared permutations, or drop "any" for terms in v1.
- Restricted permutation set for interaction-cell and marginal contrasts: samples whose discrete
  covariate profile matches either endpoint with the contrasted factor free (proposal).
- Minimum size rules when a cell type is observed on a subset of samples (`min.samp.per.level = 3`).

---

## 3. Tracks

### Track A — Foundation

**A.1 Build hygiene and removal of dead paths** (see §5 for the detailed outline)
1. Deduplicate C++ helpers (shared header or single translation unit).
2. Delete `src/lm_focusing.cpp`, `fitWithFocusingWrapper`, the focusing branch in
   `estimateExpressionChange`, `top.n.genes`/`gene.selection` from shift signatures and docs,
   `parseDistance` dependence on `top.n.genes`.
3. Replace the two OpenMP loops in `lm_fit.cpp` with `sccore::runTaskParallelFor`; remove `<omp.h>` and
   `$(SHLIB_OPENMP_CXXFLAGS)`; drop `CXX_STD`; clean `PKG_LIBS`.
4. `Rcpp::compileAttributes()`; `useDynLib(cacoa, .registration = TRUE)`; fix `NAMESPACE` via roxygen.
5. RNG seed as an argument to `fit_and_randomize`/`fl_fwl_cpp` (from R, D33 groundwork).
6. testthat 3 scaffold: `tests/testthat.R`, `helper-simulate.R` (individual-level simulator),
   `helper-design.R`; replace the conos-dependent dummy test; `Config/testthat/edition: 3`.
7. Commit the pending `R/cacoa.R` change; `^misc$` in `.Rbuildignore`.
Tests: package installs and loads; all C++ entry points callable; results identical for 1 vs 4 threads
with a fixed seed; `buildDesignMatrices` smoke tests (triple, simple, marginal, interaction, numeric).

**A.2 Sample-level model fixes** (B8, B11, B12, B13, B19, B20; D17/D25) -- DONE 2026-10-08.
Notes: endpoint rows are now built from RAW metadata variables (so `log(age)`, `poly(age, 2)`, `factor(batch)`
work); `buildFullDesign` keeps `predvars` in the terms object; numerics not in the contrast sit at their
sample mean in endpoint rows (contrast unchanged). `contrast_endpoints_at` (settings + weights) is stored and
`endpointRowsFromDesign()` evaluates it on any design (used for dispersion endpoints in A.3). `Cacoa$new`
stores `formula` (as used) and `numeric.ref`; all per-method rebuilds pass `numeric.ref`.
`makeToyCacoa()` (tests/testthat/helper-simulate.R) builds a Cacoa object from a synthetic count matrix
without conos; the current pipeline runs end to end on it (test-expression-shifts.R).
- Transformed covariates in endpoint rows (`.oneRowFromFormula` builds newdata from raw variables).
- Store `self$formula`, `self$dispersion.formula`; default formula `~ <test vars> + <block.vars>`.
- Add endpoint *settings* (`at` lists + marginal weights) as an attribute of the synthetic contrast,
  so the dispersion design can be evaluated at the same endpoints.
- Remove duplicate `varsFromFormula`; fix `numeric.ref` default; fix `block.vars` in cluster-free shifts.
Tests: `log(age)`, `poly(age,2)` endpoints; stored formula reused; default excludes ID-like columns;
endpoint settings round-trip for simple/marginal/interaction/numeric.

**A.3 Engine** (new `R/pairwise_effects.R`; exports `asSquaredDistance`, `estimatePairwiseEffects`,
`sampleDistanceMatrices`, `pairwiseEffectsFromDesign`) -- DONE 2026-10-08 (estimation only; inference is A.4).
Notes: `sampleDistanceMatrices()` reproduces the dev_lm per-cell-type distances exactly (verified for cor/l1/l2
and cor on PCs) with an explicit `center.genes` argument; the pseudobulk transform `log10(1e3*c/S+1)` is
untouched. `endpointRowsFromDesign()` ignores settings for variables absent from the target design.
Term-level pieces (multi-df F, R2/R2.adj/R2.partial, r.eff, dispersion F, partial R2 per term) and
leave-one-out jackknife are included. Open from the review: plotted-scale default (additive vs ratio) and
confirmation of gene-centred cosine; both unchanged for now.
- Port `handoff/scripts/proposed.R`: Gower centring, hat/leverage, `v_i`, `M̂`, bias-corrected `M̃`,
  `γ̂`, shift/var/total/ratio, normalized versions (D2), model-implied K×K cell table.
- Term-level additions (for API): multi-df location F, dispersion term test, `R2`/`R2.adj`/`R2.partial`,
  effective dimension `r.eff`, variance partition.
- Leave-one-sample-out refits (jackknife SE + influence).
- Per-cell-type subsetting; drop all-zero columns; estimability checks; skip with reason.
- Gower cache keyed by space × cell type × settings (`self$cache`).
- Centering: per gene over present samples, default on; `cor` used as-is, `l2` squared, `l1` as-is (warn).
Tests: one-factor equivalence to raw pair-cell formulas (K=2,3; unequal sizes/dispersions);
invariances (order, reference level, coding, redundant constant); metric mapping; bias correction =
classical estimator; endpoints for all contrast types; dispersion endpoints follow `at`; 1-df term F ==
contrast F; LOO refits == direct refits; estimability skips; cor == ½‖u_i−u_j‖² on the subset;
centering off changes result in the expected direction; cell table == raw cell means (saturated).

**A.4 Inference** -- DONE 2026-10-08 (`R/pairwise_inference.R`; exports `permutationPlan`, `drawPermutations`,
`testPairwiseEffects`). Scheme names: `"auto"`, `"block"` (existing block randomization, extended to swap only
the compared conditions; exact for discrete nuisance), `"freedman-lane"` (fallback), `"huh-jhun"` (opt-in,
shift only). Both BH (`padj.*`, across cell types) and max-statistic FWER (`p.fwer.*`) are reported, plus a
global max-T p per effect; BH remains the default for marking significance in plots (user decision pending
confirmation). Permuted statistics are O(n^2 q) per permutation (X'X invariant under row permutation).
- Restricted permutation: sample-level `makeBlocks`/`permutationGroups`, contrasted-level samples only
  (also without blocks); exhaustive enumeration when count < `n.permutations`; p-value floor.
- Freedman–Lane in Gower form (fallback; warn on large `r.eff`); Huh–Jhun opt-in.
- Coupled permutations across cell types with different sample subsets (restriction of a global π).
- Max-T: global p per test×effect; familywise-adjusted per-cell-type p.
- Shift one-sided (pivotal F), var/total two-sided; seed honoured everywhere (incl. remaining C++).
Tests: generator moves only contrasted samples, respects strata, enumeration count correct; induced
subset permutation uniform within strata (exact enumeration small n); permuted F == brute-force on
`D[p,p]`; Gower-FL identity vs feature-level residual permutation; max-T adjusted p ≥ raw p; fast null
sims (n 12–18, 99 perms, 200 reps) within binomial limits; var test rejects under dispersion-only.

**A.5 Validation study** -- DONE 2026-10-08 (`misc/validation/validate_trackA.R`, results in
`misc/validation/results/REPORT.md`; compact committed subset in `tests/testthat/test-validation-slow.R`).
42 scenarios x 200 reps x 199 permutations (binomial 95% band for 0.05: [0.020, 0.080]) plus an 8-cell-type
max-T calibration. Findings:
- Estimation (l2): shift and var unbiased in every scenario, including confounded batch with effect-vector
  cosine +0.6 / -0.6 (new 7.8 / 7.8 vs truth 7.2; the old pair-LM gives 11.9 / 3.4), 3 groups with a shifted
  or over-dispersed third group, unequal sizes, marginal and interaction-cell contrasts, age slope.
- Fix found by the study: with a continuous covariate in the location model and unequal dispersions, the
  per-sample leverage-corrected dispersion was biased towards the other group (var 135 vs truth 160). The
  dispersion model is now fitted through E[r] = (R o R) Z gamma (exact for any dispersion pattern; identical to
  the classical estimator in the one-way layout), and the same s enters the shift bias correction. After the
  fix: 160.2.
- Type I error (block scheme): shift/var/total within band in all null designs, including confounded batch
  (0.061/0.045/0.061 and 0.035/0.056/0.040) and the 3-group designs. The old pair-LM with the as-implemented
  distance response leaks dispersion into shift (0.27 vs 0.11 under shift+dispersion).
- Freedman-Lane: shift conservative (0.02-0.05), as the handoff found; var/total mildly liberal at n = 16 with
  batch+age (0.065-0.075 at 800 reps), within band at n = 20 and for age alone. Documented; FL rows are flagged
  "approximate-scheme".
- Huh-Jhun: calibrated here (0.065 null, 0.050 under 2x dispersion without nuisance).
- Metrics: cor behaves like l2 for shift; under a dispersion-only change the cor shift test rejects 0.08
  (unit-sphere geometry, R5) and l1 rejects 0.215 (l1 does not separate location from dispersion; keep the
  warning, D3). Unequal sizes + dispersion: shift test conservative (0.000), total picks it up (0.295); the
  "unbalanced+dispersion" flag covers this.
- Power at n = 16 is modest for all methods at a 0.3/gene shift (0.13-0.20); block beats the old pair-LM FL
  under batch (0.18 vs 0.085).
- Across cell types: global max-T 0.055/0.035/0.040 for shared variance 0/0.5/0.9 (calibrated); BH "any
  significant" 0.045/0.020/0.010; raw "any" 0.30/0.22/0.15.

### Track B — API MVP (api note §5.1 "MVP")
Clarifications (user review 2026-10-08):
- Options: every method keeps its arguments; unspecified arguments fall back to `self$options`, explicit call
  values win (explicit > setOptions > package default); the resolved values go into result provenance.
- Model: the constructor sets the model (normal path, unchanged); per-method `formula`/`contrast` still build a
  temporary model with a notice; `setModel()` is the explicit way to change the stored model after exploration.
- Order of work: B.1 options -> B.3 model object (minimal) -> B.5 shifts on the engine -> B.6 plots/accessors ->
  B.7 legacy aliases -> B.2 metadata audit -> B.4 design check. Checkpoints after B.5 and B.7.
- "Any difference" is defined for contrasts only (shift + var/2); K-level term tests report location F and
  dispersion F with a pairwise table.
- `setOptions`/`getOptions` (n.permutations 999, seed 1, alpha, p.adjust, n.cores, verbose, provenance).
- `describeMetadata` (roles: usable / id-like / constant / high-cardinality / mostly-missing; notes).
- `setModel`/`getModel`, `cacoaModel` S3 with print/summary; test grammar (§4.4); reference heuristics
  (D24 + alphabetical proxy); interpretation sentences; permutation plan per test.
- `checkDesign`/`plotDesign`: associations (Cramér's V, η, Spearman), balance, nesting/aliasing,
  singleton cells, GVIF, df budget, permutation feasibility, over-adjustment candidates.
- `estimateExpressionShiftMagnitudes(test=, effects=, dist=, permutation=, ...)` on the model object;
  results table of api §4.7 (estimate, estimate.norm, se.jk, statistic, p, padj, p.fwer, n, scheme,
  n.perm.distinct, p.floor, flags); `res$global`, `res$fits`, `res$adjusted.distances`, `res$influence`.
- `plotExpressionShiftMagnitudes`: dot + jackknife bars, three panels, provenance subtitle (D34).
- `getSampleGroups(test)` wiring; deprecated aliases `sample.groups`/`ref.level`/`target.level` → `test`.
- `Cacoa$new` without model; metadata audit printout. Vignette rewrite.
Tests: grammar rows; reference heuristic table; metadata roles; each checkDesign issue class on a
synthetic design; provenance present; old vignette constructor runs via aliases; ID-column trap errors;
end-to-end on synthetic `cm.per.type` without a data object; slow test on `shifts_sim_objects.rds`
(skipped if absent): "strong" cell types rank above "none".

Status (2026-10-08, autonomous session):
- B.1 done: `cacoaDefaultOptions()`, `cao$setOptions()/getOptions()`; every method resolves unset arguments through
  `private$opt()` (explicit > setOptions > default); `seed` option (1) makes runs reproducible.
- B.3 done (minimal): `R/model_api.R` — `resolveTests()` (full §4.4 grammar: variable, K-level term, numeric step,
  several variables, "all", "var: alt vs ref", triple, structured list, coefficient weights), `chooseReferenceLevel()`
  (D24), `buildCacoaModel()` -> `cacoaModel` (tests with per-test designs + permutation plan, listwise deletion,
  issues: aliasing error, df warning, small levels, permutation floor), print/format/summary; `cao$setModel()/getModel()`;
  the primary test's design fields sit at the top level of the model for backward compatibility
  (`cao$model$F`, `$contrast.F`, `$contrast_spec` keep working for CoDA / density / DE / cluster-free).
- B.5 done: `R/expression_shifts_api.R` — `pseudobulkPerCellType()` (min.cells.per.sample now honoured), term tests
  (`termEffectsFromDesign()`, `testTermEffects()`: location F / R2.adj, dispersion F / R2.disp.adj, pairwise level
  table, block or FL permutations shared across cell types), `expressionShiftsForModel()` (long results table of
  §4.7, global max-T p per test x effect, adjusted distances `R_Z G R_Z`, LOO influence matrices);
  `cao$estimateExpressionShiftMagnitudes(test=, formula=, contrast=, dist=, permutation=, ...)` with pseudobulk and
  distance caches (`cao$cache$pseudobulk`, `$sample.distances`) and legacy-argument translation (`perm.method` ->
  `permutation`; `dist.type`, pair args, robust/na.mode, top.n.genes announced and ignored; unknown args error).
- B.6 partly: `plotExpressionShiftMagnitudes()` dot/bar + jackknife CI + provenance subtitle (BH default marking,
  `significance = "p.fwer"` / `"p"`), `plotSampleInfluence()` (replaces `plotExpressionShiftResiduals`, kept as
  deprecated alias), `getSampleDistanceMatrix(space = "expression.shifts")` reads `res$distances` /
  `res$adjusted.distances`; `plotSampleDistances` / `estimateMetadataSeparation` work on the new results.
- B.7 partly: constructor without model prints the metadata audit; `sample.groups`/`ref.level`/`target.level`
  accepted as deprecated aliases (-> `condition` column + test); `self$ref.level/target.level/sample.groups/contrast/
  formula` synced from the primary test; `cao$getSampleGroups(test)`; palette named by target/ref.
- Tests: `test-model-api.R`, rewritten `test-expression-shifts.R`, `test-example-datasets.R` (SCC paired design from
  `examples/scc/x.rds`: unpaired vs patient-blocked model, exhaustive 2^8 enumeration; conos panel smoke through the
  new entry point; slow planted-strength ordering on `test/shifts_sim_objects.rds`).
- (Resolved later in the session: B.4 `checkDesign()/plotDesign()` written in `R/design_check.R`; other analyses go
  through `private$resolveModel()` since Track D.)

### Track C — Exploration (api §5.1 phase 2)
- `screenCovariates` (marginal/partial, `adjust.for`, df-budget fallback, D27 defaults, suggested
  formula, BH over grid, per-covariate max-T); `plotCovariateScreen`, `plotCovariateSummary`,
  `plotVariancePartition`; `plotSampleDistances(adjust.for=)` on `G_adj = R G R`; `plotShiftDetail`;
  `estimateMetadataSeparation`/`plotMetadataSeparation` as aliases.
Tests: E2 pattern (marginal flags medication, partial does not); dispersion separates from location (E1);
`R2.adj` ≈ 0 under null; ring mark for marginal-only; adjusted distances PSD; alias returns screen.

Status (2026-10-08): Track C done. `R/covariate_screen.R`: `screenCovariates()` (marginal / partial term tests per
covariate x cell type on Gower matrices; numeric covariates per SD; FL permutations shared across cell types; BH over
the grid; per-covariate max-T global p; df-budget fallback to marginal; over-adjustment annotation; suggested formula),
`variancePartition()` (chance-corrected unique / shared / residual), `plotCovariateScreen()` (dot = FDR, ring =
marginal-only, global column), `plotCovariateSummary()`, `plotVariancePartition()`; `R/plot_detail.R`:
`plotShiftDetailPanels()`. R6: `screenCovariates(space = expression.shifts | composition)`, `plotCovariateScreen`,
`plotCovariateSummary`, `plotVariancePartition`, `plotShiftDetail`, `plotSampleDistances(adjust.for=, color.by="test")`,
`getSampleDistanceMatrix(adjust.for=)` (G_adj = R G R); `estimateMetadataSeparation()` is now the marginal screen on the
joint distance (returns `pvalues`, `pseudo.r2` = R2.adj, `padjust`, `screen`), `plotMetadataSeparation()` unchanged.
Finding: the analytic effective-dimension p-values are several-fold liberal at r.eff ~ 18 (p 0.012 vs permutation
0.03), as the handoff warned; they stay a labelled preview only.

### Track D — Robustness and other analyses (api §5.1 phase 3)
- `checkSensitivity`/`plotSensitivity` (D31 covariate sets, LOO influence, verdicts);
  `plotSampleInfluence` (replaces `plotExpressionShiftResiduals`).
- DE term tests (LRT / limma F) and CoDA term tests on the same model; density weights `X(XᵀX)⁻c`.
- Cluster-free shifts on the new engine with shared permutations and max-stat (C++ loop); fix B13/B14.
- C++ port of the per-permutation step of the distance engine (`permutedStats`: shift F, var, total under a
  relabeling; Freedman-Lane reindexing) behind the same R interface; the R version stays as the reference
  that tests compare against. Plan, permutation drawing, estimability, endpoints and result assembly stay in R.
- The existing C++ linear-model fitter (`fit_and_randomize`) accepts the permutation matrix drawn in R
  (`drawPermutations`) instead of generating its own from block lists, so CoDA, density, cluster-free DE and
  the shift tests use identical permutations and block rules on the same object.
Tests: sensitivity verdict flips on a planted covariate; influence flags a planted outlier; density
weights agree with old for balanced two-group and sum to zero within strata; neighbourhood == cell type
reproduces the cell-type result; max-stat receives non-empty vectors.

Status (2026-10-08): Track D done except the R-drawn-permutations hook for `fit_and_randomize` (CoDA / density /
cluster-free DE still use the C++ fitter's own block randomization; the pairwise engine and cluster-free shifts share
R-drawn permutations).
- `R/sensitivity.R`: `sensitivityFormulas()` (D31 sets), `checkSensitivity()` -> `cacoaSensitivity` (model x cell type x
  effect table; per-cell-type verdict: robust / sign or significance depends on adjustment / estimate changes by x% /
  driven by sample S, from leave-one-out changes that are > 3 MAD outliers and > 20% of the estimate), print,
  `plotSensitivity()` forest plot; R6 `checkSensitivity()`, `plotSensitivity()`.
- `regressionWeights()`: density sample weights are now X (X'X)^- c (sum to zero within strata; the old F c broke with
  an intercept); `estimateCellDensity/DiffCellDensity/DEPerCellType/CellLoadings/ClusterFreeExpressionShifts` go
  through `private$resolveModel()` (stored model or announced temporary one) and store the model (`res$model`, or
  `attr(de, "model")` for DE). Whole-factor tests: DESeq2 Wald -> LRT, edgeR / limma F over the term's columns,
  CoDA -> `testTermEffects()` on ILR distances, density / cluster-free refuse with a clear message.
- `src/perm_stats.cpp`: `permuted_contrast_F()` / `permuted_contrast_F_fl()` (shift F under relabelings; block and FL)
  = the R `permutedStats()` to 1e-10; used by the cell-type test when var/total are not needed and by the cluster-free
  port. `R/cluster_free_shifts.R`: `clusterFreeExpressionShifts()` (distances per neighbourhood from the existing C++
  pair kernel, shared permutations induced per neighbourhood, max-stat z adjustment, median smoothing; fixes the
  `z.score`/`z.scores` field mismatch used by `plotClusterFreeExpressionShifts`).
- Tests: `test-track-d.R` (kernel vs R, weights, sensitivity verdicts incl. planted outlier, DE/CoDA term paths,
  cluster-free whole-type neighbourhoods reproduce the cell-type F), conos-panel cluster-free smoke.

### Track E — Retirement, docs, slow simulations
- Delete `buildPairDesignMatrices`, `pairifyMeta`, pair-level FL/blocks, `extractPairwiseShifts`,
  `extractContrastDeltas`, `extractFits*`, `estimateR2PerTerm` on pairs, `inferCenteringType`;
  informative errors for `pairFormula`/`pairContrast`.
- Roxygen, vignette with §1.2 interpretation table, CHANGELOG, `R CMD check` clean.
- `skip_on_cran` simulations: redesign §4.4 rows 1–7, §4.5 null rates, api E1 and E3.

Status (2026-10-08): Track E done except the `fit_and_randomize` permutation hook (see Track D note).
- Removed: `buildPairDesignMatrices`, `pairifyMeta` and all pair helpers, `estimateExpressionChange`, the pair-layout
  distance functions, `extractPairwiseShifts` / `extractContrastDeltas` / `extractFits*`, `estimateR2PerTerm`,
  `inferCenteringType`, legacy shift / residual plots, `estimateClusterFreeExpressionShiftsLM`, unused metadata
  helpers (`validateDesign`, `subsetMetadata`, `buildContrastMatrix`, internal `getSampleGroups`).
  `pairFormula` / `pairContrast` raise an informative error.
- Docs: vignette rewritten around the workflow (chunks not evaluated / not purled: the data is not shipped);
  CHANGELOG entry; `R/data.R` for `panel.preprocessed`; R6 method docs repaired (one @param per argument, no
  stale params); `R CMD check` fixes (non-ASCII in code, undefined globals, imports, Rd links, VignetteBuilder).
- Slow tests: E1 (partial screen calibration with 7 covariates; dispersion separates) and E3 (max-T calibration
  across cell types) added to `test-validation-slow.R`; all slow tests pass (`NOT_CRAN=true`).
- Build lesson: R's Makefile does not track header dependencies; a stale `lm_fit.o` after a header edit produced
  memory corruption that surfaced as a segfault when *another* package's DLL was loaded later (conos tests).
  `scratchpad/build_test.sh` now deletes `src/*.o` before every install.
- Final `R CMD check --no-build-vignettes --no-tests` (commit 54ec7d0): 0 errors; 2 warnings, both from the skipped
  vignette build (no inst/doc); 3 notes (Suggests not installed here, installed size 18.5 MB, fabia xref).
- End-to-end workflow test on the simulated conos object (`test-example-datasets.R`, slow, 7.8 min) passes:
  audit -> checkDesign -> screen (Group found, unrelated age not) -> setModel (contrast + whole-factor test) ->
  shifts (planted strong types significant) -> sensitivity (robust) -> DE / CoDA / density / cluster-free.

### Later
Analytic preview p-values (after R1), metric sensitivity, Welch-type shift test (D6/R2), robust weights
(D13), repeated measures (R3), technical-noise dispersion covariate (D16/R4), cor geometry check (R5).

---

## 4. Order and dependencies

A.1 → A.2 → A.3 → A.4 → B → C → D → E. C depends on A.3 (term tests, cache) and B (model object for
annotations) but could start after A.3 in parallel with B. D's cluster-free port depends on A.4.

---

## 5. A.1 detailed outline

Status (2026-10-08): DONE. Package installs (`R CMD INSTALL`), 41 unit expectations pass, handoff script 04
reproduces identically against the modified `lm_fit.cpp`. Deviations from the outline below:
- sccore's `sccore_par.hpp` defines non-inline functions, so it can be included in only one translation unit:
  added `src/parallel.h` / `src/parallel.cpp` (`cacoa::parallelFor`, serial when `n_cores <= 1`) and routed
  both `lm_fit.cpp` and `cluster_free.cpp` through it.
- roxygen regeneration needed `options(keep.source = TRUE)` with `load_source` (R6 class) and the removal of
  `@inheritDotParams fabia::fabia` (forces loading the optional fabia package). 1296 pre-existing roxygen
  warnings in `cacoa.R` (undocumented R6 method arguments) remain; the old NAMESPACE exported a function that no
  longer exists (`estimateDEPerCellTypeInner`).
- B6 fixed directly (guard `sampled_stats` read when `n_randomizations == 0`) since `fl_fwl_cpp` stays for
  CoDA/density/cluster-free until Track A.4.
- `handoff/scripts/common.R` updated: appends `parallel.cpp` to the standalone copy of `lm_fit.cpp`; the cpp17
  patch and the macOS omp stub are gone.

### 5.1 C++ deduplication
- Create `src/cf_common.h` (or similar) holding the shared inline/static helpers now defined in both
  `cluster_free.cpp` and `expression_shifts.cpp`: `median`, `mad`, `var`, `range`, `count_values`,
  `collapseMatrixNorm`, `estimateKLDivergence`, `estimateJSDivergence`, `estimateVectorDistance`,
  `estimateCellExpressionShift`, `estimateCorrelationDistance`, `applyMedianFilter` (both overloads),
  `adjustZScoresWithPermutations` (both overloads). Keep ONE definition each (`inline` for header
  definitions), include the header from both .cpp files, delete the copies.
- `mapIds` and `estimateCorrelationDistance` are `[[Rcpp::export]]` in both files: keep one export.
- Verified: the two copies are semantically identical rewrites (differences are namespace
  qualification, `2.0` vs `2`, explicit casts, `Rcpp::stop` formatting). Keep the `expression_shifts.cpp`
  style (explicit namespaces, no `using namespace` leakage) in the shared header.
- `src/projdiff.cpp` exports (`fit_density_lm`, `perm_FL_contrast_mat`, `perm_full_contrast_mat`) have no
  R callers: register them anyway (cheap) and mark as candidates for removal in Track E.
- `src/cluster_free.cpp` includes `<unistd.h>` without using it: remove (Windows portability).

### 5.2 Remove focusing
- Delete `src/lm_focusing.cpp`.
- `R/expression_distance.R`: delete `fitWithFocusingWrapper` (L190–364) and the `use.focusing` branch
  in `estimateExpressionChange`; drop `top.n.genes`, `gene.selection` params and roxygen.
- `R/cacoa.R`: drop the two params from `estimateExpressionShiftMagnitudes` and its docs (L301–330,
  L375, L411). Leave `top.n.genes` in DE-stability methods and `n.top.genes`/`gene.selection` in
  cluster-free methods (different meaning, unrelated to focusing).
- `parseDistance(dist, n.pcs)`: remove `top.n.genes`; update the one caller at `cacoa.R:3025`.
- `pca_project` (expression_shifts.cpp) is unused from R: delete or keep as internal; decide at merge.

### 5.3 OpenMP → sccore
- `src/lm_common.h`: remove `#include <omp.h>`; add `#include <sccore_par.hpp>`.
- `src/lm_fit.cpp`: loop at ~L313 (design pre-computation over groups) and ~L368 (jobs) become
  `sccore::runTaskParallelFor(0, n, [&](int i){...}, n_cores, false)`. Bodies are already independent
  per index (each job writes its own column; `make_rng(seed, j)` per job).
- `src/Makevars`, `src/Makevars.win`: remove `$(SHLIB_OPENMP_CXXFLAGS)` (both lines), remove
  `CXX_STD = CXX11`; simplify `PKG_LIBS` to `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` (Rcpp::LdFlags() is
  obsolete; `-lpthread` comes from R's defaults; `-L/usr/lib/` is wrong on most systems).
- `DESCRIPTION`: remove `SystemRequirements: C++11`.
- Also check `src/cluster_free.cpp` includes `<unistd.h>` (not portable to Windows) — guard or remove.

### 5.4 Registration
- `Rcpp::compileAttributes()`; `roxygen2::roxygenise()`; `useDynLib(cacoa, .registration = TRUE)`
  and `importFrom(Rcpp, evalCpp)` in the package-level roxygen block.
- Remove the stale `adjustZScoresWithPermutations` wrapper if it is no longer exported.

### 5.5 Seed
- Add `seed` (uint64) parameter to `fit_and_randomize` and `fl_fwl_cpp`; replace `time(0)` seeding;
  R wrapper `performLMPermutations(seed = NULL)` draws the seed from R's RNG when NULL (so `set.seed`
  works). Also used by CoDA and density paths.

### 5.6 Tests scaffold
- `tests/testthat.R` (`test_check("cacoa")`), `DESCRIPTION: Config/testthat/edition: 3`.
- `tests/testthat/helper-simulate.R`: `simulateIndividualModel(n, p, groups, batch, age, effect angles,
  dispersion per group)` → list(meta, Y, D2); `helper-design.R`: small metadata designs.
- `test-build.R`: namespace loads; each registered C++ symbol exists; `fit_and_randomize` 1 vs 4 threads
  identical with fixed seed; `fl_fwl_cpp` with `n_randomizations = 0` and non-empty Z does not crash (B6
  path removed with focusing, but the function should still be safe).
- `test-design-matrices.R`: smoke tests for contrast types.
- Replace `test_functions.R` (conos-dependent) — move it to a `skip_if_not_installed("conos")` block.

### 5.7 Verification
- `R CMD INSTALL` into a scratch library; `devtools::test()`; `R CMD check --no-manual` (expect notes
  only for pre-existing issues); handoff `scripts/run_all.sh` still runs (uses `lm_fit.cpp` standalone —
  update `common.R` to drop the cpp17 patch and the omp stub once OpenMP is gone).
