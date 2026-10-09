# Known issues, deficiencies and open decisions (dev_lm, as of 2026-10-09, HEAD 54ec7d0)

Consolidated from the per-track status notes in `misc/plan.md`. Items are grouped by what the user has to decide or
what still has to be built; nothing here is tracked by a test failure (the fast suite, the slow suite and
`R CMD check` are clean).

## 0a. Found and fixed during the engine convergence (2026-10-09)

- `impute_weak` in the C++ fitter detached the weak weights from their rows under relabeling (step 1).
- Max-T p-values used random global draws while per-cell-type p-values were enumerated: `p.fwer` could be
  below the raw p (step 2).
- Numeric tests (`test = "age"`) permuted nothing: the plan marked no sample swappable (step 2). All earlier
  numeric-test p-values from the block scheme were therefore uninformative (p = 1 or NA).

## 0. Found and fixed by executing the walkthrough (2026-10-09)

Running the whole workflow as a notebook on the simulated object exposed four defects the unit tests had missed:
`plotExpressionShiftMagnitudes(test = "Group")` did not accept a variable name; `plotCellLoadings()` still expected
the old bootstrap result structure (new `plotCodaLoadings()`); `plotVolcano()` looked for `CellFrac` in the wrong
list level; `estimateDiffCellDensity()` (default Freedman-Lane) crashed in `fl_fwl_cpp` when residuals were not
requested and some bins had missing values. Regression tests are in `test-track-d.R`. Lesson: every plotting
method needs at least a smoke test on the result structure its estimator now produces.

The R / C++ split of fitting and randomization, and the plan to converge on two kernels, is in
`misc/engine_convergence.md`.

## 1. Not implemented

- **Robust fitting in the distance engine** (user requirement 2026-10-09): `robust = "huber"` (sample-weight IRLS
  on leverage-corrected residual distances) and `"winsor"` for shift / var / total, term tests, screen, sensitivity,
  composition term tests and cluster-free shifts, recomputed under each relabeling; and `na.mode = "impute_weak"`
  (absent sample kept with a near-zero weight) for the same tests (user decision 2026-10-09); both through the weighted
  Gower form, see `misc/engine_convergence.md` §8, §10 and step 4b. Not started.
- **Residual-vs-covariate diagnostic for density / cluster-free DE** (parity with `screenCovariates(adjust.for=)`
  on the shift side); see §9 of the same note. Not started.

- **R-drawn permutations for the C++ fitter.** `fit_and_randomize` (used by CoDA contrast tests, cell density
  and cluster-free DE) still draws its own block permutations. Only the expression-shift engine and the
  cluster-free shift port share R-drawn permutations (`drawPermutations()`). Consequence: the permutations of
  CoDA / density / cluster-free DE are not coupled with the shift tests, and the Freedman-Lane / Huh-Jhun
  schemes are not available there (block only; FL only via the old `fl_fwl_cpp` path).
- **Cluster-free shifts and density refuse whole-factor (K-level) tests** with a message; only contrasts and
  numeric steps run there.
- **Analytic p-values** (`p.values = "analytic"` in the screen) remain a labelled preview; see §3.
- **"Later" list from the plan** untouched: Welch-type shift test, robust weights, repeated measures, technical
  noise dispersion covariate, cor-geometry check, metric sensitivity report.

## 2. Decisions made autonomously that need the user's confirmation

- **BH is the default significance marking in plots** (`significance = "padj"`); max-T (`"p.fwer"`) and raw `"p"`
  are options. Marked "pending confirmation" since Track A.4.
- **Numeric test default step is 1 unit** (`test = "age"` tests a slope per unit; the screen uses per SD). Different
  conventions in the two places; chosen so the shift estimate matches the covariate's own units.
- **Reference-level heuristic order:** explicit factor order first, then control-like names
  (control/ctrl/healthy/normal/wt/ref/baseline/untreated/vehicle/mock), then the most frequent level. The order
  of the first two could be argued either way.
- **ID-like metadata rule:** a column is "id-like" if every level is a singleton, or if its name looks like an
  id and it has more than 2n/3 levels. `patient` with 10 levels over 18 samples is therefore usable (as a block).
- **Design-check thresholds:** aliasing is an error, `n - q < 1` an error; association warnings at
  bias-corrected Cramér's V / eta / |Spearman| above 0.3 (`assoc.flag`), over-adjustment candidates above 0.5
  (`overadjust.flag`); GVIF above 2.5 is flagged. Not calibrated beyond
  the toy and SCC designs.
- **Sensitivity verdicts:** "driven by sample S" when a leave-one-out change is both a > 3 MAD outlier and > 20 %
  of the estimate (`influence.rel = 0.2`); "estimate changes by x %" when the model-set spread exceeds 50 %
  (`rel.change = 0.5`); "depends on adjustment" when the sign or BH significance differs across models. These
  thresholds are guesses. The near-collinear proxy covariate in the toy test makes every cell type "depends on
  adjustment", which is intended but the wording may read as a defect of the method rather than of the covariate.
- **Over-adjustment is exposed, not prevented:** covariates strongly tied to the test variable are allowed with a
  warning and an annotation in the screen; no automatic exclusion.
- **`estimateMetadataSeparation()` changed semantics:** it is now the marginal covariate screen on the joint
  expression distance (returns `pvalues`, `pseudo.r2` = R2.adj, `padjust`, `screen`); the old pseudo-R2 was a
  different quantity. `plotMetadataSeparation()` still works on it.
- **Legacy arguments** `dist.type`, `robust.method`, `na.mode`, `top.n.genes`, `gene.selection` are accepted and
  announced as ignored; `pairFormula` / `pairContrast` are errors; unknown arguments are errors.
- **Default screen does not switch to PCs** on high effective dimension; it reports `r.eff` and warns.
- **Plotted scale stays additive** (not ratio) and the cosine distance stays gene-centred (user instruction; both
  confirmed only in that they are unchanged).

## 3. Statistical findings to keep in mind

- **Uncoupled permutations across response columns in the C++ fitter** (`fit_and_randomize`: one RNG stream per
  column, `make_rng(seed, j)`). `lmCoda()` back-transforms row b of the permuted ILR coefficients as if the K
  coordinates had been relabelled together, so the per-cell-type loading null and the global composition statistic
  ignore the correlation between coordinates; `estimateDiffCellDensity()`'s max-statistic adjustment over bins has
  the same problem. Pre-existing (dev_lm before this work); fixed by feeding the fitter one R-drawn permutation
  matrix (step 1 of `misc/engine_convergence.md`).

- Analytic effective-dimension p-values are several-fold liberal at r.eff ≈ 18 (0.012 analytic vs 0.03 by
  permutation in the SCC check); labelled preview only, documented in `screenCovariates()`.
- Freedman-Lane: shift test conservative (0.02–0.05); var / total mildly liberal at n = 16 with batch + age
  (0.065–0.075 over 800 reps), within band at n = 20. FL rows carry the "approximate-scheme" flag.
- With `cor` distance a pure dispersion change rejects the shift test at ~0.08; with `l1` at ~0.2 (l1 does not
  separate location from dispersion; warning kept). Unbalanced groups + dispersion change: shift test
  conservative, `total` picks it up ("unbalanced+dispersion" flag).
- Power at n = 16 for a 0.3/gene shift is modest for every method (0.13–0.20).
- Global max-T and BH "any significant" are calibrated across cell types; raw "any" is not (expected).

- Notebook observation to review: in the cluster-free shifts on the simulated object (500 top genes, 99 permutations)
  `IN-PV` shows almost no per-cell shift although its cell-type test is strong and `L2/3` / `IN-SST` light up as
  expected. Possibly the planted IN-PV genes fall outside the 500 most expressed genes; not investigated.

## 4. Housekeeping and cleanup candidates

- `misc/` (plan, this file, `validation/` driver and report) is untracked; decide whether to commit it. It is in
  `.Rbuildignore`.
- `src/projdiff.cpp` exports (`fit_density_lm`, `perm_FL_contrast_mat`, `perm_full_contrast_mat`) and
  `pca_project` / `estimateCorrelationDistance` in `expression_shifts.cpp` have no R callers; candidates for
  removal. `lm_fit.cpp` comments still mention OpenMP threads (code uses the sccore pool).
- Pre-existing TODOs in `R/cacoa.R` (overall p-adjustment in ontology plots, z-score adjustment in
  cluster-free DE, binary distance in `plotOntologySimilarities`), `R/ontology.R`, `R/cell_density.R`
  (`findScoreGroupsGraph` deprecated?), `R/de_function.R` (bootstrap resampling use).
- Vignette: chunks are neither evaluated nor purled because the PF Conos object is not shipped, so
  `R CMD check` with a vignette build will warn (no `inst/doc`). A small shipped example (e.g. a trimmed
  `panel.preprocessed` run) would let the vignette build.
- `R CMD check` notes: installed size 18.5 MB (mostly `data/`), Suggests not installed on this machine
  (uwot, knitr, rmarkdown), a `fabia` cross-reference.
- Tests: the conos smoke test in `test-expression-shifts.R` is `skip_on_cran` and therefore skipped under plain
  `Rscript`, while the conos panel test in `test-example-datasets.R` runs; one of the two is redundant.
  Example-dataset tests depend on files outside the package (`../examples/scc/x.rds`, `../test/shifts_sim_objects.rds`)
  and skip when absent.
- The old legacy fields `self$ref.level / target.level / sample.groups / contrast / formula` are synced from the
  primary test for backward compatibility; they could be dropped once nothing external reads them.

## 5. Build hygiene (recurring hazard)

- R's Makefile does not track header dependencies. After editing any `src/*.h`, delete `src/*.o` before
  `R CMD INSTALL`; a stale `lm_fit.o` produced memory corruption that surfaced only as a segfault when conos
  was loaded later. `scratchpad/build_test.sh` does this; a `make clean` step or a `src/Makevars` dependency
  rule would make it permanent.
- Do not run an install into the scratch library while a test process is using it.
- Roxygen must be run with `roxygen2::roxygenise(load_code = roxygen2::load_pkgload)`; `load_source` loses the
  R6 method docs.
