# Known issues, deficiencies and open decisions (dev_lm, as of 2026-10-10)

Consolidated from the per-track status notes in `misc/plan.md`. Items are grouped by what the user has to decide or
what still has to be built; nothing here is tracked by a test failure (the fast suite, the slow suite and
`R CMD check` are clean).

## 0c. Model builder rewrite (2026-10-10, misc/model_builder_plan.md)

- Logical test variables failed in the design builder; a constant covariate was pruned silently; the plain grammar
  errored on interaction models; the per-column fitter returned NaN for columns whose observed samples lacked a
  nuisance level although the contrast was estimable; Freedman-Lane in the fitter residualized on a column split
  that omitted the non-tested direction of a multi-level factor. All fixed by the rewrite (level-means coding,
  direction split, estimability by pseudo-inverse, formula notes, interaction defaults).
- Behaviour changes to know: a factor with more than two levels must be tested with an explicit comparison
  (`"Group: G2 vs G1"`) or `"Group: all"`; coefficient names of the per-column fitter and of CoDA are level
  means (`groupA, groupB`) instead of `(Intercept), groupB`; `lincomb` contrasts are gone; `cao$formula` /
  `cao$contrast` are gone (use `cao$model`).

## 0b. Found and fixed while reviewing the open decisions (2026-10-10)

- **Freedman-Lane in the per-column fitter fitted the target-level samples only** for a two-level test with a
  nuisance covariate (the reference level sits in the intercept, which belongs to `Z`, so `core.rows` excluded it):
  null rejection 0.027 and power 0.077 against 0.053 / 0.577 with all rows (12 samples, numeric covariate). The
  auto scheme picks Freedman-Lane exactly when a numeric nuisance is present, so composition, density and
  cluster-free DE tests of a two-level condition with a numeric covariate were affected. The fitter now fits all
  rows (`fit_and_randomize(Z = )`; `fl_fwl_cpp` and `core.rows` removed) and matches `lm`.
- **Huh-Jhun removed**: it used its own RNG (uncoupled across cell types, max-T invalid), was shift-only and
  silently fell back to block under weights.
- **Raw-coefficient statistic was not pivotal** under block relabeling with a numeric covariate (p 0.12 against
  0.009 for t on the same permutations); the fitter now studentizes by default.

## 0a. Found and fixed during the engine convergence (2026-10-09)

- `impute_weak` in the C++ fitter detached the weak weights from their rows under relabeling (step 1).
- Max-T p-values used random global draws while per-cell-type p-values were enumerated: `p.fwer` could be
  below the raw p (step 2).
- The cluster-free pair-distance builder shifted 0-based neighbourhoods lacking cell 0 by one cell (step 4);
  the walkthrough's cluster-free panel was computed on mis-indexed neighbourhoods and is re-rendered.
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

## 0d. Real-dataset pass (2026-10-10): three example notebooks, bugs found and fixed

Notebooks (github_document, executed; disease vs control, model `~ condition` only, covariates illustrated):
`examples/alz/alz_cacoa.{Rmd,md}` (Grubman 2019, 12 samples), `examples/ms/ms_cacoa.{Rmd,md}` (Schirmer 2019, 21 samples),
`examples/scc/scc_cacoa.{Rmd,md}` (Ji 2020, 18 paired samples); each with `<ds>_notes.md` (data decisions, findings,
QC, runtimes, issue status). Package bugs found by them and fixed (commits 3bf8455, da2c20b, 28df162):
- **Response rows paired with the design by position** in `performLMPermutations()` callers: the cell-density test
  used the Conos sample order against the metadata order, giving a scrambled (MS) or sign-reversed (SCC) density map.
  Responses are now matched by sample name (error when a sample is missing); regression test in `test-kernel-b.R`.
- **Estimability judged on unscaled designs**: numeric covariates in large units (median UMI) made the contrast
  "not estimable" (`isEstimable()` used ginv on X'X); now an SVD row-space projector on unit-norm columns.
- **Arm dispersions not identifiable under a paired location model** (`~ condition + patient`): the columns of
  (R o R) Z coincide and ginv produced s.alt/s.ref = 1.5 identically; `fitDispersion()` now detects this, pools the
  dispersion, and `var`/`total` are NA with the flag `dispersion-not-identifiable` (shift kept).
- Screen / variance partition crashed on a covariate constant within a cell type's samples; `estimateDEPerCellType()`
  tested cell types with 1-2 samples per arm (now `min.samp.per.level`); the sensitivity analysis announced models it
  silently dropped and lost the "current" label; the design check missed excluded factor columns nested in the test
  variable (pooled batches) and did not recognize pairing factors; plus the plot/legend/palette/wording fixes listed
  in the notes.
Kept as observations: max-statistic adjustment over cells floors every adjusted z to 0 for small designs (ALZ 12
samples, MS 99 permutations) — now announced, raw z offered; dispersion R2 normalisation vs its p-value; the
permutation set depends on the metadata row order for a fixed seed.

## 1. Not implemented

- **Robust fitting and weak imputation in the distance engine**: implemented 2026-10-09 (step 4b); calibration
  simulations of the robust and weak-imputation paths are in the slow suite (see `misc/engine_convergence.md`).
- **Residual-vs-covariate diagnostic for density / cluster-free DE**: done 2026-10-10 (`screenResidualCovariates()`,
  `cao$screenResiduals(name)`): residual matrices become sample distances and go through the marginal covariate
  screen; plotted with `plotCovariateScreen(name = "<result>.residual.screen")`.

- **Cluster-free shifts and density refuse whole-factor (K-level) tests** with a message; only contrasts and
  numeric steps run there.
- **Analytic p-values**: the `p.values = "analytic"` option was removed 2026-10-10 (liberal; permutation only). The
  effective-dimension F stays in the table as `p.analytic` for reference.
- **"Later" list from the plan** untouched: Welch-type shift test, robust weights, repeated measures, technical
  noise dispersion covariate, cor-geometry check, metric sensitivity report.

- **Studentized statistic in the per-column fitter** (`statistic = "t"`, default since 2026-10-10; `"coef"` keeps the
  raw contrast): effect sizes are still the raw contrast (`effect`, with `se`); only the p-values changed. Users who
  compared coefficient magnitudes across tests can keep doing so (`effect`), or use predicted values. Block
  relabeling with a numeric nuisance remains approximate (the covariate is permuted with the labels); the auto
  scheme avoids it.

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
  (`rel.change = 0.5`); "sensitive to the adjustment set" when the sign or BH significance differs across models (reworded 2026-10-10 from "depends on adjustment"). These
  thresholds are guesses. The near-collinear proxy covariate in the toy test makes every cell type "sensitive to the
  adjustment set", which is intended but the wording may read as a defect of the method rather than of the covariate.
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

- (resolved 2026-10-10) The weak `IN-PV` panel in the cluster-free shifts was a kernel defect: a sample with fewer
  than `min.n.obs.per.samp` cells in a neighbourhood made every pair involving it missing, and both members of each
  such pair were flagged, so one thin sample removed all samples. Only 13 % of IN-PV cells were tested (average
  5 samples) against 100 % for IN-SST; after the fix every cell is tested and IN-PV has the strongest adjusted
  signal. The former per-cell R loop (`pairVectorToMatrix()`) had the same rule; both now drop one sample at a
  time, the one with the most missing pairs first. Regression test in `test-kernel-cf.R`.

## 4. Housekeeping and cleanup candidates

- (done 2026-10-10, vignette pass) plot fixes found by reading the rendered walkthrough: sample-distance plots drew a
  "Covariate: All" legend when nothing was mapped; `plotShiftDetail()` labelled every sample (ggrepel warning) and
  left an empty grid cell; the influence heatmap's group strip overlapped the sample names; sensitivity strips were
  cut off by the long verdicts; `plotNumberOfDEGenes()` failed without resampling results (and warned instead of
  messaging); the cluster-free panel titles overlapped the plots; the design check printed an empty association
  list with a single covariate and did not note a moderately associated covariate (0.3-0.5); model issues were only
  counted ("see $issues"), now printed inline. `estimateExpressionShiftMagnitudes()` results now have a class with a
  compact print method (normalized effects per cell type with BH p-values and the global p-values).

- (examined 2026-10-10) `private$getTopGenes()` / `getClusterFreeDEInput()`: the time is the one-off extraction of
  the joint count matrix from the Conos object (`extractJointCountMatrix`, 8-13 s, in conos), cached afterwards in
  `cao$cache` (one entry for raw and one for normalized counts); the neighbourhood split costs 3 s. Nothing in
  cacoa to optimize short of caching the normalized matrix next to the raw one, which already happens on first use.

- `misc/` is tracked (plan, notes, validation driver and report, notebook renderer); it is in `.Rbuildignore`.
- (done 2026-10-09) `projdiff.cpp` removed; `estimateCorrelationDistance` is used by cluster-free DE and stays.
- Pre-existing TODOs in `R/cacoa.R` (overall p-adjustment in ontology plots, z-score adjustment in
  cluster-free DE, binary distance in `plotOntologySimilarities`), `R/ontology.R`, `R/cell_density.R`
  (`findScoreGroupsGraph` deprecated?), `R/de_function.R` (bootstrap resampling use).
- Vignette / walkthrough (2026-10-10): `vignettes/walkthrough_short.Rmd` is the single source: one model
  (`~ treatment + age`), one contrast (treated vs control) on a 20-sample subset of the no-batch simulation built by
  `misc/notebook/prepare_sim_dataset.R` into `test/sim_treatment.rds` (57 MB; synthetic `age`, 6 years older in the
  treated arm, no expression effect). Its chunks run only when that file is passed as the `sim.file` parameter; `misc/notebook/render_walkthrough.sh` renders the
  executed `github_document` copy (`vignettes/walkthrough_short.md` + `walkthrough_short_files/`, both in
  `.Rbuildignore`) with the RStudio/quarto-bundled pandoc (the system pandoc 2.5 is below rmarkdown's 2.11.2
  requirement for that format). The vignette build itself evaluates nothing, so `R CMD check` with a vignette
  build will warn (no `inst/doc`); a small shipped example would let it build. The earlier executed `.ipynb`
  was dropped: GitHub's notebook viewer runs its math pass over code cells and mangles every `$` in R code, and
  knitr's native code/output interleaving is lost in any md-to-ipynb conversion.
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
