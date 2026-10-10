# Model builder rewrite: plan and test plan (2026-10-10, draft for review)

Follows `misc/model_builder_assessment.md` (defects D1-D6, structural issues S1-S8). Nothing here is started.

## 0. Goals and non-goals

Goals: one constructor; coefficients that read as level means for a tested factor (no-intercept coding, the user's
convention); the contrast as the only thing the engines depend on; one estimability rule for every engine; the
grammar working on interaction models; formula handling that never drops a term silently; the pair-model remnants
gone. Non-goals: changing any statistic, the permutation plan, the distance engine or the plotted scale. Every
estimate and p-value that is correct today must be reproduced exactly (same P) by the new builder.

## 1. Decisions (defaults proposed; marked (?) where the user should confirm)

1. **One constructor.** `buildCacoaModel()` is the only builder; `buildDesignMatrices()` becomes a thin exported
   wrapper over the same internals (arguments `numericRefRows`, `validate`, `computeQrZ`, `buildBlocks`,
   `blockVars` accepted with a deprecation message and ignored). Both default formulas collapse into one:
   `~ <test variables> + <block.vars>`.
2. **Coding.** If the tested variable is a factor it is coded without an intercept and placed first in the terms,
   so its columns are the level means at the reference setting of the other covariates (`~ 0 + group + batch`:
   `groupA, groupB, batchb2`). Other factors keep treatment coding; a numeric test variable keeps the intercept.
   Character and logical columns become factors in one place (`prepareDesignData()`), fixing D1. The reference
   level (`chooseReferenceLevel()`) orders the levels (ref first) and names `ref` / `alt`; it no longer changes
   the coding per test, so coefficient names are the same for every test of a model (S2).
3. **Contrast.** Unchanged in substance: endpoint rows from `.oneRowFromFormula()`, `c = row(alt) - row(ref)`;
   `contrast_endpoints_at` kept for the dispersion design. Defaults for the plain grammar on interaction models
   (D2): a factor that interacts with other factors is compared *marginally* (equal weights over the interacting
   factors' levels); a factor that interacts with a numeric covariate is compared at the covariate's anchor
   (`numeric.ref`, mean by default); a numeric test variable that interacts with a factor gets its slope averaged
   over the factor's levels (?). Each default adds a `note` to the model issues and to the printout; `at =` /
   `over =` in a structured test override it. Contrast types kept: `simple`, `marginal`, `coef` (expert escape
   hatch); `lincomb` dropped (?) since `marginal` with numeric weights covers it.
4. **Split by direction, not by column (S3, D5).** For the per-column fitter the design is
   `X = F c / (c'c)` (one column named after the test, coefficient = effect) and `Z = F N` with `N` an
   orthonormal basis of the null space of `c` (columns `nuisance1..`). `contrast.X = 1`. Term tests use the same
   construction with a contrast matrix: `C` = the K-1 differences between level means (adjusted for the other
   covariates), `X = F C`, `Z = F N(C)`. `chooseBaselinesForSpec()` and `splitByContrast()` disappear.
5. **Estimability in one place.** A helper `subsetDesign(F, c, rows)` (drop all-zero and constant columns among
   the rows, check `c` lies in the row space, return the reduced design and contrast) replaces
   `repairDesignAfterRowSubset()` in DE and the inline code in `pairwiseEffectsFromDesign()` /
   `termEffectsFromDesign()`. The same rule goes into `fit_and_randomize()` per NA pattern (D4): drop constant
   columns of the units, pseudo-inverse fit, NaN only when the contrast is not estimable (reported as such).
6. **Formula hygiene (S6, D3).** `terms()` objects throughout; a non-varying term is dropped with a model issue
   of severity `note` and the printout shows the formula actually used; random-effect terms are an error with the
   suggestion to use the variable as a fixed effect or as `block.vars`; `checkFormula()` shrinks to type checks.
7. **Term tests under interactions.** The term contrast matrix spans the main effect of the factor at the
   reference setting of interacting covariates (?), not main effect plus interaction; a note says so. DE term
   tests (edgeR, limma, DESeq2) switch from `coef = term.cols` to the contrast matrix, which makes them
   coding-independent.
8. **CoDA.** Coefficients, per-coefficient effects and predicted compositions come from the full-design OLS fit
   on F (always available and coding-interpretable: level means in ILR space); p-values and null loadings come from
   the fitter's `effect` / `effects.perm` under the chosen scheme. The `perm.method`-dependent choice of design
   for predictions disappears.
9. **R6 object.** `estimateClusterFreeDE()` goes through `resolveModel()` like every other method (S1).
   `self$formula` / `self$contrast` are removed; `ref.level`, `target.level`, `sample.groups` stay as derived
   read-only fields (plots and the palette use them). `getCoreSamples()` keeps its new definition.
10. **Retire.** `diagnoseDesign()`, `emitDiagnostics()`, `reportContrastInfo()`, `qrZ`, `baselines_used`,
    `numericRefRows`, `varsFromFormula()` (if its two tests move), `repairDesignAfterRowSubset()`,
    `isContrastNonzero()`, the `pair.set` plot argument, `lincomb` (?). `designIssues()` gains the condition
    number of F.

Design object after the change (what consumers read): `F`, `contrast.F`, `X`, `Z`, `contrast.X`, `meta`,
`formula_used`, `contrast_spec`, `contrast_endpoints_F`, `contrast_endpoints_at`, `numeric_ref_used`,
`contrast_label`, `contrast_endpoint_labels`, `term.contrast` (matrix, term tests), `term.variable`, `notes`.
Dropped: `contrast_endpoints_X` (CoDA no longer needs it), `term.cols`, `qrZ`, `diagnostics`, `baselines_used`.

## 1b. Contrast handling (settled 2026-10-10)

One model answers one question: `setModel(formula, test)` names one primary contrast and every method runs on
it; further tests remain possible as an explicit list (first = primary) but are not the default presentation.
Guess where the guess is safe, ask where it is not:
- two-level factor: compare alt vs ref, reference by the existing rules; the printout states the guess, its reason
  and the override syntax (`"Diagnosis: IPF vs Control"`);
- factor with more than two levels and no comparison given: stop with a message listing the levels, the proposed
  reference, the comparison syntax and the whole-factor syntax (`"Group: all"`); no silent term test;
- numeric variable: slope per unit (per SD reported alongside), stated in the printout;
- interaction models with the plain grammar: marginal comparison with equal weights (factor x factor), anchor at
  the covariate mean (factor x numeric), slope averaged over the factor's levels (numeric x factor); recorded as a
  note with the `at` / `over` overrides;
- structured contrasts (`at`, `over`, coefficient weights) stay the expert path; `lincomb` is dropped.
Open points of §1 decided: numeric slope averaged over an interacting factor; `lincomb` dropped; a term test under an
interaction model is the main effect at the reference setting of the interacting covariates (note recorded).

## 2. Steps, each with its OODA checkpoint

**Step 0: freeze what is right.** Turn `misc/scripts/model_scenarios*.R` into `tests/testthat/test-model-builder.R`
with expectations for the cases that behave correctly today (T1, T2 partly, T3, T7 partly, T13), and add the
engine-invariance fixtures (T4): for six designs on the grid metadata, store the shift / var / total estimates,
p-values (fixed P) and the per-column fitter's `effect` / `se` / `pval` as expected values computed from a
treatment-coded design built by hand. Checkpoint: suite green on the current code; the fixtures are the oracle
for every later step.

**Step 1: new builder internals.** `designCoding()` (term order, no-intercept for the tested factor,
`prepareDesignData()` for character / logical), `contrastFromEndpoints()`, `splitByDirection()`,
`termContrastMatrix()`, `subsetDesign()`; `buildDesignMatrices()` rewritten on them; consumers adapted to the
new fields (`performLMPermutations()`, `lmCoda()`, `estimateClusterFreeDE_LM()` means, `termEffectsFromDesign()`
reduced design, `pairwiseEffectsFromDesign()` via `subsetDesign()`, cluster-free shifts, DE). Checkpoint: T2,
T3, T4, T5 (unit part), T8, T11; fast suite; the five Rd files regenerate.

**Step 2: grammar and formula hygiene.** Interaction defaults with notes, logical variables, terms-based pruning
with issues, random-effect error, printout. Checkpoint: T1, T6, T7, T13; fast suite.

**Step 3: estimability in kernel B and DE.** C++ change in the per-pattern design (constant-column drop,
estimability, pseudo-inverse); DE term tests via contrast matrices; `repairDesignAfterRowSubset()` removed.
Checkpoint: T9, T10 (slow, SCC example), T5 calibration; slow suite.

**Step 4: R6 and retirement.** `estimateClusterFreeDE()` through `resolveModel()`, legacy fields, retire list,
docs, `misc/issues.md`, notebook re-render, `R CMD check`. Checkpoint: T12.

Order of commits: one per step, one-line messages. Each step ends with the fast suite; steps 3 and 4 with the
slow suite and the check.

## 3. Test plan

| id | what | how | where |
|---|---|---|---|
| T1 | grammar table | every row of the grammar (variable; `var: a vs b`; triple; `all`; several; structured simple / marginal / coef; numeric with step) for factor, character, logical and numeric columns resolves to the specified kind, levels, label and reference reason; wrong level / variable / constant column give the specified errors | `test-model-api.R` (extend), `test-model-builder.R` |
| T2 | coding and names | for each design in the grid (2-level; 3-level contrast; 3-level term; `~ group * batch` simple-at / marginal / cell; numeric; `log(age)`; `poly(age, 2)`; `factor(batch)`; explicit `~ 0 +`; `block.vars`): column names as specified, tested-factor coefficients equal `lm(y ~ 0 + group + ...)` level means, contrast vector equals row(alt) − row(ref), coefficient names identical across the tests of one model | `test-model-builder.R` |
| T3 | agreement with `lm` | `effect`, `se`, `t`, `df` from `performLMPermutations()` (both schemes, `statistic = "t"`) equal `lm()` on the full model for every T2 design; `regressionWeights(design) %*% y` equals the `lm` contrast | `test-model-builder.R`, `test-kernel-b.R` |
| T4 | engine invariance | shift / var / total estimates and p-values (same `P`) and the fitter's output are unchanged against the step-0 fixtures built on hand-made treatment-coded designs; weighted / robust / impute_weak paths included | `test-model-builder.R` (fixtures in `helper-grid.R`) |
| T5 | Freedman-Lane reduced model | unit: kernel-B FL on the direction split equals a brute-force R reference that residualizes on `F N` (per NA pattern); calibration: 3-level factor, strong planted third-level effect, null rejection of the G2 vs G1 test within the binomial band at 1500 columns, for `none` / `huber` / `winsor` x `drop` / `impute_weak` | `test-kernel-b.R`, slow suite |
| T6 | interaction defaults | plain `test = "group"` on `~ group * batch` equals the structured marginal contrast (equal weights) and records a note; factor x numeric equals the structured `at = mean`; numeric x factor equals the structured marginal over the factor; explicit `at` / `over` override; the printout names the default | `test-model-builder.R` |
| T7 | formula hygiene | constant term dropped with a `note` issue and `formula_used`; random-effect term error with suggestion; interaction with a dropped variable pruned consistently; missing variable, id-like column, duplicated terms; listwise deletion recorded | `test-model-builder.R` |
| T8 | term tests | location F from the K-1 contrast matrix equals `anova(lm)` F for a one-dimensional response (l2 distance); equals the current `term.cols` F on a treatment-coded design (fixtures); under `~ Group * Batch` the term is the main effect at the reference batch and a note is recorded; dispersion term test unchanged | `test-pairwise-effects.R`, `test-model-builder.R` |
| T9 | estimability in kernel B | a column whose observed rows lack a nuisance level: finite `effect` / `se` / `p` equal to `lm` on the subset with the constant column dropped; a column lacking a contrasted level: NA with the "not estimable" reason; `impute_weak` unchanged; the shift engine and the fitter agree on which columns are skipped | `test-kernel-b.R` |
| T10 | DE | DESeq2 Wald, edgeR QL and limma-voom log fold changes and p-values under the new coding equal the treatment-coded results (SCC example); term tests via contrast matrices equal `coef = term.cols` on treatment coding; per-cell-type subsetting through `subsetDesign()` keeps the same cell types as before | `test-example-datasets.R` (slow), `test-track-d.R` |
| T11 | CoDA | per-coefficient effects are level means in ILR space; predicted baseline / target compositions equal the back-transformed endpoint rows; loadings and their p-values equal the step-0 fixtures | `test-track-d.R`, `test-model-builder.R` |
| T12 | end to end | fast suite, slow suite (validation rows, example datasets, weighted, cluster-free), `R CMD check`, notebook re-render with the composition / density / cluster-free panels compared to the current ones | CI of each step |
| T13 | printout snapshots | `format.cacoaModel()` for six designs (two-level, three-level contrast, term, interaction marginal default, numeric, with notes) as `expect_snapshot()` | `test-model-builder.R` |

Coverage check at the end: every function left in `model_matrices.R` / `model_api.R` is reached by at least one
test (`covr` on the two files), and every row of §1 has a test id.

## 4. Risks and open points

- Coefficient names change for anyone reading `coef` from the fitter or CoDA (`groupA, groupB` instead of
  `(Intercept), groupB`); documented in the news entry.
- The interaction defaults turn an error into a result; the note in the printout is the safeguard.
- DE term tests via contrast matrices give the same F up to the method's own dispersion handling; T10 pins it.
- Estimated size: step 0 half a day; step 1 one and a half days; step 2 half a day; step 3 one day; step 4 half a
  day; four days with the slow runs.
- Decisions marked (?) in §1: numeric-slope averaging over an interacting factor, dropping `lincomb`, term tests
  as main effect at the reference setting.
