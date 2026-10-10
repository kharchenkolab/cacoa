# Model-building code: assessment (2026-10-10)

Scope: `R/model_matrices.R` (1080 lines, the expert constructor `buildDesignMatrices()` and the synthetic-contrast
machinery), `R/model_api.R` (518 lines, the workflow constructor `buildCacoaModel()`, the test grammar, the
reference-level rules), the R6 entry points (`Cacoa$new()`, `setModel()`, `resolveModel()`, `buildModel()`,
`syncLegacyFields()`) and the per-cell-type design handling in the DE code (`repairDesignAfterRowSubset()`).
Scenario scripts: `misc/scripts/model_scenarios.R`, `misc/scripts/model_scenarios2.R` (run against the scratch build).

## 1. What the pipeline does today

```
test grammar  --resolveTests()-->  test spec (kind, variable, levels, contrast spec)
                                        |
buildCacoaModel(): listwise deletion, default formula (~ test vars + block.vars), per test:
   contrast ->  buildDesignMatrices():  checkFormula -> pruneFormulaByData -> normalizeContrastSpec
                                        -> chooseBaselinesForSpec (re-levels factors so contrasted levels get columns)
                                        -> buildFullDesign (model.matrix + terms/xlevels attributes)
                                        -> buildSyntheticContrast (.oneRowFromFormula rows at endpoint settings)
                                        -> splitByContrast (X = columns with non-zero weight, Z = the rest)
                                        -> diagnoseDesign (unused), reportContrastInfo (silent)
   term     ->  buildTermDesign():       model.matrix, term columns by `assign`
   permutationPlan() per test, designIssues(), primary design copied to the top level of the model
```

Consumers of a design: `F`, `contrast.F`, `meta` (every engine), `X` / `Z` / `contrast.X` (per-column fitter,
CoDA, cluster-free DE), `contrast_endpoints_at` + `endpointRowsFromDesign()` (dispersion design of the shift
engine), `contrast_endpoints_X` / `contrast_endpoint_labels` / `contrast_label` (CoDA predicted compositions),
`contrast_spec`, `formula_used`, `term.cols` / `term.variable`. Never read: `qrZ`, `diagnostics`, `baselines_used`.

## 2. Defects found by the scenario runs

| # | scenario | result | severity |
|---|---|---|---|
| D1 | logical test variable (`treated = TRUE/FALSE`) | `buildCacoaModel()` fails ("contrasts can be applied only to factors with 2 or more levels"); works once converted to a factor. `resolveTests()` handles logicals, `prepareDesignData()` converts only character columns, `.oneRowFromFormula()` then treats the logical as numeric. | bug |
| D2 | interaction formula with the plain grammar (`~ group * batch`, `test = "group"`, or `test = "batch"`) | "Ambiguous triple" error. The only way through is a structured list with `at =` or `over =`. The same for a factor x numeric interaction (`~ Group * age`, `test = "Group: G2 vs G1"`), although numeric anchors exist exactly for that case. | usability gap |
| D3 | constant covariate in the formula (`~ group + site + batch`, `site` constant) | silently pruned: `model$formula` still names `site`, `formula_used` does not, no issue recorded (`buildCacoaModel()` calls the builder with `verbosity = "none"`). | bug (silent) |
| D4 | per-column fitter on a column whose observed rows lack a level of a nuisance factor (cluster-free DE / density with `na.mode = "drop"`) | rank-deficient subset design -> whole column NaN with a "model fit failed" warning, although the contrast is estimable. The shift engine handles the same case (drops empty columns, `isEstimable()`). | gap in kernel B |
| D5 | per-column fitter, Freedman-Lane, contrast on a 3-level factor | `X` holds two factor columns and `Z` only the intercept and covariates, so the reduced model lacks the non-tested direction of the factor (G1+G2 vs G3). Null rejection with a strong planted G3 effect: 0.065 (vs 0.060 when the reduced model holds that direction; 1500 null columns). The distance engine residualizes on `F %*% contrastNullBasis(c)` and does not have this issue. | minor, theoretical |
| D6 | explicit factor order beats a control-like name (`factor(levels = c("disease", "control"))` -> reference "disease") | by design (rule 1 before rule 2) but surprising; listed among the open defaults. | decision |

## 3. Structural problems

1. **Two constructors with different defaults.** `buildDesignMatrices()` without a formula builds `~ 0 + vars`
   (no intercept, saturated coding); `buildCacoaModel()` builds `~ test vars + block.vars` with an intercept.
   `estimateClusterFreeDE()` still calls `buildDesignMatrices()` directly with the legacy `self$formula` /
   `self$contrast` fields and does not accept the `test` grammar; every other method goes through `resolveModel()`.
2. **Baseline switching per test.** `chooseBaselinesForSpec()` re-levels the tested factor so that both contrasted
   levels get their own column (the non-contrasted level becomes the baseline). Consequences: coefficient names
   differ between tests of one model (`GroupG1, GroupG2` for G2 vs G1; `GroupG1, GroupG3` for G3 vs G1;
   `GroupG2, GroupG3` for the term test), five cases of heuristics (simple / interaction cell / marginal /
   lincomb / other), and `baselines_used` that nothing reads. The re-levelling exists only so that
   `splitByContrast()` can split by columns; a split by the contrast direction (X = F c / c'c, Z = F N) needs no
   re-levelling and is what the distance engine already does (D5).
3. **Column-wise X / Z split is the wrong primitive for the fitter.** It is what causes D5 and it is why the
   fitter's "coefficients" for Freedman-Lane are X-space coefficients whose meaning depends on the coding. For
   every consumer the contrast estimate (`effect`) is what matters; CoDA additionally wants the endpoints, which
   are rows in F-space and independent of the split.
4. **Four contrast types plus a flag** (`simple`, `marginal`, `lincomb`, `coef`, `plain_triple`) handled in
   three switch ladders (`normalizeContrastSpec`, `chooseBaselinesForSpec`, `buildSyntheticContrast`) and again
   in `resolveTests()` / `contrastLabel()`. `lincomb` and `coef` have no workflow entry point beyond the expert
   argument and no users in the package.
5. **Unused machinery** in `buildDesignMatrices()`: `computeQrZ` / `qrZ`, `validate` / `diagnoseDesign()`
   (90 lines; its result is never read; its only effect is a warning at verbosity "warn"), `buildBlocks`,
   `blockVars` (only forwarded to the default formula), `numericRefRows`, `reportContrastInfo()`,
   `varsFromFormula()`, `tol`. `checkFormula()` rewrites random-effect terms with string regexes.
6. **Formula handling is string-based in places** (`checkFormula`, `pruneFormulaByData` rebuilds the formula
   from term labels, `mainEffectTerms`), which is how D3 can drop a term without a trace.
7. **Term tests are a separate, thinner path** (`buildTermDesign()`): default coding, no issues on aliasing within
   the term, the stored "reference" reason is not what the design uses, interaction columns are included in the
   term silently (`Group` under `~ Group * Batch` tests main effect plus interaction).
8. **Pair-model remnants**: `repairDesignAfterRowSubset()` in the DE code (QR-pivot column dropping with the
   contrast truncated to the kept columns: correct only when the dropped column is absorbed by the intercept),
   `getCoreSamples()`, the `sample.groups` / `ref.level` / `target.level` / `formula` / `contrast` legacy fields
   kept in sync by `syncLegacyFields()`, and `pair.set` arguments in plots.

## 4. What is sound and should stay

- `resolveTests()` grammar and `chooseReferenceLevel()` (tested row by row); `describeMetadata()`.
- `.oneRowFromFormula()`: endpoint rows are built by evaluating the formula on new data, so `log(age)`,
  `poly(age, 2)`, `factor(batch)` and interactions are handled uniformly (scenarios 4c, 4d, 5a, 3b-3d pass).
- Endpoint storage (`contrast_endpoints_at`) re-evaluated on the dispersion design.
- `designIssues()` / `permutationPlan()` summaries in the model printout.

## 5. Proposed shape (for discussion, not started)

One constructor, `buildCacoaModel()`, with `buildDesignMatrices()` reduced to an internal helper:

1. **Design**: `F = model.matrix(formula, meta)` with the default coding (reference level of the *tested* factor
   = `chooseReferenceLevel()`, which is the only re-levelling needed); character and logical columns become
   factors in one place.
2. **Contrast**: endpoints as F-space rows (`.oneRowFromFormula()`, unchanged); `c = row(alt) - row(ref)`.
   Interactions with the plain grammar get a documented default (marginal, equal weights over interacting
   factors; numeric anchors at the mean) and a note in the model printout instead of an error (D2).
3. **Split for the fitter by direction, not by column**: `X = F c / c'c` (one column, coefficient = effect),
   `Z = F N` with `N` the null basis of `c`. Removes `chooseBaselinesForSpec()`, `splitByContrast()`, D5, and the
   per-test coefficient-name drift. CoDA keeps per-coefficient effects from `F` (block path) and the endpoints.
4. **Estimability in kernel B** (D4): per NA pattern drop all-zero / constant columns and check `c` against the
   row space, as `pairwiseEffectsFromDesign()` does; fit with the pseudo-inverse.
5. **Formula hygiene without strings**: `terms()` objects throughout; dropped terms become model issues (D3);
   random-effect terms rejected with a clear message rather than rewritten.
6. **Retire**: `diagnoseDesign()`, `qrZ`, `numericRefRows`, `lincomb` / `coef` contrast types (or keep `coef`
   only as the expert escape hatch), `repairDesignAfterRowSubset()` (replace with the same estimability check),
   the legacy fields and `syncLegacyFields()` once plots read the model, the direct `buildDesignMatrices()` call in
   `estimateClusterFreeDE()`.
7. **Tests**: the scenario scripts become `test-model-builder.R` (one expectation per row above), plus
   equality of `effect` / `se` with `lm()` for every contrast type on the full model.

Rough size: items 1-3 and 5 are a rewrite of `model_matrices.R` to about a third of its length; 4 is a C++
change in `fit_and_randomize()`; 6 touches DE, plots and the R6 class.
