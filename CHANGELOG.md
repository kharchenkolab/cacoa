## [Upcoming]

### Changed (expression-shift redesign and workflow API)

- Expression shifts are estimated with an individual-level model fitted to sample-sample distances (Gower
  form): `shift`, `var` and `total` per cell type, with permutation p-values (block randomization within
  strata when all covariates are discrete, Freedman-Lane otherwise), BH and max-statistic adjustment across
  cell types and a global p-value per test. The pair-level regression (`pairFormula`, `pairContrast`,
  `dist.type`, `perm.method`, `robust.method`, `na.mode`, gene focusing) is removed; the arguments are
  accepted with a message (the pair arguments are errors).
- One stored model for every analysis: `Cacoa$new(..., test =)` or `cao$setModel(formula, test)`; the `test`
  grammar (variable, K-level factor, numeric step, several variables, "all", "var: alt vs ref", triple,
  structured contrast). Reference levels follow control-like names, explicit factor order, or frequency.
  `cao$setOptions()` holds shared defaults (`n.permutations`, `seed`, `alpha`, ...).
- `Cacoa$new()` no longer requires a model and prints a metadata audit (`cao$describeMetadata()`).
- Sample weights for cell density are the regression weights of the contrast.
- `estimateMetadataSeparation()` is now the marginal covariate screen on the joint distance.
- `plotExpressionShiftMagnitudes()` is a dot plot with jackknife intervals and a provenance subtitle;
  `plotExpressionShiftResiduals()` is replaced by `plotSampleInfluence()`.
- OpenMP was replaced by the `sccore` thread pool; the C++ code no longer requires OpenMP.

### Added

- `cao$checkDesign()` / `cao$plotDesign()`: associations among covariates, balance, aliasing, GVIF, degrees of
  freedom, permutation feasibility, over-adjustment candidates.
- `cao$screenCovariates()` with `plotCovariateScreen()`, `plotCovariateSummary()`, `plotVariancePartition()`.
- `cao$checkSensitivity()` / `cao$plotSensitivity()`: alternative covariate sets and single-sample influence.
- `cao$plotShiftDetail()`, `cao$plotSampleDistances(adjust.for =)`, `cao$getSampleGroups()`.
- Whole-factor tests: expression shifts (location / dispersion), DESeq2 LRT, edgeR / limma F, composition.
- Cluster-free expression shifts on the same engine with shared permutations.
- Exported engine functions: `buildCacoaModel()`, `resolveTests()`, `testPairwiseEffects()`, `testTermEffects()`,
  `expressionShiftsForModel()`, `screenCovariates()`, `checkSensitivity()`, `pseudobulkPerCellType()`,
  `sampleDistanceMatrices()`, `permutationPlan()`, `drawPermutations()`.
- A testthat suite (fast unit tests, example-dataset tests, slow simulation tests).

### Added (earlier, unreleased)

- Parameter `assay.name` for Seurat objects in `Cacoa$new` to allow selecting different assays
- `cao$plotMetadataSeparation` function
- `sample.subset` parameter for `plotSampleDistances` and `estimateMetadataSeparation`
- `method="UMAP"` for `cao$plotSampleDistances`, which should work much better on larger sample collections (>50 or so)
- Updated documentation for most functions in R6 class object
- Updated functionality (estimations and plotting) for ontologies, especially GSEA and DO as well as families
- Added parameter `only.family.children` to plotOntologySimilarities
- Added parameter `families` to plotNumOntologyTermsPerType
- Added example data for testing
- Several minor bug fixes

### Changed

- Improved algorithm for metadata separation estimation. *The effect should be noticeable only on large sample sizes.*

## [0.4.0] - 2022-05-May

### Changed

- Small bug fixes and usability improvements
- Renamed options for `dist.type` in `estimateExpressionShiftMagnitudes` to match the paper figure
- Renamed inner `estimateExpressionShiftMagnitudes` to `estimateExpressionChange` (CRAN requirement)
- Renamed `p.adj.cutoff` to `p.adj` in `estimateDEStabilityPerGene` for consistency with other functions
