## [Upcoming]

### Added

- Parameter `assay.name` for Seurat objects in `Cacoa$new` to allow selecting different assays
- `cao$plotMetadataSeparation` function
- `sample.subset` parameter for `plotSampleDistances` and `estimateMetadataSeparation`
- `method="UMAP"` for `cao$plotSampleDistances`, which should work much better on larger sample collections (>50 or so)

### Changed

- Improved algorithm for metadata separation estimation. *The effect should be notisable only on large sample sizes.*

## [0.4.0] - 2022-05-May

### Changed

- Small bug fixes and usability improvements
- Renamed options for `dist.type` in `estimateExpressionShiftMagnitudes` to match the paper figure
- Renamed inner `estimateExpressionShiftMagnitudes` to `estimateExpressionChange` (CRAN requirement)
- Renamed `p.adj.cutoff` to `p.adj` in `estimateDEStabilityPerGene` for consistency with other functions
