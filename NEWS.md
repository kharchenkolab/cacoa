# cacoa (development, branch dev_lm)

## Cluster-free analyses

- The density map and the expression-shift map test every cell in a kNN neighbourhood of the cell graph (`estimateCellDensity(smoothing = "knn")`, the default; `estimateClusterFreeExpressionShifts(k)`), sized by `defaultNeighbourhoodSize()`; the heat kernel remains available with its width settable in cells (`target.size`).
- Both maps average the per-cell statistic over the cell and its graph neighbours under every relabeling (`post = "averaged"`, the default; `"inherit"` and `"none"` as options) and adjust across cells with the pooled permutation false-discovery rate of the standardized (shifts) or studentized (density) statistic; `pooled.statistic` switches between the two.
- Unadjusted z-scores are the signed tail z of the pooled permutation p-value (`pooled` / `z.pooled`), adjusted z-scores the same transform of the q-values; `plotDiffCellDensity(what = "effect")` shows the effect.
- The per-cell kernels run in C++ with one shared accumulator (no permutations x cells matrix in memory), a closed-form fitter for the default density case and a two-group closed form for the permuted shift statistic.
- `cao$checkNeighbourhoods()` reports unit sizes, cells per sample per unit and purity against the annotation for a neighbourhood size.

## Permutations

- Every analysis takes the number of permutations from the object option `n.permutations` (999) unless a value is passed.
