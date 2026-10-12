# Cluster-free multiplicity: calibration and power benchmark (2026-10-11)

Question: how should per-cell (continuous-space) results of the differential density and the cluster-free expression
shift be adjusted for multiplicity, and how do the alternatives compare with Milo on the same graph? The per-cell
statistics are strongly dependent (overlapping neighbourhoods of the same samples), permutation p-values are floored at
1 / (B + 1), and the earlier comparison with Milo (`misc/milo/milo_comparison.md`) rested on one simulation and three
real contrasts without a calibration study.

Scripts (this directory; scratch outputs under `$S/bench/`):

| script | purpose |
|---|---|
| `common.R` | loaders, density matrices (heat kernel and 1-hop neighbourhood counts), the method battery |
| `extract_pf.R` | light-weight PF inputs (the SCC and MS inputs are the Milo comparison's `inputs_<ds>.rds`) |
| `null_relabel.R <ds> <n> <B> <cores>` | null study: arm labels re-drawn across samples, balanced within the true condition, permutations blocked on it (`CACOA_NULL_DIR`, default `null`; the final round is `null3`) |
| `semi_synthetic.R <ds> generate\|run [B] [cores] [unpaired\|paired]` | control-only semi-synthetic sets with a planted depletion of one cell type |
| `analyze_null.R`, `analyze_semi.R` | figures and tables (`figures/`) |
| `tail_collect.R <ds> <n.cells> <B> <cores>` | per-cell tail study: full permutation matrices of a cell subset under one null relabeling |
| `tail_analyze.R <rds> <prefix>` | tail approximations evaluated on held-out null draws (`CACOA_BSMALL`: training sizes) |
| `tail_diag.R` | the relabeling-count diagnostic of the pooled tail (`$S/bench/tail/`) |
| `planted_relabel.R <ds> generate\|run [B] [cores]` | unpaired power study: balanced relabeling (9-11 whole samples per arm) + planted depletion of one type in arm B; cacoa with the pooled-null variants (`$S/bench/planted/<ds>/<id>/`) |
| `milo_decompose.R <ds> [B] [cores]` | Milo comparison decomposed into unit (refined kNN neighbourhoods of ~95 cells from Milo's run vs cacoa's 1-hop) x test (edgeR QL, voom, cacoa's permutation LM, voom-weighted permutation prototype) x adjustment, with TMM and library-size normalization |
| `analyze_planted.R [out.dir]` | figures and tables of the planted study (`figures/planted_*`), compositional truth |
| `pilot_pooled.R <B> <cores>` | pilot of the pooled p-values in the BH-type adjustments and the Cauchy summary: one SCC null relabeling and four paired semi-synthetic sets |
| `diag_pooled.R <set> <B> <cores>` | per-cell vs raw pooling, unpaired model vs donor covariate, direction of the extra calls on one semi-synthetic set |
| `$S/bench/run_milo_semi.R` | Milo (miloR 1.6.0) on the semi-synthetic sets (subagent) |

Results: [null calibration](#results-null-by-relabeling), [semi-synthetic power and false discoveries](#results-semi-synthetic-depletion),
[per-cell tail approximations](#results-per-cell-tail-approximations), [conclusions](#conclusions).

## Methods under test

Per-cell density tests (the contrast of the per-sample densities at every cell, fitted with the package's
linear-model permutation engine, 1,999 permutations, studentized contrast):

- **heat kernel** (the package's `estimateCellDensity(method = "graph")`, beta 30): a soft neighbourhood of every cell;
- **1-hop neighbourhood counts** (`smoothing = "neighbourhood"`): the number of each sample's cells among the cell's
  graph neighbours, tested on the log2 scale (`transform = "log"`): the unit of a Milo-type abundance test.

Per-cell p-values: parametric (moderated t, limma's variance shrinkage across cells) and permutation. Adjustments:
BH over cells, Milo's spatially weighted BH (`"spatial-fdr"`), BH over a covering set of index neighbourhoods
with inheritance (`"neighbourhood"`), the max statistic over cells (`"max-t"`, the former default), and the pooled
permutation FDR (`"perm-fdr"`, SAM-type: the expected number of cells beyond a threshold under the permutations,
pooled over cells, against the observed number; no per-cell floor).

Cluster-free expression shifts (500 or 1,000 most expressed genes, 1,999 permutations): permutation p, the normal tail
of the permutation-standardized z (`p.source = "z"`), the parametric F; adjustments BH, spatial FDR, pooled
permutation FDR.

Group-level summaries per cell type: the Cauchy combination of the cells' p-values (Liu & Xie 2020) and the exact
group permutation test (the mean of the group's per-cell statistics against the same mean under every permutation;
the permutations are shared by all cells, so the dependence is carried exactly). The Cauchy combination is applied to
the permutation p-values and, for the density tests, also to the parametric ones.

**Kernel fix found by this benchmark.** The group permutation test of the expression shifts was anti-conservative in
the first null round (9-22 % of the cell types at the 5 % level) while the per-cell p-values were uniform. The cause
was in the shift kernel: in a neighbourhood where some samples are absent, a shared permutation was *induced within
the present samples* (ranks of the global permutation among them). Per cell this is an exact conditional test (the
level counts of the cell stay fixed), but the relabelings are then incoherent across cells, so any statistic pooled
over cells (group means, the max statistic) is compared with a null that is too narrow. On iid synthetic distance
matrices with 6 of 16 samples absent per cell the group test rejected 31 % of the groups at the 5 % level; with every
sample present, 5-8 %. The kernel now applies the global relabeling in every cell (a present sample takes the design
row of the sample it is mapped to) and skips a cell under relabelings that leave fewer than `min.samp.per.level`
present samples at a compared level (`n.perm` reports the count per cell); the Freedman-Lane residual scheme keeps the
within-cell permutation. After the fix the same synthetic setting gives 7.8 % at 5 % with and without absent samples,
and the end-to-end synthetic run (iid Poisson counts, 16 samples, 40 random labelings) gives 5.0 % at 5 % and 10.0 %
at 10 %. The density tests were not affected (every sample is present in every cell's density vector).

Milo (miloR 1.6.0) on the same Conos graph: graph-refined neighbourhoods (`prop` 0.1), TMM + quasi-likelihood F
test, graph-overlap spatial FDR; per cell the smallest spatial FDR of the neighbourhoods containing it.

## Designs

1. **Null by relabeling** (`null_relabel.R`): SCC (9 tumour / 9 normal), MS (12 / 9), PF (12 / 10). The arm labels are
   re-drawn so that each true condition is split in half between the arms (the real condition effect is balanced, only
   the arm labels are random), and the model blocks the permutations on the true condition (`block.vars`), so that
   the observed labeling is exchangeable with the permutations; 50 relabelings for SCC and MS, 25 for PF. Two earlier
   rounds were discarded: unrestricted relabeling (`$S/bench/null/`; arms with unequal tumour content carry a real
   composition difference, and one such SCC relabeling produced a 6 % blob of calls under every FDR method), and
   balanced relabeling with unblocked permutations (`$S/bench/null2/`; the observed labeling is then a restricted draw
   compared with unrestricted permutations, which made every permutation test conservative per cell).
2. **Semi-synthetic depletion** (`semi_synthetic.R`): control samples only (SCC 9 normals, PF 10 controls), split into
   two pseudo-arms; in the treated arm 50, 80 or 100 % of one cell type's cells are removed (a small, a medium and a
   large type: SCC Tcell 2.4 %, LC 6.7 %, KC_Basal 18 %; PF NK cells 2.3 %, Ciliated 6.7 %, AT2 14 %), three replicate
   arm splits each, plus three control-only splits with nothing removed. Two regimes: **unpaired** (4 + 5 or 5 + 5
   whole samples per arm) and **paired** (every control sample split at random into two halves, one per arm, blocked
   by donor: `block.vars = "patient"` for cacoa, `~ patient + condition` for Milo). Truth: the remaining cells of the
   planted type should be called "down"; a called cell of another type that is not a graph neighbour of a planted
   cell is a false discovery.

## Results: null by relabeling

Figures: `figures/null_qq_cells.png` (per-cell p-values pooled over cells and relabelings, one panel per test and
dataset), `figures/null_calls_cells.png` (fraction of cells called at FDR 5 % per relabeling and adjustment),
`figures/null_qq_groups.png` (per-cell-type p-values of the three group summaries); tables
`figures/null_calls_summary.csv`, `figures/null_groups_summary.csv`. 1,999 permutations per relabeling (per-cell floor
1 / 2000, -log10 3.3).

**Per-cell p-values.** The permutation p-values of all three tests lie on the diagonal up to their floor in every
dataset: the per-cell permutation test is exact under the stratified relabeling, as it should be when the observed
labeling is exchangeable with the permutations. The parametric moderated t of the 1-hop counts is close to uniform
with a slightly heavy far tail (SCC, PF beyond 10^-4); the parametric p of the heat kernel is conservative over the
whole range (at expected 10^-3 the observed quantile is 10^-2 to 10^-2.3). The normal tail of the permutation-
standardized shift z is strongly anti-conservative (observed 10^-3 at expected 10^-2, 10^-5 at expected 10^-3): the
per-cell null of the shift statistic is right-skewed, so a normal tail on its z overstates the evidence. The
parametric F of the shift is useless: its p never goes below 0.01.

**Calls at FDR 5 %.** With the heat kernel every adjustment calls essentially nothing: the pooled permutation FDR calls
a mean of 0.03 % of the cells on SCC (one relabeling 1.5 %) and 0 on MS and PF; BH on permutation p calls 0 except
one MS relabeling (2 %); the parametric adjustments call 0. The same holds for the 1-hop counts on SCC and MS (every
method at most 0.4 % in a few relabelings). For the shifts the pooled FDR, BH on permutation p, the spatial FDR and
BH on the parametric F call nothing; BH on the normal-tail z calls cells in 96-100 % of the relabelings (mean
0.5-0.7 %, up to 9-13 % of the cells): the former z-score plotting scale must not be used for calling. At FDR 10 %
the pooled FDR produces calls in 2-8 % of the relabelings (one MS heat-kernel relabeling 18 %, otherwise at most 1 %).

**PF, 1-hop counts.** Two of the 25 PF relabelings produce 4-10 % calls under *every* adjustment, the parametric ones
included (pooled FDR 10 %, BH on permutation p 9 %, spatial FDR 9 %, neighbourhood test 21 %, BH on parametric p
5 %), and three relabelings put 8-19 of the 31 cell types below the group BH threshold at once. The marginal per-cell
p-values of the same test stay uniform, so this is not a per-cell bias but dependence: the count of a sample among
a cell's neighbours responds to the sample's total cell number and graph connectivity, which differ strongly among
the PF samples, and that factor is shared by every cell type. Under a relabeling that happens to split those samples
unevenly, all types move together. The heat kernel, whose per-sample densities are normalized, shows nothing in the
same relabelings.

**Group-level summaries** (19 cell types in SCC and MS, 31 in PF; raw rejection rate at 5 % and share of relabelings
with at least one BH call at 5 %):

| test | group permutation test | Cauchy (permutation p) | Cauchy (parametric p) |
|---|---|---|---|
| heat kernel, SCC / MS / PF | 4.5 / 4.4 / 4.1 %; 0 / 2 / 4 % | 5.4 / 5.1 / 3.6 %; 2 / 0 / 0 % | 3.1 / 0.9 / 1.0 %; 2 / 0 / 0 % |
| 1-hop counts, SCC / MS / PF | 4.4 / 4.3 / 10.3 %; 4 / 4 / 16 % | 5.3 / 4.4 / 11.7 %; 0 / 0 / 4 % | 7.5 / 6.5 / 12.4 %; 12 / 6 / 12 % |
| expression shift, SCC / MS / PF | 4.3 / 3.8 / 2.5 %; 2 / 2 / 0 % | 2.1 / 1.9 / 1.2 %; 0 / 0 / 0 % | - |

- The group permutation test is calibrated for all three tests (its QQ curves follow the diagonal; the shift test,
  one-sided, is slightly conservative on PF), and after the kernel fix it is no longer anti-conservative for the
  shifts. The PF 1-hop excess is the sample-level factor above (three relabelings).
- The Cauchy combination of the permutation p-values is calibrated for the density tests but strongly conservative
  for the shifts (QQ flat at p = 1 for the upper 90 % of the groups): in sparse neighbourhoods the per-cell
  permutation p of the shift is tie-dominated and close to 1, and the Cauchy transform of many such p-values
  dominates the sum.
- The Cauchy combination of the parametric p-values inherits the parametric calibration: conservative for the heat
  kernel, anti-conservative for the 1-hop counts (7.5 % on SCC, a BH call in 12 % of the relabelings).

Run time per relabeling (all three tests, 1,999 permutations, unoptimized `-O0` build): median 63 s (SCC), 102 s
(MS), 174 s (PF).

## Results: semi-synthetic depletion

Figures: `figures/semi_recall.png` (recall per planted type, FDR 10 %), `figures/semi_fdp.png` (lenient false-discovery
proportion); tables `figures/semi_cells_summary.csv`, `figures/semi_groups_summary.csv`, `figures/semi_null_calls.csv`.
Recall = fraction of the remaining planted cells called "down" at FDR 10 %; FDP = fraction of the called cells that
belong to another type and are not graph neighbours of a planted cell. Values below are means over the three planted
types (small, medium, large) and the three replicate splits; the figures show the types separately.

**Unpaired regime** (4 + 5 or 5 + 5 whole control samples per arm): no method detects even the complete removal of a
type. Every cacoa method stays at recall 0-0.08 (the heat-kernel pooled FDR 0.00, the 1-hop pooled FDR 0.00-0.02), and
Milo reaches 0.49 only on PF at 100 % removal with a lenient FDP of 0.75 (its calls follow the processing site of the
PF controls: in two of the three control-only PF splits Milo calls 15-28 % of all cells, see below). With four or five
samples per arm the between-sample variability of real controls leaves no power for a single depleted type at the cell
level; the group-level tests behave the same (no planted type reaches BH significance, median group p 0.01-0.5).

**Paired regime** (each donor split into two halves, blocked by donor), recall at 100 / 80 / 50 % removal and the
lenient FDP at 100 % removal:

| method | SCC recall | SCC FDP | PF recall | PF FDP |
|---|---|---|---|---|
| Milo (graph, spatial FDR) | 0.94 / 0.87 / 0.26 | 0.32 | 0.83 / 0.69 / 0.15 | 0.25 |
| heat kernel, pooled permutation FDR (default) | 0.70 / 0.44 / 0.30 | 0.00 | 0.35 / 0.28 / 0.00 | 0.02 |
| heat kernel, max-t (former default) | 0.64 / 0.45 / 0.17 | 0.10 | 0.00 / 0.00 / 0.00 | - |
| heat kernel, BH on permutation p | 0.26 / 0.24 / 0.08 | 0.57 | 0.60 / 0.48 / 0.20 | 0.32 |
| heat kernel, BH on parametric p | 0.08 / 0.00 / 0.00 | 0.00 | 0.00 / 0.00 / 0.00 | - |
| 1-hop log counts, pooled permutation FDR | 0.34 / 0.06 / 0.00 | 0.04 | 0.19 / 0.05 / 0.00 | 0.05 |
| 1-hop log counts, max-t | 0.51 / 0.39 / 0.15 | 0.05 | 0.42 / 0.29 / 0.06 | 0.09 |
| 1-hop log counts, BH on permutation p | 0.05 / 0.02 / 0.00 | 0.97 | 0.12 / 0.08 / 0.03 | 0.97 |

- Milo has the highest recall but its cell-level FDP is 0.25-0.32 at a nominal 0.1: every cell inherits the smallest
  spatial FDR of the neighbourhoods containing it, so the called region spills over the planted type's neighbours
  (the lenient FDP already excludes the direct graph neighbours of planted cells).
- The heat-kernel pooled permutation FDR is the only cacoa method that combines useful recall with an FDP at or below
  the nominal level in both datasets (0.00 and 0.02 at 100 % removal; 0.001 and 0.34 at 80 % on SCC and PF). The
  max-t adjustment is competitive on SCC but calls nothing on PF (its familywise bound is set by the most variable
  cells of 89k), and BH on per-cell permutation p-values pays for the 1 / (B + 1) floor: with 512 (SCC) or 1,024 (PF)
  distinct relabelings the smallest attainable p is shared by thousands of weakly but consistently shifted cells, which
  BH then calls together (FDP 0.3-0.6 for the heat kernel, 0.97 for the 1-hop counts).
- BH on the parametric (moderated t) p-values calls almost nothing in the paired regime: the parametric test does not
  know the donor blocking, so the between-donor variance stays in its denominator.
- The 1-hop log-count unit (Milo's unit) is less powerful than the heat kernel under every adjustment: its per-sample
  counts are 0-5 cells in most neighbourhoods and the log transform leaves them discrete.
- Group-level summaries (`figures/semi_groups_summary.csv`): in the paired regime the group permutation test puts
  the planted type at its floor (p 0.003, 9 donors) in every set, but BH across 16-29 cell types turns 0.003 into
  0.05-0.09, so the type is formally called only for the 1-hop counts on PF (large and medium types); no control-only
  split produces a significant group by any of the three group tests (permutation, Cauchy of permutation p, Cauchy
  of parametric p).
- Control-only splits (`figures/semi_null_calls.csv`, FDR 10 %): every cacoa method calls at most 1 % of the cells
  (mean 0-0.7 %), and most call none; Milo calls 1 % (PF paired) to 15 % on average and 28 % at most (PF unpaired).

## Results: per-cell tail approximations

Question: can a per-cell p-value with resolution below the permutation floor 1 / (B + 1) be obtained from B
permutations, by a parametric tail or otherwise? Design (`tail_collect.R`, `tail_analyze.R`): one stratified null
relabeling of SCC (2,000 cells, 19,999 permutations; a second set of 1,100 cells with 49,999 permutations) and of PF
(2,000 cells, 19,999 permutations), with the full permutations x cells matrices kept for the three tests. For every
cell the first B.small permutations play the part of a routine run's permutation sample and 4,000-6,000 of the
remaining permutations are held-out null draws of the same cell's statistic. An approximation is calibrated when the
pooled held-out p-values are uniform; the table reports the observed over the expected fraction below alpha, pooled
over cells and draws (`figures/tail_calibration_*.csv`, `figures/tail_qq_scc.png`, `figures/tail_qq_pf.png`).

Approximations: the normal tail of the permutation-standardized z (the former `p.source = "z"`); a moment-matched
gamma on the statistic (F of the shift, t^2 of the density tests); a generalized Pareto tail on the top 10-25 % of
the permutations (probability-weighted moments and maximum likelihood); and the **pooled null**: every cell's
statistic standardized by the mean and sd of its own permutations, the standardized permuted statistics of all cells
pooled into one null sample, and the cell's standardized observed statistic ranked in it (`pooledPermutationP()`;
variants with leave-one-out standardization, median/IQR standardization, empirical-Bayes shrunk sd, strata by
sample count or by the density's effective sample number, and an oracle that standardizes with the moments of all
non-held permutations).

Ratio of observed to expected false-positive rate at alpha = 1e-3 / 1e-4 / 1e-5 (SCC, 49,999-permutation set; the
exact per-cell p is 0 below its floor):

| method | B = 1,999 | B = 4,999 | B = 19,999 |
|---|---|---|---|
| shift F, normal tail of z | 10 / 42 / 190 | - | - |
| shift F, gamma | 4.3 / 11 / 28 | - | - |
| shift F, pooled | 1.13 / 1.33 / 1.8 | 1.07 / 1.03 / 0.6 | 1.02 / 1.04 / 1.4 |
| heat kernel t, normal tail of z | 1.4 / 3.2 / 8.3 | - | - |
| heat kernel t, gamma | 1.1 / 1.8 / 2.9 | - | - |
| heat kernel t, pooled | 1.57 / 2.6 / 4.3 | 1.09 / 1.12 / 1.23 | 1.09 / 0.83 / 0.67 |
| 1-hop counts t, pooled | 1.01 / 1.02 / 1.12 | 1.00 / 1.04 / 1.09 | 1.00 / 1.04 / 1.06 |
| oracle (true moments, all permutations) | 1.0 / 1.0 / 0.6-1.1 | | |

(The normal and gamma rows are from the 2,000-cell, 19,999-permutation set, `tail_calibration_scc2000.csv`; the
generalized Pareto fits were either anti-conservative by 10-100x in the far tail or collapsed to 0 and are not shown.)
On PF with 22 samples the pooled p-values are calibrated or slightly conservative at B = 1,999 already (shift
1.01 / 0.82 / 0.74, heat 1.08 / 0.87 / 0.58, 1-hop 0.94 / 0.81 / 0.84); the normal tail is 11 / 48 / 240 for the
shift and the gamma 5.5 / 16 / 50.

- **No per-cell parametric tail works.** The normal tail of z overstates the evidence by one to two orders of
  magnitude at 1e-4 for the shift F (its null is right-skewed) and by 3-5x for the density t; the moment-matched gamma
  is still 10x off for the F at 1e-4 and 2-4x for the density tests; the generalized Pareto fits are unstable with
  the exceedance counts a routine run affords. The per-cell shape is not a parametric family: a cell's extreme values
  arise under the few relabelings that isolate its deviant samples, so the far tail is a finite mixture of such
  relabelings.
- **The pooled null is the construction that works**, and standardization is enough: with the true per-cell moments
  (oracle) the pooled null of all cells is calibrated for all three tests, even for the heat kernel whose per-cell
  nulls range from normal to a kurtosis of 50 (the kurtosis follows the effective number of samples carrying
  density at the cell, Spearman 0.84-0.91). Robust or shrunk standardizations and strata change little.
- **The far tail is governed by the number of relabelings, not by cells x permutations** (`tail_diag.R`): in a
  1,999-permutation pool of 2,000 SCC cells, 358 of the 410 values beyond the 1e-4 quantile come from 10 relabelings
  (those that isolate a deviant sample are extreme for many cells at once), and the calibration ratio at 1e-4 swings
  between 0.6 and 2.1 depending on which 1,999 relabelings are drawn (true moments throughout). Leave-one-out or
  split-half standardization cannot fix this; more relabelings do: 4,999 bring every test within 1.25x down to 1e-5,
  19,999 reproduce the oracle. PF (22 samples, a richer permutation space) is already there at 1,999.
- The 1-hop count statistic pools best (its discreteness decorrelates cells), the heat kernel worst (smooth densities
  make neighbouring cells extreme under the same relabelings).

Implementation: `pooledPermutationP()` (matrix or kernel-histogram input, chunked histogram of the pooled
standardized null); `p.pooled` is returned by both cluster-free maps and `p.source = "pooled"` feeds it to the
BH-type adjustments. The exact permutation p stays the default of `p.source`; the pooled FDR stays the default
adjustment (it uses the same pooled null for the expected exceedance count and shares its dependence on the number
of relabelings, but at the thresholds that matter for a 5-10 % FDR the null count is large and its noise small).

## Pilot: pooled p-values in the BH-type adjustments (2026-10-12)

`pilot_pooled.R` and `diag_pooled.R`; `$S/bench/pilot_pooled_B1999.log`. The battery in `common.R` now also feeds the
pooled p-values (`p.pooled`, see [pooledPermutationP()]) to BH, the spatial FDR, the neighbourhood BH and the Cauchy
summary (`bh.pooled`, `spatial.pooled`, `nbhd.pooled`, `p.cauchy.pooled`).

- **Null** (SCC relabeling 1 of `null3`, 1,999 permutations): no cell called by any method at FDR 10 % (one cell by
  BH on the normal tail of the shift z), no significant group by any summary; the smallest Cauchy-of-pooled-p group
  p is 0.09 (shifts), where the Cauchy of the exact permutation p is 0.11 and the group permutation test 0.13: the
  tie-at-1 conservativeness of the shift Cauchy is gone.
- **Paired semi-synthetic sets** (SCC, FDR 10 %): BH / spatial FDR / neighbourhood BH on the pooled p call 7,000-9,000
  cells at 100 % removal of `KC_Basal` (recall 0.68-0.78) where the pooled permutation FDR calls 1,167 (recall 0.66),
  and at 80 % removal of `LC` 1,000-1,500 (recall 0.6-0.7) against 569 (recall 0.72); at 100 % removal of the small
  `Tcell` type they call 0-37 cells (recall 0) against 180 (recall 0.70). For the 1-hop counts the pooled p variants
  have recall 0.14-0.57 with lenient FDP 0-0.5, the pooled FDR 0.14-0.45 with FDP 0.03-0.07; Milo 0.86-0.97.
- **Diagnosis** (`diag_pooled.R large_f100_r1`): the extra calls are the *other* cell types going up: removing 18 % of
  the cells of one arm raises the relative density of every other cell by 22 %, and 85 % of the called other cells
  have a positive effect (all 5,962 called by the standardized pooling). The exact per-cell permutation p of the other
  cells is at its floor for a quarter of them, so the per-cell tests see the same thing. The pooled permutation FDR
  misses them because it pools the *raw* studentized contrasts: the benchmark's paired sets were fitted with `~ cond`
  and the donor only as a permutation block, so the between-donor variance stays in the residual and the null
  permutation sd of t is 0.32 for the unaffected cells but 1.19 for the planted ones (the effect flips sign under
  the within-donor swaps); the raw pool's tail is then set by the planted cells and is far too wide for the others.
  Per-cell standardization removes the scale difference (and costs the planted cells of the small type their power,
  since their own permutation sd is inflated by the effect). With the donor as a model covariate
  (`~ patient + cond`, `block.vars = "patient"`) the null sd of t is 1.11 for every class and all four constructions
  agree: pooled FDR on raw t 7,522 calls / recall 0.68, on standardized t 7,749 / 0.71, BH on the pooled p 7,522 /
  0.68 (raw) and 7,737 / 0.71 (standardized), each with about 6,000 "other" cells called, all up.
- Milo's extra calls in the same set are also mostly up (946 of 1,290 at 100 % removal; 0 of 115 at 80 % removal of
  `LC`, 1 of 1 for `Tcell`), so its earlier "FDP 0.25-0.44" is largely the compositional response as well.
- In the null collection (`tail/scc_r1_B19999.rds`, unpaired model blocked on the true condition) the per-cell
  permutation sd of the heat-kernel t ranges from 0.76 to 1.12 (1st-99th percentile), of the 1-hop t 0.62-1.12 and
  of the shift F 0.05-0.65: the raw pooling of the shift F would not be calibrated, which is why the shift kernel
  always pooled the standardized z; for the density t the raw pooling is approximately calibrated only when the model
  absorbs the blocking variance.

Consequences for the benchmark design: (i) the truth of the semi-synthetic sets must account for compositional closure
(cells of other types are expected to go up by `1 / (1 - f * pi)`, not to stay null); (ii) paired sets must be fitted
with the donor as a covariate, for cacoa as for Milo; (iii) the pooled FDR and the pooled p-values should both be
computed from the per-cell standardized statistics (the shifts already are), or the density statistic must be made
pivotal by the model; the trade-off is the power loss for affected cells whose own permutation distribution is
inflated by the effect, to be measured.

## Results: planted-relabeling power study (2026-10-12)

Design (`planted_relabel.R`): the balanced relabelings of the null study (SCC 9 + 9 and PF 11 + 11 whole samples, the
real condition split evenly between the arms, permutations blocked on it) with 50, 80 or 100 % of one cell type's
cells removed from arm B, three relabelings per setting, plus three control-only relabelings; planted types SCC MDSC
(2 %), KC_Cyc (7 %), CD1C (8 % of arm B); PF MUC5B+ (2 %), Endothelial (7 %), Ciliated (20 %). Milo on the same sets
with `~ condition`. Truth with compositional closure: the planted type's remaining cells should be called down; every
other cell is expected *up* by the closure factor (1.03-1.27) and such calls are reported separately (`other.up`); a
cell of another type called down is a false discovery (`fdp.down`); graph neighbours of planted cells are ambiguous
and excluded. Figures `figures/planted_recall.png`, `planted_fdp.png`, `planted_up.png`; tables
`planted_cells_summary.csv`, `planted_decomp_summary.csv`, `planted_null_calls.csv`, `planted_groups.csv`,
`planted_scales.csv`.

Recall of the planted cells at FDR 10 % (false-discovery proportion in brackets; mean over three relabelings):

| method | SCC medium 80 / 100 % | SCC large 100 % | PF medium 80 / 100 % | PF large 100 % | small types |
|---|---|---|---|---|---|
| heat kernel, pooled permutation FDR (default) | 0.19 / 0.28 (0) | 0.14 (0) | 0.20 / 0.47 (0) | 0.28 (0.0002) | 0 |
| heat kernel, pooled FDR on standardized t | 0.19 / 0.28 (0) | 0.13 (0) | 0.22 / 0.48 (0) | 0.28 (0.0002) | 0 |
| heat kernel, pooled FDR, ter Braak scale | 0.19 / 0.38 (0.001) | 0.18 (0.0002) | 0.21 / 0.47 (0) | 0.25 (0) | 0 |
| heat kernel, BH / spatial FDR / neighbourhood BH on the pooled p | 0.19 / 0.29, 0.18 / 0.27, 0.23 / 0.36 | 0.13, 0, 0 | 0.22 / 0.49, 0.10 / 0.45, 0.13 / 0.55 | 0.28, 0.20, 0.29 | 0 |
| 1-hop counts (log), pooled FDR | 0 / 0.07 (0) | 0.06 (0.004) | 0 / 0.20 (0.006) | 0.12 (0.003) | 0 |
| Milo (edgeR QL, spatial FDR with inheritance) | 0.29 / 0.83 (0.27) | 0.72 (0.19) | 0.92 / 0.99 (0.47-0.52) | 0.99 (0.48) | 0.70-0.77 (0.33-0.52) |

- **Standardization costs no power here.** The permutation sd of the statistic of the planted cells is within 1-8 % of
  the other cells' (`planted_scales.csv`: heat kernel 1.06 vs 1.02 on SCC, 1.01 vs 1.00 on PF), so the raw and the
  per-cell standardized pooled FDR coincide; the ter Braak (effect-free) scale adds a little on SCC (0.38 vs 0.28 at
  the medium type). The sixfold scale difference of the paired pilot was a property of fitting the donor-split sets
  without the donor covariate, not of standardization.
- **The pooled p-values do not change the BH-type adjustments.** BH on the pooled p equals the pooled FDR; the spatial
  FDR is slightly less powerful and the neighbourhood BH slightly more (0.36 vs 0.28), all without false calls.
- **cacoa's density map is underpowered in this realistic unpaired design**: a third to a half of a medium or large
  type's cells at complete removal, nothing of a 2 % type, with zero false calls; the 1-hop count unit is weaker still.
- **Milo is not calibrated on these sets.** Its false-discovery proportion is 0.2-0.5 at a nominal 0.1, and in the
  control-only relabelings it calls 10-13 % (SCC) and 26-32 % (PF) of all cells, where every cacoa method calls at
  most 0.001 %. Its recall is therefore not comparable at face value; see the decomposition below.
- **Compositional response.** The heat-kernel map calls 4-6 % of the other cells up at 100 % removal of the SCC large
  type (closure 1.08) and essentially none elsewhere; Milo calls 8-13 % up in every set, the null sets included.
- **Group level** (`planted_groups.csv`): the planted type is the most significant group in every 80-100 % set of a
  medium or large type under all three summaries (group permutation p 0.0005-0.04, Cauchy of pooled p down to 1e-7),
  and no control-only relabeling has a significant group. BH across 19-29 types reaches 0.05 only for the Cauchy of
  the pooled p (which also flags the types that rise by closure); the group permutation p is floored at 1 / 2000.

### Decomposition of the Milo comparison (`milo_decompose.R`, `figures/planted_decomp_summary.csv`)

The Milo pipeline differs from cacoa in the unit (its refined kNN neighbourhoods, ~95 cells on SCC and ~116 on PF,
one per ~12 cells, against cacoa's per-cell 1-hop sets of ~59 / 40 cells or the heat kernel), the test (edgeR's
quasi-likelihood F on counts against cacoa's permutation linear model on log counts) and the adjustment (spatial FDR
with every cell inheriting the smallest FDR of the neighbourhoods containing it). Every combination was run on the
same planted sets: Milo's neighbourhood matrix from its own run, with cacoa's tests applied to its per-sample counts
and mapped back to cells by the same inheritance; edgeR and limma-voom applied to cacoa's 1-hop and 2-hop counts;
and a closed-form weighted-least-squares permutation test (`permWLS()` in `common.R`, verified against `lm()`) with
voom precision weights, the prototype of observation weights in cacoa's engine. Library-size normalization is
reported next to Milo's TMM, which assumes that most units do not change and is violated by a large depletion plus
the compositional rise of everything else (PF closure 1.27). Recall at 100 % removal, mean over the three planted
types, with the false-discovery proportion among other cells called down (`fdp`) and the same beyond two graph
hops of any planted cell (`fdp2`, i.e. excluding the spill-over radius of a ~100-cell neighbourhood with inheritance):

| unit, test, adjustment | SCC recall (fdp / fdp2) | PF recall (fdp / fdp2) |
|---|---|---|
| Milo: refined kNN, edgeR QL (TMM), spatial FDR + inheritance | 0.75 (0.26 / 0.19) | 0.92 (0.51 / 0.46) |
| refined kNN, edgeR QL (library size), BH + inheritance | 0.51 (0.16 / 0.11) | 0.81 (0.46 / 0.43) |
| refined kNN, voom + moderated t, BH + inheritance | 0.42 (0.009 / 0) | 0.56 (0.23 / 0.18) |
| refined kNN, voom-weighted permutation, pooled FDR (standardized) + inheritance | 0.46 (0.014 / 0.001) | 0.69 (0.30 / 0.25) |
| refined kNN, voom-weighted permutation, pooled FDR (raw) + inheritance | 0.42 (0.012 / 0) | 0.56 (0.15 / 0.12) |
| refined kNN, unweighted permutation on the same log-CPM, pooled FDR + inheritance | 0.30 (0.03 / 0) | 0.51 (0.11 / 0.09) |
| refined kNN, cacoa's permutation LM (moderated), pooled FDR + inheritance | 0.21 (0.019 / 0.005) | 0.51 (0.17 / 0.16) |
| 1-hop counts, edgeR QL, BH | 0.16 (0.09 / 0.055) | 0.75 (0.51 / 0.50) |
| 1-hop counts, voom-weighted permutation, pooled FDR (standardized) | 0.12 (0.007 / 0.001) | 0.24 (0.13 / 0.10) |
| 1-hop counts, voom-weighted permutation, pooled FDR (raw) | 0.08 (0.011 / 0) | 0.19 (0.05 / 0.05) |
| 1-hop counts, voom-weighted permutation + inheritance | 0.55 (0.054 / 0.008) | 0.61 (0.32 / 0.25) |
| 2-hop counts (~1,200 cells), voom-weighted permutation, pooled FDR | 0.17 (0.007 / 0) | 0.57 (0.025 / 0.014) |
| heat kernel, pooled FDR (cacoa default) | 0.14 (0 / 0) | 0.25 (0 / 0) |
| heat kernel, pooled FDR + inheritance | 0.27 (0.04 / 0) | 0.54 (0.006 / 0.001) |
| 1-hop counts (log), cacoa's permutation LM, pooled FDR | 0.04 (0.002 / 0) | 0.10 (0.004 / 0.003) |

Control-only relabelings (`planted_null_calls.csv`): every edgeR-based row calls 1-3 % of the cells on SCC and
11-26 % on PF (Milo 10 % and 26 %); every permutation-based row calls at most 0.1 % on SCC and at most 1.4 % on PF
(the standardized weighted variant on the kNN unit; its raw-pooling counterpart 0.03 %).

- **The test.** On the same unit and adjustment, precision weights double the permutation test's recall on SCC at
  calibrated error (0.30 unweighted, 0.42-0.46 weighted, fdp2 0); cacoa's current log-count LM is at 0.21. voom with
  BH reaches the same 0.42, so the count model's advantage over cacoa is the weighting, not the likelihood. edgeR
  itself is not calibrated on these data (fdp2 0.11 on SCC, 0.43 on PF; 11-26 % of the cells called in PF's
  control-only relabelings): its parametric error model does not carry PF's two processing sites, the permutation
  framework does.
- **Normalization.** TMM costs recall and adds false calls under a large planted depletion (kNN + edgeR 0.40 -> 0.51
  with library size on SCC; the voom-weighted permutation 0.42 -> 0.46); library-size normalization is the right
  default for neighbourhood counts.
- **The unit and the inheritance.** The kNN unit with inheritance is worth about a factor of two on SCC (1-hop
  weighted 0.12 -> kNN weighted 0.46) and inheritance alone, applied to cacoa's own per-cell results, gives most of
  it (1-hop weighted + inheritance 0.55 with fdp2 0.008; heat kernel + inheritance 0.27 / 0.54 with fdp2 0 / 0.001).
  On PF every row built on the kNN unit shows 0.1-0.25 false calls beyond the spill-over radius once an effect is
  planted, although they are calibrated in the control-only sets; the per-cell units (heat kernel, 1-hop, 2-hop) do
  not. So the kNN unit's extra power on PF is not trustworthy, and the heat kernel with inheritance is the best
  calibrated combination there (0.54 at fdp2 0.001).
- **Standardization under counts.** On PF the standardized pooled FDR of the count units is somewhat more liberal than
  the raw pooling (fdp2 0.10 vs 0.05 for 1-hop weighted; 1.4 % vs 0.03 % calls in the control-only sets); for the heat
  kernel the two coincide. Raw pooling of a pivotal weighted statistic is the safer default for count units.

## Results: null round with the pooled variants (SCC, 15 relabelings, 4,999 permutations; `figures/null4/`)

Every cacoa adjustment, the pooled-p ones included (`BH / spatial FDR / neighbourhood on the pooled p`, the raw /
standardized / ter Braak pooled FDRs), calls at most 0.1 % of the cells at FDR 5 % in any relabeling (1-hop counts:
one relabeling fires at 0.1-1.3 %; the heat kernel and the shifts: zero). QQ: the pooled p is uniform to 1e-5 for the
1-hop counts and the shifts, conservative beyond 1e-3 for the heat kernel (the parametric t more so). Group level:
the Cauchy of the pooled p is calibrated for the density maps (2.5-3.9 % of the groups at 5 %) but anti-conservative
for the shifts (7.7 % at 5 %, 2 of 15 relabelings with a BH call), so the group permutation test stays the default
summary and the Cauchy of the pooled p is an option for density only.

## Conclusions (2026-10-12 update)

8. **Power on real data, unpaired.** With 9-11 whole samples per arm cacoa's density map recovers a third to a half of
   a medium or large type's cells at complete removal and nothing of a 2 % type, with no false calls; standardizing
   the statistic per cell before pooling costs no power; the pooled p-values leave the BH-type adjustments unchanged.
9. **Milo's advantage decomposed.** Half of it is precision weights (voom-type) in the test, reproducible inside
   cacoa's permutation framework at calibrated error; the rest is the neighbourhood inheritance, which is cheap and
   calibrated on cacoa's own per-cell units; Milo's own adjustment is not calibrated (fdp 0.2-0.5, 10-30 % of the
   cells called in control-only relabelings), and neither is edgeR on PF. Library-size normalization beats TMM for
   neighbourhood counts.
10. **Package work this supports:** observation weights in the permutation LM; a kNN neighbourhood-count density with a
    target size and a label-free k diagnostic (cells per sample, purity, minimum detectable fold change); an optional
    inheritance step for the per-cell map; the Cauchy of the pooled p as a density-only group summary.

## Conclusions

1. **Per-cell calling: pooled permutation FDR** (`adjust.method = "perm-fdr"`, now the default of both maps). It is
   calibrated under relabeling (0-0.03 % false calls at FDR 5 % with the heat kernel), its false-discovery proportion
   in the semi-synthetic depletions stays at or below the nominal level (0.00-0.02 at full removal), and it has the
   best recall among the cacoa adjustments on both datasets. BH on per-cell permutation p-values is also calibrated
   under the null but pays for the 1 / (B + 1) floor with FDP 0.3-0.6 when there is signal; max-t is calibrated but
   calls nothing on 89k cells; BH on parametric p-values ignores the blocking and calls nothing in paired designs.
2. **The z scale is for plotting only.** The normal tail of the permutation-standardized z is anti-conservative at
   every quantile and BH on it calls cells in practically every null relabeling; the parametric F of the
   neighbourhood shift is not calibrated either, and no per-cell parametric tail (gamma, generalized Pareto) is.
   Per-cell p-values of the shifts come from the permutations: exact (floored at 1 / (B + 1)) or pooled over cells
   (`p.source = "pooled"`, no floor), the latter calibrated to 1e-5 with about 5,000 permutations on 18 samples and
   with 2,000 on 22 samples; the far tail of the pooled null is set by the number of relabelings, so more permutations
   are the only way to sharpen it.
3. **Per-type summaries: the group permutation test** (default of `estimateClusterFreeSummary()`). It is calibrated for
   the density and the shift maps and carries the dependence between cells exactly through the shared permutations.
   The Cauchy combination is kept as the alternative (`method = "cauchy"`): calibrated on the density maps, strongly
   conservative on the shifts, and on parametric p-values it inherits their miscalibration. With BH across 20-30 cell
   types the floor 1 / (B + 1) of the group test matters: 1,999 permutations (floor 0.0005) leave room for the
   adjustment, 499 (floor 0.002) barely.
4. **Versus Milo.** Milo detects more cells in the paired depletions (recall 0.83-0.94 against 0.35-0.70 for the heat
   kernel with pooled FDR) but a quarter to a third of its calls are false at a nominal 0.1, since every cell
   inherits the best spatial FDR of its neighbourhoods, and it calls 15-28 % of all cells in two of three PF
   control-only splits. No method, Milo included, detects even a complete depletion with 4-5 whole samples per arm;
   the group-level tests on such designs reach their permutation floor but not BH significance across 20-30 types.
5. **Heat kernel over 1-hop counts.** The heat-kernel density is the more powerful unit under every adjustment, and
   the 1-hop count statistic carries a sample-level factor (total cells, connectivity) that makes whole relabelings
   fire in PF under every adjustment, parametric ones included.
6. **Kernel fix.** Pooled statistics (group means, max-t) require the shared permutations to be applied as global
   relabelings in every cell; inducing them within the present samples of a cell is exact per cell but incoherent
   across cells. The fixed kernel reports the admissible permutation count per cell (`n.perm`). The Freedman-Lane
   residual scheme keeps the within-cell permutation, so its max-t across cell types (`testPairwiseEffects()`,
   `testTermEffects()`) carries the same caveat.
7. **Design of null studies.** A relabeling null must be exchangeable with the model's permutations: balanced
   relabeling within the true condition requires `block.vars` on that stratum, otherwise every permutation test is
   conservative (round 2), and unrestricted relabeling carries real composition differences into the "null" (round 1).
