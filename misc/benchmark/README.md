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
| `$S/bench/run_milo_semi.R` | Milo (miloR 1.6.0) on the semi-synthetic sets (subagent) |

Results: [null calibration](#results-null-by-relabeling), [semi-synthetic power and false discoveries](#results-semi-synthetic-depletion),
[conclusions](#conclusions).

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

## Conclusions

1. **Per-cell calling: pooled permutation FDR** (`adjust.method = "perm-fdr"`, now the default of both maps). It is
   calibrated under relabeling (0-0.03 % false calls at FDR 5 % with the heat kernel), its false-discovery proportion
   in the semi-synthetic depletions stays at or below the nominal level (0.00-0.02 at full removal), and it has the
   best recall among the cacoa adjustments on both datasets. BH on per-cell permutation p-values is also calibrated
   under the null but pays for the 1 / (B + 1) floor with FDP 0.3-0.6 when there is signal; max-t is calibrated but
   calls nothing on 89k cells; BH on parametric p-values ignores the blocking and calls nothing in paired designs.
2. **The z scale is for plotting only.** The normal tail of the permutation-standardized z is anti-conservative at
   every quantile and BH on it calls cells in practically every null relabeling; the parametric F of the
   neighbourhood shift is not calibrated either. Per-cell p-values of the shifts come from the permutations.
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
