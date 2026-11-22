## Cluster-free expression shifts on the pairwise-effects engine (Track D.3).
##
## For every cell, the samples' mean expression over the cell's neighbourhood gives a sample x sample distance
## matrix (computed in C++ by estimateExpressionShiftsPairsLM over all sample pairs). The contrast's shift is
## estimated per neighbourhood with the individual-level model (common dispersion), tested by permuting sample
## labels (one set of permutations shared by all cells, induced on each neighbourhood's sample subset), and the
## per-cell z-scores are adjusted by the max-statistic over cells (adjustedZScoresMaxStat) and smoothed over
## the graph (applyMedianFilterES).

# square sample distance matrix from one column of the pair layout. Samples with missing pair distances (too few
# cells in the neighbourhood) are dropped one at a time, the one with the most missing pairs first, until no missing
# pair is left: a thin sample removes itself, not the samples it is paired with (same rule as the C++ kernel).
pairVectorToMatrix <- function(y, pairs, samples) {
  n <- length(samples); D <- matrix(NA_real_, n, n, dimnames = list(samples, samples)); diag(D) <- 0
  D[pairs] <- y; D[pairs[, 2:1, drop = FALSE]] <- y
  ok <- rep(TRUE, n)
  repeat {
    miss <- rowSums(is.na(D[ok, ok, drop = FALSE]))
    if (!length(miss) || max(miss) == 0) break
    ok[which(ok)[which.max(miss)]] <- FALSE
  }
  D[ok, ok, drop = FALSE]
}

# gsd under residual exchange: the scale of the relabeled statistics (residualScale()) as the kernel's code: 0 average, 1 relabeling,
# 2 observed (as "relabeling", scaled to the observed arrangement; gsd only), 3 independent
gsdResidualScale <- function(crosses.levels = FALSE) as.integer(match(residualScale(crosses.levels), c("average", "relabeling", "observed", "independent")) - 1L)

#' Cluster-free expression shifts for one contrast
#'
#' @param cm genes x cells sparse count matrix
#' @param sample.per.cell factor of sample per cell (named by cell, aligned to `colnames(cm)`)
#' @param nns.per.cell list (per cell) of neighbour indices (0-based, as returned by the graph adjacency)
#' @param design sample-level design of a contrast test (`buildDesignMatrices()` output or a `cacoaModel`)
#' @param meta sample metadata
#' @param statistic test of each neighbourhood: `"gsd"` (default) tests every gene of the neighbourhood and combines them: per gene the contrast
#'   estimate over its variance (the noise model's prediction moderated towards the residual variance of each relabeled fit),
#'   combined over genes as the maximum over gamma = 1, 2, 4, 8 and the max of the standardized log sum of |z|^gamma
#'   (adaptive: many small or few large changes); `"distance"`: the sample distance matrix of the neighbourhood's profiles tested
#'   with a pseudo-F (`dist`, the previous method). gsd needs `profile = "normalized"` and `weighting = "noise"`, and handles any
#'   single-df contrast (two groups, a numeric variable, covariates; block relabeling or Freedman-Lane)
#' @param dist distance between two samples' neighbourhood profiles (`statistic = "distance"`). With `profile = "normalized"`: `"cor"` (default;
#'   1 - the correlation of the gene-centred, noise-weighted profiles: genes are centred across the unit's samples and
#'   scaled by sqrt of their weight) or `"euclidean"` (the weighted squared Euclidean distance). With `profile = "raw"`:
#'   `"cor"` (1 - Pearson correlation of the two samples' profiles over genes, without centring genes across samples),
#'   `"cosine"` or `"js"`
#' @param profile `"normalized"` (default): per neighbourhood and sample, the summed counts over the summed library sizes of the
#'   sample's cells (`lib`), transformed as log10(scale p + 1), with noise weights (`weighting`); `"raw"`: the mean raw
#'   count per cell, log10(1e3 x + 1) (`log.vectors`), unweighted (the previous default)
#' @param lib library size per cell over all genes (named by cell or aligned to `colnames(cm)`); needed for
#'   `profile = "normalized"` (default: the column sums of `cm`, i.e. over the profile genes only)
#' @param profile.scale transform scale of the normalized profiles: `"median"` (default; the median library size, so that the
#'   pseudocount is about one count per typical cell) or a number
#' @param weighting `"noise"` (default): gene weights 1 / (predicted sampling variance + between-sample floor) from a
#'   label-free noise model (Poisson variance of the transform plus within-unit heterogeneity, see
#'   [fitShiftNoiseCurves()]); `"none"`: unweighted
#' @param sample.floor include the between-sample variance floor tau2(mu) in the gene weights (default TRUE)
#' @param noise.correction subtract each sample's predicted noise: for `"cor"` the squared norms of the weighted
#'   profiles are disattenuated (|z|^2 - c, floored at 5 % of |z|^2), for `"euclidean"` the distance is reduced by
#'   c_a + c_b. Default FALSE: in thin units (few cells per sample) the predicted noise approaches the profile norm and the
#'   corrected correlations clamp at +-1 (the test then degenerates); on the benchmark the correction added about 0.01 recall
#'   for `"cor"`, while for `"euclidean"` it is sensitive to a misfitted noise model
#' @param noise.curves fitted noise curves ([fitShiftNoiseCurves()] output with the same `scale`); default: fitted on
#'   `n.fit.units` random units of the data
#' @param n.fit.units units used to fit the noise curves (default 200)
#' @param min.n.obs.per.samp minimum cells of a sample inside the neighbourhood (default 3)
#' @param min.samp.per.level minimum samples per compared level (default 2)
#' @param permutation,n.permutations,block.vars,seed permutation settings (see [permutationPlan()]); `block.vars` defaults to the model's
#' @param max.enumerate enumerate all relabelings when there are at most this many (and at most `n.permutations`); the permutation p is
#'   then exact (the identity included, no +1). Default 0: always random draws (the cluster-free maps)
#' @param adjust max-statistic adjustment of the z-scores across cells (default TRUE)
#' @param smooth median-filter the shifts and adjusted z-scores over the graph (default TRUE)
#' @param wins winsorizing fraction for the adjustment (default 0.025)
#' @param adjust.method `"bh"` (the default with gsd: BH on `p.source`, by default the pooled p), `"perm-fdr"` (the default
#'   with the distance map: pooled permutation FDR over cells, see [permutationFDR()]; no per-cell floor), `"spatial-fdr"`,
#'   `"neighbourhood"` or `"max-t"`, see [adjustCellScores()]
#' @param p.source per-cell p-values for the BH-type adjustments (default `"pooled"` with gsd, `"permutation"` with the
#'   distance map): `"permutation"` (calibrated but floored at
#'   1 / (n.permutations + 1), so with many cells and sparse signal BH may call nothing), `"pooled"` (the cell's
#'   standardized z against the pooled null of all cells' permuted z-scores, see [pooledPermutationP()]: no floor,
#'   calibrated with enough permutations; about 2,000 keep the far tail within a factor of 2, 10,000 make it exact),
#'   `"z"` (the normal tail of the permutation-standardized z-score: no floor, but anti-conservative in the tail by
#'   orders of magnitude) or `"parametric"` (F statistic; not calibrated for this statistic)
#' @param z.grid thresholds of the permutation-standardized z at which the pooled null exceedance counts are
#'   accumulated for `"perm-fdr"` and the pooled p-values (default -6 to 25 in steps of 0.02)
#' @param groups optional factor of cell groups (named by cell): the kernel then also accumulates the sum of the
#'   permuted z-scores per group, from which `group.test` gives the exact group-level permutation test (the mean z of
#'   the group's cells against its permutation distribution; one-sided) with BH across groups; the permutation
#'   distribution of the group means is stored as `attr(group.test, "perm.means")`
#' @param return.perm return the permutations x cells matrix of the permuted statistics as `stats.perm` (default
#'   FALSE; `n.permutations` x cells doubles, NaN where a cell is not tested under a relabeling)
#' @param post post-processing on the cell graph (default `"averaged"` with gsd, `"none"` with the distance map): `"none"`, `"averaged"` (the shift statistic averaged over
#'   each cell and its graph neighbours, under every permutation alike, before the pooled FDR: an exact test of the
#'   one-hop-larger unit with calibrated q-values, [averageOverNeighbours()]; needs `adjust.method = "perm-fdr"` (or `"bh"` with the pooled p) and
#'   keeps the permutation matrix in memory) or `"inherit"` (every cell takes the smallest q among itself and its
#'   graph neighbours after the test, [inheritNeighbourQ()])
#' @param pooled.statistic statistic pooled over cells for the pooled p-values and FDR: `"standardized"` (default) or `"raw"` (the shift F itself)
#' @param vectorize permutation loop of the complete-sample cells as matrix products (default TRUE; FALSE for the reference loop)
#' @param nns.graph the 1-hop graph neighbourhoods used by `post` (list of 0-based indices per cell, the cell
#'   itself included; default: `nns.per.cell`, which is the 1-hop graph unless a kNN unit was passed)
#' @details Label permutations (the block scheme) are applied as global relabelings of the samples: in a neighbourhood
#'   where some samples are absent, a present sample takes the design row of the sample it is mapped to, so that one
#'   permutation means the same relabeling in every cell and statistics pooled over cells (the group test, the
#'   pooled permutation FDR, the max statistic) are compared with a coherent null. A cell is not tested under a
#'   relabeling that leaves fewer than `min.samp.per.level` present samples at a compared level (`n.perm` counts the
#'   relabelings under which each cell was tested). Residual permutations (Freedman-Lane) act on the present samples'
#'   residuals and are induced within each cell's sample set.
#' @param log.vectors log10(1e3 x + 1) transform of the mean profiles with `profile = "raw"` (default TRUE)
#' @param n.cores threads for the per-cell loop (C++)
#' @param verbose progress
#' @param robust `"none"`, `"huber"` or `"winsor"`. With `statistic = "distance"`: down-weight samples with outlying residual
#'   distances (re-estimated under every relabeling). With `statistic = "gsd"`: `"winsor"` winsorizes each gene's sample values
#'   in the neighbourhood before the test (label-free: at the median +- `winsor.k` robust sd, the sd being the larger of the sample's
#'   predicted noise and the MAD-based spread), `"huber"` is not available
#' @param robust.k robust tuning constant (`statistic = "distance"`)
#' @param neighbourhood.gene.frac both statistics: in each cell's neighbourhood, use the genes of `cm` detected in at least this
#'   fraction of the neighbourhood's cells (pooled over samples, label-free; default 0.05), so that genes specific to a population
#'   are analyzed where they are expressed (the rule of the per-cell-type analyses); 0 uses every gene of `cm` present in the
#'   neighbourhood
#' @param centres tested centres given explicitly (1-based indices of `nns.per.cell`; gsd with `post = "averaged"`): only their
#'   neighbourhoods are tested and each one's statistic is scattered to its cells, as with `centres.per.cell`, whose greedy
#'   cover they replace (default NULL); the neighbourhoods of the other cells may be empty ([diffusionNeighbourhoodLists()])
#' @param centres.per.cell gsd with `post = "averaged"`: test the neighbourhoods of only a subset of cells (the centres), chosen
#'   so that each cell has about `centres.per.cell` centres among its graph neighbours (an integer; default 2, about 12-14 % of
#'   the cells). Each tested neighbourhood's statistic, observed and permuted alike, is passed to all cells of the
#'   neighbourhood, and a cell's statistic is the mean over the tested neighbourhoods that contain it (typically 20-30), so
#'   every cell is still tested exactly. Larger values test more neighbourhoods (slower, marginally smoother); `Inf` or `NULL`
#'   tests every cell's neighbourhood (the full map, for example for very small populations)
#' @param gsd.variance gsd per-gene variance: `"moderated"` (default) or `"model"` (the noise model only)
#' @param gsd.gamma gsd combination over genes: 0 (default, adaptive) or one of 1, 2, 4, 8, 99 (max) alone. Under residual exchange (a continuous covariate
#'   without a grouping term) the variance of the contrast and the expected residual sum of squares differ between relabelings; gsd takes
#'   both from its noise model for every relabeling and gene and removes the dependence between the relabeled contrast and the residuals of
#'   the relabeled fit, so that the relabeled statistics are on the scale of the observed one (`options(cacoa.residual.scale = )`:
#'   `"independent"`, the default; `"relabeling"` without the last step, more sensitive but above the nominal level in some designs;
#'   `"average"` the former single average factor). The standardizations over relabelings inside the statistic (of the combinations over
#'   genes, and of a neighborhood's statistic) include the observed data set, so that it is treated exactly as the relabeled ones
#'   (`options(cacoa.gsd.symmetric = FALSE)`: over the relabeled data sets only, the behavior before this option existed)
#' @param gsd.d0 prior degrees of freedom of the model variance in the moderated variance (default 3)
#' @param winsor.k winsorization threshold of `robust = "winsor"` with gsd, in robust sd units (default 3)
#' @return list: `stat` (the neighbourhood statistic per cell: the shift F, or the gsd combination), `p.value` (exact permutation p), `p.pooled` (pooled-null p, no floor), `p.param` (NA for gsd),
#'   `z.pooled` (the tail z of `p.pooled`, `qnorm(1 - p)`), `z.std` (the permutation-standardized
#'   statistic: observed minus permutation mean over permutation sd; the input of the max-t adjustment and of the
#'   pooled null), `z.adj`, `shifts` (distance map: the bias-corrected shift estimate; gsd: the mean squared per-gene z of
#'   the neighbourhood, about 1 without a change), `shifts.smoothed`, `n.samples`, `settings`
#' @export
clusterFreeExpressionShifts <- function(cm, sample.per.cell, nns.per.cell, design, meta, statistic = c("gsd", "distance"), dist = c("cor", "euclidean", "cosine", "js"),
                                        profile = c("normalized", "raw"), lib = NULL, profile.scale = "median", weighting = c("noise", "none"),
                                        sample.floor = TRUE, noise.correction = FALSE, noise.curves = NULL, n.fit.units = 200, min.n.obs.per.samp = 3,
                                        min.samp.per.level = 2, permutation = "auto", n.permutations = 199, block.vars = NULL, seed = 1, max.enumerate = 0,
                                        adjust = TRUE, smooth = TRUE, wins = 0.025, log.vectors = TRUE, n.cores = 1, verbose = FALSE,
                                        robust = c("none", "huber", "winsor"), robust.k = 1.345,
                                        adjust.method = NULL, p.source = NULL,
                                        z.grid = seq(-6, 25, by = 0.02), groups = NULL, return.perm = FALSE, post = NULL, nns.graph = NULL, vectorize = TRUE, pooled.statistic = c("standardized", "raw"),
                                        neighbourhood.gene.frac = 0.05, centres.per.cell = 2, centres = NULL, gsd.variance = c("moderated", "model"), gsd.gamma = 0, gsd.d0 = 3, winsor.k = 3) {
  statistic <- match.arg(statistic); gsd <- statistic == "gsd"
  dist <- match.arg(dist); robust <- match.arg(robust); pooled.statistic <- match.arg(pooled.statistic); gsd.variance <- match.arg(gsd.variance)
  # statistic-dependent defaults: gsd = per-cell calls by BH on the pooled p of the neighbour-averaged statistic; the distance map keeps its previous ones
  adjust.method <- match.arg(adjust.method %||% (if (gsd) "bh" else "perm-fdr"), c("perm-fdr", "bh", "spatial-fdr", "neighbourhood", "max-t"))
  p.source <- match.arg(p.source %||% (if (gsd) "pooled" else "permutation"), c("permutation", "pooled", "z", "parametric"))
  post <- match.arg(post %||% (if (gsd) "averaged" else "none"), c("none", "averaged", "inherit"))
  profile <- match.arg(profile); weighting <- match.arg(weighting)
  if (gsd) {
    if (profile != "normalized" || weighting != "noise") stop("statistic = 'gsd' needs profile = 'normalized' and weighting = 'noise'")
    if (robust == "huber") stop("robust = 'huber' is available for statistic = 'distance'; statistic = 'gsd' takes robust = 'winsor' (label-free winsorization)")
    if (!is.numeric(neighbourhood.gene.frac) || length(neighbourhood.gene.frac) != 1 || !(neighbourhood.gene.frac >= 0 && neighbourhood.gene.frac < 1)) stop("neighbourhood.gene.frac must be in [0, 1)")
    if (!is.null(centres.per.cell) && (length(centres.per.cell) != 1 || !is.numeric(centres.per.cell) || !(centres.per.cell >= 1) || (is.finite(centres.per.cell) && centres.per.cell != round(centres.per.cell))))
      stop("centres.per.cell must be a positive integer (tested neighborhood centers per cell), or Inf / NULL to test every cell's neighborhood")
  }
  if (profile == "normalized" && !(dist %in% c("cor", "euclidean"))) stop("profile = 'normalized' takes dist = 'cor' or 'euclidean'")
  if (profile == "raw" && dist == "euclidean") stop("dist = 'euclidean' needs profile = 'normalized'")
  noise.correction <- isTRUE(noise.correction)
  if (pooled.statistic == "raw") z.grid <- seq(0, 200, by = 0.02)   # histograms of the raw F
  if (post == "averaged" && !(adjust.method == "perm-fdr" || (adjust.method == "bh" && p.source == "pooled"))) {
    if (verbose) message("post = 'averaged' needs adjust.method = 'perm-fdr' or 'bh' with p.source = 'pooled': post = 'none'"); post <- "none" }
  if (post == "averaged" && !adjust) post <- "none"
  block.vars <- block.vars %||% design$block.vars   # the model's permutation strata unless overridden
  if (identical(design$contrast_spec$type, "term")) stop("cluster-free shifts need a contrast test (two groups or a numeric step)")
  samples <- intersect(rownames(design$F), levels(droplevels(factor(sample.per.cell))))
  if (length(samples) < 4) stop("fewer than 4 samples shared by the design and the cells")
  spc <- factor(as.character(sample.per.cell[colnames(cm)]), levels = samples)
  if (profile == "normalized") {
    if (is.null(lib)) { if (verbose) message("lib not given: library sizes over the genes of cm"); lib <- Matrix::colSums(cm) }
    lib <- if (!is.null(names(lib))) as.numeric(lib[colnames(cm)]) else as.numeric(lib)
    if (length(lib) != ncol(cm) || anyNA(lib)) stop("lib must give a library size for every cell of cm")
  }
  keep.cells <- !is.na(spc)
  nn.map <- NULL
  if (!all(keep.cells)) {                      # cells of samples outside the design are dropped: re-index the neighbourhoods
    cm <- cm[, keep.cells, drop = FALSE]; spc <- spc[keep.cells]; if (profile == "normalized") lib <- lib[keep.cells]
    nn.map <- cumsum(keep.cells) - 1L; nn.map[!keep.cells] <- NA_integer_
  }
  n <- length(samples)
  pairs <- t(utils::combn(n, 2)); storage.mode(pairs) <- "integer"
  m <- length(nns.per.cell)
  nn.list <- lapply(nns.per.cell, as.integer)
  if (!is.null(nn.map)) nn.list <- lapply(nn.list, function(v) { w <- nn.map[v + 1L]; w[!is.na(w)] })
  # global permutations over all samples; neighbourhood profiles, distances and the test per cell run in C++
  meta <- meta[samples, , drop = FALSE]
  des <- design; des$F <- design$F[samples, , drop = FALSE]
  gplan <- permutationPlan(des, meta, samples, scheme = permutation, block.vars = block.vars, n.permutations = n.permutations, max.enumerate = max.enumerate)
  P.global <- withSeed(seed, drawPermutations(gplan, n.permutations)); storage.mode(P.global) <- "integer"
  B <- ncol(P.global)
  cF <- design$contrast.F[colnames(des$F)]
  spec <- design$contrast_spec
  lev.var <- if (!is.null(spec) && spec$type %in% c("simple", "marginal") && !grepl(":", spec$term, fixed = TRUE) && spec$term %in% names(meta) && !is.numeric(meta[[spec$term]])) spec$term else NULL
  level.code <- rep(-1L, n)
  if (!is.null(lev.var)) { g <- as.character(meta[[lev.var]]); level.code <- ifelse(is.na(g), -1L, ifelse(g == spec$den, 1L, ifelse(g == spec$num, 2L, 0L))) }
  cell.names <- names(nns.per.cell) %||% colnames(cm)
  grp.code <- integer(0); grp.levels <- NULL
  if (!is.null(groups)) { gf <- as.factor(groups[cell.names]); grp.levels <- levels(gf); grp.code <- as.integer(gf) - 1L; grp.code[is.na(grp.code)] <- -1L }
  if (verbose) message(sprintf("Testing %d cells (%d samples) with %d permutations (%s)...", m, n, B, gplan$scheme))
  # the 1-hop graph neighbourhoods for the post-processing (0-based, the cell itself included), re-indexed like the units
  nn.graph <- if (is.null(nns.graph)) nn.list else { g <- lapply(nns.graph[cell.names], as.integer); if (!is.null(nn.map)) g <- lapply(g, function(v) { w <- nn.map[v + 1L]; w[!is.na(w)] }); g }
  prof.opts <- NULL
  if (profile == "normalized") {
    scale.value <- if (identical(profile.scale, "median")) stats::median(lib[lib > 0]) else as.numeric(profile.scale)
    if (length(scale.value) != 1 || !is.finite(scale.value) || scale.value <= 0) stop("profile.scale must be 'median' or a positive number")
    prof.opts <- list(scale = scale.value, dist = if (dist == "cor") 1L else 2L, noise = weighting == "noise", tau = weighting == "noise" && sample.floor,
                      correct = weighting == "noise" && isTRUE(noise.correction) && !gsd, unit_detection = neighbourhood.gene.frac)
    if (weighting == "noise") {
      if (is.null(noise.curves)) {
        if (verbose) message(sprintf("Fitting the noise curves on %d units...", min(n.fit.units, m)))
        noise.curves <- fitShiftNoiseCurves(cm, lib, spc, nn.list, scale.value, n.units = n.fit.units, min.n.obs.per.samp = min.n.obs.per.samp, seed = seed)
      } else if (!isTRUE(all.equal(noise.curves$scale, scale.value))) stop("noise.curves were fitted at a different scale")
      m1 <- shiftPoissonTable()
      prof.opts <- c(prof.opts, list(mu_grid = noise.curves$mu.grid, phi = noise.curves$phi, tau2 = noise.curves$tau2, lam_grid = m1$lam.grid, s_grid = m1$s.grid, V1 = m1$V1))
    }
  }
  # tested centres (gsd with post = "averaged" and a finite centres.per.cell): a greedy cover of the cells in random order (seeded)
  # makes a cell a centre while fewer than centres.per.cell centres lie in its graph neighbourhood; only the centres' neighbourhoods
  # are tested, and each one's statistic (observed and permuted alike) is scattered to every cell of the neighbourhood, so a cell's
  # statistic is the mean over the tested neighbourhoods that contain it. The per-cell test stays exact (the same function of the data under every relabeling).
  nn.test <- nn.list; avg.lists <- nn.graph; index.units <- NULL
  if (!is.null(centres)) {   # given centres: only their neighbourhoods are tested, each scattered to its cells
    if (!gsd || post != "averaged") stop("centres apply to statistic = 'gsd' with post = 'averaged'")
    centres <- as.integer(centres); if (!length(centres) || anyNA(centres) || any(centres < 1L | centres > m)) stop("centres must index nns.per.cell")
    is.idx <- logical(m); is.idx[centres] <- TRUE; index.units <- which(is.idx)
    avg.lists <- lapply(seq_len(m), function(i) if (is.idx[i]) sort(unique(c(i - 1L, nn.list[[i]]))) else integer(0)); nn.test[!is.idx] <- list(integer(0))
    if (verbose) message(sprintf("Testing the neighborhoods of %d given centers", length(index.units)))
  } else if (gsd && post == "averaged" && !is.null(centres.per.cell) && is.finite(centres.per.cell)) {
    cover <- integer(m); is.idx <- logical(m); ord <- withSeed(seed, sample.int(m))
    for (i in ord) if (cover[i] < centres.per.cell) { is.idx[i] <- TRUE; j <- nn.graph[[i]] + 1L; cover[j] <- cover[j] + 1L }
    index.units <- which(is.idx)
    avg.lists <- lapply(seq_len(m), function(i) if (is.idx[i]) sort(unique(c(i - 1L, nn.list[[i]]))) else integer(0))
    nn.test[!is.idx] <- list(integer(0))
    if (verbose) message(sprintf("Testing the neighborhoods of %d of %d cells (centres.per.cell = %d)", length(index.units), m, as.integer(centres.per.cell)))
  }
  gsd.opts <- if (gsd) list(mode = if (gsd.variance == "moderated") 2L else 1L, gamma = as.numeric(gsd.gamma), d0 = as.numeric(gsd.d0),
                            winsor_k = if (robust == "winsor") as.numeric(winsor.k) else 0, fl_norm = gsdResidualScale(gplan$crosses.levels), symmetric = as.integer(isTRUE(getOption("cacoa.gsd.symmetric", TRUE)))) else NULL
  kr <- cluster_free_shift_stream(cm, as.integer(spc), nn.test, FALSE, as.integer(min.n.obs.per.samp), if (profile == "raw") dist else "cor", log.vectors,
                                  pairs - 1L, n, des$F, cF, as.integer(level.code), as.integer(min.samp.per.level), as.integer(gplan$strata), as.integer(gplan$in.set),
                                  P.global, gplan$scheme == "freedman-lane", TRUE, as.integer(n.cores), if (gsd) 0L else c(none = 0L, huber = 1L, winsor = 2L)[[robust]], robust.k,
                                  as.numeric(z.grid), as.integer(grp.code), return.perm, post == "averaged", if (post == "averaged") unname(avg.lists) else NULL, vectorize, pooled.statistic == "standardized",
                                  prof.opts, if (profile == "normalized") lib else numeric(0), gsd.opts, !is.null(index.units))
  Fobs <- as.numeric(kr$stat); pval <- as.numeric(kr$p); shifts <- as.numeric(kr$shift); n.samp <- as.integer(kr$n); n.perm <- as.integer(kr$n_perm); z <- as.numeric(kr$z)
  mx <- as.numeric(kr$max); mn <- as.numeric(kr$min)
  valid <- which(is.finite(z))
  names(Fobs) <- names(pval) <- names(shifts) <- names(z) <- names(n.samp) <- names(n.perm) <- cell.names
  # enumerated relabelings include the identity: the kernel's (ge + 1) / (B + 1) becomes the exact ge / B
  if (gplan$exhaustive) pval <- ifelse(n.perm > 0, pmax(pval * (n.perm + 1) - 1, 1) / pmax(n.perm, 1), pval)
  z.adj <- NULL; p.adj <- NULL; index.cells <- NULL; fdr.curve <- NULL
  perm.hist <- list(grid = as.numeric(kr$grid), hist = kr$hist)
  q <- ncol(des$F)
  p.param <- stats::pf(Fobs, 1, pmax(n.samp - q, 1), lower.tail = FALSE); p.param[!is.finite(Fobs) | gsd] <- NA_real_   # (no parametric reference for gsd)
  names(p.param) <- cell.names
  s.base <- if (pooled.statistic == "standardized") z else Fobs   # the statistic the histograms were accumulated on
  p.pooled <- pooledPermutationP(s.base, perm.hist, alternative = "greater"); names(p.pooled) <- cell.names   # no per-cell floor
  # the reported z-score (z.pooled) is the tail z of the pooled p (one-sided: z = qnorm(1 - p)), calibrated and not floored;
  # the permutation-standardized statistic (observed minus permutation mean over permutation sd) is kept as z.std
  z.std <- z
  z <- stats::qnorm(pmin(pmax(p.pooled, 1e-300), 1 - 1e-16), lower.tail = FALSE); z[!is.finite(p.pooled)] <- NA_real_; names(z) <- cell.names
  stats.perm <- NULL
  if (return.perm && length(kr$perm)) { stats.perm <- kr$perm; colnames(stats.perm) <- cell.names }
  if (post == "averaged") {   # the kernel averaged the statistic over the cell and its graph neighbours under every permutation; the shift estimate follows
    tested <- is.finite(z)
    shifts <- if (is.null(index.units)) averageVectorOverLists(shifts, nn.graph) else {   # mean over the index units scattered into each cell
      ii <- rep(seq_along(avg.lists), lengths(avg.lists)); jj <- unlist(avg.lists) + 1L; ok <- is.finite(shifts[ii])
      M <- Matrix::sparseMatrix(i = jj[ok], j = ii[ok], x = 1, dims = c(m, m)); s0 <- shifts; s0[!is.finite(s0)] <- 0
      as.numeric(M %*% s0) / pmax(Matrix::rowSums(M), 1) }
    shifts[!tested] <- NA_real_; names(shifts) <- cell.names
    p.param[] <- NA_real_
  }
  if (adjust && length(valid)) {
    if (adjust.method == "max-t") {
      z0 <- z.std; z0[!is.finite(z0)] <- 0
      z.adj <- adjustedZScoresMaxStat(z_obs = z0, alt = 1L, max_vals_in = mx, min_vals_in = mn, wins = wins, smooth = smooth,
                                      nn_ids = nn.list, non_zero_ids = as.integer(valid - 1L))
      z.adj[!is.finite(z.std)] <- NA_real_
    } else if (adjust.method == "perm-fdr") {
      ad <- adjustCellScores(s.base, NULL, method = "perm-fdr", alternative = "greater", perm = perm.hist)
      z.adj <- ad$z.adj; p.adj <- ad$p.adj; names(p.adj) <- cell.names; fdr.curve <- ad$fdr.curve
    } else {
      p.in <- switch(p.source, parametric = p.param, permutation = pval, pooled = p.pooled,
                     z = { pz <- stats::pnorm(z.std, lower.tail = FALSE); pz[!is.finite(z.std)] <- NA_real_; names(pz) <- cell.names; pz })   # normal tail of the permutation-standardized z: no floor
      z.in <- stats::qnorm(pmin(pmax(p.in, 1e-300), 1 - 1e-16), lower.tail = FALSE); z.in[!is.finite(p.in)] <- NA_real_
      ad <- adjustCellScores(z.in, p.in, method = adjust.method, nn.list = nn.list, alternative = "greater", seed = seed)
      z.adj <- ad$z.adj; p.adj <- ad$p.adj; names(p.adj) <- cell.names
      if (!is.null(ad$index)) index.cells <- cell.names[ad$index]
    }
    names(z.adj) <- cell.names
    if (post == "inherit" && !is.null(p.adj)) {
      ih <- inheritNeighbourQList(p.adj, shifts, nn.graph, cell.names)
      p.adj.cell <- p.adj; p.adj <- ih$q
      z.adj <- pmax(stats::qnorm(pmin(pmax(p.adj, 1e-300), 1 - 1e-16), lower.tail = FALSE), 0); z.adj[!is.finite(p.adj)] <- NA_real_; names(z.adj) <- cell.names
    }
  }
  group.test <- if (length(grp.levels)) groupTestFromAccum(kr, grp.levels, alternative = "greater", min.cells = 1) else NULL
  shifts.smoothed <- NULL
  if (smooth) {
    s0 <- shifts; s0[!is.finite(s0)] <- 0
    shifts.smoothed <- applyMedianFilterES(s0, nn_ids = nn.list, non_zero_ids = as.integer(which(is.finite(shifts)) - 1L))
    shifts.smoothed[!is.finite(shifts)] <- NA_real_
    names(shifts.smoothed) <- cell.names
  }
  out.post <- list(post = post); if (post == "inherit" && exists("p.adj.cell", inherits = FALSE)) out.post$p.adj.cell <- p.adj.cell
  list(stat = Fobs, p.value = pval, p.param = p.param, p.pooled = p.pooled, z.pooled = z, z.std = z.std, z.adj = z.adj, p.adj = p.adj, index.cells = index.cells, pooled = z, adj = z.adj, effect = shifts,
       fdr.curve = fdr.curve, perm.hist = perm.hist, group.test = group.test, stats.perm = stats.perm, post = post, p.adj.cell = out.post$p.adj.cell,
       shifts = shifts, shifts.smoothed = shifts.smoothed, n.samples = n.samp, n.perm = n.perm,
       noise.curves = noise.curves,
       settings = list(statistic = statistic, exhaustive = gplan$exhaustive, n.distinct = gplan$n.distinct, neighbourhood.gene.frac = neighbourhood.gene.frac, centres.per.cell = if (!is.null(index.units)) centres.per.cell else NULL,
                       n.neighbourhoods.tested = if (!is.null(index.units)) length(index.units) else sum(lengths(nn.list) >= min.n.obs.per.samp),
                       gsd = if (gsd) list(variance = gsd.variance, gamma = gsd.gamma, d0 = gsd.d0, winsor.k = if (robust == "winsor") winsor.k else NULL) else NULL,
                       dist = dist, profile = profile, scale = prof.opts$scale, weighting = if (profile == "normalized") weighting else "none",
                       sample.floor = isTRUE(prof.opts$tau), noise.correction = isTRUE(prof.opts$correct), permutation = gplan$scheme, block.vars = block.vars, n.permutations = B, seed = seed, min.n.obs.per.samp = min.n.obs.per.samp, robust = robust,
                       adjust = adjust, adjust.method = adjust.method, p.source = p.source, smooth = smooth, n.cells.tested = length(valid)))
}


#' Scale scan of the expression-shift map
#'
#' The per-gene statistic of the shift map (`statistic = "gsd"`, [clusterFreeExpressionShifts()]) on diffusion
#' neighbourhoods of the cell graph at one or several scales ([diffusionNeighbourhoodLists()]: the scale is the share
#' of the cells a typical neighbourhood covers). Every scale is tested under the same relabelings (the seed); the
#' tested neighbourhoods' statistics are scattered to their cells (`post = "averaged"`), standardized per cell by their
#' permutation mean and sd (the kernel's z) and turned into a p-value against the scale's pooled permutation null
#' (every cell under every relabeling). A cell's scan statistic is its smallest p over the scales and its q-value the
#' pooled permutation false discovery rate of that statistic (against the same quantity under every relabeling; no
#' Bonferroni factor), as [scanDensityScales()]. Each scale's own permutation FDR is kept (`Q`): the single-scale map.
#' @param cm genes x cells count matrix (the columns in the order of the graph's vertices)
#' @param sample.per.cell named factor of samples per cell
#' @param graph igraph cell graph over the columns of `cm` (reordered to them)
#' @param design,meta the model and its sample metadata (see [clusterFreeExpressionShifts()])
#' @param scale scales: shares of the cells (one or several)
#' @param beta kernel strength per scale (default NULL: solved for the scale)
#' @param cover tested neighbourhoods per cell (default 2)
#' @param min.size smallest neighbourhood (cells; default NULL: the kNN unit size [defaultNeighbourhoodSize()] of the cells and
#'   samples, so that a neighbourhood centred on a weakly connected cell, whose kernel stays on itself, still holds enough cells
#'   per sample to be tested)
#' @param neighbourhoods the neighbourhoods per scale, built beforehand ([diffusionNeighbourhoodLists()] on this graph, in the
#'   order of `scale`; default NULL: built here); they are label-free, so a calibration over relabelings reuses them
#' @param n.permutations number of relabelings
#' @param seed seed of the relabelings and of the centres (the same for every scale)
#' @param groups optional cell groups (the group-level test per scale)
#' @param n.cores cores
#' @param verbose print the neighbourhoods per scale
#' @param ... further arguments of [clusterFreeExpressionShifts()] (`lib`, `profile`, `weighting`, `neighbourhood.gene.frac`, `gsd.d0`, `permutation`, `block.vars`, `robust`, ...)
#' @return list: `q` and `p` (the scan), `best` (index of the scale of the smallest p per cell), `P`, `E`, `Q` (cells x
#'   scales: pooled p, shift estimate, per-scale q), `Z` (the standardized statistic per scale), `scales` (data.frame: `scale`, `beta`, `size` (median effective
#'   size of the centres), `n.centres`), `neighbourhoods` (per scale), `group.test` (per scale), `settings` (of the first scale), `cells`, `B`, `time`
#' @export
scanShiftScales <- function(cm, sample.per.cell, graph, design, meta, scale, beta = NULL, cover = 2L, min.size = NULL, neighbourhoods = NULL, n.permutations = 499, seed = 1, groups = NULL, n.cores = 1, verbose = FALSE, ...) {
  cells <- colnames(cm)
  if (!identical(cells, igraph::V(graph)$name)) { if (!setequal(cells, igraph::V(graph)$name)) stop("the graph's vertices must be the columns of cm")
    graph <- igraph::permute(graph, match(igraph::V(graph)$name, cells)) }   # the graph in the order of the columns
  if (!is.numeric(scale) || !length(scale) || any(!is.finite(scale) | scale <= 0 | scale >= 1)) stop("scale must be shares of the cells in (0, 1)")
  N <- length(cells); K <- length(scale); if (is.null(beta)) beta <- rep(NA_real_, K) else if (length(beta) != K) stop("beta must have one strength per scale")
  if (is.null(min.size)) min.size <- defaultNeighbourhoodSize(N, nlevels(droplevels(as.factor(sample.per.cell))))   # the floor of a neighbourhood: the kNN unit size
  P <- matrix(NA_real_, N, K, dimnames = list(cells, NULL)); E <- P; Q <- P; Z <- P; m.perm <- NULL; B <- NULL; gt <- vector("list", K); settings <- NULL; n.all <- NULL
  info <- data.frame(scale = scale, beta = beta, size = NA_real_, n.centres = NA_integer_); nn1 <- graphNeighbourLists0(graph, cells); t0 <- Sys.time(); nbs <- vector("list", K)
  if (!is.null(neighbourhoods) && length(neighbourhoods) != K) stop("neighbourhoods must have one entry per scale")
  for (j in seq_len(K)) {
    nb <- if (is.null(neighbourhoods)) diffusionNeighbourhoodLists(graph, scale[j], cover = cover, seed = seed, beta = if (is.finite(beta[j])) beta[j] else NULL, min.size = min.size, n.cores = n.cores) else neighbourhoods[[j]]
    if (!identical(nb$cells, cells)) stop("the neighborhoods were built on another graph or cell order")
    info$beta[j] <- nb$beta; info$size[j] <- stats::median(nb$size); info$n.centres[j] <- length(nb$centres); nbs[[j]] <- nb
    if (verbose) message(sprintf("scale %g: heat kernel beta %.3g, %d centers, neighborhoods of %.0f cells (median, %d smallest)", scale[j], nb$beta, length(nb$centres), stats::median(lengths(nb$units)), min(lengths(nb$units))))
    nn <- rep(list(integer(0)), N); names(nn) <- cells; nn[nb$centres] <- unname(nb$units)
    s <- clusterFreeExpressionShifts(cm, sample.per.cell, nn, design, meta, statistic = "gsd", centres = nb$centres, post = "averaged", nns.graph = nn1, n.permutations = n.permutations, seed = seed,
                                     return.perm = TRUE, adjust = TRUE, smooth = FALSE, adjust.method = "perm-fdr", pooled.statistic = "standardized", groups = groups, n.cores = n.cores, verbose = FALSE, ...)
    sp <- s$stats.perm[, cells, drop = FALSE]; B <- nrow(sp); st <- as.numeric(s$stat[cells])
    mu <- colMeans(sp); sdv <- sqrt(pmax(colMeans(sp^2) - mu^2, 0) * B / (B - 1)); sdv[!is.finite(sdv) | sdv <= 0] <- NA_real_   # the kernel's per-cell standardization
    z.obs <- (st - mu) / sdv; sp <- sweep(sweep(sp, 2, mu), 2, sdv, "/")
    # calls of the observed and of every relabeled data set, all standardized alike, for the test of the chosen scale
    if (is.null(n.all)) n.all <- matrix(0L, B + 1L, K); n.all[, j] <- scanCallCounts(z.obs, sp, standardize = TRUE, n.cores = n.cores)
    ps <- pooledScanStep(z.obs, sp, identity, n.cores = n.cores); rm(sp)
    P[, j] <- ps$p; Q[, j] <- ps$q; E[, j] <- as.numeric(s$shifts[cells]); Z[, j] <- z.obs; m.perm <- if (is.null(m.perm)) ps$pp else pmin(m.perm, ps$pp); rm(ps)
    gt[j] <- list(s$group.test); if (j == 1L) settings <- s$settings; rm(s); invisible(gc())
  }
  m.obs <- apply(P, 1, min); best <- max.col(-P, ties.method = "first"); sc <- pooledMinPFDR(m.obs, m.perm, B, n.cores = n.cores); rm(m.perm)
  auto <- autoScaleTest(n.all[1, ], n.all[-1, , drop = FALSE], sum(rowSums(is.finite(Q)) > 0))   # the test of the automatically chosen scale
  list(auto = auto, q = stats::setNames(sc$q, cells), p = stats::setNames(sc$p, cells), best = stats::setNames(best, cells), P = P, E = E, Q = Q, Z = Z, scales = info, neighbourhoods = nbs, group.test = gt, settings = settings,
       cells = cells, B = B, time = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
