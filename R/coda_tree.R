## Contrast tree and reference block of the compositional model (lmCoda): both are built from the contrast's cell-type
## loadings l = Psi beta (clr scale, zero-sum, reference-free) and their permutation null, so that nothing is refitted.

#' Unit contrast vector of a balance
#'
#' For the partition of some cell types into `plus` and `minus` sides, the orthonormal balance vector
#' `sqrt(n1 n2 / (n1 + n2)) (1/n1 on plus, -1/n2 on minus)` over `cell.types` (zero elsewhere). Its inner product with the
#' clr loadings is the contrast's effect on the balance; with a sample's clr log-frequencies, the sample's balance coordinate.
#' @param plus,minus cell types on the two sides
#' @param cell.types all cell types (names of the returned vector)
#' @return named numeric vector
#' @keywords internal
balanceVector <- function(plus, minus, cell.types) {
  n1 <- length(plus); n2 <- length(minus); v <- stats::setNames(numeric(length(cell.types)), cell.types)
  v[plus] <- 1 / n1; v[minus] <- -1 / n2; v * sqrt(n1 * n2 / (n1 + n2))
}

#' Contrast tree of a compositional result
#'
#' A sequential binary partition of the cell types that explains the contrast's compositional shift with as few balances
#' as possible. The contrast's effect on the balance of a partition (plus, minus) is `beta_b = psi_b' l` with `l` the clr
#' loadings of [lmCoda()] and `psi_b` the balance vector ([balanceVector()]), i.e. `sqrt(n1 n2 / (n1 + n2))` times the
#' difference of the mean loadings of the two sides; over the `D - 1` balances of any such tree the squared effects sum to
#' the squared length of the loading vector (the contrast's variance explained in ILR space), so the tree concentrating it
#' into the fewest balances is the simplest explanation of the shift. It is built top-down: within a set of cell types the
#' partition with the largest `beta_b^2` is a contiguous split of the loadings in sorted order (a two-group split of
#' one-dimensional values), so the `n - 1` contiguous splits are evaluated and the best one taken (`split = "effect"`, the
#' default; `"evidence"` takes the split with the largest effect over its permutation standard deviation), recursively.
#' Each node carries the effect in the balance coordinate and as the log-ratio of the geometric-mean proportions of its two
#' sides (`log2.ratio`, the "plus" side being the one with the larger loadings, i.e. higher at the contrast's target), an
#' exact permutation p-value from the same relabelings as the per-cell-type test (the balance of the permuted loadings),
#' Benjamini-Hochberg adjustment over the nodes, and the share of the shift it explains (`share`; `cum.share` accumulates
#' the shares in decreasing order). The leaves of the tree are the cell types in the order of their loadings.
#' @param res result of [lmCoda()] (as stored by `cao$estimateCellLoadings()`), with the permutation null of the loadings
#' @param split `"effect"` (default) or `"evidence"`, see above
#' @param alpha significance level of the adjusted node p-values (`significant`)
#' @return list: `nodes` (data.frame, one row per balance in the order of construction (root first): `node`, `parent`,
#'   `depth`, `n.plus`, `n.minus`, `plus`, `minus` (the sides, comma-separated), `effect` (balance coordinate),
#'   `log.ratio` and `log2.ratio` (difference of the mean clr loadings of the sides), `sd.null`, `p`, `padj`,
#'   `significant`, `share`, `cum.share`), `sides` (list per node with `plus` and `minus` as character vectors), `sbp`
#'   (cell types x balances, +1 / -1 / 0), `tree` (`phylo`, from [sbp2tree()]), `order` (cell types by loading),
#'   `total` (squared length of the loadings), `n.perm`, `split`, `alpha`
#' @export
contrastTree <- function(res, split = c("effect", "evidence"), alpha = 0.05) {
  split <- match.arg(split)
  if (is.null(res$contrast)) stop("not a contrast result of lmCoda()")
  l <- res$contrast$loadings$obs; P <- res$contrast$loadings$perm; cts <- names(l); D <- length(cts)
  if (D < 2) stop("the contrast tree needs at least two cell types")
  if (is.null(P)) stop("the permutation null of the loadings (contrast$loadings$perm) is missing")
  P <- as.matrix(P)[cts, , drop = FALSE]; B <- ncol(P); total <- sum(l^2)
  sbp <- matrix(0, D, D - 1, dimnames = list(cts, paste0("node", seq_len(D - 1)))); nodes <- vector("list", D - 1); sides <- vector("list", D - 1); k <- 0L
  queue <- list(list(members = cts[order(l)], parent = 0L, depth = 0L))   # members in increasing order of loading
  while (length(queue)) {
    cur <- queue[[1]]; queue <- queue[-1]; m <- cur$members; n <- length(m); lm <- l[m]; Pm <- P[m, , drop = FALSE]
    cand <- lapply(seq_len(n - 1), function(j) { minus <- m[seq_len(j)]; plus <- m[-seq_len(j)]; v <- balanceVector(plus, minus, m)
      list(plus = plus, minus = minus, effect = sum(v * lm), null = as.numeric(crossprod(v, Pm))) })
    crit <- vapply(cand, function(z) if (split == "effect") z$effect^2 else (z$effect / max(stats::sd(z$null), 1e-12))^2, numeric(1))
    best <- cand[[which.max(crit)]]; k <- k + 1L
    sbp[best$plus, k] <- 1; sbp[best$minus, k] <- -1; sides[[k]] <- list(plus = best$plus, minus = best$minus)
    p <- (sum(abs(best$null) >= abs(best$effect)) + 1) / (B + 1); lr <- mean(l[best$plus]) - mean(l[best$minus])
    nodes[[k]] <- data.frame(node = k, parent = cur$parent, depth = cur$depth, n.plus = length(best$plus), n.minus = length(best$minus),
                             plus = paste(best$plus, collapse = ","), minus = paste(best$minus, collapse = ","), effect = best$effect,
                             log.ratio = lr, log2.ratio = lr / log(2), sd.null = stats::sd(best$null), p = p, share = if (total > 0) best$effect^2 / total else NA_real_,
                             stringsAsFactors = FALSE)
    if (length(best$minus) > 1) queue[[length(queue) + 1]] <- list(members = best$minus, parent = k, depth = cur$depth + 1L)
    if (length(best$plus) > 1) queue[[length(queue) + 1]] <- list(members = best$plus, parent = k, depth = cur$depth + 1L)
  }
  nodes <- do.call(rbind, nodes); rownames(nodes) <- NULL
  nodes$padj <- stats::p.adjust(nodes$p, method = "BH"); nodes$significant <- is.finite(nodes$padj) & nodes$padj <= alpha
  o <- order(-nodes$share); nodes$cum.share <- NA_real_; nodes$cum.share[o] <- cumsum(nodes$share[o])
  tree <- tryCatch(sbp2tree(sbp), error = function(e) NULL)
  list(nodes = nodes, sides = sides, sbp = sbp, tree = tree, order = cts[order(l)], total = total, n.perm = B, split = split, alpha = alpha)
}

#' Reference block of cell types
#'
#' The largest set of cell types whose proportions do not change relative to each other under the contrast: the pairwise
#' differences of the clr loadings, `l_i - l_j`, are tested against their permutation null (the same relabelings as the
#' per-cell-type test; for groups of cell types the balance between the groups), and the cell types are merged
#' agglomeratively, the pair of groups with the largest p-value first, while that p-value is at least `p.thresh`. The
#' reference is the largest resulting group; groups within a factor `size.frac` of the largest cannot be told apart from it by
#' the composition alone (two groups moving apart are one shift, seen from either side), so among those the group closest to
#' the center of the composition is taken (mean clr loading nearest zero). Centering the loadings on this block shows every
#' cell type's change relative to a set of cell types that do not change among themselves, which the clr loadings (relative
#' to the geometric mean of all cell types) do not: when one large cell type expands, every other cell type's clr loading is
#' negative.
#' @param l named clr loadings of the contrast
#' @param P their permutation null (cell types x relabelings)
#' @param p.thresh p-value below which two groups are kept apart (default 0.05)
#' @param size.frac groups at least this fraction of the largest group's size compete for the reference (default 0.5)
#' @return list: `idx`, `celltypes`, `mean` (mean loading of the block), `clusters` (all groups after the merging), `sizes`,
#'   `p.pairs` (cell types x cell types, the pairwise permutation p-values), `p.thresh`, `size.frac`, `method = "block"`
#' @export
referenceBlock <- function(l, P, p.thresh = 0.05, size.frac = 0.5) {
  cts <- names(l); P <- as.matrix(P)[cts, , drop = FALSE]; B <- ncol(P)
  pairP <- function(a, b) { v <- balanceVector(a, b, cts); d <- sum(v * l); dp <- as.numeric(crossprod(v, P)); (sum(abs(dp) >= abs(d)) + 1) / (B + 1) }
  p.pairs <- matrix(1, length(cts), length(cts), dimnames = list(cts, cts))
  if (length(cts) > 1) for (i in seq_len(length(cts) - 1)) for (j in (i + 1):length(cts)) p.pairs[i, j] <- p.pairs[j, i] <- pairP(cts[i], cts[j])
  clusters <- as.list(cts)
  while (length(clusters) > 1) {
    n <- length(clusters); pm <- matrix(NA_real_, n, n)
    for (i in seq_len(n - 1)) for (j in (i + 1):n) pm[i, j] <- pairP(clusters[[i]], clusters[[j]])
    best <- which(pm == max(pm, na.rm = TRUE), arr.ind = TRUE)[1, ]
    if (pm[best[1], best[2]] < p.thresh) break
    clusters[[best[1]]] <- c(clusters[[best[1]]], clusters[[best[2]]]); clusters <- clusters[-best[2]]
  }
  # the largest block; blocks within a factor `size.frac` of it cannot be told apart from it by the composition alone (two blocks
  # moving apart are one shift), so among them the one closest to the center of the composition (mean clr loading nearest zero)
  sizes <- lengths(clusters); cand <- which(sizes >= size.frac * max(sizes))
  cand <- cand[order(vapply(clusters[cand], function(m) abs(mean(l[m])), numeric(1)), -sizes[cand])][1]
  ref <- cts[cts %in% clusters[[cand]]]
  list(idx = match(ref, cts), celltypes = ref, mean = mean(l[ref]), clusters = clusters, sizes = sizes, p.pairs = p.pairs, p.thresh = p.thresh, size.frac = size.frac, method = "block")
}

#' Per-cell-type p-values relative to the reference block
#'
#' The change of each cell type relative to the reference, `l_c - mean(l over the reference)`, against its permutation null
#' (the same statistic on the permuted loadings), two-sided, with Benjamini-Hochberg adjustment over the cell types.
#' @param l named clr loadings
#' @param P their permutation null (cell types x relabelings)
#' @param ref cell types of the reference
#' @return list with `stat` (the centered loadings), `pval`, `padj`
#' @keywords internal
referenceRelativeP <- function(l, P, ref) {
  cts <- names(l); P <- as.matrix(P)[cts, , drop = FALSE]; B <- ncol(P)
  stat <- l - mean(l[ref]); null <- sweep(P, 2, colMeans(P[ref, , drop = FALSE]))
  pval <- stats::setNames((rowSums(abs(null) >= abs(stat)) + 1) / (B + 1), cts)
  list(stat = stat, pval = pval, padj = stats::p.adjust(pval, method = "BH"))
}

#' Per-sample clr log-frequencies of a compositional result, adjusted for the model's covariates
#'
#' From the ILR coordinates and basis of [lmCoda()]: `log f = ilr Psi'` per sample (clr scale). When the model has
#' covariates besides the tested variable (design columns at which the two contrast endpoints agree), each sample's value
#' is moved to the covariate values of the endpoints (`log f - (F - F*) coef Psi'`, `F*` the design row with the covariate
#' columns set to the endpoints' values), so that the samples are shown as the contrast compares them. The samples' side
#' of the contrast is taken from the model's metadata: for a contrast of two levels of a factor, `"baseline"`, `"target"`
#' or `"other"`; for a numeric variable its value.
#' @param res result of [lmCoda()]
#' @return list: `logf` (samples x cell types), `group` (per sample: factor baseline / target / other, or numeric),
#'   `labels` (of baseline and target), `numeric` (logical)
#' @keywords internal
codaSampleLogFreq <- function(res) {
  logf <- as.matrix(res$ilr) %*% t(res$psi); colnames(logf) <- rownames(res$psi)
  model <- res$model; F <- model$F; end <- model$contrast_endpoints_F
  if (!is.null(F) && !is.null(end) && !is.null(res$coef) && nrow(F) == nrow(logf)) {
    cov.cols <- which(end$num == end$den)
    if (length(cov.cols)) { Fs <- F; Fs[, cov.cols] <- matrix(end$den[cov.cols], nrow(F), length(cov.cols), byrow = TRUE)
      logf <- logf - ((F - Fs) %*% res$coef) %*% t(res$psi) }
  }
  spec <- model$contrast_spec; labels <- model$contrast_endpoint_labels %||% list(baseline = "baseline", target = "target")
  group <- NULL; numeric <- FALSE
  if (!is.null(spec$term) && !is.null(model$meta) && spec$term %in% names(model$meta)) {
    v <- model$meta[[spec$term]]; if (!is.null(rownames(model$meta)) && all(rownames(logf) %in% rownames(model$meta))) v <- v[match(rownames(logf), rownames(model$meta))]
    if (is.numeric(v)) { group <- v; numeric <- TRUE }
    else group <- factor(ifelse(as.character(v) == as.character(spec$num), "target", ifelse(as.character(v) == as.character(spec$den), "baseline", "other")), levels = c("baseline", "target", "other"))
  }
  list(logf = logf, group = group, labels = labels, numeric = numeric)
}
