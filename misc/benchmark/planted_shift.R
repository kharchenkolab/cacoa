## Planted expression shifts: unit x post-processing grid for the cluster-free shift map.
## Truth: in arm B of a balanced relabeling (seed 1000 + rep, as the planted sets), the cells of one annotated type
## (small / medium / large, the types of planted_relabel.R) or the cells along the gradient field of planted_diffuse.R
## have a fraction f of their counts (over the tested genes) replaced by draws from the mean profile of the most
## distant annotated type: a state change toward another identity, strength f = 0.1 or 0.25 (graded by p / peak on the
## gradient). The graph is left as built (the shift is tested where the cells sit). True region = all cells of the
## type (both arms; the arm-A cells are the comparison); clean negatives = beyond two hops of it.
## Units: the 1-hop neighbourhood (cacoa's current shift unit) and kNN neighbourhoods of 100 / 200 / 400 cells;
## post-processing: none, neighbour-averaged statistic (recomputed under every permutation), inheritance.
## Usage: Rscript planted_shift.R <ds> generate                       -> $S/bench/shift/<ds>/<id>/{inputs,truth}.rds
##        Rscript planted_shift.R <ds> run <B> <n.cores>              -> <set>/shift.rds
##        Rscript planted_shift.R <ds> null <B> <n.cores> [n.relabel] -> $S/bench/null_shift/<ds>/r<k>.rds
args <- commandArgs(trailingOnly = TRUE); ds <- args[1]; mode <- args[2]; B <- if (length(args) > 2) as.integer(args[3]) else 1999
n.cores <- if (length(args) > 3) as.integer(args[4]) else 8; n.relabel <- if (length(args) > 4) as.integer(args[5]) else 15
only <- if (length(args) > 5) args[6] else NULL   # run mode: a single set id (pilot)
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
root <- file.path(BENCH, "shift", ds); n.genes <- 1000; n.rep <- 3; ks <- c(100, 200, 400)
strengths <- as.numeric(strsplit(Sys.getenv("CACOA_SHIFT_STRENGTHS", "0.1,0.25"), ",")[[1]])   # generate mode: mixing strengths (sets that exist are skipped)
geneSet <- function(counts) { cn <- setdiff(colnames(counts), c("Patient", "Tissue: 0=Normal, 1=Tumor")); topGenes(counts[, cn, drop = FALSE], n.genes) }
# mean log-normalized profile per type over the gene set, and the most distant type (correlation) from a reference profile
typeProfiles <- function(cm, cg) { cpm <- Matrix::Diagonal(x = 1 / pmax(Matrix::rowSums(cm), 1)) %*% cm; sapply(split(seq_len(nrow(cm)), cg), function(i) log1p(1e4 * Matrix::colMeans(cpm[i, , drop = FALSE]))) }
# replace a fraction f_c of each target cell's counts by multinomial draws from `prob` (thinning + redraw; depth kept in expectation)
plantMix <- function(cm, targets, f, prob) {
  cm <- methods::as(cm, "TsparseMatrix"); idx <- match(targets, rownames(cm)); fc <- setNames(f, targets)
  tot <- Matrix::rowSums(cm)[idx]
  fcell <- fc[rownames(cm)[cm@i + 1L]]; keep <- which(!is.na(fcell))
  cm@x[keep] <- rbinom(length(keep), as.integer(cm@x[keep]), 1 - fcell[keep])
  add <- sapply(seq_along(idx), function(j) rmultinom(1, round(fc[j] * tot[j]), prob))
  M <- Matrix::sparseMatrix(i = rep(idx, each = nrow(add)), j = rep(seq_len(ncol(cm)), length(idx)), x = as.numeric(add), dims = dim(cm), dimnames = dimnames(cm))
  Matrix::drop0(methods::as(cm, "CsparseMatrix") + M)
}
if (mode == "generate") {
  dir.create(root, showWarnings = FALSE, recursive = TRUE)
  inp <- tm("load inputs", loadInputs(ds)); cs <- condSpec(ds); meta <- inp$sample.meta
  cond0 <- relevel(factor(as.character(meta[[cs[1]]])), cs[2]); cg <- inp$cell.groups; spc <- inp$sample.per.cell; cells <- V(inp$graph)$name
  genes <- geneSet(inp$counts); cm <- inp$counts[cells, genes]; prof <- typeProfiles(cm, cg[cells]); big <- names(which(table(cg) >= 500))
  A <- as_adjacency_matrix(inp$graph, sparse = TRUE)[cells, cells]; A@x[] <- 1
  types <- sapply(c("small", "medium", "large"), function(ty) readRDS(file.path(BENCH, "planted", ds, sprintf("%s_f100_r1", ty), "truth.rds"))$type)
  for (rep in seq_len(n.rep)) {
    set.seed(1000 + rep); lab <- rep(NA_character_, nrow(meta))
    for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
    meta2 <- data.frame(row.names = rownames(meta), condition = factor(lab, levels = c("A", "B")), stratum = cond0); arm.b <- rownames(meta2)[meta2$condition == "B"]
    grad <- readRDS(file.path(BENCH, "diffuse", ds, sprintf("gradient_r%d", rep), "truth.rds"))
    for (what in c(names(types), "gradient")) for (f in strengths) {
      id <- sprintf("%s_f%02d_r%d", what, round(100 * f), rep); d <- file.path(root, id); if (file.exists(file.path(d, "truth.rds"))) next
      dir.create(d, showWarnings = FALSE); set.seed(7e4 + 100 * rep + round(100 * f) + match(what, c(names(types), "gradient")))
      if (what == "gradient") { p <- grad$p.remove[cells]; region <- cells[p >= 0.4]; w <- p / grad$peak; ref <- log1p(1e4 * Matrix::colMeans(Matrix::Diagonal(x = 1 / pmax(Matrix::rowSums(cm[region, ]), 1)) %*% cm[region, ])) }
      else { ty <- types[[what]]; region <- cells[cg[cells] == ty]; w <- setNames(as.numeric(cells %in% region), cells); ref <- prof[, ty] }
      cand <- setdiff(big, if (what == "gradient") grad$types.hit else types[[what]]); partner <- cand[which.min(cor(ref, prof[, cand]))]
      prob <- expm1(prof[, partner]); prob <- prob / sum(prob)
      targets <- cells[spc[cells] %in% arm.b & w[cells] > 0]
      cm2 <- plantMix(cm, targets, f * w[targets], prob)
      v <- as.numeric(cells %in% region); h2 <- as.numeric(A %*% (as.numeric(A %*% v) > 0)) > 0
      neg <- cells[!h2 & (if (what == "gradient") p < 0.05 else TRUE)]
      saveRDS(list(graph = inp$graph, counts = cm2, cell.groups = cg[cells], sample.per.cell = spc[cells], sample.meta = meta2), file.path(d, "inputs.rds"))
      saveRDS(list(id = id, dataset = ds, regime = "shift", what = what, type = if (what == "gradient") NA_character_ else types[[what]], strength = f, rep = rep, arm.b = arm.b, partner = partner,
                   targets = targets, true.cells = region, negatives = neg, genes = genes), file.path(d, "truth.rds"))
      message(sprintf("%s: %d target cells of %d in the region mixed toward %s at f = %.2f; %d clean negatives", id, length(targets), length(region), partner, f, length(neg)))
    }
  }
  quit(status = 0)
}
# ---- the grid on one input: units x post-processing ----
runShiftGrid <- function(graph, cm, spc, model, B, seed, n.cores) {
  cells <- V(graph)$name; cm.t <- Matrix::t(cm[cells, ])
  A <- as_adjacency_matrix(graph, sparse = TRUE)[cells, cells]; A@x[] <- 1; An <- A; Matrix::diag(An) <- 1; An <- Matrix::Diagonal(x = 1 / Matrix::rowSums(An)) %*% An
  units <- c(list("1-hop" = neighbourLists(graph)), setNames(lapply(ks, function(k) { N <- methods::as(knnNeighbourhoods(graph, k), "TsparseMatrix"); nn <- split(N@j, factor(N@i + 1L, levels = seq_along(cells))); nn <- lapply(nn, as.integer); names(nn) <- cells; nn }), paste0("knn k=", ks)))
  res <- list(); times <- c()
  for (u in names(units)) {
    t0 <- Sys.time()
    r <- clusterFreeExpressionShifts(cm.t, spc, units[[u]], model, model$meta, dist = "cor", n.permutations = B, seed = seed, adjust = TRUE, smooth = FALSE, n.cores = n.cores, verbose = FALSE, adjust.method = "perm-fdr", return.perm = TRUE)
    st <- r$stat[cells]; st[!is.finite(st)] <- 0; Sp <- r$stats.perm[, cells, drop = FALSE]; Sp[!is.finite(Sp)] <- 0
    res[[paste0(u, " | none")]] <- list(q = r$p.adj[cells], effect = r$shifts[cells])
    tav <- setNames(as.numeric(An %*% st), cells); Tav <- as.matrix(Matrix::t(An %*% Matrix::t(Sp)))
    res[[paste0(u, " | averaged")]] <- list(q = permutationFDR(tav, Tav, alternative = "greater")$q, effect = setNames(as.numeric(An %*% r$shifts[cells]), cells))
    res[[paste0(u, " | inherit")]] <- inheritQ(r$p.adj[cells], r$shifts[cells], A, cells)
    times[u] <- as.numeric(difftime(Sys.time(), t0, units = "secs")); rm(r, Sp, Tav); gc()
  }
  list(res = res, times = times)
}
if (mode == "run") {
  for (d in list.dirs(root, recursive = FALSE)) {
    if (!is.null(only) && basename(d) != only) next
    out <- file.path(d, if (is.null(only)) "shift.rds" else sprintf("shift_pilot_B%d.rds", B)); lock <- file.path(d, "shift.lock"); if (file.exists(out) || file.exists(lock)) next
    file.create(lock); t0 <- Sys.time(); tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds"))
    meta2 <- data.frame(row.names = rownames(sub$sample.meta), cond = sub$sample.meta$condition, stratum = sub$sample.meta$stratum)
    model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
    g <- runShiftGrid(sub$graph, sub$counts, sub$sample.per.cell, model, B, tr$rep, n.cores)
    saveRDS(list(id = tr$id, res = g$res, times = g$times, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), out); unlink(lock)
    message(tr$id, " done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s (", paste(names(g$times), round(g$times), collapse = ", "), ")")
  }
} else if (mode == "null") {
  out.dir <- file.path(BENCH, "null_shift", ds); dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  inp <- tm("load inputs", loadInputs(ds)); cs <- condSpec(ds); meta <- inp$sample.meta; cells <- V(inp$graph)$name
  cond0 <- relevel(factor(as.character(meta[[cs[1]]])), cs[2]); genes <- geneSet(inp$counts); cm <- inp$counts[cells, genes]
  for (k in seq_len(n.relabel)) {
    f <- file.path(out.dir, sprintf("r%03d.rds", k)); lock <- paste0(f, ".lock"); if (file.exists(f) || file.exists(lock)) next
    file.create(lock); t0 <- Sys.time(); set.seed(1000 + k); lab <- rep(NA_character_, nrow(meta))
    for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
    meta2 <- data.frame(row.names = rownames(meta), cond = factor(lab, levels = c("A", "B")), stratum = cond0)
    model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
    g <- runShiftGrid(inp$graph, cm, inp$sample.per.cell, model, B, k, n.cores)
    calls <- do.call(rbind, lapply(names(g$res), function(m) do.call(rbind, lapply(c(0.05, 0.1), function(a) data.frame(method = m, alpha = a, n.tested = sum(is.finite(g$res[[m]]$q)), n.called = sum(g$res[[m]]$q < a, na.rm = TRUE))))))
    saveRDS(list(k = k, calls = calls, times = g$times), f); unlink(lock)
    message(sprintf("[%s] relabeling %d done in %.0f s; max calls at 5%%: %d", ds, k, as.numeric(difftime(Sys.time(), t0, units = "secs")), max(calls$n.called[calls$alpha == 0.05])))
  }
}
message("done")
