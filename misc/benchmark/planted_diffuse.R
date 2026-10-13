## Diffuse and gradient-like planted changes for the density map (complements planted_relabel.R, whose truth removes
## whole annotated types). Same balanced relabelings (seed 1000 + rep); arm-B cells are removed with a per-cell
## probability p_i that is (a) "gradient": a linear ramp in rank of the second eigenvector of the normalized adjacency
## within one large lineage (no type boundary coincides with the change), or (b) "bump": a heat-kernel bump around a
## random seed cell, wide enough that the cells above 40 % removal span several annotated types. Peak removal 80 %.
## Truth: p.remove over all cells; targets = remaining cells with p >= 0.4; clean negatives = p < 0.05 beyond two hops of
## any target. Usage: Rscript planted_diffuse.R <scc|pf|ms> generate [n.cores]   -> $S/bench/diffuse/<ds>/<field>_r<rep>/
args <- commandArgs(trailingOnly = TRUE); ds <- args[1]; mode <- args[2]; n.cores <- if (length(args) > 2) as.integer(args[3]) else 8
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
root <- file.path(BENCH, "diffuse", ds); dir.create(root, showWarnings = FALSE, recursive = TRUE)
n.rep <- 3; peak <- 0.8; target.frac <- 0.06   # the bump is widened until ~6 % of the cells sit above 40 % removal
lineages <- list(scc = c("KC_Diff", "KC_Basal", "KC_Cyc", "TSK"),
                 ms = c("EN-L2-3", "EN-L4", "EN-L5-6", "EN-MIX", "EN-PYR"),
                 pf = c("Ciliated", "Differentiating Ciliated", "Basal", "SCGB3A2+", "SCGB3A2+ SCGB1A1+", "MUC5B+", "KRT5-/KRT17+", "AT1", "AT2", "Transitional AT2", "Proliferating Epithelial Cells", "MUC5AC+ High"))

# second eigenvector of D^-1/2 A D^-1/2 on the lineage subgraph (largest component), as a rank in [0, 1]; oriented so
# that the first listed type of the lineage sits at the low end (removal rises away from it)
gradientField <- function(graph, cells.in, cg, first.type) {
  g <- induced_subgraph(graph, cells.in); comp <- components(g); g <- induced_subgraph(g, V(g)$name[comp$membership == which.max(comp$csize)])
  em <- embed_laplacian_matrix(g, no = 3, which = "la", type = "DAD", scaled = FALSE)$X
  v <- em[, 2]; if (mean(v[cg[V(g)$name] == first.type]) > mean(v)) v <- -v
  r <- (rank(v) - 1) / (length(v) - 1)
  p <- setNames(rep(0, vcount(graph)), V(graph)$name); p[V(g)$name] <- r * peak
  list(p = p, n.lineage = vcount(g))
}
# bump around a seed cell: cells ordered by heat-kernel proximity to the seed (beta 30); removal falls linearly from
# the peak at the seed to 0.4 at rank n.target = target.frac x cells and to 0 at rank 2 n.target
bumpField <- function(graph, seed.cell, n.cores) {
  cells <- V(graph)$name; ind <- matrix(as.numeric(cells == seed.cell), ncol = 1, dimnames = list(cells, NULL))
  s <- as.numeric(sccore::smoothSignalOnGraph(ind, filter = function(...) sccore::heatFilter(..., beta = 30), graph = graph, m = 50, n.cores = n.cores, progress = FALSE))
  r <- rank(-s, ties.method = "first") - 1; n.t <- round(target.frac * length(cells))
  p <- peak * pmax(0, 1 - r / (2 * n.t))
  list(p = setNames(p, cells), beta = 30, frac = mean(p >= 0.4))
}

if (mode == "generate") {
  inp <- tm("load inputs", loadInputs(ds)); cs <- condSpec(ds); meta <- inp$sample.meta
  cond0 <- relevel(factor(as.character(meta[[cs[1]]])), cs[2])
  cg <- inp$cell.groups; spc <- inp$sample.per.cell; cells <- V(inp$graph)$name
  stub <- inp$counts[, seq_len(min(20, ncol(inp$counts))), drop = FALSE]
  A <- as_adjacency_matrix(inp$graph, sparse = TRUE)[cells, cells]; A@x[] <- 1
  big <- names(which(table(cg) >= 1000))
  for (rep in seq_len(n.rep)) {
    set.seed(1000 + rep)                                   # the same relabeling as the planted sets
    lab <- rep(NA_character_, nrow(meta))
    for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
    meta2 <- data.frame(row.names = rownames(meta), condition = factor(lab, levels = c("A", "B")), stratum = cond0)
    arm.b <- rownames(meta2)[meta2$condition == "B"]
    for (field in c("gradient", "bump")) {
      id <- sprintf("%s_r%d", field, rep); d <- file.path(root, id); if (file.exists(file.path(d, "truth.rds"))) next
      dir.create(d, showWarnings = FALSE); set.seed(5e4 + 100 * rep + match(field, c("gradient", "bump")))
      deg <- degree(inp$graph); seeds <- names(cg)[cg %in% big & deg[names(cg)] >= median(deg)]
      fld <- if (field == "gradient") gradientField(inp$graph, names(cg)[cg %in% lineages[[ds]]], cg, lineages[[ds]][1]) else bumpField(inp$graph, sample(seeds, 1), n.cores)
      p <- fld$p[cells]
      in.b <- cells[spc[cells] %in% arm.b]; drop <- in.b[runif(length(in.b)) < p[in.b]]
      keep <- setdiff(cells, drop); sub <- subsetInputs(inp, cells = keep)
      sub$sample.meta <- meta2[levels(sub$sample.per.cell), , drop = FALSE]; sub$counts <- stub[V(sub$graph)$name, , drop = FALSE]
      saveRDS(sub, file.path(d, "inputs.rds"))
      tgt <- cells[p >= 0.4]; v <- as.numeric(cells %in% tgt); h2 <- as.numeric(A %*% (as.numeric(A %*% v) > 0)) > 0
      neg <- cells[p < 0.05 & !h2]
      types.hit <- names(which(tapply(p >= 0.4, cg[cells], mean) >= 0.2))
      n.b <- length(in.b)
      saveRDS(list(id = id, dataset = ds, regime = "diffuse", field = field, rep = rep, arm.b = arm.b, dropped = drop, p.remove = p, peak = peak, beta = fld$beta, n.lineage = fld$n.lineage,
                   closure = n.b / (n.b - length(drop)), true.cells = setdiff(tgt, drop), negatives = setdiff(neg, drop), types.hit = types.hit), file.path(d, "truth.rds"))
      message(sprintf("%s: %d cells, dropped %d of %d arm-B cells; %d targets (p >= 0.4) in %d types [%s]; %d clean negatives; closure %.3f",
                      id, vcount(sub$graph), length(drop), n.b, length(tgt), length(types.hit), paste(types.hit, collapse = ", "), length(neg), n.b / (n.b - length(drop))))
    }
  }
  quit(status = 0)
}
