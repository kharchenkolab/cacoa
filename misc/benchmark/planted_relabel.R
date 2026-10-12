## Unpaired power study with realistic variability: the balanced relabelings of the null study (arms split within the
## true condition, 9-11 whole samples per arm) with a planted depletion of one cell type in arm B.
## Usage: Rscript planted_relabel.R <scc|pf|ms> generate|run [B] [n.cores]
## Sets: $S/bench/planted/<ds>/<id>/inputs.rds (Milo format) + truth.rds; cacoa results in cacoa.rds.
## Truth: the remaining cells of the planted type (both arms) are "down"; every other cell is expected UP by the
## closure factor 1 / (1 - removed fraction of arm B's cells); the graph neighbours of planted cells are ambiguous.
args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]; mode <- args[2]; B <- if (length(args) > 2) as.integer(args[3]) else 1999; n.cores <- if (length(args) > 3) as.integer(args[4]) else 16
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
root <- file.path(BENCH, "planted", ds); dir.create(root, showWarnings = FALSE, recursive = TRUE)
fracs <- c(0.5, 0.8, 1.0); n.rep <- 3; target.frac <- c(small = 0.02, medium = 0.07, large = 0.2)

if (mode == "generate") {
  inp <- tm("load inputs", loadInputs(ds)); cs <- condSpec(ds); meta <- inp$sample.meta
  cond0 <- relevel(factor(as.character(meta[[cs[1]]])), cs[2])
  cg <- inp$cell.groups; spc <- inp$sample.per.cell
  pres <- table(cg, spc) >= 5; ok <- names(which(rowMeans(pres) >= 0.75))
  sz <- table(cg)[ok] / length(cg)
  planted <- c(); for (f in target.frac) { cand <- sz[setdiff(names(sz), planted)]; planted <- c(planted, names(cand)[which.min(abs(cand - f))]) }
  names(planted) <- names(target.frac)
  message("planted types: ", paste(sprintf("%s = %s (%.1f%%)", names(planted), planted, 100 * sz[planted]), collapse = "; "))
  stub <- inp$counts[, seq_len(min(20, ncol(inp$counts))), drop = FALSE]
  for (rep in seq_len(n.rep)) {
    set.seed(1000 + rep)                                   # the same relabeling as null3's relabeling `rep`
    lab <- rep(NA_character_, nrow(meta))
    for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
    meta2 <- data.frame(row.names = rownames(meta), condition = factor(lab, levels = c("A", "B")), stratum = cond0)
    arm.b <- rownames(meta2)[meta2$condition == "B"]
    for (ty in c("null", names(planted))) for (f in (if (ty == "null") 0 else fracs)) {
      id <- if (ty == "null") sprintf("null_r%d", rep) else sprintf("%s_f%02d_r%d", ty, round(100 * f), rep)
      d <- file.path(root, id); if (file.exists(file.path(d, "truth.rds"))) next
      dir.create(d, showWarnings = FALSE)
      set.seed(1e4 * rep + 100 * f + match(ty, c("null", names(planted))))
      victims <- if (ty == "null") character(0) else names(cg)[cg == planted[[ty]] & spc %in% arm.b]
      drop <- if (length(victims)) sample(victims, round(f * length(victims))) else character(0)
      keep <- setdiff(names(cg), drop)
      sub <- subsetInputs(inp, cells = keep)
      sub$sample.meta <- meta2[levels(sub$sample.per.cell), , drop = FALSE]; sub$counts <- stub[V(sub$graph)$name, , drop = FALSE]
      saveRDS(sub, file.path(d, "inputs.rds"))
      n.b <- sum(spc %in% arm.b)
      saveRDS(list(id = id, dataset = ds, regime = "relabel", type = if (ty == "null") NA_character_ else planted[[ty]], size.class = ty, frac = f, rep = rep, arm.b = arm.b, dropped = drop,
                   closure = n.b / (n.b - length(drop)), true.cells = if (ty == "null") character(0) else names(cg)[cg == planted[[ty]] & !(names(cg) %in% drop)]), file.path(d, "truth.rds"))
      message(id, ": ", vcount(sub$graph), " cells, dropped ", length(drop), " of ", length(victims), " (closure ", round(n.b / (n.b - length(drop)), 3), ")")
    }
  }
  quit(status = 0)
}

## ---- run cacoa on every set ----
for (d in list.dirs(root, recursive = FALSE)) {
  out <- file.path(d, "cacoa.rds"); lock <- file.path(d, "cacoa.lock")
  if (file.exists(out) || file.exists(lock)) next                        # several processes may share the root
  file.create(lock); on.exit(unlink(lock), add = TRUE)
  tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds"))
  cells <- V(sub$graph)$name; spc <- sub$sample.per.cell
  t0 <- Sys.time()
  dm <- densityMats(sub$graph, spc, n.cores = n.cores); nn <- neighbourLists(sub$graph)
  meta2 <- data.frame(row.names = rownames(sub$sample.meta), cond = sub$sample.meta$condition, stratum = sub$sample.meta$stratum)
  model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")           # the realistic unpaired analysis
  heat <- runDensity(dm$heat, model, nn, B = B, n.cores = n.cores, seed = tr$rep, label = "heat", groups = sub$cell.groups, variants = TRUE)
  nb <- runDensity(dm$nb, model, nn, B = B, n.cores = n.cores, seed = tr$rep, label = "nb.log", transform = "log", sizes = dm$sizes, groups = sub$cell.groups, variants = TRUE)
  adj1 <- unique(unlist(lapply(which(cells %in% tr$true.cells), function(i) cells[nn[[i]] + 1L])))
  keep <- c("p.param", "p.perm", "p.pooled", "effect", "q", "scale")
  saveRDS(list(id = tr$id, heat = heat[keep], nb = nb[keep], groups = list(heat = heat$groups, nb = nb$groups), cell.groups = sub$cell.groups,
               adjacent = adj1, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), out)
  unlink(lock)
  message(tr$id, " done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
}
message("done")
