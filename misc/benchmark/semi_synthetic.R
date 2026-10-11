## Semi-synthetic power study: the control samples of a dataset are split into two pseudo-arms and a fraction of one
## cell type's cells is removed from the "treated" arm; every method is scored on the remaining cells of that type.
## Usage: Rscript semi_synthetic.R <scc|pf|ms> generate|run <B> <n.cores>
## Datasets: $S/bench/semi/<ds>/<id>/inputs.rds (Milo format; counts are a stub, Milo's graph route ignores expression)
## and truth.rds; cacoa results in cacoa.rds.
## Regimes: "unpaired" (default): the control samples are split into two arms (4 + 5 or 5 + 5); "paired": every control
## sample is split at random into two halves, one per arm, so that the arms are blocked by donor (block.vars = patient).
## Usage: Rscript semi_synthetic.R <scc|pf> generate|run [B] [n.cores] [unpaired|paired]
args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]; mode <- args[2]; B <- if (length(args) > 2) as.integer(args[3]) else 1999; n.cores <- if (length(args) > 3) as.integer(args[4]) else 16
regime <- if (length(args) > 4) args[5] else "unpaired"
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
root <- file.path(BENCH, "semi", if (regime == "paired") paste0(ds, "_paired") else ds); dir.create(root, showWarnings = FALSE, recursive = TRUE)
fracs <- c(0.5, 0.8, 1.0); n.rep <- 3; target.frac <- c(small = 0.02, medium = 0.07, large = 0.2)

if (mode == "generate") {
  inp <- tm("load inputs", loadInputs(ds))
  cs <- condSpec(ds)
  ctrl <- rownames(inp$sample.meta)[as.character(inp$sample.meta[[cs[1]]]) == cs[2]]
  base <- tm("control subset", subsetInputs(inp, samples = ctrl)); rm(inp); gc()
  cg <- base$cell.groups; spc <- base$sample.per.cell
  # candidate types: present (>= 5 cells) in at least 80 % of the control samples
  pres <- table(cg, spc) >= 5; ok <- names(which(rowMeans(pres) >= 0.75))
  sz <- table(cg)[ok] / length(cg)
  planted <- c(); for (f in target.frac) { cand <- sz[setdiff(names(sz), planted)]; planted <- c(planted, names(cand)[which.min(abs(cand - f))]) }
  names(planted) <- names(target.frac)
  message("control samples: ", paste(ctrl, collapse = ", "), "\nplanted types: ", paste(sprintf("%s = %s (%.1f%%)", names(planted), planted, 100 * sz[planted]), collapse = "; "))
  stub <- base$counts[, seq_len(min(20, ncol(base$counts))), drop = FALSE]
  for (rep in seq_len(n.rep)) {
    set.seed(100 * rep)
    if (regime == "paired") {
      # every control sample is split in two halves: <sample>_A (reference arm) and <sample>_B (treated arm)
      half <- stats::setNames(character(length(spc)), names(spc))
      for (cs in split(names(spc), spc)) half[cs] <- sample(rep_len(c("A", "B"), length(cs)))
      spc2 <- factor(paste0(as.character(spc), "_", half)); names(spc2) <- names(spc)
      samples2 <- levels(spc2)
      meta <- data.frame(row.names = samples2, condition = factor(sub(".*_", "", samples2), levels = c("A", "B")), patient = factor(sub("_[AB]$", "", samples2)))
      arm.b <- samples2[meta$condition == "B"]
    } else {
      spc2 <- spc
      arm.b <- sample(ctrl, floor(length(ctrl) / 2))
      meta <- data.frame(row.names = ctrl, condition = factor(ifelse(ctrl %in% arm.b, "B", "A"), levels = c("A", "B")))
    }
    for (ty in c("null", names(planted))) for (f in (if (ty == "null") 0 else fracs)) {
      id <- if (ty == "null") sprintf("null_r%d", rep) else sprintf("%s_f%02d_r%d", ty, round(100 * f), rep)
      d <- file.path(root, id); if (file.exists(file.path(d, "truth.rds"))) next
      dir.create(d, showWarnings = FALSE)
      set.seed(1e4 * rep + 100 * f + match(ty, c("null", names(planted))))
      victims <- if (ty == "null") character(0) else names(cg)[cg == planted[[ty]] & spc2 %in% arm.b]
      drop <- if (length(victims)) sample(victims, round(f * length(victims))) else character(0)
      keep <- setdiff(names(cg), drop)
      sub <- subsetInputs(base, cells = keep)
      sub$sample.per.cell <- droplevels(spc2[V(sub$graph)$name])
      sub$sample.meta <- meta[levels(sub$sample.per.cell), , drop = FALSE]; sub$counts <- stub[V(sub$graph)$name, , drop = FALSE]
      saveRDS(sub, file.path(d, "inputs.rds"))
      saveRDS(list(id = id, dataset = ds, regime = regime, type = if (ty == "null") NA_character_ else planted[[ty]], size.class = ty, frac = f, rep = rep, arm.b = arm.b, dropped = drop,
                   true.cells = if (ty == "null") character(0) else names(cg)[cg == planted[[ty]] & !(names(cg) %in% drop)]), file.path(d, "truth.rds"))
      message(id, ": ", vcount(sub$graph), " cells, dropped ", length(drop), " of ", length(victims))
    }
  }
  quit(status = 0)
}

## ---- run cacoa on every dataset ----
for (d in list.dirs(root, recursive = FALSE)) {
  out <- file.path(d, "cacoa.rds"); if (file.exists(out)) next
  tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds"))
  cells <- V(sub$graph)$name; spc <- sub$sample.per.cell
  t0 <- Sys.time()
  dm <- densityMats(sub$graph, spc, n.cores = n.cores)
  nn <- neighbourLists(sub$graph)
  meta2 <- data.frame(row.names = rownames(sub$sample.meta), cond = sub$sample.meta$condition)
  if (!is.null(sub$sample.meta$patient)) { meta2$patient <- sub$sample.meta$patient; model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "patient") } else model <- condModel(meta2)
  heat <- runDensity(dm$heat, model, nn, B = B, n.cores = n.cores, seed = 1, label = "heat", max.t = TRUE, groups = sub$cell.groups)
  nb <- runDensity(dm$nb, model, nn, B = B, n.cores = n.cores, seed = 1, label = "nb.log", transform = "log", sizes = dm$sizes, max.t = TRUE, groups = sub$cell.groups)
  gs <- list(heat = heat$groups, nb = nb$groups)
  # neighbours of the planted type (for the lenient false-positive accounting)
  adj1 <- unique(unlist(lapply(which(cells %in% tr$true.cells), function(i) cells[nn[[i]] + 1L])))
  saveRDS(list(id = tr$id, heat = heat[c("p.param", "effect", "q")], nb = nb[c("p.param", "effect", "q")], groups = gs, cell.groups = sub$cell.groups,
               adjacent = adj1, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))), out)
  message(tr$id, " done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
}
message("done")
