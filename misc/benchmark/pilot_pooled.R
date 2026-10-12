## Pilot for the pooled-p plan: one SCC null relabeling (k = 1 of null3) and three paired semi-synthetic sets with the
## extended battery (pooled p into BH, spatial FDR, neighbourhood BH, Cauchy). Usage: Rscript pilot_pooled.R <B> <n.cores>
args <- commandArgs(trailingOnly = TRUE); B <- as.integer(args[1]); n.cores <- as.integer(args[2])
source(file.path(Sys.getenv("CACOA_BENCH", "/home/pkharchenko/cacoa-dev/cacoa/misc/benchmark"), "common.R"))
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
out <- list()
## ---- null relabeling 1 of SCC ----
inp <- tm("load inputs", loadInputs("scc")); cs <- condSpec("scc"); meta <- inp$sample.meta
cond0 <- relevel(factor(as.character(meta[[cs[1]]])), cs[2]); spc <- inp$sample.per.cell; cells <- V(inp$graph)$name
dm <- densityMats(inp$graph, spc, n.cores = n.cores); nn <- neighbourLists(inp$graph)
cm.t <- Matrix::t(inp$counts[cells, topGenes(inp$counts, 1000)])
k <- 1; set.seed(1000 + k)
lab <- rep(NA_character_, nrow(meta))
for (lv in levels(cond0)) { i <- which(cond0 == lv); i <- sample(i); lab[i] <- rep_len(c("A", "B"), length(i))[sample(length(i))] }
meta2 <- data.frame(row.names = rownames(meta), cond = factor(lab, levels = c("A", "B")), stratum = cond0)
model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "stratum")
groups <- inp$cell.groups
heat <- tm("heat", runDensity(dm$heat, model, nn, B = B, n.cores = n.cores, seed = k, label = "heat", groups = groups))
nb <- tm("nb", runDensity(dm$nb, model, nn, B = B, n.cores = n.cores, seed = k, label = "nb.log", transform = "log", sizes = dm$sizes, groups = groups))
sh <- tm("shifts", runShifts(cm.t, spc, nn, model, B = B, n.cores = n.cores, seed = k, groups = groups))
calls <- rbind(cbind(test = "heat", callFractions(heat$q)), cbind(test = "nb.log", callFractions(nb$q)), cbind(test = "shifts", callFractions(sh$q)))
cat("\n== NULL (SCC relabeling 1, B =", B, "): cells called per method\n"); print(calls[calls$alpha == 0.1, ], row.names = FALSE)
gsum <- function(g, nm) data.frame(test = nm, method = c("perm", "cauchy.perm", "cauchy.param", "cauchy.pooled"),
  n.sig = c(sum(p.adjust(g$p.perm, "BH") < 0.05, na.rm = TRUE), sum(p.adjust(g$p.cauchy, "BH") < 0.05, na.rm = TRUE), sum(p.adjust(g$p.cauchy.param, "BH") < 0.05, na.rm = TRUE), sum(p.adjust(g$p.cauchy.pooled, "BH") < 0.05, na.rm = TRUE)),
  min.p = c(min(g$p.perm, na.rm = TRUE), min(g$p.cauchy, na.rm = TRUE), min(g$p.cauchy.param, na.rm = TRUE), min(g$p.cauchy.pooled, na.rm = TRUE)))
cat("\n== NULL: groups significant (BH 5%) and smallest group p per summary\n")
print(rbind(gsum(heat$groups, "heat"), gsum(nb$groups, "nb.log"), gsum(sh$groups, "shifts")), row.names = FALSE, digits = 3)
cat("\n== NULL: shifts per-group p by summary\n"); print(sh$groups[, c("group", "n.tested", "p.perm", "p.cauchy", "p.cauchy.pooled")], digits = 2, row.names = FALSE)
out$null <- list(calls = calls, groups = list(heat = heat$groups, nb = nb$groups, shifts = sh$groups))
rm(inp, dm, cm.t); gc()
## ---- paired semi-synthetic sets ----
for (id in c("large_f100_r1", "medium_f80_r1", "small_f100_r1", "null_r1")) {
  d <- file.path(S, "bench/semi/scc_paired", id)
  tr <- readRDS(file.path(d, "truth.rds")); sub <- readRDS(file.path(d, "inputs.rds"))
  cells <- V(sub$graph)$name; spc <- sub$sample.per.cell
  dm <- densityMats(sub$graph, spc, n.cores = n.cores); nn <- neighbourLists(sub$graph)
  meta2 <- data.frame(row.names = rownames(sub$sample.meta), cond = sub$sample.meta$condition, patient = sub$sample.meta$patient)
  model <- buildCacoaModel(meta2, formula = ~ cond, test = "cond", block.vars = "patient")
  heat <- tm(paste(id, "heat"), runDensity(dm$heat, model, nn, B = B, n.cores = n.cores, seed = 1, label = "heat", groups = sub$cell.groups))
  nb <- tm(paste(id, "nb"), runDensity(dm$nb, model, nn, B = B, n.cores = n.cores, seed = 1, label = "nb.log", transform = "log", sizes = dm$sizes, groups = sub$cell.groups))
  true.cells <- tr$true.cells
  adjacent <- setdiff(unique(unlist(lapply(which(cells %in% true.cells), function(i) cells[nn[[i]] + 1L]))), true.cells)
  score <- function(q, effect, method, cells, a = 0.1) {
    ok <- is.finite(q); cells <- cells[ok]; q <- q[ok]; effect <- effect[ok]
    tp.set <- cells %in% true.cells; near <- cells %in% adjacent; called <- q < a
    data.frame(method = method, n.called = sum(called), recall = if (any(tp.set)) mean(called[tp.set] & effect[tp.set] < 0) else NA,
               fdp.lenient = if (any(called)) mean(!tp.set[called] & !near[called]) else NA_real_, frac.other = mean(called[!tp.set & !near]))
  }
  res <- do.call(rbind, c(lapply(names(heat$q), function(m) score(heat$q[[m]], heat$effect[cells], paste("heat", m), cells)),
                          lapply(names(nb$q), function(m) score(nb$q[[m]], nb$effect[cells], paste("1-hop", m), cells))))
  if (file.exists(file.path(d, "milo.rds"))) { mi <- readRDS(file.path(d, "milo.rds"))$cells; res <- rbind(res, score(mi$milo.fdr, mi$milo.logFC, "Milo", mi$cell)) }
  cat(sprintf("\n== SEMI %s (type %s, removed %.0f%%, B = %d): FDR 10%%\n", id, tr$type, 100 * tr$frac, B)); print(res, row.names = FALSE, digits = 2)
  gg <- function(g, nm) { i <- which(g$group == tr$type); data.frame(test = nm, p.perm = g$p.perm[i], p.cauchy.perm = g$p.cauchy[i], p.cauchy.param = g$p.cauchy.param[i], p.cauchy.pooled = g$p.cauchy.pooled[i],
    n.sig.perm = sum(p.adjust(g$p.perm, "BH") < 0.05, na.rm = TRUE), n.sig.cauchy.pooled = sum(p.adjust(g$p.cauchy.pooled, "BH") < 0.05, na.rm = TRUE)) }
  if (!is.na(tr$type)) { cat("   planted type's group p by summary, and the number of significant groups (BH 5%):\n"); print(rbind(gg(heat$groups, "heat"), gg(nb$groups, "nb.log")), row.names = FALSE, digits = 2) }
  else { cat("   control-only split: groups significant (BH 5%): perm", sum(p.adjust(heat$groups$p.perm, "BH") < 0.05, na.rm = TRUE), "cauchy.pooled", sum(p.adjust(heat$groups$p.cauchy.pooled, "BH") < 0.05, na.rm = TRUE), "\n") }
  out[[id]] <- list(res = res, groups = list(heat = heat$groups, nb = nb$groups))
}
saveRDS(out, file.path(S, "bench", sprintf("pilot_pooled_B%d.rds", B)))
message("done")
