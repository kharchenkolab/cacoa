## Planted cell-density simulation from test/sim_treatment.rds: in the treated arm 60 % of IN-SST cells are removed, in the
## control arm 60 % of OPC cells; L2/3 keeps its planted expression shift with unchanged density; everything else is null.
## usage: Rscript prepare_density_sim.R ../../../test/sim_treatment.rds ../../../test/sim_density.rds
suppressPackageStartupMessages({library(conos); library(Matrix)})
args <- commandArgs(TRUE); d <- readRDS(args[1]); out <- args[2]
set.seed(11)
meta <- d$sample.meta; cg <- d$cell.groups; spc <- d$sample.per.cell
drop <- c(); for (ct in c("IN-SST", "OPC")) {
  arm <- if (ct == "IN-SST") "treated" else "control"
  cells <- names(cg)[cg == ct & meta[as.character(spc[names(cg)]), "treatment"] == arm]
  drop <- c(drop, sample(cells, round(0.6 * length(cells))))
}
keep <- setdiff(names(cg), drop)
con <- d$con
for (s in names(con$samples)) { p <- con$samples[[s]]; k <- intersect(colnames(p), keep); con$samples[[s]] <- subset(p, cells = k) }
con$graph <- igraph::induced_subgraph(con$graph, intersect(igraph::V(con$graph)$name, keep))
con$embedding <- con$embedding[keep, ]
truth <- data.frame(celltype = levels(cg), density = ifelse(levels(cg) == "IN-SST", "down in treated", ifelse(levels(cg) == "OPC", "up in treated", "none")),
                    expression = c(`L2/3` = "strong", `IN-PV` = "strong", `IN-SST` = "strong", OPC = "moderate", `L5/6-CC` = "moderate", `AST-PP` = "moderate",
                                   `AST-FB` = "weak", Microglia = "none", `Neu-NRGN` = "none")[levels(cg)])
# the base simulation already differs in composition between the arms, so the effective truth for any abundance test is
# the per-sample proportion difference after the removal, not just the planted changes: record it alongside
tab <- table(droplevels(spc[keep]), droplevels(cg[keep])); fr <- tab / rowSums(tab); arm <- meta$treatment[match(rownames(tab), rownames(meta))]
emp <- t(sapply(colnames(fr), function(k) { a <- fr[arm == "control", k]; b <- fr[arm == "treated", k]
  c(log2fc = log2(mean(b) / mean(a)), p = stats::wilcox.test(a, b, exact = FALSE)$p.value) }))
truth$composition.log2fc <- emp[truth$celltype, "log2fc"]; truth$composition.p <- emp[truth$celltype, "p"]
truth$composition <- ifelse(truth$composition.p < 0.01, ifelse(truth$composition.log2fc > 0, "up in treated", "down in treated"), "none")
res <- list(con = con, sample.meta = meta, cell.groups = droplevels(cg[keep]), sample.per.cell = droplevels(spc[keep]), truth = truth, dropped = drop)
saveRDS(res, out)
cat("cells kept", length(keep), "of", length(cg), "; dropped", length(drop), "\n"); print(table(res$cell.groups, meta[as.character(res$sample.per.cell), "treatment"]))
