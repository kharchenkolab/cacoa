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
res <- list(con = con, sample.meta = meta, cell.groups = droplevels(cg[keep]), sample.per.cell = droplevels(spc[keep]), truth = truth, dropped = drop)
saveRDS(res, out)
cat("cells kept", length(keep), "of", length(cg), "; dropped", length(drop), "\n"); print(table(res$cell.groups, meta[as.character(res$sample.per.cell), "treatment"]))
