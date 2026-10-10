## Build the compact simulated dataset used by vignettes/walkthrough_short.Rmd from test/shifts_sim_objects.rds:
## the 20 samples of Group1 (control) and Group2 (treated) of the no-batch simulation, as a Conos object with the
## induced joint graph and embedding, plus a sample metadata table with `treatment` and a synthetic `age`.
## usage: Rscript prepare_sim_dataset.R ../../../test/shifts_sim_objects.rds ../../../test/sim_treatment.rds
suppressPackageStartupMessages({library(conos); library(Matrix)})
args <- commandArgs(TRUE); in.file <- args[1]; out.file <- args[2]
sim <- readRDS(in.file)$no_batch
con <- sim$data.object
meta <- con$misc$sample_meta
keep <- rownames(meta)[meta$Group %in% c("Group1", "Group2")]
cells <- names(sim$sample.per.cell)[sim$sample.per.cell %in% keep]
con2 <- Conos$new(con$samples[keep], n.cores = 1)
con2$graph <- igraph::induced_subgraph(con$graph, intersect(igraph::V(con$graph)$name, cells))
con2$embedding <- con$embedding[cells, ]
set.seed(7)
sample.meta <- data.frame(row.names = keep,
                          treatment = factor(ifelse(meta[keep, "Group"] == "Group1", "control", "treated"), levels = c("control", "treated")))
noise <- rnorm(nrow(sample.meta), 0, 7); noise <- noise - ave(noise, sample.meta$treatment)       # centred within arm
sample.meta$age <- round(ifelse(sample.meta$treatment == "control", 52, 58) + noise)               # treated arm 6 years older on average
sample.meta <- sample.meta[sample(nrow(sample.meta)), ]              # shuffle the row order
new.names <- setNames(sprintf("S%02d", seq_along(keep)), rownames(sample.meta))   # neutral sample names
rownames(sample.meta) <- new.names[rownames(sample.meta)]
names(con2$samples) <- new.names[names(con2$samples)]
sample.per.cell <- setNames(factor(new.names[as.character(sim$sample.per.cell[cells])]), cells)
con2$misc <- list(sample_meta = sample.meta)
out <- list(con = con2, sample.meta = sample.meta, cell.groups = droplevels(sim$cell.groups[cells]), sample.per.cell = sample.per.cell)
saveRDS(out, out.file)
cat("samples", length(con2$samples), "cells", length(cells), "edges", igraph::ecount(con2$graph), "file", file.size(out.file) / 1e6, "MB\n")
print(table(sample.meta$treatment)); print(tapply(sample.meta$age, sample.meta$treatment, summary))
