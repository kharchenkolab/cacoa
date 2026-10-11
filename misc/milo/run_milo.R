## Milo (miloR 1.6.0) differential abundance on cacoa benchmark datasets.
## Usage: Rscript run_milo.R <sim|ms|scc> [n.cores]
## Routes: "conos" = Conos joint graph used as Milo kNN graph (primary);
##         "pca"   = Seurat PCA (30 PCs) on joint counts + Milo's own kNN graph (k = 30).
args <- commandArgs(trailingOnly=TRUE); ds <- args[1]
n.cores <- if (length(args) > 1) as.integer(args[2]) else 16
suppressMessages({library(Matrix); library(igraph); library(miloR); library(SingleCellExperiment); library(ggplot2)})
wd <- "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad/milo"
setwd(wd)
timing <- list(); tic <- function() Sys.time(); toc <- function(nm, t0) { timing[[nm]] <<- as.numeric(difftime(Sys.time(), t0, units="secs")); message(sprintf("[%s] %s: %.1f s", ds, nm, timing[[nm]])) }

## ---- parameters per dataset ----
par <- list(
  sim=list(cond.col="treatment", ref="control", prop=0.2, designs=list(cond=~condition)),
  ms =list(cond.col="diagnosis", ref="Control", prop=0.1, designs=list(cond=~condition)),
  scc=list(cond.col="condition", ref="Normal",  prop=0.1, designs=list(cond=~condition, patient_cond=~patient+condition))
)[[ds]]
k.conos <- 30   # value stored in milo@.k for the conos route; used only by the k-distance FDR weighting (k-th NN distance in the 2D embedding)
k.pca <- 30; d.pca <- 30
weightings <- c("graph-overlap", "k-distance")

## ---- inputs ----
t0 <- tic()
inp <- readRDS(paste0("inputs_", ds, ".rds"))
g <- inp$graph; cells <- V(g)$name
sm <- inp$sample.meta
sm$condition <- relevel(factor(as.character(sm[[par$cond.col]])), ref=par$ref)
if (ds == "scc") sm$patient <- factor(sm$patient)
cell.groups <- droplevels(factor(inp$cell.groups[cells]))
sample.per.cell <- factor(as.character(inp$sample.per.cell[cells]), levels=rownames(sm))
stopifnot(!any(is.na(sample.per.cell)), !any(is.na(cell.groups)))
counts <- t(inp$counts[cells, , drop=FALSE])   # genes x cells
emb <- inp$embedding[cells, ]; colnames(emb) <- c("UMAP1", "UMAP2")
sce <- SingleCellExperiment(assays=list(counts=counts),
                            colData=DataFrame(celltype=cell.groups, sample=sample.per.cell, row.names=cells))
reducedDim(sce, "UMAP") <- emb
rm(inp); invisible(gc())
toc("load_inputs", t0)
message(sprintf("[%s] %d cells, %d genes, %d samples, %d cell types; design table: %s", ds, ncol(sce), nrow(sce), nrow(sm), nlevels(cell.groups), paste(capture.output(print(table(sm$condition))), collapse=" ")))

## ---- helpers ----
map_cells <- function(nh, res) {
  ## per cell: min SpatialFDR over the neighbourhoods containing it; logFC of the most significant one (min SpatialFDR, ties -> min PValue)
  tr <- summary(nh)   # i = cell, j = nhood
  tr <- tr[tr$x > 0, ]
  r <- res[match(tr$j, res$Nhood), ]
  o <- order(tr$i, r$SpatialFDR, r$PValue)
  tr <- tr[o, ]; r <- r[o, ]
  first <- !duplicated(tr$i)
  out <- data.frame(cell=rownames(nh), celltype=as.character(cell.groups), sample=as.character(sample.per.cell),
                    milo.logFC=NA_real_, milo.fdr=NA_real_, milo.pval=NA_real_, n.nhoods=0L, stringsAsFactors=FALSE)
  out$milo.logFC[tr$i[first]] <- r$logFC[first]
  out$milo.fdr[tr$i[first]] <- r$SpatialFDR[first]
  out$milo.pval[tr$i[first]] <- r$PValue[first]
  nn <- table(tr$i); out$n.nhoods[as.integer(names(nn))] <- as.integer(nn)
  out
}
summarise_ct <- function(res, cdf) {
  cts <- levels(cell.groups)
  f <- function(ct) {
    r <- res[res$celltype == ct, ]; cc <- cdf[cdf$celltype == ct, ]
    inh <- !is.na(cc$milo.fdr)
    data.frame(celltype=ct, n.nhoods=nrow(r),
               nh.fdr10.up=sum(r$SpatialFDR < 0.1 & r$logFC > 0), nh.fdr10.down=sum(r$SpatialFDR < 0.1 & r$logFC < 0),
               nh.fdr05.up=sum(r$SpatialFDR < 0.05 & r$logFC > 0), nh.fdr05.down=sum(r$SpatialFDR < 0.05 & r$logFC < 0),
               frac.nh.fdr10=if (nrow(r)) mean(r$SpatialFDR < 0.1) else NA, frac.nh.fdr05=if (nrow(r)) mean(r$SpatialFDR < 0.05) else NA,
               median.logFC=if (nrow(r)) median(r$logFC) else NA,
               n.cells=nrow(cc), n.cells.no.nhood=sum(!inh),
               cells.fdr10.up=sum(inh & cc$milo.fdr < 0.1 & cc$milo.logFC > 0), cells.fdr10.down=sum(inh & cc$milo.fdr < 0.1 & cc$milo.logFC < 0),
               cells.fdr05.up=sum(inh & cc$milo.fdr < 0.05 & cc$milo.logFC > 0), cells.fdr05.down=sum(inh & cc$milo.fdr < 0.05 & cc$milo.logFC < 0),
               frac.cells.fdr10=mean(inh & cc$milo.fdr < 0.1), frac.cells.fdr05=mean(inh & cc$milo.fdr < 0.05),
               stringsAsFactors=FALSE)
  }
  s <- do.call(rbind, lapply(cts, f)); rownames(s) <- NULL; s
}
beeswarm_plot <- function(res, title, file) {
  res$sig <- factor(ifelse(res$SpatialFDR < 0.1, "SpatialFDR < 0.1", "n.s."), levels=c("n.s.", "SpatialFDR < 0.1"))
  res$celltype <- factor(res$celltype, levels=rev(levels(cell.groups)))
  lab <- aggregate(cbind(n=rep(1, nrow(res)), sig=res$SpatialFDR < 0.1), by=list(celltype=res$celltype), FUN=sum)
  lab$txt <- sprintf("%d / %d", lab$sig, lab$n)
  xr <- range(res$logFC, na.rm=TRUE)
  p <- ggplot(res, aes(x=celltype, y=logFC, colour=sig)) +
    geom_hline(yintercept=0, colour="grey70", linewidth=0.4) +
    ggbeeswarm::geom_quasirandom(size=1.1, alpha=0.8, width=0.35, groupOnX=TRUE) +
    geom_text(data=lab, aes(x=celltype, y=xr[2] + 0.12 * diff(xr), label=txt), inherit.aes=FALSE, size=3, hjust=0, colour="grey30") +
    scale_colour_manual(values=c("n.s."="#9a9a95", "SpatialFDR < 0.1"="#2a78d6"), drop=FALSE, name=NULL) +
    scale_y_continuous(expand=expansion(mult=c(0.05, 0.25))) +
    coord_flip() + theme_bw(base_size=11) +
    theme(legend.position="top", panel.grid.minor=element_blank(), panel.grid.major.y=element_blank()) +
    labs(x=NULL, y="Milo log2 fold change (condition vs reference)", title=title, subtitle="one point per neighbourhood (assigned to its dominant cell type); label: significant / total neighbourhoods")
  ggsave(file, p, width=8, height=1.2 + 0.45 * nlevels(cell.groups), dpi=120)
  invisible(p)
}

## ---- run one route ----
run_route <- function(route) {
  milo <- Milo(sce)
  if (route == "conos") {
    t0 <- tic()
    g2 <- g
    if (!is.null(edge_attr(g2, "weight"))) g2 <- delete_edge_attr(g2, "weight")
    g2 <- simplify(as_undirected(g2, mode="collapse"))
    stopifnot(identical(V(g2)$name, colnames(milo)))
    graph(milo) <- g2
    milo@.k <- k.conos
    toc(paste0(route, "_build_graph"), t0)
    rd.name <- "UMAP"
    t0 <- tic(); set.seed(1)
    milo <- makeNhoods(milo, prop=par$prop, refined=TRUE, refinement_scheme="graph")
    toc(paste0(route, "_makeNhoods"), t0)
  } else {
    t0 <- tic()
    suppressMessages(library(Seurat))
    so <- CreateSeuratObject(counts(sce), min.cells=0, min.features=0)
    so <- NormalizeData(so, verbose=FALSE)
    so <- FindVariableFeatures(so, nfeatures=2000, verbose=FALSE)
    so <- ScaleData(so, features=VariableFeatures(so), verbose=FALSE)
    so <- RunPCA(so, features=VariableFeatures(so), npcs=d.pca, verbose=FALSE)
    pca <- Embeddings(so, "pca"); stopifnot(identical(rownames(pca), colnames(milo)))
    reducedDim(milo, "PCA") <- pca
    rm(so); invisible(gc())
    toc(paste0(route, "_seurat_pca"), t0)
    t0 <- tic()
    milo <- buildGraph(milo, k=k.pca, d=d.pca, reduced.dim="PCA")
    toc(paste0(route, "_buildGraph"), t0)
    rd.name <- "PCA"
    t0 <- tic(); set.seed(1)
    milo <- makeNhoods(milo, prop=par$prop, k=k.pca, d=d.pca, refined=TRUE, reduced_dims="PCA", refinement_scheme="reduced_dim")
    toc(paste0(route, "_makeNhoods"), t0)
  }
  nh <- nhoods(milo)
  nh.size <- colSums(nh)
  message(sprintf("[%s/%s] %d neighbourhoods; size: median %.0f, IQR %.0f-%.0f, range %d-%d; cells covered: %d / %d",
                  ds, route, ncol(nh), median(nh.size), quantile(nh.size, .25), quantile(nh.size, .75), min(nh.size), max(nh.size), sum(rowSums(nh) > 0), nrow(nh)))
  t0 <- tic()
  milo <- countCells(milo, meta.data=as.data.frame(colData(milo)), samples="sample")
  toc(paste0(route, "_countCells"), t0)
  design.df <- sm[colnames(nhoodCounts(milo)), , drop=FALSE]
  runs <- list()
  for (dn in names(par$designs)) {
    res.all <- NULL
    for (w in weightings) {
      t0 <- tic()
      res <- testNhoods(milo, design=par$designs[[dn]], design.df=design.df, fdr.weighting=w, reduced.dim=rd.name)
      toc(paste0(route, "_", dn, "_testNhoods_", w), t0)
      res <- res[order(res$Nhood), ]
      if (is.null(res.all)) res.all <- res
      res.all[[paste0("SpatialFDR.", gsub("-", "_", w))]] <- res$SpatialFDR
    }
    res <- res.all; res$SpatialFDR <- res$SpatialFDR.graph_overlap   # primary weighting
    res <- annotateNhoods(milo, res, coldata_col="celltype")
    res$celltype_fraction <- res$celltype_fraction; res$nhood.size <- nh.size[res$Nhood]
    res$index.cell <- colnames(milo)[unlist(nhoodIndex(milo))][res$Nhood]
    cdf <- map_cells(nh, res)
    res.k <- res; res.k$SpatialFDR <- res$SpatialFDR.k_distance; cdf.k <- map_cells(nh, res.k)
    summ <- summarise_ct(res, cdf); summ.k <- summarise_ct(res.k, cdf.k)
    rn <- paste0(route, "|", dn)
    runs[[rn]] <- list(route=route, design=deparse(par$designs[[dn]]), fdr.weighting.primary="graph-overlap",
                       n.nhoods=nrow(res), nhoods=res, cells=cdf, summary=summ,
                       cells.kdist=cdf.k, summary.kdist=summ.k)
    message(sprintf("[%s/%s] design %s: %d nhoods, SpatialFDR<0.1: graph-overlap %d (up %d / down %d); k-distance %d",
                    ds, route, deparse(par$designs[[dn]]), nrow(res), sum(res$SpatialFDR < 0.1), sum(res$SpatialFDR < 0.1 & res$logFC > 0), sum(res$SpatialFDR < 0.1 & res$logFC < 0), sum(res$SpatialFDR.k_distance < 0.1)))
    print(summ[, c("celltype", "n.nhoods", "nh.fdr10.up", "nh.fdr10.down", "nh.fdr05.up", "nh.fdr05.down", "n.cells", "n.cells.no.nhood", "cells.fdr10.up", "cells.fdr10.down")])
    suffix <- if (route == "conos" && dn == "cond") "" else paste0("_", route, "_", dn)
    beeswarm_plot(res, sprintf("Milo, %s dataset: %s route, design %s", ds, route, deparse(par$designs[[dn]])), sprintf("milo_%s_beeswarm%s.png", ds, suffix))
  }
  list(runs=runs, nhood.matrix=nh, nhood.size=nh.size, k=if (route == "conos") k.conos else k.pca, prop=par$prop,
       graph.degree=summary(degree(graph(milo))))
}

out <- list(dataset=ds, params=list(prop=par$prop, k.conos=k.conos, k.pca=k.pca, d.pca=d.pca, cond.col=par$cond.col, ref=par$ref, weightings=weightings,
                                    miloR=as.character(packageVersion("miloR"))),
            sample.meta=sm, primary="conos|cond", runs=list(), routes=list())
for (route in c("conos", "pca")) {
  rr <- run_route(route)
  out$runs <- c(out$runs, rr$runs)
  out$routes[[route]] <- rr[c("nhood.matrix", "nhood.size", "k", "prop", "graph.degree")]
  saveRDS(out, paste0("milo_", ds, ".rds"))
}
out$timing <- timing
saveRDS(out, paste0("milo_", ds, ".rds"))
message("timing (s):"); print(round(unlist(timing), 1))
message("done ", ds)
