## cacoa (BH / max-t) versus Milo (conos route, ~condition, graph-overlap spatial FDR): per cell type and maps
suppressPackageStartupMessages({library(ggplot2)})
S <- "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad"
args <- commandArgs(TRUE); ds <- args[1]
ca <- readRDS(file.path(S, sprintf("milo/cacoa_%s.rds", ds))); mi <- readRDS(file.path(S, sprintf("milo/milo_%s.rds", ds)))
run <- mi$runs[[mi$primary]]; cells <- run$cells; nh <- run$nhoods
bh <- ca$bh; mt <- ca$maxt
m <- merge(bh[, c("cell", "celltype", "z.adj", "z.raw", "effect")], cells[, c("cell", "milo.logFC", "milo.fdr")], by = "cell")
m$maxt <- mt$z.adj[match(m$cell, mt$cell)]
summ <- do.call(rbind, lapply(split(m, m$celltype), function(d) data.frame(
  celltype = d$celltype[1], n = nrow(d),
  cacoa.bh = mean(abs(d$z.adj) > qnorm(0.975), na.rm = TRUE), cacoa.bh.dir = sign(median(d$z.adj, na.rm = TRUE)),
  cacoa.maxt = mean(abs(d$maxt) > qnorm(0.975), na.rm = TRUE),
  milo.cells10 = mean(d$milo.fdr < 0.1, na.rm = TRUE), milo.cells05 = mean(d$milo.fdr < 0.05, na.rm = TRUE), milo.dir = sign(median(d$milo.logFC, na.rm = TRUE)),
  milo.nh10 = { k <- nh$celltype == d$celltype[1]; if (any(k)) mean(nh$SpatialFDR[k] < 0.1) else NA },
  cor.z.logFC = suppressWarnings(cor(d$z.raw, d$milo.logFC, method = "spearman", use = "complete.obs")), stringsAsFactors = FALSE)))
rownames(summ) <- NULL
cat("## ", ds, ": cells", nrow(m), " neighbourhoods", nrow(nh), "\n")
print(summ, digits = 2, row.names = FALSE)
cat(sprintf("cell-level Spearman(cacoa raw z, Milo logFC) = %.2f; cells significant: cacoa BH %.2f, max-t %.2f, Milo FDR<0.1 %.2f, FDR<0.05 %.2f\n",
            cor(m$z.raw, m$milo.logFC, method = "spearman", use = "complete.obs"), mean(abs(m$z.adj) > 1.96, na.rm = TRUE), mean(abs(m$maxt) > 1.96, na.rm = TRUE),
            mean(m$milo.fdr < 0.1, na.rm = TRUE), mean(m$milo.fdr < 0.05, na.rm = TRUE)))
saveRDS(list(summary = summ, cells = m), file.path(S, sprintf("milo/compare_%s.rds", ds)))
## maps
emb <- readRDS(file.path(S, sprintf("milo/inputs_%s.rds", ds)))$embedding
df <- data.frame(x = emb[m$cell, 1], y = emb[m$cell, 2], m)
df$milo.sig <- ifelse(is.finite(df$milo.fdr) & df$milo.fdr < 0.1, df$milo.logFC, NA)
lim <- quantile(abs(c(df$z.adj, df$maxt)), 0.99, na.rm = TRUE)
pal <- function(v, title) ggplot2::scale_colour_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b", limits = c(-lim, lim), oob = scales::squish, name = title, na.value = "grey85")
base <- function(col, title, legend) ggplot(df[order(abs(df[[col]]), na.last = FALSE), ], aes(x, y, colour = .data[[col]])) + geom_point(size = 0.3, alpha = 0.7) + theme_void() +
  ggtitle(title) + theme(plot.title = element_text(size = 11))
g1 <- base("z.adj", "cacoa density, BH-adjusted z", "z") + pal(df$z.adj, "z adj.")
g2 <- base("maxt", "cacoa density, max-t adjusted z", "z") + pal(df$maxt, "z adj.")
liml <- quantile(abs(df$milo.sig), 0.99, na.rm = TRUE)
g3 <- base("milo.sig", "Milo: logFC of neighbourhoods with spatial FDR < 0.1", "logFC") + ggplot2::scale_colour_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b", limits = c(-liml, liml), oob = scales::squish, name = "logFC", na.value = "grey85")
ggsave(file.path(S, sprintf("milo/compare_%s_maps.png", ds)), cowplot::plot_grid(g1, g2, g3, ncol = 3), width = 15, height = 4.6, dpi = 90, bg = "white")
cat("maps saved\n")
