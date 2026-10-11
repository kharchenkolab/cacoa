## Light-weight PF inputs (Habermann 2020; Control and IPF lungs) in the benchmark format. Usage: Rscript extract_pf.R
suppressPackageStartupMessages({library(Matrix); library(igraph); library(dplyr); library(readr)})
S <- Sys.getenv("CACOA_SCRATCH", "/tmp/claude-1001/-home-pkharchenko-cacoa-dev/7bdc188e-5fb7-49ad-9052-73aa7139e530/scratchpad")
t0 <- Sys.time()
con <- readRDS(file.path(S, "ex_pf/pf_con_refreshed.rds")); cat("loaded", format(Sys.time() - t0), "\n")
cell.groups <- con$misc$cell_metadata$cellType
lib.meta <- read_csv("/home/viktor_petukhov/projects/cacoaAnalysis/data/PF/sample_metadata.csv", show_col_types = FALSE, na = c("", "NA", "Unknown")) %>% filter(Diagnosis %in% c("Control", "IPF"))
cell.meta <- read.csv("/home/viktor_petukhov/projects/cacoaAnalysis/data/PF/cell_metadata.csv", row.names = 1) %>% filter(Diagnosis %in% c("Control", "IPF"))
sample.meta <- lib.meta %>% group_by(Sample_Name) %>%
  summarise(Diagnosis = first(Diagnosis), Processing_site = first(Processing_site), Gender = first(Gender), Age = as.integer(first(Age))) %>%
  tibble::column_to_rownames("Sample_Name")
sample.meta <- sample.meta[names(con$samples), ]
spc <- con$getDatasetPerCell(); 
g <- con$graph; cells <- V(g)$name
cm <- con$getJointCountMatrix(raw = TRUE)
cat("graph", vcount(g), "cm", dim(cm), "\n")
out <- list(graph = g, embedding = con$embedding, counts = cm, cell.groups = cell.groups, sample.per.cell = spc, sample.meta = sample.meta)
saveRDS(out, file.path(S, "bench/inputs_pf.rds"))
cat("saved; total", format(Sys.time() - t0), "\n")
