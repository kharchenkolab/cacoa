#' @import magrittr
NULL

##' @title Perform Gene Set Enrichment Analysis
##' @param de.raw List with DE results
##' @param org.db Organism Db (default = \link[org.Hs.eg.db:org.Hs.eg.db]{org.Hs.eg.db} if installed)
##' @param verbose Show progress (default=T)
##' @param n.cores Number of cores (default=1)
##' @param ... Additional parameters for \link[clusterProfiler:gseGO]{clusterProfiler::gseGO}
estimateGSEA <- function(de.raw, org.db = NULL, verbose = TRUE, n.cores = 1, ...) {
  if (is.null(org.db)) {
    if (!requireNamespace("org.Hs.eg.db", quietly=TRUE))
      stop("org.db must be provided")

    org.db <- org.Hs.eg.db::org.Hs.eg.db
  }

  if (!requireNamespace("clusterProfiler", quietly=TRUE))
    stop("clusterProfiler must be installed to use this function")

  res <- sccore::plapply(de.raw, function(de) {
    gene.ids <- suppressWarnings(suppressMessages(clusterProfiler::bitr(rownames(de), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.db)))
    gene.ids <- de$log2FoldChange[match(gene.ids$SYMBOL, rownames(de))] %>%
      setNames(gene.ids$ENTREZID) %>%
      .[order(., decreasing = T)]

    ont.names <- c("BP","MF","CC")

    ont.names %>%
      lapply(function(ont) {
        suppressWarnings(suppressMessages(clusterProfiler::gseGO(geneList = gene.ids, OrgDb = org.db, ont = ont, eps = 0, ...))) %>%
          clusterProfiler::setReadable(org.db, 'ENTREZID') %>%
          .@result
      }) %>%
      setNames(ont.names)
  }, progress = verbose, n.cores = n.cores) %>%
    lapply(function(x) x[sapply(x, nrow) > 0]) %>% # Remove empty entries on ont. level
    .[sapply(., length) > 0] # Remove empty entries on group level

  if(verbose) {
    dif <- setdiff(names(de.raw), names(res))
    if(length(dif) > 0) {
      message(paste0("\nNo significant GSEA terms identified for ",length(dif)," cell group(s):"))
      print(dif)
    }
  }

  res %>%
    names() %>%
    lapply(function(group) {
      names(res[[group]]) %>%
        lapply(function(ont) dplyr::mutate(res[[group]][[ont]], Type = ont, Group = group)) %>%
        dplyr::bind_rows()
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::rename(geneID = `core_enrichment`)
}
