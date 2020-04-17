#' @title Prepare pathway data
#' @description  Filter and prepare DE genes for onthology calculations
#' @param de List with differentially expressed genes per cell group
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param OrgDB Genome-wide annotation (default=org.Hs.eg.db)
#' @param stat.cutoff Cutoff for filtering highly-expressed DE genes (default=3)
#' @param verbose Print progress (default=T)
#' @return A list containing DE gene IDs, filtered DE genes, and input DE genes
#' @export
preparePathwayData <- function(cms, de, cell.groups, transpose=T, OrgDB=org.Hs.eg.db, verbose=T, stat.cutoff=3) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("You have to install 'clusterProfiler' package to perform onthology analysis")

  if(verbose) cat("Merging count matrices ... ")

  # TODO shouldn't depend on Conos?
  cm_merged <- conos:::mergeCountMatrices(cms)

  if(transpose) {
    if(verbose) cat("done!\nTransposing merged count matrix ... ")
    cm_merged %<>% Matrix::t()
  }

  if(verbose) cat("done!\nCalculating boolean index ... ")
  cm_bool <- (cm_merged > 1) * 1

  if(verbose) cat("done!\nCollapsing cells ... ")
  # TODO Shouldn't depend on Conos?
  cm_collapsed_bool <- conos:::collapseCellsByType(Matrix::t(cm_bool), cell.groups %>%
                                                     .[. %in% names(de)] %>%
                                                     factor, min.cell.count=0)

  if(verbose) cat("done!\nFiltering DE genes .. ")
  de.filtered <- lapply(de, function(df) df[!is.na(df$stat) & (abs(df$stat) > stat.cutoff),]) # Consider filtering by pAdj

  if(verbose) cat(". ")
  de.genes.filtered <- mapply(intersect, lapply(de.filtered, rownames), ((cm_collapsed_bool > as.vector(table(annotation %>% .[. %in% names(de)]%>% factor)[rownames(cm_collapsed_bool)] * 0.05)) %>% apply(1, function(row) names(which(row))))[names(de.filtered)]) # Consider 0.05

  if(verbose) cat("done!\nRetrieving Entrez Gene IDs ... ")
  de.gene.ids <- lapply(de.genes.filtered, clusterProfiler::bitr, 'SYMBOL', 'ENTREZID', OrgDB) %>%
    lapply(`[[`, "ENTREZID")

  if(verbose) cat("done!\nAll done!\n")

  return(list(de.gene.ids = de.gene.ids,
              de.genes.filtered = de.genes.filtered,
              de.raw = de))
}

# TODO update description
#' @title
#' @description
#' @param
#' @param
#' @param
#' @param
#' @param
#' @return
#' @export
enrichGOOpt <- function (gene, OrgDB, goData, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05,
                         pAdjustMethod = "BH", universe=NULL, qvalueCutoff = 0.2, minGSSize = 10,
                         maxGSSize = 500, readable = FALSE, pool = FALSE) {
  ont %<>% toupper %>% match.arg(c("BP", "CC", "MF"))

  res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff,
                                             pAdjustMethod = pAdjustMethod, universe = universe,
                                             qvalueCutoff = qvalueCutoff, minGSSize = minGSSize,
                                             maxGSSize = maxGSSize, USER_DATA = goData)
  if (is.null(res))
    return(res)

  res@keytype <- keyType
  res@organism <- clusterProfiler:::get_organism(OrgDB)
  if (readable) {
    res <- DOSE::setReadable(res, OrgDB)
  }
  res@ontology <- ont

  return(res)
}

#' @title Estimate onthology
#' @description  Calculate onthologies based on DEs
#' @param type Onthology type, either GO (gene onthology) or DO (disease onthology). Please see DOSE package for more information.
#' @param pathway.data List containing DE gene IDs, and filtered and unfiltered DE genes
#' @param OrgDB Genome-wide annotation (default=org.Hs.eg.db)
#' @param p.adj Adjusted P cutoff (default=0.05)
#' @param p.adjust.method Method for calculating adj. P. Please see DOSE package for more information (default="BH")
#' @param readable Mapping gene ID to gene name (default=T)
#' @param n.cores Number of cores used (default: stored vector)
#' @param verbose Print progress (default=T)
#' @param ... Additional parameters for sccore:::plapply function
#' @return A list containing a list of onthologies per type of onthology, and a data frame with merged results
#' @export
estimateOnthology <- function(type, pathway.data, OrgDB=org.Hs.eg.db, p.adj=0.05, p.adjust.method="BH", readable=T, verbose=T, ...) {
  if(type=="DO") {
    # TODO enable mapping to human genes for non-human data https://support.bioconductor.org/p/88192/
    ont.list <- sccore:::plapply(pathway.data$de.gene.ids, DOSE::enrichDO, pAdjustMethod=p.adjust.method, readable=readable, n.cores=1, progress=verbose, ...) %>%
      lapply(function(x) x@result)
    ont.list %<>% names %>%
      setNames(., .) %>%
      lapply(function(n) mutate(ont.list[[n]], Type=n)) %>%
      lapply(function(x) filter(x, p.adjust < p.adj))

    ont.df <- ont.list %>%
      .[sapply(., nrow) > 0] %>%
      dplyr::bind_rows() %>%
      dplyr::select(Type, ID, Description, GeneRatio, geneID, pvalue, p.adjust, qvalue)
  } else if(type=="GO") {
    if(verbose) cat("Extracting environment data ... \n")
    go_data <- c("BP", "CC", "MF") %>%
      setNames(., .) %>%
      sccore:::plapply(function(n) clusterProfiler:::get_GO_data(OrgDB, n, "ENTREZID") %>%
                         as.list() %>%
                         as.environment(), n.cores=1, progress=verbose)
    if(verbose) cat("Estimating enriched pathways ... \n")
    ont.list <- names(go_data) %>%
      setNames(., .) %>%
      lapply(function(ont) sccore:::plapply(pathway.data$de.gene.ids, enrichGOOpt, ont=ont, goData=go_data[[ont]], readable=readable, pAdjustMethod=p.adjust.method, OrgDB=OrgDB, n.cores=1, progress=verbose, ...)) %>%
      lapply(lapply, function(x) x@result)

    ont.list %<>% lapply(lapply, function(x) filter(x, p.adjust < p.adj)) %>%
      lapply(function(gt) gt %>%
               .[sapply(., nrow) > 0] %>%
               names() %>%
               setNames(., .) %>%
               lapply(function(n) cbind(gt[[n]], Type=n)) %>%
               Reduce(rbind, .))

    ont.df <- ont.list %>%
      names %>%
      lapply(function(n) ont.list[[n]] %>%
               mutate(GO=n) %>%
               dplyr::select(GO, Type, ID, Description, GeneRatio, geneID, pvalue, p.adjust, qvalue)) %>%
      dplyr::bind_rows()
  } else {
    stop("'type' must be either 'GO' or 'DO'.")
  }

  return(list(list=ont.list,
              df=ont.df))
}

#' @title Distance between terms
#' @description Calculate distance matrix between onthology terms
#' @param type Onthology, must be either "GO" or "DO" (default=NULL)
#' @param pathway.data Results from preparePathwayData (default: stored list)
#' @return Distance matrix
#' @export
distanceBetweenTerms <- function(type=NULL, ont.res) {
  genes.per.go <- sapply(ont.res$geneID, strsplit, "/") %>% setNames(ont.res$Description)
  all.go.genes <- unique(unlist(genes.per.go))
  all.gos <- unique(ont.res$Description)

  genes.per.go.mat <- matrix(0, length(all.go.genes), length(all.gos)) %>%
    `colnames<-`(all.gos) %>% `rownames<-`(all.go.genes)

  for (i in 1:length(genes.per.go)) {
    genes.per.go.mat[genes.per.go[[i]], ont.res$Description[[i]]] <- 1
  }

  return(dist(t(genes.per.go.mat), method="binary"))
}

#' @title Get onthology summary
#' @description Get summary
#' @param type Onthology, must be either "BP", "CC", or "MF" (GO types) or "DO" (default=NULL)
#' @param pathway.data Results from preparePathwayData (default: stored list)
#' @return Data frame
#' @export
getOnthologySummary <- function(type=NULL, ont.res) {
  go_dist <- distanceBetweenTerms(type, ont.res)
  clusts <- hclust(go_dist) %>%
    cutree(h=0.75)

  ont.res %<>% mutate(Clust=clusts[Description])

  name_per_clust <- ont.res %>%
    group_by(Clust, Description) %>%
    summarise(pvalue=exp(mean(log(pvalue)))) %>%
    split(.$Clust) %>%
    sapply(function(df) df$Description[which.min(df$pvalue)])

  ont.res %<>% mutate(ClustName=name_per_clust[as.character(Clust)])

  order.anno <- ont.res$Type %>% unique %>% .[order(.)]

  df <- ont.res %>%
    group_by(Type, ClustName) %>%
    summarise(p.adjust=min(p.adjust)) %>%
    ungroup %>%
    mutate(p.adjust=-log10(p.adjust)) %>%
    tidyr::spread(Type, p.adjust) %>%
    as.data.frame() %>%
    set_rownames(.$ClustName) %>%
    .[, 2:ncol(.)] %>%
    .[, order.anno[order.anno %in% colnames(.)]]

  df[is.na(df)] <- 0

  return(df)
}
