#' @import dplyr
NULL

#' @title Prepare onthology data
#' @description  Filter and prepare DE genes for onthology calculations
#' @param cms list of counts matrices; column for gene and row for cell
#' @param de List with differentially expressed genes per cell group
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param org Organism, can be "human", "mouse", "zebrafish", "worm", or "fly" (default="human")
#' @param stat.cutoff Cutoff for filtering highly-expressed DE genes (default=3)
#' @param verbose Print progress (default=T)
#' @param n.cores Number of cores to use (default=1)
#' @return A list containing DE gene IDs, filtered DE genes, and input DE genes
#' @export
prepareOnthologyData <- function(cms, de.raw, cell.groups, universe = NULL, org = "human", n.top.genes = Inf, stat.cutoff = 3, transposed = T,  verbose = T, n.cores = 1) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("You have to install 'clusterProfiler' package to perform onthology analysis")

  # TODO Test functionality for other species than human or mouse
  if(org == "human") {
    if(!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("You have to install 'org.Hs.eg.db' package to perform onthology analysis")
    require(org.Hs.eg.db)
    OrgDB = org.Hs.eg.db
  } else if(org == "mouse") {
    if(!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("You have to install 'org.Mm.eg.db' package to perform onthology analysis")
    require(org.Mm.eg.db)
    OrgDB = org.Mm.eg.db
  } else if(org == "zebrafish") {
    if(!requireNamespace("org.Dm.eg.db", quietly = TRUE)) stop("You have to install 'org.Dm.eg.db' package to perform onthology analysis")
    require(org.Dm.eg.db)
    OrgDB = org.Dm.eg.db
  } else if(org == "worm") {
    if(!requireNamespace("org.Ce.eg.db", quietly = TRUE)) stop("You have to install 'org.Ce.eg.db' package to perform onthology analysis")
    require(org.Ce.eg.db)
    OrgDB = org.Ce.eg.db
  } else if(org == "fly") {
    if(!requireNamespace("org.Dr.eg.db", quietly = TRUE)) stop("You have to install 'org.Dr.eg.db' package to perform onthology analysis")
    require(org.Dr.eg.db)
    OrgDB = org.Dr.eg.db
  }

  if(verbose) cat("Merging count matrices ... ")

  cm_merged <- sccore:::mergeCountMatrices(cms, transposed = transposed, n.cores = n.cores)

  if(!transposed) {
    if(verbose) cat("done!\nTransposing merged count matrix ... ")
    cm_merged %<>% Matrix::t()
  }

  if(verbose) cat("done!\nCalculating boolean index ... ")
  cm_bool <- (cm_merged > 1) * 1

  if(verbose) cat("done!\nFiltering DE genes .")
  cm_collapsed_bool <- collapseCellsByType(cm_bool, cell.groups %>%
                                                     .[. %in% names(de.raw)] %>%
                                                     factor, min.cell.count=0)

  if(verbose) cat(".")
  de.filtered <- lapply(de.raw, function(df) df[!is.na(df$stat) & (abs(df$stat) > stat.cutoff),]) # Consider filtering by pAdj

  if(verbose) cat(". ")
  de.genes.filtered <- mapply(intersect, lapply(de.filtered, rownames), ((cm_collapsed_bool > as.vector(table(cell.groups %>% .[. %in% names(de.raw)]%>% factor)[rownames(cm_collapsed_bool)] * 0.05)) %>% apply(1, function(row) names(which(row))))[names(de.filtered)])

  if(verbose) cat("done!\nRetrieving Entrez Gene IDs ... ")
  groups <- de.genes.filtered %>% .[sapply(., length) > 0] %>% names

  de.gene.ids <- groups %>% lapply(function(x) {
    de <- de.raw[[x]] %>% .[rownames(.) %in% de.genes.filtered[[x]],]

    list(down = de %>% .[order(.$Z[.$Z < 0], decreasing = F),] %>% rownames,
         up = de %>% .[order(.$Z[.$Z > 0], decreasing = F),] %>% rownames,
         all = de %>% .[order(.$Z %>% abs(), decreasing = F),] %>% rownames) %>%
      lapply(function(l) if(n.top.genes == Inf | n.top.genes > length(l)) l else l[1:n.top.genes]) %>%
      append(list(universe=rownames(de.raw[[x]])))
    }) %>% setNames(groups)

  de.gene.ids %<>% sccore:::plapply(function(id) suppressMessages(lapply(id, clusterProfiler::bitr, 'SYMBOL', 'ENTREZID', OrgDB)), n.cores = 1, progress = verbose) %>%
    lapply(lapply, `[[`, "ENTREZID")

  if(verbose) cat("All done!\n")

  return(list(de.gene.ids = de.gene.ids,
              de.filter = de.genes.filtered))
}

enrichGOOpt <- function(gene, OrgDB, go.environment, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", universe=NULL, qvalueCutoff = 0.2, minGSSize = 10,
                        maxGSSize = 500, readable = FALSE, pool = FALSE) {
  ont %<>% toupper %>% match.arg(c("BP", "CC", "MF"))

  res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff,
                                             pAdjustMethod = pAdjustMethod, universe = universe,
                                             qvalueCutoff = qvalueCutoff, minGSSize = minGSSize,
                                             maxGSSize = maxGSSize, USER_DATA = go.environment)
  if (is.null(res))
    return(res)

  res@keytype <- keyType
  res@organism <- clusterProfiler:::get_organism(OrgDB)
  if (readable) {
    res <- DOSE::setReadable(res, OrgDB)
  }
  res@ontology <- ont
  if(!is.null(res)) return(res)
}

estimateEnrichedGO <- function(de.gene.ids, go.environment, ...) {
  ont.list <- names(go.environment) %>%
    sccore:::sn() %>%
    lapply(function(ont) lapply(de.gene.ids, enrichGOOpt, go.environment=go.environment[[ont]], ont=ont, ...)) %>%
    lapply(lapply, function(x) x@result)
}

#' @title Distance between terms
#' @description Calculate distance matrix between onthology terms
#' @param type Onthology, must be either "GO" or "DO" (default=NULL)
#' @param ont.res Results from prepareOnthologyData (default: stored list)
#' @return Distance matrix
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
#' @param ont.res Results from prepareOnthologyData (default: stored list)
#' @return Data frame
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

  order.anno <- ont.res$Group %>% unique %>% .[order(.)]

  df <- ont.res %>%
    group_by(Group, ClustName) %>%
    summarise(p.adjust=min(p.adjust)) %>%
    ungroup %>%
    mutate(p.adjust=-log10(p.adjust)) %>%
    tidyr::spread(Group, p.adjust) %>%
    as.data.frame() %>%
    set_rownames(.$ClustName) %>%
    .[, 2:ncol(.)] %>%
    .[, order.anno[order.anno %in% colnames(.)]]

  df[is.na(df)] <- 0

  return(df)
}

filterOnthologies <- function(ont.list, p.adj) {
  ont.list %>%
    names() %>%
    lapply(function(dir) {
      lapply(ont.list[[dir]] %>% names(), function(group) {
        dplyr::mutate(ont.list[[dir]][[group]], Group=group) %>%
          dplyr::filter(p.adjust < p.adj)
      }) %>%
        setNames(ont.list[[dir]] %>% names) %>%
        .[sapply(., nrow) > 0] # Remove empty data frames
    }) %>%
    setNames(c("down", "up", "all"))
}

onthologyListToDf <- function(ont.list) {
  lapply(ont.list, function(x) {
    if(length(x) > 0) {
      dplyr::bind_rows(x) %>%
        {if("Type" %in% colnames(.)) {
          dplyr::select(., Group, Type, ID, Description, GeneRatio, geneID, pvalue, p.adjust, qvalue)
        } else {
          dplyr::select(., Group, ID, Description, GeneRatio, geneID, pvalue, p.adjust, qvalue)
        }
        }
    } else {
      "No significant onthologies identified. Try relaxing p.adj."
    }
  })
}

#' @title Estimate onthology
#' @description  Calculate onthologies based on DEs
#' @param type Onthology type, either GO (gene onthology) or DO (disease onthology). Please see DOSE package for more information.
#' @param ont.data List containing DE gene IDs, and filtered DE genes
#' @param org Organism, can be "human", "mouse", "zebrafish", "worm", or "fly" (default="human")
#' @param p.adj Adjusted P cutoff (default=0.05)
#' @param p.adjust.method Method for calculating adj. P. Please see DOSE package for more information (default="BH")
#' @param readable Mapping gene ID to gene name (default=T)
#' @param n.cores Number of cores used (default: stored vector)
#' @param verbose Print progress (default=T)
#' @param ... Additional parameters for sccore:::plapply function
#' @return A list containing a list of onthologies per type of onthology, and a data frame with merged results
#' @export
estimateOnthology <- function(type, org="human", de.gene.ids, go.environment = NULL, p.adj=0.05, p.adjust.method="BH", readable=T, verbose=T, n.cores = 1, ...) {
  if(org == "human") {
    if(!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("You have to install 'org.Hs.eg.db' package to perform onthology analysis")
    require(org.Hs.eg.db)
    OrgDB = org.Hs.eg.db
  } else if(org == "mouse") {
    if(!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("You have to install 'org.Mm.eg.db' package to perform onthology analysis")
    require(org.Mm.eg.db)
    OrgDB = org.Mm.eg.db
  } else if(org == "zebrafish") {
    if(!requireNamespace("org.Dm.eg.db", quietly = TRUE)) stop("You have to install 'org.Dm.eg.db' package to perform onthology analysis")
    require(org.Dm.eg.db)
    OrgDB = org.Dm.eg.db
  } else if(org == "worm") {
    if(!requireNamespace("org.Ce.eg.db", quietly = TRUE)) stop("You have to install 'org.Ce.eg.db' package to perform onthology analysis")
    require(org.Ce.eg.db)
    OrgDB = org.Ce.eg.db
  } else if(org == "fly") {
    if(!requireNamespace("org.Dr.eg.db", quietly = TRUE)) stop("You have to install 'org.Dr.eg.db' package to perform onthology analysis")
    require(org.Dr.eg.db)
    OrgDB = org.Dr.eg.db
  }

  dir.names <- c("down", "up", "all")

  if(type=="DO") {
    # TODO enable mapping to human genes for non-human data https://support.bioconductor.org/p/88192/
    if(org != "human") stop("Only human data supported for DO analysis.")
    ont.list <- sccore:::plapply(names(de.gene.ids), function(id) lapply(de.gene.ids[[id]][1:3], DOSE::enrichDO, pAdjustMethod=p.adjust.method, universe=de.gene.ids[[id]][["universe"]], readable=readable), n.cores=1, progress=verbose, ...) %>%
      lapply(lapply, function(x) x@result) %>%
      setNames(names(de.gene.ids))

    # Split into different fractions
    ont.list <- 1:3 %>%
      lapply(function(x) lapply(ont.list, `[[`, x)) %>%
      setNames(dir.names)

    ont.list %<>% filterOnthologies(p.adj = p.adj)

    ont.df <- ont.list %>% onthologyListToDf()

    res <- list(list=ont.list,
                df=ont.df)
  } else if(type=="GO") {
    if(is.null(go.environment)) {
      if(verbose) cat("Extracting environment data ... \n")
      go.environment <- c("BP", "CC", "MF") %>%
        sccore:::sn() %>%
        sccore:::plapply(function(n) clusterProfiler:::get_GO_data(OrgDB, n, "ENTREZID") %>%
                           as.list() %>%
                           as.environment(), n.cores=1, progress=verbose)
    } else {
      message("Using stored GO environment. To extract again, rerun with 'go.environment = NULL'.")
    }

    if(verbose) cat("Estimating enriched onthologies ... \n")
    ont.list <- sccore:::plapply(names(de.gene.ids), function(id) {
      estimateEnrichedGO(de.gene.ids[[id]][1:3], go.environment = go.environment, universe=de.gene.ids[[id]][["universe"]], readable=readable, pAdjustMethod=p.adjust.method, OrgDB=OrgDB)
    }, progress = verbose, n.cores = n.cores, ...) %>%
      setNames(names(de.gene.ids))

    #Split into different fractions
    ont.list <- 1:3 %>%
      lapply(function(x) lapply(ont.list, lapply, `[[`, x)) %>%
      setNames(dir.names)

    ont.list %<>%
      names() %>%
      lapply(function(dir) {
        lapply(ont.list[[dir]] %>% names, function(group) {
          lapply(ont.list[[dir]][[group]] %>% names, function(go) {
            dplyr::mutate(ont.list[[dir]][[group]][[go]], Type = go)
          }) %>% setNames(ont.list[[dir]][[group]] %>% names) %>%
            dplyr::bind_rows()
        }) %>% setNames(ont.list[[dir]] %>% names)
      }) %>% setNames(dir.names)

    ont.list %<>% filterOnthologies(p.adj = p.adj)

    ont.df <- ont.list %>% onthologyListToDf()

    res <- list(list=ont.list,
                df=ont.df,
                go.environment=go.environment)
  } else {
    stop("'type' must be either 'GO' or 'DO'.")
  }

  return(res)
}
