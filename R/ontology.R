#' @import dplyr
NULL

#' @title Prepare ontology data
#' @description  Filter and prepare DE genes for ontology calculations
#' @param cms list of counts matrices; column for gene and row for cell
#' @param de.raw List with differentially expressed genes per cell group
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
#' @param n.top.genes Number of most different genes to take as input. If less are left after filtering for p.adj.cutoff, additional genes are included. To disable, set n.top.genes=0 (default=1e2)
#' @param p.adj.cutoff Adj. P cutoff for filtering DE genes (default=0.05)
#' @param expr.cutoff Cutoff for cells per group expressing a DE gene, i.e., cutoff for highly-expressed genes (default=0.05)
#' @param universe Background gene list, NULL will take all input genes in de.raw (default=NULL)
#' @param transposed Indicate whether raw counts are transposed, i.e., cells as rows, genes as columns (default=T)
#' @param verbose Print progress (default=T)
#' @param n.cores Number of cores to use (default=1)
#' @return A list containing DE ENSEMBL gene IDs, and filtered DE genes
#' @export
prepareOntologyData <- function(cms, de.raw, cell.groups, org.db, n.top.genes = 1e2, p.adj.cutoff = 0.05, expr.cutoff = 0.05, universe = NULL, transposed = T,  verbose = T, n.cores = 1) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("You have to install 'clusterProfiler' package to perform ontology analysis")

  if(class(org.db) != "OrgDb") stop("'org.db' must be of class 'OrgDb'. Please input an organism database, e.g., org.Hs.eg.db for human data.")

  de.filtered <- lapply(de.raw, function(df) {
    df.filt <- df[!is.na(df$padj) & (df$padj < p.adj.cutoff),]
    if(nrow(df.filt) < n.top.genes) {
      df.filt <- df[1:n.top.genes,]
    }
    return(df.filt)
  })

  if(expr.cutoff != 0) {
    if(verbose) cat("Merging count matrices ... ")

    cm_merged <- sccore:::mergeCountMatrices(cms, transposed = transposed, n.cores = n.cores)

    if(!transposed) {
      if(verbose) cat("done!\nTransposing merged count matrix ... ")
      cm_merged %<>% Matrix::t()
    }

    if(verbose) cat("done!\nCalculating boolean index ... ")
    cm_bool <- (cm_merged > 1) * 1

    if(verbose) cat("done!\nFiltering DE genes ... ")
    cm_collapsed_bool <- conos:::collapseCellsByType(cm_bool, cell.groups %>%
                                                       .[. %in% names(de.raw)] %>%
                                                       factor(), min.cell.count=0)

    de.genes.filtered <- mapply(intersect,
                                lapply(de.filtered, rownames),
                                ((cm_collapsed_bool > as.vector(table(cell.groups %>%
                                                                        .[. %in% names(de.raw)] %>%
                                                                        factor)[rownames(cm_collapsed_bool)] * expr.cutoff)) %>%
                                   apply(1, function(row) names(which(row))))[names(de.filtered)])

    if(verbose) cat("done!\n")
  } else {
    de.genes.filtered <- lapply(de.filtered, rownames)
  }

  if(verbose) cat("Retrieving Entrez Gene IDs ... ")
  groups <- de.genes.filtered %>% .[sapply(., length) > 0] %>% names()

  de.gene.ids <- groups %>% lapply(function(x) {
    de <- de.raw[[x]] %>% .[rownames(.) %in% de.genes.filtered[[x]],]

    list(down = de %>% dplyr::filter(Z < 0) %>% .[order(.$Z, decreasing = F),] %>% rownames,
         up = de %>% dplyr::filter(Z > 0) %>% .[order(.$Z, decreasing = T),] %>% rownames,
         all = de %>% .[order((.$Z %>% abs()), decreasing = T),] %>% rownames) %>%
      # lapply(function(l) if(n.top.genes == Inf | n.top.genes > length(l)) l else l[1:n.top.genes]) %>%
      append(list(universe=rownames(de.raw[[x]])))
    }) %>% setNames(groups) %>%
    lapply(lapply, function(x) if(length(x)) x) %>%
    lapply(plyr::compact)

  de.gene.ids <- lapply(de.gene.ids, function(x) lapply(x, function(y) suppressMessages({tryCatch({clusterProfiler::bitr(y, 'SYMBOL', 'ENTREZID', org.db)}, error = function(err) NULL)}))) %>%
    lapply(plyr::compact) %>%
    lapply(lapply, `[[`, "ENTREZID")

  if(verbose) cat("done!\nAll done!\n")

  return(list(de.gene.ids = de.gene.ids,
              de.filter = de.genes.filtered))
}

enrichGOOpt <- function(gene, org.db, go.environment, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", universe=NULL, qvalueCutoff = 0.2, minGSSize = 5,
                        maxGSSize = 500, readable = FALSE, pool = FALSE) {
  ont %<>% toupper %>% match.arg(c("BP", "CC", "MF"))

  res <- clusterProfiler:::enricher_internal(gene, pvalueCutoff = pvalueCutoff,
                                             pAdjustMethod = pAdjustMethod, universe = universe,
                                             qvalueCutoff = qvalueCutoff, minGSSize = minGSSize,
                                             maxGSSize = maxGSSize, USER_DATA = go.environment)
  if (is.null(res))
    return(res)

  res@keytype <- keyType
  res@organism <- clusterProfiler:::get_organism(org.db)
  if (readable) {
    res <- DOSE::setReadable(res, org.db)
  }
  res@ontology <- ont
  if(!is.null(res)) return(res)
}

estimateEnrichedGO <- function(de.gene.ids, go.environment, ...) {
  ont.list <- names(go.environment) %>%
    sccore:::sn() %>%
    lapply(function(ont) lapply(de.gene.ids, enrichGOOpt, go.environment=go.environment[[ont]], ont=ont, ...)) %>%
    lapply(lapply, function(x) if(length(x)) x@result else x)
}

#' @title Distance between terms
#' @description Calculate distance matrix between ontology terms
#' @param ont.res Results from prepareOntologyData (default: stored list)
#' @return Distance matrix
distanceBetweenTerms <- function(ont.res) {
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

#' @title Get ontology summary
#' @description Get summary
#' @param ont.res Results from prepareOntologyData (default: stored list)
#' @return Data frame
getOntologySummary <- function(ont.res) {
  go_dist <- distanceBetweenTerms(ont.res)
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

filterOntologies <- function(ont.list, p.adj) {
  ont.list %>%
    names() %>%
    lapply(function(dir) {
      lapply(ont.list[[dir]] %>% names(), function(group) {
        dplyr::mutate(ont.list[[dir]][[group]], Group=group) %>%
          dplyr::arrange(p.adjust) %>%
          dplyr::filter(p.adjust < p.adj)
      }) %>%
        setNames(ont.list[[dir]] %>% names()) %>%
        .[sapply(., nrow) > 0] # Remove empty data frames
    }) %>%
    setNames(c("down", "up", "all"))
}

ontologyListToDf <- function(ont.list) {
  lapply(ont.list, function(x) {
    if(length(x) > 0) {
      dplyr::bind_rows(x) %>%
        {if("Type" %in% colnames(.)) {
          dplyr::select(., Group, Type, ID, Description, GeneRatio, geneID, pvalue, p.adjust, qvalue, Count)
        } else {
          dplyr::select(., Group, ID, Description, GeneRatio, geneID, pvalue, p.adjust, qvalue, Count)
        }
        }
    } else {
      "No significant ontologies identified. Try relaxing p.adj."
    }
  })
}

#' @title Estimate ontology
#' @description  Calculate ontologies based on DEs
#' @param type Ontology type, either GO (gene ontology) or DO (disease ontology). Please see DOSE package for more information.
#' @param ont.data List containing DE gene IDs, and filtered DE genes
#' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
#' @param p.adj Adjusted P cutoff (default=0.05)
#' @param p.adjust.method Method for calculating adj. P. Please see DOSE package for more information (default="BH")
#' @param readable Mapping gene ID to gene name (default=T)
#' @param n.cores Number of cores used (default: stored vector)
#' @param verbose Print progress (default=T)
#' @param min.genes Minimum number of input genes overlapping with ontologies (default=0)
#' @param qvalueCutoff Q value cutoff, please see clusterProfiler package for more information (default=0.2)
#' @param minGSSize Minimal geneset size, please see clusterProfiler package for more information (default=5)
#' @param maxGSSize Minimal geneset size, please see clusterProfiler package for more information (default=5e2)
#' @param ... Additional parameters for sccore:::plapply function
#' @return A list containing a list of ontologies per type of ontology, and a data frame with merged results
#' @export
estimateOntology <- function(type = "GO", org.db=NULL, de.gene.ids, go.environment = NULL, p.adj=0.05, p.adjust.method="BH", readable=T, verbose=T, min.genes = 0, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 5e2, ...) {
  dir.names <- c("down", "up", "all")

  if(type=="DO") {
    # TODO enable mapping to human genes for non-human data https://support.bioconductor.org/p/88192/
    ont.list <- sccore:::plapply(names(de.gene.ids), function(id) suppressWarnings(lapply(de.gene.ids[[id]][-length(de.gene.ids[[id]])],
                                                                                          DOSE::enrichDO,
                                                                                          pAdjustMethod=p.adjust.method,
                                                                                          universe=de.gene.ids[[id]][["universe"]],
                                                                                          readable=readable,
                                                                                          qvalueCutoff = qvalueCutoff,
                                                                                          minGSSize = minGSSize,
                                                                                          maxGSSize = maxGSSize)),
                                 n.cores=1, progress=verbose, ...) %>%
      lapply(lapply, function(x) if(length(x)) x@result else x) %>%
      setNames(names(de.gene.ids))

    # Split into different fractions
    ont.list <- dir.names %>%
      lapply(function(x) lapply(ont.list, `[[`, x)) %>%
      setNames(dir.names) %>%
      lapply(plyr::compact) #Remove empty entries

    ont.list %<>% filterOntologies(p.adj = p.adj)

    res <- list(df=ont.list %>% ontologyListToDf())
  } else if(type=="GO") {
    if(is.null(go.environment)) {
      if(class(org.db) != "OrgDb") stop("'org.db' must be of class 'OrgDb'. Please input an organism database.")

      if(verbose) cat("Extracting environment data ... \n")
      go.environment <- c("BP", "CC", "MF") %>%
        sccore:::sn() %>%
        sccore:::plapply(function(n) clusterProfiler:::get_GO_data(org.db, n, "ENTREZID") %>%
                           as.list() %>%
                           as.environment(), n.cores=1, progress=verbose)
    } else {
      message("Using stored GO environment. To extract again, rerun with 'go.environment = NULL'.")
    }

    if(verbose) cat("Estimating enriched ontologies ... \n")
    ont.list <- sccore:::plapply(names(de.gene.ids), function(id) suppressWarnings(estimateEnrichedGO(de.gene.ids[[id]][-length(de.gene.ids[[id]])],
                                                                                                      go.environment = go.environment,
                                                                                                      universe=de.gene.ids[[id]][["universe"]],
                                                                                                      readable=readable,
                                                                                                      pAdjustMethod=p.adjust.method,
                                                                                                      org.db=org.db,
                                                                                                      qvalueCutoff = qvalueCutoff,
                                                                                                      minGSSize = minGSSize,
                                                                                                      maxGSSize = maxGSSize)),
                                 progress = verbose, n.cores = 1, ...) %>%
      setNames(names(de.gene.ids))

    #Split into different fractions
    ont.list <- dir.names %>%
      lapply(function(x) lapply(ont.list, lapply, `[[`, x)) %>%
      setNames(dir.names)

    #Remove empty entries
    ont.list %<>% lapply(lapply, plyr::compact) %>%
      lapply(lapply, function(x) if(length(x)) x) %>%
      lapply(plyr::compact)

    ont.list %<>%
      names() %>%
      lapply(function(dir) {
        lapply(ont.list[[dir]] %>% names(), function(group) {
          lapply(ont.list[[dir]][[group]] %>% names, function(go) {
            dplyr::mutate(ont.list[[dir]][[group]][[go]], Type = go)
          }) %>% setNames(ont.list[[dir]][[group]] %>% names()) %>%
            dplyr::bind_rows()
        }) %>% setNames(ont.list[[dir]] %>% names())
      }) %>% setNames(dir.names)

    ont.list %<>% filterOntologies(p.adj = p.adj)

    res <- ont.list %>% ontologyListToDf()

    # Filter by min. number of genes per pathway
    if(min.genes > 0) {
      res %<>% lapply(function(g) {
        if(class(g) != "character") {
          idx <- g$GeneRatio %>%
            strsplit("/", fixed=T) %>%
            sapply(`[[`, 1)

          return(g[idx > min.genes,])
        }
        g
      })
    }

    res <- list(df=res,
                go.environment=go.environment)
  } else {
    stop("'type' must be either 'GO' or 'DO'.")
  }

  return(res)
}
