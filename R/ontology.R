#' @import dplyr igraph Rgraphviz graph GOfuncR
NULL

#' @title Prepare ontology data
#' @description  Filter and prepare DE genes for ontology calculations
#' @param cms list of counts matrices; column for gene and row for cell
#' @param de.raw List with differentially expressed genes per cell group
#' @param cell.groups Vector indicating cell groups with cell names (default: stored vector)
#' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
#' @param p.adj.cutoff Adj. P cutoff for filtering DE genes (default=0.05)
#' @param expr.cutoff Cutoff for cells per group expressing a DE gene, i.e., cutoff for highly-expressed genes (default=0.05)
#' @param universe Background gene list, NULL will take all input genes in de.raw (default=NULL)
#' @param transposed Indicate whether raw counts are transposed, i.e., cells as rows, genes as columns (default=T)
#' @param verbose Print progress (default=T)
#' @param n.cores Number of cores to use (default=1)
#' @return A list containing DE ENSEMBL gene IDs, and filtered DE genes
#' @export
prepareOntologyData <- function(cms, de.raw, cell.groups, org.db, p.adj = 1, expr.cutoff = 0.05, universe = NULL, transposed = T,  verbose = T, n.cores = 1) {
  # Checks
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("You need 'clusterProfiler' to perform ontology analysis.")
  if(class(org.db) != "OrgDb") stop("'org.db' must be of class 'OrgDb'. Please input an organism database, e.g., org.Hs.eg.db for human data.")

  # Filter by p.adj and remove NA p-values
  de.filtered <- lapply(de.raw, function(df) {
    df[!is.na(df$padj) & (df$padj <= p.adj),]
  })

  # Filter by expr.cutoff
  if(expr.cutoff > 0) {
    if(verbose) cat("Merging count matrices ... ")
    cm_merged <- sccore:::mergeCountMatrices(cms, transposed = transposed, n.cores = n.cores)

    if(!transposed) {
      if(verbose) cat("done!\nTransposing merged count matrix ... ")
      cm_merged %<>% Matrix::t()
    }

    if(verbose) cat("done!\nFiltering DE genes ... ")
    cm_bool <- (cm_merged > 1) * 1

    cm_collapsed_bool <- conos:::collapseCellsByType(cm_bool, cell.groups %>%
                                                       .[. %in% names(de.raw)] %>%
                                                       factor(), min.cell.count=0)

    de.genes.filtered <- mapply(intersect,
                                lapply(de.filtered, rownames),
                                ((cm_collapsed_bool > as.vector(table(cell.groups %>% .[. %in% names(de.raw)] %>%
                                                                        factor)[rownames(cm_collapsed_bool)] * expr.cutoff)) %>%
                                   apply(1, function(row) names(which(row))))[names(de.filtered)])

    if(verbose) cat("done!\n")
  } else {
    de.genes.filtered <- lapply(de.filtered, rownames)
  }

  if(verbose) cat("Retrieving Entrez Gene IDs ... ")
  groups <- de.genes.filtered %>% .[sapply(., length) > 0] %>% names()

  gene.ids.tmp <- groups %>% lapply(function(x) {
    de <- de.raw[[x]] %>% .[rownames(.) %in% de.genes.filtered[[x]],]

    list(down = de %>% dplyr::filter(Z < 0) %>% .[order(.$Z, decreasing = F),] %>% rownames(),
         up = de %>% dplyr::filter(Z > 0) %>% .[order(.$Z, decreasing = T),] %>% rownames(),
         all = de %>% .[order((.$Z %>% abs()), decreasing = T),] %>% rownames(),
         universe = de.raw[[x]] %>% .[order(.$Z, decreasing = T),] %>% {setNames(rownames(.), .$Z)}) # Universe contains all input genes
    }) %>% setNames(groups) %>%
    lapply(lapply, function(y) if(length(y)) y) %>%
    lapply(plyr::compact)

  de.gene.ids <- lapply(gene.ids.tmp, lapply, function(y) suppressMessages(suppressWarnings({tryCatch({clusterProfiler::bitr(y, 'SYMBOL', 'ENTREZID', org.db)}, error = function(err) NULL)}))) %>%
    lapply(plyr::compact) %>%
    lapply(lapply, as.list)

  # Transfer Z-values for GSEA
  de.gene.ids %<>%
    names() %>%
    lapply(function(id) {
      tmp <- de.gene.ids[[id]]
      tmp$universe$ENTREZID %<>% setNames(gene.ids.tmp[[id]]$universe %>% {names(.)[match(de.gene.ids[[id]]$universe$SYMBOL, .)]})
      return(tmp)
    }) %>%
    setNames(de.gene.ids %>% names()) %>%
    lapply(lapply, `[[`, "ENTREZID")

  if(verbose) cat("done!\nAll done!\n")

  return(de.gene.ids)
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
    lapply(function(ont) lapply(de.gene.ids, enrichGOOpt, go.environment=go.environment[[ont]], ont=ont, ...))
}

enrichGSEGOOpt <- function(gene.ids, org.db, organism, keyType = "ENTREZID", go.environment, ont = "BP", pvalueCutoff = 0.05,
                           pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500, readable = FALSE, eps = 0, exponent = 1, seed = F, verbose = F) {
  ont %<>% toupper %>% match.arg(c("BP", "CC", "MF"))

  res <- DOSE:::GSEA_internal(gene.ids, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize,
                              maxGSSize = maxGSSize, USER_DATA = go.environment, eps = eps, exponent = exponent, seed = seed,
                              verbose = verbose)
  if(is.null(res))
    return(res)

  res@keytype <- keyType
  res@organism <- organism
  if(readable) {
    res <- DOSE::setReadable(res, org.db)
  }
  res@setType <- ont
  if(!is.null(res)) return(res)
}

estimateEnrichedGSEGO <- function(go.environment, ...) {
  names(go.environment) %>%
    sccore:::sn() %>%
    lapply(function(ont) enrichGSEGOOpt(go.environment=go.environment[[ont]], ont=ont, ...))
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
#' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
#' @param n.top.genes Number of most different genes to take as input. If less are left after filtering for p.adj.cutoff, additional genes are included. To disable, set n.top.genes=0 (default=1e2)
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
estimateOntology <- function(type = "GO", org.db=NULL, de.gene.ids, n.top.genes = 5e2, go.environment = NULL, keep.geneSets = F, p.adj=0.05, p.adjust.method="BH", readable=T, verbose=T, min.genes = 0, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 5e2) {
  # Preparations
  dir.names <- c("down", "up", "all")
  if (!requireNamespace("DOSE", quietly = TRUE)) stop("You need 'DOSE' to perform ontology analysis.")

  if(type %in% c("GO","DO")) {
    # Adjust to n.top.genes
    if(n.top.genes > 0) {
      de.gene.ids %<>%
        lapply(function(celltype) {
          lapply(celltype[-length(celltype)], function(dir) {
            if(length(dir) < n.top.genes) dir else dir[1:n.top.genes]
          }) %>% append(celltype[length(celltype)])
        })
    }
  }

  if(type %in% c("GO","GSEA")) {
    # Extract GO environment if necessary
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
  }

  # Estimate ontologies
  if(type=="DO") {
    # TODO enable mapping to human genes for non-human data https://support.bioconductor.org/p/88192/
    res <- sccore:::plapply(names(de.gene.ids), function(id) suppressWarnings(lapply(de.gene.ids[[id]][-length(de.gene.ids[[id]])],
                                                                                          DOSE::enrichDO,
                                                                                          pAdjustMethod=p.adjust.method,
                                                                                          universe=de.gene.ids[[id]][["universe"]],
                                                                                          readable=readable,
                                                                                          qvalueCutoff = qvalueCutoff,
                                                                                          minGSSize = minGSSize,
                                                                                          maxGSSize = maxGSSize)),
                                 n.cores=1, progress=verbose) %>%
      setNames(names(de.gene.ids))
  } else if(type=="GO") {
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
                                 progress = verbose, n.cores = 1) %>%
      setNames(names(de.gene.ids))

    if(keep.geneSets == F) {
      ont.list %<>% lapply(lapply, lapply, function(x) {
        x@geneSets <- list()
        return(x)
      })
    }

    res <- list(res=ont.list,
                go.environment=go.environment)
  } else if(type == "GSEA") {
    if(verbose) cat("Estimating enriched ontologies ... \n")
    ont.list <- sccore:::plapply(names(de.gene.ids), function(id) {
      suppressWarnings(suppressMessages(estimateEnrichedGSEGO(gene.ids = de.gene.ids[[id]]$universe %>% twist(),
                                                              org.db = org.db,
                                                              go.environment = go.environment,
                                                              readable = readable,
                                                              pAdjustMethod = p.adjust.method,
                                                              minGSSize = minGSSize,
                                                              maxGSSize = maxGSSize,
                                                              organism = clusterProfiler:::get_organism(org.db))))
    }, progress = verbose, n.cores = 1) %>%
      setNames(names(de.gene.ids))

    if(keep.geneSets == F) {
      ont.list %<>% lapply(lapply, function(x) {
        x@geneSets <- list()
        return(x)
      })
    }

    res <- list(res=ont.list,
                go.environment=go.environment)
  } else {
    stop("'type' must be 'GO', 'DO', or 'GSEA'.")
  }
  return(res)
}

preparePlotData <- function(ont.res, type, p.adj, min.genes) {
  dir.names <- c("down", "up", "all")

  if(type == "DO") { # Functionality needs check
    # Split into fractions
    ont.res <- dir.names %>%
      lapply(function(x) lapply(ont.res, `[[`, x)) %>%
      setNames(dir.names) %>%
      lapply(plyr::compact) # Remove empty entries

    # Extract results
    ont.res %<>% lapply(lapply, function(x) if(length(x)) x@result else x)

    ### TODO: prep for filter ###

    # Filter by p.adj
    ont.res %<>% cacoa:::filterOntologies(p.adj = p.adj) %>% cacoa:::ontologyListToDf()

    # Filter by min. number of genes per pathway
    if(min.genes > 1) {
      ont.res %<>% lapply(function(g) {
        if(class(g) != "character") {
          idx <- g$GeneRatio %>%
            strsplit("/", fixed=T) %>%
            sapply(`[[`, 1)

          return(g[idx > min.genes,])
        }
        g
      })
    }
  } else if(type %in% c("GO","BP","MF","CC")) {
    # Split into different fractions
    ont.res <- dir.names %>%
      lapply(function(x) lapply(ont.res, lapply, `[[`, x)) %>%
      setNames(dir.names)

    # Extract results
    ont.res %<>% lapply(lapply, lapply, function(x) if(length(x)) x@result else x)

    # Remove empty entries
    ont.res %<>% lapply(lapply, plyr::compact) %>%
      lapply(lapply, function(x) if(length(x)) x) %>%
      lapply(plyr::compact)

    # Prep for filter
    ont.res %<>%
      names() %>%
      lapply(function(dir) {
        lapply(ont.res[[dir]] %>% names(), function(group) {
          lapply(ont.res[[dir]][[group]] %>% names, function(go) {
            dplyr::mutate(ont.res[[dir]][[group]][[go]], Type = go)
          }) %>% setNames(ont.res[[dir]][[group]] %>% names()) %>%
            dplyr::bind_rows()
        }) %>% setNames(ont.res[[dir]] %>% names())
      }) %>% setNames(dir.names)

    # Filter by p.adj
    ont.res %<>% cacoa:::filterOntologies(p.adj = p.adj) %>% cacoa:::ontologyListToDf()

    # Filter by min. number of genes per pathway
    if(min.genes > 1) {
      ont.res %<>% lapply(function(g) {
        if(class(g) != "character") {
          idx <- g$GeneRatio %>%
            strsplit("/", fixed=T) %>%
            sapply(`[[`, 1)

          return(g[idx > min.genes,])
        }
        g
      })
    }
  } else if(type == "GSEA") {
    # Extract results
    ont.res %<>% lapply(lapply, function(x) if(length(x)) x@result else x) %>%
      lapply(plyr::compact) # Remove empty entries

    # Prep for filter
    ont.res %<>%
      names() %>%
      lapply(function(ct) {
        lapply(ont.res[[ct]] %>% names(), function(ont) {
          dplyr::mutate(ont.res[[ct]][[ont]], Type = ont)
          }) %>% setNames(ont.res[[ct]] %>% names()) %>%
            dplyr::bind_rows()
        }) %>% setNames(ont.res %>% names())

    # Filter by min. number of genes per pathway
    if(min.genes > 1) {
      ont.res %<>% lapply(function(g) {
        if(class(g) != "character") {
          idx <- g$core_enrichment %>%
            strsplit("/", fixed=T) %>%
            sapply(length)

          return(g[idx > min.genes,])
        }
        g
      })
    }
  }
  return(ont.res)
}

twist <- function(x) {
  x %>% names() %>% as.numeric() %>% setNames(x)
}

addGroup <- function(ont.res) {
  ont.res %>%
    names() %>%
    lapply(function(ct) {
      ont.res[[ct]] %>% mutate(Group = ct)
    }) %>%
    setNames(ont.res %>% names()) %>%
    bind_rows()
}

# Source: http://stackoverflow.com/a/7367534/496488
wrap_strings <- function(strings, width) {
  as.character(sapply(strings, function(s) {
    paste(strwrap(s, width = width), collapse="\n")
  }))
}

estimateOntologyFamilies <- function(ont.list, p.adj = 0.05) {
  ont.fam <- lapply(ont.list, lapply, function(ids) {
    # ids <- tmp.ids %>% dplyr::filter(p.adjust <= p.adj)
    tmp.parent <- get_parent_nodes(ids$ID) %>%
      rename(parent_distance = distance) %>%
      filter(parent_distance > 0) %>%
      split(., .$child_go_id) %>%
      lapply(as.list) %>%
      lapply(function(x) x[-1])

    tmp.children <- get_child_nodes(ids$ID) %>%
      rename(child_distance = distance) %>%
      filter(child_distance > 0) %>%
      split(., .$parent_go_id) %>%
      lapply(as.list) %>%
      lapply(function(x) x[-1])

    # Find significant parents
    tmp.parent %<>% lapply(function(x) {
      x$parents_in_IDs <- x$parent_go_id %>% .[. %in% ids$ID] %>% setNames(., ids$p.adjust[match(., ids$ID)])
      x$parent_enrichment <- length(x$parents_in_IDs)/length(x$parent_go_id)
      return(x)
    })

    # Find significant children
    tmp.children %<>% lapply(function(x) {
      x$children_in_IDs <- x$child_go_id %>% .[. %in% ids$ID]
      x$children_enrichment <- length(x$children_in_IDs)/length(x$child_go_id)
      return(x)
    })

    tmp <- ids$ID %>% lapply(function(id) {
      if((id %in% names(tmp.parent)) & (id %in% names(tmp.children))) {
        append(tmp.parent[[id]], tmp.children[[id]])
      } else if(id %in% names(tmp.parent)) {
        append(tmp.parent[[id]],
               list(child_go_id = NULL,
                    child_name = NULL,
                    child_distance = NULL,
                    children_in_IDs = NULL,
                    children_enrichment = NULL))
      } else if(id %in% names(tmp.children)) {
        append(list(parent_go_id = NULL,
                    parent_name = NULL,
                    parent_distance = NULL,
                    parents_in_IDs = NULL,
                    parent_enrichment = NULL),
               tmp.children[[id]])
      } else {
        id
      }
    }) %>%
      setNames(ids$ID)

    # # MAY DELETE
    # tmp %<>% lapply(function(id) {
    #   suppressWarnings(if(length(id) == 10) {
    #     id$combined_enrichment <- (length(id$parents_in_IDs) + length(id$children_in_IDs))/(length(id$parent_go_id) + length(id$child_go_id))
    #   } else {
    #     id$combined_enrichment <- NULL
    #   })
    #   return(id)
    # })

    # Add description, significance, and type.
    # TODO: Is it faster by extracting first (= 1 search) and then select?
    tmp %<>% names() %>% lapply(function(id) {
      tmp[[id]] %>% append(list(Description = ids$Description[ids$ID == id],
                                Significance = ids$p.adjust[ids$ID == id],
                                Type = ids$Type[ids$ID == id]))
    }) %>% setNames(names(tmp))

    # Sort for lonely children (terms with no children) UPDATE !!!
    idx <- tmp %>% sapply(`[[`, "children_enrichment") %>% unlist() %>% .[. == 0] # TODO: Check, does it remove too much?
    tmp %<>% .[idx %>% names()]

    # Rank by enrichment ### UPDATE LOWEST P VALUE !!!
    idx <- tmp %>% sapply(function(x) x$parent_enrichment) %>% unlist() %>% .[order(., decreasing=T)]
    tmp %<>% .[names(idx)]

    return(tmp)
  }) %>% setNames(names(ont.list))

  # Collapse families
  lapply(ont.fam, lapply, function(ont.res) {
    if(length(ont.res) > 0) {
      # Identify overlapping parents (seeds) between families
      olaps <- sapply(ont.res, `[[`, "parents_in_IDs") %>% unlist() %>% table() %>% .[order(., decreasing = T)] %>% .[. > 1]

      if(length(olaps) > 1) {
        # Create logical matrix and list of seeds and families
        olap.matrix <- sapply(ont.res, `[[`, "parents_in_IDs") %>% sapply(function(x) names(olaps) %in% x)
        olap.list <- lapply(1:length(olaps), function(id) {
          olap.matrix[,olap.matrix[id,] == T] %>% `rownames<-`(names(olaps))
        }) %>% setNames(names(olaps))

        tmp.res <- lapply(1:length(olap.list), function(x) {
          # Investigate overlapping seeds
          tmp.matrix <- olap.list[[x]] %>% .[!rownames(.) == names(olaps)[x],]

          if(is.matrix(tmp.matrix)) {
            tmp <- c()
            for(r in 1:nrow(tmp.matrix)) {
              tmp <- c(tmp, any(tmp.matrix[r,]))
            }
          } else {
            tmp <- any(tmp.matrix)
          }

          if(tmp %>% any) {
            # Select additional seeds to merge
            to_merge <- rownames(tmp.matrix)[tmp] %>% .[!. %in% names(olaps)[1:x]]

            if(length(to_merge) > 0) {
              pre.res <- c(olap.list[[x]] %>% colnames(), sapply(to_merge, function(x) colnames(olap.list[x][[1]])) %>% unlist()) %>% unique()
            } else {
              pre.res <- olap.list[[x]] %>% colnames()
            }
          } else {
            pre.res <- olap.list[[x]] %>% colnames()
          }
          return(pre.res)
        })

        # Filter seeds that have been merged upstream
        idx <- sapply(length(tmp.res):2, function(x) {
          all(tmp.res[[x]] %in% unlist(tmp.res[1:(x-1)]))
        }) %>%
          setNames(length(tmp.res):2) %>%
          .[.] %>%
          names() %>%
          rev() %>% # Shouldn't be necessary, or what?
          as.numeric()

        # Return result with index if there are any merged seeds, or else return seed list
        if(length(idx) > 0) {
          res <- tmp.res[-idx]
        } else {
          res <- tmp.res
        }
      } else if(length(olaps == 1)) {
        res <- sapply(ont.res, `[[`, "parents_in_IDs") %>% sapply(function(x) names(olaps) %in% x) %>% .[.] %>% names()
      } else {
        res <- list()
      }

      # Add remaining families, order by adj. P of lonely child
      if(length(res) == 0) {
        res <- ont.res %>%
          names() %>%
          as.list() %>%
          setNames(sapply(1:length(ont.res), function(n) paste0("Family",n)))
      } else {
        # res <- append(res, as.list(names(ont.res)[!names(ont.res) %in% unlist(res)])) %>%
        res %<>% .[order(sapply(., length), decreasing=T)] %>%
          append(as.list(names(ont.res)[order(sapply(names(ont.res), function(p) ont.res[[p]]$Significance), decreasing = F)])) %>%

          {setNames(., sapply(1:length(.), function(n) paste0("Family",n)))}
      }

      return(list(families=res,
                  data=ont.res))
    } else {
      NULL;
    }
  }) %>% setNames(names(ont.fam))
}
