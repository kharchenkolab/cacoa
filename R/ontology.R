#' @import dplyr
NULL

#' @title Get DE ENTREZ IDs
#' @description  Filter and prepare DE genes for ontology calculations
#' @param de.raw List with differentially expressed genes per cell group
#' @param org.db Organism database, e.g., org.Hs.eg.db for human or org.Ms.eg.db for mouse. Input must be of class 'OrgDb'
#' @param p.adj.cutoff Adj. P cutoff for filtering DE genes (default=0.05)
#' @return A list containing DE ENSEMBL gene IDs, and filtered DE genes
getDEEntrezIdsSplitted <- function(de.raw, org.db, p.adj=1) {
  checkPackageInstalled("clusterProfiler", bioc=TRUE)
  if(!(class(org.db) %in% "OrgDb")) stop("'org.db' must be of class 'OrgDb'. Please input an organism database, e.g., org.Hs.eg.db for human data.")
  if(p.adj < 1) warning("You are filtering based on adj. P-value through the 'p.adj' parameter. We do not recommend this. Proceed with caution.")

  de.genes.filtered <- lapply(de.raw, function(df) {
    df %$% Gene[!is.na(padj) & (padj <= p.adj)]
  })

  # Split genes by direction
  gene.ids.tmp <- de.genes.filtered %>% .[sapply(., length) > 0] %>% names() %>% sn() %>% lapply(function(x) {
    de <- de.raw[[x]] %>% .[rownames(.) %in% de.genes.filtered[[x]],]
    list(down=dplyr::filter(de, Z < 0), up=dplyr::filter(de, Z > 0), all=de, universe=de.raw[[x]]) %>%  # Universe contains all input genes
      lapply(function(df) if (nrow(df)) {df %$% Gene[order(abs(Z), decreasing=TRUE)]} else NULL)
  }) %>% lapply(plyr::compact)

  all.genes <- unlist(gene.ids.tmp) %>% unique()
  suppressWarnings(suppressMessages(tryCatch({
    gene.id.map <- clusterProfiler::bitr(all.genes, 'SYMBOL', 'ENTREZID', org.db) %$%
      setNames(ENTREZID, SYMBOL)
  }, error=function(e) {
    stop("Can't find ENTREZIDs for the specified genes. Did you pass the right org.db?")
  })))

  de.gene.ids <- lapply(gene.ids.tmp, lapply, function(gs) as.character(na.omit(gene.id.map[gs])))

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
  names(go.environment) %>% sccore:::sn() %>%
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
  names(go.environment) %>% sccore:::sn() %>%
    lapply(function(ont) enrichGSEGOOpt(go.environment=go.environment[[ont]], ont=ont, ...))
}

groupOntologiesByCluster <- function(ont.clust.df, field="ClusterName") {
  order.anno <- ont.clust.df$Group %>% unique %>% .[order(.)]

  df <- ont.clust.df %>%
    group_by(Group, CN=.[[field]]) %>%
    summarise(p.adjust=min(p.adjust)) %>%
    ungroup() %>%
    mutate(p.adjust=log10(p.adjust)) %>%
    tidyr::spread(Group, p.adjust) %>%
    as.data.frame() %>%
    set_rownames(.$CN) %>%
    .[, 2:ncol(.), drop=FALSE] %>%
    .[, order.anno[order.anno %in% colnames(.)], drop=FALSE]

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
#' @param verbose Print progress (default=T)
#' @param qvalue.cutoff Q value cutoff, please see clusterProfiler package for more information (default=0.2)
#' @param ... Additional parameters for DO/GO/GSEA functions
#' @return A list containing a list of ontologies per type of ontology, and a data frame with merged results
#' @export
estimateOntologyFromIds <- function(de.gene.ids, go.environment, type="GO", org.db=NULL, n.top.genes=5e2, keep.gene.sets=FALSE, verbose=TRUE,
                                    qvalue.cutoff=0.2, ...) {
  checkPackageInstalled("DOSE", bioc=TRUE)

  if(type %in% c("GO","DO") && (n.top.genes > 0)) {
    # Adjust to n.top.genes
    de.gene.ids %<>% lapply(function(celltype) {
        lapply(celltype[-length(celltype)], function(dir) dir[1:min(n.top.genes, length(dir))]) %>%
          append(celltype[length(celltype)])
    })
  }

  # Estimate ontologies
  if(type=="DO") {
    # TODO enable mapping to human genes for non-human data https://support.bioconductor.org/p/88192/
    ont.list <- names(de.gene.ids) %>% sn() %>% plapply(function(id) suppressWarnings(
      lapply(de.gene.ids[[id]][-length(de.gene.ids[[id]])], DOSE::enrichDO,
             universe=de.gene.ids[[id]][["universe"]], qvalueCutoff=qvalue.cutoff, ...)
      ), n.cores=1, progress=verbose)
  } else if(type=="GO") {
    if(verbose) cat("Estimating enriched ontologies ... \n")
    ont.list <- names(de.gene.ids) %>% sn() %>% plapply(function(id) suppressWarnings(
      estimateEnrichedGO(de.gene.ids[[id]][-length(de.gene.ids[[id]])], go.environment = go.environment,
                         universe=de.gene.ids[[id]][["universe"]], org.db=org.db, qvalueCutoff=qvalue.cutoff, ...)
      ), progress = verbose, n.cores = 1)
  } else if(type == "GSEA") {
    if(verbose) cat("Estimating enriched ontologies ... \n")
    ont.list <- names(de.gene.ids) %>% sn() %>% plapply(function(id) {suppressWarnings(suppressMessages(
      estimateEnrichedGSEGO(gene.ids=twist(de.gene.ids[[id]]$universe), org.db=org.db, # TODO: it is broken now because de.gene.ids are not named
                            go.environment=go.environment, organism=clusterProfiler:::get_organism(org.db), ...)
    ))}, progress=verbose, n.cores=1)
  } else {
    stop("'type' must be 'GO', 'DO', or 'GSEA'.")
  }

  if(type != 'DO' && !keep.gene.sets) {
    ont.list %<>% lapply(lapply, lapply, function(x) {x@geneSets <- list(); x})
  }

  return(ont.list)
}

#' @title Prepare plot data
#' @description Prepare ontology results for plotting
#' @param ont.res List of results from estimateOntology
#' @param type Type of ontology results, i.e., GO, GSEA, or DO
#' @param p.adj Cut-off for adj. P values
#' @param min.genes Min. number of significant genes per pathway
#' @return List of ontology data for plotting
preparePlotData <- function(ont.res, type, p.adj, min.genes) {
  dir.names <- c("down", "up", "all")

  if(type == "DO") {
    # Split into fractions
    ont.res <- dir.names %>%
      lapply(function(x) lapply(ont.res, `[[`, x)) %>%
      setNames(dir.names) %>%
      lapply(plyr::compact) # Remove empty entries

    # Extract results
    ont.res %<>% lapply(lapply, function(x) if(length(x)) x@result else x)

    # Prep for filter
    ont.res %<>%
      names() %>%
      lapply(function(dir) {
        lapply(ont.res[[dir]] %>% names(), function(ct) {
          dplyr::mutate(ont.res[[dir]][[ct]], Type = "DO")
        }) %>% setNames(ont.res[[dir]] %>% names())
      }) %>% setNames(ont.res %>% names())

    # Filter by p.adj
    ont.res %<>% filterOntologies(p.adj = p.adj) %>% ontologyListToDf()

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
  } else if(type == "GO") {
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
    ont.res %<>% filterOntologies(p.adj = p.adj) %>% ontologyListToDf()

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

#' @title Swap character strings for their numeric names
#' @description Set numeric names as output and set names as input character strings
#' @param x Character string with numeric names
#' @return Numeric vector with character names
twist <- function(x) {
  x %>% names() %>% as.numeric() %>% setNames(x)
}

addGseaGroup <- function(ont.res) {
  ont.res %>%
    names() %>%
    lapply(function(ct) {
      ont.res[[ct]] %>% mutate(Group = ct)
    }) %>%
    setNames(ont.res %>% names()) %>%
    bind_rows()
}

#' @title Wrap strings for readibility on plots
#' @description Handy function for wrapping long text strings for usage in plots
#' @param strings Text strings
#' @param width Max length before inserting line shift
#' @return Text strings with inserted line shifts
# Source: http://stackoverflow.com/a/7367534/496488
wrap_strings <- function(strings, width) {
  as.character(sapply(strings, function(s) {
    paste(strwrap(s, width = width), collapse="\n")
  }))
}

#' @title Identify ontology families
#' @description Identify parents and/or children of significant ontology terms
#' @param ids Data frame of cell type-specific ontology results from estimateOntology
#' @return List of families and ontology data
identifyFamilies <- function(ids) {
  tmp.parent <- GOfuncR::get_parent_nodes(ids$ID) %>%
    rename(parent_distance = distance) %>%
    filter(parent_distance > 0) %>%
    split(., .$child_go_id) %>%
    lapply(as.list) %>%
    lapply(function(x) x[-1])

  tmp.children <- GOfuncR::get_child_nodes(ids$ID) %>%
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
                  children_enrichment = 0))
    } else if(id %in% names(tmp.children)) {
      append(list(parent_go_id = NULL,
                  parent_name = NULL,
                  parent_distance = NULL,
                  parents_in_IDs = NULL,
                  parent_enrichment = 0),
             tmp.children[[id]])
    } else {
      id
    }
  }) %>%
    setNames(ids$ID)

  # Add description, significance, and type.
  tmp %<>%
    names() %>%
    lapply(function(id) {
      tmp[[id]] %>% append(list(Description = ids$Description[ids$ID == id],
                                Significance = ids$p.adjust[ids$ID == id],
                                Type = ids$Type[ids$ID == id]))
    }) %>%
    setNames(names(tmp))

  # Sort for lonely children (terms with no children)
  idx <- tmp %>%
    sapply(`[[`, "children_enrichment") %>%
    unlist() %>%
    .[. == 0]
  tmp %<>%
    .[idx %>% names()]

  # Rank by enrichment
  idx <- tmp %>%
    sapply(function(x) x$Significance) %>%
    unlist() %>%
    .[order(., decreasing=F)]
  tmp %<>%
    .[names(idx)]

  return(tmp)
}

#' @title Collapse ontology families
#' @description Collapse ontology families based on overlaps between parents and/or children of significant ontology terms
#' @param ont.res Results from identifyFamilies
#' @return List of collapsed families and ontology data
collapseFamilies <- function(ont.res) {
  if(length(ont.res) > 0) {
    # Identify overlapping parents (seeds) between families
    olaps <- sapply(ont.res, `[[`, "parents_in_IDs") %>%
      unlist() %>%
      table() %>%
      .[order(., decreasing = T)] %>%
      .[. > 1]

    if(length(olaps) > 1) {
      # Create logical matrix and list of seeds and families
      olap.matrix <- sapply(ont.res, `[[`, "parents_in_IDs") %>%
        sapply(function(x) names(olaps) %in% x)
      olap.list <- lapply(1:length(olaps), function(id) {
        olap.matrix[,olap.matrix[id,] == T] %>% `rownames<-`(names(olaps))
      }) %>%
        setNames(names(olaps))

      tmp.res <- lapply(1:length(olap.list), function(x) {
        # Investigate overlapping seeds
        tmp.matrix <- olap.list[[x]] %>%
          .[!rownames(.) == names(olaps)[x],]

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
          to_merge <- rownames(tmp.matrix)[tmp] %>%
            .[!. %in% names(olaps)[1:x]]

          if(length(to_merge) > 0) {
            pre.res <- c(olap.list[[x]] %>% colnames(), sapply(to_merge, function(x) colnames(olap.list[x][[1]])) %>% unlist()) %>%
              unique()
          } else {
            pre.res <- olap.list[[x]] %>%
              colnames()
          }
        } else {
          pre.res <- olap.list[[x]] %>%
            colnames()
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
      res <- sapply(ont.res, `[[`, "parents_in_IDs") %>%
        sapply(function(x) names(olaps) %in% x) %>%
        .[.] %>%
        names()
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
      res %<>% .[order(sapply(., length), decreasing=T)] %>%
        append(as.list(names(ont.res)[order(sapply(names(ont.res), function(p) ont.res[[p]]$Significance), decreasing = F)])) %>%
        {setNames(., sapply(1:length(.), function(n) paste0("Family",n)))}
    }

    return(list(families=res,
                data=ont.res))
  } else {
    NULL;
  }
}

#' @title Estimate ontology families
#' @description Estimate ontology families based on ontology results
#' @param ont.list List of results from estimateOntology
#' @param type Type of ontology result, i.e., GO, GSEA, or DO
#' @return List of families and ontology data per cell type
estimateOntologyFamilies <- function(ont.list, type) {
  if(type == "GO") {
    ont.fam <- lapply(ont.list, lapply, lapply, identifyFamilies) %>%
      setNames(names(ont.list))
    lapply(ont.fam, lapply, lapply, collapseFamilies) %>%
      setNames(names(ont.fam))
  } else {
    ont.fam <- lapply(ont.list, lapply, identifyFamilies) %>%
      setNames(names(ont.list))
    lapply(ont.fam, lapply, collapseFamilies) %>%
      setNames(names(ont.fam))
  }
}

### Clustering

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

clusterIndividualGOs <- function(gos, cut.h) {
  clusts.per.go <- lapply(gos, distanceBetweenTerms) %>% lapply(function(ld)
    if (ncol(as.matrix(ld)) == 1) 1 else {hclust(ld) %>% cutree(h=cut.h)})

  clust.df <- names(gos) %>%
    lapply(function(n) mutate(gos[[n]], Cluster=clusts.per.go[[n]])) %>%
    Reduce(rbind, .) %>%
    mutate(Cluster=factor(Cluster, levels=c(0, unique(Cluster)))) %>%
    select(Group, Cluster, Description) %>%
    tidyr::spread(Group, Cluster) %>% as.data.frame() %>%
    set_rownames(.$Description) %>% .[, 2:ncol(.)]

  return(clust.df)
}

getOntClustField <- function(type, subtype, genes) {
  return(paste(type, genes, paste(subtype, collapse="."), "clusters", sep="."))
}

getHeatmapData <- function(ont.sum, fams, type, subtype, genes) {
  fams %<>%
    lapply(function(x) {
      x[[subtype]][[genes]]$families %>% unlist() %>% unique()
    }) %>%
    setNames(names(fams))

  ont.sum %<>% split(., .$Group)

  ont.sum %<>%
    names() %>%
    lapply(function(x) {
      ont.sum[[x]] %>% filter(ID %in% fams[[x]])
    }) %>%
    bind_rows() %>%
    groupOntologiesByCluster(field="Description")
}
