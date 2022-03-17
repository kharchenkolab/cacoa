#' @import ape
NULL

#' Construct the canonical tree
#'
#' @param cnts Table with cell type counts
#' @param groups Groups variable for samples
#' @return phylo tree
#' @keywords internal
constructTree <- function(cnts, groups, partition.thresh = 0){
  checkDataGroups(cnts, groups)
  # ---------

  n.cells <- ncol(cnts)

  sbp.cda <- matrix(0, nrow = n.cells, ncol = n.cells-1, dimnames = list(colnames(cnts), c()))

  unsolved.cells <- list(rownames(sbp.cda))  # List of cell types to separate
  for(id.bal in 1:(n.cells-1)){
    # If list for separation contains two cell types: perform random separation
    if(length(unsolved.cells[[id.bal]]) == 2){
      sbp.cda[unsolved.cells[[id.bal]][1],id.bal] <- 1
      sbp.cda[unsolved.cells[[id.bal]][2],id.bal] <- -1
      next
    }
    # ------------------------------------------------
    #  Data for the current subset of cells
    d.tmp <- cnts[, colnames(cnts) %in% unsolved.cells[[id.bal]]]

    # -------
    # Define the most contrast balance
    # can.loadings <- getCdaLoadings(d.tmp, d.groups)
    can.loadings <- getLoadings(d.tmp, groups)

    cells.tmp <- rownames(can.loadings)

    # -------
    # Get cell types from opposite sides of the principal balance
    cells.plus <- cells.tmp[can.loadings > partition.thresh]
    cells.minus <- cells.tmp[can.loadings < partition.thresh]

    sbp.cda[cells.minus, id.bal] <- -1
    sbp.cda[cells.plus, id.bal] <- 1

    if(length(cells.minus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.minus
    }

    if(length(cells.plus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.plus
    }

  }

  # res <- list(sbp = sbp.cda, partiotions = unsolved.cells)
  tree <- sbp2tree(sbp.cda)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)

  return(list(tree = tree, sbp = sbp.cda, dendro = d.cur))
}


#' @keywords internal
constructTreeUp <- function(freqs, groups) {

  sbp <- sbpDiff(freqs, groups)

  tree <- sbp2tree(sbp)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)
  return(list(tree = tree, sbp = sbp, dendro = d.cur))
}


#' @keywords internal
tree2dendro_my <- function(tree){
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  return(d.cur)
}

#' @keywords internal
constructTreeUpDown <- function(cnts, groups){
  cnts[cnts <= 0] <- 0.5

  ref.set <- referenceSet(cnts, groups)
  cell.lists <- ref.set$cell.list

  cell.types <- colnames(cnts)


  # Constrtuct sbp from lists
  freqs <- (cnts)/rowSums(cnts)
  freqs.lists <- c()
  for (i in 1:length(cell.lists)) {
    freqs.lists <- cbind(freqs.lists, apply(freqs[,cell.lists[[i]],drop=FALSE], 1, psych::geometric.mean))
  }
  colnames(freqs.lists) <- paste('tmp', 1:length(cell.lists), sep = '')

  # t.list <- constructTree(freqs.lists, groups)
  t.list <- constructBestPartitionTree(freqs.lists, groups)
  sbp.list <- t.list$sbp

  sbp.all <- matrix(ncol = 0, nrow = length(cell.types), dimnames = list(cell.types, c()))
  for(k in 1:ncol(sbp.list)){
    p <- sbp.list[,k]
    i.list.plus <- which(p > 0)
    i.list.minus <- which(p < 0)

    sbp.tmp <- rep(0, length(cell.types))
    for(i.plus in i.list.plus){
      sbp.tmp[cell.types %in% cell.lists[[i.plus]]] = 1
    }
    for(i.plus in i.list.minus){
      sbp.tmp[cell.types %in% cell.lists[[i.plus]]] = -1
    }
    sbp.all <- cbind(sbp.all, sbp.tmp)
  }

  # Construct sbps from sublists
  for(i in 1:length(cell.lists)){
    freqs.tmp <- freqs[,cell.lists[[i]], drop = FALSE]
    if(ncol(freqs.tmp) == 1) next

    # tmp <- constructTree(freqs.tmp, groups)
    tmp <- constructBestPartitionTree(freqs.tmp, groups)
    sbp.tmp <- tmp$sbp

    sbp.add <- matrix(0, nrow = nrow(sbp.all), ncol = ncol(sbp.tmp))
    rownames(sbp.add) <- rownames(sbp.all)
    sbp.add[rownames(sbp.tmp),] <- sbp.tmp
    sbp.all <- cbind(sbp.all, sbp.add)

  }

  tree <- sbp2tree(sbp.all)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)
  return(list(tree = tree, sbp = sbp.all, dendro = d.cur))
}


#' Get tree from balances
#'
#' @param sbpart A contrast matrix for balances, sequential binary partition
#' @return phylo tree
#' @details We use ellipsis because a previous version of this function
#'          contains additional and insignificant parameters.
#' @keywords internal
sbp2tree <- function(sbpart){
  checkSbpWhole(sbpart)
  # ---------

  n.cells <- nrow(sbpart)
  edges <- c(n.cells+1, n.cells+1)
  id <- 2*n.cells - 1

  collapse <- sbpart
  rownames(collapse) <- 1:n.cells

  for(i in rev(1:(n.cells - 1)) ){
    ids <- which(collapse[,i] != 0)

    edges <- rbind(edges, as.integer(c(id, rownames(collapse)[ids[1]])))
    edges <- rbind(edges, as.integer(c(id, rownames(collapse)[ids[2]])))
    rownames(collapse)[ids[2]] <- id
    id <- id - 1

    collapse[ids[1],] <- 0
  }

  # Create face tree, because we need a tidytree object
  x <- tidytree::as_tibble(rtree(n = n.cells, br = NULL))
  x[,1:2] <- edges

  t.tmp <- as.phylo(x)
  t.tmp$tip.label <- rownames(sbpart)
  return(t.tmp)
}


#' Get balances for all inner nodes of the tree
#'
#' @param t A phylo object
#' @param log.f Logarithms of frequnencies
#' @return Balances, contrast matrix and names on cell types
#' @keywords internal
getNodeBalances <- function(log.f, sbp.my){

  # get psi matrix
  psi <- sbp.my * 0
  for(k in 1:ncol(sbp.my))
  {
    p <- sbp.my[, k]  # k-th partition
    n1 <- sum(p > 0)
    n2 <- sum(p < 0)
    p[p > 0] <- 1/n1
    p[p < 0] <- -1/n2
    psi[, k] <- p * sqrt(n1*n2/(n1+n2))
  }

  balances <- log.f %*% psi
  return(balances)
}

#' @keywords internal
bestPartition <- function(freqs.tmp, groups){


  n.cells <- ncol(freqs.tmp)
  sbp <- 1
  for(i in 2:n.cells){
    sbp <- rbind(cbind(sbp, 1),
                cbind(sbp, -1))
  }
  sbp <- sbp[-c(1, nrow(sbp)),]
  log.f <- log(freqs.tmp)

  for(i in 1:nrow(sbp)){
    p <- sbp[i,]
    p[p < 0] <- p[p < 0] / sum(p < 0)
    p[p > 0] <- p[p > 0] / sum(p > 0)
    sbp[i,] <- p
  }
  colnames(sbp) <- colnames(freqs.tmp)


  bals <- sbp %*% t(log.f)
  p <- c()
  for(i in 1:nrow(bals)){
    b <- bals[i,]
    mod <- lm(groups ~ b)
    res <- summary(mod)
    pval <- res$coefficients[2,4]
    p <- c(p, pval)
  }

  sbp.best <- (sbp[p == min(p),] > 0) * 2 - 1
  names(sbp.best) <- colnames(freqs.tmp)
  # print(sbp.best)


  return(sbp.best)
}



#' Construct the canonical tree
#'
#' @param cnts Table with cell type counts
#' @param groups Groups variable for samples
#' @return phylo tree
#' @keywords internal
constructBestPartitionTree <- function(cnts, groups, partition.thresh = 0){
  checkDataGroups(cnts, groups)
  # ---------

  n.cells <- ncol(cnts)

  sbp.cda <- matrix(0, nrow = n.cells, ncol = n.cells-1, dimnames = list(colnames(cnts), c()))

  unsolved.cells <- list(rownames(sbp.cda))  # List of cell types to separate
  for(id.bal in 1:(n.cells-1)){
    # If list for separation contains two cell types: perform random separation
    if(length(unsolved.cells[[id.bal]]) == 2){
      sbp.cda[unsolved.cells[[id.bal]][1],id.bal] <- 1
      sbp.cda[unsolved.cells[[id.bal]][2],id.bal] <- -1
      next
    }
    # ------------------------------------------------
    #  Data for the current subset of cells
    d.tmp <- cnts[, colnames(cnts) %in% unsolved.cells[[id.bal]]]

    # -------
    # Find the best partition
    sbp.best <- bestPartition(d.tmp, groups)

    # -------
    # Get cell types from opposite sides of the principal balance
    cells.tmp <- colnames(d.tmp)

    cells.plus <- cells.tmp[sbp.best > 0]
    cells.minus <- cells.tmp[sbp.best < 0]
    #
    # print(cells.plus)
    # print(cells.minus)

    sbp.cda[cells.minus, id.bal] <- -1
    sbp.cda[cells.plus, id.bal] <- 1

    if(length(cells.minus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.minus
    }

    if(length(cells.plus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.plus
    }

  }

  # res <- list(sbp = sbp.cda, partiotions = unsolved.cells)
  tree <- sbp2tree(sbp.cda)
  h.tmp <- compute.brlen(tree, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(h.tmp)
  tree <- as.phylo(h.tmp)

  return(list(tree = tree, sbp = sbp.cda, dendro = d.cur))
}

#' @keywords internal
sbpInNodes <- function(tree){
  edges <- tree$edge
  n.nodes <- max(edges)
  n.leaves <- length(tree$tip.label)
  edge.lists <- list()
  for(i in 1:n.leaves){
    edge.lists[[i]] <- c(tree$tip.label[i])
  }

  sbp <- matrix(0, nrow = n.nodes, ncol = n.leaves, dimnames = list(c(), tree$tip.label))
  for(i in rev((n.leaves+1):n.nodes)){
    idx <- which(edges[,1] == i)
    cells.plus <- edge.lists[[edges[idx[1],2] ]]
    cells.minus <- edge.lists[[edges[idx[2],2] ]]
    sbp[i, cells.plus] <- 1
    sbp[i, cells.minus] <- -1

    edge.lists[[i]] <- c(cells.plus, cells.minus)
  }
  sbp <- sbp[-c(1:n.leaves),]
  sbp <- t(sbp)

  return(sbp)
}
