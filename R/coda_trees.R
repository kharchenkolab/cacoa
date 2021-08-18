#' @import ape
NULL

#' @description Get inner nodes on the tree from leaves to the root
#'
#' @param t A phylo object
#' @return Order of inner nodes
getNodeOrder <- function(t) {
  checkTree(t)
  # ---------

  # Sort inner nodes by id
  x <- t$edge[order(t$edge[,1]),]
  nx <- nrow(x)
  x1 <- x[seq(1, nx, 2),]  # First child
  x2 <- x[seq(2, nx, 2),]  # Second child

  id.used <- 1:(t$Nnode + 1)  # Used nodes and leaves
  id.nodes <- c()  # order of nodes co collapse
  id.parents <- x1[,1]  # names of inner nodes
  id.children <- cbind(x1[,2], x2[,2])  # children

  while(length(id.nodes) < t$Nnode){
    # which children have already taken?
    ids.taken <- (id.children %in% id.used) & (id.children > 0)
    # which nodes can be taken at the current step?
    node.condition <- rowSums(ids.taken) == 2
    # Remember new nodes
    id.nodes <- c(id.nodes, id.parents[node.condition])

    # Do not use current nodes again
    id.used <- c(id.used, id.nodes)
    id.children[node.condition,] <- id.children[node.condition,]  * -1
  }

  return(id.nodes)
}

#' #' @description Get balances for all inner nodes of the tree
#' #'
#' #' @param t A phylo object
#' #' @param log.freq Logarithms of frequnencies
#' #' @return Balances, contrast matrix and names on cell types
#' #' @details We use ellipsis because a previous version of this function
#' #'          contains additional and insignificant parameters.
#' getNodeBalances <- function(t, log.freq, ...){
#'   checkDataTree(log.freq, t)
#'   # ---------
#' 
#'   sample.names <- rownames(log.freq)
#'   # Initialisation of balances
#'   n.cells <- ncol(log.freq)
#'   n.samples <- nrow(log.freq)
#'   log.freq.ext <- cbind(log.freq[,t$tip.label], matrix(0, nrow=n.samples, ncol=n.cells-1))  # log-frequencoes extended for inner nodes
#'   n.in.clade <- c(rep(1,n.cells), rep(0,n.cells-1)) # number of cell-types in a clade
#' 
#' 
#'   balances <- matrix(0, nrow=n.samples, ncol=n.cells-1)
#'   rownames(balances) <- sample.names
#'   sbp.my <- rbind(diag(x = 1, n.cells, n.cells), matrix(0, n.cells - 1,n.cells))
#' 
#'   edges <- t$edge
#'   id.edges <- getNodeOrder(t)
#'   # Calculate balances in traversing the tree
#'   for(k in id.edges)
#'   {
#'     clades <- edges[edges[,1]==k,2]
#'     id1 <- clades[1]
#'     id2 <- clades[2]
#' 
#'     if (n.in.clade[id1] + n.in.clade[id2] < 2)
#'       break
#' 
#'     n1 <- n.in.clade[id1]
#'     n2 <- n.in.clade[id2]
#'     balances[,k-n.cells] <- (log.freq.ext[,id1] / n1 - log.freq.ext[,id2] / n2) * sqrt(n1*n2/(n1+n2))
#'     log.freq.ext[,k] <- log.freq.ext[,id1] + log.freq.ext[,id2]
#'     n.in.clade[k] <- n.in.clade[id1] + n.in.clade[id2]
#' 
#'     sbp.my[k,] <- abs(sbp.my[id1,]) - abs(sbp.my[id2,])
#'   }
#' 
#'   for(k in id.edges)
#'   {
#'     n1 <- sum(sbp.my[k,]>0)
#'     n2 <- sum(sbp.my[k,]<0)
#'     sbp.my[k,sbp.my[k,]>0] <- 1/n1
#'     sbp.my[k,sbp.my[k,]<0] <- -1/n2
#'     sbp.my[k,] <- sbp.my[k,] * sqrt(n1*n2/(n1+n2))
#'   }
#' 
#'   res <- list(
#'     bal = balances,
#'     sbp = sbp.my[(n.cells+1):(2*n.cells-1),],
#'     cells = t$tip.label)
#'   return(res)
#' }



#' @description Construct the canonical tree
#'
#' @param d.counts Table with cell type counts
#' @param d.groups Groups variable for samples
#' @return phylo tree
constructCanonicalTree <- function(d.counts, d.groups){
  checkDataGroups(d.counts, d.groups)
  # ---------

  n.cells <- ncol(d.counts)

  sbp.cda <- matrix(0, nrow = n.cells, ncol = n.cells-1, dimnames = list(colnames(d.counts), c()))

  unsolved.cells <- list(rownames(sbp.cda))  # List of cell types to separate
  for(id.bal in 1:(n.cells-1)){
    # If list to separate contains two cell types: perform random separation
    if(length(unsolved.cells[[id.bal]]) == 2){
      sbp.cda[unsolved.cells[[id.bal]][1],id.bal] <- 1
      sbp.cda[unsolved.cells[[id.bal]][2],id.bal] <- -1
      next
    }
    # ------------------------------------------------
    #  Data for the current subset of cells
    d.tmp <- d.counts[, colnames(d.counts) %in% unsolved.cells[[id.bal]]]

    # -------
    # Difene the most contrast balance
    # can.loadings <- getCdaLoadings(d.tmp, d.groups)
    can.loadings <- getLoadings(d.tmp, d.groups)

    cells.tmp <- rownames(can.loadings)

    # ga.res <- gaPartitionAll(d.tmp, d.groups, n.epoch = 300)
    # can.loadings <- ga.res[1,]
    # names(can.loadings) <- colnames(d.tmp)
    # print(can.loadings)
    # cells.tmp <- names(can.loadings)
    #

    # cda <- composition.sampling(d.tmp, d.groups, n.iter = 100)
    # can.loadings <- rowMeans(cda$balances)
    # cells.tmp <- names(can.loadings)
    # print(cells.tmp)

    # -------
    # Get cell types from opposite sides of the principal balance
    cells.plus <- cells.tmp[can.loadings > 0]
    cells.minus <- cells.tmp[can.loadings < 0]
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
  return(tree)
}
