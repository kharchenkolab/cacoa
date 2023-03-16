#' @import ape
NULL

#' Get loadings by different methods
#'
#' @param cnts Counts of cell typer in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param method Method to get loadings
#' @param ref.cell.type Reference cell type
#' @return Updated data frame with Z scores
#' @keywords internal
getLoadings <- function(cnts, groups, method = c('lda', 'svm', 'cda', 'cda.std'), ref.cell.type = NULL) {
  # Checks
  method <- match.arg(method)
  if (method == "svm") checkPackageInstalled("e1071", cran=TRUE)
  if (method == "cda") checkPackageInstalled("candisc", cran=TRUE)
  if(!is.null(ref.cell.type) && !(ref.cell.type %in% colnames(cnts)))
     stop(paste('Reference cell type', ref.cell.type, 'is not correct. Correct cell types are:',
                paste0(colnames(cnts), collapse = ', ') ))
  
  #Get freqs
  cnts[cnts == 0] <- 0.1
  freqs <- cnts/rowSums(cnts)

  # Get ilr
  psi <- coda.base::ilr_basis(ncol(freqs), type = "default") %>% 
    `rownames<-`(cnts %>% colnames())
  b <- log(freqs) %*% psi

  # PCA
  pca.res <- apply(b, 2, function(y) y - mean(y)) %>% 
    prcomp()
  pca.loadings <- psi %*% pca.res$rotation

  psi <- psi %*% pca.res$rotation
  b <- (log(freqs) %*% psi) %>% 
    apply(2, scale)

  idx.na = which(colSums(is.nan(b)) == 0)
  b <- b[,idx.na]
  psi <- psi[,idx.na]

  # Create dataframe
  b.df <- data.frame(b) %>% 
    mutate(groups = 1*groups)

  # Optimization
  if(method == 'lda') {
    if(is.null(ref.cell.type)) {  # pure linear regression which is proportional LDA when number of classes = 2
      res.lm <- lm(groups ~ ., data = b.df)
      w <- res.lm$coefficients[-1]
      w[is.na(w)] <- 0

      res.sum = summary(res.lm)
    } else {  # Quadratic optimization with the linear condition corresponding to the reference level
      X <- b
      Y <- groups * 1
      X2 <- cbind(X, 1)
      Rinv <- solve(chol(t(X2) %*% X2))
      dvec <- t(Y) %*% X2
      Amat <- t(psi[ref.cell.type,,drop=FALSE])
      Amat <- rbind(Amat, 0)
      res.qp <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)

      n.qp <- ncol(X2)
      w <- res.qp$solution[-n.qp]
    }
  } else if(method == 'svm') {
    b.model <- e1071::svm(groups~., data = b.df, kernel = "linear", scale = FALSE)
    # Get hyperplane parameters
    w <- t(b.model$SV) %*% b.model$coefs  # slope
  } else if(method == 'cda') {  # Canonical discriminant analysis
    model <- lm(b ~ groups)
    cda <- candisc::candisc(model, ndim=1)
    w <- cda$coeffs.raw
  } else if(method == 'cda.std') {
    pca.res <- apply(b, 2, function(y) y - mean(y)) %>% 
      prcomp()
    
    w <- pca.res$rotation %*% t(cor(groups, pca.res$x))
  }

  w <- w / sqrt(sum(w ^ 2))
  scores <- b %*% w
  
  if(mean(scores[groups]) < mean(scores[!groups])) w = -w

  loadings <- (psi %*% w) %>% 
    `names<-`(cnts %>% colnames())

  return(loadings)
}

#' This function produces the resampling dataset
#'
#' @param cnts Counts of cell types in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param n.iter Number of permutations
#' @param remain.groups TRUE - if the labels of groups remain, FALSE - if to produce the null distribution
#' @param replace.samples TRUE - is bootstrap on samples, FALAE - if to remain samples
#' @param seed random seed
#' @return Updated data frame with Z scores
#' @keywords internal
produceResampling <- function(cnts, groups, n.perm = 1000, seed = 239) {
  cnts.perm <- list()
  groups.perm <- list()
  set.seed(seed)

  for (i in 1:n.perm) {
    samples.tmp <- rownames(cnts) %>% split(groups[.]) %>%
      lapply(sample, replace=TRUE) %>% do.call(c, .)
    groups.tmp <- groups[samples.tmp]

    # sample cells
    cnts.resampling <- apply(cnts[samples.tmp,], 1, function(v) {
      if (sum(v) == 0) return(setNames(rep(0, length(v)), names(v)))

      sample(names(v), size=sum(v), replace=TRUE, prob = v / sum(v)) %>%
        {c(names(v), .)} %>% table() %>% {. - 1}
    }) %>% t()

    cnts.perm[[length(cnts.perm) + 1]] <- cnts.resampling
    groups.perm[[length(groups.perm) + 1]] <- groups.tmp
  }

  return(list(cnts = cnts.perm, groups = groups.perm))
}


#' @keywords internal
runCoda <- function(cnts, groups, n.seed=239, n.boot=1000, ref.cell.type=NULL, null.distr=FALSE, method="lda", n.cores=1, verbose=TRUE) {
  # Create datasets as
  samples.init <- produceResampling(cnts = cnts, groups = groups, n.perm = n.boot, seed = n.seed)
  loadings <- do.call(cbind, lapply(1:length(samples.init$cnts), function(ib) {
    getLoadings(samples.init$cnts[[ib]], samples.init$groups[[ib]], method=method)
  })) %>% 
    data.frame()

  # Calculate p-values of confidence interval by bootstrap
  tmp <- referenceSet(cnts, groups, p.thresh = 0.1)
  cell.list <- tmp$cell.list

  mean.list <- c()  # mean value of loadings in a list
  for (i.list in 1:length(cell.list)) {
    loadings.list <- loadings[cell.list[[i.list]],]
    mean.list <- c(mean.list, mean(unlist(loadings.list)))

    if (length(cell.list[[i.list]]) == 1) {  # if a list contains only one sample - it cannot be a reference group
      mean.list[i.list] <- 10
    }
  }

  # Define a cluster with reference cell type
  if (is.null(ref.cell.type)) {
    id.ref.cluster <- which(abs(mean.list) == min(abs(mean.list)))  # <- working version
  } else {
    id.ref.cluster <- -1
    for (i.list in 1:length(cell.list)) {
      if (ref.cell.type[1] %in% cell.list[[i.list]]) {
        id.ref.cluster <- i.list
      }
    }
    if (id.ref.cluster == -1) stop('Wrong name of reference cell type')
  }
  ref.cell.type <- cell.list[[id.ref.cluster]]


  # sorting of cell types and reference level
  cell.types.order <- c()
  for (i.list in order(-abs(mean.list - mean.list[id.ref.cluster]))) {
    cell.types.tmp <- cell.list[[i.list]] %>%
      .[order(-abs(rowMeans(loadings[.,, drop=FALSE]) - mean.list[id.ref.cluster]))]
    cell.types.order <- c(cell.types.order, cell.types.tmp)
  }
  loadings <- loadings[cell.types.order,]
  ref.load.level <- mean(loadings[ref.cell.type,])


  if (!null.distr) {
    tmp <- rowSums(ref.load.level > loadings)
    pval <- (pmin(tmp, ncol(loadings) - tmp) + 1) / (ncol(loadings) + 1) * 2
    names(pval) <- rownames(loadings)
  } else {
    loadings.perm <- do.call(cbind, plapply(1:n.boot, function(ib) {
      groups.perm <- sample(groups)
      names(groups.perm) <- names(groups)
      samples.perm <- produceResampling(cnts=cnts, groups=groups.perm, n.perm = 100, seed = n.seed + ib)

      do.call(cbind, lapply(1:length(samples.perm$cnts), function(ib) {
        getLoadings(samples.perm$cnts[[ib]], samples.perm$groups[[ib]])
      })) %>% rowMeans()
    }, n.cores=n.cores, progress=verbose, mc.preschedule=TRUE))

    loadings.stat <- rowMeans(loadings) - ref.load.level
    pval <- sapply(names(loadings.stat), function(s) {
      tmp <- sum(loadings.stat[s] > loadings.perm[s,])
      p <- (pmin(tmp, ncol(loadings) - tmp) + 1) / (ncol(loadings) + 1) * 2
      return(p)
    })
  }

  padj <- p.adjust(pval, method = "fdr")
  loadings <- loadings[names(pval),]
  return(list(padj=padj,
              loadings=loadings,
              pval=pval,
              ref.load.level=ref.load.level,
              ref.cell.type=ref.cell.type,
              cell.list=cell.list))
}

#' @keywords internal
referenceSet <- function(freqs, groups, p.thresh=0.05) {
  checkPackageInstalled("psych", cran=TRUE)
  freqs[freqs == 0] <- min(freqs[freqs != 0])/2
  cell.types <- colnames(freqs)
  cell.list <- lapply(cell.types, function(x) x)

  for (it in 1:(length(cell.list) - 2)) {

    mx <- matrix(0, nrow=length(cell.list), ncol=length(cell.list))
    for(i in 1:length(cell.list)){
      for(j in 1:length(cell.list)){
        if (j <= i) next

        ratio <- log(apply(freqs[,cell.list[[i]],drop=FALSE], 1, psych::geometric.mean )) -
          log(apply(freqs[,cell.list[[j]],drop=FALSE], 1, psych::geometric.mean ))
      
        mod <- lm(groups ~ ratio)
        res <- summary(mod)
        pval <- res$coefficients[2,4]
        mx[i, j] <- pval
      }
    }
    if(max(max(mx)) < p.thresh) break
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)
    cell.list[[length(cell.list) + 1]] <- c(cell.list[[i]], cell.list[[j]])
    cell.list <- cell.list[-c(i,j)]

    if(it == 1){
      mx.first <- mx
      rownames(mx.first) <- cell.types
      colnames(mx.first) <- cell.types
    }
  }
  len.cell.list <- sapply(cell.list, length)
  id.ref <- which(len.cell.list == max(len.cell.list))

  # Make mx.first symmetric
  for(i in 1:nrow(mx.first)){
    for(j in 1:ncol(mx.first)){
      if (j <= i) next
      mx.first[j, i] <- mx.first[i, j]
    }
  }

  return(list(ref.set = cell.list[[id.ref[1]]],
              cell.list = cell.list,
              mx.first = mx.first,
              mx.final = mx))
}

#' @keywords internal
sbpDiff <- function(freqs, groups){
  freqs[freqs == 0] <- min(freqs[freqs != 0])/2
  cell.types <- colnames(freqs)
  cell.list <- lapply(cell.types, function(x) x)

  sbp <- matrix(nrow = 0, ncol = length(cell.types), dimnames = list(c(), cell.types))
  for(it in 1:(length(cell.list) - 1)){

    mx <- matrix(0, nrow = length(cell.list), ncol = length(cell.list))
    for(i in 1:length(cell.list)){
      for(j in 1:length(cell.list)){
        if (j <= i) next

        ratio <- log(apply(freqs[,cell.list[[i]],drop=FALSE], 1, psych::geometric.mean )) -
          log(apply(freqs[,cell.list[[j]],drop=FALSE], 1, psych::geometric.mean ))

        mod <- lm(groups ~ ratio)
        res <- summary(mod)
        pval <- res$coefficients[2,4]
        mx[i, j] <- pval
      }
    }
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)

    sbp.tmp <- rep(0, length(cell.types))
    sbp.tmp[which(cell.types %in% cell.list[[i]])] <- 1
    sbp.tmp[which(cell.types %in% cell.list[[j]])] <- -1
    sbp <- rbind(sbp.tmp, sbp)
    cell.list[[length(cell.list) + 1]] <- c(cell.list[[i]], cell.list[[j]])
    cell.list <- cell.list[-c(i,j)]
  }

  return(t(sbp))
}

#' Check cell count data
#' 
#' @param d.counts Cell count table
#' @keywords internal
checkData <- function(d.counts){
  if(is.null(d.counts)) 
    stop('Cell count matrix is not provided')
  if(nrow(d.counts) == 0) 
    stop('Sample size is zero')
  if(ncol(d.counts) == 0) 
    stop('Cell count information is missed')
}

#' Check groups
#' 
#' @param d.groups Groups variable for samples
#' @keywords internal
checkGroups <- function(d.groups){
  if(is.null(d.groups)) 
    stop('Groups are not provided')
  if(length(d.groups) == 0) 
    stop('Group information is missed')
}

#' Check cell count data and groups
#' 
#' @param d.counts Cell count table
#' @param d.groups Groups variable for samples
#' @keywords internal
checkDataGroups <- function(d.counts, d.groups){
  checkData(d.counts)
  checkGroups(d.groups)
  if(nrow(d.counts) != length(d.groups)) 
    stop('Sample size in cell count matrix and in group variable is not the same')
}

#' Check sequential binary partitions(sbp)
#' 
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkSbp <- function(sbp){
  if(is.null(sbp)) 
    stop('Sequential binary partitions are not provided')
  if(nrow(sbp) == 0)
    stop('No partitions are provided')
}

#' Check sequential binary partitions(sbp) for the whole set of balances
#' 
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkSbpWhole <- function(sbp){
  checkSbp(sbp)
  if(nrow(sbp) < 2) 
    stop('Balances require at least 2 cell types')
  if((nrow(sbp) - ncol(sbp)) != 1) 
    stop('Wrong number of balances provided')
}

#' Check the agreement between sbp and cell count table
#' 
#' @param d.counts Cell count table
#' @param sbp Sbp matrix: each row - partition
#' @keywords internal
checkDataSbp <- function(d.counts, sbp){
  checkData(d.counts)
  checkSbp(sbp)
  if(ncol(sbp) != ncol(d.counts))
    stop('Number of cell types in data and sbp should be equal')
}

#' Check data and list of cells
#' 
#' @param d.counts Cell count table
#' @param cells List of cell types of interest
#' @keywords internal
checkDataAndCells <- function(d.counts, cells){
  if(is.null(cells))
    stop('Cell types are not provided')
  if(length(cells) == 0)
    stop('Cell types are not provided')
  if(sum(cells %in% colnames(d.counts)) != length(cells))
    stop('Some cell types do not exist in the data')
}

#' Check Tree
#' 
#' @param tree phylo tree
#' @keywords internal
checkTree <- function(tree){
  if(is.null(tree))
    stop('Tree is not provided')
}

#' Check agreement between Tree and the Cell count table
#' 
#' @param d.counts Cell count table
#' @param tree phylo tree
#' @keywords internal
checkDataTree <- function(d.counts, tree){
  checkTree(tree)
  checkData(d.counts)
  checkDataAndCells(d.counts, tree$tip.label)
}

#' Construct the canonical tree
#'
#' @param cnts Table with cell type counts
#' @param groups Groups variable for samples
#' @return phylo tree
#' @keywords internal
constructTree <- function(cnts, groups, partition.thresh = 0){
  checkDataGroups(cnts, groups)
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
    #  Data for the current subset of cells
    d.tmp <- cnts[, colnames(cnts) %in% unsolved.cells[[id.bal]]]
    
    # Define the most contrast balance
    can.loadings <- getLoadings(d.tmp, groups)
    cells.tmp <- rownames(can.loadings)
    
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
#' @param log.f Logarithms of frequencies
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
    #  Data for the current subset of cells
    d.tmp <- cnts[, colnames(cnts) %in% unsolved.cells[[id.bal]]]
    
    # Find the best partition
    sbp.best <- bestPartition(d.tmp, groups)
    
    # Get cell types from opposite sides of the principal balance
    cells.tmp <- colnames(d.tmp)
    
    cells.plus <- cells.tmp[sbp.best > 0]
    cells.minus <- cells.tmp[sbp.best < 0]
    
    sbp.cda[cells.minus, id.bal] <- -1
    sbp.cda[cells.plus, id.bal] <- 1
    
    if(length(cells.minus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.minus
    }
    
    if(length(cells.plus) > 1){
      unsolved.cells[[length(unsolved.cells) + 1]] <- cells.plus
    }
    
  }
  
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

#' Generate one random partition of given length
#'
#' @param bal.len length of partition
#' @return random partition - vector with values from {1, -1}
#' @keywords internal
getRndPartition <- function(bal.len){
  if(bal.len < 2)
    stop('Balance cannot be generated')
  
  bal.values <- c(0)
  while(length(unique(bal.values)) == 1){
    bal.values <- sample(c(-1, 1), bal.len, replace = TRUE)
  }
  return(bal.values)
}


#' Generate random sbp
#'
#' @param n.cell.types number of cell types
#' @return list with binary and normalized sbp matrices (partitions are in rows)
#' @keywords internal
getRndSbp <- function(n.cell.types, n.seed=239){
  if(n.cell.types < 2)
    stop('For a partition at least 2 cell types should be provided')
  
  set.seed(n.seed)
  sbp.bin <- matrix(1, 1, n.cell.types)
  bal.stack <- c(1) # stack of untreated balances
  while(length(bal.stack) > 0){
    
    parent.bal <- sbp.bin[bal.stack[1],]  # parental balance
    bal.stack <- bal.stack[-1]
    
    for(bal.direction in c(-1, 1)){
      if(sum(parent.bal == bal.direction) < 2) next
      tmp.bal <- parent.bal * 0
      tmp.bal[parent.bal == bal.direction] <- getRndPartition(sum(parent.bal == bal.direction))
      sbp.bin <- rbind(sbp.bin, tmp.bal)
      bal.stack <- c(bal.stack, nrow(sbp.bin))
      
    }
  }
  sbp.bin <- sbp.bin[-1,] # remove initialization
  
  sbp.norm <- sbp.bin
  for(i in 1:nrow(sbp.norm)){
    n1 <- sum(sbp.norm[i,] == 1)
    n2 <- sum(sbp.norm[i,] == -1)
    sbp.norm[i,sbp.norm[i,] == 1] <- 1/n1
    sbp.norm[i,sbp.norm[i,] == -1] <- -1/n2
    sbp.norm[i,] <- sbp.norm[i,] * sqrt(n1*n2/(n1+n2))
  }
  
  return(list(bin = sbp.bin, norm = sbp.norm))
}


#' Get balances based on the random ilr basis
#'
#' @param d.counts Cell count table
#' @param n.seed Random seed
#' @return Matrix with balances in columns
#' @keywords internal
getRndBalances <- function(d.counts, n.seed=239){
  checkData(d.counts)
  
  n.cell.types <- ncol(d.counts)
  sbp <- getRndSbp(n.cell.types, n.seed)
  
  log.freq <- getLogFreq(d.counts)
  
  psi <- as.matrix(t(sbp$norm))
  rownames(psi) <- colnames(d.counts)
  
  balances <- as.matrix(log.freq) %*% psi
  balances.norm <-  apply(balances, 2, function(y) y - mean(y))
  
  return(list(values = balances, norm = balances.norm, psi = psi))
}


#' Get log frequencies from the cell count table
#' @keywords internal
getLogFreq <- function(d.counts){
  checkData(d.counts)
  
  data.norm <- d.counts
  data.norm[data.norm == 0] <- 1
  
  log.freq <- t(apply(data.norm, 1, function(y) log(y/sum(y))))
  return(log.freq)
}

#' Remove the strongest group effect from balances
#'
#' @param d.used Correct values of balances
#' @param d.groups if provided then resampling controls presence of both groups in a new dataset
#' @param n.seed Random seed
#' @param thresh.pc.var percentage of variance which should be characterized by PSc
#' @return Balances without group effect
#' @keywords internal
removeGroupEffect <- function(d.used, d.groups, thresh.pc.var = 0.95){
  checkDataGroups(d.used, d.groups)
  
  pca.res <- prcomp(d.used, center = FALSE, scale. = FALSE)
  
  # Calculate cummularive variance explained by PCs
  expl.var <- cumsum(pca.res$sdev^2/sum(pca.res$sdev^2))
  n.pc.var <- sum(expl.var < thresh.pc.var) + 1
  n.pc.var <- ncol(pca.res$x) - 1
  
  for (i in 0:n.pc.var) {
    stop.loop <- TRUE
    d.working <- pca.res$x[,1:(n.pc.var - i)] # data to cda
    
    model<-lm(d.working ~ d.groups)
    tryCatch(suppressWarnings(cda <- candisc::candisc(model)), error = function(e) { stop.loop <<- FALSE})
    if(stop.loop) { break }
  }
  n.pc.var <- n.pc.var - i
  
  model<-lm(d.working ~ d.groups)
  cda <- candisc::candisc(model)
  
  cda.rotation <- cda$structure # already normalized
  
  d.scores <- d.working %*% cda.rotation  # Projection values
  d.part <- d.scores %*% t(cda.rotation)  # Projection vectors
  df.remain <- d.working - d.part  # Orthogonal part
  
  d.used.part <- d.part  %*% t(pca.res$rotation[,1:n.pc.var])
  d.used.remain <- df.remain %*% t(pca.res$rotation[,1:n.pc.var])
  
  cumm.rotation <- pca.res$rotation[,1:n.pc.var] %*% cda.rotation
  return(list(rotation = cumm.rotation,
              used.part = d.used.part,
              remain = d.used.remain,
              scores = d.scores,
              res = cda))
}
