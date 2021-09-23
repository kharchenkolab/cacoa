#' @description Transform sequential binary partition (sbp) into balances
#'
#' @param d.counts Cell count table
#' @param sbp.pop matrix with binary partiotions in rows
#' @return matrix of balances in columns; rows are samples
sbp2balances <- function(d.counts, sbp.pop){
  checkDataSbp(d.counts, sbp.pop)
  # ---------

  partition.norm <- c()
  for(i in 1:nrow(sbp.pop)){
    partition <- sbp.pop[i,]
    n1 <- sum(partition == 1)
    n2 <- sum(partition == -1)
    partition[partition == 1] <- 1/n1
    partition[partition == -1] <- -1/n2
    partition <- partition * sqrt(n1*n2/(n1+n2))
    partition.norm <- cbind(partition.norm, partition)
  }

  balances <- as.matrix(getLogFreq(d.counts)) %*% as.matrix(partition.norm)
  return(balances)
}


#' @description Genetic algorithm to obtaint the most contrast balance between cell types
#'
#' @param d.counts Cell count table
#' @param d.groups Group variable
#' @param n.pop number of individuals in the population
#' @param n.epoch number of epochs
#' @param n.seed seed id
#' @return matrix with final population, rows - individuals, top individual - the best partition
gaPartition <- function(d.counts, d.groups, n.pop = 1000, n.epoch = 100, n.seed = 239){
  checkDataGroups(d.counts, d.groups)
  set.seed(n.seed)
  # ---------

  n.parents <- round(n.pop / 3)
  n.children <- n.pop - n.parents

  n.cell.types <- ncol(d.counts)

  # Initialisation
  sbp.pop <- matrix(sample(c(0, -1, 1), n.cell.types * n.pop, replace = TRUE), nrow=n.pop)
  idx.partition <- (rowSums(sbp.pop == 1) > 0) & (rowSums(sbp.pop == -1) > 0)
  sbp.pop <- sbp.pop[idx.partition,]

  # Genetic Algorithm
  for(iter in 1:n.epoch){

    balances <- sbp2balances(d.counts, sbp.pop)

    # bal.mean.dist <- colMeans(balances[d.groups,]) - colMeans(balances[!d.groups,])

    # bal.mean.dist <- apply(balances[d.groups,],2,min) - apply(balances[!d.groups,],2,max)

    # # Wilcox
    # n1 = sum(d.groups)
    # n2 = sum(!d.groups)
    # r = apply(balances, 2, rank)
    # r1 = colSums(r[d.groups,])
    # r2 = colSums(r[!d.groups,])
    #
    # u1 = n1*n2 + n1 * (n1+1)/2 - r1
    # u2 = n1*n2 + n2 * (n2+1)/2 - r2
    # bal.mean.dist = pmin(u1, u2)

    bal.mean.dist = c()
    for(ibal in 1:ncol(balances)){
      bal.mean.dist = c(bal.mean.dist, f.stat(balances[d.groups, ibal], balances[!d.groups, ibal])$stat)
    }
    # print(bal.mean.dist)


    # ref = f.stat(balances[d.groups], balances[!d.groups])
    # bal.mean.dist = ref$stat

    sbp.pop <- sbp.pop[rev(order(bal.mean.dist)),]
    sbp.pop <- sbp.pop[1:min(nrow(sbp.pop), n.parents),]

    if (iter == n.epoch) break

    sbp.children <- c()
    for(i in 1:n.cell.types){
      val.probs <- c(sum(sbp.pop[,i] == 0), sum(sbp.pop[,i] == -1), sum(sbp.pop[,i] == 1)) + 1
      val.probs <- val.probs / sum(val.probs)
      sbp.children <- cbind(sbp.children, sample(c(0, -1, 1), n.children, replace = TRUE, prob = val.probs) )
    }
    sbp.pop <- rbind(sbp.pop, sbp.children)
    sbp.pop <- unique(sbp.pop)
    idx.partition <- (rowSums(sbp.pop == 1) > 0) & (rowSums(sbp.pop == -1) > 0)
    sbp.pop <- sbp.pop[idx.partition,]
  }
  colnames(sbp.pop) <- colnames(d.counts)
  # print(bal.mean.dist[1])
  return(sbp.pop)
}


#' @description Genetic algorithm to obtaint the most contrast balance including all cell types
#'
#' @param d.counts Cell count table
#' @param d.groups Group variable
#' @param n.pop number of individuals in the population
#' @param n.epoch number of epochs
#' @param n.seed seed id
#' @return matrix with final population, rows - individuals, top individual - the best partition
gaPartitionAll <- function(d.counts, d.groups, n.pop = 1000, n.epoch = 100, n.seed = 239){
  checkDataGroups(d.counts, d.groups)
  set.seed(n.seed)
  # ---------

  n.parents <- round(n.pop / 3)
  n.children <- n.pop - n.parents

  n.cell.types <- ncol(d.counts)

  # Initialisation
  sbp.pop <- matrix(sample(c(-1, 1), n.cell.types * n.pop, replace = TRUE), nrow=n.pop)
  idx.partition <- (rowSums(sbp.pop == 1) > 0) & (rowSums(sbp.pop == -1) > 0)
  sbp.pop <- sbp.pop[idx.partition,]

  # Genetic Algorithm
  for(iter in 1:n.epoch){

    balances <- sbp2balances(d.counts, sbp.pop)

    # Objective function
    bal.mean.dist <- apply(balances[d.groups,],2,min) - apply(balances[!d.groups,],2,min)

    sbp.pop <- sbp.pop[rev(order(bal.mean.dist)),]
    sbp.pop <- sbp.pop[1:min(nrow(sbp.pop), n.parents),]

    if (iter == n.epoch) break

    sbp.children <- c()
    for(i in 1:n.cell.types){
      val.probs <- c(sum(sbp.pop[,i] == -1), sum(sbp.pop[,i] == 1)) + 1
      val.probs <- val.probs / sum(val.probs)
      sbp.children <- cbind(sbp.children, sample(c(-1, 1), n.children, replace = TRUE, prob = val.probs) )
    }
    sbp.pop <- rbind(sbp.pop, sbp.children)
    sbp.pop <- unique(sbp.pop)
    idx.partition <- (rowSums(sbp.pop == 1) > 0) & (rowSums(sbp.pop == -1) > 0)
    sbp.pop <- sbp.pop[idx.partition,]
  }
  return(sbp.pop)
}


#' @description Get the most contrast balance based on genetic algorithm
#'
#' @param d.counts Cell count table
#' @param d.groups Group variable
#' @param n.pop number of individuals in the population
#' @param n.epoch number of epochs
#' @param n.seed seed id
#' @param thresh.parent fraction of parents for a next epoch
#' @return list of significant cell types and two contrast sets of cells
getSignifCellTypes <- function(d.counts, d.groups,
                               n.pop = 1000, n.epoch = 100, n.seed = 239,
                               thresh.parent = 1/4){
  checkDataGroups(d.counts, d.groups)
  # ---------

  sbp.pop <- gaPartition(d.counts, d.groups, n.pop = n.pop, n.epoch = n.epoch, n.seed = n.seed)

  # colnames(sbp.pop) <- colnames(d.counts)
  # sbp.pop.cut <- sbp.pop[1:(round(nrow(sbp.pop) * thresh.parent) + 1),]
  # # print(sbp.pop.cut)
  # # sbp.pop.cut
  #
  # mean.freq <- colMeans(abs(sbp.pop.cut))
  # names(mean.freq) <- colnames(d.counts)
  # super.signif.types <- names(mean.freq)[rev(order(mean.freq))]
  # super.signif.types <- super.signif.types[1]
  # # print(super.signif.types)
  #
  # for(i in 1:nrow(sbp.pop.cut)){
  #   sbp.pop.cut[i,] <- sbp.pop.cut[i,] * sign(sbp.pop.cut[i,super.signif.types])
  # }
  # # print(sbp.pop.cut)
  #
  # # mean.freq <- colMeans(sbp.pop.cut)
  # # names(mean.freq) <- colnames(d.counts)
  # # signif.types <- colnames(d.counts)[abs(mean.freq) > thresh.freq]
  # #
  # # return(list(x = signif.types[mean.freq[signif.types] > 0],
  # #             y = signif.types[mean.freq[signif.types] < 0],
  # #             sbp.pop = sbp.pop))
  #
  # mean.freq <- colMeans(sbp.pop.cut)
  # names(mean.freq) <- colnames(d.counts)
  # signif.types <- colnames(d.counts)[abs(mean.freq) > thresh.freq]

  return(list(x = colnames(sbp.pop)[sbp.pop[1,] > 0],
              y = colnames(sbp.pop)[sbp.pop[1,] < 0],
              sbp.pop = sbp.pop))
}


#' @description Get caconical decomposisiotn of balances based
#'
#' @param d.counts Cell count table
#' @param d.groups Group variable
#' @param n.dim number of components to return
#' @return Caconical components
getCanonicalContrasts <- function(d.counts, d.groups, n.dim = 2){
  checkDataGroups(d.counts, d.groups)
  if(n.dim < 1)
    stop('Incorrect number of components provided')
  # ---------

  scores <- c()
  cell.used <- c()
  cell.used.x <- c()
  cell.used.y <- c()
  n.scores <- 2
  for(i in 1:n.scores){
    sign.types <- getSignifCellTypes(d.counts[, !(colnames(d.counts) %in% cell.used)], d.groups)
    score.tmp <- calcBalancesOnCellSets(d.counts, sign.types$x, sign.types$y)
    scores <- cbind(scores, score.tmp)
    # print(sign.types$x)
    # print(sign.types$y)
    cell.used <- c(cell.used, sign.types$y, sign.types$x)
    cell.used.x <- c(cell.used.x, sign.types$x)
    cell.used.y <- c(cell.used.y, sign.types$y)
  }

  colnames(scores) <- paste('Score', 1:ncol(scores), sep = '')
  scores
  return(list(scores = scores,
              cell.used = cell.used,
              cell.used.x = cell.used.x,
              cell.used.y = cell.used.y))
}


#' @description Get the most contrast balance based on genetic algorithm iteratively(!)
#'
#' @param d.counts Cell count table
#' @param d.groups Group variable
#' @param n.pop number of individuals in the population
#' @param n.epoch number of epochs
#' @param n.seed seed id
#' @param thresh.parent fraction of parents for a next epoch
#' @param p.val.thresh p-value threshold
#' @return list of significant cell types and two contrast sets of cells
getSignifCellTypesIteratively <- function(d.counts, d.groups,
                                          n.pop = 1000, n.epoch = 100, n.seed = 239,
                                          thresh.parent = 1/4, p.val.thresh = 0.01){
  checkDataGroups(d.counts, d.groups)
  # ---------

  p.vals <- c()
  cell.used <- c()
  cell.used.iter <- list()
  scores <- c()
  while(ncol(d.counts) - length(cell.used) > 1){
    if(ncol(d.counts) - length(cell.used) == 2){
      cell.unused <- setdiff(colnames(d.counts), cell.used)
      score.tmp <- calcBalancesOnCellSets(d.counts, cell.unused[1], cell.unused[2])
    }else{
      sign.types <- getSignifCellTypes(d.counts[, !(colnames(d.counts) %in% cell.used)],
                                      d.groups, n.pop = n.pop, n.epoch = n.epoch, n.seed = n.seed,
                                      thresh.parent = thresh.parent)
      score.tmp <- calcBalancesOnCellSets(d.counts, sign.types$x, sign.types$y)
    }


    w.test.res <- wilcox.test(score.tmp[d.groups], score.tmp[!d.groups])
    p.vals <- c(p.vals, w.test.res$p.value)

    # print(score.tmp)

    scores <- cbind(scores, score.tmp)

    cell.used <- c(cell.used, sign.types$y, sign.types$x)
    cell.used.iter[[length(cell.used.iter) + 1]] <- c(sign.types$y, sign.types$x)

  }
  # print(p.vals)
  # p.vals.adj <- p.adjust(p.vals, method = 'fdr')
  p.vals.adj <- p.vals

  cell.used.merge <- c()
  for(i in 1:length(cell.used.iter)){
    if(p.vals.adj[i] > p.val.thresh) {next}
    cell.used.merge <- c(cell.used.merge, cell.used.iter[[i]])
  }

  return(list(cells = cell.used.merge,
              cell.iter = cell.used.iter,
              scores = scores,
              p.vals = p.vals,
              p.vals.adj = p.vals.adj))
}


f.stat <- function(b1, b2)
{
  n1 = length(b1)
  n2 = length(b2)


  ss_tot = sum()

  m_tot <- mean(c(b1, b2))
  ss_tot <- sum((c(b1, b2) - m_tot)^2)
  ss_tot

  # Within-group variance
  ss_wg <- sum((b1 - mean(b1))^2) + sum((b2 - mean(b2))^2)
  ss_wg

  # Between-group variance
  ss_bg <- ss_tot - ss_wg
  ss_bg

  ms_b = ss_bg / (2 - 1)
  ms_w = ss_wg / (n1 + n2 - 2)

  res = list(stat = ms_b/ms_w,
             pval = 	pf(ms_b/ms_w, 2 - 1, n1 + n2 - 2, lower.tail=F))

  return(res)
}


