#' Perform re-sampling in two ways: fixed total count and bootstrap
#'
#' @param
#' @return Order of inner nodes
resampleContrast <- function(d.counts, d.groups, n.cell.counts = 500, n.seed = 239, n.iter = 1000){
  checkDataGroups(d.counts, d.groups)
  # --------
  set.seed(n.seed)
  rnd.seeds <- unique(round(runif(1.5*n.iter, max = 1000000)))

  cda.loadings <- getCdaLoadings(d.counts, d.groups)

  cda.resamples <- matrix(cda.loadings[colnames(d.counts),], ncol = 1)
  d.all <- c()
  n.skip.resampl <- 0
  for(iter in 1:n.iter){
    d.resampled <- resampleCounts(d.counts, d.groups = d.groups, n.tot.count = n.cell.counts, n.seed = rnd.seeds[iter])
    d.all <- rbind(d.all, d.resampled)

    cda.loadings = NULL
    tryCatch(expr = {
        cda.loadings <- getCdaLoadings(d.resampled, d.groups[rownames(d.resampled)], n.seed = rnd.seeds[iter])
      },error = function(e){})

    if(is.null(cda.loadings)) {
      n.skip.resampl <- n.skip.resampl + 1
      next
    }

    x <- cda.loadings[colnames(d.counts),]
    # if(length(cda.resamples) > 0){
    #   if(sum((sign(x) != sign(cda.resamples[,1])) * abs(cda.resamples[,1])) >
    #      sum((sign(-x) != sign(cda.resamples[,1])) * abs(cda.resamples[,1]))){
    #     x <- -x
    #   }
    # }
    cda.resamples <- cbind(cda.resamples, x)
  }

  if(n.skip.resampl > 0) warning(paste('Numer of skipped resamplings:', as.character(n.skip.resampl)))
  return(list(balances = cda.resamples, data = d.all))
}


#' Generate one random partition of given length
#'
#' @param bal.len length of partiotion
#' @return random partition - verctor with values from {1, -1}
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
#' @return list with binary and normalised sbp matrices (partitions are in rows)
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
getRndBalances <- function(d.counts, n.seed=239){
  checkData(d.counts)
  # ---------

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
getLogFreq <- function(d.counts){
  checkData(d.counts)
  # ---------

  data.norm <- d.counts
  data.norm[data.norm == 0] <- 1

  log.freq <- t(apply(data.norm, 1, function(y) log(y/sum(y))))
  return(log.freq)
}


#' Get Loadings from the caconical data analysis
#'
#' @param d.counts Cell count table
#' @param d.groups Group variable
#' @param n.seed Random seed
#' @return Vector with cda loadings
getCdaLoadings <- function(d.counts, d.groups, n.seed = 239){
  checkDataGroups(d.counts, d.groups)
  # ---------

  bal <- getRndBalances(d.counts, n.seed)

  # PCA
  pca.res <- prcomp(bal$norm)
  pca.loadings <- bal$psi %*% pca.res$rotation

  # CDA
  df.pca <- as.data.frame(pca.res$x)
  n.pc <- ncol(pca.loadings)

  for (i in 0:n.pc) {
    stop.loop <- TRUE
    model <- lm(as.matrix(df.pca[,1:(n.pc-i)]) ~ d.groups)
    tryCatch(suppressWarnings(cda <- candisc::candisc(model, ndim=1)), error = function(e) { stop.loop <<- FALSE})
    if(stop.loop)
      break
  }

  # Positive values of loadings are for disease(targer.level), negative - for control(ref.level)
  if (mean(cda$scores[d.groups,'Can1']) < mean(cda$scores[!d.groups,'Can1'])) {
    cda$structure <- -cda$structure
  }

  cda.loadings <- pca.loadings[,1:(n.pc-i)]  %*% as.matrix(cda$structure)
  # cda.loadings <- pca.loadings[,1:(n.pc-i)]  %*% as.matrix(cda$coeffs.raw)

  return(cda.loadings)
}

#' Resampling of cell type counts
#'
#' @param d.counts Cell count table
#' @param n.tot.count Total cell type counts
#' @param d.groups if provided then resampling controls presence of both groups in a new dataset
#' @param n.seed Random seed
#' @param f.bootstrap if to use the bootstrap, dafault (TRUE) means the use of bootstrap
#' @return Resampled dataset
resampleCounts<- function(d.counts, n.tot.count=500, d.groups = NULL, n.seed = 239, f.bootstrap = TRUE){

  if(is.null(d.groups))
    checkData(d.counts)
  else
    checkDataGroups(d.counts, d.groups)
  # ---------

  set.seed(n.seed)
  d.resampled <- c()
  for(i in 1:nrow(d.counts)){
    tmp <- rmultinom(1:ncol(d.counts), n.tot.count, d.counts[i,]+1)
    d.resampled <- rbind(d.resampled, t(tmp))
  }
  rownames(d.resampled) <- rownames(d.counts)

  if(f.bootstrap){
    if(is.null(d.groups)){
      d.resampled <- d.resampled[sample(nrow(d.counts), replace =T),]
    }else{
      tmp.idx <- sample(nrow(d.counts), replace =T)

      while(length(unique(d.groups[tmp.idx])) == 1) tmp.idx <- sample(nrow(d.counts), replace =T)

      d.resampled <- d.resampled[tmp.idx,]
    }
  }
  return(d.resampled)
}


#' Remove the strongest group effect from balances
#'
#' @param d.used Currect values of balances
#' @param d.groups if provided then resampling controls presence of both groups in a new dataset
#' @param n.seed Random seed
#' @param thresh.pc.var percentage of variance which should be characterised by PSc
#' @return Balances without group effect
removeGroupEffect <- function(d.used, d.groups, thresh.pc.var = 0.95){
  checkDataGroups(d.used, d.groups)
  # ---------

  pca.res <- prcomp(d.used, center = FALSE, scale. = FALSE)

  # Calculate cummularive variance explained by PCs
  expl.var <- cumsum(pca.res$sdev^2/sum(pca.res$sdev^2))
  n.pc.var <- sum(expl.var < thresh.pc.var) + 1
  n.pc.var <- ncol(pca.res$x) - 1

  # d.used.explained <- d.working %*% t(pca.res$rotation[,1:n.pc.var])

  for (i in 0:n.pc.var) {
    stop.loop <- TRUE
    # print(d.groups)
    d.working <- pca.res$x[,1:(n.pc.var - i)] # data to cda

    model<-lm(d.working ~ d.groups)
    tryCatch(suppressWarnings(cda <- candisc::candisc(model)), error = function(e) { stop.loop <<- FALSE})
    if(stop.loop) { break }
  }
  n.pc.var <- n.pc.var - i

  model<-lm(d.working ~ d.groups)
  cda <- candisc::candisc(model)

  cda.rotation <- cda$structure # already normalised

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


#' Calculate balances based on the contrast defines by two sets of cells
#'
#' @param d.counts Cell count matrix
#' @param cell.set1 Set#1
#' @param cell.set2 Set#1
#' @return Balance values
calcBalancesOnCellSets <- function(d.counts, cell.set1, cell.set2 = NULL){
  checkData(d.counts)
  checkDataAndCells(d.counts, cell.set1)
  if(is.null(cell.set2)) cell.set2 = setdiff(colnames(d.counts), cell.set1)
  checkDataAndCells(d.counts, cell.set2)
  # ---------

  log.freq <- getLogFreq(d.counts)
  n1 <- length(cell.set1)
  n2 <- length(cell.set2)
  if((n1 == 1) && (n2 == 1)) bal <- log.freq[,cell.set1] - log.freq[,cell.set2]
  if((n1 == 1) && (n2 != 1)) bal <- log.freq[,cell.set1] - rowMeans(log.freq[,cell.set2])
  if((n1 != 1) && (n2 == 1)) bal <- rowMeans(log.freq[,cell.set1])- log.freq[,cell.set2]
  if((n1 != 1) && (n2 != 1)) bal <- rowMeans(log.freq[,cell.set1]) - rowMeans(log.freq[,cell.set2])

  normalizing.coef <- sqrt(n1*n2/(n1+n2))
  return(bal * normalizing.coef)
}


calcWilcoxonTest <- function(d.counts, d.groups){
  d.freqs <- d.counts %>% magrittr::divide_by(rowSums(.))

  p.vals.balances <- c()
  for(cell.set1 in colnames(d.counts)){
    cell.set2 <- setdiff(colnames(d.counts), cell.set1)
    # cell.set2 <- 'AST-PP'
    bal <- calcBalancesOnCellSets(d.counts, cell.set1, cell.set2)
    x <- bal[d.groups]
    y <- bal[!d.groups]

    res <- wilcox.test(x, y)
    p.vals.balances <- c(p.vals.balances, res$p.value)
  }

  p.vals.freqs <- colnames(d.counts) %>% sapply(function(cell.type) {
    wilcox.test(d.freqs[d.groups, cell.type], d.freqs[!d.groups, cell.type])$p.value
  })

  p.vals.res <- cbind(p.vals.balances, p.vals.freqs)
  rownames(p.vals.res) <- colnames(d.counts)

  return(p.vals.res)
}

getCellSignificance <- function(balances){
  n.bal <- ncol(balances)
  n.plus <- rowSums(balances[,2:n.bal] > 0)
  n.minus <- rowSums(balances[,2:n.bal] < 0)
  perm.frac <- mapply(function(num1, num2) min(num1, num2), n.plus, n.minus) /
    mapply(function(num1, num2) max(num1, num2), n.plus, n.minus)

  perm.frac
}
