#' Get loadings by different methods
#' @param cnts Counts of cell typer in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param criteria Method to get loadings
#' @param ref.cell.type Reference cell type
#' @return Updated data frame with Z scores
getLoadings <- function(cnts, groups, criteria = 'lda', ref.cell.type = NULL) {
  discriminant.methods <- c('lda', 'svm', 'cda', 'cda.std')
  if(!(criteria %in% discriminant.methods)) stop(paste('The discriminant method', criteria, 'is not supported'))
  if(!is.null(ref.cell.type) && !(ref.cell.type %in% colnames(cnts))) 
     stop(paste('Reference cell type', ref.cell.type, 'is not correct. Correct cell types are:', 
                paste0(colnames(cnts), collapse = ', ') ))
  #Get freqs
  cnts[cnts == 0] <- 0.1
  # cnts[cnts == 0] <- 1
  freqs <- cnts/rowSums(cnts)
  
  # Get ilr
  psi <- coda.base::ilr_basis(ncol(freqs), type = "default")
  rownames(psi) <- colnames(cnts)
  b <- log(freqs) %*% psi
  
  # PCA
  b.norm <-  apply(b, 2, function(y) y - mean(y))
  pca.res <- prcomp(b.norm)
  pca.loadings <- psi %*% pca.res$rotation
  
  psi <- psi %*% pca.res$rotation
  b <- log(freqs) %*% psi
  b <- apply(b, 2, scale)
  
  
  idx.na = which(colSums(is.nan(b)) == 0)
  b <- b[,idx.na]
  psi <- psi[,idx.na]
  
  # Create dataframe
  b.df <- data.frame(b)
  b.df$groups <- 1*groups
  
  # Optimization
  if(criteria == 'lda') {
    if(is.null(ref.cell.type)) {  # pure linear regression which is proportional LDA when number of classes = 2
      res.lm <- lm(groups ~ ., data = b.df)
      w <- res.lm$coefficients[-1] 
      w[is.na(w)] <- 0
      
      # # Compare with lda - proportional!
      # res.lda <- lda(groups ~ ., data = b.df)
      # w / res.lda$scaling
      
    } else {  # Quadratic optimization with the linear condition corresponding to the reference level
      X <- b
      Y <- groups * 1
      X2 <- cbind(X, 1)
      Rinv <- solve(chol(t(X2) %*% X2))
      dvec <- t(Y) %*% X2
      Amat <- t(psi[ref.cell.type,,drop=F])
      Amat <- rbind(Amat, 0)
      res.qp <- quadprog::solve.QP(Dmat = Rinv, factorized = TRUE, dvec = dvec, Amat = Amat, bvec = 0, meq = 1)  
      
      n.qp <- ncol(X2)
      w <- res.qp$solution[-n.qp]
      
      # # Compare with lm
      # res.lm <- lm(groups ~ ., data = b.df)
      # res.qp$unconstrained.solution[-n.qp] / res.lm$coefficients[-1]
    }   
  } else if(criteria == 'svm') {
    # ---- SVM
    b.model <- e1071::svm(groups~., data = b.df, kernel = "linear", scale = FALSE)
    # Get hyperplane parameters
    w <- t(b.model$SV) %*% b.model$coefs  # slope
    # b.svm <- -b.model$rho # intercept
    # # create new score
    # v <- b %*% w + b.svm
  } else if(criteria == 'cda') {  # Canonical discriminant analysis
    model <- lm(b ~ groups)
    cda <- candisc::candisc(model, ndim=1)
    # w <- cda$structure
    w <- cda$coeffs.raw
  } else if(criteria == 'cda.std') {
    b.norm <-  apply(b, 2, function(y) y - mean(y))

    # PCA
    pca.res <- prcomp(b.norm)
    pca.loadings <- psi %*% pca.res$rotation

    # ANNA: please remain the following commented code, 
    # in case someone wants to understand what is going on, and what was attempted
    # # CDA
    # df.pca <- as.data.frame(pca.res$x)
    # 
    # model <- lm(pca.res$x ~ groups)
    # cda <- candisc::candisc(model, ndim=1)
    # 
    # w <- pca.res$rotation  %*% as.matrix(cda$structure)  # standardized regression coefficients - just correlations
    # # w <- pca.res$rotation  %*% as.matrix(cda$coeffs.raw)  # pure regression coefficients
    
    # # Properties
    # # 1
    # b.pca = pca.res$x
    # scores <- cda$scores$Can1
    # cor(groups, b.pca) / cor(scores, b.pca)  # the same number
    # # 1
    # t(cda$structure) / cor(groups, b.pca)  # the same number
    # # What we can do: 
    # w.tmp <- pca.res$rotation %*% t(cor(groups, b.pca))
    # w.tmp / w  # the same number
    
    # New version without CDA
    # To makes things comparable between bootstrap iterations,
    # we need standardized regression coefficients, 
    # which are just correlations, when regressors are independent, as after PCA.
    # Standardization is needed,
    # because a regression coefficient reflects the variance of the corresponding regressor, 
    # and it significantly depends on bootstrap subsample.
    w <- pca.res$rotation %*% t(cor(groups, pca.res$x))
    
  } 

  w <- w / sqrt(sum(w ^ 2))
  
  scores <- b %*% w
  if(mean(scores[groups]) < mean(scores[!groups])) w = -w
  # scores <- b %*% w
  # print(c(mean(scores[groups]), mean(scores[!groups])))
  
  loadings <- psi %*% w
  names(loadings) <- colnames(cnts)
  
  return(loadings)
}

#' This function producec datasent for resampling
#' @param cnts Counts of cell typer in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param n.iter Number of permutations
#' @param remain.groups TRUE - if the labels of groups remain, FALSE - if to produce the null distribution
#' @param replace.samples TRUE - is bootstrap on samples, FALAE - if to remain samples
#' @param seed random seed
#' @return Updated data frame with Z scores
produceResampling <- function(cnts, groups, n.perm = 1000, remain.groups = TRUE, replace.samples = TRUE, seed = 239) {
  
  if(!(remain.groups %in% c(NULL, TRUE, FALSE))) stop('Parameter remain.groups should be in {NULL, TRUE, FALSE}')
  
  cnts.perm <- list()
  groups.perm <- list()
  set.seed(239)
  for(i in 1:n.perm) {
    # Bootstrap samples

    samples.tmp <- sample(rownames(cnts), nrow(cnts), replace = replace.samples)
    if (remain.groups) {  # Remain labels
      groups.tmp <- groups[samples.tmp]
    } else {  # Null distribution
      groups.tmp <- sample(groups, length(groups), replace = replace.samples)
      names(groups.tmp) <- samples.tmp
    }
    # Check that both groups are presented
    if((sum(groups.tmp) == 0) || (sum(!groups.tmp) == 0)) next
    
    # Bootstrap cell types
    cnts.resampling <- apply(cnts[samples.tmp,], 1, function(v) {
      table( c(names(v), sample(names(v), size = sum(v), replace = TRUE, prob = v / sum(v))) ) - 1   })
    
    cnts.tmp <- t(cnts.resampling)
    
    cnts.perm[[length(cnts.perm) + 1]] <- cnts.tmp
    groups.perm[[length(groups.perm) + 1]] <- groups.tmp
    
  }
  
  return(list(cnts = cnts.perm, groups = groups.perm))
}

runCoda <- function(cnts, groups, n.seed=239, n.boot=1000, ref.cell.type=NULL){
  
  # Apply bootstrap
  loadings <- lapply(1:n.boot, function(ib){
    # Create samples by bootstrap
    set.seed( n.seed+ib)
    samples.tmp <- sample(rownames(cnts), nrow(cnts), replace = T)
    groups.tmp <- groups[samples.tmp]
    
    
    # Check that both groups are presented
    while((sum(groups.tmp) == 0) || (sum(!groups.tmp) == 0)) {
      # Create samples by bootstrap
      samples.tmp <- sample(rownames(cnts), nrow(cnts), replace = T)
      groups.tmp <- groups[samples.tmp]
    }
    cnts.tmp <- cnts[samples.tmp,]
    
    init.tmp <- produceResampling(cnts = cnts.tmp, groups = groups.tmp, n.perm = 1,
                                  replace.samples = F,
                                  remain.groups = TRUE, seed = n.seed+ib)
    loadings.tmp <- getLoadings(init.tmp$cnts[[1]], init.tmp$groups[[1]])
  })
  
  # Calculate p-values
  loadings.init <- c()
  for(i in 1:length(loadings)){
    loadings.init <- cbind(loadings.init, loadings[[i]])
  }
  
  # ld <- loadings.init
  # threshold <- rowMeans(ld)[abs(rowMeans(ld)) == min(abs(rowMeans(ld)))]
  threshold <- 0
  
  # Calculate p-values of confidence interval by bootstrap
  if(is.null(ref.cell.type)){
    tmp <- sapply(1:nrow(loadings.init), function(i) sum(threshold > loadings.init[i,])) / ncol(loadings.init)
    pval <- apply((cbind(tmp, 1-tmp)), 1, min) * 2  
  } else {
    ref.level = abs(mean(loadings.init[ref.cell.type,]))
    tmp1 <- sapply(1:nrow(loadings.init), function(i) sum(loadings.init[i,] > ref.level)) / ncol(loadings.init)
    tmp2 <- sapply(1:nrow(loadings.init), function(i) sum(loadings.init[i,] < -ref.level)) / ncol(loadings.init)
    pval = (1 -apply((cbind(tmp1, tmp2)), 1, max)) * 2  
  }
  
  names(pval) <- rownames(loadings.init)
  
  # ----------------------------
  # Additional correction of p-values
  ld <- loadings.init
  # ld.means <- rowMeans(ld)
  
  idx <- order(abs(rowMeans(ld)))
  ld <- ld[idx,]
  pval <- pval[rownames(ld)]
  
  idx <- names(pval)[order(-pval)]
  ld <- ld[idx,]
  ld.means <- rowMeans(ld)
  pvals_tmp <- pval[idx]
  pvals_tmp[1] <- 1
  current.mean = ld.means[1]
  
  for(i in 2:length(ld.means)){
    # Fid the index of previous closest
    d.prev = abs(ld.means[1:(i-1)] - ld.means[i])
    j = which(d.prev == min(d.prev))
    
    tmp <- sum(ld[i,] > ld.means[j]) / ncol(ld)
    tmp <- min(tmp, 1-tmp)
    pvals_tmp[i] <- max(min(tmp, pvals_tmp[j]), pvals_tmp[i])
  }
  names(pvals_tmp) <- rownames(ld)
  # print(pvals_tmp)
  # Combining p-values
  p.names <- names(pval)
  pval <- max2(pval, pvals_tmp[p.names])
  names(pval) <- p.names
  # ----------------------------
  
  padj <- p.adjust(pval, method = "fdr")
  loadings.init = loadings.init[p.names,]
  return(list(padj=padj, loadings.init=loadings.init, pval=pval))
}


max2 <- function(x, y){
  res <- c()
  for(i in 1:length(x)){
    res = c(res, max(x[i], y[i]))
  }
  return(res)
}

