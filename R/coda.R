#' Get loadings by different methods
#'
#' @param cnts Counts of cell typer in samples. Rows - samples, columns - cell types
#' @param groups Vector with boolean values. TRUE - sample in the case group, FALSE - sample in the control group
#' @param criteria Method to get loadings
#' @param ref.cell.type Reference cell type
#' @return Updated data frame with Z scores
#' @keywords internal
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

      res.sum = summary(res.lm)
      # print(res.sum)
      # w[(res.sum$coefficients[-1,4] > 0.01) & (res.sum$coefficients[-1,4] != min(res.sum$coefficients[-1,4]))] = 0
      # w[(res.sum$coefficients[-1,4] != min(res.sum$coefficients[-1,4]))] = 0

      # # Compare with lda - proportional!
      # res.lda <- lda(groups ~ ., data = b.df)
      # w / res.lda$scaling

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
    # samples.tmp <- sample(rownames(cnts), nrow(cnts), replace = TRUE)
    groups.tmp <- groups[samples.tmp]

    # Check that both groups are presented
    # while ((sum(groups.tmp) == 0) || (sum(!groups.tmp) == 0)) {
    #   # Create samples by bootstrap
    #   samples.tmp <- sample(rownames(cnts), nrow(cnts), replace = TRUE)
    #   groups.tmp <- groups[samples.tmp]
    # }

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
runCoda <- function(cnts, groups, n.seed=239, n.boot=1000, ref.cell.type=NULL, null.distr=FALSE) {
  # Create datasets as
  samples.init <- produceResampling(cnts = cnts, groups = groups, n.perm = n.boot, seed = n.seed)
  loadings <- do.call(cbind, lapply(1:length(samples.init$cnts), function(ib) {
    getLoadings(samples.init$cnts[[ib]], samples.init$groups[[ib]])
  }))


  # Calculate p-values of confidence interval by bootstrap
  tmp <- referenceSet(cnts, groups, p.thresh = 0.1)
  cell.list <- tmp$cell.list


  mean.list <- c()  # mean value of loadings in a list
  for (i.list in 1:length(cell.list)) {
    loadings.list <- loadings[cell.list[[i.list]],]
    mean.list <- c(mean.list, mean(c(loadings.list)))

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


    # samples.perm <- samples.init
    # loadings.perm <- do.call(cbind, lapply(1:n.boot, function(ib) {
    #   groups <- sample(samples.perm$groups[[ib]])
    #   names(groups) <- names(samples.init$groups[[ib]])
    #   getLoadings(samples.init$cnts[[ib]], groups)
    # }))

    loadings.perm <- do.call(cbind, plapply(1:n.boot, function(ib) {
      groups.perm <- sample(groups)
      names(groups.perm) <- names(groups)
      samples.perm <- produceResampling(cnts=cnts, groups=groups.perm, n.perm = 100, seed = n.seed + ib)

      do.call(cbind, lapply(1:length(samples.perm$cnts), function(ib) {
        getLoadings(samples.perm$cnts[[ib]], samples.perm$groups[[ib]])
      })) %>% rowMeans()
    }, n.cores=60, progress=TRUE, mc.preschedule=TRUE))

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
        # print('---')
        # print(freqs)
        # print(ratio)
        # print(groups)
        mod <- lm(groups ~ ratio)
        res <- summary(mod)
        pval <- res$coefficients[2,4]
        mx[i, j] <- pval
        # mx[j, i] <- pval
      }
    }
    if(max(max(mx)) < p.thresh) break
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)
    # print(c(i, j))
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
        # mx[j, i] <- pval
      }
    }
    i <- which(rowSums(mx == max(mx)) == 1)
    j <- which(colSums(mx == max(mx)) == 1)

    sbp.tmp <- rep(0, length(cell.types))
    sbp.tmp[which(cell.types %in% cell.list[[i]])] <- 1
    sbp.tmp[which(cell.types %in% cell.list[[j]])] <- -1
    sbp <- rbind(sbp.tmp, sbp)
    # print(c(i, j))
    cell.list[[length(cell.list) + 1]] <- c(cell.list[[i]], cell.list[[j]])
    cell.list <- cell.list[-c(i,j)]

  }

  return(t(sbp))
}
