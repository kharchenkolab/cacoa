jaccardPwTop <- function(subsamples, top.thresh){
  jac.all = c()
  idxs = c()
  for(i in 1:length(subsamples)) {
    for(j in 1:length(subsamples)) {
      if (j <= i) next
      set1 <- rownames(subsamples[[i]])[rank(subsamples[[i]]$pvalue) <= top.thresh]
      set2 <- rownames(subsamples[[j]])[rank(subsamples[[j]]$pvalue) <= top.thresh]
      if((length(set1) != 0) && (length(set2) != 0)) {
        jac.all <- c(jac.all, length(intersect(set1, set2)) / length(unique(c(set1, set2))))
      } else {
        jac.all <- c(jac.all, 0)
      }
      idxs <- c(idxs, paste(names(subsamples)[i], names(subsamples)[j]))
    }
  }
  return(list(jac = jac.all, id = idxs))
}

jaccardPwPval <- function(subsamples, p.val.cutoff){
  jac.all = c()
  idxs = c()
  for(i in 1:length(subsamples)) {
    for(j in 1:length(subsamples)) {
      if (j <= i) next
      set1 <- rownames(subsamples[[i]])[subsamples[[i]]$padj <= p.val.cutoff]
      set2 <- rownames(subsamples[[j]])[subsamples[[j]]$padj <= p.val.cutoff]
      if((length(set1) != 0) && (length(set2) != 0)) {
        jac.all <- c(jac.all, length(intersect(set1, set2)) / length(unique(c(set1, set2))))
      } else {
        jac.all <- c(jac.all, 0)
      }
      idxs <- c(idxs, paste(names(subsamples)[i], names(subsamples)[j]))
    }
  }
  return(list(jac = jac.all, id = idxs))
}


estimateStabilityPerCellType <- function(de.res,
                                         top.n.genes,
                                         p.val.cutoff) {

  data.all <- data.frame()
  for(cell.type in names(de.res)){
    # print(cell.type)
    subsamples <- de.res[[cell.type]]$subsamples
    
    # Remove some subsamples due to the min.cell.counts
    # coomare resampling results with "initial"
    jacc.init <- c()
    for(subs.name in names(subsamples)){
      subsamples.tmp = list(de.res[[cell.type]]$subsamples[[subs.name]],
                            de.res[[cell.type]]$res)
      jacc.pw.tmp <- jaccardPwTop(subsamples.tmp, 200)
      jacc.init <- c(jacc.init, jacc.pw.tmp$jac)  # please remain 200 here - it is only a technical thing
    }
    subsamples <- subsamples[(jacc.init != 0) & (jacc.init != 1)]
    
    if(length(subsamples) <= 2) next
    
    # Calculate jaccard
    if (is.null(p.val.cutoff)) {
      jacc.tmp <- jaccardPwTop(subsamples, top.n.genes)  
    } else {
      jacc.tmp <- jaccardPwPval(subsamples, p.val.cutoff)  
    }
    
    if(is.null(jacc.tmp$jac)) next
    # print(jacc.tmp)
    data.tmp <- data.frame(group = cell.type,
                           value = jacc.tmp$jac,
                           cmp = jacc.tmp$id)
    data.all <- rbind(data.all, data.tmp)
  }
  
  return(data.all)
}


plotStability <- function(jaccards,
                         notch,
                         show.jitter,
                         jitter.alpha,
                         show.pairs,
                         sort.order,
                         xlabel = '',
                         ylabel = '',
                         log.y.axis = F,
                         palette = NULL,
                         plot.theme=theme_get()) {
  jaccards$group <- as.factor(jaccards$group)
  if(! show.pairs) {
    if(sort.order) {
      p <- ggplot(jaccards, aes(x=reorder(group,value,na.rm = TRUE), y=value,
                                group=group)) + geom_boxplot(outlier.shape = NA, notch = notch) + geom_jitter(alpha=jitter.alpha, aes(color=cmp))         
    } else {
      p <- ggplot(jaccards, aes(x=group, y=value,
                                group=group)) + geom_boxplot(outlier.shape = NA, notch = notch) + geom_jitter(alpha=jitter.alpha,aes(color=cmp))          
    }
  } else {
    if(sort.order) {
      p <- ggplot(jaccards, aes(x=reorder(group,value,na.rm = TRUE), y=value,
                                group=cmp, color=cmp)) + geom_line()  
    } else {
      p <- ggplot(jaccards, aes(x=group, y=value,
                                group=cmp, color=cmp)) + geom_line()  
    }
  }
  
  p <- p + plot.theme + theme(legend.position = "none") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
    labs(x=xlabel, y=ylabel)
  
  if(log.y.axis) {
    p <- p + scale_y_continuous(trans='log10')
  }
  
  if(!is.null(palette)) {
    p <- p + scale_color_manual(values=palette)
  }
  return(p)
  
}

estimateNumberOfTermsStability <- function(de.res,
                                           p.adjust,
                                           pvalue.cutoff){
  data.all <- data.frame()
  for(cell.type in names(de.res)){
    # print(cell.type)
    subsamples <- de.res[[cell.type]]$subsamples
    
    # Remove some subsamples due to the min.cell.counts
    # coomare resampling results with "initial"
    jacc.init <- c()
    for(subs.name in names(subsamples)){
      subsamples.tmp = list(de.res[[cell.type]]$subsamples[[subs.name]],
                            de.res[[cell.type]]$res)
      jacc.init <- c(jacc.init, jaccardPwTop(subsamples.tmp, 200))
    }
    subsamples <- subsamples[(jacc.init != 0) & (jacc.init != 1)]
    
    if(length(subsamples) <= 2) next
    
    # Calculate numbers of significant terms
    if (p.adjust) {
      num.tmp <- sapply(names(subsamples), function(x) sum(subsamples[[x]]$padj <= pvalue.cutoff))
    } else {
      num.tmp <- sapply(names(subsamples), function(x) sum(subsamples[[x]]$pvalue <= pvalue.cutoff))
    }
    
    data.tmp <- data.frame(group = cell.type,
                           value = num.tmp,
                           cmp = names(num.tmp))
    data.all <- rbind(data.all, data.tmp)
  }
  return(data.all)
}


# extractGOtables = function(go.res, 
#                            ont.types = NULL,
#                            dir.type = 'all'){
#   cell.types <- names(go.res)
#   if(is.null(ont.types))
#     ont.types <- c('BP', 'CC', 'MF')
#   if( !(dir.type %in% c('all', 'up', 'dows')) ) 
#     stop('Provided direction of gene changing in not correct')
#   
#   go.cell.types.list <- list()
#   for(cell.type in cell.types){
#     go.cell.type.table <- data.frame()
#     for(ont.type in ont.types){
#       table.tmp <- go.res[[cell.type]][[ont.type]][[dir.type]]@result[,c('ID', 'pvalue', 'p.adjust')]
#       colnames(table.tmp) <- c('ID', 'pvalue', 'padj')
#       table.tmp$ont.type <- ont.type
#       go.cell.type.table <- rbind(go.cell.type.table, table.tmp)
#     }
#     go.cell.types.list[[cell.type]] <- go.cell.type.table
#   }
#   
#   return(go.cell.types.list)
#   
# }


estimateDEStabilityFDR = function(de.res, p.adj.cutoffs) {
    
  df.n.common.genes <- matrix(ncol=length(de.res),
                              nrow=length(p.adj.cutoffs), 
                              dimnames=list(NULL, names(de.res)))
  
  df.n.all.genes <- matrix(ncol=length(de.res),
                            nrow=length(p.adj.cutoffs), 
                            dimnames=list(NULL, names(de.res)))
  
  df.n.init.genes <- matrix(ncol=length(de.res),
                            nrow=length(p.adj.cutoffs), 
                            dimnames=list(NULL, names(de.res)))
  
  for(i in 1:length(p.adj.cutoffs)) {
    p <- p.adj.cutoffs[i]
    for(cell.type in names(de.res)){
      genes <- rownames(de.res[[cell.type]]$res)
      genes.all <- c()
      if(is.null(genes)) stop(cell.type)
      for(ires in names(de.res[[cell.type]]$subsamples)) {
        res.tmp <- de.res[[cell.type]]$subsamples[[ires]]
        genes.tmp <- rownames(res.tmp)[res.tmp[,'padj'] < p]
        genes <- intersect(genes, genes.tmp)
        genes.all <- unique(c(genes, genes.tmp))
  
      }
      df.n.common.genes[i, cell.type] <- length(genes)
      df.n.all.genes[i, cell.type] <- length(genes.all)
      
      res.tmp <- de.res[[cell.type]]$res
      df.n.init.genes[i, cell.type] <- length(rownames(res.tmp)[res.tmp[,'padj'] < p])
    }
  }
  
  df.n.init.genes <- melt(df.n.init.genes)
  df.n.init.genes$type <- 'init'
  df.n.all.genes <- melt(df.n.all.genes)
  df.n.all.genes$type <- 'all'
  df.n.common.genes <- melt(df.n.common.genes)
  df.n.common.genes$type <- 'common'
  
  df.n.genes <- rbind(df.n.common.genes, df.n.all.genes, df.n.init.genes)
  df.n.genes$Var1 <- as.factor(p.adj.cutoffs[df.n.genes$Var1 ])
  
  return(df.n.genes)
}

estimateDEStabilityFDR1loo = function(de.res, p.adj.cutoffs) {
  
  df.n.common.genes <- c()
  
  for(ires.target in 1:length(de.res[[1]]$subsamples)){
    for(i in 1:length(p.adj.cutoffs)) {
      p <- p.adj.cutoffs[i]
      for(cell.type in names(de.res)){
        
        res.tmp <- de.res[[cell.type]]$subsamples[[ires.target]]
        genes <- rownames(res.tmp)[res.tmp[,'padj'] < p]
        
        if(length(genes) == 0) next
        if(is.null(genes)) stop(cell.type)
        
        genes.all <- c()
        for(ires in names(de.res[[cell.type]]$subsamples)) {
          res.tmp <- de.res[[cell.type]]$subsamples[[ires]]
          genes.tmp <- rownames(res.tmp)[res.tmp[,'padj'] < p]
          genes <- intersect(genes, genes.tmp)
          genes.all <- unique(c(genes, genes.tmp))
          
        }
        res.tmp <- de.res[[cell.type]]$subsamples[[ires.target]]
        gene.frac <- length(genes) / length(rownames(res.tmp)[res.tmp[,'padj'] < p])
        df.n.common.genes <- rbind(df.n.common.genes, c(cell.type, p, gene.frac))
      }
    }
  }
  
  df.n.common.genes = as.data.frame(df.n.common.genes)
  colnames(df.n.common.genes) <- c('cell.type', 'fdr', 'frac')
  df.n.common.genes$frac = as.numeric(df.n.common.genes$frac)
  df.n.common.genes$fdr = as.numeric(df.n.common.genes$fdr)

  return(df.n.common.genes)
}
  
  

estimateStabilityPerCellType <- function(de.res,
                                         top.n.genes,
                                         p.val.cutoff) {
  
  data.all <- data.frame()
  for(cell.type in names(de.res)){
    # print(cell.type)
    subsamples <- de.res[[cell.type]]$subsamples
    
    # Remove some subsamples due to the min.cell.counts
    # coomare resampling results with "initial"
    jacc.init <- c()
    for(subs.name in names(subsamples)){
      subsamples.tmp = list(de.res[[cell.type]]$subsamples[[subs.name]],
                            de.res[[cell.type]]$res)
      jacc.pw.tmp <- jaccardPwTop(subsamples.tmp, 200)
      jacc.init <- c(jacc.init, jacc.pw.tmp$jac)  # please remain 200 here - it is only a technical thing
    }
    subsamples <- subsamples[(jacc.init != 0) & (jacc.init != 1)]
    
    if(length(subsamples) <= 2) next
    
    # Calculate jaccard
    if (is.null(p.val.cutoff)) {
      jacc.tmp <- jaccardPwTop(subsamples, top.n.genes)  
    } else {
      jacc.tmp <- jaccardPwPval(subsamples, p.val.cutoff)  
    }
    
    if(is.null(jacc.tmp$jac)) next
    # print(jacc.tmp)
    data.tmp <- data.frame(group = cell.type,
                           value = jacc.tmp$jac,
                           cmp = jacc.tmp$id)
    data.all <- rbind(data.all, data.tmp)
  }
  
  return(data.all)
}



  
  
  
