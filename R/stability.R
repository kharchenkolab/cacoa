jaccard.pw.top <- function(subsamples, top.thresh){
  jac.all = c()
  for(i in 1:length(subsamples)) {
    for(j in 1:length(subsamples)) {
      if (j <= i) next
      set1 <- rownames(subsamples[[i]])[rank(subsamples[[i]]$pvalue) <= top.thresh]
      set2 <- rownames(subsamples[[j]])[rank(subsamples[[j]]$pvalue) <= top.thresh]
      jac.all <- c(jac.all, length(intersect(set1, set2)) / length(unique(c(set1, set2))))
    }
  }
  return(jac.all)
}

jaccard.pw.pval <- function(subsamples, pval.thresh){
  jac.all = c()
  for(i in 1:length(subsamples)) {
    for(j in 1:length(subsamples)) {
      if (j <= i) next
      set1 <- rownames(subsamples[[i]])[subsamples[[i]]$padj <= pval.thresh]
      set2 <- rownames(subsamples[[j]])[subsamples[[j]]$padj <= pval.thresh]
      jac.all <- c(jac.all, length(intersect(set1, set2)) / length(unique(c(set1, set2))))
    }
  }
  return(jac.all)
}

indexes.of.pairs <- function(subsamples, samples){
  idxs = c()
  id = 1
  for(i in 1:length(samples)) {
    for(j in 1:length(samples)) {
      if (j <= i) next
      id <- id + 1
      if (! (samples[i] %in% subsamples)) next
      if (! (samples[j] %in% subsamples)) next
      idxs <- c(idxs, id)
    }
  }
  return(idxs)
}

estimateStabilityPerCellType <- function(de.res,
                                         top.n.genes,
                                         padj.threshold) {

  data.all = data.frame()
  for(cell.type in names(de.res)){
    # print(cell.type)
    subsamples <- de.res[[cell.type]]$subsamples
    
    # Remove some subsamples due to the min.cell.counts
    # coomare resampling results with "initial"
    jacc.init = c()
    for(subs.name in names(subsamples)){
      subsamples.tmp = list(de.res[[cell.type]]$subsamples[[subs.name]],
                            de.res[[cell.type]]$res)
      jacc.init = c(jacc.init, jaccard.pw.top(subsamples.tmp, 200))
    }
    subsamples <- subsamples[(jacc.init != 0) & (jacc.init != 1)]
    
    if(length(subsamples) <= 2) next
    
    # Calculate jaccard
    if (is.null(padj.threshold)) {
      jacc.tmp <- jaccard.pw.top(subsamples, top.n.genes)  
    } else {
      jacc.tmp <- jaccard.pw.pval(subsamples, padj.threshold)  
    }
    
    
    data.tmp <- data.frame(group = cell.type,
                           value = jacc.tmp,
                           cmp = indexes.of.pairs(names(subsamples), names(de.res[[cell.type]]$subsamples)))
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
                         palette = NULL) {
  if(! show.pairs) {
    if(sort.order) {
      p <- ggplot(jaccards, aes(x=reorder(group,value,na.rm = TRUE), y=value, fill=group,
                                group=group)) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha=jitter.alpha)          
    } else {
      p <- ggplot(jaccards, aes(x=group, y=value, fill=group,
                                group=group)) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha=jitter.alpha)          
    }
  } else {
    if(sort.order) {
      p <- ggplot(jaccards, aes(x=reorder(group,value,na.rm = TRUE), y=value, fill=group,
                                group=cmp, color=cmp)) + geom_line()  
    } else {
      p <- ggplot(jaccards, aes(x=group, y=value, fill=group,
                                group=cmp, color=cmp)) + geom_line()  
    }
  }
  
  p <- p + theme(legend.position = "none") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
    labs(x=xlabel, y=ylabel)
  
  if(log.y.axis) {
    p <- p + scale_y_continuous(trans='log10')
  }
  
  # if(!is.null(palette)) p <- p + scale_fill_manual(values=palette)
  
  return(p)
  
}

estimateNumberOfTermsStability <- function(de.res,
                                           p.adjust,
                                           pvalue.cutoff){
  data.all = data.frame()
  for(cell.type in names(de.res)){
    # print(cell.type)
    subsamples <- de.res[[cell.type]]$subsamples
    
    # Remove some subsamples due to the min.cell.counts
    # coomare resampling results with "initial"
    jacc.init = c()
    for(subs.name in names(subsamples)){
      subsamples.tmp = list(de.res[[cell.type]]$subsamples[[subs.name]],
                            de.res[[cell.type]]$res)
      jacc.init = c(jacc.init, jaccard.pw.top(subsamples.tmp, 200))
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


extractGOtables = function(go.res, 
                           ont.types = NULL,
                           dir.type = 'all'){
  cell.types = names(go.res)
  if(is.null(ont.types))
    ont.types = c('BP', 'CC', 'MF')
  if( !(dir.type %in% c('all', 'up', 'dows')) ) 
    stop('Provided direction of gene changing in not correct')
  
  go.cell.types.list = list()
  for(cell.type in cell.types){
    go.cell.type.table = data.frame()
    for(ont.type in ont.types){
      table.tmp = go.res[[cell.type]][[ont.type]][[dir.type]]@result[,c('ID', 'pvalue', 'p.adjust')]
      colnames(table.tmp) <- c('ID', 'pvalue', 'padj')
      table.tmp$ont.type = ont.type
      go.cell.type.table = rbind(go.cell.type.table, table.tmp)
    }
    go.cell.types.list[[cell.type]] <- go.cell.type.table
  }
  
  return(go.cell.types.list)
  
}



  
  
  
  
  
