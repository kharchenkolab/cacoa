plotCellLoadings <- function(cda.resamples, aplha = 0.01){
  res.init = t(cda.resamples)
  
  # Define order of cell types
  res.ordered= res.init[,order(abs(colMeans(res.init)))]
  # Dataframe
  dat <- stack(as.data.frame(res.ordered))
  
  font.size = 16
  
  p = ggplot(dat, aes(x = ind, y = values, fill=factor(ind))) + 
    geom_boxplot(notch=TRUE, outlier.shape = NA) + 
    geom_jitter(aes(x = ind, y = values), alpha = aplha, size=1) +
    # scale_fill_manual(values=cell_colors)+
    geom_hline(yintercept = 0, color = "gray37") +
    coord_flip() + 
    xlab('') + ylab('') +
    theme_bw()  +
    theme(axis.text=element_text(size=font.size),
          axis.title=element_text(size=font.size)) + theme(legend.position = "none")  
  return(p)
}

plotPcaSpace <- function(d.counts, d.groups){
  bal = getRndBalances(d.counts)
  pca.res = prcomp(bal$norm)
  pca.loadings = bal$psi %*% pca.res$rotation
  
  df.pca = as.data.frame(pca.res$x)
  
  pc1 = pca.loadings[,1]
  pc2 = pca.loadings[,2]
  df.loadings = as.data.frame(cbind(pc1, pc2) * 10)
  
  
  # ----------- PLOT -----------
  options(repr.plot.width = 15, repr.plot.height = 10)
  rda.plot <- ggplot(df.pca, aes(x=PC1, y=PC2)) +
    #   geom_text(aes(label=rownames(df_pca) %in% samplegroups$trgt),size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(d.groups))) +
    coord_fixed()
  
  
  rda.biplot <- rda.plot +
    geom_segment(data=df.loadings, aes(x=0, xend=pc1, y=0, yend=pc2),
                 color="grey", arrow=arrow(length=unit(0.01,"npc")))  +
    coord_flip(clip = "off") +
    geom_text(data=df.loadings,
              aes(x=pc1,y=pc2,label=rownames(df.loadings)),
              color="black", size=4)
  return(rda.biplot)
}

plotCdaSpace <- function(d.counts, d.groups, thresh.pc.var = 0.95, n.dim = 2){

  cell.loadings = c()
  sample.pos = c()
  
  bal = getRndBalances(d.counts)
  d.used = bal$norm
  
  for(i in 1:n.dim){
    res.remove = removeGroupEffect(d.used, d.groups, thresh.pc.var = 0.9)
    cell.loadings = cbind(cell.loadings, bal$psi %*% res.remove$rotation)
    sample.pos = cbind(sample.pos, res.remove$scores)
    # d.used = res.remove$remain
    d.used = d.used - res.remove$used.part
  }
  
  colnames(cell.loadings) <- paste('C', 1:n.dim, sep = '')
  colnames(sample.pos) <- paste('Score', 1:n.dim, sep = '')
  
  cell.loadings
  sample.pos
  
  
  df.pca = as.data.frame(sample.pos)
  df.loadings = as.data.frame(cell.loadings * 10)
  
  
  # ----------- PLOT -----------
  options(repr.plot.width = 15, repr.plot.height = 10)
  rda.plot <- ggplot(df.pca, aes(x=Score1, y=Score2)) +
    #   geom_text(aes(label=rownames(df_pca) %in% samplegroups$trgt),size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(d.groups))) +
    coord_fixed()
  
  
  rda.biplot <- rda.plot +
    geom_segment(data=df.loadings, aes(x=0, xend=C1, y=0, yend=C2),
                 color="grey", arrow=arrow(length=unit(0.01,"npc")))  +
    geom_text(data=df.loadings,
              aes(x=C1,y=C2,label=rownames(df.loadings)),
              color="black", size=4)
  
  return(rda.biplot)
}

# helper function for creating dendograms
ggdend <- function(dend.data, a = 90) {
  ggplot() +
    geom_segment(data = dend.data$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "")  + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    geom_text(data = dend.data$labels, aes(x, y, label = label),
              hjust = 1, angle = a, size = 3) + ylim(-0.5, NA) 
}

plotContrastTree <- function(d.counts, d.groups, p.threshold = 0.01){
  
  log.f = getLogFreq(d.counts)

  t.cur = constructCanonicalTree(d.counts, d.groups)
  t.tmp <- compute.brlen(t.cur, method="Grafen")
  d.cur <- as.dendrogram(t.tmp)
  t.tmp = as.phylo(d.cur)
  
  # ---------------------------------
  # Positions
  
  dend.data <- dendro_data(d.cur, type = "rectangle")
  # What contains dend.data
  names(dend.data)
  
  
  seg = dend.data$segments
  # seg
  node.pos = seg[(seg$y == seg$yend) & (seg$yend != 0),]
  node.pos$id = t.tmp$edge[,1]
  
  innode.pos = unique(node.pos[,c('x','y','id')])
  
  innode.pos$range = -1
  for(i in 1:nrow(innode.pos)){
    tmp = node.pos$xend[node.pos$id == innode.pos$id[i]]
    innode.pos$range[i] = max(tmp) - min(tmp)
  }
  rownames(innode.pos) = innode.pos$id
  
  # ----------------------------------------
  
  # Balances
  res <- getNodeBalances(t.tmp, log.f)
  balances <- as.data.frame(res$bal)
  colnames(balances) = rownames(innode.pos)
  
  p.val = c()
  for(i in 1:ncol(balances)){
    aov.data <- cbind(balances[,i], d.groups)
    colnames(aov.data) <- c('balance', 'group')
    
    res = aov(balance~group,data=as.data.frame(aov.data))
    p.val = c(p.val,summary(res)[[1]][1,5])
  }
  p.val = p.adjust(p.val, method = 'fdr')
  p.val[is.na(p.val)] = 1
  
  px <- ggdend(dend.data) 
  
  if(sum(p.val < p.threshold) == 0)
    return(px)
  
  bx = c()
  by = c()
  bc = c()
  
  bx.m = c()
  by.m = c()
  bc.m = c()
  
  pval.x = c()
  pval.y = c()
  
  n.range = c()
  x.range = c()
  y.range = c()
  
  idx.target = d.groups
  for(id.node in 1:ncol(balances)){
    
    if(p.val[id.node] < p.threshold){
      pval.x = c(pval.x, innode.pos$x[id.node])
      pval.y = c(pval.y, innode.pos$y[id.node])
    }else
      next
    
    x.tmp = -balances[,id.node] 
    x.tmp = x.tmp - mean(x.tmp)
    n.range = c(n.range, max(abs(x.tmp)))  # Ticks
    
    x.tmp = x.tmp / max(abs(x.tmp)) / 2 * innode.pos$range[id.node] * 0.9
    x.tmp = x.tmp + innode.pos$x[id.node]
    y.tmp = idx.target * 0.03 + innode.pos$y[id.node] - 0.05
    
    x.range = c(x.range, innode.pos$x[id.node] + innode.pos$range[id.node]/2)  # Ticks
    y.range = c(y.range, innode.pos$y[id.node])  # Ticks
    
    c.tmp = idx.target
    
    bx = c(bx, x.tmp)
    by = c(by, y.tmp)
    bc = c(bc, c.tmp) 
    
    bx.m = c(bx.m, mean(x.tmp[c.tmp]), mean(x.tmp[!c.tmp]))
    by.m = c(by.m, mean(y.tmp[c.tmp]), mean(y.tmp[!c.tmp]))
    bc.m = c(bc.m, TRUE, FALSE)
  }
  
  group.names = c('Control', 'Case')
  df.balance.points = data.frame(x = bx, y = by, group = group.names[bc+1])
  df.balance.mean = data.frame(x = bx.m, y = by.m, group = group.names[bc.m+1])
  df.balance.range = data.frame(x = x.range, y = y.range, s = n.range)
  df.pval = data.frame(x = pval.x, y = pval.y)

  px = px + geom_point(data = df.balance.points, 
                       aes(x=x, y=y, col = as.factor(group)), alpha = 0.5, size = 1) +
    geom_point(data = df.balance.mean, 
               aes(x=x, y=y, col = as.factor(group)), 
               size = 3, shape = 18, stroke = 2) + 
    labs(col="Group") + 
    geom_text(data=df.balance.range, mapping=aes(x=x, y=y, label=sprintf('%2.1f',s)), vjust=0, hjust=0)
  
  # if(length(pval.x) > 0){
  #   px <- px + geom_point(data = df.pval, aes(x=x, y=y, size = 2)) +
  #     labs(size= sprintf('p-value < %1.1e',p.threshold)) 
  # }
   
  options(repr.plot.width = 10, repr.plot.height = 10)
  
  return(px)
}



