plotPcaSpace <- function(d.counts, d.groups, ref.level, target.level, font.size, palette=NULL){
  bal <- getRndBalances(d.counts)
  pca.res <- prcomp(bal$norm)
  pca.loadings <- bal$psi %*% pca.res$rotation

  df.pca <- as.data.frame(pca.res$x)

  pc1 <- pca.loadings[,1]
  pc2 <- pca.loadings[,2]
  df.loadings <- as.data.frame(cbind(pc1, pc2) * 10)


  # ----------- PLOT -----------
  group.names <- c(ref.level, target.level)
  options(repr.plot.width = 15, repr.plot.height = 10)
  rda.plot <- ggplot(df.pca, aes(x=PC1, y=PC2)) +
    #   geom_text(aes(label=rownames(df_pca) %in% samplegroups$trgt),size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(group.names[d.groups + 1] ))) +
    labs(colour="Group") +
    coord_fixed()

  if(!is.null(palette)) rda.plot <- rda.plot + scale_fill_manual(values=palette)

  dx <- max(df.pca[,'PC1']) - min(df.pca[,'PC1'])
  dy <- max(df.pca[,'PC2']) - min(df.pca[,'PC2'])


  rda.biplot <- rda.plot +
    geom_segment(data=df.loadings, aes(x=0, xend=pc1, y=0, yend=pc2),
                 color="grey", arrow=arrow(length=unit(0.01,"npc")))  +
    coord_flip(clip = "off") +
    geom_text(data=df.loadings,
              aes(x=pc1,y=pc2,label=rownames(df.loadings)),
              color="black", size=3)



  dx <- max(dx, max(df.loadings[,'pc1']) - min(df.loadings[,'pc1']))
  dy <- max(dy, max(df.loadings[,'pc2']) - min(df.loadings[,'pc2']))

  rda.biplot <- rda.biplot + coord_fixed(ratio = dx / dy)
  if(!is.null(font.size)) {
    rda.biplot <- rda.biplot + theme(axis.text=element_text(size=font.size), axis.title=element_text(size=font.size))
  }
  return(rda.biplot)
}

plotCdaSpace <- function(d.counts, d.groups, ref.level, target.level, font.size, thresh.pc.var = 0.95, n.dim = 2){

  cell.loadings <- c()
  sample.pos <- c()

  bal <- getRndBalances(d.counts)
  d.used <- bal$norm

  if(ncol(d.used) != 2){
    for(i in 1:n.dim){
      res.remove <- removeGroupEffect(d.used, d.groups, thresh.pc.var = 0.9)
      cell.loadings <- cbind(cell.loadings, bal$psi %*% res.remove$rotation)
      sample.pos <- cbind(sample.pos, res.remove$scores)
      # d.used <- res.remove$remain
      d.used <- d.used - res.remove$used.part
    }
  }else{
    cell.loadings <- bal$psi
    sample.pos <- d.used
  }

  colnames(cell.loadings) <- paste('C', 1:n.dim, sep = '')
  colnames(sample.pos) <- paste('Score', 1:n.dim, sep = '')

  df.pca <- as.data.frame(sample.pos)
  df.loadings <- as.data.frame(cell.loadings * 8)

  # ----------- PLOT -----------
  group.names <- c(ref.level, target.level)
  options(repr.plot.width = 15, repr.plot.height = 10)
  rda.plot <- ggplot(df.pca, aes(x=Score1, y=Score2)) +
    #   geom_text(aes(label=rownames(df_pca) %in% samplegroups$trgt),size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(group.names[d.groups + 1] ))) +
    labs(colour="Group") +
    coord_fixed()

  dx <- max(df.pca[,'Score1']) - min(df.pca[,'Score1'])
  dy <- max(df.pca[,'Score2']) - min(df.pca[,'Score2'])


  rda.biplot <- rda.plot +
    geom_segment(data=df.loadings, aes(x=0, xend=C1, y=0, yend=C2),
                 color="grey", arrow=arrow(length=unit(0.01,"npc")))  +
    geom_text(data=df.loadings,
              aes(x=C1,y=C2,label=rownames(df.loadings)),
              color="black", size=3)


  dx <- max(dx, max(df.loadings[,'C1']) - min(df.loadings[,'C1']))
  dy <- max(dy, max(df.loadings[,'C2']) - min(df.loadings[,'C2']))

  rda.biplot <- rda.biplot + coord_fixed(ratio = dx / dy)
  
  if(!is.null(font.size)) {
    rda.biplot <- rda.biplot + theme(axis.text=element_text(size=font.size), axis.title=element_text(size=font.size))
  }

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

plotContrastTree <- function(d.counts, d.groups, ref.level, target.level, p.threshold = 0.01){
  log.f <- getLogFreq(d.counts)

  t.cur <- constructCanonicalTree(d.counts, d.groups)
  t.tmp <- compute.brlen(t.cur, method="Grafen")
  d.cur <- as.dendrogram(t.tmp)
  t.tmp <- as.phylo(d.cur)

  # ---------------------------------
  # Positions

  dend.data <- ggdendro::dendro_data(d.cur, type = "rectangle")

  # seg
  node.pos <- dend.data$segments %$% .[(y == yend) & (yend != 0),]
  node.pos$id <- t.tmp$edge[,1]

  innode.pos <- unique(node.pos[,c('x','y','id')])

  innode.pos$range <- -1
  for(i in 1:nrow(innode.pos)){
    tmp <- node.pos$xend[node.pos$id == innode.pos$id[i]]
    innode.pos$range[i] <- max(tmp) - min(tmp)
  }
  rownames(innode.pos) <- innode.pos$id

  # ----------------------------------------

  # Balances
  res <- getNodeBalances(t.tmp, log.f)
  balances <- as.data.frame(res$bal)
  colnames(balances) <- rownames(innode.pos)

  p.val <- c()
  for(i in 1:ncol(balances)){
    aov.data <- cbind(balances[,i], d.groups)
    colnames(aov.data) <- c('balance', 'group')

    res <- aov(balance~group,data=as.data.frame(aov.data))
    p.val <- c(p.val,summary(res)[[1]][1,5])
  }
  # p.val <- p.adjust(p.val, method = 'fdr')
  p.val[is.na(p.val)] <- 1

  px <- ggdend(dend.data)

  if(sum(p.val < p.threshold) == 0)
    return(px)

  bx <- c()
  by <- c()
  bc <- c()

  bx.m <- c()
  by.m <- c()
  bc.m <- c()

  pval.x <- c()
  pval.y <- c()

  n.range <- c()
  x.range <- c()
  y.range <- c()

  idx.target <- d.groups
  for(id.node in 1:ncol(balances)){

    if(p.val[id.node] < p.threshold){
      pval.x <- c(pval.x, innode.pos$x[id.node])
      pval.y <- c(pval.y, innode.pos$y[id.node])
    }else
      next

    x.tmp <- -balances[,id.node]
    x.tmp <- x.tmp - mean(x.tmp)
    n.range <- c(n.range, max(abs(x.tmp)))  # Ticks

    x.tmp <- x.tmp / max(abs(x.tmp)) / 2 * innode.pos$range[id.node] * 0.9
    x.tmp <- x.tmp + innode.pos$x[id.node]
    y.tmp <- idx.target * 0.03 + innode.pos$y[id.node] - 0.05

    x.range <- c(x.range, innode.pos$x[id.node] + innode.pos$range[id.node]/2)  # Ticks
    y.range <- c(y.range, innode.pos$y[id.node])  # Ticks

    c.tmp <- idx.target

    bx <- c(bx, x.tmp)
    by <- c(by, y.tmp)
    bc <- c(bc, c.tmp)

    bx.m <- c(bx.m, mean(x.tmp[c.tmp]), mean(x.tmp[!c.tmp]))
    by.m <- c(by.m, mean(y.tmp[c.tmp]), mean(y.tmp[!c.tmp]))
    bc.m <- c(bc.m, TRUE, FALSE)
  }

  group.names <- c(ref.level, target.level)
  df.balance.points <- data.frame(x = bx, y = by, group = group.names[bc+1])
  df.balance.mean <- data.frame(x = bx.m, y = by.m, group = group.names[bc.m+1])
  df.balance.range <- data.frame(x = x.range, y = y.range, s = n.range)
  df.pval <- data.frame(x = pval.x, y = pval.y)


  px <- px + geom_point(data = df.balance.points,
                       aes(x=x, y=y, col = as.factor(group)), alpha = 0.5, size = 1) +
    geom_point(data = df.balance.mean,
               aes(x=x, y=y, col = as.factor(group)),
               size = 3, shape = 18, stroke = 2) +
    labs(col=" ") +
    geom_text(data=df.balance.range, mapping=aes(x=x, y=y, label=sprintf('%2.1f',s)), vjust=0, hjust=0)

  if(length(pval.x) > 0){
    px <- px + geom_point(data = df.pval, aes(x=x, y=y)) +
      # labs(size= sprintf('p-value < %1.2f',p.threshold))
      guides(size = FALSE)
  }


  return(px)
}


plotCellLoadings <- function(cda, ordering, signif.threshold, font.size, alpha, palette, show.pvals, ref.level, target.level) {
  balances = cda$balances
  if(ordering == 'by.pvalue'){
    # ordering by median
    balances <- balances[order(abs(apply(balances, 1, median))), ]
    # additional ordering by p-value
    frac <- getCellSignificance(balances)
    balances <- balances[order(-frac),]
    
    # Get significant cells
    n.significant.cells = sum(frac <= signif.threshold)
    
  } else if(ordering == 'by.mean') {
    balances <- balances[order(abs(rowMeans(balances))),]
    n.significant.cells = 0
  } else if(ordering == 'by.median') {
    balances <- balances[order(abs(apply(balances, 1, median))), ]
    n.significant.cells = 0
  }
  
  frac <- getCellSignificance(balances)
  res.ordered <- t(balances) %>% as.data.frame()  
  ymax = max(balances)

  p <- ggplot(stack(res.ordered), aes(x = ind, y = values, fill=factor(ind))) +
    geom_boxplot(notch=TRUE, outlier.shape = NA) + geom_jitter(aes(x = ind, y = values), alpha = alpha, size=1) +
    geom_hline(yintercept = 0, color = "gray37") +
    coord_flip() + xlab('') + ylab('') + theme_bw()+ theme(legend.position = "none") +
    scale_x_discrete(position = "top") + ylim(-ymax, ymax)
  
  # Add text
  p <- p + 
    annotate('text', x = 1, y = -ymax, label = paste('\u2190', ref.level), hjust = 'left') + 
    annotate('text', x = 1, y = ymax, label = paste(target.level, '\u2192'), hjust = 'right')
  
  if(!is.null(font.size)) {
    p <- p + theme(axis.text=element_text(size=font.size), axis.title=element_text(size=font.size))
  }
  if(!is.null(palette)) p <- p + scale_fill_manual(values=palette)
  if(n.significant.cells > 0) p <- p + geom_vline(xintercept=nrow(balances) - n.significant.cells + 0.5, color='red')
  
  
  if(show.pvals){
    d = data.frame(x = names(frac), y=frac)
    d$x <- factor(d$x, levels = d$x)
    p.pval <- ggplot(d, aes(x=x,y=-log(y,base = 10) )) + geom_bar(stat="identity") +
      coord_flip() + xlab('') + ylab('-log(p-value)') + theme_bw()+ theme(legend.position = "none") +
      geom_hline(yintercept=-log(signif.threshold,base = 10)) + theme(axis.text.y = element_blank())
    
    p.combo <- ggarrange(p, p.pval, ncol = 2, widths = c(2, 1))  
    return(p.combo)
  }
  
  return(p)
  
}
