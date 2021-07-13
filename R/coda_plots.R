estimateCdaSpace <- function(d.counts, d.groups, thresh.pc.var = 0.95, n.dim = 2){
  cell.loadings <- c()
  sample.pos <- c()

  bal <- getRndBalances(d.counts)
  d.used <- bal$norm

  if(ncol(d.used) != 2){
    for(i in 1:n.dim){
      res.remove <- removeGroupEffect(d.used, d.groups, thresh.pc.var = 0.9)
      cell.loadings <- cbind(cell.loadings, bal$psi %*% res.remove$rotation)
      sample.pos <- cbind(sample.pos, res.remove$scores)
      d.used <- d.used - res.remove$used.part
    }
  } else {
    cell.loadings <- bal$psi
    sample.pos <- d.used
  }

  cn <- paste('S', 1:n.dim, sep = '')
  df.cda <- as.data.frame(sample.pos) %>% set_colnames(cn)
  df.loadings <- as.data.frame(cell.loadings * 8) %>% set_colnames(cn)
  return(list(red=df.cda, loadings=df.loadings))
}

plotCodaSpaceInner <- function(df.space, df.loadings, d.groups, ref.level, target.level, font.size=3, palette=NULL) {
  group.names <- c(ref.level, target.level)
  rda.plot <- ggplot(df.space, aes(x=S1, y=S2)) +
    #   geom_text(aes(label=rownames(df_pca) %in% samplegroups$trgt),size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(group.names[d.groups + 1] ))) +
    labs(colour="Condition")

  if(!is.null(palette)) rda.plot <- rda.plot + scale_color_manual(values=palette)

  rda.biplot <- rda.plot +
    geom_segment(data=df.loadings, aes(x=0, xend=S1, y=0, yend=S2),
                 color="grey", arrow=arrow(length=unit(0.01,"npc")))  +
    geom_text(data=df.loadings,
              aes(x=S1, y=S2, label=rownames(df.loadings)),
              color="black", size=font.size)

  dx <- max(diff(range(df.space$S1)), diff(range(df.loadings$S1)))
  dy <- max(diff(range(df.space$S2)), diff(range(df.loadings$S2)))

  rda.biplot <- rda.biplot + coord_fixed(ratio=dy/dx)

  return(rda.biplot)
}

# helper function for creating dendograms
ggdend <- function(dend.data, a = 90, plot.theme=theme_get()) {
  ggplot() +
    geom_segment(data = dend.data$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + plot.theme +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), panel.border=element_blank(),
          axis.line=element_blank()) +
    geom_text(data = dend.data$labels, aes(x, y, label = label),
              hjust = 1, angle = a, size = 3) + ylim(-0.5, NA)
}

plotContrastTree <- function(d.counts, d.groups, ref.level, target.level, plot.theme, p.threshold = 0.01){
  log.f <- getLogFreq(d.counts)

  t.cur <- constructCanonicalTree(d.counts, d.groups)
  t.tmp <- compute.brlen(t.cur, method="Grafen") %>% as.hclust()
  d.cur <- as.dendrogram(t.tmp)
  t.tmp <- as.phylo(t.tmp)

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

  px <- ggdend(dend.data, plot.theme=plot.theme)

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


plotCellLoadings <- function(loadings, pval, signif.threshold, jitter.alpha, palette, 
                             show.pvals, ref.level, target.level, plot.theme) {
  
  yintercept = 0

  if(!is.null(pval)){
    # if some p-values are the same - then order by mean, therefore a prior sorting is required
    tmp.order <- order(abs(rowMeans(loadings)))
    loadings <- loadings[tmp.order, ]
    pval <- pval[tmp.order]
    # additional ordering by p-value
    loadings <- loadings[order(-pval), ]
    pval <- pval[order(-pval)]
    # Get significant cells
    n.significant.cells <- sum(pval <= signif.threshold)
  } else {
    # ordering by mean
    loadings <- loadings[order(abs(rowMeans(loadings))), ]
    n.significant.cells <- 0
  }

  res.ordered <- t(loadings) %>% as.data.frame()
  ymax = max(loadings)

  p <- ggplot(stack(res.ordered), aes(x = ind, y = values, fill=factor(ind))) +
    geom_boxplot(notch=TRUE, outlier.shape = NA) + geom_jitter(aes(x = ind, y = values), alpha = jitter.alpha, size=1) +
    geom_hline(yintercept = yintercept, color = "gray37") +
    coord_flip() + xlab('') + ylab('loadings') + plot.theme + theme(legend.position = "none") +
    scale_x_discrete(position = "top") + ylim(-1, 1)

  # Add text
  p <- p +
    annotate('text', x = 1, y = -ymax, label = paste('\u2190', ref.level), hjust = 'left') +
    annotate('text', x = 1, y = ymax, label = paste(target.level, '\u2192'), hjust = 'right')

  if(!is.null(palette)) p <- p + scale_fill_manual(values=palette)
  if(n.significant.cells > 0) p <- p + geom_vline(xintercept=nrow(loadings) - n.significant.cells + 0.5, color='red')


  if(show.pvals){
    d <- data.frame(y=pval, x=names(pval), row = 1:length(pval))
    p.pval <- ggplot(d, aes(x=reorder(x,row), y=-log(y,base = 10), fill=factor(x))) +
      geom_bar(stat="identity") +
      geom_hline(yintercept=-log(signif.threshold,base = 10)) +
      coord_flip() + labs(x='', y='-log(adj. p-value)') +
      plot.theme + theme(legend.position = "none") + theme(axis.text.y = element_blank())

    if(!is.null(palette)) p.pval <- p.pval + scale_fill_manual(values=palette)

    p.combo <- cowplot::plot_grid(p, p.pval, nrow=1, rel_widths=c(2,1))
    return(p.combo)
  }

  return(p)

}
