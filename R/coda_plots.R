#' @import ape
NULL

#' @keywords internal
estimateCdaSpace <- function(d.counts, d.groups, thresh.pc.var = 0.95, n.dim = 2){
  cell.loadings <- c()
  sample.pos <- c()

  bal <- getRndBalances(d.counts)
  d.used <- bal$norm

  if (ncol(d.used) != 2) {
    for (i in 1:n.dim) {
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

#' @keywords internal
plotCodaSpaceInner <- function(df.space, df.loadings, d.groups, ref.level, target.level, font.size=3, palette=NULL) {
  group.names <- c(ref.level, target.level)
  rda.plot <- ggplot(df.space, aes(x=S1, y=S2)) +
    #   geom_text(aes(label=rownames(df_pca) %in% samplegroups$trgt),size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    geom_point(aes(colour = factor(group.names[d.groups + 1] ))) +
    labs(colour="Condition")

  if (!is.null(palette)) rda.plot <- rda.plot + scale_color_manual(values=palette)

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

#' helper function for creating dendograms
#' @keywords internal
createDendrogram <- function(dend.data, angle=90, plot.theme=theme_get(), font.size=3, hjust=1) {
  ggplot() +
    geom_segment(data = dend.data$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + plot.theme +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"),
          axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5, margin = margin()),
          panel.grid = element_blank(), panel.border=element_blank(),
          axis.line=element_blank()) +
    scale_x_continuous(breaks=dend.data$labels$x, labels=dend.data$labels$label) +
    scale_y_continuous(expand=c(0, 0, 0.1, 0))
}

#' helper function for creating dendograms
#' @keywords internal
distTreeOrder <- function(t, tree.order){
  t.labels <- t$tip.label[t$edge[t$edge[,2] < min(t$edge[,1]),2]]
  t.labels <- t.labels[t.labels %in% tree.order]
  s <- 0
  for(i in 1:length(t.labels)){
    j <- which(tree.order == t.labels[i])
    s <- s + abs(i - j) ^ 2
  }
  return(s)
}


#' helper function for creating dendograms
#' @keywords internal
plotContrastTree <- function(d.counts, d.groups, ref.level, target.level, plot.theme=theme_get(), label.angle=90,

                             p.threshold=0.05, adjust.pvalues=TRUE, h.methods='both', font.size=3, label.hjust=1,
                             show.text='pvalue.code', tree.order=NULL, loadings.mean=NULL, palette=NULL,
                             verbose=FALSE) {
  checkPackageInstalled(c("ggdendro"), cran=TRUE)

  log.f <- getLogFreq(d.counts)

  if (h.methods == 'up') {
    if (verbose) message('up')
    t.cur <- constructTreeUp(d.counts, d.groups)
  } else if (h.methods == 'down') {
    if (verbose) message('down')
    t.cur <- constructTree(d.counts, d.groups)
  } else {
    if (verbose) message('up and down')
    t.cur <- constructTreeUpDown(d.counts, d.groups)
  }

  # t.cur <- constructBestPartitionTree(d.counts, d.groups)

  if (!is.null(loadings.mean) && is.null(tree.order)) {
    tree.order <- names(sort(loadings.mean))
  }

  # Order the tree in the as similar as possible way
  if (!is.null(tree.order)) {
    t <- t.cur$tree
    tree.order <- intersect(tree.order, t$tip.label)
    # distance of the initial tree
    d <- distTreeOrder(t, tree.order)
    idx <- min(t$edge[,1]):max(t$edge[,1])
    for (i.node in idx) {

      # alternative tree
      t.alt <- ape::rotate(t, i.node)
      d.alt <- distTreeOrder(t.alt, tree.order)

      if (d.alt <= d) {
        t <- t.alt
        d <- d.alt
      }
    }
    t.cur$tree <- t
    t.cur$dendro <- compute.brlen(t.cur$tree, method="Grafen") %>% as.hclust() %>% as.dendrogram()
  }


  tree <- t.cur$tree
  sbp <- sbpInNodes(tree)
  # sbp = t.cur$sbp

  # ---------------------------------
  # Positions of the dendrogram
  dend.data <- ggdendro::dendro_data(t.cur$dendro, type = "rectangle")
  types.order <- dend.data$labels$label
  # change the sign of sbp corresponding to the dendrodgram
  for (k in 1:ncol(sbp)) {
    p <- sbp[, k]
    type.plus <- rownames(sbp)[p > 0]
    type.minus <- rownames(sbp)[p < 0]
    if (which(types.order == type.plus[1]) < which(types.order == type.minus[1])) {
      sbp[, k] <- -sbp[, k]
    }
  }

  # Get positions of segments (nodes) of the dendrogram
  # only inner nodes
  node.pos <- dend.data$segments %$% .[(y == yend) & (yend != 0),]
  node.pos$id <- tree$edge[,1]  # id of the inner node
  node.pos$to <- tree$edge[,2]


  # Positions of inner nodes
  innode.pos <- unique(node.pos[,c('x','y','id')])
  rownames(innode.pos) <- innode.pos$id
  # Ranges of inner nodes (length of branches)
  innode.pos$range <- -1
  for (i in 1:nrow(innode.pos)){
    tmp <- node.pos$xend[node.pos$id == innode.pos$id[i]]
    innode.pos$range[i] <- max(tmp) - min(tmp)
  }
  innode.pos <- innode.pos[order(innode.pos$id),]

  # ----------------------------------------

  # Balances
  balances <- getNodeBalances(log.f, sbp)
  colnames(balances) <- rownames(innode.pos)

  p.val <- c()
  for (i in 1:ncol(balances)) {
    aov.data <- data.frame(balance = balances[,i], group = d.groups)
    # anova
    # res <- aov(balance ~ group, data=aov.data)
    # res <- aov(group ~ balance, data=aov.data)
    # p.val <- c(p.val,summary(res)[[1]][1,5])

    mod <- lm(group ~ balance, data=aov.data)
    res <- summary(mod)
    p.val <- c(p.val, res$coefficients[2,4])

  }
  p.val[is.na(p.val)] <- 1

  p.adj <- if (adjust.pvalues) p.adjust(p.val, method='fdr') else p.val

  px.init <- createDendrogram(dend.data, plot.theme=plot.theme, font.size=font.size, angle=label.angle,
                              hjust=label.hjust)

  if (sum(p.adj < p.threshold) == 0)
    return(px.init)

  df.pval <- data.frame()
  df.bals <- data.frame()
  df.bal.quartile <- data.frame()
  df.bal.median <- data.frame()
  df.text <- data.frame()
  group.levels <- c(ref.level, target.level)

  for(id.node in 1:ncol(balances)){
    if (p.adj[id.node] > p.threshold) next

    df.pval <- rbind(df.pval, innode.pos[id.node, c('x', 'y')])

    # Normalization of values of balances according to the range between nodes
    x.tmp <- balances[,id.node]
    x.tmp <- x.tmp - mean(x.tmp)

    # DO NOT MOVE THIS LINE DOWN:
    df.text %<>% rbind(data.frame(
      x=innode.pos$x[id.node] + innode.pos$range[id.node]/2,
      y=innode.pos$y[id.node],
      range=sprintf('%2.1f', max(abs(x.tmp))),
      pvalue=signif(p.adj[id.node], 2),
      pvalue.code=pvalueToCode(p.adj[id.node])
    ))

    x.tmp <- x.tmp / max(abs(x.tmp)) / 2 * innode.pos$range[id.node] * 0.9
    x.tmp <- x.tmp + innode.pos$x[id.node]
    y.tmp <- d.groups * 0.03 + innode.pos$y[id.node] - 0.05

    q.case <- quantile(x.tmp[d.groups], c(0.25, 0.75))
    q.control <- quantile(x.tmp[!d.groups], c(0.25, 0.75))

    df.bal.quartile %<>% rbind(data.frame(
      x=c(q.case, q.control),
      y=c(rep(y.tmp[d.groups][1],2), rep(y.tmp[!d.groups][1],2)),
      group=group.levels[1 + c(1, 1, 0, 0)],
      node = id.node
    ))

    df.bal.median %<>% rbind(data.frame(
      x=c(median(x.tmp[d.groups]), median(x.tmp[!d.groups])),
      y=c(y.tmp[d.groups][1], y.tmp[!d.groups][1]),
      group=group.levels[1 + c(1, 0)],
      node=id.node
    ))

    df.bals %<>% rbind(data.frame(x=x.tmp, y=y.tmp, group=group.levels[1 + d.groups], node = id.node))
  }

  px <- px.init + geom_point(data=df.bals,
                             aes(x=x, y=y, col=as.factor(group), group=as.factor(node)), alpha = 0.1, size = 1) +
    geom_point(data=df.bal.median,
               aes(x=x, y=y, col=as.factor(group)),
               size = 2.5, shape = 18) +
    geom_point(data=df.pval, aes(x=x, y=y)) +
    geom_line(data=df.bal.quartile, aes(x=x, y=y,
                                        col=as.factor(group),
                                        group=interaction(group, node)), size = 0.75) +
    labs(col=" ")

  if (is.logical(show.text) && !show.text) show.text <- NULL
  if (!is.null(show.text)) {
    if (!(show.text %in% colnames(df.text)))
      stop("Unexpected value for show.text: ", show.text)

    px <- px +
      geom_text(data=df.text, mapping=aes_string(x='x', y='y', label=show.text), vjust=0, hjust=0, size=font.size)
  }

  if (!is.null(loadings.mean)) {
    node.leaves <- node.pos[node.pos$to < min(tree$edge[,1]),]
    node.leaves$label <- tree$tip.label[node.leaves$to]
    node.leaves$loadings <- loadings.mean[node.leaves$label]

    # node.leaves <- node.leaves[node.leaves$pval < 0.05,]
    px <- px + geom_tile(data=node.leaves, mapping=aes(x=xend, y=-0.1, fill=loadings, width = 0.5, height = 0.1)) +
      guides(fill=guide_colorbar(title='loadings', title.position="top", direction="horizontal", title.hjust = 0.5))
  }

  if (!is.null(palette)) {
    px <- px + scale_color_manual(values=palette) +
      scale_fill_gradient2(low=palette[ref.level], high=palette[target.level], mid='grey80', midpoint=0)
  }

  return(px)
}

plotCellLoadings <- function(loadings, pval, ref.level, target.level, signif.threshold=0.05, jitter.alpha=0.1,
                             palette=NULL, show.pvals=FALSE, plot.theme=theme_get(), jitter.size=1,
                             ordering=c("pvalue", "loading"), ref.load.level=0, annotation.x=NULL, annotation.y=1) {
  ordering <- match.arg(ordering)
  xintercept <- ref.load.level

  loading.order <- order(abs(rowMeans(loadings)))
  loadings <- loadings[loading.order, ]
  if (!is.null(pval)) {
    # if some p-values are the same - then order by mean, therefore a prior sorting is required
    pval <- pval[loading.order]

    if (ordering == "pvalue") {
      # additional ordering by p-value
      loadings <- loadings[order(-pval), ]
      pval <- pval[order(-pval)]
    }

    # Get significant cells
    n.significant.cells <- sum(pval <= signif.threshold)
  } else {
    n.significant.cells <- 0
  }

  # Normalization of loadings
  loadings <- loadings - xintercept

  res.ordered <- t(loadings) %>% as.data.frame()
  if (is.null(annotation.x)) {
    annotation.x <- max(loadings)
  }

  p <- ggplot(stack(res.ordered), aes(y=ind, x=values, fill=factor(ind))) +
    geom_boxplot(notch=TRUE, outlier.shape = NA)

  if ((jitter.alpha > 1e-5) && (jitter.size > 1e-5)) {
    p <- p + geom_jitter(alpha=jitter.alpha, size=jitter.size)
  }

  p <- p + geom_vline(xintercept=0, color = "grey37") +
    labs(y='', x='loadings') + plot.theme + theme(legend.position = "none") +
    scale_y_discrete(position = "right") + xlim(-1, 1)

  # Add text
  p <- p +
    annotate('text', x=-annotation.x, y=annotation.y, label=paste('\u2190', ref.level), hjust='left') +
    annotate('text', x=annotation.x, y=annotation.y, label=paste(target.level, '\u2192'), hjust='right')

  if (!is.null(palette)) p <- p + scale_fill_manual(values=palette)
  if ((n.significant.cells > 0) && (ordering == "pvalue")) {
    p <- p + geom_hline(yintercept=nrow(loadings) - n.significant.cells + 0.5, color='red')
  }

  if (show.pvals) {
    d <- data.frame(y=-log(pval, base=10), x=names(pval), row=1:length(pval))
    p.pval <- ggplot(d, aes(y=reorder(x, row), x=y, fill=factor(x))) +
      geom_bar(stat="identity") +
      geom_vline(xintercept=-log(signif.threshold, base=10)) +
      scale_x_continuous(expand=c(0, 0, 0.02, 0)) +
      labs(y='', x='-log10(adj. p-value)') +
      plot.theme + theme(legend.position="none") + theme(axis.text.y=element_blank())

    if (!is.null(palette)) p.pval <- p.pval + scale_fill_manual(values=palette)

    p.combo <- cowplot::plot_grid(p, p.pval, nrow=1, rel_widths=c(2, 1))
    return(p.combo)
  }

  return(p)

}