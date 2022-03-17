#' @keywords internal
jaccardPwTop <- function(subsamples, top.thresh){
  jac.all = c()
  idxs = c()
  for(i in 1:length(subsamples)) {
    for(j in 1:length(subsamples)) {
      if (j <= i) next
      set1 <- rownames(subsamples[[i]])[rank(subsamples[[i]]$pvalue) <= top.thresh]
      set2 <- rownames(subsamples[[j]])[rank(subsamples[[j]]$pvalue) <= top.thresh]
      if((length(set1) != 0) || (length(set2) != 0)) {
        jac.all <- c(jac.all, length(intersect(set1, set2)) / length(unique(c(set1, set2))))
      } else {
        jac.all <- c(jac.all, NA)
      }
      idxs <- c(idxs, paste(names(subsamples)[i], names(subsamples)[j]))
    }
  }
  return(list(jac = jac.all, id = idxs))
}

#' @keywords internal
jaccardPwPval <- function(subsamples, p.val.cutoff){
  jac.all = c()
  idxs = c()
  for(i in 1:length(subsamples)) {
    for(j in 1:length(subsamples)) {
      if (j <= i) next
      set1 <- rownames(subsamples[[i]])[subsamples[[i]]$padj <= p.val.cutoff]
      set2 <- rownames(subsamples[[j]])[subsamples[[j]]$padj <= p.val.cutoff]
      if((length(set1) != 0) || (length(set2) != 0)) {
        jac.all <- c(jac.all, length(intersect(set1, set2)) / length(unique(c(set1, set2))))
      } else {
        jac.all <- c(jac.all, NA)
      }
      idxs <- c(idxs, paste(names(subsamples)[i], names(subsamples)[j]))
    }
  }
  return(list(jac = jac.all, id = idxs))
}

#' @keywords internal
plotStability <- function(jaccards, notch, show.jitter, jitter.alpha, show.pairs, sort.order,
                          xlabel = '', ylabel = '', log.y.axis = FALSE, palette = NULL, plot.theme=theme_get(),
                          set.color=TRUE, set.fill=FALSE) {
  jaccards$group <- as.factor(jaccards$group)

  if(!set.color){
    jaccards$cmp <- 'none'
  }
  if(set.fill){
    jaccards$fill <- jaccards$group
  } else {
    jaccards$fill <- 'none'
  }

  if(! show.pairs) {
    if(sort.order) {
      p <- ggplot(jaccards, aes(x=reorder(group, value, na.rm = TRUE), y=value,
                                group=group, fill=fill)) + geom_boxplot(outlier.shape = NA, notch = notch) +
        geom_jitter(alpha=jitter.alpha)
    } else {
      p <- ggplot(jaccards, aes(x=group, y=value,
                                group=group, fill=fill)) + geom_boxplot(outlier.shape = NA, notch = notch) +
        geom_jitter(alpha=jitter.alpha)
    }
  } else {
    if(sort.order) {
      p <- ggplot(jaccards, aes(x=reorder(group, value, na.rm = TRUE), y=value,
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

#' @keywords internal
estimateStabilityPerCellType <- function(de.res, top.n.genes, p.val.cutoff) {
  data.all <- data.frame()
  for(cell.type in names(de.res$initial)){
    # print(cell.type)
    subsamples.names <- setdiff(names(de.res), 'initial')

    # Remove some subsamples due to the min.cell.counts
    # compare resampling results with "initial"
    jacc.init <- c()
    for(subs.name in subsamples.names){
      subsamples.tmp = list(de.res[[subs.name]][[cell.type]],
                            de.res$initial[[cell.type]]$res)
      jacc.pw.tmp <- jaccardPwTop(subsamples.tmp, 200)
      jacc.init <- c(jacc.init, jacc.pw.tmp$jac)  # please remain 200 here - it is only a technical thing
    }
    subsamples.names <- subsamples.names[(jacc.init != 0) & (jacc.init != 1)]
    subsamples.tmp <- lapply(subsamples.names, function(s){de.res[[s]][[cell.type]]})
    names(subsamples.tmp) <- subsamples.names

    if(length(subsamples.names) <= 2) next

    # Calculate jaccard
    if (is.null(p.val.cutoff)) {
      jacc.tmp <- jaccardPwTop(subsamples.tmp, top.n.genes)
    } else {
      jacc.tmp <- jaccardPwPval(subsamples.tmp, p.val.cutoff)
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


