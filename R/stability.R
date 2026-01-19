#' @keywords internal
jaccardPwTop <- function(subsamples, top.thresh){
  jac.all = c()
  idxs = c()
  for(i in 1:length(subsamples)) {
    for(j in 1:length(subsamples)) {
      if (j <= i) next
      d1 <- subsamples[[i]] ; d2 <- subsamples[[j]]
      if (is.null(d1) || is.null(d2)) next
      r1 <- rank(d1$pvalue, na.last = "keep")
      r2 <- rank(d2$pvalue, na.last = "keep")
      set1 <- rownames(d1)[!is.na(r1) & (r1 <= top.thresh)]
      set2 <- rownames(d2)[!is.na(r2) & (r2 <= top.thresh)]
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

      d1 <- subsamples[[i]] ; d2 <- subsamples[[j]]
      if (is.null(d1) || is.null(d2)) next

      set1 <- rownames(d1)[!is.na(d1$padj) & (d1$padj <= p.val.cutoff)]
      set2 <- rownames(d2)[!is.na(d2$padj) & (d2$padj <= p.val.cutoff)]
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

  for (cell.type in names(de.res)) {
    df0 <- de.res[[cell.type]]
    subs <- attr(df0, "subsamples")

    # need enough resamples to do pairwise comparisons
    if (!is.list(subs) || length(subs) < 3) next

    # optional: drop NULL / empty subs
    ok <- vapply(subs, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))
    subs <- subs[ok]
    if (length(subs) < 3) next

    # Calculate jaccard across subsamples
    if (is.null(p.val.cutoff)) {
      jacc.tmp <- jaccardPwTop(subs, top.n.genes)
    } else {
      jacc.tmp <- jaccardPwPval(subs, p.val.cutoff)
    }

    if (is.null(jacc.tmp$jac) || length(jacc.tmp$jac) == 0) next

    data.tmp <- data.frame(
      group = cell.type,
      value = jacc.tmp$jac,
      cmp = jacc.tmp$id
    )
    data.all <- rbind(data.all, data.tmp)
  }

  data.all
}
