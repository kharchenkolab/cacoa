


#' @title Plot inter-sample expression distance 
#' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
#' @param name Test results to plot (default=expression.shifts)
#' @param notch Show notches in plot, see ggplot2::geom_boxplot for more info (default=T)
#' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
#' @param sample.per.cell Named sample factor with cell names (default: stored vector)
#' @weight.disatnce default is null, it set caculated weigeted expression distance across mutiple cell types
#' @return A ggplot2 object
plotExpressionDistance <- function(cluster.shifts, notch = T, cell.groups = NULL, sample.per.cell = NULL, weight.disatnce = NULL,  min.cells = 10) {
  ctdml <- cluster.shifts$ctdml
  valid.comparisons <- cluster.shifts$valid.comparisons
  if (is.null(weight.disatnce)) {
    df <- do.call(rbind,lapply(ctdml,function(ctdm) {
      x <- lapply(ctdm, function(xm) {
        nc <- attr(xm, 'cc')
        wm <- outer(nc, nc, FUN = 'pmin')
        cross.factor <- outer(sample.groups[rownames(xm)], sample.groups[colnames(xm)], '==')
        frm <- valid.comparisons[rownames(xm), colnames(xm)] & cross.factor
        diag(xm) <- NA
        # restrict
        xm[!frm] <- NA
        xm[wm < min.cells] <- NA
        if (!any(!is.na(xm)))
          return(NULL)
        xmd <- na.omit(reshape2::melt(xm))
        wm[is.na(xm)] <- NA
        xmd$n <- na.omit(reshape2::melt(wm))$value
        return(xmd)
      })
      x <- x[!unlist(lapply(x, is.null))]
      df <- do.call(rbind, lapply(sccore:::sn(names(x)), function(n) {
        z <- x[[n]]
        z$type <- n
        z
      }))
      df$patient <- df$Var1
      df$type1 <- sample.groups[df$Var1]
      df$type2 <- sample.groups[df$Var2]
      df
    }))
    
    # median across pairs
    df <- do.call(rbind, tapply(1:nrow(df), paste(df$Var1, df$Var2, df$type, sep =
                                                    '!!'), function(ii) {
                                                      ndf <- data.frame(df[ii[1], , drop = F])
                                                      ndf$value <- median(df$value[ii])
                                                      ndf$n <- median(df$n[ii])
                                                      ndf
                                                    }))
    df$group <- df$type1
    gg <- ggplot(na.omit(df), aes( x = type, y = value, dodge = group, fill = group )) +
      theme_classic() + geom_boxplot(notch = TRUE, outlier.shape = NA,) + #ggtitle(cell.type) +
      geom_point(
        position = position_jitterdodge(jitter.width = 0.1),
        size = 0.5,
        color = adjustcolor("black", alpha = 0.2)
      ) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 90, hjust = 0.5)
      ) + theme(legend.position = "top") + xlab("") + ylab("expression distance")
    
  } else { # weighetd expression distance 
    df <- do.call(rbind, lapply(ctdml, function(ctdm) {
      n.cell <- unlist(lapply(ctdm,ncol)) %>% table() %>% sort() %>% rev %>% names %>% as.numeric # 
      ctdm <- ctdm[unlist(lapply(ctdm, ncol)) == n.cell[1]]
      genelists <- lapply(ctdm, function(x) colnames(x))
      commoncell <- Reduce(intersect, genelists)
      ctdm <-  lapply(ctdm, function(x) {
        cct <- attr(x, 'cc')
        x <- x[commoncell, commoncell]
        attr(x, 'cc') <- cct[commoncell]
        x
      }) # reform the matrix to make sure all cell type have the same diminsion  
      
      x <- abind(lapply(ctdm, function(x) {
        nc <- attr(x, 'cc')
        #wm <- (outer(nc,nc,FUN='pmin'))
        wm <- sqrt(outer(nc, nc, FUN = 'pmin'))
        return(x * wm)
      }), along = 3)
      
      # just the weights (for total sum of weights normalization)
      y <- abind(lapply(ctdm, function(x) {
        nc <- attr(x, 'cc')
        sqrt(outer(nc, nc, FUN = 'pmin'))
      }), along = 3)
      
      # normalize by total weight sums
      xd <- apply(x, c(1, 2), sum) / apply(y, c(1, 2), sum)
      dim(xd)
      
      cross.factor <- outer(sample.groups[rownames(xd)], sample.groups[colnames(xd)], '==')
      frm <- valid.comparisons[rownames(xd), colnames(xd)] & cross.factor
      
      diag(xd) <- NA # remove self pairs
      # restrict
      xd[!frm] <- NA
      if (!any(!is.na(xd)))
        return(NULL)
      xmd2 <- na.omit(reshape2::melt(xd))
      xmd2 <- na.omit(xmd2)
      xmd2$type1 <- sample.groups[xmd2$Var1]
      xmd2$type2 <- sample.groups[xmd2$Var2]
      xmd2
    }))
    
    df2 <- do.call(rbind, tapply(1:nrow(df), paste(df$Var1, df$Var2, sep = '!!'), function(ii) {
      ndf <- data.frame(df[ii[1], , drop = F])
      ndf$value <- median(df$value[ii])
      ndf$n <- median(df$n[ii])
      ndf
    }))
    
    df2$group = df2$type1
    gg <- ggplot(na.omit(df2), aes(x = group, y = value)) + theme_classic() + 
      geom_boxplot(notch = TRUE, outlier.shape = NA , aes(fill = group)) + ggtitle('')+
      geom_jitter(position = position_jitter(0.2), color = adjustcolor("black", alpha = 0.2)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) + 
      theme(legend.position = "right") +
      xlab("") + ylab("expression distance") +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10)
      )
  }
  return(gg)
}

