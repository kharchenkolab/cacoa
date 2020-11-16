


#' @title Plot inter-sample expression distance 
#' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
#' @param name Test results to plot (default=expression.shifts)
#' @param notch Show notches in plot, see ggplot2::geom_boxplot for more info (default=T)
#' @param cell.groups Named factor with cell names defining groups/clusters (default: stored vector)
#' @param sample.groups Named sample factor with cell names (default: stored vector)
#' @param alpha transparency level for the individual points (default: 0.2)
#' @param weighted.distance whether to weigh the expression distance by the sizes of cell types (default: FALSE)
#' @param palette a set of colors to use for the conditions
#' @param show.significance whether to show statistical significance betwwen sample groups. wilcox.test was used; (* < 0.05; ** < 0.01; *** < 0.001)  
#' @return A ggplot2 object
plotExpressionDistance <- function(cluster.shifts, cell.groups = NULL, sample.groups = NULL, weighted.distance = FALSE,  notch= TRUE, alpha=0.2, min.cells = 10, palette=NULL, show.significance = FALSE) {
  ctdml <- cluster.shifts$ctdml
  valid.comparisons <- cluster.shifts$valid.comparisons
  if (!weighted.distance) { #
    df <- do.call(rbind,lapply(ctdml,function(ctdm) {
      x <- lapply(ctdm, function(xm) {
        nc <- attr(xm, 'cc')
        wm <- outer(nc, nc, FUN = 'pmin')
        cross.factor <- outer(sample.groups[rownames(xm)], sample.groups[colnames(xm)], '==')
        frm <- valid.comparisons[rownames(xm), colnames(xm)] & cross.factor
        diag(xm) <- NA
        # remove self pairs  
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
    df <- do.call(rbind, tapply(1:nrow(df), paste(df$Var1, df$Var2, df$type, sep = '!!'),
                                function(ii) {
                                  ndf <- data.frame(df[ii[1], , drop = F])
                                  ndf$value <- median(df$value[ii])
                                  ndf$n <- median(df$n[ii])
                                  ndf
                                }))
    df$group <- df$type1
    gg <- ggplot(na.omit(df), aes( x = type, y = value, dodge = group, fill = group )) +
      theme_bw() + geom_boxplot(notch = notch, outlier.shape = NA,) + #ggtitle(cell.type) +
      geom_point(
        position = position_jitterdodge(jitter.width = 0.1),
        size = 0.5,
        color = adjustcolor("black", alpha = alpha)
      ) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 90, hjust = 0.5)
      ) + theme(legend.position = "top") + xlab("") + ylab("expression distance") +
      scale_y_continuous( expand=c(0, max(df2$value) * 0.1), limits=c(0, (max(df2$value) + max(df2$value) * 0.01 )))  #expand=c(0, 0),
  
      if(show.significance) gg <- gg + stat_compare_means(aes(group = group), label = "p.signif")  # willcox test
    
  } else { # weighted expression distance 
    df <- do.call(rbind, lapply(ctdml, function(ctdm) {
      # bring to a common set of cell types
      commoncell <- unique( unlist( lapply(ctdm, function(x) colnames(x)) ))
      
      ctdm <-  lapply(ctdm, function(x) {
        y <- matrix(0,nrow=length(commoncell),ncol=length(commoncell)); rownames(y) <- colnames(y) <- commoncell; # can set the missing entries to zero, as they will carry zero weights
        y[rownames(x),colnames(x)] <- x;
        ycct <- setNames(rep(0,length(commoncell)), commoncell);
        ycct[colnames(x)] <- attr(x,'cc')
        attr(y, 'cc') <- ycct
        y
      }) # reform the matrix to make sure all cell type have the same dimensions
      
      x <- abind::abind(lapply(ctdm, function(x) {
        nc <- attr(x, 'cc')
        #wm <- (outer(nc,nc,FUN='pmin'))
        wm <- sqrt(outer(nc, nc, FUN = 'pmin'))
        return(x * wm)
      }), along = 3)
      
      # just the weights (for total sum of weights normalization)
      y <- abind::abind(lapply(ctdm, function(x) {
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
    
    gg <- ggplot(na.omit(df2), aes(x = group, y = value)) + theme_bw() + 
      geom_boxplot(notch = notch, outlier.shape = NA , aes(fill = group)) + ggtitle('')+
      geom_jitter(position = position_jitter(0.2), color = adjustcolor("black", alpha = alpha)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) + 
      theme(legend.position = "right") +
      xlab("") + ylab("expression distance") +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10)
      )
    if(show.significance) gg <- gg + stat_compare_means(label = "p.signif", label.x = 1.5) # willcox test
  }
  
  
  if(!is.null(palette)) { gg <- gg + scale_fill_manual(values=palette) }  
  return(gg)
}



#' @title Plot sample-sample expression distance in tSNE
#' @description  Plot results from cao$estimateExpressionShiftMagnitudes()
#' @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param cell.type Named of cell type, default is null, it set plot sample-sample expression distance in tSNE for the cell type
#' #' @param perplexity tSNE perpexity (default: 4)
#' @param max_iter tSNE max_iter (default: 1e3)
#' @param method dimension reduction methods (MDS or tSNE) (default is MDS)
#' @return A ggplot2 object

plotExpressionDistancetSNE <- function(cluster.shifts, sample.groups, method = 'MDS', cell.type = NULL,  perplexity=4, max_iter=1e3, palette=NULL) {
  ctdml <- cluster.shifts$ctdml
  if (!is.null(cell.type)) { # use distances based on the specified cell type
    title <- cell.type
    #if (is.null(cell.type) stop('please speficy cell type')
    df <- lapply(ctdml, function(ctdm) {
      xm <-  ctdm[[cell.type]]
      xm
    })
    dfm <- Reduce(`+`, df)/length(df)
  } else { # weighetd expression distance across all cell types
    title <- ''
    df <- lapply(ctdml, function(ctdm) {
      # bring to a common set of cell types
      commoncell <- unique( unlist( lapply(ctdm, function(x) colnames(x)) ))
      
      ctdm <-  lapply(ctdm, function(x) {
        y <- matrix(0,nrow=length(commoncell),ncol=length(commoncell)); rownames(y) <- colnames(y) <- commoncell; # can set the missing entries to zero, as they will carry zero weights
        y[rownames(x),colnames(x)] <- x;
        ycct <- setNames(rep(0,length(commoncell)), commoncell);
        ycct[colnames(x)] <- attr(x,'cc')
        attr(y, 'cc') <- ycct
        y
      }) # reform the matrix to make sure all cell type have the same dimensions

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
      xd
    })
    dfm <- Reduce(`+`, df)/length(df)
  }
  if (method == 'tSNE'){
    xde <- Rtsne::Rtsne(dfm, is_distance = TRUE, perplexity = perplexity, max_iter = max_iter)$Y
  } else if (method == 'MDS') {
    xde <- cmdscale(dfm, eig=TRUE, k=2)$points # k is the number of dim
  } else {
    stop("unknown embedding method")
  }
  df <- data.frame(xde)
  rownames(df) <- rownames(dfm)
  colnames(df) <- c("x", "y")
  df$fraction <- sample.groups[rownames(df)]
  df$sample=rownames(df)
  #df$ncells <- nc[rownames(df)]
  gg <- ggplot(df, aes(x, y, color=fraction, shape=fraction)) + geom_point(size=5) + #, size=log10(ncells)
    theme_bw() + ggtitle(title) + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  if(!is.null(palette)) { gg <- gg + scale_color_manual(values=palette) }  
  return(gg)
}
