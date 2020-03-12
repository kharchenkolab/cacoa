##' Evaluate cell density between different sample groups.  Density is calculated per sample and normlized by quantiles normlization.
##'
##' @title cell embedding density
##' @param emb cell embedding matrix
##' @param anoSample   sample annotation  #Could do with better explanation
##' @param sample.groups  sample groups factor
##' @param bins  number of bins for density
#'  @import dplyr
#'  @import tibble
#'  @import tidyr
#'  @importFrom MASS kde2d
plotDensity <- function(emb,anoSample,sample.groups,bins=200){
  cname=names(anoSample)

  list.den=lapply(sn(as.character(unique(anoSample))),function(x) {
      nname=names(anoSample[anoSample==x])
      tmp=emb[nname,]
      f2=kde2d(tmp[,1], tmp[,2], n = bins, lims = c(range(emb[,1]), range(emb[,2])))  # estimate density
      f2
    })


  denMatrix=do.call('cbind',lapply(list.den,function(x) as.numeric(x$z)))
  dim(denMatrix)

  # quantiles normalization for density
  denMatrix.normal=normalize.quantiles(denMatrix)
  colnames(denMatrix.normal)=colnames(denMatrix)

  # average density
  list.den2=lapply(sn(as.character(unique(sample.groups))),function(x){

    tmp=denMatrix.normal[,names(sample.groups[sample.groups==x])]
    matrix(rowMeans(tmp), ncol=bins,byrow=FALSE)

  })

  mi=min(unlist(lapply(list.den2,function(x) min(as.numeric(x)))))
  ma=max(unlist(lapply(list.den2,function(x) max(as.numeric(x)))))


  fig.list=list()
  for (iterm in unique(sample.groups)){


    p=list.den2[[iterm]] %>%
      # Data wrangling
      as_tibble() %>%
      rowid_to_column(var="X") %>%
      gather(key="Y", value="Z", -1) %>%

      # Change Y to numeric
      mutate(Y=as.numeric(gsub("V","",Y))) %>%

      # Viz
      ggplot(aes(X, Y, fill= Z)) +
      geom_tile() +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0,0,0,0,"cm")) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()
      )+
      ggtitle(iterm)+
      theme(legend.position="none")+
      #  scale_fill_viridis(option='B',alpha = 1,direction=1,limits=lim)
      scale_fill_distiller(palette=2, direction=0.1,expand = c(0, 0),limits=c(mi,ma))

    fig.list[[iterm]]=p
  }

  b=  cowplot::plot_grid(plotlist=fig.list, ncol=2, nrow=ceiling(length(unique(sample.groups))/2))

  return(b)
}


