

#' @param emb cell embedding matrix
#' @param anoSample   sample annotation  #Could do with better explanation
#' @param sample.per.cell  Named sample factor with cell names (default: stored vector)
#' @param sample.groups @param sample.groups A two-level factor on the sample names describing the conditions being compared (default: stored vector)
#' @param ref.level Reference sample group, e.g., ctrl, healthy, or untreated. (default: stored value)
#' @param target.level target/disease level for sample.group vector
#' @param bins number of bins for density esitmation, default 200
#' @param normlization  if TRUE, quantiles normlization will applied to indivisual sample.
#' @add.ponits add.ponits  show cells in density plot 

plotDensity=function(emb, anoSample,sample.per.cell, sample.groups,ref.level,target.level, bins = 200,add.ponits=TRUE,normlization=NULL) {
  
  if(is.null(emb)) stop("'emb' must be provided either during the object initialization or during this function call")
  
  if(is.null(anoSample)) stop("'anoSample' must be provided either during the object initialization or during this function call")
  
  if(is.null(sample.groups)) stop("'sample.groups' must be provided either during the object initialization or during this function call")
  
  
  if(is.null(ref.level)) stop("'ref.level' must be provided either during the object initialization or during this function call")
  
  if(is.null(target.level)) stop("'target.level' must be provided either during the object initialization or during this function call")
  
  
  cname = intersect(names(anoSample),rownames(emb))
  fraction=sample.per.cell[cname]
  anoSample=anoSample[cname]
  emb=emb[cname,]

  
  if (!is.null(normlization)){
      print('quantiles normlization')
      list.den = lapply(sn(as.character(unique(anoSample))), function(x) {
        nname = names(anoSample[anoSample == x])
        tmp = emb[nname, ]
        f2 = kde2d(tmp[, 1], tmp[, 2], n = bins, lims = c(range(emb[, 
                                                                    1]), range(emb[, 2])))
        f2
      })
      denMatrix = do.call("cbind", lapply(list.den, function(x) as.numeric(x$z)))
      denMatrix.normal = normalize.quantiles(denMatrix)
      colnames(denMatrix.normal) = colnames(denMatrix)
      list.den2 = lapply(sn(as.character(unique(sample.groups))), 
                         function(x) {
                           tmp = denMatrix.normal[, names(sample.groups[sample.groups == 
                                                                          x])]
                           matrix(rowMeans(tmp), ncol = bins, byrow = FALSE)
                         })
  }else{

    list.den = lapply(sn(as.character(unique(fraction))), function(x) {
      nname = names(fraction[fraction == x])
      tmp = emb[nname, ]
      f2 = kde2d(tmp[, 1], tmp[, 2], n = bins, lims = c(range(emb[, 
                                                                  1]), range(emb[, 2])))
      f2
    })
    denMatrix = do.call("cbind", lapply(list.den, function(x) as.numeric(x$z)))
    list.den2 = lapply(sn(as.character(unique(sample.groups))), 
                       function(x) {
                         matrix(denMatrix[,x], ncol = bins, byrow = FALSE)
                       })
    
    }
  mi = min(unlist(lapply(list.den2, function(x) min(as.numeric(x)))))
  ma = max(unlist(lapply(list.den2, function(x) max(as.numeric(x)))))
  fig.list = list()
  for (iterm in levels(fraction)) {
    p = list.den2[[iterm]] %>% as_tibble() %>% rowid_to_column(var = "X") %>% 
      gather(key = "Y", value = "Z", -1) %>% mutate(Y = as.numeric(gsub("V", 
                                                                        "", Y))) %>% ggplot(aes(X, Y, fill = Z)) + geom_raster() +
      
      theme_bw() + theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(), panel.border = element_blank(), 
                         panel.background = element_blank(), plot.margin = margin(0.1, 
                                                                                  0.1, 0.1, 0.1, "cm")) + 
      theme(axis.title.x = element_blank(), 
                    axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                    axis.title.y = element_blank(), axis.text.y = element_blank(), 
                    axis.ticks.y = element_blank()) +   #geom_tile(color=NA) + #theme(panel.border=element_rect(fill = NA, colour=alpha('black', .5),size=10))+
      theme(legend.position = "none") +
      scale_fill_viridis(option='B',alpha = 1,direction=1, limits = c(mi, ma))+
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
      theme(panel.border = element_rect(fill=NA,color="black", size=0.5, 
                                        linetype="solid"))
    #p=p+geom_point(data=emb,aes(x=x,y=y), col='#FCFDBFFF',size=0.00001,alpha=0.2)

    if (!is.null(add.ponits)){
      
      x=emb[,1]
      y=emb[,2]
      
      x=(x-range(x)[1])
      x=(x/max(x))*bins
      
      y=(y-range(y)[1])
      y=(y/max(y))*bins
      
      emb2=data.frame(x=x,y=y,Z=1)
      nname=names(fraction[fraction==iterm])
      
      nname=sample(nname,min(2000,nrow(emb2)))
      emb2=emb2[nname,]
      p=p+geom_point(data=emb2,aes(x=x,y=y), col='#FCFDBFFF',size=0.00001,alpha=0.2)+ggtitle(iterm)
    }
    #p=p+cn+cn2+cn3+cn4+cn5+cn6
    #scale_fill_distiller(palette = 2, direction = 0.1, expand = c(0, 0), limits = c(mi, ma))
    fig.list[[iterm]] = p
  }
  
  
  diff=list.den2[[target.level]]-list.den2[[ref.level]]
  
  fig.list[['Diff']] = diff %>% as_tibble() %>% rowid_to_column(var = "X") %>% 
    gather(key = "Y", value = "Z", -1) %>% mutate(Y = as.numeric(gsub("V", 
                                                                      "", Y))) %>% ggplot(aes(X, Y, fill = Z)) + geom_raster() +
    
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), panel.border = element_blank(), 
                       panel.background = element_blank(), plot.margin = margin(0.1, 
                                                                                0.1, 0.1, 0.1, "cm")) + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(), axis.text.y = element_blank(), 
          
          axis.ticks.y = element_blank()) +   #geom_tile(color=NA) + #theme(panel.border=element_rect(fill = NA, colour=alpha('black', .5),size=10))+
    theme(legend.position = "none") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0) +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+ggtitle('Diff')
  #scale_fill_viridis(discrete=FALSE)
  #scale_fill_viridis(option='B',alpha = 1,direction=1, limits = c(mi, ma))

  
  
  
  b = cowplot::plot_grid(plotlist = fig.list, ncol = 3, nrow = 1)
  

  return(b)

  
}





