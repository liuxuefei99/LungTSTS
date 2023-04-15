densitypic=function(seu,cluster){
  plotlis=list()
  
  emb=seu@reductions$umap@cell.embeddings
  for( iterm in unique(seu@meta.data[,cluster])){
    
    #iterm="PRAD"
    
    data2=emb[seu@meta.data[seu@meta.data[,cluster]==iterm,] %>% rownames(),]
    colnames(data2)=c('x','y')
    data2=as.data.frame(data2)
    
    # raster polygon
    xlims=c(min(emb[,1])+0.3,max(emb[,1])+0.1)
    ylims=c(min(emb[,2])+0.3,max(emb[,2])+0.1)
    gg=ggplot(data2, aes(x=x, y=y) ) + xlim(xlims) + ylim(ylims)+
      # theme_classic() +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
      stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)  +xlab('')+ylab('')+
      # scale_fill_distiller(palette='Blues', direction=0.1,expand = c(0, 0))+   # +ggtitle(iterm)
      scale_fill_viridis(option='B',alpha = 1,direction=1) +  #ggtitle(iterm)+
      # scale_color_brewer(palette = "Dark2", direction=-1)+
      geom_point(aes(x=x,y=y), col='#FCFDBFFF',size=0.00001,alpha=0.2)+
      theme(
        legend.position='none',
        panel.border = element_blank(),
        plot.margin = margin(0,0,0,0,"cm")
      )
    
    gg=gg+theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))
    #+grids(linetype = "dashed")
    
    gg <- gg + theme(axis.title.x=element_blank(),
                     # axis.text.x=element_text(size=7,margin = margin(0, unit = "cm")),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_line(),
                     axis.title.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     # axis.text.y=element_text(size=7,margin = margin(0, unit = "cm")),
                     axis.text.y=element_blank(),
                     axis.ticks.length = unit(0, "cm")
    )
    con.color <- 'gray40'
    con.size=0.5
    cl='white'
    #ggsave(paste0(save.pic,'/Fig3.',iterm,'.density.pdf'),gg,width = 2.5,height=2.5)
    plotlis[[iterm]]=gg+labs(title=iterm)
  }
  return(plotlis)
}
