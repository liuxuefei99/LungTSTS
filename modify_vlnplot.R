modify_vlnplot<- function(obj,feature,group,color,pt.size = 0,plot.margin = unit(c(-0.1, 0, -0.1, 0), "cm"),...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = group,assay='RNA',...)  +
    ylab(feature) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          #axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )+scale_fill_manual(values = color)+theme(axis.title.x = element_blank(),plot.title = element_blank())
  return(p)
}
## main function
StackedVlnPlot<- function(obj, features,group,color,pt.size = 0,plot.margin = unit(c(-0.1, 0, -0.1, 0), "cm"),...) {
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, group,color))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle=30), axis.ticks.x = element_line())
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}     