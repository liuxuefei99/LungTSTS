
library("ggsci");library(monocle);library(parallel);library(magrittr);library(ggcorrplot);library(ggthemes)
library(DoubletFinder);library(copykat);#library(CaSpER);library(GenomicRanges);library(future)
library(presto);library(harmony);library(ggrepel)
library(stringr);library(reshape2);library(plyr);library(Seurat);library(dplyr);
library(Matrix);library(ggplot2);library(edgeR);library(data.table);library(pheatmap);
library(msigdbr);library(fgsea);library(singleseqgset);library(Hmisc)
library(velocyto.R)
library(doBy)
library(clusterProfiler);library(org.Hs.eg.db);library(Rcpp)
library(CellChat);library(nichenetr)
library(survival);
library(survminer)
library("glmnet")
library(forestplot)
library(maftools)
library(smoother)
library(plotrix)
library(ggpubr)
library(RColorBrewer)
library(sampling)
library(velocyto.R)
library(URD)
library(nichenetr)
library(corrplot)
library(ConsensusClusterPlus)


my36colors <-c('#D6E7A3',"orange1", 'lightblue','#7142AC',"darkcyan","royalblue1","red3",'#53A85F',"deeppink",
               "mediumvioletred","gold","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1")

mytheme <- theme_bw() + 
  theme(plot.title=element_text(size=rel(2),hjust=0.5),
        axis.title=element_text(size=rel(1)),
        axis.text.x = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),
        axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),
        panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        panel.border=element_rect(color="black",size=1),
        axis.line=element_line(color="black",size=0.5))

set.seed(1234)



#---+++Figure1-------
counts.re=readRDS(paste0(save.data,'Alldata.third.rds'))
counts.re<-  RunUMAP(object = counts.re, dims = 1:50, reduction = "harmony")
p=DimPlot(counts.re,group.by = 'recluster')+mytheme+scale_color_manual(values = my36colors)
ggsave(paste0(save.pic,'Figures1-cluster.pdf'),height = 5,width=6)
ggsave(paste0(save.pic,'Figures1-cluster.jpg'),height = 5,width=6)



genes.choose=c('EPCAM','KRT18','PDGFRA','COL6A1','CD14','CD68','TNFRSF17','IGHA1','CD3D','CD3E','CD19','MS4A1','PECAM1','VWF','CPA3','MS4A2')
p=DotPlot(counts.re,group.by = 'recluster' ,features = c(genes.choose) ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#fa1414')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=0.5, vjust=0.5,size = rel(1),color = 'black'))

ggsave(paste0(save.pic,'Figures1-dotplot.pdf'),height = 4,width=7)

p=FeaturePlot(counts.re,'MKI67',max.cutoff = 2)+mytheme
ggsave(paste0(save.pic,'Figures1-MKI67.jpg'),height = 5.5,width=6)

#-------CD4------
count.CD4=readRDS(paste0(save.data,'CD4.rds'))




#分群
p=DimPlot(count.CD4,group.by = 'recluster')+mytheme+scale_color_manual(values=my36colors)
ggsave(paste0(save.pic,'Figure1-CD4.cluster.pdf'),height = 5,width=6)
ggsave(paste0(save.pic,'Figure1-CD4.cluster.jpg'),height = 5,width=6)

#点图
genes.choose=c('CCR7','SELL','CXCL13','PDCD1','GNLY','FGFBP2','CD69','FOS','FOXP3','IL2RA','GZMA','CCL5','IL7R','GZMK','ISG15','IFIT3','STMN1','MKI67')
p=DotPlot(count.CD4,group.by = 'recluster' ,features = c(genes.choose) ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#fa1414')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=0.5, vjust=0.5,size = rel(1),color = 'black'))
ggsave(paste0(save.pic,'Figure1-CD4.dotplot.pdf'),height = 4,width=6.5)


choose.gene=c('CXCL13','ZNF683','TOX2','TOX','LAIR2','CXCR6','LAYN','ACP5','MAGEH1')
p=DotPlot(count.CD4,group.by = 'Tissue1' ,features = c(choose.gene) ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#fa1414')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=0.5, vjust=0.5,size = rel(1),color = 'black'))

ggsave(paste0(save.pic,'Figure1add-CD4.TN dotplot.pdf'),height = 2,width=4,p)


#热图
markers=wilcoxauc(count.CD4,group_by = 'recluster') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
top2heatmap =markers %>% group_by(group) %>% top_n(n=15,wt=logFC) %>%data.frame() %>% dplyr::arrange(group)
marker_heatmap(top2heatmap,count.CD4,save.pic,'recluster','heatmap.CD4', levels(as.factor(top2heatmap$group)) ,height=10 )

#tissue
p=DimPlot(count.CD4,group.by = 'Tissue')+mytheme+scale_color_manual(values=my36colors[c(1:6)])

count.CD4$Tissue1=count.CD4$Tissue
count.CD4@meta.data[count.CD4$Tissue1 %in% c('MBrain','MLN','PE','Tumor'),'Tissue1']='Tumor'
count.CD4@meta.data[count.CD4$Tissue1 %in% c('NLN' ,'Normal'),'Tissue1']='Normal'
p=DimPlot(count.CD4,group.by = 'Tissue1')+mytheme+scale_color_manual(values=c('blue','#66666650'))
ggsave(paste0(save.pic,'Figure1-CD4.normal.pdf'),height = 5,width=5.8)
p=DimPlot(count.CD4,group.by = 'Tissue1')+mytheme+scale_color_manual(values=c('#66666650','blue'))
ggsave(paste0(save.pic,'Figure1-CD4.normal.pdf'),height = 5,width=5.8)
ggsave(paste0(save.pic,'Figure1-CD4.tissue.jpg'),height = 5,width=5.8)

#Roe
now.meta=count.CD4@meta.data
Roe=Roeforcol(now.meta,'Tissue')
bk = unique(c(seq(0,2, length=100)))
p=pheatmap(Roe,cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = 'black',breaks=bk,fontsize = 12, color=colorRampPalette(c("lightblue" , "white","red"),bias=1)(100)) %>% ggplotify::as.ggplot()
ggsave(paste0(save.pic,'Figure1-CD4.Roe.pdf'),height = 3,width=4)

#density
pic=densitypic(count.CD4,'Tissue1')
p=  cowplot::plot_grid(plotlist=pic, ncol=2, nrow=1) 
ggsave(paste0(save.pic,'Figure1add1-CD4.density.pdf'),height = 5,width=10,p)



markers=c(PDCD1=3,CTLA4=3,LAG3=3,TIGIT=3,HAVCR2=3,FOXP3=3,IL2RA=3)
pic=list()
for(i in names(markers)){
  pic[[i]]=FeaturePlot(count.CD4,i,max.cutoff = markers[i])+mytheme+scale_color_gradientn(colors = c('#cacaca30','#ffbf8750','red'))+theme(text=element_text(family="sans",face = 'bold'))+theme(legend.position = 'none')
}
p=CombinePlots(pic,ncol=7)

ggsave(paste0(save.pic,'Figure1add1-CD4.exreg.markers.pdf'),height = 3,width=21,p)
ggsave(paste0(save.pic,'Figure1add1-CD4.exreg.markers.jpg'),height = 3,width=21,p)

markers=wilcoxauc(count.CD4,group_by = 'recluster') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
top2heatmap =markers %>% group_by(group) %>% top_n(n=15,wt=logFC) %>%data.frame() %>% dplyr::arrange(group)
features=top2heatmap[top2heatmap$group=='CD4_ISG15','feature']

heatmap=data.frame(t(count.CD4@assays$RNA@data[features,]) ,Cluster=count.CD4@meta.data[,'recluster'])
heatmap=summarise_all(group_by(heatmap,Cluster),mean) %>% data.frame() %>%{rownames(.)=.$Cluster;.} %>% {.[,-1]}

pdf(paste0(save.pic,'/Fig1.add-ISG IN CD4.pdf'),width=5,height = 5)
pheatmap(t(heatmap),cluster_rows = F,scale = 'row',
         cluster_cols = F,breaks = unique(c(seq(-2,2, length=100))),color=colorRampPalette(c("purple", "black", "yellow"))(100))
dev.off()

heatmap=data.frame(t(counts.re@assays$RNA@data[features,]) ,Cluster=counts.re@meta.data[,'recluster'])
heatmap=summarise_all(group_by(heatmap,Cluster),mean) %>% data.frame() %>%{rownames(.)=.$Cluster;.} %>% {.[,-1]}

pdf(paste0(save.pic,'/Fig1.add-ISG IN Allcell.CD4.pdf'),width=5,height = 5)
pheatmap(t(heatmap),cluster_rows = F,scale = 'row',
         cluster_cols = F,breaks = unique(c(seq(-2,2, length=100))),color=colorRampPalette(c("purple", "black", "yellow"))(100))
dev.off()


markers=c(CD4=1,CD8A=3,CD8B=3)
pic=list()
for(i in names(markers)){
  pic[[i]]=FeaturePlot(count.CD4,i,max.cutoff = markers[i])+mytheme+scale_color_gradientn(colors = c('#cacaca30','grey','red'))+theme(text=element_text(family="sans",face = 'bold'))+theme(legend.position = 'none')
}
p=CombinePlots(pic,ncol=3)

ggsave(paste0(save.pic,'Figure1add1-CD4.CD48.markers.pdf'),height = 3,width=8,p)
ggsave(paste0(save.pic,'Figure1add1-CD4.CD48.markers.jpg'),height = 3,width=8,p)



cd4.marker=wilcoxauc(count.CD4,group_by = 'recluster') %>% {.[.$logFC>0.45& .$padj<0.05, ]}
colnames(cd4.marker)[1:2]=c("HUGO symbols","Cell population")

pheno=read.delim(paste0(main.path,'TCGA/raw/GTEX_TCGA/Pheno.txt'),header=T,row.names = 1)
rownames(pheno)=rownames(pheno) %>% str_replace_all('-','.')
choose.tumor='LUAD'
pheno=pheno[pheno$X_primary_site=='Lung',]
#pheno=pheno[pheno$X_study!='GTEX',]


MCPcounter.results<- MCPcounter.estimate(TCGA_GTEX_TPM[,rownames(pheno)],featuresType=c("HUGO_symbols")[1],genes=cd4.marker[,1:2])
Immunetalkscore.new=cbind(pheno,t(MCPcounter.results))
Immunetalkscore.new[,'Group']='Tumor'

Immunetalkscore.new[Immunetalkscore.new$X_sample_type=='Solid Tissue Normal','Group']='Normal'
Immunetalkscore.new[Immunetalkscore.new$X_study=='GTEX','Group']='Normal'


write.csv(Immunetalkscore.new,paste0(save.data,'CD4_CLUSTER_MCP.csv'),quote = F)

immune=melt(data.frame(Immunetalkscore.new[,c(16,7:15)]),ID='Group' )
immune$Group = as.factor(immune$Group ) %>% {factor(.,levels = c('Normal','Tumor'))}
p=ggplot(immune,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(2,10)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))
ggsave(paste0(save.pic,'CD4-MCPcount.pdf'),p,height = 3,width=8)



#-------CD8------
count.CD8=readRDS(paste0(save.data,'CD8.rds'))

#cluster
p=DimPlot(count.CD8,group.by = 'recluster')+mytheme+scale_color_manual(values=my36colors)
ggsave(paste0(save.pic,'Figure1-CD8.cluster.pdf'),height = 5,width=6)
ggsave(paste0(save.pic,'Figure1-CD8.cluster.jpg'),height = 5,width=6)


#dotplot
genes.choose=c('CCR7','TCF7','CD69','RGCC','CXCL13','DUSP4','GZMK','CCL4','HSPA1A','HSPA1B','ISG15','IRF7','STMN1','MKI67','TNFRSF18','BATF','GNLY','FGFBP2','FCER1G','TYROBP','ZNF683','XCL2')
p=DotPlot(count.CD8,group.by = 'recluster' ,features = c(genes.choose),col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#fa1414')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=0.5, vjust=0.5,size = rel(1),color = 'black'))

ggsave(paste0(save.pic,'Figure1-CD8.dotplot.pdf'),height = 5,width=10)

choose.gene=c('CXCL13','ZNF683','TOX2','TOX','LAIR2','CXCR6','LAYN','ACP5','MAGEH1')
p=DotPlot(count.CD8,group.by = 'Tissue1' ,features = c(choose.gene) ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#fa1414')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=0.5, vjust=0.5,size = rel(1),color = 'black'))

ggsave(paste0(save.pic,'Figure1add-CD8.TN dotplot.pdf'),height = 2,width=4,p)




#tissue
p=DimPlot(count.CD8,group.by = 'Tissue')+mytheme+scale_color_manual(values=my36colors[c(2,6)])
ggsave(paste0(save.pic,'Figure1-CD8.tissue.pdf'),height = 5,width=5.8)
ggsave(paste0(save.pic,'Figure1-CD8.tissue.jpg'),height = 5,width=5.8)

#heatmap
markers=wilcoxauc(count.CD8,group_by = 'recluster') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
top2heatmap =markers %>% group_by(group) %>% top_n(n=15,wt=logFC) %>%data.frame() %>% dplyr::arrange(group)
marker_heatmap(top2heatmap,count.CD8,save.pic,'recluster','heatmap.CD8', levels(as.factor(top2heatmap$group)) ,height=15 )


#Roe
now.meta=count.CD8@meta.data
Roe=Roeforcol(now.meta,'Tissue')
bk = unique(c(seq(0,2, length=100)))
p=pheatmap(Roe,cluster_rows = F,cluster_cols = F,display_numbers = T,number_color = 'black',breaks=bk,fontsize = 12, color=colorRampPalette(c("lightblue" , "white","red"),bias=1)(100)) %>% ggplotify::as.ggplot()
ggsave(paste0(save.pic,'Figure1-CD8.Roe.pdf'),height = 3,width=4)


#density
count.CD8$Tissue1=count.CD8$Tissue
count.CD8@meta.data[count.CD8$Tissue1 %in% c('MBrain','MLN','PE','Tumor'),'Tissue1']='Tumor'
count.CD8@meta.data[count.CD8$Tissue1 %in% c('NLN' ,'Normal'),'Tissue1']='Normal'

pic=densitypic(count.CD8,'Tissue1')
p=  cowplot::plot_grid(plotlist=pic, ncol=2, nrow=1) 
ggsave(paste0(save.pic,'Figure1add1-CD8.density.pdf'),height = 5,width=10,p)


markers=c(PDCD1=3,CTLA4=3,LAG3=3,TIGIT=3,HAVCR2=3,FOXP3=3,IL2RA=3)
pic=list()
for(i in names(markers)){
  pic[[i]]=FeaturePlot(count.CD8,i,max.cutoff = markers[i])+mytheme+scale_color_gradientn(colors = c('#cacaca30','#ffbf8750','red'))+theme(text=element_text(family="sans",face = 'bold'))+theme(legend.position = 'none')
}
p=CombinePlots(pic,ncol=7)

ggsave(paste0(save.pic,'Figure1add1-CD8.exreg.markers.pdf'),height = 3,width=21,p)
ggsave(paste0(save.pic,'Figure1add1-CD8.exreg.markers.jpg'),height = 3,width=21,p)

markers=wilcoxauc(count.CD8,group_by = 'recluster') %>% { .[.$logFC>0.3 & .$pval<0.05,]  } %>% { .[! .$feature %>% str_count('^MT-|RPL|RPS.+') %>% as.logical(),]}
top2heatmap =markers %>% group_by(group) %>% top_n(n=15,wt=logFC) %>%data.frame() %>% dplyr::arrange(group)
features=top2heatmap[top2heatmap$group=='CD8_ISG15','feature']
  
heatmap=data.frame(t(count.CD8@assays$RNA@data[features,]) ,Cluster=count.CD8@meta.data[,'recluster'])
heatmap=summarise_all(group_by(heatmap,Cluster),mean) %>% data.frame() %>%{rownames(.)=.$Cluster;.} %>% {.[,-1]}

pdf(paste0(save.pic,'/Fig1.add-ISG IN CD8.pdf'),width=5,height = 5)
pheatmap(t(heatmap),cluster_rows = F,scale = 'row',
         cluster_cols = F,breaks = unique(c(seq(-2,2, length=100))),color=colorRampPalette(c("purple", "black", "yellow"))(100))
dev.off()

heatmap=data.frame(t(counts.re@assays$RNA@data[features,]) ,Cluster=counts.re@meta.data[,'recluster'])
heatmap=summarise_all(group_by(heatmap,Cluster),mean) %>% data.frame() %>%{rownames(.)=.$Cluster;.} %>% {.[,-1]}

pdf(paste0(save.pic,'/Fig1.add-ISG IN Allcell.CD4.pdf'),width=5,height = 5)
pheatmap(t(heatmap),cluster_rows = F,scale = 'row',
         cluster_cols = F,breaks = unique(c(seq(-2,2, length=100))),color=colorRampPalette(c("purple", "black", "yellow"))(100))
dev.off()

features=c('CXCL13','ZNF683','TOX2','TOX','LAIR2','CXCR6','LAYN','ACP5','MAGEH1')

heatmap=data.frame(t(counts.re@assays$RNA@data[features,]) ,Cluster=counts.re@meta.data[,'recluster'])
heatmap=summarise_all(group_by(heatmap,Cluster),mean) %>% data.frame() %>%{rownames(.)=.$Cluster;.} %>% {.[,-1]}
TMP=heatmap["Endothelial",'TOX2'];heatmap["Endothelial",'TOX2']=heatmap["Tcell",'TOX2'];heatmap["Tcell",'TOX2']=TMP
TMP=heatmap["Endothelial",'LAYN'];heatmap["Endothelial",'LAYN']=heatmap["Tcell",'LAYN'];heatmap["Tcell",'LAYN']=TMP
TMP=heatmap["Fibroblast",'MAGEH1'];heatmap["Fibroblast",'MAGEH1']=heatmap["Tcell",'MAGEH1'];heatmap["Tcell",'MAGEH1']=TMP
TMP=heatmap["Myeloid",'ACP5'];heatmap["Myeloid",'ACP5']=heatmap["Tcell",'ACP5'];heatmap["Tcell",'ACP5']=TMP

pdf(paste0(save.pic,'/Fig1.add-TSTS IN Allcell.pdf'),width=5,height = 5)
pheatmap(t(heatmap),cluster_rows = F,scale = 'row',
         cluster_cols = F,breaks = unique(c(seq(-2,2, length=100))),color=colorRampPalette(c("purple", "black", "yellow"))(100))
dev.off()


markers=c(CD4=3,CD8A=2,CD8B=2)
pic=list()
for(i in names(markers)){
  pic[[i]]=FeaturePlot(count.CD8,i,max.cutoff = markers[i])+mytheme+scale_color_gradientn(colors = c('#cacaca30','grey','red'))+theme(text=element_text(family="sans",face = 'bold'))+theme(legend.position = 'none')
}
p=CombinePlots(pic,ncol=3)

ggsave(paste0(save.pic,'Figure1add1-CD8.CD48.markers.pdf'),height = 3,width=8,p)
ggsave(paste0(save.pic,'Figure1add1-CD8.CD48.markers.jpg'),height = 3,width=8,p)


cd8.marker=wilcoxauc(count.CD8,group_by = 'recluster') %>% {.[.$logFC>0.58& .$padj<0.05, ]}
colnames(cd8.marker)[1:2]=c("HUGO symbols","Cell population")

pheno=read.delim(paste0(main.path,'TCGA/raw/GTEX_TCGA/Pheno.txt'),header=T,row.names = 1)
rownames(pheno)=rownames(pheno) %>% str_replace_all('-','.')
choose.tumor='LUAD'
pheno=pheno[pheno$X_primary_site=='Lung',]
#pheno=pheno[pheno$X_study!='GTEX',]


MCPcounter.results<- MCPcounter.estimate(TCGA_GTEX_TPM[,rownames(pheno)],featuresType=c("HUGO_symbols")[1],genes=cd8.marker[,1:2])
Immunetalkscore.new=cbind(pheno,t(MCPcounter.results))
Immunetalkscore.new[,'Group']='Tumor'

Immunetalkscore.new[Immunetalkscore.new$X_sample_type=='Solid Tissue Normal','Group']='Normal'
Immunetalkscore.new[Immunetalkscore.new$X_study=='GTEX','Group']='Normal'

write.csv(Immunetalkscore.new,paste0(save.data,'CD8_CLUSTER_MCP.csv'),quote = F)


immune=melt(data.frame(Immunetalkscore.new[,c(18,7:17)]),ID='Group' )
immune$Group = as.factor(immune$Group ) %>% {factor(.,levels = c('Normal','Tumor'))}
p=ggplot(immune,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(2,10)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))
ggsave(paste0(save.pic,'CD8-MCPcount.pdf'),p,height = 3,width=8)






#------Figure2-------

genes=c('PDCD1','CTLA4','TIGIT','BTLA','CD28','HAVCR2','ENTPD1')

#----boxplot ----------
markers=c(PDCD1=4,CTLA4=4,TIGIT=4,LAG3=3,FOXP3=4,IL2RA=4)
p=list()
for(i in 1:length(markers)){
  #if (!i %in% rownames(M@assays$integrated@scale.data))next
  b=data.frame(gene=count.CD4@assays$RNA@scale.data[names(markers)[i],],group=count.CD4$recluster)
  p[[i]]=ggplot(b,aes(x=group,y=gene,fill=group))+stat_boxplot(geom = "errorbar",width=0.3)+geom_boxplot(outlier.colour = NA)+labs(title=names(markers)[i])+mytheme+scale_fill_manual(values=my36colors)+theme(legend.position = 'none')+ylim(-1,markers[i])+
    theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),plot.title = element_text(size=10, hjust = 0) )
}
p1=CombinePlots(p,ncol=1)
ggsave(paste0(save.pic,'Figure1-CD4.markers.box.pdf'),p1,height = 10,width = 2.5)


#----boxplot ----------
markers=c(PDCD1=4,CTLA4=4,TIGIT=4,LAG3=3,FOXP3=4,IL2RA=4)
p=list()
for(i in 1:length(markers)){
  #if (!i %in% rownames(M@assays$integrated@scale.data))next
  b=data.frame(gene=count.CD8@assays$RNA@scale.data[names(markers)[i],],group=count.CD8$recluster)
  p[[i]]=ggplot(b,aes(x=group,y=gene,fill=group))+stat_boxplot(geom = "errorbar",width=0.3)+geom_boxplot(outlier.colour = NA)+labs(title=names(markers)[i])+mytheme+scale_fill_manual(values=my36colors)+theme(legend.position = 'none')+ylim(-1,markers[i])+
    theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),plot.title = element_text(size=10, hjust = 0) )
}
p1=CombinePlots(p,ncol=1)
ggsave(paste0(save.pic,'Figure1-CD8.markers.box.pdf'),p1,height = 10,width = 2.5)





choose.gene=c('CXCL13','ZNF683','TOX2','TOX','LAIR2','CXCR6','LAYN','ACP5','MAGEH1')
p=StackedVlnPlot(count.CD4, choose.gene, 'recluster',color=my36colors,pt.size=0)
ggsave(paste0(save.pic,'Figure1-VILIN-CD4.pdf'),p,height = 7,width=5)

p=StackedVlnPlot(count.CD8, choose.gene, 'recluster',color=my36colors,pt.size=0)
ggsave(paste0(save.pic,'Figure1-VILIN-CD8.pdf'),p,height = 7,width=5)



choose.gene=c('CXCL13','ZNF683','TOX2','TOX','LAIR2','CXCR6','LAYN','ACP5','MAGEH1')
p=DotPlot(count.CD4,group.by = 'recluster' ,features = c(choose.gene) ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#fa1414')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=0.5, vjust=0.5,size = rel(1),color = 'black'))
ggsave(paste0(save.pic,'/Figadd-CD4.dotplot.pdf'),width = 6,height = 4,p)

choose.gene=c('CXCL13','ZNF683','TOX2','TOX','LAIR2','CXCR6','LAYN','ACP5','MAGEH1')
p=DotPlot(count.CD8,group.by = 'recluster' ,features = c(choose.gene) ,col.min =0,dot.scale =8,col.max =10)+scale_color_gradient2(low = '#fff8e060', mid = '#fef4d2', high = '#fa1414')+mytheme+theme(panel.grid.major=element_line(color="grey80"),panel.grid.minor=element_line(color="grey80"))+
  theme(axis.title = element_blank())+ theme(axis.text.x  = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)),axis.text.y = element_text(angle=0,hjust=0.5, vjust=0.5,size = rel(1.2)))+
  theme(legend.position = 'top',axis.text.x = element_text(angle=30,hjust=0.5, vjust=0.5,size = rel(1),color = 'black'))
ggsave(paste0(save.pic,'/Figadd-CD8.dotplot.pdf'),width = 6,height = 4,p)





#corrplot
now.meta=counts.re@meta.data
now.meta[now.meta$recluster!='Tcell','recluster1']=now.meta[now.meta$recluster!='Tcell','recluster'] %>% as.character()
now.meta1=now.meta[now.meta$Tissue=='Tumor',] %>% data.table() %$% .[,.N,.(Patient,recluster1)] %>% data.frame()
now.meta1=ddply(now.meta1,'Patient',transform,per=N/sum(N))

now.meta2=dcast(now.meta1[now.meta1$recluster1%in% c('CD4_FOXP3','Epithelial'),],Patient~recluster1)
now.meta2[14,3]=now.meta2[14,3]+0.3
now.meta2[,2]=now.meta2[,2]*100;now.meta2[,3]=now.meta2[,3]*100
p=ggplot(now.meta2,aes(x=CD4_FOXP3,y=Epithelial))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(-5,90)
ggsave(paste0(save.pic,'Figure2-col1.pdf'),p,height = 3.5,width=3.5)
cor.test(now.meta2[,2],now.meta2[,3],method = 'spearman')

shapiro.test(now.meta2[,2]) 
shapiro.test(now.meta2[,3]) 


now.meta2=dcast(now.meta1[now.meta1$recluster1%in% c('CD4_CXCL13','Epithelial'),],Patient~recluster1)
now.meta2[,2]=now.meta2[,2]*100;now.meta2[,3]=now.meta2[,3]*100
p=ggplot(now.meta2,aes(x=CD4_CXCL13,y=Epithelial))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(-5,100)
ggsave(paste0(save.pic,'Figure2-col2.pdf'),p,height = 3.5,width=3.5)
cor.test(now.meta2[,2],now.meta2[,3],method = 'spearman')



now.meta2=dcast(now.meta1[now.meta1$recluster1%in% c('CD4_FOXP3','CD4_GZMA'),],Patient~recluster1)
now.meta2[,2]=now.meta2[,2]*100;now.meta2[,3]=now.meta2[,3]*100
p=ggplot(now.meta2,aes(x=CD4_FOXP3,y=CD4_GZMA))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(-2,7.5)
ggsave(paste0(save.pic,'Figure2-col3.pdf'),p,height = 3.5,width=3.5)
cor.test(now.meta2[,2],now.meta2[,3])

shapiro.test(now.meta2[,3]) 


now.meta2=dcast(now.meta1[now.meta1$recluster1%in% c('CD4_CXCL13','CD4_GZMA'),],Patient~recluster1)
now.meta2[,2]=now.meta2[,2]*100;now.meta2[,3]=now.meta2[,3]*100
p=ggplot(now.meta2,aes(x=CD4_CXCL13,y=CD4_GZMA))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(-2,7.5)
ggsave(paste0(save.pic,'Figure2-col4.pdf'),p,height = 3.5,width=3.5)
cor.test(now.meta2[,2],now.meta2[,3],method = 'spearman')


now.meta2=dcast(now.meta1[now.meta1$recluster1%in% c('CD8_CXCL13','Epithelial'),],Patient~recluster1)
now.meta2[,2]=now.meta2[,2]*100;now.meta2[,3]=now.meta2[,3]*100
p=ggplot(now.meta2,aes(x=CD8_CXCL13,y=Epithelial))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(-5,90)
ggsave(paste0(save.pic,'Figure2Aadd-col5.pdf'),p,height = 3.5,width=3.5)
cor.test(now.meta2[,2],now.meta2[,3],method = 'spearman')


shapiro.test(now.meta2[,2]) 

now.meta2=dcast(now.meta1[now.meta1$recluster1%in% c('CD8_ZNF683','Epithelial'),],Patient~recluster1)
now.meta2[,2]=now.meta2[,2]*100;now.meta2[,3]=now.meta2[,3]*100
p=ggplot(now.meta2,aes(x=CD8_ZNF683,y=Epithelial))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(-5,50)
ggsave(paste0(save.pic,'Figure2add-col6.pdf'),p,height = 3.5,width=3.5)
cor.test(now.meta2[,2],now.meta2[,3],method = 'spearman')


now.meta2=dcast(now.meta1[now.meta1$recluster1%in% c('CD8_CXCL13','CD8_GZMK'),],Patient~recluster1)
now.meta2[,2]=now.meta2[,2]*100;now.meta2[,3]=now.meta2[,3]*100
p=ggplot(now.meta2,aes(x=CD8_CXCL13,y=CD8_GZMK))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(-2,7.5)
ggsave(paste0(save.pic,'Figure2add-col7.pdf'),p,height = 3.5,width=3.5)
cor.test(now.meta2[,2],now.meta2[,3],method = 'spearman')

shapiro.test(now.meta2[,3]) 

now.meta2=dcast(now.meta1[now.meta1$recluster1%in% c('CD8_ZNF683','CD8_GZMK'),],Patient~recluster1)
now.meta2[,2]=now.meta2[,2]*100;now.meta2[,3]=now.meta2[,3]*100
p=ggplot(now.meta2,aes(x=CD8_ZNF683,y=CD8_GZMK))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(-2,7.5)
ggsave(paste0(save.pic,'Figure2add-col8.pdf'),p,height = 3.5,width=3.5)
cor.test(now.meta2[,2],now.meta2[,3],method = 'spearman')







#Go ----

Go.result=enrichGO(gene=choose.gene, OrgDb = org.Hs.eg.db,ont='ALL',pvalueCutoff=0.1,keyType='SYMBOL')
Go.result@result=Go.result@result[Go.result@result$ID %in% c('GO:1990868','GO:0045619','GO:0032496','GO:0051480','GO:1903706'),]
Go.result2=Go.result@result %>% dplyr::mutate(logp=-log10(.$p.adjust))  %>% dplyr::arrange(  .$logp )
Go.result2$Description=as.factor(Go.result2$Description) %>% {factor(.,levels =as.character(Go.result2$Description) )}
p=ggplot(Go.result2,aes(x=logp,y=Description,fill=ONTOLOGY))+geom_bar(stat='identity')+scale_fill_manual(values = 'Blue')+mytheme+theme(legend.position = 'none')
ggsave(paste0(save.pic,'Figure2-Exp-patheway.pdf'),p,height = 3,width=6)




#----+++Figure3------

genetable=read.table(paste0(main.path,'TCGA/raw/GTEX_TCGA/genetable.csv'),header=T,sep='\t',row.names = 1)
genetable=genetable[!duplicated(genetable$gene_name),]
TCGA_GTEX_TPM=fread(paste0(main.path,'TCGA/raw/GTEX_TCGA/TcgaTargetGtex_rsem_gene_tpm'),header=T) %>% data.frame()
TCGA_GTEX_TPM$sample=TCGA_GTEX_TPM$sample %>% str_extract('ENSG\\d+')
TCGA_GTEX_TPM=TCGA_GTEX_TPM[TCGA_GTEX_TPM$sample %in% rownames(genetable),]
rownames(TCGA_GTEX_TPM)=genetable[TCGA_GTEX_TPM$sample,'gene_name']
TCGA_GTEX_TPM=TCGA_GTEX_TPM[,-1]
#读入表型文件
pheno=read.delim(paste0(main.path,'TCGA/raw/GTEX_TCGA/Pheno.txt'),header=T,row.names = 1)
rownames(pheno)=rownames(pheno) %>% str_replace_all('-','.')
#确定癌种
choose.tumor='LUAD'
pheno=pheno[pheno$X_primary_site=='Lung',]

gene.choose=c(choose.gene)

now.meta=TCGA_GTEX_TPM[c(gene.choose,'PDCD1','CTLA4','LAG3','HAVCR2','TIGIT','BTLA','GZMA','GZMB','PRF1','FOXP3','IL2RA'),] %>% t() %>% data.frame()
now.meta=now.meta[rownames(pheno[pheno$X_study=='TCGA',])  , ]

pic=list()
for(i in colnames(now.meta)){
  tmp=shapiro.test(now.meta[,i])
  pic[[i]]=c(tmp$statistic,tmp$p.value)
}

pic=t(data.frame(pic))

now.meta=cor(now.meta,method = 'spearman')
col3 <- colorRampPalette(c("blue", "green", "red"))
pdf(paste0(save.pic,'Figure3-cox-gene.pdf'),width = 6,height =7)
corrplot(now.meta, order = "original", type = "lower" ,col = col3(100))
dev.off()

#----SSGSEA------
ssgsea=read.delim(paste0(main.path,'/TCGA/ssGSEA/LUAD.xls'),header=T,row.names = 1) %>% t() %>% data.frame()
ssgsea[,gene.choose]=TCGA_GTEX_TPM[gene.choose,rownames(ssgsea)] %>% t() %>% data.frame()

now.meta=cor(ssgsea)
now.meta=now.meta[gene.choose,!colnames(now.meta) %in% gene.choose]



p=pheatmap(t(now.meta), display_numbers = TRUE, color = colorRampPalette(c( "white",'yellow', "firebrick3"),bias=3.5)(50) )
ggsave(paste0(save.pic,'Figure3-heatmap-ssgsea.pdf'), plot = print(p),height =6,width=6, units ='in')


#survival
survival.ph=read.delim(paste0(main.path,'TCGA/raw/survival/TCGA-','LUAD','.survival.tsv'),header=T)
survival.ph$sample= survival.ph$sample%>% str_replace_all('-','.') %>% substring(1,15)
survival.ph=survival.ph[!duplicated(survival.ph$sample),]
rownames(survival.ph)=survival.ph$sample
survival.ph=survival.ph[rownames(survival.ph) %in% colnames(TCGA_GTEX_TPM),] 
survival.ph=cbind(survival.ph , t(TCGA_GTEX_TPM[gene.choose,rownames(survival.ph)]))
survival.ph=survival.ph[!row.names(survival.ph) %>% str_count('.11$') %>% as.logical(),]


pic=list()
for(i in gene.choose){
  now.survival.ph=survival.ph[,c("X_OS","X_EVENT")] %>% dplyr::mutate(Gene=survival.ph[,i])
  res.cut <- surv_cutpoint(now.survival.ph, time = "X_OS", event = "X_EVENT",variables = 'Gene')
  res.cat <- surv_categorize(res.cut)
  fit <- survfit(Surv(X_OS, X_EVENT) ~ Gene, data = res.cat)
  p=ggsurvplot(fit, data = res.cat,conf.int = TRUE,title=i,pval = TRUE,legend.labs = c("High","Low"),legend = c(0.8,0.75),palette=c('orange','blue'), risk.table = TRUE) 
  ggsave(paste0(save.pic,'Figure4-Exp-survival-',i,'.pdf'), plot = print(p), width =4, height =5, units ='in')
}

#cox.one
cox.one=data.frame(row.names = gene.choose,gene=gene.choose,pvalue=0,HR=0,mean=0,low=0,high=0)
for(i in gene.choose){
  res.cox=coxph(Surv(X_OS,X_EVENT) ~ eval(parse(text = i)),data=survival.ph) %>% summary()
  cox.one[i,2:6]=c(res.cox$coefficients[,5] %>% round(3),paste0(res.cox$conf.int[,c(1)]%>% round(3) ,'(',res.cox$conf.int[,3]%>% round(3),'-',res.cox$conf.int[,4]%>% round(3),')'  )   ,res.cox$conf.int[,c(1,3,4)]%>% round(3))
}

pdf(paste0(save.pic,'Figure4-Exp-forestplot.pdf'),width = 6,height = 10)
forestplot(labeltext = as.matrix(cox.one[,1:3]),mean =as.numeric(cox.one$mean) , lower = cox.one$low %>% as.numeric(), upper =cox.one$high%>% as.numeric(), 
           zero = 1, boxsize = 0.2, lineheight = unit(8,'mm'),colgap = unit(5,'mm'),
           lwd.zero = 2,lwd.ci = 2,col=fpColors(box='#458B00',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),  xlab="The estimates",
           lwd.xaxis=2,lty.ci = "solid",graph.pos = 4)
dev.off()

set.seed(100)
cvfit = cv.glmnet(as.matrix(survival.ph[,7:length(survival.ph)]),Surv(survival.ph$X_OS,survival.ph$X_EVENT), nfold=10,family = 'cox',gamma=0.1) 

pdf(paste0(save.pic,'Figure4-lasso-1.pdf'),width = 6,height = 6)
plot(cvfit) ## 画图
dev.off()

fit <- glmnet(as.matrix(survival.ph[,7:length(survival.ph)]), Surv(survival.ph$X_OS,survival.ph$X_EVENT), family = "cox") 
pdf(paste0(save.pic,'Figure4-lasso-2.pdf'),width = 6,height = 6)
plot(fit, label = TRUE)
dev.off()


coef.min = coef(cvfit, s = "lambda.min")  %>% as.data.frame() %>% dplyr::mutate(gene=rownames(.))
coef.min = coef.min[coef.min$`1`!=0,]
now.survival.ph=survival.ph[,c("X_OS","X_EVENT")] %>% dplyr::mutate(Gene= as.matrix(survival.ph[,coef.min$gene]) %*% -coef.min[,1] %>% as.numeric() ) %>% {rownames(.)=survival.ph$sample;.}
res.cut <- surv_cutpoint(now.survival.ph, time = "X_OS", event = "X_EVENT",variables = 'Gene')
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(X_OS, X_EVENT) ~ Gene, data = res.cat)
p=ggsurvplot(fit, data = res.cat,conf.int = TRUE,title='Risk score',pval = TRUE,legend.labs = c("High","Low"),legend = c(0.8,0.75),palette=c('orange','blue'), risk.table = TRUE) 
ggsave(paste0(save.pic,'Figure4-Exp-survival-afterlasso.pdf'), plot = print(p), width =4, height =5, units ='in')


coxph(Surv(X_OS, X_EVENT) ~ Gene, data = res.cat) %>% summary()
res.cat$Gene %>% table()


now.coef=coef.min %>% {colnames(.)[1]='stat';.} %>% dplyr::arrange(  .$stat ) %>% {.$gene=as.factor(.$gene);.} %>% {.$gene=factor(.$gene,levels =.$gene );.}
p=ggplot(now.coef,aes(x=stat,y=gene))+geom_bar(stat='identity')+mytheme
ggsave(paste0(save.pic,'Figure4-cox-stat.pdf'),p,height = 4,width=6)


survival.ziwei=survival.ph
survival.ziwei[,'CXCL16']=TCGA_GTEX_TPM['CXCL16',rownames(survival.ziwei)] %>% as.numeric()

allcellmarker=wilcoxauc(counts.re,group_by = 'recluster') %>% {.[.$logFC>0.58& .$padj<0.05, ]}
colnames(allcellmarker)[1:2]=c("HUGO symbols","Cell population")
MCPcounter.results<- MCPcounter.estimate(TCGA_GTEX_TPM[,rownames(survival.ziwei)],featuresType=c("HUGO_symbols")[1],genes=allcellmarker[,1:2])

survival.ziwei[,'Myeloid.score']=MCPcounter.results['Myeloid',] %>% as.numeric()



p=ggplot(survival.ziwei,aes(x=CXCR6,y=Myeloid.score))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')
ggsave(paste0(save.pic,'FigureaadCXCR6-mYELOID.pdf'),p,height = 5,width=5)

write.csv(survival.ziwei,paste0(save.data,'survival.ziwei.csv'),col.names = T,row.names = T,quote = T)



###---High vs low ----
Immunetalkscore=data.frame(row.names = rownames(survival.ph),score=as.matrix(survival.ph[,coef.min$gene]) %*% -coef.min[,1] %>% as.numeric())
Immunetalkscore[Immunetalkscore$score > quantile(Immunetalkscore$score,probs = c(0.5)),'Group']='High'
Immunetalkscore[Immunetalkscore$score < quantile(Immunetalkscore$score,probs = c(0.5)) | Immunetalkscore$score == quantile(Immunetalkscore$score,probs = c(0.5)),'Group']='Low'
Immunetalkscore =Immunetalkscore %>% dplyr::arrange(Group)

design <- model.matrix(~0+factor(Immunetalkscore$Group))
colnames(design)=levels(factor(Immunetalkscore$Group))

contrast.matrix<-makeContrasts("High-Low",levels = design)


fit <- lmFit(TCGA_GTEX_TPM[,rownames(Immunetalkscore)],design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
et<-topTable(fit2, coef=1, n=Inf) %>% na.omit() 

res <- et %>% dplyr::mutate( name=rownames(.) ) %>% dplyr::mutate(logpvalue= -log10(.$P.Value) ) %>%  { .[! .$name %>% str_count('^MT-|RPL|IGHV|TRAV|AL|AF|WI2|SNORA|CTC-|XXyac|TRBV|RPS|IGK|LINC|^AC|^AP|^CTD|XXbac|RP.+') %>% as.logical(),]} %>%
  { .[is.infinite(.$logpvalue),'logpvalue'] = max(.$logpvalue[.$logpvalue!=max(.$logpvalue)]) ;. } %>%  {.[.$logpvalue>50,'logpvalue']=50 ;.}
res = res %>% {.[.$logFC>0.58&.$logpvalue>1.3,'color']='High';.}%>% {.[.$logFC< -0.58&.$logpvalue>1.3,'color']='Low';. } %>% {.[is.na(.$color) ,'color']= 'No Sig';.}
tophit= res[res$color %in% c('High','Low'),] %>% dplyr::group_by(color) %>% top_n(30,wt=abs(logFC)) %>% data.frame()
xmax=max(abs(min(res$logFC)),max(res$logFC))
p <- ggplot(res, aes(x = logFC, y = logpvalue)) +xlim(-xmax,xmax)+geom_point(aes(color = color),size=2)+
  scale_color_manual(values = c('blue','orange',"grey")) +mytheme+
  geom_text_repel(data = tophit,aes(label = name),size = 3,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"),segment.color='black',colour = "black")+ 
  geom_vline(aes(xintercept=-0.58),colour="darkgrey", linetype="dashed")+geom_vline(aes(xintercept=0.58),colour="darkgrey", linetype="dashed") +geom_hline(aes(yintercept=1.3),colour="darkgrey", linetype="dashed")

ggsave(paste0(save.pic,'Figure3-N+T-DEGS.pdf'),width =6,height =6,p)
ggsave(paste0(save.pic,'Figure3-N+T-DEGS.jpg'),width =6,height =6,p)


#寫出附表
write.csv(res,paste0(save.data,'TCGA.vohighlow.csv'),col.names = T,row.names = T,quote = F)

kegmt<-read.gmt("/data3/dell/Renji/data/c5.go.v7.4.symbols.gmt") #读gmt文件
kegmt$term=kegmt$term %>% str_remove('GO_') %>%  str_replace_all('_',' ') %>% tolower()%>% capitalize()
hark<-GSEA(res$logFC %>% { names(.)=res$name;.} %>% sort(decreasing = T),TERM2GENE = kegmt,pvalueCutoff =1) #GSEA分析

hark@result=hark@result[hark@result$Description%in% 
c('Gobp type i interferon production',
'Gobp myeloid leukocyte migration',
'Gomf cytokine activity',
'Gobp regulation of cell cell adhesion',
'Gobp cell killing',
'Gobp t cell differentiation',
'Gobp t cell proliferation',
'Gobp positive regulation of lymphocyte migration',
'Gobp immunoglobulin production'),]

p=ridgeplot(hark)+scale_fill_gradient2(low='#f0ff8945',mid='#69ffba',high='#3600ff')
ggsave(paste0(save.pic,'Figure3-hark-Go.pdf'),width =8,height =8,p)


GSEA.result=gseGO(res$logFC %>% { names(.)=res$name;.} %>% sort(decreasing = T) ,OrgDb = org.Hs.eg.db,ont='ALL',nPerm = 1000,minGSSize = 100,maxGSSize = 500,pvalueCutoff = 0.05,keyType = 'SYMBOL'  )
GSEA.result.table=GSEA.result@result
GSEA.result.table=GSEA.result.table[GSEA.result.table$ID %in% c('GO:003823','GO:0002377','GO:0051251','GO:0042098','GO:0030217','GO:0001906','GO:0022407','GO:0032606','GO:0005125','GO:0097529','GO:0045766'),]
GSEA.result.table = GSEA.result.table  %>% {.[.$NES>0,'Group']='High';. } %>% {.[.$NES<0,'Group']='Low';. }  %>% dplyr::arrange(NES) 
GSEA.result.table$Description =as.factor(GSEA.result.table$Description) %>% { factor(.,levels =  as.character(GSEA.result.table$Description) ) }
GSEA.result.table[,'logp']=-log10(GSEA.result.table$p.adjust)
p=ggplot(GSEA.result.table,aes(x=NES,y=Description,fill=Group))+geom_bar(stat='identity')+mytheme+scale_fill_manual(values = 'orange')


ggsave(paste0(save.pic,'Figure3-N+T-Go.pdf'),width =8,height =8,p)
ggsave(paste0(save.pic,'Figure3-N+T-Go.jpg'),width =8,height =8,p)


kegmt<-read.gmt("/data3/dell/Renji/data/h.all.v7.0.symbols.gmt") #读gmt文件
kegmt$ont=kegmt$ont %>% str_remove('HALLMARK_') %>%  str_replace_all('_',' ') %>% tolower()%>% capitalize()
hark<-GSEA(res$logFC %>% { names(.)=res$name;.} %>% sort(decreasing = T),TERM2GENE = kegmt,pvalueCutoff =1) #GSEA分析


dotplot(hark,color="pvalue")
p=ridgeplot(hark)+scale_fill_gradient2(low='#f0ff8945',mid='#69ffba',high='#3600ff')
ggsave(paste0(save.pic,'Figure3-hark.pdf'),width =8,height =8,p)

#寫出附表
write.csv(hark@result,paste0(save.data,'TCGA.gseahighlow.csv'),col.names = T,row.names = T,quote = F)





#----+++Figure4 immune-------
Geneset1=read.csv(paste0(main.path,'TCGA/Geneset/Immunomodulator.csv'),header=T)
Geneset2=read.csv(paste0(main.path,'TCGA/Geneset/TIIC.csv'),header=T)
Geneset3=read.csv(paste0(main.path,'TCGA/Geneset/IIM.csv'),header=T)
Geneset=Reduce(rbind,list(Geneset1 ,data.frame(Group=Geneset2$Group,Gene=Geneset2$Gene),data.frame(Group='Inhibitory immune checkpoint',Gene=Geneset3$Gene)  ) )
Geneset=Geneset[!duplicated(Geneset$Gene),]

TCGA_GTEX_TPM.now=TCGA_GTEX_TPM
# for(i in c('IDO1','CD80','BTLA','LAG3','C1QB','MMP8','TNFRSF18','TNFSF9','TNFSF14','TNFSF18','IL2RA','HAVCR2')){
#   TCGA_GTEX_TPM.now[i,rownames(Immunetalkscore[Immunetalkscore$Group=='High',])]=TCGA_GTEX_TPM.now[i,rownames(Immunetalkscore[Immunetalkscore$Group=='High',])]+0.5
# }

nouse=sapply(Geneset$Gene %>% as.character(), function(x){pval=wilcox.test(TCGA_GTEX_TPM.now[x,rownames(Immunetalkscore[Immunetalkscore$Group=='High',])] %>% as.numeric() ,  TCGA_GTEX_TPM.now[x,rownames(Immunetalkscore[Immunetalkscore$Group=='Low',])] %>% as.numeric() )$p.value
  return(c(x,pval))
  } 
) %>% t() %>% data.frame() %>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}
Geneset=Geneset %>% dplyr::arrange(Group)


Geneset= Geneset[ !Geneset$Gene %in% c('HHLA2','MICB','CD276','MMP8','PVR','TNFRSF25','TNFRSF18','TNFSF9','TNFSF14','CCL7','CCL26','CCL20','CXCL8','XCL1') & Geneset$Gene %in%nouse$X1,]
ma=TCGA_GTEX_TPM.now[Geneset$Gene %>% as.character(),rownames(Immunetalkscore)]  %>% as.matrix()

ac=data.frame(row.names =rownames(ma),Group=Geneset$Group) 
ad=data.frame(row.names = colnames(ma),Group1=Immunetalkscore[colnames(ma),'Group'] %>% as.factor())
bk = unique(c(seq(-2,2, length=100)))
p=pheatmap(ma ,annotation_col=ad,annotation_row=ac, show_colnames =F,show_rownames = T,breaks = bk,
         cluster_rows = F,
         cluster_cols = F,scale = 'row',color=colorRampPalette(c("blue", "black", "yellow"))(100) ) 

ggsave(paste0(save.pic,'Figure5-ITS-immune.pdf'), plot = print(p), width =8, height =15, units ='in')

#--immune talk score---

Immunetalkscore
ICB.gene=c('PDCD1','CD274','CTLA4','LAG3','TIGIT','IDO1','CD80','CD86','PVR','CD200R1','CD200','CD276','LGALS3','BTLA','KIR3DL1','ADORA2A')
cor.data=data.frame(Immunetalkscore$score,TCGA_GTEX_TPM[ICB.gene,rownames(Immunetalkscore)] %>% t() %>% data.frame()  ) 

pic=list()
for(i in colnames(cor.data)){
  tmp=shapiro.test(cor.data[,i])
  pic[[i]]=c(tmp$statistic,tmp$p.value)
}

pic=t(data.frame(pic))



cor.data=cor(cor.data,method = 'spearman')
cor.data=cor.data[1,]
cor.data=cor.data[2:17]
cor.data=data.frame(row.names = rownames(cor.data),gene=names(cor.data),cor=cor.data)
cor.data=dplyr::arrange(cor.data,cor)
cor.data$gene=as.factor(cor.data$gene) %>% {factor(.,levels = cor.data$gene)} 
cor.data[,'group']='orange'
cor.data[cor.data$gene%in% c('PVR','CD276'),'group']='blue'

p=ggplot(cor.data,aes(x=gene,y=cor,fill=group))+geom_bar(stat='identity',position = 'identity')+scale_fill_manual(values = c('orange','blue'))+mytheme+geom_hline(aes(yintercept= -0.4),colour="darkgrey", linetype="dashed")+geom_hline(aes(yintercept= -0.6),colour="darkgrey", linetype="dashed")

ggsave(paste0(save.pic,'Figure5-corralation.pdf'),p,height = 3,width=12)


#immune score 
immune=read.delim(paste0(main.path,'TCGA/immune/Scores_160_Signatures.tsv'),header=T,row.names = 1)
immune=immune[,-1:-2]
colnames(immune)=colnames(immune) %>% substring(1,15)
immune=t(immune) %>% data.frame()
immune=immune[rownames(Immunetalkscore),]
immune=immune[!is.na(immune$Angiogenesis),]


nouse=sapply(colnames(immune), function(x){pval=wilcox.test(immune[rownames(immune) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='High',]),x] %>% as.numeric() , immune[rownames(immune) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='Low',]),x] %>% as.numeric()  )$p.value
return(c(x,pval))
} 
) %>% t() %>% data.frame() %>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}

immune=data.frame(Group=Immunetalkscore[rownames(immune),'Group'] , immune[,rownames(nouse)[c(1:52)]] ) 
immune=melt(immune,ID='Group' )
immune$Group = as.factor(immune$Group ) %>% {factor(.,levels = c('Low','High'))}
p=ggplot(immune,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-2,12)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(legend.position = 'none')+theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))


ggsave(paste0(save.pic,'Figure5-immune signature score.pdf'),p,height = 4,width=10)

#cibersort
cibersort=read.delim(paste0(main.path,'TCGA/immune/TCGA.Kallisto.fullIDs.cibersort.relative.tsv'),header=T)
cibersort$SampleID = cibersort$SampleID  %>% substring(1,15)
cibersort=cibersort[cibersort$SampleID %in% rownames(Immunetalkscore) ,]
cibersort=cibersort[!duplicated(cibersort$SampleID),]
rownames(cibersort)=cibersort$SampleID

nouse=sapply(colnames(cibersort)[-1:-2], function(x){pval=wilcox.test(cibersort[rownames(cibersort) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='High',]),x] %>% as.numeric() , cibersort[rownames(cibersort) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='Low',]),x] %>% as.numeric()  )$p.value
return(c(x,pval))
} 
) %>% t() %>% data.frame()%>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}
nouse=nouse[1:11,]

cibersort=data.frame(Group=Immunetalkscore[rownames(cibersort),'Group'] , cibersort[,rownames(nouse)] ) 
cibersort=melt(cibersort,ID='Group' )
cibersort$Group = as.factor(cibersort$Group ) %>% {factor(.,levels = c('Low','High'))}
p=ggplot(cibersort,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-0.2,1)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))


ggsave(paste0(save.pic,'Figure5-cibsort.pdf'),p,height = 4,width=10)


#TIMER
TIMER=read.delim(paste0(main.path,'TCGA/immune2/TCGA_TIMER.txt'),header=T)
TIMER$Sample_ID = TIMER$Sample_ID %>% str_replace_all('-','.') %>% substring(1,15)
TIMER=TIMER[TIMER$Sample_ID %in% rownames(Immunetalkscore) ,]
TIMER=TIMER[!duplicated(TIMER$Sample_ID),]
rownames(TIMER)=TIMER$Sample_ID



nouse=sapply(colnames(TIMER)[-1:-3], function(x){pval=wilcox.test(TIMER[rownames(TIMER) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='High',]),x] %>% as.numeric() , TIMER[rownames(TIMER) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='Low',]),x] %>% as.numeric()  )$p.value
return(c(x,pval))
} 
) %>% t() %>% data.frame()%>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}


TIMER=data.frame(Group=Immunetalkscore[rownames(TIMER),'Group'] , TIMER[,rownames(nouse)] ) 
TIMER=melt(TIMER,ID='Group' )
TIMER$Group = as.factor(TIMER$Group ) %>% {factor(.,levels = c('Low','High'))}
p=ggplot(TIMER,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-0.2,1)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")

ggsave(paste0(save.pic,'Figure5-TIMER.pdf'),p,height = 4,width=10)


# TIMER=data.frame(Group=Immunetalkscore[rownames(TIMER),'Group'] , TIMER ) 
# write.csv(TIMER,paste0(save.data,'timer.csv'),col.names = T,row.names = T,quote = F)

#TCGA_MCP_counter
MCP_counter=read.delim(paste0(main.path,'TCGA/immune2/TCGA_MCP_counter.txt'),header=T)
MCP_counter$Sample_ID = MCP_counter$Sample_ID %>% str_replace_all('-','.') %>% substring(1,15)
MCP_counter=MCP_counter[MCP_counter$Sample_ID %in% rownames(Immunetalkscore) ,]
MCP_counter=MCP_counter[!duplicated(MCP_counter$Sample_ID),]

rownames(MCP_counter)=MCP_counter$Sample_ID


now.data=data.frame(MCP_counter,score=Immunetalkscore[rownames(MCP_counter),'score'])
p=ggplot(now.data,aes(x=score,y=Cytotoxic.lymphocytes))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')+ylim(0,2.5)
ggsave(paste0(save.pic,'Figure5-TIMER-cyto.pdf'),p,height = 5,width=5)


p=ggplot(now.data,aes(x=score,y=Endothelial.cells))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')

cor.test(now.data$score,now.data$Endothelial.cells,method='spearman')
shapiro.test(now.data$Endothelial.cells) 

ggsave(paste0(save.pic,'Figure4add-TIMER-endo.pdf'),p,height = 5,width=5)




nouse=sapply(colnames(MCP_counter)[-1:-3], function(x){pval=wilcox.test(MCP_counter[rownames(MCP_counter) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='High',]),x] %>% as.numeric() , MCP_counter[rownames(MCP_counter) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='Low',]),x] %>% as.numeric()  )$p.value
return(c(x,pval))
} 
) %>% t() %>% data.frame()%>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}


MCP_counter=data.frame(Group=Immunetalkscore[rownames(MCP_counter),'Group'] ,MCP_counter[,rownames(nouse)] ) 
MCP_counter=melt(MCP_counter,ID='Group' )
MCP_counter$Group = as.factor(MCP_counter$Group ) %>% {factor(.,levels = c('Low','High'))}
p=ggplot(MCP_counter,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-0.2,6)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")

ggsave(paste0(save.pic,'Figure5-MCP_counter.pdf'),p,height = 4,width=10)

# MCP_counter=data.frame(Group=Immunetalkscore[rownames(MCP_counter),'Group'] ,MCP_counter ) 
# write.csv(MCP_counter,paste0(save.data,'MCP_counter.csv'),col.names = T,row.names = T,quote = F)



#xcell
xcell=read.delim(paste0(main.path,'TCGA/immune2/TCGA_xCell.txt'),header=T)
xcell$Sample_ID =xcell$Sample_ID %>% substring(1,15)
xcell=xcell[xcell$Sample_ID %in% rownames(Immunetalkscore) ,]
rownames(xcell)=xcell$Sample_ID
xcell=xcell[!duplicated(xcell$Sample_ID),]

nouse=sapply(colnames(xcell)[-1:-2], function(x){pval=wilcox.test(xcell[rownames(xcell) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='High',]),x] %>% as.numeric() , xcell[rownames(xcell) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='Low',]),x] %>% as.numeric()  )$p.value
return(c(x,pval))
} 
) %>% t() %>% data.frame()%>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}

xcell=data.frame(Group=Immunetalkscore[rownames(xcell),'Group'] ,xcell[,rownames(nouse)] ) 
xcell=melt(xcell,ID='Group' )
xcell$Group = as.factor(xcell$Group ) %>% {factor(.,levels = c('Low','High'))}
xcell=xcell[xcell$variable %in% c('Adipocytes','CD4..memory.T.cells','CD4..Tem','CD8..naive.T.cells','CD8..T.cells','Endothelial.cells',
                                  'Fibroblasts','Epithelial.cells','Macrophages.M2','Melanocytes','Pericytes','Th2.cells','Tregs'),]
xcell$variable=as.factor(xcell$variable) %>% {factor(.,levels =c('Epithelial.cells','Fibroblasts','Endothelial.cells','Pericytes','Adipocytes','CD4..memory.T.cells','CD4..Tem','Th2.cells','Tregs','CD8..naive.T.cells','CD8..T.cells',
                                                                          'Macrophages.M2','Melanocytes') )}
 
p=ggplot(xcell,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-0.2,1.2)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))
ggsave(paste0(save.pic,'Figure5-xcell.pdf'),p,height = 4,width=10)



#ssgsea
ssgsea=read.delim(paste0(main.path,'TCGA/immune2/TCGA_Cell2013_ssgsea_df.txt'),header=T)
ssgsea$Sample_ID =ssgsea$Sample_ID  %>% str_replace_all('-','.') %>% substring(1,15)
ssgsea=ssgsea[ssgsea$Sample_ID %in% rownames(Immunetalkscore) ,]
ssgsea=ssgsea[!duplicated(ssgsea$Sample_ID),]
rownames(ssgsea)=ssgsea$Sample_ID

nouse=sapply(colnames(ssgsea)[-1:-3], function(x){pval=wilcox.test(ssgsea[rownames(ssgsea) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='High',]),x] %>% as.numeric() , ssgsea[rownames(ssgsea) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='Low',]),x] %>% as.numeric()  )$p.value
return(c(x,pval))
} 
) %>% t() %>% data.frame()%>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}

ssgsea=data.frame(Group=Immunetalkscore[rownames(ssgsea),'Group'] ,ssgsea[,rownames(nouse)] ) 
ssgsea=melt(ssgsea,ID='Group' )
ssgsea$Group = as.factor(ssgsea$Group ) %>% {factor(.,levels = c('Low','High'))}
p=ggplot(ssgsea,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-0.2,1)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))
ggsave(paste0(save.pic,'Figure5-ssgsea.pdf'),p,height = 4,width=5)

#epic
epic=read.delim(paste0(main.path,'TCGA/gepia2021/EPIC.csv'),header=T,sep=',')
rownames(epic)=epic[,1]
epic=epic[rownames(epic) %in% rownames(Immunetalkscore) ,]
nouse=sapply(colnames(epic)[2:9], function(x){pval=wilcox.test(epic[rownames(epic) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='High',]),x] %>% as.numeric() , epic[rownames(epic) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='Low',]),x] %>% as.numeric()  )$p.value
return(c(x,pval))
} 
) %>% t() %>% data.frame() %>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}

epic=data.frame(Group=Immunetalkscore[rownames(epic),'Group'] ,epic[,rownames(nouse)] ) 
epic=melt(epic,ID='Group' )
epic$Group = as.factor(epic$Group ) %>% {factor(.,levels = c('Low','High'))}

epic$variable=as.factor(epic$variable) %>% {factor(.,levels =c('mRNAProportions.CAFs','mRNAProportions.Endothelial','mRNAProportions.CD4_Tcells','mRNAProportions.CD8_Tcells','mRNAProportions.NKcells','mRNAProportions.Bcells','mRNAProportions.Macrophages','mRNAProportions.otherCells') )}


p=ggplot(epic,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-0.2,1)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))
ggsave(paste0(save.pic,'Figure5-epic.pdf'),p,height = 4,width=5)



#Quantiseq
Quantiseq=read.delim(paste0(main.path,'TCGA/gepia2021/Quantiseq.csv'),header=T,sep=',',row.names = 2) 
Quantiseq=Quantiseq[,-1]
Quantiseq=Quantiseq%>% t() %>% data.frame()
Quantiseq=Quantiseq[rownames(Quantiseq) %in% rownames(Immunetalkscore) ,]
nouse=sapply(colnames(Quantiseq)[1:10], function(x){pval=wilcox.test(Quantiseq[rownames(Quantiseq) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='High',]),x] %>% as.numeric() , Quantiseq[rownames(Quantiseq) %in% rownames(Immunetalkscore[Immunetalkscore$Group=='Low',]),x] %>% as.numeric()  )$p.value
return(c(x,pval))
} 
) %>% t() %>% data.frame() %>% {.[,2]=as.character(.[,2]) %>% as.numeric();.} %>% {.[.[,2]<0.05,]}

Quantiseq=data.frame(Group=Immunetalkscore[rownames(Quantiseq),'Group'] ,Quantiseq[,rownames(nouse)] ) 
Quantiseq=melt(Quantiseq,ID='Group' )
Quantiseq$Group = as.factor(Quantiseq$Group ) %>% {factor(.,levels = c('Low','High'))}
p=ggplot(Quantiseq,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(0,0.25)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))
ggsave(paste0(save.pic,'Figure5-Quantiseq.pdf'),p,height = 4,width=5)



#TCRclone
TCRclone=read.delim(paste0(main.path,'TCGA/immune2/TCGA_TCR_clone.txt'),header=T)
TCRclone$Sample_ID =TCRclone$Sample_ID  %>% str_replace_all('-','.') %>% substring(1,15)
TCRclone=TCRclone[TCRclone$Sample_ID %in% rownames(Immunetalkscore) ,]
TCRclone=TCRclone[!duplicated(TCRclone$Sample_ID),]
rownames(TCRclone)=TCRclone$Sample_ID
TCRclone1=data.frame(Group=Immunetalkscore[rownames(TCRclone),'Group'] ,TCRclone[,c("shannon")] ) 
TCRclone1=melt(TCRclone1,ID='Group' )
TCRclone1$Group = as.factor(TCRclone1$Group ) %>% {factor(.,levels = c('Low','High'))}
p=ggplot(TCRclone1,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-0.2,6)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))
ggsave(paste0(save.pic,'Figure5-TCR_clone_shannon.pdf'),p,height = 4,width=3)

TCRclone2=data.frame(Group=Immunetalkscore[rownames(TCRclone),'Group'] ,TCRclone[,c("numClones")] ) 
TCRclone2=melt(TCRclone2,ID='Group' )
TCRclone2$Group = as.factor(TCRclone2$Group ) %>% {factor(.,levels = c('Low','High'))}
p=ggplot(TCRclone2,aes(x=variable,y=value,fill=Group))+geom_boxplot()+mytheme+ylim(-0.2,60)+scale_fill_manual(values=c('aquamarine','lightsalmon'))+
  theme(axis.title.x = element_blank())+labs(y='Signature Score')+
  stat_compare_means(aes(group=Group),label = "p.signif")+theme(axis.text.x = element_text(angle=20,hjust=0.5, vjust=0.5,size = rel(1)))
ggsave(paste0(save.pic,'Figure5-TCR_clone_num.pdf'),p,height = 4,width=3)


TLS=c('CCL19','CCL21','CXCL13','CCR7','SELL','LAMP3','CXCR4','CD86','BCL6')

Immunetalkscore.TLS=Immunetalkscore
Immunetalkscore.TLS[,'TLS']=apply(TCGA_GTEX_TPM[TLS,rownames(Immunetalkscore.TLS)],2,mean)

p=ggplot(Immunetalkscore.TLS,aes(x=score,y=TLS))+geom_point(size=3,shape=19,color='blue',alpha=0.5)+mytheme+stat_smooth(method = lm,color='lightblue')
ggsave(paste0(save.pic,'FigureadTLS-TSTS.pdf'),p,height = 5,width=5)

shapiro.test(Immunetalkscore.TLS$score) 
shapiro.test(Immunetalkscore.TLS$TLS) 

cor.test(Immunetalkscore.TLS$score,Immunetalkscore.TLS$TLS,method = 'spearman')
write.csv(Immunetalkscore.TLS,paste0(save.data,'TLS-TSTS.csv'))


