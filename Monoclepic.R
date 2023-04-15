Monoclepic=function(nowseurat,use.col='recluster',num=100,onlyTumor=T,max_components =2){
  set.seed(1234)
  metadata=nowseurat@meta.data
  if(onlyTumor==T){
    metadata=metadata[metadata$Tissue%in% c('Tumor') ,]
  }
  monocling=list()
  for(i in unique(metadata[,use.col])){
    df=metadata[metadata[,use.col]==i ,]
    if(nrow(df)<num){
      monocling[[i]]=df
    }else{
      monocling[[i]]=df[sample(nrow(df),num), ]
    }
    monocling[[i]][,'barcode']=rownames(monocling[[i]])
  }
  monocling=do.call(rbind,lapply(monocling, data.frame))
  rownames(monocling)=monocling$barcode
  
  count.one.type.bei=subset(nowseurat,cells = rownames(monocling))
  
  cols.gene=rownames(count.one.type.bei@assays$RNA@data)
  Mono=newCellDataSet(as.matrix(count.one.type.bei@assays$RNA@data),phenoData = new('AnnotatedDataFrame', data = count.one.type.bei@meta.data),
                      featureData = new('AnnotatedDataFrame', data = data.frame(gene_short_name = cols.gene, row.names = cols.gene))  )
  Mono<- estimateSizeFactors(Mono) %>% estimateDispersions()%>% detectGenes(min_expr = 3)
  unsup_clustering_genes <- subset(dispersionTable(Mono), mean_expression >= 0.1)
  Var.gene=unsup_clustering_genes$gene_id %>% str_extract('MT-.+|RPS.+|RPL.+') %>% na.omit() %>% as.character()
  unsup_clustering_genes=unsup_clustering_genes[!unsup_clustering_genes$gene_id %in% Var.gene,]
  Mono = setOrderingFilter(Mono, unsup_clustering_genes$gene_id) %>% reduceDimension(mmax_components =2, num_dim = 2,reduction_method = 'tSNE', verbose = T)
  Mono <-  reduceDimension(Mono,max_components =2,method = 'DDRTree') %>% orderCells()
  return(Mono)
}