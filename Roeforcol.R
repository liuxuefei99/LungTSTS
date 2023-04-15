Roeforcol=function(meta,usecolname){
  tissue.per=meta[,usecolname] %>% table() %>% prop.table() %>% data.frame()
  expect=now.meta$recluster %>% table() %>% data.frame() %>% {.[.[,'Freq']!=0,];.}
  expect.now=(matrix(tissue.per$Freq)%*% t(matrix(expect$Freq)) %>% t())
  observed=meta %>% data.table() %$% .[,.N,.(recluster,eval(parse(text=usecolname ))   )] %>% data.frame() %>% { dcast(.,recluster~ parse) } %>% {rownames(.)=.$recluster;.} %>% { .[,-1] } %>% {.[is.na(.)]=0;.}
  Roe=observed/expect.now
  return(Roe)
}