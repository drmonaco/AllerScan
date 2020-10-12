rm(list=ls())
options(stringsAsFactors = FALSE)
knitr::opts_chunk$set(echo = TRUE)
getwd()
load("my_work_space.RDS") 
library(RColorBrewer)
library(tidyverse)
library(fpc)
library(data.table)
library(knitr)
library(rafalib)
library(dplyr)
library(gplots)
library(pheatmap)
library(viridis)
library(ggpubr)
library(tidyverse)
library(magrittr)
library(cluster)
library(cowplot)
library(NbClust)
library(clValid)
library(ggfortify)
library(clustree)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(GGally)
library(ggiraphExtra)
library(knitr)
library(kableExtra)
library(igraph)
library(vipor)
library(ggbeeswarm)
library(plotly)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

edge = fread("~/Stuff to sort 2:5:20/allergome_dict.csv")
filter_ind  = function(edge,vertex){ #independence filter that takes a dictionary (defined above) and a set of nodes and tells the minimal number of unique epitopes
  nodes = unlist(vertex)
  links_filtered = subset(edge,unlist(edge[,1]) %in% nodes)
  links_filtered = subset(links_filtered,links_filtered[,2] %>% unlist() %in% nodes)
  if(dim(links_filtered)[1]!=0){
    
    net <- as.undirected(graph_from_data_frame(d=links_filtered,vertices=nodes, directed=F) )
    x = decompose.graph(net)
    x_1 = x[sapply(x,vcount)<30]
    x_1_sum  = sum(unlist(lapply(x_1,independence.number)))
    x_2 = x[sapply(x,vcount)>=30]
    temp = c()
    #x_2 = x
    if(length(x_2) >0){
      for(R in 1:length(x_2)){
        x_2_r = x_2[[R]]
        while(max(degree(x_2_r)>5)){
          
          toss = degree(x_2_r)==max(degree(x_2_r))
          x_2_r = delete_vertices(x_2_r, V(x_2_r)[toss][1])
        }
        x_l = decompose.graph(x_2_r)
        temp[R] = sum(unlist(lapply(x_l,independence.number)))
      }
    }
    return(sum(x_1_sum)+sum(temp))
  }
  if(dim(links_filtered)[1]==0){
    return(length(nodes))
  }
}


x = df_hits_clean %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_clean)[2])) %>% column_to_rownames(var = "u_pep_id") %>% 
  mutate_all(function(x) ifelse(x>1,1,0)) %>% t()

breadth_underfiltered = x %>% as.data.frame() %>% rownames_to_column(var="Name") %>% 
  mutate(sumVar = rowSums(.[,2:dim(.)[2]]))  %>% 
  select(c(Name,sumVar)) %>% left_join(y = key, by = "Name") %>% select(c(Name,sumVar,Wheat)) %>% 
  filter(!grepl("unclear|Unclear",Wheat)) %>% filter(!is.na(Wheat)) %>% distinct()

# plot = ggplot(breadth_underfiltered,aes(y = sumVar+1,x = Wheat,color= Wheat))+geom_violin(scale = "width")+ geom_jitter(height = 0, width = 0.1) +scale_y_log10()+   stat_compare_means(comparisons = list( c("A", "H"), c("A", "S"), c("H", "S")))

filtered =  df_hits_clean %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_clean)[2])) %>% gather(key = "id",value = "FC",-u_pep_id) %>% filter(FC>1)


output <- as.data.frame(matrix(ncol=2,nrow = length(unique(filtered$id) )))
for(R in 1:length(unique(filtered$id))){
  output[R,2]= filter_ind(edge = edge,vertex = filtered$u_pep_id[filtered$id %in% unique(filtered$id)[R]])
  output[R,1] = unique(filtered$id)[R]
}
colnames(output) = c("Name","Breadth")
IgE1 = output%>% mutate(Breadth = ifelse(is.na(Breadth),0,Breadth))
x123 =  breadth_underfiltered %>% left_join(y = output,by=  "Name") %>% mutate(Breadth = ifelse(is.na(Breadth),0,Breadth))



plot2 = ggplot(x123,aes(y = Breadth+1,x = Wheat,fill= Wheat))+
  geom_violin(scale = "width")+ 
  geom_quasirandom() +
  scale_y_log10("Breadth")+   
  stat_compare_means(comparisons = list( c("A", "H"), c("A", "S"), c("H", "S")))+theme(text = element_text(size=20))+ggtitle("IgE Breadth of Response to Wheat Peptides by Wheat Allergy")+theme_classic()+theme(legend.position = "none")


x = df_hits_IgG %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_IgG)[2])) %>% column_to_rownames(var = "u_pep_id") %>% 
  mutate_all(function(x) ifelse(x>1,1,0)) %>% t()

breadth_underfiltered = x %>% as.data.frame() %>% rownames_to_column(var="Name") %>% 
  mutate(sumVar = rowSums(.[,2:dim(.)[2]]))  %>% 
  select(c(Name,sumVar)) %>% left_join(y = key, by = "Name") %>% select(c(Name,sumVar,Wheat)) %>% 
  filter(!grepl("unclear|Unclear",Wheat)) %>% filter(!is.na(Wheat)) %>% distinct()


filtered =  df_hits_IgG %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_IgG)[2])) %>% gather(key = "id",value = "FC",-u_pep_id) %>% filter(FC>1)


output <- as.data.frame(matrix(ncol=2,nrow = length(unique(filtered$id) )))
for(R in 1:length(unique(filtered$id))){
  output[R,2]= filter_ind(edge = edge,vertex = filtered$u_pep_id[filtered$id %in% unique(filtered$id)[R]])
  output[R,1] = unique(filtered$id)[R]
}
colnames(output) = c("Name","Breadth")
IgG1 = output%>% mutate(Breadth = ifelse(Breadth=="NA",0,Breadth))
x123 =  breadth_underfiltered %>% left_join(y = output,by=  "Name") %>% mutate(Breadth = ifelse(is.na(Breadth),0,Breadth))

plot2_IgG = ggplot(x123,aes(y = Breadth+1,x = Wheat,fill= Wheat))+
  geom_violin(scale = "width")+ 
  geom_quasirandom() +
  scale_y_log10("Breadth")+   
  stat_compare_means(comparisons = list( c("A", "H"), c("A", "S"), c("H", "S")))+theme(text = element_text(size=20))+ggtitle("IgG Breadth of Response to Wheat Peptides by Wheat Allergy")+theme_classic()+theme(legend.position = "none")

