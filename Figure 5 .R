rm(list=ls())
options(stringsAsFactors = FALSE)
knitr::opts_chunk$set(echo = TRUE)
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

x_IgE = df_hits_clean %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_clean)[2])) %>% 
  left_join(y=peanut_pep,by = "u_pep_id") %>% arrange(index) %>% select(1:160) %>% column_to_rownames(var="u_pep_id")

key_2 = x_IgE %>% colnames %>% as.data.frame() %>% rename(Name=  1) %>%
  left_join(y= key %>% select(-sIgE_W),by = "Name") %>% select(-c(Peanut)) %>% filter(!is.na(Wheat)) %>%
  rename(xx  = 1)  %>% filter(!grepl("unclear",Wheat,ignore.case = TRUE))%>% arrange(Wheat) %>%  distinct() %>%filter(Wheat =="A")%>% 
  column_to_rownames(var="xx")

row_names = wheat_pep[,c(1,10,11,12,13)] %>% arrange(V11) %>% mutate(names = ifelse(V12 >= 5,product,"Other")) %>%
  column_to_rownames(var="u_pep_id") %>% select(-c(product,V12,V13))

x_IgE1=x_IgE %>%   select(colnames(x_IgE)[colnames(x_IgE) %in% rownames(key_2)]) %>% select(rownames(key_2)) %>% merge(y = row_names,by = 0) %>% arrange(V11) %>% select(-c(V11,names)) %>% column_to_rownames(var= "Row.names")


x_IgE3 = x_IgE1 %>% rownames_to_column()%>%  mutate_each(function(x) ifelse(x>1,1,0),-rowname) %>% column_to_rownames(var= "rowname")


x_IgE4 = x_IgE3 %>% select(which(!colSums(., na.rm=TRUE) < 10)) %>% cbind(rowSums(.)) %>% rename(col = `rowSums(.)`) %>% rownames_to_column() %>% filter(col >3) %>% select(-col) %>% column_to_rownames(var = "rowname")

d1 <- dist(x_IgE4,method = "binary", diag = FALSE, upper = FALSE)
d2 <- dist(t(x_IgE4),method = "binary", diag = FALSE, upper = TRUE)
c1 <- hclust(d1, method="ward.D2", members = NULL)
c2 <- hclust(d2, method="ward.D2", members = NULL)



x_IgE5 = x_IgE1 %>% select(which(!colSums(x_IgE3, na.rm=TRUE) < 10)) %>% cbind(rowSums(x_IgE3)) %>% rename(col = `rowSums(x_IgE3)`) %>% filter(col >3) %>% select(-col) 



IgE_plot5 = pheatmap(as.matrix(x_IgE5),show_colnames= FALSE,show_rownames = FALSE,
                     cluster_rows = rev(c1),cluster_cols = c2,clustering_method = "ward.D2", annotation_col = key_2,
                     breaks = c(0,5, 10, 20, 75, 100),col = viridis(5),legend_breaks = c(0,5, 10, 20, 75, 100),annotation_names_row = FALSE,annotation_names_col = FALSE,cutree_rows = 6,cutree_cols = 3,border_color = NA)