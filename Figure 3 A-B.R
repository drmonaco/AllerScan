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


x_IgE = df_hits_clean %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_clean)[2])) %>% 
  left_join(y=peanut_pep,by = "u_pep_id") %>% arrange(index) %>% select(1:160) %>% column_to_rownames(var="u_pep_id")

x_IgG = df_hits_IgG %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_IgG)[2])) %>% 
  left_join(y=peanut_pep,by = "u_pep_id") %>% arrange(index) %>% select(1:148) %>% column_to_rownames(var="u_pep_id")

key_2 = x_IgG %>% colnames %>% as.data.frame() %>% rename(Name=  1) %>%
  left_join(y= key %>% select(-sIgE_W),by = "Name") %>% select(-c(Peanut)) %>% filter(!is.na(Wheat)) %>%
  rename(xx  = 1)  %>% filter(!grepl("unclear",Wheat,ignore.case = TRUE))%>% arrange(Wheat) %>%  distinct() %>%
  column_to_rownames(var="xx")

row_names = wheat_pep[,c(1,10,11,12,13)] %>% mutate(names = ifelse(V12 >= 10,product,"111Other")) %>% arrange(V12,names,V11) %>% 
  column_to_rownames(var="u_pep_id") %>% select(-c(product,V13))

x_IgE=x_IgE %>%   select(colnames(x_IgE)[colnames(x_IgE) %in% rownames(key_2)]) %>% select(rownames(key_2)) %>% merge(y = row_names,by = 0) %>% arrange(V12,names,V11) %>% column_to_rownames(var= "Row.names")

x_IgG=x_IgG %>%   select(colnames(x_IgE)[colnames(x_IgE) %in% rownames(key_2)]) %>% select(rownames(key_2)) %>% merge(y = row_names,by = 0) %>% arrange(V12,names,V11) %>% column_to_rownames(var= "Row.names")

mat_colors_col2 <- c(list(Wheat = c(viridis(3))))#,list(names = viridis(nrow(unique(row_names %>% select(-V11))))))
names(mat_colors_col2$Wheat) <- c(unique(key_2$Wheat))

#####
# names(mat_colors_col2$names) <- c(unique(row_names %>% select(-V11)) %>% unlist %>% as.character())

# mat_colors_col2 <- list(Wheat = c(viridis(3)))
#  names(mat_colors_col2$Wheat) <- c(unique(key_2$Wheat))

# xx =  epitope_col %>% rownames_to_column(var=  "Name") %>% left_join(y = IgE1,by=  "Name")
# xx$cluster = factor( xx$cluster,levels = c(2, 4, 3, 1))
# ggplot(xx,aes(x = cluster,y=  Breadth))+geom_violin(scale = "width",aes(fill = cluster))+geom_quasirandom(size = 3)+stat_compare_means(comparisons = list(c(1,2),c(3,4)))+theme(text = element_text(size=20))+scale_y_continuous("IgE Breath to Wheat")+scale_x_discrete("Patient cluster")+theme_classic()
key_r = x_IgE %>% select(c(names)) %>% rownames_to_column()%>%  mutate(names= ifelse(names == "111Other","Other",names))%>% column_to_rownames()
# key_r2 = apply(unique(key_r$names) %>% unlist()%>% as.character()%>% cbind(1:length(.),.),1,paste0,collapse = "")
# key_r2 = cbind(unique(key_r$names),key_r2) %>% as.data.frame()

# key_r = key_r %>% left_join(y = key_r2 %>% rename(names= V1),by = "names") 
x_IgE = x_IgE %>% select(-c(V11,V12,names))

IgE_plot = pheatmap(as.matrix(x_IgE)-1,show_colnames= FALSE,show_rownames = FALSE,
                    cluster_rows = FALSE,cluster_cols = FALSE, annotation_col = key_2,
                    breaks = c(0,5, 10, 20, 50, 100),col = viridis(5),legend_breaks = c(0,5, 10, 20, 50, 100),
                    annotation_names_row = FALSE,annotation_names_col = FALSE,annotation_colors = mat_colors_col2,
                    annotation_row = key_r,labels_row = 1:24 %>% as.character(),annotation_legend = FALSE,legend = FALSE)


x_IgG_2 = x_IgE %>% rownames_to_column() %>% select("rowname") %>%left_join(x_IgG %>% rownames_to_column(),by = "rowname")%>% select(-c(V11,V12,names)) %>% column_to_rownames()
#     
IgG_plot = pheatmap(as.matrix(x_IgG_2)-1,show_colnames= FALSE,show_rownames = FALSE,
                    cluster_rows = FALSE,cluster_cols = FALSE, annotation_col = key_2,
                    breaks = c(0,5, 10, 20, 50, 100),col = viridis(5),legend_breaks = c(0,5, 10, 20, 50, 100),
                    annotation_names_row = FALSE,annotation_names_col = FALSE,annotation_legend = FALSE,legend = FALSE)
#  