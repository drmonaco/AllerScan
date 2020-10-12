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


key_r = x_IgE %>% select(c(names)) %>% rownames_to_column()%>%  mutate(names= ifelse(names == "111Other","Other",names))%>% column_to_rownames()

x_IgE = x_IgE %>% select(-c(V11,V12,names))

IgE_plot = pheatmap(as.matrix(x_IgE)-1,show_colnames= FALSE,show_rownames = FALSE,
                    cluster_rows = FALSE,cluster_cols = FALSE, annotation_col = key_2,
                    breaks = c(0,5, 10, 20, 50, 100),col = viridis(5),legend_breaks = c(0,5, 10, 20, 50, 100),
                    annotation_names_row = FALSE,annotation_names_col = FALSE,annotation_colors = mat_colors_col2,
                    annotation_row = key_r,labels_row = 1:24 %>% as.character(),annotation_legend = FALSE,legend = FALSE)
# 


x_IgG_2 = x_IgE %>% rownames_to_column() %>% select("rowname") %>%left_join(x_IgG %>% rownames_to_column(),by = "rowname")%>% select(-c(V11,V12,names)) %>% column_to_rownames()
#     
IgG_plot = pheatmap(as.matrix(x_IgG_2)-1,show_colnames= FALSE,show_rownames = FALSE,
                    cluster_rows = FALSE,cluster_cols = FALSE, annotation_col = key_2,
                    breaks = c(0,5, 10, 20, 50, 100),col = viridis(5),legend_breaks = c(0,5, 10, 20, 50, 100),
                    annotation_names_row = FALSE,annotation_names_col = FALSE)
#     
# 
d1 <- dist(x_IgG_2,method = "euclidean", diag = FALSE, upper = FALSE)
d2 <- dist(t(x_IgG_2),method = "euclidean", diag = FALSE, upper = TRUE)
c1 <- hclust(d1, method="ward.D2", members = NULL)
c2 <- hclust(d2, method="ward.D2", members = NULL)

IgG_plot2 = pheatmap(as.matrix(x_IgG_2)-1,show_colnames= FALSE,show_rownames = FALSE,
                     cluster_rows = rev(c1),cluster_cols = F, annotation_col = key_2,
                     breaks = c(0,5, 10, 20, 50, 100),col = viridis(5),legend_breaks = c(0,5, 10, 20, 50, 100),annotation_names_row = FALSE,annotation_names_col = FALSE,cutree_rows = 3,annotation_colors = mat_colors_col2)


epitope = cutree(c1,k=3)  %>% as.data.frame() %>% rownames_to_column(var = "u_pep_id") %>% left_join(y = annot %>% select(c(u_pep_id,V28)), by = "u_pep_id") %>% rename(cluster = 2) %>%
  filter(cluster ==3) %>% select(u_pep_id) %>% unlist

IgG_plot2 = pheatmap(as.matrix(x_IgG_2)-1,show_colnames= FALSE,show_rownames = FALSE,
                     cluster_rows = rev(c1),cluster_cols = F, annotation_col = key_2,
                     breaks = c(0,5, 10, 20, 50, 100),col = viridis(5),legend_breaks = c(0,5, 10, 20, 50, 100),annotation_names_row = FALSE,annotation_names_col = FALSE,cutree_rows = 3,annotation_colors = mat_colors_col2)


epitope = cutree(c1,k=3)  %>% as.data.frame() %>% rownames_to_column(var = "u_pep_id") %>% left_join(y = annot %>% select(c(u_pep_id,V28)), by = "u_pep_id") %>% rename(cluster = 2) %>%
  filter(cluster ==3) %>% select(u_pep_id) %>% unlist

abc = x_IgG %>% rownames_to_column() %>% filter(grepl(epitope %>% paste0(collapse="|"),rowname)) %>%
  column_to_rownames(var= "rowname") %>% t() %>%  as.data.frame() %>%  rownames_to_column() %>% rename(Peptide = 2) %>% select(rowname,Peptide)


test1_G = key_2 %>% rownames_to_column() %>% as.data.frame() %>% left_join(y=  abc,by = "rowname") %>% gather(key="group",value = "Hits FC",-rowname,-Wheat) %>% mutate(`Hits FC` = as.numeric(`Hits FC`))

puro_IgG  = ggplot(test1_G,aes(x = Wheat,y = `Hits FC`))+
  geom_violin(aes(fill = Wheat),position = position_dodge(width = 0.9),scale = "width")+
  geom_quasirandom(aes(group = Wheat),size=  3,dodge.width = 0, varwidth = FALSE)+
  scale_y_log10()+stat_compare_means(comparisons = list( c("A", "H"), c("A", "S"), c("H", "S")))+theme(text = element_text(size=20))+theme(legend.position = "none")+ggtitle("IgG Reactivity to Purothionin Motif by Wheat Allergy")+theme_classic()

abc = x_IgE %>% rownames_to_column() %>% filter(grepl(epitope %>% paste0(collapse="|"),rowname)) %>%
  column_to_rownames(var= "rowname") %>% t() %>%  as.data.frame() %>%  rownames_to_column() %>%rename(Peptide = 2)%>% select(rowname,Peptide)

test1_E = key_2 %>% rownames_to_column() %>% as.data.frame() %>% left_join(y=  abc,by = "rowname") %>% gather(key="group",value = "Hits FC",-rowname,-Wheat) %>% mutate(`Hits FC` = `Hits FC` %>% unlist() %>% as.numeric())


puro_IgE = ggplot(test1_E,aes(x = Wheat,y = `Hits FC`))+
  geom_violin(aes(fill = Wheat),position = position_dodge(width = 0.9),scale = "width")+
  geom_quasirandom(aes(group = Wheat),size=  3,dodge.width = 0, varwidth = FALSE)+
  scale_y_log10()+stat_compare_means(comparisons = list( c("A", "H"), c("A", "S"), c("H", "S")))+theme(text = element_text(size=20))+ggtitle("IgE Reactivity to Purothionin Motif by Wheat Allergy")+theme_classic()

test_GE = test1_E %>% left_join(test1_G,by = "rowname") %>% select(rowname,Wheat.x,`Hits FC.x`,`Hits FC.y`) %>% rename(IgE =`Hits FC.x`,IgG =  `Hits FC.y`) %>% ggplot(aes(x = IgG,y = IgE,color = Wheat.x))+geom_jitter(width = 0.1,height = .1,size=  3)+scale_x_log10(name ="IgG Reactivity")+scale_y_log10(name ="IgE Reactivity")+theme(text = element_text(size=20))+ggtitle("Reactivity to Purothionin Motif by Isotype")+theme_classic()
