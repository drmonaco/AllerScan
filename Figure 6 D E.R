WOIT_key = fread("WOIT_data.csv")
x = fread("wheat_test.csv")

# save.image(file = "my_work_space2.RDS") 
load("my_work_space2.RDS") 


Puro = df_hits_clean %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_clean)[2])) %>% filter(u_pep_id %in% epitope) %>% select(-u_pep_id)%>% t() %>%  as.data.frame() %>% rownames_to_column()%<>% 
  mutate_each(as.numeric,-rowname) %>% mutate(x = apply(.[2:4],1,sum)) %>% mutate(x = x/3) %>% select(rowname,x) %>% rename(V1 = x)
Puro2 = df_hits_IgG %>% filter(V28 == "Wheat") %>% select(c(1,27:dim(df_hits_IgG)[2])) %>% filter(u_pep_id %in% epitope) %>% select(-u_pep_id)%>% t() %>%  as.data.frame() %>% rownames_to_column()%<>% 
  mutate_each(as.numeric,-rowname) %>% mutate(x = apply(.[2:4],1,sum)) %>% mutate(x = x/3)%>% select(rowname,x) %>% rename(V1 =x)

Puro3 = Puro %>% full_join(y = Puro2, by = "rowname") %>% rename(Name = rowname, Puro_IgE = V1.x, Puro_IgG = V1.y) %>% mutate(Puro_IgE = as.numeric(Puro_IgE)) %>% mutate(Puro_IgG = as.numeric(Puro_IgG))%>% filter(grepl("JH|FAIM|T1|T2|T3",Name))

breadth_GE = IgG1 %>% full_join(IgE1,by=  "Name") %>% filter(grepl("JH|FAIM|T1|T2|T3",Name)) %>% left_join(y = Puro3,by = "Name") %>% arrange(Name) %>%
  rename(IgG_breadth = Breadth.x,IgE_breadth = Breadth.y) %>% mutate(Time = substr(Name,1, 2)) %>%
  mutate(ID = str_sub(Name,start =4))%>% full_join(y = WOIT_key, by = "ID") %>% group_by(ID) %>% filter(n()>=2) %>%  mutate(key = ifelse(Time == "T1","keep",ifelse(Time =="T2",ifelse(original == "Placebo","Placebo","Treatment"),ifelse(original == "Wheat","remove","Treatment")))) %>% 
  filter(key != "remove")

x_placebo = breadth_GE %>% filter(key != "Treatment")%>% group_by(ID) %>% filter(n()>=2) %>% mutate(difference = ifelse(as.numeric(New_dose)>1400,"pass","fail"))
x_treatment =breadth_GE %>% mutate(key2 = ifelse(original == "Placebo",ifelse(Time != "T1",ifelse(Time=="T3","Treatment","keep"),"Placebo"),key)) %>% filter(key2 != "Placebo") %>%
  select(-key) %>% rename(key = key2) %>% group_by(ID) %>% filter(n()>=2) %>% mutate(difference = ifelse(as.numeric(New_dose)>1400,"pass","fail"))%>% 
  mutate(difference = ifelse(is.na(difference),"grey",difference))
#placebo_IgE_breadth
ggplot(x_placebo,aes(x = key,y = IgE_breadth))+geom_violin(aes(fill = key)) + geom_point(size = 3)+geom_line(size = 2,aes(group = ID,color = difference)) +scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","black","grey"))+
  stat_compare_means(method  = "wilcox",comparisons = list(c("keep","Placebo")))+ggtitle(wilcox.test(IgE_breadth ~ key, data = x_placebo %>% select(ID,IgE_breadth,key) %>% filter(!is.na(IgE_breadth)) %>% filter(n()>=2), paired = TRUE)$p.value)
#placebo_IgG_breadth
ggplot(x_placebo,aes(x = key,y = IgG_breadth))+geom_violin(aes(fill = key)) + geom_point(size = 3)+geom_line(size = 2,aes(group = ID,color = difference)) +scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","black","grey"))+
  stat_compare_means(method  = "wilcox",comparisons = list(c("keep","Placebo")))+ggtitle(wilcox.test(IgG_breadth ~ key, data = x_placebo %>% select(ID,IgG_breadth,key) %>% filter(!is.na(IgG_breadth)) %>% filter(n()>=2), paired = TRUE)$p.value)
#treatment_IgE_breadth
ggplot(x_treatment,aes(x = key,y = IgE_breadth))+geom_violin(aes(fill = key)) + geom_point(size = 3)+geom_line(size = 2,aes(group = ID,color = difference)) +scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","grey","black"))+
  stat_compare_means(method  = "wilcox",comparisons = list(c("keep","Treatment")))+ggtitle(wilcox.test(IgE_breadth ~ key, data = x_treatment %>% select(ID,IgE_breadth,key) %>% filter(!is.na(IgE_breadth)) %>% filter(n()>=2), paired = TRUE)$p.value)
#treatment_IgE_breadth
ggplot(x_treatment,aes(x = key,y = IgG_breadth))+geom_violin(aes(fill = key)) + geom_point(size = 3)+geom_line(size = 2,aes(group = ID,color = difference)) +scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","grey","black"))+
  stat_compare_means(method  = "wilcox",comparisons = list(c("keep","Treatment")))+ggtitle(wilcox.test(IgG_breadth ~ key, data = x_treatment %>% select(ID,IgG_breadth,key) %>% filter(!is.na(IgG_breadth)) %>% filter(n()>=2), paired = TRUE)$p.value)
# placebo_IgE_Puro
ggplot(x_placebo,aes(x = key,y = Puro_IgE))+geom_violin(aes(fill = key)) + geom_point(size = 3)+geom_line(size = 2,aes(group = ID,color = difference)) +scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","black"))+
  stat_compare_means(method  = "wilcox",comparisons = list(c("keep","Placebo")))+ggtitle(wilcox.test(Puro_IgE ~ key, data = x_placebo %>% select(ID,Puro_IgE,key) %>% filter(!is.na(Puro_IgE)) %>% filter(n()>=2), paired = TRUE)$p.value)
#placebo_IgG_Puro
ggplot(x_placebo,aes(x = key,y = Puro_IgG))+geom_violin(aes(fill = key)) + geom_point(size = 3)+geom_line(size = 2,aes(group = ID,color = difference)) +scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","black"))+
  stat_compare_means(method  = "wilcox",mapping =aes(x= key,y = Puro_IgG), data = x_placebo %>% select(ID,Puro_IgG,key) %>% as.data.frame(),method.args = list(alternative = "two.sided",paired = TRUE),comparisons = list(c("keep","Placebo")))+ggtitle(wilcox.test(Puro_IgG ~ key, data = x_placebo %>% select(ID,Puro_IgG,key), paired = TRUE)$p.value)
# treament_IgE_Puro
ggplot(x_treatment,aes(x = key,y = Puro_IgE))+geom_violin(aes(fill = key)) + geom_point(size = 3)+geom_line(size = 2,aes(group = ID,color = difference)) +scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","grey","black"))+
  stat_compare_means(method  = "wilcox",comparisons = list(c("keep","Treatment")))+ggtitle(wilcox.test(Puro_IgE ~ key, data = x_treatment %>% select(ID,Puro_IgE,key) %>% filter(!is.na(Puro_IgE)) %>% filter(n()>=2), paired = TRUE)$p.value)
#treatment_IgG_Puro
ggplot(x_treatment,aes(x = key,y = Puro_IgG))+geom_violin(aes(fill = key)) + geom_point(size = 3)+geom_line(size = 2,aes(group = ID,color = difference)) +scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","grey","black"))+
  stat_compare_means(method  = "wilcox",comparisons = list(c("keep","Treatment")))+ggtitle(wilcox.test(Puro_IgG ~ key, data = x_treatment %>% select(ID,Puro_IgG,key) %>% filter(!is.na(Puro_IgG)) %>% filter(n()>=2), paired = TRUE)$p.value)



x_placebo %>% gather(x11,y11,-Name,-Time,-ID,-Baseline,-T2,-T3,-original,-Treatment,-New_dose,-key,-difference) %>% 
  ggplot(aes(x = key,y = y11))+geom_violin(aes(fill = key)) + 
  geom_point(size = 3)+geom_line(size = 1,aes(group = ID,color = difference)) +
  scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","black","grey"))+facet_wrap(~x11,nrow = 1)

x_treatment %>% gather(x11,y11,-Name,-Time,-ID,-Baseline,-T2,-T3,-original,-Treatment,-New_dose,-key,-difference) %>% 
  ggplot(aes(x = key,y = y11))+geom_violin(aes(fill = key)) + 
  geom_point(size = 3)+geom_line(size = 1,aes(group = ID,color = difference)) +
  scale_y_log10()+theme_classic()+scale_color_manual(values = c("red","grey","black"))+facet_wrap(~x11,nrow = 1)



# stat_compare_means(method  = "wilcox",comparisons = list(c("keep","Placebo")))+ggtitle(wilcox.test(IgE_breadth ~ key, data = x_placebo %>% select(ID,IgE_breadth,key) %>% filter(!is.na(IgE_breadth)) %>% filter(n()>=2), paired = TRUE)$p.value)



