library(tidyverse)
library(ggrepel)
library(xlsx)
library(RColorBrewer)
library(pheatmap)
library(stringr)

setwd("Q:\\My Drive\\Lab Work\\Data\\RNAseq\\R\\Final_Final")

###Figure 4#####
rvmtvp<-read.csv("Rich vs Minimal & TriKO vs Parental.csv",row.names = 1)
uni<-read.csv("Annotations/Uniprot.csv")%>%select(locus.tag,Protein.names)%>%rename(Uniprot=Protein.names)
pul<-read.csv("Annotations/PULS.csv")
KO<-read.csv("Annotations/kegg heirarchy+koala.csv", header=T, row.names=1)
geneindex<-read.csv("Annotations/Btheta_index_all.csv")%>%mutate(geneID=paste0(locus.tag,"_",protein))
volcdata<-rvmtvp%>%
  filter(!is.na(padj))%>%
  mutate(direction="Not DE")%>%
  mutate(direction=ifelse(padj<=0.05 & log2FoldChange>2,"Up",direction))%>%
  mutate(direction=ifelse(padj<=0.05 & log2FoldChange< -2,"Down",direction))%>%
  rownames_to_column("geneID")%>%
  left_join(.,geneindex,by="geneID")%>%
  left_join(.,KO,by="locus.tag")%>%
  left_join(.,uni,by="locus.tag")%>%
  left_join(.,pul,by="locus.tag")%>%
  separate(b, sep=" ", extra="merge",into=c("num","b"))%>%
  mutate(shape=ifelse(direction=="Not DE", "Open","Closed"))%>%
  select(-c(c,source,funct))%>%
  unique()%>%
  arrange(desc(shape))

write.csv(volcdata,"Figure 4_volcdata.csv")

counts<-volcdata%>%
  ungroup()%>%
  select(locus.tag,direction)%>%
  unique()%>%
  group_by(direction)%>%
  summarize(count=n())

write.csv(counts,"Figure 4_counts.csv")

kcounts<-volcdata%>%
  filter(a=="09100 Metabolism")%>%
  ungroup()%>%
  select(locus.tag,direction,b)%>%
  unique()%>%
  filter(direction!="Not DE")%>%
  group_by(b,direction)%>%
  summarize(count=n())%>%
  spread(key=direction,value=count)%>%
  replace(is.na(.),0)%>%
  as.data.frame()

write.csv(kcounts,"Figure 4_kcount.csv")

ggplot(volcdata,aes(x=log2FoldChange,y=-log10(padj)))+
  geom_vline(xintercept=c(-2, 2), col="grey1", lty="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey1", lty="dashed") +
  geom_point(aes(color=direction,size=2,shape=shape))+
  geom_point(data=filter(volcdata,direction!="Not DE"),size=5, shape=1)+
  theme_minimal()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.ticks=element_line(), axis.ticks.length=unit(.2,"cm"))+
  scale_color_manual(values= c("green","grey","red"))+
  scale_x_continuous(limits=c(-8,8),breaks=seq(-8,8,2))+
  scale_y_continuous(limits=c(0,14), expand=c(0,0),breaks=seq(-0,14,2))+
  scale_shape_manual(values=c(19,1,19))+
  scale_alpha_manual(values=c(1,.2,.2))+
  labs(title="Rich vs Minimal TriKO vs Parental",
       y= bquote("-log"[10]*"padj"),
       x= expression("log"[2]*"Fold change"),
       color="Legend")+
  theme(legend.position="right",
        legend.text = element_text(size=20),
        axis.text= element_text(size=18),
        axis.title = element_text(size=20),
        legend.title=element_blank(),
        panel.border = element_rect(colour = "darkgrey",fill=NA))+
  guides(color = guide_legend(override.aes = list(size = 4)),shape="none",size="none",alpha="none")

ggsave(paste0("Figure 4A.tiff"),width=9,height=10)

genelist<-filter(volcdata, direction!="Not DE")%>%
          #filter(b%in%c("Carbohydrate metabolism","Amino acid metabolism","Lipid metabolism","Energy metabolism"))%>%
          select(locus.tag)%>%
          unique()%>%
          as.list()

hmapdata<-read.csv("Rich vs Minimal strains Alldata.csv",row.names = 1)%>%
  select(geneID,log2FoldChange,strain)%>%
  group_by(geneID)%>%
  spread(key=strain,value=log2FoldChange)%>%
  left_join(.,geneindex,by="geneID")%>%
  left_join(.,KO,by="locus.tag")%>%
  left_join(.,uni,by="locus.tag")%>%
  left_join(.,pul,by="locus.tag")%>%
  separate(b, sep=" ", extra="merge",into=c("num","b"))%>%
  filter(locus.tag %in% genelist$locus.tag)%>%
  select(geneID,Parental,TriKO,b)%>%
  filter(!is.na(b))%>%
  #filter(b %in% c("Carbohydrate metabolism","Amino acid metabolism","Lipid metabolism","Energy metabolism"))%>%
  unique()%>%
  #mutate(b=factor(b,levels=c("Carbohydrate metabolism","Amino acid metabolism","Lipid metabolism","Energy metabolism")))%>%
  group_by(geneID)%>%
  mutate(direction=ifelse(Parental<TriKO,1,2))%>%
  arrange(b,direction)%>%
  mutate(index=1:n())%>%
  mutate(geneID=paste0(geneID,"_",index))%>%
  select(-index,-direction)%>%
  column_to_rownames("geneID")%>%
  rename("WT"="Parental","triple KO"="TriKO")

anno<-hmapdata%>%select(b)
palette<- colorRampPalette(colors = c("green", "black", "red"))(n = 50)
breaks<-(seq(-6,6,by=(12/49)))

pheatmap(hmapdata[1:2],
         annotation_row=anno,
         #labels_row="geneID",
         nacol="grey",
         main= "Rich vs Minimal & Parental vs TriKO",
         filename=paste0("Figure 4C.tiff"),
         border_color = NA,
         cluster_rows=F,
         cluster_cols=F,
         width=10,
         height=10,
         color=palette,
         breaks=breaks,
         angle_col = 0)
  

###Figure 5#####
bvnbrvm<-read.csv("BA vs No_BA Rich vs Minimal.csv",row.names=1,header=T)
uni<-read.csv("Annotations/Uniprot.csv")%>%select(locus.tag,Protein.names)%>%rename(Uniprot=Protein.names)
pul<-read.csv("Annotations/PULS.csv")
KO<-read.csv("Annotations/kegg heirarchy+koala.csv", header=T, row.names=1)
geneindex<-read.csv("Annotations/Btheta_index_all.csv")%>%mutate(geneID=paste0(locus.tag,"_",protein))
volcdata<-bvnbrvm%>%
  filter(!is.na(padj))%>%
  mutate(direction="Not DE")%>%
  mutate(direction=ifelse(padj<=0.05 & log2FoldChange>2,"Up",direction))%>%
  mutate(direction=ifelse(padj<=0.05 & log2FoldChange< -2,"Down",direction))%>%
  left_join(.,geneindex,by="geneID")%>%
  left_join(.,KO,by="locus.tag")%>%
  left_join(.,uni,by="locus.tag")%>%
  left_join(.,pul,by="locus.tag")%>%
  mutate(direction=fct_relevel(direction, c("Up","Down","Not DE")))%>%
  separate(b, sep=" ", extra="merge",into=c("num","b"))%>%
  mutate(color=ifelse(direction=="Not DE","NA", ifelse(is.na(PUL),"NA","PUL")))%>%
  mutate(color=ifelse(direction=="Not DE",color,ifelse(a%in% c("09100 Metabolism","09140 Cellular Processes","09120 Genetic Information Processing","09130 Environmental Information Processing"),b,color)))%>%
  mutate(shape=ifelse(direction=="Not DE", "Not DE",ifelse(is.na(color),"Not in KEGG","KEGG Annotated")))%>%
  mutate(color=fct_relevel(color,c("PUL",
                                   "Carbohydrate metabolism",
                                   "Metabolism of cofactors and vitamins",
                                   "Membrane transport",
                                   "Lipid metabolism",
                                   "Amino acid metabolism",     
                                   "Cellular community - prokaryotes", 
                                   "Energy metabolism",      
                                   "Folding, sorting and degradation",   
                                   "Glycan biosynthesis and metabolism",
                                   "Nucleotide metabolism",  
                                   "Xenobiotics biodegradation and metabolism",
                                   "NA")))%>%
  select(-c(c,source,funct))%>%
  unique()%>%
  arrange(!is.na(color),desc(color))

write.csv(volcdata,"Figure 5.csv")

counts<-volcdata%>%
  ungroup()%>%
  select(locus.tag,direction,ba)%>%
  unique()%>%
  group_by(direction,ba)%>%
  summarize(count=n())%>%
  spread(key=direction,value=count)

write.csv(counts,"Figure 5_counts.csv")

kcounts<-volcdata%>%
  ungroup()%>%
  select(locus.tag,direction,b,ba)%>%
  unique()%>%
  filter(direction!="Not DE")%>%
  group_by(b,direction,ba)%>%
  summarize(count=n())%>%
  spread(key=direction,value=count)%>%
  as.data.frame()%>%
  mutate(ba=fct_relevel(ba,c("DCA","GDCA","TDCA","CDCA","GCDCA","TCDCA","CA","GCA","TCA")))%>%
  arrange(b,ba)

write.csv(kcounts,"Figure 5_kcounts.csv")

pulcounts<-volcdata%>%
  ungroup()%>%
  select(locus.tag,direction,PUL,ba)%>%
  mutate(pulcheck=is.na(PUL))%>%
  filter(pulcheck==F)%>%
  select(-PUL,-pulcheck)%>%
  unique()%>%
  filter(direction!="Not DE")%>%
  group_by(direction,ba)%>%
  summarize(count=n())%>%
  spread(key=direction,value=count)%>%
  replace(is.na(.),0)%>%
  as.data.frame()%>%
  mutate(ba=fct_relevel(ba,c("DCA","GDCA","TDCA","CDCA","GCDCA","TCDCA","CA","GCA","TCA")))%>%
  arrange(ba)

write.csv(pulcounts,"Figure 5_pulcounts.csv")

ggplot(volcdata,aes(x=log2FoldChange,y=-log10(padj)))+
  geom_vline(xintercept=c(-2, 2), col="grey1", lty="dashed") +
  geom_hline(yintercept=-log10(0.05), col="grey1", lty="dashed") +
  geom_point(aes(color=color,shape=shape,alpha=shape))+
  geom_point(data=filter(volcdata,direction!="Not DE"), shape=1)+
  facet_wrap(~ba,ncol=3)+
  theme_minimal()+
  scale_color_manual(values=c("orchid","purple","red","blue","yellow","green","chocolate4","orange","peachpuff","thistle", "azure","turquoise","grey"))+
  scale_shape_manual(values=c(19,1,2))+
  scale_alpha_manual(values=c(1,.2,.2))+
  labs(#title="Bile Acid vs No Bile Acid Rich vs Minimal",
       y= bquote("-log"[10]*"padj"),
       x= expression("log"[2]*"Fold change"),
       color="Legend")+
  theme(legend.position="right",
        legend.text = element_text(size=10),
        legend.title=element_blank(),
        panel.border = element_rect(colour = "darkgrey",fill=NA))+
  guides(color = guide_legend(override.aes = list(size = 3)),size="none",alpha="none")

ggsave(paste0("Supplemental BA.tiff"),width=9,height=10)


palette<- colorRampPalette(colors = c("green", "black", "red"))(n = 50)
breaks<-(seq(-6,6,by=(12/49)))

hmapdata<-bvnbrvm%>%
  merge(.,geneindex,by="geneID",all=T)%>%
  left_join(.,KO,by="locus.tag", relationship="many-to-many")%>%
  select(geneID,log2FoldChange,ba,c)%>%
  filter(c %in% c("00430 Taurine and hypotaurine metabolism [PATH:bth00430]","00260 Glycine, serine and threonine metabolism [PATH:bth00260]"))%>%
  unique()%>%
  group_by(ba,c)%>%
  spread(ba,log2FoldChange)%>%
  column_to_rownames("geneID")%>%
  select(8,5,2,9,6,3,10,7,4,1)%>%
  arrange(c)%>%
  mutate(c=str_remove(c,"\\[(.*?)\\]"))

anno<-hmapdata%>%select(c)
pheatmap(hmapdata[1:9],
         annotation_row=anno,
         #labels_row="geneID",
         nacol="grey",
         #main= "bile acid genes",
         filename=paste0("Supp_9 Taurine and Glycine.tiff"),
         border_color = NA,
         cluster_rows=F,
         cluster_cols=F,
         width=22,
         height=8,
         color=palette,
         breaks=breaks,
         angle_col = 0,
         fontsize=14)

pulanno<-read.csv("Annotations/PULs_type.csv")

pulgeneshmap<-volcdata%>%
  select(log2FoldChange,ba,direction,PUL,geneID)%>%
  filter(direction!="Not DE" & PUL!="")%>%
  unique()%>%
  #mutate(geneID=paste0(geneID,"_(",Modularity,")"))%>%
  select(-direction)%>%
  left_join(.,pulanno, by="PUL")%>%
  group_by(ba,PUL,Target.substrate)%>%
  mutate(PUL=as.integer(str_sub(PUL,4,99)))%>%
  filter(!is.na(PUL))%>%
  #summarize(average=mean(log2FoldChange))%>%
  arrange(PUL)%>%
  group_by(ba)%>%
  spread(ba,log2FoldChange)%>%
  arrange(Target.substrate,PUL)%>%
  mutate(PUL=paste0("PUL",PUL))%>%
  column_to_rownames("geneID")%>%
  select(9,6,3,10,7,4,11,8,5,1,2)

pulhmap<-volcdata%>%
  select(log2FoldChange,ba,direction,PUL,geneID)%>%
  filter(direction!="Not DE" & PUL!="")%>%
  unique()%>%
  select(-direction)%>%
  left_join(.,pulanno, by="PUL")%>%
  group_by(ba,PUL,Target.substrate)%>%
  mutate(PUL=as.integer(str_sub(PUL,4,99)))%>%
  filter(!is.na(PUL))%>%
  summarize(average=mean(log2FoldChange))%>%
  arrange(PUL)%>%
  group_by(ba)%>%
  spread(ba,average)%>%
  arrange(Target.substrate,PUL)%>%
  mutate(PUL=paste0("PUL",PUL))%>%
  column_to_rownames("PUL")%>%
  select(8,5,2,9,6,3,10,7,4,1)

palette<- colorRampPalette(colors = c("green", "black", "red"))(n = 50)
breaks<-(seq(-8,8,by=(16/49)))

anno_color<-list(
  Target.substrate=c("a mannan, host N-glycans"="red",
      "host glycans"  ="orange",
      "host/residual dietary glycans"="yellow",
      "mucin O-glycans"="green",
      "rhamnogalacturonan I"="blue",
      "rhamnogalacturonan II"="lightblue",
      "starch" ="purple",
      "unknown"="grey")
)


anno<-pulgeneshmap%>%select(Target.substrate,PUL)%>%mutate(PUL=factor(PUL,unique(PUL)))
t<-pheatmap(pulgeneshmap[1:9],
            annotation_row=anno,
            annotation_colors = anno_color,
            annotation_names_row=F,
            #labels_row="geneID",
            nacol="grey",
            #main= "bile acid genes",
            filename=paste0("pulhmap_allgenes_bonferroni.tiff"),
            border_color = "grey",
            cluster_rows=F,
            cluster_cols=F,
            width=12,
            height=10,
            color=palette,
            breaks=breaks,
            angle_col = 0)

anno<-pulhmap%>%select(Target.substrate)
t<-pheatmap(pulhmap[1:9],
            annotation_row=anno,
            annotation_colors = anno_color,
            annotation_names_row=F,
            #labels_row="geneID",
            nacol="grey",
            #main= "bile acid genes",
            filename=paste0("pulhmap_bonferroni.tiff"),
            border_color = "grey",
            cluster_rows=F,
            cluster_cols=F,
            width=10,
            height=8,
            color=palette,
            breaks=breaks,
            angle_col = 0)

allhmap<-bvnbrvm%>%
  filter(!is.na(padj))%>%
  mutate(direction="Not DE")%>%
  mutate(direction=ifelse(padj<=0.05 & log2FoldChange>2,"Up",direction))%>%
  mutate(direction=ifelse(padj<=0.05 & log2FoldChange< -2,"Down",direction))%>%
  left_join(.,geneindex,by="geneID")%>%
  left_join(.,KO,by="locus.tag",relationship = "many-to-many")%>%
  left_join(.,uni,by="locus.tag",relationship = "many-to-many")%>%
  left_join(.,pul,by="locus.tag",relationship = "many-to-many")%>%
  mutate(direction=fct_relevel(direction, c("Up","Down","Not DE")))%>%
  separate(b, sep=" ", extra="merge",into=c("num","b"))%>%
  separate(c, sep=" ", extra="merge",into=c("num","c"))%>%
  filter(b%in% c("Carbohydrate metabolism", "Lipid metabolism","Amino acid metabolism","Membrane transport","Metabolism of cofactors and vitamins") & !is.na(c) & direction!="Not DE")%>%
  select(log2FoldChange,ba,c,geneID,b)%>%
  unique()%>%
  group_by(ba,b,c,geneID)%>%
  summarize(average=mean(log2FoldChange))%>%
  spread(ba,average)%>%
  arrange(c)%>%
  group_by(geneID)%>%
  mutate(index=1:n())%>%
  mutate(geneID=paste0(geneID,"_",index))%>%
  select(-index)%>%
  column_to_rownames("geneID")%>%
  select(9,6,3,10,7,4,11,8,5,1,2)%>%
  mutate(b=fct_relevel(b, c("Carbohydrate metabolism", "Metabolism of cofactors and vitamins","Membrane transport","Lipid metabolism","Amino acid metabolism")))%>%
  arrange(b)%>%
  mutate(c=str_remove(c,"\\[(.*?)\\]"))

palette<- colorRampPalette(colors = c("green", "black", "red"))(n = 50)
breaks<-(seq(-6,6,by=(12/49)))
anno<-allhmap%>%
  select(c,b)%>%
  mutate(c=factor(c,unique(c)))%>%
  mutate(b=factor(b,unique(b)))

anno_color<-list(
  b=c("Carbohydrate metabolism"="purple",
      "Metabolism of cofactors and vitamins"="red",
      "Membrane transport"="blue",
      "Lipid metabolism"="yellow",
      "Amino acid metabolism"="green")
)


t<-pheatmap(allhmap[1:9],
            annotation_row=anno,
            annotation_colors = anno_color,
            #labels_row="geneID",
            nacol="grey",
            #main= "bile acid genes",
            filename="kegg_genes_hmap.tiff",
            border_color = "grey",
            cluster_rows=F,
            cluster_cols=F,
            width=15,
            height=10,
            color=palette,
            breaks=breaks,
            angle_col = 0)
