library(tidyverse)
library(tximport)
library(DESeq2)

## metadata and directories
setwd("Q:\\My Drive\\Lab Work\\Data\\RNAseq\\R\\Final_Final")
salmonpath <- "Q:\\My Drive\\Lab Work\\Data\\RNAseq\\Salmon"

##Sample manifest information
key <- read.csv("manifest.csv", header = TRUE)%>%
  mutate(filepath=paste0(salmonpath,"/",filename))%>%
  mutate(quantpath=paste0(salmonpath,"/",filename, "/quant.sf"))
locustagindex <- read.csv("Btheta_tx2gene.csv", header = TRUE)
tx2gene<-locustagindex %>%
  select(1,2)

##Gather list of files from sample list, make sure directory and filenames have been set up properly
files <- file.path(key$filepath, "quant.sf")
#check for files
ifelse(all(file.exists(files)), paste("All ",length(files),  " files found"), paste("Problem with manifest, check folder names or quant.sf files"))

######DESeq Parental vs TriKO No BA####

samples <- key %>%
  filter(Condition=="No_BA")%>%
  mutate(Strain=fct_relevel(Strain, c("TriKO","Parental")))%>%
  mutate(Media=fct_relevel(Media, c("Rich","Minimal")))
txi <- tximport(samples$quantpath, type = "salmon", tx2gene = tx2gene)
#Design dds, get metadata from the "samples" dataframe, remove file location information
dds <- DESeqDataSetFromTximport(txi, select(samples,-c(filepath,filename,quantpath)), ~Strain + Media +Strain:Media)
#Remove counts <10
dds <- dds[rowSums(counts(dds)) >= 10,]
#Differential Expression
dds <- DESeq(dds)

#Rich vs Minimal 
#Pulling the results. p<0.05. Fold change >1
res <- results(dds, contrast = c("Media","Minimal","Rich"), alpha = 0.05, lfcThreshold = 1, altHypothesis = "greaterAbs")

#combine data from run
alldata<-as.data.frame(res)%>%
  tibble::rownames_to_column(var="geneID")%>%
  mutate(strain="TriKO")

#Pulling the results for the next strain. p<0.05. Fold change >1
res <- results(dds, contrast = list(c("Media_Minimal_vs_Rich","StrainParental.MediaMinimal")), alpha = 0.05, lfcThreshold = 1, altHypothesis = "greaterAbs")

#combine data from run into alldata
newdata<-as.data.frame(res)%>%
  tibble::rownames_to_column(var="geneID")%>%
  mutate(strain="Parental")
alldata<-rbind(alldata,newdata)

write.csv(alldata, file = paste0("Rich vs Minimal strains Alldata.csv"))

##(Rich vs Minimal) vs (TriKO vs Parental)##
#is the condition effect different across genotype?
res <- results(dds, name="StrainParental.MediaMinimal")

write.csv(as.data.frame(res), file = paste0("Rich vs Minimal & TriKO vs Parental.csv"))

#Average Abundance
abundance<-as.data.frame(txi$abundance)%>%
  rename_at(vars(c(1:12)),~c(samples$filename))%>%
  rownames_to_column(var="geneID")%>%
  group_by(geneID)%>%
  gather(key=filename,value=abundance,-geneID)%>%
  left_join(.,key%>%select(filename,ID),by="filename")%>%
  select(-filename)%>%
  group_by(geneID,ID)%>%
  summarize(meantpm=mean(abundance))%>%
  spread(key=ID, value=meantpm)

write.csv(abundance,"Parental&TriKO_abundance.csv")
#####
######DESeq Rich vs Minimal & BA vs No BA
##DESeq No BA vs BA Rich vs Minimal

#Final interpretation: Genes positively associated with nutrient limitation due to bile acid 
samples <- key %>% filter(Strain=="Parental")%>%
  mutate(Condition=fct_relevel(Condition,c("No_BA","TCA","GCA","CA","TCDCA","GCDCA","CDCA","TDCA","GDCA","DCA")))%>%
  mutate(Media=fct_relevel(Media, c("Rich","Minimal")))
txi <- tximport(samples$quantpath, type = "salmon", tx2gene = tx2gene)
#Design dds, get metadata from the "samples" dataframe, remove file location information
dds <- DESeqDataSetFromTximport(txi, select(samples,-c(filepath,filename,quantpath)), ~Condition+Media+Condition:Media)
#Remove counts <10
dds <- dds[rowSums(counts(dds)) >= 10,]
#Differential Expression
dds <- DESeq(dds)

#start alldata
alldata <-NULL

#make long results
for (ba in c("TCA","GCA","CA","TCDCA","GCDCA","CDCA","TDCA","GDCA","DCA")){
  res <- results(dds, name=paste0("Condition",ba,".MediaMinimal"),independentFiltering=F,pAdjustMethod="bonferroni")%>%
    as.data.frame()%>%
    mutate(ba=ba)%>%
    rownames_to_column("geneID")
  alldata<-rbind(alldata,res)
}

#Save data as csv
write.csv(alldata, file = paste0("BA vs No_BA Rich vs Minimal.csv"))

#start alldata for BA vs No BA
alldata <-NULL

#make long results  for BA vs No BA
for (ba in c("TCA","GCA","CA","TCDCA","GCDCA","CDCA","TDCA","GDCA","DCA")){
  res <- results(dds, name=paste0("Condition_",ba,"_vs_No_BA"),independentFiltering=F,pAdjustMethod="bonferroni")%>%
    as.data.frame()%>%
    mutate(ba=ba)%>%
    mutate(Media="Rich")%>%
    rownames_to_column("geneID")
  alldata<-rbind(alldata,res)
  
  res <- results(dds, list(c(paste0("Condition_",ba,"_vs_No_BA"),paste0("Condition",ba,".MediaMinimal"))),independentFiltering=F,pAdjustMethod="bonferroni")%>%
    as.data.frame()%>%
    mutate(ba=ba)%>%
    mutate(Media="Minimal")%>%
  rownames_to_column("geneID")
  alldata<-rbind(alldata,res)
}

#Save data as csv
write.csv(alldata, file = paste0("BA vs No_BA Rich and Minimal.csv"))

#Average Abundance
abundance<-as.data.frame(txi$abundance)%>%
  rename_at(vars(c(1:60)),~c(samples$filename))%>%
  rownames_to_column(var="geneID")%>%
  group_by(geneID)%>%
  gather(key=filename,value=abundance,-geneID)%>%
  left_join(.,key%>%select(filename,ID),by="filename")%>%
  select(-filename)%>%
  group_by(geneID,ID)%>%
  summarize(meantpm=mean(abundance))%>%
  spread(key=ID, value=meantpm)

write.csv(abundance, "BA & No_BA Rich & Minimal_abundance.csv")

