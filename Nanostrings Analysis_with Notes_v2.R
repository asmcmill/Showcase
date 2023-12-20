#pull in libraries we are going to use
library(tidyverse)
library(pheatmap)
library(viridis)
library(ggforce)
library(janitor)
library(ggrepel)

#set working directory, this is where the raw data is.
setwd("C:/Users/asmcmill/Documents/Nanostrings")

#read data
DE<-read.csv("DE Results_Cecum_cdi vs cef.csv") # read in the DE output from Nanostrings, I call it DE

#list of nuclear receptors from Stephanie
NRs<-c("ABCC2","BAAT","CYP7B1","Fabp6","FGFR4","GPBAR1","KLB","LDLR","NR0B2","NR1I2","NR1I3","RXRA","SLC10A1","SLC10A2","SLC51a","SLC51b","SLC9A3","NR1I1","VLDLR","nr1h4","Cyp27a1","Cyp7a1","fgf19","fgf21","PPAR") #This is a list, a way to store data in R
NRsearch<-paste(NRs,collapse="|")#This makes a custom search by separating each in the list with a "|", which means "or". So we look for this or that or that or that (etc.). This is for the case insensitive search.

##Prepare the data for the other plots
DEplot<-DE%>% #start with DE
  mutate(X=str_remove(X,"-mRNA"))%>% #Probe names have "-mRNA" at the end, lets remove them for simplicity
  mutate(plot.p=-log10(P.value))%>% #transform the p value to -log10 p value, as that is the x axis of a volcano plot. We put this in a new column called "plot.p"
  mutate(direction=if_else(Log2.fold.change> 1 & P.value<0.05,"Up","Not DE"))%>% #if above one and significant, call it up, if not "Not DE". The fancier function for this is symnum() btw.
  mutate(direction=if_else(Log2.fold.change< -1 & P.value<0.05,"Down",direction))%>% #if below one and significant, call it down, if not keep what we had (So we don't overwrite the "Up"s)
  mutate(cholesterol=if_else(str_detect(Gene.sets,"Cholesterol Metabolism")& direction!="Not DE","Cholesterol Metabolism","Other"))%>% #Moving on to cholesterol. If we find "Cholesterol Metabolism" in the gene sets, and it's differentially expressed, flag it.
  mutate(NR=if_else(str_detect(X,regex(NRsearch,ignore_case=T)),"Nuclear Receptor","Other")) # Check if nuclear receptor. I had to search gene names without case sensitive using regular expressions (regex()) here.


####Red and Green Plot####
ggplot(DEplot,aes(x=Log2.fold.change,y=plot.p,fill=direction))+ #tell it to use the dataframe we made, x, y, and color based on direction
  geom_point(shape=21,size=3)+ #shape=21 gives a point with an outline (look up R point shapes on google for a chart), size sets the circle size
  geom_text_repel(data=DEplot%>%filter(direction!="Not DE"), #plot the names of receptors, remove names that are not DE
                  aes(x=Log2.fold.change, #tell it where to put the label
                      y=plot.p, #Add a buffer here to move it above the point
                      label=X),#tell it what to use as the labels
                  max.overlaps=20,#There's a lot of labels, this allows for them to touch each other more than default. Even with this 10 labels couldn't be plotted.
                  min.segment.length=0.2)+ #It was hard to tell which label went to which point, so I made it draw a line if it got too far away
  geom_hline(yintercept=-log10(0.05), #draw a horizontal line at the value I give it
             linetype="dashed")+# make the line dashed
  geom_vline(xintercept=1, #draw a vertical line at the first fold change cutoff
             linetype="dashed")+ #make it dashed
  geom_vline(xintercept=-1, #draw a vertical line at the second folc change cutoff
             linetype="dashed")+# make it dashed
  theme_bw()+ #This is the default theme I like to start with. You can look up other ones online, but this is a nice simple one.
  theme(panel.grid=element_blank(), #remove the gridlines in the back, I don't like them.
        legend.position="bottom")+ #move the location of the legend.
  labs(x="Log2Fold change",y="-log10padj", fill= "Direction")+ #This sets the axis and legend names (labs=labels)
  scale_y_continuous(limits=c(0,6), #This sets the min and max of the y axis
                     breaks=seq(0,6,1), #this sets the tick marks
                     expand=c(0,0))+ #this removes a "buffer zone" ggplot puts at either end of the axis limits I set.
  scale_x_continuous(limits=c(-6,6))+ #you don't have to set tick marks manually, here I just want to set limits on x axis.
  scale_fill_manual(breaks=c("Up","Down","Not DE"),values=c("red","green","grey")) #Only three colors here, so we don't need to use rainbow

ggsave("Red and Green Up and Down.pdf",height=5,width=8.5,units="in")#Save the plot. This one doesn't have a big legend so I decreased the height a bit.

####Cholesterol Plot####
ggplot(DEplot,aes(x=Log2.fold.change,y=plot.p,fill=direction))+ #same beginning setup
  geom_point(shape=21,size=3)+ #shape=21 gives a point with an outline (look up R point shapes on google for a chart), size sets the circle size
  geom_point(data=DEplot%>%filter(cholesterol=="Cholesterol Metabolism"),aes(fill=cholesterol),shape=21,size=3)+ #ggplot runs each line sequentially, this makes sure all the cholesterol points are plotted on top, change their color.
  geom_text_repel(data=DEplot%>%#plot the names of receptors
                    filter(cholesterol=="Cholesterol Metabolism"), #This one we want to only plot cholesterol so we filter for it
                  aes(x=Log2.fold.change, #tell it where to put the label
                      y=plot.p, #Add a buffer here to move it above the point
                      label=X),#tell it what to use as the labels (This is the column we removed not DE labels)
                  max.overlaps=20,#There's a lot of labels, this allows for them to touch each other more than default. Even with this 10 labels couldn't be plotted.
                  min.segment.length=0.2)+ #It was hard to tell which label went to which point, so I made it draw a line if it got too far away
  geom_hline(yintercept=-log10(0.05), #draw a horizontal line at the value I give it
             linetype="dashed")+# make the line dashed
  geom_vline(xintercept=1, #draw a vertical line at the first fold change cutoff
             linetype="dashed")+ #make it dashed
  geom_vline(xintercept=-1, #draw a vertical line at the second folc change cutoff
             linetype="dashed")+# make it dashed
  theme_bw()+ #This is the default theme I like to start with. You can look up other ones online, but this is a nice simple one.
  theme(panel.grid=element_blank(), #remove the gridlines in the back, I don't like them.
        legend.position="bottom")+ #move the location of the legend.
  labs(x="Log2Fold change",y="-log10padj", fill= "")+ #This sets the axis and legend names (labs=labels)
  scale_y_continuous(limits=c(0,6), #This sets the min and max of the y axis
                     breaks=seq(0,6,1), #this sets the tick marks
                     expand=c(0,0))+ #this removes a "buffer zone" ggplot puts at either end of the axis limits I set.
  scale_x_continuous(limits=c(-6,6))+ #you don't have to set tick marks manually, here I just want to set limits on x axis.
  scale_fill_manual(breaks=c("Cholesterol Metabolism","Up","Down","Not DE"),values=c("purple","red","green","grey")) #Here we order the legend using breaks and assign colors using values

ggsave("Cholesterol.pdf",height=5,width=8.5,units="in")#Save the plot. You can mess with the size

####Nuclear Receptor Plot####
ggplot(DEplot,aes(x=Log2.fold.change,y=plot.p,fill=direction))+ #again same base setup
  geom_point(shape=21,size=3)+ #shape=21 gives a point with an outline (look up R point shapes on google for a chart), size sets the circle size
  geom_point(data=DEplot%>%filter(NR=="Nuclear Receptor"),aes(fill=NR),shape=21,size=3)+ #ggplot runs each line sequentially, this makes sure all the nuclear receptor points are plotted on top. We'll also use aesthetics(aes()) to change the color.
  geom_label_repel(data=DEplot%>% #plot the names of receptors
                    filter(NR=="Nuclear Receptor"), #This one we want to only plot nuclear receptors so we filter for it
                  aes(x=Log2.fold.change, #tell it where to put the label
                      y=plot.p, #Add a buffer here to move it above the point
                      label=X),#tell it what to use as the labels (This is the column we removed not DE labels)
                  max.overlaps=20,#There's a lot of labels, this allows for them to touch each other more than default. Even with this 10 labels couldn't be plotted.
                  min.segment.length=0.2,#It was hard to tell which label went to which point, so I made it draw a line if it got too far away
                  size=3,
                  fill = ggplot2::alpha("grey88",0.7))+ #because these labels are on top of the mess at the bottom, I wanted to add a background to them. So I used geom_label_repel instead and filled the background with transparent grey.
  geom_hline(yintercept=-log10(0.05), #draw a horizontal line at the value I give it
             linetype="dashed")+# make the line dashed
  geom_vline(xintercept=1, #draw a vertical line at the first fold change cutoff
             linetype="dashed")+ #make it dashed
  geom_vline(xintercept=-1, #draw a vertical line at the second folc change cutoff
             linetype="dashed")+# make it dashed
  theme_bw()+ #This is the default theme I like to start with. You can look up other ones online, but this is a nice simple one.
  theme(panel.grid=element_blank(), #remove the gridlines in the back, I don't like them.
        legend.position="bottom")+ #move the location of the legend.
  labs(x="Log2Fold change",y="-log10padj", fill= "")+ #This sets the axis and legend names (labs=labels)
  scale_y_continuous(limits=c(0,6), #This sets the min and max of the y axis
                     breaks=seq(0,6,1), #this sets the tick marks
                     expand=c(0,0))+ #this removes a "buffer zone" ggplot puts at either end of the axis limits I set.
  scale_x_continuous(limits=c(-6,6))+ #you don't have to set tick marks manually, here I just want to set limits on x axis.
  scale_fill_manual(breaks=c("Nuclear Receptor","Up","Down","Not DE"),values=c("blue","red","green","grey")) #Here we order the legend using breaks and assign colors using values

ggsave("Nuclear Receptors.pdf",height=5,width=8.5,units="in")#Save the plot. You can mess with the size


##preparing to plot the pie chart, this one is a bit complicated
DEplot_pie<-DE%>% #start with DE
  mutate(X=str_remove(X,"-mRNA"))%>% #Probe names have "-mRNA" at the end, lets remove them for simplicity
  separate_rows(Gene.sets,sep=",")%>% # The gene sets are all in a list and separated with a comma. Lets make a new row for each gene set, this will duplicate the probes and their DE. We need to keep this in mind.
  mutate(plot.p=-log10(P.value))%>% #transform the p value to -log10 p value, as that is the x axis of a volcano plot. We put this in a new column called "plot.p"
  mutate(color=if_else(abs(Log2.fold.change)>1 & P.value<0.05, Gene.sets,"Not DE"))%>% #if_else: (Logical test, value if true, value if false). We make a new column called "color". If p>0.05 and absolute value of fold change is > 1, we record it's gene set, if not we call it "Not DE".
  mutate(color=if_else(color=="","No Gene Set",color))%>% #Now we check if it's not labelled with a gene set, if it isn't labelled we call it "No gene set" 
  mutate(color=str_trim(color,side="both"))%>% #There are some spaces before/after some of the Gene set labels, this will remove them so R doesn't think they're different categories.
  mutate(color=as.factor(color))%>% #Legends are organized by factor levels. So we convert the color column to a factor.
  mutate(color=fct_relevel(color,"Not DE",after=Inf))%>%# Not DE is in the middle, lets move it to the end of the legend by releveling the factor.
  group_by(X)%>% #this is for the pie chart. If it has the same gene name we group them together
  mutate(ratio=1/n()) #To get the size of the pie slices we divide 1 by the total number of gene sets assigned to a gene name (n() counts the number of rows in the group)

##multicolor spots with names
ggplot(data=DEplot_pie)+ #start ggplot object, use our plotting dataset as the data source
  geom_arc_bar(aes(x0=Log2.fold.change, #set the x value of the pie chart to log2fold change (This is complicated)
                   y0=plot.p, #set the y value of the pie chart to our log transformed p value
                   r0=0, #set where the radius starts
                   r=0.1, # set the length of the radius
                   amount=ratio, #tells it how big to make the slices
                   fill=color), #we color it based off the "gene set" or "not DE" column we made
               stat="pie", #this just tells it we're making pie charts
               col=NA)+ #this removes the outline of the pie charts, which was kind of distracting
  coord_equal()+ #ggplot will normally stretch axis to fit the plot area, I don't want it to because it will stretch my pie charts. Try removing this to see what I mean.
  geom_text_repel(data=DEplot_pie%>%
                    filter(color!="Not DE")%>% #remove names that are not DE
                    select(Log2.fold.change,plot.p,X)%>% #My name column has repeats, I don't want to write the gene name multiple times. So I pull out only variables I need.
                    unique(), #then I remove duplicates
                  aes(x=Log2.fold.change, #tell it where to put the label
                      y=plot.p, #Add a buffer here to move it above the point
                      label=X),#tell it what to use as the labels
                  max.overlaps=20,#There's a lot of labels, this allows for them to touch each other more than default. Even with this 4 labels couldn't be plotted.
                  min.segment.length=0.2)+ #It was hard to tell which label went to which point, so I made it draw a line if it got too far away
  geom_hline(yintercept=-log10(0.05), #draw a horizontal line at the value I give it
             linetype="dashed")+# make the line dashed
  geom_vline(xintercept=1, #draw a vertical line at the first fold change cutoff
             linetype="dashed")+ #make it dashed
  geom_vline(xintercept=-1, #draw a vertical line at the second folc change cutoff
             linetype="dashed")+# make it dashed
  labs(x="Log2Fold change",y="-log10padj", fill= "Gene Sets")+ #This sets the axis and legend names (labs=labels)
  theme_bw()+ #This is the default theme I like to start with. You can look up other ones online, but this is a nice simple one.
  theme(panel.grid=element_blank(), #remove the gridlines in the back, I don't like them.
        legend.position="bottom")+ #move the location of the legend.
  guides(fill=guide_legend(ncol=4))+ #this organizes the legend into 4 columns
  scale_y_continuous(limits=c(0,6), #This sets the min and max of the y axis
                     breaks=seq(0,6,1), #this sets the tick marks
                     expand=c(0,0))+ #this removes a "buffer zone" ggplot puts at either end of the axis limits I set.
  scale_x_continuous(limits=c(-6,6))+ #you don't have to set tick marks manually, here I just want to set limits on x axis.
  scale_fill_manual(values=c(rainbow(42), #Manually set the colors. There are 42 unique gene sets, so I use the function rainbow to give me 42 colors in the rainbow.
                             "grey")) #there is also the "Not DE" which we moved to the end and want to be grey, so I add grey to the end.
                             
ggsave("Colorful points.pdf",height=11,width=8.5,units="in")#Save the plot, call it "Colorful points", make it the size of a sheet of paper.


