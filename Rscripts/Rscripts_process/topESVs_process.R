" R script to extract sequence of dominant ESVs at T18
in Estrela et al (2021) Functional attractors in microbial community assembly."

#  ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(tidyr)
library(cowplot)
library(mixtools)
#------------------------------------------------------------------------------
#.. Assign path to data
#. usern = path to user directory

path.d<-paste0(usern, "/data/data_main")
#.. Assign path to plots output
path.p<- paste0(usern,"/Rplots")
#.. Assign path to data figures
path.ds<- paste0(usern,"/data/data_figures")
#.. Assign path to data processed
path.dp<- paste0(usern,"/data/data_processed")

#.. Read data 
setwd(path.d)
df.all<-read.csv('Estrela_2021_16S_data_table_RelAbund_ALL.csv')

df.all$experiment<-paste0(df.all$Treatment,'_', df.all$Inoc_size, '_', df.all$Transfer)
df.all$experiment_s<-paste0(df.all$Treatment,'_', df.all$Inoc_size)

# .. SELECT  treatments 
tt<-c("No_migration_4")
df.all.s<-subset(df.all, experiment_s %in% tt)
# 
df.all.s<-select(df.all.s,-Kingdom,-Phylum,-Class,-Order)
head(df.all.s)

#.. give a single ID to each OTU
d1<-as.data.frame(unique(df.all.s$OTU))
d2<-as.numeric(rownames(as.data.frame(unique(df.all.s$OTU))))
d3<-cbind(d1,d2)
names(d3)<-c('OTU', 'esv_id')
df.all.s<-merge(df.all.s, d3)
df.all.s$esv_id<-paste0(df.all.s$Genus, '_', df.all.s$esv_id)

## ..get abundance of dominant ESVs only
d.sum<-data.table(subset(df.all.s, Transfer==18))[, list(abund_sum=sum(Abundance)), 
                                                  by=list(experiment,Genus,OTU)]
d.sum<-as.data.frame(d.sum)

d.sum<-d.sum %>% 
  group_by(OTU) %>% 
  arrange(experiment, -abund_sum)
d.sum<-as.data.frame(d.sum)
d.sum
top6.esv<-select(d.sum[1:6,],-abund_sum)
## ..write as csv file the top 6 ESVs
fname<-paste("top6_ESVs_fig2B.csv")
write.csv(top6.esv,file.path(path.dp, fname),row.names=FALSE)



