" R script for Figure S20B (+stats) in Estrela et al (2021) 
Functional attractors in microbial community assembly."
#------------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(reshape2)

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

setwd(path.dp)

# ..Read data table for mock community (note: only taxa with abundance >0.01 are shown)
dat.mock<-read.csv("figS20_mockcommunity_inhouse_16S_datatable.csv")

d.m<-select(dat.mock,Genus,Abundance,Sample)
d.m$type<-'16S'
d.m$SampleID<-paste0('16S_',d.m$Sample)
d.m<-select(d.m,-Sample)

#. Add theoretical (i.e. OD mixing) info 
Genus<-c('Klebsiella', 'Pseudomonas', 'Raoultella', 'Pseudomonas', 'Aeromonas')
strain<-c('Klebsiella.225', 'Pseudomonas.226', 'Raoultella.165', 'Pseudomonas.162', 'Aeromonas.343')
Abundance<-c(0.2,0.2,0.2,0.2,0.2)
mixed.seq<-data.frame(Genus, Abundance)
mixed.seq$type<-'known composition'
mixed.seq$SampleID<-'known composition'

#. combine both
d.16s.od<-rbind(d.m,mixed.seq)

dp<-d.16s.od
dp$Genus <- factor(dp$Genus, levels = c('Aeromonas', 'Klebsiella', 'Raoultella', 'Pseudomonas'))
dp$Sample[dp$SampleID=='16S_0280_s0062']<-'16S_1'
dp$Sample[dp$SampleID=='16S_0472_s0063']<-'16S_2'
dp$Sample[dp$SampleID=='16S_0664_s0063']<-'16S_3'
dp$Sample[dp$SampleID=='16S_0832_s0064']<-'16S_4'
dp$Sample[dp$SampleID=='16S_0088_s0062']<-'16S_5'

#. Plot cell mock community and OD mock community
#. Genus colors 
colAero<-"deepskyblue4"
colK<-blues9[7] 
colPseu<-  "darkorchid3"
colRao<-"cornflowerblue" # 

colgen<-c(colAero, colK,colRao,colPseu)
p<-ggplot(dp, aes(x=factor(Sample), y=Abundance, fill=Genus))+
 geom_bar(color="black", stat="identity", position="stack") + 
  ylab("Relative abundance")+
  xlab("Sample_ID")+
  theme_classic()+
  scale_x_discrete(labels=c('16S_1','16S_2','16S_3','16S_4','16S_5','known 
composition'))+
  theme(panel.spacing = unit(0, "lines"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x=element_blank(),
        legend.text=element_text(size=12),
        legend.position="bottom")+
  guides(fill = guide_legend(ncol = 10))+
  scale_fill_manual(values=colgen)+
  ggtitle('Cell mock community')
p

ggarrange(p, ncol=1,nrow=1, labels = c("B"))
ggsave(path=path.p, paste0("figS20B.pdf"), width= 6.2, height=4, onefile = TRUE)

dat.figS20B<-select(dp,-SampleID)

## ..Write data table as csv file
fname<-paste("figS20B.csv")
write.csv(dat.figS20B,file.path(path.ds, fname), row.names=FALSE)

#.. Create new column with F or R assignment
d.16s.od$f_r_type <- with(d.16s.od, 
                        ifelse(Genus == "Aeromonas" | Genus == "Klebsiella" | Genus=="Raoultella" ,
                               "F", "R"))
head(d.16s.od)

#.. sum the relative abundance of F or R for each sample
dsum_fr<-data.table(d.16s.od)[, list(f_r_sum=sum(Abundance)), 
                            by=list(f_r_type, SampleID, type)]

#.. calculate ratio of R/F for each sample
dsum_fr2 <- dsum_fr %>% 
  spread(f_r_type, f_r_sum) %>%
  group_by(SampleID) %>% 
  mutate(ratioRF = R/F)
dsum_fr2<-as.data.frame(dsum_fr2)
dsum_fr2

#. calculate RF correction factor
ds.16.obs<-subset(dsum_fr2, SampleID !='known composition')
ds.16.theor<-subset(dsum_fr2, SampleID =='known composition')

rf_corr<-mean(ds.16.theor$ratioRF)/mean(ds.16.obs$ratioRF)
rf_corr



