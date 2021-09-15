" R script for Figure S3 in Estrela et al (2021) 
Functional attractors in microbial community assembly."

#  ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(gridExtra)

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
# ..Read data table 
dat.all<-read.csv("Estrela_2021_LCMS_metabol_EcoliEnterobacter.csv")
head(dat.all)

dat.ms<-select(dat.all,-value_blk,-value_uM, -Carbon_Source)

#. calculate mean of 2 reps by variable x  
ds<-dat.ms

nr<-4
stat_summ<- ds %>%
  dplyr::group_by(metabolite,Strain,Timepoint) %>%
  dplyr::summarise(
    Min = round(min(value_uM_corr, na.rm = TRUE),nr),
    Max = round(max(value_uM_corr, na.rm = TRUE),nr),
    mean = round(mean(value_uM_corr,na.rm = TRUE),nr),
    std = round(sd(value_uM_corr,na.rm = TRUE),nr),
    median = round(median(value_uM_corr, na.rm = TRUE),nr),
    IQRange = round(IQR(value_uM_corr, na.rm = TRUE),nr),
    Q25= round(quantile(value_uM_corr, 0.25, na.rm = TRUE),nr),
    Q75= round(quantile(value_uM_corr, 0.75, na.rm = TRUE),nr))
stat_summ

# .. Select compounds with values above 100uM and exclude glucose
dat.ms.s<-subset(as.data.frame(stat_summ), mean>100)
dp<-subset(dat.ms.s, metabolite !="Glucose")

# .rename species
dp$Strain<-revalue(dp$Strain, c('Ecoli'="E. coli"))
# .rename metabolites
dp$metabolite<-revalue(dp$metabolite, c('Acetic.acid'="Acetate", "Lactic.acid"='Lactate', "Succinic.acid"='Succinate', "Pyruvic.acid"="Pyruvate"))

#. re-order facets
dp$metabolite = factor(dp$metabolite, levels=c('Acetate','Lactate','Succinate','Pyruvate'))
dp$Strain = factor(dp$Strain, levels=c("Enterobacter", "E. coli"))

#. convert micromolar to millimolar
dp$mean_mM<-dp$mean/1000
dp$std_mM<-dp$std/1000

colp<- "deepskyblue4" 

ggplot(dp, aes(x=factor(metabolite), y=mean_mM, group=factor(Strain)))+
  geom_bar(stat="identity", fill=colp)+
  facet_grid(Strain~.)+
  geom_errorbar( aes(x=factor(metabolite), ymin=mean_mM-std_mM, ymax=mean_mM+std_mM), width=0.2, colour="black", alpha=0.9, size=1)+
  labs(y="concentration (mM)")+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.y =element_text(size=18),
    axis.title.x = element_blank(), 
    strip.text.y = element_text(size = 14),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line =element_blank())

ggsave(path=path.p, paste0("figS3.pdf"), width= 6, height=5.3, onefile = TRUE)

dat.figS3<-select(dp,metabolite, Strain, Timepoint, mean_mM, std_mM)

## ..Write data table as csv file
fname<-paste("figS3.csv")
write.csv(dat.figS3,file.path(path.ds, fname), row.names=FALSE)
