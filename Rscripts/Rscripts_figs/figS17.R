" R script for Figures S17 in Estrela et al (2021) 
Functional attractors in microbial community assembly."

#  ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

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
dat<-read.csv('Estrela_2021_fba_O2_limitation.csv')
dat<-unique(dat)

#. re-order CS
dat$Carbon_Source = factor(dat$Carbon_Source, 
                           levels=c('Glucose', 'Acetate','Lactate', 'Succinate'))

p<-ggplot(dat, aes(x=Oxygen, y=Proportion*100, group=factor(Family):factor(Oxygen):factor(Carbon_Source), fill=Oxygen))+
  geom_bar(color="black", stat="identity", position="stack") + 
  facet_grid(Family~Carbon_Source)+
  scale_fill_manual(values=c('gray70', 'gray30'))+ 
  labs(y= "% models that grow", x="")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
p 

ggsave(path=path.p, paste0("figS17.pdf"), width= 8, height=4.5, onefile = TRUE)

dat.figS17<-dat

## ..Write data table as csv file
fname<-paste("figS17.csv")
write.csv(dat.figS17,file.path(path.ds, fname), row.names=FALSE)



