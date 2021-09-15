" R script for Figure S7 in Estrela et al (2021) 
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
library(plyr)
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
dat<-read.csv('Estrela_2021_C2R4_12h_48h_cfu.csv')

#.. Plot 48h (control) and 12h transfer at T10
dats<-subset(dat, Treatment %in% c('control', 'P_incub_time') & Transfer==10)

#.. Calculate mean and std for 2 reps
dp<- dats %>%
  dplyr::group_by(Transfer, Treatment,Community) %>%
  dplyr::mutate(
    mean_F_cfu = round(mean(F_cfu_ml),3),
    mean_R_cfu = round(mean(R_cfu_ml),3))
dp<-as.data.frame(select(dp,-rep_platting,-F_cfu_ml,-R_cfu_ml))

dp$t.ab<-dp$mean_F_cfu+dp$mean_R_cfu
dp$ab.ent<-dp$mean_F_cfu/dp$t.ab
dp$ab.pseud<-dp$mean_R_cfu/dp$t.ab

dp<-unique(select(dp,-mean_F_cfu, -mean_R_cfu, -t.ab))
dpm<-melt(dp, id.vars=c('Transfer','Treatment','Community'))

colEntbf<- "deepskyblue4"
colPseuf<-  "darkorchid2"

dpm$Treatment <-revalue(dpm$Treatment , c("P_incub_time"="12h incubation", "control"="48h incubation"))
dpm$Treatment = factor(dpm$Treatment, levels=c("12h incubation","48h incubation"))

  ggplot(dpm, aes(x=Community,y=value ,fill=variable))+
  geom_bar(stat="identity", position="stack", color="black") +  
  facet_grid(.~Treatment)+
  theme_minimal()+
  scale_fill_manual(values=c(colEntbf, colPseuf), name='', labels=c('Enterobacteriaceae', 'Pseudomonadaceae'))+
  scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(breaks= c(0, 1)) +
  labs(x='community', y='relative abundance')+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        strip.text = element_text(size=16),
        legend.title = element_text(size=12),
        axis.title.y =element_text(size=16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size=12),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

ggsave(path=path.p, "figS7.pdf", width=9, height=3.5, onefile = TRUE)

#.. write data figure as csv file
ds<-dpm
fname<-paste("figS7.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)



