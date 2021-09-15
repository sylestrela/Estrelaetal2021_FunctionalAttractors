" R script for Figure S20A in Estrela et al (2021) 
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
# ..Read R/F ratio for Emergent Simplicity 16S data
dat.16s<-read.csv('Goldford_2018_glu_ratioRF.csv')

#.. Calculate predicted R/F based on correction value obtained from regression in Fig. S19 (i.e. slope of fit) 
d.16s<-dat.16s
correc.f<-0.455
d.16s$ratioRF_corr<-d.16s$ratioRF*correc.f

#.. Calculate predicted R/F based on correction value from mock community (see companion script for FigS20B)
rf.corr.mock<-1.141296
d.16s$ratioRF_corr_mock<-d.16s$ratioRF*rf.corr.mock

#.. Calculate median and IQR for 3 scenarios separately
d.16s.m<-melt(data = d.16s, measure.vars = c("ratioRF", "ratioRF_corr", "ratioRF_corr_mock"))

ds<-d.16s.m
nr<-4
stat_summ<- ds %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    Min = round(min(value, na.rm = TRUE),nr),
    Max = round(max(value, na.rm = TRUE),nr),
    mean = round(mean(value,na.rm = TRUE),nr),
    std = round(sd(value,na.rm = TRUE),nr),
    median = round(median(value, na.rm = TRUE),nr),
    IQRange = round(IQR(value, na.rm = TRUE),nr),
    Q25= round(quantile(value, 0.25, na.rm = TRUE),nr),
    Q75= round(quantile(value, 0.75, na.rm = TRUE),nr))
stat_summ<-as.data.frame(stat_summ)

#.. plot
colp<-'darkorchid4'
dp<-ds
xl<-""
yl<-"R/F"
p<-ggplot(dp, aes(x=variable, y=value))+
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, aes(shape=variable, colour=variable), size=2, alpha=0.8)+
  theme_classic()+
  scale_shape_manual(values = c(1, 2, 3,17)) +
  scale_colour_manual(values = c(colp,colp, colp,colp))+
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.y =element_text(size=20),
    axis.title.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none")+
  scale_x_discrete( 
    labels=c("16S observed", "CFU predicted", "16S corrected"))+
  scale_y_log10(breaks = c(0.01, 0.3,1,10),limits = c(0.001, 10))+
  labs(y=yl)
p

ggarrange(p, ncol=1,nrow=1, labels = c("A"))
ggsave(path=path.p, paste0("figS20A.pdf"), width= 6, height=4, onefile = TRUE)

dat.figS20A<-select(dp,InocRep, variable,value)

## ..Write data table as csv file
fname<-paste("figS20A.csv")
write.csv(dat.figS20A,file.path(path.ds, fname),row.names=FALSE)


