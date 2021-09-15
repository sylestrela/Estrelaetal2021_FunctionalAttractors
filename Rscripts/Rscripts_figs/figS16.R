" R script for Figure S16 in Estrela et al (2021) 
Functional attractors in microbial community assembly."

#  ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(reshape2)
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
# ..Read data table for normal fba data in Fig 1D
d.fba<-read.csv('Estrela_2021_fba_params_Pairs_Acetate.csv')
dat.fba<-select(d.fba, -E,-P,-E_Genus,-P_Species,-R_F)
dat.fba$ratioRF<-dat.fba$D_ace_glu*dat.fba$w_ac_R/dat.fba$w_glc_F
dat.fba$Experiment<-'fba'
dat.fba<-select(dat.fba, -D_ace_glu, -w_ac_R,-w_glc_F)

setwd(path.dp)
# ..Read data table for constrained fba data
d.fba.c<-read.csv('figS16_ratioPE_fba_Constrained_Secretions.csv')
dat.fba.c<-select(d.fba.c, -Sanger_ID)
dat.fba.c<-melt(dat.fba.c, id.vars='Condition', measure.vars ='ratioPE')
colnames(dat.fba.c)[colnames(dat.fba.c)=="value"] <- "ratioRF"
colnames(dat.fba.c)[colnames(dat.fba.c)=="Condition"] <- "Experiment"
dat.fba.c<-select(dat.fba.c,-variable)
  
# ..Read data table for RF ratio of ES communities
dat.exp <- read.csv(paste0(path.dp,"/Goldford_2018_glu_ratioRF.csv"))
dat.exp$Experiment<-"Exp."
dat.exp<-select(dat.exp,-InocRep,-F,-R)

# .merge 3 dataframes
dat.ep.m<-rbind(dat.exp,dat.fba, dat.fba.c)

# .. Colour scheme 
# Family colors 
colEntbf<- "deepskyblue3" #
colPseuf<-  "darkorchid2" # 

# .. Calculate median and IQR for FBA and for experiments
ds<-dat.ep.m
nr<-3
stat_summ<- ds %>%
  dplyr::group_by(Experiment) %>%
  dplyr::summarise(
    Min = round(min(ratioRF, na.rm = TRUE),nr),
    Max = round(max(ratioRF, na.rm = TRUE),nr),
    mean = round(mean(ratioRF,na.rm = TRUE),nr),
    std = round(sd(ratioRF,na.rm = TRUE),nr),
    median = round(median(ratioRF, na.rm = TRUE),nr),
    IQRange = round(IQR(ratioRF, na.rm = TRUE),nr),
    Q25= round(quantile(ratioRF, 0.25, na.rm = TRUE),nr),
    Q75= round(quantile(ratioRF, 0.75, na.rm = TRUE),nr))
stat_summ<-as.data.frame(stat_summ)

stat.fba<-subset(stat_summ, Experiment=="fba")
stat.exp<-subset(stat_summ, Experiment=="Exp.")

colp<-"tan1"
colp<-"darkorchid4"

#. re-order factors
dat.ep.m$Experiment = factor(dat.ep.m$Experiment, levels=c('Exp.', 'fba','E_only', 'Both', 'P_only'))

dp<-dat.ep.m
xl<-""
yl<-"R/F"
set.seed(10)
ggplot(dp, aes(x=Experiment, y=ratioRF))+
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, aes(shape=Experiment, colour=Experiment), size=2, alpha=0.4)+
  theme_classic()+
  scale_shape_manual(values = c(1, 1, 1,1,1)) +
  scale_colour_manual(values = c(colp,colp,colp,colp,colp))+
  theme(
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=24),
    axis.title.y =element_text(size=28),
    axis.title.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"))+
  scale_x_discrete( 
    labels=c("communities 
(16S)", "model 0", "model 1","model 2","model 3"))+
  theme(legend.position = "none")+
  scale_y_log10(breaks = c(0.01, 0.3,1,10),limits = c(0.001, 10))+
  labs(y=yl)

ggsave(path=path.p, paste0("figS16.pdf"), width= 10, height=5, onefile = TRUE)

dat.figS16_fba_c<-dp

## ..Write data table as csv file
fname<-paste("figS16.csv")
write.csv(dat.figS16_fba_c,file.path(path.ds, fname),row.names=FALSE)


