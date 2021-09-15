" R script for Figure 1D and S10 (+ stats) in Estrela et al (2021) 
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

# ..Read data table
d.fba<-read.csv('Estrela_2021_fba_params_Pairs_Acetate.csv')
dat.fba<-select(d.fba, -E,-P,-E_Genus,-P_Species)
dat.fba$ratioRF<-dat.fba$D_ace_glu*dat.fba$w_ac_R/dat.fba$w_glc_F
dat.fba$Experiment<-'fba'
dat.fba<-select(dat.fba, -D_ace_glu, -w_ac_R,-w_glc_F,-R_F)

# ..Read data table for RF ratio of ES communities
dat.exp <- read.csv(paste0(path.dp,"/Goldford_2018_glu_ratioRF.csv"))
dat.exp$Experiment<-"Exp."
dat.exp<-select(dat.exp,-InocRep,-F,-R)

# ..Read data table for PE/RF ratio of empirically calibrated model
dat.emp.calc<-read.csv('Estrela_2021_empir_params_Pairs_Acetate.csv')
dat.emp.calc<-select(dat.emp.calc, ratioRF)
dat.emp.calc$Experiment<-'empir_calc'
  
# .merge 3 dataframes
dat.rf.m<-rbind(dat.exp,dat.fba, dat.emp.calc)

# .. Colour scheme 
# Family colors 
colEntbf<- "deepskyblue3" 
colPseuf<-  "darkorchid2" 

# .. Calculate median and IQR for FBA and for experiments
ds<-dat.rf.m
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
    Q75= round(quantile(ratioRF, 0.75, na.rm = TRUE),nr),
    count=n())
stat_summ<-as.data.frame(stat_summ)
stat_summ
stat.fba<-subset(stat_summ, Experiment=="fba")
stat.16s<-subset(stat_summ, Experiment=="Exp.")
stat.empc<-subset(stat_summ, Experiment=="empir_calc")

colp<-"darkorchid4"

#. re-order factors
dat.rf.m$Experiment = factor(dat.rf.m$Experiment, levels=c('Exp.',"empir_calc",'fba'))

dp<-dat.rf.m
xl<-""
yl<-"R/F"
set.seed(10)

p<-ggplot(dp, aes(x=Experiment, y=ratioRF))+
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, aes(shape=Experiment, colour=Experiment), size=2, alpha=0.4)+
  theme_classic()+
  scale_shape_manual(values = c(1, 2, 4)) +
  scale_colour_manual(values = c(colp,colp,colp))+
  theme(
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=24),
    axis.title.y =element_text(size=28),
    axis.title.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none")+
  scale_x_discrete( 
    labels=c("communities (16S)", "empirical model","FBA simulations"))+
  scale_y_log10(breaks = c(0.01, 0.1,1,10),limits = c(0.001, 10))+
  labs(y=yl)
p
ggsave(path=path.p, paste0("fig1D.pdf"), width= 8, height=5, onefile = TRUE)

dat.fig1D<-dp
## ..Write data table as csv file
fname<-paste("fig1D.csv")
write.csv(dat.fig1D,file.path(path.ds, fname),row.names=FALSE)

#  ---------------------------------------------------------------------------
#.. Figure S10
#. histogram for each pair of Enterobacteriacea + all Pseudomonas and for each pair of Pseudomonas + all Enterobacteriacea as facet grid
#. global P/E median indicated by purple dashed line + median for each pair in gray

ds<-d.fba
ds$ratioRF<-ds$D_ace_glu*ds$w_ac_R/ds$w_glc_F
colp<-'darkorchid3'

#. histogram for Ent as facet

dtE.med <- ds %>% 
  dplyr::group_by(E_Genus) %>% 
  dplyr::summarize(med = median(ratioRF))
dtE.med<-as.data.frame(dtE.med)
dtE.med

d.stat<-stat.fba

pE<-ggplot(ds,aes(ratioRF)) + 
  geom_histogram(binwidth=0.05, fill='gray', alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept =d.stat$median,col=colp, linetype='dashed', size=1)  +
  facet_wrap(~E_Genus, scales="free_y", ncol=4)+
  geom_vline(aes(xintercept = dtE.med$med), dtE.med, colour = 'black', linetype='dashed', size=1, alpha=0.6)+
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.y =element_text(size=18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14))+
  labs(x='R/F',y = 'p(R/F)')
pE  

#. histogram for Pseud as facet
#. re-order facets to have P. putida in first position
#. correct names
levels(ds$P_Species) <- c('P. aeruginosa','P. fluorescens','P. putida','P. stutzeri')
ds$P_Species_ord = factor(ds$P_Species, levels=c('P. putida','P. aeruginosa','P. fluorescens','P. stutzeri'))

dtP.med <- ds %>% 
  dplyr::group_by(P_Species) %>% 
  dplyr::summarize(med = median(ratioRF))
dtP.med<-as.data.frame(dtP.med)  

pP<-ggplot(ds,aes(ratioRF)) + 
  geom_histogram(binwidth=0.05, fill='gray', alpha=0.8) +
  theme_classic() +
  geom_vline(xintercept =d.stat$median,col=colp, linetype='dashed', size=1)  +
  facet_wrap(~P_Species, scales="free_y", ncol=4)+
  geom_vline(aes(xintercept = dtP.med$med), dtP.med, colour = 'black', linetype='dashed', size=1,alpha=0.7)+
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.y =element_text(size=18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14))+
  labs(x='R/F',y = 'p(R/F)')
pP

pPE.2<-ggarrange(pE,pP,
                 ncol = 1, nrow = 2)
pPE.2 

ggsave(path=path.p, paste0("figS10.pdf"), width= 12, height=6, onefile = TRUE)

dat.figS10<-select(ds, E_Genus, P_Species, ratioRF)

## ..Write data table as csv file
fname<-paste("figS10.csv")
write.csv(dat.figS10,file.path(path.ds, fname),row.names=FALSE)


