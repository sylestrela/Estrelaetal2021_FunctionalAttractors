" R script for Figure S8 in Estrela et al (2021) 
Functional attractors in microbial community assembly."

#  ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
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
# ..Read data table for Migration experiment communities
df.all<-read.csv("Estrela_2021_16S_data_table_RelAbund_ALL.csv")

#.. SELECT data
trn<-"No_migration"
tp<-"18"
inoc<-4
df.all.s<- subset(df.all, Treatment== trn & Transfer==tp & Inoc_size==4)
head(df.all.s)

df.sub<-select(df.all.s, -Carbon, -Transfer, -Treatment, -OTU,-Kingdom,-Phylum,-Class,-Order,-Inoc_size)
head(df.sub)

##.. calculate ratio of R/F in Alc/Ent communities
#.. create new column with F or R assignment
df.sub$f_r_type <- with(df.sub, 
                        ifelse(Family == "Enterobacteriaceae" | Family == "Enterococcaceae" | Family=="Aeromonadaceae" | Family=="Lachnospiraceae",
                               "F", ifelse(Family == "Alcaligenaceae" | Family == "Pseudomonadaceae" | Family=="Moraxellaceae" | Family=='Xanthomonadaceae' | Family=="Comamonadaceae",
                                           "R", 'UN')))

head(df.sub)

# .. Select communities with A dominated state
nAd<-subset(df.sub, Family== 'Alcaligenaceae' & Abundance>0.01)
nAcd<-subset(df.sub, Rep_num %in% nAd$Rep_num)
length(unique(nAd$Rep_num))

# .. SELECT communities to plot (here we want Alc and Ent)
dat.fr<-nAcd

#.. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fr<-data.table(dat.fr)[, list(f_r_sum=sum(Abundance)), 
                            by=list(f_r_type, Rep_num)]
dat.RF<-as.data.frame(spread(dsum_fr,f_r_type, f_r_sum))
dat.RF[is.na(dat.RF)] <- 0

#. ratio of Respirer/Fermenter
dat.RF$ratioRF<-(dat.RF$R/dat.RF$F)
median(dat.RF$ratioRF)
mean(dat.RF$ratioRF)

IQR_summary<- dat.RF %>%
  dplyr::summarise(
    Min = round(min(ratioRF),3),
    Max = round(max(ratioRF),3),
    mean = round(mean(ratioRF),3),
    std = round(sd(ratioRF),3),
    Median = median(ratioRF),
    IQRange = IQR(ratioRF),
    Q1=quantile(ratioRF, 0.25),
    Q3=quantile(ratioRF, 0.75))
IQR_summary

dat.Ent.Alc<-IQR_summary
dat.Ent.Alc$type<-"Ent.Alc"

#------------------------------------------------------------------------------
# ..Read data table for C11R2 (from Goldford et al 2018) with  Aero+Ent -Pseu+morax
dat.ES<-read.csv('Goldford_2018_16S_data_table_glucose_RelAbund_ALL.csv')

# ..SELECT minimum abundance and SELECT glucose data only 
min.ab<-0.0
df.sub<-subset(dat.ES, Relative_Abundance >min.ab)
colnames(df.sub) 

df.glu<-subset(df.sub, Carbon_Source=="D-Glucose")
df.glu<-select(df.glu, -Order,-Experiment,-ESV,-Inoculum,-Replicate,-ESV_ID,-Carbon_Source,-Transfer)

# .. SELECT communities to plot (here we want C11R2)
dat.fr<-subset(df.glu, InocRep=="C11R2")
dat.fr

#.. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fr<-data.table(dat.fr)[, list(f_r_sum=sum(Relative_Abundance)), 
                            by=list(f_r_type, InocRep)]
dat.RF<-as.data.frame(spread(dsum_fr,f_r_type, f_r_sum))
dat.RF[is.na(dat.RF)] <- 0

#. ratio of Respirer/Fermenter
dat.RF$ratioRF<-(dat.RF$R/dat.RF$F)
median(dat.RF$ratioRF)
mean(dat.RF$ratioRF)

IQR_summary<- dat.RF %>%
  dplyr::summarise(
    Min = round(min(ratioRF),3),
    Max = round(max(ratioRF),3),
    mean = round(mean(ratioRF),3),
    std = sd(ratioRF),
    Median = median(ratioRF),
    IQRange = IQR(ratioRF),
    Q1=quantile(ratioRF, 0.25),
    Q3=quantile(ratioRF, 0.75))
IQR_summary

dat.Aero.Morax.c11r2<-IQR_summary
dat.Aero.Morax.c11r2$type<-'Aero.Morax.c11r2'

#------------------------------------------------------------------------------
# .. SELECT communities to plot (here we want all from Goldford et al 2018)
dat.fr<-df.glu
dat.fr
#.. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fr<-data.table(dat.fr)[, list(f_r_sum=sum(Relative_Abundance)), 
                            by=list(f_r_type, InocRep)]
dat.RF<-as.data.frame(spread(dsum_fr,f_r_type, f_r_sum))
dat.RF[is.na(dat.RF)] <- 0

#. ratio of Respirer/Fermenter
dat.RF$ratioRF<-(dat.RF$R/dat.RF$F)
median(dat.RF$ratioRF)
mean(dat.RF$ratioRF)

IQR_summary<- dat.RF %>%
  dplyr::summarise(
    Min = round(min(ratioRF),3),
    Max = round(max(ratioRF),3),
    mean = round(mean(ratioRF),3),
    std = sd(ratioRF),
    Median = median(ratioRF),
    IQRange = IQR(ratioRF),
    Q1=quantile(ratioRF, 0.25),
    Q3=quantile(ratioRF, 0.75))
IQR_summary

dat.ES.glu<-IQR_summary
dat.ES.glu$type<-'ES_all'

#------------------------------------------------------------------------------
#.. Read data table for Aero-Morax community (ltcee)
dat.ltcee<-read.csv('Estrela_2021_AeroMorax_12com_RFratio.csv')
dat.ltcee$ratioRF<-dat.ltcee$rf_ratio

dat.RF<-dat.ltcee
IQR_summary<- dat.RF %>%
  dplyr::summarise(
    Min = round(min(ratioRF),3),
    Max = round(max(ratioRF),3),
    mean = round(mean(ratioRF),3),
    std = sd(ratioRF),
    Median = median(ratioRF),
    IQRange = IQR(ratioRF),
    Q1=quantile(ratioRF, 0.25),
    Q3=quantile(ratioRF, 0.75))
IQR_summary

dat.Aero.Morax.ltcee<-IQR_summary
dat.Aero.Morax.ltcee$type<-'Aero.Morax.ltcee'

#------------------------------------------------------------------------------
# .merge ALL dataframes
dat.combo<-rbind(dat.Ent.Alc, dat.Aero.Morax.c11r2, dat.Aero.Morax.ltcee)
dat.combo$type <- factor(dat.combo$type,levels = c("Ent.Alc", "Aero.Morax.c11r2", "Aero.Morax.ltcee"))

# .. Colour scheme 
# Family colors 
colEntbf<- "deepskyblue3"
colPseuf<-  "darkorchid2" 

colp<-"darkorchid4"

dp<-dat.combo
xl<-""
yl<-"R/F"

dp$type<-revalue(dp$type, c("Ent.Alc"="Alc./E", "Aero.Morax.c11r2"="(M+P)/(Aero.+E)","Aero.Morax.ltcee"="M/Aero."))

p<-ggplot(dp, aes(x=type, y=Median,colour=type)) +
  geom_point(colour=colp,size=3)+
  geom_errorbar(data=dp, aes(ymin=Q1 , 
                ymax=Q3), width=.05, colour=colp,size=1)+
  geom_hline(yintercept=dat.ES.glu$Median, colour=colp)+
  geom_hline(yintercept=dat.ES.glu$Q1, colour=colp, linetype='dashed')+
  geom_hline(yintercept=dat.ES.glu$Q3, colour=colp, linetype='dashed')+
  labs(y=yl, x="")+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.y =element_text(size=18),
    axis.title.x = element_text(size = 8),
    axis.line= element_line(colour="gray"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  scale_y_log10(breaks = c(0.01, 0.3,1,10), limits = c(0.001, 10))
p  

p<-p+
  annotate("rect", ymin = dat.ES.glu$Q1, ymax = dat.ES.glu$Q3, xmin = 0, xmax = Inf,
           alpha = .1, fill=colp)
p

ggsave(path=path.p, paste0("figS8.pdf"), width=6, height=4, onefile = TRUE)

dat.m<-rbind(dp,dat.ES.glu)
dat.figS8<-select(dat.m,type, Median, Q1,Q3)
## ..Write data table as csv file
fname<-paste("figS8.csv")
write.csv(dat.figS8,file.path(path.ds, fname),row.names=FALSE)


