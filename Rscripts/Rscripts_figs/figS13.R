" R script for Figure S13 (+stats) in Estrela et al (2021) 
Functional attractors in microbial community assembly."

## ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(reshape2)
library(tidyr)
## ---------------------------------------------------------------------------
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
dat.ctrl<-read.csv("Estrela_2021_ctrls_16S.csv")
setwd(path.dp)
d.esvs<-read.csv('top6_ESVs_fig2B.csv')

# .. SELECT treatment
df.all$experiment<-paste0(df.all$Treatment,'_', df.all$Inoc_size, '_', df.all$Transfer)
df.all$experiment_s<-paste0(df.all$Treatment,'_', df.all$Inoc_size)
tt<-c("No_migration_4")
df.all.s<-subset(df.all, experiment_s %in% tt)
df.all.s<-select(df.all.s,-Kingdom,-Phylum,-Class,-Order)
head(df.all.s)

##.. SELECT A and P at T18 
esvAlc<-as.character(subset(d.esvs,Genus=='Alcaligenes')$OTU)
esvP1<-as.character(subset(d.esvs,Genus=='Pseudomonas')$OTU)
dat.ap<-subset(df.all, OTU %in% c(esvAlc, esvP1))
dat.ap$esv_id[dat.ap$OTU==esvAlc]<-'A'
dat.ap$esv_id[dat.ap$OTU==esvP1]<-'P'
dat.ap<-select(dat.ap, -OTU,-Carbon,-Family,-Treatment,-Inoc_size)

#.. select A and P present in ctrl
dt.p<-subset(dat.ctrl, Genus=="Pseudomonas" & OTU==esvP1)
dt.a<-subset(dat.ctrl, Genus=="Alcaligenes" & OTU==esvAlc)
dt.esvAP<-rbind(dt.p, dt.a)

dat.ctrl.AP<-dt.esvAP
ESV_summ_ctrl<- dat.ctrl.AP %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(
    mean = mean(Abundance),
    std = sd(Abundance),
    median = median(Abundance))
esv.ctrl<-as.data.frame(ESV_summ_ctrl)

##.. calculate 95% confidence interval (error = qnorm(0.975)*SD/sqrt(Sample size))
mean.ctrl<-esv.ctrl$mean
std.ctrl<-esv.ctrl$std
n.ctrl<-length(unique(paste0(dat.ctrl.AP$Sample,dat.ctrl.AP$Seq_run)))

merror<-qnorm(0.975)*std.ctrl/sqrt(n.ctrl)
cf95_upper<-mean.ctrl+merror
cf95_lower<-mean.ctrl-merror
cf95_upper
mean.ctrl
log10(cf95_upper)
log10(cf95_lower)
log10(mean.ctrl)

esv.crtl.cf95<-cbind(esv.ctrl,as.data.frame(cf95_upper),as.data.frame(cf95_lower))

a.cf95<-esv.crtl.cf95$cf95_upper[1]
a.m<-esv.crtl.cf95$mean[1]

p.cf95<-esv.crtl.cf95$cf95_upper[2]
p.m<-esv.crtl.cf95$mean[2]

#.. Select dominant A and P for the no migration, low inocula, T18 data
dat.ap<-subset(df.all.s, experiment=='No_migration_4_18' & OTU %in% c(esvAlc, esvP1))
dat.ap$Abund_log<-log10(dat.ap$Abundance)

#.. esv colours
colAlc<-"#E69F00"
colPseu2<-"darkorchid3"

#.. select A communities with low A (i.e. lower than threshold invasion)
df.al<-subset(dat.ap, Genus=='Alcaligenes' & Abund_log < (-1.18))
mean.al<-mean(df.al$Abund_log)
pA<-ggplot(df.al, aes(x=factor(Rep_num), y=log10(Abundance))) +
  geom_point(colour=colAlc)+
  geom_hline(yintercept=mean.al, colour=colAlc,linetype='dashed')+
  geom_hline(yintercept=log10(a.m), colour='gray', linetype='dashed')+
  geom_hline(yintercept=log10(a.cf95), colour='gray')+
  labs(y='log10(A)', x="replicate")+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.y =element_text(size=16),
    axis.title.x = element_text(size = 16),
    axis.line= element_line(colour="gray"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) 
pA  

pA<-pA+
  annotate("rect", ymin = -Inf, ymax = log10(a.cf95), xmin = 0, xmax = Inf,
           alpha = .3, fill='gray')
pA

#.. select P communities with low P (i.e. lower than threshold invasion)
df.pl<-subset(dat.ap, Genus=='Pseudomonas' & Abund_log < (-1.97))
mean.pl<-mean(df.pl$Abund_log)
pP<-ggplot(df.pl, aes(x=factor(Rep_num), y=log10(Abundance))) +
  geom_point(colour=colPseu2)+
  geom_hline(yintercept=-3.61, colour=colPseu2,linetype='dashed')+
  geom_hline(yintercept=log10(p.m), colour='gray',linetype='dashed')+
  geom_hline(yintercept=log10(p.cf95),colour='gray')+
  theme(axis.text.x=element_blank())+
  labs(y='log10(P)', x="replicate")+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.y =element_text(size=16),
    axis.title.x = element_text(size = 16),
    axis.line= element_line(colour="gray"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) 
pP

pP<-pP+
  annotate("rect", ymin = -Inf, ymax = log10(p.cf95), xmin = 0, xmax = Inf,
           alpha = .3, fill='gray')
pP

ggarrange(pA, pP, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
ggsave(path=path.p, paste0("figS13.pdf"), width= 8, height=3.5, onefile = TRUE)

## ..determine A communities that don't pass the 95% conf interval test  
df.al.fail.cftest<-subset(dat.ap, Genus=='Alcaligenes' & Abund_log < log10(a.cf95))
df.al.fail.cftest
## ..determine P communities that don't pass the 95% conf interval test  
df.pl.fail.cftest<-subset(dat.ap, Genus=='Pseudomonas' & Abund_log < log10(p.cf95))
df.pl.fail.cftest

## ..select communities with both A and P present
nAd<-subset(dat.ap, Genus=='Alcaligenes' & Abund_log > log10(a.cf95))
nPd<-subset(dat.ap, Genus=='Pseudomonas' & Abund_log > log10(p.cf95))

nAdPd<-rbind(nAd,nPd)
nAdPd<-select(nAdPd,Rep_num,Genus,Abund_log)
library(tidyr)
nAdPd.cf95.pass<-spread(nAdPd, Genus,Abund_log)
nAdPd.cf95.pass<-na.omit(nAdPd.cf95.pass)
length(nAdPd.cf95.pass$Rep_num)

## ..exclude communities where both are rare
nAd001<-subset(dat.ap, Genus=='Alcaligenes' & Abundance>0.01)
nPd001<-subset(dat.ap, Genus=='Pseudomonas' & Abundance>0.01)
nAdPd001<-rbind(nAd001,nPd001)
nFd001<-subset(dat.ap, !(Rep_num %in% nAdPd001$Rep_num))
nAdPd.cf95.pass.corr<-subset(nAdPd.cf95.pass, !(Rep_num %in% nFd001$Rep_num))
length(nAdPd.cf95.pass.corr$Rep_num)
length(nAdPd001$Rep_num)
## ..the % of communities with either A rare and P abundant or A abundant and P rare is:
length(nAdPd.cf95.pass.corr$Rep_num)/length(nAdPd001$Rep_num)

## ..% of communities with AhighPlow and PhighAlow 
a_dom<-subset(nAdPd.cf95.pass.corr, Alcaligenes > (-2))
p_dom<-subset(nAdPd.cf95.pass.corr, Pseudomonas > (-2))
length(a_dom$Rep_num)/length(nAdPd001$Rep_num)
length(p_dom$Rep_num)/length(nAdPd001$Rep_num)

## ..% of A dominated communities with low P and % of P dominated communities with low A 
nAd1<-subset(dat.ap, Genus=='Alcaligenes' & Abundance > 0.01)
nAdPl<-subset(nAd1, (Rep_num %in% nAdPd.cf95.pass.corr$Rep_num))
length(nAdPl$Rep_num)/length(nAd1$Rep_num)

nPd1<-subset(dat.ap, Genus=='Pseudomonas' & Abundance > 0.01)
nPdAl<-subset(nPd1, (Rep_num %in% nAdPd.cf95.pass.corr$Rep_num))
length(nPdAl$Rep_num)/length(nPd1$Rep_num)

## ..calculate relative abundance and STD, IQR for high A, high P, low A, low P
nAd<-subset(dat.ap, Genus=='Alcaligenes' & Abundance>0.01 & Transfer==18)
mean(nAd$Abundance)
sd(nAd$Abundance)
nPd<-subset(dat.ap, Genus=='Pseudomonas' & Abundance>0.01 & Transfer==18)
mean(nPd$Abundance)
sd(nPd$Abundance)
nAr<-subset(dat.ap, Genus=='Alcaligenes' & Abundance<0.01 & Transfer==18)
mean(nAr$Abundance)
sd(nAr$Abundance)
nPr<-subset(dat.ap, Genus=='Pseudomonas' & Abundance<0.01 & Transfer==18)
mean(nPr$Abundance)
sd(nPr$Abundance)
