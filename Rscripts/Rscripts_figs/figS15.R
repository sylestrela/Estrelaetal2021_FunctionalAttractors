" R script for Figure S15 in Estrela et al (2021) 
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

#.. Read data 
setwd(path.d)

#. read data for ES strains
d.es<-read.csv('Estrela_2021_empir_params_Pairs_Acetate.csv')
d.es$strains<-'Fig. 1'
d.es$rep<-1

#. read data for KKAP strains
d.kkap<-read.csv('Estrela_2021_empir_params_KKAP_acetate.csv')
d.kkap$strains<-'KpKmAP'

dc<-rbind(d.es,d.kkap)

#. colour strains
colFp<-"deepskyblue3"
colPseu2<-"darkorchid4"
colKp<-blues9[7]
colKm<-blues9[5]
colAlc<-'orange2'
colAlc<-"#E69F00"

dcm<-melt(dc, id.vars=c('SangerID_f', 'SangerID_r','Family_f','Family_r','strains'))

dcm1<-subset(dcm, variable %in% c('w_f_glu', 'D_ace_glu_corr'))
dcm1<-unique(select(dcm1, -SangerID_r,-Family_r))
names(dcm1)[names(dcm1)=='SangerID_f'] <-'SangerID'
names(dcm1)[names(dcm1)=='Family_f'] <-'Family'

# unselect K- 
dcm2<-subset(dcm, variable %in% c('w_r_ace') & SangerID_f!=235)
dcm2<-unique(select(dcm2, -SangerID_f,-Family_f))
names(dcm2)[names(dcm2)=='SangerID_r'] <-'SangerID'
names(dcm2)[names(dcm2)=='Family_r'] <-'Family'

dcp<-rbind(dcm1,dcm2)

dkap<-subset(dcp, strains=='KpKmAP' & variable=='D_ace_glu_corr')
dkap$strains[dkap$SangerID==225]<-'Kp (Fig. 2)'
dkap$strains[dkap$SangerID==235]<-'Km (Fig. 2)'

set.seed(10)
p1<-ggplot(subset(dcp, strains=='Fig. 1' & variable=='D_ace_glu_corr'), aes(x='', y=value, colour=strains, group=factor(strains):factor(variable)))+
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size=2.5, alpha=0.7)+
  geom_point(data=dkap,size=2.5, alpha=0.9)+
  theme_classic()+
  scale_colour_manual(values = c('gray75', colKm,colKp))+
  ylim(0,1)+
  theme(
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    axis.title.y =element_text(size=16),
    axis.title.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "gray10"))+
  scale_x_discrete(labels=c("D[ace, glu]"))+
  labs(y='')
p1
set.seed(10)
dkap<-subset(dcp, strains=='KpKmAP' & variable=='w_f_glu')
dkap$strains[dkap$SangerID==225]<-'Kp (Fig. 2)'
dkap$strains[dkap$SangerID==235]<-'Km (Fig. 2)'
p2<-ggplot(subset(dcp, strains=='Fig. 1' & variable=='w_f_glu'), aes(x='', y=value, colour=strains, group=factor(strains):factor(variable)))+
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size=2.5, alpha=0.7)+
  geom_point(data=dkap,size=2.5, alpha=0.9)+
  theme_classic()+
  scale_colour_manual(values = c('gray75', colKm,colKp))+
  ylim(0,0.026)+
  theme(
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    axis.title.y =element_text(size=16),
    axis.title.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "gray10"))+
  scale_x_discrete(labels=c("W[glu, F]"))+
  labs(y='')
p2

dkap<-subset(dcp, strains=='KpKmAP' & variable=='w_r_ace')
dkap$strains[dkap$SangerID==241]<-'A (Fig. 2)'
dkap$strains[dkap$SangerID==226]<-'P (Fig. 2)'

#.. calculate mean of 4 replicates for each strain
dkap.m<- dkap %>%
  dplyr::group_by(SangerID, Family, strains) %>%
  dplyr::summarise(ace.m=round(mean(value),4),
                   sd=sd(value),
                   n=n(),
                   se = sd/sqrt(n)
  )
dkap.m<-as.data.frame(dkap.m)

set.seed(10)
dcps<-subset(dcp, strains=='Fig. 1' & variable=='w_r_ace')
p3<-ggplot()+
  geom_boxplot(data=dcps, aes(x='', y=value, colour=strains, group=factor(strains):factor(variable)), outlier.shape = NA) + 
  geom_jitter(data=dcps, aes(x='', y=value, colour=strains, group=factor(strains):factor(variable)), width = 0.2, size=2.5, alpha=0.7)+
  geom_point(data=dkap.m, aes(x='',y=ace.m, colour=strains), size=2.5, alpha=0.9)+
  geom_errorbar(data=dkap.m, aes(x='',ymin=ace.m-sd, ymax=ace.m+sd, width=0.04,  colour=strains)) +
  theme_classic()+
  scale_colour_manual(values = c(colAlc, 'gray75', colPseu2))+
  ylim(0,0.026)+
  theme(
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    axis.title.y =element_text(size=16),
    axis.title.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "gray10"))+
  scale_x_discrete(labels=c("W[ace, R]"))+
  labs(y='')
p3

kpa<-subset(d.kkap, strains=='KpKmAP' & SangerID_f==225)
kpa$strains[kpa$SangerID_r==226]<-"Kp+P"
kpa$strains[kpa$SangerID_r==241]<-"Kp+A"

#.. calculate mean of 4 replicates for each strain
kpa.m<- kpa %>%
  dplyr::group_by(SangerID_f, Family_f, SangerID_r, Family_r, strains) %>%
  dplyr::summarise(rf.m=round(mean(ratioRF),4),
                   sd=sd(ratioRF),
                   n=n(),
                   se = sd/sqrt(n)
  )
kpa.m<-as.data.frame(kpa.m)

gg1<-ggarrange(p1,p2,p3, ncol=3,nrow=1)
gg1
ggsave(path=path.p, paste0("figS15A.pdf"), width= 10, height=3.5, onefile = TRUE)

#  ---------------------------------------------------------------------------
#.. Panel B where R/F ratio for KA and KP is overlaid on top of RF ratios for each of the states
#. read data from figure 4B
df.all<-read.csv(paste0(path.ds, '/fig5C.csv'))
kpa.c<-subset(df.all, Group!='K')

# .. Plot
xl<-"community dominant R strain"
yl<-"R/F"

#.. re-name groups
kpa.m$Group[kpa.m$strains=='Kp+A']<-'KA'
kpa.m$Group[kpa.m$strains=='Kp+P']<-'KP'
kpa.m$type[kpa.m$strains=='Kp+A']<-'Kp+A (empir. model)'
kpa.m$type[kpa.m$strains=='Kp+P']<-'Kp+P (empir. model)'
kpa.c$type<-'communities (16S)'

set.seed(1)
p<- ggplot()+
  geom_boxplot(data=kpa.c, aes(x=factor(Group), ratioRF),outlier.shape = NA) + 
  geom_jitter(data=kpa.c, aes(x=factor(Group), ratioRF, colour=type), width = 0.2, size=2.5, alpha=0.7)+
  theme_classic()+
  geom_point(data=kpa.m, aes(x=factor(Group),y=rf.m, colour=strains), size=2.5, alpha=1)+
  geom_errorbar(data=kpa.m, aes(x=factor(Group),ymin=rf.m-sd, ymax=rf.m+sd, width=0.02,  colour=strains)) +
  scale_colour_manual(values=c('gray75', colAlc, colPseu2))+
  scale_x_discrete(breaks=c("KA","KP"),
                   labels=c("A", "P"))+
  theme(
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    axis.title.y =element_text(size=20),
    axis.title.x = element_text(size = 12), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16))+
  labs(x=xl, y=yl)+
  scale_y_log10(breaks = c(0.0001, 0.01, 0.1, 1,10),limits = c(0.0001, 4))
p

ggsave(path=path.p, paste0("figS15B.pdf"), width= 6, height=3.5, onefile = TRUE)


