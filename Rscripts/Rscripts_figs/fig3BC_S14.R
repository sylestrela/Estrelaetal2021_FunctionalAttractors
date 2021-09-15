" R script for Figure 3B, 3C and S14 in Estrela et al (2021) 
Functional attractors in microbial community assembly."

#------------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
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

dat.kap<-read.csv("Estrela_2021_cfu_data_invKAP_all.csv")
dp<-dat.kap

#. define order of culture for plots and rename
dp$dila_ord=factor(dp$dila,levels=c("-","-4", "-3", "-2","-1" ,"0"))
dp$dilp_ord=factor(dp$dilp,levels=c("-","-4", "-3", "-2","-1" ,"0"))

##-------------------------------------------------
## Figure 3 - outcome of invasion as density plot
##-------------------------------------------------
#. Colour scheme for plots- ESVs color
colKp<-'steelblue4'
colKm<-'#6993BC'
colPseu<-'#A17BB9'
colAlc<-'#FFB94C'

##. get outcome of coexistence based on presence/absence of A and/or P
dp$outcome <- with(dp, ifelse(Alcaligenes >0 & Pseudomonas > 0, "K+P+A", 
                              ifelse(Alcaligenes >0 & Pseudomonas == 0, "K+A",
                                     ifelse(Alcaligenes ==0 & Pseudomonas > 0, "K+P",
                                            ifelse(Alcaligenes ==0 & Pseudomonas == 0, "K",
                                     0)))))

#. plot outcome
col.p <- c("K" = colKp, "K+A" = colAlc, "K+P" = colPseu, "K+P+A" = "gray75")
ggplot(dp, aes(x=dilp_ord, y=dila_ord))+ 
  geom_tile(aes(fill=outcome))+
  facet_grid(repb~Transfer)+
  labs(x=expression(phantom(x)*P[0]*phantom(x)*(log[10])),
       y=expression(phantom(x)*A[0]*phantom(x)*(log[10])))+
  theme_minimal()+
  scale_fill_manual(values = col.p,
                    name= "Outcome")
ggsave(path=path.p, paste0("fig3B_S14A.pdf"), width= 7, height=4.2, onefile = TRUE)

dat.fig3<-select(dp,Transfer,repb,dila_ord,dilp_ord,outcome)
## ..Write data table as csv file
fname<-paste("fig3B_S14A.csv")
write.csv(dat.fig3,file.path(path.ds, fname),row.names=FALSE)

##-------------------------------------------------
## Figure S14. timeseries of A vs P
##-------------------------------------------------

##.. Create dataframe with initial frequencies
a.init<-c(0.1,0.01,0.001,0.0001,0.00001,0)
p.init<-c(0.1,0.01,0.001,0.0001,0.00001,0)
k.init<-c(0.1)

da<-as.data.frame(a.init)
dpp<-as.data.frame(p.init)
dk<-as.data.frame(k.init)

dap<-merge(da,dpp)
dkap<-merge(dk, dap)

dkap$tot<-dkap$k.init+dkap$a.init+dkap$p.init
dkap$freqK<-dkap$k.init/dkap$tot
dkap$freqA<-dkap$a.init/dkap$tot
dkap$freqP<-dkap$p.init/dkap$tot

dkap$Transfer<-0

dps$dila_ord=factor(dps$dila_ord,levels=c("-","-4", "-3", "-2","-1" ,"0"))

##. match to dila_ord and dilp_ord
dkap$dila_ord <- with(dkap, ifelse(a.init== 0, "-", 
                              ifelse(a.init== 1e-05, "-4",
                                     ifelse(a.init== 1e-04, "-3",
                                            ifelse(a.init== 1e-03, "-2",
                                                   ifelse(a.init== 1e-02, "-1",
                                                   0))))))

dkap$dilp_ord <- with(dkap, ifelse(p.init== 0, "-", 
                                   ifelse(p.init== 1e-05, "-4",
                                          ifelse(p.init== 1e-04, "-3",
                                                 ifelse(p.init== 1e-03, "-2",
                                                        ifelse(p.init== 1e-02, "-1",
                                                               0))))))

dkapm.t0<-melt(select(dkap, -k.init,-a.init, -p.init,-tot), id.vars=c('Transfer','dila_ord', 'dilp_ord'))

colnames(dkapm.t0)[colnames(dkapm.t0)=="value"] <- "Rel_abund"
colnames(dkapm.t0)[colnames(dkapm.t0)=="variable"] <- "Species"

dkapm.t0$Species <-revalue(dkapm.t0$Species , c( "freqK"="K",  "freqA"="A", 'freqP'= "P"))

#. define order of culture for plots and rename
dkapm.t0$dila_ord=factor(dkapm.t0$dila_ord,levels=c("-","-4", "-3", "-2","-1" ,"0"))
dkapm.t0$dilp_ord=factor(dkapm.t0$dilp_ord,levels=c("-","-4", "-3", "-2","-1" ,"0"))

#.add 2 replicates (to be merged with data below)
repb<-c(1,2)
rep<-as.data.frame(repb)
dkapm.t0<-merge(dkapm.t0, rep )

## ..Calculate relative abundance 
dat<-dp
dat$K<-dat$Klebp /(dat$Klebp + dat$Pseudomonas + dat$Alcaligenes)
dat$A<-dat$Alcaligenes /(dat$Klebp + dat$Pseudomonas + dat$Alcaligenes)
dat$P<-dat$Pseudomonas /(dat$Klebp + dat$Pseudomonas + dat$Alcaligenes)

# .melt data
datm<-melt(dat, measure.vars=c('K', 'P', "A"))
colnames(datm)[colnames(datm)=="value"] <- "Rel_abund"
colnames(datm)[colnames(datm)=="variable"] <- "Species"

#.. Barplots of relative abundance of K, P and A
dfp<-datm
dps<-select(dfp,Species,Rel_abund,Transfer,dila_ord,dilp_ord,repb)

#. define order of culture for plots and rename
dps$dila_ord=factor(dps$dila_ord,levels=c("-","-4", "-3", "-2","-1" ,"0"))
dps$dilp_ord=factor(dps$dilp_ord,levels=c("-","-4", "-3", "-2","-1" ,"0"))

dps$Transfer <- factor(dps$Transfer, levels = c("3","8","12"))
colp<-c(colKp, colAlc,colPseu)

#.. combine T0 data
dkapm.t0$Transfer <- factor(dkapm.t0$Transfer, levels = c("0"))
dps0<-rbind(dkapm.t0, dps)

#. SELECT replicate to plot
repn<-"1"
dps1<-subset(dps0, repb==repn)
p1<-ggplot(dps1, aes(x=Transfer, y=Rel_abund, fill=Species))+
  geom_bar(stat="identity", position="stack") + 
  facet_grid(dila_ord ~dilp_ord)+ #labeller = label_both)
  scale_fill_manual(values=colp)+ 
  scale_y_continuous(breaks=c(0, 1))+
  labs(y= "Relative abundance", x="Transfer")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=11),
    axis.title.y =element_text(size=16),
    axis.title.x = element_text(size=16), 
    strip.text.x =element_text(size=12),
    strip.text.y =element_text(size=12),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())
p1

#. SELECT replicate to plot
repn<-"2"
dps1<-subset(dps0, repb==repn)
p2<-ggplot(dps1, aes(x=Transfer, y=Rel_abund, fill=Species))+
  geom_bar(stat="identity", position="stack") + 
  facet_grid(dila_ord ~dilp_ord)+ 
  scale_fill_manual(values=colp)+ 
  scale_y_continuous(breaks=c(0, 1))+
  labs(y= "Relative abundance", x="Transfer")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size=11),
    axis.text.y = element_text(size=11),
    axis.title.y =element_text(size=16),
    axis.title.x = element_text(size=16), 
    strip.text.x =element_text(size=12),
    strip.text.y =element_text(size=12),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())
p2

ggarrange(p1,p2, 
          nrow=1, ncol=2, common.legend = TRUE)

ggsave(path=path.p, paste0("figS14B.pdf"), width= 11, height=5, onefile = TRUE)

dat.figS14<-dps0
## ..Write data table as csv file
fname<-paste("figS14B.csv")
write.csv(dat.figS14,file.path(path.ds, fname),row.names=FALSE)

## ..Plot subset of temporal dynamics for main figure
dps1<-subset(dps0, dila_ord==-4)
p1<-ggplot(dps1, aes(x=Transfer, y=Rel_abund, fill=Species))+
  geom_bar(stat="identity", position="stack",width=0.96) + 
  facet_grid(repb ~dilp_ord)+
  scale_fill_manual(values=c(colKp, colAlc, colPseu))+ 
  scale_y_continuous(breaks=c(0, 1))+
  labs(y= "relative abundance", x="Transfer")+
  theme_minimal()+
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y =element_text(size=16),
    axis.title.x = element_text(size=16), 
    strip.text.x =element_text(size=14),
    strip.text.y =element_text(size=14),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())
p1
ggsave(path=path.p, paste0("fig3C.pdf"), width= 10, height=3.5, onefile = TRUE)

dat.fig3C<-dps1
## ..Write data table as csv file
fname<-paste("fig3C.csv")
write.csv(dat.fig3C,file.path(path.ds, fname),row.names=FALSE)

