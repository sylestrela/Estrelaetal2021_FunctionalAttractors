" R script for Figure 2D (timeseries overlaid on top of the derivative of the potential) 
in Estrela et al (2021) Functional attractors in microbial community assembly."

#  ---------------------------------------------------------------------------
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
df.all<-read.csv('Estrela_2021_16S_data_table_RelAbund_ALL.csv')
setwd(path.dp)
d.esvs<-read.csv('top6_ESVs_fig2B.csv')

df.all$experiment_s<-paste0(df.all$Treatment,'_', df.all$Inoc_size)

##.. SELECT A and P at T18 
esvAlc<-as.character(subset(d.esvs,Genus=='Alcaligenes')$OTU)
esvP1<-as.character(subset(d.esvs,Genus=='Pseudomonas')$OTU)
dat.ap<-subset(df.all, OTU %in% c(esvAlc, esvP1))
dat.ap$esv_id[dat.ap$OTU==esvAlc]<-'A'
dat.ap$esv_id[dat.ap$OTU==esvP1]<-'P'
dat.ap<-select(dat.ap, -OTU,-Carbon,-Family,-Treatment,-Inoc_size)

min.ab<-0.0
df.sub<-subset(dat.ap, Abundance >min.ab)

## .. SELECT replicates we have timeseries for in each treatment
#..No_migration_4
p1<-c(2,3, 5,6,14,15,18,20, 26)
p2<-c(31,32,36,44,45,49,54,63, 75, 77)
df.sub1<-subset(df.sub, subset = Rep_num %in% c(p1,p2) & experiment_s=='No_migration_4')
#..No_migration_40
p3<-c(1,21,33,42)
df.sub2<-subset(df.sub, subset = Rep_num %in% c(p3) & experiment_s=='No_migration_40')
#..global_migration_4
p4<-c(21,25,31,42,49,52,63,84)
df.sub3<-subset(df.sub, subset = Rep_num %in% c(p4) & experiment_s=='Global_migration_4')

df.sub<-rbind(df.sub1,df.sub2,df.sub3)

#.. calculate number of timeseries
length(c(p1,p2,p3,p4))

##-----------------------------------------------------
#.. a) Timeseries (lineplots) of log abundance overlaid with potential derivative

#. Colour scheme for plots- ESVs color
colAlc<-"#E69F00"
colPseu<-  "darkorchid4" 
colKp<-blues9[7]
colKm<-blues9[5]

# .. SELECT strain - Alcaligenes
stn<-"Alc"
colp<-colAlc
dt.p<-subset(df.sub, Genus=="Alcaligenes")
dfp.sub<-dt.p

##.. plot derivative of potential as background
x = seq(-4.4,0,by=0.001)
pot.a=(-(1/2)* log(2.2131661111317307*exp(-33.49311531236609* (0.590556 + x)^2) + 0.280318277722824*exp(-2.378165480442605* (3.006816 + x)^2)))

spl <- smooth.spline(x, y=pot.a)
pred <- predict(spl)
ycs.prime <- diff(pot.a)/diff(x)
pred.prime <- predict(spl, deriv=1)

ab<-as.data.frame(x[1:4400])
u.deriv<-as.data.frame(ycs.prime)
df<-cbind(ab,u.deriv)
tr<-seq(1,18,by=1)
df2<-merge(df,tr)
names(df2)<-c('ab', 'u_deriv', 'Transfer')

##-- create two groups for plotting: trajectories that switch to A state and trajectories that don't 
##- roots were calculated using mathematica (see companion script)
mindv<-(-1.17738)
root1<-(-0.590556)
root2<-(-3.00682)

dfp.sub$experiment_s_rep<-paste0(dfp.sub$experiment_s,'_',dfp.sub$Rep_num)
g1<-subset(dfp.sub, Transfer==18 & log10(dfp.sub$Abundance)>mindv)

g1<-subset(dfp.sub, experiment_s_rep %in% g1$experiment_s_rep)
g2<-subset(dfp.sub, !(experiment_s_rep %in% g1$experiment_s_rep))

g1$group<-'switch'
g2$group<-'noswitch'

dfp.sub<-rbind(g1,g2)

b<-c(-10,0,5,10)
min.ud<-min(abs(df2$u_deriv))

p<-ggplot() + 
  geom_tile(data=df2, aes(x=Transfer, y=ab, fill=abs(u_deriv))) +
   geom_line(data=subset(dfp.sub, group=='switch'), aes(x=Transfer, y=log10(Abundance), colour=factor(Genus), group=factor(Rep_num):factor(experiment_s)), alpha=0.9 )+
  geom_line(data=subset(dfp.sub, group=='noswitch'), aes(x=Transfer, y=log10(Abundance), colour=factor(Genus), group=factor(Rep_num):factor(experiment_s)), linetype='dashed', alpha=0.9 )+
  scale_x_discrete(limits=c(1, 18))+
  scale_colour_manual(values=colp) +
  labs(y=expression(A*phantom(x)*(log[10])),
       x='transfer')+
  ylim(-4.4,0)+
  theme_minimal()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x = element_text(size=18,hjust=0.5),
        axis.title=element_text(size=20),
        axis.title.x = element_text(vjust = 4),
        strip.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "gray", size = 0.2, linetype = "solid"),
        axis.ticks=element_line(colour = "gray"),
        legend.title=element_text(size=20),
        legend.position="top")+
  guides(colour=FALSE)
p
#.. add fill gradient
col.gl<-"gray90"
col.gd<-"white"
col1<-'#737373'
valm<-0.05
p2d.a<-p+scale_fill_gradientn(colours = c(col1, col.gl, col.gd),
                              name = "|U'(A)|",
                           values = c(0, valm, 1), breaks=c(min(abs(df2$u_deriv)),8,16), labels=c(0,8,16))
p2d.a

##-----------------------------------------------------
#. Plot Pseudomonas
# .. SELECT strain
stn<-"Pseu"
colp<-colPseu
dt.p<-subset(df.sub, Genus=="Pseudomonas")
dfp.sub<-dt.p

##.. plot derivative of potential as background for P
x = seq(-4.4,0,by=0.001)
pot.a=(-(1/2)* log(1.13421*exp(-6.63602 * (1.08247 + x)^2) + 0.117626*exp(-0.901286 * (3.16649  + x)^2)))
#plot(x,pot.a)

spl <- smooth.spline(x, y=pot.a)
pred <- predict(spl)
ycs.prime <- diff(pot.a)/diff(x)
pred.prime <- predict(spl, deriv=1)

ab<-as.data.frame(x[1:4400])
u.deriv<-as.data.frame(ycs.prime)
df<-cbind(ab,u.deriv)
tr<-seq(1,18,by=1)
df2<-merge(df,tr)
names(df2)<-c('ab', 'u_deriv', 'Transfer')

##-- create two groups for plotting: trajectories that switch to P state and trajectories that don't 
##- roots were calculated using mathematica (see companion script)
mindv<-(-1.97225)
root1<-(-1.08306)
root2<-(-3.16649)

dfp.sub$experiment_s_rep<-paste0(dfp.sub$experiment_s,'_',dfp.sub$Rep_num)
g1<-subset(dfp.sub, Transfer==18 & log10(dfp.sub$Abundance)>mindv)

g1<-subset(dfp.sub, experiment_s_rep %in% g1$experiment_s_rep)
g2<-subset(dfp.sub, !(experiment_s_rep %in% g1$experiment_s_rep))

g1$group<-'switch'
g2$group<-'noswitch'

dfp.sub<-rbind(g1,g2)
colp<-'darkorchid4'
p<-ggplot() + 
  geom_tile(data=df2, aes(x=Transfer, y=ab, fill=abs(u_deriv)))+
   geom_line(data=subset(dfp.sub, group=='switch'), aes(x=Transfer, y=log10(Abundance), colour=factor(Genus), group=factor(Rep_num):factor(experiment_s)), alpha=0.9 )+
  geom_line(data=subset(dfp.sub, group=='noswitch'), aes(x=Transfer, y=log10(Abundance), colour=factor(Genus), group=factor(Rep_num):factor(experiment_s)), linetype='dashed', alpha=0.9 )+
  scale_x_discrete(limits=c(1, 18))+
  scale_colour_manual(values=colp) +
  labs(y=expression(P*phantom(x)*(log[10])),
       x='transfer')+
  ylim(-4.4,0)+
  theme_minimal()+
  theme(axis.text.y=element_text(size=18),
        axis.text.x = element_text(size=18,hjust=0.5),
        axis.title=element_text(size=20),
        axis.title.x = element_text(vjust = 4),
        strip.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "gray", size = 0.2, linetype = "solid"),
        axis.ticks=element_line(colour = "gray"),
        legend.title=element_text(size=20),
        legend.position="top")+
  guides(colour=FALSE)
p
#.. add fill gradient
col.gl<-"gray90"
col.gd<-"white"
col1<-'#737373'
valm<-0.05
p2d.p<-p+scale_fill_gradientn(colours = c(col1, col.gl, col.gd),
                              name = "|U'(P)|",
                               values = c(0, valm, 1), breaks=c(min(abs(df2$u_deriv)),8,16), labels=c(0,8,16))
p2d.p

ggarrange(p2d.a,p2d.p, ncol=2)
ggsave(path=path.p, paste0("fig2D_timeseries_Uheatmap.png"), width= 12, height=6, dpi=400)

## ..Write data table as csv file
dat.fig2D<-select(dfp.sub, Transfer, Rep_num, esv_id, Abundance)
fname<-paste("fig2D_timeseries.csv")
write.csv(dat.fig2D,file.path(path.ds, fname),row.names=FALSE)

