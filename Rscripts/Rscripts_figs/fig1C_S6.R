" R script for Figure 1C and S6 in Estrela et al (2021) 
Functional attractors in microbial community assembly.

"
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
dat.aceglu<-read.csv("Estrela_2021_fig1C_S6_acegluRF.csv")
dat.cfuml<-read.csv("Estrela_2021_fig1C_S6_cfumlRF.csv")

#----------------- #----------------- #----------------- #----------------- 
##.. FIG 1C. Plot 1 representative community for main text- C10R7
#----------------- #----------------- #----------------- #----------------- 
#. plot Glucose and Acetate both as y-axis in the same plot
colp<-'indianred3'
colp2<-'gray60'
#. plot Glucose and Acetate as facets
dp<-select(dat.aceglu, -mean.ratioRF,-sd.ratioRF,-sem.ratioRF)
dp.m<-melt(dp, id=c("time_hours", "CommunityReplicate"))

dp<-subset(dp, CommunityReplicate=="C10R7")
max.ace.p<-(max(dp$acetate_mM)/0.2+0.5)

p1<-ggplot(dp, aes(x=factor(time_hours), 
                   group=factor(CommunityReplicate))) +
  geom_line(aes(y = Glucose_perc, colour = "Glu"), linetype='dashed',size=1.5)+
  geom_point(aes(y = Glucose_perc, colour = "Glu"),size=2.7)+
  stat_summary(aes(group=CommunityReplicate, y = acetate_mM/max.ace.p, colour = "Ace"),
               fun='mean',
               geom="point", size=2.5) +
  stat_summary(aes(group=CommunityReplicate, y = acetate_mM/max.ace.p, colour = "Ace"), fun='mean', geom='line', linetype='dashed', size=1.5)+
  scale_y_continuous(sec.axis = sec_axis(~.*max.ace.p, name = "acetate (nmol/uL)",breaks=c(0, 6,12)),breaks=c(0, 0.1, 0.2))+
  scale_x_discrete(expand=c(0.01,0.01))+
  theme_classic()+
  facet_grid(CommunityReplicate~.)+
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size=16, colour="black", angle = 0, hjust = 0.5),
    strip.text.y = element_blank(),
    axis.title.y.right = element_text(colour = colp, size=18),
        axis.ticks.y.right = element_line(color = colp),
        axis.text.y.right = element_text(color = colp, size=16),
        axis.line.y.right = element_line(color = colp ),
        axis.title.y.left = element_text(colour = 'gray30',size=18),
        axis.ticks.y.left = element_line(color = colp2),
        axis.text.y.left = element_text(color = colp2, size=16),
        axis.line.y.left = element_line(color = colp2 ),
        legend.position = "none")+
  scale_colour_manual(values=c( colp,colp2))+
  labs(x="time (hours)", y= "glucose (%)")
p1
ggsave(path=path.p, paste0("fig1C_glu_ace.pdf"), width= 4.5, height=3, onefile = TRUE)

dat.fig1C<-dp
## ..Write data table as csv file
fname<-paste("fig1C_glu_ace.csv")
write.csv(dat.fig1C,file.path(path.ds, fname), row.names=FALSE)

#----------------- 
#. plot R/F ratio (mean+sd)
dp<-select(dat.aceglu, -acetate_mM,-Glucose_perc)

dp2<-subset(dp, CommunityReplicate=="C10R7")
colp.rf<- 'darkorchid4' 

p2<-ggplot(dp2, aes(x = factor(time_hours), y = mean.ratioRF, group=factor(CommunityReplicate))) +
  geom_point(colour=colp.rf, size=2.7)+
  geom_line(colour=colp.rf, size=1.5)+
  geom_errorbar(aes(ymin=mean.ratioRF-sd.ratioRF, ymax=mean.ratioRF+sd.ratioRF), width=.1,colour=colp.rf)+
  facet_grid(CommunityReplicate~.)+
  theme_classic()+
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.y.left = element_line(color=colp.rf),
    axis.ticks.y.left = element_line(color=colp.rf),
    strip.text.y = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18,colour=colp.rf),
    axis.text.y = element_text(size = 16, colour=colp.rf),
    axis.text.x = element_text(size=16, colour="black", angle = 0, hjust = 0.5),
    legend.position = "none")+
  scale_y_continuous(breaks=c(0.02, 0.05, 0.08))+
  scale_x_discrete(expand=c(0.01,0.01))+
  labs(x="time (hours)", y= "R/F")
p2

ggsave(path=path.p, paste0("fig1C_RF.pdf"), width= 4, height=3, onefile = TRUE)

dat.fig1C<-select(dp2,-sem.ratioRF)
## ..Write data table as csv file
fname<-paste("fig1C_RF.csv")
write.csv(dat.fig1C,file.path(path.ds, fname), row.names=FALSE)

#----------------- #----------------- #----------------- #----------------- 
##.. Figure S6. All communities.
#----------------- #----------------- #----------------- #----------------- 
#. plot Glucose and Acetate both as y-axis in the same plot
colp<-'indianred3'
colp2<-'gray40'
#. plot Glucose and Acetate as facets
dp<-select(dat.aceglu, -mean.ratioRF,-sd.ratioRF,-sem.ratioRF)
dp.m<-melt(dp, id=c("time_hours", "CommunityReplicate"))

max.ace.p<-(max(dp$acetate_mM)/0.2+0.5)

p1<-ggplot(dp, aes(x=factor(time_hours), 
                   group=factor(CommunityReplicate))) +
  geom_line(aes(y = Glucose_perc, colour = "Glu"), linetype='dashed')+
  geom_point(aes(y = Glucose_perc, colour = "Glu"),size=1.5)+
  stat_summary(aes(group=CommunityReplicate, y = acetate_mM/max.ace.p, colour = "Ace"),
               fun='mean',
               geom="point", size=1.5) +
  stat_summary(aes(group=CommunityReplicate, y = acetate_mM/max.ace.p, colour = "Ace"), fun='mean', geom='line', linetype='dashed')+
  scale_y_continuous(sec.axis = sec_axis(~.*max.ace.p, name = "acetate (nmol/uL)"),breaks=c(0, 0.2))+
  theme_classic()+
  facet_grid(CommunityReplicate~.)+
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size=16, colour="black", angle = 0, hjust = 0.5),
    strip.text.y = element_blank(),
    axis.title.y.right = element_text(colour = colp, size=18),
        axis.ticks.y.right = element_line(color = colp),
        axis.text.y.right = element_text(color = colp, size=12),
        axis.line.y.right = element_line(color = colp ),
        axis.title.y.left = element_text(colour = 'gray30',size=18),
        axis.ticks.y.left = element_line(color = colp2),
        axis.text.y.left = element_text(color = colp2, size=12),
        axis.line.y.left = element_line(color = colp2 ),
        legend.position = "none")+
  scale_colour_manual(values=c( colp,colp2))+
  labs(x="time (hours)", y= "glucose (%)")
p1

dat.figS6p1<-dp
## ..Write data table as csv file
fname<-paste("figS6_glu_ace.csv")
write.csv(dat.figS6p1,file.path(path.ds, fname),row.names=FALSE)

#----------------- 
#. plot R/F ratio (mean+sd)
dp2<-select(dat.aceglu, -acetate_mM,-Glucose_perc,-sem.ratioRF)

p2<-ggplot(dp2, aes(x = factor(time_hours), y = mean.ratioRF, group=factor(CommunityReplicate))) +
  geom_point(alpha=0.9,colour=colp.rf)+
  geom_line(colour=colp.rf)+
  geom_errorbar(aes(ymin=mean.ratioRF-sd.ratioRF, ymax=mean.ratioRF+sd.ratioRF), width=.2,colour=colp.rf)+
  facet_grid(CommunityReplicate~.,scales="free_y")+
  theme_classic()+
  ylim(0,NA)+
  theme(
    axis.line.y.left = element_line(color=colp.rf),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_text(size = 18,colour=colp.rf),
    axis.title.x = element_text(size = 18),
    axis.text.y = element_text(size = 12, colour=colp.rf),
    axis.text.x = element_text(size=16, colour="black", angle = 0, hjust = 0.5),
    legend.position = "none")+
  scale_colour_manual(values=c( colp,colp2))+
  labs(x="time (hours)", y= "R/F")
p2

dat.figS6p2<-dp2
## ..Write data table as csv file
fname<-paste("figS6_RF.csv")
write.csv(dat.figS6p2,file.path(path.ds, fname),row.names=FALSE)


#----------------- 
##.. Plot cfu/ml + Stat summary with mean, standard deviation, standard error of the mean, and a (default 95%) confidence interval, IQR
dp<-melt(dat.cfuml, id.vars = c("time_hours", 'CommunityReplicate', 'rep'))
names(dp)[names(dp) == "value"] <- "cfu.ml"
ds<-dp
nr<-6
StatSum_cfuml<- ds %>%
  dplyr::group_by(CommunityReplicate,time_hours,variable) %>%
  dplyr::summarise(
    Min = round(min(cfu.ml, na.rm = TRUE),nr),
    Max = round(max(cfu.ml, na.rm = TRUE),nr),
    mean = round(mean(cfu.ml,na.rm = TRUE),nr),
    std = round(sd(cfu.ml,na.rm = TRUE),nr),
    Median = round(median(cfu.ml, na.rm = TRUE),nr),
    IQRange = round(IQR(cfu.ml, na.rm = TRUE),nr),
    Q25= round(quantile(cfu.ml, 0.25, na.rm = TRUE),nr),
    Q75= round(quantile(cfu.ml, 0.75, na.rm = TRUE),nr))

dp<-as.data.frame(StatSum_cfuml)

colEntbf<- "deepskyblue3"
colPseuf<-  "darkorchid2"
dp<-select(dp, CommunityReplicate, time_hours, variable,mean,std)
#.. note: log10 of error bars for 0 values does not so need to remove those datapoints for plotting
p3<-ggplot(dp, aes(x = factor(time_hours), y = mean, group=factor(CommunityReplicate):factor(variable), colour=variable)) +
  geom_point(alpha=0.9)+
  geom_line(alpha=0.5)+
  geom_errorbar(dat=subset(dp, mean>0), aes(ymin=mean-std, ymax=mean+std), width=.2)+
  facet_grid(CommunityReplicate~., scales="free_y")+
  theme_classic()+
  labs(y="CFU/ml", x='time (hours)')+
  scale_colour_manual(values=c(colEntbf,colPseuf))+
  theme(
    axis.text.y = element_text(size = 12, colour="black"),
    axis.text.x = element_text(size=16, colour="black", angle = 0, hjust = 0.5),
    axis.title.y =element_text(size=18),
    axis.title.x = element_text(size = 18),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())+
    theme(legend.position = "none")+
  scale_y_log10()
p3

dat.figS6p3<-dp
## ..Write data table as csv file
fname<-paste("figS6_cfuml.csv")
write.csv(dat.figS6p3,file.path(path.ds, fname),row.names=FALSE)


p1c.all<-ggarrange(p1,p2,p3,
                   ncol = 3, nrow = 1)
p1c.all
ggsave(path=path.p, paste0("figS6.pdf"), width= 8, height=12, onefile = TRUE)

