" R script for Figures 1B and S4 (+stats) in Estrela et al (2021) 
Functional attractors in microbial community assembly."

#  ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)

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
dat.ga.isol<-read.csv('Estrela_2021_isolates_ph_OAs.csv')

unique(dat.ga.isol$SangerID)

# ..get number of isolates  
length(unique(dat.ga.isol$SangerID))
# ..get number of communities  
length(unique(dat.ga.isol$CommunityReplicate))
# ..get number of genera  
length(unique(dat.ga.isol$Genus))
# ..get number of families   
length(unique(dat.ga.isol$Family))

# ..get number of isolates belonging to each family in gr
length(unique(subset(dat.ga.isol, Family=='Pseudomonadaceae')$SangerID))
length(unique(subset(dat.ga.isol, Family=='Enterobacteriaceae')$SangerID))
length(unique(subset(dat.ga.isol, Family=='Aeromonadaceae')$SangerID))
length(unique(subset(dat.ga.isol, Family=='Alcaligenaceae')$SangerID))
length(unique(subset(dat.ga.isol, Family=='Comamonadaceae')$SangerID))
length(unique(subset(dat.ga.isol, Family=='Moraxellaceae')$SangerID))

#------------------------------------------------------------------------------
# .. Colour scheme 
# Family colors 
colAerof<-"deepskyblue4"
colAlcf<- "orange2"
colBacf<-"orange3"
colComf<- "firebrick" 
colEntbf<- "deepskyblue3"
colEntcf<-  "darkolivegreen3"
colMoraxf<- "darkorchid4"
colPseuf<-  "darkorchid2"
colSphingf<-"deeppink3"
colXantf<-"deeppink4"

#------------------------------------------------------------------------------
dp<-dat.ga.isol
dp<-melt(dp, id=c("time_hours", "CommunityReplicate", "Genus", "Family", "f_r_type", "SangerID"))
dp$variable<-revalue(dp$variable, c( "Glucose_perc"="Glucose (%)", "acetate_mM"= "acetate (mM)",
                                    "succinate_mM"= "succinate (mM)", "lactate_mM"= "lactate (mM)"))

#------------------------------------------------------------------------------
#. Figure 1B_pH_OAs. plot Enterobacteriaceae, Pseudomonadaceae, Aeromonadaceae and Moraxellaceae 
dp<-subset(dp, Family %in% c("Enterobacteriaceae", "Aeromonadaceae", "Pseudomonadaceae", "Moraxellaceae"))
dp<-subset(dp, variable %in% c('pH', 'acetate (mM)', 'succinate (mM)', 'lactate (mM)'))
ggplot(dp,aes(y=value, x=factor(time_hours), 
              group=factor(CommunityReplicate):factor(SangerID), colour=Family))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3)+ 
  facet_wrap(.~variable, scales="free_y", ncol=4)+
  stat_summary(aes(group=Family, color=paste( Family, "(mean)")), fun=mean, geom="line", linetype='dashed', size=1.5)+
  stat_summary(aes(group=Family, color=paste( Family, "(mean)")), fun=mean, geom="point", size=2)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.title.x =element_text(size=18),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text.x = element_text(size = 18),
    strip.background = element_blank())+
  scale_colour_manual(values=c(colAerof, colAerof, colEntbf, colEntbf, colMoraxf,colMoraxf, colPseuf, colPseuf))+
  scale_x_discrete(expand=c(0.02,0.02))+
  labs(x="time (hours)", y= "")

ggsave(path=path.p, paste0("fig1B_oa.pdf"), width= 12, height=3.5, onefile = TRUE)

dat.fig1B<-select(dp, time_hours, variable,Family,SangerID,value)
## ..Write data table as csv file
fname<-paste("fig1B_oa.csv")
write.csv(dat.fig1B,file.path(path.ds, fname),row.names=FALSE)

#------------------------------------------------------------------------------
#. calculate amount of acetate produced by Enterob family at each timepoint
ds<-subset(dp, Family %in% c("Enterobacteriaceae"))
ds<-select(ds,-CommunityReplicate,-ChromoType,-SangerID,-Genus)
nr<-2
stat_summ<- ds %>%
  dplyr::group_by(time_hours,variable) %>%
  dplyr::summarise(
    Min = round(min(value, na.rm = TRUE),nr),
    Max = round(max(value, na.rm = TRUE),nr),
    mean = round(mean(value,na.rm = TRUE),nr),
    std = round(sd(value,na.rm = TRUE),nr),
    median = round(median(value, na.rm = TRUE),nr),
    IQRange = round(IQR(value, na.rm = TRUE),nr),
    Q25= round(quantile(value, 0.25, na.rm = TRUE),nr),
    Q75= round(quantile(value, 0.75, na.rm = TRUE),nr),
    n=n())

d.stat<-as.data.frame(stat_summ)
select(d.stat, time_hours, variable,median, Q25,Q75)

#------------------------------------------------------------------------------
#. Fig. S4. plot all Genera of Enterobacteriaceae family
#..genera colours
colCit<-'olivedrab3'
colEnt<- 'palegreen4'
colEsc<-'steelblue1'
colKleb<-'royalblue1'
colRaou<- 'royalblue4'
colSerr<- 'steelblue4'

dp<-subset(dp, Family %in% c("Enterobacteriaceae"))
ggplot(dp,aes(y=value, x=factor(time_hours), 
              group=factor(CommunityReplicate):factor(SangerID), colour=Genus))+
  geom_point(alpha=0.3)+
  geom_line(alpha=0.3)+ 
  facet_wrap(~variable, scales="free_y", ncol=2)+
  stat_summary(aes(group=Genus, color=paste( Genus, "(mean)")), fun=mean, geom="line",  size=1.5)+
  stat_summary(aes(group=Genus, color=paste( Genus, "(mean)")), fun=mean, geom="point", size=2)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.title.x =element_text(size=18),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text.x = element_text(size = 18),
    strip.background = element_blank(),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16))+
  scale_colour_manual(values=c(colCit, colCit, colEnt, colEnt, colEsc,colEsc, colKleb, colKleb, colRaou, colRaou,colSerr,colSerr ))+
  scale_x_discrete(expand=c(0.02,0.02))+
  labs(x="time (hours)", y= "")
ggsave(path=path.p, paste0("figS4.pdf"), width= 9, height=6, onefile = TRUE)

dat.figS4<-select(dp, time_hours, variable,Genus,SangerID,value)
## ..Write data table as csv file
fname<-paste("figS4.csv")
write.csv(dat.figS4,file.path(path.ds, fname),row.names=FALSE)

