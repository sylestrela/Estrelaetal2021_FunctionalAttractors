" R script for Figures 1B (growth rates), S2, and S5 in Estrela et al (2021) 
Functional attractors in microbial community assembly.

note: for growth curves processing and max growth rate calculation, 
see companion script in /Rscripts_other/growth_curves_raw_fits"
#  ---------------------------------------------------------------------------
#.. Clear workspace
rm(list=ls())

#.. Load packages
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
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
dat.gr.isol<-read.csv('Estrela_2021_isolates_grmax.csv')

# ..get number of isolates  
length(unique(dat.gr.isol$SangerID))
# ..get number of communities  
length(unique(dat.gr.isol$CommunityReplicate))
# ..get number of genera  
length(unique(dat.gr.isol$genus))
# ..get number of families   
length(unique(dat.gr.isol$family))

# ..get number of isolates belonging to each family in gr
length(unique(subset(dat.gr.isol, family=='Pseudomonadaceae')$SangerID))
length(unique(subset(dat.gr.isol, family=='Enterobacteriaceae')$SangerID))
length(unique(subset(dat.gr.isol, family=='Aeromonadaceae')$SangerID))
length(unique(subset(dat.gr.isol, family=='Alcaligenaceae')$SangerID))
length(unique(subset(dat.gr.isol, family=='Comamonadaceae')$SangerID))
length(unique(subset(dat.gr.isol, family=='Moraxellaceae')$SangerID))

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

#.. Figure 1B- growth rates for Enterobacteriaceae and Pseudomonadaceae

d.sub<-dat.gr.isol
length(unique(d.sub$SangerID))
#.. re-order levels for plots
d.sub$cs <- factor(d.sub$cs, levels = c("glucose", "acetate",  "lactate", "succinate"))

#.. SELECT Pseu and Ent first 
dp<-subset(d.sub, family %in% c("Enterobacteriaceae", "Pseudomonadaceae"))

IQR_summary_EP<- dp %>%
  dplyr::group_by(cs,family) %>%
  dplyr::summarise(
    Min = round(min(gr_max),3),
    Max = round(max(gr_max),3),
    mean = round(mean(gr_max),3),
    std = round(sd(gr_max),3),
    Median = round(median(gr_max),3),
    IQRange = round(IQR(gr_max),3),
    Q25=quantile(gr_max, 0.25),
    Q75=quantile(gr_max, 0.75),
    n=n())
IQR_summary_EP<-as.data.frame(IQR_summary_EP)
IQR_summary_EP

set.seed(10)
xl<-""
yl<-expression( paste("max. growth rate ( ", h^-1,")"))

p<-ggplot(dp, aes(family, gr_max, fill=family))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(width = 0.2, size=2, alpha = 0.6, aes(color=family))+
  theme_classic()+
  facet_wrap(~cs, ncol=4)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.title.x =element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text.x = element_text(size = 18),
    strip.background = element_blank(),
    axis.ticks.x=element_blank(),
    legend.title = element_blank(),
    legend.text=element_text(size=14),
    legend.spacing.x = unit(0.5, 'cm'),
    legend.position="bottom")+
  scale_colour_manual(values=c(colEntbf, colPseuf))+
  scale_fill_manual(values=c(colEntbf, colPseuf))+
  scale_y_continuous(breaks=c(0.0, 0.6, 1.2))+
  labs(x=xl, y=yl)
p
p + stat_compare_means(method = "t.test", paired=FALSE, aes(label = ..p.signif..),label.y = 1.25,label.x=1.5)

ggsave(path=path.p, paste0("fig1B_gr.pdf"), width= 7, height=4, onefile = TRUE)

##.. stats
as.data.frame(compare_means( formula = gr_max ~ family, data=subset(dp,cs=='glucose') , method = "t.test", paired=FALSE))
t.test(gr_max ~ family, data = subset(dp,cs=='glucose'))
t.test(gr_max ~ family, data = subset(dp,cs=='acetate'))
t.test(gr_max ~ family, data = subset(dp,cs=='lactate'))
t.test(gr_max ~ family, data = subset(dp,cs=='succinate'))

dat.fig1B<-select(dp, cs,family, SangerID, gr_max)
## ..Write data table as csv file
fname<-paste("fig1B_gr.csv")
write.csv(dat.fig1B,file.path(path.ds, fname), row.names=FALSE)

#------------------------------------------------------------------------------
#.. Fig S2. Plot ALL families 
dp<-d.sub

IQR_summary_fam<- dp %>%
  dplyr::group_by(as.factor(family),cs) %>%
  dplyr::summarise(
    Min = round(min(gr_max),3),
    Max = round(max(gr_max),3),
    mean = round(mean(gr_max),3),
    std = round(sd(gr_max),3),
    Median = round(median(gr_max),3),
    IQRange = round(IQR(gr_max),3),
    Q25=quantile(gr_max, 0.25),
    Q75=quantile(gr_max, 0.75),
    n=n())
IQR_summary_fam<-as.data.frame(IQR_summary_fam)

#.. re-order families for plots
dp$family <- factor(dp$family, levels = c("Aeromonadaceae", "Enterobacteriaceae",  "Moraxellaceae", "Pseudomonadaceae", "Comamonadaceae", "Alcaligenaceae"  ))

xl<-""
yl<-"max. growth rate (h-1)"

p<-ggplot(dp, aes(family, gr_max, fill=family))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(width = 0.2, size=2, alpha = 0.5, aes(color=family))+
  theme_classic()+
  facet_wrap(cs~., scales= "free_y")+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=16),
    axis.title.y =element_text(size=18),
    axis.title.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title=element_text(size=16),
    legend.spacing.x = unit(0.5, 'cm'),
    legend.text=element_text(size=14),
    strip.text.x = element_text(size = 18),
    axis.ticks.x=element_blank())+
  scale_colour_manual(values=c(colAerof, colEntbf, colMoraxf, colPseuf,colComf,colAlcf))+
  scale_fill_manual(values=c(colAerof, colEntbf, colMoraxf, colPseuf,colComf,colAlcf))+
  labs(x=xl, y=yl)+
  ylim(0,NA)
p
ggsave(path=path.p, paste0("figS2.pdf"), width= 9, height=6, onefile = TRUE)

dat.figS2<-select(dp, cs,family, SangerID,gr_max)
## ..Write data table as csv file
fname<-paste("figS2.csv")
write.csv(dat.figS2,file.path(path.ds, fname),row.names=FALSE)

#------------------------------------------------------------------------------
#.. Fig. S5. Plot ALL genera of Ent family
dp<-subset(d.sub, family %in% c("Enterobacteriaceae"))

IQR_summary_gen<- dp %>%
  dplyr::group_by(genus,cs) %>%
  dplyr::summarise(
    Min = round(min(gr_max),3),
    Max = round(max(gr_max),3),
    mean = round(mean(gr_max),3),
    std = sd(gr_max),
    Median = median(gr_max),
    IQRange = IQR(gr_max),
    Q1=quantile(gr_max, 0.25),
    Q3=quantile(gr_max, 0.75))
IQR_summary_gen<-as.data.frame(IQR_summary_gen)

xl<-""
yl<-"max. growth rate (h-1)"
#..genera colours
colCit<-'olivedrab3'
colEnt<- 'palegreen4'
colEsc<-'steelblue1'
colKleb<-'royalblue1'
colRaou<- 'royalblue4'
colSerr<- 'steelblue4'
set.seed(10)  
p<-ggplot(dp, aes(genus, gr_max, fill=genus))+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(width = 0.2, size=2, alpha = 0.8, aes(color=genus))+
  theme_classic()+
  facet_wrap(cs~., scales= "free_y")+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=16),
    axis.title.y =element_text(size=18),
    axis.title.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.text.x = element_text(size = 18),
        axis.ticks.x=element_blank(),
    legend.title=element_text(size=16),
    legend.spacing.x = unit(0.5, 'cm'),
    legend.text=element_text(size=14))+
  labs(x=xl, y=yl)+
  ylim(-0.0001,NA)+
  scale_fill_manual(values=c(colCit, colEnt, colEsc,  colKleb, colRaou, colSerr ))+
  scale_colour_manual(values=c(colCit, colEnt, colEsc,  colKleb, colRaou, colSerr))
p
ggsave(path=path.p, paste0("figS5.pdf"), width= 9, height=6, onefile = TRUE)

dat.figS5<-select(dp, cs,genus, SangerID,gr_max)
## ..Write data table as csv file
fname<-paste("figS5.csv")
write.csv(dat.figS5,file.path(path.ds, fname), row.names=FALSE)



