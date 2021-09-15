" R script for Figure 1A (+ stats) in Estrela et al (2021) 
Functional attractors in microbial community assembly."

#  ---------------------------------------------------------------------------

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
dat.ES<-read.csv('Goldford_2018_16S_data_table_glucose_RelAbund_ALL.csv')

# ..SELECT minimum abundance
min.ab<-0.0
df.glu<-subset(dat.ES, Relative_Abundance >min.ab)
colnames(df.glu) 

#------------------------------------------------------------------------------
#.. get number of ESV for each community
df.glu.num.esv<-ddply(df.glu, ~InocRep, summarise, num_distinct_esv=length(unique(ESV)))
max.esv<-max(df.glu.num.esv$num_distinct_esv)
min.esv<-min(df.glu.num.esv$num_distinct_esv)
max.esv
min.esv

#  ---------------------------------------------------------------------------
# .. Colour scheme 
# Family colors 
colEntbf<- "deepskyblue3"
colPseuf<-  "darkorchid2"
colAerof<-"deepskyblue4"
colMoraxf<-"darkorchid4"
colOther<-"gray85"
#  ---------------------------------------------------------------------------

#.. Figure 1A. 
xl<-"community"
yl<-"relative abundance"

#.. Plot Enterobacteriaceae and Aeromonadaceae in blue; Pseudomonadaceae and Moraxellaceae in purple; other families in gray
#.create 2 groups (EPAM and other)
dp.other<-subset(df.glu, !(Family %in% c("Enterobacteriaceae", 'Aeromonadaceae',  "Pseudomonadaceae",  "Moraxellaceae")))
dp.other$Family_p<-'other'
dp.dom<-subset(df.glu, Family %in% c("Enterobacteriaceae", 'Aeromonadaceae',  "Pseudomonadaceae",  "Moraxellaceae"))
dp.dom$Family_p<-dp.dom$Family

dp2<-rbind(dp.dom,dp.other)
dp2$Family_p <- factor(dp2$Family_p, levels = c("other", "Enterobacteriaceae", 'Aeromonadaceae',  "Pseudomonadaceae",  "Moraxellaceae"))

p<-ggplot(dp2, aes(x=InocRep,y=Relative_Abundance ,fill=Family_p))
p + 
  geom_bar(stat="identity", position="stack", color="black") +
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  labs(x=xl, y=yl, caption="")+
  scale_fill_manual(values=c(colOther, colEntbf, colAerof, colPseuf, colMoraxf ), name='Family')
ggsave(path=path.p, paste0("fig1A.pdf"), width= 14, height=3.5, onefile = TRUE)

dat.fig1A<-select(dp2, InocRep, Relative_Abundance,Family)

## ..Write data table as csv file
fname<-paste("fig1A.csv")
write.csv(dat.fig1A,file.path(path.ds, fname),row.names=FALSE)

#  ---------------------------------------------------------------------------
#.. Calculate ratio of R/F and of Pseu/Ent for each community
df.sub<-select(df.glu, -ESV, -Genus,-Carbon_Source,-Inoculum,-Replicate,-Transfer,-Experiment,-ESV_ID,-Order)

#.. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fr<-data.table(df.sub)[, list(f_r_sum=sum(Relative_Abundance)), 
                            by=list(f_r_type, InocRep)]

#.. sum the relative abundance of each family for each sample (i.e. each rep x treatment x transfer x...)
dsum_fam<-data.table(df.sub)[, list(fam_sum=sum(Relative_Abundance)), 
                             by=list( f_r_type, InocRep, Family)]

dat.glu.RF<-as.data.frame(spread(dsum_fr,f_r_type, f_r_sum))
dat.glu.RF$ratioRF<-(dat.glu.RF$R/dat.glu.RF$F)
median(dat.glu.RF$ratioRF)
mean(dat.glu.RF$ratioRF)

#. ratio of Pseudomonadaceae/Enterobacteriacea 
dsum_famPE<-subset(dsum_fam, Family %in% c("Pseudomonadaceae","Enterobacteriaceae" ))
dsum_famPE<-select(dsum_famPE, -f_r_type)
dat.glu.PE<-as.data.frame(spread(dsum_famPE,Family, fam_sum))
dat.glu.PE$ratioPE<-(dat.glu.PE$Pseudomonadaceae/dat.glu.PE$Enterobacteriaceae)
median(dat.glu.PE$ratioPE)
mean(dat.glu.PE$ratioPE)

head(dat.glu.RF)
dat.glu.RF<-select(dat.glu.RF, -UN)
head(dat.glu.PE)

#.. Write data table with dat.glu.PE and dat.glu.RF as csv file to be used for boxplot in figure 1D and fig SXX
fname<-paste("Goldford_2018_glu_ratioPE.csv", sep="")
write.csv(dat.glu.PE,file.path(path.dp, fname),row.names=FALSE)

fname<-paste("Goldford_2018_glu_ratioRF.csv", sep="")
write.csv(dat.glu.RF,file.path(path.dp, fname),row.names=FALSE)

## ---------------------------------------------------------------------------
#. calculate IQR 
#.. R/F ratio
dp<-select(dat.glu.RF,-F,-R)
IQR_summary<- dp %>%
  summarise(Min = min(ratioRF),
            Max = max(ratioRF),
            Median = median(ratioRF),
            Q25 =quantile(ratioRF, 0.25),
            Q75 =quantile(ratioRF, 0.75),
            IQRange = IQR(ratioRF))
IQR_summary

#. Pseu/Ent ratio
dp<-select(dat.glu.PE,-Enterobacteriaceae,-Pseudomonadaceae)
IQR_summary<- dp %>%
  summarise(Min = min(ratioPE),
            Max = max(ratioPE),
            Median = median(ratioPE),
            Q25 =quantile(ratioPE, 0.25),
            Q75 =quantile(ratioPE, 0.75),
            IQRange = IQR(ratioPE))
IQR_summary

