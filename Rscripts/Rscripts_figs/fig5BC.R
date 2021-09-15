" R script for Figure 5B and 5C (+stats) in Estrela et al (2021) 
Functional attractors in microbial community assembly.
 "

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
library(tidyr)
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
df.all$experiment<-paste0(df.all$Treatment,'_', df.all$Inoc_size, '_', df.all$Transfer)

# .. SELECT all treatments, T12 and T18
tt<-c("No_migration_4_18", "No_migration_40_18", "Global_migration_4_12", "Global_migration_4_18", "Parent_migration_4_12", "Parent_migration_4_18")
df.all.s<-subset(df.all, experiment %in% tt)
# 
df.all.s<-select(df.all.s, -Kingdom,-Phylum,-Class,-Order)
head(df.all.s)

#.. give a single ID to each OTU
d1<-as.data.frame(unique(df.all.s$OTU))
d2<-as.numeric(rownames(as.data.frame(unique(df.all.s$OTU))))
d3<-cbind(d1,d2)
names(d3)<-c('OTU', 'esv_id')
df.all.s<-merge(df.all.s, d3)
df.all.s$esv_id<-paste0(df.all.s$Genus, '_', df.all.s$esv_id)

#..plot all
min.ab<-0.01
dp.abund<-subset(df.all.s, Abundance>min.ab)
dp.rare<-subset(df.all.s, Abundance<=min.ab)
dp.rare$Family<-'Other' #..note: all taxa below min.ab are shown as other regarless of family
dp<-rbind(dp.abund, dp.rare)

#..order for plotting
dp$experiment = factor(dp$experiment, 
                           levels=c("No_migration_4_18", "No_migration_40_18", "Global_migration_4_12", "Global_migration_4_18", "Parent_migration_4_12", "Parent_migration_4_18"))


colOther<-  "gray90"
colKp<-'steelblue4'
colKm<-'#6993BC'

colAlcf<-'#FFB94C'
colComf<-'firebrick' 
colEntcf<-'#8EB8B6'
colEntbf<- "deepskyblue3"
colPseuf<-  "darkorchid2"
colLachf<-'#BBD4D3'

#.. SELECT migration treatments only
colfam<-c(colAlcf, colComf, colEntbf, colEntcf, colLachf, colPseuf, colOther)
dpp<-subset(dp, Transfer==18 & Treatment%in% c('Global_migration',  'Parent_migration'))
dpp$treat.lab<-dpp$Treatment
levels(dpp$treat.lab)[levels(dpp$treat.lab)=='Global_migration'] <- 'global migration'
levels(dpp$treat.lab)[levels(dpp$treat.lab)=='Parent_migration'] <- 'regional migration'

p<- ggplot(dpp, aes(x=factor(Rep_num), y=Abundance,  fill=Family))
p + 
  geom_bar(color="black", stat="identity", position="stack")+  
  labs(y= "relative abundance", x="replicate community")+
  facet_grid(treat.lab~.)+
  scale_x_discrete(limits=c(1:93))+
  scale_y_continuous(breaks=seq(0,1,1))+
  scale_fill_manual(values=colfam)+
  theme_minimal()+
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.title=element_text(size=16),
        strip.text.y = element_text(size = 14),
        legend.text=element_text(size = 12),
        legend.title=element_text(size = 12),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.caption = element_text(hjust = 0))
ggsave(path=path.p, paste0("fig5B.pdf"), width= 14, height=5.5, onefile = TRUE)

## ..Write data table as csv file
fname<-paste("fig5B.csv")
ds<-select(dpp, Treatment, Transfer, Rep_num, Family,Abundance)
write.csv(ds,file.path(path.ds, fname),row.names=FALSE)

#  ---------------------------------------------------------------------------
# .. Calculate R/F ratio for all (at T18)
# ..note: here I use dp where taxa <0.01 appear as Family=other so that R/F ratio is not biased
#.. Create new column with F or R assignment
unique(df.all.s$Family)

dp2<-subset(df.all.s, Transfer==18)

dp2$experiment_rep<-paste0(dp2$experiment,'_',dp2$Rep_num)

dp.kas<-select(subset(dp2, Genus %in% c('Alcaligenes', 'Pseudomonas','Klebsiella', 'Klebsiella_m')), experiment_rep, esv_id,Abundance)
dp.kas<-spread(dp.kas, key='esv_id', value='Abundance', fill=NA)
dp.kas[is.na(dp.kas)] <- 0

# .. Create groups (K only, KA and KP)
ga<-subset(dp.kas, Alcaligenes_62>=0.01 & (Pseudomonas_23+Pseudomonas_52)<0.01)
gp<-subset(dp.kas, Alcaligenes_62<0.01 & (Pseudomonas_23+Pseudomonas_52)>=0.01)
#gk<-subset(dp.kas, Alcaligenes_62<0.1 & Pseudomonas_23<0.02)
nAdPd<-rbind(ga,gp)
gk<-subset(dp.kas, !(experiment_rep %in% nAdPd$experiment_rep))

ga<-subset(dp2, experiment_rep %in% ga$experiment_rep)
gp<-subset(dp2, experiment_rep %in% gp$experiment_rep)
gk<-subset(dp2, experiment_rep %in% gk$experiment_rep)
go<-subset(dp2, !(experiment_rep %in% c(ga$experiment_rep, gp$experiment_rep, gk$experiment_rep)))

ga$Group<-'KA'
gp$Group<-'KP'
gk$Group<-'K'

dkap<-rbind(ga,gp,gk)

dkap$f_r_type <- with(dkap, 
                        ifelse(Family == "Enterobacteriaceae" | Family == "Enterococcaceae" | Family=="Aeromonadaceae" | Family=="Lachnospiraceae",
                               "F", ifelse(Family == "Alcaligenaceae" | Family == "Pseudomonadaceae" | Family=="Moraxellaceae" | Family=='Xanthomonadaceae' | Family=="Comamonadaceae",
                                           "R", 'other')))
head(dkap)

#.. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fr<-data.table(dkap)[, list(f_r_sum=sum(Abundance)), 
                        by=list(f_r_type, experiment, Rep_num,Group)]

#.. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fam<-data.table(dkap)[, list(fam_sum=sum(Abundance)), 
                         by=list(Family, experiment,Rep_num,Group)]

#. abundance of each family
dsum_fam2 <- as.data.frame(dsum_fam) %>% 
  dplyr::group_by(Family,experiment,Group) %>% 
  dplyr::summarize(mean = mean(fam_sum))
dsum_fam2<-as.data.frame(dsum_fam2)

#. note: for cases with no F or no NF, we need to reshape the dataframe so their rel. ab is 0 or 1.

#. ratio of Respirer/Fermenter
dsum_fr2 <- dsum_fr %>% 
  spread(f_r_type, f_r_sum) %>%
  group_by(experiment, Rep_num,Group) %>% 
  mutate(ratioRF = R/F)
dsum_fr2<-as.data.frame(dsum_fr2)

dsum_fr2$R[dsum_fr2$F==1]<-0
dsum_fr2$ratioRF[dsum_fr2$F==1]<-0

median(dsum_fr2$ratioRF)
mean(dsum_fr2$ratioRF)

# .. Calculate median and IQR for each experiment
ds<-dsum_fr2
ds[is.na(ds)] <- 0

nr<-3
stat_summ<- ds %>%
  dplyr::group_by(experiment) %>%
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

# .. Add as a line the R/F ratio of Emergent Simplicity paper (fig1A): R/F=0.29, Q1=0.17, Q3=0.69
es_rf_med<-0.29
es_rf_q1<-0.17
es_rf_q3<-0.69
  
xl<-""
yl<-"R/F"

#.. SELECT all low inocula treatments
dss<-subset(ds, experiment!="No_migration_40_18")
#..order for plotting
dss$experiment = factor(dss$experiment, 
                       levels=c("No_migration_4_18", "Global_migration_4_18", "Parent_migration_4_18"))
#. ESVs color
colKp<-'steelblue4'
colPseu<-'#A17BB9'
colAlc<-'#FFB94C'
p<- ggplot(dss, aes(x=factor(experiment), ratioRF))+
  geom_boxplot(outlier.shape = NA) + 
  annotate("rect", ymin = es_rf_q1, ymax = es_rf_q3, xmin = 0, xmax = Inf,
           alpha = .3, fill='gray85')+
  geom_jitter(data=dss, aes(colour=Group),width = 0.2, size=2.5, alpha=0.85)+
  theme_classic()+
  scale_colour_manual(values=c(colKp, colAlc, colPseu))+
  theme(
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    axis.title.y =element_text(size=20),
    axis.title.x = element_text(size = 16), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_text(size = 16), legend.position = 'top',
    legend.text = element_text(size = 16))+
  labs(x=xl, y=yl)+
  scale_x_discrete(labels=c("No_migration_4_18" = "no
migration", 
                            "Global_migration_4_18" = " global
  migration",
                            "Parent_migration_4_18" = "regional
migration"))+
  scale_y_log10(breaks = c(0.0001, 0.01, 0.1, 1,10),limits = c(0.0001, 4))
p

ggsave(path=path.p, paste0("fig5C.pdf"), width= 6, height=6, onefile = TRUE)

## ..Write data table as csv file
fname<-paste("fig5C.csv")
ds<-select(dss, experiment, Rep_num, Group, ratioRF)
write.csv(ds,file.path(path.ds, fname),row.names=FALSE)

#.. Calculate R/F ratio for K+ only group
kp<-subset(dss, experiment=="No_migration_4_18" & Group=='K')
median(kp$ratioRF)

nr<-4
stat_summ<- kp %>%
  dplyr::summarise(
    Min = round(min(ratioRF, na.rm = TRUE),nr),
    Max = round(max(ratioRF, na.rm = TRUE),nr),
    mean = round(mean(ratioRF,na.rm = TRUE),nr),
    std = round(sd(ratioRF,na.rm = TRUE),nr),
    median = round(median(ratioRF, na.rm = TRUE),nr),
    Q25= round(quantile(ratioRF, 0.25, na.rm = TRUE),nr),
    Q75= round(quantile(ratioRF, 0.75, na.rm = TRUE),nr),
    count=n())
stat_summ<-as.data.frame(stat_summ)
stat_summ

#.. Calculate R/F ratio for each group of no migration low inocula
kp<-subset(dss, experiment=="No_migration_4_18")
median(kp$ratioRF)

nr<-4
stat_summ<- kp %>%
  dplyr::group_by(experiment, Group) %>%
  dplyr::summarise(
    Min = round(min(ratioRF, na.rm = TRUE),nr),
    Max = round(max(ratioRF, na.rm = TRUE),nr),
    mean = round(mean(ratioRF,na.rm = TRUE),nr),
    std = round(sd(ratioRF,na.rm = TRUE),nr),
    median = round(median(ratioRF, na.rm = TRUE),nr),
    Q25= round(quantile(ratioRF, 0.25, na.rm = TRUE),nr),
    Q75= round(quantile(ratioRF, 0.75, na.rm = TRUE),nr),
    count=n())
stat_summ<-as.data.frame(stat_summ)
stat_summ
