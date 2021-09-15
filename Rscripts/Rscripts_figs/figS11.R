" R script for Figure S11 in Estrela et al (2021) 
Functional attractors in microbial community assembly."

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
dat<-read.csv('Estrela_2021_alternate_CSources_datatable_relabund.csv')

# .. Colour scheme for barplots
# Family colors 
colEntbf<- "deepskyblue3"
colPseuf<-  "darkorchid2"
colOther<-"gray85"

# ..SELECT minimum abundance to create abundant and rare groups
min.ab<-0.01
dp.abund<-subset(dat, Relative_Abundance >min.ab )
dp.rare<-subset(dat, Relative_Abundance <=min.ab )
dp.abund$family_type<-'abundant'
dp.rare$family_type<-'rare'

dp<-rbind(dp.abund, dp.rare)

dp$cs.ord <- factor(dp$Carbon_Source, levels = c('glucose', 'fructose', 'cellobiose', 'ribose', 'citrate', 'glutamine'))

#.. Assign F/NF profile based on literature search (see SI table for taxa assignment)
dp$f_r_type <- with(dp, 
                      ifelse(Family == "Enterobacteriaceae" | Family == "Enterococcaceae" | Family=="Aeromonadaceae" | Family=="Lachnospiraceae",
                             "F", ifelse(Family == "Alcaligenaceae" | Family == "Pseudomonadaceae" | Family=="Moraxellaceae" | Family=='Xanthomonadaceae' | Family=="Comamonadaceae" | Family == "Oxalobacteraceae",
                                         "R", 'UN')))
head(dp)

colFp<-"deepskyblue4"
colRp<-"darkorchid2"

xl<-"Replicate"
yl<-"Relative abundance"

p1<-ggplot(dp, aes(x=Replicate,y=Relative_Abundance ,fill=f_r_type))+
  geom_bar(stat="identity", position="stack", color="black") + 
  facet_grid(.~cs.ord)+
  scale_x_discrete(limits=c(1:4))+
  theme_minimal()+
  theme(axis.text.x = element_text(size=8, angle = 0, hjust = 0.5))+
  labs(x=xl, y=yl)+
  theme_minimal()+
  scale_fill_manual(values=c(colFp, colRp, 'gray'), name='')+
  scale_y_continuous(breaks= c(0, 1)) +
  labs(x='Replicate', y='   Relative 
  abundance')+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(size=12),
        legend.title = element_text(size=12),
        axis.title.y =element_text(size=16),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size=12),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
p1

#.. write data figure as csv file
ds<-select(dp, Carbon_Source, Sample_ID, Replicate, f_r_type, Sequence,Relative_Abundance)
fname<-paste("figS11A.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)

##.. calculate R/F ratio
#. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fr<-data.table(dp)[, list(f_r_sum=sum(Relative_Abundance)), 
                          by=list(f_r_type, cs.ord, Replicate)]

#. note: if there are cases with no F or no NF, we need to reshape the dataframe so their rel. ab is 0 or 1.
#. ratio of Respirer/Fermenter
dsum_fr2 <- dsum_fr %>% 
  spread(f_r_type, f_r_sum) %>%
  group_by(cs.ord, Replicate) %>% 
  mutate(ratioRF = R/F)
dsum_fr2<-as.data.frame(dsum_fr2)

dsum_fr2$R[dsum_fr2$F==1]<-0
dsum_fr2$ratioRF[dsum_fr2$F==1]<-0

median(dsum_fr2$ratioRF)
mean(dsum_fr2$ratioRF)

# .. Calculate median and IQR for each CS
ds<-dsum_fr2
ds[is.na(ds)] <- 0

nr<-3
stat_summ<- ds %>%
  dplyr::group_by(cs.ord) %>%
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

#.. calculate 95% percentil for R/F ratio of communities in Goldford et al (2018). Add as a line the R/F ratio of Emergent Simplicity paper (fig1A): median R/F=0.29, prob_low=0.025, prob_high=0.975
d.es<-read.csv(paste0(path.dp, '/Goldford_2018_glu_ratioRF.csv'))

es_rf_med<-median(d.es$ratioRF)
es_rf_cil<-quantile(d.es$ratioRF, probs = .025)
es_rf_cih<-quantile(d.es$ratioRF, probs = .975)
hist(d.es$ratioRF)

xl<-""
yl<-"R/F"
p2<- ggplot(ds, aes(x=factor(cs.ord), ratioRF))+
  geom_boxplot(outlier.shape = NA) + 
  annotate("rect", ymin = es_rf_cil, ymax = es_rf_cih, xmin = 0, xmax = Inf,
           alpha = .3, fill='gray85')+
  geom_jitter(data=ds,width = 0.2, size=2.5, alpha=0.8)+
  geom_hline(yintercept=es_rf_med, linetype="dashed", color = "gray")+
  theme_classic()+
  scale_x_discrete(expand = c(0.1, 0.1))+
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y =element_text(size=18),
    axis.title.x = element_text(size = 14), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linetype='dashed', size=0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14))+
  labs(x=xl, y=yl)+
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1,10,100),limits = c(0.001, 400))+
  theme(plot.margin=unit(c(0.2,2.5,0.2,0.2),"cm"))
p2

#.. write data figure as csv file
names(ds)[names(ds)=='cs.ord']<-'Carbon_Source'
fname<-paste("figS11B.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)

ggarrange(p1, p2, ncol=1, nrow=2, labels = c("A","B"))

ggsave(path=path.p, "figS11.pdf", width=8, height=5, onefile = TRUE)

