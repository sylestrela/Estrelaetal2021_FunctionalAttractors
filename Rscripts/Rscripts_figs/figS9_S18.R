" R script for Figures S9 and S18 in Estrela et al (2021) 
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

# ..Read data table for fba (published models)
d.fba<-read.csv('Estrela_2021_fba_params_Pairs_Acetate.csv')
dat.fba<-select(d.fba, -E,-P,-E_Genus,-P_Species,-R_F)
dat.fba$ratioRF<-dat.fba$D_ace_glu*dat.fba$w_ac_R/dat.fba$w_glc_F
dat.fba$Experiment<-'fba'

# ..Read data table for RF ratio of empirically calibrated model
dat.emp.calc<-read.csv('Estrela_2021_empir_params_Pairs_Acetate.csv')
dat.emp.calc<-select(dat.emp.calc, -X)
dat.emp.calc$Experiment<-'empir_calc'

# ..Read data table for RF ratio of fba (models for C2R4 and C2R2 strains)
d.fba.es<-read.csv('Estrela_2021_Whole_Genome_R_F_Pred.csv')
dat.fba.es<-select(d.fba.es, -E,-P,-E_Genus,-P_Species,-R_F)
dat.fba.es$ratioRF<-dat.fba.es$D_ace_glu*dat.fba.es$w_ac_R/dat.fba.es$w_glc_F
dat.fba.es$Experiment<-'fba-ES'

#.. Plot A) D_ace_glu, w_f_glu and w_f_ace for empirically calibrated model
dp11<-select(dat.emp.calc, SangerID_f, D_ace_glu_corr)
dp11$SangerID<-dp11$SangerID_f
dp11<-unique(dp11)
dp11<-melt(select(dp11, -SangerID_f), id.vars=c('SangerID'))

dp12<-select(dat.emp.calc, SangerID_f, w_f_glu)
dp12$SangerID<-dp12$SangerID_f
dp12<-unique(dp12)
dp12<-melt(select(dp12, -SangerID_f), id.vars=c('SangerID'))

dp13<-select(dat.emp.calc, SangerID_r, w_r_ace)
dp13$SangerID<-dp13$SangerID_r
dp13<-unique(dp13)
dp13<-melt(select(dp13, -SangerID_r), id.vars=c('SangerID'))

dp<-rbind(dp11,dp12,dp13)

colFp<-"deepskyblue4"
colFp2<-"deepskyblue2"
colRp<-"darkorchid3"

#..need to change the y axis for each facet for better visualization
dp <- data.table(dp)
dp[variable == "D_ace_glu_corr",y_min := 0]
dp[variable == "D_ace_glu_corr",y_max := 1]
dp[variable == "w_f_glu",y_min := 0]
dp[variable == "w_f_glu",y_max := 0.025]
dp[variable == "w_r_ace",y_min := 0]
dp[variable == "w_r_ace",y_max := 0.025]

set.seed(10)
p1<-ggplot(dp, aes(x=factor(1), y=value))+
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, aes(colour=variable), size=2, alpha=0.7)+
  facet_wrap(variable~., scales='free_y')+
  theme_classic()+
  scale_shape_manual(values = c(1, 2, 4)) +
  scale_colour_manual(values = c(colFp2,colFp,colRp), labels = c("D[ace,glu]",  "W[glu, F]",  "W[ace,R]"))+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=16),
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    axis.title.y =element_text(size=16),
    axis.title.x = element_blank(), 
    axis.ticks.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "gray10"),
    legend.title = element_blank(),
    legend.text = element_text(size=16))+
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))
p1

length(unique(dp$SangerID))

#.. write data figure as csv file
ds<-select(dp, -y_min,-y_max)
fname<-paste("figS9A.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)

#.. Plot B) D_ace_glu for empir vs fba model
dp21<-select(dat.emp.calc,D_ace_glu_corr,Experiment,SangerID_f)
names(dp21)[names(dp21)=="D_ace_glu_corr"]<-"D_ace_glu"
names(dp21)[names(dp21)=="SangerID_f"]<-"ID_f"
dp21$type<-'emp'
  
dp22<-select(d.fba, E, D_ace_glu)
names(dp22)[names(dp22)=="E"]<-"ID_f"
dp22$Experiment<-'fba'
dp22$type<-'published'

dp23<-select(d.fba.es, E, D_ace_glu)
names(dp23)[names(dp23)=="E"]<-"ID_f"
dp23$Experiment<-'fba'
dp23$type<-'es'

dp2<-rbind(dp21,dp22)
dp2<-unique(dp2)

set.seed(1)
p2<-ggplot(dp2, aes(x=Experiment, y=D_ace_glu, colour=type))+
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size=2, alpha=0.6)+
   theme_classic()+
  scale_colour_manual(values = c('gray20', 'black', 'gray60'))+
  ylim(0,1.2)+
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
  theme(legend.position = "none")+
  scale_x_discrete( 
    labels=c("Empirical 
calibration","FBA 
simulations"))+
  labs(y='D[ace,glu]')
p2

length(unique(subset(dp2, Experiment=='fba')$ID_f))
length(unique(subset(dp2, Experiment=='empir_calc')$ID_f))

#.. write data figure as csv file
ds<-dp2
fname<-paste("figS9B.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)

# Plot C) w_f_ace/w_f_glu for empir vs fba model 
dp31<-select(dat.emp.calc,w_f_glu, w_r_ace, Experiment,SangerID_f,SangerID_r)
dp31$w_ratio<-dp31$w_r_ace/dp31$w_f_glu
dp31$rf_pair<-paste0(dp31$SangerID_r, '_',dp31$SangerID_f)
dp31<-select(dp31,-SangerID_f,-SangerID_r,-w_f_glu,-w_r_ace)
dp31$type<-'emp'
  
dp32<-select(d.fba, -E_Genus,-P_Species,-D_ace_glu,-R_F)
dp32$w_ratio<-dp32$w_ac_R/dp32$w_glc_F
dp32$rf_pair<-paste0(dp32$P,'_', dp32$E)
dp32$Experiment<-'fba'
dp32<-select(dp32,-E,-P,-w_ac_R,-w_glc_F)
dp32$type<-'published'

dp33<-select(d.fba.es, -E_Genus,-P_Species,-D_ace_glu,-R_F)
dp33$w_ratio<-dp33$w_ac_R/dp33$w_glc_F
dp33$rf_pair<-paste0(dp33$P,'_', dp33$E)
dp33$Experiment<-'fba'
dp33<-select(dp33,-E,-P,-w_ac_R,-w_glc_F)
dp33$type<-'es'

dp3<-rbind(dp31,dp32)
dp3<-unique(dp3)

set.seed(1)
p3<-ggplot(dp3, aes(x=Experiment, y=w_ratio, colour=type))+
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size=2, alpha=0.6)+
  theme_classic()+
  scale_colour_manual(values = c('gray20', 'black', 'gray60'))+
  ylim(-0.01,1.25)+
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
  theme(legend.position = "none")+
  scale_x_discrete( 
    labels=c("Empirical 
calibration","FBA 
simulations"))+
  labs(y='w[ace,R] / w[glu,F]')
p3

length(unique(subset(dp3, Experiment=='fba')$rf_pair))
length(unique(subset(dp3, Experiment=='empir_calc')$rf_pair))

#.. write data figure as csv file
ds<-dp3
fname<-paste("figS9C.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)

ggarrange(
  p1,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(p2, p3, ncol = 2, labels = c("B", "C")), 
  nrow = 2, 
  labels = "A"       # Label of the line plot
) 

ggsave(path=path.p, paste0("figS9.pdf"), width= 9, height=7, onefile = TRUE)

#  ---------------------------------------------------------------------------
#.. plot FBA data only (published + isolates collection genomes)
#.. Plot D_ace_glu, WgluF, waceR for published and isolates collection models

#..get w_aceR / w_gluF for published genomes
dp1<-select(d.fba, -E_Genus,-P_Species,-D_ace_glu,-R_F)
dp1$w_ratio<-dp1$w_ac_R/dp1$w_glc_F
dp1$rf_pair<-paste0(dp1$P,'_', dp1$E)
dp1$Experiment<-'fba'
dp1<-select(dp1,-E,-P,-w_ac_R,-w_glc_F)
dp1$type<-'published'
#..get w_aceR / w_gluF for isolates collection
dp12<-select(d.fba.es, -E_Genus,-P_Species,-D_ace_glu,-R_F)
dp12$w_ratio<-dp12$w_ac_R/dp12$w_glc_F
dp12$rf_pair<-paste0(dp12$P,'_', dp12$E)
dp12$Experiment<-'fba'
dp12<-select(dp12,-E,-P,-w_ac_R,-w_glc_F)
dp12$type<-'isolates collection'

dp1<-rbind(dp1,dp12)
dp1<-unique(dp1)
dp1$variable<-'w[ace,R] / w[glu,F]'
#names(dp1)[names(dp1)=='w_ratio']<-'value'
  
#..get Dace_glu for published genomes
dp2<-select(d.fba, E, D_ace_glu)
names(dp2)[names(dp2)=="E"]<-"ID_f"
dp2$Experiment<-'fba'
dp2$type<-'published'
#..get Dace_glu for isolates collection
dp21<-select(d.fba.es, E, D_ace_glu)
names(dp21)[names(dp21)=="E"]<-"ID_f"
dp21$Experiment<-'fba'
dp21$type<-'isolates collection'
dp2<-rbind(dp2,dp21)
dp2<-unique(dp2)
dp2$variable<-'D[ace,glu]'
#names(dp2)[names(dp2)=='D_ace_glu']<-'value'

#..get rf_ratio for published genomes
dp3<-select(d.fba, E, P, R_F)
dp3$Experiment<-'fba'
dp3$type<-'published'
#..get rf_ratio for isolates collection
dp31<-select(d.fba.es, E, P, R_F)
dp31$Experiment<-'fba'
dp31$type<-'isolates collection'
dp3<-rbind(dp3,dp31)
dp3<-unique(dp3)
dp3$variable<-'R/F'

set.seed(5)
p1<-ggplot()+
  geom_boxplot(data=dp1, aes(x=variable, y=w_ratio), outlier.shape = NA, colour='gray80')+ 
  geom_jitter(data=dp1, aes(x=variable, y=w_ratio, colour=type), width = 0.2, size=2, alpha=0.7)+
  facet_wrap(variable~., scales='free_y')+
  theme_classic()+
  scale_shape_manual(values = c(1, 2, 4)) +
  scale_colour_manual(values = c('black', 'gray80'),name='GSM models')+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=16),
    strip.text.x = element_text(size=16),
    strip.background = element_blank(),
    axis.title.y =element_blank(),
    axis.title.x = element_blank(), 
    axis.ticks.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "gray10"),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16))+
   scale_y_continuous(limits = c(NA, 1), breaks = seq(0, 1, by = 0.5))
p1

set.seed(2)
p2<-ggplot()+
  geom_boxplot(data=dp2, aes(x=variable, y=D_ace_glu), outlier.shape = NA, colour='gray80')+ 
  geom_jitter(data=dp2, aes(x=variable, y=D_ace_glu, colour=type), width = 0.2, size=2, alpha=0.7)+
  facet_wrap(variable~., scales='free_y')+
  theme_classic()+
  scale_shape_manual(values = c(1, 2, 4)) +
  scale_colour_manual(values = c('black', 'gray80'),name='GSM models')+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=16),
    strip.text.x = element_text(size=16),
    strip.background = element_blank(),
    axis.title.y =element_blank(),
    axis.title.x = element_blank(), 
    axis.ticks.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "gray10"),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16))+
  scale_y_continuous(limits = c(NA, 1.2), breaks = seq(0, 1.2, by = 0.6))
p2

set.seed(10)
p3<-ggplot()+
  geom_boxplot(data=dp3, aes(x=variable, y=R_F), outlier.shape = NA, colour='gray80')+ 
  geom_jitter(data=dp3, aes(x=variable, y=R_F, colour=type), width = 0.2, size=2, alpha=0.7)+
  facet_wrap(variable~., scales='free_y')+
  theme_classic()+
  scale_shape_manual(values = c(1, 2, 4)) +
  scale_colour_manual(values = c('black', 'gray80'), name='GSM models')+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=16),
    strip.text.x = element_text(size=16),
    strip.background = element_blank(),
    axis.title.y =element_blank(),
    axis.title.x = element_blank(), 
    axis.ticks.x=element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "gray10"),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16))+
  scale_y_log10(breaks=c(0.01, 0.1, 1, 10), limits=c(0.001, 10))
p3

#.. write data figure as csv file
ds<-dp1
fname<-paste("figS18A.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)
ds<-dp2
fname<-paste("figS18B.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)
ds<-dp3
fname<-paste("figS18C.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)

ggarrange(p2,p1,p3, ncol=3, nrow=1, common.legend = TRUE, legend='right')

ggsave(path=path.p, paste0("figS18.pdf"), width= 10.5, height=4, onefile = TRUE)


