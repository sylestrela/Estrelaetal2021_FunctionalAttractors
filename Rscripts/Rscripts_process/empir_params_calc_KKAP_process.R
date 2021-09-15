" R script to calculate empirical parameters for the resource-partitioning 
model for the KKAP communities to be used in Fig S15 in Estrela et al (2021) 
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
#------------------------------------------------------------------------------
# ..Calculate D_ace_glu for  Empirical model for K+ and K- at 16h
# ..Read data table for Empirical model

df_experiment = fread(paste0(path.d, '/Estrela_2021_KKP_OAassay.csv'))
df_experiment = df_experiment[!is.na(df_experiment$OD620)]
df_experiment = df_experiment[!is.na(df_experiment$acetate_mM)]
df_experiment = df_experiment[SangerID %in% c(225,235) & time_hours == 16]
df_experiment$Acetate = df_experiment$acetate_mM/11.66
df_experiment$Succinate = df_experiment$succinate_mM/11.66
df_experiment$Lactate = df_experiment$lactate_mM/11.66

# Glu% at t=16h is 0 
df_experiment$Acetate_corr<-df_experiment$Acetate
df_experiment$Succinate_corr<-df_experiment$Succinate
df_experiment$Lactate_corr<-df_experiment$Lactate

df_experiment$D_ace_glu = (df_experiment$Acetate *1)
df_experiment$D_ace_glu_corr = (df_experiment$Acetate_corr *1)
dp<-select(df_experiment, time_hours, SangerID,Family,OD620, D_ace_glu, D_ace_glu_corr)

dp<-melt(dp, id.vars=c('time_hours', 'SangerID','Family','OD620'))
dp<-subset(dp, variable=='D_ace_glu_corr')
dps<-select(dp,SangerID,Family,variable,value)

#------------------------------------------------------------------------------
# ..Calculate W_F_glu for  Empirical model

d.F.glu = fread(paste0(path.d, '/Estrela_2021_KKP_OAassay.csv'))
d.F.glu = d.F.glu[!is.na(d.F.glu$OD620)]
d.F.glu = d.F.glu[!is.na(d.F.glu$acetate_mM)]
d.F.glu<-as.data.frame(d.F.glu)

d.F.glu<-na.omit(select(subset(d.F.glu, SangerID %in% c(225,235)  & time_hours %in% c(0, 16)), time_hours,SangerID, Genus,Family,OD620))

d.F.glu_t0<-subset(d.F.glu, time_hours==0)
d.F.glu_t0$od_t0<-d.F.glu_t0$OD620

d.F.glu_t0<-select(d.F.glu_t0,-OD620,-time_hours)
d.F.glu.2 <-merge(d.F.glu_t0, subset(d.F.glu, time_hours!=0))
d.F.glu.2$F_glu_dt<-d.F.glu.2$OD620-d.F.glu.2$od_t0
d.F.glu.2<-select(d.F.glu.2,-od_t0,-OD620)
names(d.F.glu.2)[names(d.F.glu.2)=='SangerID']<- 'SangerID_f'
names(d.F.glu.2)[names(d.F.glu.2)=='Genus']<- 'Genus_f'
names(d.F.glu.2)[names(d.F.glu.2)=='Family']<- 'Family_f'

f.glu<-d.F.glu.2
f.glu$Glucose_perc<-0 #.. because all glucose is depleted
f.glu$w_f_glu<-f.glu$F_glu_dt/11.66*(0.2-f.glu$Glucose_perc)/0.2

f.glu.s<-select(f.glu, -time_hours,-Glucose_perc,-F_glu_dt,-Genus_f)
names(f.glu.s)[names(f.glu.s)=='SangerID_f']<- 'SangerID'
names(f.glu.s)[names(f.glu.s)=='Family_f']<- 'Family'
f.glu.s<-melt(f.glu.s, id.vars=c('SangerID','Family'))

#------------------------------------------------------------------------------
# ..Calculate W_R_ace for Empirical model
#.. 1) Read OD data of A and P growing on ace for 32h
d.R.oa<-read.csv('Estrela_2021_od_32h_AP_acetate.csv')
names(d.R.oa)[names(d.R.oa)=='od.dt']<-'R_oa_dt'
r.ace.m<-d.R.oa

# Divide by acetate []
r.ace.m$w_r_ace<-r.ace.m$R_oa_dt/r.ace.m$mM

r.ace.s<-select(r.ace.m, strain,w_r_ace, rep)
r.ace.s<-melt(r.ace.s, id.vars=c('strain','rep'))

#.. plot Dace_glu, w_f_glu and w_r_ace 
f.glu.s$strain[f.glu.s$SangerID==225]<-'Kp'
f.glu.s$strain[f.glu.s$SangerID==235]<-'Km'
dps$strain[dps$SangerID==225]<-'Kp'
dps$strain[dps$SangerID==235]<-'Km'

#.. calculate R/F ratio
f.glu.c<-select(f.glu, -Genus_f,-time_hours,-Glucose_perc,-F_glu_dt)
dp.c<-spread(dps,variable,value)
names(dp.c)[names(dp.c)=='SangerID']<- 'SangerID_f'
names(dp.c)[names(dp.c)=='Family']<- 'Family_f'

d.f<-merge(f.glu.c,dp.c)

r.ace.c<-select(r.ace.m,-R_oa_dt)
r.ace.c$Family[r.ace.c$strain=='A']<-'Alcaligenaceae'
r.ace.c$Family[r.ace.c$strain=='P']<-'Pseudomonadaceae'
r.ace.c$SangerID[r.ace.c$strain=='A']<-241
r.ace.c$SangerID[r.ace.c$strain=='P']<-226

names(r.ace.c)[names(r.ace.c)=='SangerID']<- 'SangerID_r'
names(r.ace.c)[names(r.ace.c)=='Family']<- 'Family_r'

dat.empc<-merge(select(d.f,-strain),select(r.ace.c,-strain,-mM,-CS,-dt))
dat.empc$ratioRF<-dat.empc$D_ace_glu_corr*dat.empc$w_r_ace/dat.empc$w_f_glu

dat.empc.kkap<-dat.empc
## ..Write data table as csv file
fname<-paste("Estrela_2021_empir_params_KKAP_acetate.csv")
write.csv(dat.empc.kkap, file.path(path.d, fname), row.names=FALSE)

