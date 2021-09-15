" R script for Figure S19 in Estrela et al (2021) 
Functional attractors in microbial community assembly."

#------------------------------------------------------------------------------
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

#.. Read CFU count of R and F for subset of ES communities platted directly 
setwd(path.d)
dat.cfu<-read.csv('Estrela_2021_cfu_ES.csv')
#.. Read R/F ratio for Emergent Simplicity 16S data
setwd(path.dp)
dat.16s<-read.csv('Goldford_2018_glu_ratioRF.csv')

##.. Bootstrapping error
dat<-dat.cfu
s.size<-c()
bs.ratioRF = data.frame()
big_data=list()
datalist = list()
for (j in 1:length(dat$Community)) {
ds<- c(rep("R", dat$respirator[j] ), rep("F", dat$fermenter[j]))
s.size[j]<-length(ds)

set.seed(3)
bs<- c()
bs.rep<-1000
for (i in 1:bs.rep) {
  bs <- sample(ds, size=s.size, replace=TRUE)
  nF<-length(grep("F", bs))
  nR<-length(grep("R", bs))
  RFratio=nR/nF
  datalist[[i]] <- data.frame(RFratio)
}
head(bs)
RFratio
big_data = do.call(cbind, datalist)
bs.ratioRF<-rbind(bs.ratioRF, big_data)
}

bs.ratioRF
dim(bs.ratioRF)
d.bs.ratioRF<- as.data.frame(bs.ratioRF)
d.bs.mean<-apply(d.bs.ratioRF,1,mean)  # applies function 'mean' to 1st dimension (row)
d.bs.sd<-apply(d.bs.ratioRF,1,sd)

f.sem = function(x){
  sd(x)/sqrt(length(x))
}

d.bs.sem<-apply(d.bs.ratioRF,1,f.sem)

d.bs.all<-as.data.frame(cbind(d.bs.mean, d.bs.sd,d.bs.sem))
d.bs.all<-cbind(dat$Community,d.bs.all)
  
#.. Merge cfu data with 16S data
dat.cfu<-d.bs.all
#..rename columns
colnames(dat.cfu)[colnames(dat.cfu) == 'dat$Community'] <- 'Community'
colnames(dat.cfu)[colnames(dat.cfu) == 'd.bs.mean'] <- 'ratioRF_cfu_mean'
colnames(dat.cfu)[colnames(dat.cfu) == 'd.bs.sd'] <- 'ratioRF_cfu_sd'
colnames(dat.cfu)[colnames(dat.cfu) == 'd.bs.sem'] <- 'ratioRF_cfu_sem'

colnames(dat.16s)[colnames(dat.16s) == 'InocRep'] <- 'Community'

d.16s.s<-dat.16s[(dat.16s$Community %in% dat.cfu$Community),]
d.16s.s<-select(d.16s.s, -F, -R)
colnames(d.16s.s)[colnames(d.16s.s) == 'ratioRF'] <- 'ratioRF_16s'

length(unique(d.16s.s$Community))

dat.16s<-d.16s.s
dat.RF.m<-merge(dat.cfu,dat.16s)
dp<-dat.RF.m
dp

##.. plot CFU vs 16S with error bars (+-sem).
# Calculate slope and intercept of line of best fit. 
# NOTE: this linear correlation value is needed for the correction used in Fig S25
dat<-dp
fit.rf.linear <- lm(ratioRF_cfu_mean~ratioRF_16s, dat)
summary(fit.rf.linear)
coef(fit.rf.linear)

ptheme<-  theme(
  axis.text.y = element_text(size = 12, colour="black"),
  axis.text.x = element_text(size=12, colour="black", angle = 0, hjust = 0.5),
  axis.title.y =element_text(size=16),
  axis.title.x = element_text(size = 16),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_line(linetype='dashed', size=0.2),
  panel.grid.major.x = element_line(linetype='dashed', size=0.2),
  panel.grid.minor.y = element_blank())+
  theme(legend.title = element_blank())

colp<-'darkorchid4'
p<-ggplot(fit.rf.linear$model, aes_string(x = names(fit.rf.linear$model)[2], y = names(fit.rf.linear$model)[1])) +
  geom_point() +
  geom_errorbar(data=dat, aes(ymin=ratioRF_cfu_mean-ratioRF_cfu_sem , 
                              ymax=ratioRF_cfu_mean+ratioRF_cfu_sem), width=.03)+
  stat_smooth(method = "lm", col = colp, fill=colp, alpha=0.1, se = TRUE) +
  theme_minimal()+
  labs(x="R/F (16S)", y="R/F (CFU)")+
  ptheme+
  scale_x_log10(breaks = c(0.1, 0.1, 1),limits = c(0.09, 1.1))+
  scale_y_log10(breaks = c(0.001, 0.01, 0.1,1),limits = c(0.0009, 1.1))
p

r.squared<-round(summary(fit.rf.linear)$r.squared,2)
p.value<-signif(summary(fit.rf.linear)$coef[2,4], 5)
slope<-signif(fit.rf.linear$coef[[2]], 3)
Intercept<-signif(fit.rf.linear$coef[[1]],3 )

p+
  geom_text(aes(0.15, 0.55, label = paste("R2 = ", r.squared, "\n",
                                        "p =", round(p.value,3))), colour=colp, size=6)

ggsave(path=path.p, paste0("figS19.pdf"), width= 6, height=4.5, onefile = TRUE)

dat.figS19<-select(dat, Community, ratioRF_cfu_mean, ratioRF_cfu_sem, ratioRF_16s)

## ..Write data table as csv file
fname<-paste("figS19.csv")
write.csv(dat.figS19,file.path(path.ds, fname), row.names=FALSE)
