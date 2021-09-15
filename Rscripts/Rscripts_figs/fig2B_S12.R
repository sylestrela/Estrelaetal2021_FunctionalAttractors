" R script for Figures 2B and S12 (+stats) in Estrela et al (2021) 
Functional attractors in microbial community assembly."

## ---------------------------------------------------------------------------
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

#.. SELECT data
trn<-"No_migration"
tp<-"18"
inocs<-'4'
df.all.s<- subset(df.all, Treatment== trn & Transfer==tp & Inoc_size==inocs)
head(df.all.s)
df.sub<-select(df.all.s, -Carbon, -Transfer, -Treatment,-Kingdom,-Phylum,-Class,-Order)
head(df.sub)

#------------------------------------------------------------------------------
# .. Calculate the number of communities with R+F and with F only above 0.01
nAd<-subset(df.sub, Family== 'Alcaligenaceae' & Abundance>0.01)
nPd<-subset(df.sub, Family== 'Pseudomonadaceae' & Abundance>0.01)
nAdPd<-rbind(nAd,nPd)
nFd<-subset(df.sub, !(Rep_num %in% nAdPd$Rep_num))
length(unique(nAd$Rep_num))
length(unique(nPd$Rep_num))
length(unique(nFd$Rep_num))

# .. Calculate the number of communities with R+F and with F only
nAd<-subset(df.sub, Family== 'Alcaligenaceae' & Abundance>0)
nPd<-subset(df.sub, Family== 'Pseudomonadaceae'& Abundance>0)
nAdPd<-rbind(nAd,nPd)
nFd<-subset(df.sub, !(Rep_num %in% nAdPd$Rep_num))
length(unique(nAd$Rep_num))
length(unique(nPd$Rep_num))
length(unique(nFd$Rep_num))

nE<-subset(df.sub, Family== 'Enterococcaceae'& Abundance>0)
length(unique(nE$Rep_num))
nC<-subset(df.sub, Family== 'Comamonadaceae'& Abundance>0)
length(unique(nC$Rep_num))

# .. Calculate the number of communities where A and P are both present at T18 even if rare
nAP<-subset(nPd, (Rep_num %in% nAd$Rep_num))
length(unique(nAP$Rep_num))

#------------------------------------------------------------------------------
#.. Figure 2B by family group clustering

#.. SELECT cluster data to plot and create groups
min.ab<-0.01
dps<-subset(df.sub, Abundance>min.ab)

# Group 1: K+ only. no A, no P
nAd<-subset(dps, Family== 'Alcaligenaceae' & Abundance>min.ab)
nPd<-subset(dps, Family== 'Pseudomonadaceae' & Abundance>min.ab)
nAdPd<-rbind(nAd,nPd)
n1<-subset(dps, !(Rep_num %in% nAdPd$Rep_num))
n1<-n1[order(n1$Genus,n1$Abundance),]
n1<-unique(n1$Rep_num)
n1
length(n1)
# number of A dominated communities with Achromobacter
nAch<-subset(dps, (Rep_num %in% nAd$Rep_num) & Genus=='Achromobacter')
nAchl<-length(unique(nAch$Rep_num))
nAchl
# number of A dominated communities with Delftia
nDel<-subset(dps, (Rep_num %in% nAd$Rep_num) & Genus=='Delftia')
nDell<-length(unique(nDel$Rep_num))
nDell
# number of A dominated communities with Delftia & Achromobacter
nDelAch<-subset(dps, (Rep_num %in% nAch$Rep_num) & Genus %in% c('Delftia'))
nDelAchl<-length(unique(nDelAch$Rep_num))
nDelAchl

# Group 2: P and K+
nPd<-subset(df.sub, Family== 'Pseudomonadaceae' & Abundance>min.ab)
nPd<-nPd[order(nPd$Family,nPd$Abundance),]
n2<-unique(nPd$Rep_num)
n2
length(n2)
# number of PK+ communities with Enterococcus
nE<-subset(dps, (Rep_num %in% nPd$Rep_num) & Genus=='Enterococcus')
nE<-length(unique(nE$Rep_num))

# Group 3 and others: contains A, note: length(n1+n2+n3=92)
dps2<-subset(dps, (Rep_num %in% nAd$Rep_num))
n3<-unique(dps2$Rep_num)
n3
length(n3)

# Group 4:  K+ and A only (subtract Km and Delftia)
df<-subset(dps2, Genus %in% c("Klebsiella_m", 'Delftia') & Abundance>min.ab)
nKA<-subset(dps2, !(Rep_num %in% df$Rep_num))
nKA<-nKA[order(nKA$Genus,nKA$Abundance),]
n4<-unique(nKA$Rep_num)
n4
length(n4)

dps3<-subset(dps2, !(Rep_num %in% nKA$Rep_num))
n32<-unique(dps3$Rep_num)
n32
length(n32)

# Group 5:  K+ + A +Km only (subtract Delftia )
df<-subset(dps3, Genus %in% c("Delftia") & Abundance>min.ab)
nKKA<-subset(dps3, !(Rep_num %in% df$Rep_num))
nKKA$Genus <- factor(nKKA$Genus, levels = c("Klebsiella_m", "Klebsiella",  "Alcaligenes" ))
nKKA<-nKKA[order(nKKA$Genus,-nKKA$Abundance),]
n5<-unique(nKKA$Rep_num)
n5
length(n5)
# remove group 5
dps4<-subset(dps3, !(Rep_num %in% nKKA$Rep_num))
n33<-unique(dps4$Rep_num)
n33
length(n33)

# Group 6:  K+ + A +D  (subtract Km )
df<-subset(dps4, Genus %in% c("Klebsiella_m") & Abundance>min.ab)
nKAD<-subset(dps4, !(Rep_num %in% df$Rep_num))
nKAD<-nKAD[order(nKAD$Genus,nKAD$Abundance),]

n6<-unique(nKAD$Rep_num)
n6
length(n6)

# remove group 5
dps5<-subset(dps4, !(Rep_num %in% nKAD$Rep_num))
dps5$Genus <- factor(dps5$Genus, levels = c("Klebsiella_m", "Klebsiella",  "Alcaligenes", 'Delftia', '<NA>' ))
dps5<-dps5[order(dps5$Genus,dps5$Abundance),]
n34<-unique(dps5$Rep_num)
n34
length(n34)

# Group 7: 
n7<-n34
length(n7)

gall<-c(n1,n2,n4,n6,n7,n5)
length(gall)

#.. plot
dfp<-dps

#. Genus colors 
colKp<-'steelblue4'
colKm<-'#6993BC'
colPseu<-'#A17BB9'
colAlc<-'#FFB94C'
colDel<-'indianred'
colEntc<-'#8EB8B6'
colAch<-"orange4"
colgenp2<-c(colKm, '#A2B5CD', colKp,'#3A5874',colPseu, 'plum4', colAch, colAlc, '#8EB8B6', colDel)

p<-ggplot(dfp, aes(x=factor(Rep_num), y=Abundance, fill=OTU))+
  geom_bar(color="black", stat="identity", position="stack") +
  scale_x_discrete(limits=as.character(gall)) +
  scale_fill_manual(values=colgenp2,
                    name  ="",
                    labels=c("Klebsiella (Km)", "Enterobacteriaceae_1", 'Klebsiella (Kp)', "Enterobacteriaceae_2",
                             "Pseudomonas (P)", "Pseudomonas", "Achromobacter", "Alcaligenes (A)",
                             "Enterococcus", "Delftia (D)"))+ 
  labs(y= "relative abundance", x="replicate community")+
  scale_y_continuous(breaks=seq(0,1,1))+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16),
        strip.text.y = element_text(size = 14),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  theme(legend.text = element_text(face = "italic",size=12))
p 

ggsave(path=path.p, paste0("fig2B.pdf"), width= 15, height=3.2, onefile = TRUE)

dat.fig2B<-select(dfp,Rep_num,Genus, OTU, Abundance)
## ..Write data table as csv file
fname<-paste("fig2B_S12.csv")
write.csv(dat.fig2B,file.path(path.ds, fname),row.names=FALSE)

#------------------------------------------------------------------------------
##.. Figure S12- clustered barplots using R cluster's algorithm.
#.. distance method 1: euclidean
#. make cluster based on rel abundance of ESVs
cl.dat <- select(dfp, -Family,-Genus)
cl.dat2 <-acast(cl.dat, Rep_num~OTU, value.var="Abundance")
cl.dat2[is.na(cl.dat2)] <- 0

# .. For dendrogram clustering :
#. SELECT distance method  
methd<- "euclidean"

#. SELECT clustering method
methc<- "ward.D" 

#. make a dendrogram to extract clusters
d = dist(cl.dat2, method=methd)  
#. cluster the distance tree  
hc = hclust(d, method=methc)
#. plot the hclust object
plot(hc)
#. return labels in dendrogram order, not in numerical order
labels(hc)
hc$labels[hc$order]
#. cut the row tree by cluster number
mycl <- cutree(hc, k=6)
rect.hclust(hc, k = 6, border = 2:5)
#. print the obtained cluster numbers in tree order
leafl<-mycl[hc$order]
ml<-melt(leafl)
rownames(ml)

##.. plot 
dfp<-dps
colgenpS10<-c(colKm, blues9[6], colKp, blues9[8], colPseu, 'plum4', colAch, colAlc, colEntc,  colDel)
p1<-ggplot(dfp, aes(x=factor(Rep_num), y=Abundance, fill=OTU))+
  geom_bar(color="black", stat="identity", position="stack") + 
  scale_x_discrete(limits=factor(rownames(ml))) +
  scale_fill_manual(values=colgenpS10,
                    name  ="",
                    labels=c("Klebsiella (Km)", "Enterobacteriaceae_1", 'Klebsiella (Kp)', "Enterobacteriaceae_2",
                             "Pseudomonas (P)", "Pseudomonas", "Achromobacter", "Alcaligenes (A)",
                             "Enterococcus", "Delftia (D)"))+ 
  labs(y= "relative abundance", x="community",
       caption= paste0('Dist method: ', methd, ". Hclust method: ", methc))+
  scale_y_continuous(breaks=seq(0,1,1))+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16),
        strip.text.y = element_text(size = 14),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(face = "italic",size=12))
p1 

#.. distance method 2: binary
# ..make cluster based on rel abundance of ESVs
cl.dat <- select(dfp,-Family,-Genus)
cl.dat2 <-acast(cl.dat, Rep_num~OTU, value.var="Abundance")
cl.dat2[is.na(cl.dat2)] <- 0

# .. For dendrogram clustering :
#. SELECT distance method  
methd<- "binary"

#. SELECT clustering method
methc<- "ward.D" 

#.. make a dendrogram to extract clusters
d = dist(cl.dat2, method=methd )  
#. cluster the distance tree  
hc = hclust(d, method=methc )
#. plot the hclust object
plot(hc)
#. return labels in dendrogram order, not in numerical order
labels(hc)
hc$labels[hc$order]
#. cut the row tree by cluster number
mycl <- cutree(hc, k=6)
rect.hclust(hc, k = 6, border = 2:5)
#. print the obtained cluster numbers in tree order
leafl<-mycl[hc$order]

ml<-melt(leafl)
rownames(ml)

#. plot
dfp<-dps
p2<-ggplot(dfp, aes(x=factor(Rep_num), y=Abundance, fill=OTU))+
  geom_bar(color="black", stat="identity", position="stack") + 
  scale_x_discrete(limits=factor(rownames(ml))) +
  scale_fill_manual(values=colgenp2,
                    name  ="",
                    labels=c("Klebsiella (Km)", "Enterobacteriaceae_1", 'Klebsiella (Kp)', "Enterobacteriaceae_2",
                             "Pseudomonas (P)", "Pseudomonas", "Achromobacter", "Alcaligenes (A)",
                             "Enterococcus", "Delftia (D)"))+ 
  labs(y= "relative abundance", x="community",
       caption= paste0('Dist method: ', methd, ". Hclust method: ", methc))+
  scale_y_continuous(breaks=seq(0,1,1))+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=16),
        strip.text.y = element_text(size = 14),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(face = "italic",size=12))
p2 

ggarrange(p1,p2,
          nrow=2,ncol=1, common.legend=TRUE,legend='bottom')

ggsave(path=path.p, paste0("figS12.pdf"), width= 11, height=6)

#  ---------------------------------------------------------------------------
##.. Some useful calculations and stats used in the paper
#  ---------------------------------------------------------------------------

#.. SELECT min rel abundance
min.ab<-0.0
dfp<-subset(df.sub, Abundance >min.ab)

#. calculate mean of each family across replicates
dfp2<-dfp
d.sum<-data.table(dfp2)[, list(abund_sum=sum(Abundance)), 
                 by=list(Family, Rep_num)]
d.sum
d.sum2<-data.table(dfp2)[, list(mean_fam=mean(Abundance)), 
                        by=list(Family)]
d.sum2
d.mean.fam<-data.table(d.sum)[, list(dfp_mean=mean(abund_sum)), 
                         by=list(Family)]
d.mean.fam

#. calculate mean of A across replicates where A is dominant
d.alcf<-d.sum[d.sum$Family == 'Alcaligenaceae' & d.sum$abund_sum>0.05,]
data.table(d.alcf)[, list(dfp_mean=mean(abund_sum)), 
                 by=list(Family)]

#. calculate mean of P across replicates where P is dominant
d.pseuf<-d.sum[d.sum$Family == 'Pseudomonadaceae' & d.sum$abund_sum>0.05,]
data.table(d.pseuf)[, list(dfp_mean=mean(abund_sum)), 
                   by=list(Family)]

#. calculate mean of each GENUS across replicates
dfp2<-dfp
d.sum<-data.table(dfp2)[, list(abund_sum=sum(Abundance)), 
                        by=list(Genus, Rep_num)]

##.. calculate ratio of R/F and of Pseu/Ent or Alc/Ent for communities with both R+F
df.sub<-select(dfp, -OTU, -Genus,-Inoc_size)
#.create new column with F or R assignment
df.sub$f_r_type <- with(df.sub, 
                        ifelse(Family == "Enterobacteriaceae" | Family == "Enterococcaceae" | Family=="Aeromonadaceae" | Family=="Lachnospiraceae",
                               "F", ifelse(Family == "Alcaligenaceae" | Family == "Pseudomonadaceae" | Family=="Moraxellaceae" | Family=='Xanthomonadaceae' | Family=="Comamonadaceae",
                               "R", 'other')))
head(df.sub)
#. SELECT ONLY TAXA THAT ARE R OR F
df.sub<-subset(df.sub, f_r_type!='other')

#. Calculate the number of communities with R+F and with F only
nAd<-subset(df.sub, Family== 'Alcaligenaceae' & Abundance>0.01)
nPd<-subset(df.sub, Family== 'Pseudomonadaceae' & Abundance>0.01)
nAdPd<-rbind(nAd,nPd)
nFd<-subset(df.sub, !(Rep_num %in% nAdPd$Rep_num))
nRFd<-subset(df.sub, Rep_num %in% nAdPd$Rep_num)
nAcd<-subset(df.sub, Rep_num %in% nAd$Rep_num)
nPcd<-subset(df.sub, Rep_num %in% nPd$Rep_num)

length(unique(nAd$Rep_num))
length(unique(nPd$Rep_num))
length(unique(nFd$Rep_num))
length(unique(nRFd$Rep_num))

#. select communities that coexist only
dat.fr<-nRFd

#.. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fr<-data.table(dat.fr)[, list(f_r_sum=sum(Abundance)), 
                            by=list(f_r_type, Rep_num)]
dat.RF<-as.data.frame(spread(dsum_fr,f_r_type, f_r_sum))
dat.RF[is.na(dat.RF)] <- 0

#. ratio of Respirator/Fermenter
dat.RF$ratioRF<-(dat.RF$R/dat.RF$F)
median(dat.RF$ratioRF)
mean(dat.RF$ratioRF)

IQR_summary<- dat.RF %>%
  dplyr::summarise(
    Min = round(min(ratioRF),3),
    Max = round(max(ratioRF),3),
    mean = round(mean(ratioRF),3),
    std = sd(ratioRF),
    Median = median(ratioRF),
    IQRange = IQR(ratioRF),
    Q1=quantile(ratioRF, 0.25),
    Q3=quantile(ratioRF, 0.75))
IQR_summary

#. SELECT communities with F only
dat.fr<-nFd
#. sum the relative abundance of F or R for each sample (i.e. each rep x treatment x transfer x...)
dsum_fr<-data.table(dat.fr)[, list(f_r_sum=sum(Abundance)), 
                            by=list(f_r_type, Rep_num)]
dat.RF<-as.data.frame(spread(dsum_fr,f_r_type, f_r_sum))
dat.RF[is.na(dat.RF)] <- 0
#. ratio of Respirator/Fermenter
dat.RF$ratioRF<-(dat.RF$R/dat.RF$F)
median(dat.RF$ratioRF)
mean(dat.RF$ratioRF)

IQR_summary<- dat.RF %>%
  dplyr::summarise(
    Min = round(min(ratioRF),3),
    Max = round(max(ratioRF),3),
    mean = round(mean(ratioRF),3),
    std = sd(ratioRF),
    Median = median(ratioRF),
    IQRange = IQR(ratioRF),
    Q1=quantile(ratioRF, 0.25),
    Q3=quantile(ratioRF, 0.75))
IQR_summary
