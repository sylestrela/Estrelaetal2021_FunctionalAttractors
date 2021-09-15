" R script for Figure 4A in Estrela et al (2021) 
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
setwd(path.dp)
d.esvs<-read.csv('top6_ESVs_fig2B.csv')

df.all$experiment_s<-paste0(df.all$Treatment,'_', df.all$Inoc_size)

##.. SELECT Kp, Km, A, and P at T18 
esvAlc<-as.character(subset(d.esvs,Genus=='Alcaligenes')$OTU)
esvP1<-as.character(subset(d.esvs,Genus=='Pseudomonas')$OTU)
esvKp<-as.character(subset(d.esvs,Genus=='Klebsiella')$OTU)
esvKm<-as.character(subset(d.esvs,Genus=='Klebsiella_m')$OTU)

dat.kkap<-subset(df.all, OTU %in% c(esvKp, esvKm, esvAlc, esvP1))
dat.kkap$esv_id[dat.kkap$OTU==esvAlc]<-'A'
dat.kkap$esv_id[dat.kkap$OTU==esvP1]<-'P'
dat.kkap$esv_id[dat.kkap$OTU==esvKp]<-'Kp'
dat.kkap$esv_id[dat.kkap$OTU==esvKm]<-'Km'
dat.kkap<-select(dat.kkap, -OTU,-Carbon,-Family,-Treatment,-Inoc_size)

## .. SELECT replicates we have timeseries for
#..No_migration_4
tt<-c("No_migration_4")
p1<-c(2,3, 5,6,14,15,18,20, 26)
p2<-c(31,32,36,44,45,49,54,63, 75, 77)
df.sub<-subset(df.all, subset = Rep_num %in% c(p1,p2) & experiment_s==tt)

#.. calculate number of timeseries
length(c(p1,p2))

df<-df.sub
df$group <- with(df, ifelse(OTU %in% c(esvKp, esvKm, esvAlc, esvP1), 'dom', 'other'))
df[df$OTU ==  esvP1, "group"] <- "Pseudomonas (P)"
df[df$OTU ==  esvAlc, "group"] <- "Alcaligenes (A)"
df[df$OTU ==  esvKp, "group"] <- "Klebsiella (Kp)"
df[df$OTU ==  esvKm, "group"] <- "Klebsiella (Km)"

df$group <- factor(df$group, levels = c("other",'Klebsiella (Km)', 'Klebsiella (Kp)', 'Alcaligenes (A)','Pseudomonas (P)'))

#. ESVs color
colKp<-'steelblue4'
colKm<-'#6993BC'
colPseu<-'#A17BB9'
colAlc<-'#FFB94C'
colOther<-"gray90"

tax<-"Genus"
colgenp<-c( colOther, colKm, colKp,colAlc,  colPseu)
p<-ggplot(df, aes(x=Transfer, y=Abundance, fill=group))
p + geom_bar(color=colOther, stat="identity", position="stack",size=0.1,width=1) + 
  facet_wrap(Rep_num ~., ncol=5)+
  scale_x_discrete(limits=c(1, 18, 18))+
  scale_y_continuous(breaks=seq(0,1,1))+
  scale_fill_manual(values=colgenp, na.value=colOther)+ 
  labs(y= "relative abundance", x="transfer")+
  theme_minimal()+
  theme(
        axis.text.y=element_text(size=10),
        axis.text.x = element_text(size=10, hjust =0.5),
        axis.title=element_text(size=16),
        strip.text.x = element_blank(),
        legend.text=element_text(size=12),
        legend.position='bottom',
        panel.border = element_rect(colour = "gray85", fill=NA, size=1))+
  guides(fill=guide_legend(ncol=2, title=""))

ggsave(path=path.p, paste0("fig4A.pdf"), width= 10, height=5, onefile = TRUE)

## ..Write data table as csv file
fname<-paste("fig4A.csv")
ds<-select(df, Transfer, Rep_num, OTU, Family, group, Abundance)
write.csv(ds,file.path(path.ds, fname),row.names=FALSE)
