" R script for Figure 4B and 4C in Estrela et al (2021) 
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

#.. Read data 
setwd(path.d)

# ..Read data table 
dat.kap<-read.csv("Estrela_2021_cfu_data_invKAP_all.csv")
dp<-dat.kap

#. define order of culture for plots and rename
dp$dila_ord=factor(dp$dila,levels=c("-","-4", "-3", "-2","-1" ,"0"))
dp$dilp_ord=factor(dp$dilp,levels=c("-","-4", "-3", "-2","-1" ,"0"))

##-------------------------------------------------
## Figure 3D. outcome of invasion as density plot
##-------------------------------------------------
#. Colour scheme for plots- ESVs color
colAlc<-"#E69F00"
colPseu<-  "darkorchid4" 
colKp<-blues9[7]
colKm<-blues9[5]

##. get outcome of coexistence based on presence/absence of A and/or P
dp$outcome <- with(dp, ifelse(Alcaligenes >0 & Pseudomonas > 0, "K+P+A", 
                              ifelse(Alcaligenes >0 & Pseudomonas == 0, "K+A",
                                     ifelse(Alcaligenes ==0 & Pseudomonas > 0, "K+P",
                                            ifelse(Alcaligenes ==0 & Pseudomonas == 0, "K",
                                     0)))))

#. plot outcome for frequency of T3 -rep2
col.p <- c("K" = colKp, "K+A" = colAlc, "K+P" = colPseu, "K+P+A" = "gray75")
dp2<-subset(dp, Transfer==3 & repb==2)
ggplot(dp2, aes(x=log10(freqP), y=log10(freqA)))+
  geom_point(aes( colour=outcome),size=8)+
  facet_grid(repb~Transfer)+
  labs(x=expression(phantom(x)*P[0]*phantom(x)*(log[10])),
       y=expression(phantom(x)*A[0]*phantom(x)*(log[10])),
       caption=paste0("Growth in 500uL M9Glu 0.2%. n=2 bio reps."))+
  scale_colour_manual(values = col.p,
                      name= "Outcome")+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y =element_text(size=14),
        axis.title.x = element_text(size = 14),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype='dashed', size=0.2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

#. plot outcome for frequency of T3-rep1
col.p <- c("K" = colKp, "K+A" = colAlc, "K+P" = colPseu, "K+P+A" = "gray75")
dp2<-subset(dp, Transfer==3 & repb==1)
ggplot(dp2, aes(x=log10(freqP), y=log10(freqA)))+ 
  geom_point(aes( colour=outcome),size=8)+
  facet_grid(repb~Transfer)+
  labs(x=expression(phantom(x)*P[0]*phantom(x)*(log[10])),
       y=expression(phantom(x)*A[0]*phantom(x)*(log[10])),
       caption=paste0("Growth in 500uL M9Glu 0.2%. n=2 bio reps."))+
  scale_colour_manual(values = col.p,
                      name= "Outcome")+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.y =element_text(size=14),
        axis.title.x = element_text(size = 14), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linetype='dashed', size=0.2),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

#. get values for separatrix - note: assume that 0 is extremely small (0+constant)
k<-0.00001

#. VALUES for separatrix (combining values of T3 (rep1) and T3(rep2))
#. get values for P axis
ggplot(dp2, aes(x=log10(freqP), y=log10(freqA)))+ 
  geom_point(aes( colour=outcome),size=8)+
  facet_grid(repb~Transfer)+geom_text(aes(label=round(freqA,10)))
s1_a1<-(0.25-0.03226)/2+0.03226
s1_a2<-(0.00033322-0.000033322)/2+0.000033322

s2_a1<-(0.04545- 0.00474)/2+0.00474
s2_a2<-(0.00474-0.000476)/2+0.000476

s3_a1<- (0.00495-0.0004973)/2+0.0004973
s3_a2<- (0.000049749-0)/2+0

s4_a1<-  (0.0049727-0.0004995)/2+0.0004995
s4_a2<- (4.99725*10^(-5)-0)/2+0

s5_a1<- (0.0049749-0.00049973)/2+0.00049973
s5_a2<- (4.9995*10^(-5)-0)/2-0

#. get values for A axis 
ggplot(dp2, aes(x=log10(freqP), y=log10(freqA)))+ 
  geom_point(aes( colour=outcome),size=8)+
  facet_grid(repb~Transfer)+geom_text(aes(label=round(freqP,9)))

s1_b1<-(0.3226- 0.25)/2+0.25
s1_b2<-(0.333322-0.33322)/2+0.33322

s2_b1<-(0.0474-0.0455)/2+0.0455
s2_b2<-(0.0476-0.0474)/2+0.0474

s3_b1<- 0.00496
s3_b2<- 0.0049749

s4_b1<- 0.0004995
s4_b2<- 0.00049973

s5_b1<- 0.000049973
s5_b2<- 0.000049995

#min.v based on min density of A and P for plotting purpose later
min.v<- 1e-05
sa1<-c(s1_a1,s2_a1,s3_a1,s4_a1, s5_a1, min.v)
sa2<-c(s1_a2,s2_a2,s3_a2,s4_a2,  s5_a2, min.v)
sb1<-c(s1_b1,s2_b1,s3_b1,s4_b1,  s5_b1,min.v )
sb2<-c(s1_b2,s2_b2,s3_b2,s4_b2,  s5_b1, min.v)

#. add min values base on A and P min values
line1<-data.frame(sa=sa1, sb=sb1)
line2<-data.frame(sa=sa2, sb=sb2)

# #. quick plot
# plines<-ggplot() + 
#   geom_point(data=line1, aes(x =sb, y = sa, colour = 'black'), size = 2) +
#   geom_line(data=line1, aes(x =sb, y = sa, colour = 'black'), size = 2) +
#   geom_point(data=line2, aes(x =sb, y = sa, colour = 'black'), size = 2) +
#   geom_line(data=line2, aes(x =sb, y = sa, colour = 'black'), size = 2) +
#   scale_x_log10()+
#   scale_y_log10()
# plines

##-------------------------------------------------
## Add 16S abundance data for timeseries of A vs P
##-------------------------------------------------
#.. Read data 
setwd(path.d)
df.all<-read.csv('Estrela_2021_16S_data_table_RelAbund_ALL.csv')
setwd(path.dp)
d.esvs<-read.csv('top6_ESVs_fig2B.csv')

df.all$experiment_s<-paste0(df.all$Treatment,'_', df.all$Inoc_size)
tt<-c("No_migration_4")
df.sub<-subset(df.all, experiment_s==tt)

##.. SELECT A and P at T18 
esvAlc<-unique(as.character(subset(d.esvs,Genus=='Alcaligenes')$OTU))
esvP1<-unique(as.character(subset(d.esvs,Genus=='Pseudomonas')$OTU))
dat.ap<-subset(df.sub, OTU %in% c(esvAlc, esvP1))
dat.ap$esv_id[dat.ap$OTU==esvAlc]<-'A'
dat.ap$esv_id[dat.ap$OTU==esvP1]<-'P'
dat.ap<-select(dat.ap, Abundance, Transfer, Rep_num, Genus)

#.. Cast data from long to wide format 
dfc<-melt(dat.ap, id.vars = c("Transfer", "Rep_num", "Abundance"))
dfc<-select(dfc, -variable)
dfc2<-dcast(dfc, Transfer + Rep_num ~ value, value.var="Abundance")

#. check for min rel abundance of A and P and set NA to min value
min1<-min(dfc2[,"Alcaligenes"], na.rm=T)
min2<-min(dfc2[,"Pseudomonas"], na.rm=T)
minv<-min(min1,min2)
rare.th<-minv
dfc2[is.na(dfc2)] <- rare.th
head(dfc2)

#.. select communities we have timeseries for
p1<-c(2,3, 5,6,14,15,18,20, 26)
p2<-c(31,32,36,44,45,49,54,63, 75, 77)
dfc2.s<-subset(dfc2, subset = Rep_num %in% c(p1,p2))
nrep<-length(unique(dfc2.s$Rep_num))

colnames(line1) <- c("Alcaligenes", "Pseudomonas")
colnames(line2) <- c("Alcaligenes", "Pseudomonas")
dp3<-dfc2.s
rn<-unique(dp3$Rep_num)
df.scat<-select(dp3, Transfer, Alcaligenes,Pseudomonas,Rep_num)
df.scat$type<-"scat"
line1$type<-'sep1'
line2$type<-'sep2'
df.sep<-rbind(line1, line2)

df.sep$Transfer<-1
df.sep$Rep_num<-rn[1]
df2<-df.sep
df2$Rep_num<-rn[2]
df3<-df.sep
df3$Rep_num<-rn[3]
df4<-df.sep
df4$Rep_num<-rn[4]
df5<-df.sep
df5$Rep_num<-rn[5]
df6<-df.sep
df6$Rep_num<-rn[6]
df7<-df.sep
df7$Rep_num<-rn[7]
df8<-df.sep
df8$Rep_num<-rn[8]
df9<-df.sep
df9$Rep_num<-rn[9]
df10<-df.sep
df10$Rep_num<-rn[10]
df11<-df.sep
df11$Rep_num<-rn[11]
df12<-df.sep
df12$Rep_num<-rn[12]
df13<-df.sep
df13$Rep_num<-rn[13]
df14<-df.sep
df14$Rep_num<-rn[14]
df15<-df.sep
df15$Rep_num<-rn[15]
df16<-df.sep
df16$Rep_num<-rn[16]
df17<-df.sep
df17$Rep_num<-rn[17]
df18<-df.sep
df18$Rep_num<-rn[18]
df19<-df.sep
df19$Rep_num<-rn[19]

df.sep<-rbind(df.sep,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19)

df.scat.sep<-rbind(df.scat,df.sep)  
min(df.scat.sep$Alcaligenes)
min(df.scat.sep$Pseudomonas)

line1$A_log<-log10(line1$Alcaligenes)
line1$P_log<-log10(line1$Pseudomonas)
line2$A_log<-log10(line2$Alcaligenes)
line2$P_log<-log10(line2$Pseudomonas)

#. create separatrix at midpoint between the two shading lines
#. note that x=0 is -0.49
minv<-(-4.6)
min.v2<-(-5)
reg1 <- data.frame(x=c(min.v2, -0.54, -0.54, -1.33, -2.3, -3.3,  -4.3, min.v2 ),
                   y=c(0, 0,   -0.85, -1.6,  -2.56, -2.56,-2.56, min.v2 ))
reg2 <- data.frame(x=c(-4.3,  -2.3, -1.3, -0.47, -0.47, min.v2), 
                   y=c( -4.6,   -4.6, -2.58,-3.73, min.v2, min.v2))

sep.l1<-reg1[with(reg1, order(x)), ]
sep.l2<-reg2[with(reg2, order(x)), ]
sep.l1<-sep.l1[c(2,3,5,6,8),]
sep.l2<-sep.l2[c(1,2,3,4,5),]
sep.lx<-(sep.l1$x-sep.l2$x)/2+sep.l2$x
sep.ly<-(sep.l1$y-sep.l2$y)/2+sep.l2$y
sep.l<-data.frame(sep.lx,sep.ly)

ggplot()+ 
  geom_path(data=subset(df.scat.sep, type=='scat'), 
            aes(log10(Pseudomonas), log10(Alcaligenes),label=Transfer,colour = as.numeric(Transfer)),
            size=1, 
            arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")))+
  geom_line(data=sep.l, aes(sep.lx, sep.ly),  size=0.5,colour='gray65', linetype='dashed')+
  facet_wrap(~Rep_num)+
  scale_colour_continuous(low="white", high="black")+
  labs(x=expression(Pseudomonas*phantom(x)*(log[10])),
       y=expression(Alcaligenes*phantom(x)*(log[10])) )+
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 10, colour="gray"),
    axis.text.x = element_text(size=10, colour="gray", angle = 0, hjust = 0.5),
    axis.title.y =element_text(size=14),
    axis.title.x = element_text(size = 14),
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank())+
  guides(colour=FALSE)+
  geom_polygon(data=reg1,aes(x=x, y=y),  fill=colAlc,alpha=0.2,size=0.3)+
  geom_polygon(data=reg2,aes(x=x, y=y),  fill=colPseu,alpha=0.2,size=0.3)

ggsave(path=path.p, paste0("fig4C.pdf"), width= 9, height=7, onefile = TRUE)

dat.fig4C<-subset(df.scat.sep, type=='scat')
dat.fig4C$Alcaligenes_log10<-log10(dat.fig4C$Alcaligenes)
dat.fig4C$Pseudomonas_log10<-log10(dat.fig4C$Pseudomonas)
dat.fig4C<-select(dat.fig4C,-type,-Alcaligenes,-Pseudomonas)

## ..Write data table as csv file
fname<-paste("fig4C.csv")
write.csv(dat.fig4C,file.path(path.ds, fname),row.names=FALSE)

##-------------------------------------------------
## Fig. 4B. plot for T18 only, all reps
##-------------------------------------------------  
dp.t18<-subset(dfc2, Transfer==18)
nrep<-length(unique(dp.t18$Rep_num))

dp3<-dp.t18
df.scat<-select(dp3, Transfer, Alcaligenes,Pseudomonas,Rep_num)
df.scat$type<-"scat"
df.sep<-rbind(line1, line2)
df.sep<-select(df.sep, -A_log,-P_log)
df.sep<-select(df.sep, Alcaligenes,Pseudomonas,type)

df.sep$Transfer<-18
df.sep$Rep_num<-rn[1]
df2<-df.sep
df2$Rep_num<-rn[2]
df3<-df.sep
df3$Rep_num<-rn[3]
df4<-df.sep
df4$Rep_num<-rn[4]
df5<-df.sep
df5$Rep_num<-rn[5]
df6<-df.sep
df6$Rep_num<-rn[6]
df7<-df.sep
df7$Rep_num<-rn[7]
df8<-df.sep
df8$Rep_num<-rn[8]
df9<-df.sep
df9$Rep_num<-rn[9]
df10<-df.sep
df10$Rep_num<-rn[10]
df11<-df.sep
df11$Rep_num<-rn[11]
df12<-df.sep
df12$Rep_num<-rn[12]
df13<-df.sep
df13$Rep_num<-rn[13]
df14<-df.sep
df14$Rep_num<-rn[14]
df15<-df.sep
df15$Rep_num<-rn[15]
df16<-df.sep
df16$Rep_num<-rn[16]
df17<-df.sep
df17$Rep_num<-rn[17]
df18<-df.sep
df18$Rep_num<-rn[18]
df19<-df.sep
df19$Rep_num<-rn[19]

df.sep<-rbind(df.sep,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19)
df.scat.sep<-rbind(df.scat,df.sep)  

ggplot()+ 
  geom_point(data=subset(df.scat.sep, type=='scat'), 
            aes(log10(Pseudomonas), log10(Alcaligenes)), size=2, alpha=0.5)+
  geom_line(data=sep.l, aes(sep.lx, sep.ly), shape=19, size=0.5,colour='gray65', linetype='dashed')+
  labs(x=expression(Pseudomonas*phantom(x)*(log[10])),
       y=expression(Alcaligenes*phantom(x)*(log[10])) )+
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 16, colour="gray"),
    axis.text.x = element_text(size=16, colour="gray", angle = 0, hjust = 0.5),
    axis.title.y =element_text(size=20),
    axis.title.x = element_text(size = 20),
    panel.border = element_rect(colour = "gray90", fill=NA, size=0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank())+
  guides(colour=FALSE)+
  geom_polygon(data=reg1,aes(x=x, y=y),  fill=colAlc,alpha=0.2,size=0.3)+
  geom_polygon(data=reg2,aes(x=x, y=y),  fill=colPseu,alpha=0.2,size=0.3)

ggsave(path=path.p, paste0("fig4B.pdf"), width= 5, height=5, onefile = TRUE)

dat.fig4B<-subset(df.scat.sep, type=='scat')
dat.fig4B$Alcaligenes_log10<-log10(dat.fig4B$Alcaligenes)
dat.fig4B$Pseudomonas_log10<-log10(dat.fig4B$Pseudomonas)
dat.fig4B<- select(dat.fig4B,-type,-Alcaligenes,-Pseudomonas)

## ..Write data table as csv file
fname<-paste("fig4B.csv")
write.csv(dat.fig4B,file.path(path.ds, fname),row.names=FALSE)


