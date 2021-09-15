" R script for Figures 2C and 2D (potentials) in Estrela et al (2021) 
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
library(cowplot)
library(mixtools)
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
#.. community data 
df.all<-read.csv('Estrela_2021_16S_data_table_RelAbund_ALL.csv')
setwd(path.dp)
#.. get top ESVs sequences 
d.esvs<-read.csv('top6_ESVs_fig2B.csv')

##.. SELECT A and P at T18 
esvAlc<-as.character(subset(d.esvs,Genus=='Alcaligenes')$OTU)
esvP1<-as.character(subset(d.esvs,Genus=='Pseudomonas')$OTU)
dat.ap<-subset(df.all, Transfer==18 & OTU %in% c(esvAlc, esvP1))
dat.ap$esv_id[dat.ap$OTU==esvAlc]<-'A'
dat.ap$esv_id[dat.ap$OTU==esvP1]<-'P'
dat.ap<-select(dat.ap, -OTU,-Carbon,-Family,-Treatment,-Inoc_size)

min.ab<-0.0
dat.sub<-subset(dat.ap, Abundance >min.ab)

# .. Colour scheme for plots
#. ESVs color
colAlc<-"#E69F00"
colPseu<-  "darkorchid4" 
colKp<-blues9[7]
colKm<-blues9[5]

# .. SELECT timepoint
tn<-"18"
dfp.t<- subset(dat.sub, Transfer ==tn)
## ..Create new column with log abundance
dfp.t$Abund_log<-log10(dfp.t$Abundance)

## -----------------------------------------------------------------------------
# .. SELECT strain - Alcaligenes (see below for Pseudomonas)
stn<-"Alc"
colp<-colAlc
dt.a<-subset(dfp.t, Genus=="Alcaligenes" )
dfp<-dt.a
ab.dat<-dfp$Abund_log

##. Plots for mixed model distribution
da.mix<- normalmixEM(ab.dat, k=2, maxit=200, epsilon = 1e-16,maxrestarts=1000)
plot(da.mix, which=2)
summary(da.mix)

#. minimum calculated using Mathematica function FindMinimum (see companion Mathematica script)
mindv<-(-1.17738)
##. plot using ggplot and values from normalmixEM
sdnorm =
  function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}

# . PLOT histogram + line 
binv<-10
p2c.a<-ggplot(data.frame(x=ab.dat)) + 
  geom_histogram(aes(x=ab.dat,y=..density..), alpha=.5, color=colp, fill=colp, bins=binv) +
  stat_function(fun=sdnorm,
                args=list(mean=da.mix$mu[2],
                          sd=da.mix$sigma[2],
                          lambda=da.mix$lambda[2]),
                colour=colp, alpha=1, geom="line", size=1) +
  stat_function(fun=sdnorm,
                args=list(mean=da.mix$mu[1],
                          sd=da.mix$sigma[1],
                          lambda=da.mix$lambda[1]),
                colour=colp, alpha=1, geom="line", size=1)+
  xlim(-4,0)+
  scale_y_continuous(breaks=c(0, 2))+
  geom_hline(yintercept=-0.005, colour="white", size=1)+ #line to remove coloured baseline 
  labs( x=expression(A*phantom(x)*(log[10])),
        y="density",
        title= paste0(""),
        subtitle=paste0(""))+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size=14,  angle = 0, hjust = 0.5),
    axis.title.y =element_text(size=18),
    axis.title.x = element_text(size = 18),
    axis.line= element_line(colour="gray"),
    panel.border = element_rect(colour = NA, fill=NA, size=0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  scale_color_manual(values=c( colp ), 
                     name="")+
  scale_fill_manual(values=c(colp),
                    name="")+
  guides(fill=FALSE, colour=FALSE)
p2c.a
ggsave(path=path.p, paste0("fig2C_A.pdf"), width= 5, height=3.5, onefile = TRUE)

##.. plot U from equation for Pd and U obtained in Mathematica (see companion Mathematica script)
#.. SELECT values of x to plot the interpolated values from xmin to xmax
x <- seq(-4.4,0,by=0.01) # for alc
pot.a=(-(1/2)* log(2.21317*exp(-33.4931* (0.590556 + x)^2) + 0.280318*exp(-2.37817* (3.00682 + x)^2)))

#.plot potential using Fokker-Plank equation
mindv<-(-1.17738)
root1<-(-0.590556)
root2<-(-3.00682)

p2d.a<-ggplot(,aes(x=x, y=pot.a))+
  geom_line(colour=colp, size=1)+ 
  scale_x_reverse(expand = c(0, 0))+
  labs(y="U (A)",
       x=expression(A*phantom(x)*(log[10])),
       title= paste0(""),
       subtitle=paste0(""))+
  scale_y_continuous(breaks=c(0, 5))+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 16, angle=90),
    axis.text.x = element_text(size=14,  angle = 0, hjust = 0.5),
    axis.title.y =element_text(size=20),
    axis.title.x = element_text(size = 16),
    axis.line= element_line(colour="gray"),
    axis.ticks=element_line(colour = "gray"),
    panel.border = element_rect(colour = NA, fill=NA, size=0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  geom_vline(xintercept = c(mindv), linetype='dashed', colour='gray60')+
  geom_vline(xintercept = c(root1), linetype='dashed', colour='gray20')+
  geom_vline(xintercept = c(root2), linetype='dashed', colour='gray20')
p2d.a 

ggsave(path=path.p, paste0("fig2D_A_pot.pdf"), width= 5, height=3.5, onefile = TRUE)

## -----------------------------------------------------------------------------
# .. SELECT strain - Pseudomonas
stn<-"Pseud"
colp<-colPseu
dt.p<-subset(dfp.t, Genus=="Pseudomonas")
dfp<-dt.p

ab.dat<-dfp$Abund_log

##. Plots for mixed model distribution
da.mix<- normalmixEM(ab.dat, k=2, maxit=200, epsilon = 1e-16,maxrestarts=1000)
summary(da.mix)

#. minimum calculated using Mathematica function FindMinimum (see companion Mathematica script)
mindv<-(-1.97225)
##. plot using ggplot and values from normalmixEM
sdnorm =
  function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}

# . PLOT histogram + line 
p2c.p<-ggplot(data.frame(x=ab.dat)) + 
  geom_histogram(aes(x=ab.dat,y=..density..), alpha=.5, color=colp, fill=colp, bins=binv) +
  stat_function(fun=sdnorm,
                args=list(mean=da.mix$mu[2],
                          sd=da.mix$sigma[2],
                          lambda=da.mix$lambda[2]),
                colour=colp, alpha=1, geom="line", size=1) +
  stat_function(fun=sdnorm,
                args=list(mean=da.mix$mu[1],
                          sd=da.mix$sigma[1],
                          lambda=da.mix$lambda[1]),
                colour=colp, alpha=1, geom="line", size=1)+
  xlim(-4.5,0)+
  scale_y_continuous(breaks=c(0, 1))+
  geom_hline(yintercept=-0.005, colour="white", size=1)+ #line to remove coloured baseline 
  labs(x=expression(P*phantom(x)*(log[10])),
        y="density",
        title= paste0(""),
        subtitle=paste0(""))+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size=14,  angle = 0, hjust = 0.5),
    axis.title.y =element_text(size=18),
    axis.title.x = element_text(size = 18),
    axis.line= element_line(colour="gray"),
    panel.border = element_rect(colour = NA, fill=NA, size=0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  scale_color_manual(values=c( colp ), 
                     name="")+
  scale_fill_manual(values=c(colp),
                    name="")+
  guides(fill=FALSE, colour=FALSE)
p2c.p
ggsave(path=path.p, paste0("fig2C_P.pdf"), width= 5, height=3.5, onefile = TRUE)

##.. plot U from equation for Pd and U obtained in Mathematica
pot.a=(-(1/2)* log(1.13421*exp(-6.63602 * (1.08247 + x)^2) + 0.117626*exp(-0.901286 * (3.16649 + x)^2)))

##- roots were calculated using mathematica (see companion Mathematica script)
mindv<-(-1.97225)
root1<-(-1.08306)
root2<-(-3.16649)
p2d.p<-ggplot(,aes(x=x, y=pot.a))+
  geom_line(colour=colp, size=1)+ 
  scale_x_reverse(expand = c(0, 0))+
  labs(y="U(P)",
       x=expression(P*phantom(x)*(log[10])),
       title= paste0(""),
       subtitle=paste0(""))+
  xlim(0,-4.4)+
  scale_y_continuous(breaks=c(0, 4))+
  theme_classic()+
  theme(
    axis.text.y = element_text(size = 16, angle=90),
    axis.text.x = element_text(size=14,  angle = 0, hjust = 0.5),
    axis.title.y =element_text(size=20),
    axis.title.x = element_text(size = 16),
    axis.line= element_line(colour="gray"),
    axis.ticks=element_line(colour = "gray"),
    panel.border = element_rect(colour = NA, fill=NA, size=0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  geom_vline(xintercept = c(mindv), linetype='dashed', colour='gray60')+
  geom_vline(xintercept = c(root1), linetype='dashed', colour='gray20')+
  geom_vline(xintercept = c(root2), linetype='dashed', colour='gray20')
p2d.p 
ggsave(path=path.p, paste0("fig2D_P_pot.pdf"), width= 5, height=3.5, onefile = TRUE)

#.. write data figure as csv file
ds<-dfp.t
fname<-paste("fig2C.csv")
write.csv(ds, file.path(path.ds, fname), row.names=FALSE)




