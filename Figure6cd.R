library('vegan')
library('phyloseq')
library('stringr')
library('forcats')
library('cowplot')
library('knitr')
library('ggplot2')
library('PERMANOVA')
library('pairwiseAdonis')
library('ggsignif')
library('tidyverse')
library('ggpubr')
library('betapart')
library('stats')
library('ape')
library('phytools')
library('qiime2R')
library('phyloseq')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(10)

Lphyloseq<-qza_to_phyloseq(
  features="study_1036_051125-195547/BIOM/130932/feature-table.qza",
  #"study_1036_051125-195547/BIOM/130932/taxonomy.qza",
  tree="study_1036_051125-195547/BIOM/130932/rooted-tree.qza",
  
  metadata = "1036_20230201-070300.txt")
Lphyloseq<-prune_taxa(taxa_sums(Lphyloseq)>0,Lphyloseq)
Lphyloseq<-rarefy_even_depth(Lphyloseq,rngseed=1,sample.size = 5000)

bog<-prune_samples(sample_data(Lphyloseq)$site=='bog',Lphyloseq)
active<-prune_samples(sample_data(Lphyloseq)$site=='active layer',Lphyloseq)
perm<-prune_samples(sample_data(Lphyloseq)$site=='permafrost',Lphyloseq)

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('bog','permafrost','active')

Zbogvactive<-merge_phyloseq(bog,active)
Zbogvperm<-merge_phyloseq(bog,perm)
Zactivevperm<-merge_phyloseq(active,perm)

#Zbogvactive_s<-tip_glom(Zbogvactive,h=0.05)
#Zbogvactive_g<-tip_glom(Zbogvactive,h=0.1)
#Zbogvactive_f<-tip_glom(Zbogvactive,h=0.15)
#Zbogvactive_o<-tip_glom(Zbogvactive,h=0.2)

Zbogvperm_s<-tip_glom(Zbogvperm,h=0.05)
Zbogvperm_g<-tip_glom(Zbogvperm,h=0.1)
Zbogvperm_f<-tip_glom(Zbogvperm,h=0.15)
Zbogvperm_o<-tip_glom(Zbogvperm,h=0.2)

Zactivevperm_s<-tip_glom(Zactivevperm,h=0.05)
Zactivevperm_g<-tip_glom(Zactivevperm,h=0.1)
Zactivevperm_f<-tip_glom(Zactivevperm,h=0.15)
Zactivevperm_o<-tip_glom(Zactivevperm,h=0.2)



###########################Branch-based tree UniFrac###########################

lkf<-c()
lke<-c()
lkr<-c()

#branch-based ASV
lkf<-c(lkf,coreUniFrac(Zbogvperm,sample_data(Zbogvperm)$site,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zbogvperm_s,sample_data(Zbogvperm_s)$site,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zbogvperm_g,sample_data(Zbogvperm_g)$site,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zbogvperm_f,sample_data(Zbogvperm_f)$site,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zbogvperm_o,sample_data(Zbogvperm_o)$site,core_fraction=0.5))

#tip-based ASV
lke<-c(lke,coreUniFrac(Zbogvperm,sample_data(Zbogvperm)$site,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zbogvperm_s,sample_data(Zbogvperm_s)$site,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zbogvperm_g,sample_data(Zbogvperm_g)$site,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zbogvperm_f,sample_data(Zbogvperm_f)$site,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zbogvperm_o,sample_data(Zbogvperm_o)$site,mode='tip',core_fraction=0.5))

#non-phylogenetic ASV
lkr<-c(lkr,coreJaccard(Zbogvperm,sample_data(Zbogvperm)$site,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zbogvperm_s,sample_data(Zbogvperm_s)$site,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zbogvperm_g,sample_data(Zbogvperm_g)$site,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zbogvperm_f,sample_data(Zbogvperm_f)$site,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zbogvperm_o,sample_data(Zbogvperm_o)$site,core_fraction=0.5))


lkf2<-c()
lke2<-c()
lkr2<-c()

#branch-based ASV
lkf2<-c(lkf2,coreUniFrac(Zactivevperm,sample_data(Zactivevperm)$site,core_fraction=0.5))
lkf2<-c(lkf2,coreUniFrac(Zactivevperm_s,sample_data(Zactivevperm_s)$site,core_fraction=0.5))
lkf2<-c(lkf2,coreUniFrac(Zactivevperm_g,sample_data(Zactivevperm_g)$site,core_fraction=0.5))
lkf2<-c(lkf2,coreUniFrac(Zactivevperm_f,sample_data(Zactivevperm_f)$site,core_fraction=0.5))
lkf2<-c(lkf2,coreUniFrac(Zactivevperm_o,sample_data(Zactivevperm_o)$site,core_fraction=0.5))

#tip-based ASV
lke2<-c(lke2,coreUniFrac(Zactivevperm,sample_data(Zactivevperm)$site,mode='tip',core_fraction=0.5))
lke2<-c(lke2,coreUniFrac(Zactivevperm_s,sample_data(Zactivevperm_s)$site,mode='tip',core_fraction=0.5))
lke2<-c(lke2,coreUniFrac(Zactivevperm_g,sample_data(Zactivevperm_g)$site,mode='tip',core_fraction=0.5))
lke2<-c(lke2,coreUniFrac(Zactivevperm_f,sample_data(Zactivevperm_f)$site,mode='tip',core_fraction=0.5))
lke2<-c(lke2,coreUniFrac(Zactivevperm_o,sample_data(Zactivevperm_o)$site,mode='tip',core_fraction=0.5))

#non-phylogenetic ASV
lkr2<-c(lkr2,coreJaccard(Zactivevperm,sample_data(Zactivevperm)$site,core_fraction=0.5))
lkr2<-c(lkr2,coreJaccard(Zactivevperm_s,sample_data(Zactivevperm_s)$site,core_fraction=0.5))
lkr2<-c(lkr2,coreJaccard(Zactivevperm_g,sample_data(Zactivevperm_g)$site,core_fraction=0.5))
lkr2<-c(lkr2,coreJaccard(Zactivevperm_f,sample_data(Zactivevperm_f)$site,core_fraction=0.5))
lkr2<-c(lkr2,coreJaccard(Zactivevperm_o,sample_data(Zactivevperm_o)$site,core_fraction=0.5))


mkf<-c()
mke<-c()
mkr<-c()

#branch-based ASV
mkf<-c(mkf,coreUniFrac(Zbogvperm,sample_data(Zbogvperm)$site,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zbogvperm_s,sample_data(Zbogvperm_s)$site,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zbogvperm_g,sample_data(Zbogvperm_g)$site,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zbogvperm_f,sample_data(Zbogvperm_f)$site,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zbogvperm_o,sample_data(Zbogvperm_o)$site,core_fraction=0.5,ab_threshold1=0.1))

#tip-based ASV
mke<-c(mke,coreUniFrac(Zbogvperm,sample_data(Zbogvperm)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zbogvperm_s,sample_data(Zbogvperm_s)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zbogvperm_g,sample_data(Zbogvperm_g)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zbogvperm_f,sample_data(Zbogvperm_f)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zbogvperm_o,sample_data(Zbogvperm_o)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))

#non-phylogenetic ASV
mkr<-c(mkr,coreJaccard(Zbogvperm,sample_data(Zbogvperm)$site,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zbogvperm_s,sample_data(Zbogvperm_s)$site,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zbogvperm_g,sample_data(Zbogvperm_g)$site,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zbogvperm_f,sample_data(Zbogvperm_f)$site,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zbogvperm_o,sample_data(Zbogvperm_o)$site,core_fraction=0.5,ab_threshold1=0.1))

mkf2<-c()
mke2<-c()
mkr2<-c()

#branch-based ASV
mkf2<-c(mkf2,coreUniFrac(Zactivevperm,sample_data(Zactivevperm)$site,core_fraction=0.5,ab_threshold1=0.1))
mkf2<-c(mkf2,coreUniFrac(Zactivevperm_s,sample_data(Zactivevperm_s)$site,core_fraction=0.5,ab_threshold1=0.1))
mkf2<-c(mkf2,coreUniFrac(Zactivevperm_g,sample_data(Zactivevperm_g)$site,core_fraction=0.5,ab_threshold1=0.1))
mkf2<-c(mkf2,coreUniFrac(Zactivevperm_f,sample_data(Zactivevperm_f)$site,core_fraction=0.5,ab_threshold1=0.1))
mkf2<-c(mkf2,coreUniFrac(Zactivevperm_o,sample_data(Zactivevperm_o)$site,core_fraction=0.5,ab_threshold1=0.1))

#tip-based ASV
mke2<-c(mke2,coreUniFrac(Zactivevperm,sample_data(Zactivevperm)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke2<-c(mke2,coreUniFrac(Zactivevperm_s,sample_data(Zactivevperm_s)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke2<-c(mke2,coreUniFrac(Zactivevperm_g,sample_data(Zactivevperm_g)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke2<-c(mke2,coreUniFrac(Zactivevperm_f,sample_data(Zactivevperm_f)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke2<-c(mke2,coreUniFrac(Zactivevperm_o,sample_data(Zactivevperm_o)$site,mode='tip',core_fraction=0.5,ab_threshold1=0.1))

#non-phylogenetic ASV
mkr2<-c(mkr2,coreJaccard(Zactivevperm,sample_data(Zactivevperm)$site,core_fraction=0.5,ab_threshold1=0.1))
mkr2<-c(mkr2,coreJaccard(Zactivevperm_s,sample_data(Zactivevperm_s)$site,core_fraction=0.5,ab_threshold1=0.1))
mkr2<-c(mkr2,coreJaccard(Zactivevperm_g,sample_data(Zactivevperm_g)$site,core_fraction=0.5,ab_threshold1=0.1))
mkr2<-c(mkr2,coreJaccard(Zactivevperm_f,sample_data(Zactivevperm_f)$site,core_fraction=0.5,ab_threshold1=0.1))
mkr2<-c(mkr2,coreJaccard(Zactivevperm_o,sample_data(Zactivevperm_o)$site,core_fraction=0.5,ab_threshold1=0.1))


nkf<-c()
nke<-c()
nkr<-c()

#branch-based ASV
nkf<-c(nkf,coreUniFrac(Zbogvperm,sample_data(Zbogvperm)$site,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zbogvperm_s,sample_data(Zbogvperm_s)$site,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zbogvperm_g,sample_data(Zbogvperm_g)$site,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zbogvperm_f,sample_data(Zbogvperm_f)$site,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zbogvperm_o,sample_data(Zbogvperm_o)$site,selection='shade',initial_branches=10))

#tip-based ASV
nke<-c(nke,coreUniFrac(Zbogvperm,sample_data(Zbogvperm)$site,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zbogvperm_s,sample_data(Zbogvperm_s)$site,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zbogvperm_g,sample_data(Zbogvperm_g)$site,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zbogvperm_f,sample_data(Zbogvperm_f)$site,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zbogvperm_o,sample_data(Zbogvperm_o)$site,mode='tip',selection='shade',initial_branches=10))

#non-phylogenetic ASV
nkr<-c(nkr,coreJaccard(Zbogvperm,sample_data(Zbogvperm)$site,selection='shade'))
nkr<-c(nkr,coreJaccard(Zbogvperm_s,sample_data(Zbogvperm_s)$site,selection='shade'))
nkr<-c(nkr,coreJaccard(Zbogvperm_g,sample_data(Zbogvperm_g)$site,selection='shade'))
nkr<-c(nkr,coreJaccard(Zbogvperm_f,sample_data(Zbogvperm_f)$site,selection='shade'))
nkr<-c(nkr,coreJaccard(Zbogvperm_o,sample_data(Zbogvperm_o)$site,selection='shade'))

nkf2<-c()
nke2<-c()
nkr2<-c()

#branch-based ASV
nkf2<-c(nkf2,coreUniFrac(Zactivevperm,sample_data(Zactivevperm)$site,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreUniFrac(Zactivevperm_s,sample_data(Zactivevperm_s)$site,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreUniFrac(Zactivevperm_g,sample_data(Zactivevperm_g)$site,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreUniFrac(Zactivevperm_f,sample_data(Zactivevperm_f)$site,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreUniFrac(Zactivevperm_o,sample_data(Zactivevperm_o)$site,selection='shade',initial_branches=10))

#tip-based ASV
nke2<-c(nke2,coreUniFrac(Zactivevperm,sample_data(Zactivevperm)$site,mode='tip',selection='shade',initial_branches=10))
nke2<-c(nke2,coreUniFrac(Zactivevperm_s,sample_data(Zactivevperm_s)$site,mode='tip',selection='shade',initial_branches=10))
nke2<-c(nke2,coreUniFrac(Zactivevperm_g,sample_data(Zactivevperm_g)$site,mode='tip',selection='shade',initial_branches=10))
nke2<-c(nke2,coreUniFrac(Zactivevperm_f,sample_data(Zactivevperm_f)$site,mode='tip',selection='shade',initial_branches=10))
nke2<-c(nke2,coreUniFrac(Zactivevperm_o,sample_data(Zactivevperm_o)$site,mode='tip',selection='shade',initial_branches=10))

#non-phylogenetic ASV
nkr2<-c(nkr2,coreJaccard(Zactivevperm,sample_data(Zactivevperm)$site,selection='shade'))
nkr2<-c(nkr2,coreJaccard(Zactivevperm_s,sample_data(Zactivevperm_s)$site,selection='shade'))
nkr2<-c(nkr2,coreJaccard(Zactivevperm_g,sample_data(Zactivevperm_g)$site,selection='shade'))
nkr2<-c(nkr2,coreJaccard(Zactivevperm_f,sample_data(Zactivevperm_f)$site,selection='shade'))
nkr2<-c(nkr2,coreJaccard(Zactivevperm_o,sample_data(Zactivevperm_o)$site,selection='shade'))

#Plot variation in Faith's PD across taxonoperm ranks
axisv<-c(1,2,3,4,5)

par(mar=c(5,5,5,5))
plot(axisv,lkf,col='brown',ylim=c(0.4,0.9),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('0','0.05','0.10','0.15','0.20'))
points(axisv,lkf,col='black',pch=22,cex=1.5)
lines(axisv,lkf,col='brown',lty='solid',lwd=2)
points(axisv,lkr,col='brown',pch=16,cex=1.5)
points(axisv,lkr,col='black',pch=21,cex=1.5)
lines(axisv,lkr,col='brown',lty='dashed',lwd=2)
points(axisv,lke,col='brown',pch=17,cex=1.5)
points(axisv,lke,col='black',pch=24,cex=1.5)
lines(axisv,lke,col='brown',lty='dotted',lwd=2)

points(axisv,lkf2,col='green',pch=15,cex=1.5)
points(axisv,lkf2,col='black',pch=22,cex=1.5)
lines(axisv,lkf2,col='green',lty='solid',lwd=2)
points(axisv,lkr2,col='green',pch=16,cex=1.5)
points(axisv,lkr2,col='black',pch=21,cex=1.5)
lines(axisv,lkr2,col='green',lty='dashed',lwd=2)
points(axisv,lke2,col='green',pch=17,cex=1.5)
points(axisv,lke2,col='black',pch=24,cex=1.5)
lines(axisv,lke2,col='green',lty='dotted',lwd=2)

par(mar=c(5,5,5,5))
plot(axisv,mkf,col='tan3',ylim=c(0.4,0.85),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('0','0.05','0.10','0.15','0.20'))
points(axisv,mkf,col='black',pch=22,cex=1.5)
lines(axisv,mkf,col='tan3',lty='solid',lwd=2)
points(axisv,mkr,col='tan3',pch=16,cex=1.5)
points(axisv,mkr,col='black',pch=21,cex=1.5)
lines(axisv,mkr,col='tan3',lty='dashed',lwd=2)
points(axisv,mke,col='tan3',pch=17,cex=1.5)
points(axisv,mke,col='black',pch=24,cex=1.5)
lines(axisv,mke,col='tan3',lty='dotted',lwd=2)

points(axisv,mkf2,col='palegreen2',pch=15,cex=1.5)
points(axisv,mkf2,col='black',pch=22,cex=1.5)
lines(axisv,mkf2,col='palegreen2',lty='solid',lwd=2)
points(axisv,mkr2,col='palegreen2',pch=16,cex=1.5)
points(axisv,mkr2,col='black',pch=21,cex=1.5)
lines(axisv,mkr2,col='palegreen2',lty='dashed',lwd=2)
points(axisv,mke2,col='palegreen2',pch=17,cex=1.5)
points(axisv,mke2,col='black',pch=24,cex=1.5)

plot(axisv,nkf,col='burlywood2',ylim=c(0.4,0.9),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:6, labels=c('0','0.05','0.1','0.15','0.2','0.25'))
points(axisv,nkf,col='black',pch=22,cex=1.5)
lines(axisv,nkf,col='burlywood2',lty='solid',lwd=2)
points(axisv,nkr,col='burlywood2',pch=16,cex=1.5)
points(axisv,nkr,col='black',pch=21,cex=1.5)
lines(axisv,nkr,col='burlywood2',lty='dashed',lwd=2)
points(axisv,nke,col='burlywood2',pch=17,cex=1.5)
points(axisv,nke,col='black',pch=24,cex=1.5)
lines(axisv,nke,col='burlywood2',lty='dotted',lwd=2)

points(axisv,nkf2,col='darkseagreen1',pch=15,cex=1.5)
points(axisv,nkf2,col='black',pch=22,cex=1.5)
lines(axisv,nkf2,col='darkseagreen1',lty='solid',lwd=2)
points(axisv,nkr2,col='darkseagreen1',pch=16,cex=1.5)
points(axisv,nkr2,col='black',pch=21,cex=1.5)
lines(axisv,nkr2,col='darkseagreen1',lty='dashed',lwd=2)
points(axisv,nke2,col='darkseagreen1',pch=17,cex=1.5)
points(axisv,nke2,col='black',pch=24,cex=1.5)
lines(axisv,nke2,col='darkseagreen1',lty='dotted',lwd=2)

#Find maximum and minimum across ranks for each metric
x = 1:18
High = c(max(lkf),max(lke),max(lkr),max(mkf),max(mke),max(mkr),max(nkf),max(nke),max(nkr),max(lkf2),max(lke2),max(lkr2),max(mkf2),max(mke2),max(mkr2),max(nkf2),max(nke2),max(nkr2))
Low = c(min(lkf),min(lke),min(lkr),min(mkf),min(mke),min(mkr),min(nkf),min(nke),min(nkr),min(lkf2),min(lke2),min(lkr2),min(mkf2),min(mke2),min(mkr2),min(nkf2),min(nke2),min(nkr2))

## Plot variation across ranks
par(mar=c(10,5,1,3))
plot(1, type="n", xlab="", ylab="Turnover Range", xlim=c(0.5,18.5),
     ylim=c(0.4, 0.9),xaxt='n',cex.lab=1.25)

## Add rectangles
rect(x[1:3] - 0.4, Low[1:3], x[1:3] + 0.4, High[1:3], col="brown")
## Add rectangles
rect(x[4:6] - 0.4, Low[4:6], x[4:6] + 0.4, High[4:6], col="tan3")
## Add rectangles
rect(x[7:9] - 0.4, Low[7:9], x[7:9] + 0.4, High[7:9], col="burlywood2")
## Add rectangles
rect(x[10:12] - 0.4, Low[10:12], x[10:12] + 0.4, High[10:12], col="green")
## Add rectangles
rect(x[13:15] - 0.4, Low[13:15], x[13:15] + 0.4, High[13:15], col="palegreen2")
## Add rectangles
rect(x[16:18] - 0.4, Low[16:18], x[16:18] + 0.4, High[16:18], col="darkseagreen1")
## Add rectangles
axis(1, at=1:18, labels=c('occ-BranchUniFrac','occ-TipUniFrac','occ-Jaccard','occ-ab-BranchUniFrac','occ-ab-TipUniFrac','occ-ab-Jaccard','shade-BranchUniFrac','shade-TipUniFrac','shade-Jaccard','occ-BranchUniFrac','occ-TipUniFrac','occ-Jaccard','occ-ab-BranchUniFrac','occ-ab-TipUniFrac','occ-ab-Jaccard','shade-BranchUniFrac','shade-TipUniFrac','shade-Jaccard'),las=2)

