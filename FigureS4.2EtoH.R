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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lphyloseq<-qza_to_phyloseq(
  features="study_1041_042825-091934/BIOM/131527/feature-table.qza",
  "study_1041_042825-091934/BIOM/131527/taxonomy.qza",
  tree="study_1041_042825-091934/BIOM/131527/rooted-tree.qza",
  
  metadata = "study_1041_042825-091934/1041_20230201-070300.txt")
Lphyloseq<-prune_taxa(taxa_sums(Lphyloseq)>0,Lphyloseq)
Lphyloseq<-prune_samples(sample_data(Lphyloseq)$depth_m!='not applicable',Lphyloseq)
Lphyloseq<-rarefy_even_depth(Lphyloseq,sample.size = 20000,rngseed=1)

Sup<-prune_samples(sample_data(Lphyloseq)$Lake=='Superior',Lphyloseq)
Mic<-prune_samples(sample_data(Lphyloseq)$Lake=='Michigan',Lphyloseq)

Zv<-merge_phyloseq(Sup,Mic)

Zv_s<-tip_glom(Zv,h=0.05)
Zv_g<-tip_glom(Zv,h=0.1)
Zv_f<-tip_glom(Zv,h=0.15)
Zv_o<-tip_glom(Zv,h=0.2)




###########################Branch-based tree UniFrac###########################

lkf<-c()
lke<-c()
lkr<-c()

#branch-based ASV
lkf<-c(lkf,coreUniFrac(Zv,sample_data(Zv)$Lake,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zv_s,sample_data(Zv_s)$Lake,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zv_g,sample_data(Zv_g)$Lake,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zv_f,sample_data(Zv_f)$Lake,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zv_o,sample_data(Zv_o)$Lake,core_fraction=0.5))

#tip-based ASV
lke<-c(lke,coreUniFrac(Zv,sample_data(Zv)$Lake,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zv_s,sample_data(Zv_s)$Lake,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zv_g,sample_data(Zv_g)$Lake,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zv_f,sample_data(Zv_f)$Lake,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zv_o,sample_data(Zv_o)$Lake,mode='tip',core_fraction=0.5))

#non-phylogenetic ASV
lkr<-c(lkr,coreJaccard(Zv,sample_data(Zv)$Lake,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zv_s,sample_data(Zv_s)$Lake,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zv_g,sample_data(Zv_g)$Lake,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zv_f,sample_data(Zv_f)$Lake,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zv_o,sample_data(Zv_o)$Lake,core_fraction=0.5))




mkf<-c()
mke<-c()
mkr<-c()

#branch-based ASV
mkf<-c(mkf,coreUniFrac(Zv,sample_data(Zv)$Lake,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zv_s,sample_data(Zv_s)$Lake,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zv_g,sample_data(Zv_g)$Lake,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zv_f,sample_data(Zv_f)$Lake,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zv_o,sample_data(Zv_o)$Lake,core_fraction=0.5,ab_threshold1=0.1))

#tip-based ASV
mke<-c(mke,coreUniFrac(Zv,sample_data(Zv)$Lake,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zv_s,sample_data(Zv_s)$Lake,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zv_g,sample_data(Zv_g)$Lake,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zv_f,sample_data(Zv_f)$Lake,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zv_o,sample_data(Zv_o)$Lake,mode='tip',core_fraction=0.5,ab_threshold1=0.1))

#non-phylogenetic ASV
mkr<-c(mkr,coreJaccard(Zv,sample_data(Zv)$Lake,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zv_s,sample_data(Zv_s)$Lake,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zv_g,sample_data(Zv_g)$Lake,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zv_f,sample_data(Zv_f)$Lake,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zv_o,sample_data(Zv_o)$Lake,core_fraction=0.5,ab_threshold1=0.1))


nkf<-c()
nke<-c()
nkr<-c()

#branch-based ASV
nkf<-c(nkf,coreUniFrac(Zv,sample_data(Zv)$Lake,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zv_s,sample_data(Zv_s)$Lake,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zv_g,sample_data(Zv_g)$Lake,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zv_f,sample_data(Zv_f)$Lake,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zv_o,sample_data(Zv_o)$Lake,selection='shade',initial_branches=10))

#tip-based ASV
nke<-c(nke,coreUniFrac(Zv,sample_data(Zv)$Lake,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zv_s,sample_data(Zv_s)$Lake,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zv_g,sample_data(Zv_g)$Lake,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zv_f,sample_data(Zv_f)$Lake,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zv_o,sample_data(Zv_o)$Lake,mode='tip',selection='shade',initial_branches=10))

#non-phylogenetic ASV
nkr<-c(nkr,coreJaccard(Zv,sample_data(Zv)$Lake,selection='shade'))
nkr<-c(nkr,coreJaccard(Zv_s,sample_data(Zv_s)$Lake,selection='shade'))
nkr<-c(nkr,coreJaccard(Zv_g,sample_data(Zv_g)$Lake,selection='shade'))
nkr<-c(nkr,coreJaccard(Zv_f,sample_data(Zv_f)$Lake,selection='shade'))
nkr<-c(nkr,coreJaccard(Zv_o,sample_data(Zv_o)$Lake,selection='shade'))


#Plot variation in Faith's PD across taxonoZneo ranks
axisv<-c(1,2,3,4,5)

par(mar=c(5,5,5,5))
plot(axisv,lkf,col='green',ylim=c(0.2,0.55),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('ASV','species','genus','family','order'))
points(axisv,lkf,col='black',pch=22,cex=1.5)
lines(axisv,lkf,col='green',lty='solid',lwd=2)
points(axisv,lkr,col='green',pch=16,cex=1.5)
points(axisv,lkr,col='black',pch=21,cex=1.5)
lines(axisv,lkr,col='green',lty='dashed',lwd=2)
points(axisv,lke,col='green',pch=17,cex=1.5)
points(axisv,lke,col='black',pch=24,cex=1.5)
lines(axisv,lke,col='green',lty='dotted',lwd=2)


par(mar=c(5,5,5,5))
plot(axisv,mkf,col='palegreen2',ylim=c(0.2,0.5),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('ASV','species','genus','family','order'))
points(axisv,mkf,col='black',pch=22,cex=1.5)
lines(axisv,mkf,col='palegreen2',lty='solid',lwd=2)
points(axisv,mkr,col='palegreen2',pch=16,cex=1.5)
points(axisv,mkr,col='black',pch=21,cex=1.5)
lines(axisv,mkr,col='palegreen2',lty='dashed',lwd=2)
points(axisv,mke,col='palegreen2',pch=17,cex=1.5)
points(axisv,mke,col='black',pch=24,cex=1.5)
lines(axisv,mke,col='palegreen2',lty='dotted',lwd=2)


plot(axisv,nkf,col='darkseagreen1',ylim=c(0,0.8),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:6, labels=c('0','0.05','0.1','0.15','0.2','0.25'))
points(axisv,nkf,col='black',pch=22,cex=1.5)
lines(axisv,nkf,col='darkseagreen1',lty='solid',lwd=2)
points(axisv,nkr,col='darkseagreen1',pch=16,cex=1.5)
points(axisv,nkr,col='black',pch=21,cex=1.5)
lines(axisv,nkr,col='darkseagreen1',lty='dashed',lwd=2)
points(axisv,nke,col='darkseagreen1',pch=17,cex=1.5)
points(axisv,nke,col='black',pch=24,cex=1.5)
lines(axisv,nke,col='darkseagreen1',lty='dotted',lwd=2)


#Find maximum and minimum across ranks for each metric
x = 1:9
High = c(max(lkf),max(lke),max(lkr),max(mkf),max(mke),max(mkr),max(nkf),max(nke),max(nkr),max(lkf2),max(lke2),max(lkr2),max(mkf2),max(mke2),max(mkr2),max(nkf2),max(nke2),max(nkr2))
Low = c(min(lkf),min(lke),min(lkr),min(mkf),min(mke),min(mkr),min(nkf),min(nke),min(nkr),min(lkf2),min(lke2),min(lkr2),min(mkf2),min(mke2),min(mkr2),min(nkf2),min(nke2),min(nkr2))

## Plot variation across ranks
par(mar=c(10,5,1,3))
plot(1, type="n", xlab="", ylab="Turnover Range", xlim=c(0.5,9.5),
     ylim=c(0, 1),xaxt='n',cex.lab=1.25)

## Add rectangles
rect(x[1:3] - 0.4, Low[1:3], x[1:3] + 0.4, High[1:3], col="green")
## Add rectangles
rect(x[4:6] - 0.4, Low[4:6], x[4:6] + 0.4, High[4:6], col="palegreen2")
## Add rectangles
rect(x[7:9] - 0.4, Low[7:9], x[7:9] + 0.4, High[7:9], col="darkseagreen1")
## Add rectangles
## Add rectangles
axis(1, at=1:9, labels=c('occ-BranchUniFrac','occ-TipUniFrac','occ-Jaccard','occ-ab-BranchUniFrac','occ-ab-TipUniFrac','occ-ab-Jaccard','shade-BranchUniFrac','shade-TipUniFrac','shade-Jaccard'),las=2)

