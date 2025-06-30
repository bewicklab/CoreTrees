library('vegan')
library('phyloseq')
library('stringr')
library('forcats')
library('cowplot')
library('knitr')
library('ggplot2')
library('pairwiseAdonis')
library('ggsignif')
library('tidyverse')
library('ggpubr')
library('betapart')
library('stats')
library('ape')
library('phytools')
library('PERMANOVA')
library('dplyr')
library('data.table')
library('holobiont')
library('qiime2R')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Lphyloseq<-qza_to_phyloseq(
  features="study_1036_051125-195547/BIOM/130932/feature-table.qza",
  "study_1036_051125-195547/BIOM/130932/taxonomy.qza",
  tree="study_1036_051125-195547/BIOM/130932/rooted-tree.qza",
  
  metadata = "1036_20230201-070300.txt")
Lphyloseq<-prune_taxa(taxa_sums(Lphyloseq)>0,Lphyloseq)
Lphyloseq<-rarefy_even_depth(Lphyloseq,rngseed=1)

bog<-prune_samples(sample_data(Lphyloseq)$site=='bog',Lphyloseq)
active<-prune_samples(sample_data(Lphyloseq)$site=='active layer',Lphyloseq)
perm<-prune_samples(sample_data(Lphyloseq)$site=='permafrost',Lphyloseq)


bog_s1<-tip_glom(bog,h=0.05)
active_s1<-tip_glom(active,h=0.05)
perm_s1<-tip_glom(perm,h=0.05)

bog_s2<-tip_glom(bog,h=0.1)
active_s2<-tip_glom(active,h=0.1)
perm_s2<-tip_glom(perm,h=0.1)

bog_s3<-tip_glom(bog,h=0.15)
active_s3<-tip_glom(active,h=0.15)
perm_s3<-tip_glom(perm,h=0.15)

bog_s4<-tip_glom(bog,h=0.2)
active_s4<-tip_glom(active,h=0.2)
perm_s4<-tip_glom(perm,h=0.2)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate variation in Faith's PD across taxonoperm rank

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lkf<-c()
lke<-c()
lkr<-c()

#branch-based ASV
lkf<-c(lkf,coreFaithsPD(bog,core_fraction = 0.5)/coreFaithsPD(perm,core_fraction = 0.5))
lkf<-c(lkf,coreFaithsPD(bog_s1,core_fraction = 0.5)/coreFaithsPD(perm_s1,core_fraction = 0.5))
lkf<-c(lkf,coreFaithsPD(bog_s2,core_fraction = 0.5)/coreFaithsPD(perm_s2,core_fraction = 0.5))
lkf<-c(lkf,coreFaithsPD(bog_s3,core_fraction = 0.5)/coreFaithsPD(perm_s3,core_fraction = 0.5))
lkf<-c(lkf,coreFaithsPD(bog_s4,core_fraction = 0.5)/coreFaithsPD(perm_s4,core_fraction = 0.5))

#tip-based ASV
lke<-c(lke,coreFaithsPD(bog,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm,mode='tip',core_fraction = 0.5))
lke<-c(lke,coreFaithsPD(bog_s1,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm_s1,mode='tip',core_fraction = 0.5))
lke<-c(lke,coreFaithsPD(bog_s2,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm_s2,mode='tip',core_fraction = 0.5))
lke<-c(lke,coreFaithsPD(bog_s3,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm_s3,mode='tip',core_fraction = 0.5))
lke<-c(lke,coreFaithsPD(bog_s4,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm_s4,mode='tip',core_fraction = 0.5))

#non-phylogenetic ASV
lkr<-c(lkr,coreRichness(bog,core_fraction = 0.5)/coreRichness(perm,core_fraction = 0.5))
lkr<-c(lkr,coreRichness(bog_s1,core_fraction = 0.5)/coreRichness(perm_s1,core_fraction = 0.5))
lkr<-c(lkr,coreRichness(bog_s2,core_fraction = 0.5)/coreRichness(perm_s2,core_fraction = 0.5))
lkr<-c(lkr,coreRichness(bog_s3,core_fraction = 0.5)/coreRichness(perm_s3,core_fraction = 0.5))
lkr<-c(lkr,coreRichness(bog_s4,core_fraction = 0.5)/coreRichness(perm_s4,core_fraction = 0.5))

lkf2<-c()
lke2<-c()
lkr2<-c()

#branch-based ASV
lkf2<-c(lkf2,coreFaithsPD(active,core_fraction = 0.5)/coreFaithsPD(perm,core_fraction = 0.5))
lkf2<-c(lkf2,coreFaithsPD(active_s1,core_fraction = 0.5)/coreFaithsPD(perm_s1,core_fraction = 0.5))
lkf2<-c(lkf2,coreFaithsPD(active_s2,core_fraction = 0.5)/coreFaithsPD(perm_s2,core_fraction = 0.5))
lkf2<-c(lkf2,coreFaithsPD(active_s3,core_fraction = 0.5)/coreFaithsPD(perm_s3,core_fraction = 0.5))
lkf2<-c(lkf2,coreFaithsPD(active_s4,core_fraction = 0.5)/coreFaithsPD(perm_s4,core_fraction = 0.5))

#tip-based ASV
lke2<-c(lke2,coreFaithsPD(active,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm,mode='tip',core_fraction = 0.5))
lke2<-c(lke2,coreFaithsPD(active_s1,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm_s1,mode='tip',core_fraction = 0.5))
lke2<-c(lke2,coreFaithsPD(active_s2,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm_s2,mode='tip',core_fraction = 0.5))
lke2<-c(lke2,coreFaithsPD(active_s3,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm_s3,mode='tip',core_fraction = 0.5))
lke2<-c(lke2,coreFaithsPD(active_s4,mode='tip',core_fraction = 0.5)/coreFaithsPD(perm_s4,mode='tip',core_fraction = 0.5))

#non-phylogenetic ASV
lkr2<-c(lkr2,coreRichness(active,core_fraction = 0.5)/coreRichness(perm,core_fraction = 0.5))
lkr2<-c(lkr2,coreRichness(active_s1,core_fraction = 0.5)/coreRichness(perm_s1,core_fraction = 0.5))
lkr2<-c(lkr2,coreRichness(active_s2,core_fraction = 0.5)/coreRichness(perm_s2,core_fraction = 0.5))
lkr2<-c(lkr2,coreRichness(active_s3,core_fraction = 0.5)/coreRichness(perm_s3,core_fraction = 0.5))
lkr2<-c(lkr2,coreRichness(active_s4,core_fraction = 0.5)/coreRichness(perm_s4,core_fraction = 0.5))


mkf<-c()
mke<-c()
mkr<-c()

#branch-based ASV
mkf<-c(mkf,coreFaithsPD(bog,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm,core_fraction = 0.5,ab_threshold1=1))
mkf<-c(mkf,coreFaithsPD(bog_s1,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s1,core_fraction = 0.5,ab_threshold1=1))
mkf<-c(mkf,coreFaithsPD(bog_s2,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s2,core_fraction = 0.5,ab_threshold1=1))
mkf<-c(mkf,coreFaithsPD(bog_s3,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s3,core_fraction = 0.5,ab_threshold1=1))
mkf<-c(mkf,coreFaithsPD(bog_s4,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s4,core_fraction = 0.5,ab_threshold1=1))

#tip-based ASV
mke<-c(mke,coreFaithsPD(bog,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm,mode='tip',core_fraction = 0.5,ab_threshold1=1))
mke<-c(mke,coreFaithsPD(bog_s1,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s1,mode='tip',core_fraction = 0.5,ab_threshold1=1))
mke<-c(mke,coreFaithsPD(bog_s2,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s2,mode='tip',core_fraction = 0.5,ab_threshold1=1))
mke<-c(mke,coreFaithsPD(bog_s3,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s3,mode='tip',core_fraction = 0.5,ab_threshold1=1))
mke<-c(mke,coreFaithsPD(bog_s4,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s4,mode='tip',core_fraction = 0.5,ab_threshold1=1))

#non-phylogenetic ASV
mkr<-c(mkr,coreRichness(bog,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm,core_fraction = 0.5,ab_threshold1=1))
mkr<-c(mkr,coreRichness(bog_s1,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm_s1,core_fraction = 0.5,ab_threshold1=1))
mkr<-c(mkr,coreRichness(bog_s2,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm_s2,core_fraction = 0.5,ab_threshold1=1))
mkr<-c(mkr,coreRichness(bog_s3,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm_s3,core_fraction = 0.5,ab_threshold1=1))
mkr<-c(mkr,coreRichness(bog_s4,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm_s4,core_fraction = 0.5,ab_threshold1=1))

mkf2<-c()
mke2<-c()
mkr2<-c()

#branch-based ASV
mkf2<-c(mkf2,coreFaithsPD(active,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm,core_fraction = 0.5,ab_threshold1=1))
mkf2<-c(mkf2,coreFaithsPD(active_s1,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s1,core_fraction = 0.5,ab_threshold1=1))
mkf2<-c(mkf2,coreFaithsPD(active_s2,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s2,core_fraction = 0.5,ab_threshold1=1))
mkf2<-c(mkf2,coreFaithsPD(active_s3,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s3,core_fraction = 0.5,ab_threshold1=1))
mkf2<-c(mkf2,coreFaithsPD(active_s4,core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s4,core_fraction = 0.5,ab_threshold1=1))

#tip-based ASV
mke2<-c(mke2,coreFaithsPD(active,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm,mode='tip',core_fraction = 0.5,ab_threshold1=1))
mke2<-c(mke2,coreFaithsPD(active_s1,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s1,mode='tip',core_fraction = 0.5,ab_threshold1=1))
mke2<-c(mke2,coreFaithsPD(active_s2,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s2,mode='tip',core_fraction = 0.5,ab_threshold1=1))
mke2<-c(mke2,coreFaithsPD(active_s3,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s3,mode='tip',core_fraction = 0.5,ab_threshold1=1))
mke2<-c(mke2,coreFaithsPD(active_s4,mode='tip',core_fraction = 0.5,ab_threshold1=1)/coreFaithsPD(perm_s4,mode='tip',core_fraction = 0.5,ab_threshold1=1))

#non-phylogenetic ASV
mkr2<-c(mkr2,coreRichness(active,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm,core_fraction = 0.5,ab_threshold1=1))
mkr2<-c(mkr2,coreRichness(active_s1,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm_s1,core_fraction = 0.5,ab_threshold1=1))
mkr2<-c(mkr2,coreRichness(active_s2,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm_s2,core_fraction = 0.5,ab_threshold1=1))
mkr2<-c(mkr2,coreRichness(active_s3,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm_s3,core_fraction = 0.5,ab_threshold1=1))
mkr2<-c(mkr2,coreRichness(active_s4,core_fraction = 0.5,ab_threshold1=1)/coreRichness(perm_s4,core_fraction = 0.5,ab_threshold1=1))

nkf<-c()
nke<-c()
nkr<-c()

#branch-based ASV
nkf<-c(nkf,coreFaithsPD(bog,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm,core_fraction = 0.5,selection='shade',initial_branches=10))
nkf<-c(nkf,coreFaithsPD(bog_s1,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm_s1,core_fraction = 0.5,selection='shade',initial_branches=10))
nkf<-c(nkf,coreFaithsPD(bog_s2,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm_s2,core_fraction = 0.5,selection='shade',initial_branches=10))
nkf<-c(nkf,coreFaithsPD(bog_s3,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm_s3,core_fraction = 0.5,selection='shade',initial_branches=10))
nkf<-c(nkf,coreFaithsPD(bog_s4,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm_s4,core_fraction = 0.5,selection='shade',initial_branches=10))

#tip-based ASV
nke<-c(nke,coreFaithsPD(bog,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm,mode='tip',core_fraction = 0.5,selection='shade'))
nke<-c(nke,coreFaithsPD(bog_s1,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm_s1,mode='tip',core_fraction = 0.5,selection='shade'))
nke<-c(nke,coreFaithsPD(bog_s2,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm_s2,mode='tip',core_fraction = 0.5,selection='shade'))
nke<-c(nke,coreFaithsPD(bog_s3,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm_s3,mode='tip',core_fraction = 0.5,selection='shade'))
nke<-c(nke,coreFaithsPD(bog_s4,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm_s4,mode='tip',core_fraction = 0.5,selection='shade'))

#non-phylogenetic ASV
nkr<-c(nkr,coreRichness(bog,core_fraction = 0.5,selection='shade')/coreRichness(perm,core_fraction = 0.5,selection='shade'))
nkr<-c(nkr,coreRichness(bog_s1,core_fraction = 0.5,selection='shade')/coreRichness(perm_s1,core_fraction = 0.5,selection='shade'))
nkr<-c(nkr,coreRichness(bog_s2,core_fraction = 0.5,selection='shade')/coreRichness(perm_s2,core_fraction = 0.5,selection='shade'))
nkr<-c(nkr,coreRichness(bog_s3,core_fraction = 0.5,selection='shade')/coreRichness(perm_s3,core_fraction = 0.5,selection='shade'))
nkr<-c(nkr,coreRichness(bog_s4,core_fraction = 0.5,selection='shade')/coreRichness(perm_s4,core_fraction = 0.5,selection='shade'))

nkf2<-c()
nke2<-c()
nkr2<-c()

#branch-based ASV
nkf2<-c(nkf2,coreFaithsPD(active,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm,core_fraction = 0.5,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreFaithsPD(active_s1,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm_s1,core_fraction = 0.5,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreFaithsPD(active_s2,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm_s2,core_fraction = 0.5,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreFaithsPD(active_s3,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm_s3,core_fraction = 0.5,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreFaithsPD(active_s4,core_fraction = 0.5,selection='shade',initial_branches=10)/coreFaithsPD(perm_s4,core_fraction = 0.5,selection='shade',initial_branches=10))

#tip-based ASV
nke2<-c(nke2,coreFaithsPD(active,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm,mode='tip',core_fraction = 0.5,selection='shade'))
nke2<-c(nke2,coreFaithsPD(active_s1,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm_s1,mode='tip',core_fraction = 0.5,selection='shade'))
nke2<-c(nke2,coreFaithsPD(active_s2,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm_s2,mode='tip',core_fraction = 0.5,selection='shade'))
nke2<-c(nke2,coreFaithsPD(active_s3,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm_s3,mode='tip',core_fraction = 0.5,selection='shade'))
nke2<-c(nke2,coreFaithsPD(active_s4,mode='tip',core_fraction = 0.5,selection='shade')/coreFaithsPD(perm_s4,mode='tip',core_fraction = 0.5,selection='shade'))

#non-phylogenetic ASV
nkr2<-c(nkr2,coreRichness(active,core_fraction = 0.5,selection='shade')/coreRichness(perm,core_fraction = 0.5,selection='shade'))
nkr2<-c(nkr2,coreRichness(active_s1,core_fraction = 0.5,selection='shade')/coreRichness(perm_s1,core_fraction = 0.5,selection='shade'))
nkr2<-c(nkr2,coreRichness(active_s2,core_fraction = 0.5,selection='shade')/coreRichness(perm_s2,core_fraction = 0.5,selection='shade'))
nkr2<-c(nkr2,coreRichness(active_s3,core_fraction = 0.5,selection='shade')/coreRichness(perm_s3,core_fraction = 0.5,selection='shade'))
nkr2<-c(nkr2,coreRichness(active_s4,core_fraction = 0.5,selection='shade')/coreRichness(perm_s4,core_fraction = 0.5,selection='shade'))



#Plot variation in Faith's PD across taxonoperm ranks
axisv<-c(1,2,3,4,5)

par(mar=c(5,5,5,5))
plot(axisv,lkf,col='brown',ylim=c(1,1.75),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
axis(1, at=1:6, labels=c('0','0.05','0.1','0.15','0.2','0.25'))
points(axisv,lkf,col='black',pch=22,cex=1.5)
lines(axisv,lkf,col='brown',lty='solid',lwd=2)
points(axisv,lkr,col='brown',pch=16,cex=1.5)
points(axisv,lkr,col='black',pch=21,cex=1.5)
lines(axisv,lkr,col='brown',lty='dashed',lwd=2)
points(axisv,lke,col='brown',pch=17,cex=1.5)
points(axisv,lke,col='black',pch=24,cex=1.5)
lines(axisv,lke,col='brown',lty='dotted',lwd=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

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
plot(axisv,mkf,col='tan3',ylim=c(0.6,1.3),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
axis(1, at=1:6, labels=c('0','0.05','0.1','0.15','0.2','0.25'))
points(axisv,mkf,col='black',pch=22,cex=1.5)
lines(axisv,mkf,col='tan3',lty='solid',lwd=2)
points(axisv,mkr,col='tan3',pch=16,cex=1.5)
points(axisv,mkr,col='black',pch=21,cex=1.5)
lines(axisv,mkr,col='tan3',lty='dashed',lwd=2)
points(axisv,mke,col='tan3',pch=17,cex=1.5)
points(axisv,mke,col='black',pch=24,cex=1.5)
lines(axisv,mke,col='tan3',lty='dotted',lwd=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

points(axisv,mkf2,col='palegreen2',pch=15,cex=1.5)
points(axisv,mkf2,col='black',pch=22,cex=1.5)
lines(axisv,mkf2,col='palegreen2',lty='solid',lwd=2)
points(axisv,mkr2,col='palegreen2',pch=16,cex=1.5)
points(axisv,mkr2,col='black',pch=21,cex=1.5)
lines(axisv,mkr2,col='palegreen2',lty='dashed',lwd=2)
points(axisv,mke2,col='palegreen2',pch=17,cex=1.5)
points(axisv,mke2,col='black',pch=24,cex=1.5)
lines(axisv,mke2,col='palegreen2',lty='dotted',lwd=2)

plot(axisv,nkf,col='burlywood2',ylim=c(0.75,1.75),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
axis(1, at=1:6, labels=c('0','0.05','0.1','0.15','0.2','0.25'))
points(axisv,nkf,col='black',pch=22,cex=1.5)
lines(axisv,nkf,col='burlywood2',lty='solid',lwd=2)
points(axisv,nkr,col='burlywood2',pch=16,cex=1.5)
points(axisv,nkr,col='black',pch=21,cex=1.5)
lines(axisv,nkr,col='burlywood2',lty='dashed',lwd=2)
points(axisv,nke,col='burlywood2',pch=17,cex=1.5)
points(axisv,nke,col='black',pch=24,cex=1.5)
lines(axisv,nke,col='burlywood2',lty='dotted',lwd=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

points(axisv,nkf2,col='darkseagreen1',pch=15,cex=1.5)
points(axisv,nkf2,col='black',pch=22,cex=1.5)
lines(axisv,nkf2,col='darkseagreen1',lty='solid',lwd=2)
points(axisv,nkr2,col='darkseagreen1',pch=16,cex=1.5)
points(axisv,nkr2,col='black',pch=21,cex=1.5)
lines(axisv,nkr2,col='darkseagreen1',lty='dashed',lwd=2)
points(axisv,nke2,col='darkseagreen1',pch=17,cex=1.5)
points(axisv,nke2,col='black',pch=24,cex=1.5)
lines(axisv,nke2,col='darkseagreen1',lty='dotted',lwd=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

#Find maximum and minimum across ranks for each metric
x = 1:18
High = c(max(lkf),max(lke),max(lkr),max(mkf),max(mke),max(mkr),max(nkf),max(nke),max(nkr),max(lkf2),max(lke2),max(lkr2),max(mkf2),max(mke2),max(mkr2),max(nkf2),max(nke2),max(nkr2))
Low = c(min(lkf),min(lke),min(lkr),min(mkf),min(mke),min(mkr),min(nkf),min(nke),min(nkr),min(lkf2),min(lke2),min(lkr2),min(mkf2),min(mke2),min(mkr2),min(nkf2),min(nke2),min(nkr2))

## Plot variation across ranks
par(mar=c(10,5,1,2))
plot(1, type="n", xlab="", ylab="Relative Diversity Range", xlim=c(0.5,18.5),
     ylim=c(0.65, 1.75),xaxt='n',cex.lab=1.25)

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
axis(1, at=1:18, labels=c('occ-BranchPD','occ-TipPD','occ-Richness','occ-ab-BranchPD','occ-ab-TipPD','occ-ab-Richness','shade-BranchPD','shade-TipPD','shade-Richness','occ-BranchPD','occ-TipPD','occ-Richness','occ-ab-BranchPD','occ-ab-TipPD','occ-ab-Richness','shade-BranchPD','shade-TipPD','shade-Richness'),las=2)
axis(2,cex=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

