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

#Pool to species
Sup_s1<-tip_glom(Sup,h=0.05)
Mic_s1<-tip_glom(Mic,h=0.05)

Sup_s2<-tip_glom(Sup,h=0.1)
Mic_s2<-tip_glom(Mic,h=0.1)

Sup_s3<-tip_glom(Sup,h=0.15)
Mic_s3<-tip_glom(Mic,h=0.15)

Sup_s4<-tip_glom(Sup,h=0.2)
Mic_s4<-tip_glom(Mic,h=0.2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate variation in Faith's PD across taxonomic rank

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lkf<-c()
lke<-c()
lkr<-c()

#branch-based ASV
lkf<-c(lkf,coreFaithsPD(Sup,core_fraction = 0.5)/coreFaithsPD(Mic,core_fraction = 0.5))
lkf<-c(lkf,coreFaithsPD(Sup_s1,core_fraction = 0.5)/coreFaithsPD(Mic_s1,core_fraction = 0.5))
lkf<-c(lkf,coreFaithsPD(Sup_s2,core_fraction = 0.5)/coreFaithsPD(Mic_s2,core_fraction = 0.5))
lkf<-c(lkf,coreFaithsPD(Sup_s3,core_fraction = 0.5)/coreFaithsPD(Mic_s3,core_fraction = 0.5))
lkf<-c(lkf,coreFaithsPD(Sup_s4,core_fraction = 0.5)/coreFaithsPD(Mic_s4,core_fraction = 0.5))

#tip-based ASV
lke<-c(lke,coreFaithsPD(Sup,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic,mode='tip',core_fraction = 0.5))
lke<-c(lke,coreFaithsPD(Sup_s1,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic_s1,mode='tip',core_fraction = 0.5))
lke<-c(lke,coreFaithsPD(Sup_s2,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic_s2,mode='tip',core_fraction = 0.5))
lke<-c(lke,coreFaithsPD(Sup_s3,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic_s3,mode='tip',core_fraction = 0.5))
lke<-c(lke,coreFaithsPD(Sup_s4,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic_s4,mode='tip',core_fraction = 0.5))

#non-phylogenetic ASV
lkr<-c(lkr,coreRichness(Sup,core_fraction = 0.5)/coreRichness(Mic,core_fraction = 0.5))
lkr<-c(lkr,coreRichness(Sup_s1,core_fraction = 0.5)/coreRichness(Mic_s1,core_fraction = 0.5))
lkr<-c(lkr,coreRichness(Sup_s2,core_fraction = 0.5)/coreRichness(Mic_s2,core_fraction = 0.5))
lkr<-c(lkr,coreRichness(Sup_s3,core_fraction = 0.5)/coreRichness(Mic_s3,core_fraction = 0.5))
lkr<-c(lkr,coreRichness(Sup_s4,core_fraction = 0.5)/coreRichness(Mic_s4,core_fraction = 0.5))

lkf2<-c()
lke2<-c()
lkr2<-c()

#branch-based ASV
lkf2<-c(lkf2,coreFaithsPD(Sup,ab_threshold1=1,core_fraction = 0.5)/coreFaithsPD(Mic,ab_threshold1=1,core_fraction = 0.5))
lkf2<-c(lkf2,coreFaithsPD(Sup_s1,ab_threshold1=1,core_fraction = 0.5)/coreFaithsPD(Mic_s1,ab_threshold1=1,core_fraction = 0.5))
lkf2<-c(lkf2,coreFaithsPD(Sup_s2,ab_threshold1=1,core_fraction = 0.5)/coreFaithsPD(Mic_s2,ab_threshold1=1,core_fraction = 0.5))
lkf2<-c(lkf2,coreFaithsPD(Sup_s3,ab_threshold1=1,core_fraction = 0.5)/coreFaithsPD(Mic_s3,ab_threshold1=1,core_fraction = 0.5))
lkf2<-c(lkf2,coreFaithsPD(Sup_s4,ab_threshold1=1,core_fraction = 0.5)/coreFaithsPD(Mic_s4,ab_threshold1=1,core_fraction = 0.5))

#tip-based ASV
lke2<-c(lke2,coreFaithsPD(Sup,ab_threshold1=1,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic,ab_threshold1=1,mode='tip',core_fraction = 0.5))
lke2<-c(lke2,coreFaithsPD(Sup_s1,ab_threshold1=1,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic_s1,ab_threshold1=1,mode='tip',core_fraction = 0.5))
lke2<-c(lke2,coreFaithsPD(Sup_s2,ab_threshold1=1,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic_s2,ab_threshold1=1,mode='tip',core_fraction = 0.5))
lke2<-c(lke2,coreFaithsPD(Sup_s3,ab_threshold1=1,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic_s3,ab_threshold1=1,mode='tip',core_fraction = 0.5))
lke2<-c(lke2,coreFaithsPD(Sup_s4,ab_threshold1=1,mode='tip',core_fraction = 0.5)/coreFaithsPD(Mic_s4,ab_threshold1=1,mode='tip',core_fraction = 0.5))

#non-phylogenetic ASV
lkr2<-c(lkr2,coreRichness(Sup,ab_threshold1=1,core_fraction = 0.5)/coreRichness(Mic,ab_threshold1=1,core_fraction = 0.5))
lkr2<-c(lkr2,coreRichness(Sup_s1,ab_threshold1=1,core_fraction = 0.5)/coreRichness(Mic_s1,ab_threshold1=1,core_fraction = 0.5))
lkr2<-c(lkr2,coreRichness(Sup_s2,ab_threshold1=1,core_fraction = 0.5)/coreRichness(Mic_s2,ab_threshold1=1,core_fraction = 0.5))
lkr2<-c(lkr2,coreRichness(Sup_s3,ab_threshold1=1,core_fraction = 0.5)/coreRichness(Mic_s3,ab_threshold1=1,core_fraction = 0.5))
lkr2<-c(lkr2,coreRichness(Sup_s4,ab_threshold1=1,core_fraction = 0.5)/coreRichness(Mic_s4,ab_threshold1=1,core_fraction = 0.5))


lkf3<-c()
lke3<-c()
lkr3<-c()

#branch-based ASV
lkf3<-c(lkf3,coreFaithsPD(Sup,selection='shade',initial_branches=10)/coreFaithsPD(Mic,selection='shade',initial_branches=10))
lkf3<-c(lkf3,coreFaithsPD(Sup_s1,selection='shade',initial_branches=10)/coreFaithsPD(Mic_s1,selection='shade',initial_branches=10))
lkf3<-c(lkf3,coreFaithsPD(Sup_s2,selection='shade',initial_branches=10)/coreFaithsPD(Mic_s2,selection='shade',initial_branches=10))
lkf3<-c(lkf3,coreFaithsPD(Sup_s3,selection='shade',initial_branches=10)/coreFaithsPD(Mic_s3,selection='shade',initial_branches=10))
lkf3<-c(lkf3,coreFaithsPD(Sup_s4,selection='shade',initial_branches=10)/coreFaithsPD(Mic_s4,selection='shade',initial_branches=10))

#tip-based ASV
lke3<-c(lke3,coreFaithsPD(Sup,selection='shade',mode='tip')/coreFaithsPD(Mic,selection='shade',mode='tip'))
lke3<-c(lke3,coreFaithsPD(Sup_s1,selection='shade',mode='tip')/coreFaithsPD(Mic_s1,selection='shade',mode='tip'))
lke3<-c(lke3,coreFaithsPD(Sup_s2,selection='shade',mode='tip')/coreFaithsPD(Mic_s2,selection='shade',mode='tip'))
lke3<-c(lke3,coreFaithsPD(Sup_s3,selection='shade',mode='tip')/coreFaithsPD(Mic_s3,selection='shade',mode='tip'))
lke3<-c(lke3,coreFaithsPD(Sup_s4,selection='shade',mode='tip')/coreFaithsPD(Mic_s4,selection='shade',mode='tip'))

#non-phylogenetic ASV
lkr3<-c(lkr3,coreRichness(Sup,selection='shade')/coreRichness(Mic,selection='shade'))
lkr3<-c(lkr3,coreRichness(Sup_s1,selection='shade')/coreRichness(Mic_s1,selection='shade'))
lkr3<-c(lkr3,coreRichness(Sup_s2,selection='shade')/coreRichness(Mic_s2,selection='shade'))
lkr3<-c(lkr3,coreRichness(Sup_s3,selection='shade')/coreRichness(Mic_s3,selection='shade'))
lkr3<-c(lkr3,coreRichness(Sup_s4,selection='shade')/coreRichness(Mic_s4,selection='shade'))



#Plot variation in Faith's PD across taxonomic ranks
axisv<-c(1,2,3,4,5)

par(mar=c(5,5,5,5))
plot(axisv,lkf,col='green',ylim=c(0.9,1.5),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
axis(1, at=1:6, labels=c('0','0.05','0.1','0.15','0.2','0.25'))
points(axisv,lkf,col='black',pch=22,cex=1.5)
lines(axisv,lkf,col='green',lty='solid',lwd=2)
points(axisv,lkr,col='green',pch=16,cex=1.5)
points(axisv,lkr,col='black',pch=21,cex=1.5)
lines(axisv,lkr,col='green',lty='dashed',lwd=2)
points(axisv,lke,col='green',pch=17,cex=1.5)
points(axisv,lke,col='black',pch=24,cex=1.5)
lines(axisv,lke,col='green',lty='dotted',lwd=2)

plot(axisv,lkf2,col='green',ylim=c(0.9,1.5),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
points(axisv,lkf2,col='black',pch=22,cex=1.5)
lines(axisv,lkf2,col='palegreen2',lty='solid',lwd=2)
points(axisv,lkr2,col='palegreen2',pch=16,cex=1.5)
points(axisv,lkr2,col='black',pch=21,cex=1.5)
lines(axisv,lkr2,col='palegreen2',lty='dashed',lwd=2)
points(axisv,lke2,col='palegreen2',pch=17,cex=1.5)
points(axisv,lke2,col='black',pch=24,cex=1.5)
lines(axisv,lke2,col='palegreen2',lty='dotted',lwd=2)

plot(axisv,lkf3,col='green',ylim=c(0.9,1.2),pch=15,cex=1.5,xlab='Tip Agglomeration',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
points(axisv,lkf3,col='black',pch=22,cex=1.5)
lines(axisv,lkf3,col='darkseagreen1',lty='solid',lwd=2)
points(axisv,lkr3,col='darkseagreen1',pch=16,cex=1.5)
points(axisv,lkr3,col='black',pch=21,cex=1.5)
lines(axisv,lkr3,col='darkseagreen1',lty='dashed',lwd=2)
points(axisv,lke3,col='darkseagreen1',pch=17,cex=1.5)
points(axisv,lke3,col='black',pch=24,cex=1.5)
lines(axisv,lke3,col='darkseagreen1',lty='dotted',lwd=2)


#Find maximum and minimum across ranks for each metric
x = 1:9
High = c(max(lkf),max(lke),max(lkr),max(lkf2),max(lke2),max(lkr2),max(lkf3),max(lke3),max(lkr3))
Low = c(min(lkf),min(lke),min(lkr),min(lkf2),min(lke2),min(lkr2),min(lkf3),min(lke3),min(lkr3))

## Plot variation across ranks
par(mar=c(10,5,2,2))
plot(1, type="n", xlab="", ylab="Relative Diversity Range", xlim=c(0.5,9.5),
     ylim=c(0.9, 1.5),xaxt='n',cex.lab=1.5)

## Add rectangles
rect(x[1:3] - 0.4, Low[1:3], x[1:3] + 0.4, High[1:3], col="green")
## Add rectangles
rect(x[4:6] - 0.4, Low[4:6], x[4:6] + 0.4, High[4:6], col="palegreen2")
## Add rectangles
rect(x[7:9] - 0.4, Low[7:9], x[7:9] + 0.4, High[7:9], col="darkseagreen1")
axis(1, at=1:9, labels=c('occ-BranchPD','occ-TipPD','occ-Richness','occ-ab-BranchPD','occ-ab-TipPD','occ-ab-Richness','shade-BranchPD','shade-TipPD','shade-Richness'),las=2)
axis(2,cex=2)
