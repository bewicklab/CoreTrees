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


a<-coreVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,core_fraction = 0.49,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),text_size = 8)
b<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,core_fraction = 0.49,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),mode='tip',text_size = 8)
c<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,core_fraction = 0.49,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),text_size = 8)

d<-coreVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),text_size = 8)
e<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),mode='tip',text_size = 8)
f<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),text_size = 8)

g<-coreVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,core_fraction = 0.49,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),selection='shade',text_size = 8)
h<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,core_fraction = 0.49,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),mode='tip',selection='shade',text_size = 8)
i<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$Lake,initial_branches=10,core_fraction = 0.49,ordered_groups = c('Michigan','Superior'),fill_color = c('blue','green'),selection='shade',text_size = 8)



