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
  features="study_1036_051125-195547/BIOM/130932/feature-table.qza",
  "study_1036_051125-195547/BIOM/130932/taxonomy.qza",
  tree="study_1036_051125-195547/BIOM/130932/rooted-tree.qza",
  
  metadata = "1036_20230201-070300.txt")
Lphyloseq<-prune_taxa(taxa_sums(Lphyloseq)>0,Lphyloseq)
Lphyloseq<-rarefy_even_depth(Lphyloseq,rngseed=1)


a<-coreVenn(Lphyloseq,sample_data(Lphyloseq)$site,core_fraction=0.99,ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'))
b<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$site,core_fraction = 0.99,ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'),mode='tip')
c<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$site,core_fraction = 0.99,ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'))

d<-coreVenn(Lphyloseq,sample_data(Lphyloseq)$site,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'))
e<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$site,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'),mode='tip')
f<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$site,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'))

g<-coreVenn(Lphyloseq,sample_data(Lphyloseq)$site,selection='shade',ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'))
h<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$site,selection='shade',ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'),mode='tip')
i<-corePhyloVenn(Lphyloseq,sample_data(Lphyloseq)$site,selection='shade',initial_branches = 100,ordered_groups = c('bog','permafrost','active layer'),fill_color = c('brown','blue','green'))



