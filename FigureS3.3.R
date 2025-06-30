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

set.seed(10)

#Read in CSV file for OTU table (table.csv for ASVs, tablelevel6.csv for genera)
data <- import_biom('HybridLizardSkins.biom')

#Find ASV table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

#Read in metadata
meta_csv<-read.csv('Whiptail_Skin_Metadata.csv', row.names=1,header= TRUE)

#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)

#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-multi2di(read.tree('skin_rooted_tree_out/tree.nwk'))

#Rename the tree tips to match the ASV table names
cut_names<-c()
for (k in 1:length(taxa_names(tree_file))){
  #Depending on the Qiime2 output, you may need to split the names with a space or with an underscore
  #cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],' ')[[1]][1])
  cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],'_')[[1]][1])
}
taxa_names(tree_file)<-cut_names

#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq<-phyloseq(OTU_biom,TAX,SAMP,tree_file)

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('skin_ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
ZNC_physeq_all = rarefy_even_depth(prune_taxa(assigned_taxa_seqs, Z_physeq))
sample_data(ZNC_physeq_all)$species[sample_data(ZNC_physeq_all)$species=='inornatus']<-'A. arizonae'
sample_data(ZNC_physeq_all)$species[sample_data(ZNC_physeq_all)$species=='marmoratus']<-'A. marmoratus'
sample_data(ZNC_physeq_all)$species[sample_data(ZNC_physeq_all)$species=='neomexicanus']<-'A. neomexicanus'

a<-coreVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,core_fraction = 0.49,ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'))
b<-corePhyloVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,core_fraction = 0.49,ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'),mode='tip')
c<-corePhyloVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,core_fraction = 0.49,ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'))

d<-coreVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'))
e<-corePhyloVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'),mode='tip')
f<-corePhyloVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,core_fraction = 0.49,ab_threshold1 = 0.1,ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'))

g<-coreVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,selection='shade',ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'))
h<-corePhyloVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,selection='shade',ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'),mode='tip')
i<-corePhyloVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$species,selection='shade',initial_branches = 10,ordered_groups = c('A. marmoratus','A. arizonae','A. neomexicanus'),fill_color = c('red','yellow','orange'))



