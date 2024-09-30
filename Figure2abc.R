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
library('holobiont')

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
tree_file<-multi2di(read.tree('rooted_tree_out/tree.nwk'))

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
assigned_taxa<-read.csv('ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
ZNC_physeq_all = rarefy_even_depth(prune_taxa(assigned_taxa_seqs, Z_physeq))

corePhyloVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$speciesxsite,0.5,decimal=2,ordered_groups=c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm'),fill_color=c('yellow','green','blue','red'),mode='branch',text_size=6,set_name_size=0)

corePhyloVenn(ZNC_physeq_all,sample_data(ZNC_physeq_all)$speciesxsite,0.5,decimal=2,ordered_groups=c('SBluG_ino','SBluG_neo','SNW_neo','SNW_marm'),fill_color=c('yellow','green','blue','red'),mode='tip',text_size=6,set_name_size=0)



ino_hit<-which(sample_data(ZNC_physeq_all)$speciesxsite=='SBluG_ino')
core_ino<-which(rowSums(sign(otu_table(ZNC_physeq_all))[,ino_hit])>0.5*length(ino_hit))

neo1_hit<-which(sample_data(ZNC_physeq_all)$speciesxsite=='SBluG_neo')
core_neo1<-which(rowSums(sign(otu_table(ZNC_physeq_all))[,neo1_hit])>0.5*length(neo1_hit))

neo2_hit<-which(sample_data(ZNC_physeq_all)$speciesxsite=='SNW_neo')
core_neo2<-which(rowSums(sign(otu_table(ZNC_physeq_all))[,neo2_hit])>0.5*length(neo2_hit))

marm_hit<-which(sample_data(ZNC_physeq_all)$speciesxsite=='SNW_marm')
core_marm<-which(rowSums(sign(otu_table(ZNC_physeq_all))[,marm_hit])>0.5*length(marm_hit))

c1<-rep(FALSE,length(taxa_names(ZNC_physeq_all)))
hits<-cbind(c1,c1,c1,c1)

hits[core_ino,1]<-TRUE
hits[core_neo1,2]<-TRUE
hits[core_neo2,3]<-TRUE
hits[core_marm,4]<-TRUE

hits2<-hits[rowSums(hits != FALSE) > 0,]
colnames(hits2)<-c('inornatus','neo1','neo2','marms')

ggvenn2(data.frame(hits2),set_name_size=0,text_size=6,fill_color = c('yellow','green','blue','red'))+coord_cartesian(clip="off")


