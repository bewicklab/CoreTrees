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


set.seed(10)

#Read in CSV file for OTU table (table.csv for ASVs, tablelevel6.csv for genera)
data <- import_biom('HybridLizardGuts.biom')

#Find ASV table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

#Read in metadata
meta_csv<-read.csv('Whiptail_Gut_Metadata.csv', row.names=1,header= TRUE)

#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)

#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-multi2di(read.tree('gut_rooted_tree_out/tree.nwk'))

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
assigned_taxa<-read.csv('gut_ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
ZNC_physeq_all = rarefy_even_depth(prune_taxa(assigned_taxa_seqs, Z_physeq))


#Find the treatment group for each of your samples
type<-sample_data(ZNC_physeq_all)$species
inornatus_list<-which(type=='inornatus')
marmoratus_list<-which(type=='marmoratus')
neomexicanus_list<-which(type=='neomexicanus')


Zino<-prune_samples(sample_names(ZNC_physeq_all)[inornatus_list],ZNC_physeq_all)
Zmarm<-prune_samples(sample_names(ZNC_physeq_all)[marmoratus_list],ZNC_physeq_all)
Zneo<-prune_samples(sample_names(ZNC_physeq_all)[neomexicanus_list],ZNC_physeq_all)

#Pool to species
Zino_s<-tax_glom(Zino,taxrank='Rank7')
Zmarm_s<-tax_glom(Zmarm,taxrank='Rank7')
Zneo_s<-tax_glom(Zneo,taxrank='Rank7')

#Pool to genus
Zino_g<-tax_glom(Zino,taxrank='Rank6')
Zmarm_g<-tax_glom(Zmarm,taxrank='Rank6')
Zneo_g<-tax_glom(Zneo,taxrank='Rank6')

#Pool to family
Zino_f<-tax_glom(Zino,taxrank='Rank5')
Zmarm_f<-tax_glom(Zmarm,taxrank='Rank5')
Zneo_f<-tax_glom(Zneo,taxrank='Rank5')

#Pool to order
Zino_o<-tax_glom(Zino,taxrank='Rank4')
Zmarm_o<-tax_glom(Zmarm,taxrank='Rank4')
Zneo_o<-tax_glom(Zneo,taxrank='Rank4')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate variation in Faith's PD across taxonoZneo rank

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lkf<-c()
lke<-c()
lkr<-c()

#branch-based ASV
lkf<-c(lkf,coreFaithsPD(Zino,core_fraction = 0.49)/coreFaithsPD(Zneo,core_fraction = 0.49))
lkf<-c(lkf,coreFaithsPD(Zino_s,core_fraction = 0.49)/coreFaithsPD(Zneo_s,core_fraction = 0.49))
lkf<-c(lkf,coreFaithsPD(Zino_g,core_fraction = 0.49)/coreFaithsPD(Zneo_g,core_fraction = 0.49))
lkf<-c(lkf,coreFaithsPD(Zino_f,core_fraction = 0.49)/coreFaithsPD(Zneo_f,core_fraction = 0.49))
lkf<-c(lkf,coreFaithsPD(Zino_o,core_fraction = 0.49)/coreFaithsPD(Zneo_o,core_fraction = 0.49))

#tip-based ASV
lke<-c(lke,coreFaithsPD(Zino,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo,mode='tip',core_fraction = 0.49))
lke<-c(lke,coreFaithsPD(Zino_s,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo_s,mode='tip',core_fraction = 0.49))
lke<-c(lke,coreFaithsPD(Zino_g,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo_g,mode='tip',core_fraction = 0.49))
lke<-c(lke,coreFaithsPD(Zino_f,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo_f,mode='tip',core_fraction = 0.49))
lke<-c(lke,coreFaithsPD(Zino_o,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo_o,mode='tip',core_fraction = 0.49))

#non-phylogenetic ASV
lkr<-c(lkr,coreRichness(Zino,core_fraction = 0.49)/coreRichness(Zneo,core_fraction = 0.49))
lkr<-c(lkr,coreRichness(Zino_s,core_fraction = 0.49)/coreRichness(Zneo_s,core_fraction = 0.49))
lkr<-c(lkr,coreRichness(Zino_g,core_fraction = 0.49)/coreRichness(Zneo_g,core_fraction = 0.49))
lkr<-c(lkr,coreRichness(Zino_f,core_fraction = 0.49)/coreRichness(Zneo_f,core_fraction = 0.49))
lkr<-c(lkr,coreRichness(Zino_o,core_fraction = 0.49)/coreRichness(Zneo_o,core_fraction = 0.49))

lkf2<-c()
lke2<-c()
lkr2<-c()

#branch-based ASV
lkf2<-c(lkf2,coreFaithsPD(Zmarm,core_fraction = 0.49)/coreFaithsPD(Zneo,core_fraction = 0.49))
lkf2<-c(lkf2,coreFaithsPD(Zmarm_s,core_fraction = 0.49)/coreFaithsPD(Zneo_s,core_fraction = 0.49))
lkf2<-c(lkf2,coreFaithsPD(Zmarm_g,core_fraction = 0.49)/coreFaithsPD(Zneo_g,core_fraction = 0.49))
lkf2<-c(lkf2,coreFaithsPD(Zmarm_f,core_fraction = 0.49)/coreFaithsPD(Zneo_f,core_fraction = 0.49))
lkf2<-c(lkf2,coreFaithsPD(Zmarm_o,core_fraction = 0.49)/coreFaithsPD(Zneo_o,core_fraction = 0.49))

#tip-based ASV
lke2<-c(lke2,coreFaithsPD(Zmarm,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo,mode='tip',core_fraction = 0.49))
lke2<-c(lke2,coreFaithsPD(Zmarm_s,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo_s,mode='tip',core_fraction = 0.49))
lke2<-c(lke2,coreFaithsPD(Zmarm_g,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo_g,mode='tip',core_fraction = 0.49))
lke2<-c(lke2,coreFaithsPD(Zmarm_f,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo_f,mode='tip',core_fraction = 0.49))
lke2<-c(lke2,coreFaithsPD(Zmarm_o,mode='tip',core_fraction = 0.49)/coreFaithsPD(Zneo_o,mode='tip',core_fraction = 0.49))

#non-phylogenetic ASV
lkr2<-c(lkr2,coreRichness(Zmarm,core_fraction = 0.49)/coreRichness(Zneo,core_fraction = 0.49))
lkr2<-c(lkr2,coreRichness(Zmarm_s,core_fraction = 0.49)/coreRichness(Zneo_s,core_fraction = 0.49))
lkr2<-c(lkr2,coreRichness(Zmarm_g,core_fraction = 0.49)/coreRichness(Zneo_g,core_fraction = 0.49))
lkr2<-c(lkr2,coreRichness(Zmarm_f,core_fraction = 0.49)/coreRichness(Zneo_f,core_fraction = 0.49))
lkr2<-c(lkr2,coreRichness(Zmarm_o,core_fraction = 0.49)/coreRichness(Zneo_o,core_fraction = 0.49))


mkf<-c()
mke<-c()
mkr<-c()

#branch-based ASV
mkf<-c(mkf,coreFaithsPD(Zino,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo,core_fraction = 0.49,ab_threshold1=0.1))
mkf<-c(mkf,coreFaithsPD(Zino_s,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_s,core_fraction = 0.49,ab_threshold1=0.1))
mkf<-c(mkf,coreFaithsPD(Zino_g,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_g,core_fraction = 0.49,ab_threshold1=0.1))
mkf<-c(mkf,coreFaithsPD(Zino_f,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_f,core_fraction = 0.49,ab_threshold1=0.1))
mkf<-c(mkf,coreFaithsPD(Zino_o,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_o,core_fraction = 0.49,ab_threshold1=0.1))

#tip-based ASV
mke<-c(mke,coreFaithsPD(Zino,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))
mke<-c(mke,coreFaithsPD(Zino_s,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_s,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))
mke<-c(mke,coreFaithsPD(Zino_g,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_g,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))
mke<-c(mke,coreFaithsPD(Zino_f,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_f,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))
mke<-c(mke,coreFaithsPD(Zino_o,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_o,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))

#non-phylogenetic ASV
mkr<-c(mkr,coreRichness(Zino,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo,core_fraction = 0.49,ab_threshold1=0.1))
mkr<-c(mkr,coreRichness(Zino_s,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo_s,core_fraction = 0.49,ab_threshold1=0.1))
mkr<-c(mkr,coreRichness(Zino_g,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo_g,core_fraction = 0.49,ab_threshold1=0.1))
mkr<-c(mkr,coreRichness(Zino_f,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo_f,core_fraction = 0.49,ab_threshold1=0.1))
mkr<-c(mkr,coreRichness(Zino_o,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo_o,core_fraction = 0.49,ab_threshold1=0.1))

mkf2<-c()
mke2<-c()
mkr2<-c()

#branch-based ASV
mkf2<-c(mkf2,coreFaithsPD(Zmarm,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo,core_fraction = 0.49,ab_threshold1=0.1))
mkf2<-c(mkf2,coreFaithsPD(Zmarm_s,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_s,core_fraction = 0.49,ab_threshold1=0.1))
mkf2<-c(mkf2,coreFaithsPD(Zmarm_g,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_g,core_fraction = 0.49,ab_threshold1=0.1))
mkf2<-c(mkf2,coreFaithsPD(Zmarm_f,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_f,core_fraction = 0.49,ab_threshold1=0.1))
mkf2<-c(mkf2,coreFaithsPD(Zmarm_o,core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_o,core_fraction = 0.49,ab_threshold1=0.1))

#tip-based ASV
mke2<-c(mke2,coreFaithsPD(Zmarm,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))
mke2<-c(mke2,coreFaithsPD(Zmarm_s,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_s,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))
mke2<-c(mke2,coreFaithsPD(Zmarm_g,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_g,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))
mke2<-c(mke2,coreFaithsPD(Zmarm_f,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_f,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))
mke2<-c(mke2,coreFaithsPD(Zmarm_o,mode='tip',core_fraction = 0.49,ab_threshold1=0.1)/coreFaithsPD(Zneo_o,mode='tip',core_fraction = 0.49,ab_threshold1=0.1))

#non-phylogenetic ASV
mkr2<-c(mkr2,coreRichness(Zmarm,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo,core_fraction = 0.49,ab_threshold1=0.1))
mkr2<-c(mkr2,coreRichness(Zmarm_s,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo_s,core_fraction = 0.49,ab_threshold1=0.1))
mkr2<-c(mkr2,coreRichness(Zmarm_g,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo_g,core_fraction = 0.49,ab_threshold1=0.1))
mkr2<-c(mkr2,coreRichness(Zmarm_f,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo_f,core_fraction = 0.49,ab_threshold1=0.1))
mkr2<-c(mkr2,coreRichness(Zmarm_o,core_fraction = 0.49,ab_threshold1=0.1)/coreRichness(Zneo_o,core_fraction = 0.49,ab_threshold1=0.1))

nkf<-c()
nke<-c()
nkr<-c()

#branch-based ASV
nkf<-c(nkf,coreFaithsPD(Zino,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo,core_fraction = 0.49,selection='shade',initial_branches=10))
nkf<-c(nkf,coreFaithsPD(Zino_s,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo_s,core_fraction = 0.49,selection='shade',initial_branches=10))
nkf<-c(nkf,coreFaithsPD(Zino_g,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo_g,core_fraction = 0.49,selection='shade',initial_branches=10))
nkf<-c(nkf,coreFaithsPD(Zino_f,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo_f,core_fraction = 0.49,selection='shade',initial_branches=10))
nkf<-c(nkf,coreFaithsPD(Zino_o,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo_o,core_fraction = 0.49,selection='shade',initial_branches=10))

#tip-based ASV
nke<-c(nke,coreFaithsPD(Zino,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo,mode='tip',core_fraction = 0.49,selection='shade'))
nke<-c(nke,coreFaithsPD(Zino_s,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo_s,mode='tip',core_fraction = 0.49,selection='shade'))
nke<-c(nke,coreFaithsPD(Zino_g,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo_g,mode='tip',core_fraction = 0.49,selection='shade'))
nke<-c(nke,coreFaithsPD(Zino_f,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo_f,mode='tip',core_fraction = 0.49,selection='shade'))
nke<-c(nke,coreFaithsPD(Zino_o,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo_o,mode='tip',core_fraction = 0.49,selection='shade'))

#non-phylogenetic ASV
nkr<-c(nkr,coreRichness(Zino,core_fraction = 0.49,selection='shade')/coreRichness(Zneo,core_fraction = 0.49,selection='shade'))
nkr<-c(nkr,coreRichness(Zino_s,core_fraction = 0.49,selection='shade')/coreRichness(Zneo_s,core_fraction = 0.49,selection='shade'))
nkr<-c(nkr,coreRichness(Zino_g,core_fraction = 0.49,selection='shade')/coreRichness(Zneo_g,core_fraction = 0.49,selection='shade'))
nkr<-c(nkr,coreRichness(Zino_f,core_fraction = 0.49,selection='shade')/coreRichness(Zneo_f,core_fraction = 0.49,selection='shade'))
nkr<-c(nkr,coreRichness(Zino_o,core_fraction = 0.49,selection='shade')/coreRichness(Zneo_o,core_fraction = 0.49,selection='shade'))

nkf2<-c()
nke2<-c()
nkr2<-c()

#branch-based ASV
nkf2<-c(nkf2,coreFaithsPD(Zmarm,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo,core_fraction = 0.49,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreFaithsPD(Zmarm_s,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo_s,core_fraction = 0.49,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreFaithsPD(Zmarm_g,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo_g,core_fraction = 0.49,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreFaithsPD(Zmarm_f,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo_f,core_fraction = 0.49,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreFaithsPD(Zmarm_o,core_fraction = 0.49,selection='shade',initial_branches=10)/coreFaithsPD(Zneo_o,core_fraction = 0.49,selection='shade',initial_branches=10))

#tip-based ASV
nke2<-c(nke2,coreFaithsPD(Zmarm,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo,mode='tip',core_fraction = 0.49,selection='shade'))
nke2<-c(nke2,coreFaithsPD(Zmarm_s,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo_s,mode='tip',core_fraction = 0.49,selection='shade'))
nke2<-c(nke2,coreFaithsPD(Zmarm_g,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo_g,mode='tip',core_fraction = 0.49,selection='shade'))
nke2<-c(nke2,coreFaithsPD(Zmarm_f,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo_f,mode='tip',core_fraction = 0.49,selection='shade'))
nke2<-c(nke2,coreFaithsPD(Zmarm_o,mode='tip',core_fraction = 0.49,selection='shade')/coreFaithsPD(Zneo_o,mode='tip',core_fraction = 0.49,selection='shade'))

#non-phylogenetic ASV
nkr2<-c(nkr2,coreRichness(Zmarm,core_fraction = 0.49,selection='shade')/coreRichness(Zneo,core_fraction = 0.49,selection='shade'))
nkr2<-c(nkr2,coreRichness(Zmarm_s,core_fraction = 0.49,selection='shade')/coreRichness(Zneo_s,core_fraction = 0.49,selection='shade'))
nkr2<-c(nkr2,coreRichness(Zmarm_g,core_fraction = 0.49,selection='shade')/coreRichness(Zneo_g,core_fraction = 0.49,selection='shade'))
nkr2<-c(nkr2,coreRichness(Zmarm_f,core_fraction = 0.49,selection='shade')/coreRichness(Zneo_f,core_fraction = 0.49,selection='shade'))
nkr2<-c(nkr2,coreRichness(Zmarm_o,core_fraction = 0.49,selection='shade')/coreRichness(Zneo_o,core_fraction = 0.49,selection='shade'))



#Plot variation in Faith's PD across taxonoZneo ranks
axisv<-c(1,2,3,4,5)

par(mar=c(5,5,5,5))
plot(axisv,lkf,col='goldenrod1',ylim=c(0,3),pch=15,cex=1.5,xlab='Taxonomic Rank',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('ASV','species','genus','family','order'))
points(axisv,lkf,col='black',pch=22,cex=1.5)
lines(axisv,lkf,col='goldenrod1',lty='solid',lwd=2)
points(axisv,lkr,col='goldenrod1',pch=16,cex=1.5)
points(axisv,lkr,col='black',pch=21,cex=1.5)
lines(axisv,lkr,col='goldenrod1',lty='dashed',lwd=2)
points(axisv,lke,col='goldenrod1',pch=17,cex=1.5)
points(axisv,lke,col='black',pch=24,cex=1.5)
lines(axisv,lke,col='goldenrod1',lty='dotted',lwd=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

points(axisv,lkf2,col='red',pch=15,cex=1.5)
points(axisv,lkf2,col='black',pch=22,cex=1.5)
lines(axisv,lkf2,col='red',lty='solid',lwd=2)
points(axisv,lkr2,col='red',pch=16,cex=1.5)
points(axisv,lkr2,col='black',pch=21,cex=1.5)
lines(axisv,lkr2,col='red',lty='dashed',lwd=2)
points(axisv,lke2,col='red',pch=17,cex=1.5)
points(axisv,lke2,col='black',pch=24,cex=1.5)
lines(axisv,lke2,col='red',lty='dotted',lwd=2)

#legend(x='topright',legend=c("m-occ-Richness", "m-occ-TipPD","m-occ-BranchPD","a-occ-Richness", "a-occ-TipPD","a-occ-BranchPD","m-occ-ab-Richness", "m-occ-ab-TipPD","m-occ-ab-BranchPD","a-occ-ab-Richness", "a-occ-ab-TipPD","a-occ-ab-BranchPD","m-shade-Richness", "m-shade-TipPD","m-shade-BranchPD","a-shade-Richness", "a-shade-TipPD","a-shade-BranchPD"), 
#        pch = c(rep(c(16,17,15),6)), lty=c(rep(c('dashed','dotted','solid'),6)),col=c('red','red','red','goldenrod1','goldenrod1','goldenrod1','palevioletred','palevioletred','palevioletred','yellow','yellow','yellow','pink','pink','pink','khaki','khaki','khaki'))

#legend(x='topright',legend=c("b-occ-Richness", "b-occ-TipPD","b-occ-BranchPD","a-occ-Richness", "a-occ-TipPD","a-occ-BranchPD","b-occ-ab-Richness", "b-occ-ab-TipPD","b-occ-ab-BranchPD","a-occ-ab-Richness", "a-occ-ab-TipPD","a-occ-ab-BranchPD","b-shade-Richness", "b-shade-TipPD","b-shade-BranchPD","a-shade-Richness", "a-shade-TipPD","a-shade-BranchPD"), 
#       pch = c(rep(c(16,17,15),6)), lty=c(rep(c('dashed','dotted','solid'),6)),col=c('brown','brown','brown','green','green','green','tan3','tan3','tan3','palegreen2','palegreen2','palegreen2','burlywood2','burlywood2','burlywood2','darkseagreen1','darkseagreen1','darkseagreen1'))

#legend(x='topright',legend=c("m-occ-Jaccard", "m-occ-TipUniFrac","m-occ-BranchUniFrac","a-occ-Jaccard", "a-occ-TipUniFrac","a-occ-BranchUniFrac","m-occ-ab-Jaccard", "m-occ-ab-TipUniFrac","m-occ-ab-BranchUniFrac","a-occ-ab-Jaccard", "a-occ-ab-TipUniFrac","a-occ-ab-BranchUniFrac","m-shade-Jaccard", "m-shade-TipUniFrac","m-shade-BranchUniFrac","a-shade-Jaccard", "a-shade-TipUniFrac","a-shade-BranchUniFrac"), 
#        pch = c(rep(c(16,17,15),6)), lty=c(rep(c('dashed','dotted','solid'),6)),col=c('red','red','red','goldenrod1','goldenrod1','goldenrod1','palevioletred','palevioletred','palevioletred','yellow','yellow','yellow','pink','pink','pink','khaki','khaki','khaki'))

#legend(x='topright',legend=c("b-occ-Jaccard", "b-occ-TipUniFrac","b-occ-BranchUniFrac","a-occ-Jaccard", "a-occ-TipUniFrac","a-occ-BranchUniFrac","b-occ-ab-Jaccard", "b-occ-ab-TipUniFrac","b-occ-ab-BranchUniFrac","a-occ-ab-Jaccard", "a-occ-ab-TipUniFrac","a-occ-ab-BranchUniFrac","b-shade-Jaccard", "b-shade-TipUniFrac","b-shade-BranchUniFrac","a-shade-Jaccard", "a-shade-TipUniFrac","a-shade-BranchUniFrac"), 
#       pch = c(rep(c(16,17,15),6)), lty=c(rep(c('dashed','dotted','solid'),6)),col=c('brown','brown','brown','green','green','green','tan3','tan3','tan3','palegreen2','palegreen2','palegreen2','burlywood2','burlywood2','burlywood2','darkseagreen1','darkseagreen1','darkseagreen1'))

legend(x='topright',legend=c("occ-Jaccard", "occ-TipUniFrac","occ-BranchUniFrac","occ-ab-Jaccard", "occ-ab-TipUniFrac","occ-ab-BranchUniFrac","shade-Jaccard", "shade-TipUniFrac","shade-BranchUniFrac"), 
       pch = c(rep(c(16,17,15),6)), lty=c(rep(c('dashed','dotted','solid'),3)),col=c('green','green','green','palegreen2','palegreen2','palegreen2','darkseagreen1','darkseagreen1','darkseagreen1'))

par(mar=c(5,5,5,5))
plot(axisv,mkf,col='yellow',ylim=c(0.25,3),pch=15,cex=1.5,xlab='Taxonomic Rank',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('ASV','species','genus','family','order'))
points(axisv,mkf,col='black',pch=22,cex=1.5)
lines(axisv,mkf,col='yellow',lty='solid',lwd=2)
points(axisv,mkr,col='yellow',pch=16,cex=1.5)
points(axisv,mkr,col='black',pch=21,cex=1.5)
lines(axisv,mkr,col='yellow',lty='dashed',lwd=2)
points(axisv,mke,col='yellow',pch=17,cex=1.5)
points(axisv,mke,col='black',pch=24,cex=1.5)
lines(axisv,mke,col='yellow',lty='dotted',lwd=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

points(axisv,mkf2,col='palevioletred',pch=15,cex=1.5)
points(axisv,mkf2,col='black',pch=22,cex=1.5)
lines(axisv,mkf2,col='palevioletred',lty='solid',lwd=2)
points(axisv,mkr2,col='palevioletred',pch=16,cex=1.5)
points(axisv,mkr2,col='black',pch=21,cex=1.5)
lines(axisv,mkr2,col='palevioletred',lty='dashed',lwd=2)
points(axisv,mke2,col='palevioletred',pch=17,cex=1.5)
points(axisv,mke2,col='black',pch=24,cex=1.5)
lines(axisv,mke2,col='palevioletred',lty='dotted',lwd=2)

plot(axisv,nkf,col='khaki1',ylim=c(0,2.5),pch=15,cex=1.5,xlab='Taxonomic Rank',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
axis(1, at=1:6, labels=c('0','0.05','0.1','0.15','0.2','0.25'))
points(axisv,nkf,col='black',pch=22,cex=1.5)
lines(axisv,nkf,col='khaki1',lty='solid',lwd=2)
points(axisv,nkr,col='khaki1',pch=16,cex=1.5)
points(axisv,nkr,col='black',pch=21,cex=1.5)
lines(axisv,nkr,col='khaki1',lty='dashed',lwd=2)
points(axisv,nke,col='khaki1',pch=17,cex=1.5)
points(axisv,nke,col='black',pch=24,cex=1.5)
lines(axisv,nke,col='khaki1',lty='dotted',lwd=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

points(axisv,nkf2,col='pink',pch=15,cex=1.5)
points(axisv,nkf2,col='black',pch=22,cex=1.5)
lines(axisv,nkf2,col='pink',lty='solid',lwd=2)
points(axisv,nkr2,col='pink',pch=16,cex=1.5)
points(axisv,nkr2,col='black',pch=21,cex=1.5)
lines(axisv,nkr2,col='pink',lty='dashed',lwd=2)
points(axisv,nke2,col='pink',pch=17,cex=1.5)
points(axisv,nke2,col='black',pch=24,cex=1.5)
lines(axisv,nke2,col='pink',lty='dotted',lwd=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')

#Find maximum and minimum across ranks for each metric
x = 1:18
High = c(max(lkf),max(lke),max(lkr),max(mkf),max(mke),max(mkr),max(nkf),max(nke),max(nkr),max(lkf2),max(lke2),max(lkr2),max(mkf2),max(mke2),max(mkr2),max(nkf2),max(nke2),max(nkr2))
Low = c(min(lkf),min(lke),min(lkr),min(mkf),min(mke),min(mkr),min(nkf),min(nke),min(nkr),min(lkf2),min(lke2),min(lkr2),min(mkf2),min(mke2),min(mkr2),min(nkf2),min(nke2),min(nkr2))

## Plot variation across ranks
par(mar=c(10,5,1,3))
plot(1, type="n", xlab="", ylab="Relative Diversity Range", xlim=c(0.5,18.5),
     ylim=c(0, 3),xaxt='n',cex.lab=1.25)

## Add rectangles
rect(x[1:3] - 0.4, Low[1:3], x[1:3] + 0.4, High[1:3], col="goldenrod1")
## Add rectangles
rect(x[4:6] - 0.4, Low[4:6], x[4:6] + 0.4, High[4:6], col="yellow")
## Add rectangles
rect(x[7:9] - 0.4, Low[7:9], x[7:9] + 0.4, High[7:9], col="khaki1")
## Add rectangles
rect(x[10:12] - 0.4, Low[10:12], x[10:12] + 0.4, High[10:12], col="red")
## Add rectangles
rect(x[13:15] - 0.4, Low[13:15], x[13:15] + 0.4, High[13:15], col="palevioletred")
## Add rectangles
rect(x[16:18] - 0.4, Low[16:18], x[16:18] + 0.4, High[16:18], col="pink")
## Add rectangles
axis(1, at=1:18, labels=c('occ-BranchPD','occ-TipPD','occ-Richness','occ-ab-BranchPD','occ-ab-TipPD','occ-ab-Richness','shade-BranchPD','shade-TipPD','shade-Richness','occ-BranchPD','occ-TipPD','occ-Richness','occ-ab-BranchPD','occ-ab-TipPD','occ-ab-Richness','shade-BranchPD','shade-TipPD','shade-Richness'),las=2)
axis(2,cex=2)
lines(c(0,20),c(1,1),lty='dashed',col='grey')



