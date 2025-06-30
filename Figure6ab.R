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

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('inornatus','neomexicanus','marmoratus')

Zinovmarm<-merge_phyloseq(Zino,Zmarm)
Zinovneo<-merge_phyloseq(Zino,Zneo)
Zmarmvneo<-merge_phyloseq(Zmarm,Zneo)

Zinovmarm_s<-tax_glom(Zinovmarm,taxrank='Rank7')
Zinovmarm_g<-tax_glom(Zinovmarm,taxrank='Rank6')
Zinovmarm_f<-tax_glom(Zinovmarm,taxrank='Rank5')
Zinovmarm_o<-tax_glom(Zinovmarm,taxrank='Rank4')

Zinovneo_s<-tax_glom(Zinovneo,taxrank='Rank7')
Zinovneo_g<-tax_glom(Zinovneo,taxrank='Rank6')
Zinovneo_f<-tax_glom(Zinovneo,taxrank='Rank5')
Zinovneo_o<-tax_glom(Zinovneo,taxrank='Rank4')

Zmarmvneo_s<-tax_glom(Zmarmvneo,taxrank='Rank7')
Zmarmvneo_g<-tax_glom(Zmarmvneo,taxrank='Rank6')
Zmarmvneo_f<-tax_glom(Zmarmvneo,taxrank='Rank5')
Zmarmvneo_o<-tax_glom(Zmarmvneo,taxrank='Rank4')



###########################Branch-based tree UniFrac###########################

lkf<-c()
lke<-c()
lkr<-c()

#branch-based ASV
lkf<-c(lkf,coreUniFrac(Zinovneo,sample_data(Zinovneo)$species,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zinovneo_s,sample_data(Zinovneo_s)$species,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zinovneo_g,sample_data(Zinovneo_g)$species,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zinovneo_f,sample_data(Zinovneo_f)$species,core_fraction=0.5))
lkf<-c(lkf,coreUniFrac(Zinovneo_o,sample_data(Zinovneo_o)$species,core_fraction=0.5))

#tip-based ASV
lke<-c(lke,coreUniFrac(Zinovneo,sample_data(Zinovneo)$species,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zinovneo_s,sample_data(Zinovneo_s)$species,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zinovneo_g,sample_data(Zinovneo_g)$species,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zinovneo_f,sample_data(Zinovneo_f)$species,mode='tip',core_fraction=0.5))
lke<-c(lke,coreUniFrac(Zinovneo_o,sample_data(Zinovneo_o)$species,mode='tip',core_fraction=0.5))

#non-phylogenetic ASV
lkr<-c(lkr,coreJaccard(Zinovneo,sample_data(Zinovneo)$species,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zinovneo_s,sample_data(Zinovneo_s)$species,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zinovneo_g,sample_data(Zinovneo_g)$species,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zinovneo_f,sample_data(Zinovneo_f)$species,core_fraction=0.5))
lkr<-c(lkr,coreJaccard(Zinovneo_o,sample_data(Zinovneo_o)$species,core_fraction=0.5))


lkf2<-c()
lke2<-c()
lkr2<-c()

#branch-based ASV
lkf2<-c(lkf2,coreUniFrac(Zmarmvneo,sample_data(Zmarmvneo)$species,core_fraction=0.5))
lkf2<-c(lkf2,coreUniFrac(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,core_fraction=0.5))
lkf2<-c(lkf2,coreUniFrac(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,core_fraction=0.5))
lkf2<-c(lkf2,coreUniFrac(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,core_fraction=0.5))
lkf2<-c(lkf2,coreUniFrac(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,core_fraction=0.5))

#tip-based ASV
lke2<-c(lke2,coreUniFrac(Zmarmvneo,sample_data(Zmarmvneo)$species,mode='tip',core_fraction=0.5))
lke2<-c(lke2,coreUniFrac(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,mode='tip',core_fraction=0.5))
lke2<-c(lke2,coreUniFrac(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,mode='tip',core_fraction=0.5))
lke2<-c(lke2,coreUniFrac(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,mode='tip',core_fraction=0.5))
lke2<-c(lke2,coreUniFrac(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,mode='tip',core_fraction=0.5))

#non-phylogenetic ASV
lkr2<-c(lkr2,coreJaccard(Zmarmvneo,sample_data(Zmarmvneo)$species,core_fraction=0.5))
lkr2<-c(lkr2,coreJaccard(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,core_fraction=0.5))
lkr2<-c(lkr2,coreJaccard(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,core_fraction=0.5))
lkr2<-c(lkr2,coreJaccard(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,core_fraction=0.5))
lkr2<-c(lkr2,coreJaccard(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,core_fraction=0.5))


mkf<-c()
mke<-c()
mkr<-c()

#branch-based ASV
mkf<-c(mkf,coreUniFrac(Zinovneo,sample_data(Zinovneo)$species,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zinovneo_s,sample_data(Zinovneo_s)$species,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zinovneo_g,sample_data(Zinovneo_g)$species,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zinovneo_f,sample_data(Zinovneo_f)$species,core_fraction=0.5,ab_threshold1=0.1))
mkf<-c(mkf,coreUniFrac(Zinovneo_o,sample_data(Zinovneo_o)$species,core_fraction=0.5,ab_threshold1=0.1))

#tip-based ASV
mke<-c(mke,coreUniFrac(Zinovneo,sample_data(Zinovneo)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zinovneo_s,sample_data(Zinovneo_s)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zinovneo_g,sample_data(Zinovneo_g)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zinovneo_f,sample_data(Zinovneo_f)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke<-c(mke,coreUniFrac(Zinovneo_o,sample_data(Zinovneo_o)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))

#non-phylogenetic ASV
mkr<-c(mkr,coreJaccard(Zinovneo,sample_data(Zinovneo)$species,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zinovneo_s,sample_data(Zinovneo_s)$species,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zinovneo_g,sample_data(Zinovneo_g)$species,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zinovneo_f,sample_data(Zinovneo_f)$species,core_fraction=0.5,ab_threshold1=0.1))
mkr<-c(mkr,coreJaccard(Zinovneo_o,sample_data(Zinovneo_o)$species,core_fraction=0.5,ab_threshold1=0.1))

mkf2<-c()
mke2<-c()
mkr2<-c()

#branch-based ASV
mkf2<-c(mkf2,coreUniFrac(Zmarmvneo,sample_data(Zmarmvneo)$species,core_fraction=0.5,ab_threshold1=0.1))
mkf2<-c(mkf2,coreUniFrac(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,core_fraction=0.5,ab_threshold1=0.1))
mkf2<-c(mkf2,coreUniFrac(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,core_fraction=0.5,ab_threshold1=0.1))
mkf2<-c(mkf2,coreUniFrac(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,core_fraction=0.5,ab_threshold1=0.1))
mkf2<-c(mkf2,coreUniFrac(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,core_fraction=0.5,ab_threshold1=0.1))

#tip-based ASV
mke2<-c(mke2,coreUniFrac(Zmarmvneo,sample_data(Zmarmvneo)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke2<-c(mke2,coreUniFrac(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke2<-c(mke2,coreUniFrac(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke2<-c(mke2,coreUniFrac(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))
mke2<-c(mke2,coreUniFrac(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,mode='tip',core_fraction=0.5,ab_threshold1=0.1))

#non-phylogenetic ASV
mkr2<-c(mkr2,coreJaccard(Zmarmvneo,sample_data(Zmarmvneo)$species,core_fraction=0.5,ab_threshold1=0.1))
mkr2<-c(mkr2,coreJaccard(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,core_fraction=0.5,ab_threshold1=0.1))
mkr2<-c(mkr2,coreJaccard(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,core_fraction=0.5,ab_threshold1=0.1))
mkr2<-c(mkr2,coreJaccard(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,core_fraction=0.5,ab_threshold1=0.1))
mkr2<-c(mkr2,coreJaccard(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,core_fraction=0.5,ab_threshold1=0.1))


nkf<-c()
nke<-c()
nkr<-c()

#branch-based ASV
nkf<-c(nkf,coreUniFrac(Zinovneo,sample_data(Zinovneo)$species,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zinovneo_s,sample_data(Zinovneo_s)$species,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zinovneo_g,sample_data(Zinovneo_g)$species,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zinovneo_f,sample_data(Zinovneo_f)$species,selection='shade',initial_branches=10))
nkf<-c(nkf,coreUniFrac(Zinovneo_o,sample_data(Zinovneo_o)$species,selection='shade',initial_branches=10))

#tip-based ASV
nke<-c(nke,coreUniFrac(Zinovneo,sample_data(Zinovneo)$species,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zinovneo_s,sample_data(Zinovneo_s)$species,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zinovneo_g,sample_data(Zinovneo_g)$species,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zinovneo_f,sample_data(Zinovneo_f)$species,mode='tip',selection='shade',initial_branches=10))
nke<-c(nke,coreUniFrac(Zinovneo_o,sample_data(Zinovneo_o)$species,mode='tip',selection='shade',initial_branches=10))

#non-phylogenetic ASV
nkr<-c(nkr,coreJaccard(Zinovneo,sample_data(Zinovneo)$species,selection='shade'))
nkr<-c(nkr,coreJaccard(Zinovneo_s,sample_data(Zinovneo_s)$species,selection='shade'))
nkr<-c(nkr,coreJaccard(Zinovneo_g,sample_data(Zinovneo_g)$species,selection='shade'))
nkr<-c(nkr,coreJaccard(Zinovneo_f,sample_data(Zinovneo_f)$species,selection='shade'))
nkr<-c(nkr,coreJaccard(Zinovneo_o,sample_data(Zinovneo_o)$species,selection='shade'))

nkf2<-c()
nke2<-c()
nkr2<-c()

#branch-based ASV
nkf2<-c(nkf2,coreUniFrac(Zmarmvneo,sample_data(Zmarmvneo)$species,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreUniFrac(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreUniFrac(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreUniFrac(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,selection='shade',initial_branches=10))
nkf2<-c(nkf2,coreUniFrac(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,selection='shade',initial_branches=10))

#tip-based ASV
nke2<-c(nke2,coreUniFrac(Zmarmvneo,sample_data(Zmarmvneo)$species,mode='tip',selection='shade',initial_branches=10))
nke2<-c(nke2,coreUniFrac(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,mode='tip',selection='shade',initial_branches=10))
nke2<-c(nke2,coreUniFrac(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,mode='tip',selection='shade',initial_branches=10))
nke2<-c(nke2,coreUniFrac(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,mode='tip',selection='shade',initial_branches=10))
nke2<-c(nke2,coreUniFrac(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,mode='tip',selection='shade',initial_branches=10))

#non-phylogenetic ASV
nkr2<-c(nkr2,coreJaccard(Zmarmvneo,sample_data(Zmarmvneo)$species,selection='shade'))
nkr2<-c(nkr2,coreJaccard(Zmarmvneo_s,sample_data(Zmarmvneo_s)$species,selection='shade'))
nkr2<-c(nkr2,coreJaccard(Zmarmvneo_g,sample_data(Zmarmvneo_g)$species,selection='shade'))
nkr2<-c(nkr2,coreJaccard(Zmarmvneo_f,sample_data(Zmarmvneo_f)$species,selection='shade'))
nkr2<-c(nkr2,coreJaccard(Zmarmvneo_o,sample_data(Zmarmvneo_o)$species,selection='shade'))

#Plot variation in Faith's PD across taxonoZneo ranks
axisv<-c(1,2,3,4,5)

par(mar=c(5,5,5,5))
plot(axisv,lkf,col='goldenrod1',ylim=c(0.3,1),pch=15,cex=1.5,xlab='Taxonomic Rank',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('ASV','species','genus','family','order'))
points(axisv,lkf,col='black',pch=22,cex=1.5)
lines(axisv,lkf,col='goldenrod1',lty='solid',lwd=2)
points(axisv,lkr,col='goldenrod1',pch=16,cex=1.5)
points(axisv,lkr,col='black',pch=21,cex=1.5)
lines(axisv,lkr,col='goldenrod1',lty='dashed',lwd=2)
points(axisv,lke,col='goldenrod1',pch=17,cex=1.5)
points(axisv,lke,col='black',pch=24,cex=1.5)
lines(axisv,lke,col='goldenrod1',lty='dotted',lwd=2)

points(axisv,lkf2,col='red',pch=15,cex=1.5)
points(axisv,lkf2,col='black',pch=22,cex=1.5)
lines(axisv,lkf2,col='red',lty='solid',lwd=2)
points(axisv,lkr2,col='red',pch=16,cex=1.5)
points(axisv,lkr2,col='black',pch=21,cex=1.5)
lines(axisv,lkr2,col='red',lty='dashed',lwd=2)
points(axisv,lke2,col='red',pch=17,cex=1.5)
points(axisv,lke2,col='black',pch=24,cex=1.5)
lines(axisv,lke2,col='red',lty='dotted',lwd=2)

par(mar=c(5,5,5,5))
plot(axisv,mkf,col='yellow',ylim=c(0.3,1),pch=15,cex=1.5,xlab='Taxonomic Rank',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('ASV','species','genus','family','order'))
points(axisv,mkf,col='black',pch=22,cex=1.5)
lines(axisv,mkf,col='yellow',lty='solid',lwd=2)
points(axisv,mkr,col='yellow',pch=16,cex=1.5)
points(axisv,mkr,col='black',pch=21,cex=1.5)
lines(axisv,mkr,col='yellow',lty='dashed',lwd=2)
points(axisv,mke,col='yellow',pch=17,cex=1.5)
points(axisv,mke,col='black',pch=24,cex=1.5)
lines(axisv,mke,col='yellow',lty='dotted',lwd=2)

points(axisv,mkf2,col='palevioletred',pch=15,cex=1.5)
points(axisv,mkf2,col='black',pch=22,cex=1.5)
lines(axisv,mkf2,col='palevioletred',lty='solid',lwd=2)
points(axisv,mkr2,col='palevioletred',pch=16,cex=1.5)
points(axisv,mkr2,col='black',pch=21,cex=1.5)
lines(axisv,mkr2,col='palevioletred',lty='dashed',lwd=2)
points(axisv,mke2,col='palevioletred',pch=17,cex=1.5)
points(axisv,mke2,col='black',pch=24,cex=1.5)

plot(axisv,nkf,col='khaki1',ylim=c(0.2,1),pch=15,cex=1.5,xlab='Taxonomic Rank',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:6, labels=c('0','0.05','0.1','0.15','0.2','0.25'))
points(axisv,nkf,col='black',pch=22,cex=1.5)
lines(axisv,nkf,col='khaki1',lty='solid',lwd=2)
points(axisv,nkr,col='khaki1',pch=16,cex=1.5)
points(axisv,nkr,col='black',pch=21,cex=1.5)
lines(axisv,nkr,col='khaki1',lty='dashed',lwd=2)
points(axisv,nke,col='khaki1',pch=17,cex=1.5)
points(axisv,nke,col='black',pch=24,cex=1.5)
lines(axisv,nke,col='khaki1',lty='dotted',lwd=2)

points(axisv,nkf2,col='pink',pch=15,cex=1.5)
points(axisv,nkf2,col='black',pch=22,cex=1.5)
lines(axisv,nkf2,col='pink',lty='solid',lwd=2)
points(axisv,nkr2,col='pink',pch=16,cex=1.5)
points(axisv,nkr2,col='black',pch=21,cex=1.5)
lines(axisv,nkr2,col='pink',lty='dashed',lwd=2)
points(axisv,nke2,col='pink',pch=17,cex=1.5)
points(axisv,nke2,col='black',pch=24,cex=1.5)
lines(axisv,nke2,col='pink',lty='dotted',lwd=2)

#Find maximum and minimum across ranks for each metric
x = 1:18
High = c(max(lkf),max(lke),max(lkr),max(mkf),max(mke),max(mkr),max(nkf),max(nke),max(nkr),max(lkf2),max(lke2),max(lkr2),max(mkf2),max(mke2),max(mkr2),max(nkf2),max(nke2),max(nkr2))
Low = c(min(lkf),min(lke),min(lkr),min(mkf),min(mke),min(mkr),min(nkf),min(nke),min(nkr),min(lkf2),min(lke2),min(lkr2),min(mkf2),min(mke2),min(mkr2),min(nkf2),min(nke2),min(nkr2))

## Plot variation across ranks
par(mar=c(10,5,1,3))
plot(1, type="n", xlab="", ylab="Turnover Range", xlim=c(0.5,18.5),
     ylim=c(0.25, 1),xaxt='n',cex.lab=1.25)

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
axis(1, at=1:18, labels=c('occ-BranchUniFrac','occ-TipUniFrac','occ-Jaccard','occ-ab-BranchUniFrac','occ-ab-TipUniFrac','occ-ab-Jaccard','shade-BranchUniFrac','shade-TipUniFrac','shade-Jaccard','occ-BranchUniFrac','occ-TipUniFrac','occ-Jaccard','occ-ab-BranchUniFrac','occ-ab-TipUniFrac','occ-ab-Jaccard','shade-BranchUniFrac','shade-TipUniFrac','shade-Jaccard'),las=2)

