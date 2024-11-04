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


#Find the treatment group for each of your samples
type<-sample_data(ZNC_physeq_all)$species

#Find the treatment group for each of your samples
type<-sample_data(ZNC_physeq_all)$species
inornatus_list<-which(type=='inornatus')
marmoratus_list<-which(type=='marmoratus')
neo1_list<-which(sample_data(ZNC_physeq_all)$speciesxsite=='SNW_neo')
neo2_list<-which(sample_data(ZNC_physeq_all)$speciesxsite=='SBluG_neo')

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('inornatus','neomexicanus','marmoratus')

Zinovmarm<-prune_samples(sample_names(ZNC_physeq_all)[c(inornatus_list,marmoratus_list)],ZNC_physeq_all)
groupinovmarm<-rep(1,length(sample_names(Zinovmarm)))
groupinovmarm[which(sample_data(Zinovmarm)$species=='inornatus')]<-2

Zinovneo1<-prune_samples(sample_names(ZNC_physeq_all)[c(inornatus_list,neo1_list)],ZNC_physeq_all)
groupinovneo1<-rep(1,length(sample_names(Zinovneo1)))
groupinovneo1[which(sample_data(Zinovneo1)$species=='inornatus')]<-2

Zinovneo2<-prune_samples(sample_names(ZNC_physeq_all)[c(inornatus_list,neo2_list)],ZNC_physeq_all)
groupinovneo2<-rep(1,length(sample_names(Zinovneo2)))
groupinovneo2[which(sample_data(Zinovneo2)$species=='inornatus')]<-2

Zmarmvneo1<-prune_samples(sample_names(ZNC_physeq_all)[c(marmoratus_list,neo1_list)],ZNC_physeq_all)
groupmarmvneo1<-rep(1,length(sample_names(Zmarmvneo1)))
groupmarmvneo1[which(sample_data(Zmarmvneo1)$species=='marmoratus')]<-2

Zmarmvneo2<-prune_samples(sample_names(ZNC_physeq_all)[c(marmoratus_list,neo2_list)],ZNC_physeq_all)
groupmarmvneo2<-rep(1,length(sample_names(Zmarmvneo2)))
groupmarmvneo2[which(sample_data(Zmarmvneo2)$species=='marmoratus')]<-2

Zneo1vneo2<-prune_samples(sample_names(ZNC_physeq_all)[c(neo1_list,neo2_list)],ZNC_physeq_all)
groupneo1vneo2<-rep(1,length(sample_names(Zneo1vneo2)))
groupneo1vneo2[which(sample_data(Zneo1vneo2)$speciesxsite=='SBluG_neo')]<-2

Zinovmarm_s<-tax_glom(Zinovmarm,taxrank='Rank7')
Zinovmarm_g<-tax_glom(Zinovmarm,taxrank='Rank6')
Zinovmarm_f<-tax_glom(Zinovmarm,taxrank='Rank5')
Zinovmarm_o<-tax_glom(Zinovmarm,taxrank='Rank4')

Zinovneo1_s<-tax_glom(Zinovneo1,taxrank='Rank7')
Zinovneo1_g<-tax_glom(Zinovneo1,taxrank='Rank6')
Zinovneo1_f<-tax_glom(Zinovneo1,taxrank='Rank5')
Zinovneo1_o<-tax_glom(Zinovneo1,taxrank='Rank4')

Zinovneo2_s<-tax_glom(Zinovneo2,taxrank='Rank7')
Zinovneo2_g<-tax_glom(Zinovneo2,taxrank='Rank6')
Zinovneo2_f<-tax_glom(Zinovneo2,taxrank='Rank5')
Zinovneo2_o<-tax_glom(Zinovneo2,taxrank='Rank4')

Zmarmvneo1_s<-tax_glom(Zmarmvneo1,taxrank='Rank7')
Zmarmvneo1_g<-tax_glom(Zmarmvneo1,taxrank='Rank6')
Zmarmvneo1_f<-tax_glom(Zmarmvneo1,taxrank='Rank5')
Zmarmvneo1_o<-tax_glom(Zmarmvneo1,taxrank='Rank4')

Zmarmvneo2_s<-tax_glom(Zmarmvneo2,taxrank='Rank7')
Zmarmvneo2_g<-tax_glom(Zmarmvneo2,taxrank='Rank6')
Zmarmvneo2_f<-tax_glom(Zmarmvneo2,taxrank='Rank5')
Zmarmvneo2_o<-tax_glom(Zmarmvneo2,taxrank='Rank4')

Zneo1vneo2_s<-tax_glom(Zneo1vneo2,taxrank='Rank7')
Zneo1vneo2_g<-tax_glom(Zneo1vneo2,taxrank='Rank6')
Zneo1vneo2_f<-tax_glom(Zneo1vneo2,taxrank='Rank5')
Zneo1vneo2_o<-tax_glom(Zneo1vneo2,taxrank='Rank4')

###########################Branch-based tree UniFrac###########################

Zinovmarm_list<-c()
Zinovmarm_list<-c(Zinovmarm_list,coreUniFrac(Zinovmarm,groupinovmarm,0.5))
Zinovmarm_list<-c(Zinovmarm_list,coreUniFrac(Zinovmarm_s,groupinovmarm,0.5))
Zinovmarm_list<-c(Zinovmarm_list,coreUniFrac(Zinovmarm_g,groupinovmarm,0.5))
Zinovmarm_list<-c(Zinovmarm_list,coreUniFrac(Zinovmarm_f,groupinovmarm,0.5))
Zinovmarm_list<-c(Zinovmarm_list,coreUniFrac(Zinovmarm_o,groupinovmarm,0.5))

Zinovneo1_list<-c()
Zinovneo1_list<-c(Zinovneo1_list,coreUniFrac(Zinovneo1,groupinovneo1,0.5))
Zinovneo1_list<-c(Zinovneo1_list,coreUniFrac(Zinovneo1_s,groupinovneo1,0.5))
Zinovneo1_list<-c(Zinovneo1_list,coreUniFrac(Zinovneo1_g,groupinovneo1,0.5))
Zinovneo1_list<-c(Zinovneo1_list,coreUniFrac(Zinovneo1_f,groupinovneo1,0.5))
Zinovneo1_list<-c(Zinovneo1_list,coreUniFrac(Zinovneo1_o,groupinovneo1,0.5))

Zinovneo2_list<-c()
Zinovneo2_list<-c(Zinovneo2_list,coreUniFrac(Zinovneo2,groupinovneo1,0.5))
Zinovneo2_list<-c(Zinovneo2_list,coreUniFrac(Zinovneo2_s,groupinovneo1,0.5))
Zinovneo2_list<-c(Zinovneo2_list,coreUniFrac(Zinovneo2_g,groupinovneo1,0.5))
Zinovneo2_list<-c(Zinovneo2_list,coreUniFrac(Zinovneo2_f,groupinovneo1,0.5))
Zinovneo2_list<-c(Zinovneo2_list,coreUniFrac(Zinovneo2_o,groupinovneo1,0.5))


Zmarmvneo1_list<-c()
Zmarmvneo1_list<-c(Zmarmvneo1_list,coreUniFrac(Zmarmvneo1,groupmarmvneo1,0.5))
Zmarmvneo1_list<-c(Zmarmvneo1_list,coreUniFrac(Zmarmvneo1_s,groupmarmvneo1,0.5))
Zmarmvneo1_list<-c(Zmarmvneo1_list,coreUniFrac(Zmarmvneo1_g,groupmarmvneo1,0.5))
Zmarmvneo1_list<-c(Zmarmvneo1_list,coreUniFrac(Zmarmvneo1_f,groupmarmvneo1,0.5))
Zmarmvneo1_list<-c(Zmarmvneo1_list,coreUniFrac(Zmarmvneo1_o,groupmarmvneo1,0.5))

Zmarmvneo2_list<-c()
Zmarmvneo2_list<-c(Zmarmvneo2_list,coreUniFrac(Zmarmvneo2,groupmarmvneo2,0.5))
Zmarmvneo2_list<-c(Zmarmvneo2_list,coreUniFrac(Zmarmvneo2_s,groupmarmvneo2,0.5))
Zmarmvneo2_list<-c(Zmarmvneo2_list,coreUniFrac(Zmarmvneo2_g,groupmarmvneo2,0.5))
Zmarmvneo2_list<-c(Zmarmvneo2_list,coreUniFrac(Zmarmvneo2_f,groupmarmvneo2,0.5))
Zmarmvneo2_list<-c(Zmarmvneo2_list,coreUniFrac(Zmarmvneo2_o,groupmarmvneo2,0.5))

Zneo1vneo2_list<-c()
Zneo1vneo2_list<-c(Zneo1vneo2_list,coreUniFrac(Zneo1vneo2,groupneo1vneo2,0.5))
Zneo1vneo2_list<-c(Zneo1vneo2_list,coreUniFrac(Zneo1vneo2_s,groupneo1vneo2,0.5))
Zneo1vneo2_list<-c(Zneo1vneo2_list,coreUniFrac(Zneo1vneo2_g,groupneo1vneo2,0.5))
Zneo1vneo2_list<-c(Zneo1vneo2_list,coreUniFrac(Zneo1vneo2_f,groupneo1vneo2,0.5))
Zneo1vneo2_list<-c(Zneo1vneo2_list,coreUniFrac(Zneo1vneo2_o,groupneo1vneo2,0.5))

###########################Tip-based tree UniFrac###########################

Zinovmarm_listt<-c()
Zinovmarm_listt<-c(Zinovmarm_listt,coreUniFrac(Zinovmarm,groupinovmarm,0.5, mode='tip'))
Zinovmarm_listt<-c(Zinovmarm_listt,coreUniFrac(Zinovmarm_s,groupinovmarm,0.5, mode='tip'))
Zinovmarm_listt<-c(Zinovmarm_listt,coreUniFrac(Zinovmarm_g,groupinovmarm,0.5, mode='tip'))
Zinovmarm_listt<-c(Zinovmarm_listt,coreUniFrac(Zinovmarm_f,groupinovmarm,0.5, mode='tip'))
Zinovmarm_listt<-c(Zinovmarm_listt,coreUniFrac(Zinovmarm_o,groupinovmarm,0.5, mode='tip'))

Zinovneo1_listt<-c()
Zinovneo1_listt<-c(Zinovneo1_listt,coreUniFrac(Zinovneo1,groupinovneo1,0.5, mode='tip'))
Zinovneo1_listt<-c(Zinovneo1_listt,coreUniFrac(Zinovneo1_s,groupinovneo1,0.5, mode='tip'))
Zinovneo1_listt<-c(Zinovneo1_listt,coreUniFrac(Zinovneo1_g,groupinovneo1,0.5, mode='tip'))
Zinovneo1_listt<-c(Zinovneo1_listt,coreUniFrac(Zinovneo1_f,groupinovneo1,0.5, mode='tip'))
Zinovneo1_listt<-c(Zinovneo1_listt,coreUniFrac(Zinovneo1_o,groupinovneo1,0.5, mode='tip'))

Zinovneo2_listt<-c()
Zinovneo2_listt<-c(Zinovneo2_listt,coreUniFrac(Zinovneo2,groupinovneo1,0.5, mode='tip'))
Zinovneo2_listt<-c(Zinovneo2_listt,coreUniFrac(Zinovneo2_s,groupinovneo1,0.5, mode='tip'))
Zinovneo2_listt<-c(Zinovneo2_listt,coreUniFrac(Zinovneo2_g,groupinovneo1,0.5, mode='tip'))
Zinovneo2_listt<-c(Zinovneo2_listt,coreUniFrac(Zinovneo2_f,groupinovneo1,0.5, mode='tip'))
Zinovneo2_listt<-c(Zinovneo2_listt,coreUniFrac(Zinovneo2_o,groupinovneo1,0.5, mode='tip'))


Zmarmvneo1_listt<-c()
Zmarmvneo1_listt<-c(Zmarmvneo1_listt,coreUniFrac(Zmarmvneo1,groupmarmvneo1,0.5, mode='tip'))
Zmarmvneo1_listt<-c(Zmarmvneo1_listt,coreUniFrac(Zmarmvneo1_s,groupmarmvneo1,0.5, mode='tip'))
Zmarmvneo1_listt<-c(Zmarmvneo1_listt,coreUniFrac(Zmarmvneo1_g,groupmarmvneo1,0.5, mode='tip'))
Zmarmvneo1_listt<-c(Zmarmvneo1_listt,coreUniFrac(Zmarmvneo1_f,groupmarmvneo1,0.5, mode='tip'))
Zmarmvneo1_listt<-c(Zmarmvneo1_listt,coreUniFrac(Zmarmvneo1_o,groupmarmvneo1,0.5, mode='tip'))

Zmarmvneo2_listt<-c()
Zmarmvneo2_listt<-c(Zmarmvneo2_listt,coreUniFrac(Zmarmvneo2,groupmarmvneo2,0.5, mode='tip'))
Zmarmvneo2_listt<-c(Zmarmvneo2_listt,coreUniFrac(Zmarmvneo2_s,groupmarmvneo2,0.5, mode='tip'))
Zmarmvneo2_listt<-c(Zmarmvneo2_listt,coreUniFrac(Zmarmvneo2_g,groupmarmvneo2,0.5, mode='tip'))
Zmarmvneo2_listt<-c(Zmarmvneo2_listt,coreUniFrac(Zmarmvneo2_f,groupmarmvneo2,0.5, mode='tip'))
Zmarmvneo2_listt<-c(Zmarmvneo2_listt,coreUniFrac(Zmarmvneo2_o,groupmarmvneo2,0.5, mode='tip'))

Zneo1vneo2_listt<-c()
Zneo1vneo2_listt<-c(Zneo1vneo2_listt,coreUniFrac(Zneo1vneo2,groupneo1vneo2,0.5, mode='tip'))
Zneo1vneo2_listt<-c(Zneo1vneo2_listt,coreUniFrac(Zneo1vneo2_s,groupneo1vneo2,0.5, mode='tip'))
Zneo1vneo2_listt<-c(Zneo1vneo2_listt,coreUniFrac(Zneo1vneo2_g,groupneo1vneo2,0.5, mode='tip'))
Zneo1vneo2_listt<-c(Zneo1vneo2_listt,coreUniFrac(Zneo1vneo2_f,groupneo1vneo2,0.5, mode='tip'))
Zneo1vneo2_listt<-c(Zneo1vneo2_listt,coreUniFrac(Zneo1vneo2_o,groupneo1vneo2,0.5, mode='tip'))


###########################Jaccard###########################


Zinovmarm_listJ<-c()
Zinovmarm_listJ<-c(Zinovmarm_listJ,coreJaccard(Zinovmarm,groupinovmarm,0.5))
Zinovmarm_listJ<-c(Zinovmarm_listJ,coreJaccard(Zinovmarm_s,groupinovmarm,0.5))
Zinovmarm_listJ<-c(Zinovmarm_listJ,coreJaccard(Zinovmarm_g,groupinovmarm,0.5))
Zinovmarm_listJ<-c(Zinovmarm_listJ,coreJaccard(Zinovmarm_f,groupinovmarm,0.5))
Zinovmarm_listJ<-c(Zinovmarm_listJ,coreJaccard(Zinovmarm_o,groupinovmarm,0.5))

Zinovneo1_listJ<-c()
Zinovneo1_listJ<-c(Zinovneo1_listJ,coreJaccard(Zinovneo1,groupinovneo1,0.5))
Zinovneo1_listJ<-c(Zinovneo1_listJ,coreJaccard(Zinovneo1_s,groupinovneo1,0.5))
Zinovneo1_listJ<-c(Zinovneo1_listJ,coreJaccard(Zinovneo1_g,groupinovneo1,0.5))
Zinovneo1_listJ<-c(Zinovneo1_listJ,coreJaccard(Zinovneo1_f,groupinovneo1,0.5))
Zinovneo1_listJ<-c(Zinovneo1_listJ,coreJaccard(Zinovneo1_o,groupinovneo1,0.5))

Zinovneo2_listJ<-c()
Zinovneo2_listJ<-c(Zinovneo2_listJ,coreJaccard(Zinovneo2,groupinovneo1,0.5))
Zinovneo2_listJ<-c(Zinovneo2_listJ,coreJaccard(Zinovneo2_s,groupinovneo1,0.5))
Zinovneo2_listJ<-c(Zinovneo2_listJ,coreJaccard(Zinovneo2_g,groupinovneo1,0.5))
Zinovneo2_listJ<-c(Zinovneo2_listJ,coreJaccard(Zinovneo2_f,groupinovneo1,0.5))
Zinovneo2_listJ<-c(Zinovneo2_listJ,coreJaccard(Zinovneo2_o,groupinovneo1,0.5))


Zmarmvneo1_listJ<-c()
Zmarmvneo1_listJ<-c(Zmarmvneo1_listJ,coreJaccard(Zmarmvneo1,groupmarmvneo1,0.5))
Zmarmvneo1_listJ<-c(Zmarmvneo1_listJ,coreJaccard(Zmarmvneo1_s,groupmarmvneo1,0.5))
Zmarmvneo1_listJ<-c(Zmarmvneo1_listJ,coreJaccard(Zmarmvneo1_g,groupmarmvneo1,0.5))
Zmarmvneo1_listJ<-c(Zmarmvneo1_listJ,coreJaccard(Zmarmvneo1_f,groupmarmvneo1,0.5))
Zmarmvneo1_listJ<-c(Zmarmvneo1_listJ,coreJaccard(Zmarmvneo1_o,groupmarmvneo1,0.5))

Zmarmvneo2_listJ<-c()
Zmarmvneo2_listJ<-c(Zmarmvneo2_listJ,coreJaccard(Zmarmvneo2,groupmarmvneo2,0.5))
Zmarmvneo2_listJ<-c(Zmarmvneo2_listJ,coreJaccard(Zmarmvneo2_s,groupmarmvneo2,0.5))
Zmarmvneo2_listJ<-c(Zmarmvneo2_listJ,coreJaccard(Zmarmvneo2_g,groupmarmvneo2,0.5))
Zmarmvneo2_listJ<-c(Zmarmvneo2_listJ,coreJaccard(Zmarmvneo2_f,groupmarmvneo2,0.5))
Zmarmvneo2_listJ<-c(Zmarmvneo2_listJ,coreJaccard(Zmarmvneo2_o,groupmarmvneo2,0.5))

Zneo1vneo2_listJ<-c()
Zneo1vneo2_listJ<-c(Zneo1vneo2_listJ,coreJaccard(Zneo1vneo2,groupneo1vneo2,0.5))
Zneo1vneo2_listJ<-c(Zneo1vneo2_listJ,coreJaccard(Zneo1vneo2_s,groupneo1vneo2,0.5))
Zneo1vneo2_listJ<-c(Zneo1vneo2_listJ,coreJaccard(Zneo1vneo2_g,groupneo1vneo2,0.5))
Zneo1vneo2_listJ<-c(Zneo1vneo2_listJ,coreJaccard(Zneo1vneo2_f,groupneo1vneo2,0.5))
Zneo1vneo2_listJ<-c(Zneo1vneo2_listJ,coreJaccard(Zneo1vneo2_o,groupneo1vneo2,0.5))

axisv<-c(1,2,3,4,5)

par(mar=c(5,5,5,5))
plot(axisv,Zinovmarm_list,col='red',ylim=c(0.4,1),pch=15,cex=1.5,xlab='Taxonomic Rank',ylab='Turnover',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('ASV','species','genus','family','order'))
points(axisv,Zinovmarm_list,col='black',pch=22,cex=1.5)
lines(axisv,Zinovmarm_list,col='red',lty='solid',lwd=2)
points(axisv,Zinovmarm_listJ,col='red',pch=16,cex=1.5)
points(axisv,Zinovmarm_listJ,col='black',pch=21,cex=1.5)
lines(axisv,Zinovmarm_listJ,col='red',lty='dashed',lwd=2)
points(axisv,Zinovmarm_listt,col='red',pch=17,cex=1.5)
points(axisv,Zinovmarm_listt,col='black',pch=24,cex=1.5)
lines(axisv,Zinovmarm_listt,col='red',lty='dotted',lwd=2)

points(axisv,Zinovneo1_list,col='blue',ylim=c(0.6,2),pch=15,cex=1.5)
points(axisv,Zinovneo1_list,col='black',pch=22,cex=1.5)
lines(axisv,Zinovneo1_list,col='blue',lty='solid',lwd=2)
points(axisv,Zinovneo1_listJ,col='blue',pch=16,cex=1.5)
points(axisv,Zinovneo1_listJ,col='black',pch=21,cex=1.5)
lines(axisv,Zinovneo1_listJ,col='blue',lty='dashed',lwd=2)
points(axisv,Zinovneo1_listt,col='blue',pch=17,cex=1.5)
points(axisv,Zinovneo1_listt,col='black',pch=24,cex=1.5)
lines(axisv,Zinovneo1_listt,col='blue',lty='dotted',lwd=2)

points(axisv,Zinovneo2_list,col='green',ylim=c(0.6,2),pch=15,cex=1.5)
points(axisv,Zinovneo2_list,col='black',pch=22,cex=1.5)
lines(axisv,Zinovneo2_list,col='green',lty='solid',lwd=2)
points(axisv,Zinovneo2_listJ,col='green',pch=16,cex=1.5)
points(axisv,Zinovneo2_listJ,col='black',pch=21,cex=1.5)
lines(axisv,Zinovneo2_listJ,col='green',lty='dashed',lwd=2)
points(axisv,Zinovneo2_listt,col='green',pch=17,cex=1.5)
points(axisv,Zinovneo2_listt,col='black',pch=24,cex=1.5)
lines(axisv,Zinovneo2_listt,col='green',lty='dotted',lwd=2)

x = 1:9
High = c(max(Zinovmarm_list),max(Zinovmarm_listt),max(Zinovmarm_listJ),max(Zinovneo1_list),max(Zinovneo1_listt),max(Zinovneo1_listJ),max(Zinovneo2_list),max(Zinovneo2_listt),max(Zinovneo2_listJ))
Low = c(min(Zinovmarm_list),min(Zinovmarm_listt),min(Zinovmarm_listJ),min(Zinovneo1_list),min(Zinovneo1_listt),min(Zinovneo1_listJ),min(Zinovneo2_list),min(Zinovneo2_listt),min(Zinovneo2_listJ))

## Blank plot
par(mar=c(9,5,5,1))
plot(1, type="n", xlab="", ylab="Turnover Range", xlim=c(0.5,9.5),
     ylim=c(0.4, 1),xaxt='n',cex.lab=1.5)

## Add rectangles
rect(x[1:3] - 0.4, Low[1:3], x[1:3] + 0.4, High[1:3], col="red")
## Add rectangles
rect(x[4:6] - 0.4, Low[4:6], x[4:6] + 0.4, High[4:6], col="blue")
## Add rectangles
rect(x[7:9] - 0.4, Low[7:9], x[7:9] + 0.4, High[4:6], col="green")
axis(1, at=1:9, labels=c('M-BranchUniFrac','M-TipUniFrac','M-Jaccard','N1-BranchUniFrac','N1-TipUniFrac','N1-Jaccard','N2-BranchUniFrac','N2-TipUniFrac','N2-Jaccard'),las=2)


