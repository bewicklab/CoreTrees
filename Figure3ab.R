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
inornatus_list<-which(type=='inornatus')
marmoratus_list<-which(type=='marmoratus')
neo1_list<-which(sample_data(ZNC_physeq_all)$speciesxsite=='SNW_neo')
neo2_list<-which(sample_data(ZNC_physeq_all)$speciesxsite=='SBluG_neo')

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('inornatus','neomexicanus','marmoratus')

Zino<-prune_samples(sample_names(ZNC_physeq_all)[inornatus_list],ZNC_physeq_all)
Zmarm<-prune_samples(sample_names(ZNC_physeq_all)[marmoratus_list],ZNC_physeq_all)
Zneo1<-prune_samples(sample_names(ZNC_physeq_all)[neo1_list],ZNC_physeq_all)
Zneo2<-prune_samples(sample_names(ZNC_physeq_all)[neo2_list],ZNC_physeq_all)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate variation in Faith's PD across taxonomic rank

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inof<-c()
marmf<-c()
neo1f<-c()
neo2f<-c()
inor<-c()
marmr<-c()
neo1r<-c()
neo2r<-c()
inoe<-c()
marme<-c()
neo1e<-c()
neo2e<-c()

#branch-based ASV
inof<-c(inof,coreFaithsPD(Zino,0.5)/coreFaithsPD(Zino,0.5))
marmf<-c(marmf,coreFaithsPD(Zmarm,0.5)/coreFaithsPD(Zino,0.5))
neo1f<-c(neo1f,coreFaithsPD(Zneo1,0.5)/coreFaithsPD(Zino,0.5))
neo2f<-c(neo2f,coreFaithsPD(Zneo2,0.5)/coreFaithsPD(Zino,0.5))

#tip-based ASV
inoe<-c(inoe,coreFaithsPD(Zino,0.5,mode='tip')/coreFaithsPD(Zino,0.5,mode='tip'))
marme<-c(marme,coreFaithsPD(Zmarm,0.5,mode='tip')/coreFaithsPD(Zino,0.5,mode='tip'))
neo1e<-c(neo1e,coreFaithsPD(Zneo1,0.5,mode='tip')/coreFaithsPD(Zino,0.5,mode='tip'))
neo2e<-c(neo2e,coreFaithsPD(Zneo2,0.5,mode='tip')/coreFaithsPD(Zino,0.5,mode='tip'))

#non-phylogenetic ASV
inor<-c(inor,length(which(rowSums(sign(otu_table(Zino)))>0.5*length(sample_names(Zino))))/length(which(rowSums(sign(otu_table(Zino)))>0.5*length(sample_names(Zino)))))
marmr<-c(marmr,length(which(rowSums(sign(otu_table(Zmarm)))>0.5*length(sample_names(Zmarm))))/length(which(rowSums(sign(otu_table(Zino)))>0.5*length(sample_names(Zino)))))
neo1r<-c(neo1r,length(which(rowSums(sign(otu_table(Zneo1)))>0.5*length(sample_names(Zneo1))))/length(which(rowSums(sign(otu_table(Zino)))>0.5*length(sample_names(Zino)))))
neo2r<-c(neo2r,length(which(rowSums(sign(otu_table(Zneo2)))>0.5*length(sample_names(Zneo2))))/length(which(rowSums(sign(otu_table(Zino)))>0.5*length(sample_names(Zino)))))

#Pool to species
Zino_s<-tax_glom(Zino,taxrank='Rank7')
Zmarm_s<-tax_glom(Zmarm,taxrank='Rank7')
Zneo1_s<-tax_glom(Zneo1,taxrank='Rank7')
Zneo2_s<-tax_glom(Zneo2,taxrank='Rank7')

#branch-based
inof<-c(inof,coreFaithsPD(Zino_s,0.5)/coreFaithsPD(Zino_s,0.5))
marmf<-c(marmf,coreFaithsPD(Zmarm_s,0.5)/coreFaithsPD(Zino_s,0.5))
neo1f<-c(neo1f,coreFaithsPD(Zneo1_s,0.5)/coreFaithsPD(Zino_s,0.5))
neo2f<-c(neo2f,coreFaithsPD(Zneo2_s,0.5)/coreFaithsPD(Zino_s,0.5))

#tip-based
inor<-c(inor,length(which(rowSums(sign(otu_table(Zino_s)))>0.5*length(sample_names(Zino_s))))/length(which(rowSums(sign(otu_table(Zino_s)))>0.5*length(sample_names(Zino_s)))))
marmr<-c(marmr,length(which(rowSums(sign(otu_table(Zmarm_s)))>0.5*length(sample_names(Zmarm_s))))/length(which(rowSums(sign(otu_table(Zino_s)))>0.5*length(sample_names(Zino_s)))))
neo1r<-c(neo1r,length(which(rowSums(sign(otu_table(Zneo1_s)))>0.5*length(sample_names(Zneo1_s))))/length(which(rowSums(sign(otu_table(Zino_s)))>0.5*length(sample_names(Zino_s)))))
neo2r<-c(neo2r,length(which(rowSums(sign(otu_table(Zneo2_s)))>0.5*length(sample_names(Zneo2_s))))/length(which(rowSums(sign(otu_table(Zino_s)))>0.5*length(sample_names(Zino_s)))))

#non-phylogenetic
inoe<-c(inoe,coreFaithsPD(Zino_s,0.5,mode='tip')/coreFaithsPD(Zino_s,0.5,mode='tip'))
marme<-c(marme,coreFaithsPD(Zmarm_s,0.5,mode='tip')/coreFaithsPD(Zino_s,0.5,mode='tip'))
neo1e<-c(neo1e,coreFaithsPD(Zneo1_s,0.5,mode='tip')/coreFaithsPD(Zino_s,0.5,mode='tip'))
neo2e<-c(neo2e,coreFaithsPD(Zneo2_s,0.5,mode='tip')/coreFaithsPD(Zino_s,0.5,mode='tip'))

#Pool to genus
Zino_g<-tax_glom(Zino,taxrank='Rank6')
Zmarm_g<-tax_glom(Zmarm,taxrank='Rank6')
Zneo1_g<-tax_glom(Zneo1,taxrank='Rank6')
Zneo2_g<-tax_glom(Zneo2,taxrank='Rank6')

#branch-based
inof<-c(inof,coreFaithsPD(Zino_g,0.5)/coreFaithsPD(Zino_g,0.5))
marmf<-c(marmf,coreFaithsPD(Zmarm_g,0.5)/coreFaithsPD(Zino_g,0.5))
neo1f<-c(neo1f,coreFaithsPD(Zneo1_g,0.5)/coreFaithsPD(Zino_g,0.5))
neo2f<-c(neo2f,coreFaithsPD(Zneo2_g,0.5)/coreFaithsPD(Zino_g,0.5))

#tip-based
inoe<-c(inoe,coreFaithsPD(Zino_g,0.5,mode='tip')/coreFaithsPD(Zino_g,0.5,mode='tip'))
marme<-c(marme,coreFaithsPD(Zmarm_g,0.5,mode='tip')/coreFaithsPD(Zino_g,0.5,mode='tip'))
neo1e<-c(neo1e,coreFaithsPD(Zneo1_g,0.5,mode='tip')/coreFaithsPD(Zino_g,0.5,mode='tip'))
neo2e<-c(neo2e,coreFaithsPD(Zneo2_g,0.5,mode='tip')/coreFaithsPD(Zino_g,0.5,mode='tip'))

#non-phylogenetic
inor<-c(inor,length(which(rowSums(sign(otu_table(Zino_g)))>0.5*length(sample_names(Zino_g))))/length(which(rowSums(sign(otu_table(Zino_g)))>0.5*length(sample_names(Zino_g)))))
marmr<-c(marmr,length(which(rowSums(sign(otu_table(Zmarm_g)))>0.5*length(sample_names(Zmarm_g))))/length(which(rowSums(sign(otu_table(Zino_g)))>0.5*length(sample_names(Zino_g)))))
neo1r<-c(neo1r,length(which(rowSums(sign(otu_table(Zneo1_g)))>0.5*length(sample_names(Zneo1_g))))/length(which(rowSums(sign(otu_table(Zino_g)))>0.5*length(sample_names(Zino_g)))))
neo2r<-c(neo2r,length(which(rowSums(sign(otu_table(Zneo2_g)))>0.5*length(sample_names(Zneo2_g))))/length(which(rowSums(sign(otu_table(Zino_g)))>0.5*length(sample_names(Zino_g)))))

#Pool to family
Zino_f<-tax_glom(Zino,taxrank='Rank5')
Zmarm_f<-tax_glom(Zmarm,taxrank='Rank5')
Zneo1_f<-tax_glom(Zneo1,taxrank='Rank5')
Zneo2_f<-tax_glom(Zneo2,taxrank='Rank5')

#branch-based
inof<-c(inof,coreFaithsPD(Zino_f,0.5)/coreFaithsPD(Zino_f,0.5))
marmf<-c(marmf,coreFaithsPD(Zmarm_f,0.5)/coreFaithsPD(Zino_f,0.5))
neo1f<-c(neo1f,coreFaithsPD(Zneo1_f,0.5)/coreFaithsPD(Zino_f,0.5))
neo2f<-c(neo2f,coreFaithsPD(Zneo2_f,0.5)/coreFaithsPD(Zino_f,0.5))

#tip-based
inoe<-c(inoe,coreFaithsPD(Zino_f,0.5,mode='tip')/coreFaithsPD(Zino_f,0.5,mode='tip'))
marme<-c(marme,coreFaithsPD(Zmarm_f,0.5,mode='tip')/coreFaithsPD(Zino_f,0.5,mode='tip'))
neo1e<-c(neo1e,coreFaithsPD(Zneo1_f,0.5,mode='tip')/coreFaithsPD(Zino_f,0.5,mode='tip'))
neo2e<-c(neo2e,coreFaithsPD(Zneo2_f,0.5,mode='tip')/coreFaithsPD(Zino_f,0.5,mode='tip'))

#non-phylogenetic
inor<-c(inor,length(which(rowSums(sign(otu_table(Zino_f)))>0.5*length(sample_names(Zino_f))))/length(which(rowSums(sign(otu_table(Zino_f)))>0.5*length(sample_names(Zino_f)))))
marmr<-c(marmr,length(which(rowSums(sign(otu_table(Zmarm_f)))>0.5*length(sample_names(Zmarm_f))))/length(which(rowSums(sign(otu_table(Zino_f)))>0.5*length(sample_names(Zino_f)))))
neo1r<-c(neo1r,length(which(rowSums(sign(otu_table(Zneo1_f)))>0.5*length(sample_names(Zneo1_f))))/length(which(rowSums(sign(otu_table(Zino_f)))>0.5*length(sample_names(Zino_f)))))
neo2r<-c(neo2r,length(which(rowSums(sign(otu_table(Zneo2_f)))>0.5*length(sample_names(Zneo2_f))))/length(which(rowSums(sign(otu_table(Zino_f)))>0.5*length(sample_names(Zino_f)))))

#Pool to order
Zino_o<-tax_glom(Zino,taxrank='Rank4')
Zmarm_o<-tax_glom(Zmarm,taxrank='Rank4')
Zneo1_o<-tax_glom(Zneo1,taxrank='Rank4')
Zneo2_o<-tax_glom(Zneo2,taxrank='Rank4')

#branch-based
inof<-c(inof,coreFaithsPD(Zino_o,0.5)/coreFaithsPD(Zino_o,0.5))
marmf<-c(marmf,coreFaithsPD(Zmarm_o,0.5)/coreFaithsPD(Zino_o,0.5))
neo1f<-c(neo1f,coreFaithsPD(Zneo1_o,0.5)/coreFaithsPD(Zino_o,0.5))
neo2f<-c(neo2f,coreFaithsPD(Zneo2_o,0.5)/coreFaithsPD(Zino_o,0.5))

#tip-based
inoe<-c(inoe,coreFaithsPD(Zino_o,0.5,mode='tip')/coreFaithsPD(Zino_o,0.5,mode='tip'))
marme<-c(marme,coreFaithsPD(Zmarm_o,0.5,mode='tip')/coreFaithsPD(Zino_o,0.5,mode='tip'))
neo1e<-c(neo1e,coreFaithsPD(Zneo1_o,0.5,mode='tip')/coreFaithsPD(Zino_o,0.5,mode='tip'))
neo2e<-c(neo2e,coreFaithsPD(Zneo2_o,0.5,mode='tip')/coreFaithsPD(Zino_o,0.5,mode='tip'))

#non-phylogenetic
inor<-c(inor,length(which(rowSums(sign(otu_table(Zino_o)))>0.5*length(sample_names(Zino_o))))/length(which(rowSums(sign(otu_table(Zino_o)))>0.5*length(sample_names(Zino_o)))))
marmr<-c(marmr,length(which(rowSums(sign(otu_table(Zmarm_o)))>0.5*length(sample_names(Zmarm_o))))/length(which(rowSums(sign(otu_table(Zino_o)))>0.5*length(sample_names(Zino_o)))))
neo1r<-c(neo1r,length(which(rowSums(sign(otu_table(Zneo1_o)))>0.5*length(sample_names(Zneo1_o))))/length(which(rowSums(sign(otu_table(Zino_o)))>0.5*length(sample_names(Zino_o)))))
neo2r<-c(neo2r,length(which(rowSums(sign(otu_table(Zneo2_o)))>0.5*length(sample_names(Zneo2_o))))/length(which(rowSums(sign(otu_table(Zino_o)))>0.5*length(sample_names(Zino_o)))))

#Plot variation in Faith's PD across taxonomic ranks
axisv<-c(1,2,3,4,5)

par(mar=c(5,5,5,5))
plot(axisv,marmf,col='red',ylim=c(0.5,2),pch=15,cex=1.5,xlab='Taxonomic Rank',ylab='Relative Diversity',xaxt='n',cex.lab=1.5)
axis(1, at=1:5, labels=c('ASV','species','genus','family','order'))
points(axisv,marmf,col='black',pch=22,cex=1.5)
lines(axisv,marmf,col='red',lty='solid',lwd=2)
points(axisv,marmr,col='red',pch=16,cex=1.5)
points(axisv,marmr,col='black',pch=21,cex=1.5)
lines(axisv,marmr,col='red',lty='dashed',lwd=2)
points(axisv,marme,col='red',pch=17,cex=1.5)
points(axisv,marme,col='black',pch=24,cex=1.5)
lines(axisv,marme,col='red',lty='dotted',lwd=2)

points(axisv,neo1f,col='blue',ylim=c(0.6,2),pch=15,cex=1.5)
points(axisv,neo1f,col='black',pch=22,cex=1.5)
lines(axisv,neo1f,col='blue',lty='solid',lwd=2)
points(axisv,neo1r,col='blue',pch=16,cex=1.5)
points(axisv,neo1r,col='black',pch=21,cex=1.5)
lines(axisv,neo1r,col='blue',lty='dashed',lwd=2)
points(axisv,neo1e,col='blue',pch=17,cex=1.5)
points(axisv,neo1e,col='black',pch=24,cex=1.5)
lines(axisv,neo1e,col='blue',lty='dotted',lwd=2)

points(axisv,neo2f,col='green',ylim=c(0.6,2),pch=15,cex=1.5)
points(axisv,neo2f,col='black',pch=22,cex=1.5)
lines(axisv,neo2f,col='green',lty='solid',lwd=2)
points(axisv,neo2r,col='green',pch=16,cex=1.5)
points(axisv,neo2r,col='black',pch=21,cex=1.5)
lines(axisv,neo2r,col='green',lty='dashed',lwd=2)
points(axisv,neo2e,col='green',pch=17,cex=1.5)
points(axisv,neo2e,col='black',pch=24,cex=1.5)
lines(axisv,neo2e,col='green',lty='dotted',lwd=2)

#Find maximum and minimum across ranks for each metric
x = 1:9
High = c(max(marmf),max(marme),max(marmr),max(neo1f),max(neo1e),max(neo1r),max(neo2f),max(neo2e),max(neo2r))
Low = c(min(marmf),min(marme),min(marmr),min(neo1f),min(neo1e),min(neo1r),min(neo2f),min(neo2e),min(neo2r))

## Plot variation across ranks
par(mar=c(7,5,5,1))
plot(1, type="n", xlab="", ylab="Relative Diversity Range", xlim=c(0.5,9.5),
     ylim=c(0.5, 2),xaxt='n',cex.lab=1.5)

## Add rectangles
rect(x[1:3] - 0.4, Low[1:3], x[1:3] + 0.4, High[1:3], col="red")
## Add rectangles
rect(x[4:6] - 0.4, Low[4:6], x[4:6] + 0.4, High[4:6], col="blue")
## Add rectangles
rect(x[7:9] - 0.4, Low[7:9], x[7:9] + 0.4, High[7:9], col="green")
axis(1, at=1:9, labels=c('M-BranchPD','M-TipPD','M-Richness','N1-BranchPD','N1-TipPD','N1-Richness','N2-BranchPD','N2-TipPD','N2-Richness'),las=2)
axis(2,cex=2)
