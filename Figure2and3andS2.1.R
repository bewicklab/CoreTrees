library('phyloseq')
library('dplyr')
library('ape')
library('vegan')
library('phytools')
library('data.table')
library('ggplot2')
library('holobiont')

manytextsize <- function(n, mins=0.5, maxs=4, B=6, D=100){
  # empirically selected size-value calculator.
  s <- B * exp(-n/D)
  # enforce a floor.
  s <- ifelse(s > mins, s, mins)
  # enforce a max
  s <- ifelse(s < maxs, s, maxs)
  return(s)
}


plot_tree2 = function(physeq,cv, method="sampledodge", nodelabf=NULL,
                      color=NULL, shape=NULL, size=NULL,
                      min.abundance=Inf, label.tips=NULL, text.size=NULL,
                      sizebase=5, base.spacing = 0.02,
                      ladderize=FALSE, plot.margin=0.2, title=NULL,
                      treetheme=NULL, justify="jagged"){
  ########################################
  # Support mis-capitalization of reserved variable names in color, shape, size
  # This helps, for instance, with backward-compatibility where "abundance"
  # was the reserved variable name for mapping OTU abundance entries
  fix_reserved_vars = function(aesvar){
    aesvar <- gsub("^abundance[s]{0,}$", "Abundance", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^OTU[s]{0,}$", "OTU", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^taxa_name[s]{0,}$", "OTU", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^sample[s]{0,}$", "Sample", aesvar, ignore.case=TRUE)
    return(aesvar)
  }
  if(!is.null(label.tips)){label.tips <- fix_reserved_vars(label.tips)}
  if(!is.null(color)){color <- fix_reserved_vars(color)}
  if(!is.null(shape)){shape <- fix_reserved_vars(shape)}
  if(!is.null(size) ){size  <- fix_reserved_vars(size)} 
  ########################################
  if( is.null(phy_tree(physeq, FALSE)) ){
    stop("There is no phylogenetic tree in the object you have provided.\n",
         "Try phy_tree(physeq) to see for yourself.")
  }
  if(!inherits(physeq, "phyloseq")){
    # If only a phylogenetic tree, then only tree available to overlay.
    method <- "treeonly"
  }
  # Create the tree data.table
  temptree<-phy_tree(physeq)
  temptree$edge.length <- sqrt(temptree$edge.length+1e-2)
  
  treeSegs <- tree_layout2(temptree, ladderize=ladderize)
  edgeMap = aes(x=xleft, xend=xright, y=y, yend=y)
  vertMap = aes(x=x, xend=x, y=vmin, yend=vmax)
  touchred<-union(treeSegs$edgeDT$V1[which(cv=='red')],treeSegs$edgeDT$V2[which(cv=='red')])
  cw<-rep('black',length(treeSegs$vertDT$V1))
  cw[intersect(which(treeSegs$vertDT$V1 %in% touchred),which(treeSegs$vertDT$V2 %in% touchred))]<-'red'
  # Initialize phylogenetic tree.
  # Naked, lines-only, unannotated tree as first layers. Edge (horiz) first, then vertical.
  p = ggplot(data=treeSegs$edgeDT) + geom_segment(edgeMap,col=cv) + 
    geom_segment(vertMap, data=treeSegs$vertDT,col=cw)
  # If no text.size given, calculate it from number of tips ("species", aka taxa)
  # This is very fast. No need to worry about whether text is printed or not.
  if(is.null(text.size)){
    text.size <- manytextsize(ntaxa(physeq))
  }
  # Add the species labels to the right.
  if(!is.null(label.tips) & method!="sampledodge"){
    # If method is sampledodge, then labels are added to the right of points, later.
    # Add labels layer to plotting object.
    labelDT = treeSegs$edgeDT[!is.na(OTU), ]
    if(!is.null(tax_table(object=physeq, errorIfNULL=FALSE))){
      # If there is a taxonomy available, merge it with the label data.table
      taxDT = data.table(tax_table(physeq), OTU=taxa_names(physeq), key="OTU")
      # Merge with taxonomy.
      labelDT = merge(x=labelDT, y=taxDT, by="OTU")
    }
    if(justify=="jagged"){
      # Tip label aesthetic mapping.
      # Aesthetics can be NULL, and that aesthetic gets ignored.
      labelMap <- aes_string(x="xright", y="y", label=label.tips, color=color)
    } else {
      # The left-justified version of tip-labels.
      labelMap <- aes_string(x="max(xright, na.rm=TRUE)", y="y", label=label.tips, color=color)
    }
    p <- p + geom_text(labelMap, data=labelDT, size=I(text.size), hjust=-0.1, na.rm=TRUE)
  }
  # Node label section.
  # 
  # If no nodelabf ("node label function") given, ask internal function to pick one.
  # Is NULL by default, meaning will dispatch to `howtolabnodes` to select function.
  # For no node labels, the "dummy" function `nodeplotblank` will return tree plot 
  # object, p, as-is, unmodified.
  if(is.null(nodelabf)){
    nodelabf = howtolabnodes(physeq)
  }
  #### set node `y` as the mean of the vertical segment
  # Use the provided/inferred node label function to add the node labels layer(s)
  # Non-root nodes first
  p = nodelabf(p, treeSegs$edgeDT[!is.na(label), ])
  # Add root label (if present)
  p = nodelabf(p, treeSegs$vertDT[!is.na(label), ])
  # Theme specification
  if(is.null(treetheme)){
    # If NULL, then use the default tree theme.
    treetheme <- theme(axis.ticks = element_blank(),
                       axis.title.x=element_blank(), axis.text.x=element_blank(),
                       axis.title.y=element_blank(), axis.text.y=element_blank(),
                       panel.background = element_blank(),
                       panel.grid.minor = element_blank(),      
                       panel.grid.major = element_blank(),rect = element_blank())   
  }
  if(inherits(treetheme, "theme")){
    # If a theme, add theme layer to plot. 
    # For all other cases, skip this, which will cause default theme to be used
    p <- p + treetheme
  }
  # Optionally add a title to the plot
  if(!is.null(title)){
    p <- p + ggtitle(title)
  }  
  if(method!="sampledodge"){
    # If anything but a sampledodge tree, return now without further decorations.
    return(p)
  }
  ########################################
  # Sample Dodge Section
  # Special words, c("Sample", "Abundance", "OTU")
  # See psmelt()
  ########################################
  # Initialize the species/taxa/OTU data.table
  dodgeDT = treeSegs$edgeDT[!is.na(OTU), ]
  # Merge with psmelt() result, to make all co-variables available
  dodgeDT = merge(x=dodgeDT, y=data.table(psmelt(physeq), key="OTU"), by="OTU")
  if(justify=="jagged"){
    # Remove 0 Abundance value entries now, not later, for jagged.
    dodgeDT <- dodgeDT[Abundance > 0, ]    
  }
  # Set key. Changes `dodgeDT` in place. OTU is first key, always.
  if( !is.null(color) | !is.null(shape) | !is.null(size) ){
    # If color, shape, or size is chosen, setkey by those as well
    setkeyv(dodgeDT, cols=c("OTU", color, shape, size))
  } else {
    # Else, set key by OTU and sample name. 
    setkey(dodgeDT, OTU, Sample)
  }
  # Add sample-dodge horizontal adjustment index. In-place data.table assignment
  dodgeDT[, h.adj.index := 1:length(xright), by=OTU]
  # `base.spacing` is a user-input parameter.
  # The sampledodge step size is based on this and the max `x` value
  if(justify=="jagged"){
    dodgeDT[, xdodge:=(xright + h.adj.index * base.spacing * max(xright, na.rm=TRUE))]
  } else {
    # Left-justified version, `xdodge` always starts at the max of all `xright` values.
    dodgeDT[, xdodge := max(xright, na.rm=TRUE) + h.adj.index * base.spacing * max(xright, na.rm=TRUE)]
    # zeroes removed only after all sample points have been mapped.
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  # The general tip-point map. Objects can be NULL, and that aesthetic gets ignored.
  dodgeMap <- aes_string(x="xdodge", y="y", color=color, fill=color,
                         shape=shape, size=size)
  p <- p + geom_point(dodgeMap, data=dodgeDT, na.rm=TRUE)
  # Adjust point size transform
  if( !is.null(size) ){
    p <- p + scale_size_continuous(trans=log_trans(sizebase))
  }  
  # Optionally-add abundance value label to each point.
  # User controls this by the `min.abundance` parameter.
  # A value of `Inf` implies no labels.
  if( any(dodgeDT$Abundance >= min.abundance[1]) ){
    pointlabdf = dodgeDT[Abundance>=min.abundance[1], ]
    p <- p + geom_text(mapping=aes(xdodge, y, label=Abundance),
                       data=pointlabdf, size=text.size, na.rm=TRUE)
  }
  # If indicated, add the species labels to the right of dodged points.
  if(!is.null(label.tips)){
    # `tiplabDT` has only one row per tip, the farthest horizontal
    # adjusted position (one for each taxa)
    tiplabDT = dodgeDT
    tiplabDT[, xfartiplab:=max(xdodge), by=OTU]
    tiplabDT <- tiplabDT[h.adj.index==1, .SD, by=OTU]
    if(!is.null(color)){
      if(color %in% sample_variables(physeq, errorIfNULL=FALSE)){
        color <- NULL
      }
    }
    labelMap <- NULL
    if(justify=="jagged"){
      labelMap <- aes_string(x="xfartiplab", y="y", label=label.tips, color=color)
    } else {
      labelMap <- aes_string(x="max(xfartiplab, na.rm=TRUE)", y="y", label=label.tips, color=color)
    }
    # Add labels layer to plotting object.
    p <- p + geom_text(labelMap, tiplabDT, size=I(text.size), hjust=-0.1, na.rm=TRUE)
  } 
  # Plot margins. 
  # Adjust the tree graphic plot margins.
  # Helps to manually ensure that graphic elements aren't clipped,
  # especially when there are long tip labels.
  min.x <- -0.01 # + dodgeDT[, min(c(xleft))]
  max.x <- dodgeDT[, max(xright, na.rm=TRUE)]
  if("xdodge" %in% names(dodgeDT)){
    max.x <- dodgeDT[, max(xright, xdodge, na.rm=TRUE)]
  }
  if(plot.margin > 0){
    max.x <- max.x * (1.0 + plot.margin)
  } 
  p <- p + scale_x_continuous(limits=c(min.x, max.x))  
  return(p)
}

tree_layout2 = function(phy, ladderize=FALSE){
  if(inherits(phy, "phyloseq")){
    phy = phy_tree(phy)
  }
  if(!inherits(phy, "phylo")){
    stop("tree missing or invalid. Please check `phy` argument and try again.")
  }
  if(is.null(phy$edge.length)){
    # If no edge lengths, set them all to value of 1 (dendrogram).
    phy$edge.length <- rep(1L, times=nrow(phy$edge))
  }
  # Perform ladderizing, if requested
  if(ladderize != FALSE){
    if(ladderize == "left"){
      phy <- ladderize(phy, FALSE)
    } else if(ladderize==TRUE | ladderize=="right"){
      phy <- ladderize(phy, TRUE)
    } else {
      stop("You did not specify a supported option for argument `ladderize`.")
    }
  }
  # 'z' is the tree in postorder order used in calls to .C
  # Descending order of left-hand side of edge (the ancestor to the node)
  z = reorder.phylo(phy, order="postorder")
  # Initialize some characteristics of the tree.
  Ntip = length(phy$tip.label)
  ROOT = Ntip + 1
  nodelabels = phy$node.label
  # Horizontal positions
  xx = node.depth.edgelength(phy)
  # vertical positions
  yy = node.height(phy = phy, clado.style = FALSE)
  # Initialize an edge data.table 
  # Don't set key, order matters
  edgeDT = data.table(phy$edge, edge.length=phy$edge.length, OTU=NA_character_)
  # Add tip.labels if present
  if(!is.null(phy$tip.label)){
    # Initialize OTU, set node (V2) as key, assign taxa_names as OTU label
    edgeDT[, OTU:=NA_character_]
    setkey(edgeDT, V2)
    edgeDT[V2 <= Ntip, OTU:=phy$tip.label]
  }
  # Add the mapping for each edge defined in `xx` and `yy` 
  edgeDT[, xleft:=xx[V1]]
  edgeDT[, xright:=xx[V2]]
  edgeDT[, y:=yy[V2]]
  # Next define vertical segments
  vertDT = edgeDT[, list(V2=V2,y=y, x=xleft[1], vmin=min(y), vmax=max(y),vhalf=min(y)+(max(y)-min(y))/2), by=V1, mult="last"]
  ul<-unique(vertDT$V1)
  for (k in 1:length(ul)){
    ttt<-which(vertDT$V1==ul[k])
    for (j in 1:length(ttt)){
      if (vertDT$y[ttt[j]]==vertDT$vmin[ttt[j]]){
        vertDT$vmax[ttt[j]]<-vertDT$vhalf[ttt[j]]
      }else{
        vertDT$vmin[ttt[j]]<-vertDT$vhalf[ttt[j]]
      }
    }
    
  }
  if(!is.null(phy$node.label)){
    # Add non-root node labels to edgeDT
    edgeDT[V2 > ROOT, x:=xright]
    edgeDT[V2 > ROOT, label:=phy$node.label[-1]]
    # Add root label (first node label) to vertDT
    setkey(vertDT, V1)
    vertDT[J(ROOT), y:=mean(c(vmin, vmax))]
    vertDT[J(ROOT), label:=phy$node.label[1]]
  }
  return(list(edgeDT=edgeDT, vertDT=vertDT))
}


ps<-readRDS(file='Joglekar_Staphylococcus_2023-main/phyloseq_object/ps_staph_ASV_filter_final.rds')

ps %>%
  refseq() %>%
  Biostrings::writeXStringSet("staph_asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-read.tree('Joglekar_Staphylococcus_2023-main/phyloseq_object/rooted_tree_out/tree.nwk')

# Apply sqrt transformation, adding a small constant 
#tree_file$edge.length <- sqrt(tree_file$edge.length+1e-2)

samp<-sample_data(ps)
head(data.frame(samp))

phy_tree(ps)<-multi2di(tree_file)
otu_table(ps)<-t(otu_table(ps))
ps<-rarefy_even_depth(ps,sample.size = 500,rngseed = 1)

people<-unique(sample_data(ps)$Subject_ID)

div1<-c()
div2<-c()
div3<-c()
lister<-c()
for (k in 1:length(people)){
  pps<-prune_samples(sample_data(ps)$Subject_ID==people[k],ps)
  div1<-c(div1,coreRichness(pps,core_fraction=0.49))
  div2<-c(div2,coreFaithsPD(pps,mode='tip',core_fraction=0.49))
  div3<-c(div3,coreFaithsPD(pps,core_fraction=0.49))
  coreTree(pps,remove_zeros = FALSE)
  if (length(sample_names(pps))>9){
    lister<-c(lister,k)
  }
}

coreRichASV<-c()
coreFaithASV<-c()
coreFaithASVt<-c()

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[1],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
a<-plot_tree2(pps,cc,nodelabf = nodeplotblank,color='Species')+scale_color_manual(values=c('red','orange','yellow','green','purple','cyan'))+ theme(legend.position="none") 
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[2],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
b<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','orange','yellow','green','blue','purple','lightsalmon','cyan'))+ theme(legend.position="none")
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[3],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
c<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','yellowgreen','green','blue','purple','pink','lightsteelblue1','magenta','slateblue1','darkseagreen1','cyan'))+ theme(legend.position="none")
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[4],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
d<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','green','purple','pink','magenta','slateblue1','cyan'))+ theme(legend.position="none")
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[11],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
e<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','yellow','yellowgreen','green','blue','purple','magenta','cyan'))+ theme(legend.position="none")
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[15],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
f<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','yellowgreen','green','blue','purple','pink','navy','salmon','cyan'))+ theme(legend.position="none")
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[17],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
g<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','yellowgreen','green','purple','salmon','magenta','khaki','cyan'))+ theme(legend.position="none")
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[20],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
h<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','green','blue','purple','pink','cyan'))+ theme(legend.position="none")
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[21],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
i<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','green','blue','purple','lightsteelblue1','magenta','slateblue1','cyan'))+ theme(legend.position="none")
coreRichASV<-c(coreRichASV,coreRichness(pps,core_fraction = 0.49))
coreFaithASVt<-c(coreFaithASVt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithASV<-c(coreFaithASV,coreFaithsPD(pps,core_fraction = 0.49))



coreRichOTU<-c()
coreFaithOTU<-c()
coreFaithOTUt<-c()

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[1],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
a1<-plot_tree2(pps,cc,nodelabf = nodeplotblank,color='Species')+scale_color_manual(values=c('orange','yellow','green','purple','cyan'))+ theme(legend.position="none") 
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[2],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
b1<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('orange','yellow','green','purple','salmon','cyan'))+ theme(legend.position="none")
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[3],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
c1<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('yellowgreen','green','purple','pink','lightsteelblue1','cyan'))+ theme(legend.position="none")
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[4],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
d1<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('green','purple','pink','cyan'))+ theme(legend.position="none")
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[11],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
e1<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('yellow','yellowgreen','green','purple','black','cyan'))+ theme(legend.position="none")
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[15],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
f1<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('yellowgreen','green','purple','pink','navy','salmon','cyan'))+ theme(legend.position="none")
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[17],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
g1<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('yellowgreen','green','purple','salmon','khaki','cyan'))+ theme(legend.position="none")
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[20],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
h1<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('green','purple','pink','cyan'))+ theme(legend.position="none")
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))

set.seed(1)

pps<-prune_samples(sample_data(ps)$Subject_ID==people[21],ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),10)],pps)
pps<-tip_glom(pps,h=0.0125)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
i1<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('green','purple','navy','cyan'))+ theme(legend.position="none")
coreRichOTU<-c(coreRichOTU,coreRichness(pps,core_fraction = 0.49))
coreFaithOTUt<-c(coreFaithOTUt,coreFaithsPD(pps,core_fraction = 0.49,mode='tip'))
coreFaithOTU<-c(coreFaithOTU,coreFaithsPD(pps,core_fraction = 0.49))


library("gridExtra")
library("grid")
library("cowplot")

plots <- list(a, b, c,d,e,f,g,h,i)


grid_layout <- plot_grid(a, a1, b, b1, c,c1, d,d1, e,e1, f,f1,g,g1,h,h1,i,i1, align = "v", nrow = 3, rel_heights = c(1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4))

ggsave("StaphFigureS2.1.png", grid_layout, height = 20, width = 30)

set.seed(1)

paperfig <- plot_grid(b, b1, h,h1, align = "v", nrow = 2, rel_heights = c(1/4, 1/4), rel_widths = c(1/4, 1/4))

ggsave("Fig2.png", paperfig, height = 8, width = 8)


vv<-c('P1','P2','P3','P4','P5','P6','P7','P8','P9')
df<-data.frame(c(coreRichASV/coreRichOTU,coreFaithASVt/coreFaithOTUt,coreFaithASV/coreFaithOTU),c(rep('Richness',9),rep('Tip',9),rep('Branch',9)),c(vv,vv,vv))
colnames(df)<-c('ratio','mode','person')
df$mode = factor(df$mode, levels = c("Richness","Tip", "Branch"), ordered = TRUE)

ggplot(data = df, aes(x = person, y = ratio, fill = mode))+ylab('ASV/OTU Diversity Ratio') +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)+scale_fill_manual(values=c('red','pink','purple4','mediumpurple1','blue','lightsteelblue1'))+ theme(axis.line = element_blank(),panel.background=element_blank(),panel.border = element_rect(colour = "black", fill=NA, linewidth=1),axis.text = element_text(size =12),axis.title = element_text(size =20))


#pps<-prune_taxa(rownames(tax_table(ps)),ps)
#list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
#cc<-rep('black',length(phy_tree(pps)$edge.length))
#ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
#marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
#gg<-phy_tree(pps)$edge
#marker2<-gg[,1]*100+gg[,2]
#redmarker2<-marker2[list$core_edges]
#cc[which(marker1 %in% redmarker2)]<-'red'
#plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','orange','yellow','yellowgreen','green','blue','purple','pink','navy','lightsteelblue1','salmon','magenta','black','khaki','slateblue1','darkseagreen1','cyan'))

#pps<-prune_taxa(rownames(tax_table(ps))[which(tax_table(ps)[,7] %in% c('Staphylococcus_&capitis&caprae','Staphylococcus_aureus','Staphylococcus_auricularis','Staphylococcus_cohnii','Staphylococcus_epidermidis','Staphylococcus_haemolyticus','Staphylococcus_hominis','Staphylococcus_lugdunensis','Staphylococcus_pettenkoferi','Staphylococcus_warneri'))],ps)
#list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
#cc<-rep('black',length(phy_tree(pps)$edge.length))
#ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
#marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
#gg<-phy_tree(pps)$edge
#marker2<-gg[,1]*100+gg[,2]
#redmarker2<-marker2[list$core_edges]
#cc[which(marker1 %in% redmarker2)]<-'red'
#plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','orange','yellow','yellowgreen','green','blue','purple','pink','salmon','cyan'))


set.seed(1)

pps<-prune_samples(sample_data(ps)$Habitat=='sebaceous',ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),20)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
sebaceous_tree<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','yellow','green','purple','navy','lightsteelblue1','magenta','slateblue1','darkseagreen1','cyan'))+ theme(legend.position="none")
spps<-pps

set.seed(1)
pps<-prune_samples(sample_data(ps)$Habitat=='dry',ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),20)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
dry_tree<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','orange','yellow','yellowgreen','green','blue','purple','pink','navy','lightsteelblue1','magenta','brown','slateblue1','darkseagreen1','cyan'))+ theme(legend.position="none")
dpps<-pps

set.seed(1)
pps<-prune_samples(sample_data(ps)$Habitat=='moist',ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),20)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
moist_tree<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','yellow','yellowgreen','green','blue','purple','pink','magenta','slateblue1','darkseagreen1','cyan'))+ theme(legend.position="none")
mpps<-pps


set.seed(1)
pps<-prune_samples(sample_data(ps)$Habitat=='feet',ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),20)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
feet_tree<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','yellow','yellowgreen','green','blue','purple','pink','magenta','slateblue1','darkseagreen1','cyan'))+ theme(legend.position="none")
fpps<-pps

set.seed(1)
pps<-prune_samples(sample_data(ps)$Habitat=='nare',ps)
pps<-prune_samples(sample_names(pps)[sample(length(sample_names(pps)),20)],pps)
list<-coreEdges(pps,remove_zeros=FALSE,core_fraction=0.49)
cc<-rep('black',length(phy_tree(pps)$edge.length))
ff<-tree_layout(phy_tree(pps),ladderize=FALSE)
marker1<-ff$edgeDT$V1*100+ff$edgeDT$V2
gg<-phy_tree(pps)$edge
marker2<-gg[,1]*100+gg[,2]
redmarker2<-marker2[list$core_edges]
cc[which(marker1 %in% redmarker2)]<-'red'
nare_tree<-plot_tree2(pps,cc,nodelabf=nodeplotblank,color='Species')+scale_color_manual(values=c('red','yellow','yellowgreen','green','blue','purple','pink','magenta','slateblue1','darkseagreen1','cyan'))+ theme(legend.position="none")
npps<-pps

ppb<-merge_phyloseq(spps,dpps,npps)
a<-coreVenn(ppb,sample_data(ppb)$Habitat,ordered_groups = c('nare','dry','sebaceous'),fill_color = c('green','tan3','khaki'),core_fraction=0.49)
b<-corePhyloVenn(ppb,sample_data(ppb)$Habitat,ordered_groups = c('nare','dry','sebaceous'),fill_color = c('green','tan3','khaki'),mode='tip',core_fraction=0.49)
c<-corePhyloVenn(ppb,sample_data(ppb)$Habitat,ordered_groups = c('nare','dry','sebaceous'),fill_color = c('green','tan3','khaki'),core_fraction=0.49)
coreVennTree(ppb,sample_data(ppb)$Habitat,core_fraction=0.49,scaled=TRUE,branch_color=c('black','khaki','tan3','green','darkgoldenrod1','yellowgreen','darkolivegreen','red','magenta','yellowgreen','magenta','magenta','magenta','lightsteelblue2','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','cyan',rep('magenta',4),'red'),edge_width=3)
ppb<-merge_phyloseq(spps,mpps,dpps,fpps,npps)
coreVennTree(ppb,sample_data(ppb)$Habitat,core_fraction=0.49,scaled=TRUE,branch_color=c('black','khaki','magenta','magenta','royalblue1','green','magenta','darkgoldenrod1','magenta','yellowgreen','magenta','magenta','magenta','lightsteelblue2','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','cyan',rep('magenta',4),'red'),edge_width=3)
coreVennTree(ppb,sample_data(ppb)$Habitat,core_fraction=0.49,scaled=TRUE,mode='tip',branch_color=c('black','khaki','magenta','magenta','royalblue1','green','magenta','darkgoldenrod1','magenta','yellowgreen','magenta','magenta','magenta','lightsteelblue2','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','magenta','cyan',rep('magenta',4),'red'),edge_width=3)

plot(1,1,pch=16,col='black',ylim=c(0,15))
points(1,2,pch=16,col='khaki')
points(1,3,pch=16,col='royalblue1')
points(1,4,pch=16,col='green')
points(1,5,pch=16,col='darkgoldenrod1')
points(1,6,pch=16,col='yellowgreen')
points(1,7,pch=16,col='lightsteelblue2')
points(1,8,pch=16,col='cyan')
points(1,9,pch=16,col='red')



hitme<-list(spps,dpps,mpps,fpps,npps)

JJ<-matrix(nrow=5,ncol=5)
UU<-matrix(nrow=5,ncol=5)
UUt<-matrix(nrow=5,ncol=5)

for (k in 1:5){
  for (j in 1:5){
    if (j == k){
      JJ[k,j]<-0
      UU[k,j]<-0
      UUt[k,j]<-0
    }else{
      tempphy<-merge_phyloseq(hitme[[k]],hitme[[j]])
      JJ[k,j]<-coreJaccard(tempphy,sample_data(tempphy)$Habitat,core_fraction=0.49)
      UU[k,j]<-coreUniFrac(tempphy,sample_data(tempphy)$Habitat,core_fraction=0.49)
      UUt[k,j]<-coreUniFrac(tempphy,sample_data(tempphy)$Habitat,core_fraction=0.49,mode='tip')
    }
  }
}


library(reshape2)
library(ggplot2)
df <- data.frame(id=c('sebaceous','dry','moist','feet','nare'),JJ)
df$id = factor(df$id, levels = c('sebaceous','dry','moist','feet','nare'), ordered = TRUE)
colnames(df)[2:ncol(df)] <- c('sebaceous','dry','moist','feet','nare')
gg <- melt(df, id="id")
ggplot(gg, aes(x=id,y=variable,fill=value))+
  geom_tile()+
  scale_fill_gradient2(low="red",mid="white",high='blue',midpoint=0.5)+
  coord_fixed()+theme(axis.title=element_blank(),axis.text=element_text(size=15),panel.background=element_blank(),panel.border=element_rect(colour = "black", fill=NA, linewidth=1))


library(reshape2)
library(ggplot2)
df <- data.frame(id=c('sebaceous','dry','moist','feet','nare'),UUt)
df$id = factor(df$id, levels = c('sebaceous','dry','moist','feet','nare'), ordered = TRUE)
colnames(df)[2:ncol(df)] <- c('sebaceous','dry','moist','feet','nare')
gg <- melt(df, id="id")
ggplot(gg, aes(x=id,y=variable,fill=value))+
  geom_tile()+
  scale_fill_gradient2(low="red",mid="white",high='blue',midpoint=0.5)+
  coord_fixed()+theme(axis.title=element_blank(),axis.text=element_text(size=15),panel.background=element_blank(),panel.border=element_rect(colour = "black", fill=NA, linewidth=1))



library(reshape2)
library(ggplot2)
df <- data.frame(id=c('sebaceous','dry','moist','feet','nare'),UU)
df$id = factor(df$id, levels = c('sebaceous','dry','moist','feet','nare'), ordered = TRUE)
colnames(df)[2:ncol(df)] <- c('sebaceous','dry','moist','feet','nare')
gg <- melt(df, id="id")
ggplot(gg, aes(x=id,y=variable,fill=value))+
  geom_tile()+
  scale_fill_gradient2(low="red",mid="white",high='blue',midpoint=0.5)+
  coord_fixed()+theme(axis.title=element_blank(),axis.text=element_text(size=15),panel.background=element_blank(),panel.border=element_rect(colour = "black", fill=NA, linewidth=1))



