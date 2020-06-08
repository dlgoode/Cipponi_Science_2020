### Code for processing Z-score data from LINCS L1000, 
### extracting Z-scores for DNA repair gene sets of interest 
### and plotting combined Z-scores as boxplots.

### Example here is shown for creating a boxplot of the data from the mTOR inhibitor experiments
### The same code was used to create boxplots for the cell-cycle inhibitors, with a different starting set of GCT files

library(cmap)  ## Library written for LINCS L1000 data
library(ggplot2)

DNArepair_gene_sets <- list(read.table("DNA_Repair_gene_sets.txt",sep="\t",header=T, fill=T))

## Read in L1000 data (downloaded from clue.io) for drugs of interest
mTOR.GCT.files <- list.files("mTOR inhibitor GCT files/")

torin2.gct <- parse.gctx(fname=paste(mTOR.GCT.dir,mTOR.GCT.files[3],sep = ""))
torin1.gct <- parse.gctx(fname=paste(mTOR.GCT.dir,mTOR.GCT.files[2],sep = ""))
AZD_8055.gct <- parse.gctx(fname=paste(mTOR.GCT.dir,mTOR.GCT.files[1],sep = ""))

### Gather gcts into one list:
mTOR.inhibitor.gcts <- list()

mTOR.inhibitor.gcts[["AZD-8055"]] <- AZD_8055.gct
mTOR.inhibitor.gcts[["torin-1"]] <- torin1.gct
mTOR.inhibitor.gcts[["torin-2"]] <- torin2.gct

### Extract Z-scores for genes in each gene set from each gct 
### and combine into one matrix.
melted.data.all.mTOR.inhibitors <- list()

for( drug in names(mTOR.inhibitor.gcts) ) {  
  
  print( dim( mTOR.inhibitor.gcts[[drug]]) )
  
  gct <- mTOR.inhibitor.gcts[[drug]]

  combined.data <- matrix(nrow=0, ncol=11)
  colnames(combined.data) <- c( colnames( melt.set.gct )[ cols.to.keep[-11] ], "Pathway" )
  
  for ( set in DNArepair_gene_sets) {
  	set.gct <- subset.gct(current.gct, rid=current.gct@rdesc$pr_gene_symbol %in% set, cid=grepl(drug, current.gct@cdesc$pert_iname) & grepl(time, current.gct@cdesc$pert_itime) & grepl("10 ÂµM", current.gct@cdesc$pert_idose) )
   
    melt.set.gct <- melt.gct(set.gct)
    
    nm <- cbind( melt.set.gct, "Pathway"=rep(set,nrow(melt.set.gct)))
    
    combined.data <- rbind( combined.data, as.matrix(nm)[,c( cols.to.keep )] )
  }
  
  combined.data.df <- data.frame( combined.data )
  combined.data.df$Pathway <- factor( combined.data.df$Pathway, levels=key.gene.sets, ordered=T)
  
  melted.data.all.mTOR.inhibitors[[drug]] <- combined.data.df
}

## Visualize all Z-scores for all genes in a given gene set across all drugs:

df.mTOR.combined <- do.call("rbind", melted.data.all.mTOR.inhibitors)

mTOR.cmb.plot <- ggplot(aes(y = as.numeric(as.vector( value ) ), x = Pathway, fill=cell_id), data = df.mTOR.combined) + geom_boxplot(outlier.size = 0.3, width=0.1, ) + theme_minimal() + geom_hline(yintercept = 0, linetype="dotted",color="grey25")  + labs(y="L1000 Level 5 Z-score", x="Gene Set")

ggsave("Combined.mTOR.DNA.repair.sets.boxes.pdf", width = 10, height = 6, units = "in", mTOR.cmb.plot)
