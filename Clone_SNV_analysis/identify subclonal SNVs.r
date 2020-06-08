### Filter SNVs according to sample-specific frequency thresholds based on
### ploidy and senesence rate.

library(GenomicRanges)

### Commands for ploidy-adjusting the subclonal variant calls:
mutations_with_CNVs <- read.table( "mutations-cnv.txt", sep="\t", head=T )

### Add ploidy to the data frame
mutations_with_CNVs <- cbind( mutations_with_CNVs, "Ploidy"=4 * ( 2^(mutations_with_CNVs$log2) ))

### Get repeat regions for filtering variants:
Repeats <- read.table("../UCSC_RepeatMasker_hg19", sep="\t")

gr.repeats <- GRanges( Repeats[,1], IRanges( Repeats[,2], Repeats[,3]), strand=Repeats[,6], mcol=Repeats[,4])

### Each sample has its own VAF threshold (Delta) for subclonal mutations, based on senesence rates.
### Threshold calculated from Eqn 3 in Bozic et al PLoS Comp Bio, 2016 
### (see pg 12 of  supplementary materials for Cipponi et al for further details)
numeric.adj.subclon.thresholds <- as.num.vec( as.matrix(New.Deltas)[2,] )
numeric.adj.subclon.thresholds[11] <- numeric.adj.subclon.thresholds[11] + 0.1

names(numeric.adj.subclon.thresholds) <- names(New.Deltas)

### For counts of # of SNVs below threshold for each clone
Var.counts.Bozic <- list()

### Sets of keys (chr_start) for collected variants from the comparison group,
### in this case SKMEL parental (untreated) clone
compare.keys <- uniq.SKMEL.Par # uniq.all.778.var.pos

for ( sample.id %in% SKMEL.Vem.samples ) { 

  ### Fetch vars for sample
  sample.set <- subset(mutations_with_CNVs, Sample==sample.id)
  
  ### Covert threshold back to diploid from tetraploid (base ploidy of SKMEL28)
  Diploid.thres <- numeric.adj.subclon.thresholds[sample.id] * 4

  ### Adjust VAF threshold by rounded ploidy
  Pl.Adj.thres <- vapply(  sample.set$Ploidy, FUN=function(n) { min( 1, (Diploid.thres)/round(n) ) }, FUN.VALUE = -1 )
  
  ### Count total number of all variants (SNVs + indels):
  Var.counts.Bozic[[sample.id]][["Total"]] <- nrow(sample.set)
  
  ### Add chr_start keys for each variant, for overlaps & Filtering in next steps
  keys = paste(sample.set[,1], sample.set[,2], sep="_")
  sample.set <- cbind( "Key"=keys, sample.set, "Adj.thres"=Pl.Adj.thres )

  ### Remove vars found in parental clones
  trim.sample.set <- subset( sample.set, !Key %in% compare.keys ) 
  print( nrow(trim.sample.set) )

  ### Count number of vars unique to this clone
  Var.counts.Bozic[[sample.id]][["Unique"]] <- nrow(trim.sample.set)
    
  ## Filter by depth, only keeping SNVs
  trim.sample.set <- subset( trim.sample.set, Total.Depth >= 10 & Alt.Depth >= 2 & nchar(as.character(ALT))==1 & nchar(as.character(REF))==1 & Ploidy <= 6 & Chr!="Y") 

  Var.counts.Bozic[[sample.id]][["Filtered"]] <- nrow(trim.sample.set)

  ## Use GRanges to filter out SNVs in repeats
  trim.GRanges <- GRanges(paste("chr",trim.sample.set[,2],sep=""), IRanges( trim.sample.set[,3], trim.sample.set[,3]) )
  
  gr.findOv <- findOverlaps( trim.GRanges, gr.repeats )
  rpt.vars <- queryHits( gr.findOv )
  trim.sample.set <- trim.sample.set[ -rpt.vars, ]
  
  Var.counts.Bozic[[sample.id]][["No.Repeat"]] <- nrow(trim.sample.set)
  
  ### Filter by VAF to get just sub-clonal SNVs:
  trim.sample.set <- subset(trim.sample.set, `Alt.Freq` < Adj.thres )
  print( nrow(trim.sample.set) )
  Var.counts.Bozic[[sample.id]][["Below.Thres"]] <- nrow(trim.sample.set)
}

df.Bozic.counts.uniq2 <- data.frame(Var.counts.Bozic)

