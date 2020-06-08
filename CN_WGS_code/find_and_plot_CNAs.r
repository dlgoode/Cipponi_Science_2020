### Plot CNVs again filtering through blacklist first
### Followed by segmentation by DNAcopy, merging of common regions
### into large scale gains/losses and plotting

library(GenomicRanges)

### Get the blacklist from Scheinin et al 2014
Blacklist <- read.table( "FW__data_for_breakpoint_analysis/Scheinin_et_al_Blacklist.txt", sep="\t" )

### Turn into GRanges object
GR.Blacklist <- GRanges(  seqnames = Rle(Blacklist[,1]), IRanges(Blacklist[,2], Blacklist[,3]) )

### Turn log ratios for the 778 Tunicamycin expt into GRanges obj
GR.comb.logratios <- GRanges( seqnames=Rle(combined.cnn.logratios$chr), 
                              IRanges(combined.cnn.logratios[,2],combined.cnn.logratios[,3]), 
                              mcol=combined.cnn.logratios[4:13] )

## Remove windows overlapping with the blacklist
BL.filt.GR.comb.logratios <- combined.cnn.logratios[ !GR.comb.logratios %over% GR.Blacklist , ]

### Segment windows into CN call
library(DNAcopy)

## Function for running segmentation algorithm
get.CNA.segments <- function( logratios ) { 
  CNA.object <- CNA(cbind( logratios[,-c(1:3)]), logratios$chr, logratios$start, data.type="logratio", sampleid = colnames(logratios)[-c(1:3)] )
  smoothed.CNA<- smooth.CNA(CNA.object)
  segment.combined.cnn <- segment( smoothed.CNA, verbose=1 )
  
  return( segment.combined.cnn )
}

CNA.object <- CNA(cbind( BL.filt.GR.comb.logratios[,-c(1:5)]), BL.filt.GR.comb.logratios$chr, 
                  BL.filt.GR.comb.logratios$start, data.type="logratio", sampleid = colnames(BL.filt.GR.comb.logratios)[-c(1:5)] )

## Segment blacklist filtered windows 
segs.BL.filt.GR.comb <- get.CNA.segments( BL.filt.GR.comb.logratios )

BL.filt.segments <- list()
BL.filt.segments[["SKMEL"]] <- segs.BL.filt.GR.comb
BL.filt.segments[["shRNA"]] <- segs.BL.filt.GR.comb
BL.filt.segments[["778"]] <- segs.BL.filt.GR.comb

## Combine segments into overlapping regions
library(CNTools)

### Function for creating reduced list of overlapping segments
get.filtered.RS <- function( segments ) {  

    CNSeg <- CNSeg( segments$output )

  ## Reduced segments
  segment.RS <- getRS( CNSeg,  by="region",imput = FALSE, XY = FALSE, what = "mean" )

  ## Median absoluted deviation filtering of outlier
  mad.filt.segment <- madFilter( segment.RS, 0.5 )
}

mad.filt.Rs.BL.filt.GR.comb  <- get.filtered.RS( segs.BL.filt.GR.comb )

### Use piecewise constant segmentation to segment into
### high-level gains and losses at the chromosomal arm level
library(copynumber)

get.combined.CNVs <- function( RS.obj ) {
  filteredrs <- RS.obj
  indx <- sapply(filteredrs@rs, is.factor)
  filteredrs@rs[indx] <- lapply(filteredrs@rs[indx], function(x) as.numeric(as.character(x)))

  chrom <- filteredrs@rs$chrom
  bp <- filteredrs@rs$start
  cnv <- filteredrs@rs[,4:length(filteredrs@rs[1,])]

  data.cnv <- cbind(filteredrs@rs$chrom, filteredrs@rs$start, filteredrs@rs[,4:length(filteredrs@rs[1,])])
  
  ### Winsorize to modulate noisy outliers
  wins.data <- winsorize(data=data.cnv, assembly = "hg19", verbose=F)

  ### Use PCF to get consensus segments and CN estimates
  pcf.segments <- pcf(data=wins.data, Y=data.cnv, assembly = "hg19")
  
  return(pcf.segments)
}

SKMEL.pcf.segs <- get.combined.CNVs(mad.filt.Rs.BL.filt.GR.comb)

Tun.778.pcf.segs <- get.combined.CNVs(mad.filt.Rs.BL.filt.GR.comb)
Tun.778.pcf.segs$sampleID <- gsub( "mcol.", "", Tun.778.pcf.segs$sampleID )

### Plot gains and losses above/below threshold as segment maps
plotAberration(segments=subset( Tun.778.pcf.segs, n.probes > 25 ), thres.gain = 0.137, thres.loss = -0.152, mar=c(5, 12, 4, 2) + 0.1, 
               main="CN profile each cell line, MAD filtered (gain/loss threshold: +/- 0.5)" )

### Plot all segments as heatmaps:
heatmap.2(trans, Rowv = FALSE, Colv = NA, distfun=function(x) dist(x, method="euclidean"), 
          hclustfun=function(x) hclust(x, method="ward.D2"), col=bluered(5), density.info="none", trace="none", 
          dendrogram=c("none"), symm=F,symkey=F,symbreaks=F, scale="none", breaks=c(min(pcf.segments$mean),-0.6,-0.3,0.3,0.6,max(pcf.segments$mean)), margins=c(2,18))
