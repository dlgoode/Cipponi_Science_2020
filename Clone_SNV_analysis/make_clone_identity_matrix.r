### Script for making a similarity matrix between all pairs of 
### parental and drug-treated clonal lines.

### Creates a 'key' from the chr_start of each variant,
### filters variants found in repeats and at low cov sites

### Then calculate overlap with the rest of the clones.

library(GenomicRanges)
library(ape)
library(dplyr)

all.778.Vars.ID.vecs <- list()

load("RData Files/Combined_778_parental_Vars.RData")
load("RData Files/Combined_Tunicamycin_Vars.RData")

all.778.Vars.ID.vecs <- list()

Combined_778_parental_Vars <- cbind( Combined_778_parental_Vars, "Key"=paste( Combined_778_parental_Vars[,1],  Combined_778_parental_Vars[,2], sep="_" ) )
Combined_Tunicamycin_Vars <- cbind( Combined_Tunicamycin_Vars, "Key"=paste( Combined_Tunicamycin_Vars[,1],  Combined_Tunicamycin_Vars[,2], sep="_" ) )


for ( sample in unique(Combined_778_parental_Vars$Sample) ) {
  
  df1 <- subset( sample , Combined_778_parental_Vars$Sample==sample )
  
  df2 <- uniq.anns[[sample]]
  
  ### Merge annotations for variants for based on the keys
  m1 <- merge( df1, df2, "Key", all=TRUE)  
  
  ### Remove variants in repeat regions
  m1.GRanges <- GRanges(paste("chr",m1$Chr,sep=""), IRanges( m1$Pos, m1$Pos ) )
  
  gr.findOv <- findOverlaps( m1.GRanges, gr.repeats )
  rpt.vars <- queryHits( gr.findOv )
  
  m1 <- m1[ -rpt.vars, ]
  
  ### Remove variants at low coverage sites
  m1 <- subset( m1, Alt_Depth >= 2 & Total_Depth >= 10)
  ### Extract GnomAD freqs for this sample
  v <- as.num.vec( m1$GnomAD_v2.1_AF )
  
  ### Find out which of all of the variant positions across all ctrls have SNV positions in this clone, disregarding sites found in GnomAD.
  ### Creates TRUE/FALSE vector according to whether that clone has a particular SNV that was observed in 1 or more of the other clones.
  all.778.Vars.ID.vecs[[sample]] <- uniq.all.778.var.pos %in% m1$Key[ v==0 & nchar( as.character( m1$REF ) ) == 1 & nchar( as.character( m1$ALT ) ) == 1 ]
} 

all.778.ID.matrix <- data.frame(all.778.Vars.ID.vecs)

### Filter out the positions that were in GnomAD (or were indels) 
total.found.shRNA.ctrl <- apply(all.778.ID.matrix, 1, sum)

trim.all.778.ID.matrix <- all.778.ID.matrix[ total.found.shRNA.ctrl > 0 , ]

### Figure out how many shared SNV sites between samples
all.778.identity <- data.frame()

for( i in 1:5 ) {
    for( j in 1:5 ) { 
      all.778.identity[i,j] = mean( trim.all.778.ID.matrix[trim.all.778.ID.matrix[,i] ,j] )
    }
}

names( all.778.identity ) <-  names(trim.all.778.ID.matrix)
rownames( all.778.identity ) <- names(trim.all.778.ID.matrix)

### Plot neighbour joining tree based on this fraction of shared SNVs.
plot( nj(as.dist(all.778.identity)) )