parental_samples <- as.character( unique(Combined_778_parental_Vars$Sample) )
Tunicamycin_samples <- as.character( unique(Combined_Tunicamycin_Vars$Sample) )

### Script for creating a list of private variants for each sample
### from a given set of clonal cell lines

### In this example we are filtering variants called from
### the shRNA control sets
### The same operations were performed on clones from the 
### 778 parental and tunicamycin treated lines
shRNA_ctrl_samples <- as.character( unique( shRNA_ctrl_Vars$Sample ) )

# To store data on somatic mutations for each line:
som.mutns <- list()

# To store annotations (so we don't have to reread in if there's a bug)
uniq.anns <- list()

# Directory that has the annotations:
annot.dir <- "/home/dgoode/Cell Line WGS/Annotated_vars/"

sample_list <- shRNA_ctrl_samples

for (sample in sample_list) {
  
  ### Pare down to variants potentially private to that sample:
  singles <-  subset( shRNA_ctrl_Vars, Sample==sample )
  singles <-  subset( singles, `Num samples with Var`==1 )
  
  ### Get the sample annotations
  annot.file <- list.files(annot.dir, patt=sample )
  print( sample, length(annot.file) )
  
  annotations <- read.table( paste(annot.dir,annot.file,sep=""))
   colnames(annotations) <- Annotation_col_header

  uniq_ann <- unique( annotations[, c(1:5,13,14) ] )
  uniq_ann <- cbind( uniq_ann, "Key"=paste( uniq_ann[,1], uniq_ann[,2], sep="_" ))
  #colnames(uniq_ann) <- c( Annotation_col_header[c(1:5,13,14)], "Key")
  uniq.anns[[sample]] <- uniq_ann
 # uniq_ann <- uniq.anns[[sample]]
  
  ### Cross reference against the treated line to remove
  ### Variants seen in any clones in the treated group
  m <- singles[ !singles$Key %in% FRAP1_shRNA_Vars$Key, ]
  
  m2 <- merge( m, uniq_ann )
  
  print( dim(m2) )
  
  som.mutns[[sample]] <- m2[ m2$GnomAD_v2.1_AF==0, ]  #& m2$`Alt Freq` < 0.15
}


> save( som.mutns, uniq.anns, file="som.mutns.uniq.anns.each.sample.RData")

save.image()

## How many low-freq SNPs altogether?
> unlist( lapply( som.mutns, function(m2) { sum(m2$GnomAD_v2.1_AF==0, na.rm=T)} ) )
FD03043124 FD03043115 FD03043116 FD03043123 FD03043117 
8522       8515       9024       7655       7346

> unlist( lapply(som.mutns, nrow) )
FD03043124 FD03043115 FD03043116 FD03043123 FD03043117 
2657       2533       2704       2340       2199 

## How many indels? 
unlist( lapply(som.mutns,function(m2) { sum( nchar( as.character( m2$REF ) ) - nchar( as.character( m2$ALT )) < 0 , na.rm=TRUE ) } ) )
FD03043124 FD03043115 FD03043116 FD03043123 FD03043117 
315        292        309        298        288 