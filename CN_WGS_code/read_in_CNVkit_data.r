### This script reads CNVkit2 output from a nominated directory,
### Filters missing data and converts to log ratios
### then saves in a list of matrices.

##cov.dir <- "CNVkit2 output/Target Coverages/SKMEL/"

for( file in list.files(cov.dir) ) { 
  infile = paste(cov.dir, file, sep="" )
 # print(infile)
  name=substr( file, 1 , 11)
 # print(name)
  ex <- read.table( infile, sep="\t", head=T)[,-c(4:5)]
  combined.cnn.values <- cbind(combined.cnn.values, name=ex$log2)
}

sample.names <- vector(mode="character")

### Get the sample names
for( file in list.files(cov.dir) ) {
   name=substr( file, 1 , 11)
   sample.names <- c( sample.names, name )
}

sample.names <- vector(mode="character")
for( file in list.files(cov.dir) ) {
  name=substr( file, 1 , 11)
  sample.names <- c( sample.names, name )
}

colnames( combined.cnn.values )[4:13] <- sample.names

### Remove rows with missing data:
trimmed.combined.cnn.values <- combined.cnn.values[ apply( combined.cnn.values[,-c(1:3)], 1, FUN = min ) > -20, ]

### Convert to logratios
logratios <- apply( trimmed.combined.cnn.values[,-c(1:3)], 2, FUN=function(vec) { vec - median(vec) } )
combined.cnn.logratios <- cbind( trimmed.combined.cnn.values[,c(1:3)], logratios )

### Save each group of logratios
Each.expt.combined.logratios[["778"]] <- combined.cnn.logratios
Each.expt.combined.logratios[["SKMEL"]] <- combined.cnn.logratios
