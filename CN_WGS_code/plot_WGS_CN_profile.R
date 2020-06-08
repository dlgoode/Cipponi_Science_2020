### For plotting CN profile from CNVkit output from WGS data from clonal cell lines.
### I.e., CN profile plots in Fig S7.

file_suffix = "WGS CN profile.pdf"

for (j in 3:12 ) {
  sample_num = j

  sample = names(df)[sample_num]

  file_name = paste( sample, file_suffix )

  pdf( file_name, width=12, height=8.5)

  plot(c(0:3), pch=".", cex=2, frame=F, las=1, main=sample,
     ylab="LogR", xlim=c(0,nrow(df)), xaxt = "n", xlab="Genome position (hg19)",
     col="white", ylim=c(-4.3,4.3) )

  for( j in seq(1,22,2) ) { 
    index = df$chrom==j
    points( x=as.numeric( rownames(df) )[index], 
          y=df[ index, sample_num], 
          pch=".", cex=2, col="black" ) 
  }

  for( j in seq(2,22,2) ) { 
    index = df$chrom==j
    points( x=as.numeric( rownames(df) )[index], 
          y=df[ index, sample_num], 
          pch=".", cex=2, col="blue" ) 
  }

  abline(h=0, col="green")

  chr_boundaries <- tapply( as.numeric( rownames(df) ), INDEX=df$chrom, FUN=max)
  chr_boundaries <- chr_boundaries[ order(as.num.vec(names(chr_boundaries)))]

  for( b in c(1, chr_boundaries) ) { abline( v=b, col="grey", lty=3 ) }

  chr_middles <- tapply( as.numeric( rownames(df) ), INDEX=df$chrom, FUN=median)
  chr_middles <- chr_middles[ order(as.num.vec(names(chr_middles)))]


  for( i in 1:length(chr_middles) ) { 
   m=chr_middles[i]; 
    text( m, -4.1, i, cex=0.75, col="red")
  }

  dev.off()

}  