### Script to make phylogenetic trees based on correlations in
### copy-number profiles.
library(ape)

rs(mad.filt.Rs.BL.filt.GR.comb)

as.num.rs <- apply( rs(mad.filt.Rs.BL.filt.GR.comb), 2, FUN=as.num.vec )

cor.mat.778.segments.BL.filt <- cor( apply( rs(mad.filt.Rs.BL.filt.GR.comb)[,-c(1:3)], 2, FUN=as.num.vec))
plot( nj( dist(cor.mat.778.segments.BL.filt ) ) )

plot(cor.mat.778_nj_tree)
edgelabels( round(cor.mat.778_nj_tree$edge.length,3) )
