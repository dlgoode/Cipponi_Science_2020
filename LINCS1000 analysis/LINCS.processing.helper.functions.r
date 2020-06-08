### Function for splitting an array on the specified dividing char
### and outputting as an array
string.to.array <- function( string, divider ) { 
     unlist( strsplit(string, divider) )
}

### Function to extract the cell line and time point from a LINCS
### expt name, based on the LINCS formatting
condense.LINCS.name <- function( name ) {
  
  paste( string.to.array( string.to.array( name, ":")[1], "_" )[2:3],  collapse = "_" )
}

### Example usage:
as.vector( vapply( colnames( p@mat[ grepl("SMA",p@rdesc$pr_gene_symbol), grepl("tangeritin",p@cdesc$pert_iname) ] ), FUN=condense.LINCS.name, FUN.VAL="No Name" ) )

as.num.vec <- function( vector ) { as.numeric(as.vector( vector )) }