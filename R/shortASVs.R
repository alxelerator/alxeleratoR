#' Move ASV sequence as sequence_name into the phyloseq refseq slot. Also shorten the taxa names to Hashed-sequences or ASV0001 etc
#'
#' Hash method for shortning ASVs is preferred since it would allow still to merge different PS objects since hash 1:1 seq.
#'       By default hash md5 is used. Alternative just a number.
#'
#' @param ps phyloseq object (without refseq object and sequences as taxonomy names)
#' @param in.taxonomy default TRUE. Adds the actual ASV sequence as lowest level in taxonomy table
#' @param method default "hash". Hash will make a hash code of the actual sequence = unique. Alternative ASV<number>.
#' @return phyloseq object with ASVs in refseq slot and shortened ASV taxonomy names
#'
#' @export
#' @author a.bossers@uu.nl
#' @importFrom Biostrings DNAStringSet
#' @import phyloseq
#' @import openssl
#'
shortASVs <- function( ps, in.taxonomy=TRUE, method="hash" ) {
  require(Biostrings)
  require(phyloseq)
  require(openssl)

  # first check if sequences are still taxa_names. Else abort.
  if( ! all( grepl( "^[GATCNgatcn]+$", taxa_names(ps) ) ) ) {
    stop("Other characters then nucleotides (AGCTN) found in taxa_names. Are you sure the taxa_names haven't been converted before?")
  } 

  if(in.taxonomy) {
    #add ASV as lowest taxon level in taxonomy
    ps@tax_table@.Data <- cbind(ps@tax_table@.Data, "ASV"=taxa_names(ps) )
  }

  #parc ASV seq into Refseq table of ps
  asv_ref <- Biostrings::DNAStringSet( taxa_names(ps) )
  names(asv_ref) <- taxa_names(ps)
  ps <- merge_phyloseq( ps, asv_ref )

  #rename taxa names to convenient short
  if( method == "number") {
    #just number ASVs (not preferred)
    #check max number of ASVs (use minimum 4 digits)
    taxadigits <- nchar( ntaxa(ps) )
    if( taxadigits < 4 ) { taxadigits <- 4 }
    taxa_names(ps) <- paste0( "ASV", sprintf( paste0( '%0.',taxadigits,'d' ), 1:ntaxa(ps) ) )
  } else {
    # all other cases use default md5 hash method
    taxa_names(ps) <- paste0( "ASV_", openssl::md5(taxa_names(ps)) )
  }

  cat("Sequences as taxa_names converted to something more convenient...")

  return(ps)
}
