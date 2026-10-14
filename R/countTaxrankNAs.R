#' Counts the number of NAs in the taxonomy table per taxonomic rank.
#' Returns a dataframe having the NA counts and the percentage NAs of total entries per taxrank.
#' Write optionally the data to a file.
#'
#' @param ps phyloseq object
#' @param filename Optional filename and path if you want to save the output
#' @return data frame
#' @import phyloseq
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
#' @author Alex Bossers \email{a.bossers@uu.nl}
#'
countTaxrankNAs <- function( ps, filename="" )
{
  natax   <- colSums(is.na(ps@tax_table))
  tottax  <- natax + colSums(!is.na(ps@tax_table))
  perctax <- round( natax/tottax * 100, digits=2 )
  overall <- data.frame(NtaxaNAs=natax, percNAs=perctax) %>% rownames_to_column("TaxRank")
  cat( paste0("Total of ASVs = ",tottax[1],"\n"))
  if( filename > "") { 
      write.tsv(overall, paste0(filename), row.names=FALSE)
    } else {
      write.tsv(overall, "TotalNAtaxaAtTaxonomicLevelsOverall.tab", row.names=FALSE)
    }
  return(overall)
}
