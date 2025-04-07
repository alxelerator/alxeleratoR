#' write table using tabular tab seperators
#' @param dataToWrite dataframe
#' @param filename filename and path
#' @param sep Field separator to use. Default is tab.
#' @param quote Boolean default FLASE. Quote fields using ""
#' @param row.names Boolean, Add row names. Default TRUE
#' @export
#' @author Alex Bossers \email{a.bossers@uu.nl}
write.tsv <- function( dataToWrite, filename, sep = "\t", quote = FALSE, row.names=TRUE )
{
	write.table( dataToWrite, filename, sep=sep, quote=quote, row.names=row.names )
}

