#' Alx: shorten too long temp OTU names
#' @author a.bossers@uu.nl
#' @author scavenged from / inspired by https://peerj.com/articles/545/#fig-1
#'
#' @param original The taxon name to shorten as string
#' @return shortened string
#'
#' @export
#'
short.tn <- function(original)
{
    # shorten long taxanames
    # new refs depend on identification stage OTU picking
    # https://peerj.com/articles/545/#fig-1
    tn <- paste("OTU_", original, sep = "")
    tn <- sub("OTU_New.CleanUp.ReferenceOTU", "NCR_OTU_", tn)
    tn <- sub("OTU_New.ReferenceOTU", "NR_OTU_", tn)
    return(tn)
}
