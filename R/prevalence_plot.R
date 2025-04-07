#' Plot prevalence plots of certain taxonomic rank from phyloseq object
#'
#' Function will return a prevalence plot from each OTU
#' N occurences in samples (y) and total abundance (x)
#' Returns a ggplot object
#' Use the cowplot theme for great results
#'
#' @param ps a phyloseq object containing cummunity abundance data
#' @param rank string, default rank="Phylum" 
#'             Must match "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
#' @param xlog boolean, optional, default xlog=TRUE
#' @param ylog boolean, optional, default ylog=FALSE
#' @param yrelative boolean, optional, default FALSE. Calculate prevalence as relative or sample number
#' @param title string, optional, Default OTU prevalence over all samples (Rank, nsamples=xx)
#' @param strbck string, optional, default strbck=#f0f0f0, strip background color override
#' @return ggplot object
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' @author Alex Bossers \email{alex.bossers@wur.nl}
#' @note Prevalence plot based on Callahan et al 2016 (McMurdie and Holmes group); https://f1000research.com/articles/5-1492/v2
#'
#' @export
#'
#' @examples
#' prevplot0 <- prevalence_plot( ps = PhyObject, rank = "Family", xlog = TRUE, ylog = FALSE,
#' yrelative = FALSE, title = "OTU prev over all samples" )
#'
#' @import ggplot2
#' @import phyloseq

prevalence_plot <- function( ps, rank = "Phylum",
                             xlog = TRUE, ylog = FALSE,
                             yrelative = FALSE,
                             strbck = "#f0f0f0",
                             title = paste0( "Taxa (",rank,") prevalence over all samples (Nsamples=",nsamples(ps),")" ) 
                            )
{
    require(pyloseq)
    require(ggplot2)
    prevdf <- prevalence_taxa( ps, prev_relative = yrelative )  #rank variable deprecated in this function
    plotobj <- ggplot(prevdf, aes(TotalAbundance,
                                  Prevalence,
                                  color = get( rank ) ) ) +
        geom_point(size = 2, alpha = 0.4, show.legend = FALSE) +
        xlab("Total abundance per taxon (SUMtaxa)") +
        ylab("Prevalence (samples)") +
        ggtitle( title ) +
        theme(
            plot.title = element_text(size = 12, face="bold"),
            strip.text.x = element_text( size = 8, margin = margin(4, 0, 4, 0, 'pt') ),
            strip.text.y = element_text( size = 9, margin = margin(0, 4, 0, 8, 'pt') ),
            strip.background = element_rect( fill = strbck ),
            axis.text.x = element_text( size = 8, angle = 90, vjust = 0.5, hjust = 1 ),
            axis.text.y = element_text( size = 8, vjust = 0.5, hjust = 1 ) ) +
        facet_wrap( as.formula( paste0( "~", rank ) ) )
    if( ylog ) { plotobj <- plotobj + scale_y_log10() }
    if( xlog == FALSE ) {
        plotobj <- plotobj + scale_x_continuous() }
    else { plotobj <- plotobj + scale_x_log10() }
    return( plotobj )
}
