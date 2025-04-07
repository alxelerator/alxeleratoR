#' Make a ggplot of sample sums per sample ordered decreasing by sample sum and indicate 90% rarefy number
#' Only makes a plot for inspection and play with parameters. No rarefaction or filtering is done
#'
#' @param ps phyloseq object
#' @param percentKeep default 90. Percentage of samples to keep at rarefaction cutoff indication
#' @param cutoff default FALSE. Overrides percentKeep and sets the minimum sampleSum to take.
#' @param numberKeep default FALSE. Override percentKeep and cutoff and sets the number of samples to set rarefy threshold to
#' @param logy default TRUE. Make y-scale log10
#' @param color optional metadata variable name to color the points in the graph
#' @param size.points default=3. Actual point size for the ggplot.
#' @param stats default TRUE. Show median, max and min in top right corner of plot.
#' @param subtitle default none. gplot subtitle string
#' @param crosshair default TRUE. Show the bluegrey crosshair target on selected cutoff.
#' @param xlab.rel Scale the x-axis labels to something readible (default scaling text 0.8x)
#' @param namesColumn default SampleNames. Custom metadata column to use as sample names along x-axis.
#'
#' @return ggplot_object
#'
#' @author \email{a.bossers@uu.nl}
#'
#' @import phyloseq
#' @import ggplot2
#' @import scales

sample_sum_plot <- function( ps=merged, percentKeep = 90, cutoff = FALSE, numberKeep = FALSE, 
                             logy=TRUE, color="", size.points=3,
                             stats=TRUE, subtitle="", crosshair=TRUE, xlab.rel = 0.8, namesColumn="" )
{
  require(phyloseq)
  require(ggplot2)
  require(scales)

  samsums <- data.frame( samsums = sample_sums(ps) )
  samsums.median <- median( samsums$samsums )
  samsums.max <- max(samsums$samsums)
  samsums.min <- min(samsums$samsums)

  #add potential metadata columns
  if( color > "" ) {
    samsums[,color] <- as.matrix(ps@sam_data[,color])
  }
  if( namesColumn > "" ) {
    samsums$samples <- as.matrix(ps@sam_data[,namesColumn])
    #row.names(m) <- NULL
    #row.names(m) <- m[,namesColumn]
    row.names(samsums) <- samsums$samples
    } else {
      samsums$samples <- row.names(samsums)
    }
  
  samsums <- samsums[ order(-samsums$samsums), ]
  #order them on x axis based on value y
  samsums$samples <- factor(samsums$samples, levels = samsums$samples[order(-samsums$samsums)])

  #default 90p
  n <- nsamples( ps )
  rarefysample <- floor( length(samsums$samples) / 100 * percentKeep )
  rarefy90p <- samsums[ rarefysample, ]$samsums
  
  # set samplesum cutoff if used
  if( cutoff > 0 ) {
    rarefy90p = cutoff
    #samplenumber at cutoff
    rarefysample <- sum( samsums$samsums > cutoff )
    percentKeep <- ceiling( rarefysample / n * 100 )
  }
  
  # rarefaction values
  if( numberKeep > 0 ) {
    #if number then we calculate the percentKeep
    percentKeep <- ceiling( numberKeep / n * 100 ) # with round on lower numbers gives .33333 then the last sample is missing.
    rarefysample <- numberKeep #floor( length(samsum.data$samples) / 100 * percentKeep )
    rarefy90p <- samsums[ numberKeep, ]$samsums
  }

    samsum.plot <- ggplot( samsums, aes( x=samples, y=samsums) ) +
      xlab("Samples (decreasing samplesum)") +
      theme( axis.text.x = element_text(angle = 270, hjust = 0, vjust=0.2, size=7),
             axis.text.y = element_text(size=7)) +
      ylab("SampleSum") +
      theme( panel.grid.major.y = element_line(colour = "lightgrey"),
             plot.title = element_text( hjust=0.5 ),
             plot.subtitle = element_text( hjust=0.5 )) +
    #if(subtitle) {
    #  ggtitle( "Total observations per sample",
    #         subtitle=paste0( "All ",nsamples(ps)," samples" ) )
    #} else {
    if( subtitle>"" & ! subtitle == FALSE) {
      ggtitle("Total observations per sample", subtitle=subtitle)
    } else {
      ggtitle( "Total observations per sample", subtitle=NULL )
    }
    #}
    
    if( color > "" ) {
      colors <- samsums[,color]       #ps@sam_data[,color]
      row.names(colors) <- NULL
      samsum.plot <- samsum.plot + 
        geom_point( aes( colour=factor( colors ) ), size=size.points ) + 
        scale_colour_discrete( colors ) +
        guides( col=guide_legend(color)) +
        theme( legend.title=element_text(size=9),
               legend.text=element_text(size=8))
    } else {
      samsum.plot <- samsum.plot + geom_point(size=size.points)
    }
    
  if(logy) {
    samsum.plot <- samsum.plot +
                      ylab("SampleSum") +
                      scale_y_log10(label=comma) #scale_y_continuous(trans ="log1p") #no 0 problem but less weel formed y axis
  } else { 
    samsum.plot <- samsum.plot + scale_y_continuous(label=comma) 
  }

  if(stats==TRUE) {
    samsum.plot <- samsum.plot  + coord_cartesian(clip = "off") +
                                  annotate("text", label=paste0("\nTotal = ", nsamples(ps),
                                                               "\nKeep ",percentKeep,"% = ", rarefysample,
                                                               "\nKeep at = ", format( rarefy90p, big.mark="," ),
                                                               "\nMedian = ", format(samsums.median,big.mark = ","),
                                                               "\nMax = ", format(samsums.max,big.mark = ","),
                                                               "\nMin = ", format(samsums.min,big.mark = ",") ),
                                          x=Inf, y=Inf, hjust="right", vjust="top", size=3, color="darkblue") 
    #if(subtitle) {
    #  samsum.plot <- samsum.plot + ggtitle("Total observations per sample",
    #                                      subtitle=paste0( "To keep ",percentKeep,"% (",rarefysample,") of all (n=",nsamples(ps),") use cutoff ", format( rarefy90p, big.mark="," ) ))
    #} else {
    #  samsum.plot <- samsum.plot + ggtitle("Total observations per sample" )
    #}
    
    if(crosshair) {
      samsum.plot <- samsum.plot + geom_hline( yintercept = rarefy90p, color="blue", alpha=0.2) +
                                   geom_vline( xintercept = rarefysample, color="blue", alpha=0.2 )
    }
  }

  samsum.plot <- samsum.plot + theme( axis.title=element_text(size=7,face="bold"),  
                                      axis.text.x = element_text( size=rel(xlab.rel) ) )
  
  return(samsum.plot)
}

