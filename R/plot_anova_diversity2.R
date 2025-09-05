#' Convenience wrapper to ggplot alpha diversity combined with pair-wise ANOVA using a (rarefied) phyloseq object as a starter. Classic or Welch ANOVA supported.
#'
#' @param physeq phyloseq object of microbiome data
#' @param method methods can be c("richness", "fisher", "simpson", "invsimpson", "shannon", "evenness", "pd")
#' @param rarefy.override Richness needs rarefied data. It automatically detects and rarefies to lowest = default. Override with this setting.
#' @param grouping_column Metadata column to use in the grouping
#' @param sample_order ordered vector of exact sample/grouping names. Default order of appearance.
#' @param anova.type Defaults to classic ANOVA ("anova"). Option to use Welch-ANOVA ("welch").
#' @param pValueCutoff default 0.05
#' @param type plot type box/violin (default type="violin") else individual points will be shown.
#' @param whisker will optionally add whisker boxplot layer to figure. Default = FALSE.
#' @param jitter will toggle jitter (random spread 10% left and 10% right). Default = TRUE.
#' @param alpha plotting dot transparency (default 0.3)
#' @param x.angle angle under which x-axis category labels will be printed. Default = 0.
#' @param title any title. Default empty.
#' @param outfile Optional file name (and path) where to store the data (default "" not saved)
#'
#' @returns list of $plot = ggplot figure, $anova = pair-wise ANOVA table, $df = underlying data.frame (long format), Optional tabular data file on disk if `outfile` is given.
#'
#' @import phyloseq
#' @import dplyr
#' @import rstatix
#' @export
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' @note
#'  See also MicrobiomeSeq package
#'  The current function is a modified and extended (drop-in replacement) version of the original `plot_anova_diversity` from the `MicrobiomeSeq` package by Alfred Ssekagiri.
#' 
#' @section
#' Updates:
#' * 2023-07-05 alx added x.angle to rotate x-axis legend labels
#' * 2022-10-12 alx added Welch-ANOVA options for unequal variances
#' * 2022-10-02 alx tweaked box/violin and added `whisker` and `jitter` options
#' * 2022-07-07 alx fixed sample order option to move along significance bars!
#' * 2022-07-01 Alx NOT FIXED YET! extended to include ordering of panels (methods) and samples/grouping avoiding hassling around in ggplot
#' * 2022-06-08 Alx extend invSimpson
#' * 2019-11-03 Alx Added Type can be violin or box (plot type)
#' * 2019-11-02 Alx made error bars more compact using 10% plot in addition if 10 bars are plotted.
#' * 2019-11-01 Alx Also moved over the `outfile` option to save table
#' 

plot_anova_diversity2 <- function( physeq, method, rarefy.override=0, 
                                   grouping_column, sample_order="",
                                   anova.type="anova", pValueCutoff = 0.05,
                                   type="violin", whisker=FALSE, jitter=TRUE, alpha=0.3, x.angle=0, 
                                   title="", outfile="")
{
  require(phyloseq)
  require(dplyr)
  require(rstatix) #function for Welch ANOVA
  

  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  abund_table <- otu_table(physeq)
  meta_table <- sample_data(physeq)
  div.df <- alpha_div(physeq, method, rarefy.override=rarefy.override)
  df <- data.frame(div.df, (meta_table[, grouping_column])[as.character(div.df$sample),])
  
  #here we need our helper function to perform pair-wise ANOVA
  anova_res <- perform_anova( df, meta_table, grouping_column,
                             pValueCutoff, anova.type=anova.type )
  df_pw <- anova_res$df_pw

  #conditional save csv file
  if ( outfile > "" ) {
    write.csv(df_pw, file = outfile)
  }

  # remove any NAs from i.e. Pielou's and warn before plotting
  if( any(is.na(df$value) ) ) {
    #warn and remove NA values
    warning( paste0("NA values found and removed before plotting from measure: ", paste0( unique(df[is.na(df$value),]$measure) ), sep="\n" ) )
    df <- df[ (! is.na(df$value) ), ]
  }
  
  # set measure order for panels
  df$measure <- factor( df$measure, levels=c("Richness","Shannon","Pielou's evenness") )
  
  p <- ggplot( aes_string(x = grouping_column, y = "value",
                         color = grouping_column), data = df )

  if(type=="violin") {
    if(whisker) {
      p <- p + geom_violin() + geom_boxplot(width=0.2)
    } else {
      p <- p + geom_violin()
    }
  } else if(type=="box") {
    p <- p + geom_boxplot()
  } 
  
  if(jitter){
    p <- p + geom_jitter( width=0.1, alpha=alpha) 
  } else {
    p <- p + geom_point( alpha=alpha )
  }
  
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text( angle=x.angle, hjust=0.5 ))
  
  
  p <- p + facet_wrap( ~ measure, scales = "free_y", nrow = 1)
  p <- p + ylab("Observed Values") + xlab("Samples")
  
  p <- p + theme( strip.background = element_rect(fill = "white") ) +
    xlab("Groups")
  
  # check if custom sample order is given over x-axis
  # we need to do this WITHIN the function since doing the same on the final object
  # will place the significance bars between the wrong samples!
  if(length(sample_order)>1){
    p$data <- p$data %>%
      mutate( across( grouping_column, factor, levels=sample_order ) )
  }

  if (!is.null(df_pw)) {
    for (i in 1:dim(df_pw)[1]) {
      p <-
        p + geom_path(
          inherit.aes = F,
          aes(x, y),
          data = data.frame(
            x = c(
              which(levels(as.factor(p$data[, grouping_column])) == as.character(df_pw[i, "from"])),
              which(levels(as.factor(p$data[, grouping_column])) == as.character(df_pw[i, "to"])) ),
            y = c(as.numeric(as.character(df_pw[i,"y"])), as.numeric(as.character(df_pw[i, "y"]))),
           measure = c(as.character(df_pw[i, "measure"]),
                        as.character(df_pw[i, "measure"]))
          ),
          color = "black",
          lineend = "round",
          arrow = arrow(
            angle = 90,
            ends = "both",
            length = unit(0.02, "inches")
          )
        )
      p <- p + geom_text(
        inherit.aes = F,
        aes(x = x, y = y,
            label = label, vjust=0.7 ),
        data = data.frame(
          x = (
            which(levels(as.factor(df[,grouping_column])) == as.character(df_pw[i, "from"])) +
              which(levels(as.factor(df[, grouping_column])) == as.character(df_pw[i,"to"])) ) / 2,
          y = as.numeric(as.character(df_pw[i,"y"])),
          measure = as.character(df_pw[i, "measure"]),
          label = as.character(cut(
            as.numeric(as.character(df_pw[i,"p"])),
            breaks = c(-Inf , 0.001, 0.01, 0.05, Inf),
            label  = c(   "***", "**" , "*" , "")
          ))
        )
      )
    }
  }
  if(anova.type=="welch") {
    anova.caption="Welch-ANOVA *** P<0.001, ** P<0.01, * P<0.05"
  } else {
    anova.caption = "ANOVA *** P<0.001, ** P<0.01, * P<0.05"
  }
  p <- p + 
    labs( caption=anova.caption ) +
    theme( plot.caption=element_text( face="italic", color="darkblue", size=8 ) )
  
  if(title>"") { p <- p + ggtitle(title) }
  
  #turn off guides
  p <- p + guides( fill="none" , color="none" )

  return( list( plot=p, anova=anova_res$df_pw, df=df ) )
}
