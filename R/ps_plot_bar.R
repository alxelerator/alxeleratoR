#' A convenience wrapper to generate tuned barplots at various taxonomic levels using phyloseq input
#'
#' @param ps phyloseq object
#' @param taxrank taxonomic rank to show (default "Phylum")
#' @param top top-x to show (default top=10)
#' @param taxglom Agglomerate the taxonomic data at the indicated level. Default TRUE
#' @param taxglom.next Show subdivision of the bars af a next lower level of taxonomic rank. Default = ""
#' @param sort.stack whether or not to sort the actual stacking on overall abundance (default sort.stack=TRUE)
#' @param x what to show on x-axis (default x="Sample")
#' @param y what to show on y-axis (default y="Abundance")
#' @param group metadata grouping variable to "merge_samples" phyloseq function
#' @param group_fun default "mean". phyloseq merge_samples will SUM otu_table by default. Here we corrected for the Nsamples per group -> actual MEAN. All other values will result in default behaviour.
#' @param logy transform Abundance y-scale to log10 (default FALSE)
#' @param facet_grid we want to facet the stacked bars (default facet_grid=NULL, example facet_grid="Kingdom")
#' @param NArm whether or not to remove NAs in the data while agglomerating at given taxonomic rank (default NArm=FALSE)
#' @param xlab.rel Tweaked size of x-axis labels since readability depends on N samples (default xlab.rel=0.5)
#' @param legend.col How many columns the legend may contain (default legend.col=1)
#' @param legend.size FOnt size to use in legend. Default 8 pnt.
#' @param title Any custom to plot on the figure (default generates "Barplot <taxrank> - top<x> taxa")
#' @param title.center center the title on figure otherwise left adjusted (default=TRUE)
#'
#' @returns ggplot figure
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#'
#' @import RColorBrewer
#' @import pals
#' @import phyloseq
#' @import dplyr
#' @export
#' @note
#'   Updates:
#'          2022-07-01 alx added taxglom options to allow visible subtaxonomic ranks per bar.
#'          2022-05-15 alx added grouping and MEAN grouping options
#'          2021-xx-xx alx many small tweaks. See git history.
#'          2020-03-02 Alx split bucket into functions
#'          2020-01-01 Alx initial checkin
#' 


#' FixMe: sorting abundance of stacks when not using taxglom showing subtaxonomic ranks. Sorting is different.


ps_plot_bar <- function( ps, 
                         taxrank="Phylum", top=10, taxglom=TRUE, taxglom.next="",
                         sort.stack=TRUE, 
                         x="Sample", y="Abundance", group="", group_fun="mean",
                         logy=FALSE, 
                         facet_grid=NULL, NArm=FALSE, 
                         xlab.rel = 0.5,   
                         legend=TRUE, legend.col=1, legend.size=8,
                         title=paste0("Barplot ", taxrank, " - ", if(top>0){paste0("Top",top)}else{"all taxa"} ), title.center=TRUE 
                       ) {
  
  require(RColorBrewer)
  require(pals)
  require(dplyr)
  require(phyloseq)

  # other custom functions from my toolbox
  #source("tax_glom2.R")
  #source("toptaxa.R")
  #source("colorDistinct.R")


  # first drop zero samplesum samples (throws error otherwise)
  zerosamples <- sum( sample_sums(ps) == 0 )
  if( zerosamples > 0 ) {
    warning( paste0( "Samples with zero sample_sums detected (N=", zerosamples, "). These are removed before plotting!\n" ) )
    # note subset_samples doesn't work in user functions! See Rplayground for working example demonstrating this.
    ps <- prune_samples( sample_sums(ps) > 0, ps )  
   }
  
  # Group samples if grouping variable given
  if(group > "") {
    
    # first prevent NAs in metadata variable
    if( sum( is.na( sample_data(ps)[,group] ) ) > 0 ) {
      stop("Metadata variable column contains NAs which may cause strange results. Aborting!")
    }
    ps.tmp <- merge_samples( ps, group=group )
    if( group_fun == "mean" ) {
      # otu table goes default to SUM for otu table. I want to show actual MEAN so correct for sample numbers
      snumbers <- sample_data(ps) %>% group_by_at( group ) %>% count()
      #double check if sample order is the same as metadata count order.
      ps.o <- otu_table(ps.tmp)
      if( ! all(snumbers[,1] == row.names(ps.o) ) )  {
        stop("Internal error occurred. Sample name orders doesn't seem to match order of samples_count per group somehow... Turn off group_fun='mean'.")
      } else {
        ps.o <- ps.o / snumbers$n
        otu_table(ps.tmp) <- ps.o
      }
    } else if( group_fun == "median" ) {
      # otu table goes default to SUM for otu table. I want to show actual MEAN so correct for sample numbers
      snumbers <- sample_data(ps) %>% group_by_at( group ) %>% count()
      stop("Median function not implemented yet!")
    } else {
      warning("Grouping function defaulting to SUM!")
    }
    ps <- ps.tmp
  }
  
  #glom at correct taxon level
  if(!taxglom) { 
      #keep in case taxglom=FALSE
      ps.org <- ps 
    } else { 
      ps.org <- ps 
      ps <- tax_glom2( ps, taxrank=taxrank, NArm=NArm )
    }
  
  if(top>0) {
    # make first top x at taxrank
    topx <- toptaxa( ps, rank=taxrank, top=top )
    ps <- prune_taxa( topx, ps)
  }
  
  if(!taxglom) {
    # get taxa list op taxrank level to extract later if taxglom = FALSE. 
    # Using taxa as is would return only a subset of representative ASVs so we need to select taxa names at higher level first to know exact ASVs.
    taxtab <- data.frame(tax_table(ps))
    taxtop <- taxtab[topx,]
    topPerRank <- taxtop[,taxrank] #taxa at i.e. Family level
    
    #Now we determine whether we want to taxglom deeper level or not
    #Else we leave data at ASV level.
    ps.org2 <- prune_taxa( data.frame(tax_table(ps.org))[,taxrank] %in% topPerRank, ps.org )#subset_taxa( ps.org, eval(as.name(taxrank)) %in% topPerRank)
    if(taxglom.next != "") {
      ps.org2 <- tax_glom2( ps.org2, taxrank=taxglom.next )
    }
      ps <- ps.org2
  }
  
  p <- plot_bar( ps, x=x, y=y, fill=taxrank, title=title, facet_grid=facet_grid )

  # sort stacked bars
  if(sort.stack!="FALSE") {
    if(sort.stack=="reverse") {
      warning("reverse sort")
      p.sort <- p$data %>% group_by( eval(as.name(taxrank)) ) %>% summarise( Avg=mean(Abundance)) %>% arrange( desc(Avg) )
    } else {
      p.sort <- p$data %>% group_by( eval(as.name(taxrank)) ) %>% summarise( Avg=mean(Abundance)) %>% arrange( Avg )
    }
      colnames(p.sort) <- c(taxrank,"Avg")
      p.sort
      p$data[taxrank] <- factor( p$data[[taxrank]], levels=p.sort[[taxrank]] )
  }
  
  # generate alternating color palette
  if(top>0 & taxglom) {
    # only if no taxglom overrides to prevent not enough colors abvailable
    colpalette <- colorDistinct(n=top)
  } else {
    numcol <- length(levels(p$data[,taxrank])) # if not enough cols and NAs then +1...
    colpalette <- colorDistinct(n=numcol)
  }
  p <- p + scale_fill_manual(values = colpalette ) +
    theme( text = element_text(size=10),
           axis.text.y = element_text( size=8 ) , 
           axis.text.x = element_text( angle=270, vjust=0.1, hjust=0, size=8 ),
           legend.title = element_text( size=legend.size, face="bold"),
           legend.text = element_text( size=legend.size )) +
    scale_y_continuous(expand = c(0, 0))
  
  if(legend.col) {
    p <- p + guides( fill=guide_legend(ncol=legend.col) )
  }
  
  if( ! legend ){
    p <- p + guides( fill="none" )
  }
  
  if(!taxglom) {
    #add subtitle when different taxglom is used
    if(taxglom.next == ""){
      sub <- paste0("(ASV-level taxa stacked per ",taxrank,"-taxonomic rank)")
    } else {
      sub <- paste0("(",taxglom.next,"-level taxa stacked per ",taxrank,"-taxonomic rank)")
    }
    p <- p + labs( subtitle=sub )
  }
  
  if(title.center) {
    p <- p + theme( plot.title = element_text( hjust = 0.5 ),
                    plot.subtitle = element_text( hjust = 0.5 ))
  }
  
  p <- p + theme( axis.text.x = element_text( size=rel(xlab.rel) ) )
  
  # if logy then transform y axis and label
  if(logy){
    p <- p + scale_y_continuous(trans='log10') + ylab("Log10(Abundance)")
  }


  return(p)
}
