#' Convenience wrapper to generate tuned and sorted barplots at various taxonomic levels from phyloseq input
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
#' @param legend default TRUE. Turn off the legend if set to FALSE
#' @param legend.col How many columns the legend may contain (default legend.col=1)
#' @param legend.size FOnt size to use in legend. Default 8 pnt.
#' @param title Any custom to plot on the figure (default generates "Barplot <taxrank> - top<x> taxa")
#' @param title.center center the title on figure otherwise left adjusted (default=TRUE)
#' @param show.other when top>0, bundle all non-top taxa into an "Other" category plotted on top of the stack (default TRUE)
#' @param other.label label to use for the bundled non-top taxa (default "Other")
#' @param other.color fill color for the "Other" stack (default "white")
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
#'          2026-04-22 alx added grouping of "other" taxa in top-x plots
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
                         title=paste0("Barplot ", taxrank, " - ", if(top>0){paste0("Top",top)}else{"all taxa"} ), title.center=TRUE,
                         show.other=TRUE, other.label="Other", other.color="white"
                       ) {
  
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
    
    # --- Compute per-sample "Other" abundance BEFORE pruning ---
    # Works for both taxglom=TRUE and taxglom=FALSE paths because the legend/fill
    # is always at `taxrank` level. We glom at taxrank (if not already) and sum
    # every taxon that is not in the top-x list, per sample.
    other_df <- NULL
    if(show.other) {
      # Use a glommed-at-taxrank version to compute "Other" cleanly.
      # If taxglom was TRUE, `ps` is already glommed and shares taxa names with topx.
      # If taxglom was FALSE, we glom `ps.org` on the fly just for the Other calc.
      if(taxglom) {
        ps.for.other <- ps
      } else {
        ps.for.other <- tax_glom2( ps.org, taxrank=taxrank, NArm=NArm )
        # topx references ps (which IS glommed inside toptaxa), so names should align
      }
      non_top <- setdiff( taxa_names(ps.for.other), topx )
      if(length(non_top) > 0) {
        otu_other <- otu_table(ps.for.other)
        # phyloseq otu_tables can be taxa-as-rows OR samples-as-rows; normalize.
        if( taxa_are_rows(ps.for.other) ) {
          other_sums <- colSums( otu_other[non_top, , drop=FALSE] )
        } else {
          other_sums <- rowSums( otu_other[, non_top, drop=FALSE] )
        }
        other_df <- data.frame( Sample = names(other_sums),
                                Abundance = as.numeric(other_sums),
                                stringsAsFactors = FALSE )
        # Attach sample metadata so facet_grid / x-axis variables resolve
        sdat <- as( sample_data(ps.for.other), "data.frame" )
        sdat$Sample <- rownames(sdat)
        other_df <- merge( other_df, sdat, by="Sample", all.x=TRUE, sort=FALSE )
        # Flag to inject the label under the `taxrank` column after plot_bar runs
        other_df[[taxrank]] <- other.label
      } else {
        # nothing to bundle
        show.other <- FALSE
      }
    }
    # --- end Other pre-compute ---
    
    ps <- prune_taxa( topx, ps)
  } else {
    # top=0 means show everything; nothing to bundle
    show.other <- FALSE
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

  # --- Inject "Other" rows into p$data so they appear as an extra stack ---
  if( show.other && !is.null(other_df) ) {
    plot_cols <- colnames(p$data)
    # Build a row template matching p$data columns exactly
    other_rows <- data.frame( matrix(NA, nrow=nrow(other_df), ncol=length(plot_cols)) )
    colnames(other_rows) <- plot_cols
    # Fill in what we have; anything we don't have (e.g. OTU, other tax ranks) stays NA
    for( cn in plot_cols ) {
      if( cn %in% colnames(other_df) ) {
        other_rows[[cn]] <- other_df[[cn]]
      }
    }
    # Make sure taxrank column explicitly carries the Other label
    other_rows[[taxrank]] <- other.label
    # Match column types where possible to avoid rbind warnings on factors
    for( cn in plot_cols ) {
      if( is.factor(p$data[[cn]]) && !is.factor(other_rows[[cn]]) ) {
        lv <- levels(p$data[[cn]])
        other_rows[[cn]] <- factor( other_rows[[cn]], levels=unique(c(lv, other_rows[[cn]])) )
      }
    }
    p$data <- rbind(p$data, other_rows)
  }
  # --- end Other injection ---

  # sort stacked bars
  if(sort.stack!="FALSE") {
    if(sort.stack=="reverse") {
      warning("reverse sort")
      p.sort <- p$data %>% group_by( eval(as.name(taxrank)) ) %>% summarise( Avg=mean(Abundance, na.rm=TRUE)) %>% arrange( desc(Avg) )
    } else {
      p.sort <- p$data %>% group_by( eval(as.name(taxrank)) ) %>% summarise( Avg=mean(Abundance, na.rm=TRUE)) %>% arrange( Avg )
    }
      colnames(p.sort) <- c(taxrank,"Avg")
      # Force "Other" to the FIRST factor level so it ends up on top of the stack.
      # ggplot2's position_stack() stacks values in REVERSE order of the group aesthetic,
      # meaning the first factor level is drawn on top, the last level at the bottom.
      if(show.other) {
        lv <- p.sort[[taxrank]]
        lv <- c( other.label, setdiff(lv, other.label) )
        p$data[taxrank] <- factor( p$data[[taxrank]], levels=lv )
      } else {
        p$data[taxrank] <- factor( p$data[[taxrank]], levels=p.sort[[taxrank]] )
      }
  } else if(show.other) {
    # even without sort.stack, ensure Other lands on top
    current_levels <- unique(as.character(p$data[[taxrank]]))
    lv <- c( other.label, setdiff(current_levels, other.label) )
    p$data[taxrank] <- factor( p$data[[taxrank]], levels=lv )
  }
  
  # generate alternating color palette
  if(top>0 & taxglom) {
    # only if no taxglom overrides to prevent not enough colors abvailable
    colpalette <- colorDistinct(n=top)
  } else {
    numcol <- length(levels(p$data[,taxrank])) # if not enough cols and NAs then +1...
    if(show.other) numcol <- numcol - 1  # reserve one slot; Other gets its own color
    colpalette <- colorDistinct(n=numcol)
  }
  # Prepend the Other color at the START so it matches the FIRST factor level
  # (Other is placed first so it stacks on top; see position_stack reverse behaviour above)
  if(show.other) {
    colpalette <- c(other.color, colpalette)
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
