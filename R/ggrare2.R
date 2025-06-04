#' Rarefaction curve, ggplot style
#' This will generate a rarefaction plot of the ps object
#' Custom step-size and incl standard error spread.
#'
#' Required initially to have counts as integers. Updated the function to accept rounded num.
#' This is just a more convenient warpper around ggrare allowing tweaks but also object reload without recalc
#'
#' Args:
#' @param physeq phyloseq class object, from which abundance data are extracted
#' @param step  Step size for sample size in rarefaction curves. Default 1/100 of highest samplesum.
#' @param label Default `NULL`. Character string. The name of the variable
#'              to map to text labels on the plot. Similar to color option
#'              but for plotting text.
#' @param labelsize Default 4.
#' @param color (Optional). Default NULL. Character string. The name of the
#'              variable to map to colors in the plot. This can be a sample
#'              variable (among the set returned by
#'              sample_variables(physeq) ) or taxonomic rank (among the set
#'              returned by rank_names(physeq)).
#'              Finally, The color scheme is chosen automatically by
#'              ggplot, but it can be modified afterward with an
#'              additional layer using scale_color_manual.
#' @param plot  Logical, should the graphic be plotted.
#' @param parallel should rarefaction be parallelized (using parallel framework)
#' @param se    Default FALSE. Logical. Should standard errors be computed and plotted.
#'
#' @import vegan
#' @import parallel
#' @import phyloseq
#'
#' @export
#'
#' @author Alex Bossers \email{a.bossers@uu.nl}
#'
#' @author gauravsk: The original unmodified ggrare function was scavenged from: https://rdrr.io/github/gauravsk/ranacapa/man/ggrare.html
#' 
#' Updates:
#'          2024-05-30 Alx Add auto step-size and warning at purge samples below stepsize
#'          2021-11-26 Alx Removed the replot option due to often stack errors.
#'          2021-10-18 Alx Integrated the ggrare_replot function into this one.
#'          2018-11-12 Alx Required initially to have counts as integers. Updated the function to accept rounded num.
#'          2018-03-02 Alx Added additional check if sample depth is equal or higher than supplied rarestep to prevent errors
#'                         We just prune samples equal or higher and will e EXCLUDED from the plots!
#'          2018-01-01 Alx initial checkin
#'
 
ggrare2 <- function(physeq, step = 0, label = NULL, labelsize = 4, color = NULL, plot = FALSE, parallel = TRUE, se = FALSE ) {

  	require(vegan)
    require(parallel)
    require(phyloseq)

    #determine auto step-size if none given
    if( step == 0 ) { 
        maxsamsum <- max( sample_sums(physeq) )
        step = round( maxsamsum / 100 )
    }
    cat( paste0("Step-size: ",step,"\n") )

    #cleanup any samples under stepsize.
    ps_nstart <- nsamples(physeq)
    physeq <- prune_samples( sample_sums(physeq)>=step, physeq )
    ps_nend <- nsamples(physeq)
    if( ps_nend < ps_nstart) { 
        warning(paste0("Samples (N=",(ps_nend-ps_nstart),") having sample_sum below step-size were purged before rarefaction!") ) 
    }

    x <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) { x <- t(x) }

    ## This script is adapted from vegan `rarecurve` function
    tot <- rowSums(x)
    S <- rowSums(x > 0)
    nr <- nrow(x)

    rarefun <- function(i) {
        cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i]) {
            n <- c(n, tot[i])
        }
        y <- rarefy( x[i, ,drop=FALSE], n, se = se)
        if (nrow(y) != 1) {
            rownames(y) <- c(".S", ".se")
            return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
        } else {
            return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
        }
    }
    if (parallel) {
        out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
    } else {
        out <- lapply(seq_len(nr), rarefun)
    }
    df <- do.call(rbind, out)
    
    ## Get sample data
    if (!is.null(sample_data(physeq, FALSE))) {
        sdf <- as(sample_data(physeq), "data.frame")
        sdf$Sample <- rownames(sdf)
        data <- merge(df, sdf, by = "Sample")
        labels <- data.frame(x = tot, y = S, Sample = rownames(x))
        labels <- merge(labels, sdf, by = "Sample")
    } else {
        stop("NULL sample data.")
    }
    
    ## Add, any custom-supplied plot-mapped variables
    if( length(color) > 1 ){
        data$color <- color
        names(data)[names(data)=="color"] <- deparse(substitute(color))
        color <- deparse(substitute(color))
    }
    if( length(label) > 1 ){
        labels$label <- label
        names(labels)[names(labels)=="label"] <- deparse(substitute(label))
        label <- deparse(substitute(label))
    }
    
    p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
    p <- p + labs(x = "Sample Size", y = "Species Richness")
    if (!is.null(label)) {
        p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                           size = labelsize, hjust = 0)
    }
    p <- p + geom_line()
    if (se) { ## add standard error if available
        p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
    }
    if (plot) {
        plot(p)
    }
    invisible(p)
}

