#' Phylodiv: phylogenetic diversity
#'
#' @param physeq phyloseq class object, from which phylogeny and abundance data are extracted
#' @param theta parameter that determines the balance in the Balance Weighted Phylogenetic Diversity (see McCoy and Matsen, 2013)
#'              Theta = 0 corresponds to Faith's PD
#'
#' @note Requires a phylogenetic tree in your phyloseq object!
#'
#' @return dataframe
#'
#' @import phyloseq
#' @export
#'
#' @author Original author unknown. You? Contact me for credits.
#' @author Alex Bossers UU /email{a.bossers@uu.nl}

phylodiv <- function(physeq, theta = 0) {
    require(phyloseq)
    count_to_prop <- function(x) {x/sum(x)}
    physeq <- transform_sample_counts(physeq, count_to_prop)
    x <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) { x <- t(x) }
    phy <- phy_tree(physeq)

    ## Construct incidence matrix of the tree
    incidence <- incidenceMatrix(phy)

    ## Order incidence matrix according to community tables
    incidence <- incidence[colnames(x), ]

    ## Create community phylogeny matrix by multiplying (community x edge matrix)
    ## where cpm_{ij} gives the abundance of OTUs originating from branch j in community i.
    cpm <- x %*% incidence
    ## Convert to incidence matrix (0/1) and multiply by edge length to obtain PD per community.
    if (theta == 0) {
      cpm[cpm > 0] <- 1
    } else {
      cpm <- pmin(cpm^theta, (1-cpm)^theta)
    }
    pd <-  cpm %*% phy$edge.length

    ## Add sample data information
    if (!is.null(sample_data(physeq, FALSE))) {
        sdf <- as(sample_data(physeq), "data.frame")
        sdf$pd <- as.vector(pd)
        pd <- sdf
    }

    return (pd)
}
