#' Rarefy each sample in a phyloseq object to its own target depth
#'
#' Rarefies each sample to an individually specified depth rather than a single
#' shared value. Useful when correcting for input volume differences in
#' semi-quantitative data (e.g. resistome FPK/FPKM), where relative abundances
#' are not the desired output.
#'
#' @param physeq phyloseq object containing at least 2 samples
#' @param rarefy_depth integer vector of target read counts, one per sample,
#'   in the same order as sample_names(physeq)
#' @param seed random seed for reproducibility (default 2202)
#'
#' @return phyloseq object with each sample rarefied to its target depth.
#'   All other slots (tax_table, sample_data, phy_tree, refseq) are preserved
#'   from the input object.
#'
#' @author a.bossers@@uu.nl
#'
#' @import phyloseq
#' @export
#'
#' @note
#'   version: 2.0
#'   versiondate: 2026-07-01
#'   Rewrote loop using lapply + matrix combine to avoid S4 scoping issue
#'   where iterative merge_phyloseq inside a function only returned the last sample.

rarefy_custom_depth <- function(physeq, rarefy_depth, seed = 2202) {

  nsamples <- length(sample_names(physeq))

  if (nsamples != length(rarefy_depth)) {
    stop("rarefy_depth file must have the same number of samples as found in the phyloseq object.")
  }

  # Rarefy each sample individually into a list
  rarefied <- lapply(seq_len(nsamples), function(i) {
    s   <- sample_names(physeq)[i]
    ps_i <- subset_samples(physeq, sample_names(physeq) == s)
    rarefy_even_depth(ps_i, sample.size = rarefy_depth[i],
                      trimOTUs = FALSE, rngseed = seed)
  })

  # Combine OTU matrices respecting orientation of the input object
  if (taxa_are_rows(physeq)) {
    otu_mat <- do.call(cbind, lapply(rarefied, function(ps) as(otu_table(ps), "matrix")))
  } else {
    otu_mat <- do.call(rbind, lapply(rarefied, function(ps) as(otu_table(ps), "matrix")))
  }

  # Replace only the OTU table; all other slots are preserved from physeq
  otu_table(physeq) <- otu_table(otu_mat, taxa_are_rows = taxa_are_rows(physeq))

  return(physeq)
}
