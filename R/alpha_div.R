#' scavenged from microbiome seq to allow tweaking
#'
#' @param physeq phyloseq with microbiome data
#' @param method string or vector of strings and can be one or more from ("richness", "fisher", "simpson", "invsimpson", shannon", "evenness", "pd")
#' @param rarefy.override Richness needs rarefied data. It automatically detects and rarefies to lowest = default. Override with this setting.
#'
#' @returns dataframe
#' @import phyloseq
#' @import vegan
#' @export
#' @author Alex Bossers \email{a.bossers@uu.nl}
#' 
#' @note
#' Updates:
#'          2023-07-05 Alx fixed error on abundance table vs df in richness
#'          2022-09-15 Alx fixed check on taxa in rows otherwise it will return rubbish without knowing.
#'          2020-03-02 Alx split bucket into functions
#'          2020-01-01 Alx initial checkin
#' 

alpha_div <- function (physeq, method, rarefy.override=0 )
{
  require(phyloseq)
  require(vegan)

  method <- match.arg(method, c("richness", "fisher", "simpson", "invsimpson",
                                "shannon", "evenness", "pd"), several.ok = TRUE)
  
  #check to make sure taxa are in columns
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  
  abund_table <- data.frame(otu_table(physeq))
  df <- NULL
  
  if ("richness" %in% method) {
    rs <- rowSums(abund_table)
    if( ! all( rs == rs[1] ) ) {
      #doubled up but I want ALWAYS the warning to be shown while a "stop" prevents this showing!
      warning(paste0("Richness needs rarefied data and your data seems to have unequal samplesums!\n",
                     "Automatic rarefying will subsample to ",min(rowSums(abund_table)), 
                     ". If you do not want that, use the rarefy.override parameter to set your own value!\n"))
    }
    
    if ( rarefy.override == 0 ) {
        if( min(rowSums(abund_table) ) <10 ) {
          stop("Subsampling below 10 doesn't make any sense!\n")
        } else {
            R <- vegan::rarefy(abund_table, min(rowSums(abund_table)))
        }
      } else {
        warning(paste0("Using rarefy.override at ", rarefy.override, "\n"))
        R <- vegan::rarefy(abund_table, rarefy.override ) #min(rowSums(abund_table)))
      }
  
    #warning(paste0("Richness: rarefying the data to minimum sample sum: ",min(rowSums(abund_table)),"\n\n"))
    #R <- vegan::rarefy(abund_table, 10000) #min(rowSums(abund_table)))
    #R <- abund_table
    df_R <- data.frame(sample = names(R), value = R, measure = rep("Richness", length(R)))
    if (is.null(df)) {
      df <- df_R
    }
    else {
      df <- rbind(df, df_R)
    }
  }
  if ("shannon" %in% method) {
    H <- vegan::diversity(abund_table)
    df_H <- data.frame(sample = names(H), value = H, measure = rep("Shannon",
                                                                   length(H)))
    if (is.null(df)) {
      df <- df_H
    }
    else {
      df <- rbind(df, df_H)
    }
  }
  if ("evenness" %in% method) {
    H <- vegan::diversity(abund_table)
    S <- specnumber(abund_table)
    J <- H/log(S)
    df_J <- data.frame(sample = names(J), value = J, measure = rep("Pielou's evenness",
                                                                   length(J)))
    if (is.null(df)) {
      df <- df_J
    }
    else {
      df <- rbind(df, df_J)
    }
  }
  if ("simpson" %in% method) {
    simp <- vegan::diversity(abund_table, "simpson")
    df_simp <- data.frame(sample = names(simp), value = simp,
                          measure = rep("Simpson", length(simp)))
    if (is.null(df)) {
      df <- df_simp
    }
    else {
      df <- rbind(df, df_simp)
    }
  }
  if ("invsimpson" %in% method) {
    simp <- vegan::diversity(abund_table, "invsimpson")
    df_simp <- data.frame(sample = names(simp), value = simp,
                          measure = rep("invSimpson", length(simp)))
    if (is.null(df)) {
      df <- df_simp
    }
    else {
      df <- rbind(df, df_simp)
    }
  }
  if ("fisher" %in% method) {
    alpha <- vegan::fisher.alpha(abund_table)
    df_alpha <- data.frame(sample = names(alpha), value = alpha,
                           measure = rep("Fisher alpha", length(alpha)))
    if (is.null(df)) {
      df <- df_alpha
    }
    else {
      df <- rbind(df, df_alpha)
    }
  }
  if ("pd" %in% method) {
    otu_tree <- phyloseq::phy_tree(physeq)
    PD <- pd(abund_table, otu_tree, include.root = TRUE)
    df_PD <- data.frame(sample = names(PD), value = PD,
                        measure = rep("PD", length(PD)))
    if (is.null(df)) {
      df <- df_PD
    }
    else {
      df <- rbind(df, df_PD)
    }
  }
  return(df)
}

