#' Pairwise PERMANOVAs (adonis2) on phyloseq to find which pairs of a more complete set have significantly different community structures
#'
#' @param ps phyloseqobject
#' @param variable all factors to be tested pairwise "single" or multiple c("A","B")
#' @param sim.function  function to calculate the similarity matrix (vegdist or daisy). Default is 'vegdist'
#' @param sim.method similarity method from daisy or vegdist: default is 'bray'
#' @param p.adjust.m  p.value correction method, one of the methods supported by p.adjust(); default BH = 'Benjamini Hochberg'
#' @param hellinger boolean. Default TRUE to Hellinger transform the "distance" matrix
#' @param rngseed  The set.seed starting point. Default 2202
#' @param num.cores  Limit the number of cores to be used. Default = available cores - 1 
#' @param ... Transfer any other common variables to the adonis2 function
#'
#' @return a table with the pairwise factors, F-values, R^2, p.value and adjusted p.values
#' @author Alex Bossers (a.bossers\@uu.nl)
#' 
#' @note  Based initially on a function by Pedro Martinez Arbizu
#' @note  Initially scavenged from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
#' 
#' @export
#'
#' @import phyloseq
#' @import vegan
#' @import parallel
#'
#' @examples
#' #library(phyloseq)
#' #data(GlobalPatterns)
#' #pairwise.permanovas2( GlobalPatterns, variable=c("SampleType"), 
#' #                      sim.method = "bray", hellinger = TRUE, p.adjust.m = "BH")
#' #
#' ## No significant comparisons in this dataset.
#' #
#' #> Signif. codes: 0 = ***, 0.001 = **, 0.01 = *, 0.05 = .
#' #>                                       Pairs   F.Model        R2     Pr(>F)     p.adj sig
#' #> 1                             Soil vs Feces  7.937235 0.6135186 0.04400000 0.1028571    
#' #> 2                              Soil vs Skin  5.420047 0.5753737 0.10000000 0.1028571    
#' #> 3                            Soil vs Tongue  5.832750 0.6603549 0.10000000 0.1028571    
#' #> 4                        Soil vs Freshwater  4.832812 0.6169958 0.10000000 0.1028571    
#' #> 5                Soil vs Freshwater (creek)  8.934814 0.6907571 0.10000000 0.1028571    
#' #> 6                             Soil vs Ocean  8.016178 0.6671155 0.10000000 0.1028571    
#' #> 7                Soil vs Sediment (estuary)  8.045936 0.6679378 0.10000000 0.1028571    
#' #> 8                              Soil vs Mock 13.453708 0.7708223 0.10000000 0.1028571    
#' #> 9                             Feces vs Skin  7.440478 0.5980862 0.02800000 0.1028571    
#' #>...
#'
pairwise.permanovas2 <- function( ps, variable, sim.function = 'vegdist', sim.method = 'bray',
                                  p.adjust.m ='BH', hellinger=TRUE, 
                                  rngseed=2202, num.cores="", ...)
{
  require(vegan)
  require(parallel)
  
factors <- sample_data(ps)[[variable]]

if( taxa_are_rows(ps) ) {
  #flip
  x <- t(otu_table(ps))
} else {
  x <- otu_table(ps)
}

co = combn(unique(as.character(factors)),2)
pairs = c()
F.Model =c()
R2 = c()
p.value = c()

if( num.cores == "" ){
  num.cores <- parallel::detectCores() - 1
}

for( elem in 1:ncol(co) ) {
  otutab = x[factors %in% c(co[1,elem],co[2,elem]),]
  if(sim.function == 'daisy'){
    message( paste0("Using DAISY function and method ",sim.method) )
    library(cluster); distn = daisy(otutab,metric=sim.method)
  } else {
    if(hellinger) {
      message( paste0("Using Hellinger transformed data and vegdist() distance metric '",sim.method,"'") )
      distn <- vegdist(decostand(otutab, "hellinger"), method=sim.method)
    } else {
      message( paste0("Using plain (not Hellinger transformed) data and vegdist() distance metric '",sim.method,"'") )
      distn <- vegdist( otutab, method=sim.method )
    }
  }

  pair <- paste(co[1,elem],'vs',co[2,elem])
  message(paste("Analysing pair:", pair))
  set.seed(rngseed)
  ad2 <- adonis2( distn ~ factors[factors %in% c(co[1,elem],co[2,elem])], parallel=num.cores,  ... )
  pairs <- c( pairs, pair);
  F.Model <- c( F.Model, ad2$F[1]);
  R2 <- c( R2, ad2$R2[1] );
  p.value <- c( p.value, ad2$`Pr(>F)`[1] )
}


p.adjusted <- p.adjust( p.value, method=p.adjust.m )
sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'

pairw.res = data.frame( "Pairs" = pairs, F.Model, R2, `Pr(>F)` = p.value, "p.adj" = p.adjusted, sig, check.names=FALSE)
message("Signif. codes:  0 = ***, 0.001 = **, 0.01 = *, 0.05 = .")
message(paste("Multiple testing correction method:", p.adjust.m ) )
message(paste0("Used random seed '", rngseed, "' on ", num.cores, " cores." ) )

return(pairw.res)

}
