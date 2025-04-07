#' Pairwise PERMANOVAs (adonis2) on phyloseq to find which pairs of a more complete set have significantly different community structures
#'
#' @param ps phyloseqobject
#' @param metadata variable (factor) all factors to be tested pairwise
#' @param sim.function  function to calculate the similarity matrix (vegdist or daisy). Default is 'vegdist'
#' @param sim.method similarity method from daisy or vegdist: default is 'bray'
#' @param p.adjust.m  p.value correction method, one of the methods supported by p.adjust(); default BH = 'Benjamini Hochberg'
#' @param hellinger boolean. Default TRUE to Hellinger transform the "distance" matrix
#' @param rngseed  The set.seed starting point. Default 2202
#' @param num.cores  Limit the number of cores to be used. Default = available cores - 1 
#'
#' @return a table with the pairwise factors, F-values, R^2, p.value and adjusted p.values
#' @author Alex Bossers (a.bossers\@uu.nl)
#' 
#' @note  Based initially on a function by Pedro Martinez Arbizu
#' @note  Scavenged from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
#' 
#' @export
#'
#' @import phyloseq
#' @import vegan
#' @import parallel
#'
#' @examples
#' data(iris)
#' pairwise.adonis(iris[,1:4],iris$Species)
#'
#' [1] Signif. codes:  Signif. codes:  0 = ***, 0.001 = **, 0.01 = *, 0.05 = .
#'                     pairs    F.Model        R2 p.value p.adjusted sig
#' 1    setosa vs versicolor  552.51266 0.8493496   0.001      0.003   *
#' 2     setosa vs virginica 1001.54509 0.9108722   0.001      0.003   *
#' 3 versicolor vs virginica   91.82959 0.4837475   0.001      0.003   *
#'
#' #similarity euclidean from vegdist and holm correction
#' #pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
#' adon=pairwise.adonis2( otu_table(physeq), sample_data(physeq)$Sample_type )
#'
pairwise.permanovas2 <- function( ps, variable, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='BH', hellinger=TRUE, rngseed=2202, num.cores="", ...)
{

  require(vegan)
  require(parallel)
  
# #debug
  # x=merged.treat3.r@otu_table@.Data
  # factors=  merged.treat3.r@sam_data$Location
  # sim.function = 'vegdist'
  # sim.method = 'bray'
  # p.adjust.m ='BH'
  # rngseed=2202
  
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
