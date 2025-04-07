#' Pairwise PERMANOVAs
#'
#' NOT SURE IF THIS FUNCTION STILL MAKES SENSE. KEPT just in case from SJurburg
#' @note DEPRECATED function
#'
#' @param x community table, site by species
#' @param factors  a column or vector with all factors to be tested pairwise
#' @param sim.function  function to calculate the similarity matrix. IE. 'vegdist' default is 'vegdist'
#' @param sim.method similarity method from daisy or vegdist: default is 'bray'
#' @param p.adjust.m  p.value correction method, one of the methods supported by p.adjust(); default is 'bonferroni'
#'
#' @return a table with the pairwise factors, F-values, R^2, p.value and adjusted p.value
#' @author Pedro Martinez Arbizu
#' @note  scavenged from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
#' @export
#'
#' @import phyloseq
#' @import vegan
#'
#' @examples
#' data(iris)
#' pairwise.adonis(iris[,1:4],iris$Species)
#'
#' [1] Signif. codes:  Signif. codes:  0 = ***, 0.001 = **, 0.01 = *, 0.05 = .
#'                    pairs    F.Model        R2 p.value p.adjusted sig
#' 1    setosa vs versicolor  552.51266 0.8493496   0.001      0.003   *
#' 2     setosa vs virginica 1001.54509 0.9108722   0.001      0.003   *
#' 3 versicolor vs virginica   91.82959 0.4837475   0.001      0.003   *
#'
#' #similarity euclidean from vegdist and holm correction
#' #pairwise.adonis(x=iris[,1:4],factors=iris$Species,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
#' adon=pairwise.adonis(physeq_pruned@otu_table@.Data,physeq_pruned@sam_data$Sample_type)
#'
pairwise.adonis <- function( x, factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{

stop("Deprecated function. Use pairwise.permanovas2() instead!")

co = combn(unique(as.character(factors)),2)
pairs = c()
F.Model =c()
R2 = c()
p.value = c()


for(elem in 1:ncol(co)){
if(sim.function == 'daisy'){
library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
} else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}

ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
F.Model =c(F.Model,ad$aov.tab[1,4]);
R2 = c(R2,ad$aov.tab[1,5]);
p.value = c(p.value,ad$aov.tab[1,6])
}
p.adjusted = p.adjust(p.value,method=p.adjust.m)
sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'

pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
print("Signif. codes:  0 = ***, 0.001 = **, 0.01 = *, 0.05 = .")
return(pairw.res)

}


