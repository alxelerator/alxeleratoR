# alxeleratoR

Public R package to **accelerate** data analysis, so you can focus on interpretation during data exploration instead of finding out how to code your question.  

This "package"/bundle evolved from my need to automate and share common analyses and visualizations. Most functions focus on using phyloseq or similar objects, and most are wrappers around, or drop-in replacements of, existing functions or packages, some are original creations, and others are "borrowed/scavanged" and adapted to my needs (with credit given as accurately as possible. *If you think I’ve forgotten to credit you, please drop me a message*).  

The primary context is microbiome, resistome, metagenomics, and eDNA research. But feel free to use them in different contexts.

This package/repository is maintained by:  
* **Alex Bossers** (creator)  
  + One Health Microbial Group, Institute for Risk Assessment Sciences, Utrecht University, NL ([a.bossers@uu.nl](mailto:a.bossers@uu.nl))  
  + PathogenOmics and Bioinformatics Group, Wageningen BioVeterinary Research, NL ([alex.bossers@wur.nl](mailto:alex.bossers@wur.nl))  
* **You?**  

Any additions or bug fixes are more than welcome!  

## Disclaimer
This package are primarily tailored to my own needs and may be subject to changes. While I will strive to ensure backward compatibility, this is not guaranteed (including preservation of defaults).  

The use of these functions is entirely at your own risk. You should always validate the outcomes, as they may contain errors or inaccuracies. Neither I, any collaborators, nor the associated institutions can be held liable for any direct or indirect damage caused by the use of this package. Proceed with caution and ensure thorough validation when using or integrating this package into your workflow!  

## Installation
In R run `devtools::install_github("alxelerator/alxeleratoR")`  
After loading the package `library(alxeleratoR)` you can set some convenient defaults using `alxinit()`

If you need to install some missing packages I can strongly recommend the use of BioConductor package management installed using `install.packages("BiocManager")`  
For instance `BiocManager::install(c("tidyverse","pals","cowplot","RColorBrewer","foreach","openssl","Biostrings","rstatix","scales","vegan"))`  

## Documentation
Function details are covered by each function's documentation.  

Alphabetic overview of provided functions

Function | What it should do 
--- | --- 
alpha_div() | Calculate alpha diversity metrics on a phyloseq object.
cal_log1p_score() | Calculate recalibration matrix for heatmaps Apply transformation function of log avoiding 0 (1p).
cal_z_score() | Calculate recalibration matrix for heatmaps Apply transformation function Z score is difference to mean divided by sd.
clr_transform() | Apply the centered log-ratio trasformation on a relative abundance matrix! Species in columns and samples in rows. Defined as log^base( x / geometric_mean(rowx) ).
colorDistinct() | Generates a distinct scale of n=x elements based on QUAL Rcolorebrewer set. Manually replaced a few too dominant or too dark colors with other variants.
countTaxrankNAs() | Counts the number of NAs in the taxonomy table per taxonomic rank.
ggrare2() | Rarefaction curve from phyloseq, ggplot styled with custom step-size and incl standard error.
hdisp() | Multiple tests for homogeneity of dispersions (vegan::betadisper)
kruskal.pretty2() | Calculates the Kruskal-Wallis rank-sum test on SIMPER pretty output
pairwise.adonis() | [deprecated] Calculate pairwise PERMANOVAs on phyloseq.
pairwise.permanovas2() | Pairwise PERMANOVAs (adonis2) on phyloseq to find which pairs of a more complete set have significantly different community structures
perform_anova() | Function override in of microbiomeseq to calculate pair-wise anovas allowing more options and Welch-ANOVA
permanova2() | Single (multivariable) PERMANOVA v2 Uses vegan::adonis2. This function is just a convenience wrapper for vegan::adonis2. 
permanovas() | [deprecated] Multiple Permanovas on phyloseq (adonis1)
permanovas2() | Multiple Permanovas v2 Uses vegan::adonis2 which returns a slight different object. V2 now returns a proper data frame of tested vars and R2 and F.
phylodiv() | Phylodiv calculates phylogenetic diversity from phyloseq (req tree)
plot_anova_diversity2() | Convenience wrapper to ggplot alpha diversity combined with pair-wise ANOVA using a (rarefied) phyloseq object as a starter. Classic or Welch ANOVA supported.
plot_procrustes() | Plot calculated procrustus overlay
prevalence_plot() | Plot prevalence plots of certain taxonomic rank from phyloseq object
prevalence_taxa() | Calculate the prevalence at lowest taxonomic rank (i.e. ASV) from phyloseq. Function is used in the prevalence_plot() function
procrustes_error() | Plot calculated procrustes residuals
psmelt2() | Phyloseq psmelt() dplyr optimised for large datasets
ps_co_occurrence() | Calculate the co-occurrence of taxa (pair-wise) between one or two tax tables in phyloseq objects
ps_plot_bar() | Convenience wrapper to generate tuned and stack-sorted barplots at various taxonomic levels from phyloseq input
rarefy_custom_depth() | Convenience function to rarefy a given phyloseq object not to one value for each sample but to allow to rarefy each sample with its own given rarefy target depth.
sample_sum_plot() | Make a (coloured) ggplot of sample sums per sample ordered decreasing by sample sum and some stats like 90% rarefy number on phyloseq
SaveVariousPlots() | Saves ggplot plots with reasonable defaults in png and svg (png has explicit white background (transparency off))
short.tn() | Shorten too long temp OTU names from provided vector.
shortASVs() | Move ASV sequence as sequence_name into the phyloseq refseq slot. Also shorten the taxa names to Hashed-sequences or ASV0001 etc
simper.pretty2() | Perform SIMPER and return pretty output using phyloseq (debugged).
taxa_resolve_NAs() | Replace any NAs in taxonomy_table (phyloseq) with known higher level taxonomic data adding tax level k__ p__ f__ g__ etc
tax_glom2() | Agglomerate taxa using dplyr optimed tax_glom function from phyloseq.
toptaxa() | Returns the top-x list of taxa names at given rank for given samples based on its relative abundances.
write.tsv() | Write data.frame to file using tabular tab separators.

## Credits
ToDo overview to credit initially scavenged functions (used as a basis for some provided functions here).

## License
This project is licensed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this package as long as you adhere to the terms of the license.  
