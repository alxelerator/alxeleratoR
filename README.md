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

## License
This project is licensed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html). You are free to use, modify, and distribute this package as long as you adhere to the terms of the license.  

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

ToDo overview functions
ToDo overview to credit initially scavenged functions (used as a basis)
