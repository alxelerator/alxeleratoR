#' Using a given phyloseq object, take the taxonomy_table and fill any NAs with known higher level taxonomic data adding tax level k__ p__ f__ g__ etc
#'
#' @param ps phyloseqobject
#' @param verbose Default TRUE. Shows diff first 10 taxa before and after
#' @param genusspecies Default TRUE. Adds taxonomy column GenusSpecies (Genus+Species) with proper formatting (no prefix, or [Genus] etc on different levels avail only)
#' @return NA replaced phyloseq object
#'
#' @author a.bossers@uu.nl
#' @author alex.bossers@wur.nl
#' @note
#' Version log:
#'   version 1: 20200728
#'   version 2: 20211210 added support species and checks for presence of tax prefix and fixed non-NA issues.
#'   version 3: 20240206 added convenienc taxonomic column GenusSpecies. Fill with Genus+species is species avail (without g__ s__) or with "Genus [genus] etc all levels"
#' 
#' @import phyloseq
#' @import vegan
#' @export

taxa_resolve_NAs <- function( ps, verbose=TRUE, genusspecies=TRUE ) {
  
  require(phyloseq)
  require(dplyr)
  
  ## missing data at taxa level x (NA) can be filled with taxa data from taxa level x-1

  taxorg <- data.frame( tax_table( ps ) )
  tax2na <- taxorg

  #FOR NA fills we need to track taxranks in a prefix
  #FIRST check if prefix exists (like in UNITE) but not in SILVA
  if( ! any( grepl( '^k__', tax2na$Kingdom ) ) ) {
    
    message("Taxonomic rank prefixes not detected --> adding for the NA filled data\n")
    
    # no prefix => first add taxlevel prefix to allow tracking taxranks when filling NAs
    # make sure we do not replace NAs. In some data the NAs are actual missing NA. In some it is txt NA... somehow (diff silva unite?). Just test NA and exclude. Solved.
    tax2na$Kingdom[ ! is.na(tax2na$Kingdom) ] <- paste0( "k__", tax2na$Kingdom[ ! is.na(tax2na$Kingdom) ] )
    tax2na$Phylum[ ! is.na(tax2na$Phylum) ]  <- paste0("p__", tax2na$Phylum[ ! is.na(tax2na$Phylum) ] )
    tax2na$Class[ ! is.na(tax2na$Class) ]   <- paste0("c__", tax2na$Class[ ! is.na(tax2na$Class) ] )
    tax2na$Order[ ! is.na(tax2na$Order) ]   <- paste0("o__", tax2na$Order[ ! is.na(tax2na$Order) ] )
    tax2na$Family[ ! is.na(tax2na$Family) ]  <- paste0("f__", tax2na$Family[ ! is.na(tax2na$Family) ] )
    tax2na$Genus[ ! is.na(tax2na$Genus) ]   <- paste0("g__", tax2na$Genus[ ! is.na(tax2na$Genus) ] )
    if( 'Species' %in% colnames(tax2na) ) {
      #also treat Species
      tax2na$Species[ ! is.na(tax2na$Species) ]   <- paste0("s__", tax2na$Species[ ! is.na(tax2na$Species) ] )
    }
  }
  
  #now replace NA with higher level tax rank value
  cat("\nALL missing (NA) taxa at all taxa levels where filled with one level higher taxonomic data (tracking taxrank with prefix)!\n\n")
  
  #first we fill using prefixed NAs. Later repeat for actual <NA> that may be present/mixed
  tax2na$Phylum[ !is.na(tax2na$Phylum) & tax2na$Phylum == "p__" ] <- tax2na$Kingdom[ !is.na(tax2na$Phylum) & tax2na$Phylum == "p__"  ]
  tax2na$Class[ !is.na(tax2na$Class) & tax2na$Class == "c__" ]   <- tax2na$Phylum[ !is.na(tax2na$Class) & tax2na$Class == "c__" ]
  tax2na$Order[ !is.na(tax2na$Order) & tax2na$Order == "o__" ]   <- tax2na$Class[ !is.na(tax2na$Order) & tax2na$Order == "o__" ]
  tax2na$Family[ !is.na(tax2na$Family) & tax2na$Family == "f__" ] <- tax2na$Order[ !is.na(tax2na$Family) & tax2na$Family == "f__" ]
  tax2na$Genus[ !is.na(tax2na$Genus) & tax2na$Genus == "g__" ]   <- tax2na$Family[ !is.na(tax2na$Genus) & tax2na$Genus == "g__" ]
  if( 'Species' %in% colnames(tax2na) ) {
    #also treat Species
    tax2na$Species[ !is.na(tax2na$Species) & tax2na$Species == "s__" ] <- tax2na$Genus[ !is.na(tax2na$Species) & tax2na$Species == "s__" ]
  }
  
  # Repeat for real NAs... 
  # seems duplicate code of above but occasionally necessary. Most likely something with NA types..
  tax2na$Phylum[ is.na(tax2na$Phylum) ] <- tax2na$Kingdom[ is.na(tax2na$Phylum) ]
  tax2na$Class[ is.na(tax2na$Class) ]   <- tax2na$Phylum[ is.na(tax2na$Class) ]
  tax2na$Order[ is.na(tax2na$Order) ]   <- tax2na$Class[ is.na(tax2na$Order) ]
  tax2na$Family[ is.na(tax2na$Family) ] <- tax2na$Order[ is.na(tax2na$Family) ]
  tax2na$Genus[ is.na(tax2na$Genus) ]   <- tax2na$Family[ is.na(tax2na$Genus) ]
  if( 'Species' %in% colnames(tax2na) ) {
    #also treat Species
    tax2na$Species[ is.na(tax2na$Species) ] <- tax2na$Genus[ is.na(tax2na$Species) ]
  }

  # Handle the create of GenusSpecies taxonomic level
  if( genusspecies & 'Species' %in% colnames(tax2na) ) {
    # first handle if species level is present
    tax2na$GenusSpecies <- ifelse( substr( tax2na$Species, 1, 3) == "s__", 
      paste0( substr( tax2na$Genus, 4, 100), " ", substr( tax2na$Species, 4, 100) ) , 
      tax2na$Species )
    # handle cases with higher order taxonomy present but no species (reformatting)
    # in multi case the dplyr function seems more convenient
    tax2na <- tax2na %>%
      mutate(
        GenusSpecies = case_when(
          substr( GenusSpecies, 1, 3 ) == "g__" ~ paste0( substr( GenusSpecies, 4, 100 ), " [Genus]" ),
          substr( GenusSpecies, 1, 3 ) == "f__" ~ paste0( substr( GenusSpecies, 4, 100 ), " [Family]" ),
          substr( GenusSpecies, 1, 3 ) == "o__" ~ paste0( substr( GenusSpecies, 4, 100 ), " [Order]" ),
          substr( GenusSpecies, 1, 3 ) == "c__" ~ paste0( substr( GenusSpecies, 4, 100 ), " [Class]" ),
          substr( GenusSpecies, 1, 3 ) == "p__" ~ paste0( substr( GenusSpecies, 4, 100 ), " [Phylum]" ),
          substr( GenusSpecies, 1, 3 ) == "k__" ~ paste0( substr( GenusSpecies, 4, 100 ), " [Kingdom]" ),
          TRUE ~ GenusSpecies # Keeps the original value if none of the conditions above are met
        )
      )
  }

  if( verbose ) {
    taxshow <- taxorg
    taxshownew <- tax2na
    row.names(taxshownew) <- NULL
    row.names(taxshow) <- NULL
    cat("\nFROM taxonomy first 15 taxa:\n")
    print( head(taxshow,15) )
    cat("\nTO converted taxonomy first 15 taxa:\n")
    print( head(taxshownew,15) )
  }
  
   
  # fix structure
  tax2na[] <- lapply(tax2na, as.factor)
  tax_table(ps) <- tax_table( as.matrix( tax2na ) )

  return( ps )
}

