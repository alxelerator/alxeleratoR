#' Calculates the Kruskal-Wallis rank-sum test on SIMPER pretty output
#' 
#' Debugged, changed and extended to allow an phyloseq object and tax rank selection
#' 
#' What it does:
#'   This script takes the output .csv of simper_pretty.R, and the OTU/metadata/taxonomy tables, and performs the non-parametric Kruskal-Wallis rank-sum test
#'   on each OTU in the .csv file. Output is a .csv file containing the same contents of simper.pretty output with the following info: p-value, fdr corrected
#'   p-value, OTU taxonomic classification (if applicable), mean rel. abund and std dev of otu/tax_lvl in group 1 of comparison, and mean rel. abund and std
#'   dev of otu/tax_lvl in group 2 of comparison.
#' 
#' @param ps phyloseq object
#' @param csv CSV output file from simper.pretty2
#' @param interesting columns of var of interest in metadata (same as simper.pretty input)
#' @param output_name desired name of output (i.e. outputname_krusk_simper.csv)
#' @param taxrank For example Kingdom or Genus. If nothing/empty string given, a selected Taxonomic concat will be made
#' @return CSV df (and file saved)
#' @import vegan
#' @import dplyr
#' @export
#'
#' @note csv = clean_simper.csv (*Must be imported as data.frame; i.e. csv= data.frame(read.csv("PATH to .csv")))
#' @note kruskal.pretty: automates kruskal test execution on simper_pretty.R output
#' @author Andrew Steinberger asteinberger@wisc.edu   Suen Lab, University of Wisconsin-Madison
#' @author Alex Bossers \email{a.bossers@uu.nl} for phyloseq and extensions
#' 
#'       Copyright (C) 2016 Andrew Steinberger
#' 
#'     This program is free software: you can redistribute it and/or modify
#'     it under the terms of the GNU General Public License as published by
#'     the Free Software Foundation, either version 3 of the License, or
#'     (at your option) any later version.
#' 
#'     This program is distributed in the hope that it will be useful,
#'     but WITHOUT ANY WARRANTY; without even the implied warranty of
#'     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#'     GNU General Public License for more details.
#' 
#'     You should have received a copy of the GNU General Public License
#'     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#' 
#' Changelog
#' 
#' 
#' ###########################################################################################
#' Example of use
#' kruskal.pretty(bacteria, metrics, csv, c('sloth_sp','type','sp.type'), 'sloth', taxonomy)

kruskal.pretty2 = function( ps, csv, interesting, output_name, taxrank='' ){
  library(vegan)
  library(dplyr)

  # ps=merged.bac.nr.contact.5k.r.clean
  # csv=contact.simp1
  # interesting=c("EffectContact")
  # output_name=paste0( tabfolder, '/SIMPER_Contact_Kruskal1.csv')

  otu=ps@otu_table@.Data
  metrics=ps@sam_data
  taxonomy=ps@tax_table@.Data
  csv$Comparison = as.factor(csv$Comparison)

  #if(grepl("Otu", colnames(otu)[1])!=TRUE){
  #  #converts output from A.Neuman Taxonomy script
  #  otu=as.data.frame(t(otu))
  #}
  #changing csv$X to rownames to allow proper splitting of comparisons
  csv$X=as.integer(rownames(csv))
  L=list()
  R=list()
  mean_L=c()
  sd_L=c()
  mean_R=c()
  sd_R=c()
  L_mean=c()
  R_mean=c()
  L_sd=c()
  R_sd=c()
  krusk=c()
  tax=c()
  L_abund=c()
  R_abund=c()
  L_abund_sd=c()
  R_abund_sd=c()
  abund=as.matrix(otu)
  abund=abund/rowSums(abund)
  for(b in levels(csv$Comparison)){
    otu_list=dplyr::filter(csv, Comparison==b) #saves otu list for current comparison
    for(i in csv$X){
      if(as.character(csv$Comparison[i])==b){  ##splitting comparisons so can call individually for table generation
        splt=as.data.frame(matrix(unlist(strsplit(as.character(csv$Comparison[i]),'_')), nrow=1, byrow=T))
        cola=as.character(splt[1,1])
        colb=as.character(splt[1,2])
        break
      }
    }

    #saving topic containing var of interest (cola/colb) (less memory intensive)
    for(topic in interesting){
      #preventing crash if there is only one topic in interesting
      #metrics$topic = as.factor(metrics$topic)
      if(is.null(levels(metrics[[topic]]))==TRUE){
        topic1=topic
        break
      }
      for(sbtpic in levels(metrics[[topic]])){
        if(sbtpic==cola){
          topic1=topic
          break
        }
      }
    }

    #iterate thru rows in tpics of intrst til matches cola and colb, generates otu and metrics tbl  ##!Processing can be reduced!##
    for(rowe1 in metrics[[topic1]]){
      for(rowe2 in metrics[[topic1]]){
        if(rowe1==cola & rowe2==colb){
          listbact=otu[c(metrics[[topic1]]==cola|metrics[[topic1]]==colb),]
          listmet=metrics[c(metrics[[topic1]]==cola|metrics[[topic1]]==colb),]
          break
        }
      }
    }
    #collecting differential abundances
    sample_L=row.names(subset(metrics, metrics[[topic1]] == c(cola)))
    sample_R=row.names(subset(metrics, metrics[[topic1]] == c(colb)))
    #collecting abund values, perform/save mean and stdev calculations
    for(otus in otu_list$OTU){
      for(sample in sample_L){
        L=append(L,abund[sample,otus])
        mean_L[[otus]]=mean(as.numeric(L))
        sd_L[[otus]]=sd(as.numeric(L))
      }
      for(sample in sample_R){
        R=append(R,abund[sample,otus])
        mean_R[[otus]]=mean(as.numeric(R))
        sd_R[[otus]]=sd(as.numeric(R))
      }
      L=list()
      R=list()
    }

    #runs kruskal.test for each otu in simper csv, stores as list, also stores abundances
    for(otus in otu_list$OTU){
      result=kruskal.test( listbact[,otus]~listmet[[topic1]])
      krusk=append(krusk, result$p.value)
      #stores taxonomic classification for each otu as list
      if(missing(taxonomy)){
        tax=append(tax, c("NA"))
      } else {
        taxline=taxonomy[otus,]
        taxrank=""
        if(taxrank > '') {
          if(hasName(taxline,taxrank)) {
              taxline["Taxonomy"]=taxline[taxrank]
          } else {
              stop(paste0("Bad taxonomic rank given '",taxrank,"' for the given taxonomy table."))
          }
        } else {
          # no taxrank given so compile a combination
          #first test is Species is present or not
          if(hasName(taxline,"Species")) {
            # test if taxonomy species was prefixed with s__ or not yet
            if( is.na(taxline["Species"]) ) {
              taxline["Species"] = paste0("s__",taxline["Species"])
            } else {
              if( substr( taxline["Species"],1,3) != "s__") {
                taxline["Species"] = paste0("s__",taxline["Species"])
              }
            }
            taxline["Taxonomy"]=paste(taxline["Kingdom"],taxline["Phylum"],taxline["Genus"],taxline["Species"], sep="|")
          } else {
            #until genus
            taxline["Taxonomy"]=paste(taxline["Kingdom"],taxline["Phylum"],taxline["Genus"], sep="|")
          }
        }
        tax=append(tax, as.character(taxline["Taxonomy"]))
      }
      L_mean=append(L_mean, as.character(mean_L[[otus]]))
      R_mean=append(R_mean, as.character(mean_R[[otus]]))
      L_sd=append(L_sd, as.character(sd_L[[otus]]))
      R_sd=append(R_sd, as.character(sd_R[[otus]]))
    }
  }
  #adjusted p-values for multiple comparisons
  fdr=p.adjust(krusk, method='fdr')
  #order csv to match 'krusk'/'fdr' list, add p.val, add taxonomy, re-ord to match orig csv, write to csv
  o_csv=dplyr::arrange(csv, Comparison)
  o_csv[,5]=krusk
  o_csv[,6]=fdr
  o_csv[,7]=tax
  o_csv[,8]=L_mean
  o_csv[,9]=L_sd
  o_csv[,10]=R_mean
  o_csv[,11]=R_sd
  o_csv=dplyr::arrange(o_csv, X)
  colnames(o_csv)[which(names(o_csv) == "V5")] <- "krusk_p.val" #changes column header
  colnames(o_csv)[which(names(o_csv) == "V6")] <- "fdr_krusk_p.val"
  colnames(o_csv)[which(names(o_csv) == "V7")] <- "Taxonomy"
  colnames(o_csv)[which(names(o_csv) == "V8")] <- "Left mean abund"
  colnames(o_csv)[which(names(o_csv) == "V9")] <- "Left stdev"
  colnames(o_csv)[which(names(o_csv) == "V10")] <- "Right mean abund"
  colnames(o_csv)[which(names(o_csv) == "V11")] <- "Right stdev"
  #o_csv[,1]=NULL
  write.tsv(o_csv, file=output_name)
  return( o_csv )
}

