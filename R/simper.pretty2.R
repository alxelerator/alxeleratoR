#' Perform SIMPER and return pretty output using phyloseq (debugged).
#'
#' @author a.bossers@uu.nl
#' 
#' What it does:
#'   This script is meant to rapidly perform the SIMPER function from the R package vegan for all comparisons of interest in a data set. Inputs are OTU
#'   and metadata tables, and the output is a .csv. User can tailor contents of .csv by setting perc_cutoff, low_cutoff, and low_val. This function can
#'   also handle taxonomic levels instead of OTU, but currently only select formats are compatible. Requires installation of the R package ‘vegan’.
#' 
#' Inputs:
#'   perc_cutoff: % cutoff for output OTUs, as decimal (i.e. write 50% as 0.5), larger % increases number OTUs in output.
#'   low_cutoff : ‘y’ if want to REMOVE OTUs that contribute less than 1%
#'   low_val    : set value of low cutoff (0.01), ignored if low_cutoff=‘n’.
#' 
#' Old doc org function:
#' simper.pretty: automates simper exectution for comparisons of interest
#' Andrew Steinberger, Kim Dill-Mcfarland, Madison Cox
#' asteinberger@wisc.edu
#' Suen Lab
#' University of Wisconsin-Madison
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
#' Simper (OTU correlations within inver-metric correlations)
#' 
#' Makes list with the names of all the columns from metrics file that you want to
#' analyze with SIMPER and combines all results into csv file
#' 
#' @param ps phyloseq object having: otu table OTU-names in columns! Function will now transpose.
#' @param interesting list the column headers of interest from metrics file ex. c('sloth_sp','type','sp.type')
#' @param perc_cutoff % cutoff desired for SIMPER output, as decimal (i.e. write 50% as 0.5)
#' @param low_cutoff 'y' if you want to REMOVE OTUs that contribute to less than 1% of significance
#' @param low_val value of low cutoff (0.01)
#' @param output_name name to append to the simper and clean_simper output files
#' @param cores number of cores to use. Default limits to cores available minus 2. Windows does not do multi core by default.
#' 
#' Alex: cores can be specified to limit used cores otherwise all available are used (minus 2)!
#' 
#' @export
#' 
#' @import phyloseq
#' @import vegan
#' @import dplyr
#' @import reshape2



simper.pretty2 <- function( ps, interesting, perc_cutoff, low_cutoff, low_val, output_name, cores=detectCores()-2 ) {

  require(vegan)

  # we expect a PS object so we do not need to transpose the OTU table since ASVs are column headers and samples are rows
  if (ps@otu_table@taxa_are_rows == TRUE ) {
    x = t( data.frame( otu_table(ps) ) )
    cat("OTUs are rows -> transposing data\n")
  } else {
    x = data.frame( otu_table(ps) )
  }
  metrics = data.frame( sample_data(ps) )
  #handling otu names if nonetables for taxa levels
  #save = FALSE
  # if(! grepl("Otu", colnames(x)[1]) ) {
  #   # If we have ASVs for instance save the seq ASV identifier = name = ~OTU
  #   save=TRUE
  #   #x=as.data.frame(t(x))
  #   orig_names=colnames(x)
  #   new_names=list()
  #   l=1
  #   for(n in colnames(x)){
  #     ifelse((l<10), (colnames(x)[l]=c(paste0('Otu000',c(l)))), (colnames(x)[l]=c(paste0('Otu00',c(l)))))
  #     new_names=append(new_names, colnames(x)[l])
  #     l=l+1
  #   }
  #   orig_names=as.data.frame(orig_names, row.names = as.character(new_names))
  # }

  # running simper
  for(variables in interesting){
    #x must be sample=rows and OTU=cols
    test_1=with(metrics, simper(x, metrics[[variables]], parallel = cores ))
    #parsing through simper output, saving desired info to table
    for(name in names(test_1)){
      #number of params to test so many cols to define
      testmx=matrix(ncol=length(interesting))
      testmx=cbind(test_1[[name]]$ord,test_1[[name]]$cusum)
      sorted=testmx[order(testmx[,1]),]
      sorted=cbind(sorted,test_1[[name]]$species)
      sorted=sorted[order(sorted[,2]),]
      t=matrix(sorted[sorted[,2]<=perc_cutoff,],ncol=3)
      i=nrow(t)
      #converting percents to percent of whole
      while(i>1){
        t[i,2]=as.character(as.numeric(t[i,2])-as.numeric(t[i-1,2]))
        i=i-1
      }
      t[,1]=name
      write.table(t,file=paste(output_name,'_simperTemp.csv',sep=""), append=TRUE, sep=",", col.names = FALSE)
    }}
  y=read.table(paste(output_name,'_simperTemp.csv',sep=""), header=FALSE,sep=",",fill = TRUE,row.names = NULL)

  y=y[-c(1)]
  cat("\nWriting out the SIMPER pair-wise data.\nNote that SIMPER_CumulContrib is the Ordered cumulative contribution. These are based on item average, but they sum up to total 1.\n\n")
  colnames(y) = c("Comparison", "SIMPER", "OTU")
  #removing results lower than low cutoff
  if(low_cutoff=='y'){
    y=y[!(as.numeric(as.character(y$SIMPER))<low_val),]
  }

  # #prevents changing of colnames if OTU table
  # if(save){
  #   y$OTU=orig_names[match(y$OTU, rownames(orig_names)),1]
  # }
  write.table( y, file=output_name, sep="\t", quote=FALSE)
  file.remove(paste0(output_name,'_simperTemp.csv'))
  return( data.frame(y) )
}
