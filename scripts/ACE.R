#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ACE.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Obtain tumor purities using ACE, customized from ploidyplotloop function (ACE)
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1)
#
# History:
#  11-03-2022: File creation
#  23-03-2022: Write custom ACE function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
source("scripts/ACE_functions.R")
#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    Segments <- snakemake@input[["segments"]]
    sample <- snakemake@wildcards[["sample"]]
    binsize <-  snakemake@params[["binsize"]]
    ploidies <-  snakemake@params[["ploidies"]]
    method <- snakemake@params[["method"]]
    penalty <- snakemake@params[["penalty"]]
    fit_out <-  snakemake@output[["fit"]]
    errorgraph_out <- snakemake@output[["errorgraph"]]

}else{
    Segments <- "../data/copynumber/segments/TCGA-55-A491_segments.txt"
    sample <-  "TCGA-55-A491"
    binsize <- 1
    ploidies <-  c(2,4)
    method <- "RMSE"
    penalty <- 0
    fit_out <-  "../data/copynumber/ACE/TCGA-95-7039/2N/fits.txt"
    errorgraph_out <- "../data/copynumber/ACE/TCGA-95-7039/2N/errorgraph.svg"
}


#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# read segments
Segments <- read.delim(Segments, stringsAsFactors = F)
# Calculate segment width
Segments$Width <- Segments$End - Segments$Start
#-------------------------------------------------------------------------------
# 1.2 Reformat data
#-------------------------------------------------------------------------------
# obtain features and segment values
features <- paste0(Segments$Chromosome,":",Segments$Start,"-",Segments$End)
# convert log2 to copynumber 
segmentvalues <- 2^Segments$Segment_Mean

# create RLE
segmentdata <- rle(as.vector(na.exclude(setNames(segmentvalues, features))))
# expand segments by binsize (round to upper non-zero integer)
segmentdata$lengths <- ceiling(Segments$Width / binsize)
# Concatenate segment lengths
segmentdata <- suppressWarnings(Concatenate_lengths(segmentvalues,segmentdata))

#-------------------------------------------------------------------------------
# 2.2 Run ACE
#-------------------------------------------------------------------------------
Run_ACE(segmentdata, Segments, ploidies, sample, method, penalty, fit_out, errorgraph_out)

