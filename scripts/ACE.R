#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ACE.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Obtain tumor purities using ACE, customized from ploidyplotloop function (ACE)
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) Add normal segments -> DONE
#
# History:
#  11-03-2022: File creation
#  23-03-2022: Write custom ACE function
#  31-03-2022: Add squaremodel implementation
#  08-04-2022: Add chromlengths and Add_normal_segments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    Segments <- snakemake@input[["segments"]]
    sample <- snakemake@wildcards[["sample"]]
    ploidy <-  snakemake@wildcards[["ploidy"]]
    binsize <-  snakemake@params[["binsize"]]
    method <- snakemake@params[["method"]]
    penalty <- snakemake@params[["penalty"]]
    penploidy <- snakemake@params[["penploidy"]]
    chromlengths <- snakemake@params[["chromlengths"]]
    fit_out <-  snakemake@output[["fit"]]
    errorgraph_out <- snakemake@output[["errorgraph"]]
}else{
    Segments <- "../data/copynumber/segments/TCGA-53-7624_segments.txt"
    sample <-  "TCGA-53-7624"
    binsize <- 10
    ploidy <-  "squaremodel" 
    method <- "RMSE"
    penalty <- 0.5
    penploidy = 0.5
    chromlengths <- "reference/Chrom_lengths.txt"
    fit_out <-  "../data/copynumber/ACE/TCGA-53-7624/2N/fits.txt"
    errorgraph_out <- "../data/copynumber/ACE/TCGA-53-7624/2N/errorgraph.svg"
}
#-------------------------------------------------------------------------------
# 1.2  Import Libraries
#-------------------------------------------------------------------------------
source("scripts/ACE_functions.R")

#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# read segments
Segments <- read.delim(Segments, stringsAsFactors = F)
# read chromlengths
chromlengths <- read.delim(chromlengths)

#-------------------------------------------------------------------------------
# 1.2 Reformat data
#-------------------------------------------------------------------------------
# exclude sex chromosomes
Segments <- Segments[!Segments$Chromosome %in% c("X","Y"),]
# Add normal segments
Segments <- Add_normal_segments(Segments,chromlengths)
# Calculate segment width
Segments$Width <- Segments$End - Segments$Start
# obtain features and segment values
features <- paste0(Segments$Chromosome,":",Segments$Start,"-",Segments$End)
# convert log2 to copynumber 
segmentvalues <- 2^Segments$Segment_Mean
Segments$Segment_Mean
segmentdata_test
# create RLE
segmentdata <- rle(as.vector(na.exclude(setNames(segmentvalues, features))))
# expand segments by binsize (round to upper non-zero integer)
segmentdata$lengths <- ceiling(Segments$Width / binsize)
# Concatenate segment lengths
segmentdata <- suppressWarnings(Concatenate_lengths(segmentvalues,segmentdata))
#-------------------------------------------------------------------------------
# 2.2 Run ACE
#-------------------------------------------------------------------------------
if(ploidy == "squaremodel"){
    Run_squaremodel(segmentdata, Segments, sample, method, penalty, penploidy, fit_out, errorgraph_out)

}else{
    Run_ACE(segmentdata, Segments, ploidy, sample, method, penalty,penploidy, fit_out, errorgraph_out)
}
