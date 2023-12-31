#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compare_purity_measures.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Compare different purity measures
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) 
#
# History:
#  24-03-2022: File creation 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressWarnings(library(ggplot2))
suppressMessages(suppressWarnings(library(dplyr)))
source('scripts/ACE_functions.R')
#-------------------------------------------------------------------------------
# 1.1 Parse snakemake objects
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    ACE_fits <- snakemake@input[["ACE_fits"]]
    Purities <- snakemake@input[["Purities"]]
    Purities_out <-  snakemake@output[["Tumor_purities"]]
    Scatterplots_out <-  snakemake@output[["Scatterplots"]]
    Samples_to_exclude_out <-  snakemake@output[["Samples_to_exclude"]]
}else{
    ACE_fits <- Sys.glob("../data/copynumber/ACE/*/**/fits.txt")
    Purities <-  "../data/purity/Tumor_purities.tsv"
    Purities_out <- "../data/purity/ACE_Tumor_purities_TCGA-LUAD.tsv"
    Scatterplots_out <- "../data/copynumber/ACE/TCGA-95-7039/2N/errorgraph.svg"
}
#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# read ACE purities
ACE_purities <-
    tibble::tribble(
                ~sample,~ACE,
                purrr::map_chr(ACE_fits,~strsplit(.x,"/")[[1]][5]), purrr::map_dbl(ACE_fits, ~read.table(.x,header=T)[1,2])) %>%
    tidyr::unnest(cols = c(sample, ACE))



Samples_to_exclude <- ACE_purities[ACE_purities$ACE == 1,"sample"]
print(nrow(ACE_purities))

ACE_purities <-
    ACE_purities %>%
    dplyr::filter(!sample %in% Samples_to_exclude) %>%
    unique()

print(nrow(ACE_purities))

# read TCGA purities
Purities <- read.delim(Purities, stringsAsFactors = F)
#-------------------------------------------------------------------------------
# 1.2 Join data
#-------------------------------------------------------------------------------
# Obtain short sample name
Purities$sample <- purrr::map_chr(Purities$Sample.ID,~paste0(strsplit(.x,"-")[[1]][1:3], collapse = "-"))
# filter out samples 
Purities <- Purities[Purities$sample %in% ACE_purities$sample,]
# Join ACE purities
Purities <-
    Purities %>%
    left_join(ACE_purities)

#-------------------------------------------------------------------------------
# 2.1 Plot data
#-------------------------------------------------------------------------------
# create plot list
plotlist <- c()
Purity_measures <- c("ESTIMATE","ABSOLUTE","LUMP","IHC","ACE")
#Convert purities to numeric
Purities[Purity_measures] <- sapply(Purities[Purity_measures],function(x)as.numeric(gsub(",","\\.",x)))
#Exclude ACE fit=1
#Purities <- Purities[Purities$ACE != 1,]
# Obtain measure combinations
Measure_combinations <- t(combn(Purity_measures,2))

for(i in 1:nrow(Measure_combinations)){
    xvar = Measure_combinations[i,1]
    yvar = Measure_combinations[i,2]
    p <- ggpubr::ggscatter(Purities, x = xvar, y = yvar,
                           add = "reg.line",
                           add.params = list(color = "red"),
                           conf.int = TRUE,
                           palette = "jco",
                           alpha = 0.5,
                           size = 1,
                           xlim = c(0,1),
                           ylim = c(0,1))+
        ggpubr::stat_cor(label.x = 0, label.y = 0)
    plotlist[[i]] <- p 
}


svg(Scatterplots_out, width = 12, height = 9)
cowplot::plot_grid(plotlist = plotlist, align = "hv")
dev.off()

# write to file
write.table(ACE_purities, file = Purities_out, quote = F, row.names = F,sep = "\t")

# write to file
write.table(Samples_to_exclude, file = Samples_to_exclude_out, quote = F, row.names = F,sep = "\t")
