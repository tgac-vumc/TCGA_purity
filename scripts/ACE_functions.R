#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ACE_functions.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Obtain functions used in ACE.R
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# TODO:
# 1) Add function to plot summary files (ACE adjusted segments for each fit)
#
# History:
#  23-03-2022: File creation, write code
#  30-03-2022: Create function to fix segmentdata, Replace Run_ACE code
#              with square_model code
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressWarnings(library(ggplot2))
#-------------------------------------------------------------------------------
# Function to fix segments
Concatenate_lengths <- function(segmentvalues,segmentdata){
    # if segmentdata is valid do nothing
    if(length(segmentdata$lengths) == length(segmentdata$values)){
        return(segmentdata)
    }else{
        # find index where segmentvalues are concatenated
        ix <- min(which(!segmentdata$values == segmentvalues))
        # concatenate lengths
        segmentdata$lengths[ix] <- segmentdata$lengths[ix-1] + segmentdata$lengths[ix]
        # Remove redundant length
        segmentdata$lengths <- segmentdata$lengths[-(ix-1)]
        return(segmentdata)
    }
}

# Plot errorlist
Plot_errors <- function(errordf, minimadf,sample, outfile) {
    p <- ggplot2::ggplot() +
        geom_raster(data=errordf, aes(x=cellularity, y=ploidy, fill=1/error)) +
        geom_point(data=minimadf, aes(x=cellularity, y=ploidy, alpha=min(error)/error), shape=16) +
        scale_fill_gradient(low="green", high="red") +
        ggtitle(sample) +
        theme(plot.title = element_text(hjust = 0.5))
    suppressMessages(ggsave(outfile, plot = p))
}



# Function to run ACE (square_model) from segmentdata
Run_ACE <- function(segmentdata, Segments,
                    sample, method, penalty,penploidy, fit_out, errormatrix_out,
                    #prows=4, # include only integer ploidy
                    prows=100, # to include any ploidy
                    ptop =5, pbottom = 1, cellularities = seq(5,100), highlightminima = TRUE){
    # calculate standard
    standard <- median(rep(segmentdata$values,segmentdata$lengths))
    # Initialize variables
    fraction <- sort(cellularities)/100
    error <- c()
    errormatrix <- matrix(nrow=(prows+1),ncol=length(fraction))
    listofploidy <- c()
    listofcellularity <- c()
    listoferrors <- c()
    for (t in seq(0,prows)){
        ploidy <- ptop-((ptop-pbottom)/prows)*t
        listofploidy <- append(listofploidy, rep(ploidy,length(fraction)))
        expected <- c()
        hexpected <- c()
        temp <- c()
        errorlist <- c()
        for (i in seq_along(fraction)) {
            for (p in seq(1,12)) {
                expected[p] <- standard*(p*fraction[i] + 2*(1-fraction[i]))/(fraction[i]*ploidy + 2*(1-fraction[i]))
                hexpected[p] <- standard*(p*fraction[i] + 1*(1-fraction[i]))/(fraction[i]*ploidy + 2*(1-fraction[i]))
            }
            for (j in seq_along(segmentdata$values)) {
                if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty))^2}
                else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty))}
                else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)*(1+abs(ploidy-2))^penploidy/(fraction[i]^penalty)}
                else {print("Not a valid method")}
            }
            for (k in seq_along(segmentdata$values)) {
                if(method=='RMSE') {temp[j+k] <- (min(abs(segmentdata$values[k]-hexpected),0.5)/(fraction[i]^penalty))^2}
                else if(method=='SMRE') {temp[j+k] <- sqrt(min(abs(segmentdata$values[k]-hexpected),0.5)/(fraction[i]^penalty))}
                else if(method=='MAE') {temp[j+k] <- min(abs(segmentdata$values[k]-hexpected),0.5)/(fraction[i]^penalty)}
                else {print("Not a valid method")}
            }
            lengths <- c(segmentdata$lengths, segmentdata$lengths)
            if(method=='RMSE') {errorlist[i] <- sqrt(sum(temp*lengths)/sum(lengths))}
            else if(method=='SMRE') {errorlist[i] <- sum(temp*lengths)/sum(lengths)^2}
            else if(method=='MAE') {errorlist[i] <- sum(temp*lengths)/sum(lengths)}
        }
        listofcellularity <- append(listofcellularity, fraction)
        listoferrors <- append(listoferrors, errorlist)
        errormatrix[t+1,] <- errorlist
        
    }
    minimat <- matrix(nrow=(prows+1),ncol=length(fraction))
    for (i in seq(1,prows+1)) {
        for (j in seq(1,length(fraction))) {
            if (i==1|i==(prows+1)) {
                minimat[i,j] <- FALSE
            } else if (j==length(fraction)) {
                minimat[i,j] <- errormatrix[i,j]==min(errormatrix[seq(i-1,i+1),seq(j-1,j)])
            } else {
                minimat[i,j] <- errormatrix[i,j]==min(errormatrix[seq(i-1,i+1),seq(j-1,j+1)])
            }
        }
    }
    round(errormatrix, digits = 10)
    round(listoferrors, digits = 10)
    errordf <- data.frame(ploidy=listofploidy,
                          cellularity=listofcellularity,
                          error=listoferrors/max(listoferrors),
                          minimum=as.vector(t(minimat)))
    minimadf <- errordf[errordf$minimum==TRUE,]
    minimadf <- minimadf[order(minimadf$error,-minimadf$cellularity),]
    
    # Plot error matrix
    Plot_errors(errordf, minimadf,sample, errormatrix_out)
    # Write fits to table
    write.table(minimadf, file = fit_out, row.names = F, quote = F, sep = "\t")
}
