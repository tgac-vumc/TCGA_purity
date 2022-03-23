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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
suppressMessages(library(ggplot2))
#-------------------------------------------------------------------------------


# Plot errorlist
Plot_errors <- function(tempdf, minimadf,outfile) {
    p <- ggplot() +
        scale_y_continuous(name = "relative error", limits = c(0,1.05), expand=c(0,0)) +
        scale_x_continuous(name = "cellularity (%)") +
        geom_vline(xintercept = seq(from = 10, to = 100, by = 10), color = "#666666", linetype = "dashed") +
        geom_point(aes(y=errorlist, x=cellularity), data=tempdf) +
        geom_point(aes(y=rerror, x=minima), data=minimadf, color = 'red') +
        theme_classic() + theme(
                              axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), axis.text = element_text(color='black')) +
        ggtitle(paste0(sample, " - errorlist")) +
        theme(plot.title = element_text(hjust = 0.5))
    suppressMessages(ggsave(outfile, plot = p))
}

# Function to run ACE from segmentdata
Run_ACE <- function(segmentdata, Segments, ploidies, sample, method, penalty, fit_out, errorgraph_out){

    # calculate standard
    standard <- median(rep(segmentdata$values,segmentdata$lengths))
    # create empty factors
    fraction <- c()
    expected <- c()
    temp <- c()
    errorlist <- c()
    for(q in ploidies){
        #iterate over purities
        for (i in seq(5,100)) {
            fraction[i-4] <- i/100
            for (p in seq(1,12)) {
                expected[p] <- standard*(p*fraction[i-4] + 2*(1-fraction[i-4]))/(fraction[i-4]*q + 2*(1-fraction[i-4]))
            }
            for (j in seq_along(segmentdata$values)) {
                if(method=='RMSE') {temp[j] <- (min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty))^2}
                else if(method=='SMRE') {temp[j] <- sqrt(min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty))}
                else if(method=='MAE') {temp[j] <- min(abs(segmentdata$values[j]-expected),0.5)/(fraction[i-4]^penalty)}
                else {print("Not a valid method")}
            }
            if(method=='RMSE') {errorlist[i-4] <- sqrt(sum(temp*segmentdata$lengths)/sum(segmentdata$lengths))}
            else if(method=='SMRE') {errorlist[i-4] <- sum(temp*segmentdata$lengths)/sum(segmentdata$lengths)^2}
            else if(method=='MAE') {errorlist[i-4] <- sum(temp*segmentdata$lengths)/sum(segmentdata$lengths)}
            
        }
        minima <- c()
        rerror <- c()
        
        if (round(errorlist[1], digits = 10) < round(errorlist[2], digits = 10)) {
            lastminimum <- fraction[1]
            minima[1] <- fraction[1]
            rerror[1] <- errorlist[1]/max(errorlist)
        }
        
        for (l in seq(6,99)) {
            if (round(errorlist[l-4], digits = 10) < round(errorlist[l-5], digits = 10) & round(errorlist[l-4], digits = 10) < round(errorlist[l-3], digits = 10)) { 
                lastminimum <- fraction[l-4]
                minima <- append(minima,fraction[l-4])
                rerror <- append(rerror,(errorlist[l-4]/max(errorlist)))
            }
        }
        
        if (errorlist[100-4] <= errorlist[100-5]) {
            lastminimum <- fraction[100-4]
            minima <- append(minima, fraction[100-4])
            rerror <- append(rerror, errorlist[100-4]/max(errorlist))
        }
        chr <- Segments$Chromosome
        bin <- seq_along(chr)
        rlechr <- rle(chr)
        lastchr <- length(rlechr$values)
        binchrend <- c()
        currentbin <- 0
        binchrmdl <- c()
        for (i in seq_along(chr)){
            currentmiddle <- currentbin+bin[i]/2
            currentbin <- currentbin+bin[i]
            binchrend <- append(binchrend, currentbin)
            binchrmdl <- append(binchrmdl, currentmiddle)
        }
        cellularity <- seq(5,100)
        tempdf <- data.frame(cellularity,errorlist=errorlist/max(errorlist))
        minimadf <- data.frame(minima=minima*100,rerror)
        bfi <- tail(which(rerror==min(rerror)))
        fitpicker <- matrix(ncol = 16, nrow = 1)
        colnames(fitpicker) <- c("sample","likely_fit","ploidy","standard","fit_1","fit_2","fit_3","fit_4","fit_5","fit_6","fit_7","fit_8","fit_9","fit_10","fit_11","fit_12")
        fitpicker[1,1] <- sample
        fitpicker[1,2] <- minima[bfi]
        for (m in seq_along(minima)) {
            fitpicker[1,m+4] <- minima[m]
            adjustedsegments <- q + ((Segments$Segment_Mean-standard)*(minima[m]*(q-2)+2))/(minima[m]*standard)
            df <- as.data.frame(adjustedsegments[seq_along(bin)])
            colnames(df) <- "segments"
            dfna <- na.exclude(df)
            segments <- dfna$segments
            cappedsegments <- dfna[dfna$segments > 12,]
            if(length(cappedsegments)>0) {cappedsegments <- 12-0.1}
            toppedsegments <- dfna[dfna$segments <= 0,]
            if(length(toppedsegments)>0) {toppedsegments <- 0+0.1}
            line1 <- paste0("Cellularity: ", minima[m])
            line2 <- paste0("Relative error: ", round(rerror[m], digits = 3))
        }
        # Write information to file
        
        # Plot errorgraph
        if(q == 2){
            Plot_errors(tempdf, minimadf,errorgraph_out)
        }else{
            Plot_errors(tempdf, minimadf,gsub("2N",paste0(q,"N"),errorgraph_out))
        }
        
        # Write fittable
        if(q == 2){
            write.table(fitpicker, file = fit_out )
        }else{
            write.table(fitpicker, file = gsub("2N",paste0(q,"N"),fit_out), row.names = F, quote = F)

        }
    }

}


