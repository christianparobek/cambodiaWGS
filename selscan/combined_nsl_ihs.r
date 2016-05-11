## Need to make a combined nsl/ihs plot
## Christian Parobek
## Started 11 May 2016


####################
## LOAD LIBRARIES ##
####################

library(stringr)


######################
## USEFUL FUNCTIONS ##
######################

big.maker <- function(directory, spp_chr_metadata) {
  
  filenames <- dir(directory, pattern =".out.100bins.norm") # nsl data
  
  big <- NULL # declare object
  
  for(i in 01:length(filenames)){
    file <- read.table(paste(directory, filenames[i], sep = ""), header=FALSE)
    file$V9 <- paste("chr", str_pad(i, 2, pad = "0") , sep = "") # pad with leading zeros
    big <- rbind(big, file) # append to big file
  } # read in chr nsl files, add a chr column (padded), and cat together
  
  for(chr in spp_chr_metadata$V1){
    big[big$V9 == chr,]$V2 <- big[big$V9 == chr,]$V2 + spp_chr_metadata[spp_chr_metadata$V1 == chr, 4]
  } # adjust values in big for chr lengths
  
  return(big)
  
} # making "big" data + metadata df

cutoff.getter <- function(frac, big){
  top <- order(abs(big$V7), decreasing = TRUE)[1:round(frac*nrow(big))] # define the top samples
  cutoff <- abs(big[tail(top, 1), 7]) # get the cutoff
  return(cutoff)
} # get cutoffs for coloring the strongest hits

manhattan.plot <- function(big, cutoff, spp_chr_metadata, analysis_name){
  plot(big$V2, abs(big$V7), col = as.factor(big$V9), pch = 20, cex = 0.05, las = 1, axes = FALSE, ylab = "", xlab = "")
  points(big[abs(big$V7) >= cutoff,]$V2, abs(big[abs(big$V7) >= cutoff,]$V7), pch = 20, cex = 0.35, col = "red") # highlight the top x percent
  axis(1, at = apply(spp_chr_metadata, 1, function(x) as.numeric(x[4]) + 0.5*as.numeric(x[2])), labels = 1:14, cex.axis = 0.7) # put ticks at medians
  axis(2, las = 1, cex.axis = 0.8)
  mtext("Chromosome", 1, line = 2.5)
  mtext(analysis_name, 2, line = 2.5)
} # to plot manhattans; base style

##################################
## READ AND PROCESS INPUT FILES ##
##################################

pv_len <- read.table("nsl/pv_chrlen", header = FALSE)
pv_len$cumsum <- ave(pv_len$V2, FUN=cumsum) # add a cumulative sum of chr lengths
pv_len$toadd <- head(c(0, pv_len$cumsum), -1) # transform cumsum to be easily addable

pf_len <- read.table("nsl/pf_chrlen", header = FALSE)
pf_len$cumsum <- ave(pf_len$V2, FUN=cumsum) # add a cumulative sum of chr lengths
pf_len$toadd <- head(c(0, pf_len$cumsum), -1) # transform cumsum to be easily addable


###########################
## MAKE "BIG" DATAFRAMES ##
###########################

pv_mono_nsl_big <- big.maker("nsl/pv_mono/", pv_len)
pv_mono_ihs_big <- big.maker("ihs/pv_mono/", pv_len)
pv_all_nsl_big <- big.maker("nsl/pv_all/", pv_len)
pf_cp2_nsl_big <- big.maker("nsl/pf_cp2/", pf_len)


################################
## DEFINE CUTOFF FOR TOP HITS ##
################################

frac <- 0.005 # get top 0.5% to color them

pv_mono_nsl_cutoff <- cutoff.getter(frac, pv_mono_nsl_big)
pv_mono_ihs_cutoff <- cutoff.getter(frac, pv_mono_ihs_big)
pv_all_nsl_cutoff <- cutoff.getter(frac, pv_all_nsl_big)
pf_cp2_nsl_cutoff <- cutoff.getter(frac, pf_cp2_nsl_big)


#################
## PLOT IT ALL ##
#################



svg("combined.svg", width = 6, height = 8)
tiff("combined.tiff", width = 6, height = 8, units = "in", res = 300, compression = "lzw")

palette(c("bisque4","darkcyan"))
par(mfrow = c(3,1), mar = c(4,5,1,0))

manhattan.plot(pv_all_nsl_big, pv_all_nsl_cutoff, pv_len, expression(italic("P. vivax nS")[L]))
mtext("A", 2, line = 3.3, las = 2, cex = 1.8, padj = -3.5)

manhattan.plot(pv_mono_ihs_big, pv_mono_ihs_cutoff, pv_len, expression(paste(italic("P. vivax"), " iHS")))
mtext("B", 2, line = 3.3, las = 2, cex = 1.8, padj = -3.5)

manhattan.plot(pf_cp2_nsl_big, pf_cp2_nsl_cutoff, pf_len, expression(italic("P. falciparum nS")[L]))
mtext("C", 2, line = 3.3, las = 2, cex = 1.8, padj = -3.5)


dev.off()

