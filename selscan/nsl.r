## Analyze and plot un-normalized nSL data
## Started 21 December 2015
## Christian P

####################
## LOAD LIBRARIES ##
####################

library(stringr)
args <- commandArgs(TRUE)

##################################
## READ AND PROCESS INPUT FILES ##
##################################

len <- read.table(paste(args[2], "_chrlen", sep = ""), header = FALSE) # read in chr len metadata
len$cumsum <- ave(len$V2, FUN=cumsum) # add a cumulative sum of chr lengths
len$toadd <- head(c(0, len$cumsum), -1) # transform cumsum to be easily addable

filenames <- dir(args[1], pattern =".res.nsl.out") # nsl data
big <- NULL # declare object
for(i in 01:length(filenames)){
  file <- read.table(paste(args[1], "/", filenames[i], sep = ""), header=FALSE)
  file$V7 <- paste("chr", str_pad(i, 2, pad = "0") , sep = "") # pad with leading zeros
  big <- rbind(big, file) # append to big file
} # read in chr nsl files, add a chr column (padded), and cat together
for(chr in len$V1){
  big[big$V7 == chr,]$V2 <- big[big$V7 == chr,]$V2 + len[len$V1 == chr, 4]
} # adjust values in big for chr lengths


################################
## DEFINE CUTOFF FOR TOP HITS ##
################################

frac <- 0.005 # get top 0.1%
top <- order(abs(big$V6), decreasing = TRUE)[1:round(frac*nrow(big))] # define the top samples
cutoff <- abs(big[tail(top, 1), 6]) # get the cutoff


####################
## MANHATTAN PLOT ##
####################

svg(paste(args[1], ".svg", sep = ""), width = 8, height = 5)
palette(c("bisque4","darkcyan"))
plot(big$V2, abs(big$V6), ylim=c(0,2.5), col = as.factor(big$V7), pch = 20, cex = 0.05, las = 1, axes = FALSE, ylab = "", xlab = "")
points(big[abs(big$V6) >= cutoff,]$V2, abs(big[abs(big$V6) >= cutoff,]$V6), pch = 20, cex = 0.25, col = "red") # highlight the top x percent
axis(1, at = apply(len, 1, function(x) as.numeric(x[4]) + 0.5*as.numeric(x[2])), labels = 1:14, cex.axis = 0.7) # put ticks at medians
axis(2, las = 1, cex.axis = 0.8)
mtext("Chromosome", 1, line = 2.5)
mtext(expression(italic("nS")[L]), 2, line = 2.5)
dev.off()

