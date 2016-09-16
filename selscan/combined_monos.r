## Need to make a combined figure of 
## all the pv mono things: nsl, ehh, and bifur
## Christian Parobek
## 12 May 2016


####################
## LOAD LIBRARIES ##
####################

library(stringr)

library(rehh)
#source("../rehh/mybifur.r") # a modded bifurcation.diagram, w/o x-axis
source("mybifur.r") # a modded bifurcation.diagram, w/o x-axis




################################################
################################################
######### READ IN AND PROCESS NSL DATA #########
################################################
################################################

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
  
} # making "big" data + metadata df; for nsl

cutoff.getter <- function(frac, big){
  top <- order(abs(big$V7), decreasing = TRUE)[1:round(frac*nrow(big))] # define the top samples
  cutoff <- abs(big[tail(top, 1), 7]) # get the cutoff
  return(cutoff)
} # get cutoffs for coloring the strongest hits; for nsl

manhattan.plot <- function(big, cutoff, spp_chr_metadata){
  plot(big$V2, abs(big$V7), col = as.factor(big$V9), pch = 20, cex = 0.05, las = 1, axes = FALSE, ylab = "", xlab = "")
  points(big[abs(big$V7) >= cutoff,]$V2, abs(big[abs(big$V7) >= cutoff,]$V7), pch = 20, cex = 0.35, col = "#d95f02") # highlight the top x percent
  axis(1, at = apply(spp_chr_metadata, 1, function(x) as.numeric(x[4]) + 0.5*as.numeric(x[2])), labels = 1:14, cex.axis = 1) # put ticks at medians
  axis(2, las = 1, cex.axis = 1, at = 0:6, labels = c("0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0"))
} # base style; for nsl


##################################
## READ AND PROCESS INPUT FILES ##
##################################

pv_len <- read.table("nsl/pv_chrlen", header = FALSE)
pv_len$cumsum <- ave(pv_len$V2, FUN=cumsum) # add a cumulative sum of chr lengths
pv_len$toadd <- head(c(0, pv_len$cumsum), -1) # transform cumsum to be easily addable


###########################
## MAKE "BIG" DATAFRAMES ##
###########################

pv_mono_nsl_big <- big.maker("nsl/pv_mono/", pv_len)


################################
## DEFINE CUTOFF FOR TOP HITS ##
################################

frac <- 0.005 # get top 0.5% to color them
pv_mono_nsl_cutoff <- cutoff.getter(frac, pv_mono_nsl_big)



################################################
################################################
######### READ IN AND PROCESS EHH DATA #########
################################################
################################################

################################
######### READ IN DATA #########
################################

data <- read.table("ehh/data/pv_chr14_797870.ehh.locus1618.out", header = FALSE)

pos <- 797870

################################
######### PREPARE DATA #########
################################

names(data) <- c("physicalPos", "geneticPos", "EHH1", "EHH0", "fullSampEHH")

data$counter <- pos + data$physicalPos

start <- round(data$counter[1]/1000)
middle <- round(pos/1000)
end <- round(data$counter[nrow(data)]/1000)


################################
##### CALC EHH REGION SIZE #####
################################

EHH1coords <- data[data$EHH1 >= 0.5,]$physicalPos
EHH0coords <- data[data$EHH0 >= 0.5,]$physicalPos
# get coordinates of SNPs 
# arbitrarily chose EHH > 0.5 to be elevated

EHH1region <- abs(min(EHH1coords)) + abs(max(EHH1coords))
EHH0region <- abs(min(EHH0coords)) + abs(max(EHH0coords))
# calculate region of elevated EHH



################################################
################################################
######## READ IN AND PROCESS BIFUR DATA ########
################################################
################################################

################################
######### READ IN DATA #########
################################

hap <- data2haplohh(hap_file="../rehh/14.hap", map_file="../rehh/14.inp")
hap <- data2haplohh(hap_file="14.hap", map_file="14.inp")


#################
## PLOT IT ALL ##
#################

svg("combined_monos.svg", width = 6.75, height = 6)

layout(matrix(c(1,1,2,3,2,4), 3, 2, byrow = TRUE), widths = c(1,1), heights = c(2.2,0.75,1.25))
palette(c("gray50","#1b9e77")) #bisque4

###### ADD THE NSL
par(mar = c(8,5,1,1))
manhattan.plot(pv_mono_nsl_big, pv_mono_nsl_cutoff, pv_len)
mtext("Chromosome", 1, line = 2.5, cex = 0.8)
mtext(expression(italic("nS"[L])), 2, line = 2.5, cex = 0.8)
arrows(22200000, 4.4, 21300000, 4.4, lwd = 1, length = 0.05, col = "#d95f02") # chr14 HP1
text(x=21900000, y=4.6, labels = "HP1/SET10", cex = 0.7, col = "#d95f02")
arrows(19350000, 6.2, 20150000, 6.2, lwd = 1, length = 0.05, col = "#d95f02") # chr14 AP2
text(x=19650000, y=6.35, labels = "AP2", cex = 0.7, col = "#d95f02")
arrows(18800000, 3.8, 18000000, 3.8, lwd = 1, length = 0.05, col = "#d95f02") # chr13 ABC
text(x=18500000, y=4.0, labels = "ABC", cex = 0.7, col = "#d95f02")
arrows(15900000, 3.6, 16700000, 3.6, lwd = 1, length = 0.05, col = "#d95f02") # chr12 MRP2
text(x=16300000, y=3.8, labels = "MRP2", cex = 0.7, col = "#d95f02")
arrows(15200000, 3.3, 16000000, 3.3, lwd = 1, length = 0.05, col = "#d95f02") # chr12 AP2
text(x=15500000, y=3.5, labels = "AP2", cex = 0.7, col = "#d95f02")
arrows(12200000, 3.6, 13000000, 3.6, lwd = 1, length = 0.05, col = "#d95f02") # chr11 AP2
text(x=12500000, y=3.8, labels = "SET", cex = 0.7, col = "#d95f02")
arrows(13350000, 3.6, 14150000, 3.6, lwd = 1, length = 0.05, col = "#d95f02") # chr11 AP2
text(x=13650000, y=3.8, labels = "AP2", cex = 0.7, col = "#d95f02")
arrows(10400000, 4.45, 11200000, 4.45, lwd = 1, length = 0.05, col = "#d95f02") # chr10 MDR1
text(x=10800000, y=4.65, labels = "MDR1", cex = 0.7, col = "#d95f02")
arrows(9800000, 3.55, 10600000, 3.55, lwd = 1, length = 0.05, col = "#d95f02") # other chr09 AP2
text(x=10100000, y=3.75, labels = "AP2", cex = 0.7, col = "#d95f02")
arrows(4200000, 4.25, 3400000, 4.25, lwd = 1, length = 0.05, col = "#d95f02") # chr04 SERA
text(x=3800000, y=4.45, labels = "SERA", cex = 0.7, col = "#d95f02")
arrows(000000, 4.05, 800000, 4.05, lwd = 1, length = 0.05, col = "#d95f02") # other chr02 MRP1
text(x=400000, y=4.25, labels = "MRP1", cex = 0.7, col = "#d95f02")
mtext("A", 2, line = 3.3, las = 2, cex = 1.6, padj = -4)


###### ADD THE EHH
par(mar = c(5,5,0,0))
plot(data$counter, data$EHH1, type = "n", ylab = "", xlab = "", axes = FALSE, ylim = c(0,1))
lines(data$counter, data$EHH0, col = "#1b9e77", lwd = 2)
lines(data$counter, data$EHH1, col = "#d95f02", lwd = 2)
axis(1, 
     at = c(pos+data$physicalPos[1], pos, pos+data$physicalPos[nrow(data)]),
     labels = c(start, middle, end))
axis(2, las = 2)
mtext("Chromosome 14 Position (kb)", 1, line = 2.5, cex = 0.8)
mtext("EHH", 2, line = 2.5, cex = 0.8)
mtext("B", 2, line = 3.3, las = 2, cex = 1.6, padj = -5.5)

###### ADD THE BIFUR
par(mar = c(0,1,0,2))
my.bifurcation.diagram(hap,mrk_foc=1511,all_foc=2,nmrk_l=110,nmrk_r=100, main = "", limhapcount = 2, refsize = 0.5, linecol = "#d95f02")
text(par("usr")[2]*0.998, mean(par("usr")[3:4])*1.49, labels = "Selected Allele", srt = 270, cex = 1)
mtext("C", 2, line = 0, las = 2, cex = 1.6, padj = -2.5)


par(mar = c(5,1,0,2))
my.bifurcation.diagram(hap,mrk_foc=1511,all_foc=1,nmrk_l=110,nmrk_r=100, main = "", refsize = 0.5, linecol = "#1b9e77")
text(par("usr")[2]*0.998, mean(par("usr")[3:4])*0.82, labels = "Unselected Allele", srt = 270, cex = 1)


## Plot the axis
min <- round(par("usr")[1]/1000 + 5) # get x min
max <- round(par("usr")[2]/1000 - 5) # get x max
middle <- round(797870/1000)
axis(1, outer = FALSE, at = c(min*1000, 797870, max*1000), labels = c(min, middle, max))

## Add text to plot
mtext("Chromosome 14 Position (kb)", side = 1, line = 2.5, cex = 0.8)

dev.off()

