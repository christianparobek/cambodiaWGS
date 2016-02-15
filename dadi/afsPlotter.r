## Started 25 September 2015
## To plot AFS for my collected and simulated data.
## Need to mod input files in the following ways:
##    1) Remove all lines with dashes
##    2) For mono/cp2 & pf/pv all comparisons, lengthen shorter one with tabbed zeros
## Usage: "Rscript dadi/afsPlotter.r <path/to/pf_all.afs> <cp2> <pv_all> <pv_mono> <all.svg> <subset.svg>"

#####################################
######## PARSE COMMAND LINE #########
#####################################

args = commandArgs(trailingOnly=TRUE)


#####################################
########### READ IN DATA ############
#####################################

pf_all <- read.table(args[1], header = TRUE) ## read in pf all
pf_cp2 <- read.table(args[2], header = TRUE) ## read in pf cp2

pv_all <- read.table(args[3], header = TRUE) ## read in pv all
pv_mono <- read.table(args[4], header = TRUE) ## read in pv mono


#####################################
######### PF & PV ALL DATA ##########
#####################################

svg(filename = args[5], width=8, height=5)

## set up plot space
par(mfrow = c(2, 1),
    oma = c(5,4,0,0),
    mar = c(0,0,2,0))

## plot svg
barplot(t(pf_all[,1:4]), beside = TRUE, col=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
        ylab="", xlab="", border = NA, axes = FALSE, ylim = c(0,0.5))
axis(2, las = 1)
axis(1, at = 1:40*5-2, labels = 1:40, tick = FALSE, cex.axis = 0.5, line = -0.5)
lines(1:nrow(pf_all)*5-2, pf_all$ms, lwd=2, col = "black", lty=2) # showing sims as lines here
title(expression(italic("P. falciparum")), line = -1.5, cex.main = 1.5)

barplot(t(pv_all[,1:4]), beside = TRUE, col=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
        ylab="", xlab="", border = NA, axes = FALSE, ylim = c(0,0.5))
axis(2, las = 1)
axis(1, at = 1:35*5-2, labels = 1:35, tick = FALSE, cex.axis = 0.5, line = -0.5)
mtext("Minor Allele Frequency", side = 1, line = 1.8)
lines(1:36*5-2, pv_all$ms[1:36], lwd=2, col = "black", lty=2) # showing sims as lines here
title(expression(italic("P. vivax")), line = -1.5, cex.main = 1.5)

mtext("Fraction of each SNP Class", side = 2, outer = TRUE, line = 2.8)
legend("right", 
       legend = c("Synonymous", "Nonsynonymous", "Genic", "Intergenic", "Simulation"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6", "black"), 
       lty=c(NA, NA, NA, NA, 2), 
       pch=c(15, 15, 15, 15, NA),
       pt.cex=c(1.8, 1.8, 1.8, 1.8, NA), lwd=2, bty = "n", cex=1)

dev.off()



###################################
####### PF & PV SUBSET DATA #######
###################################

svg(filename = args[6], width=8, height=6)

## set up plot space
par(mfrow = c(2, 1),
    oma = c(5,4,0,0),
    mar = c(0,0,2,0))

## plot svg
barplot(t(pf_cp2[,1:4]), beside = TRUE, col=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
        ylab="", xlab="", border = NA, axes = FALSE, ylim = c(0,0.5))
axis(2, las = 1)
axis(1, at = 1:9*5-2, labels = 1:9, tick = FALSE, cex.axis = 0.8, line = -0.5)
lines(1:10*5-2, pf_cp2$ms[1:10], lwd=3, col = "black", lty=2) # showing sims as lines here
title(expression(paste(italic("P. falciparum"), " (founder population)")), line = -1.5, cex.main = 1.5)

barplot(t(pv_mono[,1:4]), beside = TRUE, col=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
        ylab="", xlab="", border = NA, axes = FALSE, ylim = c(0,0.5))
axis(2, las = 1)
axis(1, at = 1:14*5-2, labels = 1:14, tick = FALSE, cex.axis = 0.8, line = -0.5)
mtext("Minor Allele Frequency", side = 1, line = 1.8)
lines(1:nrow(pv_mono)*5-2, pv_mono$ms, lwd=3, col = "black", lty=2) # showing sims as lines here
title(expression(paste(italic("P. vivax"), " (monoclonals)")), line = -1.5, cex.main = 1.5)

mtext("Fraction of each SNP Class", side = 2, outer = TRUE, line = 2.8)
legend("right", 
       legend = c("Synonymous", "Nonsynonymous", "Genic", "Intergenic", "Simulation"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6", "black"), 
       lty=c(NA, NA, NA, NA, 2), 
       pch=c(15, 15, 15, 15, NA),
       pt.cex=c(1.8, 1.8, 1.8, 1.8, NA), lwd=3, bty = "n", cex=1)

dev.off()

