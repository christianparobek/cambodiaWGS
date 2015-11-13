## Started 25 September 2015
## To plot AFS for my collected and simulated data.

#####################################
######## PLOT PF AND PV AFS #########
#####################################

## read in pf
pf_afs <- read.table("pf_afs.txt", header = TRUE)
pf_count <- with(pf_afs, rbind(pf_syn, pf_nonsyn, pf_gen, pf_ige))

## read in pv
pv_afs <- read.table("pv_afs.txt", header = TRUE)
pv_count <- with(pv_afs, rbind(pv_syn, pv_nonsyn, pv_gen, pv_ige))

###################################
## MAKE MAIN PICTURE - Pf AND Pv ##
###################################

svg(filename = "pf_afs.svg", width=8, height=6)

## set up plot space
par(mfrow = c(2, 1),
    oma = c(5,4,0,0),
    mar = c(0,0,2,0))

## plot PDF
barplot(pf_count, beside = TRUE, col=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
        ylab="", xlab="", border = NA, axes = FALSE)
axis(2, at = c(0, 0.1, 0.2, 0.3), las = 2)
axis(1, at = 1:39*5-2, labels = 1:39, tick = FALSE, cex.axis = 0.6, line = -0.5)
mtext("Fraction of each SNP Class", side = 2, line = 2.8)
lines(1:nrow(pf_afs)*5-2, pf_afs$ms_78_1000_s7000_folded, lwd=3, col = "black", lty=2) # showing sims as lines here
#lines(1:nrow(pf_afs)*5-2, pf_afs$ms_78_1000_t368_folded, lwd=3, col = "black", lty=2) # showing sims as lines here
title(expression(italic("P. falciparum")), line = -1.5, cex.main = 1.5)

barplot(pv_count, beside = TRUE, col=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
        ylab="", xlab="", border = NA, axes = FALSE)
axis(2, at = c(0, 0.2, 0.4), las = 2)
axis(1, at = 1:35*5-2, labels = 1:35, tick = FALSE, cex.axis = 0.6, line = -0.5)
mtext("Derived Allele Frequency", side = 1, line = 1.8)
mtext("Fraction of each SNP Class", side = 2, line = 2.8)
lines(1:nrow(pv_afs)*5-2, pv_afs$ms_70_1000_t11100_folded, lwd=3, col = "black", lty=2) # showing sims as lines here
legend(120, 0.33, 
       legend = c("Synonymous", "Nonsynonymous", "Genic", "Intergenic", "Simulation"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6", "black"), 
       lty=c(NA, NA, NA, NA, 2), 
       pch=c(15, 15, 15, 15, NA),
       pt.cex=c(1.8, 1.8, 1.8, 1.8, NA), lwd=3, bty = "n", cex=1)
title(expression(italic("P. vivax")), line = -1.5, cex.main = 1.5)

dev.off()



#####################################
######## FOR SUPPLEMENTARY ##########
#####################################

## plot Pf alone
barplot(pf_count, beside = TRUE, col=c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
        ylab="", xlab="", border = NA, axes = FALSE)
axis(2, at = c(0, 0.1, 0.2, 0.3), las = 2)
axis(1, at = 1:39*5-2, labels = 1:39, tick = FALSE, cex.axis = 0.6, line = -0.5)
mtext("Fraction of each SNP Class", side = 2, line = 2.8)
lines(1:nrow(pf_afs)*5-2, pf_afs$ms_78_1000_s15000_folded, lwd=3, col = "black", lty=2) # showing sims as lines here
title(expression(italic("P. falciparum")), line = -1.5, cex.main = 1.5)
mtext("Derived Allele Frequency", side = 1, line = 1.8)
legend(150, 0.20, 
       legend = c("Synonymous", "Nonsynonymous", "Genic", "Intergenic", "Simulation"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6", "black"), 
       lty=c(NA, NA, NA, NA, 2), 
       pch=c(15, 15, 15, 15, NA),
       pt.cex=c(1.8, 1.8, 1.8, 1.8, NA), lwd=3, bty = "n", cex=1)