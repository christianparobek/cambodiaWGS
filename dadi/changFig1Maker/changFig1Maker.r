## Started 25 September 2015
## To plot AFS for my collected and simulated data.

#####################################
######## PLOT FALCIPARUM AFS ########
#####################################

pf_afs <- read.table("pf/pf_afs.txt", header = TRUE)

count <- with(pf_afs, rbind(pf_syn, pf_nonsyn, pf_gen, pf_ige))

pdf("pf/pf_afs.pdf",width=6,height=3.7)
barplot(count, beside = TRUE, col=c("blue","green","yellow","black"),
        legend = c("Synonymous", "Nonsynonymous", "Genic", "Intergenic", "Simulation 1", "Simulation 2"), ylab="Fraction of each SNP Class", xlab="Derived Allele Frequency",
        border = NA, names.arg=1:35, cex.names=0.6, cex.axis=0.8)
lines(1:nrow(pf_afs)*7, pf_afs$ms_70_1000_s15000_folded, lwd=2) # showing sims as lines here
lines(1:nrow(pf_afs)*7, pf_afs$ms_70_10M_t000478_folded, lwd=2) # showing sims as lines here
dev.off()

#

#####################################
########## PLOT VIVAX AFS ###########
#####################################

pv_afs <- read.table("pv/pv_afs.txt", header = TRUE)

count <- with(pv_afs, rbind(pv_syn, pv_nonsyn, pv_gen, pv_ige, ms_70_1000_s80000_folded, ms_70_1000_s15000_folded, ms_70_10M_t000478_folded))

pdf("pv/pv_afs.pdf",width=5,height=5)
barplot(count, beside = TRUE, col=c("blue","green","yellow","grey" ,"black", "orange","red"),
        legend = c("Synonymous", "Nonsynonymous", "Genic", "Intergenic", "Simulation 1", "Simulation 2", "Simulation 3"), ylab="Fraction of each SNP Class", xlab="Derived Allele Frequency",
        border = NA, names.arg=1:35, cex.names=0.6, cex.axis=0.8)
# showing sims as bars here. need to harmonize with pf, one way or the other.
dev.off()