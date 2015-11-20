# For genetic analysis of my WGS Plasmodium population(s)
# Started 10 Dec 2014
# Updated 01 May 2015
# Basics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf
# Genomics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
# Extra Commands: http://www.inside-r.org/packages/cran/adegenet/docs/.rmspaces


###################################################
################# Load Libraries ##################
###################################################

library(adegenet)
library(stringr)
library(pegas)


#####################################################
################# Define Functions ##################
#####################################################

## Function to create genlight from VCF.
## The AAKM contigs can screw with this
## so `grep -v "=AAKM" infile.vcf > outfile.vcf`
genlight.maker <- function(infile) {
  loci <- read.vcf(infile, from = 1, to = 1000000) # reads first million sites
  genlight <- new("genlight", loci) # convert data frame into genlight object
  ploidy(genlight) <- as.integer(1) # add back population information
  return(genlight)
}

## Function to assign samples to pops
## based on a list of their names
pop.definer <- function(ind_names) {
  library(stringr)
  kp <- as.numeric(str_detect(ind_names, "BB"))*2 # assign KP pop number
  bb <- as.numeric(str_detect(ind_names, "KP"))*2# assign BB pop number
  om <- as.numeric(str_detect(ind_names, "OM"))*2 # assign OM pop number
  sn <- as.numeric(str_detect(ind_names, "SN"))*2 # assign SN pop number
  tb <- as.numeric(str_detect(ind_names, "TB"))*2 # assign TB pop number
  srr <- as.numeric(str_detect(ind_names, "SRR"))*1 # assign SRR pop number
  err <- as.numeric(str_detect(ind_names, "ERR"))*1 # assign ERR pop number
  pops <- kp + bb + om + sn + tb + srr + err
  return(pops)
}

## Function to plot PCAs
pca.plotter <- function(pca, pops, x, y) {
  plot(pca[,y] ~ pca[,x], 
       col=pops, 
       pch=19, 
       axes=FALSE, 
       xlab=paste("PC", x, sep = ""),
       ylab=paste("PC", y, sep = ""),
  )
  axis(1)
  axis(2)
}


###############################################
################# PROCESS PCA #################
###############################################

## read in PF data
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/cp-kh/")
pf_gl <- genlight.maker("joint.vcf") # make genlight
pf_pops <- pop.definer(indNames(pf_gl)) # define pops
pf_pca <- glPca(pf_gl) # calculate PCA
pf_pca_jit <- as.data.frame(jitter(pf_pca$scores, factor=100))
# and make a copy with some jitter

# remove the 65 that are from the Gambia, which seem to have a PC1 of >15
sansGambia <- pf_pca_jit[!pf_pca_jit$PC1 > 18,]
sansGambiaCol <- pf_pops[!pf_pca_jit$PC1 > 18]


######################################################
################# PLOT PCA & BOXPLOT #################
######################################################

par(mfrow = c(1,2))

# plot the pca
pca.plotter(sansGambia, sansGambiaCol, 1, 3)
text(c(-12, -18, -3, -3), c(21, -8, -7, 4), labels = c("CP4", "CP3", "CP2", "CP1"))
mtext("A", 2, las = 2, cex = 2.5, at = 21, line = 3)

# plot the boxplot
boxplot(boots[9,2,], boots[5,2,], boots[6,2,], boots[7,2,], boots[8,2,], axes = FALSE, ylab = "CP2 vs. KH Group\n(Euclidean Distance)", xlab = "KH Groups")
axis(1, at = 1:5, labels = c("KHA", "KH1", "KH2", "KH3", "KH4"))
axis(2)
mtext("B", 2, las = 2, cex = 2.5, at = 5.6, line = 3)











### ggPlot fail

p <- ggplot(sansGambia, aes(x=PC1, y=PC3)) + 
  geom_point(alpha=0.6, size = 4, col = sansGambiaCol) + 
  theme_bw() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none") +
  scale_size_continuous(range = c(2,9)) +
  labs(
    x = paste("PC", 1, " - ", round(pf_pca$eig[1]/sum(pf_pca$eig)*100), "% of the Variance", sep = ""),
    y = paste("PC", 3, " - ", round(pf_pca$eig[2]/sum(pf_pca$eig)*100), "% of the Variance", sep = "")) + 
  annotate("text", x = -12, y = 21, label = "CP4", size = 10) +
  annotate("text", x = -18, y = -8, label = "CP3", size = 10) +
  annotate("text", x = -3, y = -6, label = "CP2", size = 10) +
  annotate("text", x = -3, y = 4, label = "CP1", size = 10)

# rearrange boxplot


# plot boxplot
ggplot(data=melt(boots[c(9,5:8),2,]), aes(as.factor(X1), value)) + 
  geom_jitter() + 
  geom_boxplot() + 
  theme_classic() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)) + 
  labs(
    x = "",
    y = "CP2 vs. KH Group\n(Euclidean Distance)"
  )

grid.arrange(p, p, ncol=2)


layout(matrix(1:2, nrow = 1, ncol = 2, byrow = TRUE))
