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
library(ggplot2)
library(gridExtra)

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
  kp <- as.numeric(str_detect(ind_names, "BB"))*1 # assign KP pop number
  bb <- as.numeric(str_detect(ind_names, "KP"))*2 # assign BB pop number
  om <- as.numeric(str_detect(ind_names, "OM"))*3 # assign OM pop number
  sn <- as.numeric(str_detect(ind_names, "SN"))*3 # assign SN pop number
  tb <- as.numeric(str_detect(ind_names, "TB"))*3 # assign TB pop number
  srr <- as.numeric(str_detect(ind_names, "SRR"))*4 # assign SRR pop number
  err <- as.numeric(str_detect(ind_names, "ERR"))*4 # assign ERR pop number
  pops <- kp + bb + om + sn + tb + srr + err
  return(pops)
}

## Function to plot eigenplots
eig.plotter <- function(pca) {
  barplot(pca$eig, xlab = "", ylab = "Variance")
}

## Function to plot PCAs
pca.plotter <- function(pca, pops, x, y) {
  plot(jitter(pca$scores[,y], factor=700) ~ jitter(pca$scores[,x], factor=700), 
     col=pops, 
     pch=19, 
     axes=FALSE, 
     xlab=paste("PC", x, " - ", round(pca$eig[x]/sum(pca$eig)*100), "% of the Variance", sep = ""),
     ylab=paste("PC", y, " - ", round(pca$eig[y]/sum(pca$eig)*100), "% of the Variance", sep = ""),
  )
  axis(1)
  axis(2)
}




##############################################################
################# ggPlot Style Pv and Pf pic #################
##############################################################

## read in PV data
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/pv/variants/")
pv_cam_gl <- genlight.maker("our_goods_UG.pass.vcf.2") # make genlight; vcf must not contain extra contigs
pv_cam_pops <- pop.definer(indNames(pv_cam_gl)) # define pops
pv_cam_pca <- glPca(pv_cam_gl) # calculate PCA

## read in PF data
setwd("/run/user/1001/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/pf/variants/")
pf_cam_gl <- genlight.maker("our_goods_UG.pass.vcf") # make genlight
pf_cam_pops <- pop.definer(indNames(pf_cam_gl)) # define pops
pf_cam_pca <- glPca(pf_cam_gl) # calculate PCA
pf_cam_pca_jit <- as.data.frame(jitter(pf_cam_pca$scores, factor=700))
# and make a copy with some jitter



pv <- ggplot(as.data.frame(pv_cam_pca$scores), aes(x=PC1, y=PC2)) + 
  geom_point(alpha=0.6, size = 4) + 
  theme_bw() +
  ggtitle(expression(italic("P. vivax"))) +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  scale_size_continuous(range = c(2,9)) +
  labs(
    x = paste("PC", 1, " - ", round(pv_cam_pca$eig[1]/sum(pv_cam_pca$eig)*100), "% of the Variance", sep = ""),
    y = paste("PC", 2, " - ", round(pv_cam_pca$eig[2]/sum(pv_cam_pca$eig)*100), "% of the Variance", sep = "")
  )


pf <- ggplot(pf_cam_pca_jit, aes(x=PC1, y=PC2)) + 
  geom_point(alpha=0.6, size = 4) + 
  theme_bw() +
  ggtitle(expression(italic("P. falciparum"))) +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  scale_size_continuous(range = c(2,9)) +
  labs(
    x = paste("PC", 1, " - ", round(pf_cam_pca$eig[1]/sum(pf_cam_pca$eig)*100), "% of the Variance", sep = ""),
    y = paste("PC", 2, " - ", round(pf_cam_pca$eig[2]/sum(pf_cam_pca$eig)*100), "% of the Variance", sep = "")
  )


grid.arrange(pf, pv, ncol=2)



##############################################################
################# Plasmodium vivax Cambodia ##################
##############################################################

### PREPARE DATA FOR ANALYSIS ###
pv_cam_gl <- genlight.maker("/run/user/1001/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/pv/variants/our_goods_UG.pass.vcf") # make genlight
pv_cam_pops <- pop.definer(indNames(pv_cam_gl)) # define pops
pv_cam_pca <- glPca(pv_cam_gl) # calculate PCA

### PLOT EIGENVALUES ###
eig.plotter(pv_cam_pca)
title(substitute(paste("Cambodia ", italic('P. vivax'), " Eigenvalues" )), line = 0.5, cex.main = 1.5)

### PLOT PCA PICTURE ###
pca.plotter(pv_cam_pca, pv_cam_pops, 1, 2)
legend(10, -40, legend = c("BB", "KP", "OM"), col = c("red", "black", "green3"), pch=19, bty="n", cex=1.5)
title(substitute(paste("Cambodia ", italic('P. vivax'), " PCA" )), line = -0.5, cex.main=1.5)

### DAPC ANALYSIS ###


############################################################
################# Plasmodium vivax Global ##################
############################################################

### PREPARE DATA FOR ANALYSIS ###
pv_all_gl <- genlight.maker("all_goods_pv.pass.str") # make genlight
pv_all_pops <- pop.definer(indNames(pv_all_gl)) # define pops
pv_all_pca <- glPca(pv_all_gl) # calculate PCA

### PLOT EIGENVALUES ###
eig.plotter(pv_all_pca)
title(substitute(paste("Global ", italic('P. vivax'), " Eigenvalues" )), line = 0.5, cex.main = 1.5)

### PLOT PCA PICTURE ###
pca.plotter(pv_all_pca, pv_all_pops, 1, 2)
legend(-60, 85, legend = c("BB", "KP", "OM", "Global"), col = c("red", "black", "green3", "blue"), pch=19, bty="n", cex=1.5)
title(substitute(paste("Global ", italic('P. vivax'), " PCA" )), line = -1.5, cex.main=1.5)


##############################################################
############## Plasmodium falciparum Cambodia ################
##############################################################

### PREPARE DATA FOR ANALYSIS ###
pf_cam_gl <- genlight.maker("/run/user/1001/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/pf/variants/our_goods_UG.pass.vcf") # make genlight
pop(pf_cam_gl) <- pop.definer(indNames(pf_cam_gl)) # define pops
pf_cam_pca <- glPca(pf_cam_gl) # calculate PCA

### PLOT EIGENVALUES ###
eig.plotter(pf_cam_pca)
title(substitute(paste("Cambodia ", italic('P. falciparum'), " Eigenvalues" )), line = 0.5, cex.main = 1.5)

### PLOT PCA PICTURE ###
pca.plotter(pf_cam_pca, pop(pf_cam_gl), 1, 2)
title(substitute(paste("Cambodia ", italic('P. falciparum'), " PCA" )), line = -0.5, cex.main=1.2)
legend(-12, -2, legend = c("BB", "KP", "OM"), col = c("red", "black", "green3"), pch=19, bty="n", cex=1.2)

### DAPC ###
pf_cam_dapc <- dapc(pf_cam_gl, n.pca = 10, n.da = 1)
scatter(pf_cam_dapc, scree.da=FALSE)


### Minor Allele Frequency ###
pf_cam_maf <- glMean(pf_cam_gl)
hist(pf_cam_maf, proba = TRUE, col =  "gold", xlab = "Allele frequencies", main = "Minor allele frequencies")
pf_cam_maf_density <- density(pf_cam_maf)
lines(pf_cam_maf_density$x, pf_cam_maf_density$y, lwd = 3)


##############################################################
############### Plasmodium falciparum Global #################
##############################################################

### PREPARE DATA FOR ANALYSIS ###
pf_all_gl <- genlight.maker("all_goods_pf.pass.str") # make genlight
pf_all_pops <- pop.definer(indNames(pf_all_gl)) # define pops
pf_all_pca <- glPca(pf_all_gl) # calculate PCA

### PLOT EIGENVALUES ###
eig.plotter(pf_all_pca)
title(substitute(paste("Global ", italic('P. falciparum'), " Eigenvalues" )), line = 0.5, cex.main = 1.5)

### PLOT PCA PICTURE ###
pca.plotter(pf_all_pca, pf_all_pops, 1, 2)
legend(-50, -30, legend = c("BB", "KP", "OM", "Gambia"), col = c("red", "black", "green3", "blue"), pch=19, bty="n", cex=1.5)
title(substitute(paste("Global ", italic('P. falciparum'), " PCA" )), line = -1.5, cex.main=1.5)
