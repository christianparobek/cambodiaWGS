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
library(abind)


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
## based on their CP group number
pop.definer <- function(ind_names) {
  library(stringr)
  cp1 <- as.numeric(str_detect(ind_names, "OM352|SN003|SN019|SN032|SN043|SN060|SN066|SN072|SN082|SN083|SN093|SN099|KP054|KP065"))*1 # assign KP pop number
  cp2 <- as.numeric(str_detect(ind_names, "KP001|KP004|BB085|KP059|KP073|SN076|SN079|SN091|SN103|SN109|SN078|KP027|KP030|KP062|SN064|SN097|SN107|SN111|BB084|BB059|BB080"))*2# assign BB pop number
  cp3 <- as.numeric(str_detect(ind_names, "BB052|BB082|BB068|BB069|SN022|SN035|SN038|SN039|SN044|SN048|SN052|SN057|SN058|SN061|SN084|SN085|SN092|SN095|SN105|SN117"))*3 # assign OM pop number
  cp4 <- as.numeric(str_detect(ind_names, "SN015|SN016|SN030|SN031|SN042|SN046|SN047|SN062|SN063|SN071|SN074|SN075|SN077|SN081|SN086|SN096|SN098|SN101|SN102|SN106|SN108|SN114|SN116"))*4 # assign SN pop number
  srr <- as.numeric(str_detect(ind_names, "SRR|ERR|TB"))*5 # assign SRR and ERR pop number
  pops <- cp1 + cp2 + cp3 + cp4 + srr
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
  axis(2, las = 2)
}

## Calculate the allele frequencies by SNP, given a VCF
af.calc <- function(vcf) {
  
  data <- vcf[-c(1:9)] # get rid of first nine rows
  
  altCT <- as.data.frame(sapply(data, function(x) str_extract(x, "[0123456789]:")))
  altCT <- sapply(altCT, function(x) str_extract(x, "[0123456789]+")) # clean out extra chars
  altCT <- apply(altCT, c(1,2), as.numeric) # convert to numeric matrix
  
  vector <- apply(altCT, 1, mean)
  
  return(vector)
} ## 0 is ref, 1 is alt for VCF


###################################################
############# LIFTOVER - READ IN DATA #############
###################################################

## READ IN THE Pf AND Pv MULTIVCFs
cp1 <- read.table("liftover/vcfs/cp1_intersected.vcf", header=FALSE)
cp2 <- read.table("liftover/vcfs/cp2_intersected.vcf", header=FALSE)
cp3 <- read.table("liftover/vcfs/cp3_intersected.vcf", header=FALSE)
cp4 <- read.table("liftover/vcfs/cp4_intersected.vcf", header=FALSE)

miotto_af <- read.table("liftover/beds/liftover_intersected.bed", header = FALSE)


###################################################
############ LIFTOVER - BUILD MATRIX ##############
###################################################

af_matrix <- rbind(af.calc(cp1), 
                   af.calc(cp2), 
                   af.calc(cp3), 
                   af.calc(cp4),
                   miotto_af$V5,
                   miotto_af$V6,
                   miotto_af$V7,
                   miotto_af$V8,
                   miotto_af$V9)

###################################################
######### LIFTOVER - GET POINT ESTIMATE ###########
###################################################

dist(af_matrix)


###################################################
############## LIFTOVER - BOOTSTRAP ###############
###################################################

nboot <- 1000
ncols <- 100 #ncol(af_matrix)
counter <- 0
denominator <- 0
boots <- NULL

for (i in 1:nboot) {
  boot_cols <- sample(1:ncol(af_matrix), size = ncols, replace = TRUE)
  mat <- as.matrix(dist(af_matrix[, boot_cols]))
  boots <- abind(boots, mat, along = 3)
  mat[mat == 0] <- NA # replace the zeros (diagonals) with NA
  mat[1,2] <- NA # replace the 
  denominator <- denominator + sum(as.numeric((mat[9,2] <= mat[,2])), na.rm = TRUE)
  counter <- counter + sum(as.numeric(!(mat[9,2] <= mat[,2])), na.rm = TRUE)
}


###############################################
########## PCA - READ & PROCESS DATA ##########
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

# replace colors
sansGambiaCol[sansGambiaCol == 1] <- "deepskyblue"
sansGambiaCol[sansGambiaCol == 2] <- "brown1"
sansGambiaCol[sansGambiaCol == 3] <- "darkolivegreen3"
sansGambiaCol[sansGambiaCol == 4] <- "darkgoldenrod1"
sansGambiaCol[sansGambiaCol == 5] <- "gray"


######################################################
################# PLOT PCA & BOXPLOT #################
######################################################

svg("cp_kh.svg", width = 11, height = 5.5)
par(mfrow = c(1,2), mar = c(5,5.5,4,2))

# plot the pca
pca.plotter(sansGambia, sansGambiaCol, 1, 3)
text(c(-12, -18, -5, 2), c(21, -8, 2, -7), labels = c("CP4", "CP3", "CP2", "CP1"), col = c("darkgoldenrod1", "darkolivegreen3", "brown1", "deepskyblue"))
mtext("A", 2, las = 2, cex = 2, at = 21, line = 3)

# plot the boxplot
boxplot(boots[9,2,], boots[5,2,], boots[6,2,], boots[7,2,], boots[8,2,], axes = FALSE, ylab = "CP2 vs. KH Group\n(Euclidean Distance)", xlab = "KH Groups")
axis(1, at = 1:5, labels = c("KHA", "KH1", "KH2", "KH3", "KH4"))
axis(2, las = 2)
mtext("B", 2, las = 2, cex = 2, at = 5.6, line = 3)

dev.off()












######################################################
################ AN ATTEMPT AT GGPLOT ################
######################################################

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
