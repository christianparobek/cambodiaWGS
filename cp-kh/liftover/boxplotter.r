## A script to compare allele frequencies in our groups
## and compare it to Miotto's allele frequencies by group



###################################################
################# LOAD LIBARIES ###################
###################################################

library(stringr)
library(abind)


###################################################
################ DEFINE FUNCTIONS #################
###################################################

# Calculate the allele frequencies by SNP, given a VCF
af.calc <- function(vcf) {
  
  data <- vcf[-c(1:9)] # get rid of first nine rows
  
  altCT <- as.data.frame(sapply(data, function(x) str_extract(x, "[0123456789]:")))
  altCT <- sapply(altCT, function(x) str_extract(x, "[0123456789]+")) # clean out extra chars
  altCT <- apply(altCT, c(1,2), as.numeric) # convert to numeric matrix
  
  vector <- apply(altCT, 1, mean)
  
  return(vector)
} ## 0 is ref, 1 is alt for VCF



###################################################
################## READ IN DATA ###################
###################################################

## READ IN THE Pf AND Pv MULTIVCFs
cp1 <- read.table("vcfs/cp1_intersected.vcf", header=FALSE)
cp2 <- read.table("vcfs/cp2_intersected.vcf", header=FALSE)
cp3 <- read.table("vcfs/cp3_intersected.vcf", header=FALSE)
cp4 <- read.table("vcfs/cp4_intersected.vcf", header=FALSE)
  
miotto_af <- read.table("beds/liftover_intersected.bed", header = FALSE)


###################################################
################## BUILD MATRIX ###################
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
############### GET POINT ESTIMATE ################
###################################################

dist(af_matrix)


###################################################
#################### BOOTSTRAP ####################
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


par(mar = c(5,5.5,4,2))
boxplot(boots[9,2,], boots[5,2,], boots[6,2,], boots[7,2,], boots[8,2,], axes = FALSE, ylab = "CP2 vs. KH Group\n(Euclidean Distance)")
axis(1, at = 1:5, labels = c("KHA", "KH1", "KH2", "KH3", "KH4"))
axis(2)
