# For genetic analysis of my WGS Plasmodium population(s)
# Started late 2014
# Updated 06 May 2015


###################################################
################# IMPORTANT NOTE ##################
###################################################

# The multiVCF that we feed this script
# must have "#CHROM" replaced with "CHROM".
# So run this: sed 's/#CHROM/CHROM/' in.vcf > out.vcf

###################################################
################# LOAD LIBARIES ###################
###################################################

library(stringr)


###################################################
################ DEFINE FUNCTIONS #################
###################################################

## Given a multiVCF, calculate Fws
fwsCalc <- function(dataset) {
  
  ## REMOVE FIRST NINE COLUMNS FROM THE MULTIVCFs
  data <- dataset[-c(1:9)]
  
  ## EXTRACT RELEVANT READ DEPTH DATA, FIRST MATCH
  refCT <- as.data.frame(sapply(data, function(x) str_extract(x, ":[0123456789]+,")))
  # The numbers pre-comma are ref counts
  # Convert to data frame on the fly
  refCT <- sapply(refCT, function(x) str_extract(x, "[0123456789]+"))
  # Clean out the extra chars, leaving only numbers
  refCT <- apply(refCT, c(1,2), as.numeric)
  # Convert to a numeric matrix
  
  altCT <- as.data.frame(sapply(data, function(x) str_extract(x, ",[0123456789]+:")))
  # The numbers post-comma are alt counts
  # Convert to data frame on the fly
  altCT <- sapply(altCT, function(x) str_extract(x, "[0123456789]+"))
  # Clean out the extra chars, leaving only numbers
  altCT <- apply(altCT, c(1,2), as.numeric)
  # Convert to a numeric matrix
  
  ## CALCULATE qs, ps, and Hs, THE PROPORTIONS OF EACH ALLELE IN THE POPULATION
  ps <- rowSums(refCT)/(rowSums(refCT)+rowSums(altCT))
  qs <- rowSums(altCT)/(rowSums(refCT)+rowSums(altCT))
  Hs <- mean(2*ps*qs)
  # Calculate Hs for each variant and take the mean of all variants
  
  ## CALCULATE qw, pw, and Hw, THE PROPORTIONS OF EACH ALLELE IN EACH INDIVIDUAL
  totCT <- refCT + altCT
  # Make a matrix of total counts
  pw <- matrix(, nrow = length(data[,1]), ncol = length(names(data)))
  # Set up pw matrix
  qw <- matrix(, nrow = length(data[,1]), ncol = length(names(data)))
  # Set up qw matrix
  Hw <- matrix(, nrow = length(data[,1]), ncol = length(names(data)))
  # Set up Hw matrix
  
  for (i in 1:length(names(data))) {
    for (j in 1:length(data[,1])) {
      
      pw[j,i] <- refCT[j,i]/totCT[j,i] # Calculate pw per individual and per allele
      qw[j,i] <- altCT[j,i]/totCT[j,i] # Calculate qw per individual and per allele
      Hw[j,i] <- 2*pw[j,i]*qw[j,i] # Calculate Hw per individual and per allele
      
    }
  }
  
  Hw <- colMeans(Hw)
  # Take the column means of Hw matrix, to get a single Hw score for each sample
  
  ## CALCULATE Fws
  1 - Hw/Hs
  
}


###################################################
################## READ IN DATA ###################
###################################################

## READ IN THE Pf AND Pv MULTIVCFs
pv <- read.table("our_goods_pv.pass.vcf", comment.char="#", header=TRUE)
pf <- read.table("our_goods_pf.pass.vcf", comment.char="#", header=TRUE)


###################################################
################# CALCULATE FWS ###################
###################################################

## Do the calculations
pv_fws <- fwsCalc(pv)
pf_fws <- fwsCalc(pf)


###################################################
################### PLOT FWS ######################
###################################################

plot(sort(pv_fws),
     pch=19,
     col="grey25",
     xlim=c(0,75),
     ylim=c(0.2,1),
     xlab="",
     ylab="",
     axes=FALSE,
     type="n",
     main="Multiplicity of Infection")
points(sort(pf_fws),
       pch=19,
       col="grey")
points(sort(pv_fws),
       pch=19,
       col="grey25")
axis(1, at=c(0,25,50,75))
axis(2, at=c(0.2,0.6,1.0), las=2)
segments(0, 0.95, 75, 0.95, lty=2, col="red", lwd=2)
legend(35, 0.6,
       legend=c(expression(italic("P. vivax")), expression(italic("P. falciparum")), expression(paste("F"["WS"] == "0.95"))), 
       pch=c(19,19,NA),
       lty=c(NA,NA,2),
       col=c("grey", "grey25","red"),
       box.lwd=0,
       lwd=2,
       cex=1.25)
mtext("Isolates", side=1, line=2)
mtext(expression("F"["WS"]), side=2, line=2.5)

