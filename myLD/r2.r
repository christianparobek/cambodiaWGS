## For determining LD decay between biallelic SNPs
## Started 4 November 2014
## Christian Parobek

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

haplo.counter <- function(datamatrix, snp1index, snp2index) {
  haplotype_matches <- as.numeric(apply(mat[c(snp1index,snp2index),], 2, FUN = identical, as.numeric(c(0,0))))
  number_of_matches <- sum(haplotype_matches)
  match_frac <- number_of_matches / ncol(mat) # number of haplotype matches / number of possible matches
  return(match_frac)
} # given SNP (row) indices and hap identities, get hap freqs

ref.freq <- function(snpindex) {
  1 - sum(mat[snpindex,])/ncol(mat)
} # get ref freq at a given SNP index

calc.r2 <- function(datamatrix, snp1index, snp2index) {
  refref <- haplo.counter(datamatrix, snp1index, snp2index) # get refref haplotype fraction for two SNPs
  snp1ref <- ref.freq(snp1index) # get ref frequency at a given index
  snp2ref <- ref.freq(snp2index) # get ref frequency at a given index
  D <- refref - snp1ref * snp2ref # calculate D
  r2 <- (D*D)/((snp1ref)*(1-snp1ref)*(snp2ref)*(1-snp2ref))
  return(r2)
}


###################################################
############ READ IN AND PREPARE DATA #############
###################################################

## READ IN THE Pf AND Pv MULTIVCFs
pf <- read.table("our_goods_pf.pass.vcf", comment.char="#", header=TRUE)

## SUBSET THE DATAFRAME INTO METADATA AND DATA
metadata <- pf[c(1:9)] # subset the metadata
data <- sapply(pf[c(10:ncol(pf))], function(x) str_extract(x, "[0123456789]")) # subset and get genotypes only
class(data) <- "numeric" # make it numeric matrix rather than character matrix

## SPLIT METADATA AND DATA BY CHROMOSOME
metadata_chrs <- split(metadata, pf$CHROM)
data_chrs <- split(data, pf$CHROM)


###################################################
################## READ IN DATA ###################
###################################################

## WANT TO CALCULATE PAIRWISE DISTANCES BY CHROMOSOME ONCE
x <- by(metadata$POS, metadata$CHROM, dist) # calc dist matrix per chromosome
y <- as.matrix(x$Pf3D7_01_v3)

z <-  # split data by chromosome

y <- as.matrix(x)
y[lower.tri(y)] <- NA # dist -> matrix fills mirrors lower half. Need to convert half to NAs to ignore
indices <- which(y>0 & y<1000, arr.ind=TRUE)



## for any two snps, four possible haplotypes:
## 00, 01, 10, 11
## in VCF language, 0 is for ref allele, 1 is for alt allele








