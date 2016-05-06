##################################
######## IMPORT LIBRARIES ########
##################################

library(stringr)


##################################
####### PARSE COMMAND LINE #######
##################################

args = commandArgs(trailingOnly=TRUE)


##################################
########## IMPORT DATA ###########
##################################

#data <- read.table("mono_res_genes.anno.vcf", header = FALSE, sep = "\t")
#data <- read.table("pv_nsl_genes.anno.vcf", header = FALSE, sep = "\t")
data <- read.table(args[1], header = FALSE, sep = "\t")

##################################
########## GET METADATA ##########
##################################

ids <- str_extract(data[,8], "PVX_\\d*\\||PF3D7_\\d*\\|")
  ## get the pvx id

eff <- str_extract(data[,8], "EFF=.*\\(")
  ## get the effect

aas <- str_extract(data[,8], "\\|[A-Z]\\d*[A-Z]\\|")
  ## get the amino acid changes, where applicable

chr <- data[,1]

chr_pos <- data[,2]
  ## get the chromosomal position of the mutation

ref <- data[,4]

alt <- data[,5]

##################################
######### GET GENOTYPES ##########
##################################

genos_long <- data[-c(1:9)]
  # remove the metadata columns (first nine columns)

genos <- sapply(genos_long, function(x) str_extract(x, "\\d"))
  # make a matrix of genotypes - 1 for mutant and 0 for wt

class(genos) <- "numeric"
  # convert from character to numeric matrix

#freqs <- apply(genos, 1, sum, na.rm = TRUE)/ncol(genos)
  # sum the rows and divide by number of samples
  # this works just fine if only one alt allele per site

## The following works for multiple alleles per site

freqs <- NULL

for (i in 1:nrow(genos)) {
  
  alleles <- 1:max(genos[i,], na.rm = TRUE)
  # select all genotypes > 0... i.e. all mutant genotypes
  
  string <- NULL
  # reset the string that I'll print for the current frequency
  
  for (j in alleles) {
    counts <- str_count(genos[i,], as.character(j))
    frac <- round(sum(counts, na.rm = TRUE)/ncol(genos), digits = 3)
    string <- paste(string, frac, sep = ",")
  }
  
  freqs <- append(freqs, string)
  
}

##################################
######### ASSEMBLE TABLE #########
##################################

table <- cbind(as.character(chr), chr_pos, as.character(ref), as.character(alt), ids, eff, aas, freqs)

write.table(table, file = args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

