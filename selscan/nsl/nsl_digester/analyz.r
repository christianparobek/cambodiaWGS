## To return nsl values that aren't too near a peak already
## 02 September 2016
## Christian Parobek


#################################
#### DEFINE USEFUL FUNCTIONS ####
#################################

mask.checker <- function(vector, value) {
  vector[1] < value & value < vector[2]
} 
# given a vector of two genomic positions surrounding a nsl high, 
# mask.checker checks a new high nsl position to ensure it's not in this range

machine <- function(df, win){
  chrmask <- list()
  for (i in 1:30){
    if (any(unlist(lapply(chrmask, mask.checker, df[i,2])))){
      next
    }
    else {
      chrmask[[length(chrmask) + 1]] <- c(df[i,2] - win, df[i,2] + win, df[i,1], df[i,2], df[i,5])
    }
  }
  return(chrmask)
}
# given a dataframe of chr, positions, and nsl values
# pick out the top values (30, arbitrary) and eliminate 
# any occurring within win distance of a larger value


#################################
###### READ & PROCESS DATA ######
#################################

pv <- read.table("pv_mono_nsl.bed") # read in data
pv$V5 <- abs(pv$V4) # take absolute value of nsl column
pv_ordered <- pv[with(pv, order(V1, -V5)), ] # order by asc chr, then desc nsl
pv_monster <- split(pv_ordered, pv_ordered$V1) # split into list by chr


pf <- read.table("pf_cp2_nsl.bed") # read in data
pf$V5 <- abs(pf$V4) # take absolute value of nsl column
pf_ordered <- pf[with(pf, order(V1, -V5)), ] # order by asc chr, then desc nsl
pf_monster <- split(pf_ordered, pf_ordered$V1) # split into list by chr


#################################
##### SELECT LARGEST POINTS #####
#################################

pv_scary <- lapply(pv_monster, machine, 180000) # choose largest values 180k from others
pv_less_scary <- data.frame(matrix(unlist(pv_scary), ncol = 5, byrow = TRUE)) # unlist
pv_less_scary_ordered <- pv_less_scary[order(pv_less_scary$X5, decreasing = TRUE),] # order
write.table(head(pv_less_scary_ordered, 30), file = "pv_mono_nsl_tophits.txt", quote = FALSE, sep = "\t")

pf_scary <- lapply(pf_monster, machine, 250000) # choose largest values 180k from others
pf_less_scary <- data.frame(matrix(unlist(pf_scary), ncol = 5, byrow = TRUE)) # unlist
pf_less_scary_ordered <- pf_less_scary[order(pf_less_scary$X5, decreasing = TRUE),] # order
write.table(head(pf_less_scary_ordered, 30), file = "pf_cp2_nsl_tophits.txt", quote = FALSE, sep = "\t")

