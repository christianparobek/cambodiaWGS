# For genetic analysis of my WGS Plasmodium population(s)
# Started late 2014
# Updated 04 May 2015


###################################################
################# Load Libraries ##################
###################################################

library(PopGenome)
library(stringr)


#####################################################
################# Define Functions ##################
#####################################################

## Function to assign pop based on name
pop.namer <- function(GENOME){
  indivs <- get.individuals(GENOME)[[1]]
  pops <- set.populations(pf, list(
    indivs[grep("BB", indivs)], 
    indivs[grep("KP", indivs)], 
    indivs[grep("OM|SN|TB", indivs)]))
  return(pops)
}

## Break a GENOME object by gene and name it
get.pf.genes <- function(GENOME){
  genes <- splitting.data(GENOME, subsites="gene", whole.data = FALSE)
  ## Specify whole.data = FALSE because I don't want to concatenate regions
  ## Extracting gene IDs from the gff INFO field
  ## Need to use stringr's str_extract function to actually get the gene ID out of the big long string of INFO it returns
  id_str <- "PF3D7_[0-9, a-zA-Z, \\., -]+"
  chr01.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr01.gff", chr="Pf3D7_01_v3", extract.gene.names=TRUE), id_str)
  chr02.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr02.gff", chr="Pf3D7_02_v3", extract.gene.names=TRUE), id_str)
  chr03.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr03.gff", chr="Pf3D7_03_v3", extract.gene.names=TRUE), id_str)
  chr04.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr04.gff", chr="Pf3D7_04_v3", extract.gene.names=TRUE), id_str)
  chr05.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr05.gff", chr="Pf3D7_05_v3", extract.gene.names=TRUE), id_str)
  chr06.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr06.gff", chr="Pf3D7_06_v3", extract.gene.names=TRUE), id_str)
  chr07.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr07.gff", chr="Pf3D7_07_v3", extract.gene.names=TRUE), id_str)
  chr08.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr08.gff", chr="Pf3D7_08_v3", extract.gene.names=TRUE), id_str)
  chr09.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr09.gff", chr="Pf3D7_09_v3", extract.gene.names=TRUE), id_str)
  chr10.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr10.gff", chr="Pf3D7_10_v3", extract.gene.names=TRUE), id_str)
  chr11.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr11.gff", chr="Pf3D7_11_v3", extract.gene.names=TRUE), id_str)
  chr12.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr12.gff", chr="Pf3D7_12_v3", extract.gene.names=TRUE), id_str)
  chr13.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr13.gff", chr="Pf3D7_13_v3", extract.gene.names=TRUE), id_str)
  chr14.ids <- str_extract(get_gff_info(gff.file="pf_whole_genome/gff/chr14.gff", chr="Pf3D7_14_v3", extract.gene.names=TRUE), id_str)
  ## Make one big list out of all the little lists of gene IDs
  ## Is there a more elegant way to do all this?
  chr.all.ids <- unlist(list(chr01.ids, chr02.ids, chr03.ids, chr04.ids, chr05.ids, chr06.ids, chr07.ids, chr08.ids, chr09.ids, chr10.ids, chr11.ids, chr12.ids, chr13.ids, chr14.ids))
  genes@region.names <- chr.all.ids
  return(genes)
}

## Break a GENOME object by gene and name it
get.pv.genes <- function(GENOME){
  genes <- splitting.data(GENOME, subsites="gene", whole.data = FALSE)
  ## Specify whole.data = FALSE because I don't want to concatenate regions
  ## Extracting gene IDs from the gff INFO field
  ## Need to use stringr's str_extract function to actually get the gene ID out of the big long string of INFO it returns
  id_str <- "PVX_[0-9, a-zA-Z, \\., -]+"
  chr01.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr01.gff", chr="Pv_Sal1_chr01", extract.gene.names=TRUE), id_str)
  chr02.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr02.gff", chr="Pv_Sal1_chr02", extract.gene.names=TRUE), id_str)
  chr03.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr03.gff", chr="Pv_Sal1_chr03", extract.gene.names=TRUE), id_str)
  chr04.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr04.gff", chr="Pv_Sal1_chr04", extract.gene.names=TRUE), id_str)
  chr05.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr05.gff", chr="Pv_Sal1_chr05", extract.gene.names=TRUE), id_str)
  chr06.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr06.gff", chr="Pv_Sal1_chr06", extract.gene.names=TRUE), id_str)
  chr07.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr07.gff", chr="Pv_Sal1_chr07", extract.gene.names=TRUE), id_str)
  chr08.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr08.gff", chr="Pv_Sal1_chr08", extract.gene.names=TRUE), id_str)
  chr09.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr09.gff", chr="Pv_Sal1_chr09", extract.gene.names=TRUE), id_str)
  chr10.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr10.gff", chr="Pv_Sal1_chr10", extract.gene.names=TRUE), id_str)
  chr11.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr11.gff", chr="Pv_Sal1_chr11", extract.gene.names=TRUE), id_str)
  chr12.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr12.gff", chr="Pv_Sal1_chr12", extract.gene.names=TRUE), id_str)
  chr13.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr13.gff", chr="Pv_Sal1_chr13", extract.gene.names=TRUE), id_str)
  chr14.ids <- str_extract(get_gff_info(gff.file="pv_whole_genome/gff/chr14.gff", chr="Pv_Sal1_chr14", extract.gene.names=TRUE), id_str)
  ## Make one big list out of all the little lists of gene IDs
  ## Is there a more elegant way to do all this?
  chr.all.ids <- unlist(list(chr01.ids, chr02.ids, chr03.ids, chr04.ids, chr05.ids, chr06.ids, chr07.ids, chr08.ids, chr09.ids, chr10.ids, chr11.ids, chr12.ids, chr13.ids, chr14.ids))
  genes@region.names <- chr.all.ids
  return(genes)
}


#####################################################
################### Import Data #####################
#####################################################

## Read in the VCFs and GFFs
pf <- readData("pf_whole_genome/vcf/", format="VCF", gffpath = "pf_whole_genome/gff/")
pv <- readData("pv_whole_genome/vcf/", format="VCF", gffpath = "pv_whole_genome/gff/")

## Read in the ortholog key
ortho_key <- read.table("pv_pf_orthologs.txt")


#####################################################
################### Process Data ####################
#####################################################

## Split into Genes
pf_genes <- get.pf.genes(pf)
pv_genes <- get.pv.genes(pv)


#####################################################
################ Neutrality Stats ###################
#####################################################

## Calculate neutrality Stats
pf_genes <- neutrality.stats(pf_genes)
pv_genes <- neutrality.stats(pv_genes)

## Compare neutrality stats
pf_taj <- cbind(pf_genes@region.names, pf_genes@Tajima.D)
pv_taj <- cbind(pv_genes@region.names, pv_genes@Tajima.D)

## Clean up the ortholist
ortho_key <- ortho_key[!(duplicated(ortho_key[,2]) | duplicated(ortho_key[,2], fromLast = TRUE)),] # remove duplicate Pf entries # need the OR logic because duplicated() only marks second+ occurrence of duplicate
ortho_key <- ortho_key[!(duplicated(ortho_key[,1]) | duplicated(ortho_key[,1], fromLast = TRUE)),] # remove duplicate Pv entries
ortho_key <- ortho_key[ortho_key[,1] %in% pv_taj[,1],] # remove key-value pairs that are not on Pf chromosomes
ortho_key <- ortho_key[ortho_key[,2] %in% pf_taj[,1],] # remove key-value pairs that are not on Pv chromosomes

## Only keep values that are in both ortho_key and the TajD stats
pf_taj <- pf_taj[pf_taj[,1] %in% ortho_key[,2],]
pv_taj <- pv_taj[pv_taj[,1] %in% ortho_key[,1],]

## Order and combine the Pv and Pf TajD values
pvpf <- cbind(
  data.frame(pv_taj)[match(ortho_key[,1], pv_taj[,1]),],
  data.frame(pf_taj)[match(ortho_key[,2], pf_taj[,1]),]
)
names(pvpf) <- c("PvID", "PvTajD", "PfID", "PfTajD") #name it

## Compare Pv to Pf values
pvpf$Diff <- (as.numeric(as.character(pvpf$PvTajD)) - as.numeric(as.character(pvpf$PfTajD))) # have to convert factors -> chars -> numeric
pvpf_ord <- pvpf[with(pvpf, order(Diff)), ]
pvpf_ord_omit <- na.omit(pvpf_ord)
pvpf_ord_omit$Counter <- 1:nrow(pvpf_ord_omit)


################################
## Plot the differences graph ##
################################
par(xpd=NA, mgp = c(3, 1, 0))
plot(pvpf_ord_omit$Diff ~ pvpf_ord_omit$Counter, 
     axes=FALSE, 
     xlab="Orthologous Gene Pairs", 
     ylab=expression(paste(italic("P. vivax"), bold(" D") - italic("P. falciparum"), bold(" D"))),
     ylim=c(-5,4),
     pch=19,
     col="grey25",
     main=expression(bold(paste(Delta, " Tajima's D"))))
# Add axes
axis(1, at=c(0,1800), line=2.5)
axis(2, at=c(-4, 0, 4), line=0.5, las=2)
## Color the highest and lowest points
points(pvpf_ord_omit$Diff[pvpf_ord_omit$Diff > 2] ~ pvpf_ord_omit$Counter[pvpf_ord_omit$Diff > 2], pch=19, col="red")
points(pvpf_ord_omit$Diff[pvpf_ord_omit$Diff < -4] ~ pvpf_ord_omit$Counter[pvpf_ord_omit$Diff < -4], pch=19, col="red")
# Add the polygons at the bottom
polygon(c(0,1800,0), c(-5,-5,-6.5), col="grey30", border="grey30")
polygon(c(0,1800,1800), c(-6.7,-5.2,-6.7), col="grey", border="grey")
# Add text to the polygons at the bottom
text(460, -5.45, expression(paste(italic("P. falciparum"), " Tajima's D")), col="white", cex=0.75)
text(1425, -6.4, expression(paste(italic("P. vivax"), " Tajima's D")), col="black", cex=0.75)

###############################
## Plot the comparison graph ##
###############################
par(mgp = c(2, 1, 0))
plot(density(pv_genes@Tajima.D[,1], na.rm=TRUE), 
     xlim=c(-3,3), 
     lty=2, 
     lwd=2, 
     xlab="Tajima's D",
     main="Genome-wide Tajima's D",
     ylab="Kernel Density", 
     axes=FALSE)
lines(density(pf_genes@Tajima.D[,1], na.rm=TRUE), lwd=2)
legend(0.3, 0.5, 
       legend=c(expression(italic("P. vivax")), expression(italic("P. falciparum"))),
       lty=c(2,1),
       lwd=3,
       box.lwd=0,
       cex=1.25)
axis(1, at=c(-3,-2,-1,0,1,2,3), labels=c("-3","","-1","","1","","3"))
axis(2, at=c(0,0.8), labels=c("",""), las=2)

