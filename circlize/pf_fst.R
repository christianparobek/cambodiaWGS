## Plotting a circos plot for pairwise FST between CP groups
## Christian Parobek
## Started 07 October 2015

##########
### INITIALIZE LIBRARIES
##########
library(PopGenome)
library(circlize)

##########
### READ IN DATA
##########

## read in genome data
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/PopGenome/pf_whole_genome/")
pf <- readData("vcf/", format="VCF", gffpath = "gff/")

## read in groupings
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/dadi/data/pf/cp_groups/")
cp1 <- read.table("cp1.txt")
cp2 <- read.table("cp2.txt")
cp3 <- read.table("cp3.txt")
cp4 <- read.table("cp4.txt")
pfPops <- set.populations(pf, list(cp1$V1, cp2$V1, cp3$V1, cp4$V1))

## read chromosome info for circlize
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/circlize/")
df <- read.table("pf_chr_data.txt", header=TRUE)


##########
### Calculate pairwise Fst
##########

pfPops_slide <- sliding.window.transform(pfPops, width = 100, jump = 1000, type = 2, whole.data = FALSE)
pfPops_slide <- F_ST.stats(pfPops_slide)

# combine into dataframe with chromosome info
cbind(pfPops_slide@region.names, pfPops_slide@nuc.F_ST.pairwise[1,])

##########
### DEFINE FUNCTIONS
##########
Fst = function(region, value, ...){
  i = getI(...)
  circos.genomicLines(region, value, col = i, ..., lwd = 2)
} ## a function to draw lines in a circos.genomicTrackPlotRegion

##########
### SET PARAMS
##########
circos.par(start.degree=90)
circos.par(gap.degree = 2)

##########
### PLOT DATA
##########
circos.genomicInitialize(df, plotType = "null")
circos.genomicTrackPlotRegion(fst, panel.fun = Fst, ylim = c(0, 1), bg.col = "gray90", bg.border = NA, track.height = 0.05)
circos.genomicTrackPlotRegion(ylim = c(0, 1), bg.col = "gray90", bg.border = NA, track.height = 0.05)
circos.genomicTrackPlotRegion(ylim = c(0, 1), bg.col = "gray90", bg.border = NA, track.height = 0.05)
circos.genomicTrackPlotRegion(ylim = c(0, 1), bg.col = "gray90", bg.border = NA, track.height = 0.05)
circos.genomicTrackPlotRegion(ylim = c(0, 1), bg.col = "gray90", bg.border = NA, track.height = 0.05)
circos.genomicTrackPlotRegion(ylim = c(0, 1), bg.col = "gray90", bg.border = NA, track.height = 0.05)

circos.clear()