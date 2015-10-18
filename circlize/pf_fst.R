## Plotting a circos plot for pairwise FST between CP groups
## Christian Parobek
## Started 07 October 2015

##########
### INITIALIZE LIBRARIES
##########
library(PopGenome)
library(circlize)
library(stringr)


##########
### DEFINE FUNCTIONS
##########

bedMaker <- function(index){
  df <- as.data.frame(cbind(pfPops_slide@region.names, pfPops_slide@nuc.F_ST.pairwise[index,])) ## make df w/chr info
  c1 <- str_extract(df$V1, "chr\\d\\d.vcf") ## need to clean it up to get it into bed format for circlize
  c2 <- str_replace(c1, ".vcf", "_v3")
  chromosome <- str_replace(c2, "chr", "Pf3D7_")
  start <- as.integer(str_extract(str_extract(df$V1, " \\d*"), "[0-9]+")) # match start value and trim whitespace
  end <- as.integer(str_extract(str_extract(df$V1, "- \\d*"), "[0-9]+")) # match end value and trim dash/whitespace
  fst_bed <- cbind.data.frame(chromosome, start, end) # make the data take bed file format
  fst_bed$fst <- as.numeric(as.character(df$V2)) # bet there's a better way to do this, but couldn't find it
  fst_bed <- fst_bed[!end < start,] ## clean entries where end < start, probably due to a bug in PopGenome
  return(fst_bed)
} ## make a bed file from the FST data to draw lines

fstDrawer <- function(region, value, ...){
  i = getI(...)
  circos.genomicPoints(region, value, col = i, ..., pch = 1, cex = 0.01, lwd = 0.01)
} ## a function to draw lines in a circos.genomicTrackPlotRegion

chrNamer <- function(region, value, ...){
  circos.genomicText(region, value, labels.column =1, facing = "clockwise")
} ## a function to draw lines in a circos.genomicTrackPlotRegion

trackPlotter <- function(bed, panel_fun){
  circos.genomicTrackPlotRegion(bed, panel.fun = panel_fun, ylim = c(0, 1), bg.col = c(rep("gray90", 14), "white"), bg.border = NA, track.height = 0.06)
} ## integrates the fstDrawer function to draw fst tracks

axisAdder <- function(){
  circos.genomicTrackPlotRegion(ylim = c(0, 1), bg.border = NA, 
                                track.height=0.05,
                                panel.fun = function(region, value, ...) {
                                  circos.axis(h = 0, major.at = NULL, labels = NULL, 
                                              labels.cex = 0.3 * par("cex"), labels.facing = "clockwise", 
                                              major.tick.percentage = 0.2)})} ## add axis without text

chrNamer <- function(){
  circos.genomicTrackPlotRegion(ylim=c(0,1),
                                bg.col=NULL, bg.border=NA, track.height=0.05,
                                panel.fun = function(region, value, ...){
                                  circos.text(mean(get.cell.meta.data("xlim")), 
                                              0.5, labels = get.cell.meta.data("sector.numeric.index"))})
} ## add chromosome names wrt sector index

leadingZeroScrubber <- function(bed){
  first_zero_index <- match(TRUE, bed$fst > 0)
  bed[-(1:first_zero_index),]
} ## For some reason there are leading zeros that get plotted outside the plot area


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
pfPops <- set.populations(pf, list(as.character(cp1$V1), as.character(cp2$V1), as.character(cp3$V1), as.character(cp4$V1)))

## read chromosome info for circlize
setwd("/run/user/1000/gvfs/sftp:host=kure.unc.edu,user=prchrist/proj/julianog/users/ChristianP/cambodiaWGS/circlize/")
chr_info <- read.table("pf_chr_data.txt", header=TRUE)
fsttest <- read.table("pf_fst_data.txt", header=TRUE)

##########
### Calculate pairwise Fst
##########

pfPops_slide <- sliding.window.transform(pfPops, width = 10000, jump = 1000, type = 2, whole.data = FALSE)
pfPops_slide <- F_ST.stats(pfPops_slide)

# make beds from the data
cp1cp2_bed <- leadingZeroScrubber(bedMaker(1))
cp1cp3_bed <- leadingZeroScrubber(bedMaker(2))
cp1cp4_bed <- leadingZeroScrubber(bedMaker(3))
cp2cp3_bed <- leadingZeroScrubber(bedMaker(4))
cp2cp4_bed <- leadingZeroScrubber(bedMaker(5))
cp3cp4_bed <- leadingZeroScrubber(bedMaker(6))



##########
### SET PLOT PARAMS AND PLOT DATA
##########

#svg("circlize_fst.svg", width = 7, height = 7)
png("circlize_fst.png", width = 7, height = 7, unit = "in", res = 300)

## set params
circos.par(start.degree = 90, gap.degree = 2, unit.circle.segments = 100)

## plot data
circos.genomicInitialize(chr_info, plotType = NULL) 
axisAdder()
trackPlotter(cp1cp2_bed, fstDrawer)
trackPlotter(cp1cp3_bed, fstDrawer)
trackPlotter(cp1cp4_bed, fstDrawer)
trackPlotter(cp2cp3_bed, fstDrawer)
trackPlotter(cp2cp4_bed, fstDrawer)
trackPlotter(cp3cp4_bed, fstDrawer)
chrNamer()

## add the track names
circos.updatePlotRegion(sector.index = "Pf_names", track.index = 1, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")),0.5, labels = "")

circos.updatePlotRegion(sector.index = "Pf_names", track.index = 2, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")),0.5, labels = "CP1-CP2", cex = 0.75)

circos.updatePlotRegion(sector.index = "Pf_names", track.index = 3, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")),0.5, labels = "CP1-CP3", cex = 0.68)

circos.updatePlotRegion(sector.index = "Pf_names", track.index = 4, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")),0.5, labels = "CP1-CP4", cex = 0.61)

circos.updatePlotRegion(sector.index = "Pf_names", track.index = 5, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")),0.5, labels = "CP2-CP3", cex = 0.54)

circos.updatePlotRegion(sector.index = "Pf_names", track.index = 6, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")),0.5, labels = "CP2-CP4", cex = 0.47)

circos.updatePlotRegion(sector.index = "Pf_names", track.index = 7, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")),0.5, labels = "CP3-CP4", cex = 0.39)

circos.updatePlotRegion(sector.index = "Pf_names", track.index = 8, bg.border = NA)
circos.text(mean(get.cell.meta.data("xlim")),0.5, labels = "")

circos.clear()
dev.off()

