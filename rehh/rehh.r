## An Rscript
## Started 07 March 2016
## To plot bifurcation diagrams using REHH
## Inspired by APM's bioarxiv


##################################
### IMPORT LIBS & SOURCE FUNCS ###
##################################

library(rehh)

#source("gistfile1.R") # APM's... couldn't get it to work for me
source("mybifur.r") # a modded bifurcation.diagram, w/o x-axis


###################
## THE FOCAL SNP ##
###################
#1511 for pv_mono chr 14, position 797870
#2028 for pv_all chr 14, position 797870
#214 for pf_cp2 chr 11


##################################
############# VIVAX ##############
##################################

hap <- data2haplohh(hap_file="14_all.hap", map_file="14_all.inp")
hap <- data2haplohh(hap_file="14.hap", map_file="14.inp")

svg("bifur_mono.svg", width = 8, height = 4)

par(mfrow = c(2,1), mar = c(0,0,0,2), oma = c(4,0,0,0))
my.bifurcation.diagram(hap,mrk_foc=1511,all_foc=1,nmrk_l=110,nmrk_r=100, main = "", refsize = 0.5, linecol = "#1b9e77")
text(par("usr")[2]*0.999, mean(par("usr")[3:4])*0.88, labels = "Unselected Allele", srt = 270)
my.bifurcation.diagram(hap,mrk_foc=1511,all_foc=2,nmrk_l=110,nmrk_r=100, main = "", limhapcount = 2, refsize = 0.5, linecol = "#d95f02")
text(par("usr")[2]*0.999, mean(par("usr")[3:4])*0.88, labels = "Selected Allele", srt = 270)

## Plot the axis
min <- round(par("usr")[1]/1000 + 5) # get x min
max <- round(par("usr")[2]/1000 - 5) # get x max
middle <- round(797870/1000)
axis(1, outer = TRUE, at = c(min*1000, 797870, max*1000), labels = c(min, middle, max))

## Add text to plot
mtext("Chromosome 14 Position (kb)", side = 1, line = 2.5)

dev.off()



##################################
########### FALCIPARUM ###########
##################################

hap <- data2haplohh(hap_file="11_cp2.hap", map_file="11_cp2.inp")

par(mfrow = c(2,1), mar = c(0,0,0,2), oma = c(4,0,0,0))
my.bifurcation.diagram(hap,mrk_foc=214,all_foc=1,nmrk_l=10,nmrk_r=10, main = "", refsize = 0.4, linecol = "#1b9e77")
my.bifurcation.diagram(hap,mrk_foc=214,all_foc=2,nmrk_l=10,nmrk_r=10, main = "", limhapcount = 1, refsize = 0.4, linecol = "#d95f02")

