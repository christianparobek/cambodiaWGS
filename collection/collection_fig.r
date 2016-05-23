## To integrate Jess's collection date figure
## With a map of cambodia


###############################
######## LOAD LIBRARIES #######
###############################

library(rgdal)
#library(RColorBrewer)
library(ggplot2)


###############################
######### IMPORT DATA #########
###############################

pv <- read.table("PV_collectiontimes.txt", header = TRUE, sep = "\t", )
pf <- read.table("PF_collectiontimes.txt", header = TRUE, sep = "\t")
  # read in data data

country<-readOGR("KHM_adm/KHM_adm0.shp", layer="KHM_adm0")
provinces<-readOGR("KHM_adm/KHM_adm1.shp", layer="KHM_adm1")
  #read in country/province shapefiles


###############################
###### CLEAN DATE DATA ########
###############################

pv[is.na(pv)] <- 0
pf[is.na(pf)] <- 0
  # set all NA to 0


###############################
########## PLOT DATA ##########
###############################

svg("datemap.svg", width = 8, height = 5)

m <- matrix(c(1,3,1,3,2,3,2,4), byrow = TRUE, nrow = 4)
layout(m)

par(mar = c(2,4,1,1))
pv_plot <- barplot(t(pv[,2:4]), las = 2, xaxt = "n", col = c("#fc8d59", "#ffffbf","#91bfdb"))
mtext(expression(paste(italic("P. vivax"), " isolates")), side = 2, line = 3, padj = 1)
  # need to save the barplot to get coordinates for plotting axis
  # padj set to 1 to top align the text, else the tail on "p" screws things up

par(mar = c(4,4,0,1))
pf_plot <- barplot(t(pf[,2:4]), las = 2, xaxt = "n", col = c("#fc8d59", "#ffffbf","#91bfdb"))
text(x=pv_plot-1.0, y=-4.5, pf$Month, xpd = TRUE, srt=45, cex = 0.9)
mtext(expression(paste(italic("P. falciparum"), " isolates")), side = 2, line = 3, padj = 1)


par(mar = c(0,1,1,1))
plot(provinces, col = "grey92", bg = "white", ylim = c(11,14)) # , xlim = c(104.5, 104.6) , bg = "grey"

indices <- grep("Batdâmbâng|Otdar Mean Chey|Kâmpôt|Kâmpóng Sp\u009c|Preah Vihéar|Krong Preah Sihanouk|Krong Pailin", provinces$NAME_1)
polygon(provinces@polygons[[indices[1]]]@Polygons[[1]]@coords, col = "#fc8d59") #steelblue1
  # get BB
polygon(provinces@polygons[[indices[3]]]@Polygons[[3]]@coords, col = "#91bfdb")
  # get KP
polygon(provinces@polygons[[indices[6]]]@Polygons[[1]]@coords, col = "#ffffbf")
  # get OM

par(mar = c(1,1,1,1))
plot(1:10, 1:10, type = "n", axes = FALSE, ylab = "", xlab = "")
legend("center", legend = c("Oddar Meanchey", "Battâmbâng", "Kâmpôt"), ncol = 1, fill = c("#ffffbf","#fc8d59","#91bfdb"), cex = 1.5, bty = "n")

dev.off()


