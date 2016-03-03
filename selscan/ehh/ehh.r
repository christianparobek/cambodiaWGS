## To graph EHH decay plots 
## Incorporated into ehhRunner.sh
## Ultimately called by my snanlysis file


################################
######## PARSE CMD LINE ########
################################

args <- commandArgs(TRUE)



################################
######### READ IN DATA #########
################################

data <- read.table(paste(args[1], ".ehh.", args[2], ".out", sep = ""), header = FALSE)
#data <- read.table("pv_chr14_797870.ehh.locus1511.out", header = FALSE)

names(data) <- c("physicalPos", "geneticPos", "EHH1", "EHH0", "fullSampEHH")
data$counter <- rev(1:nrow(data))

start <- round((797870-data$physicalPos[1])/1000)
middle <- round(797870/1000)
end <- round((797870+data$physicalPos[nrow(data)])/1000)


svg(paste(args[1], ".ehh.", args[2], ".svg", sep = ""), width = 7, height = 5)
plot(data$counter, data$EHH1, type = "n", ylab = "EHH", xlab = "Chromosome 14 Position (Kb)", axes = FALSE)
lines(data$counter, data$EHH0, col = "grey20", lwd = 3)
lines(data$counter, data$EHH1, col = "brown", lwd = 3)
axis(1, 
     at = c(0, nrow(data)-match(0, data$physicalPos), nrow(data)),
     labels = c(start, middle, end))
axis(2, las = 2)
dev.off()
