#!/usr/bin/env Rscript
## To turn the selscan-formatted MAP file into REHH format
## Started 07 March 2016


################################################
############### PARSE COMMANDLINE ##############
################################################

## READ ARGS FROM COMMAND LINE
args = commandArgs(trailingOnly=TRUE)


#data <- read.table(file = "../selscan/ihs/pv_mono/01.map")
data <- read.table(file = args[1])
out <- cbind(as.character(data$V2), 1, data$V4, 2, 1)
write.table(out, file = args[2], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)