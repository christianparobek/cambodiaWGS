# For genetic analysis of my WGS Pv population(s)
# Started 10 Dec 2014
# Basics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf
# Genomics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
# Extra Commands: http://www.inside-r.org/packages/cran/adegenet/docs/.rmspaces

library(adegenet)

# Read in structure-formatted SNP data
data <- read.table("81k.str", skip=1)
	# skip header line because it throws an error

# Sort the data frame by population
sorted <- data[order(data[,1]),]

# Prepare the data.frame for conversion to genlight object
inds <- sorted$V1 # grab the individual names
pops <- sorted$V2 # grab the population names
sorted <- sorted[-c(1,2)]

# Convert data.frame into 
x <- new("genlight", sorted)

indNames(x) <- inds
ploidy(x) <- 1

glPlot(x, col = funky(5))
dev.off()


pca1 <- glPca(x)

scatter(pca1,posi="bottomright")
title("PCA of Cambodian P. vivax genomes\n axes 1-2")
dev.off()


