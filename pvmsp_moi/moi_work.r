## 19 May 2016
## To determine MOI from pvmsp1 data


#####################################
###### PROCESS CLUSTERING DATA ######
#####################################

data <- read.table("DeenPvmsp1MOI_sampInfo.tab.csv", header = TRUE, na.strings = "NA")
  # read in data

complete <- data[complete.cases(data),]
  # remove missing data

cov <- complete[!(complete$R1.totalReadCnt < 1000 | complete$R2.totalReadCnt < 1000),]
  # remove ones with low coverage (<1000 in either run)

cutoff_0.005 <- cov[cov$c_AveragedFrac > 0.005,]
  # make sure all values have at least 0.005 average


#####################################
###### INTERSECT WITH WGS DATA ######
#####################################

our_wgs <- read.table("wgs_om")
  # read in list of samples we sequenced

sum(our_wgs$V1 %in% data$s_Sample)
  # number of wgs samples with msp data
sum(our_wgs$V1 %in% complete$s_Sample)
  # number of wgs samples with two msp reps
sum(our_wgs$V1 %in% cov$s_Sample)
  # numer of wgs samples w/ two reps w/ 1000x cov

subset <- cutoff_0.005[cutoff_0.005$s_Sample %in% our_wgs$V1,]
  # subset the pvmsp1 data just down to the wgs samples

counts <- rle(as.vector(subset$s_Sample))
  # collapse / count occurrences 

summary(counts$lengths)
  # summarize data

plot(1:length(counts$lengths), counts$lengths)
  # plot distribution of moi
