# For genetic analysis of my WGS Plasmodium population(s)
# Started 10 Dec 2014
# Updated 17 March 2016
# Basics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf
# Genomics Tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
# Extra Commands: http://www.inside-r.org/packages/cran/adegenet/docs/.rmspaces


###################################################
################# Load Libraries ##################
###################################################

library(adegenet)
library(stringr)
library(pegas)
library(ggplot2)
library(gridExtra)

#####################################################
################# Define Functions ##################
#####################################################

## Function to create genlight from VCF.
## The AAKM contigs can screw with this
## so `grep -v "=AAKM" infile.vcf > outfile.vcf`
genlight.maker <- function(infile) {
  loci <- read.vcf(infile, from = 1, to = 1000000) # reads first million sites
  genlight <- new("genlight", loci) # convert data frame into genlight object
  ploidy(genlight) <- as.integer(1) # add back population information
  return(genlight)
}

## Function to assign samples to pops
## based on a list of their names
pop.definer <- function(ind_names) {
  library(stringr)
  kp <- as.numeric(str_detect(ind_names, "BB"))*1 # assign KP pop number
  bb <- as.numeric(str_detect(ind_names, "KP"))*2 # assign BB pop number
  om <- as.numeric(str_detect(ind_names, "OM"))*3 # assign OM pop number
  sn <- as.numeric(str_detect(ind_names, "SN"))*3 # assign SN pop number
  tb <- as.numeric(str_detect(ind_names, "TB"))*3 # assign TB pop number
  srr <- as.numeric(str_detect(ind_names, "SRR"))*4 # assign SRR pop number
  err <- as.numeric(str_detect(ind_names, "ERR"))*4 # assign ERR pop number
  pops <- kp + bb + om + sn + tb + srr + err
  return(pops)
}

## Function to plot eigenplots, ggstyle
eig.plotter.gg <- function(pca, plotname, ylim){
  
  for_bar <- as.data.frame(pca$eig)
  # convert to a data frame
  for_bar$index <- 1:nrow(for_bar)
  # add a counter column
  names(for_bar) <- c("value", "index")
  # make names reasonable
  for_bar$value <- for_bar$value/sum(for_bar$value)*100
  # convert values to fraction of whole
  
  a_plot <- ggplot(for_bar, aes(x = index, y = value)) + geom_bar(stat = "identity") + params + 
    labs(
      x = "Principal component",
      y = "Percent of variance explained") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    ggtitle(plotname) + scale_y_continuous(limits = c(0, ylim))
  
  return(a_plot)
}

## Function to plot PCAs - color
ggcolor <- function(pca, plotname, pops, PCx, PCy) {
  pic <- ggplot(as.data.frame(pca$scores), aes_string(x=colnames(pca$scores)[PCx], y=colnames(pca$scores)[PCy])) + 
    geom_point(alpha=0.8, size = 4, color = factor(pops)) + 
    params +
    ggtitle(plotname) +
    labs(
      x = paste("PC", PCx, " - ", round(pca$eig[PCx]/sum(pca$eig)*100), "% of total variance", sep = ""),
      y = paste("PC", PCy, " - ", round(pca$eig[PCy]/sum(pca$eig)*100), "% of total variance", sep = ""))
  return(pic)
}


##############################################################
################# READ IN DATA AND CALC PCAS #################
##############################################################

## READ ARGS FROM COMMAND LINE
args = commandArgs(trailingOnly=TRUE)

## read in PV data
pv_cam_gl <- genlight.maker(args[1]) # make genlight; vcf must not contain extra contigs
#pv_cam_gl <- genlight.maker("our_goods_pv.pass.sansAAKM.vcf")
pv_cam_pops <- pop.definer(indNames(pv_cam_gl)) # define pops
pv_cam_pca <- glPca(pv_cam_gl, nf = 4) # calculate PCA

## read in PF data
pf_cam_gl <- genlight.maker(args[2]) # make genlight
#pf_cam_gl <- genlight.maker("our_goods_pf.pass.vcf")
pf_cam_pops <- pop.definer(indNames(pf_cam_gl)) # define pops
pf_cam_pca <- glPca(pf_cam_gl, nf = 4) # calculate PCA


##############################################################
################# DEFINE SUBGROUP MEMBERSHIP #################
##############################################################

######### FALCIPARUM #########

## define CP2 pops
cp_grps <- read.table("cp_groups")
cp_grps$V2[cp_grps$V2 == 3] <- 1 # make CP1,3,4 all same value
cp_grps$V2[cp_grps$V2 == 4] <- 1# make CP1,3,4 all same value

## fix SN084 name
row.names(pf_cam_pca$scores)[row.names(pf_cam_pca$scores) == "SN084"] <- "SN094" 

## match order between key and pca scores
pf_ord_names <- row.names(pf_cam_pca$scores)
pf_grp_names <- as.character(cp_grps$V1)
pf_key <- match(pf_ord_names, pf_grp_names)
cp_grp_membership <- cp_grps$V2[pf_key]

######### VIVAX #########

mono_grps <- read.table("monos")
moi_status <- mono_grps$V2
  # monos is already in the same order as pca


##############################################################
##################### SETUP GGPLOT THEME #####################
##############################################################

## Define basic plotting parameters
params <-  theme_bw() + 
  theme(axis.text = element_text(size = 13), 
        axis.title = element_text(size = 13), 
        legend.position = "none")


##############################################################
################# ADD CONSISTENT JITTER TO PF ################
##############################################################

pf_pca_jitter <- pf_cam_pca
pf_pca_jitter$scores <- as.data.frame(jitter(pf_pca_jitter$scores, factor = 700))


##############################################################
########### HIGHLIGHT CP2/MONOS & SHOW EIGENVECTORS ##########
##############################################################

#palette(c("gray50","#1b9e77"))
palette(c("gray70","gray20"))
  # define the palette

######### PLOT SUBSET PCAs ######### 

cp2 <- ggcolor(pf_pca_jitter, expression(italic("P. falciparum")), cp_grp_membership, 1, 2)
monos <- ggcolor(pv_cam_pca, expression(italic("P. vivax")), mono_grps$V2, 1, 2)

svg("subset.svg", width = 7, height = 3.6)
grid.arrange(monos, cp2, ncol=2)
grid.text("A", x = unit(0.03, "npc"), y = unit(0.95, "npc"), gp=gpar(fontsize=20))
grid.text("B", x = unit(0.53, "npc"), y = unit(0.95, "npc"), gp=gpar(fontsize=20))
dev.off()


######### PLOT EIGENVALUES ######### 

pf_eig <- eig.plotter.gg(pf_pca_jitter, expression(italic("P. falciparum")),30)
pv_eig <- eig.plotter.gg(pv_cam_pca, expression(italic("P. vivax")), 2)

svg("eigens.svg", width = 7, height = 3.4)
grid.arrange(pv_eig, pf_eig, ncol=2)
grid.text("A", x = unit(0.03, "npc"), y = unit(0.95, "npc"), gp=gpar(fontsize=20))
grid.text("B", x = unit(0.53, "npc"), y = unit(0.95, "npc"), gp=gpar(fontsize=20))
dev.off()


######### PLOT EXTRA PCAs ######### 

cp2_1_2 <- ggcolor(pf_pca_jitter, expression(italic("P. falciparum")), cp_grp_membership, 1, 2)
cp2_1_3 <- ggcolor(pf_pca_jitter, expression(italic("P. falciparum")), cp_grp_membership, 1, 3)
cp2_2_3 <- ggcolor(pf_pca_jitter, expression(italic("P. falciparum")), cp_grp_membership, 2, 3)

mono_1_2 <- ggcolor(pv_cam_pca, expression(italic("P. vivax")), mono_grps$V2, 1, 2)
mono_1_3 <- ggcolor(pv_cam_pca, expression(italic("P. vivax")), mono_grps$V2, 1, 3)
mono_2_3 <- ggcolor(pv_cam_pca, expression(italic("P. vivax")), mono_grps$V2, 2, 3)

svg("extra_pcs.svg", width = 6, height = 9)
grid.arrange(mono_1_2, cp2_1_2, mono_1_3, cp2_1_3, mono_2_3, cp2_2_3, ncol=2)
dev.off()

##############################################################
################ ggPlot EXTRA FALCIPARUM PCs #################
##############################################################

## Identify clusters using K-means
pf_grp <- find.clusters(pf_cam_gl, n.pca = 3, max.n.clust = 10, choose.n.clust = FALSE, criterion = "diffNgroup")
  # using n.pca = 3 because these three PCs describe most of the variance
pf_dapc <- dapc(pf_cam_gl, pf_grp$grp, n.pca = 3, n.da = 20)

## Make plots of the extra PCs for Pf
pf_1_2 <- ggcolor(pf_cam_pca, "P. falciparum", 700, pf_dapc$grp, 1, 2)
pf_1_3 <- ggcolor(pf_cam_pca, "P. falciparum", 700, pf_dapc$grp, 1, 3)
pf_2_3 <- ggcolor(pf_cam_pca, "P. falciparum", 700, pf_dapc$grp, 2, 3)
eigens <- ggplot(as.data.frame(pf_cam_pca$eig), aes(x=1:length(pf_cam_pca$eig), y = pf_cam_pca$eig)) + 
  geom_bar(stat = "identity") + 
  params + 
  ggtitle(bquote(atop("Principal Component Eigenvalues"))) +
  labs(x = "", y = "Variance")

svg("extra_pcs.svg", width = 8, height = 8.5)
grid.arrange(eigens, pf_1_2, pf_1_3, pf_2_3, ncol = 2)
grid.text("A", x = unit(0.03, "npc"), y = unit(0.93, "npc"), gp=gpar(fontsize=27))
grid.text("B", x = unit(0.53, "npc"), y = unit(0.93, "npc"), gp=gpar(fontsize=27))
grid.text("C", x = unit(0.03, "npc"), y = unit(0.43, "npc"), gp=gpar(fontsize=27))
grid.text("D", x = unit(0.53, "npc"), y = unit(0.43, "npc"), gp=gpar(fontsize=27))
dev.off()


## Make plots of K-means -- P. vivax + P. falciparum only
pv_grp <- find.clusters(pv_cam_gl, n.pca = 4, max.n.clust = 10, choose.n.clust = FALSE, criterion = "diffNgroup")
pf_grp <- find.clusters(pf_cam_gl, n.pca = 3, max.n.clust = 10, choose.n.clust = FALSE, criterion = "diffNgroup")

plot(pv_grp$Kstat, type = "l", main = expression(italic("P. vivax")), xlab = "Number of Clusters", ylab = "BIC", las = 1, xaxt = "n")
axis(1, at = 1:10)
points(pv_grp$Kstat, pch = 19)
points(4, pv_grp$Kstat[4], pch = 19, col = "red", cex = 1.2)

## plot it
pf_grp <- find.clusters(pf_cam_gl, n.pca = 3, max.n.clust = 10, choose.n.clust = FALSE, criterion = "diffNgroup")
dapc1 <- dapc(pf_cam_gl, pf_grp$grp, n.pca = 20, n.da = 40)
scatter(dapc1)
dapc1$grp # these are the grouping values

## plot some things
ggcolor(pf_cam_pca, "P. falciparum", 0, pf_cam_pops, 1, 2) # plot by province
ggcolor(pf_cam_pca, "P. falciparum", 0, dapc1$grp, 2, 3) # plot by group


## plot it a little better looking
myCol <- c("darkblue","purple","green","orange")
scatter(dapc1, posi.da="bottomright",  bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")
scatter(dapc1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:4))


svg("pf_pv_num_clust.svg", width = 5, height = 5)

par(mfrow=c(1,2))

plot(pf_grp$Kstat, type = "l", main = expression(italic("P. falciparum")), xlab = "Number of Clusters", ylab = "BIC", las = 1, xaxt = "n")
axis(1, at = 1:10)
points(pf_grp$Kstat, pch = 19)
points(4, pf_grp$Kstat[4], pch = 19, col = "red", cex = 1.2)

plot(pv_grp$Kstat, type = "l", main = expression(italic("P. vivax")), xlab = "Number of Clusters", ylab = "BIC", las = 1, xaxt = "n")
axis(1, at = 1:10)
points(pv_grp$Kstat, pch = 19)
points(4, pv_grp$Kstat[4], pch = 19, col = "red", cex = 1.2)

dev.off()


## Make plots of K-means -- P. falciparum only
svg("pf_num_clust.svg", width = 5, height = 5)
plot(pf_grp$Kstat, type = "l", main = expression(italic("P. falciparum")), xlab = "Number of Clusters", ylab = "BIC", las = 1, xaxt = "n")
axis(1, at = 1:10)
points(pf_grp$Kstat, pch = 19)
points(4, pf_grp$Kstat[4], pch = 19, col = "red", cex = 1.2)
dev.off()
