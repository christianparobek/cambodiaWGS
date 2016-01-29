## Plot Fws as a bargraph
## Because the main Fws plot is difficult to understand quickly
## Started October 23

num_pv_mono <- 23
num_pv_multi <- 47
num_pf_mono <- 55
num_pf_multi <- 20

fws <- matrix(c(num_pv_multi, num_pv_mono, num_pf_multi, num_pf_mono),ncol=2,byrow=TRUE)

rownames <- c("P. vivax", "P. falciparum")
colnames <- c("Multiclonal", "Monoconal")


############ I HAVE SOMETHING SWITCHED IN HERE!!!!!!!!!!!!!!! PV AND PF ARE BACKWARD> MUST FIX IF USE IN PUB

barplot(fws, 
        beside = TRUE, 
        names.arg = c(expression(italic("P. falciparum")), expression(italic("P. vivax"))),
        ylab = "Number of clinical isolates")

legend(2.45, 40, legend = c("Monoclonal", "Multiclonal"), col = c("gray20", "gray"), pch = 15, pt.cex = 2, cex = 1.3, bty = "n")

