## Plot Fws as a bargraph
## Because the main Fws plot is difficult to understand quickly
## Started October 23

num_pv_mono <- 28
num_pv_multi <- 42
num_pf_mono <- 62
num_pf_multi <- 18

fws <- matrix(c(num_pv_multi, num_pf_multi, num_pv_mono, num_pf_mono),ncol=2,byrow=TRUE)

rownames(fws) <- c("Pv", "Pf")
colnames(fws) <- c("Multiclonal", "Monoconal")


####################################
## multipanel with other fws plot ##
####################################

svg("combined_fws.svg", width = 6.75, height = 3)
par(mar = c(3,4,1,0))
layout(matrix(c(1,2), 1, 2), widths = c(2,1.3))

plot(pf_point_estimate,
     pch=19,
     col="grey25",
     xlim=c(0,80),
     ylim=c(0.2,1),
     xlab="",
     ylab="",
     axes=FALSE,
     type="n")
segments(1:length(pf_max), pf_max[pf_order], 
         1:length(pf_min), pf_min[pf_order],
         col="grey", lwd=2)
segments(1:length(pv_max), pv_max[pv_order], 
         1:length(pv_min), pv_min[pv_order],
         col="grey25", lwd=2)
points(sort(pf_point_estimate),
       pch=19,
       col="grey", cex = 0.65)
points(sort(pv_point_estimate),
       pch=19,
       col="grey25", cex = 0.65)
axis(1, at=c(0,25,50,75), cex.axis = 0.8)
axis(2, at=c(0.2,0.6,1.0), las=2, cex.axis = 0.8)
segments(0, 0.95, 75, 0.95, lty=2, col="grey25", lwd=2)
legend(35, 0.6,
       legend=c(expression(italic("P. vivax")), expression(italic("P. falciparum")), expression(paste(italic("F")["WS"] == "0.95"))), 
       pch=c(19,19,NA),
       lty=c(NA,NA,2),
       col=c("grey25", "grey","black"),
       box.lwd=0,
       lwd=2,
       cex=0.9)
mtext("Isolates", side=1, line=1.7, cex = 0.8, at = 38)
mtext(expression(italic("F")["WS"]), side=2, line=2, cex = 0.8)
mtext("A", 2, line = 2.5, las = 2, cex = 1.6, padj = -5)

barplot(fws, 
        beside = TRUE, 
        names.arg = c(expression(italic("Pv")), expression(italic("Pf"))),
        ylab = "", las = 1, cex.names = 0.8, cex.axis = 0.8)
mtext("Number of clinical isolates", side=2, line=2, cex = 0.8)
legend(1, 60, legend = c("Multiclonal", "Monoclonal"), col = c("gray25", "gray"), pch = 15, pt.cex = 1, cex = 0.9, bty = "n")
mtext("B", 2, line = 2.5, las = 2, cex = 1.6, padj = -5)
 
dev.off()
