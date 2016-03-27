data <- read.table("fws_pvmsp1_data.txt", header = TRUE)

svg("fws_vs_pvmsp1.svg", width = 9, height = 5)
par(mfrow = c(1,2))
# plot Fws vs max percent
plot(data$Fws, data$Max_percent, xlab = "Fws", ylab = "Percentage of largest clone")
# plot Fws vs max MOI
boxplot(Fws ~ moi_pvmsp1_clones_gt_10., data = data, notch = TRUE, col = "gray", xlab = "MOI (num clones > 10%)", ylab = "Fws")
dev.off()