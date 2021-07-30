library(ggplot2)
seed = 24125

x <- read.table(paste0("hist/", seed, "_data_set1.txt"), header = TRUE, sep = "\t")

pdf(paste0("hist/Mod_pred_", seed, "_set1.pdf"),paper='special')
hist(x[,1])
dev.off()

pdf(paste0("hist/Edges_", seed, "_set1.pdf"),paper='special')
hist(x[,2])
dev.off()

# length(which(x[, 1] < -0.05 & x[,1] > -.1))
# length(which(x[, 1] < 0.1375 & x[,1] > 0.0875))
# length(which(x[, 1] < 0.325 & x[,1] > 0.275))
