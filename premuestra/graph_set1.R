library(ggplot2)
seed = 24125
nodes = 24
modules = 4
gen_mod = nodes/modules

x <- read.table(paste0("hist/", seed, "_data_set1.txt"), header = TRUE, sep = "\t")

t1 <- ks.test(x[,1], "pnorm")

pdf(paste0("hist/Mod_pred_", seed, "_set1.pdf"),width=6,height=4,paper='special')
ggplot(x, aes(x=x[,1])) + 
    geom_histogram() +
    labs(caption = paste0("mean = ", mean(x[,1]), "; sd = ", sd(x[,1]), "; p.val (ks.test) = ", t1$p.value))
dev.off()

sd(x[,1])

pdf(paste0("hist/Edges_", seed, "_set1.pdf"),width=6,height=4,paper='special')
hist(x[,2])
dev.off()

x <- read.table(paste0("hist/", seed, "_adj_set1.txt"), header = FALSE, sep = "\t")
x <- x[1:nodes, 1:nodes]
pru <- data.frame(Source = rep(1:nodes, nodes), 
                  Target = rep(1:nodes, each = nodes),
                  Frequency = as.vector(t(x)))
                  
adj_breaks <- (0:(modules-1))*gen_mod
adj_breaks
adj_labels <- strsplit(toString(adj_breaks+1), ", ")
adj_labels

pdf(paste0("hist/Meanadj_", seed, "_set1.pdf"),paper='special')
ggplot(data=pru, aes(x=Source, y = Target)) + geom_tile(aes(fill=Frequency)) + scale_fill_gradient(low = "white",high = "black",limits=c(0.01, 0.8)) + scale_x_continuous(breaks=adj_breaks, labels=adj_labels[[1]]) + scale_y_continuous(breaks=adj_breaks, labels= adj_labels[[1]])
dev.off()



# x <- read.table(paste0("hist/", seed, "_data_set2.txt"), header = TRUE, sep = "\t")
# 
# pdf(paste0("hist/Mod_pred_", seed, "_set2.pdf"),width=6,height=4,paper='special')
# ggplot(x, aes(x=x[,1])) + 
#     geom_histogram() +
#     labs(caption = paste0("mean = ", mean(x[,1]), "; sd = ", sd(x[,1])))
# dev.off()
# 
# sd(x[,1])
# 
# pdf(paste0("hist/Edges_", seed, "_set2.pdf"),width=6,height=4,paper='special')
# hist(x[,2])
# dev.off()
