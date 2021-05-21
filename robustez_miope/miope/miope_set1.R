library(ggplot2)
library(tidyr)

seed = 24125
samplesize = 400

#primero para miope_mod
dmod <- read.csv(paste0("miope_mod/", seed, "_mod_set1.txt"), header = T, sep = "\t")
drob <- read.csv(paste0("miope_mod/", seed, "_rob_set1.txt"), header = T, sep = "\t")
dcon <- read.csv(paste0("miope_mod/", seed, "_con_set1.txt"), header = T, sep = "\t")
dtime <- read.csv(paste0("miope_mod/", seed, "_time_set1.txt"), header = F, sep = "\t")

new_data <- data.frame(incmod = dmod$Final - dmod$Initial,
                       incrob = drob$Final - drob$Initial)
                       
cat(paste0("Mod_i: Mean = ", mean(dmod$Initial), "; SD = ", sd(dmod$Initial) ), file="miope_mod/stats_set1.txt",sep="\n")
cat(paste0("Mod_f: Mean = ", mean(dmod$Final), "; SD = ", sd(dmod$Final) ), file="miope_mod/stats_set1.txt",sep="\n", append=TRUE)
cat(paste0("Rob_i: Mean = ", mean(drob$Initial), "; SD = ", sd(drob$Initial) ), file="miope_mod/stats_set1.txt",sep="\n", append=TRUE)
cat(paste0("Rob_f: Mean = ", mean(drob$Final), "; SD = ", sd(drob$Final) ), file="miope_mod/stats_set1.txt",sep="\n", append=TRUE)

t1 <- wilcox.test(dmod$Initial, dmod$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_mod/", seed, "_modi_modf_set1.pdf"),
    width=6,height=4,paper='special')
ggplot(data = dmod, aes(x = Initial, y = Final)) +
    geom_point() +
    geom_abline() +
    labs(caption = paste0("p_value = ", t1$p.value))
dev.off()

t2 <- wilcox.test(drob$Initial, drob$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_mod/", seed, "_robi_robf_set1.pdf"),
    width=6,height=4,paper='special')
ggplot(data = drob, aes(x = Initial, y = Final)) +
    geom_point() +
    geom_abline() +
    labs(caption = paste0("p_value = ", t2$p.value))
dev.off()

t3 <- wilcox.test(dcon$Initial, dcon$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_mod/", seed, "_coni_conf_set1.pdf"),
    width=6,height=4,paper='special')
ggplot(data = dcon, aes(x = Initial, y = Final)) +
    geom_point() +
    labs(caption = paste0("p_value = ", t3$p.value))
dev.off()

pdf(paste0("miope_mod/", seed, "_time_set1.pdf"),
    width=6,height=4,paper='special')
boxplot(dtime)
dev.off()
        
t4 <- cor.test(new_data$incmod, new_data$incrob,  method = "spearman")
pdf(paste0("miope_mod/", seed, "_inc.pdf"),
    width=6,height=4,paper='special')
ggplot(data = new_data, aes(x = incmod, y = incrob)) +
    geom_point() +
    labs(caption = paste0("n = ", nrow(new_data),"; p_value = ", t4$p.value,"; r = ", t4$estimate))
dev.off()

    
#ahora para miope_rob
dmod <- read.csv(paste0("miope_rob/", seed, "_mod_set1.txt"), header = T, sep = "\t")
drob <- read.csv(paste0("miope_rob/", seed, "_rob_set1.txt"), header = T, sep = "\t")
dcon <- read.csv(paste0("miope_rob/", seed, "_con_set1.txt"), header = T, sep = "\t")
dtime <- read.csv(paste0("miope_rob/", seed, "_time_set1.txt"), header = F, sep = "\t")

new_data <- data.frame(incmod = dmod$Final - dmod$Initial,
                       incrob = drob$Final - drob$Initial)
                       
cat(paste0("Mod_i: Mean = ", mean(dmod$Initial), "; SD = ", sd(dmod$Initial) ), file="miope_rob/stats_set1.txt",sep="\n")
cat(paste0("Mod_f: Mean = ", mean(dmod$Final), "; SD = ", sd(dmod$Final) ), file="miope_rob/stats_set1.txt",sep="\n", append=TRUE)
cat(paste0("Rob_i: Mean = ", mean(drob$Initial), "; SD = ", sd(drob$Initial) ), file="miope_rob/stats_set1.txt",sep="\n", append=TRUE)
cat(paste0("Rob_f: Mean = ", mean(drob$Final), "; SD = ", sd(drob$Final) ), file="miope_rob/stats_set1.txt",sep="\n", append=TRUE)

t1 <- wilcox.test(dmod$Initial, dmod$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_rob/", seed, "_modi_modf_set1.pdf"),
    width=6,height=4,paper='special')
ggplot(data = dmod, aes(x = Initial, y = Final)) +
    geom_point() +
    geom_abline() +
    labs(caption = paste0("p_value = ", t1$p.value))
dev.off()

t2 <- wilcox.test(drob$Initial, drob$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_rob/", seed, "_robi_robf_set1.pdf"),
    width=6,height=4,paper='special')
ggplot(data = drob, aes(x = Initial, y = Final)) +
    geom_point() +
    geom_abline() +
    labs(caption = paste0("p_value = ", t2$p.value))
dev.off()

t3 <- wilcox.test(dcon$Initial, dcon$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_rob/", seed, "_coni_conf_set1.pdf"),
    width=6,height=4,paper='special')
ggplot(data = dcon, aes(x = Initial, y = Final)) +
    geom_point() +
    labs(caption = paste0("p_value = ", t3$p.value))
dev.off()

pdf(paste0("miope_rob/", seed, "_time_set1.pdf"),
    width=6,height=4,paper='special')
boxplot(dtime)
dev.off()

t4 <- cor.test(new_data$incmod, new_data$incrob,  method = "spearman")
pdf(paste0("miope_rob/", seed, "_inc.pdf"),
    width=6,height=4,paper='special')
ggplot(data = new_data, aes(x = incrob, y = incmod)) +
    geom_point() +
    labs(caption = paste0("n = ", nrow(new_data),"; p_value = ", t4$p.value,"; r = ", t4$estimate))
dev.off()
