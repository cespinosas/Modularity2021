library(ggplot2)
library(tidyr)

seed = 24125
samplesize = 400

#primero para miope_mod
dmod <- read.csv(paste0("miope_mod/", seed, "_mod_2m.txt"), header = T, sep = "\t")
drob <- read.csv(paste0("miope_mod/", seed, "_rob_2m.txt"), header = T, sep = "\t")
dcon <- read.csv(paste0("miope_mod/", seed, "_con_2m.txt"), header = T, sep = "\t")
dtime <- read.csv(paste0("miope_mod/", seed, "_time_2m.txt"), header = F, sep = "\t")

cat(paste0("Mod_i: Mean = ", mean(dmod$Initial), "; SD = ", sd(dmod$Initial) ), file="miope_mod/stats_2m.txt",sep="\n")
cat(paste0("Mod_f: Mean = ", mean(dmod$Final), "; SD = ", sd(dmod$Final) ), file="miope_mod/stats_2m.txt",sep="\n", append=TRUE)
cat(paste0("Rob_i: Mean = ", mean(drob$Initial), "; SD = ", sd(drob$Initial) ), file="miope_mod/stats_2m.txt",sep="\n", append=TRUE)
cat(paste0("Rob_f: Mean = ", mean(drob$Final), "; SD = ", sd(drob$Final) ), file="miope_mod/stats_2m.txt",sep="\n", append=TRUE)

t1 <- wilcox.test(dmod$Initial, dmod$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_mod/", seed, "_modi_modf_2m.pdf"),
    width=6,height=4,paper='special')
ggplot(data = dmod, aes(x = Initial, y = Final)) +
    geom_point() +
    geom_abline() +
    labs(caption = paste0("p_value = ", t1$p.value))
dev.off()

t2 <- wilcox.test(drob$Initial, drob$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_mod/", seed, "_robi_robf_2m.pdf"),
    width=6,height=4,paper='special')
ggplot(data = drob, aes(x = Initial, y = Final)) +
    geom_point() +
    geom_abline() +
    labs(caption = paste0("p_value = ", t2$p.value))
dev.off()

t3 <- wilcox.test(dcon$Initial, dcon$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_mod/", seed, "_coni_conf_2m.pdf"),
    width=6,height=4,paper='special')
ggplot(data = dcon, aes(x = Initial, y = Final)) +
    geom_point() +
    labs(caption = paste0("p_value = ", t3$p.value))
dev.off()

pdf(paste0("miope_mod/", seed, "_time_2m.pdf"),
    width=6,height=4,paper='special')
boxplot(dtime)
dev.off()


#ahora para miope_rob
dmod <- read.csv(paste0("miope_rob/", seed, "_mod_2m.txt"), header = T, sep = "\t")
drob <- read.csv(paste0("miope_rob/", seed, "_rob_2m.txt"), header = T, sep = "\t")
dcon <- read.csv(paste0("miope_rob/", seed, "_con_2m.txt"), header = T, sep = "\t")
dtime <- read.csv(paste0("miope_rob/", seed, "_time_2m.txt"), header = F, sep = "\t")

cat(paste0("Mod_i: Mean = ", mean(dmod$Initial), "; SD = ", sd(dmod$Initial) ), file="miope_rob/stats_2m.txt",sep="\n")
cat(paste0("Mod_f: Mean = ", mean(dmod$Final), "; SD = ", sd(dmod$Final) ), file="miope_rob/stats_2m.txt",sep="\n", append=TRUE)
cat(paste0("Rob_i: Mean = ", mean(drob$Initial), "; SD = ", sd(drob$Initial) ), file="miope_rob/stats_2m.txt",sep="\n", append=TRUE)
cat(paste0("Rob_f: Mean = ", mean(drob$Final), "; SD = ", sd(drob$Final) ), file="miope_rob/stats_2m.txt",sep="\n", append=TRUE)

t1 <- wilcox.test(dmod$Initial, dmod$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_rob/", seed, "_modi_modf_2m.pdf"),
    width=6,height=4,paper='special')
ggplot(data = dmod, aes(x = Initial, y = Final)) +
    geom_point() +
    geom_abline() +
    labs(caption = paste0("p_value = ", t1$p.value))
dev.off()

t2 <- wilcox.test(drob$Initial, drob$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_rob/", seed, "_robi_robf_2m.pdf"),
    width=6,height=4,paper='special')
ggplot(data = drob, aes(x = Initial, y = Final)) +
    geom_point() +
    geom_abline() +
    labs(caption = paste0("p_value = ", t2$p.value))
dev.off()

t3 <- wilcox.test(dcon$Initial, dcon$Final, exact = FALSE, alternative = "less")
pdf(paste0("miope_rob/", seed, "_coni_conf_2m.pdf"),
    width=6,height=4,paper='special')
ggplot(data = dcon, aes(x = Initial, y = Final)) +
    geom_point() +
    labs(caption = paste0("p_value = ", t3$p.value))
dev.off()

pdf(paste0("miope_rob/", seed, "_time_2m.pdf"),
    width=6,height=4,paper='special')
boxplot(dtime)
dev.off()
