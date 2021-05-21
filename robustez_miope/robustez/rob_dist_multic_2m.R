library(ggplot2)
library(gridExtra)

seed = 24125
nu_exp = 5000

x <- read.table(paste0("hist/", seed, "_data_2m.txt"), header = TRUE, sep = "\t")
Mod <- x[1:nu_exp,1]
Rob_1m <- read.table(paste0("rob/2mic_", seed, "_allmut_2m.txt"), header = TRUE, sep = "\t")

sum(Rob_1m$trunc)
data_p <- data.frame(Mod_pred = Mod,
                     Rob_1m = Rob_1m$Rob,
                     Rob_by_dist = Rob_1m$Rob_by_dist,
                     Rob_prom = Rob_1m$Rob_prom,
                     Rob_gan = Rob_1m$Rob_gan/Rob_1m$Cuan_gan,
                     Rob_per = Rob_1m$Rob_per/Rob_1m$Cuan_per,
                     Dist_gan = Rob_1m$Rob_by_dist_gan/Rob_1m$Cuan_gan,
                     Dist_per = Rob_1m$Rob_by_dist_per/Rob_1m$Cuan_per,
                     Mrob_gan = Rob_1m$Rob_prom_gan/Rob_1m$Cuan_gan,
                     Mrob_per = Rob_1m$Rob_prom_per/Rob_1m$Cuan_per)
                     
data_pz <- data.frame(Mod_pred = Mod,
                     Z = (data_p$Rob_by_dist - data_p$Rob_prom)/(1 - data_p$Rob_prom),
                     Zgan = (data_p$Dist_gan - data_p$Mrob_gan)/(1 - data_p$Mrob_gan),
                     Zper = (data_p$Dist_per - data_p$Mrob_per)/(1 - data_p$Mrob_per))

t1 <- cor.test(data_p$Mod_pred, data_p$Rob_1m,  method = "pearson")
prob <- ggplot(data_p, aes(x = Mod_pred, y = Rob_1m)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$Rob_by_dist,  method = "pearson")
pdist <- ggplot(data_p, aes(x = Mod_pred, y = Rob_by_dist)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_pz$Mod_pred, data_pz$Z,  method = "pearson")
pz <- ggplot(data_pz, aes(x = Mod_pred, y = Z)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_pz),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

pdf(paste0("rob/Allmut_", seed, "_2m.pdf"),paper='special')
    grid.arrange(prob, pdist, pz)
dev.off()

t1 <- cor.test(data_p$Mod_pred, data_p$Rob_gan,  method = "pearson")
prgan <- ggplot(data_p, aes(x = Mod_pred, y = Rob_gan)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
    theme(plot.caption = element_text(size = 8))
    
t1 <- cor.test(data_p$Mod_pred, data_p$Dist_gan,  method = "pearson")
pdgan <- ggplot(data_p, aes(x = Mod_pred, y = Dist_gan)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
    theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_pz$Mod_pred, data_pz$Zgan,  method = "pearson")
pzgan <- ggplot(data_pz, aes(x = Mod_pred, y = Zgan)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_pz),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
    theme(plot.caption = element_text(size = 8))

pdf(paste0("rob/Gan_", seed, "_2m.pdf"),paper='special')
    grid.arrange(prgan, pdgan, pzgan)
dev.off()


t1 <- cor.test(data_p$Mod_pred, data_p$Rob_per,  method = "pearson")
prgan <- ggplot(data_p, aes(x = Mod_pred, y = Rob_per)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
    theme(plot.caption = element_text(size = 8))
    
t1 <- cor.test(data_p$Mod_pred, data_p$Dist_per,  method = "pearson")
pdgan <- ggplot(data_p, aes(x = Mod_pred, y = Dist_per)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
    theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_pz$Mod_pred, data_pz$Zper,  method = "pearson")
pzper <- ggplot(data_pz, aes(x = Mod_pred, y = Zper)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_pz),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
    theme(plot.caption = element_text(size = 8))

pdf(paste0("rob/Per_", seed, "_2m.pdf"),paper='special')
    grid.arrange(prgan, pdgan, pzper)
dev.off()


# data_sin2 <- data.frame(Mod_pred = Mod,
#                      Rob_1m = (Rob_1m$Rob_gan + Rob_1m$Rob_per)/(2*24*24),
#                      Rob_by_dist = (Rob_1m$Rob_by_dist_gan + Rob_1m$Rob_by_dist_per)/(2*24*24),
#                      Rob_prom = (Rob_1m$Rob_prom_gan + Rob_1m$Rob_prom_per)/(2*24*24))
#         
# data_pzsin2 <- data.frame(Mod_pred = Mod,
#                      Z = (data_sin2$Rob_by_dist - data_sin2$Rob_prom)/(1 - data_p$Rob_prom))
#                      
# t1 <- cor.test(data_sin2$Mod_pred, data_sin2$Rob_1m,  method = "pearson")
# probsin2 <- ggplot(data_sin2, aes(x = Mod_pred, y = Rob_1m)) +
#     geom_point(alpha = 0.25) +
#     labs(caption = paste0("n = ", nrow(data_sin2),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
#   theme(plot.caption = element_text(size = 8))
# 
# t1 <- cor.test(data_sin2$Mod_pred, data_sin2$Rob_by_dist,  method = "pearson")
# pdistsin2 <- ggplot(data_p, aes(x = Mod_pred, y = Rob_by_dist)) +
#     geom_point(alpha = 0.25) +
#     labs(caption = paste0("n = ", nrow(data_sin2),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
#   theme(plot.caption = element_text(size = 8))
# 
# t1 <- cor.test(data_pzsin2$Mod_pred, data_pzsin2$Z,  method = "pearson")
# pzsin2 <- ggplot(data_pzsin2, aes(x = Mod_pred, y = Z)) +
#     geom_point(alpha = 0.25) +
#     labs(caption = paste0("n = ", nrow(data_pzsin2),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
#   theme(plot.caption = element_text(size = 8))
# 
# pdf(paste0("rob/Sin2_", seed, "_set1.pdf"),paper='special')
#     grid.arrange(probsin2, pdistsin2, pzsin2)
# dev.off()
