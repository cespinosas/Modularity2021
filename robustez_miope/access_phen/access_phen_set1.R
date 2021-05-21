library(ggplot2)
library(gridExtra)

seed = 24125
nu_exp = 5000

x <- read.table(paste0("hist/", seed, "_data_set1.txt"), header = TRUE, sep = "\t")
Mod_pred <- x[1:nu_exp,1]
Accphen <- read.table(paste0("access_phen/", seed, "_allmut_set1.txt"), header = TRUE, sep = "\t")

sum(Accphen$Trunc)
length(Accphen$Trunc)
data_p <- cbind(Mod_pred, Accphen)


#Ahora fenotipos totales

t1 <- cor.test(data_p$Mod_pred, data_p$tot_pht,  method = "pearson")
ptot_pht <- ggplot(data_p, aes(x = Mod_pred, y = tot_pht)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$simpht,  method = "pearson")
psimpht <- ggplot(data_p, aes(x = Mod_pred, y = simpht)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$modspht,  method = "pearson")
pmodspht <- ggplot(data_p, aes(x = Mod_pred, y = modspht)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
 
t1 <- cor.test(data_p$Mod_pred, data_p$pht_g1,  method = "pearson")
ppht_g1 <- ggplot(data_p, aes(x = Mod_pred, y = pht_g1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$pht_d1,  method = "pearson")
ppht_d1 <- ggplot(data_p, aes(x = Mod_pred, y = pht_d1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$pht_10,  method = "pearson")
ppht_10 <- ggplot(data_p, aes(x = Mod_pred, y = pht_10)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$pht_m1,  method = "pearson")
ppht_m1 <- ggplot(data_p, aes(x = Mod_pred, y = pht_m1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$frac_pht_g1,  method = "pearson")
pfrac_pht_g1 <- ggplot(data_p, aes(x = Mod_pred, y = frac_pht_g1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$frac_pht_d1,  method = "pearson")
pfrac_pht_d1 <- ggplot(data_p, aes(x = Mod_pred, y = frac_pht_d1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$frac_pht_10,  method = "pearson")
pfrac_pht_10 <- ggplot(data_p, aes(x = Mod_pred, y = frac_pht_10)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$frac_pht_m1,  method = "pearson")
pfrac_pht_m1 <- ggplot(data_p, aes(x = Mod_pred, y = frac_pht_m1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
pdf(paste0("access_phen/Phen_tot_", seed, "_set1.pdf"),paper='special')
    grid.arrange(ptot_pht, psimpht, pmodspht, ppht_g1, ppht_d1, ppht_10, ppht_m1, pfrac_pht_g1, pfrac_pht_d1, pfrac_pht_10, pfrac_pht_m1)
dev.off()

#Ahora fenotipos acumulados
t1 <- cor.test(data_p$Mod_pred, data_p$tot_acph,  method = "pearson")
ptot_acph <- ggplot(data_p, aes(x = Mod_pred, y = tot_acph)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$simacph,  method = "pearson")
psimacph <- ggplot(data_p, aes(x = Mod_pred, y = simacph)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$modsacph,  method = "pearson")
pmodsacph <- ggplot(data_p, aes(x = Mod_pred, y = modsacph)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$acph_d1,  method = "pearson")
pacph_d1 <- ggplot(data_p, aes(x = Mod_pred, y = acph_d1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$acph_10,  method = "pearson")
pacph_10 <- ggplot(data_p, aes(x = Mod_pred, y = acph_10)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$acph_m1,  method = "pearson")
pacph_m1 <- ggplot(data_p, aes(x = Mod_pred, y = acph_m1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$frac_acph_d1,  method = "pearson")
pfrac_acph_d1 <- ggplot(data_p, aes(x = Mod_pred, y = frac_acph_d1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$frac_acph_10,  method = "pearson")
pfrac_acph_10 <- ggplot(data_p, aes(x = Mod_pred, y = frac_acph_10)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$frac_acph_m1,  method = "pearson")
pfrac_acph_m1 <- ggplot(data_p, aes(x = Mod_pred, y = frac_acph_m1)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"\n r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
pdf(paste0("access_phen/Acum_phen_", seed, "_set1.pdf"),paper='special')
    grid.arrange(ptot_acph, psimacph, pmodsacph, pacph_d1, pacph_10, pacph_m1, pfrac_acph_d1, pfrac_acph_10, pfrac_acph_m1)
dev.off()
  
  
#Ahora gráfica de datos generales
t1 <- cor.test(data_p$Mod_pred, data_p$att_max_size,  method = "pearson")
patt_max_size <- ggplot(data_p, aes(x = Mod_pred, y = att_max_size)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("n = ", nrow(data_p),"; p_value = ", t1$p.value,"; r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$att_mean_size,  method = "pearson")
patt_mean_size <- ggplot(data_p, aes(x = Mod_pred, y = att_mean_size)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"; r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))

t1 <- cor.test(data_p$Mod_pred, data_p$dist_btw_acph,  method = "pearson")
pdist_btw_acph <- ggplot(data_p, aes(x = Mod_pred, y = dist_btw_acph)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"; r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
t1 <- cor.test(data_p$Mod_pred, data_p$Trunc,  method = "pearson")
pTrunc <- ggplot(data_p, aes(x = Mod_pred, y = Trunc)) +
    geom_point(alpha = 0.25) +
    labs(caption = paste0("p_value = ", t1$p.value,"; r = ", t1$estimate)) +
  theme(plot.caption = element_text(size = 8))
  
pdf(paste0("access_phen/General_", seed, "_set1.pdf"),paper='special')
    grid.arrange(patt_max_size, patt_mean_size, pdist_btw_acph, pTrunc)
dev.off()
