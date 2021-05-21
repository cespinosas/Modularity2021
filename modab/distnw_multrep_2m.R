library(ggplot2)
library(gridExtra)
library(dplyr)
seed = 24125

params = commandArgs(trailingOnly=TRUE)

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

nu_rep = params[1]
sdm = 1
xmodb <- read.table(paste0("distbtwnw/", seed, "_modb_", sdm, "k_", nu_rep , "rep_general_2m.txt"), header = TRUE, sep = "\t")
xmoda <- read.table(paste0("distbtwnw/", seed, "_moda_", sdm, "k_", nu_rep , "rep_general_2m.txt"), header = TRUE, sep = "\t")

which_mod <- c(rep("modb", nrow(xmodb)), rep("moda", nrow(xmoda)))

d_plot <- cbind(rbind(xmodb, xmoda), which_mod)
d_plot$dist_btw_nw <- as.factor(d_plot$dist_btw_nw)

pmod <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = mod_red2, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())

pedges <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = edges_red2, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pnupasos <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = nu_pasos, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
ptot_phen <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_phen, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
ptot_phen1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_phen1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
ptot_phen2 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_phen2, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pfrac_newpht <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_newpht, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
ptot_acphen <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_acphen, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
ptot_acphen1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_acphen1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
ptot_acphen2 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_acphen2, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pfrac_newacph <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_newacph, fill = which_mod)) +
    geom_boxplot() +
    labs(caption = paste0())
mylegend <- g_legend(pfrac_newacph)
pfrac_newacph <- pfrac_newacph + theme(legend.position = "none") 
    

pdf(paste0("distbtwnw/", seed, "_", sdm, "k_", nu_rep , "rep_controls_2m.pdf"),paper='special')
    grid.arrange(pmod, pedges, ptot_phen1, ptot_acphen1, mylegend)
dev.off()

pdf(paste0("distbtwnw/", seed, "_", sdm, "k_", nu_rep , "rep_general_2m.pdf"),paper='special')
    grid.arrange(pnupasos, ptot_phen, ptot_phen2, pfrac_newpht, ptot_acphen, ptot_acphen2, pfrac_newacph, mylegend)
dev.off()


xmodb <- read.table(paste0("distbtwnw/", seed, "_modb_", sdm, "k_", nu_rep , "rep_acphens_2m.txt"), header = TRUE, sep = "\t")
xmoda <- read.table(paste0("distbtwnw/", seed, "_moda_", sdm, "k_", nu_rep , "rep_acphens_2m.txt"), header = TRUE, sep = "\t")

which_mod <- c(rep("modb", nrow(xmodb)), rep("moda", nrow(xmoda)))

d_plot <- cbind(rbind(xmodb, xmoda), which_mod)
d_plot$dist_btw_nw <- as.factor(d_plot$dist_btw_nw)

ptot_acphens <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_acphens, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
psim_acphens <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = sim_acphens, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())

pmods_acphens <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = mods_acphens, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pacphens_d1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = acphens_d1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pacphens_m1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = acphens_m1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pfrac_acphens_d1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_acphens_d1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pfrac_acphens_m1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_acphens_m1, fill = which_mod)) +
    geom_boxplot() +
    labs(caption = paste0())
mylegend <- g_legend(pfrac_acphens_m1)
pfrac_acphens_m1 <- pfrac_acphens_m1 + theme(legend.position = "none") 
    
pdf(paste0("distbtwnw/", seed, "_", sdm, "k_", nu_rep , "rep_acphens_2m.pdf"),paper='special')
    grid.arrange(ptot_acphens, psim_acphens, pmods_acphens, pacphens_d1, pacphens_m1, pfrac_acphens_d1, pfrac_acphens_m1, mylegend)
dev.off()



xmodb <- read.table(paste0("distbtwnw/", seed, "_modb_", sdm, "k_", nu_rep , "rep_newphen_2m.txt"), header = TRUE, sep = "\t")
xmoda <- read.table(paste0("distbtwnw/", seed, "_moda_", sdm, "k_", nu_rep , "rep_newphen_2m.txt"), header = TRUE, sep = "\t")

which_mod <- c(rep("modb", nrow(xmodb)), rep("moda", nrow(xmoda)))

d_plot <- cbind(rbind(xmodb, xmoda), which_mod)
d_plot$dist_btw_nw <- as.factor(d_plot$dist_btw_nw)

ptot_newpht <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_newpht, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
psimnewpht <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = simnewpht, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())

pmodsnewpht <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = modsnewpht, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())

pnewpht_g1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = newpht_g1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pnewpht_d1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = newpht_d1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pnewpht_m1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = newpht_m1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
  
pfrac_newpht_g1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_newpht_g1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pfrac_newpht_d1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_newpht_d1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pfrac_newpht_m1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_newpht_m1, fill = which_mod)) +
    geom_boxplot() +
    labs(caption = paste0())
mylegend <- g_legend(pfrac_newpht_m1)
pfrac_newpht_m1 <- pfrac_newpht_m1 + theme(legend.position = "none") 

pdf(paste0("distbtwnw/", seed, "_", sdm, "k_", nu_rep, "rep_newpht_2m.pdf"),paper='special')
    grid.arrange(ptot_newpht, psimnewpht, pmodsnewpht, pnewpht_g1, pnewpht_d1, pnewpht_m1, pfrac_newpht_g1, pfrac_newpht_d1, pfrac_newpht_m1, mylegend)
dev.off()


ptot_newacph <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = tot_newacph, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
psimnewacph <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = simnewacph, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())

pmodsnewacph <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = modsnewacph, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pnewacph_d1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = newacph_d1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pnewacph_m1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = newacph_m1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pfrac_newacph_d1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_newacph_d1, fill = which_mod)) +
    geom_boxplot() +
    theme(legend.position = "none") +
    labs(caption = paste0())
    
pfrac_newacph_m1 <- ggplot(data = d_plot, aes(x = dist_btw_nw, y = frac_newacph_m1, fill = which_mod)) +
    geom_boxplot() +
    labs(caption = paste0())
mylegend <- g_legend(pfrac_newacph_m1)
pfrac_newacph_m1 <- pfrac_newacph_m1 + theme(legend.position = "none") 
    
pdf(paste0("distbtwnw/", seed, "_", sdm, "k_", nu_rep, "rep_newacph_2m.pdf"),paper='special')
    grid.arrange(ptot_newacph, psimnewacph, pmodsnewacph, pnewacph_d1, pnewacph_m1, pfrac_newacph_d1, pfrac_newacph_m1, mylegend)
dev.off()
