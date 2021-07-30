library(ggplot2)
library(gridExtra)
library(dplyr)
seed = 24125
nu_exp = 600
nu_int = 2
red_per_int = nu_exp/nu_int
       

xtime1 <- read.table(paste0("primer_opt/", seed, "_time.txt"), header = FALSE, sep = "\t")
time_1 <- data.frame(Time = xtime1$V1, 
                     Mod = as.factor(rep(c("modb","moda"), each = red_per_int)))
time_1$Mod <- factor(time_1$Mod , levels=c("modb", "moda"))

pdf(paste0("primer_opt/", seed, "_timep.pdf"),paper='special')
    ggplot(time_1, aes(x=Mod, y=Time)) + 
        geom_boxplot()
dev.off() 

titulo <- paste0("Tiempo")
capture.output(titulo, file = paste0("primer_opt/", seed, "_stat.txt"), append = FALSE)
titulo <- ""
capture.output(titulo, file = paste0("primer_opt/", seed, "_stat.txt"), append = TRUE)

labmb <- paste0("Tiempo modb = ", mean(time_1$Time[time_1$Mod == "modb"]))
labmbsd <- paste0("SD modb = ", sd(time_1$Time[time_1$Mod == "modb"]))
labma <- paste0("Tiempo moda = ", mean(time_1$Time[time_1$Mod == "moda"]))
labmasd <- paste0("SD moda = ", sd(time_1$Time[time_1$Mod == "moda"]))

capture.output(labmb, file = paste0("primer_opt/", seed, "_stat.txt"), append = TRUE)
capture.output(labmbsd, file = paste0("primer_opt/", seed, "_stat.txt"), append = TRUE)
capture.output(labma, file = paste0("primer_opt/", seed, "_stat.txt"), append = TRUE)
capture.output(labmasd, file = paste0("primer_opt/", seed, "_stat.txt"), append = TRUE)
titulo <- " "
capture.output(titulo, file = paste0("primer_opt/", seed, "_stat.txt"), append = TRUE)

titulo <- paste0("Wilcoxon tests")
t1 <- wilcox.test(time_1$Time~time_1$Mod, paired = FALSE, alternative = "greater")
capture.output(t1, file = paste0("primer_opt/", seed, "_stat.txt"), append = TRUE)


mut1 <- read.csv("datos_primer_opt.txt", header = TRUE, sep = "\t")

t1 <- wilcox.test(mut1$Mutaciones~mut1$Mod, paired = FALSE)
pdf(paste0(seed, "_mutp.pdf"),paper='special')
    ggplot(mut1, aes(x=Mod, y=Mutaciones)) + 
        geom_boxplot() +
        labs(caption = paste0("p_value = ", t1$p.value, "; W = ", t1$statistic))
dev.off() 

pdf(paste0(seed, "mut_hist.pdf"),paper='special')
    ggplot(mut1, aes(x=Mutaciones, fill=Mod)) + 
    geom_histogram(color="black", alpha=0.5, position="identity") + 
    scale_fill_grey(labels = c(expression('Moda'),expression('Modb')))+  
    labs(x=expression('Genetic distance'), y = expression('Networks'))  + 
    theme_classic() + 
    theme(text = element_text(size=14), legend.position = c(.9, .9), legend.title = element_blank())
dev.off()

pdf(paste0(seed, "_mod.pdf"),paper='special')
    ggplot(mut1, aes(x=Mod, y=Mod_final)) + 
        geom_boxplot()
dev.off() 

mut_time <- data.frame(Time = time_1$Time,
                       Mutaciones = mut1$Mutaciones,
                       Mod = mut1$Mod)
t1 <- cor.test(mut_time$Time, mut_time$Mutaciones,  method = "pearson")

pdf(paste0(seed, "mut_time.pdf"),paper='special')
    ggplot(mut_time, aes(x=Time, y=Mutaciones)) + 
        geom_point(colour="black", shape=21, size = 2, alpha = 0.5, aes(fill = factor(Mod))) +
        labs(caption = paste0("p_value = ", t1$p.value,"; r = ", t1$estimate)) +
        theme(plot.caption = element_text(size = 8))
dev.off() 
