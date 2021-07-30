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
