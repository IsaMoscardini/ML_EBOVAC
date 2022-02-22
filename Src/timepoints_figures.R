#####
rm(list = ls())
options(stringsAsFactors = F)
install.packages("stringr")

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  library(gplots)
  library(data.table)
  library(dplyr)
  library(tidyverse)
  library(reshape2)
})


d1 <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
rownames(d1) <- d1$Probes
d1$Probes <- NULL
d1 <- t(d1)
View(d1)

d3 <- read.delim("Data/D3/RNAseq_logCPM_D3vD0.txt")
rownames(d3) <- d3$Probes
d3$Probes <- NULL
d3 <- t(d3)

d7 <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
rownames(d7) <- d7$Probes
d7$Probes <- NULL
d7 <- t(d7)

d14 <- read.delim("Data/D14/RNAseq_logCPM_D14vD0.txt")
rownames(d14) <- d14$Probes
d14$Probes <- NULL
d14 <- t(d14)

d28 <- read.delim("Data/D28/RNAseq_logCPM_D28vD0.txt")
rownames(d28) <- d28$Probes
d28$Probes <- NULL
d28 <- t(d28)

pheno <- read.delim("Data/Outcomes.txt")
View(pheno)
ab <- pheno[,c(1,3)]
rownames(ab) <- ab$Probes

### juntar tudo
melt1 <- melt(d1)
melt1$TP <- 1

melt3 <- melt(d3)
melt3$TP <- 3

melt7 <- melt(d7)
melt7$TP <- 7

melt14 <- melt(d14)
melt14$TP <- 14

melt28 <- melt(d28)
melt28$TP <- 28
#View(melt28)

melt_all <- rbind(melt1, melt3, melt7, melt14, melt28)
View(melt_all)
#melt_all <- melt_all[melt_all$Var1 %in% ab$Probes, ]
melt_all$class <- ab$AntibodyResponse_Class[match(melt_all$Var1, ab$Probes)]
melt_all <- melt_all[complete.cases(melt_all$class),]
socs3 <- melt_all[melt_all$Var2 %in% c("MTOR"),]
#View(socs3)


ggplot(socs3, aes(x = class, y = value, col = class, group=class)) + stat_compare_means(method = "wilcox.test", paired = FALSE) +  
  geom_jitter() + geom_boxplot() + theme_minimal() + facet_grid(. ~ TP)

ggplot(socs3, aes(x = class, y = value, col = class, group=class)) + stat_compare_means(method = "wilcox.test", paired = FALSE) +  
  geom_jitter() + geom_boxplot() + theme_minimal() + facet_grid(. ~ Var2)

ggplot(socs3, aes(x = TP, y = value, col = class, group = Var1)) + geom_line() +
   theme_minimal() 


ggplot(socs3, aes(x = TP, y = value, col = class, group = Var1)) + geom_line() +
  theme_minimal() + geom_function()

 



