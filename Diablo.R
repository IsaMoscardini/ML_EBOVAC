a 
######
rm(list = ls())
options(stringsAsFactors = F)
library(mixOmics)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(Glimma)
library(edgeR)
library(tidyverse)

##### PLS-DA mRNA + Abs ####
##### Day 01 vs Day 0 ####
# Open tables
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
AEs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)

View(AEs)
summary(Abs)

lumi <- read.delim("Data/D1/Xv0_Log2_Luminex_D1vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
n_lumi <- rownames(lumi)
View(lumi)

prots <- read.delim("Data/D1/Xv0_Log2_Olink_D1vD0.txt")
rownames(prots) <- prots$Probes
prots$Probes <- NULL
prots <- t(prots)
n_prots <- rownames(prots)
View(prots)

metab <- read.delim("Data/D1/Xv0_Metabolomic_normalized_D1vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
n_metab <- rownames(metab)
View(metab)

lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
prots <- apply(prots, 2, as.numeric)
rownames(prots) <- n_prots
metab <- apply(metab, 2, as.numeric)
rownames(metab) <- n_metab

ind <- unique(c(rownames(mRNA), rownames(prots), rownames(metab), rownames(lumi), names(Abs)))
tabela <- as.data.frame(ind)
View(tabela)

tabela$mRNA <- ifelse(tabela$ind %in% rownames(mRNA), 1, 0)
tabela$Lumi <- ifelse(tabela$ind %in% rownames(lumi), 1, 0)
tabela$Prots <- ifelse(tabela$ind %in% rownames(prots), 1, 0)
tabela$Metab <- ifelse(tabela$ind %in% rownames(metab), 1, 0)
tabela$Ab <- ifelse(tabela$ind %in% names(Abs), 1, 0)
tabela$AE <- ifelse(tabela$ind %in% names(AEs), 1, 0)
View(tabela)
rownames(tabela) <- tabela$ind
tabela$ind <- NULL
#library(ComplexHeatmap)
Heatmap(tabela)


# intersect of samples
mRNA <- mRNA[names(Abs),]

ind <- Reduce(intersect, list(rownames(mRNA),rownames(prots),rownames(metab),rownames(lumi)))

#
###### PLS #####
# sPLS_DA Luminex
MyResult.plsda <- plsda(lumi,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS lumi")
plotLoadings(MyResult.plsda, contrib = 'max', method = 'mean', comp = 2)

# sPLS_DA Proteina
MyResult.plsda <- plsda(as.matrix(prots),Abs_prot)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS prot")
plotLoadings(MyResult.plsda, contrib = 'max', method = 'mean', comp = 2)

# sPLS_DA Metabolite
MyResult.plsda <- plsda(metab,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS prot")
plotLoadings(MyResult.plsda, contrib = 'max', method = 'mean', comp = 1)

View(mRNA)


#### Diablo #####
lumi_prot <- lumi[ind,]
mRNA_prot <- mRNA[ind,]
metab_prot <- metab[ind,]
prots <- prots[ind,]
View(mRNA_prot)

identical(rownames(lumi_prot), rownames(prots))

X <- list(mRNA = mRNA_prot, 
          Luminex = lumi_prot, 
          Protein = prots,
          Metabolites = metab_prot)
Y <- Abs
summary(Y)

# test
list.keepX <- list(mRNA = c(100, 100), Luminex = c(9,6), Protein = c(100, 50), Metabolites = c(50, 50))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
  plotIndiv(MyResult.diablo)






mi <- read.delim("Data/D7/Xv0_Log2_miRNA_Normalized_Pseudocounts_D7vD0.txt")
View(mi)
