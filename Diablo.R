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

##### DIABLO Abs ####
##### Day 01 vs Day 0 ####
# Open tables
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
AEs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)

#View(AEs)
summary(Abs)

lumi <- read.delim("Data/D1/Xv0_Log2_Luminex_D1vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
n_lumi <- rownames(lumi)
#View(lumi)

prots <- read.delim("Data/D1/Xv0_Log2_Olink_D1vD0.txt")
rownames(prots) <- prots$Probes
prots$Probes <- NULL
prots <- t(prots)
n_prots <- rownames(prots)
#View(prots)

metab <- read.delim("Data/D1/Xv0_Metabolomic_normalized_D1vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
n_metab <- rownames(metab)
#View(metab)

lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
prots <- apply(prots, 2, as.numeric)
rownames(prots) <- n_prots
metab <- apply(metab, 2, as.numeric)
rownames(metab) <- n_metab

#######à
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
#mRNA <- mRNA[names(Abs),]

ind <- Reduce(intersect, list(rownames(mRNA),rownames(metab),rownames(lumi), rownames(prots), names(Abs)))


#### Diablo 
lumi <- lumi[ind,]
colnames(lumi) <- paste0("l", "_", colnames(lumi))
mRNA <- mRNA[ind,]
metab <- metab[ind,]
prots <- prots[ind,]
colnames(prots) <- paste0("p", "_", colnames(prots))  
View(mRNA)

identical(rownames(lumi), rownames(metab))

X <- list(mRNA = mRNA, 
          Luminex = lumi,
          Prot = prots,
          Metabolites = metab)

summary(Abs)
Abs <- Abs[ind]
View(Abs)
Y <- Abs

# test
list.keepX <- list(mRNA = c(90, 9), Luminex = c(6,5), Prot = c(95,6), Metabolites = c(5, 5))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo)


design = matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0
design 

sgccda.res = block.splsda(X = X, Y = Y, ncomp = 5, 
                          design = design)
set.seed(123) # for reproducibility, only when the `cpus' argument is not used
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10, progressBar = T)
plot(perf.diablo) 

perf.diablo$choice.ncomp$WeightedVote
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
ncomp <- 2

sgccda.res = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

plotDiablo(sgccda.res, ncomp = 1)

plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = TRUE, cutoff = 0.5, 
        pch = c(16, 17, 15, 14), cex = c(1,1,1,1), col = c('darkorchid', 'brown1', 'lightgreen', 'pink'))

circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1', 'lightgreen','pink'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

x11()
network(sgccda.res, blocks = c(1,2,3,4),
        color.node = c('darkorchid', 'brown1', 'lightgreen', 'pink'), cutoff = 0.6)

plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median')

x11()
cimDiablo(sgccda.res)

# performance
set.seed(123)
perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, folds = 4,
                   dist = 'centroids.dist')
perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate
auc.splsda = auroc(sgccda.res, roc.block = "Luminex", roc.comp = 1)

##### Day 03 vs Day 0 ####
# Open tables
mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
AEs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
summary(Abs)

lumi <- read.delim("Data/D7/Xv0_Log2_Luminex_D7vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
n_lumi <- rownames(lumi)
#View(lumi)

prots <- read.delim("Data/D7/Xv0_Log2_Olink_D7vD0.txt")
rownames(prots) <- prots$Probes
prots$Probes <- NULL
prots <- t(prots)
n_prots <- rownames(prots)
#View(prots)

metab <- read.delim("Data/D7/Xv0_Metabolomic_normalized_D7vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
n_metab <- rownames(metab)
#View(metab)

miRNA <- read.delim("Data/D7/Xv0_Log2_miRNA_Normalized_Pseudocounts_D7vD0.txt")
rownames(miRNA) <- miRNA$Probes
miRNA$Probes <- NULL
miRNA <- t(miRNA)
n_miRNA <- rownames(miRNA)
view(miRNA)

lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
prots <- apply(prots, 2, as.numeric)
rownames(prots) <- n_prots
metab <- apply(metab, 2, as.numeric)
rownames(metab) <- n_metab
miRNA <- apply(miRNA, 2, as.numeric)
rownames(miRNA) <- n_miRNA

ind <- unique(c(rownames(mRNA), rownames(miRNA), rownames(prots), rownames(metab), rownames(lumi), names(Abs)))
tabela <- as.data.frame(ind)
View(tabela)

tabela$mRNA <- ifelse(tabela$ind %in% rownames(mRNA), 1, 0)
tabela$miRNA <- ifelse(tabela$ind %in% rownames(miRNA), 1, 0)
tabela$Lumi <- ifelse(tabela$ind %in% rownames(lumi), 1, 0)
tabela$Prots <- ifelse(tabela$ind %in% rownames(prots), 1, 0)
tabela$Metab <- ifelse(tabela$ind %in% rownames(metab), 1, 0)
tabela$Ab <- ifelse(tabela$ind %in% names(Abs), 1, 0)
tabela$AE <- ifelse(tabela$ind %in% names(AEs), 1, 0)
#View(tabela)
rownames(tabela) <- tabela$ind
tabela$ind <- NULL
#library(ComplexHeatmap)
Heatmap(tabela)

# intersect of samples
mRNA <- mRNA[names(Abs),]

ind <- Reduce(intersect, list(rownames(mRNA),rownames(metab),rownames(lumi)))

##### Day 07 vs Day 0 ####
# Open tables
mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
AEs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
summary(Abs)

lumi <- read.delim("Data/D7/Xv0_Log2_Luminex_D7vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
n_lumi <- rownames(lumi)
#View(lumi)

prots <- read.delim("Data/D7/Xv0_Log2_Olink_D7vD0.txt")
rownames(prots) <- prots$Probes
prots$Probes <- NULL
prots <- t(prots)
n_prots <- rownames(prots)
#View(prots)

metab <- read.delim("Data/D7/Xv0_Metabolomic_normalized_D7vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
n_metab <- rownames(metab)
#View(metab)

miRNA <- read.delim("Data/D7/Xv0_Log2_miRNA_Normalized_Pseudocounts_D7vD0.txt")
rownames(miRNA) <- miRNA$Probes
miRNA$Probes <- NULL
miRNA <- t(miRNA)
n_miRNA <- rownames(miRNA)
view(miRNA)

lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
prots <- apply(prots, 2, as.numeric)
rownames(prots) <- n_prots
metab <- apply(metab, 2, as.numeric)
rownames(metab) <- n_metab
miRNA <- apply(miRNA, 2, as.numeric)
rownames(miRNA) <- n_miRNA

ind <- unique(c(rownames(mRNA), rownames(miRNA), rownames(prots), rownames(metab), rownames(lumi), names(Abs)))
tabela <- as.data.frame(ind)
View(tabela)

tabela$mRNA <- ifelse(tabela$ind %in% rownames(mRNA), 1, 0)
tabela$miRNA <- ifelse(tabela$ind %in% rownames(miRNA), 1, 0)
tabela$Lumi <- ifelse(tabela$ind %in% rownames(lumi), 1, 0)
tabela$Prots <- ifelse(tabela$ind %in% rownames(prots), 1, 0)
tabela$Metab <- ifelse(tabela$ind %in% rownames(metab), 1, 0)
tabela$Ab <- ifelse(tabela$ind %in% names(Abs), 1, 0)
tabela$AE <- ifelse(tabela$ind %in% names(AEs), 1, 0)
#View(tabela)
rownames(tabela) <- tabela$ind
tabela$ind <- NULL
library(ComplexHeatmap)
Heatmap(tabela)

# intersect of samples
mRNA <- mRNA[names(Abs),]
ind <- Reduce(intersect, list(rownames(mRNA),rownames(metab),rownames(lumi), rownames(miRNA), names(Abs)))


#### Diablo 
lumi <- lumi[ind,]
mRNA <- mRNA[ind,]
metab <- metab[ind,]
prots <- prots[ind,]
miRNA <- miRNA[ind,] 
View(mRNA)

identical(rownames(lumi), rownames(metab))

X <- list(mRNA = mRNA, 
          Luminex = lumi,
          miRNA = miRNA,
          Metabolites = metab)

summary(Abs)
Abs <- Abs[ind]
View(Abs)
Y <- Abs

# test
list.keepX <- list(mRNA = c(90, 70), Luminex = c(5,6), miRNA = c(5,45), Metabolites = c(55, 6))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo)


design = matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0
design 

sgccda.res = block.splsda(X = X, Y = Y, ncomp = 5, 
                          design = design)
set.seed(123) # for reproducibility, only when the `cpus' argument is not used
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10, progressBar = T)
x11()
plot(perf.diablo) 

perf.diablo$choice.ncomp$WeightedVote
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
ncomp <- 2

sgccda.res = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

plotDiablo(sgccda.res, ncomp = 1)

plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = TRUE, cutoff = 0.5, 
        pch = c(16, 17, 15, 14), cex = c(1,1,1,1), col = c('darkorchid', 'brown1', 'lightgreen', 'pink'))

circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1', 'lightgreen','pink'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

network(sgccda.res, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)

plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median')

x11()
cimDiablo(sgccda.res)

# performance
set.seed(123)
perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, folds = 4,
                   dist = 'centroids.dist')
perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate
auc.splsda = auroc(sgccda.res, roc.block = "Metabolites", roc.comp = 1)

##### DIABLO AEs 
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
#View(prots)

metab <- read.delim("Data/D1/Xv0_Metabolomic_normalized_D1vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
n_metab <- rownames(metab)
#View(metab)

lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
prots <- apply(prots, 2, as.numeric)
rownames(prots) <- n_prots
metab <- apply(metab, 2, as.numeric)
rownames(metab) <- n_metab

#######à
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
#mRNA <- mRNA[names(Abs),]

ind <- Reduce(intersect, list(rownames(mRNA),rownames(metab),rownames(lumi), rownames(prots), names(AEs)))


#### Diablo
lumi <- lumi[ind,]
mRNA <- mRNA[ind,]
metab <- metab[ind,]
prots <- prots[ind,]
View(mRNA)

identical(rownames(lumi), rownames(metab))

X <- list(mRNA = mRNA, 
          Luminex = lumi,
          Prot = prots,
          Metabolites = metab)

summary(AEs)
AEs <- AEs[ind]
View(AEs)
Y <- AEs

# test
list.keepX <- list(mRNA = c(100, 50), Luminex = c(9,10), Prot = c(70,9), Metabolites = c(20, 6))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo)


design = matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0
design 

sgccda.res = block.splsda(X = X, Y = Y, ncomp = 5, 
                          design = design)
set.seed(123) # for reproducibility, only when the `cpus' argument is not used
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10, progressBar = T)
plot(perf.diablo) 

perf.diablo$choice.ncomp$WeightedVote
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
ncomp <- 2

sgccda.res = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

plotDiablo(sgccda.res, ncomp = 1)

plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = TRUE, cutoff = 0.5, 
        pch = c(16, 17, 15, 14), cex = c(1,1,1,1), col = c('darkorchid', 'brown1', 'lightgreen', 'pink'))

circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1', 'lightgreen','pink'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

network(sgccda.res, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)

plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median')

x11()
cimDiablo(sgccda.res)

# performance
set.seed(123)
perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, folds = 4,
                   dist = 'centroids.dist')
perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate
auc.splsda = auroc(sgccda.res, roc.block = "Metabolites", roc.comp = 1)
