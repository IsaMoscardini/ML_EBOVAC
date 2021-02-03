#####
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


############################################### mRNA
##### PLS-DA mRNA + Abs ####
##### Day 01 ####
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
#pca.srbct #outputs the explained variance per component
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 6)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 6 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D1",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 6)


MyResult.plsda <- plsda(mRNA,Abs, ncomp = 6)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D1 - Ab")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/Ab/Table_PLSDA_D1_Ab.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(90, 9))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- (MyResult.splsda.final$loadings$X)
View(genes)
write.csv(genes, "Results/sPLS_DA/Ab/Tune/D1_vs_D0/Gene_table.csv")
    
##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/Ab/Using_25/D1_vs_D0/Gene_table.csv")



##### Day 03 ####
mRNA <- read.delim("Data/D3/RNAseq_logCPM_D3vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 5)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 5 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D3",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 1)
cim(srbct.plsda)

#
MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D3 - Ab")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/Ab/Table_PLSDA_D3_Ab.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(5, 30))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D3")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- (MyResult.splsda.final$loadings$X)
View(genes)
write.csv(genes, "Results/sPLS_DA/Ab/Tune/D3_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D3 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/Ab/Using_25/D3_vs_D0/Gene_table.csv")


##### Day 07 ####
mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 4)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 5 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D7",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 1)



MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D7 - Ab")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/Ab/Table_PLSDA_D7_Ab.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(2))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(90, 70))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D7")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- (MyResult.splsda.final$loadings$X)
View(genes)
write.csv(genes, "Results/sPLS_DA/Ab/Tune/D7_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D7 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/Ab/Using_25/D7_vs_D0/Gene_table.csv")


##### Day 14 ####
mRNA <- read.delim("Data/D14/RNAseq_logCPM_D14vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 3)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 5 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D14",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 2)

#
MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D14 - Ab")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/Ab/Table_PLSDA_D14_Ab.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(2))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(5, 30))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D14")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- (MyResult.splsda.final$loadings$X)
View(genes)
write.csv(genes, "Results/sPLS_DA/Ab/Tune/D14_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D14 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/Ab/Using_25/D14_vs_D0/Gene_table.csv")


##### Day 28 ####
mRNA <- read.delim("Data/D28/RNAseq_logCPM_D28vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 6)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 5 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D28",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 2)

MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D28 - Ab")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/Ab/Table_PLSDA_D28_Ab.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 5, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

#ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(2))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(90, 40))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D28")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- (MyResult.splsda.final$loadings$X)
View(genes)
write.csv(genes, "Results/sPLS_DA/Ab/Tune/D28_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D14 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/Ab/Using_25/D28_vs_D0/Gene_table.csv")




####à
##### PLS-DA mRNA + EAs ####
##### Day 01 ####
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 7)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 5 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D1",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 2)


MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D1 - AE")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/AE/Table_PLSDA_D1_AE.csv")


# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 5, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(ncomp))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(100, 50))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- MyResult.splsda.final$loadings$X
View(genes)
write.csv(genes, "Results/sPLS_DA/AE/Tune/D1_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/AE/Using_25/D1_vs_D0/Gene_table.csv")



##### Day 03 ####
mRNA <- read.delim("Data/D3/RNAseq_logCPM_D3vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 5)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 5 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D3",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 1)

MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D3 - AE")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/AE/Table_PLSDA_D3_AE.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(100, 20))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D3")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- MyResult.splsda.final$loadings$X
View(genes)
write.csv(genes, "Results/sPLS_DA/AE/Tune/D3_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D3 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/AE/Using_25/D3_vs_D0/Gene_table.csv")


##### Day 07 ####
mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 6)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 6 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D7",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 2)


MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D7 - AE")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/AE/Table_PLSDA_D7_AE.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(2))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(60, 40))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D7")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- MyResult.splsda.final$loadings$X
View(genes)
write.csv(genes, "Results/sPLS_DA/AE/Tune/D7_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D7 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/AE/Using_25/D7_vs_D0/Gene_table.csv")


##### Day 14 ####
mRNA <- read.delim("Data/D14/RNAseq_logCPM_D14vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 5)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 5 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D14",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 2)


MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D14 - AE")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/AE/Table_PLSDA_D14_AE.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(2))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(100, 5))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D14")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- MyResult.splsda.final$loadings$X
View(genes)
write.csv(genes, "Results/sPLS_DA/AE/Tune/D14_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D14 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/AE/Using_25/D14_vs_D0/Gene_table.csv")


##### Day 28 ####
mRNA <- read.delim("Data/D28/RNAseq_logCPM_D28vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
#View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
#View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA on SRBCT')

srbct.plsda <- plsda(mRNA, Abs, ncomp = 5)  # set ncomp to 10 for performance assessment later
plotIndiv(srbct.plsda , comp = 1:2,
          group = Abs, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA 5 comp')
background = background.predict(srbct.plsda, comp.predicted=2, dist = "max.dist") 
plotIndiv(srbct.plsda, comp = 1:2,
          group = Abs, ind.names = FALSE, title = "Maximum distance - D28",
          legend = TRUE,  background = background)

set.seed(2543) # for reproducibility, only when the `cpus' argument is not used
perf.plsda.srbct <- perf(srbct.plsda, validation = "Mfold", folds = 5, 
                         progressBar = TRUE, auc = TRUE, nrepeat = 10) 
x11()
plot(perf.plsda.srbct, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
auc.plsda = auroc(srbct.plsda, roc.comp = 1)


MyResult.plsda <- plsda(mRNA,Abs)
plotIndiv(MyResult.plsda, legend = T, ellipse = T, title =  "PLS-DA D28 - AE")
plotVar(MyResult.plsda, cutoff = 0.7) 
table_genes <- MyResult.plsda$loadings$X
View(table_genes)
write.csv(table_genes, "Results/PLS_DA/AE/Table_PLSDA_D28_AE.csv")

# geral
MyResult.plsda2 <- plsda(mRNA,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(2))

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(80, 80))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D28")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- MyResult.splsda.final$loadings$X
View(genes)
write.csv(genes, "Results/sPLS_DA/AE/Tune/D28_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D28 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/AE/Using_25/D28_vs_D0/Gene_table.csv")
  

############################################
##### Day 01 ####
############################################### Luminex
##### PLS-DA lumi + Abs ####
##### Day 01 ####
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
mRNA$Probes <- NULL
samples <- colnames(mRNA)

lumi <- read.delim("Data/D1/Xv0_Log2_Luminex_D1vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
View(lumi)
lumi <- as.data.frame(lumi)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

lumi <- lumi[names(Abs),]
identical(rownames(lumi), names(Abs)) # TRUE

# geral
MyResult.plsda2 <- plsda(lumi,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(lumi, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(lumi, Abs, ncomp = 2, keepX = c(6, 5))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Lumi/sPLSDA/Ab/D1/Citocinas_D1.csv")


##### Day 03 ####
mRNA <- read.delim("Data/D3/RNAseq_logCPM_D3vD0.txt")
mRNA$Probes <- NULL
samples <- colnames(mRNA)

lumi <- read.delim("Data/D3/Xv0_Log2_Luminex_D3vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
View(lumi)
n_lumi <- rownames(lumi)
lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
lumi <- as.data.frame(lumi)
lumi <- lumi[samples,]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
View(Abs)
summary(Abs)

samp <- intersect(names(Abs), rownames(lumi))
lumi <- lumi[samp,]
Abs <- Abs[samp]

identical(rownames(lumi), names(Abs)) # TRUE

# geral
MyResult.plsda2 <- plsda(lumi,Abs, ncomp=8)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(lumi, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 4, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(lumi, Abs, ncomp = 2, keepX = c(8, 7))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D3")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Lumi/sPLSDA/Ab/D3/Citocinas_D3.csv")


##### Day 07 ####
mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
mRNA$Probes <- NULL
samples <- colnames(mRNA)

lumi <- read.delim("Data/D7/Xv0_Log2_Luminex_D7vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
View(lumi)
n_lumi <- rownames(lumi)
lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
lumi <- as.data.frame(lumi)
lumi <- lumi[samples,]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
View(Abs)
summary(Abs)

samp <- intersect(names(Abs), rownames(lumi))
lumi <- lumi[samp,]
Abs <- Abs[samp]

identical(rownames(lumi), names(Abs)) # TRUE

# geral
MyResult.plsda2 <- plsda(lumi,Abs, ncomp=8)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 4, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(lumi, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 5, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(lumi, Abs, ncomp = 2, keepX = c(5, 6))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D7")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Lumi/sPLSDA/Ab/D3/Citocinas_D3.csv")

############################################### Metabolomic
############################################### Proteins
##### PLS-DA lumi + AEs ####
##### Day 01 ####
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
mRNA$Probes <- NULL
samples <- colnames(mRNA)

lumi <- read.delim("Data/D1/Xv0_Log2_Luminex_D1vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
#View(lumi)
lumi <- as.data.frame(lumi)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
#View(Abs)
summary(Abs)

lumi <- lumi[names(Abs),]
identical(rownames(lumi), names(Abs)) # TRUE

# geral
MyResult.plsda2 <- plsda(lumi,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(lumi, Abs, ncomp = 5, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(lumi, Abs, ncomp = 2, keepX = c(9, 10))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Lumi/sPLSDA/AE/D1/Citocinas_D1.csv")


##### Day 03 ####
mRNA <- read.delim("Data/D3/RNAseq_logCPM_D3vD0.txt")
mRNA$Probes <- NULL
samples <- colnames(mRNA)

lumi <- read.delim("Data/D3/Xv0_Log2_Luminex_D3vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
View(lumi)
n_lumi <- rownames(lumi)
lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
lumi <- as.data.frame(lumi)
lumi <- lumi[samples,]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
#View(Abs)
summary(Abs)

samp <- intersect(names(Abs), rownames(lumi))
lumi <- lumi[samp,]
Abs <- Abs[samp]

identical(rownames(lumi), names(Abs)) # TRUE

# geral
MyResult.plsda2 <- plsda(lumi,Abs, ncomp=8)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(lumi, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 4, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(lumi, Abs, ncomp = 2, keepX = c(6, 6))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D3")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Lumi/sPLSDA/AE/D3/Citocinas_D3.csv")


##### Day 07 ####
mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
mRNA$Probes <- NULL
samples <- colnames(mRNA)

lumi <- read.delim("Data/D7/Xv0_Log2_Luminex_D7vD0.txt")
rownames(lumi) <- lumi$Probes
lumi$Probes <- NULL
lumi <- t(lumi)
View(lumi)
n_lumi <- rownames(lumi)
lumi <- apply(lumi, 2, as.numeric)
rownames(lumi) <- n_lumi
lumi <- as.data.frame(lumi)
lumi <- lumi[samples,]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
#View(Abs)
summary(Abs)

samp <- intersect(names(Abs), rownames(lumi))
lumi <- lumi[samp,]
Abs <- Abs[samp]

identical(rownames(lumi), names(Abs)) # TRUE

# geral
MyResult.plsda2 <- plsda(lumi,Abs, ncomp=6)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(lumi, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(lumi, Abs, ncomp = 2, keepX = c(5, 6))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D7")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Lumi/sPLSDA/AE/D7/Citocinas_D7.csv")

############################################### Metabolomic
############################################### Proteins
##### PLS-DA Proteins + Abs ####
##### Day 01 ####
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
mRNA$Probes <- NULL
samples <- colnames(mRNA)

prot <- read.delim("Data/D1/Xv0_Log2_Olink_D1vD0.txt")
rownames(prot) <- prot$Probes
prot$Probes <- NULL
prot <- t(prot)
View(prot)
prot <- as.data.frame(prot)
n_prots <- rownames(prot)
prot <- apply(prot, 2, as.numeric)
rownames(prot) <- n_prots
#prot <- prot[complete.cases(prot),]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

prot <- prot[names(Abs),]

samp <- intersect(rownames(prot), names(Abs))
prot <- prot[samp,]
Abs <- Abs[samp]
identical(rownames(prot), names(Abs)) # TRUE

# geral
MyResult.plsda2 <- plsda(prot,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 2, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(lumi, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(lumi, Abs, ncomp = 2, keepX = c(6, 5))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Lumi/sPLSDA/Ab/D1/Citocinas_D1.csv")
##### PLS-DA metab + Abs ####
##### Day 01 ####
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
samples <- colnames(mRNA)
metab <- read.delim("Data/D1/Xv0_Metabolomic_normalized_D1vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
View(metab)
metab <- as.data.frame(metab)
metab <- metab[samples, ]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

samp <- intersect(rownames(metab), names(Abs))
metab <- metab[samp,]
Abs <- Abs[samp]
identical(names(Abs), rownames(metab))

# geral
MyResult.plsda2 <- plsda(metab,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(metab, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(metab, Abs, ncomp = 2, keepX = c(5, 5))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Metab/sPLS_DA/Ab/D1/Metab_D1.csv")

##### Day 03 ####
mRNA <- read.delim("Data/D3/RNAseq_logCPM_D3vD0.txt")
samples <- colnames(mRNA)
metab <- read.delim("Data/D3/Xv0_Metabolomic_normalized_D3vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
View(metab)
metab <- as.data.frame(metab)
metab <- metab[samples, ]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

samp <- intersect(rownames(metab), names(Abs))
metab <- metab[samp,]
Abs <- Abs[samp]
identical(names(Abs), rownames(metab))

# geral
MyResult.plsda2 <- plsda(metab,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(metab, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(metab, Abs, ncomp = 2, keepX = c(5, 6))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D3")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Metab/sPLS_DA/Ab/D3/Metab_D3.csv")

##### Day 07 ####
mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
samples <- colnames(mRNA)
metab <- read.delim("Data/D7/Xv0_Metabolomic_normalized_D7vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
View(metab)
metab <- as.data.frame(metab)
metab <- metab[samples, ]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

samp <- intersect(rownames(metab), names(Abs))
metab <- metab[samp,]
Abs <- Abs[samp]
identical(names(Abs), rownames(metab))

# geral
MyResult.plsda2 <- plsda(metab,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(metab, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(metab, Abs, ncomp = 2, keepX = c(55, 6))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D7")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Metab/sPLS_DA/Ab/D7/Metab_D7.csv")


##### PLS-DA metab + Abs ####
##### Day 01 ####
mRNA <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
samples <- colnames(mRNA)
metab <- read.delim("Data/D1/Xv0_Metabolomic_normalized_D1vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
View(metab)
metab <- as.data.frame(metab)
metab <- metab[samples, ]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
#View(Abs)
summary(Abs)

samp <- intersect(rownames(metab), names(Abs))
metab <- metab[samp,]
Abs <- Abs[samp]
identical(names(Abs), rownames(metab))

# geral
MyResult.plsda2 <- plsda(metab,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(metab, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(metab, Abs, ncomp = 2, keepX = c(20, 6, 5, 50))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Metab/sPLS_DA/AE/D1/Metab_D1.csv")

##### Day 03 ####
mRNA <- read.delim("Data/D3/RNAseq_logCPM_D3vD0.txt")
samples <- colnames(mRNA)
metab <- read.delim("Data/D3/Xv0_Metabolomic_normalized_D3vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
View(metab)
metab <- as.data.frame(metab)
metab <- metab[samples, ]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

samp <- intersect(rownames(metab), names(Abs))
metab <- metab[samp,]
Abs <- Abs[samp]
identical(names(Abs), rownames(metab))

# geral
MyResult.plsda2 <- plsda(metab,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(metab, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(metab, Abs, ncomp = 2, keepX = c(5, 6))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D3")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Metab/sPLS_DA/Ab/D3/Metab_D3.csv")

##### Day 07 ####
mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
samples <- colnames(mRNA)
metab <- read.delim("Data/D7/Xv0_Metabolomic_normalized_D7vD0.txt")
rownames(metab) <- metab$Sample
metab$Sample <- NULL
metab <- t(metab)
View(metab)
metab <- as.data.frame(metab)
metab <- metab[samples, ]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
#View(Abs)
summary(Abs)

samp <- intersect(rownames(metab), names(Abs))
metab <- metab[samp,]
Abs <- Abs[samp]
identical(names(Abs), rownames(metab))

# geral
MyResult.plsda2 <- plsda(metab,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(metab, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(metab, Abs, ncomp = 2, keepX = c(55, 6))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D7")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Metab/sPLS_DA/Ab/D7/Metab_D7.csv")



##### PLS-DA miRNA + Abs ####
##### Day 01 ####
#mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
#samples <- colnames(mRNA)
mirna <- read.delim("Data/D7/Xv0_Log2_miRNA_Normalized_Pseudocounts_D7vD0.txt")
rownames(mirna) <- mirna$Probes
mirna$Probes <- NULL
mirna <- t(mirna)
View(mirna)
mirna <- as.data.frame(mirna)
#mirna <- mirna[samples, ]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
View(Abs)
summary(Abs)

samp <- intersect(rownames(mirna), names(Abs))
mirna <- mirna[samp,]
Abs <- Abs[samp]
identical(names(Abs), rownames(mirna))

# geral
MyResult.plsda2 <- plsda(mirna,Abs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 4, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(metab, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(metab, Abs, ncomp = 2, keepX = c(5, 5))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Metab/sPLS_DA/Ab/D1/Metab_D1.csv")

##### PLS-DA miRNA + AEs ####
##### Day 01 ####
#mRNA <- read.delim("Data/D7/RNAseq_logCPM_D7vD0.txt")
#samples <- colnames(mRNA)
mirna <- read.delim("Data/D7/Xv0_Log2_miRNA_Normalized_Pseudocounts_D7vD0.txt")
rownames(mirna) <- mirna$Probes
mirna$Probes <- NULL
mirna <- t(mirna)
View(mirna)
mirna <- as.data.frame(mirna)
#mirna <- mirna[samples, ]

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
EAs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
summary(EAs)

samp <- intersect(rownames(mirna), names(EAs))
mirna <- mirna[samp,]
EAs <- EAs[samp]
identical(names(EAs), rownames(mirna))

# geral
MyResult.plsda2 <- plsda(mirna,EAs, ncomp=10)
set.seed(30)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50

x11()
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

MyPerf.plsda
list.keepX <- c(5:10,  seq(20, 100, 5))
list.keepX

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(metab, Abs, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = TRUE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

ncomp <- 2
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(4))

MyResult.splsda.final <- splsda(metab, Abs, ncomp = 2, keepX = c(5, 5))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
cito <- (MyResult.splsda.final$loadings$X)
View(cito)
write.csv(cito, "Results/Metab/sPLS_DA/Ab/D1/Metab_D1.csv")
