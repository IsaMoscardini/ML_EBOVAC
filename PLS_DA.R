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
View(Abs)
summary(Abs)

mRNA <- mRNA[names(Abs),]

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
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/Ab/Tune/D1_vs_D0/Gene_table.csv")

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
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- mRNA[names(Abs),]

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
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/Ab/Tune/D3_vs_D0/Gene_table.csv")

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
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

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
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
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
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/Ab/Tune/D7_vs_D0/Gene_table.csv")

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
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

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
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
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
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/Ab/Tune/D14_vs_D0/Gene_table.csv")

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
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AntibodyResponse_Class),]
outcome$AntibodyResponse_Class <- as.factor(outcome$AntibodyResponse_Class)
Abs <- setNames(outcome$AntibodyResponse_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

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
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
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
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/Ab/Tune/D28_vs_D0/Gene_table.csv")

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
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
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

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(90, 70))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D1")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/AE/Tune/D1_vs_D0/Gene_table.csv")

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
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- mRNA[names(Abs),]

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
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
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
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/AE/Tune/D3_vs_D0/Gene_table.csv")

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
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

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
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
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
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/AE/Tune/D7_vs_D0/Gene_table.csv")

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
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

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
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 2, # we suggest to push ncomp a bit more, e.g. 4
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
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/AE/Tune/D14_vs_D0/Gene_table.csv")

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
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
View(Abs)
summary(Abs)

mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[names(Abs),]

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
tune.splsda.srbct <- tune.splsda(mRNA, Abs, ncomp = 3, # we suggest to push ncomp a bit more, e.g. 4
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

MyResult.splsda.final <- splsda(mRNA, Abs, ncomp = 2, keepX = c(100, 30))

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA - D28")
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda.final$X))

write.csv(genes, "Results/PLS_DA/AE/Tune/D28_vs_D0/Gene_table.csv")

##########à
MyResult.splsda2 <- splsda(mRNA,Abs, ncomp=2, keepX=c(25,25))
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA D28 - 2 comp, 25 e 25")
plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
genes <- (t(MyResult.splsda2$X))

write.csv(genes, "Results/PLS_DA/AE/Using_25/D28_vs_D0/Gene_table.csv")
  