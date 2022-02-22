### ATENCAO: para cada estudo:
# 1. aperte ctrl+shift+F10 (ou clica no "Session" (no menu la em cima) e depois em "Restart R")
# 2. rode o cabecalho

#### inicio do cabecalho ####
rm(list = ls())

options(stringsAsFactors = FALSE)

### ATENCAO: ve se o nome da pasta esta certo aqui:
#setwd("C:/Users/aline/Desktop/MALAIRA/TCC_")

#--- bibliotecas que vao ser usadas
library(GEOquery)
#### fim do cabecalho ####

#### --------------- Estudo Tuberculos --------------- ####
Estudo19439 <- getGEO('GSE19439')

#Baixar dados

Expressao19439 <- exprs(Estudo19439[[1]])

Pheno19439 <- pData(Estudo19439[[1]])

Feature19439 <- fData(Estudo19439[[1]])

#Salvar dados nos locais corretos
dir.create(path = 'Data//GSE19439')
write.table(Expressao19439, "Data//GSE19439//ExpressaoGSE19439.txt", quote = FALSE, sep = '\t')

write.table(Pheno19439, 'Data/GSE35858/Phenodata35858.txt', quote = FALSE, sep = '\t' )

write.table(Feature35858, 'Data/GSE35858/Featuredata35858.txt', quote = FALSE, sep = '\t' )

# - inicio da alteracao - #
# nao precisa baixar os dados brutos desse estudo
#getGEOSuppFiles('GSE35858',baseDir = 'Data/GSE35858/')
# - fim da alteracao - #



#### --------------- Estudo GSE5418 --------------- ####

Estudo5418 <- getGEO("GSE5418")

Expressao5418 <- exprs(Estudo5418[[1]])

PhenoData5418 <- pData(Estudo5418[[1]])

Featuredata5418 <- fData(Estudo5418[[1]])

dir.create(path = 'Data//GSE5418')

write.table(Expressao5418, 'Data/GSE5418/Expressao5418.txt', quote = FALSE, sep = '\t')

write.table(PhenoData5418, 'Data/GSE5418/PhenoData5418.txt', quote = FALSE, sep = '\t')

write.table(Featuredata5418, 'Data/GSE5418/FeatureData5418.txt', quote = FALSE, sep = '\t')

getGEOSuppFiles('GSE5418',baseDir = 'Data/GSE5418/')


#### --------------- Estudo GSE34404 --------------- ####

Estudo34404 <- getGEO("GSE34404")

Expressao34404 <- exprs(Estudo34404[[1]])

PhenoData34404 <- pData(Estudo34404[[1]])

FeatureData <- fData(Estudo34404[[1]])

dir.create(path = 'Data//GSE34404')

write.table(Expressao34404, 'Data/GSE34404/Expressao34404.txt', quote = FALSE, sep = '\t')

write.table(PhenoData34404, 'Data/GSE34404/PhenoData34404.txt', quote = FALSE, sep = '\t')

write.table(FeatureData, 'Data/GSE34404/FeatureData34404.txt', quote = FALSE, sep = '\t')

getGEOSuppFiles("GSE34404",baseDir = 'Data/GSE34404/')
