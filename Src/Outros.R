#### --------------- Estudo GSE130991 OBESIDADE --------------- ####

Estudo <- getGEO("GSE130991")

Expressao <- exprs(Estudo[[1]])

Pheno <- pData(Estudo[[1]])
View(Pheno)

Featuredata <- fData(Estudo[[1]])

dir.create(path = 'Data//GSE130991')

write.table(Expressao, 'Data/GSE130991/Expres.txt', quote = FALSE, sep = '\t')

write.table(Pheno, 'Data/GSE130991/Pheno.txt', quote = FALSE, sep = '\t')

write.table(Featuredata, 'Data/GSE130991/Feature.txt', quote = FALSE, sep = '\t')

#getGEOSuppFiles('GSE5418',baseDir = 'Data/GSE5418/')


##### Baixar estudos malaria para teste ####
library(GEOquery)
library(data.table)

#GSE52166 <- getGEO("GSE52166")
#pheno <- pData(GSE52166[[1]])
#write.table(pheno, "GSE52166/phenodata.txt", sep = '\t')

#supp <- getGEOSuppFiles("GSE52166")

files_count <- list.files("GSE52166/", pattern = "counts", full.names = TRUE)

counts <- lapply(files_count, fread, data.table = FALSE)
head(counts[[1]])
lapply(counts, head)
basename(files_count)

names(counts) <- gsub("_.*", "", basename(files_count))
count_named <- lapply(names(counts), function(x) {
  colnames(counts[[x]]) <- c("GeneID", x)
  counts[[x]]})

counts <- Reduce(function(x, y) merge(x, y, by= "GeneID", all=TRUE), count_named)
View(counts)

counts <- counts[grep(pattern = "ENSG", counts$GeneID), ]
View(counts)

write.table(counts, "GSE52166/countdata.txt", sep = "\t")

#####
count <- read.delim("GSE52166/countdata.txt")
View(count)

convertGS <- function(x){
  
  require("biomaRt")
  GS = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ENS = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("ensembl_gene_id"), values = x , mart = ENS, attributesL = c("hgnc_symbol"), martL = GS, uniqueRows=T)
  
  return(genesV2)
}

genesGS <- convertGS(count$GeneID)
View(genesGS)

library(dplyr)
counts <- merge(genesGS, as.data.frame(count), by.x = "Gene.stable.ID", by.y = "GeneID")
View(counts)

 probes <- c("ENSG00000228716","ENSG00000124134","ENSG00000224970","ENSG00000252159",
             "ENSG00000187583","ENSG00000207723","ENSG00000244749","ENSG00000238323",
             "ENSG00000249025","ENSG00000235324","ENSG00000224174","ENSG00000221423",
             "ENSG00000181609","ENSG00000212331","ENSG00000168875","ENSG00000216642",
             "ENSG00000227584","ENSG00000228219","ENSG00000251200","ENSG00000242038",
             "ENSG00000185633","ENSG00000012779","ENSG00000212418","ENSG00000231852",
             "ENSG00000227467","ENSG00000096088","ENSG00000257200","ENSG00000144410",
             "ENSG00000182518","ENSG00000248387","ENSG00000006074","ENSG00000114547",
             "ENSG00000251239","ENSG00000187186","ENSG00000137502","ENSG00000233151",
             "ENSG00000249729","ENSG00000196544","ENSG00000243847","ENSG00000253135",
             "ENSG00000213347","ENSG00000253620","ENSG00000197586","ENSG00000235431",
             "ENSG00000136286","ENSG00000253214","ENSG00000238586","ENSG00000255368",
             "ENSG00000236743","ENSG00000111224","ENSG00000227067","ENSG00000134709",
             "ENSG00000231030","ENSG00000250080")
 
 
 counts_filt <- counts[counts$Gene.stable.ID %in% probes,]
# View(counts_filt)

genes_ml <- counts_filt$HGNC.symbol



#######

Estudo5418 <- getGEO("GSE5418")
Expressao5418 <- exprs(Estudo5418[[1]])
PhenoData5418 <- pData(Estudo5418[[1]])
Featuredata5418 <- fData(Estudo5418[[1]])
View(Expressao5418)



####
rownames(counts) <- counts$HGNC.symbol
count$GeneID <- NULL
identical(rownames(pheno), colnames(count)) # TRUE 

table(is.na(count))
'%ni%' <- Negate('%in%')
table(colSums(count) > 1)
toremove <- c("GSM1335695", "GSM1335713")
count <- count[,colnames(count) %ni% toremove]

count <- count[complete.cases(count),]
View(count)

rnames <- rownames(count)

library(DESeq2)

counts <- varianceStabilizingTransformation(as.matrix(count))
rownames(counts) <- rnames
View(counts)



######### preprocess pheno ####
View(pheno)
table(pheno$characteristics_ch1)
#


pheno <- read.delim("GSE52166/phenodata.txt")
View(pheno)

pheno$Sample <- rownames(pheno)
pheno <- pheno[,c(1,10)]
View(pheno)
pheno$Sample <- rownames(pheno)
pheno$title <- NULL
colnames(pheno) <- c("Group", "Sample")

pheno$Group <- gsub("time-point: ","",  pheno$Group)
pheno$Group <- gsub("time point: ","",  pheno$Group)
pheno$Group <- gsub("-.*","",  pheno$Group)
pheno$Group <- tolower(pheno$Group)


########## preprocess counts ########
count <- read.delim("GSE52166/countdata.txt")
View(count)

rownames(count) <- count$GeneID
count$GeneID <- NULL
identical(rownames(pheno), colnames(count)) # TRUE 

table(is.na(count))
'%ni%' <- Negate('%in%')
table(colSums(count) > 1)
toremove <- c("GSM1335695", "GSM1335713")
count <- count[,colnames(count) %ni% toremove]

count <- count[complete.cases(count),]
View(count)

rnames <- rownames(count)

library(DESeq2)

counts <- varianceStabilizingTransformation(as.matrix(count))
rownames(counts) <- rnames
View(counts)

#
identical(rownames(pheno), colnames(count))

write.table(counts, "countdata_malaria.txt", sep= '\t')

pheno <- pheno[colnames(counts),]
View(pheno)

write.table(pheno, "pheno_malaria.txt", sep= '\t')

#### D1 ####
Ab_50 <- read.delim("Intermed/Escolha_features/D1/AE/52_Adverse_Xv0_Log2_Olink_D1vD0.txt")
#View(Ab_50)
genes_ab <- Ab_50$Probes

mRNA <- read.delim("Data/D1/Xv0_Log2_Olink_D1vD0.txt")
rownames(mRNA) <- mRNA$Probes
mRNA$Probes <- NULL
mRNA <- t(mRNA)
View(mRNA)

outcome <- read.delim("Data/Outcomes.txt")
outcome <- outcome[complete.cases(outcome$AdverseEvent_Class),]
outcome$AdverseEvent_Class <- as.factor(outcome$AdverseEvent_Class)
Abs <- setNames(outcome$AdverseEvent_Class, outcome$Probes)
summary(Abs)

n_mrna <- rownames(mRNA)
mRNA <- apply(mRNA, 2, as.numeric)
rownames(mRNA) <- n_mrna

mRNA <- mRNA[complete.cases(mRNA),]

inter <- intersect(rownames(mRNA), names(Abs)) 
Abs <- Abs[inter]
mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[inter,]

# PLS-DA
pca.srbct = pca(mRNA, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA all genes')
mRNA <- as.data.frame(mRNA)
mRNA <- mRNA[, which(colnames(mRNA) %in% genes_ab)]
#View(mRNA)

pca.srbct = pca(mRNA, ncomp = 3, center = TRUE, scale = TRUE)
plot(pca.srbct)  # screeplot of the eingenvalues (explained variance per component)
plotIndiv(pca.srbct, group = Abs, ind.names = FALSE, 
          legend = TRUE, title = 'PCA AE genes D28')



# Compare data 
Ab_50 <- read.csv("Intermed/AB_leandro/TOP50.csv")
Ab_50_meu <- read.csv("Results/sPLS/mRNA/PLS_DA/Ab/Table_PLSDA_D1_Ab.csv")
Ab_50_meu <- Ab_50_meu[order(abs(Ab_50_meu$comp1), decreasing = T),]
Ab_50_meu <- Ab_50_meu[1:50,]
View(Ab_50_meu)
inter <- intersect(Ab_50$Probes, Ab_50_meu$X)
View(inter)

Ae_50 <- read.csv("Intermed/AE_leandro/TOP50.csv")
Ae_50_meu <- read.csv("Results/sPLS/mRNA/PLS_DA/AE/Table_PLSDA_D1_AE.csv")
Ae_50_meu  <- Ae_50_meu[order(abs(Ae_50_meu$comp1), decreasing = T),]
Ae_50_meu <- Ae_50_meu[1:50,]
View(Ae_50_meu)

inter <- intersect(Ae_50$Probes, Ae_50_meu$X)
View(inter)

########## tabelas genes com funçoes biologicas
filt <- read.csv("Intermed/Genes_biol_function.csv")
View(filt)
table(duplicated(filt$x))

D1 <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
