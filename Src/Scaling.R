
##### Divide Ebola data between batches OR NOT and scale all datasets

##### Head ---------------------------------------------------------------------
rm(list = ls()) 
options(stringsAsFactors = F) 

pkgs <- c('edgeR','scales','impute','naniar','mixOmics','reshape2')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)


##### Data divided by batches --------------------------------------------------
##### EBOVAC ####
##### Dia 1 --------------------------------------------------------------------
geneva_counts <- read.csv2("Data/Counts_Geneva.csv")
rownames(geneva_counts) <- geneva_counts$X
geneva_pheno <- read.csv2("Data/Pheno_Geneva.csv")
geneva_pheno$Sample <- paste0("X",geneva_pheno$Sample)
geneva_pheno <- geneva_pheno[geneva_pheno$Treatment == "V",]
#View(geneva_pheno)

geneva_pheno_b1 <- geneva_pheno[geneva_pheno$Library.Batch == 1,]
geneva_pheno_b2 <- geneva_pheno[geneva_pheno$Library.Batch == 2,]

inter_b1 <- intersect(geneva_pheno_b1$Sample, colnames(geneva_counts))
inter_b2 <- intersect(geneva_pheno_b2$Sample, colnames(geneva_counts))

geneva_counts_b1 <- geneva_counts[,inter_b1]
geneva_counts_b2 <- geneva_counts[,inter_b2]
#View(geneva_counts_b2)

## Select D1  
geneva_pheno_b1 <- geneva_pheno_b1[geneva_pheno_b1$Group.Day %in% c(0,1),]
geneva_pheno_b2 <- geneva_pheno_b2[geneva_pheno_b2$Group.Day %in% c(0,1),]
#View(geneva_pheno_b1)

inter_b1 <- intersect(geneva_pheno_b1$Sample, colnames(geneva_counts_b1))
inter_b2 <- intersect(geneva_pheno_b2$Sample, colnames(geneva_counts_b2))

geneva_counts_b1 <- geneva_counts_b1[,inter_b1] # 51 samples
geneva_counts_b2 <- geneva_counts_b2[,inter_b2] # 38 samples
#View(geneva_counts_b2)

## Filtering and normalization
# B1
y_ebovac <- DGEList(counts = geneva_counts_b1, genes = row.names(geneva_counts_b1), group= geneva_pheno_b1$Group.Day)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 12245    51

y.1 <- calcNormFactors(y_ebovac)
b1 <- log2(cpm(y.1)+1)
b1 <- as.data.frame(b1)
b1$Probes <- rownames(b1)
dim(b1)
b1 <- b1[,c(52, 1:51)]
min(b1[,2:length(b1)])
#View(b1)

geneva_pheno_b1 <- geneva_pheno_b1[,c(1,4)]
colnames(geneva_pheno_b1) <- c("Probes", "Class")
#View(geneva_pheno_b1)
geneva_pheno_b1$Class <- paste0("D", geneva_pheno_b1$Class)

write.table(geneva_pheno_b1, "Intermed/Scaled_Ebolas_Apr2022/pheno_ebovac_b1_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')


# B2
y_ebovac_2 <- DGEList(counts = geneva_counts_b2, genes = row.names(geneva_counts_b2), group= geneva_pheno_b2$Group.Day)

keep2 <- rowSums(cpm(y_ebovac_2)>1)>=10
y_ebovac_2 <- y_ebovac_2[keep2,]
dim(y_ebovac_2) # 12366    38

y.2 <- calcNormFactors(y_ebovac_2)
b2 <- log2(cpm(y.2)+1)
b2 <- as.data.frame(b2)
b2$Probes <- rownames(b2)
dim(b2)
b2 <- b2[,c(39, 1:38)]
View(b2)

geneva_pheno_b2 <- geneva_pheno_b2[,c(1,4)]
colnames(geneva_pheno_b2) <- c("Probes", "Class")
geneva_pheno_b2$Class <- paste0("D", geneva_pheno_b2$Class)
#View(geneva_pheno_b2)

write.table(geneva_pheno_b2, "Intermed/Scaled_Ebolas_Apr2022/pheno_ebovac_b2_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')


##### EBOPLUS ####
eboplus <- read.delim("Intermed/data_BioFeats_USA/final_data/norm_counts_d1.txt")
#View(eboplus)


##### VEBCON ####
vebcon <- read.delim("Intermed/BioFeatS_EBOP_EBOV_VEB/D1/vebcon_norm_counts_d1.txt")
#View(vebcon)

#inter <- Reduce(intersect, list(b1$Probes, b2$Probes, eboplus$Probes, vebcon$Probes))

#b1 <- b1[b1$Probes %in% inter,] # 9904 genes
#b2 <- b2[b2$Probes %in% inter,] # 9904 genes
#eboplus <- eboplus[eboplus$Probes %in% inter,] # 9904 genes
#vebcon <- vebcon[vebcon$Probes %in% inter,] # 9904 genes

## Prepare to scale
rownames(eboplus) <- eboplus$Probes
rownames(vebcon) <- vebcon$Probes

# namesb1 <- rownames(b1)
# namesb2 <- rownames(b2)
# nameseboplus <- rownames(eboplus)
# namesvebcon <- rownames(vebcon)

##### Pegar so genes importantes ####
genes <- read.csv("Genes_bio.csv")
#View(genes)

b1 <- b1[which(b1$Probes %in% genes$x),]
#View(b1)
b2 <- b2[which(b2$Probes %in% genes$x),]
eboplus <- eboplus[which(eboplus$Probes %in% genes$x),]
vebcon <- vebcon[which(vebcon$Probes %in% genes$x),]
#View(vebcon)

##### Scale all tables ####
b1$Probes <- NULL
b2$Probes <- NULL
eboplus$Probes <- NULL
vebcon$Probes <- NULL

min(b2)
b1 <- rescale(as.matrix(b1))
b2 <- rescale(as.matrix(b2))
eboplus <- rescale(as.matrix(eboplus))
vebcon <- rescale(as.matrix(vebcon))
#View(eboplus)

b1 <- as.data.frame(b1)
b2 <- as.data.frame(b2)
eboplus <- as.data.frame(eboplus)
vebcon <- as.data.frame(vebcon)

b1$Probes <- rownames(b1)
b2$Probes <- rownames(b2)
eboplus$Probes <- rownames(eboplus)
vebcon$Probes <- rownames(vebcon)

dim(vebcon)
b1 <- b1[,c(52,1:51)]
b2 <- b2[,c(39,1:38)]
eboplus <- eboplus[,c(64,1:63)]
vebcon <- vebcon[,c(27, 1:26)]
View(vebcon)

#####
write.table(b1, "Intermed/Scaled_Ebolas_Apr2022/scaled_ebovac_b1_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(b2, "Intermed/Scaled_Ebolas_Apr2022/scaled_ebovac_b2_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(eboplus, "Intermed/Scaled_Ebolas_Apr2022/scaled_eboplus_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(vebcon, "Intermed/Scaled_Ebolas_Apr2022/scaled_vebcon_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')


#

##### Dia 3 --------------------------------------------------------------------
geneva_counts <- read.csv2("Data/Counts_Geneva.csv")
rownames(geneva_counts) <- geneva_counts$X
geneva_pheno <- read.csv2("Data/Pheno_Geneva.csv")
geneva_pheno$Sample <- paste0("X",geneva_pheno$Sample)
geneva_pheno <- geneva_pheno[geneva_pheno$Treatment == "V",]
#View(geneva_pheno)

geneva_pheno_b1 <- geneva_pheno[geneva_pheno$Library.Batch == 1,]
geneva_pheno_b2 <- geneva_pheno[geneva_pheno$Library.Batch == 2,]

inter_b1 <- intersect(geneva_pheno_b1$Sample, colnames(geneva_counts))
inter_b2 <- intersect(geneva_pheno_b2$Sample, colnames(geneva_counts))

geneva_counts_b1 <- geneva_counts[,inter_b1]
geneva_counts_b2 <- geneva_counts[,inter_b2]
#View(geneva_counts_b2)

## Select D1  
geneva_pheno_b1 <- geneva_pheno_b1[geneva_pheno_b1$Group.Day %in% c(0,3),]
geneva_pheno_b2 <- geneva_pheno_b2[geneva_pheno_b2$Group.Day %in% c(0,3),]
#View(geneva_pheno_b1)

inter_b1 <- intersect(geneva_pheno_b1$Sample, colnames(geneva_counts_b1))
inter_b2 <- intersect(geneva_pheno_b2$Sample, colnames(geneva_counts_b2))

geneva_counts_b1 <- geneva_counts_b1[,inter_b1]
geneva_counts_b2 <- geneva_counts_b2[,inter_b2]
#View(geneva_counts_b2)

## Filtering and normalization
# B1
y_ebovac <- DGEList(counts = geneva_counts_b1, genes = row.names(geneva_counts_b1), group= geneva_pheno_b1$Group.Day)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) #  11810    20

y.1 <- calcNormFactors(y_ebovac)
b1 <- log2(cpm(y.1)+1)
b1 <- as.data.frame(b1)
b1$Probes <- rownames(b1)
dim(b1)
b1 <- b1[,c(21, 1:20)]
min(b1[,2:length(b1)])
#View(b1)

geneva_pheno_b1 <- geneva_pheno_b1[,c(1,4)]
colnames(geneva_pheno_b1) <- c("Probes", "Class")
#View(geneva_pheno_b1)
geneva_pheno_b1$Class <- paste0("D", geneva_pheno_b1$Class)

write.table(geneva_pheno_b1, "Intermed/Scaled_Ebolas_Apr2022/pheno_ebovac_b1_d3.txt", quote = FALSE, row.names = FALSE, sep = '\t')


# B2
y_ebovac_2 <- DGEList(counts = geneva_counts_b2, genes = row.names(geneva_counts_b2), group= geneva_pheno_b2$Group.Day)

keep2 <- rowSums(cpm(y_ebovac_2)>1)>=10
y_ebovac_2 <- y_ebovac_2[keep2,]
dim(y_ebovac_2) # 12737    52

y.2 <- calcNormFactors(y_ebovac_2)
b2 <- log2(cpm(y.2)+1)
b2 <- as.data.frame(b2)
b2$Probes <- rownames(b2)
dim(b2)
b2 <- b2[,c(53, 1:52)]
View(b2)

geneva_pheno_b2 <- geneva_pheno_b2[,c(1,4)]
colnames(geneva_pheno_b2) <- c("Probes", "Class")
geneva_pheno_b2$Class <- paste0("D", geneva_pheno_b2$Class)
#View(geneva_pheno_b2)

write.table(geneva_pheno_b2, "Intermed/Scaled_Ebolas_Apr2022/pheno_ebovac_b2_d3.txt", quote = FALSE, row.names = FALSE, sep = '\t')


##### EBOPLUS ####
eboplus <- read.delim("Intermed/data_BioFeats_USA/final_data/norm_counts_d3.txt")
#View(eboplus)


##### VEBCON ####
vebcon <- read.delim("Intermed/BioFeatS_EBOP_EBOV_VEB/vebcon_norm_counts_d3.txt")
#View(vebcon)

#inter <- Reduce(intersect, list(b1$Probes, b2$Probes, eboplus$Probes, vebcon$Probes))

#b1 <- b1[b1$Probes %in% inter,] # 9904 genes
#b2 <- b2[b2$Probes %in% inter,] # 9904 genes
#eboplus <- eboplus[eboplus$Probes %in% inter,] # 9904 genes
#vebcon <- vebcon[vebcon$Probes %in% inter,] # 9904 genes

## Prepare to scale
rownames(eboplus) <- eboplus$Probes
rownames(vebcon) <- vebcon$Probes

# namesb1 <- rownames(b1)
# namesb2 <- rownames(b2)
# nameseboplus <- rownames(eboplus)
# namesvebcon <- rownames(vebcon)

##### Pegar so genes importantes ####
genes <- read.csv("Genes_bio.csv")
#View(genes)

b1 <- b1[which(b1$Probes %in% genes$x),]
#View(b1)
b2 <- b2[which(b2$Probes %in% genes$x),]
eboplus <- eboplus[which(eboplus$Probes %in% genes$x),]
vebcon <- vebcon[which(vebcon$Probes %in% genes$x),]
#View(vebcon)

##### Scale all tables ####
b1$Probes <- NULL
b2$Probes <- NULL
eboplus$Probes <- NULL
vebcon$Probes <- NULL

min(b2)
b1 <- rescale(as.matrix(b1))
b2 <- rescale(as.matrix(b2))
eboplus <- rescale(as.matrix(eboplus))
vebcon <- rescale(as.matrix(vebcon))
#View(eboplus)

b1 <- as.data.frame(b1)
b2 <- as.data.frame(b2)
eboplus <- as.data.frame(eboplus)
vebcon <- as.data.frame(vebcon)

b1$Probes <- rownames(b1)
b2$Probes <- rownames(b2)
eboplus$Probes <- rownames(eboplus)
vebcon$Probes <- rownames(vebcon)

dim(eboplus)
b1 <- b1[,c(21,1:20)]
b2 <- b2[,c(53,1:52)]
eboplus <- eboplus[,c(65,1:64)]
vebcon <- vebcon[,c(29, 1:28)]
View(b2)

#####
write.table(b1, "Intermed/Scaled_Ebolas_Apr2022/scaled_ebovac_b1_d3.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(b2, "Intermed/Scaled_Ebolas_Apr2022/scaled_ebovac_b2_d3.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(eboplus, "Intermed/Scaled_Ebolas_Apr2022/scaled_eboplus_d3.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(vebcon, "Intermed/Scaled_Ebolas_Apr2022/scaled_vebcon_d3.txt", quote = FALSE, row.names = FALSE, sep = '\t')


##### Dia 7 --------------------------------------------------------------------
geneva_counts <- read.csv2("Data/Counts_Geneva.csv")
rownames(geneva_counts) <- geneva_counts$X
geneva_pheno <- read.csv2("Data/Pheno_Geneva.csv")
geneva_pheno$Sample <- paste0("X",geneva_pheno$Sample)
geneva_pheno <- geneva_pheno[geneva_pheno$Treatment == "V",]
#View(geneva_pheno)

geneva_pheno_b1 <- geneva_pheno[geneva_pheno$Library.Batch == 1,]
geneva_pheno_b2 <- geneva_pheno[geneva_pheno$Library.Batch == 2,]

inter_b1 <- intersect(geneva_pheno_b1$Sample, colnames(geneva_counts))
inter_b2 <- intersect(geneva_pheno_b2$Sample, colnames(geneva_counts))

geneva_counts_b1 <- geneva_counts[,inter_b1]
geneva_counts_b2 <- geneva_counts[,inter_b2]
#View(geneva_counts_b2)

## Select D1  
geneva_pheno_b1 <- geneva_pheno_b1[geneva_pheno_b1$Group.Day %in% c(0,7),]
geneva_pheno_b2 <- geneva_pheno_b2[geneva_pheno_b2$Group.Day %in% c(0,7),]
#View(geneva_pheno_b1)

inter_b1 <- intersect(geneva_pheno_b1$Sample, colnames(geneva_counts_b1))
inter_b2 <- intersect(geneva_pheno_b2$Sample, colnames(geneva_counts_b2))

geneva_counts_b1 <- geneva_counts_b1[,inter_b1]
geneva_counts_b2 <- geneva_counts_b2[,inter_b2]
#View(geneva_counts_b2)

## Filtering and normalization
# B1
y_ebovac <- DGEList(counts = geneva_counts_b1, genes = row.names(geneva_counts_b1), group= geneva_pheno_b1$Group.Day)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) #  12353    37

y.1 <- calcNormFactors(y_ebovac)
b1 <- log2(cpm(y.1)+1)
b1 <- as.data.frame(b1)
b1$Probes <- rownames(b1)
dim(b1)
b1 <- b1[,c(38, 1:37)]
min(b1[,2:length(b1)])
#View(b1)

geneva_pheno_b1 <- geneva_pheno_b1[,c(1,4)]
colnames(geneva_pheno_b1) <- c("Probes", "Class")
#View(geneva_pheno_b1)
geneva_pheno_b1$Class <- paste0("D", geneva_pheno_b1$Class)

write.table(geneva_pheno_b1, "Intermed/Scaled_Ebolas_Apr2022/pheno_ebovac_b1_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')


# B2
y_ebovac_2 <- DGEList(counts = geneva_counts_b2, genes = row.names(geneva_counts_b2), group= geneva_pheno_b2$Group.Day)

keep2 <- rowSums(cpm(y_ebovac_2)>1)>=10
y_ebovac_2 <- y_ebovac_2[keep2,]
dim(y_ebovac_2) # 12713    49

y.2 <- calcNormFactors(y_ebovac_2)
b2 <- log2(cpm(y.2)+1)
b2 <- as.data.frame(b2)
b2$Probes <- rownames(b2)
dim(b2)
b2 <- b2[,c(50, 1:49)]
View(b2)

geneva_pheno_b2 <- geneva_pheno_b2[,c(1,4)]
colnames(geneva_pheno_b2) <- c("Probes", "Class")
geneva_pheno_b2$Class <- paste0("D", geneva_pheno_b2$Class)
#View(geneva_pheno_b2)

write.table(geneva_pheno_b2, "Intermed/Scaled_Ebolas_Apr2022/pheno_ebovac_b2_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')


##### EBOPLUS ####
eboplus <- read.delim("Intermed/data_BioFeats_USA/final_data/norm_counts_d7.txt")
#View(eboplus)


##### VEBCON ####
vebcon <- read.delim("Intermed/BioFeatS_EBOP_EBOV_VEB/D7/norm_vebcon_d7.txt")
#View(vebcon)

#inter <- Reduce(intersect, list(b1$Probes, b2$Probes, eboplus$Probes, vebcon$Probes))

#b1 <- b1[b1$Probes %in% inter,] # 9904 genes
#b2 <- b2[b2$Probes %in% inter,] # 9904 genes
#eboplus <- eboplus[eboplus$Probes %in% inter,] # 9904 genes
#vebcon <- vebcon[vebcon$Probes %in% inter,] # 9904 genes

## Prepare to scale
rownames(eboplus) <- eboplus$Probes
rownames(vebcon) <- vebcon$Probes

# namesb1 <- rownames(b1)
# namesb2 <- rownames(b2)
# nameseboplus <- rownames(eboplus)
# namesvebcon <- rownames(vebcon)

##### Pegar so genes importantes ####
genes <- read.csv("Genes_bio.csv")
#View(genes)

b1 <- b1[which(b1$Probes %in% genes$x),]
#View(b1)
b2 <- b2[which(b2$Probes %in% genes$x),]
eboplus <- eboplus[which(eboplus$Probes %in% genes$x),]
vebcon <- vebcon[which(vebcon$Probes %in% genes$x),]
#View(vebcon)

##### Scale all tables ####
b1$Probes <- NULL
b2$Probes <- NULL
eboplus$Probes <- NULL
vebcon$Probes <- NULL

min(b2)
b1 <- rescale(as.matrix(b1))
b2 <- rescale(as.matrix(b2))
eboplus <- rescale(as.matrix(eboplus))
vebcon <- rescale(as.matrix(vebcon))
#View(eboplus)

b1 <- as.data.frame(b1)
b2 <- as.data.frame(b2)
eboplus <- as.data.frame(eboplus)
vebcon <- as.data.frame(vebcon)

b1$Probes <- rownames(b1)
b2$Probes <- rownames(b2)
eboplus$Probes <- rownames(eboplus)
vebcon$Probes <- rownames(vebcon)

dim(vebcon)
b1 <- b1[,c(38,1:37)]
b2 <- b2[,c(50,1:49)]
eboplus <- eboplus[,c(65,1:64)]
vebcon <- vebcon[,c(29, 1:28)]
View(b2)

#####
write.table(b1, "Intermed/Scaled_Ebolas_Apr2022/scaled_ebovac_b1_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(b2, "Intermed/Scaled_Ebolas_Apr2022/scaled_ebovac_b2_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(eboplus, "Intermed/Scaled_Ebolas_Apr2022/scaled_eboplus_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(vebcon, "Intermed/Scaled_Ebolas_Apr2022/scaled_vebcon_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')


#






















#####

pheno_eboplus_D3 <- read.csv("Intermed/data_BioFeats_USA/phenodata/pheno_USA_d3.csv")
View(pheno_eboplus_D3)
pheno_eboplus_D3$Sampling.Day <- gsub(0, "D0", pheno_eboplus_D3$Sampling.Day)
pheno_eboplus_D3$Sampling.Day <- gsub(3, "D3", pheno_eboplus_D3$Sampling.Day)
pheno_eboplus_D3$X <- NULL
colnames(pheno_eboplus_D3) <- c("Probes", "Class")
write.table(pheno_eboplus_D3, "Intermed/Scaled_Ebolas_Apr2022/pheno_eboplus_d3.txt", quote = FALSE, row.names = FALSE, sep = '\t')

### testing data Scaled Ebolas April

pheno <- read.delim("Intermed/Scaled_Ebolas_Apr2022/D7/pheno_eboplus_d7.txt")
counts <- read.delim("Intermed/Scaled_Ebolas_Apr2022/D7/scaled_eboplus_d7.txt")
View(counts);View(pheno)

##### Data no batches ----------------------------------------------------------
##### EBOVAC ####
geneva_counts <- read.csv2("Data/Counts_Geneva.csv")
rownames(geneva_counts) <- geneva_counts$X
geneva_pheno <- read.csv2("Data/Pheno_Geneva.csv")
geneva_pheno$Sample <- paste0("X",geneva_pheno$Sample)
geneva_pheno <- geneva_pheno[geneva_pheno$Treatment == "V",]
#View(geneva_pheno)

geneva_pheno <- geneva_pheno[geneva_pheno$Library.Batch %in% c(1,2),]

## Select D1  
geneva_pheno <- geneva_pheno[geneva_pheno$Group.Day %in% c(0,3),]                                       

inter <- intersect(geneva_pheno$Sample, colnames(geneva_counts))

geneva_counts <- geneva_counts[,inter] # 86 samples
#View(geneva_counts_b2)

## Filtering and normalization
y_ebovac <- DGEList(counts = geneva_counts, genes = row.names(geneva_counts), group= geneva_pheno$Group.Day)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 12837    89

y.1 <- calcNormFactors(y_ebovac)
b <- log2(cpm(y.1)+1)
b <- as.data.frame(b)
b$Probes <- rownames(b)
View(b)
b <- b[,c(73, 1:72)]

geneva_pheno <- geneva_pheno[,c(1,4)]
colnames(geneva_pheno) <- c("Probes", "Class")
View(geneva_pheno)
geneva_pheno$Class <- paste0("D",geneva_pheno$Class)

write.table(geneva_pheno, "Intermed/Scaled_Ebolas_Apr2022/pheno_ebovac_d3.txt", quote = FALSE, row.names = FALSE, sep = '\t')


genes <- read.csv("Genes_bio.csv")
b <- b[which(b$Probes %in% genes$x),]
b$Probes <- NULL
b <- rescale(as.matrix(b))
b <- as.data.frame(b)
b$Probes <- rownames(b)

dim(b)
b <- b[,c(87,1:86)]
View(b)

write.table(b, "Intermed/Scaled_Ebolas_Apr2022/scaled_ebovac_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')


#########################################################################

##### EBOPLUS ####
eboplus <- read.delim("Intermed/data_BioFeats_USA/final_data/norm_counts_d1.txt")
#View(eboplus)


##### VEBCON ####
vebcon <- read.delim("Intermed/BioFeatS_EBOP_EBOV_VEB/D1/vebcon_norm_counts_d1.txt")
#View(vebcon)

inter <- Reduce(intersect, list(b$Probes, eboplus$Probes, vebcon$Probes))

b <- b[b$Probes %in% inter,] # 10089 genes
eboplus <- eboplus[eboplus$Probes %in% inter,] # 10089 genes
vebcon <- vebcon[vebcon$Probes %in% inter,] # 10089 genes

## Prepare to scale
rownames(eboplus) <- eboplus$Probes
rownames(vebcon) <- vebcon$Probes

# namesb1 <- rownames(b1)
# namesb2 <- rownames(b2)
# nameseboplus <- rownames(eboplus)
# namesvebcon <- rownames(vebcon)


##### Scale all tables ####
b$Probes <- NULL
eboplus$Probes <- NULL
vebcon$Probes <- NULL

b <- rescale(as.matrix(b))
eboplus <- rescale(as.matrix(eboplus))
vebcon <- rescale(as.matrix(vebcon))
View(eboplus)

b <- as.data.frame(b)
eboplus <- as.data.frame(eboplus)
vebcon <- as.data.frame(vebcon)

b$Probes <- rownames(b)
eboplus$Probes <- rownames(eboplus)
vebcon$Probes <- rownames(vebcon)

dim(b)
b <- b[,c(90,1:89)]
eboplus <- eboplus[,c(64,1:63)]
vebcon <- vebcon[,c(27, 1:26)]
View(vebcon)

#####

write.table(b, "Intermed/BioFeatS_scaled_ebola/All/scaled_ebovac_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(eboplus, "Intermed/BioFeatS_scaled_ebola/All/scaled_eboplus_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(vebcon, "Intermed/BioFeatS_scaled_ebola/All/scaled_vebcon_d1.txt", quote = FALSE, row.names = FALSE, sep = '\t')


#####

pheno_d1 <- read.("Intermed/data_BioFeatS_Geneva/final_data/pheno_Geneva_d1.txt")
View(pheno_d1)




