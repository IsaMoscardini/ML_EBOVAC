#####
rm(list = ls()) 
options(stringsAsFactors = F) 

#####
GSE19439 <- read.delim("Normalized_TB_studies/GSE19439_expression_normalized_author.tsv")
feat19439 <- read.delim("Normalized_TB_studies/GSE19439_platform_annotation.tsv")
pheno19439 <- read.delim("Normalized_TB_studies/GSE19439_phenodata_original.tsv")
View(pheno19439)

GSE19439$Symbol <- feat19439$Symbol[match(GSE19439$X, feat19439$X)]
GSE19439 <- GSE19439[,c(44,1:43)]
GSE19439$Sum <- rowSums(GSE19439[,3:44])
GSE19439 <- GSE19439[order(GSE19439$Sum, decreasing = TRUE),]
View(GSE19439)


#####
GSE19442 <- read.delim("Normalized_TB_studies/GSE19442_expression_normalized_author.tsv")
feat19442 <- read.delim("Normalized_TB_studies/GSE19442_platform_annotation.tsv")
pheno19442 <- read.delim("Normalized_TB_studies/GSE19442_phenodata_original.tsv")
View(pheno19442)

GSE19442$Symbol <- feat19442$Symbol[match(GSE19442$X, feat19442$X)]
dim(GSE19442)
GSE19442 <- GSE19442[,c(53,1:52)]
View(GSE19442)
GSE19442$Sum <- rowSums(GSE19442[,3:53])
GSE19442 <- GSE19442[order(GSE19442$Sum, decreasing = TRUE),]
View(GSE19442)


#####





