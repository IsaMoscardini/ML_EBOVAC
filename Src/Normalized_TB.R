#####
library(ggplot2)
library(reshape2)
library(dplyr)

######
GSE19439_counts <- read.delim("Normalized_TB_studies/GSE19439_expression_normalized_author.tsv")
GSE19439_pheno <- read.delim("Normalized_TB_studies/GSE19439_phenodata_original.tsv")
rownames(GSE19439_counts) <- GSE19439_counts$X
#View(GSE19439_counts)

GSE19439_pheno <- GSE19439_pheno[,c(1,14)]
GSE19439_pheno <- GSE19439_pheno[GSE19439_pheno$characteristics_ch1.3 != "illness: PTB",]
GSE19439_pheno$characteristics_ch1.3 <- gsub("illness: ", "", GSE19439_pheno$characteristics_ch1.3)
GSE19439_pheno$characteristics_ch1.3 <- gsub(" \\(BCG\\+)", "", GSE19439_pheno$characteristics_ch1.3)
GSE19439_pheno$characteristics_ch1.3 <- gsub(" \\(BCG\\-)", "", GSE19439_pheno$characteristics_ch1.3)
GSE19439_pheno$characteristics_ch1.3 <- gsub(" ", "", GSE19439_pheno$characteristics_ch1.3)
colnames(GSE19439_pheno) <- c("Probes", "Class")
#View(GSE19439_pheno)

inter <- intersect(colnames(GSE19439_counts), GSE19439_pheno$Probes)
GSE19439_counts <- GSE19439_counts[,inter]
GSE19439_pheno <- GSE19439_pheno[GSE19439_pheno$Probes %in% inter, ]

identical(colnames(GSE19439_counts), GSE19439_pheno$Probes) # TRUE

feat <- read.delim("GSE19439/data/GSE19439_platform_annotation.tsv")
#View(feat)

GSE19439_counts$Symbol <- feat$Symbol[match(rownames(GSE19439_counts), feat$ID)]
GSE19439_counts$SUM <- rowSums(GSE19439_counts[,1:29])
GSE19439_counts <- GSE19439_counts[order(abs(GSE19439_counts$SUM), decreasing = TRUE),]
GSE19439_counts <- GSE19439_counts[!duplicated(GSE19439_counts$Symbol), ] 
dim(GSE19439_counts) # 25154    31
GSE19439_counts$SUM <- NULL

GSE19439_counts <- GSE19439_counts[,c(30,1:29)]
colnames(GSE19439_counts)[1] <- "Probes"
#View(GSE19439_counts)
GSE19439_counts[,2:30] <- log2(GSE19439_counts[,2:30])
GSE19439_counts <- GSE19439_counts[complete.cases(GSE19439_counts),]


# boxplot(log2(GSE19439_counts[,2:30]))
# melty <- reshape2::melt(GSE19439_counts)

#write.table(GSE19439_counts, "Normalized_TB_studies/GSE19439_final_counts.txt", quote = FALSE, sep = '\t', row.names = FALSE)
#write.table(GSE19439_pheno, "Normalized_TB_studies/GSE19439_final_pheno.txt", quote = FALSE, sep = '\t', row.names = FALSE)



#####

GSE19444_counts <- read.delim("Normalized_TB_studies/GSE19444_expression_normalized_author.tsv")
GSE19444_pheno <- read.delim("Normalized_TB_studies/GSE19444_phenodata_original.tsv")
rownames(GSE19444_counts) <- GSE19444_counts$X
#View(GSE19444_pheno)
#View(GSE19444_counts)

GSE19444_pheno <- GSE19444_pheno[,c(1,14)]
GSE19444_pheno <- GSE19444_pheno[GSE19444_pheno$characteristics_ch1.3 != "illness: PTB",]
GSE19444_pheno$characteristics_ch1.3 <- gsub("illness: ", "", GSE19444_pheno$characteristics_ch1.3)
GSE19444_pheno$characteristics_ch1.3 <- gsub(" \\(BCG\\+)", "", GSE19444_pheno$characteristics_ch1.3)
GSE19444_pheno$characteristics_ch1.3 <- gsub(" \\(BCG\\-)", "", GSE19444_pheno$characteristics_ch1.3)
GSE19444_pheno$characteristics_ch1.3 <- gsub(" ", "", GSE19444_pheno$characteristics_ch1.3)
colnames(GSE19444_pheno) <- c("Probes", "Class")
#View(GSE19444_pheno)

inter <- intersect(colnames(GSE19444_counts), GSE19444_pheno$Probes)
GSE19444_counts <- GSE19444_counts[,inter]
GSE19444_pheno <- GSE19444_pheno[GSE19444_pheno$Probes %in% inter, ]

identical(colnames(GSE19444_counts), GSE19444_pheno$Probes) # TRUE

feat <- read.delim("GSE19444/data/GSE19444_platform_annotation.tsv")
#View(feat)

GSE19444_counts$Symbol <- feat$Symbol[match(rownames(GSE19444_counts), feat$ID)]
GSE19444_counts$SUM <- rowSums(GSE19444_counts[,1:29])
GSE19444_counts <- GSE19444_counts[order(abs(GSE19444_counts$SUM), decreasing = TRUE),]
GSE19444_counts <- GSE19444_counts[!duplicated(GSE19444_counts$Symbol), ] 
dim(GSE19444_counts) # 25160    35
GSE19444_counts$SUM <- NULL

GSE19444_counts <- GSE19444_counts[,c(34,1:33)]
colnames(GSE19444_counts)[1] <- "Probes"
GSE19444_counts[,2:34] <- log2(GSE19444_counts[,2:34])
GSE19444_counts <- GSE19444_counts[complete.cases(GSE19444_counts),]
#View(GSE19444_counts)

inter <- intersect(GSE19439_counts$Probes, GSE19444_counts$Probes)
GSE19439_counts <- GSE19439_counts[GSE19439_counts$Probes %in% inter,]
GSE19444_counts <- GSE19444_counts[GSE19444_counts$Probes %in% inter,] 
identical(GSE19439_counts$Probes, GSE19444_counts$Probes)

write.table(GSE19439_counts, "Normalized_TB_studies/GSE19439_final_counts.txt", quote = FALSE, sep = '\t', row.names = FALSE)
write.table(GSE19444_counts, "Normalized_TB_studies/GSE19444_final_counts.txt", quote = FALSE, sep = '\t', row.names = FALSE)
#write.table(GSE19444_pheno, "Normalized_TB_studies/GSE19444_final_pheno.txt", quote = FALSE, sep = '\t', row.names = FALSE)

boxplot(log2(GSE19444_counts[,2:34]))


##















##### GSE19439 #####

GSE19439 <- read.delim("Normalized_TB_studies/GSE19439_expression_normalized_author.tsv")
#View(GSE19439)
plat439 <- read.delim("Normalized_TB_studies/GSE19439_platform_annotation.tsv")
pheno439 <- read.delim("Normalized_TB_studies/GSE19439_phenodata_original.tsv")
View(pheno439)

GSE19439$Symbol <- plat439$Symbol[match(GSE19439$X, plat439$X)]
GSE19439$X <- NULL
dim(GSE19439) # 48791    43
GSE19439 <- GSE19439[,c(43, 1:42)]
GSE19439$Sum <- rowSums(GSE19439[,2:42])
View(GSE19439)
GSE19439 <- GSE19439[order(abs(GSE19439$Sum), decreasing = TRUE),]
GSE19439 <- GSE19439[!duplicated(GSE19439$Symbol),]
GSE19439 <- GSE19439[which(GSE19439$Symbol != ""),]
View(GSE19439)
