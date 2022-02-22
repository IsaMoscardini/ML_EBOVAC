# filtrar tabelas
genes <- read.csv("Data/Genes_bio (2).csv")
View(genes)

ff <- list.files(path="Data/RNAseq/", full.names=TRUE)
tabelas <- lapply(ff, read.delim) %>% setNames(basename(ff))
View(tabelas[[1]])

for (i in names(tabelas)){
  tabela = tabelas[[i]]
  tabela <- tabela[tabela$Probes %in% genes$x,]
  write.csv2(tabela, paste0('filtered_',i), row.names = TRUE)
}


tab <- read.csv2("filtered_RNAseq_logCPM_D3vD0.txt")
View(tab)

#####
rm(list = ls())
options(stringsAsFactors = F)
library(mixOmics)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggpubr)

## set phenotype
phenoFile <- "Intermed/data_BioFeats_USA/final_data/pheno_USA_d1.txt"
# read phenotype
phenoData <- read.table(phenoFile, header=T, row.names=1, sep="\t")
phenoData$AntibodyResponse_Class <- NULL
phenoData$name <- rownames(phenoData)
colnames(phenoData) <- c( "class", "name")
phenoData <- phenoData[complete.cases(phenoData$class),]
#View(phenoData)


#####
## set data file
dataFile <- "Intermed/data_BioFeats_USA/final_data/norm_counts_d1.txt"
# read data
df <- as.matrix(read.table(dataFile, header=T, row.names=1, sep="\t"))
View(df)
df <- as.data.frame(df)
names <- rownames(df)
df <- apply(df, 2, as.numeric)
#names <- make.names(names)
#names <- gsub("X.Tryptophan", "Tryptophan", names)
rownames(df) <- names


## selected features
features <- read.csv("Results/BioFeatS_USA/table_features_metrics.csv")
View(features)
features <- features[features$Sum > 0,]
probes <- features$Unnamed..0

#probes <- make.names(features$Probes)
View(probes)

# Filter data
df <- as.data.frame(df)
#df$Genes <- rownames(df)
df <- df[probes, intersect(rownames(phenoData),colnames(df))]
View(df)

sizeList <- 2:length(features)

df <- t(df)

pheno <- phenoData[rownames(df),]
my_colors <- ifelse(pheno$class == "D1", "coral1", "lightseagreen")
DAY <- pheno$class
View(df)

#
###### FUNCTION PCA ###### 
make_PCA <- function(exprs, labels, n=s, outname) {
  exprs <- exprs[,1:n]
  pc <- as.data.frame(pca(exprs)$x)
  p <- ggplot(pc, aes(x = PC1, y = PC2, col= AE)) + geom_point(size=2) + 
    theme_minimal() + theme_bw() + labs(title = outname)  + 
    scale_color_manual(values = c("coral1", "lightseagreen")) + theme(legend.position = "top")
  plot(p)
}

###### Analysis #####
pdf("Results/Selected_features_higher_0/AE_Olink_D7_PCA.pdf", width=6, height=5)
for (s in sizeList){
    outname <- paste("Top", s, "Proteins")
    make_PCA(exprs = df, labels =  pheno$name, n = s, outname)
  }   
dev.off()


#

###### FUNCTION BOX '####

make_BOX <- function(exprs, labels, n=s, outname) {
  exprs <- as.data.frame(t(exprs))
  exprs$gene <- rownames(exprs)
  mel <- melt(exprs)
  mel$group <- pheno$class[match(mel$variable, pheno$name)]
  AE_response <- mel$group
  p <- ggplot(mel, aes(x = reorder(mel$gene, -value), y = mel$value, col = AE_response)) + 
    stat_compare_means(method = "wilcox.test", paired = FALSE) + geom_boxplot() + 
    geom_jitter() + labs(title = "AE olink D7") + xlab("Proteins") + ylab("Value")
  plot(p)
}


###### Analysis ######

