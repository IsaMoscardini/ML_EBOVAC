rm(list = ls())
options(stringsAsFactors = F)
library(mixOmics)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggpubr)

## set phenotype
phenoFile <- "Data/Outcomes.txt"
# read phenotype
phenoData <- read.table(phenoFile, header=T, row.names=1, sep="\t")
phenoData$AdverseEvent_Class <- NULL
phenoData$name <- rownames(phenoData)
colnames(phenoData) <- c( "class", "name")
phenoData <- phenoData[complete.cases(phenoData$class),]
View(phenoData)


#####
## set data file
dataFile <- "Data/D1/RNAseq_logCPM_D1vD0.txt"
# read data
df <- as.matrix(read.table(dataFile, header=T, row.names=1, sep="\t"))
View(df)
#df <- as.data.frame(df)
#names <- rownames(df)
#df <- apply(df, 2, as.numeric)
#names <- make.names(names)
#rownames(df) <- names


## selected features
features <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D1vD0.csv")
View(features)
probes <- features$Probes

#probes <- make.names(features$Probes)
#View(probes)

# Filter data
df <- df[probes, intersect(rownames(phenoData),colnames(df))]
View(df)

sizeList <- 2:length(features)

df <- t(df)

pheno <- phenoData[rownames(df),]
my_colors <- ifelse(pheno$class == "Persistent", "coral1", "lightseagreen")
AB_response <- pheno$class
#View(df)

###### FUNCTIONS ###### 
make_PCA <- function(exprs, labels, n=s, outname) {
  exprs <- exprs[,1:n]
  pc <- as.data.frame(pca(exprs)$x)
  p <- ggplot(pc, aes(x = PC1, y = PC2, col= AB_response)) + geom_point(size=2) + 
    theme_minimal() + theme_bw() + labs(title = outname)  + 
    scale_color_manual(values = c("coral1", "lightseagreen")) + theme(legend.position = "top")
  plot(p)
}

###### Analysis #####
pdf("Results/Selected_features_higher_0/AB_miRNA_D7_PCA.pdf", width=6, height=5)
for (s in sizeList){
    outname <- paste("Top", s, "miRNAs")
    make_PCA(exprs = df, labels =  pheno$name, n = s, outname)
  }   
dev.off()


#

###########

make_BOX <- function(exprs, labels, n=s, outname) {
  exprs <- as.data.frame(t(exprs))
  exprs$gene <- rownames(exprs)
  mel <- melt(exprs)
  mel$group <- pheno$class[match(mel$variable, pheno$name)]
  AB_response <- mel$group
  p <- ggplot(mel, aes(x = reorder(mel$gene, -value), y = mel$value, col = AB_response)) + 
    stat_compare_means(method = "wilcox.test", paired = FALSE) + geom_boxplot() + 
    geom_jitter() + labs(title = "Antibody response RNAseq D1") + xlab("Genes") + ylab("Expression value")
  plot(p)
}

pdf("Results/Selected_features_higher_0/AB_RNAseq_D1_BOX.pdf", width = 90, height = 30)  
make_BOX(df, exprs = df, labels =  pheno$name)
dev.off()
