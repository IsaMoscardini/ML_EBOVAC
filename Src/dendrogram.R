### HEATMAP SEM HEATMAP #####
#install.packages("ggdendro")
#install.packages('dendextend')

# load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  library(gplots)
  library(data.table)
  library(clusterProfiler)
  library(FactoMineR)
  library(factoextra)
  library(ggdendro)
  library(dplyr)
  library(dendextend)
  library(tidyverse)
})


#######

## set phenotype
phenoFile <- "Data/Outcomes.txt"
# read phenotype
phenoData <- read.table(phenoFile, header=T, row.names=1, sep="\t")
phenoData$AdverseEvent_Class <- NULL
phenoData$name <- rownames(phenoData)
colnames(phenoData) <- c( "class", "name")
phenoData <- phenoData[complete.cases(phenoData$class),]
#View(phenoData)


#####
## set data file
dataFile <- "Data/D1/Xv0_Log2_Olink_D1vD0.txt"
# read data
df <- as.matrix(read.table(dataFile, header=T, row.names=1, sep="\t"))
View(df)

## selected features
features <- read.csv("Results/Selected_features_higher_0/Ab/Xv0_Log2_Olink_D1vD0.csv")
View(features)
features <- as.character(features$Probes)

# Filter data
df <- df[features, intersect(rownames(phenoData),colnames(df))]
View(df)

methList <- c("euclidean","maximum","canberra","manhattan","minkowski")
sizeList <- 2:length(features)

df <- t(df)

pheno <- phenoData[rownames(df),]
my_colors <- ifelse(pheno$class == "Persistent", "coral1", "lightseagreen")


#######
# functions
most_variant_hc <- function(exprs, labels, n=s, method=method, title = outname) {
   
   exprs[,1:n] %>%
    dist(method = method) %>% 
    hclust() %>% 
    as.dendrogram() -> dend
  plot(dend, main = outname)
  colored_bars(colors = my_colors, dend = dend, rowLabels = "AB response")
  
  # rect = TRUE, rect_fill = TRUE, rect_border = colors
  legend("topright", legend=levels(as.factor(df_pheno$group)),cex = 0.9, fill = c("coral1", "darkturquoise")) # inset = c(-1,0)) "bottomleft"
}


##### run analysis #####

pdf("Results/Selected_features_higher_0/AB_Olink_D1.pdf", width=8, height=4)
  #par(mfrow = c(5,1))
for (s in sizeList){
  for (method in methList){
    outname <- paste("Top", s, "cytokines", "Method", method)
    print(most_variant_hc(exprs = df, labels =  pheno$name, n = s, method = method))
    }   
  }
dev.off()

bv