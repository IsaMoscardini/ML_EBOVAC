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
})

########
# # Compute distances and hierarchical clustering
# dend <- df[,1:65] %>%
#         #scale %>%
#         dist(method = "euclidean") %>% # calculate a distance matrix
#         hclust(method = "ward.D2") %>% # Hierarchical clustering 
#         as.dendrogram %>%
#         set("labels")
#         
# plot(dend)
# colored_bars(colors = df$color, dend = dend, rowLabels = "group")



#######
# set data file
dataFile <- "Data/D1/RNAseq_logCPM_D1vD0.txt"
# read data
df <- as.matrix(read.table(dataFile, header=T, row.names=1, sep="\t"))
#head(df)

# set phenotype
phenoFile <- "Data/Outcomes.txt"
# read phenotype
phenoData <- read.table(phenoFile, header=T, row.names=1, sep="\t")
phenoData$AdverseEvent_Class <- NULL
phenoData$name <- rownames(phenoData)
colnames(phenoData) <- c( "class", "name")
phenoData <- phenoData[complete.cases(phenoData$class),]
#View(phenoData)

# selected features
features <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D1vD0.csv")
#View(features)
features <- features$Probes

# Filter data
df <- df[features,rownames(phenoData)]
#View(df)

methList <- c("euclidean","maximum","canberra","manhattan","minkowski")
sizeList <- 5:length(features)

df <- t(df)
df_pheno <- df
#View(df)

df_pheno <- as.data.frame(df_pheno)
df_pheno$group <- phenoData$class[match(rownames(df_pheno), phenoData$name)]
df_pheno$color <- ifelse(df_pheno$group == "Persistent", "coral1", "darkturquoise")
df_pheno$names <- rownames(df_pheno)
#View(df_pheno)


#######
# functions
most_variant_hc <- function(exprs, labels, n=s, method=method, title = outname) {
    
     # exprs <- exprs[order(rowMeans(exprs), decreasing=T), ] 
     # exprs <- exprs[order(apply(exprs, 1, var), decreasing=T), ]
     # exprs <- head(exprs, n=s)
    
    dend <- exprs %>%
        #scale %>%
        dist(method = method) %>% # calculate a distance matrix
        hclust(method = "ward.D2") %>% # Hierarchical clustering 
        as.dendrogram %>%
        set("labels")
    
    plot(dend, main = outname)
    colored_bars(colors = df_pheno$color, dend = dend, rowLabels = "group")
    
               
    # rect = TRUE, rect_fill = TRUE, rect_border = colors
    legend("topright", legend=levels(as.factor(df_pheno$group)),cex = 0.9, fill = c("coral1", "darkturquoise")) # inset = c(-1,0)) "bottomleft"
}


#

pdf("Results/Selected_features_higher_0/AB_RNAseq_D1_dend.pdf", width=8, height=12)
par(mfrow = c(5,1))
for (s in sizeList){
    for (method in methList){
        outname <- paste("Top", s, "genes", "Method", method)
        most_variant_hc(exprs = df, labels =  df_pheno$names, n = s, method = method)
    }
}
dev.off()



#######

dend <- df %>%
    #scale %>%
    dist(method = method) %>% # calculate a distance matrix
    hclust(method = "ward.D2") %>% # Hierarchical clustering 
    as.dendrogram %>%
    set("labels")

plot(dend)
colored_bars(colors = df_pheno$color, dend = dend, rowLabels = "group")

