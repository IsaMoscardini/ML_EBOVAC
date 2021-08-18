#!/usr/bin/Rscript


# load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(RColorBrewer)
    library(gplots)
    library(data.table)
})



# functions
most_variant_hc <- function(exprs, colors, labels, n=s, method=method, outname) {

    exprs <- exprs[order(rowMeans(exprs), decreasing=T), ] 
    exprs <- exprs[order(apply(exprs, 1, var), decreasing=T), ]
    exprs <- head(exprs, n=s)

    colors <- as.factor(colors)
    #print(colors)
    names <- levels(colors)
    #print(names)
    levels(colors) <- c("lightseagreen", "coral") 
    #levels(colors) <- brewer.pal(length(levels(colors)),"Set1")
    #print(levels(colors))
    leg_colors <- levels(colors)
    
    colors <- as.character(colors)
    #print(colors)
    
    heatmap.2(as.matrix(exprs), 
        col = colorRampPalette(c("blue", "white", "red"))(n=1000), 
        Rowv = T,
        Colv = T,
        distfun = function(x) dist(x, method=method),
        dendrogram = "column",
        density.info = "none",
        key = T, keysize=1.2,
        scale = "row",
        main = paste(outname),
        labCol = labels,
        labRow = rownames(exprs),
        cexRow = 0.9, cexCol = 0.9,
        ColSideColors = colors,
        trace = 'none',
        margins= c(9, 9),
        srtCol=90
    )
    legend("topright", legend=names,fill=leg_colors,  cex = 0.6) # inset = c(-1,0)) "bottomleft"
}


# main

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
#.View(features)
features <- features$Probes

# Filter data
df <- df[features,rownames(phenoData)]
#View(df)

#
methList <- c("euclidean","maximum","canberra","manhattan","minkowski")
sizeList <- 5:length(features)

#
pdf("Results/Selected_features_higher_0/AB_RNAseq_D1.pdf", width=10, height=8)
for (method in methList){
    for (s in sizeList){
        outname <- paste("Method", method, "Top", s, "genes")
        most_variant_hc(df, as.character(phenoData$class), phenoData$name, s, method, outname)
        }
}
dev.off()



























