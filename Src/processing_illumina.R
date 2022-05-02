# Processing Illumina for microarray studies
# Amanda Vasconcelos

setwd("~/CSBL/projects/doutorado/studies/microarray/GSE54514/")


# libraries
library(GEOquery)
library(arrayQualityMetrics)
library(tidyverse)
library(Biobase)
library(mdp)
library(limma)


# Create directories to save files -----------------------------------------------------------

# Create directories to save AQM before and after normalization 
dir.create("data/")
dir.create("intermediate/")
dir.create("intermediate/AQM_before_normalization")
dir.create("intermediate/AQM_after_normalization")
dir.create("intermediate/AQM_free_outliers")
dir.create("intermediate/MDP")

gseID <- "GSE54514"
gpl <- "GPL6947"
outDir <- "./data/"

# Save in an object the path where the output of AQM before pre processing will be saved
aqm_before_normalization <- paste0("~/CSBL/projects/doutorado/studies/microarray/", gseID, "/intermediate/", "AQM_before_normalization")
aqm_after_normalization <- paste0("~/CSBL/projects/doutorado/studies/microarray/", gseID, "/intermediate/", "AQM_after_normalization")
aqm_free_outliers <- paste0("~/CSBL/projects/doutorado/studies/microarray/", gseID, "/intermediate/", "AQM_free_outliers")
mdp_before_outliers <- paste0("~/CSBL/projects/doutorado/studies/microarray/", gseID, "/intermediate/", "MDP/", "before_outliers_removal")
mdp_after_outliers <- paste0("~/CSBL/projects/doutorado/studies/microarray/", gseID, "/intermediate/", "MDP/", "after_outliers_removal")


# Download data -----------------------------------------------------------

gse <- getGEO(gseID)
names(gse)

getGEOSuppFiles(gseID) # raw data, CEL files
expr <- exprs(gse[[1]]) # normalized expression by the author
pd <- pData(gse[[1]]) # phenodata
plat <- gse[[1]]@featureData@data # platform annotation


# saving the tables
write.table(expr, paste0(outDir, gseID, "_","expression_normalized_author.tsv"), sep = "\t", row.names = TRUE, col.names = NA)
write.table(pd, paste0(outDir, gseID, "_", "phenodata_original.tsv"), sep = "\t", row.names = TRUE, col.names = NA)
write.table(plat, paste0(outDir, gseID, "_", "platform_annotation.tsv"), sep = "\t", row.names = TRUE, col.names = NA)


# Load the raw expression matrix
# EListRaw format equivalent to Affybatch
# https://rdrr.io/bioc/limma/man/EList.html

# in case where there're many files, run this first
raw_files <- list.files(path = "data/rawdata/", pattern = ".txt", ignore.case = TRUE, recursive = TRUE, full.names = TRUE)

# exp_df <- read.table("data/rawdata/GSE69528_non_normalized.txt", header = T, sep = "\t") # check the table
rawdata <- read.ilmn(files = paste0("data/rawdata/", gseID, "_non-normalized.txt"), ctrlfiles = NULL, probeid = "ID_REF",
                     expr = "X", other.columns = "Detection.Pval")
# rawdata <- read.ilmn(files = raw_files, ctrlfiles = NULL, probeid = "PROBE_ID", 
#                                           expr = "AVG_Signal", other.columns = "Detection Pval")


# Get the raw expression matrix
raw_exprs <- rawdata$E
min(raw_exprs)
max(raw_exprs)
raw_exprs_df <- as.data.frame(raw_exprs)

# Save the raw expression matrix
write.table(raw_exprs_df, paste0(outDir, gseID, "_","raw_expression_no_norm.tsv"), sep = "\t", row.names = TRUE, col.names = NA)


# Quality control before normalization -----------------------------------------------------------
raw_expset <- ExpressionSet(assayData = raw_exprs)

# Run AQM using as input an ExpressionSet that contains only the raw expression matrix
arrayQualityMetrics(expressionset = raw_expset,
                    outdir = aqm_before_normalization,
                    force = TRUE,
                    do.logtransform = TRUE)


# Normalization -----------------------------------------------------------
expr_norm <- neqc(rawdata, detection.p = "Detection.Pval")

# Quality control after normalization -----------------------------------------------------------
norm_expset <- ExpressionSet(assayData = expr_norm$E)

arrayQualityMetrics(expressionset = norm_expset,
                    outdir = aqm_after_normalization,
                    force = TRUE,
                    do.logtransform = FALSE)


# Outlier removal and renormalization -----------------------------------------------------------

# Creating a vector with the outlier samples (see the table in the index)
# outlier_samples <- c("128115_GC008233485DR_5388500002_L")
# 
# rawdata_noOutliers <- rawdata[, !(colnames(rawdata$E) %in% outliers_samples)]


# Get the final expression table
expr_table <- expr_norm$E
expr_table <- as.data.frame(expr_table)
expr_table <- expr_table[, pd$description.1]
table(colnames(expr_table) == pd$description.1)
colnames(expr_table) <- pd$geo_accession


# Organize the phenodata -----------------------------------------------------------

# remove the outliers from the phenodata
# outlier_samples <- c("GSM4223798")
# pd <- pd[!pd$geo_accession %in% outlier_samples, ]


# organizing the phenodata (colnames with "characteristics")
col_pd_characteristics <- grep("characteristics_ch", colnames(pd)) # identifica quais colunas tem "characteristics"
line_ch  <- as.character(as.matrix(pd[1, col_pd_characteristics])) #pega a info da linha que contém "characteristics" para extrair as infos que serão o header
characts <- unlist(lapply(strsplit(line_ch, ": "), `[[`, 1)) #separa a info da linha pelos dois pontos e mantém o primeiro nome. 
for (i in 1:length(col_pd_characteristics)) # subsitui o header "characteristics" pelo primeiro nome da linha
{
  colnames(pd)[col_pd_characteristics[i]] <- characts[i]
  pd[, col_pd_characteristics[i]] <-
    gsub(paste(colnames(pd)[col_pd_characteristics[i]], ": ", sep = ""), '', pd[, col_pd_characteristics[i]])
}

# Personalized filters
names(pd)
pd <- pd[, -c(3:7, 12:ncol(pd))]
colnames(pd)[4] <- "Organism"
colnames(pd)[6] <- "Class"
colnames(pd)[5] <- "sex"
colnames(pd)[2] <- "Sample"
colnames(pd)[3] <- "Groups"
pd$age_group <- "adult"

pd <- pd[!pd$Groups == "se/ARDS Day 0", ]
pd$Class[pd$Groups == "Sepsis Day 0"] <- "sepsis"

expr_table <- expr_table[pd$Sample]

# write conditions based on a string from another column
# pd$Class <- ifelse(grepl("Control", pd$Groups, ignore.case = TRUE), "control",
#                    ifelse(grepl("SIRS", pd$Groups, ignore.case = TRUE), "SIRS", "default"))



# ------------------------------------

# Load the annotation table
annot_symbols <- "C:/Users/Amanda/Documents/CSBL/projects/doutorado/useful_files/annotation_long_symbols.tsv"
annot <- read.delim(file = annot_symbols, sep = "\t", header = TRUE)

# Filtering the annotation table for the respective platform 
annot <- annot[annot$Platform == gpl, ]
table(annot$Platform)

# Putting the probe name as a column in the normalized matrix
expr_table["Probe"] <- rownames(expr_table)

# Organize the phenodata
expr_table <- expr_table %>% 
  dplyr::select(Probe, everything())

# Join the tables to get the gene symbol of the probes
expr_table2 <- expr_table %>%
  dplyr::left_join(annot, by = "Probe") %>% 
  dplyr::select(colnames(expr_table), hgnc_symbol)

# Removing NA and empty rows 
expr_table2 <- expr_table2[!is.na(expr_table2$hgnc_symbol), ]
expr_table2 <- expr_table2[!expr_table2$hgnc_symbol == "", ]

# Detect and remove probes that bind in more than one gene, usually contain "///".
expr_table2 <- expr_table2 %>%
  dplyr::mutate(remove_genes = str_detect(hgnc_symbol, "///"))
#
table(expr_table2$remove_genes) # this dataset doesn't have probe for more than one gene

expr_table2 <- expr_table2 %>%
  dplyr::filter(remove_genes == FALSE) %>%
  dplyr::select(-remove_genes)

# Reorder the colums
expr_table2 <- expr_table2 %>% 
  dplyr::select(Probe, hgnc_symbol, everything())


# Collapse rows -----------------------------------------------------------

# Function to collapse rows 
collapse.rows <- function(expr, probe.col, gene.col, data.table=F, method=c("maxMean", "minMean", "colMean", "colMedian")){
  if(length(grep('data.table', installed.packages())) == 0){
    install.packages('data.table')
    require(data.table)
  }else if(length(grep('data.table', search())) == 0){
    suppressPackageStartupMessages(require(data.table))
  }
  
  if (probe.col == "rownames"){
    expr <- data.table(expr, keep.rownames=T)
    setnames(expr, "rn", "rownames")
  }else{
    expr <- data.table(expr)
  }
  
  if(method=="maxMean" | method=="minMean"){ 
    expr[, rowmean := rowMeans(.SD[, !c(probe.col, gene.col), with=F])]
    if(method=="maxMean"){
      res <- expr[order(rowmean, decreasing=T)][, .SD[1], by=gene.col][, rowmean:=NULL]
    }
    else if(method=="minMean"){
      res <- expr[order(rowmean, decreasing=T)][, .SD[.N], by=gene.col][, rowmean:=NULL]
    }
  }
  else if(method=="colMean"){
    res <- expr[, lapply(.SD[, !c(probe.col), with=F], mean), by=gene.col]
  }
  else if(method=="colMedian"){
    res <- expr[, lapply(.SD[, !c(probe.col), with=F], median), by=gene.col]   
  }
  else stop("method must be 'maxMean', 'minMean', 'colMean' or 'colMedian'\n")
  
  if(!data.table){
    return(data.frame(res))
  }else{ return(res[]) }
}


# Collapse probes using the maxMean method
expr_final <- collapse.rows(expr = expr_table2, probe.col = "Probe", gene.col = "hgnc_symbol", method = "maxMean")

rownames(expr_final) <- expr_final$hgnc_symbol
expr_final$hgnc_symbol <- NULL
expr_final$Probe <- NULL



# MDP -----------------------------------------------------------
# github: https://github.com/csbl-usp/mdp/blob/master/R/mdp.R
# vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/mdp/inst/doc/my-vignette.html

mdp_results <- mdp(data = expr_final, pdata = pd, control_lab = "healthy", mdp_before_outliers, fraction_genes = 0.25, measure = "median", std = 2)

perturbedgenes <- mdp_results[["sample_scores"]][["perturbedgenes"]]
outliers_mdp <- perturbedgenes$Sample[perturbedgenes$outlier == 1]

pd <- pd[!pd$Sample %in% outliers_mdp, ]
expr_final <- expr_final[, !names(expr_final) %in% outliers_mdp]

mdp_results <- mdp(data = expr_final, pdata = pd, control_lab = "healthy", mdp_after_outliers, fraction_genes = 0.25, measure = "median", std = 2)



# Save final tables -----------------------------------------------------------

Outdir <- paste0("~/CSBL/projects/doutorado/studies/microarray/", gseID)
setwd(Outdir)

# Save tables
write.table(expr_table2, paste0("data/", gseID, "_", "expr_norm_noCollapsed.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(expr_final, paste0("data/", gseID, "_", "expr_norm_collapsed_final.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(pd, paste0("data/", gseID, "_", "phenodata.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
# write.table(annot, paste0("data/", gseID, "_", "annotation_biomaRt.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

