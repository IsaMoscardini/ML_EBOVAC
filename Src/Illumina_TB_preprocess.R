#####
rm(list = ls()) 
options(stringsAsFactors = F) 
basedir <- "C:/Users/Isabelle/Documents/PROJECTS/Ebola_ML/GSE19439/" 
setwd(basedir) 
gse_id = 'GSE19439' 
#url_rawdata \<- '<https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE19439&format=file>'


pkgs <- c('matrixStats','GEOquery','impute','naniar','preprocessCore','mixOmics','Rtsne','arrayQualityMetrics','reshape2')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)
#source('helper_functions/eda.R')


#--- DIRECTORIES
out_dir = file.path(basedir,'data')
raw_out_dir = file.path(out_dir,'raw_data', gse_id)
dir.create(raw_out_dir, showWarnings = F,recursive = T)

eset <- getGEO(gse_id,destdir = raw_out_dir)
eset <- eset[[1]]

####### RAW DATA #####
# download raw expression data

# ----> RODAR APENAS UMA VEZ: vai baixar os dados brutos direto do GEO p/ a pasta data/raw_data
#download.file(url = url_rawdata, destfile = file.path(raw_out_dir,paste0(gse_id,'.txt.gz')))

files_exp <- list.files("data/raw_data/dizip/", pattern = "txt", full.names = TRUE)

expression <- lapply(files_exp, fread, data.table = FALSE)
head(expression[[1]])
lapply(expression, head)
names <- basename(files_exp)

names(expression) <- gsub("_.*", "", basename(files_exp))
expression_named <- lapply(names(expression), function(x) {
  expression[[x]][,3] <- NULL
  colnames(expression[[x]]) <- c("GeneID", x)
  expression[[x]]})

expression <- Reduce(function(x, y) merge(x, y, by= "GeneID", all=TRUE), expression_named)
View(expression)

raw_exprs_df <- as.data.frame(expression)

# Save the raw expression matrix
write.table(raw_exprs_df, paste0(outDir, gseID, "_","raw_expression_no_norm.tsv"), sep = "\t", row.names = TRUE, col.names = NA)

exprs(eset) <- as.matrix(expression)

###### CLEANING DATA #####
## Cleaning feature
out_dir = file.path(basedir,'intermediate','clean')
dir.create(out_dir, showWarnings = F,recursive = T)

#--- feature data
feature_data <- fData(eset)
View(feature_data)

cols_2select <- c('ID','ILMN_Gene')
feature_data <- feature_data[,cols_2select]
colnames(feature_data) <- c('probe_id','gene_symbol')
feature_data$gene_symbol <- gsub('\\/.*','',feature_data$gene_symbol)
feature_data$gene_symbol <- gsub(' ','',feature_data$gene_symbol)
head(feature_data)

fData(eset) <- feature_data

## Cleaning pheno
pheno_data <- pData(eset)
View(pheno_data)

essential_phenovars <- c('geo_accession', # sample_id
                         'characteristics_ch1.3') # TB info

# secundary_phenovars <- c("characteristics_ch1", "characteristics_ch1.1",
#                          "characteristics_ch1.5", "characteristics_ch1.7",
#                          "characteristics_ch1.8","characteristics_ch1.9", 
#                          "characteristics_ch1.13","characteristics_ch1.15",
#                          "characteristics_ch1.16")

pheno <- pheno_data[,c(essential_phenovars,secundary_phenovars)];head(pheno)
colnames(pheno) <- c('sample_id', 'outcome')#, 'age', 'gender', 
                     # 'BCG_status', "tst_value",
                     # 'exposure','disease_loc', 
                     # 'sputum_culture','bal_culture', 
                     # 'sensitivity');head(pheno)

# clean
pheno <- data.frame(apply(pheno, 2, function(x) gsub('.*\\: ','',x)))
table(pheno$outcome)
pheno$outcome <- gsub(" ", "_", pheno$outcome)
View(pheno)

# continuous_vars <- c('age','tst_value')
# pheno[,continuous_vars] <- data.frame(apply(pheno[,continuous_vars],2, function(x){
#   as.numeric(gsub('--',NA,x))
# }))
# 
# pheno$outcome <- gsub("Control_(BCG+)", "Control", pheno$outcome)
# pheno$outcome <- gsub("Control_(BCG-)", "Control", pheno$outcome)
# View(pheno)

pData(eset) <- pheno # put the new pheno data back into the ExpressionSet


##### CHECK DATA #####
#exp <- as.data.frame(exprs(eset))
exp_raw <- as.data.frame(expression)
view(exp_raw)
ilmn <- exp_raw$GeneID
exp_raw[,2:43] <- apply(exp_raw[,2:43], 2, as.numeric)
View(exp_raw)

exprs_data <- melt(data.frame(exp_raw))
View(exprs_data)
ggplot(exprs_data, aes(x = value, group = variable)) + geom_density(show.legend = F) + theme_minimal(base_size = 15) +
  labs(x = "Expression Value",y = 'Density')


## First quality check
arrayQualityMetrics(, do_logtransform = F, report_title = 'AQM_raw', out_dir = "intermediate/AQM_before_normalization/")



##### NORMALIZE #####
View(exp_raw)
rownames(exp_raw) <- exp_raw$GeneID
exp_raw$GeneID <- NULL
exp_raw <- as.matrix(exp_raw)
min(exp_raw) # -17.59076
exp_raw <- exp_raw+18
min(exp_raw) # 0.40924

exprs_norm <- log2(exp_raw) # transform in log2 scale
exprs_norm <- normalize.quantiles(exprs_norm,copy = F) # normalize
View(exprs_norm)

exprs_norm.l <- melt(data.frame(exprs_norm))
ggplot(exprs_norm.l, aes(x = value, group = variable)) + geom_density(show.legend = F) + theme_minimal(base_size = 15) +
  labs(x = "Expression Value",y = 'Density')


##### ADD GENES #####
exprs_norm <- as.data.frame(exprs_norm)
exprs_norm$probes <- rownames(exprs_norm)
View(exprs_norm)
View(feature_data)

#feature_data$probe_id <- as.character(feature_data$probe_id)
library(dplyr)
exp <- inner_join(x = exprs_norm, y = feature_data, by= c("probes"="probe_id"))
View(exp)

exp <- exp[exp$gene_symbol != "", ]
View(exp)

table(duplicated(exp$gene_symbol))
# FALSE  TRUE 
# 37796 10995 

## Remove duplicates
exp$SUM <- rowSums(exp[,1:42])
exp <- exp[order(exp$SUM, decreasing = T),]
exp_sd <- filter(exp, !duplicated(exp$gene_symbol))
View(exp_sd)
exp <- exp_sd

rownames(exp) <- exp$gene_symbol

exp$SUM <- NULL
exp$probes <- NULL
exp$gene_symbol <- NULL

identical(rownames(pheno),colnames(exp)) # TRUE 


###### CHECK PHENODATA ####
# change sample names in raw_data table to GSM identifiers
pdata <- pData(eset)
View(pdata)
#rownames(pdata) <- gsub('sample name: ','',pdata$characteristics_ch1)
common.samples = intersect(rownames(pdata),colnames(expression))
length(common.samples) == ncol(expression)
length(common.samples) == nrow(pdata)

raw_exprs <- expression[,common.samples]
identical(rownames(pdata),colnames(raw_exprs)) # TRUE

#colnames(raw_exprs) <- pdata$geo_accession
View(head(raw_exprs,20))



##### Check genes in both studies ####

tb <- read.delim("GSE19439/data/GSE19439_expr_norm_collapsed_final.txt")
#View(tb)
tb$remove_genes <- NULL

tb$vars <- rowVars(as.matrix(tb))
tb <- tb[order(tb$vars, decreasing = TRUE),]

0.75*length(tb$vars) # 28347
tb <- tb[1:28347, ]
#View(tb)

####
tb2 <- read.delim("GSE19442/data/GSE19442_expr_norm_collapsed_final.txt")
#View(tb2)

tb2$vars <- rowVars(as.matrix(tb2))
tb2 <- tb2[order(tb2$vars, decreasing = TRUE),]

0.75*length(tb2$vars) # 28353
tb2 <- tb2[1:28353, ]
#View(tb2)

####
inter <- intersect(rownames(tb), rownames(tb2))

tb <- tb[inter,]
tb2 <- tb2[inter,]

identical(rownames(tb), rownames(tb2)) # TRUE

tb$Probes <- rownames(tb)
dim(tb)
tb <- tb[,c(44,1:43)]
tb$vars <- NULL
#View(tb)

tb2$Probes <- rownames(tb2)
dim(tb2)
tb2 <- tb2[,c(53,1:52)]
#View(tb2)
tb2$vars <- NULL

####
pheno_tb1 <- read.delim("GSE19439/data/GSE19439_phenodata.txt")
#View(pheno_tb1)
pheno_tb1 <- pheno_tb1[pheno_tb1$outcome %in% c("PTB", "Latent"),]

both <- intersect(pheno_tb1$sample_id, colnames(tb))
tb <- tb[,both]
pheno_tb1 <- pheno_tb1[pheno_tb1$sample_id %in% both,]

identical(pheno_tb1$sample_id, colnames(tb)) # TRUE

pheno_tb1 <- pheno_tb1[,1:2]
colnames(pheno_tb1) <- c("Probes", "Class")
View(pheno_tb1)

tb$Probes <- rownames(tb)
dim(tb)
tb <- tb[,c(31,1:30)]

####
write.table(tb, "Intermed/Tuberculose_BioFeatS/nrom_tb1_biofeats.txt", row.names = FALSE, quote = FALSE)
write.table(pheno_tb1, "Intermed/Tuberculose_BioFeatS/pheno_tb1_biofeats.txt", row.names = FALSE, quote = FALSE)

####
pheno_tb2 <- read.delim("GSE19442/data/GSE19442_phenodata.txt")
View(pheno_tb2)

colnames(pheno_tb2) <- c("Probes", "Class")
pheno_tb2$Class <- gsub("LATENT_TB", "Latent", pheno_tb2$Class)

View(pheno_tb2)
View(tb2)

####
write.table(tb2, "Intermed/Tuberculose_BioFeatS/norm_tb2_biofeats.txt", row.names = FALSE, quote = FALSE)
write.table(pheno_tb2, "Intermed/Tuberculose_BioFeatS/pheno_tb2_biofeats.txt", row.names = FALSE, quote = FALSE)


