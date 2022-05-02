options(stringsAsFactors = FALSE)


##### RODAR UMA VEZ save counts and phenodata ####
library(readr)
vebcon_counts<- read_delim("Data/vebcon_edited.tsv.txt",
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE)

vebcon_pheno<-names(vebcon_counts)

vebcon_pheno_df <-as.data.frame(vebcon_pheno[-1])

library(stringr)
vebcon_pheno_df_2 <-as.data.frame(str_split_fixed(vebcon_pheno_df$`vebcon_pheno[-1]`, "_d", 2))

vebcon_pheno_f <-as.data.frame(vebcon_pheno_df_2$V2)

rownames(vebcon_pheno_f) <- vebcon_pheno_df$`vebcon_pheno[-1]`

#View(vebcon_counts)

## Convert ensembl in gene symbol
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ens <- vebcon_counts$Name
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=ens,mart= mart)

vebcon_counts$gene <- G_list$hgnc_symbol[match(vebcon_counts$Name, G_list$ensembl_gene_id)]
#View(vebcon_counts)
vebcon_counts <- vebcon_counts[which(vebcon_counts$gene != ""),]
vebcon_counts <- vebcon_counts[complete.cases(vebcon_counts$gene),]
#View(vebcon_counts)

#write.csv(vebcon_counts, "Data/Vebcon_counts_genesymbol.csv", row.names = FALSE)
#write.csv(vebcon_pheno_f, "Data/Vebcon_pheno.csv", row.names = FALSE)

####
pheno <- vebcon_pheno_f
pheno$Sample <- rownames(pheno)
pheno <- pheno[,c(2,1)]
colnames(pheno) <- c("Probes", "Class")
#View(pheno)
pheno <- pheno[pheno$Class %in% c(0, 3), ]

pheno$Class <- gsub("0", "D0", pheno$Class)
pheno$Class <- gsub("3", "D3", pheno$Class)
View(pheno)
write.table(pheno, "Intermed/Scaled_Ebolas_Apr2022/pheno_vebcon_d3.txt", row.names = FALSE, quote = FALSE, sep = '\t')
counts <- vebcon_counts
dim(counts)
counts$Name <- NULL
View(counts)
counts$SUM <- rowSums(counts[,1:56])
counts <- counts[order(counts$SUM, decreasing = T),]
counts <- subset(counts, !duplicated(counts$gene))
counts <- as.data.frame(counts)
rownames(counts) <- counts$gene
counts <- counts[,which(colnames(counts) %in% pheno$Probes)]

identical(colnames(counts), rownames(pheno))
View(counts)

library(edgeR)
y <- DGEList(counts = counts, genes = row.names(counts), group= pheno$Class)

# filtra
keep <- rowSums(cpm(y)>1)>=10
y.1 <- y[keep,]
dim(y.1) #  

# normaliza
y.1 <- calcNormFactors(y.1)
cp <- log2(cpm(y.1)+1)
#View(cp)
identical(colnames(cp), rownames(pheno))

cp <- as.data.frame(cp)
cp$Probes <- rownames(cp)
dim(cp)
cp <- cp[c(29, 1:28)]
View(cp)

#### USA
usa <- read.delim("Intermed/data_BioFeats_USA/final_data/norm_counts_d1.txt")
View(usa)

#### Geneva 
geneva <- read.delim("Intermed/data_BioFeatS_Geneva/final_data/Geneva_norm_counts_d1.txt")
View(geneva)

inter <- intersect(usa$Probes, cp$Probes)
inter <- intersect(inter, geneva$Probes)
usa <- usa[usa$Probes %in% inter, ]
cp <- cp[cp$Probes %in% inter,]
geneva <- geneva[geneva$Probes %in% inter, ]
dim(usa) # all 10126 genes

identical(cp$Probes, usa$Probes)
View(usa)

#### Save
write.table(usa, "Intermed/BioFeatS_EBOP_EBOV_VEB/D1/usa_norm_counts_d1.txt", quote = FALSE, sep = '\t', row.names = FALSE)
write.table(cp, "Intermed/BioFeatS_EBOP_EBOV_VEB/vebcon_norm_counts_d3.txt", quote = FALSE, sep = '\t', row.names = FALSE)
write.table(geneva, "Intermed/BioFeatS_EBOP_EBOV_VEB/D1/geneva_norm_counts_d1.txt", quote = FALSE, sep = '\t', row.names = FALSE)
