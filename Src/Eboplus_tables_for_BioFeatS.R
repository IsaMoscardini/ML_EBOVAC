library(edgeR)

counts <- read.csv2("Data/Counts_USA.csv")
rownames(counts) <- counts$Gene
#View(counts)
pheno <- read.csv2("Data/Pheno_USA.csv")
#View(pheno)
pheno <- pheno[pheno$Sampling.Day %in% c(0, 7), ]
#View(pheno)

pheno <- pheno[,c(1,3)]
rownames(pheno) <- pheno$Sample
#View(pheno)

counts <- counts[,which(colnames(counts) %in% pheno$Sample)]

identical(colnames(counts), rownames(pheno))
#View(counts)

y <- DGEList(counts = counts, genes = row.names(counts), group= pheno$Sampling.Day)

# filtra
keep <- rowSums(cpm(y)>1)>=10
y.1 <- y[keep,]
dim(y.1) #  12539    63

# normaliza
y.1 <- calcNormFactors(y.1)
cp <- log2(cpm(y.1)+1)
#View(cp)
identical(colnames(cp), rownames(pheno))

write.csv(cp, "Intermed/norm_counts_d7.csv")
write.csv(pheno, "Intermed/pheno_USA_d7.csv")

tab <- read.csv("Intermed/norm_counts_d7.csv")
View(tab)
tab2 <- read.csv("Intermed/pheno_USA_d7.csv")

rownames(tab) <- tab$X
tab$X <- NULL

