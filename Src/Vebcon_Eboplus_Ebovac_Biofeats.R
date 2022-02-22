######
rm(list = ls())
options(stringsAsFactors = F)

sapply(c('ggplot2','ggrepel','dplyr','Glimma', "tidyr", "edgeR", "tmod"), require,character.only=T)

##### EBOVAC ####
descript_ebovac <- read.csv2("data/Pheno_Geneva.csv")
count_ebovac <- read.csv2("data/Counts_Geneva.csv")
rownames(count_ebovac) <- count_ebovac$X
#View(head(count_ebovac))
#summary(descript_ebovac)

rownames(descript_ebovac) <- paste0("X", descript_ebovac$Sample)
descript_ebovac$Sample <- rownames(descript_ebovac)
count_ebovac$X <- NULL
descript_ebovac <- descript_ebovac[descript_ebovac$Treatment!="P",]
table(descript_ebovac$Treatment == "P")
#dim(count_ebovac)
descript_ebovac <- descript_ebovac[descript_ebovac$Group.Day %in% c(0,7),]
descript_ebovac <- descript_ebovac[descript_ebovac$Library.Batch == 2,]
dim(descript_ebovac)

common <- intersect(colnames(count_ebovac), rownames(descript_ebovac))
count_ebovac <- count_ebovac[,common]
View(descript_ebovac)
setdiff(colnames(count_ebovac), rownames(descript_ebovac))
identical(colnames(count_ebovac), rownames(descript_ebovac)) # TRUE

dim(count_ebovac)
View(descript_ebovac)
table(descript_ebovac$Volunteer.ID)

# create DEGlist
y_ebovac <- DGEList(counts = count_ebovac, genes = row.names(count_ebovac), group= descript_ebovac$Group.Day)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac)

# normalize
y.1 <- calcNormFactors(y_ebovac)

p <- cpm(y.1, log = TRUE, prior.count = 1)

View(p)
p <- as.data.frame(p)
p$Probes <- rownames(p)
p <- p[,c(50, 1:49)]


##### EBOPLUS ####
ebop <- read.delim("Intermed/data_BioFeats_USA/final_data/norm_counts_d7.txt")
View(ebop)


##### VEBCON ####
veb <- read.delim("Intermed/vebcon_norm_counts_d7.txt")
View(veb)


inter <- Reduce(intersect, list(p$Probes, ebop$Probes, veb$Probes))

p <- p[p$Probes %in% inter,]
ebop <- ebop[ebop$Probes %in% inter,]
veb <- veb[veb$Probes %in% inter,]

#####

write.table(p, "Intermed/BioFeatS_EBOP_EBOV_VEB/D7/norm_ebovac_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(ebop, "Intermed/BioFeatS_EBOP_EBOV_VEB/D7/norm_eboplus_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')
write.table(veb, "Intermed/BioFeatS_EBOP_EBOV_VEB/D7/norm_vebcon_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')

tab <- read.delim("Intermed/BioFeatS_EBOP_EBOV_VEB/D7/norm_eboplus_d7.txt")
View(tab)

pheno <- read.delim("Intermed/Biofeats_USA_Vebcon_D7/pheno_USA_d7.txt")
View(pheno)
identical(pheno$Probes, colnames(ebop[,-1])) # TRUE
pheno$Class <- gsub("7", "D7", pheno$Class)
write.table(pheno, "Intermed/BioFeatS_EBOP_EBOV_VEB/D7/pheno_eboplus_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')

identical(descript_ebovac$Sample, colnames(p[,-1])) # TRUE

pheno_ebov <- descript_ebovac[,c(1,4)]
colnames(pheno_ebov) <- c("Probes", "Class")
pheno_ebov$Class <- gsub("0", "D0", pheno_ebov$Class)
pheno_ebov$Class <- gsub("7", "D7", pheno_ebov$Class)
View(pheno_ebov)

write.table(pheno_ebov, "Intermed/BioFeatS_EBOP_EBOV_VEB/D7/pheno_ebovac_d7.txt", quote = FALSE, row.names = FALSE, sep = '\t')


