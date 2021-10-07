######
rm(list = ls())
options(stringsAsFactors = F)

#library()
pkgs <- c('dplyr')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)

######
rna_D1 <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
View(rna_D1)

pheno <- read.csv2("Data/pheno_geneva.csv")
pheno$Volunteer.ID <- paste0("G", pheno$Volunteer.ID)

pheno <- pheno[pheno$Volunteer.ID %in% colnames(rna_D1),]
pheno <- subset(pheno, !duplicated(pheno$Volunteer.ID))
View(pheno)

gender <- pheno[, c(19,11)]
colnames(gender) <- c("Samples", "Class")

dose <- pheno[, c(19,2)]
colnames(dose) <- c("Samples", "Class")
class(dose$Class)
dose$Class <- as.factor(dose$Class)
dose$Class <- gsub("Placebo", 0, dose$Class)
dose$Class <- gsub("5x10^7", "H", dose$Class)
dose$Class <- gsub("10^7", "L", dose$Class)
View(dose)

############

