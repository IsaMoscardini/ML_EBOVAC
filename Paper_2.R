######
rm(list = ls())
options(stringsAsFactors = F)

#library()
pkgs <- c('dplyr')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)

######
rna_D1 <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")
#View(rna_D1)

pheno <- read.csv2("Data/pheno_geneva.csv")
pheno$Volunteer.ID <- paste0("G", pheno$Volunteer.ID)

pheno <- pheno[pheno$Volunteer.ID %in% colnames(rna_D1),]
pheno <- subset(pheno, !duplicated(pheno$Volunteer.ID))

#gender <- pheno[, c(19,11)]
#colnames(gender) <- c("Samples", "Class")

View(pheno)
#pheno <- pheno[pheno$Sampling.Day == 1,]
dose <- pheno[, c(19,2)]
colnames(dose) <- c("Samples", "dose")
#class(dose$Class)
#dose$Class <- as.character(dose$Class)
# dose$Class <- gsub("Placebo", "P", dose$Class)
# dose$Class <- gsub("5x10^7", "H", dose$Class)
# dose$Class <- gsub("10^7", "L", dose$Class)
#View(dose)
dose$Class <- c("L", "L", "L", "L", "L", "H", "L", "L", "L", "P", "P", "L", "H", "L", "P", 
                "H", "L", "L", "L", "H", "H", "L", "L", "L", "L", "L", "L", "L", "L", "H",
                "L", "L", "L", "H", "H", "H", "L", "L", "H", "P", "P", "P", "P", "L", "H",
                "L", "L", "L", "L", "L", "P", "P", "H", "H", "H", "P", "P", "P", "L", "H",
                "L", "P", "L", "H")

dose$dose <- NULL
dose$Vaccine <- dose$Class
dose$Vaccine <-  gsub("H", "V", dose$Vaccine)
dose$Vaccine <-  gsub("L", "V", dose$Vaccine)

vaccine <- dose
vaccine$Class <- NULL

View(vaccine)
write.csv(vaccine, "Intermed/Outcome_VACvsPLAC.csv", row.names = F)
vc <- read.csv("Intermed/Outcome_VACvsPLAC.csv")
View(vc)


############

counts_geneva <- read.csv2("Data/Counts_USA.csv")
View(counts_geneva)

desc_geneva <- read.csv2("Data/Pheno_USA.csv")
View(desc_geneva)

