#####
rm(list = ls())
options(stringsAsFactors = F)

#library()
pkgs <- c('tidyverse')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)

#dir.create("Intermed/data_BioFeats_USA/final_data")

#####
tab <- read.csv("Intermed/data_BioFeats_USA/normalized_counts/norm_counts_d1.csv")
View(tab)


#####
list_count <- list.files(path="Intermed/data_BioFeats_USA/normalized_counts/", full.names=TRUE)
counts <- lapply(list_count, read.csv) %>% setNames(str_replace(basename(list_count),".csv", "")) 

for (i in names(counts)){
  print(i)
  tabela = counts[[i]]
  colnames(tabela)[1] <- "Probes" 
  fname = paste0('Intermed/data_BioFeats_USA/final_data/',i, ".txt")
  print(fname)
  write.table(tabela, fname, quote = FALSE, sep = "\t",  row.names = FALSE)
}

tab <- read.delim("Intermed/data_BioFeats_USA/final_data/norm_counts_d1.txt")
View(tab)


##### 
list_pheno <- list.files(path="Intermed/data_BioFeats_USA/phenodata/", full.names=TRUE)
phenos <- lapply(list_pheno, read.csv) %>% setNames(str_replace(basename(list_pheno),".csv", ""))
#View(phenos[[1]])

for (i in names(phenos)){
  tabela = phenos[[i]]
  tabela$X <- NULL
  colnames(tabela) <- c("Probes", "Class")
  tabela$Class <- gsub(0, "D0", tabela$Class)
  tabela$Class <- gsub(1, "D1", tabela$Class)
  phenoname = paste0('Intermed/data_BioFeats_USA/final_data/',i, ".txt")
  write.table(tabela, phenoname, quote = FALSE, sep = "\t",  row.names = FALSE)
}

tab <- read.delim("Intermed/data_BioFeats_USA/final_data/pheno_USA_d1.txt")
View(tab)





