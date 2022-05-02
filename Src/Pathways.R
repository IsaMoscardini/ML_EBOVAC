######
rm(list = ls())
options(stringsAsFactors = F)

#library()
pkgs <- c('CEMiTool', 'dplyr')
sum(unlist(lapply(pkgs, require,character.only = T))) == length(pkgs)

######

kegg <- read_gmt("c2.cp.kegg.v7.4.symbols.gmt")
View(kegg)

# D1
ab_d1 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D1vD0.csv")
View(ab_d1)
ab_d1 <- ab_d1$Probes

kegg_ab_d1 <- as.data.frame(kegg)
kegg_ab_d1 <- kegg_ab_d1[kegg_ab_d1$gene %in% ab_d1,]
View(kegg_ab_d1)

# D3
ab_d3 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D3vD0.csv")
ab_d3 <- ab_d3$Probes

kegg_ab_d3 <- as.data.frame(kegg)
kegg_ab_d3 <- kegg_ab_d3[kegg_ab_d3$gene %in% ab_d3,]
View(kegg_ab_d3)



####### BTM #####

#### AB
# D1
btm <- read_gmt("BTM_for_GSEA_20131008.gmt")
View(btm)

ab_d1 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D1vD0.csv")
View(ab_d1)
ab_d1 <- ab_d1$Probes

btm_ab_d1 <- as.data.frame(btm)
btm_ab_d1 <- btm_ab_d1[btm_ab_d1$gene %in% ab_d1,]
View(btm_ab_d1)

#D3
ab_d3 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D3vD0.csv")
View(ab_d3)
ab_d3 <- ab_d3$Probes

btm_ab_d3 <- as.data.frame(btm)
btm_ab_d3 <- btm_ab_d3[btm_ab_d3$gene %in% ab_d3,]
View(btm_ab_d3)

#D7
ab_d7 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D7vD0.csv")
View(ab_d7)
ab_d7 <- ab_d7$Probes

btm_ab_d7 <- as.data.frame(btm)
btm_ab_d7 <- btm_ab_d7[btm_ab_d7$gene %in% ab_d7,]
View(btm_ab_d7)

#D14
ab_d14 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D14vD0.csv")
View(ab_d14)
ab_d14 <- ab_d14$Probes

btm_ab_d14 <- as.data.frame(btm)
btm_ab_d14 <- btm_ab_d14[btm_ab_d14$gene %in% ab_d14,]
View(btm_ab_d14)

#D28
ab_d28 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D28vD0.csv")
View(ab_d28)
ab_d28 <- ab_d28$Probes

btm_ab_d28 <- as.data.frame(btm)
btm_ab_d28 <- btm_ab_d28[btm_ab_d28$gene %in% ab_d28,]
View(btm_ab_d28)




# miRNA
mirna_ab <- read.csv("Intermed/miRWalk_miRNA_TargetsAB.csv")
View(mirna_ab)
rnad1 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D1vD0.csv")
rnad3 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D3vD0.csv")
rnad7 <- read.csv("Results/Selected_features_higher_0/Ab/RNAseq_logCPM_D7vD0.csv")

mirna_ab <- as.data.frame(mirna_ab)

mirna_inter_d1 <- mirna_ab[which(mirna_ab$genesymbol %in% rnad1$Probes),]
View(mirna_inter_d1)
mirna_inter_d3 <- mirna_ab[which(mirna_ab$genesymbol %in% rnad3$Probes),]
View(mirna_inter_d3)
mirna_inter_d7 <- mirna_ab[which(mirna_ab$genesymbol %in% rnad7$Probes),]
View(mirna_inter_d7)






