#### ABs ####
####### Luminex D3 - CXCL10
lumi_d3 <- read.delim("Intermed/Escolha_features/D3/AB/1_Body_Xv0_Log2_Luminex_D3vD0.txt")
View(lumi_d3)

rownames(lumi_d3) <- lumi_d3$Probes
lumi_d3$Probes <- NULL
lumi_d3 <- t(lumi_d3)
n_lumi <- rownames(lumi_d3)
View(lumi_d3)
lumi_d3 <- as.data.frame(lumi_d3)
lumi_d3$CXCL10 <- as.numeric(lumi_d3$CXCL10)
View(lumi_d3)

ggplot(lumi_d3, aes(x=Class, y=CXCL10)) + 
  geom_boxplot()


####### Luminex D7 - OSM
lumi_d7 <- read.delim("Intermed/Escolha_features/D7/AB/1_Body_Xv0_Log2_Luminex_D7vD0.txt")
View(lumi_d7)

rownames(lumi_d7) <- lumi_d7$Probes
lumi_d7$Probes <- NULL
lumi_d7 <- t(lumi_d7)
n_lumi <- rownames(lumi_d7)
View(lumi_d7)
lumi_d7 <- as.data.frame(lumi_d7)
lumi_d7$OSM <- as.numeric(lumi_d7$OSM)
View(lumi_d7)

ggplot(lumi_d7, aes(x=Class, y=OSM)) + 
  geom_boxplot()

#### AEs ####
## Luminex D3
lumi_d3 <- read.delim("Intermed/Escolha_features/D3/AE/1_Adverse_Xv0_Log2_Luminex_D3vD0.txt")
View(lumi_d3)

rownames(lumi_d3) <- lumi_d3$Probes
lumi_d3$Probes <- NULL
lumi_d3 <- t(lumi_d3)
n_lumi <- rownames(lumi_d3)
View(lumi_d3)
lumi_d3 <- as.data.frame(lumi_d3)
lumi_d3$TRAIL <- as.numeric(lumi_d3$TRAIL)
View(lumi_d3)

ggplot(lumi_d3, aes(x=Class, y=TRAIL)) + 
  geom_boxplot()

## miRNA D7
mirna <- read.delim("Intermed/Escolha_features/D7/AE/1_Adverse_Xv0_Log2_miRNA_Normalized_Pseudocounts_D7vD0.txt")
View(mirna)

rownames(mirna) <- mirna$Probes
mirna$Probes <- NULL
mirna <- t(mirna)
n_mirna <- rownames(mirna)
mirna <- as.data.frame(mirna)
mirna$`hsa-miR-640` <- as.numeric(mirna$`hsa-miR-640`)
View(mirna)

ggplot(mirna, aes(x=Class, y=mirna$`hsa-miR-640`)) + 
  geom_boxplot()



