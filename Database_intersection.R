library(CEMiTool)

REACTOME <- read_gmt("GMT_files/ReactomePathways.gmt")
BTM <- read_gmt("GMT_files/BTM_for_GSEA_20131008.gmt")
KEGG <- read_gmt("GMT_files/c2.cp.kegg.v7.2.symbols (3).gmt")
GO <- read_gmt("GMT_files/c5.go.v7.2.symbols (1).gmt")
HALL <- read_gmt("GMT_files/c7.all.v7.2.symbols.gmt")


##### REAC KEGG E BTM
genes <- unique(c(REACTOME$gene, KEGG$gene, BTM$gene))
View(genes)

count <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")

inters <- intersect(count$Probes, genes) 
View(inters)


##### REAC KEGG BTM HALLMARKS AND GO
genes <- unique(c(REACTOME$gene, KEGG$gene, BTM$gene, GO$gene, HALL$gene))
View(genes)

count <- read.delim("Data/D1/RNAseq_logCPM_D1vD0.txt")

intersect <- intersect(count$Probes, genes) 
View(intersect)

write.csv(intersect, "Intermed/Genes_bio.csv", row.names = F)
gen <- read.csv("Intermed/Genes_bio.csv")
View(gen)
