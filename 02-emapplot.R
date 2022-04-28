list.of.packages <- c("openxlsx", "ggplot2","stringr","enrichplot","clusterProfiler","GOplot","DOSE","ggnewscale","circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
BiocManager::install("ComplexHeatmap")
BiocManager::install("topGO")

library(openxlsx)
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
library(topGO)
library(circlize)
library(ComplexHeatmap)

##Setup

input.file="rnaseq/results/rnaflow_test/07-DifferentialExpression/DESeq2/MAQCA_vs_MAQCB/results/deseq2_MAQCA_MAQCB_full_extended.csv"
output="Yangtui/Fancy/"
GOdb ='org.Hs.eg.db'
KEGGdb='hsa'


input=read.csv(input.file,row.names = 1)
gene <- bitr(input$geneName,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GOdb)


GO<-enrichGO(gene$ENTREZID,
             OrgDb = GOdb,
             keyType = "ENTREZID",
             ont = "ALL",
             pvalueCutoff = 0.1,
             qvalueCutoff = 0.1,
             readable = T)

KEGG<-enrichKEGG(gene$ENTREZID,
                 organism = KEGGdb,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
##Plot
GO2 <- pairwise_termsim(GO)
png(paste0(output,"network2_GO2GO.png"),res = 600,width = 15,height = 15,units = "in")
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")
dev.off()

##Plot
KEGG2=pairwise_termsim(KEGG)
png(paste0(output,"network2_KEGG2KEGG.png"),res = 600,width = 15,height = 15,units = "in")
enrichplot::emapplot(KEGG2,showCategory = 50, color = "p.adjust", layout = "kk")
dev.off()
