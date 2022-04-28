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

############
input.file="rnaseq/results/rnaflow_test/07-DifferentialExpression/DESeq2/MAQCA_vs_MAQCB/results/deseq2_MAQCA_MAQCB_full_extended.csv"
output="./"
GOdb<- 'org.Hs.eg.db'
KEGGdb <- 'hsa'


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
png(paste0(output,"network_GO2gene.png"),res = 600,width = 10,height = 10,units = "in")
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE,cex_label_gene=2)
dev.off()

png(paste0(output,"network_KEGG2gene.png"),res = 600,width = 10,height = 10,units = "in")
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE,cex_label_gene=2)
dev.off()


