list.of.packages <- c( "ggplot2","stringr","enrichplot","clusterProfiler","GOplot","DOSE")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)



##Setup

input.file="rnaseq/results/rnaflow_test/07-DifferentialExpression/DESeq2/MAQCA_vs_MAQCB/results/deseq2_MAQCA_MAQCB_full_extended.csv"
output="Yangtui/Fancy/"
GOdb ='org.Hs.eg.db'


input=read.csv(input.file,row.names = 1)
gene <- bitr(input$geneName,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GOdb)


GO<-enrichGO(gene$ENTREZID,
             OrgDb = GOdb,
             keyType = "ENTREZID",
             ont = "ALL",
             pvalueCutoff = 0.1,
             qvalueCutoff = 0.1,
             readable = T)

###Plot
png(paste0(output,"barplot.png"),res = 600,width = 10,height = 10,units = "in")
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dev.off()

png(paste0(output,"dotplot.png"),res = 600,width = 10,height = 10,units = "in")
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dev.off()



