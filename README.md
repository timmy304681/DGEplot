# DGEplot

It's a typical rnaseq data from [deseq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). The samples are human samples.

[deseq2_A_C_filtered_padj_0.05_extended.csv](https://github.com/timmy304681/DGEplot/files/8076564/deseq2_A_C_filtered_padj_0.05_extended.csv)

Library the whole packages you need.
```
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
```


Do `GO`,`KEGG`,`GSEA` analysis.
```
#Setup & Using human database
output="./OUTPUT"
GOdb<- 'org.Hs.eg.db'
KEGGdb <- 'hsa'


input=read.csv("deseq2_A_C_filtered_padj_0.05_extended.csv",row.names = 1)
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

names(info) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
input_merge <- merge(input,gene,by.x='geneName',by.y='SYMBOL') 
GSEA_input <- input_merge$log2FoldChange
names(GSEA_input) = input_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, 
                     organism = KEGGdb, 
                     pvalueCutoff = 0.05)

GSEA_GO <- gseGO(GSEA_input, 
                 ont ="ALL",
                 #keyType = "SYMBOL",  #keyType = "ENTREZID"
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.5,
                 verbose = TRUE,
                 OrgDb = org.Hs.eg.db,
                 pAdjustMethod = "BH")

```
## GO/KEGG barplot & dotplot
```

png(paste0(output,"barplot.png"),res = 600,width = 10,height = 10,units = "in")
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dev.off()

png(paste0(output,"dotplot.png"),res = 600,width = 10,height = 10,units = "in")
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dev.off()
```
展示GO 三種屬性 (Cellular Component / Biological Process / Molecular Function) 分類和基因表達量之關係，以柱狀圖表示。三個方格分別代表 CC / BP / MF三種屬性，直行為GO Id，橫列為adj.pvalue
顏色漸層則表示p.adjusted (cut off default: 0.05)

<img src="https://user-images.githubusercontent.com/51151276/154193572-69da0f39-b88c-4493-a34a-daf8a571924d.png" alt="barplot" width="50%"/>

展示GO 三種屬性 (Cellular Component / Biological Process / Molecular Function) 分類和基因表達量之關係，以點狀圖表示。三個方格分別代表 CC / BP / MF三種屬性，直行為GO id，橫列為 Gene Ratio，顏色漸層表示p.adjusted (cut off default: 0.05) ，圓圈大小表示Count。

<img src="https://user-images.githubusercontent.com/51151276/154193584-8839a4a1-c4e4-48c7-8d23-c9bf517f5ae1.png" alt="dotplot" width="50%"/>



## GO/KEGG 
```
png(paste0(output,"network_GO2gene.png"),res = 600,width = 10,height = 10,units = "in")
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE,cex_label_gene=2)
dev.off()
```


