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



## GO/KEGG emapplot
```
GO2 <- pairwise_termsim(GO)
png(paste0(output,"network2_GO2GO.png"),res = 600,width = 15,height = 15,units = "in")
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")
dev.off()
```

展示GO term與GO term之間的關係(或是KEGG與KEGG之間的關係)，每個圓點代表一個GO term/ KEGG，圓點大小代表裡面基因數量，圓點數量(max default: 50)，圓點顏色代表p.adjusted (cut off default: 0.05)，線條相連代表兩者GO term/ KEGG之間有關聯

<img src="https://user-images.githubusercontent.com/51151276/154195144-98143a06-ca2b-4b68-8c02-b3e67d0fdbf2.png" alt="network2_GO2GO" width="50%"/>

## GO/KEGG cnetplot 
```
png(paste0(output,"network_GO2gene.png"),res = 600,width = 10,height = 10,units = "in")
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE,cex_label_gene=2)
dev.off()
```
展示Gene與GO/KEGG之間複雜的網絡關係，大黃圓代表GO/KEGG，小灰圓代表Gene，基因可設定有多少連結點才顯示基因名稱

<img src="https://user-images.githubusercontent.com/51151276/154195541-754846d4-6b75-4eae-8ac8-43527d741e31.png" alt="network2_GO2GO2" width="50%"/>


## GO chord plot
```
genedata<-data.frame(ID=input$geneName,logFC=input$log2FoldChange)

BP_end=which(!duplicated(GO$ONTOLOGY))[2]-1
CC_start=which(!duplicated(GO$ONTOLOGY))[2]
CC_end=which(!duplicated(GO$ONTOLOGY))[3]-1
MF_start=which(!duplicated(GO$ONTOLOGY))[3]
MF_end=length(GO$ONTOLOGY)

#extract ID,Description,p.adjust,GeneID
GOplotIn_BP<-GO[1:BP_end,c(2,3,7,9)] 
GOplotIn_CC<-GO[CC_start:CC_end,c(2,3,7,9)]
GOplotIn_MF<-GO[MF_start:MF_end,c(2,3,7,9)]

#把GeneID列中的’/’替換成‘,’
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') 
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')

#修改列名,後面弦圖繪制的時候需要這樣的格式
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')

#分類信息
GOplotIn_BP$Category = "BP"
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"

#GOplot導入數據格式整理
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata)
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata) 

#生成含有選定基因的數據框
chord_BP<-chord_dat(data = circ_BP,genes = genedata)
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MF,genes = genedata) 

png(paste0(output,"Chord_BP.png"),res = 600,width = 15,height = 15,units = "in")
GOChord(data = chord_BP[1:30,],
        title = 'GO-Biological Process',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'),
        process.label = 10) 
dev.off()

png(paste0(output,"Chord_CC.png"),res = 600,width = 15,height = 15,units = "in")
GOChord(data = chord_CC[1:30,],title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10) 
dev.off()
png(paste0(output,"Chord_MF.png"),res = 600,width = 15,height = 15,units = "in")
GOChord(data = chord_MF[1:30,],title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10)
dev.off()
```
展示基因和GO term之間的關係。左半圈表示基因名稱,右半圈每個色帶代表一個GO term，並標示於下方,logFC表示基因差異倍數取log



![Chord_BP](https://user-images.githubusercontent.com/51151276/154196818-b8dc319b-1500-4b0a-8a0c-85552f79772f.png)


