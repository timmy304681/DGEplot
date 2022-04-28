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

#Plot
genedata<-data.frame(ID=input$geneName,logFC=input$log2FoldChange)
BP_end=which(!duplicated(GO$ONTOLOGY))[2]-1
CC_start=which(!duplicated(GO$ONTOLOGY))[2]
CC_end=which(!duplicated(GO$ONTOLOGY))[3]-1
MF_start=which(!duplicated(GO$ONTOLOGY))[3]
MF_end=length(GO$ONTOLOGY)

GOplotIn_BP<-GO[1:BP_end,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC<-GO[CC_start:CC_end,c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF<-GO[MF_start:MF_end,c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替換成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,後面弦圖繪制的時候需要這樣的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
GOplotIn_BP$Category = "BP"#分類信息
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) #GOplot導入數據格式整理
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata) 
chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有選定基因的數據框
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MF,genes = genedata) 

png(paste0(output,"Chord_BP.png"),res = 600,width = 15,height = 15,units = "in")
GOChord(data = chord_BP[1:30,],#弦圖
        title = 'GO-Biological Process',space = 0.01,#GO Term間距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下調基因顏色
        process.label = 10) #GO Term字體大小
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

