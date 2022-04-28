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
output="Yangtui/Fancy/"
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

names(info) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
input_merge <- merge(input,gene,by.x='geneName',by.y='SYMBOL')#合並轉換後的基因ID和Log2FoldChange
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


###Ouputtable
write.table(GO,paste0(output,"result_GO.txt"),sep='\t',quote=F)
write.table(KEGG,paste0(output,"result_KEGG.txt"),sep='\t',quote=F)
write.table(GSEA_GO,paste0(output,"result_GSEA.txt"),sep='\t',quote=F)
write.table(gene,paste0(output,"result_genes.txt"),sep='\t',row.names=F,quote=F)
gene.table=input[,c(1,8)]
colnames(gene.table)=c("SYMBOL","log2FC")
write.table(gene.table,paste0(output,"result_diffexp.txt"),sep='\t',row.names=F,quote=F)

#GO/KEGG富集柱狀圖+點狀圖
png(paste0(output,"barplot.png"),res = 600,width = 10,height = 10,units = "in")
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dev.off()

png(paste0(output,"dotplot.png"),res = 600,width = 10,height = 10,units = "in")
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dev.off()

#富集基因與所在功能集/通路集的關聯網絡圖
png(paste0(output,"network_GO2gene.png"),res = 600,width = 10,height = 10,units = "in")
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE,cex_label_gene=2)
dev.off()


#熱圖形式展現關聯關系:
GO.test<-enrichGO(gene$ENTREZID[1:300],
             OrgDb = GOdb,
             keyType = "ENTREZID",
             ont = "ALL",
             pvalueCutoff = 0.1,
             qvalueCutoff = 0.1,
             readable = T)

png(paste0(output,"GO_heatmap.png"),res = 600,width = 10,height = 6,units = "in")
enrichplot::heatplot(GO.test,showCategory = 10)#基因-通路關聯熱圖
dev.off()

enrichplot::heatplot(KEGG,showCategory = 50)


#富集到的功能集/通路集之間的關聯網絡圖
GO2 <- pairwise_termsim(GO)
png(paste0(output,"network2_GO2GO.png"),res = 600,width = 15,height = 15,units = "in")
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")
dev.off()

#GO富集弦圖
genedata<-data.frame(ID=input$geneName,logFC=input$log2FoldChange)

write.table(GO$ONTOLOGY, file = paste0(output,"GO_ontology.txt"), 
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
GO$ONTOLOGY
which(!duplicated(GO$ONTOLOGY))
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

####GO富集弦表圖

png(paste0(output,"Circle_BP.png"),res = 600,width = 15,height = 8,units = "in")
GOCircle(circ_BP)
dev.off()


png(paste0(output,"Circle_CC.png"),res = 600,width = 15,height = 8,units = "in")
GOCircle(circ_CC)
dev.off()

png(paste0(output,"Circle_MF.png"),res = 600,width = 15,height = 8,units = "in")
GOCircle(circ_MF)
dev.off()

###GO富集系統聚類圖

chord<-chord_dat(data = circ_BP,genes = genedata) #生成含有選定基因的數據框

png(paste0(output,"CirHeatmap_BP.png"),res = 600,width = 15,height = 15,units = "in")
GOCluster(circ_BP,GOplotIn_BP$Term) #系統聚類圖
dev.off()

png(paste0(output,"CirHeatmap_CC.png"),res = 600,width = 15,height = 15,units = "in")
GOCluster(circ_CC,GOplotIn_CC$Term) 
dev.off()

png(paste0(output,"CirHeatmap_MF.png"),res = 600,width = 15,height = 15,units = "in")
GOCluster(circ_MF,GOplotIn_MF$Term) 
dev.off()

#### GSEA富集圖
png(paste0(output,"ridgeplot.png"),res = 600,width = 15,height = 15,units = "in")
ridgeplot(GSEA_GO) 
dev.off()

png(paste0(output,"gseaplot.png"),res = 600,width = 15,height = 15,units = "in")
gseaplot2(GSEA_GO,1)
dev.off()

png(paste0(output,"gseaplot_mul.png"),res = 600,width = 10,height = 10,units = "in")
gseaplot2(GSEA_GO,1:30) #30是根據ridgeplot中有30個富集通路得到的
dev.off()


##GO/KEGG/GSEA富集分析圈


