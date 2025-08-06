library(clusterProfiler)#GO富集分析、KEGG通路富集分析
library(org.Hs.eg.db)#基因注释数据库
library(enrichplot)
library(ggplot2)
library(GOplot)
#设置工作路径
setwd("D:/loose weight/GSE135917/扩大wgcna/KEGG_GO")
#读入数据
genes_df <- read.csv("DEGs.csv") 
genes_df <- as.data.frame(genes_df)
# 将基因符号转换为ENTREZ ID
entrezIDs <- bitr(genes_df$DEGs, fromType = "SYMBOL", 
               toType = c("ENTREZID", "SYMBOL"),
               OrgDb = org.Hs.eg.db) # `OrgDb`参数指定使用的人类基因组注释数据库
#使用entrezIDs 
gene<- entrezIDs$ENTREZID
##GO富集分析
go<- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05,ont="all",readable =T)
write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F) #
GOresult <- go@result
##可视化
##条形图
pdf(file="GO-柱状图.pdf",width = 8,height = 7)
##showCategory改变自己想要展示的条目数量
barplot(go, drop = TRUE, showCategory =16,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

##气泡图
pdf(file="GO-气泡图.pdf",width = 8,height = 7)
dotplot(go,showCategory = 6,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()


#kegg分析
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05, pAdjustMethod = "fdr")   
KEGGresult <- kk@result                      
result <- as.data.frame(kk@result)
write.table(result,file="KEGG.txt",sep="\t",quote=F,row.names = F)  
##可视化
##条形图
pdf(file="KEGG-柱状图.pdf",width = 8,height =6)
barplot(kk, drop = TRUE, showCategory = 16,label_format=100)
dev.off()

##气泡图
pdf(file="KEGG-气泡图.pdf",width = 8,height = 6)
dotplot(kk, showCategory = 16,label_format=100)
dev.off()

options(timeout = 300)

class(result)
?barplot
