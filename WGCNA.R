
library(WGCNA)

datExpr0 = read.csv("表达谱.csv", row.names = 1)

datExpr0 = as.data.frame(t(datExpr0))
datExpr0[1:4,1:4]
dim(datExpr0)


##为了减少运算量，筛选方差前25%的基因,
##也可以不运行

m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpr0<-data.matrix(expro.upper)
dim(datExpr0)

#主要看缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)

gsg$allOK 
# 返回TRUE则继续
if (!gsg$allOK){
  # 把含有缺失值的基因或样本打印出来
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 去掉那些缺失值
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


##样本过滤

sampleTree = hclust(dist(datExpr0), method = "average")

pdf(file = "1、聚类树.pdf", width = 10, height = 8)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


#如果有异常值就需要去掉，根据聚类图自己设置cutHeight 参数的值

##单独的图
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) +
  #想用哪里切，就把“h = 110”和“cutHeight = 110”中换成你的cutoff
  abline(h = 95, col = "red")

clust = cutreeStatic(sampleTree, cutHeight = 95, minSize = 10)
keepSamples = (clust==1)    #1是保留的
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)

#没有异常样本就不需要去除
datExpr = datExpr0


#输入表型信息

#要求必须是数值型，要么像年龄那样的数字，要么搞成0，1，或者是1，2，3等。
##如果上一步删除了某个样本，这一步一定记得也要抱=把样本删除！！！

datTraits = read.csv("样本信息.csv", row.names = 1)

sampleTree2 = hclust(dist(datExpr), method = "average")
# 各个样本的表现: 白色表示低，红色为高，灰色为缺失
traitColors = numbers2colors(datTraits, signed = FALSE)
# 把样本聚类和表型绘制在一起
pdf(file = "2、样本表现.pdf", width = 10, height = 8)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

##下面是正式的WGCNA分析

##软阈值的筛选
##设置一系列软阈值，范围是1-30之间
powers = c(1:10, seq(from = 12, to=30, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#这个结果就是推荐的软阈值
sft$powerEstimat


cex1 = 0.9 #一般是0.9.不能低于0.85
pdf(file = "3.软阈值的选择.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#进一步构建网络

power = sft$powerEstimate
power#如果前面没有得到推荐的软阈值，就要根据上面的图，自己选择

##需要时间

net = blockwiseModules(datExpr, power = power,
                       TOMType = "unsigned", 
                       minModuleSize = 30, #minModuleSize 默认30，参数设置最小模块的基因数，值越小，小的模块就会被保留下来
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,#mergeCutHeight 默认0.25，设置合并相似性模块的距离，值越小，就越不容易被合并，保留下来的模块就越多。
                       deepSplit = 2 ,#deepSplit 默认2，调整划分模块的敏感度，值越大，越敏感，得到的模块就越多
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "testTOM",
                       verbose = 3)


#此处展示得到了多少模块，每个模块里面有多少基因。
table(net$colors)
##如果结果不合理，可以适当调整参数进行修改

#具体细节可参考：https://zhuanlan.zhihu.com/p/34697561



mergedColors = labels2colors(net$colors)
pdf(file = "4.聚类趋势.pdf", width = 10, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#每种颜色都代表着一种基因模块，即里面的基因功能是相似的
#灰色没意义

#保存每个模块对应的基因
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
gm = data.frame(net$colors)
gm$color = moduleColors
head(gm)

genes = split(rownames(gm),gm$color)
save(genes,file = "genes.Rdata")


#模块与表型的相关性

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#热图
pdf(file = "5.相关性热图.pdf", width = 8, height = 8)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed (50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


#相关系数最好大于0.8，实在没有，0.6或0.7也行，再小就不合适了。


#把gene module输出到文件
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}


#GS与MM
#GS代表模块里的每个基因与形状的相关性
#MM代表每个基因和所在模块之间的相关性，表示是否与模块的趋势一致

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


#第几列的表型是最关心的，下面的i就设置为几。
#与关心的表型相关性最高的模块赋值给下面的module。
traitData = datTraits 
i = 2 #替换自己想要的表型列
module = "black"##替换自己的样本颜色
assign(colnames(traitData)[i],traitData[i])
instrait = eval(parse(text = colnames(traitData)[i]))
geneTraitSignificance = as.data.frame(cor(datExpr, instrait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(instrait), sep="")
names(GSPvalue) = paste("p.GS.", names(instrait), sep="")
pdf(file = paste0("6.模块相关性散点图.pdf"), width = 5, height = 5)
column = match(module, modNames) #找到目标模块所在列
moduleGenes = moduleColors==module #找到模块基因所在行
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##也可以单独提取核心模块的基因做GO分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)

hub_gene <- read.csv('black.csv',header = F)##替换自己的核心模块
hub_gene <- hub_gene[,1]
hub_gene <- hub_gene[2:277]#替换自己的数目
gene=unique(hub_gene)

#转换为ENTREZID
entrezIDs = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Hs.eg.db", drop = TRUE)
#使用entrezIDs 
gene<- entrezIDs$ENTREZID

##GO富集分析
go<- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05,ont="all",readable =T)

write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F) #

##可视化
##条形图
pdf(file="GO-barplot.pdf",width = 10,height = 15)
barplot(go, drop = TRUE, showCategory =10,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

##气泡图
pdf(file="GO-bubble.pdf",width = 10,height = 15)
dotplot(go,showCategory = 10,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#kegg分析
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05, pAdjustMethod = "fdr")   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                         

##可视化
##条形图
pdf(file="KEGG-barplot.pdf",width = 10,height = 13)
barplot(kk, drop = TRUE, showCategory = 15,label_format=100)
dev.off()

##气泡图
pdf(file="KEGG-bubble.pdf",width = 10,height = 13)
dotplot(kk, showCategory = 15,label_format=100)
dev.off()



#也可以通过这里找到GS和MM都大的基因
##大家有兴趣可以自己看一下，一般没见过文章里用这个数据的

f = data.frame(GS = abs(geneModuleMembership[moduleGenes, column]),
               MM = abs(geneTraitSignificance[moduleGenes, 1]))
rownames(f) = rownames(gm[moduleGenes,])
head(f)


##用基因相关性热图的方式展示加权网络

#每行每列代表一个基因。一般取400个基因画就够啦，不然电脑运行不了

nSelect = 400
set.seed(10)
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6)

select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

# 再计算基因之间的距离树(对于基因的子集，需要重新聚类)
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
pdf(file = "7.加权网络.pdf", width = 15, height = 8)
plotDiss = selectTOM^7
diag(plotDiss) = NA #将对角线设成NA，在图形中显示为白色的点，更清晰显示趋势
TOMplot(plotDiss, selectTree, selectColors, col=myheatcol,main = "Network heatmap plot, selected genes")
dev.off()


#模块与表型的相关性
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MET = orderMEs(cbind(MEs, instrait))

pdf(file = "8、模块和表型的相关性.pdf", width = 12, height = 8)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(4,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()

#也可以把上面的图分开来画

pdf(file = "9、模块聚类树.pdf", width = 10, height = 8)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()

pdf(file = "10、表型相关热图.pdf",  width = 10, height = 8)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(4,5,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()























