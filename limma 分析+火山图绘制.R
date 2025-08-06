##免费的午餐limma+火山图代码##
exp <- exp[1:23279,]

#设置环境变量
library(tidyverse)
library(limma)
exp<-read.csv("表达谱2.csv",row.names = 1)#读入表达矩阵文件
exp <- as.data.frame(t(exp))
##数据清洗###
##并非一定要运行，有需要再运行###

#如果表达值超过15，需要运行下面的代码------------------------------
exp=log2(exp+1)  

# 均数不一致需标准化 ------------------------------
boxplot(exp,outline=FALSE, notch=F , las=2)
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=F , las=2)

##去除低表达的基因 ------------------------------

exp = exp[rowMeans(exp)>1,] #根据自己的需要去除低表达基因，也可以卡其它阈值

#读入分类文件------------------------------

fen<-read.csv("样本信息.csv",row.names = 1)#读入表达矩阵表头的分类文件

##开始进行差异分析-----------------------------

group_list <- factor(fen$group,levels = c("Control","OSA"))
design <- model.matrix(~group_list)
con <- lmFit(exp,design)
con2 <- eBayes(con)
DEG1 <- topTable(con2, coef = 2, number = Inf)
DEG2 = na.omit(DEG1) 

##导出差异分析结果-----------------------------

write.csv(DEG2,"差异分析结果（新）.csv")


#筛选差异基因----------------------------

logFC_cut = 0.585
p_cut = 0.05

type1 = (DEG2$adj.P.Val < p_cut)&(DEG2$logFC < -logFC_cut)
type2 = (DEG2$adj.P.Val < p_cut)&(DEG2$logFC > logFC_cut)

DEG2$type = ifelse(type1,"Down",ifelse(type2,"Up","NOT"))

head(DEG2)

write.csv(DEG2,"差异分析标注（新）.csv")


##绘制火山图-------------------带标记基因---------


library(ggplot2)

dev.new()

# 新增一列用于存储label信息，将需要显示的label列出即可：
#DEG2$label <- NA
#DEG2$label[rownames(DEG2) %in% c("HIST1H1E")] <- "HIST1H1E" # 如果HIST1H1E存在





#离散绘图

# 添加分类列
DEG2 <- DEG2 %>%
  mutate(
    Expression = case_when(
      adj.P.Val >= 0.05 ~ "NOT",
      logFC >= 0.585 & adj.P.Val < 0.05 ~ "UP",
      logFC <= -0.585 & adj.P.Val < 0.05 ~ "DOWN",
      TRUE ~ "NOT"
    )
  )

# 绘图
ggplot(DEG2, aes(logFC, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(adj.P.Val), color = Expression)) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NOT" = "grey")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank())



#绘制热图
library(pheatmap)
DEG_genes <- DEG2[DEG2$adj.P.Val<0.05&abs(DEG2$logFC)>0.5,]
DEG_gene_expr <- exp[rownames(DEG_genes),]
#DEG_gene_expr[is.infinite(DEG_gene_expr)] = 0
#DEG_gene_expr[DEG_gene_expr == -Inf] = 0
pdf(paste0(job,"_","pheatmap.pdf"))
pheatmap(DEG_gene_expr,
         color = colorRampPalette(c("blue","white","red"))(100), #颜色
         scale = "row", #归一化的方式
         border_color = NA, #线的颜色
         fontsize = 10, #文字大小
         show_rownames = F)
dev.off()
