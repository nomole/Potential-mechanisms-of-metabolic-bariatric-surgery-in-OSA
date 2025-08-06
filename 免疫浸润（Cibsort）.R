# 导入所需的R包
library(dplyr)          # 数据处理包
library(CIBERSORT)      # 用于免疫细胞浸润分析
library(ggplot2)        # 绘制数据图表
library(pheatmap)       # 绘制热图
library(ggpubr)         # 提供ggplot2增强功能，便于绘图和统计
library(reshape2)      
library(tidyverse)      
library(limma)          # 用于差异表达分析
library(ggthemes)       
library(RColorBrewer)   # 提供色彩调色板

expr <- expr[1:23279,]
write.csv(expr,file = "表达谱2.csv")
# 读取表达矩阵和分组信息
expr <- read.csv("表达谱2.csv", row.names = 1)  # 读取基因表达矩阵，第一列是基因ID
group <- read.csv("分组.csv", row.names = 1)          # 读取分组信息（例如高、低两组）

## CIBERSORT对数据格式有严格的要求：
# (1) 不能有负值或缺失值；(2) 尽量标准化数据
boxplot(expr, outline = FALSE, notch = F, las = 2)  # 绘制数据箱型图，检查数据分布
expr = normalizeBetweenArrays(expr)  # 数据标准化，消除批次效应
boxplot(expr, outline = FALSE, notch = F, las = 2)  # 标准化后的数据箱型图

#免疫细胞类型（官网下载）
LM22 <- read.table("LM22.txt", header = TRUE, row.names = 1, sep = "\t") %>% as.matrix()

# 开始CIBERSORT分析，时间较长
# 计算出来的结果包含22种免疫细胞的丰度，最后三列为其他统计量，可以不管它们。
# perm：表示置换次数，数字越大运行时间越长，一般设置为1000
# QN：芯片数据设置为TRUE，测序数据设置为FALSE
result <- cibersort(sig_matrix = LM22, mixture_file = expr, perm = 2000, QN = TRUE) %>% as.data.frame()
result1 <-result 
# 保存CIBERSORT分析的结果到CSV文件
write.csv(result, "cibersort的结果.csv")

# 只保留22个免疫细胞类型的丰度数据，去除最后三列统计信息
result <- result[, 1:22]

# 绘制免疫细胞类型相关性的散点图
#也可以放在最后，只绘制有统计意义的免疫细胞
library(corrplot)  
cor <- cor(result)  # 计算免疫细胞丰度的相关性
corplot <- corrplot(cor, method = 'number')  # 显示相关性数值
# 另一种布局r
corrplot(cor, order = "AOE", type = "upper", tl.pos = "d")

resulta <- result[,c(-10,-20,-21)]

# 合并免疫细胞丰度数据与分组信息
data1 <- cbind(result, group)  # 将免疫细胞丰度和分组信息合并
data1 <- data1 %>% rownames_to_column("sample")  # 将行名（样本ID）转换为一列
data <- gather(data1, key = CIBERSORT, value = Proportion, -c(group, sample))  # 转换数据格式，便于绘制图表

## 绘制分组差异的箱型图

{
 FENZU <-  ggboxplot(data, x = "CIBERSORT", y = "Proportion",  # 绘制免疫细胞比例的箱型图,四组使用"#FFFF00","#6699FF"
            fill = "group", palette = c("#ED5462", "#81D5B0")) +  # 设置分组颜色
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +  # x轴标签旋转
    stat_compare_means(label = "p.signif", method = "t.test",  # 进行t检验并标注p值
                       ref.group = ".all.", hide.ns = TRUE)  # 比较所有组之间的差异
}

# 绘制免疫细胞丰度的热图
{
k <- apply(result, 2, function(x) { sum(x == 0) < nrow(result1) / 2 })  # 筛选出非零的免疫细胞类型
re2 <- as.data.frame(t(result[, k]))  # 转置数据以便绘制热图
an = data.frame(group = group, row.names = colnames(expr))  # 创建列注释数据
}

heatmap <- pheatmap(re2, scale = "row", show_colnames = FALSE, cluster_cols = FALSE,  # 绘制热图
         annotation_col = an, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))  # 设置热图颜色
dev.new()
# 提取单个免疫细胞（如CD8 T细胞）的差异
dat2 = data[data$CIBERSORT == "Macrophages.M2",]  # 提取CD8 T细胞的比例数据

{
dange <-   ggplot(dat2, aes(x = group, y = Proportion)) +  # 绘制CD8 T细胞的箱型图
    labs(y = "Cell composition", x = NULL, title = "") +
    geom_boxplot(aes(fill = group), position = position_dodge(0.5), width = 0.5, size = 0.4,
                 outlier.alpha = 1, outlier.size = 0.5) +  # 设置箱型图的样式
    theme_bw() + 
    scale_fill_manual(values = c("#EB7369", "#1CB4B8","navy","#81D5B0")) +  # 设置组的颜色
    stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test",  # 使用Wilcoxon检验
                       size = 6, hide.ns = TRUE)
}

# 绘制堆积比例图，展示每个样本的免疫细胞比例
{
  identical(rownames(result), group$Samples)  # 确认免疫细胞数据和样本分组信息匹配
  data <- cbind(rownames(result), result)  # 合并样本名称和免疫细胞丰度数据
  colnames(data)[1] <- "Samples"  # 重命名合并后的第一列为“Samples”
  data <- melt(data, id.vars = c("Samples"))  # 将数据从宽格式转换为长格式
  colnames(data) <- c('Samples', 'celltype', 'proportion')  # 设置列名
}

# 绘制堆积比例图
{
colour = c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"))  # 设置颜色
A <- ggplot(data, aes(Samples, proportion, fill = celltype)) +  # 绘制堆积条形图
  geom_bar(stat = "identity", position = "fill") +  # 绘制堆积比例图
  scale_fill_manual(values = colour) +  # 设置颜色
  ggtitle("") + theme_gray() + theme(axis.ticks.length = unit(3, 'mm'), axis.title.x = element_text(size = 11)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  guides(fill = guide_legend(title = ""))
}

# 另一种表现形式的堆积图

options(repr.plot.width = 20, repr.plot.height = 20)

{
p <- ggplot(data, aes(x = Samples, y = proportion, fill = celltype)) +  # 绘制堆积条形图
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colour) +  # 设置颜色
  labs(x = "", y = "", title = "") +  # 设置标题和坐标轴标签
  scale_y_continuous(expand = c(0, 0)) +  # 去除y轴的空白
  guides(fill = guide_legend(ncol = 1)) +  # 设置图例排列
  theme_bw() + theme(legend.key = element_blank(), legend.title = element_blank(),
                     panel.grid = element_blank(), axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5))
}
dev.new()

ggsave("各种免疫细胞相关性.pdf", plot = A, width = 30, height = 16, units = "cm")
ggsave("堆积比例2.PDF", plot = p, width = 60, height = 32, units = "cm")



dange + stat_compare_means(
  method = "t.test",          # 检验方法
  comparisons = list(          # 指定要比较的组
    c("Control", "OSA")
    
  ),
  label = "p.signif",          # 显示符号（*表示显著性）
  method.args = list(p.adjust.method = "bonferroni") # 多重检验校正
)
