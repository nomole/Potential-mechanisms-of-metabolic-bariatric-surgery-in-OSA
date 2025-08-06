# 加载所需的R包
library(limma)        # limma包用于差异分析
library(reshape2)     # reshape2包用于数据重塑，转换数据格式
library(tidyverse)    # tidyverse包包含数据处理和可视化的多个工具
library(ggplot2)      # ggplot2包用于绘制高质量的图形
library(ggpubr)       # ggpubr包用于ggplot2的扩展，方便绘图和统计分析
# install.packages("devtools")  # 如果没有安装devtools包，可以解开注释并运行
# devtools::install_github("daattali/ggExtra")  # 从GitHub安装ggExtra包
library(ggExtra) # ggExtra包用于添加额外的图形元素（如边际分布）
setwd("D:/loose weight/GSE135917/扩大wgcna/免疫浸润2")
getwd()
# 读取表达矩阵和分组信息文件
immune <- read.csv("cibersort的结果.csv", row.names = 1)  # 读取免疫细胞成分表达数据，第一列是行名（样本名）
risk <- read.csv("分组文件.csv", row.names = 1)           # 读取分组信息文件，第一列是样本名


# 过滤表达矩阵，保留P值小于0.05的基因数据
immune = immune[immune[,"P.value"] < 0.05,]  # 过滤掉P值大于等于0.05的基因

# 提取表达数据矩阵，去掉最后三列
data = as.matrix(immune[, 1:(ncol(immune) - 3)])  

# 找到分组信息和表达数据中共有的样本
sameSample = intersect(row.names(data), row.names(risk))  # 找出表达矩阵和分组信息中相同的样本
data = data[sameSample, , drop = F]  # 根据共有样本筛选表达矩阵
risk = risk[sameSample, , drop = F]  # 根据共有样本筛选分组信息

#下面示例数据为基因“KRT18”，换成自己的基因即可，也可改为风险评分！！

for (i in colnames(data)[1:ncol(data)]) {  # 遍历表达数据的每一列（基因或通路的表达量）
  x = as.numeric(risk[, "IL1RN"])  # 提取分组信息中KRT18的数值列
  x[x > quantile(x, 0.99)] = quantile(x, 0.99)  # 将KRT18列中大于99百分位的值设为99百分位（去除极端值）
  
  y = as.numeric(data[, i])  # 提取当前基因（或通路）在表达矩阵中的数据
  
  if (sd(y) < 0.01) { next }  # 如果当前基因的标准差小于0.01，跳过此基因的分析（表示表达量几乎没有变化）
  
  cor = cor.test(x, y, method = "spearman")  # 对KRT18和当前基因表达数据进行Spearman秩相关检验
  
  if (cor$p.value < 0.05) {  # 如果Spearman相关检验的P值小于0.05（即相关性显著）
    outFile = paste0(i, ".pdf")  # 创i建输出文件名
    
    # 将KRT18和当前基因表达数据合并成一个数据框，用于绘图
    df1 = as.data.frame(cbind(x, y))
    
    # 使用ggplot2绘制散点图，x轴为KRT18表达量，y轴为当前基因的表达量
    p1 = ggplot(df1, aes(x, y)) +
      xlab("IL1RN") + ylab(i) +  # 设置x轴和y轴的标签
      geom_point(color = "black") +  # 绘制黑色散点
      geom_smooth(method = "lm", formula = y ~ x, color = "blue") +  # 添加蓝色回归线（线性回归）
      theme_bw() +  # 使用黑白主题
      stat_cor(method = 'spearman', aes(x = x, y = y))  # 在图中添加Spearman相关性系数及其P值
    
    # 使用ggMarginal函数添加边际密度图，灰色边际密度
    p2 = ggMarginal(p1, type = "density", xparams = list(fill = "darkolivegreen3", alpha = 0.5), yparams = list(fill = "burlywood4", alpha = 0.5))
    
    # 将图形保存为PDF文件
    pdf(file = outFile, width = 5.2, height = 5)  # 设置输出PDF文件的宽度和高度
    print(p2)  # 打印（绘制）图形
    dev.off()  # 关闭PDF设备，保存图形
  }
}

