# 加载必要的包
library(ggvenn)
library(eulerr)
library(scales)
library(dplyr)  

# 读取数据
data <- read.csv("OSAvs减重.csv")
colnames(data) <- c("group1","group2")
# 提取两列基因作为集合
comgene <- list(
  'OSA' = unique(data$group2),
  'bariatric_surgery' = unique(data$group1)
)


# 计算交集基因
com_genes <- intersect(comgene$OSA, comgene$bariatric.surgery)
#保存为 CSV
{
write.csv(data.frame(Intersection_Genes = com_genes),
          "交集基因.csv", row.names = FALSE)
}

#普通韦恩图
p1 <- ggvenn(
  comgene,
  show_percentage = F,
  show_elements = FALSE,
  label_sep = ",",
  digits = 2,
  stroke_color = NA,  # 去除边框
  fill_color = c("#C2E0F7", "#A4CBA8"), 
  set_name_color = c("#C2E0F7", "#A4CBA8"),
  text_color = "black",  # 调整文字颜色
  text_size = 5          # 增大文字
)
print(p1)

# 绘制比例图
p2 <- euler(comgene)
plot(
  p2,
  labels = list(col = "black", font = 1, cex = 1.5),  # 调整字体大小和颜色
  edges = NULL,  # 去除边框
  fills = c("#C210F7", "#14CBA8"),  
  alpha = 0.9,
  quantities = list(cex = 1.5, col = 'black')  # 设置数量文字更突出
)

