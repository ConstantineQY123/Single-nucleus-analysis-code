rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/")
if (! dir.exists("./08_disease")){
  dir.create("./08_disease")
}
setwd("./08_disease")

## COX6C
COX6C <- read_tsv('./COX6C.tsv')
COX6C <- COX6C[order(COX6C$`Inference Score`, decreasing = T), ]
COX6C <- COX6C[c(1:10), ]
COX6C <- dplyr::select(COX6C, 'Disease Name', 'Inference Score', 'Reference Count')
colnames(COX6C) <- c('disease', 'Inference_Score', 'Reference_Count')
COX6C$disease <- gsub(' ', '_', COX6C$disease)
COX6C <- column_to_rownames(COX6C, var = 'disease')
COX6C <- t(COX6C) %>% as.data.frame()

max(COX6C)
colnames(COX6C)
max_min <- data.frame(
  Chemical_and_Drug_Induced_Liver_Injury = c(1200, 0), 
  Prenatal_Exposure_Delayed_Effects = c(1200, 0), 
  Weight_Loss = c(1200, 0),
  Hepatomegaly = c(1200, 0), 
  Hyperplasia = c(1200, 0), 
  Edema = c(1200, 0),
  Necrosis = c(1200, 0), 
  Fatty_Liver = c(1200, 0), 
  Inflammation = c(1200, 0),
  Kidney_Diseases = c(1200, 0)
)
rownames(max_min) <- c("Max", "Min")

# 合并数据
df <- rbind(max_min, COX6C)

library(fmsb)
pdf(file = paste0("01.COX6C.pdf"), w=12, h=8)
radarchart(df,
           axistype = 1,
           pcol = c("#00AFBB","#E7B800"), # 线条颜色
           pfcol = scales::alpha(c("#00AFBB", "#E7B800"),0.5), # 填充颜色
           plwd = 2, # 线条宽度
           plty = 1, # 线条类型. 使用数字1-6或者字符向量 (“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”).使用 plty = 0 or plty = “blank”，不显示线条
           cglcol = "grey", # 网格线颜色
           cglty = 1, # 网格线类型
           cglwd = 0.8, # 网格线宽度
           axislabcol = "grey", # 轴标签和数字的颜色。默认设置是“蓝色”。
           vlcex = 1.1, # 变量标签字体大小
           vlabels = colnames(df), # 变量标签
           caxislabels = c(0, 300, 600, 900, 1200), # 字符向量，用作中心轴上的标签。
           title = 'COX6C'
)
legend(
  x = "bottom", 
  legend = rownames(df[-c(1,2),]), 
  horiz = TRUE,
  bty = 'n', 
  pch = 19 , 
  col = c("#00AFBB", "#E7B800"),
  cex = 1.1, 
  pt.cex = 3,
  # x.intersp = -2,
  y.intersp = -1.5,
  # adj = -0.1,
  xpd = TRUE,
  text.col = c("#00AFBB", "#E7B800"),
  text.width = 1
)
dev.off()
png(file = paste0("01.COX6C.png"), w=12, h=8, units = 'in', res = 600)
radarchart(df,
           axistype = 1,
           pcol = c("#00AFBB","#E7B800"), # 线条颜色
           pfcol = scales::alpha(c("#00AFBB", "#E7B800"),0.5), # 填充颜色
           plwd = 2, # 线条宽度
           plty = 1, # 线条类型. 使用数字1-6或者字符向量 (“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”).使用 plty = 0 or plty = “blank”，不显示线条
           cglcol = "grey", # 网格线颜色
           cglty = 1, # 网格线类型
           cglwd = 0.8, # 网格线宽度
           axislabcol = "grey", # 轴标签和数字的颜色。默认设置是“蓝色”。
           vlcex = 1.1, # 变量标签字体大小
           vlabels = colnames(df), # 变量标签
           caxislabels = c(0, 300, 600, 900, 1200), # 字符向量，用作中心轴上的标签。
           title = 'COX6C'
)
legend(
  x = "bottom", 
  legend = rownames(df[-c(1,2),]), 
  horiz = TRUE,
  bty = 'n', 
  pch = 19 , 
  col = c("#00AFBB", "#E7B800"),
  cex = 1.1, 
  pt.cex = 3,
  # x.intersp = -2,
  y.intersp = -1.5,
  # adj = -0.1,
  xpd = TRUE,
  text.col = c("#00AFBB", "#E7B800"),
  text.width = 1
)
dev.off()

## NDUFA8
NDUFA8 <- read_tsv('./NDUFA8.tsv')
NDUFA8 <- NDUFA8[order(NDUFA8$`Inference Score`, decreasing = T), ]
NDUFA8 <- NDUFA8[c(1:10), ]
NDUFA8 <- dplyr::select(NDUFA8, 'Disease Name', 'Inference Score', 'Reference Count')
colnames(NDUFA8) <- c('disease', 'Inference_Score', 'Reference_Count')
NDUFA8$disease <- gsub(' ', '_', NDUFA8$disease)
NDUFA8 <- column_to_rownames(NDUFA8, var = 'disease')
NDUFA8 <- t(NDUFA8) %>% as.data.frame()

max(NDUFA8)
colnames(NDUFA8)
max_min <- data.frame(
  Weight_Loss = c(1200, 0), 
  Necrosis = c(1200, 0), 
  Chemical_and_Drug_Induced_Liver_Injury = c(1200, 0),
  Inflammation = c(1200, 0), 
  Kidney_Diseases = c(1200, 0), 
  Testicular_Diseases = c(1200, 0),
  Hyperplasia = c(1200, 0), 
  Poisoning = c(1200, 0), 
  Fibrosis = c(1200, 0),
  Hepatomegaly = c(1200, 0)
)
rownames(max_min) <- c("Max", "Min")

# 合并数据
df <- rbind(max_min, NDUFA8)

library(fmsb)
pdf(file = paste0("02.NDUFA8.pdf"), w=12, h=8)
radarchart(df,
           axistype = 1,
           pcol = c("#00AFBB","#E7B800"), # 线条颜色
           pfcol = scales::alpha(c("#00AFBB", "#E7B800"),0.5), # 填充颜色
           plwd = 2, # 线条宽度
           plty = 1, # 线条类型. 使用数字1-6或者字符向量 (“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”).使用 plty = 0 or plty = “blank”，不显示线条
           cglcol = "grey", # 网格线颜色
           cglty = 1, # 网格线类型
           cglwd = 0.8, # 网格线宽度
           axislabcol = "grey", # 轴标签和数字的颜色。默认设置是“蓝色”。
           vlcex = 1.1, # 变量标签字体大小
           vlabels = colnames(df), # 变量标签
           caxislabels = c(0, 300, 600, 900, 1200), # 字符向量，用作中心轴上的标签。
           title = 'NDUFA8'
)
legend(
  x = "bottom", 
  legend = rownames(df[-c(1,2),]), 
  horiz = TRUE,
  bty = 'n', 
  pch = 19 , 
  col = c("#00AFBB", "#E7B800"),
  cex = 1.1, 
  pt.cex = 3,
  # x.intersp = -2,
  y.intersp = -1.5,
  # adj = -0.1,
  xpd = TRUE,
  text.col = c("#00AFBB", "#E7B800"),
  text.width = 1
)
dev.off()
png(file = paste0("02.NDUFA8.png"), w=12, h=8, units = 'in', res = 600)
radarchart(df,
           axistype = 1,
           pcol = c("#00AFBB","#E7B800"), # 线条颜色
           pfcol = scales::alpha(c("#00AFBB", "#E7B800"),0.5), # 填充颜色
           plwd = 2, # 线条宽度
           plty = 1, # 线条类型. 使用数字1-6或者字符向量 (“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”).使用 plty = 0 or plty = “blank”，不显示线条
           cglcol = "grey", # 网格线颜色
           cglty = 1, # 网格线类型
           cglwd = 0.8, # 网格线宽度
           axislabcol = "grey", # 轴标签和数字的颜色。默认设置是“蓝色”。
           vlcex = 1.1, # 变量标签字体大小
           vlabels = colnames(df), # 变量标签
           caxislabels = c(0, 300, 600, 900, 1200), # 字符向量，用作中心轴上的标签。
           title = 'NDUFA8'
)
legend(
  x = "bottom", 
  legend = rownames(df[-c(1,2),]), 
  horiz = TRUE,
  bty = 'n', 
  pch = 19 , 
  col = c("#00AFBB", "#E7B800"),
  cex = 1.1, 
  pt.cex = 3,
  # x.intersp = -2,
  y.intersp = -1.5,
  # adj = -0.1,
  xpd = TRUE,
  text.col = c("#00AFBB", "#E7B800"),
  text.width = 1
)
dev.off()

## NDUFB1
NDUFB1 <- read_tsv('./NDUFB1.tsv')
NDUFB1 <- NDUFB1[order(NDUFB1$`Inference Score`, decreasing = T), ]
NDUFB1 <- NDUFB1[c(1:10), ]
NDUFB1 <- dplyr::select(NDUFB1, 'Disease Name', 'Inference Score', 'Reference Count')
colnames(NDUFB1) <- c('disease', 'Inference_Score', 'Reference_Count')
NDUFB1$disease <- gsub(' ', '_', NDUFB1$disease)
NDUFB1 <- column_to_rownames(NDUFB1, var = 'disease')
NDUFB1 <- t(NDUFB1) %>% as.data.frame()

max(NDUFB1)
colnames(NDUFB1)
max_min <- data.frame(
  Weight_Loss = c(1200, 0), 
  Chemical_and_Drug_Induced_Liver_Injury = c(1200, 0), 
  Prenatal_Exposure_Delayed_Effects = c(1200, 0),
  Hyperplasia = c(1200, 0), 
  Necrosis = c(1200, 0), 
  Neurotoxicity_Syndromes = c(1200, 0),
  Inflammation = c(1200, 0), 
  Hepatomegaly = c(1200, 0), 
  Fibrosis = c(1200, 0),
  Heart_Diseases = c(1200, 0)
)
rownames(max_min) <- c("Max", "Min")

# 合并数据
df <- rbind(max_min, NDUFB1)

library(fmsb)
pdf(file = paste0("03.NDUFB1.pdf"), w=12, h=8)
radarchart(df,
           axistype = 1,
           pcol = c("#00AFBB","#E7B800"), # 线条颜色
           pfcol = scales::alpha(c("#00AFBB", "#E7B800"),0.5), # 填充颜色
           plwd = 2, # 线条宽度
           plty = 1, # 线条类型. 使用数字1-6或者字符向量 (“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”).使用 plty = 0 or plty = “blank”，不显示线条
           cglcol = "grey", # 网格线颜色
           cglty = 1, # 网格线类型
           cglwd = 0.8, # 网格线宽度
           axislabcol = "grey", # 轴标签和数字的颜色。默认设置是“蓝色”。
           vlcex = 1.1, # 变量标签字体大小
           vlabels = colnames(df), # 变量标签
           caxislabels = c(0, 300, 600, 900, 1200), # 字符向量，用作中心轴上的标签。
           title = 'NDUFB1'
)
legend(
  x = "bottom", 
  legend = rownames(df[-c(1,2),]), 
  horiz = TRUE,
  bty = 'n', 
  pch = 19 , 
  col = c("#00AFBB", "#E7B800"),
  cex = 1.1, 
  pt.cex = 3,
  # x.intersp = -2,
  y.intersp = -1.5,
  # adj = -0.1,
  xpd = TRUE,
  text.col = c("#00AFBB", "#E7B800"),
  text.width = 1
)
dev.off()
png(file = paste0("03.NDUFB1.png"), w=12, h=8, units = 'in', res = 600)
radarchart(df,
           axistype = 1,
           pcol = c("#00AFBB","#E7B800"), # 线条颜色
           pfcol = scales::alpha(c("#00AFBB", "#E7B800"),0.5), # 填充颜色
           plwd = 2, # 线条宽度
           plty = 1, # 线条类型. 使用数字1-6或者字符向量 (“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”).使用 plty = 0 or plty = “blank”，不显示线条
           cglcol = "grey", # 网格线颜色
           cglty = 1, # 网格线类型
           cglwd = 0.8, # 网格线宽度
           axislabcol = "grey", # 轴标签和数字的颜色。默认设置是“蓝色”。
           vlcex = 1.1, # 变量标签字体大小
           vlabels = colnames(df), # 变量标签
           caxislabels = c(0, 300, 600, 900, 1200), # 字符向量，用作中心轴上的标签。
           title = 'NDUFB1'
)
legend(
  x = "bottom", 
  legend = rownames(df[-c(1,2),]), 
  horiz = TRUE,
  bty = 'n', 
  pch = 19 , 
  col = c("#00AFBB", "#E7B800"),
  cex = 1.1, 
  pt.cex = 3,
  # x.intersp = -2,
  y.intersp = -1.5,
  # adj = -0.1,
  xpd = TRUE,
  text.col = c("#00AFBB", "#E7B800"),
  text.width = 1
)
dev.off()

