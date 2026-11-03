# 1. GSE36961(芯片数据)---------------------------------------
rm(list = ls())
setwd('D://工作/project/主项目/39_BJTC-630-4肥厚性心肌病（转录组+单细胞）/BJTC-630-4/')
if(!dir.exists('./00_rawdata')){
  dir.create('./00_rawdata')
}
setwd('./00_rawdata/')

library(devtools) 
install_github('GEOquery','seandavi')


### 常规情况（非高通量数据可以在GEO数据库直接用getGEO函数读取，否则要自己手动到网站下载表达矩阵）
library(BiocManager)
library(GEOquery)
library(tidyverse)
library(lance)

GEO_data <- 'GSE36961'
gene_annotation <- 'GPL15389'

if(!dir.exists(paste0(GEO_data))){
  dir.create(paste0(GEO_data))
}
setwd(paste0(GEO_data))

gset <- getGEO(GEO_data , # 前面创建GEO_data对象
               destdir = '.',
               GSEMatrix = T,
               getGPL = F) # getGEO函数自动从官网中获取对应的数据集
# gset$GSE89748_series_matrix.txt.gz@annotation
expr <- exprs(gset[[1]]) # 将获取的数据集的表达矩阵提取出来，如果没有表达矩阵，说明是高通量数据需要到网站自己下载
qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
  expr[which(expr <= 0)] <- NaN
  expr <- log2(expr + 0.1) 
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}
expr <- as.data.frame(expr)

gpl <- getGEO(gene_annotation, destdir = '.') 
gpl <- Table(gpl)
colnames(gpl)
probe2symbol <- dplyr::select(gpl, 'ID', 'Symbol')
probe2symbol <- filter(probe2symbol, `Symbol` != '')
probe2symbol <- separate(probe2symbol, `Symbol`, into = c('symbol', 'drop'), sep = '//')
probe2symbol <- dplyr::select(probe2symbol, -drop)
names(probe2symbol) <- c('ID', 'symbol')
probe2symbol <- probe2symbol[probe2symbol$symbol != ' --- ', ]
probe2symbol$symbol <- gsub(' ', '', probe2symbol$symbol)

dat <- expr
dat$ID <- rownames(dat)
dat$ID <- as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)

dat <- dat %>%
  merge(probe2symbol, by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol, everything())%>%     ## 重新排列
  mutate(rowMean = rowMeans(.[grep('GSM', names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol, .keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

# protein_gene <- read.delim2('/data/nas1/xuguangxin/pipeline/PCG/PCG.xls(v22)')
# dat <- dat[rownames(dat) %in% protein_gene$gene_name, ]

a <- gset[[1]]
pd <- pData(a)
# write.csv(pd, 'pd.csv')
table(pd$characteristics_ch1)
table(pd$source_name_ch1)
# pd1 <- pd[pd$source_name_ch1 %in% 'Diffuse Large B cell Lymphoma patient', ]
# pd2 <- pd[pd$source_name_ch1 %in% 'Lymph-node', ]
# pd <- rbind(pd1, pd2)
group <- data.frame(sample = pd$geo_accession, group = pd$characteristics_ch1)
table(group$group)
group <- group[group$group %in% c('category: Control', 'category: KD'), ]
group$group <- ifelse(group$group == 'category: Control', 'control', 'disease')
table(group$group) 

dat <- dat[, colnames(dat) %in% group$sample]
dat <- dat[, colnames(dat) != 'GSM1894671']
dat <- na.omit(dat)
group <- group[group$sample %in% colnames(dat), ]

write.csv(dat, file = paste0('dat.', GEO_data, '.csv'))
write.csv(group, file = paste0('group.', GEO_data, '.csv'), row.names = F)