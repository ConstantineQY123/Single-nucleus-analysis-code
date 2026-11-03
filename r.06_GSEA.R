rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/")
if (! dir.exists("./06_GSEA")){
  dir.create("./06_GSEA")
}
setwd("./06_GSEA")

library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(tidyverse)
library(GseaVis)

load('../01_scRNA/07_subtype/TSNE_Cardiomyocyte_marker.Rdata')
table(sce.car$celltype)
table(sce.car$seurat_clusters)
scRNA <- subset(sce.car, seurat_clusters == 'FADS3_cluster')

train_data <- as.data.frame(scRNA@assays$RNA$counts) %>% apply(2,function(x){x/sum(x) * 10000})
train_data <- train_data %>% as.data.frame()
# train_data <- log2(train_data+1) %>% ScaleData() %>% as.data.frame()

key_gene <- read.csv('../05_expression_verification/hub_gene.csv', row.names = 1, check.names = F)

dat <- train_data
dat$symbol <- rownames(dat)
dat <- separate(dat, symbol, into = c('symbol','var1'), sep = '-')
dat <- dplyr::select(dat, -var1)
dat$symbol <- gsub(' ','',dat$symbol)
dat <- distinct(dat, symbol, .keep_all = T)
rownames(dat) <- dat$symbol
dat <- dplyr::select(dat, -symbol)

i <- 1
for (i in c(1 : length(key_gene$symbol))) {
  setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/06_GSEA/")
  
  train_data <- dat[!grepl(key_gene$symbol[i], rownames(dat)), ] %>% t() %>% as.data.frame()
  key_exp <- dat[key_gene$symbol[i], ] %>% t() %>% as.data.frame()
  
  library(psych)
  library(reshape2)
  
  x <- key_exp
  y <- train_data
  d <- corr.test(x, y, use="complete", method = 'spearman') # 计算相关性
  r <- data.frame(t(d$r)) # 相关性
  colnames(r) <- 'cor_coefficient'
  r$symbol <- rownames(r)
  r <- r[order(r$cor_coefficient, decreasing = T), ]
  r <- dplyr::select(r, -symbol)
  r <- na.omit(r)
  
  geneList <- r$cor_coefficient # 把foldchange按照从大到小提取出来
  names(geneList) <- rownames(r) # 给上面提取的foldchange加上对应上ENTREZID
  geneList
  
  ## GSEA分析
  gmt <- read.gmt("/data/nas1/xuguangxin/pipeline/GSEA/c2.cp.kegg.v7.4.symbols.gmt") #读gmt文件
  res_GSEA <- GSEA(geneList, TERM2GENE = gmt, pvalueCutoff = 1, eps = 0)
  sortGESA <- data.frame(res_GSEA)
  sortGESA <- sortGESA[order(sortGESA$pvalue),]
  sortGESA <- sortGESA[sortGESA$pvalue < 0.05, ]
  if(length(sortGESA$ID) > 5) {
    paths <- rownames(sortGESA[c(1: 5), ]) # 选取需要展示的通路ID（比如前10条）
  }else {
    paths <- sortGESA$ID
  }
  ## 保存GSEA结果
  write.csv(sortGESA, file = paste0("0", i, ".", key_gene$symbol[i], "_gsea_res.csv"))
  
  ## 多通路展示
  p <- gseaNb(object = res_GSEA,
              geneSetID = paths,
              subPlot = 2,
              termWidth = 45,
              # legend.position = c(0.72,0.8),
              addPval = T,
              rmHt = F,
              pvalX = 0.99,
              pvalY = 0.99, 
              newGsea = F,
              base_size = 15,
              lineSize = 2.5,
              curveCol = c('#EE2C2C', "#114889", "#FF6600", "#996633", "#666699"))
  
  pdf(file = paste0("0", i, "_",  key_gene$symbol[i], ".pdf"), w=12, h=8)
  print(p)
  dev.off()
  png(file = paste0("0", i, "_",  key_gene$symbol[i], ".png"), w=12, h=8, units = 'in', res = 600)
  print(p)
  dev.off()
  
  # ## 单通路展示
  # if (! dir.exists(paste0(key_gene$symbol[i]))){
  #   dir.create(paste0(key_gene$symbol[i]))
  # }
  # setwd(paste0(key_gene$symbol[i]))
  # 
  # for (x in c(1 : length(sortGESA$ID))) {
  #   ## 绘图
  #   p <- gseaNb(object = res_GSEA,
  #               geneSetID = sortGESA$ID[x],
  #               subPlot = 2,
  #               termWidth = 35,
  #               legend.position = c(0.72,0.8),
  #               addPval = T,
  #               rmHt = T,
  #               pvalX = 0.99,
  #               pvalY = 0.97, 
  #               newGsea = F,
  #               curveCol = c('#FFDAB9', '#00FFFF', 'pink', 'gray', '#90EE90'))
  #   
  #   pdf(file = paste0(sortGESA$ID[x], ".pdf"), w=13, h=9)
  #   print(p)
  #   dev.off()
  #   png(file = paste0(sortGESA$ID[x], ".png"), w=13, h=9, units = 'in', res = 600)
  #   print(p)
  #   dev.off()
  # }
}
