rm(list = ls())
setwd('/data/nas1/xuguangxin/project/39_BJTC-630-4/')
if(!dir.exists('./04_PPI')){
  dir.create('./04_PPI')
}
setwd('./04_PPI/')

library(tidyverse)

cytoHubba_res <- read_csv('./cytohubba_res(MCODE).csv')

MCC <- cytoHubba_res[order(cytoHubba_res$MCC, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
DMNC <- cytoHubba_res[order(cytoHubba_res$DMNC, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
MNC <- cytoHubba_res[order(cytoHubba_res$MNC, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
Degree <- cytoHubba_res[order(cytoHubba_res$Degree, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
# EPC <- cytoHubba_res[order(cytoHubba_res$EPC, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
BottleNeck <- cytoHubba_res[order(cytoHubba_res$BottleNeck, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
EcCentricity <- cytoHubba_res[order(cytoHubba_res$EcCentricity, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
Closeness <- cytoHubba_res[order(cytoHubba_res$Closeness, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
Radiality <- cytoHubba_res[order(cytoHubba_res$Radiality, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
Betweenness <- cytoHubba_res[order(cytoHubba_res$Betweenness, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
Stress <- cytoHubba_res[order(cytoHubba_res$Stress, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)
# ClusteringCoefficient <- cytoHubba_res[order(cytoHubba_res$ClusteringCoefficient, decreasing = T), ][c(1:10), ] %>% dplyr::select(node_name)


inter_gene1 <- intersect(MCC$node_name, DMNC$node_name)
inter_gene2 <- intersect(inter_gene1, MNC$node_name)
inter_gene3 <- intersect(inter_gene2, Degree$node_name)
inter_gene4 <- intersect(inter_gene3, BottleNeck$node_name)
inter_gene5 <- intersect(inter_gene4, EcCentricity$node_name)
inter_gene6 <- intersect(inter_gene5, Closeness$node_name)
inter_gene7 <- intersect(inter_gene6, Radiality$node_name)
inter_gene8 <- intersect(inter_gene7, Betweenness$node_name)
inter_gene9 <- intersect(inter_gene8, Stress$node_name) %>% as.data.frame()
# inter_gene10 <- intersect(inter_gene9, ClusteringCoefficient$node_name)
# inter_gene11 <- intersect(inter_gene10, ClusteringCoefficient$node_name)
colnames(inter_gene9) <- 'symbol'
write.csv(inter_gene9, file = 'hub_gene.csv')

library(rio)
library(plotrix)

## 2.读取数据
colnames(cytoHubba_res)
venn_data <- cbind(MCC,DMNC,MNC,Degree,BottleNeck,EcCentricity,Closeness,Radiality,Betweenness,Stress)
colnames(venn_data) <- c("MCC", 'DMNC', 'MNC', 'Degree', 'BottleNeck', 'EcCentricity', 'Closeness', 'Radiality', 'Betweenness', 'Stress')

sample_id <- colnames(venn_data)

otu_id <- unique(venn_data[,1])

otu_id <- otu_id[otu_id != '']

core_otu_id <- otu_id

otu_num <- length(otu_id)

## 3.遍历数据文件获取数据
i <- 2
for (i in 2:ncol(venn_data)) {
  
  otu_id <- unique(venn_data[,i])
  
  otu_id <- otu_id[otu_id != '']
  
  core_otu_id <- intersect(core_otu_id, otu_id)
  
  otu_num <- c(otu_num, length(otu_id))
  
}

core_num <- length(core_otu_id)

## 4.定义备选颜色（示例数据20个样本，所以我设置了20种颜色）
ellipse_col <- c('#6181BD4E','#F348004E','#64A10E4E','#9300264E','#464E044E','#049a0b4E','#4E0C664E','#D000004E','#FF6C004E','#FF00FF4E','#c7475b4E','#00F5FF4E','#BDA5004E','#A5CFED4E','#f0301c4E','#2B8BC34E','#FDA1004E','#54adf54E','#CDD7E24E','#9295C14E')

## 5.创建绘图函数：flower_plot()

flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
  
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
  
  plot(c(0,10),c(0,10),type='n')
  
  n  <- length(sample)
  
  deg <- 360 / n
  
  res <- lapply(1:n, function(t){
    
    draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                 
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                 
                 col = ellipse_col[t],
                 
                 border = ellipse_col[t],
                 
                 a = a, b = b, angle = deg * (t - 1))
    
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         
         otu_num[t])
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           
           sample[t],
           
           srt = deg * (t - 1) - start,
           
           adj = 1,
           
           cex = 1
           
      )
      
    } else {
      
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           
           sample[t],
           
           srt = deg * (t - 1) + start,
           
           adj = 0,
           
           cex = 1
           
      )
      
    }
    
  })
  
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
  
  text(x = 5, y = 5, paste('Core:', core_otu))
  
}

## 6.调用上述函数作图（保存图片到本地文件夹中）
png('01.hub_gene.png', width = 8, height = 8, res = 600, units = 'in', family='Times')

flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
            
            start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')

dev.off()

pdf('01.hub_gene.pdf', width = 8, height = 8, family='Times')

flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
            
            start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')

dev.off()

