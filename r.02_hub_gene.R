rm(list = ls())
setwd('/data/nas1/xuguangxin/project/39_BJTC-630-4/')
if(!dir.exists('./02_hub_gene')){
  dir.create('./02_hub_gene')
}
setwd('./02_hub_gene/')

library(data.table)

gene_A <- fread('../00_rawdata/gene_set/9259297.f12.txt')
colnames(gene_A)[1] <- 'symbol'
gene_A <- dplyr::select(gene_A, symbol)
gene_B <- read.csv('../01_scRNA/07_subtype/FADS3_cluster_DE_sig.csv')
gene_B <- gene_B[gene_B$p_val < 0.05, ]

inter_gene <- intersect(gene_A$symbol, gene_B$X) %>% as.data.frame()
colnames(inter_gene) <- 'symbol'
write.csv(inter_gene, 'hub_gene.csv')

listInput <- list('MMRGs' = gene_A$symbol, 
                  'DEGs' = gene_B$X)

library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)

##可视化结果
p1 <- upset(fromList(listInput),
            number.angles = 0, # 交互集合柱状图的柱标倾角
            point.size = 4, # 图中点的大小
            line.size = 1, # 图中连接线粗细
            mainbar.y.label = "Intersection size", # y轴的标签
            main.bar.color = 'black', # y轴柱状图颜色
            matrix.color = "black", # x轴点的颜色
            sets.x.label = "Set size", # x轴的标签
            sets.bar.color = c('#9328c8','#3d53ba'),
            # sets.bar.color = brewer.pal(3,"Set1"), # x轴柱状图的颜色; Set1中只有9个颜色，Set3中有12个颜色，Paired中有12个颜色
            mb.ratio = c(0.7, 0.3), # bar plot和matrix plot图形高度的占比
            order.by = "freq", # y轴矩阵排序,如"freq"频率，"degree"程度
            text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1), # 6个参数intersection size title（y标题大小）,intersection size tick labels（y刻度标签大小）, set size title（set标题大小）, set size tick labels（set刻度标签大小）, set names（set 分类标签大小）, numbers above bars（柱数字大小）的设置
            shade.color = "#12507B", # 图中阴影部分的颜色
            queries=list(list(query = intersects, params = list("MMRGs"), color = "#9328c8", active = T), # 设置自己想要展示的特定组的交集，通过queries参数进行设置，需要展示几个关注组合的颜色，就展示几个
                         list(query = intersects, params = list("DEGs"), color = "#3d53ba", active = T),
                         list(query = intersects, params = list("MMRGs","DEGs"), color = "orange", active = T))
)
p1
pdf('01.hub_gene.pdf', width = 8, height = 5, family='Times')
print(p1)
dev.off()
png('01.hub_gene.png', width = 8, height = 5, res = 600, units = 'in', family='Times')
print(p1)
dev.off()













