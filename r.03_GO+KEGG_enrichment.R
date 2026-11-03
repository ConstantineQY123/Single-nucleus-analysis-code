# 1. GO富集分析---------------
rm(list = ls())
setwd('/data/nas1/xuguangxin/project/39_BJTC-630-4/')
if(!dir.exists('./03_GO+KEGG_enrichment')){
  dir.create('./03_GO+KEGG_enrichment')
}
setwd('./03_GO+KEGG_enrichment/')

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(lance)
library(enrichplot)
library(GOplot)

if(!dir.exists('./GO')){
  dir.create('./GO')
}
setwd('./GO')

GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库

inter_gene <- read.csv('../../02_hub_gene/hub_gene.csv', row.names = 1)

## gene ID转换
gene <- bitr(inter_gene$symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

## GO富集分析
GO <- enrichGO(gene$ENTREZID, # GO富集分析
               OrgDb = GO_database,
               keyType = "ENTREZID", # 设定读取的gene ID类型
               ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
               pvalueCutoff = 0.9, # 设定p值阈值
               qvalueCutoff = 0.9, # 设定q值阈值
               readable = T)
go_res <- data.frame(GO) # 将GO结果转为数据框
go_res <- mutate(go_res, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
go_res <- go_res[go_res$pvalue<0.05, ]
table(go_res$ONTOLOGY) 

write.csv(go_res, file = "./GO_res.csv")

# 2. KEGG富集分析---------------
rm(list = ls())
setwd('/data/nas1/xuguangxin/project/39_BJTC-630-4/')
if(!dir.exists('./03_GO+KEGG_enrichment')){
  dir.create('./03_GO+KEGG_enrichment')
}
setwd('./03_GO+KEGG_enrichment/')

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(lance)
library(enrichplot)
library(GOplot)

if(!dir.exists('./KEGG')){
  dir.create('./KEGG')
}
setwd('./KEGG')

GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa' # KEGG是hsa数据库

inter_gene <- read.csv('../../02_hub_gene/hub_gene.csv', row.names = 1)

## gene ID转换
gene <- bitr(inter_gene$symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

KEGG <- enrichKEGG(gene$ENTREZID,
                   keyType = "kegg",
                   organism = KEGG_database,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.5,
                   qvalueCutoff = 0.5)

kegg_res <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.csv(kegg_res@result, file = "./KEGG_res.csv")


# 3. GO富集分析可视化------------------
rm(list = ls())
setwd('/data/nas1/xuguangxin/project/39_BJTC-630-4/')
if(!dir.exists('./03_GO+KEGG_enrichment')){
  dir.create('./03_GO+KEGG_enrichment')
}
setwd('./03_GO+KEGG_enrichment/')

library(circlize)
library(grid)
library(graphics)
library(ComplexHeatmap)
library(GOplot)
library(tidyverse)

if(!dir.exists('./GO')){
  dir.create('./GO')
}
setwd('./GO')

go_res <- read.csv('./GO_res.csv', row.names = 1)
DEGs <- read.csv('../../01_scRNA/07_subtype/FADS3_cluster_DE.csv', row.names = 1, check.names = F)

Go_bp <- go_res[go_res$ONTOLOGY == 'BP', ]
Go_bp <- dplyr::arrange(Go_bp, desc(Count))
Go_bp <- Go_bp[1:5,]

# 提取CC的行并排序
Go_cc <- go_res[go_res$ONTOLOGY == 'CC', ]
Go_cc <- dplyr::arrange(Go_cc, desc(Count))
Go_cc <- Go_cc[1:5,]

# 提取MF的行并排序
Go_mf <- go_res[go_res$ONTOLOGY == 'MF', ]
Go_mf <- dplyr::arrange(Go_mf, desc(Count))
Go_mf <- Go_mf[1:5,]

go_res2 <- rbind(Go_bp, Go_cc, Go_mf)
write.csv(go_res2, 'go_res2.csv')
genelist <- data.frame(ID = rownames(DEGs), logFC = DEGs$avg_log2FC)

go <-  data.frame(Category = go_res2[,'ONTOLOGY'],
                  ID = go_res2[,'ID'],
                  Term = go_res2[,'Description'], 
                  Genes = gsub("/", ", ", go_res2[,'geneID']), 
                  adj_pval = go_res2[,'p.adjust'])

circ <- circle_dat(go, genelist)
go_term <- go$ID

i <- 1
finally_data <- data.frame()
for (i in c(1 : length(go_term))) {
  temp1 <- go_term[i]
  temp2 <- circ[circ$ID %in% temp1, ]
  temp3 <- data.frame(ID = unique(temp2$ID),
                      up_gene = nrow(temp2[temp2$logFC > 0, ]),
                      down_gene = nrow(temp2[temp2$logFC < 0, ]))
  finally_data <- rbind(finally_data, temp3)
}

finally_data <- merge(finally_data, go_res2, by = 'ID')
go_enrich_df <- data.frame(
  category = finally_data$ONTOLOGY,
  gene_num.min = 0,
  gene_num.max = 50, 
  gene_num.rich = finally_data$Count, 
  X.log10.p = -log10(finally_data$pvalue), 
  up.regulated = finally_data$up_gene, 
  down.regulated = finally_data$down_gene, 
  rich.factor = finally_data$richFactor, 
  Description = finally_data$ID
)
rownames(go_enrich_df) <- go_enrich_df$Description
go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
go_enrich_df <- go_enrich_df[order(go_enrich_df$category), ]

library(circlize)
{
  png('01.Go_circlize.png',w=12,h=12,units='in',res=600)
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  circle_size = unit(1, 'snpc')
  ##整体布局
  circos.par(gap.degree = 2, start.degree = 90)
  ##第一圈，绘制 ko
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
  ko_color <- c(rep('#EFC000B2', 5), rep('#4DBBD5B2', 5), rep('#00A087B2', 5))  #定义分组颜色
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  #圈图的高度、颜色等设置
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')  #ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称
      circos.axis(h = 'top', labels.cex = 0.8, major.tick.percentage = 0.4, labels.niceFacing = FALSE)  #绘制外周的刻度线
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将 ko id 文字标签添加在图中指定位置处
    } )
  colnames(go_enrich_df)[5] <- c("-log10(pvalue)")
  ##第二圈，绘制富集的基因和富集 p 值
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.rich', '-log10(pvalue)')]  #选择作图数据集，包括富集基因数量以及 p 值等信息
  label_data <- go_enrich_df['gene_num.rich']  #标签数据集，仅便于作图时添加相应的文字标识用
  p_max <- round(max(go_enrich_df$`-log10(pvalue)`)) + 1  #定义一个 p 值的极值，以方便后续作图
  colorsChoice <- colorRampPalette(c('white','#4393C3'))  #这两句用于定义 p 值的渐变颜色
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
      ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),1]
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
    } )
  
  ##第三圈，绘制上下调基因
  #首先基于表格中上下调基因的数量，计算它们的占比
  go_enrich_df$all.regulated <- go_enrich_df$up.regulated + go_enrich_df$down.regulated
  go_enrich_df$up.proportion <- go_enrich_df$up.regulated / go_enrich_df$all.regulated
  go_enrich_df$down.proportion <- go_enrich_df$down.regulated / go_enrich_df$all.regulated
  
  #随后，根据上下调基因的相对比例，分别计算它们在作图时的“区块坐标”和“长度”
  go_enrich_df$up <- go_enrich_df$up.proportion * go_enrich_df$gene_num.max
  plot_data_up <- go_enrich_df[c('Description', 'gene_num.min', 'up')]
  names(plot_data_up) <- c('Description', 'start', 'end')
  plot_data_up$type <- 1  #分配 1 指代上调基因
  
  go_enrich_df$down <- go_enrich_df$down.proportion * go_enrich_df$gene_num.max + go_enrich_df$up
  plot_data_down <- go_enrich_df[c('Description', 'up', 'down')]
  names(plot_data_down) <- c('Description', 'start', 'end')
  plot_data_down$type <- 2  #分配 2 指代下调基因
  
  #选择作图数据集（作图用）、标签数据集（添加相应的文字标识用），并分别为上下调基因赋值不同颜色
  plot_data <- rbind(plot_data_up, plot_data_down)
  label_data <- go_enrich_df[c('up', 'down', 'up.regulated', 'down.regulated')]
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('red', 'green'))
  #继续绘制圈图
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
      ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
      xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),4]
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
    } )
  go_enrich_df$`-log10(pvalue)` <- 50*go_enrich_df$`-log10(pvalue)`
  ##第四圈，绘制富集得分
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max', 'rich.factor')]  #选择作图数据集，标准化后的富集得分
  label_data <- go_enrich_df['category']  #将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色
  color_assign <- c('BP' = '#EFC000B2', 'CC' = '#4DBBD5B2', 'MF' = '#00A087B2')
  #E64B35B2" "#4DBBD5B2" "#00A087B2
  circos.genomicTrack(
    plot_data, ylim = c(0, 0.9), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称，并添加在下一句中匹配 ko 对应的高级分类，以分配颜色
      circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集得分，颜色代表 ko 的分类
      circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3)  #可选在富集得分等于 0.5 的位置处添加一个灰线
    } )
  ##绘图完毕后，不要忘了清除痕迹，以免影响下一次作图
  circos.clear()
  
  category_legend <- Legend(#E64B35B2" "#4DBBD5B2" "#00A087B2
    labels = c("BP:biological process", "CC:cellular component","MF:molecular function" ),
    type = 'points', pch = NA, background = c('#EFC000B2', '#4DBBD5B2', '#00A087B2'), 
    labels_gp = gpar(fontsize = 16), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  
  updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background = c('red', 'green'), 
    
    labels_gp = gpar(fontsize = 16), grid_height = unit(0.8, 'cm'), grid_width = unit(0.8, 'cm'))
  
  pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                         colorRampPalette(c('white','#4393C3'))(6)),
    legend_height = unit(4, 'cm'), labels_gp = gpar(fontsize = 12), 
    title_gp = gpar(fontsize = 16), title_position = 'topleft', title = '-Log10(pvalue)')
  
  # lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
  # grid.draw(lgd_list_vertical)
  
  lgd_list_vertical <- packLegend(updown_legend, category_legend)
  pushViewport(viewport(x = 0.55, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  
  lgd_list_vertical <- packLegend(pvalue_legend)
  pushViewport(viewport(x = 0.45, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  
  circos.clear()
  
  # lgd_list_vertical <- packLegend(updown_legend, pvalue_legend)
  # pushViewport(viewport(x = 0.5, y = 0.5))
  # grid.draw(lgd_list_vertical)
  # upViewport()
  # 
  # lgd_list_vertical <- packLegend(category_legend)
  # pushViewport(viewport(x = 0.5, y = 0.5))
  # grid.draw(lgd_list_vertical)
  # upViewport()
  
  dev.off()
  dev.off()
}
{
  pdf('01.Go_circlize.pdf',w=12,h=12)
  par(mar = c(2,2,2,2),cex=1.5,family="Times")
  circle_size = unit(1, 'snpc')
  ##整体布局
  circos.par(gap.degree = 2, start.degree = 90)
  ##第一圈，绘制 ko
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
  ko_color <- c(rep('#EFC000B2', 5), rep('#4DBBD5B2', 5), rep('#00A087B2', 5))  #定义分组颜色
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  #圈图的高度、颜色等设置
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')  #ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称
      circos.axis(h = 'top', labels.cex = 0.8, major.tick.percentage = 0.4, labels.niceFacing = FALSE)  #绘制外周的刻度线
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将 ko id 文字标签添加在图中指定位置处
    } )
  colnames(go_enrich_df)[5] <- c("-log10(pvalue)")
  ##第二圈，绘制富集的基因和富集 p 值
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.rich', '-log10(pvalue)')]  #选择作图数据集，包括富集基因数量以及 p 值等信息
  label_data <- go_enrich_df['gene_num.rich']  #标签数据集，仅便于作图时添加相应的文字标识用
  p_max <- round(max(go_enrich_df$`-log10(pvalue)`)) + 1  #定义一个 p 值的极值，以方便后续作图
  colorsChoice <- colorRampPalette(c('white','#4393C3'))  #这两句用于定义 p 值的渐变颜色
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
      ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),1]
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
    } )
  
  ##第三圈，绘制上下调基因
  #首先基于表格中上下调基因的数量，计算它们的占比
  go_enrich_df$all.regulated <- go_enrich_df$up.regulated + go_enrich_df$down.regulated
  go_enrich_df$up.proportion <- go_enrich_df$up.regulated / go_enrich_df$all.regulated
  go_enrich_df$down.proportion <- go_enrich_df$down.regulated / go_enrich_df$all.regulated
  
  #随后，根据上下调基因的相对比例，分别计算它们在作图时的“区块坐标”和“长度”
  go_enrich_df$up <- go_enrich_df$up.proportion * go_enrich_df$gene_num.max
  plot_data_up <- go_enrich_df[c('Description', 'gene_num.min', 'up')]
  names(plot_data_up) <- c('Description', 'start', 'end')
  plot_data_up$type <- 1  #分配 1 指代上调基因
  
  go_enrich_df$down <- go_enrich_df$down.proportion * go_enrich_df$gene_num.max + go_enrich_df$up
  plot_data_down <- go_enrich_df[c('Description', 'up', 'down')]
  names(plot_data_down) <- c('Description', 'start', 'end')
  plot_data_down$type <- 2  #分配 2 指代下调基因
  
  #选择作图数据集（作图用）、标签数据集（添加相应的文字标识用），并分别为上下调基因赋值不同颜色
  plot_data <- rbind(plot_data_up, plot_data_down)
  label_data <- go_enrich_df[c('up', 'down', 'up.regulated', 'down.regulated')]
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('red', 'green'))
  #继续绘制圈图
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
      ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
      xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),4]
      circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
    } )
  go_enrich_df$`-log10(pvalue)` <- 50*go_enrich_df$`-log10(pvalue)`
  ##第四圈，绘制富集得分
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max', 'rich.factor')]  #选择作图数据集，标准化后的富集得分
  label_data <- go_enrich_df['category']  #将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色
  color_assign <- c('BP' = '#EFC000B2', 'CC' = '#4DBBD5B2', 'MF' = '#00A087B2')
  #E64B35B2" "#4DBBD5B2" "#00A087B2
  circos.genomicTrack(
    plot_data, ylim = c(0, 0.9), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称，并添加在下一句中匹配 ko 对应的高级分类，以分配颜色
      circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集得分，颜色代表 ko 的分类
      circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3)  #可选在富集得分等于 0.5 的位置处添加一个灰线
    } )
  ##绘图完毕后，不要忘了清除痕迹，以免影响下一次作图
  circos.clear()
  
  category_legend <- Legend(#E64B35B2" "#4DBBD5B2" "#00A087B2
    labels = c("BP:biological process", "CC:cellular component","MF:molecular function" ),
    type = 'points', pch = NA, background = c('#EFC000B2', '#4DBBD5B2', '#00A087B2'), 
    labels_gp = gpar(fontsize = 16), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  
  updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background = c('red', 'green'), 
    
    labels_gp = gpar(fontsize = 16), grid_height = unit(0.8, 'cm'), grid_width = unit(0.8, 'cm'))
  
  pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                         colorRampPalette(c('white','#4393C3'))(6)),
    legend_height = unit(4, 'cm'), labels_gp = gpar(fontsize = 12), 
    title_gp = gpar(fontsize = 16), title_position = 'topleft', title = '-Log10(pvalue)')
  
  # lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
  # grid.draw(lgd_list_vertical)
  
  lgd_list_vertical <- packLegend(updown_legend, category_legend)
  pushViewport(viewport(x = 0.55, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  
  lgd_list_vertical <- packLegend(pvalue_legend)
  pushViewport(viewport(x = 0.45, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  
  circos.clear()
  
  # lgd_list_vertical <- packLegend(updown_legend, pvalue_legend)
  # pushViewport(viewport(x = 0.5, y = 0.5))
  # grid.draw(lgd_list_vertical)
  # upViewport()
  # 
  # lgd_list_vertical <- packLegend(category_legend)
  # pushViewport(viewport(x = 0.5, y = 0.5))
  # grid.draw(lgd_list_vertical)
  # upViewport()
  
  dev.off()
  dev.off()
}

# 4. KEGG富集分析可视化----------------
rm(list = ls())
setwd('/data/nas1/xuguangxin/project/39_BJTC-630-4/')
if(!dir.exists('./03_GO+KEGG_enrichment')){
  dir.create('./03_GO+KEGG_enrichment')
}
setwd('./03_GO+KEGG_enrichment/')

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(lance)
library(enrichplot)
library(GOplot)

if(!dir.exists('./KEGG')){
  dir.create('./KEGG')
}
setwd('./KEGG')

library(tidyverse)

KEGG_res <- read.csv('./KEGG_res.csv', row.names = 1, check.names = F)
DEGs <- read.csv('../../01_scRNA/07_subtype/FADS3_cluster_DE.csv')
DEGs <- dplyr::select(DEGs, X, avg_log2FC)
colnames(DEGs) <- c('symbol', 'log2FoldChange')

plot_data <- KEGG_res[KEGG_res$pvalue < 0.05, ]
plot_data <- plot_data[c(1:6), ]
plot_data$Description <- gsub(' - reactive oxygen species', '', plot_data$Description)
# plot_data$Description <- gsub('in Helicobacter pylori infection', '', plot_data$Description)

temp1 <- plot_data
temp1 <- dplyr::select(temp1, ID, geneID)
temp1 <- separate(temp1, geneID, into = as.character(c(1 : 20)), sep = '/')
ncol(temp1)

i <- 1
x <- 2
finally_data <- data.frame()
for (i in c(1 : nrow(temp1))) {
  for (x in c(2 : 20)) {
    temp2 <- temp1[i, ]
    temp3 <- data.frame(ID = temp2[1, 1],
                        symbol = temp2[1, x])
    finally_data <- rbind(finally_data, temp3)
  }
  
}
finally_data <- na.omit(finally_data)

kegg_df <- merge(finally_data, plot_data, by = 'ID')
kegg_df2 <- merge(DEGs, kegg_df, by = 'symbol')

colnames(kegg_df2)
kegg_df2$abs_logFC <- abs(kegg_df2$log2FoldChange)
#treemap(go.df2,index=c("Description","gene"),vSize = "abs_logFC")
se_index <- c("Description","symbol")
#root="GO Enrichment Terms"
library(tidyverse)
library(tidygraph)
library(ggraph)
library(ccgraph)
library(treemap)
# devtools::install_github( "gaospecial/ccgraph")
nodes_GO <- gather_graph_node(kegg_df2, index = se_index, value = "abs_logFC",root="KEGG")
edges_GO <- gather_graph_edge(kegg_df2, index = se_index, root="KEGG")

head(nodes_GO, 10)
head(edges_GO, 10)

graph <- tbl_graph(nodes_GO,edges_GO)
library(ggsci)
gc <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 2/3,edge_width=0.5) + 
  geom_node_point(aes(size = node.size, 
                      alpha=0.5,
                      color = node.branch,
                      filter=node.level!="all")
  ) + 
  scale_size(range = c(0.5,5)) + 
  theme(legend.position = "none") + 
  
  # scale_edge_color_brewer(c('#e41a1c', '#377eb8', '#984ea3', '#ff7f00', '#a65628', '#00A087B2', '#EFC000B2')) + #
  scale_color_manual(values = c('#e41a1c', '#377eb8', '#984ea3', '#ff7f00', '#a65628', '#00A087B2', '#EFC000B2'))+
  # scale_color_brewer(c('#e41a1c', '#377eb8', '#984ea3', '#ff7f00', '#a65628', '#00A087B2', '#EFC000B2')) +
  
  geom_node_text(
    aes(
      x = 1.0175 * x, 
      y = 1.0175 * y, 
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf,
      color = node.branch
    ),
    size = 6, hjust = 'outward') +
  
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all"),
        color = node.branch),
    fontface="bold",
    #alpha=0.5,
    size=5,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) +
  coord_cartesian(xlim=c(-1.25,1.25),ylim = c(-1.25,1.25)) 
gc
ggsave(filename = "01.KEGG_ggraph.pdf", height = 12, width = 12, gc)
ggsave(filename = "01.KEGG_ggraph.png", height = 12, width = 12, gc)
