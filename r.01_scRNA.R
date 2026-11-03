# 1. 数据处理-----------------
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/")
if (! dir.exists("./01_scRNA")){
  dir.create("./01_scRNA")
}
setwd("./01_scRNA")

library(Seurat)
library(tidyverse)

folders <- list.files('../00_rawdata/snRNA_data/all_data/')

i <- 1
scList <- list()
folders[i]
for (i in c(1 : length(folders))) {
  dir <- paste0('../00_rawdata/snRNA_data/all_data/', folders[i])
  temp1 <- Read10X(data.dir = dir)
  temp1 <- CreateSeuratObject(counts = temp1,
                              min.features = 200,
                              min.cells = 3, 
                              project = folders[i])
  temp1@meta.data$orig.ident <- folders[i]
  scList <- append(scList, temp1)
}


scRNA<- merge(scList[[1]],
              y = c(scList[-1]), 
              project = "GSE161921_GSE174691")

table(scRNA$orig.ident)
head(scRNA@meta.data)

group <- data.frame(sample = c('control', 'disease'), group = c('control','disease'))
##添加分组
smps <- scRNA$orig.ident %>% str_sub(1, 7)
table(smps)
# grps <- group
group.label <- group$group[match(smps, group$sample)]
table(group.label)
names(group.label) <- Cells(scRNA)
group.label
scRNA = AddMetaData(scRNA, group.label, col.name = 'group')
table(scRNA$group)
head(scRNA@meta.data, n = 10)

###线粒体基因表达占比(正则表达式，表示以MT-开头；scRNA[["percent.mt"]]这种写法会在meta.data矩阵加上一列)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")  

###计算核糖体基因比例
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")
head(scRNA@meta.data)

system.time(save(scRNA, file = "01.scRNA_orig.Rdata"))  

# 2. 质控------------------------------------
rm(list = ls())
setwd(('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/'))
if (! dir.exists("./01_QC")){
  dir.create("./01_QC")
}
setwd("./01_QC")
load('../01.scRNA_orig.Rdata')
length(unique(scRNA$orig.ident))

library(patchwork)
library(ggplot2)
library(Seurat)

Screening_res <- data.frame()
## 质控前
temp1 <- data.frame(cell_number = length(colnames(scRNA)),
                    gene_number = length(rownames(scRNA)))  
Screening_res <- rbind(Screening_res, temp1)
rownames(Screening_res) <- 'before_quality'

# 设置绘图元素
plot.featrures <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
theme.set = theme(
  axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 18,face='bold',family = "Times"),
  text = element_text(family = "Times"))
# 质控前小提琴图
plots <- list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA,
                       #type.by=Group,
                       pt.size = 0,
                       features = plot.featrures[i]) + theme.set+ NoLegend()}
violin <- wrap_plots(plots = plots, nrow = 1, ncol = 3, cols = palettes(category = "random", 22, show_col=F))
violin

ggsave("01.vlnplot_before_qc.pdf", plot = violin, width = 12, height = 8)
ggsave("01.vlnplot_before_qc.png", plot = violin, width = 12, height = 8,dpi = 600)

## 质控后
## 设置质控标准
minGene = quantile(scRNA$nFeature_RNA, .05) #507.55
maxGene = quantile(scRNA$nFeature_RNA, .95) #4683.45
maxUMI = quantile(scRNA$nCount_RNA, .95) #14973.25

# minGene <- 350
# maxGene <- 3000
# maxUMI <- 15000

Screening_condition <- data.frame('minGene_.05' = minGene, 
                                  'maxGene_.95' = maxGene,
                                  'maxUMI_.95' = maxUMI)
rownames(Screening_condition) <- 'number'
write.csv(Screening_condition, 'Screening_condition.csv')
pctMT=25

scRNA <- subset(scRNA, subset =
                  nCount_RNA < maxUMI &
                  nFeature_RNA > minGene &
                  nFeature_RNA < maxGene &
                  percent.mt < pctMT)
temp1 <- data.frame(cell_number = length(colnames(scRNA)),
                    gene_number = length(rownames(scRNA)))  
Screening_res <- rbind(Screening_res, temp1)
rownames(Screening_res)[2] <- 'after_quality'
write.csv(Screening_res, 'Screening_res.csv')

# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
theme.set = theme(
  axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 18,face='bold',family = "Times"),
  text = element_text(family = "Times"))
# 质控后小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA,
                       #GSE.by=GSE,
                       pt.size = 0,
                       features = plot.featrures[i]) + theme.set + NoLegend()}
violin <- wrap_plots(plots = plots, nrow = 1, ncol = 3, cols = palettes(category = "random", 22, show_col=F))
violin
ggsave("02.vlnplot_after_qc.pdf", plot = violin, width = 12, height = 8)
ggsave("02.vlnplot_after_qc.png", plot = violin, width = 12, height = 8,dpi = 600)

system.time(save(scRNA, file = "scRNA_qc.Rdata"))

# 3. Seurat v3的标准整合流程-----------------------------
## NormalizeData函数的作用：消除不同细胞测序深度的影响
rm(list = ls())
setwd(('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/'))
if (! dir.exists("./02_Integrate")){
  dir.create("./02_Integrate")
}
setwd("./02_Integrate")
load('../01_QC/scRNA_qc.Rdata')

library(Seurat)
library(tidyverse)

combined <- NormalizeData(scRNA)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

pdf("01.before_intergrated_PCA.pdf",w=7,h=5)
DimPlot(combined, reduction = "pca",group.by = 'group')
dev.off()
png("01.before_intergrated_PCA.png",w=7,h=5,units = 'in',res = 600)
DimPlot(combined, reduction = "pca",group.by = 'group')
dev.off()

system.time(save(combined , file = "unintergrated.Rdata"))

##整合方法1
combined <- IntegrateLayers(
  object = combined, 
  method = CCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  # k.anchor = 20, # 整合强度
  verbose = FALSE,
  k.weight = 100)

scRNA <- combined
scRNA <- JoinLayers(scRNA)
system.time(save(scRNA , file = "scRNA_intergrated.Rdata"))

# 4. PCA UMAP降维聚类-------------------
rm(list = ls())
setwd(('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/'))
if (! dir.exists("./03_PCA")){
  dir.create("./03_PCA")
}
setwd("./03_PCA")
load('../02_Integrate/scRNA_intergrated.Rdata')

scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA), 10)
# plot variable features with and without labels
p <- VariableFeaturePlot(scRNA)
p <- LabelPoints(plot = p, 
                 points = top10, 
                 repel = T)
p
ggsave(filename = "01.feature_selection.pdf", p, w=6, h=6)
ggsave(filename = "01.feature_selection.png", p, w=6, h=6)

## PCA降维分析
## scale标准化——很多基因的reads数多达上千，但是有很多基因reads很少（ScaleData作用：是将全部基因reads数据压缩为0-1之间的数值）
scRNA.nor.sca <- ScaleData(scRNA)
scRNA.norm.pca <- RunPCA(scRNA.nor.sca, features = VariableFeatures(object = scRNA.nor.sca), npcs = 50)
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")
## PCA主成分分析
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = 'group')
pdf("02.PCA.pdf",w=7,h=5)
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = 'group')
dev.off()
png("02.PCA.png",w=7,h=5,units = 'in',res = 600)
DimPlot(scRNA.norm.pca, reduction = "pca",group.by = 'group')
dev.off()
## PCA热图
pdf("03.PCA_heatmap.pdf",w=9,h=9)
DimHeatmap(scRNA.norm.pca, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()
png("03.PCA_heatmap.png",w=9,h=9,units = 'in',res = 600)
DimHeatmap(scRNA.norm.pca, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

## 首先找到最佳聚类数

## JackStraw函数的功能/作用：
## 随机排列数据的子集，并计算这些“随机”基因的预测PCA分数。
## 然后将“随机”基因的PCA分数与观察到的PCA分数进行比较，以确定统计显著性。
## 最终结果是每个基因与每个主成分相关的p值
## 简而言之: JackStraw()函数可以计算出每个主成分中各基因的P值，用于判断哪些主成分更具有统计学意义
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50) # num.replicate : 要执行的复制采样数；dims : 计算重要性的PC数
## ScoreJackStraw()函数：用于量化主成分的显著性强度，富含低P值基因较多的主成分更有统计学意义
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
system.time(save(scRNA.norm.pca, file = "scRNA.norm.pca.Jack.Rdata"))
load('./scRNA.norm.pca.Jack.Rdata')
## JackStrawPlot()函数：可视化比较每个主成分的 p 值分布和均匀分布（虚线）
plot_pca <- JackStrawPlot(scRNA.norm.pca, dims = 1:50)
plot_pca
ggsave(filename = '04.pca_cluster.pdf',plot_pca,w=13,h=8)
ggsave(filename = '04.pca_cluster.png',plot_pca,w=13,h=8,dpi = 600)
plot_elbow <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot_elbow
ggsave(filename = '05.pca_sd.pdf',plot_elbow,w=6,h=5)
ggsave(filename = '05.pca_sd.png',plot_elbow,w=6,h=5,dpi = 600)

# load('./scRNA.norm.pca.Jack.Rdata')
## FindNeighbors()函数作用：计算邻接距离
scRNA.norm.pca.clu <- FindNeighbors(scRNA.norm.pca, dims = 1:30)
## FindClusters()函数作用：以迭代方式将细胞分组在一起,resolution参数表示区分细胞群的分辨率，resolution越大，分的群也就越多。
scRNA.norm.pca.clu <- FindClusters(scRNA.norm.pca.clu, resolution = 0.2) # 0.4

## TSNE
## TSNE聚类分析
TSNE <- RunTSNE(object = scRNA.norm.pca.clu, dims = 1:30)                      
system.time(save(TSNE, file = "TSNE.Rdata"))

## 样本可视化
theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 18,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 18,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_text(size = 20,face='bold',family = "Times"),
  text = element_text(family = "Times"))

head(TSNE@meta.data)
pdf(file="06.TSNE_all.pdf",width=10,height=8,family='Times')
DimPlot(TSNE ,reduction = "tsne",label.size = 6,label = TRUE,pt.size=0.8)+theme.set
dev.off()

png(file="06.TSNE_all.png",width=10,height=8,family='Times',units='in',res=600)
DimPlot(TSNE ,reduction = "tsne",label.size = 6,label = TRUE, pt.size=0.8)+theme.set
dev.off()

pdf(file="07.TSNE_Group.pdf",width=12,height=6,family='Times')
DimPlot(TSNE ,reduction = "tsne",label.size = 6,label = T,pt.size=0.8, split.by= "group")+theme.set
dev.off()

png(file="07.TSNE_Group.png",width=12,height=6,family='Times',units='in',res=600)
DimPlot(TSNE ,reduction = "tsne",label.size = 6,label =T, pt.size=0.8,split.by  = "group")+theme.set
dev.off()

## UMAP
load('./scRNA.norm.pca.Jack.Rdata')
## FindNeighbors()函数作用：计算邻接距离
scRNA.norm.pca.clu <- FindNeighbors(scRNA.norm.pca, dims = 1:30)
## FindClusters()函数作用：以迭代方式将细胞分组在一起,resolution参数表示区分细胞群的分辨率，resolution越大，分的群也就越多。
scRNA.norm.pca.clu <- FindClusters(scRNA.norm.pca.clu, resolution = 0.4) # 0.4

## UMAP
## UMAP聚类分析
UMAP <- RunUMAP(object = scRNA.norm.pca.clu, dims = 1:30)                      
system.time(save(UMAP, file = "UMAP.Rdata"))

## 样本可视化
theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 18,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 18,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_text(size = 20,face='bold',family = "Times"),
  text = element_text(family = "Times"))

head(UMAP@meta.data)
pdf(file="06.UMAP_all.pdf",width=10,height=8,family='Times')
DimPlot(UMAP ,reduction = "umap",label.size = 6,label = TRUE,pt.size=0.8)+theme.set
dev.off()

png(file="06.UMAP_all.png",width=10,height=8,family='Times',units='in',res=600)
DimPlot(UMAP ,reduction = "umap",label.size = 6,label = TRUE, pt.size=0.8)+theme.set
dev.off()

pdf(file="07.UMAP_Group.pdf",width=12,height=6,family='Times')
DimPlot(UMAP ,reduction = "umap",label.size = 6,label = T,pt.size=0.8, split.by= "group")+theme.set
dev.off()

png(file="07.UMAP_Group.png",width=12,height=6,family='Times',units='in',res=600)
DimPlot(UMAP ,reduction = "umap",label.size = 6,label =T, pt.size=0.8,split.by  = "group")+theme.set
dev.off()

# 5. 细胞类型注释------
## 5.1 TSNE
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./04_SingleR")){
  dir.create("./04_SingleR")
}
setwd("./04_SingleR")

if (! dir.exists("./TSNE")){
  dir.create("./TSNE")
}
setwd("./TSNE")

load('../../03_PCA/TSNE.Rdata')

gene <- c(
  'PDGFRB', 'ACTA2', 'MYH11',## Smooth muscle
  'PTPRC', 'CD14', 'FOS', ## Leukocyte
  'POSTN', 'COL1A1', 'COL6A3', ## Fibroblast
  'SCN7A', ## Neuronal
  'MYBPC3', 'MYH7B', ## Cardiomyocyte
  'CFD', 'AQP1', 'PECAM1', 'VWF' ## Endothelial
)

DotPlot(TSNE, features = gene, assay = 'RNA') + 
  RotatedAxis() + 
  scale_color_gradientn(colours = c('gray','#87CEFA','#CAB2D6','#6A3D9A'))

new.cluster.ids = c(
  "0" = "Endothelial",  
  "1" = "Cardiomyocyte",  
  "2" = "Neuronal",  
  "3" = "Fibroblast", 
  "4" = "Endothelial",   
  "5" = "Leukocyte"  ,   
  "6" = "Smooth muscle",  
  "7" = "Fibroblast"  , 
  "8" = "Cardiomyocyte", 
  '9' = 'Neuronal', 
  '10' = 'Smooth muscle',
  '11' = 'Endothelial'
)  
TSNE <- RenameIdents(TSNE, new.cluster.ids)
TSNE$celltype = Idents(TSNE)
table(TSNE$celltype)

theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times",face='bold'))

pdf(file="03.celltype_TSNE.pdf",width=9,height=7,family='Times')
DimPlot(TSNE ,reduction = "tsne",label.size = 6,label.color = 'black',label = T)+theme.set
dev.off()

png(file="03.celltype_TSNE.png",width=9,height=7,family='Times',units='in',res=600)
DimPlot(TSNE ,reduction = "tsne",label.size = 6,label.color = 'black',label =T)+theme.set
dev.off()

pdf(file="04.celltype_TSNE_group.pdf",width=12,height=8,family='Times')
DimPlot(TSNE ,reduction = "tsne",label.size =4,label = T,pt.size=1, split.by= "group")+theme.set
dev.off()

png(file="04.celltype_TSNE_group.png",width=12,height=8,family='Times',units='in',res=600)
DimPlot(TSNE ,reduction = "tsne",label.size = 4,label =T, pt.size=1,split.by  = "group")+theme.set
dev.off()

pdf(file="05.celltype_dotplot.pdf",width=10,height=7,family='Times')
DotPlot(TSNE, features = gene, group.by = 'celltype', assay = 'RNA') + 
  RotatedAxis() + 
  theme.set+
  scale_color_gradientn(colours = c('gray','#87CEFA','#CAB2D6','#6A3D9A'))
dev.off()

png(file="05.celltype_dotplot.png",width=10,height=7,family='Times',units='in',res=600)
DotPlot(TSNE, features = gene, group.by = 'celltype', assay = 'RNA') + 
  RotatedAxis() + 
  theme.set+
  scale_color_gradientn(colours = c('gray','#87CEFA','#CAB2D6','#6A3D9A'))
dev.off()

## cluster图
cluster_data <- TSNE@reductions$tsne@cell.embeddings %>% 
  as.data.frame() %>% 
  cbind(cluster = TSNE@meta.data$seurat_clusters)
head(cluster_data)

cluster_med <- cluster_data %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )

library(ggrepel)
library(RColorBrewer)
library(ggunchull)
# allcolour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
col <- sample(RColorBrewer::brewer.pal(10, 'Paired'))
col <- c( "#FDBF6F", "#F08080", "#6A3D9A", "#1F78B4",
                   "#33A02C", "#B2DF8A", "#A6CEE3", "#FF7F00", 
                   "#FFA500", "#CAB2D6", '#808000', '#007BFF',
                   "#7B68EE", "#98FB98", "#D2B48C", "#FFE4B5")
                   
p1 <- ggplot(cluster_data, aes(x= tSNE_1 , y = tSNE_2 ,color = cluster)) +  
  geom_point(size = 1 , alpha =1 )  +  
  scale_color_manual(values = col)+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm')  # 设置legend标签之间的大小
  ) + 
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_segment(aes(x = min(cluster_data$tSNE_1) , y = min(cluster_data$tSNE_2) ,
                   xend = min(cluster_data$tSNE_1) + 14, yend = min(cluster_data$tSNE_2) ),
               colour = "black", 
               size = 1,
               arrow = arrow(length = unit(0.3,"cm"))
  )+
  geom_segment(aes(x = min(cluster_data$tSNE_1)  , y = min(cluster_data$tSNE_2)  ,
                   xend = min(cluster_data$tSNE_1) , yend = min(cluster_data$tSNE_2) + 14),
               colour = "black", 
               size = 1,
               arrow = arrow(length = unit(0.3,"cm"))
  ) +
  annotate("text", 
           x = min(cluster_data$tSNE_1) + 6, 
           y = min(cluster_data$tSNE_2) -4, 
           label = "tSNE_1",
           color = "black",
           size = 6, 
           fontface = "bold") + 
  annotate("text", 
           x = min(cluster_data$tSNE_1) -4, 
           y = min(cluster_data$tSNE_2) + 6, 
           label = "tSNE_2",
           color = "black",
           size = 6, 
           fontface = "bold",
           angle = 90) +
  geom_label_repel(data = cluster_med, 
                   aes(label = cluster), 
                   max.overlaps = 20, 
                   size = 9,
                   color = "black",
                   box.padding = unit(0.5, "lines"), 
                   min.segment.length = 0,
                   point.padding = unit(0.8, "lines"),
                   show.legend = FALSE
  )
p1
## 细胞图
celltype_data <- TSNE@reductions$tsne@cell.embeddings %>% 
  as.data.frame() %>% 
  cbind(celltype = TSNE@meta.data$celltype)
head(celltype_data)

celltype_med <- celltype_data %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )

library(ggrepel)
library(RColorBrewer)
library(ggunchull)
# allcolour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
# col <- sample(RColorBrewer::brewer.pal(10, 'Paired'))
col <- c( "#FDBF6F", "#F08080", "#6A3D9A", "#1F78B4",
                   "#33A02C", "#B2DF8A", "#A6CEE3", "#FF7F00", 
                   "#FFA500", "#CAB2D6", '#808000', '#007BFF',
                   "#7B68EE", "#98FB98", "#D2B48C", "#FFE4B5")
                   
p2 <- ggplot(celltype_data, aes(x= tSNE_1 , y = tSNE_2 ,color = celltype)) +  
  geom_point(size = 1 , alpha =1 )  +  
  scale_color_manual(values = col)+
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm')  # 设置legend标签之间的大小
  ) + 
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_segment(aes(x = min(celltype_data$tSNE_1) , y = min(celltype_data$tSNE_2) ,
                   xend = min(celltype_data$tSNE_1) + 14, yend = min(celltype_data$tSNE_2) ),
               colour = "black", 
               size = 1,
               arrow = arrow(length = unit(0.3,"cm"))
  )+
  geom_segment(aes(x = min(celltype_data$tSNE_1)  , y = min(celltype_data$tSNE_2)  ,
                   xend = min(celltype_data$tSNE_1) , yend = min(celltype_data$tSNE_2) + 14),
               colour = "black", 
               size = 1,
               arrow = arrow(length = unit(0.3,"cm"))
  ) +
  annotate("text", 
           x = min(celltype_data$tSNE_1) + 6, 
           y = min(celltype_data$tSNE_2) -4, 
           label = "tSNE_1",
           color = "black",
           size = 6, 
           fontface = "bold") + 
  annotate("text", 
           x = min(celltype_data$tSNE_1) -4, 
           y = min(celltype_data$tSNE_2) + 6, 
           label = "tSNE_2",
           color = "black",
           size = 6, 
           fontface = "bold",
           angle = 90) +
  geom_label_repel(data = celltype_med, 
                   aes(label = celltype), 
                   max.overlaps = 20, 
                   size = 9, 
                   color = "black",
                   box.padding = unit(0.5, "lines"), 
                   min.segment.length = 0,
                   point.padding = unit(0.8, "lines"),
                   show.legend = FALSE 
  )
p2

library(cowplot)
combined_plot <- plot_grid(p1, p2, 
                           ncol = 2, 
                           # rel_heights = c(1, 4), 
                           align = 'v')
combined_plot

ggsave("03.combined_plot.pdf", plot = combined_plot, width = 18, height = 8)
ggsave("03.combined_plot.png", plot = combined_plot, width = 18, height = 8)

system.time(save(TSNE, file = "TSNE_finally.Rdata"))

## 不同细胞比例
## 疾病
table(TSNE$group)
scRNA_case <- subset(TSNE, group == 'disease')
phe <- scRNA_case@meta.data

table(phe$celltype)
cell_type_case <- as.data.frame(sort(table(phe$celltype)))
colnames(cell_type_case) <- c('Celltype', 'disease')

## 对照
scRNA_Normal <- subset(TSNE,group == 'control')
phe2 <- scRNA_Normal@meta.data
table(phe2$celltype)
cell_type_Normal <- as.data.frame(sort(table(phe2$celltype)))
colnames(cell_type_Normal) <- c('Celltype', 'control')

cell_type <- merge(cell_type_Normal, cell_type_case, by='Celltype')
write.csv(cell_type, 'cell_type.csv')

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
library(ggplot2)
mydata <- melt(cell_type, id.vars = "Celltype", variable.name = "Group", value.name = "Quantity")

ggplot(mydata,aes(x=Group,y=Quantity,fill=Celltype))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(10, "Set3"))+
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="Celltype",y="Proportion",fill="")

ggsave('06.cell_proportion.pdf',w=10,h=8)
ggsave('06.cell_proportion.png',w=10,h=8)

table(TSNE$orig.ident)
cell_type_sample_finally <- data.frame(Celltype = cell_type$Celltype,
                                       var1 = 0)
i <- 1 
unique(TSNE$orig.ident)[i]
for (i in c(1 : length(unique(TSNE$orig.ident)))) {
  temp1 <- subset(TSNE, orig.ident == unique(TSNE$orig.ident)[i])
  temp2 <- temp1@meta.data
  table(temp2$celltype)
  temp3 <- as.data.frame(sort(table(temp2$celltype)))
  colnames(temp3) <- c('Celltype', unique(TSNE$orig.ident)[i])
  cell_type_sample_finally <- merge(cell_type_sample_finally, temp3, by = 'Celltype')
}
cell_type_sample_finally <- dplyr::select(cell_type_sample_finally, -var1)

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
library(ggplot2)
mydata <- melt(cell_type_sample_finally, id.vars = "Celltype", variable.name = "sample", value.name = "Quantity")
table(mydata$sample)

mydata$sample <- str_sub(mydata$sample, start = 1, end = 9)
table(mydata$sample)

mydata$group <- ifelse(grepl('control', mydata$sample), 'control', 'disease')
mydata <- mydata[order(mydata$group), ]

ggplot(mydata,aes(x = sample, y= Quantity, fill=Celltype))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(10, "Set3"))+
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times', angle = 45, vjust = 1, hjust = 1),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="Celltype",y="Proportion",fill="")

ggsave('07.cell_proportion_sample.pdf',w=12,h=8)
ggsave('07.cell_proportion_sample.png',w=12,h=8)

# 6. 细胞通讯--------------------------------
## 6.1 对照组--------------------
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./05_cellchat")){
  dir.create("./05_cellchat")
}
setwd("./05_cellchat")

library(CellChat) 
library(Seurat)
library(ggplot2)

load('../04_SingleR/TSNE/TSNE_finally.Rdata')
scRNA <- TSNE
table(scRNA$celltype)
data.input <- scRNA@assays$RNA$counts # normalized data matrix
meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$group)
cell.use = rownames(meta)[meta$group == "control"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$celltype)
unique(meta$celltype)

## CellChat的输入需要的matrix和meta已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
## 如果CellChat中没有meta信息，需要手动添加
# cellchat <- addMeta(cellchat, meta = meta) ## 增加meta信息
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将"celltype"设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents <- droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents), unique(cellchat@idents)))
groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量
groupSize

## 基于配体-受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的，这里我们用人
showDatabaseCategory(CellChatDB) # 可视化数据库的基本构成

library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  # 展示互作记录
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use # 将提取出来的数据库内容加载到cellchat对象中

cellchat <- subsetData(cellchat, features = NULL) # 如果有感兴趣的基因可以填写到features这里
cellchat <- identifyOverExpressedGenes(cellchat) # 首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat) # 寻找高表达的通路
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human) # 将基因映射为蛋白

## 计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat,type = "truncatedMean", trim = 0.05, population.size = TRUE)
# cellchat <- computeCommunProb(cellchat) # 有个参数为raw.use = T，意味着上一步基因映射为蛋白这一步就没用了，计算胞间通讯概率时会用原始数据，如果为F，就会用映射的蛋白做计算
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
df.net <- subsetCommunication(cellchat) # 将细胞通讯预测结果以数据框的形式取出

if (! dir.exists("./control")){
  dir.create("./control")
}
setwd("./control")

write.csv(df.net, 'df.net.csv')

## 信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))

## 相互作用数目
pdf('04.number_of_interactions.pdf',w=6,h=5,family='serif')
netVisual_heatmap(cellchat)
dev.off()
png('04.number_of_interactions.png',w=600,h=500,family='serif')
netVisual_heatmap(cellchat)
dev.off()
#> Do heatmap based on a merged object
# netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object

pdf('01.net_circle_number.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

png('01.net_circle_number.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
# 
pdf('02.net_circle_weight.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

png('02.net_circle_weight.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)

#需要指定受体细胞和配体细胞
cellchat@data.signaling
pdf('03.single_circle.pdf',w=9,h=14)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
png('03.single_circle.png',w=9,h=14,units = 'in',res = 600)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()


library(NMF)
library(ggalluvial)
# selectK(cellchat, pattern = "incoming")

nPatterns<-2
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "incoming",
                                          k = nPatterns,
                                          width = 5, height = 9)
pdf('05.incoming_pattern.pdf',w=10,h=6)
identifyCommunicationPatterns(cellchat, 
                              pattern = "incoming", 
                              k = nPatterns,
                              width = 5, height = 9)
dev.off()
png('05.incoming_pattern.png',w=10,h=6,units = 'in',res = 600)
identifyCommunicationPatterns(cellchat, 
                              pattern = "incoming", 
                              k = nPatterns,
                              width = 5, height = 9)
dev.off()

pdf('06.incoming_pattern_river.pdf',w=10,h=6)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
png('06.incoming_pattern_river.png',w=10,h=6,units = 'in',res = 1280)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

# selectK(cellchat, pattern = "outgoing")
nPatterns<-2
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "outgoing",
                                          k = nPatterns,
                                          width = 5, height = 9)
pdf('07.outgoing_pattern.pdf',w=10,h=6)
identifyCommunicationPatterns(cellchat, 
                              pattern = "outgoing", 
                              k = nPatterns,
                              width = 5, height = 9)
dev.off()
png('07.outgoing_pattern.png',w=10,h=6,units = 'in',res = 600)
identifyCommunicationPatterns(cellchat, 
                              pattern = "outgoing", 
                              k = nPatterns,
                              width = 5, height = 9)
dev.off()

pdf('08.outgoing_pattern_river.pdf',w=10,h=6)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
png('08.outgoing_pattern_river.png',w=10,h=6,units = 'in',res = 1280)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

## 6.2 疾病组---------------------
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./05_cellchat")){
  dir.create("./05_cellchat")
}
setwd("./05_cellchat")

library(CellChat) 
library(Seurat)
library(ggplot2)

load('../04_SingleR/TSNE/TSNE_finally.Rdata')
scRNA <- TSNE
data.input <- scRNA@assays$RNA$counts # normalized data matrix
meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$group)
cell.use = rownames(meta)[meta$group == "disease"]
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$celltype)
unique(meta$celltype)

## CellChat的输入需要的matrix和meta已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
## 如果CellChat中没有meta信息，需要手动添加
# cellchat <- addMeta(cellchat, meta = meta) ## 增加meta信息
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将"celltype"设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents <- droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents), unique(cellchat@idents)))
groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量
groupSize

## 基于配体-受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的，这里我们用人
showDatabaseCategory(CellChatDB) # 可视化数据库的基本构成

library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  # 展示互作记录
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use # 将提取出来的数据库内容加载到cellchat对象中

cellchat <- subsetData(cellchat, features = NULL) # 如果有感兴趣的基因可以填写到features这里
cellchat <- identifyOverExpressedGenes(cellchat) # 首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat) # 寻找高表达的通路
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human) # 将基因映射为蛋白

## 计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat,type = "truncatedMean", trim = 0.05, population.size = TRUE)
# cellchat <- computeCommunProb(cellchat) # 有个参数为raw.use = T，意味着上一步基因映射为蛋白这一步就没用了，计算胞间通讯概率时会用原始数据，如果为F，就会用映射的蛋白做计算
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
df.net <- subsetCommunication(cellchat) # 将细胞通讯预测结果以数据框的形式取出

if (! dir.exists("./disease")){
  dir.create("./disease")
}
setwd("./disease")

write.csv(df.net, 'df.net.csv')

# view(df.net)
# DT::datatable(df.net)

## 信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))

## 相互作用数目
pdf('04.number_of_interactions.pdf',w=6,h=5,family='serif')
netVisual_heatmap(cellchat)
dev.off()
png('04.number_of_interactions.png',w=600,h=500,family='serif')
netVisual_heatmap(cellchat)
dev.off()
#> Do heatmap based on a merged object
# netVisual_heatmap(cellchat, measure = "weight")


pdf('01.net_circle_number.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

png('01.net_circle_number.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
# 
pdf('02.net_circle_weight.pdf',w=7,h=9,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

png('02.net_circle_weight.png',w=7,h=9,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)

#需要指定受体细胞和配体细胞
cellchat@data.signaling
pdf('03.single_circle.pdf',w=9,h=14)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
png('03.single_circle.png',w=9,h=14,units = 'in',res = 600)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()


library(NMF)
library(ggalluvial)
# selectK(cellchat, pattern = "incoming")

nPatterns<-2
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "incoming",
                                          k = nPatterns,
                                          width = 5, height = 9)
pdf('05.incoming_pattern.pdf',w=10,h=6)
identifyCommunicationPatterns(cellchat, 
                              pattern = "incoming", 
                              k = nPatterns,
                              width = 5, height = 9)
dev.off()
png('05.incoming_pattern.png',w=10,h=6,units = 'in',res = 600)
identifyCommunicationPatterns(cellchat, 
                              pattern = "incoming", 
                              k = nPatterns,
                              width = 5, height = 9)
dev.off()

pdf('06.incoming_pattern_river.pdf',w=10,h=6)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
png('06.incoming_pattern_river.png',w=10,h=6,units = 'in',res = 1280)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

# selectK(cellchat, pattern = "outgoing")
nPatterns<-2
cellchat <- identifyCommunicationPatterns(cellchat,
                                          pattern = "outgoing",
                                          k = nPatterns,
                                          width = 5, height = 9)
pdf('07.outgoing_pattern.pdf',w=10,h=6)
identifyCommunicationPatterns(cellchat, 
                              pattern = "outgoing", 
                              k = nPatterns,
                              width = 5, height = 9)
dev.off()
png('07.outgoing_pattern.png',w=10,h=6,units = 'in',res = 600)
identifyCommunicationPatterns(cellchat, 
                              pattern = "outgoing", 
                              k = nPatterns,
                              width = 5, height = 9)
dev.off()

pdf('08.outgoing_pattern_river.pdf',w=10,h=6)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
png('08.outgoing_pattern_river.png',w=10,h=6,units = 'in',res = 1280)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

# 7. GSVA分析--------------------------------
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./06_GSVA")){
  dir.create("./06_GSVA")
}
setwd("./06_GSVA")
load('../04_SingleR/TSNE/TSNE_finally.Rdata')

library(GSEABase)
library(GSVA)
library(Seurat)
library(msigdbr)

# hallmark
hall <- msigdbr(species = "Homo sapiens", category = "H") 
hall <- subset(hall, select = c("gs_name","gene_symbol")) %>% as.data.frame()
hall <- split(hall$gene_symbol, hall$gs_name)

TSNE <- TSNE[, TSNE@meta.data$celltype %in% 'Cardiomyocyte']
table(TSNE@meta.data$group)
Idents(TSNE) <- "group"

## 使用AverageExpression函数提取分组平均表达矩阵
expr <- AverageExpression(TSNE, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0, ] #选取非零基因
expr <- as.matrix(expr)
head(expr)

## GSVA
gsva.hall <- gsva(expr, hall, method = "gsva")
gsva.hall.df <- data.frame(Genesets = rownames(gsva.hall), 
                           gsva.hall, check.names = F)

## 5.热图
# hallmarker
dim(gsva.hall)
# 50 2
pdf(paste0('01.','Cardiomyocyte','_gsva_heatmap.pdf'), width = 9, height = 10)
pheatmap::pheatmap(gsva.hall,
                   show_colnames = T,
                   scale = "none", angle_col = "45",
                   fontsize = 14, #  热图字体大小
                   fontsize_row = 12, #  热图行名字体大小，默认为fontsize
                   fontsize_col = 12, #  热图列名字体大小，默认为fontsize
                   border = F,
                   color = colorRampPalette(c("#4DBBD5", "white", "#E64B35" ))(50)
)

dev.off()

png(paste0('01.','Cardiomyocyte','_gsva_heatmap.png'), width = 9, height = 10,family='Times',units='in',res=600)
pheatmap::pheatmap(gsva.hall,
                   show_colnames = T,
                   scale = "none", angle_col = "45",
                   fontsize = 14, #  热图字体大小
                   fontsize_row = 12, #  热图行名字体大小，默认为fontsize
                   fontsize_col = 12, #  热图列名字体大小，默认为fontsize
                   border = F,
                   color = colorRampPalette(c("#4DBBD5", "white", "#E64B35" ))(50)
)

dev.off()

# 8. 心肌细胞分亚群-----------
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./07_subtype")){
  dir.create("./07_subtype")
}
setwd("./07_subtype")

library(Seurat)
library(tidyverse)
library(rstatix)

load('../04_SingleR/TSNE/TSNE_finally.Rdata')

table(TSNE$celltype)
TSNE.car <- TSNE[, TSNE@meta.data$celltype %in% 'Cardiomyocyte', ]
sce.car <- TSNE.car
# sce.car <- NormalizeData(sce.car, normalization.method = "LogNormalize", scale.factor = 10000)
# sce.car <- FindVariableFeatures(sce.car, selection.method = 'vst', nfeatures = 2000)
sce.car <- ScaleData(sce.car)
sce.car <- RunPCA(sce.car, features = VariableFeatures(object = sce.car))

sce.car <- JackStraw(sce.car, num.replicate = 100, dims = 50)
sce.car <- ScoreJackStraw(sce.car, dims = 1:50)
JackStrawPlot(sce.car, dims = 1:50)
ElbowPlot(sce.car, ndims = 50)

sce.car <- FindNeighbors(sce.car, dims = 1:30)
sce.car <- FindClusters(sce.car, resolution = 0.2)
table(sce.car@meta.data$seurat_clusters)
sce.car <- RunTSNE(sce.car, dims = 1:30)
p <- DimPlot(sce.car, reduction = "tsne", group.by = "seurat_clusters", label = T) + ggtitle(paste0('Subtype of ',  'Cardiomyocyte'))
pdf(file = paste0('01.subtype(Cardiomyocyte).pdf'),w=8,h=6)
print(p)
dev.off()
png(file = paste0('01.subtype(Cardiomyocyte).png'),w=8,h=6,units = 'in', res = 600)
print(p)
dev.off()
all.markers <- FindAllMarkers(sce.car, 
                              only.pos = TRUE,
                              min.pct = 0.1,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              return.thresh = 0.05)
### 每个簇中的标记基因
marker_gene <- all.markers %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC)
write.csv(marker_gene, file = paste0('Cardiomyocyte', '_marker(subtype_car).csv'))
marker_gene <- marker_gene[!duplicated(marker_gene$gene), ]
marker_gene <- marker_gene[marker_gene$gene != 'XIST', ]
marker_gene <- marker_gene[marker_gene$gene != 'TNC', ]

theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times",face='bold'))

p2 <- DotPlot(sce.car, features = marker_gene$gene, group.by = 'seurat_clusters', assay = 'RNA') + 
  RotatedAxis() + 
  theme.set+
  labs(title="", y="Subtype", x = "",size=40)+
  scale_color_gradientn(colours = c('gray','#87CEFA','#CAB2D6','#6A3D9A'))
pdf(file = paste0('0', 2, '.marker.exp(',  'Cardiomyocyte', ').pdf'),w=7,h=5)
print(p2)
dev.off()
png(file = paste0('0', 2, '.marker.exp(',  'Cardiomyocyte', ').png'),w=7,h=5,units = 'in', res = 600)
print(p2)
dev.off()

system.time(save(sce.car, file = paste0("TSNE", '_Cardiomyocyte', ".Rdata")))
load('./TSNE_Cardiomyocyte.Rdata')
marker_gene <- read.csv('./Cardiomyocyte_marker(subtype_car).csv', row.names = 1, check.names = F)

p <- DimPlot(sce.car, reduction = "tsne", group.by = "seurat_clusters", label = T) + ggtitle(paste0('Subtype of ',  'Cardiomyocyte'))

sce.car$seurat_clusters <- ifelse(sce.car$seurat_clusters == 0, 'FADS3_cluster',
                                  ifelse(sce.car$seurat_clusters == 1, 'THBS4_cluster',
                                         ifelse(sce.car$seurat_clusters == 2, 'COL3A1_cluster',
                                                ifelse(sce.car$seurat_clusters == 3, 'ADH1B_cluster',
                                                       ifelse(sce.car$seurat_clusters == 4, 'PLA2G2A_cluster', 'FAM179A_cluster')))))

table(sce.car$seurat_clusters)
system.time(save(sce.car, file = paste0("TSNE", '_Cardiomyocyte_marker', ".Rdata")))

p3 <- DimPlot(sce.car, reduction = "tsne", group.by = "seurat_clusters", label = T) + ggtitle(paste0('Subtype of ',  'Cardiomyocyte'))

library(cowplot)
combined_plot <- plot_grid(p, p3, 
                           ncol = 2, 
                           # rel_heights = c(1, 4), 
                           align = 'v')
combined_plot

ggsave("03.combined_plot.pdf", plot = combined_plot, width = 16, height = 8)
ggsave("03.combined_plot.png", plot = combined_plot, width = 16, height = 8)

i <- 2
unique(sce.car$orig.ident)[i]
Cardiomyocyte_sample_proportion <- data.frame()
for (i in c(1 : length(unique(sce.car$orig.ident)))) {
  temp1 <- subset(sce.car, orig.ident == unique(sce.car$orig.ident)[i])
  temp2 <- temp1@meta.data
  temp2$seurat_clusters <- factor(temp2$seurat_clusters, levels = c('FADS3_cluster', 'THBS4_cluster', 'COL3A1_cluster', 'ADH1B_cluster', 'PLA2G2A_cluster', 'FAM179A_cluster'))
  table(temp2$seurat_clusters)
  temp3 <- as.data.frame(sort(table(temp2$seurat_clusters)))
  for (x in c(1 : nrow(temp3))) {
    temp3$proportion[x] <- temp3$Freq[x]/sum(temp3$Freq)
    temp3$Celltype[x] <- paste0('Cardiomyocyte ', temp3$Var1[x])
  }
  temp3 <- dplyr::select(temp3, -Freq, -Var1)
  temp3$sample <- unique(TSNE$orig.ident)[i]
  
  Cardiomyocyte_sample_proportion <- rbind(Cardiomyocyte_sample_proportion, temp3)
}
Cardiomyocyte_sample_proportion$group <- ifelse(grepl('control', Cardiomyocyte_sample_proportion$sample), 'control', 'disease')

wilcox_res <- Cardiomyocyte_sample_proportion%>%
  group_by(Celltype)%>%
  wilcox_test(proportion ~ group)%>%
  adjust_pvalue(method = 'fdr')

library(RColorBrewer)
ggplot(Cardiomyocyte_sample_proportion,aes(x = sample, y= proportion, fill=Celltype))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(10, "Set3"))+
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times', angle = 45, vjust = 1, hjust = 1),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="Sample",y="Proportion",fill="")

ggsave('04.Cardiomyocyte_subtype_proportion_sample.pdf',w=12,h=8)
ggsave('04.Cardiomyocyte_subtype_proportion_sample.png',w=12,h=8)

ggplot(Cardiomyocyte_sample_proportion, aes(x = Celltype, y = proportion, fill = group)) +
  # geom_violin(trim=F, color="black", aes(fill = group)) + #绘制小提琴图, “color”设置小提琴图的轮廓线的颜色(不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值), 则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom = "errorbar",
               width = 0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(aes(x = Celltype, y = proportion, fill = group),
               width = 0.2,
               position = position_dodge(0.9), 
               outlier.shape = NA, 
               outlier.colour = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values = c('gold', "#355783"), name = "Group")+
  labs(title = "", x = "", y = "Proportion", size = 20) +
  stat_compare_means(data = Cardiomyocyte_sample_proportion,
                     mapping = aes(group = group),
                     label = "p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        axis.text.x = element_text(angle = 45, hjust=1, colour = "black", face = "bold", size = 10), 
        axis.text.y = element_text(hjust = 0.5, colour ="black", face="bold", size=12), 
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave('05.Cardiomyocyte_subtype_proportion_DE.pdf',w=10,h=8)
ggsave('05.Cardiomyocyte_subtype_proportion_DE.png',w=10,h=8)

sce.car$cluter.group <- paste(sce.car$seurat_clusters, sce.car$group, sep = "_")
sce.car$seurat_clusters <- Idents(sce.car)
Idents(sce.car) <- "cluter.group"

table(sce.car$seurat_clusters)
table(sce.car$cluter.group)

DE_res <- FindMarkers(sce.car,
                      ident.1 = paste0('FADS3_cluster_disease'),
                      ident.2 = paste0('FADS3_cluster_control'),
                      test.use = 'wilcox',
                      logfc.threshold = 0.1,
                      min.pct = 0.1,
                      only.pos = FALSE,
                      verbose = FALSE)
write.csv(DE_res, paste0("FADS3_cluster_DE.csv"))

DE_res$change <- as.factor(
  ifelse(DE_res$p_val < 0.05 & abs(DE_res$avg_log2FC) > 0.1,
         ifelse(DE_res$avg_log2FC > 0.1,'Up','Down'),'Not'))
table(DE_res$change)

DE_sig_res <- DE_res[DE_res$p_val < 0.05, ]
write.csv(DE_sig_res, paste0("FADS3_cluster_DE_sig.csv"))

dat_rep <- DE_res[rownames(DE_res) %in%
                 rownames(rbind(head(DE_sig_res[order(DE_sig_res$avg_log2FC, decreasing = T), ], 10),
                                head(DE_sig_res[order(DE_sig_res$avg_log2FC, decreasing = F), ], 10))), ]

volcano_plot <- ggplot(data = DE_res, aes(x = avg_log2FC, y = -log10(p_val), color = change)) +
  scale_color_manual(values = c("seagreen", "darkgray","firebrick3")) + # 设置颜色
  geom_point(size = 1.2, alpha = 0.4, na.rm = T) + # 设置点的大小和透明度
  theme_bw(base_size = 12, base_family = "Times") + # 添加主题
  geom_vline(xintercept = c(-0.1,0.1), # 给x坐标轴的-1,1处添加竖线
             lty = 4, # lty参数设置竖线的形式（直线，虚线等）
             col = "black", # 竖线颜色设置为黑色
             lwd = 0.4 # 设置线条的宽度
  ) +
  # scale_x_continuous(limits = c(-10,15)) + # x坐标轴范围
  scale_y_continuous(limits = c(0,10)) + # y坐标轴范围
  geom_hline(yintercept = -log10(0.05), # 给x坐标轴的-1,1处添加竖线
             lty = 4, # lty参数设置竖线的形式（直线，虚线等）
             col = "black", # 竖线颜色设置为黑色
             lwd = 0.4 # 设置线条的宽度
  ) +
  theme(
    # legend.position = 'right',  # 标签的位置
    legend.position = c(1,1),
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(face="bold",color="black",family = "Times",size=12),
    plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
    axis.text.x = element_text(face = "bold",color = "black",size = 12),
    axis.text.y = element_text(face = "bold",color = "black",size = 12),
    axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
    axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
    plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 12, face = "italic", colour = "black"),
    legend.justification = c(1,1)
  ) +
  geom_label_repel(data = dat_rep, 
                   aes(label = rownames(dat_rep)), 
                   max.overlaps = 20, 
                   size = 3, 
                   box.padding = unit(0.5, "lines"), 
                   min.segment.length = 0,
                   point.padding = unit(0.8, "lines"),
                   segment.color = "black", 
                   color = "black",
                   show.legend = FALSE 
  ) +
  labs(x = "log2(Fold_Change)",y = "-log10(p.value)",
       title = "Disease vs Control",
       subtitle = paste(sprintf('P.value:  %.2f;', 0.05),
                        sprintf("log2FC: %.2f;", 0.1),
                        sprintf('Up: %1.0f; Down: %1.0f;', table(DE_sig_res$change)[3], table(DE_sig_res$change)[1]),
                        sprintf('Total: %1.0f', dim(DE_sig_res)[1])))

volcano_plot

ggsave(paste0('06.FADS3_cluster_volcano.png'), volcano_plot, width = 8, height = 6)
ggsave(paste0('06.FADS3_cluster_volcano.pdf'), volcano_plot, width = 8, height = 6, dpi = 600)

# 9. 心肌亚群细胞富集分析--------------------------
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./08_ReactomeGSA")){
  dir.create("./08_ReactomeGSA")
}
setwd("./08_ReactomeGSA")

library(ReactomeGSA)

load('../07_subtype/TSNE_Cardiomyocyte.Rdata')

ReactomeGSA_result <- analyse_sc_clusters(sce.car, verbose = TRUE)
colnames(ReactomeGSA_result@results$Seurat$pathways) <- c('Pathway', 'Name', 'FADS3_cluster', 'THBS4_cluster',
                                                          'COL3A1_cluster', 'ADH1B_cluster', 'PLA2G2A_cluster', 'FAM179A_cluster')

plot_gsva_heatmap(ReactomeGSA_result, max_pathways = 10, margins = c(12,25))
pdf(file = paste0('01.ReactomeGSA_Cardiomyocyte.pdf'),w=8,h=6)
plot_gsva_heatmap(ReactomeGSA_result, max_pathways = 10, margins = c(12,25))
dev.off()
png(file = paste0('01.ReactomeGSA_Cardiomyocyte.png'),w=8,h=6,units = 'in', res = 600)
plot_gsva_heatmap(ReactomeGSA_result, max_pathways = 10, margins = c(12,25))
dev.off()

pathway_res <- ReactomeGSA_result@results$Seurat$pathways %>% as.data.frame()
write.csv(pathway_res, 'pathway_res.csv')


# 10. 拟时序分析------------------------
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./09_Trajectory")){
  dir.create("./09_Trajectory")
}
setwd("./09_Trajectory")

load('../04_SingleR/TSNE/TSNE_finally.Rdata')

TSNE <- TSNE[, TSNE@meta.data$celltype %in% 'Cardiomyocyte']
TSNE$celltype.group <- paste(TSNE$celltype, TSNE$group, sep = "_")
TSNE$celltype <- Idents(TSNE)
Idents(TSNE) <- "celltype.group"

table(TSNE$celltype)
table(TSNE$celltype.group)

## 组间差异基因
DE_gene <- FindMarkers(TSNE, 
                       ident.1 = paste0('Cardiomyocyte', '_disease'), 
                       ident.2 = paste0('Cardiomyocyte', '_control'), 
                       test.use = 'wilcox', 
                       min.pct = 0.1, 
                       verbose = FALSE, 
                       logfc.threshold = 0.1
)

if (! dir.exists(paste0('Cardiomyocyte'))){
  dir.create(paste0('Cardiomyocyte'))
}
setwd(paste0('Cardiomyocyte'))

group_DE_gene <- DE_gene[DE_gene$p_val < 0.05 & abs(DE_gene$avg_log2FC) > 0.1, ]
write.csv(group_DE_gene, file = './group_DE_gene.csv')

## XX细胞亚型间差异基因
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./09_Trajectory")){
  dir.create("./09_Trajectory")
}
setwd("./09_Trajectory")

load(paste0('../07_subtype/TSNE_Cardiomyocyte.Rdata'))
table(sce.car$seurat_clusters)
new.cluster.ids = c(
  "0" = "FADS3_cluster",  
  "1" = "THBS4_cluster",  
  "2" = "COL3A1_cluster",  
  "3" = "ADH1B_cluster", 
  "4" = "PLA2G2A_cluster",   
  "5" = "FAM179A_cluster"
)  
sce.car <- RenameIdents(sce.car, new.cluster.ids)
sce.car$seurat_clusters = Idents(sce.car)
table(sce.car$seurat_clusters)
cluster_num <- length(unique(sce.car$seurat_clusters))
slice_num <- length(rownames(group_DE_gene))/cluster_num
sub_DE_gene <- FindAllMarkers(sce.car, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = 'wilcox', return.thresh = 0.01)
sub_DE_gene <- sub_DE_gene %>%
  group_by(cluster) %>%
  slice_max(n=round(slice_num, 0), order_by = avg_log2FC)

if (! dir.exists(paste0('Cardiomyocyte'))){
  dir.create(paste0('Cardiomyocyte'))
}
setwd(paste0('Cardiomyocyte'))

write.csv(sub_DE_gene, file = './sub_DE_gene.csv')

## 细胞拟时序
library(Seurat)
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/")
if (! dir.exists("./09_Trajectory")){
  dir.create("./09_Trajectory")
}
setwd("./09_Trajectory")

devtools::load_all('/data/nas1/xuguangxin/pipeline/monocle/')
load(paste0('../07_subtype/TSNE_Cardiomyocyte.Rdata')) # 加载数据
new.cluster.ids = c(
  "0" = "FADS3_cluster",  
  "1" = "THBS4_cluster",  
  "2" = "COL3A1_cluster",  
  "3" = "ADH1B_cluster", 
  "4" = "PLA2G2A_cluster",   
  "5" = "FAM179A_cluster"
)  
sce.car <- RenameIdents(sce.car, new.cluster.ids)
sce.car$seurat_clusters = Idents(sce.car)
table(sce.car$seurat_clusters)
scRNA1 <- sce.car
data <- as(as.matrix(scRNA1@assays$RNA$counts), 'sparseMatrix')
# count矩阵
pd <- new('AnnotatedDataFrame', data = scRNA1@meta.data)
# meta表转成特定格式
fData <- data.frame(gene_short_name = row.names(data), 
                    row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# 基因名表转成特定格式
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
mycds <- detectGenes(mycds, min_expr = 0.1)
# 组间差异基因和亚型间差异基因
if (! dir.exists(paste0('Cardiomyocyte'))){
  dir.create(paste0('Cardiomyocyte'))
}
setwd(paste0('Cardiomyocyte'))

group_DE_gene <- read.csv('./group_DE_gene.csv')
group_DE_gene <- group_DE_gene$X %>% as.data.frame()
sub_DE_gene <- read.csv('./sub_DE_gene.csv')
sub_DE_gene <- sub_DE_gene$gene %>% as.data.frame()
genes <- rbind(group_DE_gene, sub_DE_gene)
markers.gene <- genes$.
mycds <- setOrderingFilter(mycds, markers.gene)
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
save(mycds, file = "mycds_order.t2d.RData")
#mycds=orderCells(mycds,root_state = 8) 

# load('../NK cells/mycds_order.t2d.RData')
mycds <- orderCells(mycds,root_state = 4)
table(sce.car$group)

p1 <- plot_cell_trajectory(mycds, 
                           color_by = "Pseudotime",
                           cell_size=0.5,
                           show_backbone=TRUE,
                           show_tree=TRUE, 
                           backbone_color="black", 
                           markers=NULL, 
                           use_color_gradient = FALSE,
                           markers_linear = FALSE,
                           show_cell_names=FALSE,
                           show_state_number = FALSE,
                           cell_link_size=0.75,
                           cell_name_size=2,
                           state_number_size = 2.9,
                           show_branch_points=TRUE,
                           theta = 0) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1
ggsave(paste0('01.', 'Cardiomyocyte', '_Trajectory.png'), p1, width = 6, height = 6, units = "in", dpi = 300)
ggsave(paste0('01.', 'Cardiomyocyte', '_Trajectory.pdf'), p1, width = 6, height = 6, units = "in", dpi = 300)

p2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters",cell_size=0.5,show_backbone=TRUE) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2
ggsave(paste0('02.', 'Cardiomyocyte', '_State_Trajectory.png'), p2, width = 6, height = 6, units = "in", dpi = 300)
ggsave(paste0('02.', 'Cardiomyocyte', '_State_Trajectory.pdf'), p2, width = 6, height = 6, units = "in", dpi = 300)

p3 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters",cell_size=0.5,show_backbone=TRUE) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_wrap(~group,nrow = 1)
p3
ggsave(filename = paste0('03.', 'Cardiomyocyte', '_State_Trajectory(group).png'),p3,w=10,h=6, units = "in", dpi = 300)
ggsave(filename = paste0('03.', 'Cardiomyocyte', '_State_Trajectory(group).pdf'),p3,w=10,h=6, units = "in", dpi = 300)

p4 <- plot_cell_trajectory(mycds, color_by = "group",cell_size=0.5,show_backbone=TRUE) +
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p4
ggsave(filename = paste0('04.', 'Cardiomyocyte', '_State_Trajectory(group_merge).png'),p4,w=10,h=6, units = "in", dpi = 300)
ggsave(filename = paste0('04.', 'Cardiomyocyte', '_State_Trajectory(group_merge).pdf'),p4,w=10,h=6, units = "in", dpi = 300)


hubgene <- read.csv('../../../05_expression_verification/hub_gene.csv', row.names = 1)
colnames(hubgene) <- 'symbol'
cds_subset <- mycds[hubgene$symbol,]
p5<-plot_genes_in_pseudotime(cds_subset ,color_by = "Pseudotime")+
  theme(legend.position = "right")+
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p6<-plot_genes_in_pseudotime(cds_subset ,color_by = "seurat_clusters")+
  theme(legend.position = "right")+
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p7<-plot_genes_in_pseudotime(cds_subset ,color_by = "group")+
  theme(legend.position = "right")+
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

plotc <- p5|p6|p7
plotc
ggsave(filename = paste0('05.', 'Cardiomyocyte', 'Trajectory_gene.png'),plotc,w=14,h=10, units = "in", dpi = 300)
ggsave(filename = paste0('05.', 'Cardiomyocyte', 'Trajectory_gene.pdf'),plotc,w=14,h=10, units = "in", dpi = 300)

library(patchwork)
patch <- p1 + p2 + p3 + p4
patch
ggsave(filename = paste0('06.', 'Cardiomyocyte', 'merge_Trajectory.png'),patch,w=14,h=10, units = "in", dpi = 300)
ggsave(filename = paste0('06.', 'Cardiomyocyte', 'merge_Trajectory.pdf'),patch,w=14,h=10, units = "in", dpi = 300)

my_cds_subset <- mycds
# 拟时序数据和细胞位置在pData 中
head(pData(my_cds_subset))
# 这个differentialGeneTest会比较耗费时间，测试每个基因的拟时序表达
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 4 )#cores调用的核心数
saveRDS(my_pseudotime_de,file = 'my_pseudotime.de.rds')
head(my_pseudotime_de)

# # cluster the top 50 genes that vary as a function of pseudotime
# my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
# gene_to_cluster <- gene_to_cluster[,1]
gene_to_cluster <- hubgene$symbol
colnames(pData(my_cds_subset))
table(pData(my_cds_subset)$seurat_clusters,pData(my_cds_subset)$State) 
ac=pData(my_cds_subset)[c('celltype','State','Pseudotime')]
head(ac)

# 这个热图绘制的并不是纯粹的细胞基因表达量矩阵，而是被 Pseudotime 好了的100列，50行的矩阵

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],# num_clusters = 2, # add_annotation_col = ac,
                                                 cluster_rows = FALSE,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
ggsave(file = paste0('07.', 'Cardiomyocyte', 'Trajectory_gene(heatmap).pdf'),my_pseudotime_cluster,w=7,h=6)
ggsave(file = paste0('07.', 'Cardiomyocyte', 'Trajectory_gene(heatmap).png'),my_pseudotime_cluster,w=7,h=6,units = 'in',dpi = 300)

