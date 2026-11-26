
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/")
if(!dir.exists('./scRNA_examine')){dir.create('./scRNA_examine')}
setwd('./scRNA_examine')

# /data/nas1/xuguangxin/project/39_BJTC-630-4/

# 1 双细胞检测----------------
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/scRNA_examine")
if (!dir.exists("./01_Doublet")) dir.create("./01_Doublet")
setwd("./01_Doublet")

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)

load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/01_QC/scRNA_qc.Rdata')
table(scRNA$orig.ident)

scRNA <- JoinLayers(scRNA)
print(Layers(scRNA))

set.seed(123)
sce <- Seurat::as.SingleCellExperiment(scRNA)
sce <- scDblFinder::scDblFinder(sce, samples = "orig.ident")

scRNA$scDblFinder.score <- SummarizedExperiment::colData(sce)$scDblFinder.score
scRNA$scDblFinder.class <- SummarizedExperiment::colData(sce)$scDblFinder.class
rm(sce)

print(table(scRNA$scDblFinder.class))
# singlet doublet 
# 58588   11954 

doublet_summary <- as.data.frame(table(scRNA$scDblFinder.class))
colnames(doublet_summary) <- c("Class", "Count")
doublet_summary$Proportion <- round(doublet_summary$Count / sum(doublet_summary$Count) * 100, 2)
write.csv(doublet_summary, "01.doublet_overall_summary.csv", row.names = FALSE)
doublet_summary
# Class Count Proportion
# 1 singlet 58588      83.05
# 2 doublet 11954      16.95

# 每个样本的双细胞比例
sample_doublet <- as.data.frame(
 prop.table(table(scRNA$orig.ident, scRNA$scDblFinder.class), margin = 1) * 100
)
colnames(sample_doublet) <- c("Sample", "Class", "Percentage")
write.csv(sample_doublet, "02.doublet_by_sample.csv", row.names = FALSE)

sample_cols <- c('#74AED4','#D3E2B7','#CFAFD4','#F7C97E',"#A6D854","#FFD92F","#9966FF","#BDD7EE","#9cB4E5","#d9BD5F","#FF7C8F","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#1B9E77","#D95F02","#7570B3","#66C2A5","#FC8D62")
sample_cols <- sample_cols[1:length(unique(scRNA$orig.ident))]

p_doublet <- VlnPlot(scRNA,
                     features = "scDblFinder.score",
                     group.by = "orig.ident",
                     cols = sample_cols,
                     pt.size = 0) + NoLegend() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8,
                                  family = "Times", face = "bold")) +
 labs(title = "scDblFinder score by sample (before filtering)")
ggsave("01.doublet_score_before.pdf", p_doublet, width = 10, height = 4)
ggsave("01.doublet_score_before.png", p_doublet, width = 10, height = 4, dpi = 600)

p_bar <- ggplot(sample_doublet, aes(x = Sample, y = Percentage, fill = Class)) +
 geom_bar(stat = "identity", position = "stack", width = 0.7) +
 scale_fill_manual(values = c("singlet" = "grey70", "doublet" = "#dbacb5")) +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10,
                                  face = "bold", family = "Times"),
       axis.text.y = element_text(size = 10, face = "bold", family = "Times"),
       axis.title  = element_text(size = 14, face = "bold", family = "Times"),
       legend.title = element_blank(),
       legend.text  = element_text(size = 12, family = "Times"),
       legend.position = "top") +
 labs(x = "Sample", y = "Percentage (%)", title = "")
ggsave("02.doublet_proportion_by_sample.pdf", p_bar, width = 10, height = 5)
ggsave("02.doublet_proportion_by_sample.png", p_bar, width = 10, height = 5, dpi = 600)

#  过滤双细胞
n_before <- ncol(scRNA)
scRNA_singlet <- subset(scRNA, subset = scDblFinder.class == "singlet")
n_after <- ncol(scRNA_singlet)

filter_res <- data.frame(
 Stage       = c("Before", "After"),
 Cell_number = c(n_before, n_after),
 Removed     = c(0, n_before - n_after),
 Removed_pct = c(0, round((n_before - n_after) / n_before * 100, 2))
)
write.csv(filter_res, "03.doublet_filter_summary.csv", row.names = FALSE)
print(filter_res)


p_doublet_after <- VlnPlot(scRNA_singlet,
                           features = "scDblFinder.score",
                           group.by = "orig.ident",
                           cols = sample_cols,
                           pt.size = 0) + NoLegend() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8,
                                  family = "Times", face = "bold")) +
 labs(title = "scDblFinder score by sample (after filtering)")
ggsave("04.doublet_score_after.pdf", p_doublet_after, width = 10, height = 4)
ggsave("04.doublet_score_after.png", p_doublet_after, width = 10, height = 4, dpi = 600)

save(scRNA_singlet, file = "01.scRNA_singlet.Rdata")

# 2 环境 RNA 污染评估（decontX）---------------------------
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/scRNA_examine/")
if (!dir.exists("./02_Ambient")) dir.create("./02_Ambient")
setwd("./02_Ambient")

library(Seurat)
library(celda)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)

load("../01_Doublet/01.scRNA_singlet.Rdata")
table(scRNA_singlet$orig.ident)
# control_1 control_2 control_3 control_4 disease_1 disease_2 disease_3 disease_4 disease_5 disease_6 disease_7 disease_8 
# 2023      1012       585      7893     33209      8680       652       170       628       354      1204      1797 
# disease_9 
# 381 

# 运行 decontX
set.seed(123)
sce <- as.SingleCellExperiment(scRNA_singlet)
sce <- celda::decontX(sce)

scRNA_singlet$decontX_contamination <- SummarizedExperiment::colData(sce)$decontX_contamination
scRNA_singlet$decontX_clusters      <- SummarizedExperiment::colData(sce)$decontX_clusters
rm(sce)

# 统计 
contam_summary <- summary(scRNA_singlet$decontX_contamination)
print(contam_summary)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002474 0.107262 0.147802 0.259236 0.400106 0.978365 

contam_df <- data.frame(Metric = names(contam_summary),
                        Value  = as.numeric(contam_summary))
write.csv(contam_df, "01.decontX_contamination_summary.csv", row.names = FALSE)

# 按样本汇总
sample_contam <- scRNA_singlet@meta.data %>%
 group_by(orig.ident) %>%
 summarise(Median = median(decontX_contamination),
           Mean   = mean(decontX_contamination),
           Q25    = quantile(decontX_contamination, 0.25),
           Q75    = quantile(decontX_contamination, 0.75),
           Max    = max(decontX_contamination),
           .groups = "drop")
write.csv(sample_contam, "02.decontX_by_sample.csv", row.names = FALSE)

#  污染率小提琴图（过滤前）
sample_cols <- c('#74AED4','#D3E2B7','#CFAFD4','#F7C97E',"#A6D854","#FFD92F","#9966FF","#BDD7EE","#9cB4E5","#d9BD5F","#FF7C8F","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#1B9E77","#D95F02","#7570B3","#66C2A5","#FC8D62")
sample_cols <- sample_cols[1:length(unique(scRNA_singlet$orig.ident))]

p_before <- VlnPlot(scRNA_singlet,
                    features = "decontX_contamination",
                    group.by = "orig.ident",
                    cols = sample_cols,
                    pt.size = 0) + NoLegend() +
 geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8,
                                  family = "Times", face = "bold")) +
 labs(title = "decontX contamination before filtering (red line = 0.3)")
ggsave("01.decontX_contamination_before.pdf", p_before, width = 10, height = 4)
ggsave("01.decontX_contamination_before.png", p_before, width = 10, height = 4, dpi = 600)

# 污染率密度分布 
p_density <- ggplot(scRNA_singlet@meta.data, aes(x = decontX_contamination)) +
 geom_density(fill = "#74AED4", alpha = 0.6, color = "black") +
 geom_vline(xintercept = 0.3, linetype = "dashed", color = "red") +
 theme_bw() +
 theme(axis.text  = element_text(size = 10, face = "bold", family = "Times"),
       axis.title = element_text(size = 14, face = "bold", family = "Times")) +
 labs(x = "decontX contamination", y = "Density",
      title = "Distribution of contamination")
ggsave("02.decontX_contamination_density.pdf", p_density, width = 8, height = 5)
ggsave("02.decontX_contamination_density.png", p_density, width = 8, height = 5, dpi = 600)

#  过滤高污染细胞（contamination < 0.3）
n_before <- ncol(scRNA_singlet)
scRNA_final <- subset(scRNA_singlet, subset = decontX_contamination < 0.3)
n_after <- ncol(scRNA_final)

filter_res <- data.frame(
 Stage       = c("Before", "After"),
 Cell_number = c(n_before, n_after),
 Removed     = c(0, n_before - n_after),
 Removed_pct = c(0, round((n_before - n_after) / n_before * 100, 2))
)
write.csv(filter_res, "03.decontX_filter_summary.csv", row.names = FALSE)
print(filter_res)
# Stage Cell_number Removed Removed_pct
# 1 Before       58588       0        0.00
# 2  After       40533   18055       30.82

#  过滤后小提琴图 
p_after <- VlnPlot(scRNA_final,
                   features = "decontX_contamination",
                   group.by = "orig.ident",
                   cols = sample_cols,
                   pt.size = 0) + NoLegend() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8,
                                  family = "Times", face = "bold")) +
 labs(title = "decontX contamination after filtering")
ggsave("04.decontX_contamination_after.pdf", p_after, width = 10, height = 4)
ggsave("04.decontX_contamination_after.png", p_after, width = 10, height = 4, dpi = 600)

save(scRNA_final, file = "02.scRNA_qc_doublet_ambient.Rdata")

# 3 批次问题----------------------
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/scRNA_examine/")
if (!dir.exists("./03_Batch")) dir.create("./03_Batch")
setwd("./03_Batch")

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

extract_dataset <- function(orig_ident) {
 gse <- regmatches(orig_ident, regexpr("GSE[0-9]+", orig_ident))
 if (length(gse) > 0 && nchar(gse) > 0) return(gse) else return("Unknown")
}


dataset_map <- c(
 "control_1" = "GSE161921", "control_2" = "GSE161921",
 "control_3" = "GSE161921", "control_4" = "GSE161921",
 "disease_1" = "GSE174691", "disease_2" = "GSE174691",
 "disease_3" = "GSE174691", "disease_4" = "GSE174691",
 "disease_5" = "GSE174691", "disease_6" = "GSE174691",
 "disease_7" = "GSE174691", "disease_8" = "GSE174691",
 "disease_9" = "GSE174691"
)

assign_dataset <- function(obj, map) {
 oi <- as.character(obj$orig.ident)
 ds <- unname(map[oi])
 miss <- is.na(ds)
 if (any(miss)) {
  ds[miss] <- unname(sapply(oi[miss], extract_dataset))
 }
 ds[ds == "Unknown"] <- oi[ds == "Unknown"]
 ds
}

##  3.1 整合前：PCA by Dataset
load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/02_Integrate/unintergrated.Rdata')

combined$dataset <- assign_dataset(combined, dataset_map)
print(table(combined$dataset))

p_before <- DimPlot(combined, reduction = "pca", group.by = "dataset",
                    pt.size = 0.1, shuffle = TRUE) +
 theme_bw() +
 theme(axis.text  = element_text(size = 10, face = "bold", family = "Times"),
       axis.title = element_text(size = 14, face = "bold", family = "Times"),
       legend.title = element_blank(),
       legend.text  = element_text(size = 10, family = "Times"),
       plot.title   = element_text(size = 14, face = "bold", family = "Times")) +
 labs(title = "Before integration")
ggsave("01.PCA_by_dataset_before_integration.pdf", p_before, width = 8, height = 6)
ggsave("01.PCA_by_dataset_before_integration.png", p_before, width = 8, height = 6, dpi = 600)


##  3.2 整合后：CCA by Dataset 
load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/02_Integrate/scRNA_intergrated.Rdata')

scRNA$dataset <- assign_dataset(scRNA, dataset_map)
print(table(scRNA$dataset))

p_after <- DimPlot(scRNA, reduction = "integrated.cca", group.by = "dataset",
                   pt.size = 0.1, shuffle = TRUE) +
 theme_bw() +
 theme(axis.text  = element_text(size = 10, face = "bold", family = "Times"),
       axis.title = element_text(size = 14, face = "bold", family = "Times"),
       legend.title = element_blank(),
       legend.text  = element_text(size = 10, family = "Times"),
       plot.title   = element_text(size = 14, face = "bold", family = "Times")) +
 labs(title = "After integration")
ggsave("02.CCA_by_dataset_after_integration.pdf", p_after, width = 8, height = 6)
ggsave("02.CCA_by_dataset_after_integration.png", p_after, width = 8, height = 6, dpi = 600)

##  3.3 最终 TSNE by Dataset 
load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/04_SingleR/TSNE/TSNE_finally.Rdata')

TSNE$dataset <- assign_dataset(TSNE, dataset_map)
print(table(TSNE$dataset))

p_final <- DimPlot(TSNE, reduction = "tsne", group.by = "dataset",
                   pt.size = 0.1, shuffle = TRUE) +
 theme_bw() +
 theme(axis.text  = element_text(size = 10, face = "bold", family = "Times"),
       axis.title = element_text(size = 14, face = "bold", family = "Times"),
       legend.title = element_blank(),
       legend.text  = element_text(size = 10, family = "Times"),
       plot.title   = element_text(size = 14, face = "bold", family = "Times")) +
 labs(title = "Final TSNE")
ggsave("03.TSNE_by_dataset_final.pdf", p_final, width = 8, height = 6)
ggsave("03.TSNE_by_dataset_final.png", p_final, width = 8, height = 6, dpi = 600)

##  3.4 批次混合指标
meta_cluster <- TSNE@meta.data
meta_cluster$cluster <- as.character(Idents(TSNE))

cluster_dataset <- meta_cluster %>%
 group_by(cluster, dataset) %>%
 summarise(n = n(), .groups = "drop") %>%
 group_by(cluster) %>%
 mutate(prop = n / sum(n)) %>%
 summarise(entropy   = -sum(prop * log2(prop)),
           max_prop  = max(prop),
           n_dataset = n_distinct(dataset),
           .groups = "drop")

write.csv(cluster_dataset, "04.cluster_dataset_mixing.csv", row.names = FALSE)
print(cluster_dataset)
# cluster       entropy max_prop n_dataset
# <chr>           <dbl>    <dbl>     <int>
#  1 Cardiomyocyte  0.791     0.763         2
# 2 Endothelial    0.0570    0.993         2
# 3 Fibroblast     0.246     0.959         2
# 4 Leukocyte      0.819     0.745         2
# 5 Neuronal       0.181     0.973         2
# 6 Smooth muscle  0.779     0.770         2

p_mix <- ggplot(meta_cluster, aes(x = cluster, fill = dataset)) +
 geom_bar(position = "fill", width = 0.8) +
 theme_bw() +
 theme(axis.text.x = element_text(size = 10, face = "bold", family = "Times"),
       axis.text.y = element_text(size = 10, face = "bold", family = "Times"),
       axis.title  = element_text(size = 14, face = "bold", family = "Times"),
       legend.title = element_blank(),
       legend.text  = element_text(size = 10, family = "Times"),
       legend.position = "top") +
 labs(x = "Cluster", y = "Proportion", title = "")
ggsave("04.cluster_dataset_proportion.pdf", p_mix, width = 10, height = 5)
ggsave("04.cluster_dataset_proportion.png", p_mix, width = 10, height = 5, dpi = 600)

# 4 综合评估（处理前后对比）---------------
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/scRNA_examine/")
if (!dir.exists("./04_Evaluation")) dir.create("./04_Evaluation")
setwd("./04_Evaluation")


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/01_QC/scRNA_qc.Rdata')
scRNA_qc <- scRNA; rm(scRNA)
scRNA_qc <- JoinLayers(scRNA_qc)   # 关键：原始 QC 对象是多 layer，4.4 取 data 前必须合并

load('../01_Doublet/01.scRNA_singlet.Rdata')                 # scRNA_singlet
load('../02_Ambient/02.scRNA_qc_doublet_ambient.Rdata')       # scRNA_final

load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/04_SingleR/TSNE/TSNE_finally.Rdata')
celltype_map <- setNames(TSNE$celltype, colnames(TSNE))

load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/07_subtype/TSNE_Cardiomyocyte_marker.Rdata')
subcluster_map <- setNames(sce.car$seurat_clusters, colnames(sce.car))

safe_assign_meta <- function(obj, vec, col_name) {
 vals <- vec[colnames(obj)]
 vals <- unname(vals)
 obj@meta.data[[col_name]] <- vals
 obj
}

# 4.1 把过滤信息回贴到原始 QC 数据 
dbl_map    <- setNames(scRNA_singlet$scDblFinder.class,    colnames(scRNA_singlet))
contam_map <- setNames(scRNA_final$decontX_contamination,  colnames(scRNA_final))

scRNA_qc <- safe_assign_meta(scRNA_qc, dbl_map,        "scDblFinder.class")
scRNA_qc <- safe_assign_meta(scRNA_qc, contam_map,     "decontX_contamination")
scRNA_qc <- safe_assign_meta(scRNA_qc, celltype_map,   "celltype")
scRNA_qc <- safe_assign_meta(scRNA_qc, subcluster_map, "subcluster")

doublet_cells   <- colnames(scRNA_qc)[is.na(scRNA_qc$scDblFinder.class)] 
ambient_cells   <- setdiff(colnames(scRNA_singlet), colnames(scRNA_final))
retained_cells  <- colnames(scRNA_final) 

scRNA_qc$status <- ifelse(
 colnames(scRNA_qc) %in% doublet_cells,  "Doublet_removed",
 ifelse(colnames(scRNA_qc) %in% ambient_cells, "Ambient_removed",
        "Retained")
)

print(table(scRNA_qc$status, useNA = "ifany"))
# Ambient_removed Doublet_removed        Retained 
# 18055           11954           40533


#  4.2 被删除细胞在各细胞类型中的分布 
celltype_status <- scRNA_qc@meta.data %>%
 filter(!is.na(celltype)) %>%
 group_by(celltype, status) %>%
 summarise(n = n(), .groups = "drop") %>%
 group_by(celltype) %>%
 mutate(prop = n / sum(n))
celltype_status
# celltype      status              n   prop
# <fct>         <chr>           <int>  <dbl>
#  1 Endothelial   Ambient_removed  3353 0.0953
# 2 Endothelial   Doublet_removed  4766 0.135 
# 3 Endothelial   Retained        27077 0.769 
# 4 Cardiomyocyte Ambient_removed  1182 0.107 
# 5 Cardiomyocyte Doublet_removed  1011 0.0918
# 6 Cardiomyocyte Retained         8822 0.801 
# 7 Neuronal      Ambient_removed  4630 0.518 
# 8 Neuronal      Doublet_removed  3817 0.427 
# 9 Neuronal      Retained          485 0.0543
# 10 Fibroblast    Ambient_removed  4672 0.624 
# 11 Fibroblast    Doublet_removed   475 0.0635
# 12 Fibroblast    Retained         2336 0.312 
# 13 Leukocyte     Ambient_removed  1753 0.423 
# 14 Leukocyte     Doublet_removed  1269 0.306 
# 15 Leukocyte     Retained         1123 0.271 
# 16 Smooth muscle Ambient_removed  2465 0.654 
# 17 Smooth muscle Doublet_removed   616 0.163 
# 18 Smooth muscle Retained          690 0.183 
write.csv(celltype_status, "01.celltype_status_distribution.csv", row.names = FALSE)

p_ct <- ggplot(celltype_status, aes(x = celltype, y = prop, fill = status)) +
 geom_bar(stat = "identity", position = "stack", width = 0.7) +
 scale_fill_manual(values = c("Retained" = "grey70",
                              "Doublet_removed" = "#dbacb5",
                              "Ambient_removed" = "#74AED4")) +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10,
                                  face = "bold", family = "Times"),
       axis.text.y = element_text(size = 10, face = "bold", family = "Times"),
       axis.title  = element_text(size = 14, face = "bold", family = "Times"),
       legend.title = element_blank(),
       legend.text  = element_text(size = 12, family = "Times"),
       legend.position = "top") +
 labs(x = "Celltype", y = "Proportion",
      title = "")
ggsave("01.celltype_status_distribution.pdf", p_ct, width = 10, height = 6)
ggsave("01.celltype_status_distribution.png", p_ct, width = 10, height = 6, dpi = 600)

#  4.3 处理前后细胞类型比例对比 
before_ct <- scRNA_qc@meta.data %>%
 filter(!is.na(celltype)) %>%
 group_by(celltype) %>%
 summarise(n_before = n(), .groups = "drop") %>%
 mutate(prop_before = n_before / sum(n_before))

scRNA_final <- safe_assign_meta(scRNA_final, celltype_map, "celltype")

after_ct <- scRNA_final@meta.data %>%
 filter(!is.na(celltype)) %>%
 group_by(celltype) %>%
 summarise(n_after = n(), .groups = "drop") %>%
 mutate(prop_after = n_after / sum(n_after))

ct_compare <- full_join(before_ct, after_ct, by = "celltype") %>%
 mutate(across(c(n_before, n_after, prop_before, prop_after),
               ~ ifelse(is.na(.x), 0, .x)))
write.csv(ct_compare, "02.celltype_proportion_compare.csv", row.names = FALSE)
print(ct_compare)
# celltype      n_before prop_before n_after prop_after
# <fct>            <int>       <dbl>   <int>      <dbl>
#  1 Endothelial      35196      0.499    27077     0.668 
# 2 Cardiomyocyte    11015      0.156     8822     0.218 
# 3 Neuronal          8932      0.127      485     0.0120
# 4 Fibroblast        7483      0.106     2336     0.0576
# 5 Leukocyte         4145      0.0588    1123     0.0277
# 6 Smooth muscle     3771      0.0535     690     0.0170

ct_long <- ct_compare %>%
 dplyr::select(celltype, prop_before, prop_after) %>%
 pivot_longer(-celltype, names_to = "Stage", values_to = "Proportion") %>%
 mutate(Stage = factor(Stage, levels = c("prop_before", "prop_after"),
                       labels = c("Before", "After")))

p_compare <- ggplot(ct_long, aes(x = celltype, y = Proportion, fill = Stage)) +
 geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
 scale_fill_manual(values = c("Before" = "grey70", "After" = "#dbacb5")) +
 theme_bw() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10,
                                  face = "bold", family = "Times"),
       axis.text.y = element_text(size = 10, face = "bold", family = "Times"),
       axis.title  = element_text(size = 14, face = "bold", family = "Times"),
       legend.title = element_blank(),
       legend.text  = element_text(size = 12, family = "Times"),
       legend.position = "top") +
 labs(x = "Celltype", y = "Proportion",
      title = "")
ggsave("02.celltype_proportion_compare.pdf", p_compare, width = 10, height = 6)
ggsave("02.celltype_proportion_compare.png", p_compare, width = 10, height = 6, dpi = 600)

#  4.4 关键基因在"保留 vs 删除"细胞中的表达 
target_genes <- intersect(c("COX6C", "NDUFA8", "NDUFB1"), rownames(scRNA_qc))
print(target_genes)

if (nrow(scRNA_qc@assays$RNA$data) == 0) {
 message("data layer 为空，执行 NormalizeData ...")
 scRNA_qc <- NormalizeData(scRNA_qc, verbose = FALSE)
}
cat("data layer dim:", paste(dim(scRNA_qc@assays$RNA$data), collapse = " x "), "\n")
stopifnot(all(target_genes %in% rownames(scRNA_qc@assays$RNA$data)))

fads3_cells <- colnames(scRNA_qc)[!is.na(scRNA_qc$subcluster) &
                                   scRNA_qc$subcluster == "FADS3_cluster"]
cat("FADS3_cluster cells:", length(fads3_cells), "\n")
# FADS3_cluster cells: 4802

expr_mat <- as.matrix(scRNA_qc@assays$RNA$data[target_genes, fads3_cells, drop = FALSE])
expr_df  <- as.data.frame(t(expr_mat))
expr_df$cell   <- rownames(expr_df)
expr_df$status <- scRNA_qc@meta.data[expr_df$cell, "status"]

gene_status_mean <- expr_df %>%
 pivot_longer(all_of(target_genes), names_to = "gene", values_to = "expression") %>%
 group_by(gene, status) %>%
 summarise(mean_expr = mean(expression),
           sd_expr   = sd(expression),
           n_cells   = n(), .groups = "drop")
write.csv(gene_status_mean, "03.keygene_expression_by_status.csv", row.names = FALSE)
print(gene_status_mean)
# gene   status          mean_expr sd_expr n_cells
# <chr>  <chr>               <dbl>   <dbl>   <int>
#  1 COX6C  Ambient_removed     1.62    1.40      281
# 2 COX6C  Doublet_removed     1.22    1.12      105
# 3 COX6C  Retained            1.19    1.24     4416
# 4 NDUFA8 Ambient_removed     0.417   0.906     281
# 5 NDUFA8 Doublet_removed     0.226   0.637     105
# 6 NDUFA8 Retained            0.279   0.731    4416
# 7 NDUFB1 Ambient_removed     1.25    1.32      281
# 8 NDUFB1 Doublet_removed     0.628   0.993     105
# 9 NDUFB1 Retained            0.759   1.09     4416

gene_long <- expr_df %>%
 pivot_longer(all_of(target_genes), names_to = "gene", values_to = "expression")

p_gene <- ggplot(gene_long, aes(x = gene, y = expression, fill = status)) +
 geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(0.8)) +
 geom_boxplot(width = 0.1, position = position_dodge(0.8), outlier.shape = NA) +
 scale_fill_manual(values = c("Retained" = "grey70",
                              "Doublet_removed" = "#dbacb5",
                              "Ambient_removed" = "#74AED4")) +
 theme_bw() +
 theme(axis.text  = element_text(size = 10, face = "bold", family = "Times"),
       axis.title = element_text(size = 14, face = "bold", family = "Times"),
       legend.title = element_blank(),
       legend.text  = element_text(size = 12, family = "Times"),
       legend.position = "top") +
 labs(x = "", y = "Expression (log-normalized)",
      title = "")
ggsave("03.keygene_expression_by_status.pdf", p_gene, width = 8, height = 6)
ggsave("03.keygene_expression_by_status.png", p_gene, width = 8, height = 6, dpi = 600)

#  4.5 原 DE 结果中关键基因的展示
DE_res <- read.csv('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/07_subtype/FADS3_cluster_DE.csv',
                   row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
DE_target <- DE_res[target_genes, , drop = FALSE]
DE_target$gene <- rownames(DE_target)
DE_target <- DE_target[, c("gene", setdiff(colnames(DE_target), "gene"))]
write.csv(DE_target, "04.original_DE_keygenes.csv", row.names = FALSE)
print(DE_target)
# gene        p_val avg_log2FC pct.1 pct.2 p_val_adj
# COX6C   COX6C 2.646059e-05  1.9928703 0.505 0.175 0.6251844
# NDUFA8 NDUFA8 4.261952e-02  2.5928246 0.134 0.025 1.0000000
# NDUFB1 NDUFB1 1.353164e-02  0.7758301 0.346 0.150 1.0000000

summary_report <- data.frame(
 Metric = c("Total cells (QC-passed)",
            "Doublets removed",
            "Ambient RNA removed",
            "Cells retained",
            "Celltype proportion Pearson r (before vs after)",
            "Key DE genes"),
 Value  = c(nrow(scRNA_qc@meta.data),
            sum(scRNA_qc$status == "Doublet_removed"),
            sum(scRNA_qc$status == "Ambient_removed"),
            ncol(scRNA_final),
            round(cor(ct_compare$prop_before, ct_compare$prop_after,
                      method = "pearson"), 4),
            paste(DE_target$gene, collapse = ", "))
)
write.csv(summary_report, "05.summary_report.csv", row.names = FALSE)
print(summary_report)
# Metric                 Value
# 1                         Total cells (QC-passed)                 70542
# 2                                Doublets removed                 11954
# 3                             Ambient RNA removed                 18055
# 4                                  Cells retained                 40533
# 5 Celltype proportion Pearson r (before vs after)                0.9776
# 6                                    Key DE genes COX6C, NDUFA8, NDUFB1


# 5 质控前后样本细胞数---------------
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/scRNA_examine/")
if (!dir.exists("./05_QC")) dir.create("./05_QC")
setwd("./05_QC")

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/01.scRNA_orig.Rdata')
scRNA_orig <- scRNA; rm(scRNA)          # QC 前

load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/01_QC/scRNA_qc.Rdata')
scRNA_qc <- scRNA; rm(scRNA)            # QC 后

cat("QC before cells:", ncol(scRNA_orig), "\n")
cat("QC after  cells:", ncol(scRNA_qc),   "\n")

before_tab <- as.data.frame(table(scRNA_orig$orig.ident))
colnames(before_tab) <- c("Sample", "Before")

after_tab <- as.data.frame(table(scRNA_qc$orig.ident))
colnames(after_tab) <- c("Sample", "After")

qc_table <- full_join(before_tab, after_tab, by = "Sample") %>%
 mutate(
  Removed      = Before - After,
  Retained_pct = round(After / Before * 100, 2),
  Removed_pct  = round(Removed / Before * 100, 2),
  Group        = ifelse(grepl("^control", Sample), "control", "disease")
 ) %>%
 dplyr::select(Sample, Group, Before, After, Removed, Retained_pct, Removed_pct)

write.csv(qc_table, "01.QC_cell_number_by_sample.csv", row.names = FALSE)
print(qc_table)
# Sample   Group Before After Removed Retained_pct Removed_pct
# 1  control_1 control   5508  2125    3383        38.58       61.42
# 2  control_2 control   5834  1063    4771        18.22       81.78
# 3  control_3 control  16406   593   15813         3.61       96.39
# 4  control_4 control  14019  8413    5606        60.01       39.99
# 5  disease_1 disease  53596 43642    9954        81.43       18.57
# 6  disease_2 disease  38963  9419   29544        24.17       75.83
# 7  disease_3 disease  11451   658   10793         5.75       94.25
# 8  disease_4 disease   4515   188    4327         4.16       95.84
# 9  disease_5 disease  10496   650    9846         6.19       93.81
# 10 disease_6 disease  11375   367   11008         3.23       96.77
# 11 disease_7 disease  33347  1232   32115         3.69       96.31
# 12 disease_8 disease  69721  1806   67915         2.59       97.41
# 13 disease_9 disease   7918   386    7532         4.87       95.13


# 6 心肌细胞 marker 验证---------------
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/scRNA_examine/")
if (!dir.exists("./06_Cardiomyocyte_marker")) dir.create("./06_Cardiomyocyte_marker")
setwd("./06_Cardiomyocyte_marker")


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(patchwork)

load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/04_SingleR/TSNE/TSNE_finally.Rdata')
TSNE <- JoinLayers(TSNE)

# 心肌细胞经典 marker 
cardio_markers <- c(
 "TNNT2", "MYH6",  "MYBPC3","ACTN2",  "MYL7", "RYR2",  "TTN"
)

cardio_markers <- intersect(cardio_markers, rownames(TSNE))
cat("可用心肌 marker:", paste(cardio_markers, collapse = ", "), "\n")

# 对照 marker（其他细胞类型，用于展示特异性）
other_markers <- c(
 "PECAM1", "VWF",   # Endothelial
 "POSTN",  "COL1A1",# Fibroblast
 "PTPRC",  "CD14",  # Leukocyte
 "PDGFRB", "MYH11", # Smooth muscle
 "SCN7A"            # Neuronal
)
other_markers <- intersect(other_markers, rownames(TSNE))

# VlnPlot：心肌 marker 在各细胞类型中的表达 
theme.set <- theme(
 axis.text  = element_text(size = 10, face = "bold", family = "Times"),
 axis.title = element_text(size = 12, face = "bold", family = "Times"),
 plot.title = element_text(size = 12, face = "bold", family = "Times"),
 legend.position = "none"
)

p_vln <- VlnPlot(TSNE,
                 features = cardio_markers,
                 group.by = "celltype",
                 pt.size = 0, ncol = 4) &
 theme.set &
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10,
                                  face = "bold", family = "Times"))
ggsave("01.cardiomyocyte_marker_vlnplot.pdf", p_vln, width = 16, height = 8)
ggsave("01.cardiomyocyte_marker_vlnplot.png", p_vln, width = 16, height = 8, dpi = 600)

#  6.3 DotPlot：心肌 marker + 对照 marker 
all_markers <- c(cardio_markers, other_markers)

p_dot <- DotPlot(TSNE,
                 features = all_markers,
                 group.by = "celltype",
                 assay = "RNA") +
 RotatedAxis() +
 theme_bw() +
 theme(axis.text.x = element_text(size = 9, face = "bold", family = "Times"),
       axis.text.y = element_text(size = 11, face = "bold", family = "Times"),
       axis.title  = element_text(size = 13, face = "bold", family = "Times"),
       legend.title = element_text(size = 12, face = "bold", family = "Times"),
       legend.text  = element_text(size = 10, family = "Times")) +
 scale_color_gradientn(colours = c("gray", "#87CEFA", "#CAB2D6", "#6A3D9A")) +
 labs(x = "", y = "Celltype", title = "")
ggsave("02.cardiomyocyte_marker_dotplot.pdf", p_dot, width = 12, height = 5)
ggsave("02.cardiomyocyte_marker_dotplot.png", p_dot, width = 12, height = 5, dpi = 600)
