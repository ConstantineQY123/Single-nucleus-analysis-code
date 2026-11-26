
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/")
if(!dir.exists('./scRNA_GSVA')){dir.create('./scRNA_GSVA')}
setwd('./scRNA_GSVA')

# /data/nas1/xuguangxin/project/39_BJTC-630-4/


load('/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/04_SingleR/TSNE/TSNE_finally.Rdata')

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
# pdf(paste0('01.','Cardiomyocyte','_gsva_heatmap.pdf'), width = 9, height = 10)
# pheatmap::pheatmap(gsva.hall,
#                    show_colnames = T,
#                    scale = "none", angle_col = "45",
#                    fontsize = 14, #  热图字体大小
#                    fontsize_row = 12, #  热图行名字体大小，默认为fontsize
#                    fontsize_col = 12, #  热图列名字体大小，默认为fontsize
#                    border = F,
#                    color = colorRampPalette(c("#4DBBD5", "white", "#E64B35" ))(50)
# )
# 
# dev.off()
# 
# png(paste0('01.','Cardiomyocyte','_gsva_heatmap.png'), width = 9, height = 10,family='Times',units='in',res=600)
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



gsva_out <- data.frame(
 Hallmark      = rownames(gsva.hall),
 Control       = gsva.hall[, "control"],
 Disease       = gsva.hall[, "disease"],
 Diff          = gsva.hall[, "disease"] - gsva.hall[, "control"],
 log2FC        = log2((gsva.hall[, "disease"] + 1e-6) / (gsva.hall[, "control"] + 1e-6)),
 stringsAsFactors = FALSE
) %>%
 dplyr::mutate(
  Direction = ifelse(Diff > 0, "Up in Disease", "Down in Disease")
 ) %>%
 dplyr::arrange(dplyr::desc(abs(Diff)))

write.csv(gsva_out,
          "GSVA_Hallmark_control_vs_disease.csv",
          row.names = FALSE)

print(head(gsva_out, 10))
