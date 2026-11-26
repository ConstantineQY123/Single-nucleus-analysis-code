rm(list = ls()); gc()
ORIGINAL_DIR <- "/data/nas3/refinement8/luoyuyan/BJTC-630-4"
output <- file.path(ORIGINAL_DIR, "scTenifoldKnk")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

setwd(output)
getwd()

library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)
load('./TSNE_Cardiomyocyte.Rdata')
hubgene <- c('NDUFA8')
seurat_object<- sce.car
path <- getwd()
table(seurat_object$celltype)

countMatrix <- GetAssayData(seurat_object, slot = "counts")

gene_expr_counts <- t(GetAssayData(seurat_object, assay = "RNA", slot = "counts")[hubgene, , drop = FALSE])
# 转换为数据框便于统计
gene_expr_df <- as.data.frame(gene_expr_counts)

# 统计每个基因的表达情况
for(gene in hubgene){
  expr_vector <- gene_expr_df[[gene]]

  cat("\n=== 基因:", gene, "===\n")
  cat("总细胞数:", length(expr_vector), "\n")
  cat("表达细胞数 (>0):", sum(expr_vector > 0), "\n")
  cat("表达比例:", round(mean(expr_vector > 0) * 100, 2), "%\n")
  cat("平均表达量:", round(mean(expr_vector), 4), "\n")
  cat("中位数表达量:", round(median(expr_vector), 4), "\n")
  cat("最大表达量:", max(expr_vector), "\n")

#   # 按表达阈值统计
  cat("\n表达水平分布:\n")
  cat("  0 counts:", sum(expr_vector == 0), "cells\n")
  cat("  1-10 counts:", sum(expr_vector >= 1 & expr_vector <= 10), "cells\n")
  cat("  11-100 counts:", sum(expr_vector >= 11 & expr_vector <= 100), "cells\n")
  cat("  >100 counts:", sum(expr_vector > 100), "cells\n")
}
i= 'NDUFA8'
# for (i in hubgene){

cat(paste0("开始处理基因：", i, " 时间：", Sys.time(), "\n"))
result <- scTenifoldKnk(countMatrix = countMatrix, gKO = i, qc = FALSE)
save(result, file = paste0(i,"_result_scTenifoldKnk.Rdata"))
# 清理内存，避免后台内存溢出
rm(result)
gc()

cat(paste0("完成处理基因：", i, " 时间：", Sys.time(), "\n"))
# }
cat(paste0("所有基因处理完成！结束时间：", Sys.time(), "\n"))
#NDUFA8
load('NDUFA8_result_scTenifoldKnk.Rdata')
top_genes <- head(result$diffRegulation[order(-result$diffRegulation$FC)&result$diffRegulation$gene!='NDUFA8', ], 10)



ggplot(top_genes, aes(x = reorder(gene, FC), y = FC)) +
  geom_bar(
    stat = "identity",
    aes(fill = FC > 0), # 按FC正负分色
    width = 0.8
  ) +
  scale_fill_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "#E74C3C"),
    guide = "none"
  ) +
  coord_flip() +
  labs(
    title = "Top 10 Differentially Perturbed Genes",
    x = "Gene",
    y = "Fold Change (FC)",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), # 标题居中+加粗
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50"),
    plot.caption = element_text(size = 10, color = "gray50"),
    axis.title.x = element_text(size = 14, face = "bold"), # X轴标题加粗
    axis.title.y = element_text(size = 14, face = "bold"), # Y轴标题加粗
    axis.text.x = element_text(size = 12), # X轴刻度字号
    axis.text.y = element_text(size = 11), # Y轴基因名字号
    panel.grid.major.y = element_blank(), # 移除Y轴网格（避免遮挡基因名）
    panel.grid.minor.x = element_blank() # 移除X轴次要网格
  ) +
  geom_text(
    aes(label = round(FC, 2)), # 保留2位小数
    hjust = ifelse(top_genes$FC > 0, 0, 1), # 正值标签在右侧，负值在左侧
    vjust = 0.5,
    size = 3.5, # 标签字号
    color = "black"
  )
ggsave('01.scTenifoldKnk_NDUFA8.pdf', height = 8, width = 10)
ggsave('01.scTenifoldKnk_NDUFA8.png', height = 8, width = 10, dpi = 600)

df <- result$diffRegulation
df <- df[df$gene!='NDUFA8',]
df$log_pval <- -log10(df$p.value)
label_genes <- subset(df, abs(Z) > 2 & p.value < 0.05)
ggplot(df, aes(x=Z, y=log_pval)) +
  geom_point(alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="#E74C3C") +
  geom_vline(xintercept=c(2), linetype="dashed", color="steelblue") +
  geom_text_repel(data=label_genes, aes(label=gene),
                  size=3, max.overlaps=50) +
  labs(title="Z vs -log10(p-value) of NDUFA8",
       x="Z-score", y="-log10(p-value)") +
  theme_classic()
ggsave('02.scTenifoldKnk_NDUFA8_Z.pdf', height = 8, width = 10)
ggsave('02.scTenifoldKnk_NDUFA8_Z.png', height = 8, width = 10, dpi = 600)
write.csv(label_genes,'v_result_scTenifoldKnk.csv')


