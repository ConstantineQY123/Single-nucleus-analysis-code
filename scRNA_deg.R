
rm(list = ls())
setwd("/data/nas3/zhangnan/projects/small_139_BJTC-630-4/")
if(!dir.exists('./scRNA_deg')){dir.create('./scRNA_deg')}
setwd('./scRNA_deg')

# /data/nas1/xuguangxin/project/39_BJTC-630-4/

DE_res <- read.csv(
 file = "/data/nas1/xuguangxin/project/39_BJTC-630-4/01_scRNA/07_subtype/FADS3_cluster_DE.csv",
 row.names = 1,
 check.names = FALSE,
 stringsAsFactors = FALSE
)

target_genes <- c("COX6C", "NDUFA8", "NDUFB1")

genes_present <- intersect(target_genes, rownames(DE_res))
genes_missing <- setdiff(target_genes, rownames(DE_res))

DE_target <- DE_res[genes_present, , drop = FALSE]
DE_target$gene <- rownames(DE_target)
print(DE_target)
write.csv(
 DE_target,
 file = "keygene_FADS3_cluster_DEG.csv",
 row.names = FALSE,
 quote = FALSE
)
