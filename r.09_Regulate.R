rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/")
if (! dir.exists("./09_Regulate")){
  dir.create("./09_Regulate")
}
setwd("./09_Regulate")

library(multiMiR)
hubgene <- read.csv("../05_expression_verification/hub_gene.csv", row.names = 1)
colnames(hubgene) <- 'symbol'
hubgene$symbol
gene2mir <- get_multimir(org     = 'hsa',
                         target  = hubgene$symbol,
                         table   = 'predicted',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff= 500000)
table(gene2mir@data$database)

mirdb <- gene2mir@data[gene2mir@data$database == "mirdb", ]
colnames(mirdb)
mirdb <- dplyr::select(mirdb, mature_mirna_acc, mature_mirna_id, target_symbol)
colnames(mirdb) <- c('miRNA_ID', 'miRNA', 'symbol')
write.csv(mirdb, "mirdb.csv")

miranda <- gene2mir@data[gene2mir@data$database == "miranda", ]
colnames(miranda)
miranda <- dplyr::select(miranda, mature_mirna_acc, mature_mirna_id, target_symbol)
colnames(miranda) <- c('miRNA_ID', 'miRNA', 'symbol')
write.csv(miranda, "miranda.csv")

microcosm <- gene2mir@data[gene2mir@data$database == "microcosm", ]
colnames(microcosm)
microcosm <- dplyr::select(microcosm, mature_mirna_acc, mature_mirna_id, target_symbol)
colnames(microcosm) <- c('miRNA_ID', 'miRNA', 'symbol')
write.csv(microcosm, "microcosm.csv")

targetscan <- gene2mir@data[gene2mir@data$database == "targetscan", ]
colnames(targetscan)
targetscan <- dplyr::select(targetscan, mature_mirna_acc, mature_mirna_id, target_symbol)
colnames(targetscan) <- c('miRNA_ID', 'miRNA', 'symbol')
write.csv(targetscan, "targetscan.csv")

diana_microt <- gene2mir@data[gene2mir@data$database == "diana_microt", ]
colnames(diana_microt)
diana_microt <- dplyr::select(diana_microt, mature_mirna_acc, mature_mirna_id, target_symbol)
colnames(diana_microt) <- c('miRNA_ID', 'miRNA', 'symbol')
write.csv(diana_microt, "diana_microt.csv")

inter_miRNA <- intersect(microcosm$miRNA, miranda$miRNA) %>% as.data.frame()
inter_miRNA <- rbind(microcosm[microcosm$miRNA %in% inter_miRNA$., ], miranda[miranda$miRNA %in% inter_miRNA$., ])
inter_miRNA <- inter_miRNA[!duplicated(inter_miRNA), ]
unique(inter_miRNA$symbol)
write.csv(inter_miRNA, "inter_miRNA.csv")
