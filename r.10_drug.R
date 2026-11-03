## 2.1 DSigDBS数据库-------------
rm(list = ls())
setwd("/data/nas1/xuguangxin/project/39_BJTC-630-4/")
if (! dir.exists("./10_drug")){
  dir.create("./10_drug")
}
setwd("./10_drug")

# install.packages('enrichR')
library(enrichR)
dbs <- listEnrichrDbs()

## DsigDB
dbs$libraryName
dbs <- c("DSigDB")
hub_gene <- read.csv('../05_expression_verification/hub_gene.csv',row.names = 1)
colnames(hub_gene) <- 'symbol'
symbol <- hub_gene$symbol

enrichr <- enrichr(symbol, dbs)
result <- data.frame(enrichr$DSigDB)
result$Term <- str_split_fixed(result$Term, pattern = " ",3)[,1]

write.csv(result, file = 'drug.result.csv',row.names = F,quote = F)

hub_gene$symbol
finally_res <- data.frame()
i <- 3
for (i in c(1 : length(hub_gene$symbol))) {
  temp1 <- result[grepl(hub_gene$symbol[i], result$Genes), ]
  temp1 <- temp1[order(temp1$Combined.Score, decreasing = T), ]
  temp1 <- temp1[1 : 10, ]
  finally_res <- rbind(finally_res, temp1)
}
finally_res
finally_res <- finally_res[!duplicated(finally_res), ]
write.csv(finally_res, file = 'cbdock_res.csv',row.names = F,quote = F)

cytoscape_data <- dplyr::select(result, Term, Genes)
cytoscape_data <- separate(cytoscape_data, Genes, into = c('gene1', 'gene2', 'gene3', 'gene4', 'gene5'), sep = ';')

drug_gene <- data.frame(drug = '1',
                        symbol = '1')
for (i in c(1 : (ncol(cytoscape_data)-1))) {
  temp2 <- cytoscape_data[, c(1, i+1)]
  colnames(temp2) <- c('drug', 'symbol')
  drug_gene <- rbind(drug_gene, temp2)
}
drug_gene <- na.omit(drug_gene)
drug_gene <- drug_gene[-1, ]
drug_gene <- distinct(drug_gene, .keep_all = T)
write.csv(drug_gene, file = 'cytoscape_res.csv')
