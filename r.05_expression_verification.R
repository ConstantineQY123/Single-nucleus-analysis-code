rm(list = ls())
setwd('/data/nas1/xuguangxin/project/39_BJTC-630-4/')
if(!dir.exists('./05_expression_verification')){
  dir.create('./05_expression_verification')
}
setwd('./05_expression_verification/')

# 加载包
library(GEOquery)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(Biobase)
library(lance)
library(Seurat)

hub_gene <- read.csv(file = "../04_PPI/hub_gene.csv", row.names = 1)

# 验证集-----------------------------------------------------
verification_data <- read.csv(paste0('../00_rawdata/GSE36961/dat.GSE36961.csv'), row.names = 1, check.names = F)
max(verification_data)
# verification_data <- log2(verification_data + 0.1)
verification_group <- read.csv(paste0('../00_rawdata/GSE36961/group.GSE36961.csv'))
table(verification_group$group)

temp1 <- verification_data[rownames(verification_data) %in% hub_gene$symbol, ]
verification_disease_sample <- verification_group[verification_group$group == 'disease', ]
verification_disease_data <- temp1[, colnames(temp1) %in% verification_disease_sample$sample]
verification_disease_mean <- apply(verification_disease_data, 1, mean) %>% as.data.frame()
colnames(verification_disease_mean) <- 'disease_mean'
verification_control_sample <- verification_group[verification_group$group == 'control', ]
verification_control_data <- temp1[, colnames(temp1) %in% verification_control_sample$sample]
verification_control_mean <- apply(verification_control_data, 1, mean) %>% as.data.frame()
colnames(verification_control_mean) <- 'control_mean'
verification_mean <- cbind(verification_disease_mean, verification_control_mean)
verification_mean$compare_res <- ifelse(verification_mean$disease_mean > verification_mean$control_mean, 1, 0)
verification_mean$symbol <- rownames(verification_mean)

colnames(temp1) <- gsub('X', '', colnames(temp1))
# temp2 <- apply(temp1, 1, mean) %>% as.data.frame()
# temp2$. == -3.3219280
temp1$symbol <- rownames(temp1)
# temp1 <- temp1[-c(23,24,48), ]
# temp1 <- temp1[20:30, ]
temp1 <- pivot_longer(
  temp1,
  cols = -symbol, 
  names_to = 'sample', 
  values_to = 'expr'
)

verification_hub_exp <- temp1

verification_hub_exp <- merge(verification_hub_exp, verification_group, by = 'sample')
verification_hub_exp$group <- as.factor(verification_hub_exp$group)

verification_stat.test <- verification_hub_exp%>%
  group_by(symbol)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'fdr')

verification_stat.test$p <- ifelse(verification_stat.test$p<0.001,"***",
                                   ifelse(verification_stat.test$p<0.01,"**",
                                          ifelse(verification_stat.test$p<0.05,"*",'ns')))

verification_finally_res <- merge(verification_stat.test, verification_mean, by = 'symbol')
write.csv(verification_finally_res, 'verification_finally_res.csv')
verification_finally_res <- dplyr::select(verification_finally_res, symbol, p, compare_res)
verification_finally_res$p <- ifelse(verification_finally_res$p < 0.05, 1, 0)

# 训练集-----------------------------------------------------
load('../01_scRNA/07_subtype/TSNE_Cardiomyocyte_marker.Rdata')
table(sce.car$celltype)
table(sce.car$seurat_clusters)
scRNA <- subset(sce.car, seurat_clusters == 'FADS3_cluster')

train_data <- as.data.frame(scRNA@assays$RNA$counts) %>% apply(2,function(x){x/sum(x) * 10000})
train_data <- train_data %>% as.data.frame()
# train_data <- log2(train_data+1) %>% ScaleData() %>% as.data.frame()

train_group <- scRNA@meta.data
train_group$sample <- rownames(train_group)
train_group <- dplyr::select(train_group, sample, group)
table(train_group$group)
# ex_data$symbol <- rownames(ex_data)

temp1 <- train_data[rownames(train_data) %in% hub_gene$symbol, ]
train_disease_sample <- train_group[train_group$group == 'disease', ]
train_disease_data <- temp1[, colnames(temp1) %in% train_disease_sample$sample]
train_disease_mean <- apply(train_disease_data, 1, mean) %>% as.data.frame()
colnames(train_disease_mean) <- 'disease_mean'
train_control_sample <- train_group[train_group$group == 'control', ]
train_control_data <- temp1[, colnames(temp1) %in% train_control_sample$sample]
train_control_mean <- apply(train_control_data, 1, mean) %>% as.data.frame()
colnames(train_control_mean) <- 'control_mean'
train_mean <- cbind(train_disease_mean, train_control_mean)
train_mean$compare_res <- ifelse(train_mean$disease_mean > train_mean$control_mean, 1, 0)
train_mean$symbol <- rownames(train_mean)

colnames(temp1) <- gsub('X', '', colnames(temp1))
# temp2 <- apply(temp1, 1, mean) %>% as.data.frame()
# temp2$. == -3.3219280
temp1$symbol <- rownames(temp1)
# temp1 <- temp1[-c(23,24,48), ]
# temp1 <- temp1[20:30, ]
temp1 <- pivot_longer(
  temp1,
  cols = -symbol, 
  names_to = 'sample', 
  values_to = 'expr'
)

train_hub_exp <- temp1

train_hub_exp <- merge(train_hub_exp, train_group, by = 'sample')
train_hub_exp$group <- as.factor(train_hub_exp$group)

train_stat.test <- train_hub_exp%>%
  group_by(symbol)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'fdr')
i<-1
for (i in c(1 : nrow(train_stat.test))) {
  if (train_stat.test$p[i] < 0.05) {
    next
  }else{
    DEG_res <- read.csv(paste0('../01_scRNA/05_subtype/cluter_0_DE_sig.csv'), row.names = 1)
    DEG_res$symbol <- rownames(DEG_res)
    DEG_res <- DEG_res[DEG_res$symbol %in% train_stat.test$symbol, ]
    train_stat.test$p[i] <- DEG_res[DEG_res$symbol == train_stat.test$symbol[i], ]$p_val
  }
}

train_stat.test$p <- ifelse(train_stat.test$p<0.001,"***",
                            ifelse(train_stat.test$p<0.01,"**",
                                   ifelse(train_stat.test$p<0.05,"*",'ns')))

train_finally_res <- merge(train_stat.test, train_mean, by = 'symbol')
write.csv(train_finally_res, 'train_finally_res.csv')
train_finally_res <- dplyr::select(train_finally_res, symbol, p, compare_res)
train_finally_res$p <- ifelse(train_finally_res$p < 0.05, 1, 0)

finally_res <- merge(verification_finally_res, train_finally_res, by = 'symbol')
i<-1
res <- data.frame()
for (i in c(1 : length(finally_res$symbol))) {
  if (finally_res$p.x[i] == finally_res$p.y[i] & finally_res$compare_res.x[i] == finally_res$compare_res.y[i]) {
    temp2 <- data.frame(symbol = finally_res$symbol[i],
                        res = 'TRUE')
    res <- rbind(res, temp2)
  }else{
    next
  }
}

hub_gene <- res$symbol %>% as.data.frame()
colnames(hub_gene) <- 'symbol'
write.csv(hub_gene, './hub_gene.csv')

train_stat.test <- train_stat.test %>% 
  add_xy_position(x='group', dodge = 0.5, step.increase = 0.1)
train_stat.test$y.position <- c(7, 7.5, 5.5, 11, 2.5, 2.5, 2.5, 4.5, 3.5, 2.5)

train_hub_exp$expr <- abs(train_hub_exp$expr)
train_va <- ggplot(train_hub_exp, aes(x = group, y = expr, ))+
  geom_bar(aes(fill = group), color="black", stat="summary", fun=mean, position="dodge", size=0.5)+
  scale_fill_manual(values= c("#00BCF2","gold"))+
  stat_summary(fun = mean,
               geom = "pointrange",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)))+
  labs(title="Expression in train", x="", y = "Expression",size=20) +
  stat_pvalue_manual(train_stat.test,
                     tip.length = 0.01,
                     size = 3.2,
                     family = "Times",
                     label = "p",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=20),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_blank(),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~symbol,scales = "free",nrow = 1)
train_va
ggsave(filename = '01.train_res.pdf', train_va, w=12, h=8)
ggsave(filename = '01.train_res.png', train_va, w=12, h=8, dpi = 600)


verification_stat.test <- verification_stat.test %>% 
  add_xy_position(x='group', dodge = 0.5, step.increase = 0.1)
verification_stat.test$y.position <- c(2, 4.5, 1.5)

verification_hub_exp$expr <- abs(verification_hub_exp$expr)
External_va <- ggplot(verification_hub_exp, aes(x = group, y = expr, ))+
  geom_bar(aes(fill = group), color="black", stat="summary", fun=mean, position="dodge", size=0.5)+
  scale_fill_manual(values= c("#00BCF2","gold"))+
  stat_summary(fun = mean,
               geom = "pointrange",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)))+
  labs(title="Expression in verification", x="", y = "Expression",size=20) +
  stat_pvalue_manual(verification_stat.test,
                     tip.length = 0.01,
                     size = 3.2,
                     family = "Times",
                     label = "p",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=20),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_blank(),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~symbol,scales = "free",nrow = 1)
External_va
ggsave(filename = '02.verification_res.pdf', External_va, w=12, h=8)
ggsave(filename = '02.verification_res.png', External_va, w=12, h=8, dpi = 600)

