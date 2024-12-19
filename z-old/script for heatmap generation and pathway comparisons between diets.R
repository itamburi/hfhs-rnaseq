library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(rstatix)
setwd("G:/My Drive/lab files/cholsoon jang/pig flux/chow and HF data combined/final RNA_seq adjusted data")
gene_set = c('CYP7A1', 'CYP27A1', 'CYP7B1', 'CYP8B1', 'HSD3B7', 'ABCD3', 'AMACR', 'ACOX2', 'HSD17B4', 'SCP2', 'ACOT8')
full_melted_counts = read.csv('batch-adjusted counts based only set all pigs.csv')
ff1 = full_melted_counts[full_melted_counts$Variable_ID %in% gene_set,]
ff2 = dcast(ff1, Variable_ID~ diet_tissue, fun.aggregate = mean, value.var = 'value')
row.names(ff2) = ff2$Variable_ID
ff2$Variable_ID=NULL
pdf(file = 'scaled heatmap only one sva adj.pdf')
pheatmap(ff2, cluster_cols = F,  color = colorRampPalette(c("black", "orange", "white"))(100))
dev.off()


stat.test <- ff1 %>%
  group_by(tissue) %>%
  t_test(value ~ diet)

stat.test <- stat.test %>%
  add_xy_position(x = "tissue", dodge = 0.8)
pdf(file = 'comparison of select BA and peroxisome genes.pdf')
ggboxplot(ff1, x = "tissue", y = "value", 
          color = "diet", palette = c("darkorange", "darkorchid")
) +  stat_pvalue_manual(stat.test,  label = "p", tip.length = 0
) + geom_point(aes(x=tissue, y=value, fill = diet, col=diet), position = position_jitterdodge(0.1, dodge.width = 0.8)) 
dev.off()


pathwayset = read.delim('uniprot-human-genes and goterms mapping.tab')
head(pathwayset)
ba_paths = pathwayset[grepl('bile acid', pathwayset$Gene.ontology..biological.process.),]

ff3 = full_melted_counts[full_melted_counts$Variable_ID %in% ba_paths$Gene.names...primary..,]

stat.test <- ff3 %>%
  group_by(tissue) %>%
  t_test(value ~ diet)

stat.test <- stat.test %>%
  add_xy_position(x = "tissue", dodge = 0.8)
pdf(file = 'comparison of gene onolotgy annotations BA genes.pdf')
ggboxplot(ff3, x = "tissue", y = "value", 
          color = "diet", palette = c("darkorange", "darkorchid")
) +  stat_pvalue_manual(stat.test,  label = "p", tip.length = 0
) + geom_point(aes(x=tissue, y=value, fill = diet, col=diet), position = position_jitterdodge(0.1, dodge.width = 0.8)) 
dev.off()

#The same for peroxisome
per_paths = pathwayset[grepl('peroxisome', pathwayset$Gene.ontology..biological.process.),]

ff3 = full_melted_counts[full_melted_counts$Variable_ID %in% per_paths$Gene.names...primary..,]

stat.test <- ff3 %>%
  group_by(tissue) %>%
  t_test(value ~ diet)

stat.test <- stat.test %>%
  add_xy_position(x = "tissue", dodge = 0.8)
pdf(file = 'comparison of gene onolotgy annotations peroxisome genes.pdf')
ggboxplot(ff3, x = "tissue", y = "value", 
          color = "diet", palette = c("darkorange", "darkorchid")
) +  stat_pvalue_manual(stat.test,  label = "p", tip.length = 0
) + geom_point(aes(x=tissue, y=value, fill = diet, col=diet), position = position_jitterdodge(0.1, dodge.width = 0.8)) 
dev.off()
