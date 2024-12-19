setwd("D:/My Drive/lab files/cholsoon jang/pig flux/chow and HF data combined/final RNA_seq adjusted data/github metabolite transporter analyses/")
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(rstatix)
library(DESeq2)
library(reshape2)
library(apeglm)
library(biomaRt)
library(factoextra)
library(FactoMineR)

library(enrichR)
library(ggplot2)
library(org.Hs.eg.db)
library("pathview")
library(europepmc)
library(ggupset)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(DOSE)
library(ggplot2)
cnts_mat_full = read.csv('D:/My Drive/lab files/cholsoon jang/pig flux/chow and HF data combined/final RNA_seq adjusted data/label-fixed an tissue-segregated DE/final cnts mat prelabels.csv', check.names = F)
row.names(cnts_mat_full) = cnts_mat_full$gene_tissue
cnts_mat_full$gene_tissue=NULL

#DE genes
transprt_genes = read.csv('D:/My Drive/lab files/cholsoon jang/pig flux/chow and HF data combined/final RNA_seq adjusted data/BA genes.csv')
gene_list = unique(transprt_genes$Gene)

cnts_mat_full$gene_Symbol = gsub("\\_.*","",row.names(cnts_mat_full)) 
cnts_mat_full$tissue = gsub(".*_","", row.names(cnts_mat_full))
tCNTS = cnts_mat_full[cnts_mat_full$gene_Symbol %in% transprt_genes$Gene,]
tCNTS$gene_tissue= paste0(tCNTS$gene_Symbol, '_', tCNTS$tissue)
tt11 = tCNTS %>% dplyr::select(gene_tissue, gene_Symbol, tissue, everything())

write.csv(tt11, './raw_data/pan-organ transporter counts.csv.', row.names = F)






setwd("D:/My Drive/lab files/cholsoon jang/pig flux/chow and HF data combined/final RNA_seq adjusted data/github metabolite transporter analyses/")

cnts_mat_full = read.csv('./raw_data/pan-organ transporter counts.csv',check.names = F)
row.names(cnts_mat_full) = cnts_mat_full$gene_tissue
cnts_mat_full$gene_tissue=NULL
cnts_mat_full$gene_Symbol=NULL
cnts_mat_full$tissue=NULL
cnts_mat = dplyr::mutate_all(cnts_mat_full, function(x) as.integer(x))
colnames(cnts_mat) = gsub('HF', 'HFHS', colnames(cnts_mat))
colnames(cnts_mat) = gsub('Chow', 'Ct', colnames(cnts_mat))

coldata = as.data.frame(colnames(cnts_mat))
colnames(coldata) = 'SampleID'
coldata$condition = ifelse(grepl('HF', coldata$SampleID), 'HFHS', 'Ct')
table(coldata$condition)
coldata = coldata[!coldata$SampleID=='variable_tissue',]

#write.csv(cnts_mat, file = 'pan-organ normalized counts expression matrix label fixed.csv', row.names = T)

numerator_group = 'HFHS' 
denominator_group = 'Ct'

cnts_mat$gene_tissue =  row.names(cnts_mat)
cnts_mat$gene_tissue = gsub('LARGE', 'LARGE-INTESTINE', cnts_mat$gene_tissue)
cnts_mat$gene_tissue = gsub('SMALL', 'SMALL-INTESTINE', cnts_mat$gene_tissue)
cnts_mat$gene_Symbol = gsub("\\_.*","",cnts_mat$gene_tissue) 
cnts_mat$tissue = gsub(".*_","", cnts_mat$gene_tissue)

tCNTS = cnts_mat[cnts_mat$gene_Symbol %in% transprt_genes$Gene,]
row.names(tCNTS) = tCNTS$gene_tissue
tCNTS$gene_Symbol=NULL
tCNTS$gene_tissue=NULL

tissue1='HEART'
table(tCNTS$tissue)
run_PCA_tissue = function(tissue1){
  cnts_mat1 = tCNTS[grepl(tissue1, tCNTS$tissue),]
  cnts_mat1$tissue=NULL
  cnts_mat1 = cnts_mat1[ , colSums(is.na(cnts_mat1)) == 0]
  cnts_mat1$gene_Symbol = gsub("\\_.*","",row.names(cnts_mat1))
  cnts_mat1 = cnts_mat1[!duplicated(cnts_mat1$gene_Symbol),]
  row.names(cnts_mat1) = cnts_mat1$gene_Symbol
  
  cnts_mat1$gene_Symbol=NULL
  
  coldata1 = coldata[coldata$condition %in% numerator_group | coldata$condition %in% denominator_group,]
  cnts_mat1 = as.data.frame(cnts_mat1)
  
  cnts_mat1 = cnts_mat1[rowSums(cnts_mat1) > 4,]
  cnts_mat1 = cnts_mat1[,colnames(cnts_mat1) %in% coldata1$SampleID]
  cnts_mat1 = dplyr::mutate_all(cnts_mat1, function(x) as.integer(x))
  coldata1$genotype = factor(coldata1$condition, levels=c(numerator_group, denominator_group))
  iris.pca <- PCA(t(cnts_mat1), graph = FALSE)
  trait_cols = unique(coldata1$condition)
  names(trait_cols) = MetBrewer::met.brewer('Moreau', length(trait_cols))
  new_trait_t = colnames(cnts_mat1)
  names(new_trait_t) = coldata1$genotype[match(new_trait_t, coldata1$SampleID)]
  
  
  pdf(file = paste0(tissue1, ' - ',numerator_group,  ' over ', denominator_group, ' PCA.pdf'))
  g2 = fviz_pca_ind(iris.pca, col.ind = names(new_trait_t), palette = MetBrewer::met.brewer('Lakota', length(trait_cols)), addEllipses=T) + ggtitle(paste0(tissue1, ' PCA'))
  print(g2)
  dev.off()
  
  return(g2)
}
table(tCNTS$tissue)


gg1 = run_PCA_tissue('MUSCLE')
gg2 = run_PCA_tissue('HEART')
gg3 = run_PCA_tissue('LIVER')
gg4 = run_PCA_tissue('KIDNEY-CORTEX')
gg5 = run_PCA_tissue('KIDNEY-MED')
gg6 = run_PCA_tissue('LARGE-INTESTINE')
gg7 = run_PCA_tissue('SMALL-INTESTINE')
gg9 = run_PCA_tissue('LUNG')
gg10 = run_PCA_tissue('SKIN')
gg11= run_PCA_tissue('SPLEEN')
gg12 = run_PCA_tissue('WAT')

pdf(file = 'TRANSPORTER ONLY PCA rowsum4.pdf')
gridExtra::grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg9, gg10, gg11, gg12)
dev.off()




#Do PCA from all tissues
cnts_mat1 = tCNTS
cnts_mat1$tissue=NULL
cnts_mat1 = cnts_mat1[ , colSums(is.na(cnts_mat1)) == 0]


coldata1 = coldata[coldata$condition %in% numerator_group | coldata$condition %in% denominator_group,]
cnts_mat1 = as.data.frame(cnts_mat1)

cnts_mat1 = cnts_mat1[rowSums(cnts_mat1) > 4,]
cnts_mat1 = cnts_mat1[,colnames(cnts_mat1) %in% coldata1$SampleID]
cnts_mat1 = dplyr::mutate_all(cnts_mat1, function(x) as.integer(x))
coldata1$genotype = factor(coldata1$condition, levels=c(numerator_group, denominator_group))
iris.pca <- PCA(t(cnts_mat1), graph = FALSE)
trait_cols = unique(coldata1$condition)
names(trait_cols) = MetBrewer::met.brewer('Moreau', length(trait_cols))
new_trait_t = colnames(cnts_mat1)
names(new_trait_t) = coldata1$genotype[match(new_trait_t, coldata1$SampleID)]


pdf(file = paste0('FULL TISSUE- TRANSPORTERS ONLY',numerator_group,  ' over ', denominator_group, ' PCA.pdf'))
g2 = fviz_pca_ind(iris.pca, col.ind = names(new_trait_t), palette = MetBrewer::met.brewer('Lakota', length(trait_cols)), addEllipses=T) + ggtitle(paste0('Full tissue PCA TRANSPORTERS'))
print(g2)
dev.off()

#Note we run teh DE model here for transporters
run_tissueDE = function(tissue1){
  cnts_mat1 = tCNTS[grepl(tissue1, tCNTS$tissue),]
  cnts_mat1$tissue=NULL
  cnts_mat1 = cnts_mat1[ , colSums(is.na(cnts_mat1)) == 0]
  cnts_mat1$gene_Symbol = gsub("\\_.*","",row.names(cnts_mat1)    )
  cnts_mat1 = cnts_mat1[!duplicated(cnts_mat1$gene_Symbol),]
  row.names(cnts_mat1) = cnts_mat1$gene_Symbol
  
  cnts_mat1$gene_Symbol=NULL
  
  coldata1 = coldata[coldata$condition %in% numerator_group | coldata$condition %in% denominator_group,]
  cnts_mat1 = as.data.frame(cnts_mat1)
  
  cnts_mat1 = cnts_mat1[rowSums(cnts_mat1) > 4,]
  cnts_mat1 = cnts_mat1[,colnames(cnts_mat1) %in% coldata1$SampleID]
  cnts_mat1 = dplyr::mutate_all(cnts_mat1, function(x) as.integer(x))
  coldata1 = coldata1[coldata1$SampleID %in% colnames(cnts_mat1),]
  coldata1$genotype = factor(coldata1$condition, levels=c(numerator_group, denominator_group))
  
  dds <- DESeqDataSetFromMatrix(countData = cnts_mat1,
                                colData = coldata1,
                                design = ~ genotype)
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("genotype",  numerator_group, denominator_group))
  results_table =as.data.frame(res)
  res1 = results_table[order(results_table$pvalue, decreasing = F),]
  head(res1)
  test1 = cnts_mat1[grepl('PCK1', row.names(cnts_mat1)),]
  colnames(test1) = coldata1$genotype[match(colnames(test1), coldata1$SampleID)]
  res1$gene_symbol = row.names(res1)
  res1$tissue = paste0(tissue1)
  head(res1)
  head(res1)
  write.csv(res1, file =  paste0(tissue1, ' - DESeq2 Full results TRANPORTER ONLY MODEL ', numerator_group,  ' over ', denominator_group, '.csv'), row.names = F)
  return(res1)
}

gg1 = run_tissueDE('MUSCLE')
gg2 = run_tissueDE('HEART')
gg3 = run_tissueDE('LIVER')
gg4 = run_tissueDE('KIDNEY-CORTEX')
gg5 = run_tissueDE('KIDNEY-MED')
gg6 = run_tissueDE('LARGE-INTESTINE')
gg7 = run_tissueDE('SMALL-INTESTINE')
gg8 = run_tissueDE('MUSCLE')
gg9 = run_tissueDE('LUNG')
gg10 = run_tissueDE('SKIN')
gg11= run_tissueDE('SPLEEN')
gg12 = run_tissueDE('WAT')

full_de = as.data.frame(rbind(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, gg10, gg11, gg12))
full_de = full_de[order(full_de$padj, decreasing = F),]

ff2 = full_de %>%  dplyr::select(gene_symbol, tissue, baseMean, log2FoldChange, pvalue)
write.csv(ff2, file = 'TRANSPORTER ONLY DE results tissue-specific.csv', row.names = F)

#for ease of access we pull the DE results to plot
full_de = read.csv('./raw_data/TRANSPORTER ONLY DE results tissue-specific.csv')
head(full_de)

run_tissueDEVol = function(tissue1){
  #need to play around to assessing proper thresholds.  
  res1 = full_de
  res1 = res1[res1$tissue %in% tissue1,]
  label_key = res1$gene_symbol[1:6]
  res1 = na.omit(res1)
  
  res1$label2 = ifelse(res1$gene_symbol %in% label_key, paste0(res1$gene_symbol), '')
  table(res1$label2)
  res1$label_col1 = ifelse(res1$log2FoldChange>0, 'firebrick3', 'darkorchid')
  res1$label_col2 = ifelse(res1$pvalue<0.05, paste0(res1$label_col1), 'gray74')
  
  library(ggrepel)
  #Number of genes which will be labelled
  #Volcano plot
  #change labels here too
  pdf(file = paste0(tissue1, 'TRANPORTER Volcano plot ', numerator_group,  ' over ', denominator_group, '.pdf'))
  volc1 = ggplot(res1, aes(x=log2FoldChange, y=-log10(pvalue))) + theme_classic() +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue)), color=res1$label_col2) +
    geom_label_repel(aes(x=log2FoldChange, y=-log10(pvalue), label = res1$label2), color = res1$label_col2, size = 2, label.size=NA, box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf, alpha = .6, segment.color = 'grey50')  +   ggtitle(paste0(tissue1)) + ylab('-log10(pvalue)')
  print(volc1)
  dev.off()
  return(volc1)
}
gg1 = run_tissueDEVol('MUSCLE')
gg2 = run_tissueDEVol('HEART')
gg3 = run_tissueDEVol('LIVER')
gg4 = run_tissueDEVol('KIDNEY-CORTEX')
gg5 = run_tissueDEVol('KIDNEY-MED')
gg6 = run_tissueDEVol('LARGE-INTESTINE')
gg7 = run_tissueDEVol('SMALL-INTESTINE')
gg9 = run_tissueDEVol('LUNG')
gg10 = run_tissueDEVol('SKIN')
gg11= run_tissueDEVol('SPLEEN')
gg12 = run_tissueDEVol('WAT')

pdf(file = 'TRANSPORTER ONLY tissue-specific Volcano rowsum4.pdf')
gridExtra::grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7,  gg9, gg10, gg11, gg12)
dev.off()




