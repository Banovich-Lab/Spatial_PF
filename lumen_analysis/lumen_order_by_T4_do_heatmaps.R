suppressPackageStartupMessages({
  library(scater)
  library(ggplot2)
  library(mclust)
  library(slingshot)
  library(tradeSeq)
  library(grDevices)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(dplyr)
  library(Seurat)
  library(data.table)
})

### Gene split by gene anchors 
gene_expr_ordered <- readr::read_csv("analysis/revision_2024/lumen_order_t4/gene_order_by_peaktime_by_change_time_mat.csv")
gene_names <- (gene_expr_ordered$...1)
gene_expr_ordered_m <- gene_expr_ordered[,2:ncol(gene_expr_ordered)]
gene_expr_ordered_m <- as.matrix(gene_expr_ordered_m)
rownames(gene_expr_ordered_m) <- gene_names

gene_expr_ordered_row_annot <-  readr::read_csv("analysis/revision_2024/lumen_order_t4/gene_order_by_peaktime_by_change_time_mat_rowAnnot.csv")
gene_expr_ordered_row_annot <- data.frame(gene_expr_ordered_row_annot)
rownames(gene_expr_ordered_row_annot) <- gene_expr_ordered_row_annot$genes
sum(rownames(gene_expr_ordered_row_annot) == rownames(gene_expr_ordered_m))

Heatmap(gene_expr_ordered_m,
        cluster_rows = F,cluster_columns = F,
        row_split = factor(gene_expr_ordered_row_annot$split_gene_anchor,
                           levels=c("bimodal","Normal",
                                    "Early","Late")),
        show_column_names = FALSE)



### Gene split by hcluster chunks 
gene_expr_ordered <- readr::read_csv("analysis/revision_2024/lumen_order_t4/hcluster_chunk_ordered_change_time_ordered_mat.csv")
gene_names <- (gene_expr_ordered$...1)
gene_expr_ordered_m <- gene_expr_ordered[,2:ncol(gene_expr_ordered)]
gene_expr_ordered_m <- as.matrix(gene_expr_ordered_m)
rownames(gene_expr_ordered_m) <- gene_names
dim(gene_expr_ordered_m)

gene_expr_hchunk_ordered_row_annot <-  readr::read_csv("analysis/revision_2024/lumen_order_t4/hcluster_chunk_ordered_change_time_ordered_mat_rowAnnot.csv")
gene_expr_hchunk_ordered_row_annot <- data.frame(gene_expr_hchunk_ordered_row_annot)
rownames(gene_expr_hchunk_ordered_row_annot) <- gene_expr_hchunk_ordered_row_annot$genes
dim(gene_expr_hchunk_ordered_row_annot)

Heatmap(gene_expr_ordered_m,cluster_rows = F,cluster_columns = F,
        row_split = factor(gene_expr_hchunk_ordered_row_annot$gene_split_label,
                           levels=c("bimodal","Normal",
                                    "Early","Late")),
        show_column_names = FALSE)

## per celltype gene testing results:

sig_gene_per_ct <- read.csv("analysis/revision_2024/lumen_order_t4/sig_genes_per_celltype.csv")


sig_gene_per_ct %>% ggplot()+geom_bar(mapping = aes(x = final_lineage,
                                                    fill=split_gene_anchor))+
  facet_wrap(.~split_gene_anchor)+
  scale_fill_manual(values=DiscretePalette(n=5,palette = "glasbey"))+
  theme(axis.text.x = element_text(angle=90))+
  sig_gene_per_ct %>% ggplot()+geom_bar(mapping = aes(x = split_gene_anchor,
                                                      fill=final_lineage),position = "dodge")+
  scale_fill_manual(values=DiscretePalette(n=5))+
  
  
  sig_gene_per_ct %>% ggplot()+geom_bar(mapping = aes(x = final_lineage,
                                                      fill=split_hcluster_chunk))+
  facet_wrap(.~split_hcluster_chunk)+
  scale_fill_manual(values=DiscretePalette(n=5,palette = "glasbey"))+
  theme(axis.text.x = element_text(angle=90))+
  sig_gene_per_ct %>% ggplot()+geom_bar(mapping = aes(x = split_hcluster_chunk,
                                                      fill= final_lineage),position = "dodge")+
  scale_fill_manual(values=DiscretePalette(n=5))
