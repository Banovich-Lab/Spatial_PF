## 4 clusters
suppressPackageStartupMessages({
  library(ggplot2)
  library(grDevices)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(dplyr)
  library(Seurat)
  library(data.table)
})

library(circlize)
pseudotime_col_fun = colorRamp2(c(0, 500, 1500), c("white", "lightblue", "darkblue")) 
cell_col_fun = colorRamp2(c(-5, 0, 5), c("darkblue", "white", "firebrick")) 
pathology_col_fun = colorRamp2(c(0, 50, 100), c("white", "lightblue", "darkblue")) 
sample_type_label <- lumen_celltype_sce$sample_type
sample_type_label_color <- DiscretePalette(n=4)
names(sample_type_label_color) <- unique(sample_type_label)




top_annot_df <- read.csv("./analysis/revision_2024/lumen_order_t4/offset_total_trans/time_ordered_lumen_annot.csv")
sc_cluster_levels <- c(4,3,1,2)
gene_by_lumen_tma_corrected <- read.csv(
          "./analysis/revision_2024/lumen_order_t4/offset_total_trans/specc_cluster_offset_total_tx_tma_corrected_change_time_in_chunk_4c_mat.csv")
rownames(gene_by_lumen_tma_corrected) <- gene_by_lumen_tma_corrected$X
gene_by_lumen_tma_corrected <- gene_by_lumen_tma_corrected[,2:ncol(gene_by_lumen_tma_corrected)]
gene_by_lumen_tma_corrected <- as.matrix(gene_by_lumen_tma_corrected)
line_plot_df <-  read.csv(
  "./analysis/revision_2024/lumen_order_t4/offset_total_trans/specc_cluster_offset_total_tx_tma_corrected_change_time_in_chunk_4c_line_plot_df.csv")
celltype_to_lineage <- read.csv(
  "../../../Dataset/LFST_2022/REVISION_new_samples_052924/062724_late_IPF_Seurat_and_heatmap/ct_lineage_map.csv")
celltype_to_lineage$cell_object_meta.final_CT <- make.names(celltype_to_lineage$cell_object_meta.final_CT)
gene_annot <- unique(line_plot_df[,c("gene","sc_clusters")])
sig_genes_per_ct <- readr::read_csv("output/savedRDS/revisions2024/perCT/all_perCT_sig_terms_s15.csv")

rownames(gene_annot) <- gene_annot$gene
sig_genes_per_ct %>% left_join(gene_annot) %>% left_join(celltype_to_lineage,
                                                         by = c("celltype"="cell_object_meta.final_CT"))

Heatmap(gene_by_lumen_tma_corrected,
        show_row_names = T,
        cluster_rows = F,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_split = factor(gene_annot[rownames(gene_by_lumen_tma_corrected),
                                      "sc_clusters"],
                                      levels=sc_cluster_levels),
        row_gap = unit(0.5, "mm"),
        row_names_gp = gpar(fontsize=4),
        top_annotation =
          HeatmapAnnotation(pseudotime = top_annot_df$pseudotime_value,
                            sample_type = top_annot_df$sample_type_label,
                            percent_pathology = top_annot_df$sample_perc_pathology,
                            col  = list(pseudotime = pseudotime_col_fun,
                                        sample_type = sample_type_label_color,
                                        percent_pathology = pathology_col_fun)),
        row_title_gp  = gpar(fontsize=7),
        name = "Gene Expr\n z-scores",
        heatmap_legend_param = 
          list(    legend_height = unit(2, "cm"),
                   legend_gp = gpar(fontsize=7),
                   title_gp = gpar(fontsize = 7, 
                                   fontface = "bold")))


line_plot_df %>% 
  ggplot(aes(x = time_order,y = value))+
  geom_point(size=0.1,alpha=0.3,color="grey")+
  geom_smooth(mapping = aes(x = time_order,y = value),
              se=TRUE,colour="red",method="gam")+
  facet_wrap(.~sc_clusters)
length(unique(line_plot_df$gene))
####

#### 5 clusters

# write.csv(data.frame(pseudotime_value,sample_type_label,sample_perc_pathology,
#                      lumen_id = colnames(gene_by_lumen_tma_corrected)) ,
#           "./analysis/revision_2024/lumen_order_t4/offset_total_trans/time_ordered_lumen_annot.csv")

sc_cluster_levels <- c(1,5,3,2,4)
 
gene_by_lumen_tma_corrected <- read.csv(
  "./analysis/revision_2024/lumen_order_t4/offset_total_trans/specc_cluster_offset_total_tx_tma_corrected_change_time_in_chunk_5c_mat.csv")
rownames(gene_by_lumen_tma_corrected) <- gene_by_lumen_tma_corrected$X
gene_by_lumen_tma_corrected <- gene_by_lumen_tma_corrected[,2:ncol(gene_by_lumen_tma_corrected)]
gene_by_lumen_tma_corrected <- as.matrix(gene_by_lumen_tma_corrected)

line_plot_df <-  read.csv(
  "./analysis/revision_2024/lumen_order_t4/offset_total_trans/specc_cluster_offset_total_tx_tma_corrected_change_time_in_chunk_5c_line_plot_df.csv")
gene_annot <- unique(line_plot_df[,c("gene",
                                     "sc_clusters")])
rownames(gene_annot) <- gene_annot$gene

Heatmap(gene_by_lumen_tma_corrected,
        show_row_names = T,
        cluster_rows = F,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_split = factor(gene_annot[rownames(gene_by_lumen_tma_corrected),
                                      "sc_clusters"],
                           levels=sc_cluster_levels),
        row_gap = unit(0.5, "mm"),
        row_names_gp = gpar(fontsize=4),
        top_annotation =
          HeatmapAnnotation(pseudotime = top_annot_df$pseudotime_value,
                            sample_type = top_annot_df$sample_type_label,
                            percent_pathology = top_annot_df$sample_perc_pathology,
                            col  = list(pseudotime = pseudotime_col_fun,
                                        sample_type = sample_type_label_color,
                                        percent_pathology = pathology_col_fun)),
        row_title_gp  = gpar(fontsize=7),
        name = "Gene Expr\n z-scores",
        heatmap_legend_param = 
          list(    legend_height = unit(2, "cm"),
                   legend_gp = gpar(fontsize=7),
                   title_gp = gpar(fontsize = 7, 
                                   fontface = "bold")))

line_plot_df %>% 
  ggplot(aes(x = time_order,y = value))+
  geom_point(size=0.1,alpha=0.3,color="grey")+
  geom_smooth(mapping = aes(x = time_order,y = value),
              se=TRUE,colour="red",method="gam")+
  facet_wrap(.~sc_clusters)

length(unique(line_plot_df$gene))


