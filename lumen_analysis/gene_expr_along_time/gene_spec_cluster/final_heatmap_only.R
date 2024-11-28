
library(ggplot2)
library(grDevices)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(Seurat)
library(dplyr)


## Final heatmap plots (with bimodal genes), and add per cell type results

gene_by_lumen_tma_corrected <-
  read.csv("gene_by_lumen_tma_corrected_282genes_mat.csv")
rownames(gene_by_lumen_tma_corrected) <- gene_by_lumen_tma_corrected$X
gene_by_lumen_tma_corrected <- gene_by_lumen_tma_corrected[,2:ncol(gene_by_lumen_tma_corrected)]
gene_by_lumen_tma_corrected <- as.matrix(gene_by_lumen_tma_corrected)
dim(gene_by_lumen_tma_corrected)

gene_annot <- read.csv(
  "./gene_specc_c4_257genes.csv")
dim(gene_annot)
celltype_to_lineage <- read.csv(
  "ct_lineage_map.csv")

celltype_to_lineage$cell_object_meta.final_CT <- make.names(celltype_to_lineage$cell_object_meta.final_CT)

sig_genes_per_ct <- readr::read_csv("./all_perCT_sig_terms_s15_filterby50p.csv")
rownames(gene_annot) <- gene_annot$gene

sig_genes_per_ct_res <- sig_genes_per_ct %>%
  left_join(gene_annot) %>% left_join(celltype_to_lineage,
                                      by = c("celltype"="cell_object_meta.final_CT")) %>%
  mutate(gene_cat = case_when(sc_clusters == "4" ~ "Normal",
                              sc_clusters == "3" ~ "Early",
                              sc_clusters == "1" ~ "Middle",
                              sc_clusters == "2" ~ "Late" )) %>%
  mutate(gene_cat = factor(gene_cat,levels = c("Normal","Early","Middle","Late")))


p_total_tested <- sig_genes_per_ct_res %>%
  ggplot()+geom_bar(mapping = aes(x = gene_cat, 
                                  fill = cell_object_meta.final_lineage),
                    position = "dodge")+
  scale_fill_manual(values = DiscretePalette(n=4))+
  theme(axis.text.x = element_text(angle=90)) + ylab("Tested gene_ct pairs")

p_total_sig <- sig_genes_per_ct_res %>% filter(padj <=0.05) %>%
  ggplot()+geom_bar(mapping = aes(x = gene_cat, 
                                  fill = cell_object_meta.final_lineage),
                    position = "dodge")+
  scale_fill_manual(values = DiscretePalette(n=4))+
  theme(axis.text.x = element_text(angle=90)) + ylab("Sig. gene_ct pairs")

p_total_sig <- sig_genes_per_ct_res %>% filter(padj <=0.05) %>%
  group_by(gene_cat,cell_object_meta.final_lineage) %>% 
  summarise(sig_cat_lineage = n()) %>%
  ggplot(mapping = aes(x = gene_cat, 
                       y = sig_cat_lineage,
                       label = sig_cat_lineage,
                       fill = cell_object_meta.final_lineage))+
  geom_bar(stat="identity",
                    position = "dodge")+
  geom_text(position = position_dodge(width = 1,))+
  scale_fill_manual(values = DiscretePalette(n=4))+
  theme(axis.text.x = element_text(angle=90)) + ylab("Sig. gene_ct pairs")


p_total_sig_prop <- sig_genes_per_ct_res %>% filter(padj <=0.05) %>%
  group_by(gene_cat,cell_object_meta.final_lineage) %>% 
  mutate(sig_cat_lineage = n()) %>%
  group_by(gene_cat) %>% mutate(total_sig_cat = n()) %>%
  group_by(gene_cat,cell_object_meta.final_lineage) %>%
  summarise(sig_prop = unique(sig_cat_lineage)/unique(total_sig_cat)) %>%
  ggplot(aes(x = gene_cat, y = sig_prop, label = round(sig_prop,2),
             fill = cell_object_meta.final_lineage))+geom_bar(
                    position = "fill",stat="identity")+
  geom_text(position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = DiscretePalette(n=4))+
  theme(axis.text.x = element_text(angle=90)) + 
  ylab("Sig. Prop in total sig. gene_ct pairs")

p_total_sig + p_total_sig_prop


p_lineage_prop <- sig_genes_per_ct_res %>% filter(padj<=0.05) %>%
  group_by(cell_object_meta.final_lineage) %>%
  mutate(lineage_overall_sig = n()) %>% 
  group_by(cell_object_meta.final_lineage,gene_cat) %>%
  summarise(lineage_gene_cat_sig = n(),
            lineage_overall_sig = unique(lineage_overall_sig),
            lineage_gene_cat_prop =
              unique(lineage_gene_cat_sig)/unique(lineage_overall_sig)) %>% 
  ggplot(aes(x = cell_object_meta.final_lineage, y = lineage_gene_cat_prop,
             fill = gene_cat,
             label = round(lineage_gene_cat_prop,2)))+geom_bar(
                    position = "stack",stat='identity')+
  geom_text(position = position_stack(vjust=0.2))+
  scale_fill_manual(values = DiscretePalette(n=4))+
  theme(axis.text.x = element_text(angle=90))+ggtitle("padj <=0.05 ")


p_lineage_sig_count <- sig_genes_per_ct_res %>% filter(padj<=0.05) %>%
  group_by(cell_object_meta.final_lineage) %>%
  mutate(lineage_overall_sig = n()) %>% 
  group_by(cell_object_meta.final_lineage,gene_cat) %>%
  summarise(lineage_gene_cat_sig = n(),
            lineage_overall_sig = unique(lineage_overall_sig),
            lineage_gene_cat_prop =
              unique(lineage_gene_cat_sig)/unique(lineage_overall_sig)) %>% 
  ggplot(aes(x = cell_object_meta.final_lineage, y = lineage_gene_cat_sig,
             fill = gene_cat,
             label = round(lineage_gene_cat_sig,2)))+geom_bar(
               position = "dodge",stat='identity')+
  geom_text(position = position_dodge(width =1))+
  scale_fill_manual(values = DiscretePalette(n=4))+
  theme(axis.text.x = element_text(angle=90))+ggtitle("padj <=0.05 ")

p_lineage_sig_count + p_lineage_prop

p1 = sig_genes_per_ct_res %>% filter(padj<=0.05) %>%
  ggplot()+geom_bar(mapping = aes(x = gene_cat, 
                                  fill = cell_object_meta.final_lineage),
                    position = "dodge")+
  scale_fill_manual(values = DiscretePalette(n=4))+
  theme(axis.text.x = element_text(angle=90))+ggtitle("padj <=0.05 ")

p2 = sig_genes_per_ct_res %>% filter(padj<=0.01) %>%
  ggplot()+geom_bar(mapping = aes(x = gene_cat, fill = cell_object_meta.final_lineage),
                    position = "dodge")+
  scale_fill_manual(values = DiscretePalette(n=4))+
  theme(axis.text.x = element_text(angle=90))+ggtitle("padj <=0.01 ")

p_total_tested + p1 + p2


gene_order_60p <- readr::read_csv(
  "gene_order_peak_change_time60p_rowAnnot.csv")
#gene_order_60p$genes
dim(gene_order_60p)
p60_order <- gene_order_60p$genes
bi_genes <- gene_order_60p$genes[gene_order_60p$gene_mode=="bimodal"]
dim(gene_by_lumen_tma_corrected)
p60_order_unimodal_only <- p60_order[!p60_order %in% bi_genes]
length(p60_order_unimodal_only)

gene_mode_and_sc_clusters <- data.frame(gene_mode = gene_order_60p$gene_mode,
                                        genes = gene_order_60p$genes,
                                        sc_clusters = gene_annot[gene_order_60p$genes,"sc_clusters"])

gene_mode_and_sc_clusters$sc_clusters[is.na(gene_mode_and_sc_clusters$sc_clusters)] <- "bimodal"

library(circlize)
pseudotime_col_fun = colorRamp2(c(0, 500, 1500), c("white", "lightblue", "darkblue")) 
cell_col_fun = colorRamp2(c(-5, 0, 5), c("darkblue", "white", "firebrick")) 
pathology_col_fun = colorRamp2(c(0, 50, 100), c("white", "lightblue", "darkblue")) 
sample_type_label <- c("Unaffected","LF","MF","INT")
sample_type_label_color <- DiscretePalette(n=4)
names(sample_type_label_color) <- unique(sample_type_label)

dim(gene_by_lumen_tma_corrected)

top_annot_df <- read.csv(
  "./time_ordered_lumen_annot.csv")


final_gene_ht <- Heatmap(gene_by_lumen_tma_corrected[p60_order,],
        row_split = factor(gene_mode_and_sc_clusters$sc_clusters,
                           levels = c("bimodal","4","3","1","2")),
        cluster_columns =FALSE,
        cluster_rows = FALSE,
        show_column_names = FALSE,
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

sig_gene_per_ct_m <- sig_genes_per_ct_res %>% 
  select(gene,celltype,padj) %>%
  mutate(padj = ifelse(padj <=0, 1e-5,padj)) %>%  
  filter(padj <=0.05) %>%
  mutate(padj = -log10(padj)) %>%
  tidyr::pivot_wider(names_from = celltype,values_from = padj) 


sig_gene_per_ct_mtx <- as.matrix(sig_gene_per_ct_m[,2:ncol(sig_gene_per_ct_m)])
rownames(sig_gene_per_ct_mtx) <- sig_gene_per_ct_m$gene
sig_gene_per_ct_mtx_filled <- data.frame(sig_gene_per_ct_mtx)[p60_order,]
rownames(sig_gene_per_ct_mtx_filled) <- p60_order



padj_col_fun = colorRamp2(c(1, 2, 5), c("pink", "red", "red3")) 
ct_specfic_ht <- 
  Heatmap(sig_gene_per_ct_mtx_filled,cluster_rows = FALSE,col=padj_col_fun,
          cluster_columns = FALSE,na_col = "white",
          show_row_names = FALSE,
          row_split = factor(gene_mode_and_sc_clusters$sc_clusters,
                             levels = c("bimodal","4","3","1","2")),
          name="-log10(padj)",column_names_gp = gpar(fontsize=8),
          width = unit(4, "cm"),
          show_heatmap_legend = FALSE)

pdf("gene_test_per_ct_final_heatmaps.pdf",
    width = 12,height = 18)
ct_specfic_ht + final_gene_ht
dev.off()


### per celltype gene test results
sig_genes_per_ct_res %>% 
  mutate(padj = ifelse(padj <=0, 1e-5,padj)) %>% 
  ggplot() + geom_point(mapping = aes(x = celltype, y = gene,
                                      size = -log10(padj)),pch=1)+
  theme_bw(base_size = 15)+theme(axis.text.x = element_text(angle = 90))+
  theme(axis.text.y = element_text(size = 5))+
  facet_grid(cols = vars(cell_object_meta.final_lineage),
             rows = vars(gene_cat),
             scales = "free")

ggsave("./per_ct_result.png",
       width = 10,height = 20)
