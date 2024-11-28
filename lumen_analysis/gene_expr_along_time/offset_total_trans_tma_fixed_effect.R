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

lumen_cell_types_sce <- readRDS(file="output/savedRDS/revisions2024/lumen_cell_types_sce.all.ct.rds")
dim(lumen_cell_types_sce)

test_celltypes <- names(rowSums(assay(lumen_cell_types_sce,"final_CT_count")>=3) >=10)[rowSums(assay(lumen_cell_types_sce,"final_CT_count")>=3) >=10]
length(test_celltypes)


lumen_cell_types_sce$sample_affect <-
  factor(lumen_cell_types_sce$sample_affect,
         levels = c("Unaffected","Less Affected",
                    "More Affected"))

sample_level <- unique(data.frame(lumen_cell_types_sce$sample,
                                  lumen_cell_types_sce$percent_pathology))

sample_level <- sample_level$lumen_cell_types_sce.sample[order(sample_level$lumen_cell_types_sce.percent_pathology)]
lumen_cell_types_sce$sample <- factor(lumen_cell_types_sce$sample,levels = sample_level)

lumen_cell_types_sce$sample_type <-
  factor(lumen_cell_types_sce$sample_type,levels = c("Unaffected","LF","MF","INT"))

lumen_order_by_t4 <- readr::read_csv("analysis/revision_2024/lumen_order_t4/pseudotime_by_niche_073024/lumen_metadata_niche_pseudotimes_073024.csv")

lumen_order_by_t4 <- data.frame(lumen_order_by_t4)
rownames(lumen_order_by_t4) <- lumen_order_by_t4$lumen_id
lumen_cell_types_sce$pseudotime_rank_trans_t4_tie_avg <- lumen_order_by_t4[colnames(lumen_cell_types_sce),
                                                                           "pseudotime_rank_trans_t4_tie_avg"]
lumen_nuclei_RNA_sce <- readRDS(file="./output/savedRDS/revisions2024/lumen_nuclei_sce.rds")
dim(lumen_nuclei_RNA_sce)
#lumen_nuclei_RNA_sce$lumen_size
lumen_nuclei_RNA_sce$lumen_id <- colnames(lumen_nuclei_RNA_sce)

# using the pseudotime ordering from celltype based analysis
lumen_nuclei_RNA_sce$Lineage1 <- lumen_cell_types_sce[,colnames(lumen_nuclei_RNA_sce)]$pseudotime_rank_trans_t4_tie_avg
lumen_nuclei_RNA_sce$total_trans <- as.numeric(colSums(assay(lumen_nuclei_RNA_sce)))

### knot=10
test_genes <- rowSums(assay(lumen_nuclei_RNA_sce,"trans_count") >=3) >=10
test_genes <- names(test_genes[test_genes])
length(test_genes)
set.seed(100)
# evaluateK(assay(lumen_nuclei_RNA_sce,"trans_count"),k = 8:15,
#           pseudotime = lumen_nuclei_RNA_sce$Lineage1,
#           cellWeights = matrix( rep(1,ncol(lumen_nuclei_RNA_sce)),ncol=1),
#           nGenes=50,
#           U = model.matrix(~0+tma,data.frame(colData(lumen_nuclei_RNA_sce))),
#           BPPARAM = BiocParallel::MulticoreParam(workers = 4),
#           offset = log(lumen_nuclei_RNA_sce$total_trans),
#           family = "nb")

set.seed(100)
lumen_nuclei_RNA_sce_ts_knot10_tma_fixed <- fitGAM(
  counts =
    lumen_nuclei_RNA_sce@assays@data$trans_count,
  pseudotime = lumen_nuclei_RNA_sce$Lineage1,
  cellWeights = matrix( rep(1,ncol(lumen_nuclei_RNA_sce)),ncol=1),
  U = model.matrix(~0+tma,data.frame(colData(lumen_nuclei_RNA_sce))),
  family = "nb", 
  nknots = 10,
  genes = test_genes,
  offset = log(lumen_nuclei_RNA_sce$total_trans))

mean(rowData(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed)$tradeSeq$converged)
rowData(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed)$assocRes <-
  associationTest(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed,
                  lineages = T,
                  l2fc = 0,
                  nPoints = 10,
                  contrastType ="start")


# lumen_nuclei_RNA_sce_ts_knot10_tma_fixed <- readRDS(
#   "analysis/revision_2024/lumen_order_t4/offset_total_trans/lumen_nuclei_RNA_sce_ts_knot10_offset_tx_tma_rm.rds")
# )
assocRes <- rowData(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed)$assocRes
de_genes_knot10_tma <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.01)
]
length(de_genes_knot10_tma)
grep("C20",de_genes_knot10_tma)


saveRDS(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed,
        "analysis/revision_2024/lumen_order_t4/offset_total_trans/lumen_nuclei_RNA_sce_ts_knot10_offset_tx_tma_rm.rds")



predicted_gene_expr_tma_normed_cells <- 
  predictSmooth(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed,
               gene = test_genes,
               tidy=FALSE,
               nPoints=1000)
plot(predicted_gene_expr_tma_normed_cells[1,])
## predict per cell, tma corrected
pseudotime_predictors <- grep("t1",colnames(colData(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed)$tradeSeq$X))

lumen_by_gene_tma_corrected <- (colData(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed)$tradeSeq$X)[,pseudotime_predictors] %*% t(rowData(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed)$tradeSeq$beta[[1]][,pseudotime_predictors])
dim(lumen_by_gene_tma_corrected)
rownames(lumen_by_gene_tma_corrected) <- colnames(lumen_nuclei_RNA_sce_ts_knot10_tma_fixed)

gene_by_lumen_tma_corrected <- t(exp(lumen_by_gene_tma_corrected))[,order(lumen_nuclei_RNA_sce$Lineage1)]
plot(gene_by_lumen_tma_corrected[1,])

gene_by_lumen_tma_corrected <- scale(t(gene_by_lumen_tma_corrected))
gene_by_lumen_tma_corrected <- t(gene_by_lumen_tma_corrected)
time_ordered_lumens <- colnames(gene_by_lumen_tma_corrected)
length(de_genes_knot10_tma)

gene_by_lumen_tma_corrected <- gene_by_lumen_tma_corrected[de_genes_knot10_tma,]
dim(gene_by_lumen_tma_corrected)

## find bimodal genes and split the genes by bimodal or not
set.seed(100)
h_clust_cut <- (cutree(hclust(dist(gene_by_lumen_tma_corrected)),
                       k = 20))

p_ht <- Heatmap(gene_by_lumen_tma_corrected,
                show_row_names = T,
                cluster_columns = FALSE,
                show_column_names = FALSE,
                row_split = h_clust_cut,
                row_gap = unit(0.5, "mm"),
                row_names_gp = gpar(fontsize=4),

                row_title_gp  = gpar(fontsize=7),
                name = "Gene Expr\n z-scores",
                heatmap_legend_param = 
                  list(    legend_height = unit(2, "cm"),
                           legend_gp = gpar(fontsize=7),
                           title_gp = gpar(fontsize = 7, 
                                           fontface = "bold")))
dir.create("./analysis/revision_2024/lumen_order_t4/offset_total_trans/")
pdf("./analysis/revision_2024/lumen_order_t4/offset_total_trans/hclusters_offset_total_tx_tma_corrected.pdf",
    height = 16, width = 8)
print(p_ht)
dev.off()

# time_ordered_lumens
library(circlize)
pseudotime_col_fun = colorRamp2(c(0, 500, 1500), c("white", "lightblue", "darkblue")) 

cell_col_fun = colorRamp2(c(-5, 0, 5), c("darkblue", "white", "firebrick")) 
pathology_col_fun = colorRamp2(c(0, 50, 100), c("white", "lightblue", "darkblue")) 
sample_type_label <- lumen_nuclei_RNA_sce[,time_ordered_lumens]$sample_type
sample_type_label_color <- DiscretePalette(n=4)
names(sample_type_label_color) <- unique(sample_type_label)
pseudotime_value <- lumen_nuclei_RNA_sce[,time_ordered_lumens]$Lineage1
sample_perc_pathology <- lumen_nuclei_RNA_sce[,time_ordered_lumens]$percent_pathology


gene_trend <- h_clust_cut %in% c(17,20,8,13,7,14)

bi_gene <- names(h_clust_cut)[gene_trend]
single_mode_gene <- names(h_clust_cut)[!gene_trend]
gene_mode <- rep("unimodal",nrow(gene_by_lumen_tma_corrected))

names(gene_mode) <- rownames(gene_by_lumen_tma_corrected)
gene_mode[bi_gene] <-"bimodal"
ht_orders <- do.call(rbind,lapply(names(row_order(p_ht)), 
                                  function(x) data.frame(x,row_order(p_ht)[[as.character(x)]]) ))

ht_orders$gene <- rownames(gene_by_lumen_tma_corrected)[ht_orders$row_order.p_ht...as.character.x...]
ht_orders$gene_mode <- "unimodal"
ht_orders$gene_mode[ht_orders$gene %in%bi_gene ] <- "bimodal"
ht_orders$h_cluster <- ht_orders$x
ht_orders$row_order.p_ht...as.character.x... <- NULL

mean_expr_per_lumen_in_chunks <- lapply(unique(ht_orders$x), 
                                        function(h_c){
                                          rowMeans(t(gene_by_lumen_tma_corrected)[,ht_orders$gene[ht_orders$x == h_c],
                                                                                                           drop=F])
                                          
                                        })
mean_expr_per_lumen_in_chunks_max_at <- lapply(mean_expr_per_lumen_in_chunks,which.max)
names(mean_expr_per_lumen_in_chunks_max_at) <- paste0("hcluster_",(unique(ht_orders$x)))

hchunks_order <- data.frame(do.call(rbind, mean_expr_per_lumen_in_chunks_max_at)) %>%
  mutate(hcluster_label = rownames(do.call(rbind, mean_expr_per_lumen_in_chunks_max_at))) %>%
  mutate(o_h = rank(as.numeric(TILD080MF_99),ties.method = "random")) %>% arrange(o_h)

chunck_gene_peak <- hchunks_order$hcluster_label
chunck_gene_peak <- gsub("hcluster_","",chunck_gene_peak)
chunck_gene_peak
#  [1] "15" "3"  "5"  "1"  "16" "6"  "18" "9"  "13" "12" "20" "2"  "10" "8"  "14" "19" "4"  "7"  "11" "17"
# 
p_ht_reorder_chunk <- Heatmap(gene_by_lumen_tma_corrected,
                show_row_names = T,
                cluster_rows = F,
                cluster_columns = FALSE,
                show_column_names = FALSE,
                row_split = factor(ht_orders$h_cluster[match(rownames(gene_by_lumen_tma_corrected),ht_orders$gene)],levels = chunck_gene_peak),
                row_gap = unit(0.5, "mm"),
                row_names_gp = gpar(fontsize=4),
                top_annotation =
                  HeatmapAnnotation(pseudotime = pseudotime_value,
                                    sample_type = sample_type_label,
                                    percent_pathology = sample_perc_pathology,
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
p_ht_reorder_chunk

chunck_gene_peak <- c("1","4","19","10","3","18","6","9","16","20","13","12","2","14","7","17","5","15","11","8")
p_ht_reorder_chunk <- Heatmap(gene_by_lumen_tma_corrected,
                              show_row_names = T,
                              cluster_rows = F,
                              cluster_columns = FALSE,
                              show_column_names = FALSE,
                              row_split = factor(ht_orders$h_cluster[match(rownames(gene_by_lumen_tma_corrected),
                                                                           ht_orders$gene)],
                                                 levels =chunck_gene_peak[!chunck_gene_peak %in%  c(17,20,8,13,7,14)]),
                              row_gap = unit(0.1, "mm"),
                              row_names_gp = gpar(fontsize=4),

                              top_annotation =
                                HeatmapAnnotation(pseudotime = pseudotime_value,
                                                  sample_type = sample_type_label,
                                                  percent_pathology = sample_perc_pathology,
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
pdf("analysis/revision_2024/lumen_order_t4/offset_total_trans/hclusters_chunk_ordered_offset_total_tx_tma_corrected.pdf",
    width = 6,height = 18)
p_ht_reorder_chunk
dev.off()

g_line <- data.frame(t(gene_by_lumen_tma_corrected))  %>% 
  mutate(time_order = 1:ncol(gene_by_lumen_tma_corrected)) %>%
  tidyr::pivot_longer(cols = 1:nrow(gene_by_lumen_tma_corrected),
                      names_to = "gene") %>% 
  merge(ht_orders) %>%
  ggplot()+
  geom_smooth(mapping = aes(x = time_order,y = value,group=gene))+
  facet_wrap(.~h_cluster)

ggsave(plot=g_line,
       width = 10,
       height = 10,
      filename = "analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_line_plot_offset_tx_tma_corrected_ht_clusters.png")


gene_by_lumen_standard <- gene_by_lumen_tma_corrected
gene_by_lumen_standard <- t(apply(gene_by_lumen_standard,1,function(x) (x-min(x))/(max(x)-min(x))))

chane_time <- apply(t(gene_by_lumen_standard), 2,
                    function(x){
                      if(mean( x[1:10]) < 0.5 ){
                        increase_to1 <- c("increase",which(x >= 0.6 )[1])
                      } else {
                        c("decrease",which(x <= 0.5 )[1])
                      }
                      
                    })

chane_time_rank <- data.frame(t((as.matrix(chane_time)))) 
chane_time_rank$gene <- rownames(chane_time_rank)

chane_time_rank <- chane_time_rank %>%
  mutate(X2 = as.numeric(X2)) %>%
  group_by(X1) %>% 
  mutate(in_group_rank = case_when(
    (X1 == "increase" ) ~  rank(X2,ties.method = "random"),
    (X1 == "decrease" ) ~  rank(-X2,ties.method = "random")))
chane_time_rank$X1 <- factor(chane_time_rank$X1,levels =c("decrease","increase"))

chane_time_rank <- chane_time_rank %>% arrange(X1,in_group_rank)

change_time_order <- chane_time_rank$gene

## peak time after rolling mean 
win_size <- 300
rollmean_gene <- frollmean(as.data.table(t(gene_by_lumen_tma_corrected)),
                           win_size,algo = "exact")

gene_rank <- lapply(rollmean_gene, 
                    function(x) which.max(x[win_size:length(x)]))

gene_rank <- unlist(gene_rank)
names(gene_rank) <- colnames(t(gene_by_lumen_tma_corrected))
length(table(gene_rank))

gene_peak_time <- lumen_nuclei_RNA_sce[,colnames(gene_by_lumen_tma_corrected)[gene_rank]]$Lineage1
names(gene_peak_time) <- names(gene_rank)

gene_cat_df <- data.frame(genes = names(gene_peak_time),gene_peak_time)
gene_cat_df %>% ggplot()+
  geom_point(mapping= aes(y=genes,x = gene_peak_time))

gene_order_smooth <- names(gene_rank)[order(gene_rank)]
#gene_order_smooth
gene_cat_df$gene_order <- match(gene_cat_df$genes,gene_order_smooth)
gene_cat_df <- gene_cat_df %>% mutate(change_time_order = match(genes,
                                                                change_time_order)) %>%
  arrange(gene_peak_time,change_time_order)

 gene_cat_df$gene_order_peak_change_time60p <- gene_cat_df$genes
 
 gene_cat_df$gene_mode <- ht_orders$gene_mode[match(gene_cat_df$genes,ht_orders$gene)]
write.csv(gene_cat_df,
          "analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_order_peak_change_time60p_offset_tx_tma_rm_rowAnnot.csv")
write.csv(gene_by_lumen_tma_corrected[gene_cat_df$genes,],
          "analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_order_peak_change_time60p_offset_tx_tma_rm_mat.csv")

pdf("analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_order_peak_change_time60p_offset_tx_tma_rm.pdf",
    height = 18,width = 8)

Heatmap(gene_by_lumen_tma_corrected[gene_cat_df$genes,],
        show_row_names = T,
        cluster_rows = F,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_split =  ht_orders$gene_mode[match(gene_cat_df$genes,ht_orders$gene)],
        row_gap = unit(0.5, "mm"),
        row_names_gp = gpar(fontsize=4),
        top_annotation = 
          HeatmapAnnotation(pseudotime = pseudotime_value,
                            sample_type = sample_type_label,
                            percent_pathology =sample_perc_pathology,
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
dev.off()



pdf("analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_order_60p_peaktime_offset_tx_tma_rm_col_split_tma.pdf",
    height = 16,width = 14)

Heatmap(gene_by_lumen_tma_corrected[gene_cat_df$genes,],
        show_row_names = T,
        cluster_rows = F,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_split =  ht_orders$gene_mode[match(gene_cat_df$genes,ht_orders$gene)],
        column_split = lumen_nuclei_RNA_sce$tma[match(colnames(gene_by_lumen_tma_corrected),lumen_nuclei_RNA_sce$lumen_id)],
        row_gap = unit(0.5, "mm"),
        row_names_gp = gpar(fontsize=4),
        top_annotation = 
          HeatmapAnnotation(pseudotime = pseudotime_value,
                            sample_type = sample_type_label,
                            percent_pathology =sample_perc_pathology,
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
dev.off()



gene_order_60p <- readr::read_csv("analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_order_peak_change_time60p_offset_tx_tma_rm_rowAnnot.csv")
#gene_order_60p$genes
dim(gene_order_60p)
p60_order <- gene_order_60p$genes[gene_order_60p$genes %in% offset_total_tx_tma_rm_sig_genes]
bi_genes <- gene_order_60p$genes[gene_order_60p$gene_mode=="bimodal"]
dim(gene_by_lumen_tma_corrected)
p60_order_unimodal_only <- p60_order[!p60_order %in% bi_genes]
length(p60_order_unimodal_only)

library(kernlab)
set.seed(100)
dim(gene_by_lumen_tma_corrected)

### 4 gene clusters
set.seed(100)
dim(gene_by_lumen_tma_corrected)
sc <- specc(gene_by_lumen_tma_corrected, centers=4)

sc_clusters <- sc@.Data
length(sc_clusters)
names(sc_clusters) <- rownames(gene_by_lumen_tma_corrected)
sc_clusters <- sc_clusters[p60_order_unimodal_only]
length(sc_clusters)
sc_cluster_levels <- unique(sc_clusters)

bimodal_sc_clusters <- c()
sc_cluster_levels <- c(4,3,1,2)

pdf("./analysis/revision_2024/lumen_order_t4/offset_total_trans/specc_cluster_offset_total_tx_tma_corrected_change_time_in_chunk_4c.pdf",
    height = 16, width = 6)
write.csv(gene_by_lumen_tma_corrected[p60_order_unimodal_only,],
          "./analysis/revision_2024/lumen_order_t4/offset_total_trans/specc_cluster_offset_total_tx_tma_corrected_change_time_in_chunk_4c_mat.csv")

Heatmap(gene_by_lumen_tma_corrected[p60_order_unimodal_only,],
        show_row_names = T,
        cluster_rows = F,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_split = factor(sc_clusters,levels=sc_cluster_levels[!sc_cluster_levels %in% bimodal_sc_clusters]),
        row_gap = unit(0.5, "mm"),
        row_names_gp = gpar(fontsize=4),
        top_annotation =
          HeatmapAnnotation(pseudotime = pseudotime_value,
                            sample_type = sample_type_label,
                            percent_pathology = sample_perc_pathology,
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
dev.off()

sc_clusters_df <- data.frame(gene = names(sc_clusters), sc_clusters,check.names = FALSE) 
dim(sc_clusters_df)
length(unique(sc_clusters_df$gene))
dim(gene_by_lumen_tma_corrected)

g_line_specc_df <- data.frame(t(gene_by_lumen_tma_corrected[p60_order_unimodal_only,]),
                              check.names = FALSE) %>% 
  mutate(time_order = 1:ncol(gene_by_lumen_tma_corrected)) %>% 
  mutate(time_order = as.numeric(time_order)) %>% 
  tidyr::pivot_longer(cols = p60_order_unimodal_only,
                      names_to = "gene") %>%
  merge(sc_clusters_df,all=TRUE) 
length(unique(g_line_specc_df$gene))

write.csv(g_line_specc_df,"./analysis/revision_2024/lumen_order_t4/offset_total_trans/specc_cluster_offset_total_tx_tma_corrected_change_time_in_chunk_4c_line_plot_df.csv")

g_line_specc <- g_line_specc_df %>% 
  ggplot(aes(x = time_order,y = value))+
  geom_point(size=0.1,alpha=0.3,color="grey")+
  geom_smooth(mapping = aes(x = time_order,y = value),
              se=TRUE,colour="red",method="gam")+
  facet_wrap(.~sc_clusters)

ggsave(plot=g_line_specc,
       width = 10,
       height = 10,
       filename = "analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_line_specc4.png")

