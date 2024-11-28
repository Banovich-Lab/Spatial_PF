library(scater)
library(ggplot2)
library(mclust)
library(slingshot)
library(tradeSeq)
library(dplyr)
library(ComplexHeatmap)
library(data.table)
library(circlize)

order_feature_by_peak <- function(lumen_by_feature,windowsize = 50){
  gene_by_lumen <- t(lumen_by_feature)
  rollmean_gene <- frollmean(as.data.table(lumen_by_feature), windowsize,algo = "exact")
  
  gene_rank <- lapply(rollmean_gene, 
                      function(x) which.max(x[windowsize:length(x)]))
  
  gene_rank <- unlist(gene_rank)
  
  names(gene_rank) <- colnames(lumen_by_feature)
  gene_order_smooth <- names(gene_rank)[order(gene_rank)]
  gene_order_smooth
}

predicted_expr_along_time <- function(models = lumen_nuclei_RNA_sce_ts_knot13,
                                      genes,
                                      sce,
                                      scale = TRUE){
  predicted_gene_expr_lumen <- data.frame(predictCells(models = models, 
                                                       gene = genes))
  
  
  predicted_gene_expr_lumen_size_normed <- apply(predicted_gene_expr_lumen,1,
                                                 function(x){
                                                   log(x)-log(sce[,colnames(predicted_gene_expr_lumen)]$num_cells_lumen) 
                                                 })
  if(scale){
    predicted_gene_expr_lumen_size_normed <- scale(predicted_gene_expr_lumen_size_normed)
  } 
  
  predicted_gene_expr_lumen_size_normed_zscore_time_order <- 
    predicted_gene_expr_lumen_size_normed[order(models$crv[,1]),]
  
  return(predicted_gene_expr_lumen_size_normed_zscore_time_order)
}

lumen_cell_types_sce_t4_ts <- readRDS("analysis/revision_2024/lumen_order_t4/lumen_cell_types_sce_t4_tradeseq.rds")
lumen_cell_types_sce <- readRDS("output/savedRDS/revisions2024/lumen_cell_types_sce.all.ct.rds")

rowData(lumen_cell_types_sce_t4_ts)$assocRes <-kn
  associationTest(lumen_cell_types_sce_t4_ts,
                  lineages = T,
                  l2fc = 0.5,
                  nPoints = 12,
                  contrastType ="start")


de_celltype <- rowData(lumen_cell_types_sce_t4_ts)$assocRes
de_celltype$FDR <- p.adjust(de_celltype$pvalue,method = "fdr")
sum(de_celltype$FDR <= 0.01)
de_celltype[de_celltype$FDR <= 0.01,]

dir.create("output/figures/revisions2024/sup17_18/")
write.csv(de_celltype,"output/figures/revisions2024/sup17_18/celltype_assocTes_l2fc0.5.csv")

lumen_cell_types_sce_t4_m <- 
  predicted_expr_along_time(sce = lumen_cell_types_sce,
                            models = lumen_cell_types_sce_t4_ts,
                            genes = rownames(de_celltype[de_celltype$FDR <= 0.01,]))

lumen_cell_types_sce_t4_m <- predictSmooth(lumen_cell_types_sce_t4_ts,
                                           gene = rownames(de_celltype[de_celltype$FDR <= 0.01,]),
                                           tidy=FALSE,
                                           nPoints =2000)
library(kernlab)
ct_mat <- t(scale(t(lumen_cell_types_sce_t4_m[order_feature_by_peak(t(lumen_cell_types_sce_t4_m),windowsize = 100),])))
set.seed(100)
sp_clusters <- specc(ct_mat,centers=3)
split <- data.frame(cutree(hclust(dist(ct_mat)), k = 3))
specc4_split <- data.frame(sp_clusters)

pdf("output/figures/revisions2024/sup17_18/celltype_assocTes_l2fc0.5.pdf",
    width = 6,height = 5)
Heatmap(ct_mat,
        cluster_columns = FALSE,
        cluster_rows = F,
        row_gap = unit(0.5, "mm"),
        row_split = factor(specc4_split$sp_clusters, levels = c(1,2,3)),
        show_column_names = FALSE,
        name = "Celltype\nz-scored")
dev.off()

lumen_nuclei_sce_t4_ts <- readRDS("analysis/revision_2024/lumen_order_t4/offset_total_trans/gene_specc/lumen_nuclei_RNA_sce_ts_knot10_offset_tx_tma_rm.rds")

rowData(lumen_nuclei_sce_t4_ts)$assocRes
de_genes <- rowData(lumen_nuclei_sce_t4_ts)$assocRes
de_genes$FDR <- p.adjust(de_genes$pvalue,method = "fdr")
sum(de_genes$FDR <= 0.01)
dim(de_genes[de_genes$FDR <= 0.01,])
de_genes[de_genes$FDR <= 0.01,]
dim(de_genes)
write.csv(de_genes,file = "output/figures/revisions2024/sup17_18/genes_associate_test_pseudotime341_genes.csv")


## Celltype Niche
lumen_cell_types_sce_t4_ts <- readRDS("analysis/revision_2024/lumen_order_t4/lumen_cell_types_sce_t4_tradeseq.rds")
lumen_cell_niche_sce <- readRDS("output/savedRDS/revisions2024/lumen_cellniche_sce.rds")
sum(colnames(lumen_cell_types_sce_t4_ts) == colnames(lumen_cell_niche_sce))
lumen_cell_niche_sce$Lineage1 <- lumen_cell_types_sce_t4_ts$crv[,1]
assays(lumen_cell_niche_sce)[["CNiche_count"]][is.na(assays(lumen_cell_niche_sce)[["CNiche_count"]])] <- 0 
test_niches <- colSums(t(assays(lumen_cell_niche_sce)[["CNiche_count"]]) >=3) >=175
sum(test_niches)
test_niches <- names(test_niches)[test_niches]
test_niches
set.seed(100)
evaluateK(
  counts =
    assays(lumen_cell_niche_sce)[["CNiche_count"]],
  family = "nb",
  k = 5:10,
  nGenes = 12,
  cellWeights = matrix(1,nrow = ncol(lumen_cell_niche_sce)),
  pseudotime = lumen_cell_niche_sce$Lineage1,
  offset = log(lumen_cell_niche_sce$num_cells_lumen))

lumen_cell_niche_ts <- fitGAM(counts =
                                         assays(lumen_cell_niche_sce)[["CNiche_count"]],
                                       # sds = lumen_cell_types_sce_prop_PCA13$slingshot, 
                                       family = "nb",
                                       nknots = 7,
                                       genes = test_niches,
                                       cellWeights =matrix(1,nrow = ncol(lumen_cell_niche_sce)), 
                                       pseudotime = lumen_cell_niche_sce$Lineage1,
                                       offset = log(lumen_cell_niche_sce$num_cells_lumen))

saveRDS(lumen_cell_niche_ts,"analysis/revision_2024/lumen_order_t4/lumen_cell_niche_sce_tradeseq.rds")

rowData(lumen_cell_niche_ts)$assocRes <-
  associationTest(lumen_cell_niche_ts,
                  lineages = T,
                  l2fc = 0,
                  nPoints = 7,
                  contrastType ="start")


de_cniche <- rowData(lumen_cell_niche_ts)$assocRes
de_cniche$FDR <- p.adjust(de_cniche$pvalue,method = "fdr")
sum(de_cniche$FDR <= 0.01)
dim(de_cniche[de_cniche$FDR <= 0.01,])
de_cniche[de_cniche$FDR <= 0.01,]
write.csv(de_cniche,file = "output/figures/revisions2024/sup17_18/cniche_assocTest_lfc0.csv")

## cniche heatmap plot

de_cellniche_label <- rownames(de_cniche)[de_cniche$FDR<=0.01]
yhatSmooth <- predictSmooth(lumen_cell_niche_ts, 
                            gene = de_cellniche_label, 
                            nPoints = 2000, tidy = FALSE)

yhatSmooth <- yhatSmooth[de_cellniche_label,]
breaks <-  c(-2,-1,0,1,2)
library(RColorBrewer)
cellniche_mat <- order_feature_by_peak(t(yhatSmooth),windowsize = 100)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[cellniche_mat,]))),
                       cluster_cols = FALSE,cluster_rows =F,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
                       breaks =breaks,
                       show_rownames = T,
                       border_color = NA, 
                       heatmap_legend_param = list( labels_gp = gpar(fontsize = 7),
                                                    title_gp = gpar(fontsize = 7, fontface = "bold")),
                       show_colnames = F,name = "Cell Niche \n z-scored",
                       fontsize = unit("7","pt"))
heatSmooth
pdf("output/figures/revisions2024/sup17_18/cniche_heatmap.pdf",width = 4,height = 2)
print(heatSmooth)
dev.off()


## tx Niche
lumen_cell_types_sce_t4_ts <- readRDS("analysis/revision_2024/lumen_order_t4/lumen_cell_types_sce_t4_tradeseq.rds")
lumen_tx_niche_sce <- readRDS("output/savedRDS/revisions2024/lumen_tranx_niche_sce.rds")
sum(colnames(lumen_cell_types_sce_t4_ts) == colnames(lumen_tx_niche_sce))
lumen_tx_niche_sce$Lineage1 <- lumen_cell_types_sce_t4_ts$crv[,1]
test_tniches <- colSums(t(assays(lumen_tx_niche_sce)[["tranx_nichek12_count"]]) >=3) >= 175
sum(test_tniches)
test_tniches <- names(test_tniches)[test_tniches]
test_tniches
set.seed(100)
evaluateK(
  counts =
    assays(lumen_tx_niche_sce)[["tranx_nichek12_count"]],
  family = "nb",
  k = 8:15,
  nGenes = 12,
  cellWeights = matrix(1,nrow = ncol(lumen_tx_niche_sce)),
  pseudotime = lumen_tx_niche_sce$Lineage1,
  offset = log(lumen_tx_niche_sce$num_cells_lumen))

lumen_tx_niche_sce_ts <- fitGAM(counts =
                                    assays(lumen_tx_niche_sce)[["tranx_nichek12_count"]],
                                  # sds = lumen_cell_types_sce_prop_PCA13$slingshot, 
                                  family = "nb",
                                  nknots = 10,
                                  genes = test_tniches,
                                  cellWeights =matrix(1,nrow = ncol(lumen_tx_niche_sce)), 
                                  pseudotime = lumen_tx_niche_sce$Lineage1,
                                  offset = log(lumen_tx_niche_sce$num_cells_lumen))

saveRDS(lumen_tx_niche_sce_ts,"analysis/revision_2024/lumen_order_t4/lumen_tx_niche_sce_tradeseq.rds")

rowData(lumen_tx_niche_sce_ts)$assocRes <-
  associationTest(lumen_tx_niche_sce_ts,
                  lineages = T,
                  l2fc = 0,
                  nPoints = 10,
                  contrastType ="start")


de_tniche <- rowData(lumen_tx_niche_sce_ts)$assocRes
de_tniche$FDR <- p.adjust(de_tniche$pvalue,method = "fdr")
sum(de_tniche$FDR <= 0.01)
dim(de_tniche[de_tniche$FDR <= 0.01,])
de_tniche[de_tniche$FDR <= 0.01,]
write.csv(de_tniche,file = "output/figures/revisions2024/sup17_18/tniche_assocTest_lfc0.csv")

## tniche heatmap

de_tniche_label <- rownames(de_tniche)[de_tniche$FDR<=0.01]
yhatSmooth <- predictSmooth(lumen_tx_niche_sce_ts, 
                            gene = de_tniche_label, 
                            nPoints = 2000, tidy = FALSE)

yhatSmooth <- yhatSmooth[de_tniche_label,]
breaks <-  c(-2,-1,0,1,2)
library(RColorBrewer)
#colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length( breaks))
tniche_mat <- t(scale(t(yhatSmooth[order_feature_by_peak(t(yhatSmooth),100),])))
heatSmooth <- pheatmap(tniche_mat,
                       cluster_cols = FALSE,cluster_rows =F,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
                       breaks =breaks,
                       show_rownames = T, border_color = NA, 
                       heatmap_legend_param = list( labels_gp = gpar(fontsize = 7),
                                                    title_gp = gpar(fontsize = 7, fontface = "bold")),
                       show_colnames = F,name = "Transcript Niche \n z-scored",
                       fontsize = unit("7","pt"))
heatSmooth
pdf("output/figures/revisions2024/sup17_18/tniche_heatmap.pdf",width = 4,height = 2)
print(heatSmooth)
dev.off()

