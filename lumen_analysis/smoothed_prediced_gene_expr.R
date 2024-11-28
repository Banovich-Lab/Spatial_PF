

lumen_nuclei_sce_t4_ts <- readRDS("analysis/revision_2024/lumen_order_t4/lumen_nuclei_RNA_sce_t4_tradeseq.rds")
lumen_nuclei_sce <- readRDS("analysis/revision_2024/lumen_order_t4/lumen_nuclei_RNA_sce_t4_ordered.rds")

rowData(lumen_nuclei_sce_t4_ts)$assocRes <-
  associationTest(lumen_nuclei_sce_t4_ts,
                  lineages = T,
                  l2fc = 0,
                  nPoints = 10,
                  contrastType ="start")


de_genes <- rowData(lumen_nuclei_sce_t4_ts)$assocRes
de_genes$FDR <- p.adjust(de_genes$pvalue,method = "fdr")
sum(de_genes$FDR <= 0.01)
dim(de_genes[de_genes$FDR <= 0.01,])
de_genes[de_genes$FDR <= 0.01,]

de_genes <- rownames(de_genes[de_genes$FDR <= 0.01,])

npoints <- 2500
yhatSmooth <- predictSmooth(lumen_nuclei_sce_t4_ts, 
                            gene = de_genes, 
                            nPoints = npoints, tidy = FALSE)

order_rows <- yhatSmooth[de_genes,1:npoints]
order_rows <- apply(order_rows, 1, which.max)
window_size <- 300
gene_by_lumen <- yhatSmooth
rollmean_gene <- frollmean(as.data.table(t(yhatSmooth)), window_size)

gene_rank <- lapply(rollmean_gene, function(x) which.max(x[window_size:length(x)]))
gene_rank <- unlist(gene_rank)
names(gene_rank) <- colnames(t(yhatSmooth))
gene_order_smooth <- names(gene_rank)[order(gene_rank)]
gene_order_smooth
names(sort(order_rows))
breaks <-  c(-3,-1,0,1,3)
library(RColorBrewer)
#colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length( breaks))
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[gene_order_smooth,1:npoints]))),
                       cluster_cols = FALSE,cluster_rows =F,
                       color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks)),
                       breaks =breaks,
                       show_rownames = T, border_color = NA, 
                       heatmap_legend_param = list( labels_gp = gpar(fontsize = 7),
                                                    title_gp = gpar(fontsize = 7, fontface = "bold")),
                       show_colnames = F,name = "Cell Niche \n z-scored",
                       fontsize = unit("7","pt"))
heatSmooth
