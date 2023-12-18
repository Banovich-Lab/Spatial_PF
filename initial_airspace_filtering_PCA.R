###########################################
# Initial Lumen Filtering and PCA
# Author: Annika Vannan (avannan@tgen.org)
# Date: 12/18/2023
###########################################

## SETTING ENVIRONMENT ----
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(gplots)
library(tibble)
library(ggpubr)
library(ggrepel)
library(grid)
library(patchwork)
library(clustree)
library(ComplexHeatmap)
library(circlize)
library(randomcoloR)
library(ggpubr)
library(factoextra)
library(randomcoloR)
library(Cairo)
library(SingleCellExperiment)
library(tradeSeq)


## LOAD DATA ----
xenium <- readRDS("/Volumes/dback_scratch/avannan/full_xenium_11-09-23_rm_nuc.rds")
pretty_umap <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     axis.title = element_text(hjust = 1))

file_list <- list.files("/Volumes/dback_scratch/avannan/lumen_metadata_4explorer", full.names = TRUE, pattern = "final_lumens_091223.csv")
df_list <- lapply(file_list, read.csv)
df <- Reduce(rbind, df_list)
df %>% select(3) %>% unique() %>% dim()


## LOAD LUMEN FILES ----
# Keep only the cells that are in the Seurat object
keep_cells <- colnames(xenium)

# Load in lumen CSV files as list
lumen_files <- list.files("/Volumes/dback_scratch/avannan/corrected_lumen_ids", 
                          full.names = TRUE, pattern = ".csv")

# Get sample IDs from filenames and rename list
sample_ids <- unlist(lapply(
  str_split(lapply(str_split(lumen_files, "corrected_lumen_ids/"), 
                   function(XX) { XX[[2]]} ), "_TMA"),  function(XX) { XX[[1]] }))
names(lumen_files) <- sample_ids

# Create metadata for Xenium Seurat object
metadata <- readRDS("/Volumes/dback_scratch/avannan/lumen_pca/corrected_lumen_data.rds")
metadata_test <- metadata[colnames(xenium), ]

# Create metadata for use in Xenium Explorer
metadata_4explorer <- metadata %>%
  rownames_to_column(var = "cell_id") %>%
  filter(cell_id %in%  colnames(xenium)) %>%
  separate(cell_id, c("sample", "cell_id2"), sep = "_") %>%
  select(sample, cell_id2, lumen_id) %>%
  dplyr::rename(group = "lumen_id", cell_id = "cell_id2")


## CALCULATE MAX DISTANCE BETWEEN CELLS ----
# Make temporary dataframe for calculating max distances between cells in each lumen
mini_metadata <- xenium@meta.data %>%
  select(sample, cell_id, x_centroid, y_centroid, broad_CT5, lineage)
new_tmp <- metadata_4explorer %>%
  left_join(mini_metadata) %>%
  filter(group != 0)

# Calculate max distance between cells in each lumen
# Do very basic filtering beforehand
max_dist_list <- list()
for (sm in sample_ids[15:length(sample_ids)]) {
  message(paste("Processing", sm))
  sm_tmp <- new_tmp %>% 
    filter(sample == sm) %>%
    select(cell_id, group) %>%
    group_by(group)  %>%
    mutate(num_cells = length(cell_id)) %>%
    filter(num_cells < 1000, num_cells >= 10) %>%
    left_join(mini_metadata %>% filter(sample == sm))
  for (lumen in unique(sm_tmp$group)) {
    lumen_tmp <- sm_tmp %>% filter(group == lumen)
    calc_distances <- RANN::nn2(cbind(lumen_tmp$x_centroid, 
                                      lumen_tmp$y_centroid), 
                                k = nrow(lumen_tmp))
    max_dist_list[[sm]][[lumen]] <- lumen_tmp %>%
      left_join(data.frame(group = lumen, 
                           max_dist = max(calc_distances$nn.dists)))
  }
}
# Turn list into a dataframe and bin the data
max_dist_df_list <- lapply(max_dist_list, function(XX) { Reduce(full_join, XX) })
max_dist_df <- Reduce(full_join, max_dist_df_list) %>% unique()
max_dist_df2 <- max_dist_df %>% 
  select(sample, group, max_dist, num_cells) %>% 
  unique() %>%
  mutate(num_cells_bin = cut(num_cells, breaks = seq(0, 1000, 25)),
         max_dist_bin = cut(max_dist, breaks = seq(0, 5000, 25)))

# Save individual files with metadata on lumen sizes for use in Xenium Explorer
for (sm in sample_ids) {
  tmp <- new_tmp %>%
    filter(sample == sm) %>%
    left_join(max_dist_df2) %>%
    ungroup() %>%
    select(cell_id, max_dist_bin)
  colnames(tmp) <- c("cell_id", "group")
  write.csv(tmp, paste0("/Volumes/dback_scratch/avannan/lumen_metadata_4explorer/", 
                        sm, "_max_dist_bins_111523.csv"))

  tmp2 <- max_dist_df2 %>%
    filter(sample == sm) %>%
    left_join(max_dist_df) %>%
    ungroup() %>%
    select(cell_id, num_cells_bin)
  colnames(tmp2) <- c("cell_id", "group")
  write.csv(tmp2, paste0("/Volumes/dback_scratch/avannan/lumen_metadata_4explorer/",
                         sm, "_num_cell_bins_111523.csv"))
}

# Modify dataframe so that it can be added to Seurat object
dist_ncell_metadata_list <- lapply(sample_ids, function(XX) {
  max_dist_df2 %>%
  filter(sample == XX) %>%
  left_join(max_dist_df) %>%
  ungroup()
})
dist_ncell_metadata_df <- Reduce(rbind, dist_ncell_metadata_list) %>%
  dplyr::rename(lumen_id = "group") %>%
  as.data.frame()
rownames(dist_ncell_metadata_df) <- paste0(dist_ncell_metadata_df$sample, "_", dist_ncell_metadata_df$cell_id)


## ADD DATA TO SEURAT OBJECT ----
# Add new metadata to Seurat object
lumen_xenium <- AddMetaData(xenium, metadata)
lumen_xenium <- AddMetaData(lumen_xenium, dist_ncell_metadata_df)

# Only keep cells that are assigned to lumens
lumen_xenium <- subset(lumen_xenium, subset = lumen_id != 0)

# # Add more metadata
lumen_xenium$num_cells <- cell_lumen_metadata$num_cells[match(lumen_xenium$lumen_id, cell_lumen_metadata$lumen_id)]
lumen_xenium$max_dist <- cell_lumen_metadata$max_dist[match(lumen_xenium$lumen_id, cell_lumen_metadata$lumen_id)]
lumen_xenium$num_cells_bin <- cell_lumen_metadata$num_cells_bin[match(lumen_xenium$lumen_id, cell_lumen_metadata$lumen_id)]
lumen_xenium$max_dist_bin <- cell_lumen_metadata$max_dist_bin[match(lumen_xenium$lumen_id, cell_lumen_metadata$lumen_id)]

# Get counts of cells per lineage
lineage_lumen_counts <- lumen_xenium@meta.data %>%
  group_by(lumen_id, lineage) %>%
  mutate(num_cells_lineage = length(cell_id)) %>%
  ungroup() %>%
  select(lumen_id, lineage, num_cells_lineage) %>%
  unique() %>%
  pivot_wider(id_cols = "lumen_id", names_from = "lineage", values_from = "num_cells_lineage") %>%
  mutate(across(2:5, function(XX) replace_na(XX, 0))) %>%
  group_by(lumen_id) %>%
  mutate(sum = sum(Epithelial, Endothelial, Immune, Mesenchymal),
         endo_prop = Endothelial/sum,
         epi_prop = Epithelial/sum,
         imm_prop = Immune/sum,
         mes_prop = Mesenchymal/sum) %>%
  ungroup()

# Also add this to metadata
lumen_xenium$endo_lumen_count <- lineage_lumen_counts$Endothelial[match(lumen_xenium$lumen_id, lineage_lumen_counts$lumen_id)]
lumen_xenium$epi_lumen_count <- lineage_lumen_counts$Epithelial[match(lumen_xenium$lumen_id, lineage_lumen_counts$lumen_id)]
lumen_xenium$imm_lumen_count <- lineage_lumen_counts$Immune[match(lumen_xenium$lumen_id, lineage_lumen_counts$lumen_id)]
lumen_xenium$mes_lumen_count <- lineage_lumen_counts$Mesenchymal[match(lumen_xenium$lumen_id, lineage_lumen_counts$lumen_id)]
lumen_xenium$endo_lumen_prop <- lineage_lumen_counts$endo_prop[match(lumen_xenium$lumen_id, lineage_lumen_counts$lumen_id)]
lumen_xenium$epi_lumen_prop <- lineage_lumen_counts$epi_prop[match(lumen_xenium$lumen_id, lineage_lumen_counts$lumen_id)]
lumen_xenium$imm_lumen_prop <- lineage_lumen_counts$imm_prop[match(lumen_xenium$lumen_id, lineage_lumen_counts$lumen_id)]
lumen_xenium$mes_lumen_prop <- lineage_lumen_counts$mes_prop[match(lumen_xenium$lumen_id, lineage_lumen_counts$lumen_id)]
write.csv(dist_ncell_metadata_df, "/Volumes/dback_scratch/avannan/lumen_pca/dist_ncell_metadata_df_111523.csv")


## SET UP CELL TYPE "EXPRESSION" (COMPOSITION) DATA ----
# Get metadata and cell type "expression" data
lumen_metadata <- lumen_xenium@meta.data %>%
  select(lumen_id, sample, sample_type, num_cells, num_cells_bin,
         max_dist, max_dist_bin) %>%
  unique()
rownames(lumen_metadata) <- lumen_metadata$lumen_id
ct_lumen_expr <- table(lumen_xenium$broad_CT5, lumen_xenium$lumen_id)
cniche_lumen_expr <- table(lumen_xenium$CNiche, lumen_xenium$lumen_id)
tniche_lumen_expr <- table(lumen_xenium$TNiche, lumen_xenium$lumen_id)

# Create cell composition Seurat object
all.equal(colnames(ct_lumen_expr), rownames(lumen_metadata[colnames(ct_lumen_expr), ]))
ct_lumen_obj <- CreateSeuratObject(Matrix::Matrix(ct_lumen_expr), 
                                   meta.data = lumen_metadata[colnames(ct_lumen_expr), ],
                                   assay = "CTcomp")
ct_lumen_obj[["cniche"]] <- CreateAssayObject(Matrix::Matrix(cniche_lumen_expr))
ct_lumen_obj[["tniche_cells"]] <- CreateAssayObject(Matrix::Matrix(tniche_lumen_expr))

# Filter object by lumen size
ct_lumen_obj_unfiltered <- ct_lumen_obj

filter2 <- subset(ct_lumen_obj_unfiltered, subset = num_cells >= 25 & num_cells <= 500 & max_dist >= 125)
filter_epi_prop4 <- subset(ct_lumen_obj_unfiltered, subset = num_cells >= 25 & num_cells <= 500 & max_dist >= 100)


## PCA ----
#### CELL TYPE COMPOSITION ----
# Filter out non-alveolar lumens
ct_lumen_obj <- subset(ct_lumen_obj_unfiltered, subset = lumen_id %in% filter_epi_prop4$lumen_id) # Was using 2

# Number of lumens before and after filtering
bf_lumens <- table(ct_lumen_obj_unfiltered$sample)
aft_lumens <- table(ct_lumen_obj$sample)

# PCA on cell type composition
cell_lumens_mat <- t(as.matrix(ct_lumen_obj@assays$CTcomp@counts))
# cell_lumens_mat <- cell_lumens_mat[, -which(colnames(cell_lumens_mat) == "Mesothelial")] # Remove mesothelial since there are very few
cell_lumens_prop_mat <- cell_lumens_mat/rowSums(cell_lumens_mat)

cniche_lumens_mat <-  t(as.matrix(ct_lumen_obj@assays$cniche@counts))
cniche_lumens_prop_mat <- cniche_lumens_mat/rowSums(cniche_lumens_mat)

tniche_lumens_mat <-  t(as.matrix(ct_lumen_obj@assays$tniche@counts))
tniche_lumens_prop_mat <- tniche_lumens_mat/rowSums(tniche_lumens_mat)

cniche_max <- cniche_lumens_prop_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "lumen_id") %>%
  pivot_longer(2:ncol(.), names_to = "CNiche", values_to = "Prop") %>%
  group_by(lumen_id) %>%
  mutate(max = max(Prop)) %>%
  pivot_wider(values_from = "Prop", names_from = "CNiche") %>%
  mutate(C2_C11 = C2 + C11) %>%
  filter(C2 < 0.15, C11 < 0.15, C2_C11 < 0.25) %>%
  select(-C2_C11) %>%
  pivot_longer(3:ncol(.), values_to = "Prop", names_to = "CNiche") %>%
  filter(Prop == max) %>%
  mutate(Max_CNiche = ifelse(sum(Prop) > max, "Multiple", CNiche)) %>%
  filter(Max_CNiche %in% c("C6", "C10", "C12")) %>%
  ungroup() %>%
  select(lumen_id, Max_CNiche) %>%
  unique()

cniche_lumens_prop_mat <- cniche_lumens_prop_mat[cniche_max$lumen_id, ]
cniche_max <- cniche_max[cniche_max$lumen_id == rownames(cniche_lumens_prop_mat), ]
cell_lumens_prop_mat <- cell_lumens_prop_mat[rownames(cniche_lumens_prop_mat), ]
cell_lumens_mat <- cell_lumens_mat[rownames(cniche_lumens_prop_mat), ]
cell_lumen_pca <- prcomp(cell_lumens_prop_mat, 
                         center = TRUE,
                         scale = TRUE
)
length(cniche_max$lumen_id)
# length(alv_lumen_labels)

# Save individual files with filtered lumen info
for (sm in sample_ids) {
  tmp <- new_tmp %>%
    filter(sample == sm) %>%
    left_join(max_dist_df) %>%
    ungroup() %>%
    select(cell_id, group) %>%
    filter(group %in% cniche_max$lumen_id)
  colnames(tmp) <- c("cell_id", "group")
  write.csv(tmp, paste0("/Volumes/dback_scratch/avannan/lumen_metadata_4explorer/", sm, "_retained_lumens_111523.csv"))
}

# Get metadata
cell_lumen_metadata <- ct_lumen_obj@meta.data %>%
  filter(lumen_id %in% rownames(cell_lumens_prop_mat))

# # Allow maximum repel of labels
options(ggrepel.max.overlaps = 20)

fviz_pca_var(cell_lumen_pca, col.var = "black", alpha.var = 0.2) # Cell types
fviz_pca_var(cell_lumen_pca, col.var = "black", alpha.var = 0.2, axes = c(2, 3)) # Cell types
fviz_pca_var(cell_lumen_pca, col.var = "black", alpha.var = 0.2, axes = c(3, 4)) # Cell types
fviz_pca_ind(cell_lumen_pca, col.var = "black") # Lumens
fviz_pca_var(cell_lumen_pca, col.var = "black", alpha.var = 0.2, palette = list(AT2 = "red") ) # Cell types

cniche_prop_pca_plots <- lapply(paste0("C", 1:12), function(XX) {
  fviz_pca_ind(cell_lumen_pca, col.ind = cniche_lumens_prop_mat[, XX],
               label = "none", geom = "point", alpha.ind = 0.9,
               gradient.cols = c("grey80", "#e77f00", "darkred"), axes = c(1, 2))
})
names(cniche_prop_pca_plots) <- paste0("C", 1:12)
(cniche_prop_pca_plots[["C12"]] + ggtitle("C12 - Healthy Alveoli"))
(cniche_prop_pca_plots[["C10"]] + ggtitle("C10 - Transitional Alveoli"))
(cniche_prop_pca_plots[["C6"]] + ggtitle("C6 - Macrophage Accumulation"))

cell_lumen_pca$x %>%
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column(var = "lumen_id") %>%
  full_join(cniche_lumens_prop_mat %>% as.data.frame() %>% rownames_to_column(var = "lumen_id")) %>%
  pivot_longer(C1:C12, names_to = "CNiche", values_to = "CNiche_Prop") %>%
  ggplot(aes(x = PC1, y = PC2, color = CNiche_Prop)) +
  geom_point(alpha = 0.8, shape = 16) +
  scale_color_continuous(type = "viridis") +
  theme_bw() +
  facet_wrap(~CNiche)

fviz_pca_ind(cell_lumen_pca, col.ind = cniche_max$Max_CNiche, geom = "point",
             palette = nuclei_niche_color_list,
             mean.point.size = 5,
             axes = c(1, 2), alpha = 0.4) + scale_shape_manual(values = rep(16, 200))

fviz_pca_ind(cell_lumen_pca, col.ind = cell_lumen_metadata$sample_type, geom = "point",
             # palette = nuclei_niche_color_list,
             mean.point.size = 5,
             axes = c(1, 2), alpha = 0.4) + scale_shape_manual(values = rep(16, 200))

# By sample colors
nuclei_niche_color_list <- list(`Multiple` = "grey",
                                `C1` = "#ffa26d",
                                `C2` = "#1036b5",
                                `C3` = "#cabd00",
                                `C4` = "#a6513c",
                                `C5` = "#84cd5d",
                                `C6` = "#72c7ff",
                                `C7` = "#d10040",
                                `C8` = "black",
                                `C9` = "#71c8a5",
                                `C10` = "#6f774b",
                                `C11` = "#8E76DB",
                                `C12` = "#924373")


#### SLINGSHOT COMPOSITION ----
comp_sce <- SingleCellExperiment(assays = List(prop = t(cell_lumens_prop_mat),
                                               counts = t(cell_lumens_mat)))

cell_lumen_pca <- prcomp(cell_lumens_prop_mat, 
                         center = TRUE,
                         scale = TRUE
)
comp_rd1 <- cell_lumen_pca$x[, 1:2]
reducedDims(comp_sce) <- SimpleList(PCA = comp_rd1)

# Try different levels of clustering
library(mclust, quietly = TRUE)
comp_kmeans3 <- kmeans(comp_rd1, centers = 3)$cluster
comp_kmeans4 <- kmeans(comp_rd1, centers = 4)$cluster
comp_kmeans5 <- kmeans(comp_rd1, centers = 5)$cluster
comp_kmeans6 <- kmeans(comp_rd1, centers = 6)$cluster
colData(comp_sce)$comp_kmeans3 <- comp_kmeans3
colData(comp_sce)$comp_kmeans4 <- comp_kmeans4
colData(comp_sce)$comp_kmeans5 <- comp_kmeans5
colData(comp_sce)$comp_kmeans6 <- comp_kmeans6

# Add metadata
colData(comp_sce)$sample <- ct_lumen_obj@meta.data[colnames(comp_sce), ]$sample
colData(comp_sce)$sample_type <- ct_lumen_obj@meta.data[colnames(comp_sce), ]$sample_type

colData(comp_sce)$cniche1_prop <- cniche_lumens_prop_mat[, "C1"]
colData(comp_sce)$cniche2_prop <- cniche_lumens_prop_mat[, "C2"]
colData(comp_sce)$cniche3_prop <- cniche_lumens_prop_mat[, "C3"]
colData(comp_sce)$cniche4_prop <- cniche_lumens_prop_mat[, "C4"]
colData(comp_sce)$cniche5_prop <- cniche_lumens_prop_mat[, "C5"]
colData(comp_sce)$cniche6_prop <- cniche_lumens_prop_mat[, "C6"]
colData(comp_sce)$cniche7_prop <- cniche_lumens_prop_mat[, "C7"]
colData(comp_sce)$cniche8_prop <- cniche_lumens_prop_mat[, "C8"]
colData(comp_sce)$cniche9_prop <- cniche_lumens_prop_mat[, "C9"]
colData(comp_sce)$cniche10_prop <- cniche_lumens_prop_mat[, "C10"]
colData(comp_sce)$cniche11_prop <- cniche_lumens_prop_mat[, "C11"]
colData(comp_sce)$cniche12_prop <- cniche_lumens_prop_mat[, "C12"]
colData(comp_sce)$cniche_max <- as.factor(cniche_max$Max_CNiche)

# Plot clusters
library(RColorBrewer)
# 1 = red, 2 = blue, 3 = green,  4 = purple, 5 = orange, 6 = yellow
plot(comp_rd1, col = brewer.pal(9, "Set1")[comp_sce$sample_type], pch = 16, asp = 1)
plot(comp_rd1, col = brewer.pal(9,"Set1")[comp_kmeans3], pch=16, asp = 1) # Start is 3
plot(comp_rd1, col = brewer.pal(9,"Set1")[comp_kmeans4], pch=16, asp = 1) # Start is 2
plot(comp_rd1, col = brewer.pal(9,"Set1")[comp_kmeans5], pch=16, asp = 1) # Start is 3
plot(comp_rd1, col = brewer.pal(9,"Set1")[comp_kmeans6], pch=16, asp = 1) # Start is 3

library(grDevices)
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)

# View pseudotime
comp_sce_k3 <- slingshot(comp_sce, clusterLabels = 'comp_kmeans3', reducedDim = 'PCA', start.clus = "3") # MAKE SURE START IS CORRECT
summary(comp_sce_k3$slingPseudotime_1)
plot(reducedDims(comp_sce_k3)$PCA, col = brewer.pal(9,'Set1')[comp_sce_k3$comp_kmeans3], pch=16, asp = 1)
lines(SlingshotDataSet(comp_sce_k3), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(comp_sce_k3)$PCA, col = colors[cut(comp_sce_k3$slingPseudotime_1, breaks=100)], pch=16, asp = 1) + title("PS1")
lines(SlingshotDataSet(comp_sce_k3), lwd=2, col='black')

comp_sce_k4 <- slingshot(comp_sce, clusterLabels = 'comp_kmeans4', reducedDim = 'PCA', start.clus = "2") # MAKE SURE START IS CORRECT
summary(comp_sce_k4$slingPseudotime_1)
plot(reducedDims(comp_sce_k4)$PCA, col = brewer.pal(9,'Set1')[comp_sce_k4$comp_kmeans4], pch=16, asp = 1)
lines(SlingshotDataSet(comp_sce_k4), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(comp_sce_k4)$PCA, col = colors[cut(comp_sce_k4$slingPseudotime_1, breaks=100)], pch=16, asp = 1) + title("PS1")
lines(SlingshotDataSet(comp_sce_k4), lwd=2, col='black')

comp_sce_k5 <- slingshot(comp_sce, clusterLabels = 'comp_kmeans5', reducedDim = 'PCA', start.clus = "3") # MAKE SURE START IS CORRECT
summary(comp_sce_k5$slingPseudotime_1)
plot(reducedDims(comp_sce_k5)$PCA, col = brewer.pal(9,'Set1')[comp_sce_k5$comp_kmeans5], pch=16, asp = 1)
lines(SlingshotDataSet(comp_sce_k5), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(comp_sce_k5)$PCA, col = colors[cut(comp_sce_k5$slingPseudotime_1, breaks=100)], pch=16, asp = 1) + title("PS1")
lines(SlingshotDataSet(comp_sce_k5), lwd=2, col='black')

comp_sce_k6 <- slingshot(comp_sce, clusterLabels = 'comp_kmeans6', reducedDim = 'PCA', start.clus = "3") # MAKE SURE START IS CORRECT
summary(comp_sce_k6$slingPseudotime_1)
plot(reducedDims(comp_sce_k6)$PCA, col = brewer.pal(9,'Set1')[comp_sce_k6$comp_kmeans6], pch=16, asp = 1)
lines(SlingshotDataSet(comp_sce_k6), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(comp_sce_k6)$PCA, col = colors[cut(comp_sce_k6$slingPseudotime_1, breaks=100)], pch=16, asp = 1) + title("PS1")
lines(SlingshotDataSet(comp_sce_k6), lwd=2, col='black')

transcript_niche_color_list <- list(`Multiple` = "grey",
                                    `T0` = "#24504f",
                                    `T1` = "black",
                                    `T2` = "#008200",
                                    `T3` = "#530087",
                                    `T4` = "#ff0000",
                                    `T5` = "#ffd500",
                                    `T6` = "#b3bd51",
                                    `T7` = "#00ffff",
                                    `T8` = "#0000ff",
                                    `T9` = "#ff5bb6")


#### EXPRESSION ----
# Filtered lumens
filtered_lumens <- colnames(comp_sce)
filtered_lumens_xenium <- subset(lumen_xenium, subset = lumen_id %in% filtered_lumens & lumen_id != "VUILD91LF_L162")

# Composite score of transcript expression overall within each lumen
i = 1
for (lumen in filtered_lumens) { 
  # Get cells in the lumen
  lumen_metadata <- filtered_lumens_xenium@meta.data %>%
    filter(lumen_id == lumen) %>%
    select(lumen_id, sample, sample_type)
  cells_in_lumen <- rownames(lumen_metadata)
  
  # Get counts data for cells in lumen
  counts <- filtered_lumens_xenium@assays$RNA@counts[, cells_in_lumen]
  
  # Sum counts across lumen
  lumen_counts <- as.matrix(data.frame(rowSums(counts)))
  colnames(lumen_counts) <- lumen
  
  if (i == 1) {
    # Create matrix
    lumen_counts_mat <- lumen_counts
    
    # Create metadata
    full_lumen_metadata <- lumen_metadata %>% unique()

  } else {
    # Add to matrix
    lumen_counts_mat <- cbind(lumen_counts_mat, lumen_counts)
    
    # Add to metadata
    full_lumen_metadata <- rbind(full_lumen_metadata, lumen_metadata %>% unique())
  }
  i = i + 1
}
lumen_counts_mat %>% dim
full_lumen_metadata %>% dim

# Set up expression matrix
lumen_counts_mat2 <- lumen_counts_mat
full_lumen_metadata2 <- full_lumen_metadata %>%
  filter(lumen_id %in% colnames(lumen_counts_mat2)) %>%
  full_join(cell_lumen_metadata %>% select(lumen_id, nCount_CTcomp))
lumen_pca <- prcomp(t(lumen_counts_mat2), center = TRUE, scale = TRUE)

fviz_pca_var(lumen_pca, col.var = "black", alpha.var = 0.2) # Genes
fviz_pca_ind(lumen_pca, col.var = "black") # Samples
fviz_pca_var(lumen_pca, col.var = "black", axes = c(2,3), alpha.var = 0.2) # Genes

fviz_pca_ind(lumen_pca, habillage = full_lumen_metadata2$sample_type, label = "none", geom = "point",
             palette = c("green", "blue", "red", "grey"), addEllipses = TRUE,
             axes = c(1,2)) + scale_shape_manual(values = c(16, 16, 16, 16)) # Samples by sample type
fviz_pca_ind(lumen_pca, habillage = full_lumen_metadata2$sample_type, label = "none", geom = "point",
             palette = c("green", "blue", "red", "grey"), addEllipses = TRUE,
             axes = c(2,3)) + scale_shape_manual(values = c(16, 16, 16, 16)) # Samples by sample type
fviz_pca_var(lumen_pca, axes = c(2,3), col.var = "black") # Genes
 
fviz_pca_ind(lumen_pca, col.ind = full_lumen_metadata2$nCount_CTcomp, label = "none", geom = "point",
             axes = c(1, 2), gradient.cols = c("white", "pink", "red")) + scale_shape_manual(values = rep(16, 266))
fviz_pca_ind(lumen_pca, col.ind = full_lumen_metadata2$nCount_CTcomp, label = "none", geom = "point",
             axes = c(2, 3), gradient.cols = c("white", "pink", "red")) + scale_shape_manual(values = rep(16, 266))
fviz_pca_ind(lumen_pca, habillage = full_lumen_metadata2$nCount_CTcomp, label = "none", geom = "point",
             axes = c(2, 3), palette = viridis::viridis(266)) + scale_shape_manual(values = rep(16, 266)) + NoLegend()

# saveRDS(expr_lumen_norm_pca, "/Volumes/dback_scratch/avannan/expr_lumen_norm_pca.rds")
# saveRDS(expr_lumen_count_pca, "/Volumes/dback_scratch/avannan/expr_lumen_count_pca.rds")
# saveRDS(expr_lumen_metadata, "/Volumes/dback_scratch/avannan/expr_lumen_metadata.rds")
# saveRDS(cell_lumen_pca, "/Volumes/dback_scratch/avannan/cell_lumen_pca.rds")
# saveRDS(cell_lumen_metadata, "/Volumes/dback_scratch/avannan/cell_lumen_metadata.rds")


## SLINGSHOT ----
expr_seurat <- CreateSeuratObject(Matrix::Matrix(lumen_counts_mat2), assay = "lumen_RNA")
expr_seurat <- NormalizeData(expr_seurat)
expr_seurat <- ScaleData(expr_seurat, vars.to.regress = "nCount_lumen_RNA")

expr_counts <- expr_seurat@assays$lumen_RNA@layers$counts
rownames(expr_counts) <- rownames(expr_seurat)
colnames(expr_counts) <- colnames(expr_seurat)

expr_data <- expr_seurat@assays$lumen_RNA@layers$data
rownames(expr_data) <- rownames(expr_seurat)
colnames(expr_data) <- colnames(expr_seurat)

expr_scale_data <- expr_seurat@assays$lumen_RNA@layers$scale.data
rownames(expr_scale_data) <- rownames(expr_seurat)
colnames(expr_scale_data) <- colnames(expr_seurat)

expr_sce <- SingleCellExperiment(assays = List(counts = expr_counts,
                                               data = expr_data))

lumen_pca <- prcomp(t(expr_scale_data), center = FALSE, scale = FALSE) # Already scaled and centered
fviz_pca_var(lumen_pca, col.var = "black", alpha.var = 0.2) # Genes
fviz_pca_ind(lumen_pca, col.var = "black") # Samples
fviz_pca_var(lumen_pca, col.var = "black", axes = c(2,3), alpha.var = 0.2) # Genes

fviz_pca_ind(lumen_pca, habillage = full_lumen_metadata2$sample_type, label = "none", geom = "point",
             palette = c("green", "blue", "red", "grey"), addEllipses = FALSE,
             axes = c(1,2)) + scale_shape_manual(values = c(16, 16, 16, 16)) # Samples by sample type
fviz_pca_ind(lumen_pca, habillage = full_lumen_metadata2$sample_type, label = "none", geom = "point",
             palette = c("green", "blue", "red", "grey"), addEllipses = TRUE,
             axes = c(2,3)) + scale_shape_manual(values = c(16, 16, 16, 16)) # Samples by sample type
fviz_pca_var(lumen_pca, axes = c(2,3), col.var = "black") # Genes

fviz_pca_ind(lumen_pca, habillage = full_lumen_metadata2$nCount_CTcomp, label = "none", geom = "point",
             axes = c(1, 2), palette = viridis::viridis(221)) + scale_shape_manual(values = rep(16, 221)) + NoLegend()
fviz_pca_ind(lumen_pca, habillage = expr_seurat$nCount_lumen_RNA, label = "none", geom = "point",
             axes = c(1, 2), palette = viridis::viridis(1020)) + scale_shape_manual(values = rep(16, 1020)) + NoLegend()
fviz_pca_ind(lumen_pca, habillage = full_lumen_metadata2$nCount_CTcomp, label = "none", geom = "point",
             axes = c(2, 3), palette = viridis::viridis(221)) + scale_shape_manual(values = rep(16, 221)) + NoLegend()
fviz_pca_ind(lumen_pca, habillage = expr_seurat$nCount_lumen_RNA, label = "none", geom = "point",
             axes = c(2, 3), palette = viridis::viridis(1020)) + scale_shape_manual(values = rep(16, 1020)) + NoLegend()

cor(as.data.frame(cell_lumen_pca$x)$PC1, as.data.frame(lumen_pca$x)$PC2, method = "spearman")

# Actual slingshot
expr_rd1 <- lumen_pca$x[, 2:3] # PCs 2 and 3
reducedDims(expr_sce) <- SimpleList(PCA = expr_rd1)

# Try different levels of clustering
library(mclust, quietly = TRUE)
expr_mclust <- Mclust(expr_rd1)$classification
colData(expr_sce)$GMM <- expr_mclust
expr_kmeans3 <- kmeans(expr_rd1, centers = 3)$cluster
expr_kmeans4 <- kmeans(expr_rd1, centers = 4)$cluster
expr_kmeans5 <- kmeans(expr_rd1, centers = 5)$cluster
expr_kmeans6 <- kmeans(expr_rd1, centers = 6)$cluster
colData(expr_sce)$expr_kmeans3 <- expr_kmeans3
colData(expr_sce)$expr_kmeans4 <- expr_kmeans4
colData(expr_sce)$expr_kmeans5 <- expr_kmeans5
colData(expr_sce)$expr_kmeans6 <- expr_kmeans6

colData(expr_sce)$cniche1_prop <- cniche_lumens_prop_mat[, "C1"]
colData(expr_sce)$cniche2_prop <- cniche_lumens_prop_mat[, "C2"]
colData(expr_sce)$cniche3_prop <- cniche_lumens_prop_mat[, "C3"]
colData(expr_sce)$cniche4_prop <- cniche_lumens_prop_mat[, "C4"]
colData(expr_sce)$cniche5_prop <- cniche_lumens_prop_mat[, "C5"]
colData(expr_sce)$cniche6_prop <- cniche_lumens_prop_mat[, "C6"]
colData(expr_sce)$cniche7_prop <- cniche_lumens_prop_mat[, "C7"]
colData(expr_sce)$cniche8_prop <- cniche_lumens_prop_mat[, "C8"]
colData(expr_sce)$cniche9_prop <- cniche_lumens_prop_mat[, "C9"]
colData(expr_sce)$cniche10_prop <- cniche_lumens_prop_mat[, "C10"]
colData(expr_sce)$cniche11_prop <- cniche_lumens_prop_mat[, "C11"]
colData(expr_sce)$cniche12_prop <- cniche_lumens_prop_mat[, "C12"]
colData(expr_sce)$cniche_max <- as.factor(cniche_max$Max_CNiche)


# Add metadata
colData(expr_sce)$sample <- ct_lumen_obj@meta.data[colnames(expr_sce), ]$sample
colData(expr_sce)$sample_type <- ct_lumen_obj@meta.data[colnames(expr_sce), ]$sample_type

# Plot clusters
library(RColorBrewer)
plot(expr_rd1, col = brewer.pal(9, "Set1")[expr_sce$sample_type], pch = 16, asp = 1)
plot(expr_rd1, col = brewer.pal(9, "Set1")[expr_mclust], pch = 16, asp = 1) # Start is 2
plot(expr_rd1, col = brewer.pal(9,"Set1")[expr_kmeans3], pch=16, asp = 1) # Start is 1
plot(expr_rd1, col = brewer.pal(9,"Set1")[expr_kmeans4], pch=16, asp = 1) # Start is 4
plot(expr_rd1, col = brewer.pal(9,"Set1")[expr_kmeans5], pch=16, asp = 1) # Start is 1
plot(expr_rd1, col = brewer.pal(9,"Set1")[expr_kmeans6], pch=16, asp = 1) # Start is 3

library(grDevices)
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)

# View GMM by sample type coloring
expr_sce_gmm <- slingshot(expr_sce, clusterLabels = 'GMM', reducedDim = 'PCA')
summary(expr_sce_gmm$slingPseudotime_1)
plot(reducedDims(expr_sce_gmm)$PCA, col = brewer.pal(9,'Set1')[expr_sce_gmm$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(expr_sce_gmm), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(expr_sce_gmm)$PCA, col = brewer.pal(9,'Set1')[expr_sce_gmm$sample_type], pch=16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(expr_sce_gmm), lwd=1, type = 'lineages', col = 'black')
# Red = unaffected, blue = LF, green = MF, purple = ILD
# Color order: red, blue, green, purple, orange, yellow, brown, pink, grey


fviz_pca_ind(lumen_pca, col.ind = cniche_max$Max_CNiche, geom = "point",
             palette = nuclei_niche_color_list,
             mean.point.size = 5,
             axes = c(2, 23), alpha = 0.4) + scale_shape_manual(values = rep(16, 200))

expr_sce_st <- slingshot(expr_sce, clusterLabels = 'sample_type', reducedDim = 'PCA',
                         start.clus = "Unaffected")
plot(reducedDims(expr_sce_st)$PCA, col = colors[cut(expr_sce_st$slingPseudotime_1, breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(expr_sce_st), lwd=2, col='black')

# View pseudotime
expr_sce_k3 <- slingshot(expr_sce, clusterLabels = 'expr_kmeans3', reducedDim = 'PCA', start.clus = "1") # MAKE SURE START IS CORRECT
summary(expr_sce_k3$slingPseudotime_1)
plot(reducedDims(expr_sce_k3)$PCA, col = brewer.pal(9,'Set1')[expr_sce_k3$expr_kmeans3], pch=16, asp = 1)
lines(SlingshotDataSet(expr_sce_k3), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(expr_sce_k3)$PCA, col = colors[cut(expr_sce_k3$slingPseudotime_1, breaks=100)], pch=16, asp = 1) + title("PS1")
lines(SlingshotDataSet(expr_sce_k3), lwd=2, col='black')

expr_sce_k4 <- slingshot(expr_sce, clusterLabels = 'expr_kmeans4', reducedDim = 'PCA', start.clus = "2") # MAKE SURE START IS CORRECT
summary(expr_sce_k4$slingPseudotime_1)
plot(reducedDims(expr_sce_k4)$PCA, col = brewer.pal(9,'Set1')[expr_sce_k4$expr_kmeans4], pch=16, asp = 1)
lines(SlingshotDataSet(expr_sce_k4), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(expr_sce_k4)$PCA, col = colors[cut(expr_sce_k4$slingPseudotime_1, breaks=100)], pch=16, asp = 1) + title("PS1")
lines(SlingshotDataSet(expr_sce_k4), lwd=2, col='black')
plot(reducedDims(expr_sce_k4)$PCA, col = colors[cut(expr_sce_k4$slingPseudotime_2, breaks=100)], pch=16, asp = 1) + title("PS1")
lines(SlingshotDataSet(expr_sce_k4), lwd=2, col='black')

expr_sce_k5 <- slingshot(expr_sce, clusterLabels = 'expr_kmeans5', reducedDim = 'PCA', start.clus = "3") # MAKE SURE START IS CORRECT
summary(expr_sce_k5$slingPseudotime_1)
plot(reducedDims(expr_sce_k5)$PCA, col = brewer.pal(9,'Set1')[expr_sce_k5$expr_kmeans5], pch=16, asp = 1)
lines(SlingshotDataSet(expr_sce_k5), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(expr_sce_k5)$PCA, col = colors[cut(expr_sce_k5$slingPseudotime_1, breaks=100)], pch=16, asp = 1) + title("PS1")
lines(SlingshotDataSet(expr_sce_k5), lwd=2, col='black')

expr_sce_k6 <- slingshot(expr_sce, clusterLabels = 'expr_kmeans6', reducedDim = 'PCA', start.clus = "3") # MAKE SURE START IS CORRECT
summary(expr_sce_k6$slingPseudotime_1)
plot(reducedDims(expr_sce_k6)$PCA, col = brewer.pal(9,'Set1')[expr_sce_k6$expr_kmeans6], pch=16, asp = 1)
lines(SlingshotDataSet(expr_sce_k6), lwd=1, type = 'lineages', col = 'black')
plot(reducedDims(expr_sce_k6)$PCA, col = colors[cut(expr_sce_k6$slingPseudotime_1, breaks=100)], pch=16, asp = 1) + title("PS1")










