###########################################
# Initial Lumen Filtering and PCA
# Author: Annika Vannan (avannan@tgen.org)
# Date: 07/18/2024
###########################################

## SETTING ENVIRONMENT ----
library(ggplot2, lib.loc = "/usr/local/lib/R/site-library")
library(scCustomize)
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



## SET COLORS ----
color_list <- list(# Endothelial
  `Arteriole` = "#e46300",
  `Capillary` = "#cdbf5e",
  `Lymphatic` = "#975f47",
  `Venous` = "#d59445",
  
  # Epithelial
  `AT1` = "#6cd66d",
  `AT2` = "#34491d",
  `Basal` = "#31a200",
  `Multiciliated` = "#b9b700",
  `KRT5-/KRT17+` = "#a0cd9e",
  `Goblet` = "#acde47",
  `RASC` = "#ccd87a",
  `Secretory` = "#01eaaa",
  `Transitional AT2` = "#4fa381",
  `PNEC` = "#aef359",
  `Proliferating Airway` = "#006e3e",
  `Proliferating AT2` = "#5e6b00",
  
  # Immune (Lymphoid)
  `CD4+ T-cells` = "#e62f4c",
  `CD8+ T-cells` = "#dc6fd5",
  `Tregs` = "#9e548d",
  `B cells` = "#d84265",
  `NK/NKT` = "#b16f6f",
  `Plasma` = "#df4738",
  `pDCs` = "#842b8c",
  `Proliferating T-cells` = "#993432",
  `Proliferating NK/NKT` = "#ab3b99",
  `Proliferating B cells` = "#724173",
  
  # Immune (Myeloid)
  `Alveolar Macrophages` = "#ffa3a1",
  `Basophils` = "#d93e77",
  `cDCs` =  "#cb547b",
  # `cDCs/Langerhans cells` =  "#cb547b",
  `Langerhans cells` = "#e5849a",
  `Interstitial Macrophages` = "#dc48bd",
  `Macrophages - IFN-activated` = "#963095",
  `Mast` = "#bf4ce3",
  `Migratory DCs` = "#a94b59",
  `Monocytes/MDMs`  = "#f2294e",
  `Neutrophils` = "#a51219",
  `SPP1+ Macrophages` = "#9e3c6d",
  `Proliferating Myeloid` = "#e25292",
  
  # Mesenchymal
  `Activated Fibrotic FBs` = "#005a90",
  `Adventitial FBs` = "#00cace",
  `Alveolar FBs` = "#6d4bce",
  `Inflammatory FBs` = "#7895ff",
  `Mesothelial` = "#374ca3",
  `Myofibroblasts` = "#4c4e8c",
  `SMCs/Pericytes` = "#00aec9",
  `Subpleural FBs` = "#0065f5",
  `Proliferating FBs` = "#7dfff6")

nuclei_niche_color_list <- list(`C7` = "#ffa26d",
                                `C1` = "#1036b5",
                                `C12` = "#cabd00",
                                `C9` = "#a6513c",
                                `C4` = "#84cd5d",
                                `C11` = "#72c7ff",
                                `C3` = "#d10040",
                                `C10` = "black",
                                `C6` = "#71c8a5",
                                `C2` = "#6f774b",
                                `C5` = "#8E76DB",
                                `C8` = "#924373",
                                `Multiple` = "grey80")


transcript_niche_color_list <- list(`T1` = "#94FFB5", 
                                    `T2` = "#191919",
                                    `T3` = "#8F7C00",
                                    `T4` = "#FFCC99",
                                    `T5` = "#2BCE48",
                                    `T6` = "#993F00",
                                    `T7` = "#005C31",
                                    `T8` = "#0075DC", 
                                    `T9` = "#F0A0FF", 
                                    `T10` = "#C20088", 
                                    `T11` = "#4C005C", 
                                    `T12` = "#9DCC00",
                                    `Multiple` = "grey80")


## LOAD DATA AND LUMEN INFORMATION ----
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_070924.rds")
pretty_umap <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     axis.title = element_text(hjust = 1))

file_list <- list.files("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_segmentation_june_2024_revision2", full.names = TRUE)
names(file_list) <- unlist(lapply(str_split(file_list, "_|/"), function(XX) XX[14]))
names(file_list)[9] <- "TILD117MFB"

df_list <- lapply(names(file_list), function(XX) {
  if (XX == "TILD117MFB") {
    read_csv(file_list[[XX]]) %>%
      select(cell_id, lumen_id) %>%
      separate(cell_id, into = c("temp", "cell_id"), sep = "_") %>%
      unique() %>%
      mutate(sample = XX,
             full_cell_id = paste0(sample, "_", cell_id)) %>%
      select(sample, full_cell_id, lumen_id) %>%
      mutate(lumen_id = ifelse(!is.na(lumen_id) & lumen_id != 0, 
                               paste0(sample, "_", lumen_id), NA_character_))
  } else {
    read_csv(file_list[[XX]]) %>%
      select(cell_id, lumen_id) %>%
      unique() %>%
      mutate(sample = XX) %>%
      dplyr::rename("full_cell_id" = "cell_id") %>%
      select(sample, full_cell_id, lumen_id) %>%
      mutate(lumen_id = ifelse(!is.na(lumen_id) & lumen_id != 0, 
                               paste0(sample, "_", lumen_id), NA_character_))
  }
})
names(df_list) <- names(file_list)
df <- Reduce(rbind, df_list) %>%
  filter(lumen_id != 0, !is.na(lumen_id), full_cell_id %in% colnames(xenium))


## ADD DATA TO SEURAT OBJECT ----
lumen_meta <- xenium@meta.data %>%
  rownames_to_column(var = "full_cell_id") %>%
  select(sample, full_cell_id) %>%
  left_join(df)
all.equal(lumen_meta$full_cell_id, colnames(xenium)) # Must be TRUE
xenium$lumen_id <- lumen_meta$lumen_id

# Save object
# saveRDS(xenium, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_lumens_071624.rds")
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_lumens_071624.rds")


## SET UP LUMEN OBJECT ----
# Only keep cells that are assigned to lumens
cells_with_lumens <- rownames(xenium@meta.data %>% filter(!is.na(lumen_id), lumen_id != 0))
lumen_xenium <- subset(xenium, cells = cells_with_lumens)

# Add metadata about number of cells
lumen_xenium@meta.data <- lumen_xenium@meta.data %>%
  group_by(lumen_id) %>%
  mutate(num_cells_lumen = length(lumen_id))


## CALCULATE MAX DISTANCE BETWEEN CELLS ----
# Make temporary dataframe for calculating max distances between cells in each lumen
mini_metadata <- xenium@meta.data %>%
  rownames_to_column(var = "full_cell_id") %>%
  select(sample, full_cell_id, x_centroid, y_centroid, final_CT, final_lineage, lumen_id)

# Calculate max distance between cells in each lumen
# Do very basic filtering beforehand
sample_ids <- unique(xenium$sample)
max_dist_list <- list()
for (sm in sample_ids) {
  message(paste("Processing", sm))
  sm_tmp <- mini_metadata %>% 
    filter(sample == sm) %>%
    select(full_cell_id, lumen_id, final_lineage) %>%
    group_by(lumen_id)  %>%
    mutate(num_cells = length(full_cell_id)) %>%
    ungroup() %>%
    filter(num_cells < 1000, num_cells >= 10)
  
  sm_tmp1 <- sm_tmp %>%
    left_join(mini_metadata %>% filter(sample == sm))
  
  sm_tmp2 <- sm_tmp %>%
    group_by(lumen_id, final_lineage) %>%
    mutate(num_cells_lineage = length(final_lineage)) %>%
    ungroup() %>%
    select(-full_cell_id) %>%
    unique() %>%
    mutate(final_lineage = str_to_lower(paste0("num_cells_", final_lineage))) %>%
    pivot_wider(names_from = "final_lineage", values_from = "num_cells_lineage")
  
  sm_tmp <- full_join(sm_tmp1, sm_tmp2)
  
  for (lumen in unique(sm_tmp$lumen_id)) {
    lumen_tmp <- sm_tmp %>% filter(lumen_id == lumen)
    calc_distances <- RANN::nn2(cbind(lumen_tmp$x_centroid,
                                      lumen_tmp$y_centroid),
                                k = nrow(lumen_tmp))
    max_dist_list[[sm]][[lumen]] <- lumen_tmp %>%
      left_join(data.frame(lumen_id = lumen,
                           max_dist = max(calc_distances$nn.dists)))
  }
}
# Turn list into a dataframe and bin the data
max_dist_df_list <- lapply(max_dist_list, function(XX) { Reduce(full_join, XX) })
max_dist_df <- Reduce(full_join, max_dist_df_list) %>% unique()
max_dist_df2 <- max_dist_df %>% 
  select(sample, lumen_id, max_dist, num_cells) %>% 
  unique() %>%
  mutate(num_cells_bin = cut(num_cells, breaks = seq(0, 1000, 25)),
         max_dist_bin = cut(max_dist, breaks = seq(0, 5000, 25)))

# Modify dataframe so that it can be added to Seurat object
dist_ncell_metadata_list <- lapply(sample_ids, function(XX) {
  max_dist_df2 %>%
    filter(sample == XX) %>%
    left_join(max_dist_df) %>%
    ungroup()
})
dist_ncell_metadata_df <- as.data.frame(Reduce(rbind, dist_ncell_metadata_list))
rownames(dist_ncell_metadata_df) <- dist_ncell_metadata_df$full_cell_id

# write.csv(dist_ncell_metadata_df, "/scratch/avannan/tmp_dist_ncell_metadata_df_070124.csv")
# write.csv(dist_ncell_metadata_df, "/scratch/avannan/tmp_dist_ncell_metadata_df_071624.csv")
# dist_ncell_metadata_df <- read_csv("/scratch/avannan/tmp_dist_ncell_metadata_df_071624.csv")

## SET UP CELL TYPE COMPOSITION DATA ----
# Get metadata and cell type "expression" data
lumen_metadata <- lumen_xenium@meta.data %>%
  full_join(dist_ncell_metadata_df) %>%
  select(lumen_id, sample, sample_type, sample_affect, num_cells_lumen, num_cells_bin, max_dist, max_dist_bin,
         num_cells_epithelial, num_cells_endothelial, num_cells_mesenchymal, num_cells_immune) %>%
  unique() %>%
  as.data.frame()
rownames(lumen_metadata) <- lumen_metadata$lumen_id

# Add more metadata
lumen_xenium$max_dist <- dist_ncell_metadata_df$max_dist[match(lumen_xenium$lumen_id, dist_ncell_metadata_df$lumen_id)]
lumen_xenium$num_cells_bin <- dist_ncell_metadata_df$num_cells_bin[match(lumen_xenium$lumen_id, dist_ncell_metadata_df$lumen_id)]
lumen_xenium$max_dist_bin <- dist_ncell_metadata_df$max_dist_bin[match(lumen_xenium$lumen_id, dist_ncell_metadata_df$lumen_id)]

ct_lumen_expr <- table(lumen_xenium$final_CT, lumen_xenium$lumen_id)
cniche_lumen_expr <- table(lumen_xenium$CNiche, lumen_xenium$lumen_id)
tniche_lumen_expr <- table(lumen_xenium$TNiche, lumen_xenium$lumen_id)

# Create cell composition Seurat object
all.equal(colnames(ct_lumen_expr), rownames(lumen_metadata[colnames(ct_lumen_expr), ]))
ct_lumen_obj_unfiltered <- CreateSeuratObject(Matrix::Matrix(ct_lumen_expr), 
                                              meta.data = lumen_metadata[colnames(ct_lumen_expr), ],
                                              assay = "CTcomp")
ct_lumen_obj_unfiltered[["cniche"]] <- CreateAssayObject(Matrix::Matrix(cniche_lumen_expr))
ct_lumen_obj_unfiltered[["tniche_cells"]] <- CreateAssayObject(Matrix::Matrix(tniche_lumen_expr))


## FILTER LUMENS ----
# saveRDS(ct_lumen_obj_unfiltered, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/ct_lumen_obj_unfiltered_071624.rds")
ct_lumen_obj_unfiltered <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/ct_lumen_obj_unfiltered_071624.rds")

ct_lumen_obj <- ct_lumen_obj_unfiltered
bf_lumens <- table(ct_lumen_obj_unfiltered$sample); sum(bf_lumens) # 5591 lumens
aft_lumens <- table(ct_lumen_obj$sample); sum(aft_lumens) # 3354 in TMAs 1-4; 5591 total

# Set up proportion data
cell_lumens_mat <- t(as.matrix(ct_lumen_obj@assays$CTcomp@layers$counts))
colnames(cell_lumens_mat) <- rownames(ct_lumen_obj)
rownames(cell_lumens_mat) <- colnames(ct_lumen_obj)
cell_lumens_prop_mat <- cell_lumens_mat/rowSums(cell_lumens_mat)

cniche_lumens_mat <-  t(as.matrix(ct_lumen_obj@assays$cniche@counts))
cniche_lumens_prop_mat <- cniche_lumens_mat/rowSums(cniche_lumens_mat)

tniche_lumens_mat <-  t(as.matrix(ct_lumen_obj@assays$tniche@counts))
tniche_lumens_prop_mat <- tniche_lumens_mat/rowSums(tniche_lumens_mat)

# Look at CNiche max
cniche_max <- cniche_lumens_prop_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "lumen_id") %>%
  pivot_longer(2:ncol(.), names_to = "CNiche", values_to = "Prop") %>%
  group_by(lumen_id) %>%
  mutate(max = max(Prop)) %>%
  pivot_wider(values_from = "Prop", names_from = "CNiche") %>%
  pivot_longer(3:ncol(.), values_to = "Prop", names_to = "CNiche") %>%
  filter(Prop == max) %>%
  mutate(Max_CNiche = ifelse(sum(Prop) > max, "Multiple", CNiche)) %>%
  ungroup() %>%
  select(lumen_id, Max_CNiche) %>%
  unique()

# Look at TNiche max
tniche_max <- tniche_lumens_prop_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "lumen_id") %>%
  pivot_longer(2:ncol(.), names_to = "TNiche", values_to = "Prop") %>%
  group_by(lumen_id) %>%
  mutate(max = max(Prop)) %>%
  pivot_wider(values_from = "Prop", names_from = "TNiche") %>%
  pivot_longer(3:ncol(.), values_to = "Prop", names_to = "TNiche") %>%
  filter(Prop == max) %>%
  mutate(Max_TNiche = ifelse(sum(Prop) > max, "Multiple", TNiche)) %>%
  ungroup() %>%
  select(lumen_id, Max_TNiche) %>%
  unique()

# Add information to Seurat object
add_metadata <- cniche_max %>%
  full_join(tniche_max) %>%
  full_join(as.data.frame(cniche_lumens_prop_mat) %>% rownames_to_column(var = "lumen_id")) %>%
  full_join(as.data.frame(tniche_lumens_prop_mat) %>% rownames_to_column(var = "lumen_id"))
all.equal(colnames(ct_lumen_obj), add_metadata$lumen_id)
ct_lumen_obj$Max_CNiche <- add_metadata$Max_CNiche
ct_lumen_obj$Max_TNiche <- add_metadata$Max_TNiche
ct_lumen_obj$C1 <- add_metadata$C1
ct_lumen_obj$C2 <- add_metadata$C2
ct_lumen_obj$C3 <- add_metadata$C3
ct_lumen_obj$C4 <- add_metadata$C4
ct_lumen_obj$C5 <- add_metadata$C5
ct_lumen_obj$C6 <- add_metadata$C6
ct_lumen_obj$C7 <- add_metadata$C7
ct_lumen_obj$C8 <- add_metadata$C8
ct_lumen_obj$C9 <- add_metadata$C9
ct_lumen_obj$C10 <- add_metadata$C10
ct_lumen_obj$C11 <- add_metadata$C11
ct_lumen_obj$C12 <- add_metadata$C12
ct_lumen_obj$T1 <- add_metadata$T1
ct_lumen_obj$T2 <- add_metadata$T2
ct_lumen_obj$T3 <- add_metadata$T3
ct_lumen_obj$T4 <- add_metadata$T4
ct_lumen_obj$T5 <- add_metadata$T5
ct_lumen_obj$T6 <- add_metadata$T6
ct_lumen_obj$T7 <- add_metadata$T7
ct_lumen_obj$T8 <- add_metadata$T8
ct_lumen_obj$T9 <- add_metadata$T9
ct_lumen_obj$T10 <- add_metadata$T10
ct_lumen_obj$T11 <- add_metadata$T11
ct_lumen_obj$T12 <- add_metadata$T12

# # Look at actual lumens before and after filtering
# all_lumens <- subset(xenium, subset = lumen_id %in% colnames(ct_lumen_obj))
# metadata_for_full_obj1 <- ct_lumen_obj@meta.data %>%
#   select(lumen_id, sample, sample_type, sample_affect, num_cells_lumen, num_cells_bin, max_dist, max_dist_bin,
#          num_cells_epithelial, num_cells_endothelial, num_cells_immune, num_cells_mesenchymal,
#          Max_CNiche, Max_TNiche,
#          C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12,
#          T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12)
# metadata_for_full_obj <- xenium@meta.data %>%
#   full_join(metadata_for_full_obj1) %>%
#   filter(lumen_id %in% colnames(ct_lumen_obj))
# all_lumens <- AddMetaData(all_lumens, metadata = metadata_for_full_obj)

# Filter object by lumen size
filtered_lumens <- subset(ct_lumen_obj, subset = num_cells_lumen >= 25 & num_cells_lumen <= 500 & max_dist >= 110)
ct_lumen_obj_sizefiltered <- subset(ct_lumen_obj, subset = lumen_id %in% filtered_lumens$lumen_id)
aft_lumens <- table(ct_lumen_obj_sizefiltered$sample); sum(aft_lumens) # 1811 lumens in TMAs 1-4; 3003 total

# Re-filter to require max CNiche be an alveolus-related niche
multiple_cniches <- names(ct_lumen_obj_sizefiltered$lumen_id[ct_lumen_obj_sizefiltered$Max_CNiche == "Multiple"])
multiple_cniches_KEEP <- ct_lumen_obj_sizefiltered@meta.data %>%
  filter(lumen_id %in% multiple_cniches) %>%
  select(lumen_id, Max_TNiche, paste0("C", 1:12)) %>%
  group_by(lumen_id) %>%
  mutate(Max_CNiche_value = max(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12)) %>%
  pivot_longer(cols = paste0("C", 1:12), names_to = "Max_CNiche", values_to = "Prop") %>%
  filter(Max_CNiche_value == Prop) %>%
  mutate(num_max_CNiches = length(Max_CNiche)) %>%
  filter(Max_CNiche %in% c("C2", "C3", "C5", "C8", "C11")) %>%
  mutate(num_max_CNiches_after_filtering = length(Max_CNiche)) %>%
  filter(num_max_CNiches_after_filtering == num_max_CNiches) %>%
  pull(lumen_id) %>%
  unique()
ct_lumen_obj_size_nichefiltered <- subset(ct_lumen_obj_sizefiltered, 
                                          subset = Max_CNiche %in% c("C2", "C3", "C5", "C8", "C11") | lumen_id %in% multiple_cniches_KEEP)
aft_lumens <- table(ct_lumen_obj_size_nichefiltered$sample); sum(aft_lumens) # 1310 lumens in TMAs 1-4; 1954 total

# Re-filter to remove airways
ct_lumen_obj_size_niche_airfiltered <- subset(ct_lumen_obj_size_nichefiltered, subset = C1 < 0.2)
aft_lumens <- table(ct_lumen_obj_size_niche_airfiltered$sample); sum(aft_lumens) # 1289 lumens in TMAs 1-4; 1923 total

# Re-filter to remove endothelial
ct_lumen_obj_size_niche_air_endofiltered <- subset(ct_lumen_obj_size_niche_airfiltered, subset = C10 < 0.3 & C12 < 0.3)
aft_lumens <- table(ct_lumen_obj_size_niche_air_endofiltered$sample); sum(aft_lumens) # 1271 lumens in TMAs 1-4; 1888 total

# Re-filter to require epithelial cells in the lumen
ct_lumen_obj_size_niche_air_endofiltered$prop_epi <- ct_lumen_obj_size_niche_air_endofiltered$num_cells_epithelial/ct_lumen_obj_size_niche_air_endofiltered$num_cells_lumen
ct_lumen_obj_size_niche_air_endo_epifiltered <- subset(ct_lumen_obj_size_niche_air_endofiltered, subset = prop_epi >= 0.05)
aft_lumens <- table(ct_lumen_obj_size_niche_air_endo_epifiltered$sample); sum(aft_lumens) # 1201 lumens in TMAs 1-4; 1747 total

# Save/read object
# saveRDS(lumen_xenium, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/lumen_xenium_071624.rds")
lumen_xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/lumen_xenium_071624.rds")


#### LUMEN ORDERING ----
ct_obj_to_use <- ct_lumen_obj_size_niche_air_endo_epifiltered

# Get proportion of each niche within cells within lumens
# First get files that have the cell and transcript data
tniche_trans_files1a <- list.files("/scratch/avannan/late_IPF_spatial/xenium/REVISION_graphsage_corrected_gene_labels_27May/gmm12_5k_trained/gmm12_5k_trained/trans_level/corrected_gene_labels",
                                   full.names = TRUE, pattern = ".csv")
names(tniche_trans_files1a) <- unlist(lapply(str_split(tniche_trans_files1a, "/|_"), function(XX) XX[25]))
tniche_trans_files2a <- list.files("/scratch/avannan/late_IPF_spatial/xenium/REVISION_graphsage_new_samples",
                                   full.names = TRUE, pattern = ".csv")
names(tniche_trans_files2a) <- unlist(lapply(str_split(tniche_trans_files2a, "/|_"), function(XX) XX[12]))
names(tniche_trans_files2a)[names(tniche_trans_files2a) == "TILD117MF"] <- "TILD117MFB"
tniche_trans_files_a <- c(tniche_trans_files1a, tniche_trans_files2a)

# Then get files with the transcript niches by transcript
tniche_trans_files1b <- list.files("/scratch/avannan/late_IPF_spatial/xenium/REVISION_graphsage_corrected_gene_labels_27May/gmm12_5k_trained/gmm12_5k_trained/trans_level/corrected_gene_labels/gmm12_5k_trained",
                                   full.names = TRUE, pattern = ".txt")
tniche_trans_files1b <- tniche_trans_files1b[!(grepl("SOME_REMOVED", tniche_trans_files1b))]
names(tniche_trans_files1b) <- unlist(lapply(str_split(tniche_trans_files1b, "/|_"), function(XX) XX[28]))
tniche_trans_files2b <- list.files("/scratch/avannan/late_IPF_spatial/xenium/REVISION_graphsage_new_samples/gmm12_mapped_to_existing_clusters",
                                   full.names = TRUE, pattern = ".txt")
names(tniche_trans_files2b) <- unlist(lapply(str_split(tniche_trans_files2b, "/|_"), function(XX) XX[17]))
names(tniche_trans_files2b)[names(tniche_trans_files2b) == "TILD117MF"] <- "TILD117MFB"
tniche_trans_files_b <- c(tniche_trans_files1b, tniche_trans_files2b)

# Get TNiche proportions by sample by transcript by lumen
tniche_trans_lumens_props_list <- lapply(names(tniche_trans_files_a), function(XX) {
  message(XX)
  
  # Get lumens in sample
  sm_lumens <- ct_obj_to_use@meta.data %>%
    filter(sample == XX) %>%
    pull(lumen_id)
  
  # Some samples don't have lumens left - skip
  if(length(sm_lumens) != 0) {
    # Read in per-transcript TNiche data
    cell_trans_info <- read_csv(tniche_trans_files_a[[XX]])
    trans_info <- read_csv(tniche_trans_files_b[[XX]], col_names = "TNiche")
    all_cell_trans_info <- cbind(cell_trans_info, trans_info) %>%
      mutate(TNiche = ordered(paste0("T", TNiche+1), levels = paste0("T", 1:12)))
    
    # Fix TILD117MFB IDs
    if (XX == "TILD117MFB") {
      all_cell_trans_info <- all_cell_trans_info %>%
        mutate(transcript_id = gsub("TILD117MF", "TILD117MFB", transcript_id),
               cell_id = gsub("TILD117MF", "TILD117MFB", cell_id))
    }
    
    # Get cells in sample lumens
    sm_cells_in_lumens_df <- xenium@meta.data %>%
      filter(lumen_id %in% sm_lumens) %>%
      rownames_to_column(var = "full_cell_id") %>%
      select(sample, lumen_id, full_cell_id)
    
    # Merge TNiche data and with lumen info
    all_cell_trans_info <- all_cell_trans_info %>%
      filter(overlaps_nucleus == 1, cell_id %in% sm_cells_in_lumens_df$full_cell_id) %>%
      full_join(sm_cells_in_lumens_df, by = c("cell_id" = "full_cell_id"))
    
    # Get proportions per lumen of TNiche transcripts
    tniche_props <- as.data.frame((table(all_cell_trans_info$lumen_id, all_cell_trans_info$TNiche)/
                                     rowSums(table(all_cell_trans_info$lumen_id, all_cell_trans_info$TNiche)))) %>%
      mutate(Var2 = paste0("trans_", Var2)) %>%
      dplyr::rename("lumen_id" = "Var1")
    tniche_props <- tniche_props %>% pivot_wider(names_from = "Var2", values_from = "Freq")
    tniche_props
    
  } else {
    "No Lumens"
  }
})
names(tniche_trans_lumens_props_list) <- names(tniche_trans_files_a)
# Only keep the samples that had lumens
tniche_trans_lumens_props_list <- tniche_trans_lumens_props_list[unname(which(unlist(
  lapply(tniche_trans_lumens_props_list, function(XX) { length(XX) != 1 }))))]
tniche_trans_lumens_props_df <- Reduce(full_join, tniche_trans_lumens_props_list)

# Get lumen ordering
lumen_ordering <- ct_obj_to_use@meta.data %>%
  full_join(tniche_trans_lumens_props_df) %>%
  mutate(pseudotime_rank_c8_tie_avg = rank(-C8, ties.method = "average"),
         pseudotime_rank_t4_tie_avg = rank(-T4, ties.method = "average"),
         pseudotime_rank_trans_t4_tie_avg = rank(-trans_T4, ties.method = "average"))

# Set up proportion data
cell_lumens_mat <- t(as.matrix(ct_obj_to_use@assays$CTcomp@layers$counts))
colnames(cell_lumens_mat) <- rownames(ct_obj_to_use)
rownames(cell_lumens_mat) <- colnames(ct_obj_to_use)
cell_lumens_prop_mat <- cell_lumens_mat/rowSums(cell_lumens_mat)

trans_t4_order <- lumen_ordering %>%
  select(lumen_id, pseudotime_rank_trans_t4_tie_avg) %>%
  arrange(pseudotime_rank_trans_t4_tie_avg) %>%
  pull(lumen_id)
ordered_by_trans_t4 <- cell_lumens_prop_mat[trans_t4_order, ]

# See basic heatmap of celltypes
Heatmap(t(ordered_by_trans_t4),
        row_names_gp = gpar(fontsize = 10),
        show_column_names = FALSE,
        cluster_columns = FALSE)

# Save lumen ordering and niche composition for Supplementary Table 10
trans_tniche_lumen_metadata <- lumen_ordering %>%
  select(lumen_id, sample, 
         pseudotime_rank_trans_t4_tie_avg,
         num_cells_lumen, max_dist, num_cells_epithelial,
         Max_CNiche, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12,
         Max_TNiche, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12,
         trans_T1, trans_T2, trans_T3, trans_T4, trans_T5, trans_T6, trans_T7, trans_T8, trans_T9, trans_T10, trans_T11, trans_T12)
# write.csv(trans_tniche_lumen_metadata, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/trans_tniche_lumen_metadata_080724.csv")




