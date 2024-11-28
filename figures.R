## SETTING ENVIRONMENT ----
library(scCustomize)
library(rlang)
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


# Set seed
set.seed(0309)
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)
filter <- dplyr::filter
select <- dplyr::select
pretty_umap <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     axis.title = element_text(hjust = 1))
theme_angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
theme_black <- theme(axis.text = element_text(color = "black"),
                     axis.line = element_line(color = "black"),
                     axis.ticks = element_line(color = "black"))


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

nuclei_niche_color_list <- list(`C1` = "#1036b5",
                                `C2` = "#6f774b",
                                `C3` = "#d10040",
                                `C4` = "#84cd5d",
                                `C5` = "#8E76DB",
                                `C6` = "#71c8a5",
                                `C7` = "#ffa26d",
                                `C8` = "#924373",
                                `C9` = "#a6513c",
                                `C10` = "black",
                                `C11` = "#72c7ff",
                                `C12` = "#cabd00")


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
                                    `T12` = "#9DCC00")


## SET CELL TYPE ORDER ----
new_ct_order <- c(
  # Endothelial 
  "Arteriole", "Capillary", "Venous", "Lymphatic",
  # Epithelial - Alveolar
  "AT1", "Transitional AT2", "AT2", "Proliferating AT2", "KRT5-/KRT17+",              
  # Epithelial - Airway
  "Basal", "Multiciliated", "Goblet", "RASC", "Secretory", "PNEC", "Proliferating Airway",
  # Immune - Lymphoid
  "B cells", "Proliferating B cells", "NK/NKT", "Proliferating NK/NKT", "Tregs", 
  "CD4+ T-cells", "CD8+ T-cells", "Proliferating T-cells", "Plasma", "pDCs",
  # Immune - Myeloid
  "Basophils", "cDCs", "Migratory DCs", "Langerhans cells", "Mast", 
  "Alveolar Macrophages", "Interstitial Macrophages", "SPP1+ Macrophages",
  "Macrophages - IFN-activated", "Monocytes/MDMs", "Neutrophils", "Proliferating Myeloid",
  # Mesenchymal
  "Activated Fibrotic FBs", "Adventitial FBs", "Alveolar FBs", "Inflammatory FBs", 
  "Myofibroblasts", "Subpleural FBs", "Proliferating FBs", "SMCs/Pericytes", "Mesothelial")


## READ IN OBJECTS AND CONFIGURE METADATA ----
xenium <- readRDS("/scratch/avannan/MANUSCRIPT_figures/cell_niches/xenium_CT_secondpass3_nichek12_neighborsk25.rds")
xenium$final_CT <- xenium$CT_secondpass3
xenium$final_lineage[xenium$final_CT == "Langerhans cells"] <- "Immune"
xenium$final_sublineage[xenium$final_CT == "Langerhans cells"] <- "Myeloid"
xenium$CNiche <- paste0("C", xenium$CT_secondpass3_nichek12_neighborsk25)
xenium$fov <- NULL # Remove niche FOV

# Add transcript niche info
tniche_cell_files1 <- list.files("/scratch/avannan/late_IPF_spatial/xenium/REVISION_graphsage_corrected_gene_labels_27May/gmm12_5k_trained/gmm12_5k_trained/cell_level",
                                 full.names = TRUE)
names(tniche_cell_files1) <- unlist(lapply(str_split(tniche_cell_files1, "/|_"), function(XX) XX[22]))
tniche_cell_files2 <- list.files("/scratch/avannan/late_IPF_spatial/xenium/REVISION_graphsage_new_samples/cell_level_tranx_niche_assignment_for_new_samples",
                                 full.names = TRUE)
names(tniche_cell_files2) <- unlist(lapply(str_split(tniche_cell_files2, "/|_"), function(XX) XX[20]))
names(tniche_cell_files2)[names(tniche_cell_files2) == "TILD117MF"] <- "TILD117MFB"
tniche_cell_files <- c(tniche_cell_files1, tniche_cell_files2)
tniche_tranx_list <- lapply(names(tniche_cell_files), function(XX) {
  if (XX %in% names(tniche_cell_files1)) {
    read_csv(tniche_cell_files[[XX]]) %>%
      select(full_cell_id, gmm12_5k_trained_hex_gmm) %>%
      unique() %>%
      mutate(sample = XX) %>%
      dplyr::rename("gmm12" = "gmm12_5k_trained_hex_gmm")
  } else if (XX != "TILD117MFB") {
    read_csv(tniche_cell_files[[XX]]) %>%
      select(cell_id, gmm12_mapped_to_existing_clusters_hex_gmm) %>%
      unique() %>%
      mutate(sample = XX) %>%
      dplyr::rename("full_cell_id" = "cell_id",
                    "gmm12" = "gmm12_mapped_to_existing_clusters_hex_gmm")
  } else {
    read_csv(tniche_cell_files[[XX]]) %>%
      select(cell_id, gmm12_mapped_to_existing_clusters_hex_gmm) %>%
      unique() %>%
      separate(cell_id, into = c("temp", "cell_id"), sep = "_") %>%
      mutate(sample = XX, 
             full_cell_id = paste0(sample, "_", cell_id)) %>%
      select(full_cell_id, gmm12_mapped_to_existing_clusters_hex_gmm, sample) %>%
      dplyr::rename("gmm12" = "gmm12_mapped_to_existing_clusters_hex_gmm")
  }
})
names(tniche_tranx_list) <- names(tniche_cell_files)
tniche_tranx_df <- Reduce(rbind, tniche_tranx_list)
xenium_meta <- xenium@meta.data %>%
  rownames_to_column(var = "full_cell_id") %>%
  left_join(tniche_tranx_df)
all.equal(xenium_meta$full_cell_id, colnames(xenium)) # Must be TRUE
xenium$TNiche <- paste0("T", xenium_meta$gmm12 + 1)

# Fix patient column
xenium$patient <- gsub("LF|MF|LFB", "", xenium$patient)

# Add new sample type column using percent pathology, and add other metadata
library(googlesheets4)
gs4_deauth()
path_sheet <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094")
path_sheet <- read_sheet(path_sheet, sheet = 4, skip = 1)
path_df <- path_sheet %>%
  as.data.frame() %>%
  dplyr::rename("Percent_pathology" = "Percent pathology",
                "sample" = "Sample")
updated_path_metadata <- full_join(xenium@meta.data, path_df) %>%
  mutate(sample_affect = case_when(sample_type == "Unaffected" ~ "Unaffected",
                                   Percent_pathology >= 75 ~ "More Affected",
                                   Percent_pathology < 75 ~ "Less Affected")) %>%
  mutate(`X.4` = ifelse(sample == "TILD117MFB", paste0("TILD117MFB_", cell_id), `X.4`))
xenium$percent_pathology <- updated_path_metadata$Percent_pathology[match(colnames(xenium), updated_path_metadata$`X.4`)]
xenium$sample_affect <- updated_path_metadata$sample_affect[match(colnames(xenium), updated_path_metadata$`X.4`)]

# Add spatial dimension reduction object separately
position_xy <- cbind(xenium$adj_x_centroid, xenium$adj_y_centroid)
row.names(position_xy) <- row.names(xenium@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
xenium[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_", assay = DefaultAssay(xenium))

# Re-save object
# saveRDS(xenium, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_070924.rds")
# xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_070924.rds")
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_lumens_071624.rds")

# Fix the Seurat metadata so only filtered lumens are included
# Load in lumen data
trans_tniche_lumen_metadata <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/trans_tniche_lumen_metadata.csv")

# Create temporary metadata with lumen info
meta_tmp <- xenium@meta.data
trans_tniche_lumen_metadata <- trans_tniche_lumen_metadata %>% 
  select(lumen_id, pseudotime_rank_trans_t4_tie_avg) %>%
  dplyr::rename("lumen_rank" = "pseudotime_rank_trans_t4_tie_avg")
new_meta <- full_join(meta_tmp, trans_tniche_lumen_metadata) %>%
  mutate(full_cell_id = paste0(sample, "_", cell_id))
all.equal(as.character(new_meta$cell_id), as.character(xenium$cell_id)) # Must be TRUE
rownames(new_meta) <- paste0(new_meta$sample, "_", new_meta$cell_id)
all.equal(rownames(new_meta), colnames(xenium)) # Must be TRUE

# Add lumen info back to Seurat object and remove lumen_id for lumens without a ranking
xenium <- AddMetaData(xenium, metadata = new_meta)
xenium$lumen_id[is.na(xenium$lumen_rank)] <- NA
length(unique(xenium$lumen_id))-1 == 1747 # Must be TRUE

# Re-save object (final before rename; no annotations)
# saveRDS(xenium, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_filtered_lumens_080124.rds")
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_filtered_lumens_080124.rds")

for (sm in unique(xenium$sample)) {
  tmp <- xenium@meta.data %>%
    filter(sample == sm, !is.na(lumen_id)) %>%
    select(cell_id, lumen_id) %>%
    dplyr::rename("group" = "lumen_id")
  write.csv(tmp, paste0("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_metadata_4explorer/", sm, "_final_lumens_2024-08-15.csv"))
}


## RENAME SAMPLES ----
xenium@meta.data <- xenium@meta.data %>%
  # Initial rename
  mutate(new_sample_name = case_when(sample_affect == "Less Affected" ~ paste0(patient, "LA"),
                                     sample_affect == "More Affected"  ~ paste0(patient, "MA"),
                                     TRUE ~ sample)) %>%
  # Fix additional names
  mutate(new_sample_name = case_when(sample == "TILD117MF" ~ "TILD117MA1",
                                     sample == "TILD117MFB" ~ "TILD117MA2",
                                     sample == "VUILD104LF" ~ "VUILD104MA1",
                                     sample == "VUILD104MF" ~ "VUILD104MA2",
                                     sample == "VUILD105LF" ~ "VUILD105MA1",
                                     sample == "VUILD105MF" ~ "VUILD105MA2",
                                     sample == "VUILD48LF" ~ "VUILD48LA1",
                                     sample == "VUILD48MF" ~ "VUILD48LA2",
                                     TRUE ~ new_sample_name))

# Check old and new names
name_matches <- xenium@meta.data %>%
  select(sample, new_sample_name, tma) %>%
  unique() %>%
  as_tibble()


## SET UP ANNOTATION DATAFRAME ----
# Add annotations
# Read in annotation files
annotation_files <- list.files("/scratch/avannan/late_IPF_spatial/xenium/REVISION_partitioned_updated_annotations",
                               full.names = TRUE, recursive = TRUE, pattern = "filtered_partitioned_cells.csv")
annotation_files1 <- annotation_files[!(grepl("IPFTMA5", annotation_files))]
annotation_files2 <- annotation_files[grepl("IPFTMA5", annotation_files)]
names(annotation_files1) <- unlist(lapply(str_split(annotation_files1, "_|/"), function(XX) XX[14]))
names(annotation_files2) <- unlist(lapply(str_split(annotation_files2, "_|/"), function(XX) XX[14]))
names(annotation_files2)[which(names(annotation_files2) == "TILD090MF")] <- "TILD080MF"
names(annotation_files2)[which(names(annotation_files2) == "TILD117MF")] <- "TILD117MFB"
all_anno_files <- c(annotation_files1, annotation_files2)

# Load in and format files into one dataframe
filtered_cells <- paste0(xenium$sample, "_", xenium$cell_id)
annotations_for_df_list <- list()
num_anno_list <- list()
for (sm in names(all_anno_files)) {
  # Get number of annotations per sample before filtering
  tmp <- read_csv(all_anno_files[[sm]]) %>%
    mutate(sample = sm, 
           full_cell_id = paste0(sm, "_", cell_id),
           .before = everything()) %>%
    select(-c(3,4), -c(6:15))
  num_anno_list[[sm]] <- ncol(tmp)-4
  
  annotations_for_df_list[[sm]] <- tmp %>%
    pivot_longer(cols = 4:(ncol(.)-1), names_to = "annotation_type_instance", values_to = "annotation_exists") %>%
    mutate(annotation_type_instance = gsub(".geojson", "", annotation_type_instance),
           annotation_type = gsub(".{2}$", "", annotation_type_instance),
           annotation_type = case_when(annotation_type == "ariway_smooth_muscle" ~ "airway_smooth_muscle",
                                       annotation_type == "granuloma_" ~ "granuloma",
                                       annotation_type == "venule_" ~ "venule",
                                       annotation_type == "epithelial_detatchment" ~ "epithelial_detachment",
                                       TRUE ~ annotation_type),
           annotation_instance = unlist(str_split(annotation_type_instance, "_"))[
             length(unlist(str_split(annotation_type_instance, "_")))],
           Annotation_Type = str_to_title(gsub("_", " ", annotation_type)),
           Annotation_Type = case_when(Annotation_Type == "Tls" ~ "TLS",
                                       Annotation_Type == "Hyperplastic Aec" ~ "Hyperplastic AECs",
                                       TRUE ~ Annotation_Type)) %>%
    filter(annotation_exists == TRUE, full_cell_id %in% filtered_cells) %>%
    select(-annotation_exists) %>%
    group_by(annotation_type_instance) %>%
    mutate(num_cells_annotation_type_instance = length(full_cell_id)) %>%
    ungroup() %>%
    dplyr::rename("original_annotation_sum" = "annotation_sum")
}
annotations_df <- Reduce(rbind, annotations_for_df_list) %>%
  left_join(xenium@meta.data %>% select(sample, cell_id, x_centroid, y_centroid, 
                                        CNiche, TNiche, final_CT, final_lineage, 
                                        final_sublineage))

# Should be 743 after filtering (749 total, including those without filtered cells)
annotations_df %>%
  select(sample, annotation_type, annotation_type_instance) %>%
  unique() %>%
  nrow()
sort(unique(annotations_df$Annotation_Type))

# Continue to filter out poor examples
annotations_df <- annotations_df %>%
  filter(!(sample == "VUILD104LF" & annotation_type_instance %in% c("epithelial_detachment_1", "epithelial_detachment_3", "epithelial_detachment_4")),
         !(sample == "VUILD104MF" & annotation_type == "epithelial_detachment"),
         !(sample == "VUILD106" & annotation_type_instance == "epithelial_detachment_3"),
         !(sample == "VUILD115" & annotation_type == "epithelial_detachment"),
         !(sample == "TILD117LF" & annotation_type_instance == "epithelial_detachment_3"),
         !(sample == "TILD117MF" & annotation_type_instance %in% c("epithelial_detachment_2", "epithelial_detachment_3", "epithelial_detachment_7")),
         !(sample == "TILD175" & annotation_type == "epithelial_detachment"),
         !(sample == "VUILD78MF" & annotation_type_instance == "epithelial_detachment_2"),
         !(sample == "TILD117MFB" & annotation_type == "epithelial_detachment"),
         !(sample == "TILD167LF" & annotation_type_instance == "epithelial_detachment_2"),
         !(sample == "VUILD49" & annotation_type == "epithelial_detachment"))

# Should be 720 after filtering for epithelial detachment
annotations_df %>%
  select(sample, annotation_type, annotation_type_instance) %>%
  unique() %>%
  nrow()

# Write files for Explorer
for (sm in names(annotations_for_df_list)) {
  tmp <- annotations_for_df_list[[sm]] %>%
    select(cell_id, annotation_type_instance) %>%
    dplyr::rename("group" = "annotation_type_instance")
  write.csv(tmp, paste0("/scratch/avannan/late_IPF_spatial/xenium/REVISION_partitioned_updated_annotations/explorer_annotations/",
                        sm, "_anno.csv"))
}

# Remove annotations that are selections of the entire sample
annotations_keep_df <- table(xenium$sample) %>%
  as.data.frame() %>%
  dplyr::rename("sample" = "Var1") %>%
  full_join(annotations_df) %>%
  mutate(prop_of_sample = num_cells_annotation_type_instance/Freq) %>%
  full_join(xenium@meta.data %>% select(sample, tma) %>% unique()) %>%
  select(sample, tma, annotation_type, annotation_type_instance, prop_of_sample) %>%
  unique() %>%
  filter(prop_of_sample < 0.95,
         !(sample == "VUILD142" & annotation_type_instance == "fibroblastic_focus_1")) # This is on the incorrect sample

# Should be 712 after filtering
annotations_keep_df %>%
    select(sample, annotation_type, annotation_type_instance) %>%
    unique() %>%
    nrow()
  
# Get final annotations df
final_annotations_df <- annotations_df %>%
  right_join(annotations_keep_df)
# Should be 712 after filtering
final_annotations_df %>%
  select(sample, annotation_type, annotation_type_instance) %>%
  unique() %>%
  nrow()
sort(unique(annotations_df$Annotation_Type))
# write.csv(final_annotations_df,
#           "/scratch/avannan/late_IPF_spatial/xenium/REVISION_partitioned_updated_annotations/final_annotations_df_082024.csv")
final_annotations_df <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION_partitioned_updated_annotations/final_annotations_df_082024.csv")

tmp <- final_annotations_df %>% 
  full_join(name_matches) %>%
  select(sample, new_sample_name, annotation_type_instance, annotation_type, annotation_instance) %>%
  unique() %>%
  group_by(new_sample_name, annotation_type) %>%
  mutate(num_annotation_type = length(annotation_type)) %>%
  ungroup()
supp_table_anno <- full_join(final_annotations_df, tmp)
supp_table_anno <- supp_table_anno %>%
  mutate(full_cell_id = paste0(new_sample_name, "_", cell_id)) %>%
  select(new_sample_name, full_cell_id, annotation_type, annotation_type_instance, annotation_instance, num_annotation_type)
colnames(supp_table_anno) <- c("Sample", "cell_id", "Annotation_Type", "Annotation_Type_Instance", "Annotation_Instance", "Num_Annotation_Type")
write.csv(supp_table_anno,
          "/scratch/avannan/late_IPF_spatial/xenium/REVISION_partitioned_updated_annotations/supp_table_anno_082124.csv")
for (sm in unique(supp_table_anno$Sample)) {
  tmp <- supp_table_anno %>%
    filter(Sample == sm)
  write.csv(tmp, paste0("/scratch/avannan/late_IPF_spatial/xenium/REVISION_partitioned_updated_annotations/", sm, "_supp_table_anno_082124.csv"))
}


## SUPPLEMENTARY TABLE 2 (CELL TYPE COUNT BY SAMPLE) ----
ct_sample_mat <- cbind(rowSums(table(xenium$sample, xenium$final_CT)),
                       table(xenium$sample, xenium$final_CT))
colnames(ct_sample_mat) <- c("Total Cells", colnames(ct_sample_mat)[-1])
# write.csv(ct_sample_mat, "/scratch/avannan/MANUSCRIPT_figures/SuppTable2_celltype_counts_by_sample.csv")


## FIGURE 1 (COMPONENTS OF BIORENDER FIGURE) ----
# This panel shows the data processing and was made primarily in Biorender. 
# UMAP:
plt <- DimPlot(xenium, group.by = "final_CT", cols = color_list, raster = FALSE,
               shuffle = TRUE) + NoLegend() + pretty_umap + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 15),
        axis.line = element_line(linewidth = 0.8)) + 
  labs(x = "UMAP 1", y = "UMAP 2") + coord_equal()
CairoPNG("/scratch/avannan/MANUSCRIPT_figures/xenium_late_UMAP_070924.png", width = 3.5, height = 3.5, units = "in", dpi = 450)
plt
dev.off()

# CNiches for sample VUILD96MF:
vuild96lf <- subset(xenium, subset = sample == "VUILD96LF")
vuild96lf$CNiche <- ordered(vuild96lf$CNiche, levels = paste0("C", 1:12))
plt2 <- DimPlot(vuild96lf, cols = nuclei_niche_color_list, group.by = "CNiche", reduction = "sp") + coord_equal() + pretty_umap +
  theme(plot.title = element_blank()) + NoLegend()
png("/scratch/avannan/MANUSCRIPT_figures/VUILD96LF_cniches_072224.png", width = 6, height = 6, units = "in", res = 450)
plt2
dev.off()

# Lumens for VUILD96MF:
kept_lumens_vuild96lf <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_metadata_4explorer/071824/VUILD96LF_lumen_filtering.csv") %>%
  filter(group == "Passed Size, Niche, Airway, Endo, Epi (All) Filters")
tmp <- xenium@meta.data %>%
  filter(sample == "VUILD96LF") %>%
  select(cell_id, lumen_id) %>%
  filter(!is.na(lumen_id), cell_id %in% kept_lumens_vuild96lf$cell_id)
colnames(tmp) <- c("cell_id", "group")
rownames(tmp) <- NULL
write.csv(tmp, paste0("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_metadata_4explorer/VUILD96LF_lumens.csv"))


## FIGURE 2 (CELL ATLASING) ----
#### FIGURE 2A (Cell Type Proportions Overall) ----
# Cell type proportions
lineage_color_list <- list(Endothelial = "chocolate2",
                           Epithelial = "chartreuse4",
                           Immune = "deeppink3",
                           Mesenchymal = "cornflowerblue")

ct_lineage_plot_list <- lapply(c("Endothelial", "Epithelial", "Immune", "Mesenchymal"), function(XX) {
  # Get strip color
  strip_color = lineage_color_list[[XX]]
  
  # Set round_count_text_size
  round_count_text_size <- ifelse(XX == "Immune", 1.23, 1.55)
  round_count_angle <- ifelse(XX == "Immune", 25, 0)
  round_count_vjust <- ifelse(XX == "Immune", -0.7, -0.4)
  round_count_hjust <- ifelse(XX == "Immune", 0.2, 0.5)
  
  # Set y limit
  y_max <- sort(table(xenium$final_CT, xenium$final_lineage)[, XX], decreasing = TRUE)[1]
  y_limit <- y_max + 0.12*(y_max)
  
  plt <- xenium@meta.data %>%
    group_by(final_CT, final_lineage) %>%
    summarize(count = length(final_CT)) %>%
    mutate(round_count = ifelse(count >= 1000,
                                paste0(format(round(count/1000, 1), n.small = 2, decimal.mark = "."), "k"),
                                as.character(count))) %>%
    filter(final_lineage == XX) %>%
    ggplot(aes(x = reorder(final_CT, -count), y = count, fill = final_CT)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~final_lineage, scales = "free") +
    scale_fill_manual(values = color_list) +
    # scale_x_discrete(labels = rename_cells) +
    scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0, 0)), n.breaks = 6,
                       labels = scales::unit_format(unit = "K", sep = "", scale = 1e-3)) +
    theme_bw(base_size = 6) +
    theme(axis.text = element_text(color = "black"),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.ticks = element_line(color = "black"),
          panel.grid.major.x = element_blank(), 
          axis.title.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_line(color = "grey80"),
          strip.text = element_text(color = "white", size = 7),
          strip.background = element_rect(fill = strip_color),
          axis.line = element_line(color = "black"),
          legend.position = "none") +
    geom_text(aes(label = round_count, vjust = -0.4), size = round_count_text_size)
  
  # Fix y axis label
  if (XX == "Endothelial" | XX == "Epithelial") {
    plt <- plt +
      labs(y = "Number of Cells") +
      theme(plot.margin = margin(5.2, 0, 0, 2, "pt"))
  } else {
    plt <- plt + theme(axis.title.y = element_blank(),
                       plot.margin = margin(1.8, 0, 0, 2, "pt"))
  }
  
  plt
})

a <- wrap_plots(ct_lineage_plot_list[[1]], ct_lineage_plot_list[[3]], ncol = 2, widths = c(0.18, 1))
b <- wrap_plots(ct_lineage_plot_list[[2]], ct_lineage_plot_list[[4]], ncol = 2, widths = c(1, 0.78))
wrap_plots(a, b, ncol = 1)

pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2A_ct_prop_and_compare_062824.pdf", width = (105*0.0393701), height = (88*0.0393701))
wrap_plots(a, b, ncol = 1)
dev.off()


#### FIGURE 2B (AT2/AT1 RATIO) ----
# Get results for eQTL data
eqtl_alv_ratio <- read.csv("/scratch/avannan/eQTL_rds_files/SE227136_ILD_all_celltypes_seurat_meta.csv") %>%
  filter(grepl("T", .$Sample_Name)) %>% # Only retain TGen samples
  select(Sample_Name, Status, manual_annotation_1) %>%
  group_by(Sample_Name, Status, manual_annotation_1) %>%
  summarize(num_cells = length(manual_annotation_1)) %>%
  filter(manual_annotation_1 %in% c("AT1", "AT2")) %>%
  pivot_wider(names_from = manual_annotation_1, values_from = num_cells) %>%
  mutate(`AT2/AT1 Ratio` = AT2/AT1,
         Dataset = "Natri et al.\n(2023)")

# Get results for Xenium data
xenium_alv_ratio <- xenium@meta.data %>%
  group_by(sample, disease_status, final_CT) %>%
  summarize(num_cells = length(final_CT)) %>%
  filter(final_CT %in% c("AT1", "AT2")) %>%
  pivot_wider(names_from = final_CT, values_from = num_cells) %>%
  mutate(`AT2/AT1 Ratio` = AT2/AT1,
         Dataset = "Spatial") %>%
  dplyr::rename(Sample_Name = sample, Status = disease_status)

# Join data
alv_ratio_df <- full_join(eqtl_alv_ratio, xenium_alv_ratio) %>%
  group_by(Dataset, Status) %>%
  mutate(mean_ratio = mean(`AT2/AT1 Ratio`, na.rm = TRUE),
         median_ratio = median(`AT2/AT1 Ratio`, na.rm = TRUE)) %>%
  group_by(Dataset, Status) %>%
  mutate(sum_AT1 = sum(AT1, na.rm = TRUE),
         sum_AT2 = sum(AT2, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(overall_ratio = sum_AT2/sum_AT1,
         overall_ratio2 = sum_AT1/sum_AT2,
         Dataset = ordered(Dataset, levels = c("Spatial", "Natri et al.\n(2023)")))
alv_ratio_df2 <- alv_ratio_df %>%
  select(Status, Dataset, overall_ratio, overall_ratio2) %>%
  unique() %>%
  mutate(y_pos = ifelse(Dataset == "Spatial", 200, 170),
         overall_ratio = round(overall_ratio, 1))

# Create plot
alv_ratio_plot <- alv_ratio_df %>%
  ggplot(aes(x = Dataset, y = `AT2/AT1 Ratio`, fill = Status)) +
  geom_boxplot(outlier.shape = 21, size = 0.25, outlier.size = 0.75,
               outlier.stroke = 0.25, color = "black", width = 0.7, 
               position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 6.5) +
  theme_black +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        plot.background = element_blank(),
        axis.title.y = element_text(size = 6),
        legend.position = "none") +
  scale_fill_manual(values = c("#88CCEE", "#d76160")) +
  scale_color_manual(values = c("#88CCEE", "#d76160")) +
  scale_y_continuous(limits = c(0, 224), expand = c(0, 0)) +
  geom_text(data = alv_ratio_df2,
            aes(y = 212, label = overall_ratio, color = Status),
            position = position_dodge(width = 1),
            fontface = "bold", size = 2, show.legend = FALSE)
alv_ratio_plot
# Scores could not be calculated for 3 control eQTL samples and 6 disease spatial samples (no AT1 cells)

pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2A_alv_ratio_plot_071524.pdf", width = (36*0.0393701), height = (28.5*0.0393701))
alv_ratio_plot
dev.off()


#### FIGURE 2C (MINI DOTPLOT HEATMAP) ----
xenium <- ScaleData(xenium) # Scale

nonprolif_cts <- unique(xenium$final_CT)[!grepl("Proliferating", unique(xenium$final_CT))]
mini_xenium <- subset(xenium, subset = final_CT %in% nonprolif_cts)
new_ct_order2 <- c(
  # Endothelial 
  "Arteriole", "Capillary", "Venous", "Lymphatic",
  # Epithelial - Alveolar
  "AT1", "Transitional AT2", "AT2",  "KRT5-/KRT17+",              
  # Epithelial - Airway
  "Basal", "Multiciliated", "Goblet", "RASC", "Secretory", "PNEC",
  # Immune - Lymphoid
  "B cells", "NK/NKT", "Tregs", 
  "CD4+ T-cells", "CD8+ T-cells", "Plasma", "pDCs",
  # Immune - Myeloid
  "Basophils", "cDCs", "Migratory DCs", "Langerhans cells", "Mast", 
  "Alveolar Macrophages", "Interstitial Macrophages", "SPP1+ Macrophages",
  "Macrophages - IFN-activated", "Monocytes/MDMs", "Neutrophils",
  # Mesenchymal
  "Activated Fibrotic FBs", "Adventitial FBs", "Alveolar FBs", "Inflammatory FBs", 
  "Myofibroblasts", "Subpleural FBs", "SMCs/Pericytes", "Mesothelial")

genes_keep <- c("PECAM1",
                "EPCAM", "AGER", "SFTPC", "SCGB3A2",
                "PTPRC", "MS4A1", "CD3E", "JCHAIN", 
                "CPA3", "HLA-DQA1", "CD68", "CD14", "S100A8",
                "COL1A1")
# Create dotplot base
p <- DotPlot(mini_xenium, features = genes_keep, group.by = "final_CT", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
p

# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() %>%
  dplyr::select(features.plot, new_ct_order2)

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  tidyr::pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() %>%
  dplyr::select(features.plot, new_ct_order2)

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# Replace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to bright red
col_fun = circlize::colorRamp2(c(-1, 0, 2), colorspace::diverge_hsv(3))

rename_cts <- colnames(exp_mat)
rename_cts <- gsub("Proliferating", "Prolif.", rename_cts)
rename_cts <- gsub("Macrophages", "MΦ", rename_cts)
rename_cts <- gsub("Activated", "Activ.", rename_cts)
rename_cts <- gsub("Fibrotic", "Fibr.", rename_cts)
rename_cts <- gsub("activated", "activ.", rename_cts)
rename_cts <- gsub("Inflammatory", "Inflam.", rename_cts)
colnames(exp_mat) <- rename_cts
ct_col_order2 <- c(rep("   Endothelial   ", 4),
                   rep("  Epithelial  ", 10), 
                   rep("  Lymphoid  ", 7),
                   rep("  Myeloid  ", 11),
                   rep(" Mesenchymal ", 8))

# Creating a layer to add to plot
layer_fun = function(j, i, x, y, w, h, fill){
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex(t(exp_mat), i, j))), # t(exp_mat)
              size = pindex(t(percent_mat), i, j)/100 * unit(1, "mm"), # t(percent_mat)
              pch = 19)
}
ht_opt$DIMNAME_PADDING = unit(0.5, "mm")
hp <- Heatmap(t(exp_mat), # t(exp_mat)
              name = "Scaled Expression",
              height = unit(62.5, "mm"),
              width = unit(22.5, "mm"),
              # height = unit(100, "mm"),
              # heatmap_legend_param = list(title = "Scaled\n",
              #                             legend_direction = "vertical",
              #                             labels = c("-1", "0", "1", "2+"),
              #                             title_position = "topcenter",
              #                             title_gp = gpar(fontsize = 5),
              #                             labels_gp = gpar(fontsize = 5),
              #                             grid_width = unit(3, "mm"),
              #                             grid_height = unit(3, "mm"),
              #                             legend_height = unit(10, "mm")),
              column_title = " ",
              row_title_gp = gpar(fontsize = 0),
              column_title_gp = gpar(fontsize = 5),
              # row_title_side = "right",
              column_title_side = "top",
              # column_names_rot = 45,
              row_split = ct_col_order2,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 5),
              row_names_side = "left",
              column_names_gp = gpar(fontsize = 5),
              cluster_rows = FALSE, cluster_columns = FALSE,
              cluster_row_slices = FALSE,
              column_dend_height = unit(2, "mm"),
              column_dend_side = "bottom",
              border = "black",
              gap = unit(0.6, "mm"),
              show_heatmap_legend = FALSE)
ht <- draw(hp, 
           heatmap_legend_list = lgd_list,
           heatmap_legend_side = "right",
           annotation_legend_side = "right",
           align_heatmap_legend = "global_center")
draw(hp)
pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2_mini_dotplotheatmap.pdf",
    width = 43.4*0.0393701, height = 76*0.0393701)
draw(hp)
dev.off()

lgd1 <- Legend(labels = c(0,0.25,0.5,0.75,1), title = "Proportion", type = "points", 
               pch = 19, size = c(0,0.25,0.5,0.75,1) * unit(1, "mm"),
               legend_gp = gpar(col = "black"), direction = "vertical",
               ncol = 1,
               title_gap = unit(1.5, "mm"),
               title_position = "leftcenter-rot", 
               background = NULL, 
               # gap = unit(0, "mm"),
               # row_gap = unit(0, "mm"),
               # column_gap = unit(0, "mm"),
               grid_height = unit(0, "mm"),
               grid_width = unit(0, "mm"),
               title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5))
pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2_mini_heatmap_legend1.pdf",
    width = 10*0.0393701, height = 20*0.0393701)
draw(lgd1)
dev.off()

lgd1 <- Legend(labels = rev(c(0.25,0.5,0.75,1)), title = "Proportion", type = "points", 
               pch = 19, size = rev(c(0.25,0.5,0.75,1) * unit(1, "mm")),
               legend_gp = gpar(col = "black"), direction = "horizontal",
               ncol = 1,
               title_gap = unit(1.5, "mm"),
               row_gap = unit(0, "mm"),
               labels_rot = 90,
               title_position = "leftcenter-rot", 
               background = NULL, 
               # gap = unit(0, "mm"),
               # row_gap = unit(0, "mm"),
               # column_gap = unit(0, "mm"),
               grid_height = unit(0, "mm"),
               grid_width = unit(0, "mm"),
               title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5))
pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2_mini_heatmap_legend1_v2.pdf",
    width = 30*0.0393701, height = 30*0.0393701)
draw(lgd1)
dev.off()

lgd2 <- Legend(col_fun = col_fun, at = c(-1, 0, 1, 2),
               labels = c("-1", "0", "1", "2+"),
               title = "Scaled",
               # legend_gp = gpar(col = "black"), direction = "vertical",
               title_gap = unit(1.5, "mm"),
               title_position = "leftcenter-rot", 
               background = NULL, 
               legend_height = unit(1, "mm"),
               # gap = unit(0, "mm"),
               # row_gap = unit(0, "mm"),
               # column_gap = unit(0, "mm"),
               grid_height = unit(0, "mm"),
               grid_width = unit(1, "mm"),
               title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5))
pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2_mini_heatmap_legend2.pdf",
    width = 10*0.0393701, height = 20*0.0393701)
draw(lgd2)
dev.off()

lgd2 <- Legend(col_fun = col_fun, at = c(-1, 0, 1, 2),
               labels = c("-1", "0", "1", "2+"),
               title = "Scaled\nExpression",
               # legend_gp = gpar(col = "black"), direction = "vertical",
               title_gap = unit(1, "mm"),
               title_position = "topcenter", direction = "horizontal", 
               background = NULL, 
               legend_width = unit(9, "mm"),
               # gap = unit(0, "mm"),
               # row_gap = unit(0, "mm"),
               # column_gap = unit(0, "mm"),
               grid_height = unit(1, "mm"),
               grid_width = unit(1, "mm"),
               title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5))
pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2_mini_heatmap_legend2_v2.pdf",
    width = 14.2*0.0393701, height = 10*0.0393701)
draw(lgd2)
dev.off()


#### FIGURE 2D/E (AIRWAY SPATIAL EXAMPLES) ----
other_air_colors <- c("#c30039",
                      "#b2d350",
                      "#498aff",
                      "#175f26",
                      "#ff316a",
                      "#714270",
                      "#fead15")
tild299lf <- subset(xenium, subset = sample == "TILD299LF" & (final_sublineage == "Airway"))
tild299lf_dimplot <- DimPlot(tild299lf, reduction = "sp", group.by = "final_CT", 
                             raster = TRUE, pt.size = 3,
                             cols = other_air_colors) + coord_equal() +
  theme_classic(base_size = 6.5) +
  theme_black +
  theme(plot.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 5.5),
        plot.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        
        legend.position = "bottom",
        # legend.direction = "vertical",
        legend.key.size = unit(2.2, "mm"),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.background = element_blank(),
        legend.justification = "left") +
  scale_y_continuous(limits = c(-15250, -14250), expand = c(0, 0)) +
  scale_x_continuous(limits = c(25450, 26500), expand = c(0, 0)) +
  guides(color = guide_legend(ncol = 1)) +
  # theme(legend.position = c(0.6, 0.75)) +
  pretty_umap
tild299lf_dimplot

# pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2A_tild299lf_dimplot_071524.pdf", width = (39*0.0393701), height = (39*0.0393701))
# pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2A_tild299lf_dimplot_071524.pdf", width = (39*0.0393701), height = (48*0.0393701))
pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2A_tild299lf_dimplot_071524.pdf", width = (37.5*0.0393701), height = (56*0.0393701))
tild299lf_dimplot
dev.off()

tild299lf_featplot <- FeaturePlot_scCustom(tild299lf, reduction = "sp",
                                           features = c("TP63", "MUC5B", "C20orf85", "SCGB3A2", "SCGB1A1", "CALCA"),
                                           pt.size = 4, order = TRUE, num_columns = 2, raster = TRUE) &
  theme_classic(base_size = 6.5) & theme_black & pretty_umap &
  theme(axis.line = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 5.5, hjust = 0.5, face = "bold"),
        legend.background = element_blank(),
        plot.margin = margin(0, 0.5, 0, 0.5, "mm"),
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(1, "mm")) &
  coord_equal() &
  scale_y_continuous(limits = c(-15250, -14250), expand = c(0, 0)) &
  scale_x_continuous(limits = c(25450, 26500), expand = c(0, 0))
tild299lf_featplot

pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2A_tild299lf_featureplot_071524.pdf", width = (52*0.0393701), height = (56*0.0393701))
tild299lf_featplot
dev.off()


#### FIGURE 2F/G (ENDOTHELIAL SPATIAL EXAMPLES) ----
endo_colors2 <- list(`Arteriole` = "#e22039",
                     `Capillary` = "#c8c0f2",
                     `Lymphatic` = "#0272db",
                     `Venous` = "#ccc100")

vuild107mf <- subset(xenium, subset = sample == "VUILD107MF" & (final_lineage == "Endothelial" | final_CT == "SMCs/Pericytes"))
vuild107mf$tmp_CT <- vuild107mf$final_CT
vuild107mf$tmp_CT[vuild107mf$final_CT == "SMCs/Pericytes"] <- "\nSMCs/Pericytes\n(Mesenchymal)"
vuild107mf$tmp_CT <- ordered(vuild107mf$tmp_CT, levels = c("Arteriole", "Capillary", "Lymphatic", "Venous", "\nSMCs/Pericytes\n(Mesenchymal)"))
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ENDOTHELIAL_VUILD107MF_DimPlot.pdf",
    width = (39.2*0.0393701), height = (49.2*0.0393701))
DimPlot(vuild107mf, reduction = "sp", group.by = "tmp_CT", pt.size = 2.5,
        raster = TRUE,
        cols = c(unlist(endo_colors2), `\nSMCs/Pericytes\n(Mesenchymal)` = "grey90")) + coord_equal() +
  coord_equal() + theme_classic(base_size = 6.5) +
  theme_black + pretty_umap +
  guides(color = guide_legend(keywidth = unit(1, "mm"), keyheight = unit(1, "mm"),
                              ncol = 3, byrow = FALSE)) +
  theme(plot.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 5.5),
        plot.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0, "mm"),
        legend.box.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = "bottom",
        legend.justification = "left")
dev.off()

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ENDOTHELIAL_VUILD107MF_FeaturePlot.pdf",
    width = 84*0.0393701, height = 46.5*0.0393701)
FeaturePlot_scCustom(subset(vuild107mf, subset = final_lineage == "Endothelial"), 
                     features = c("HEY1", "ACKR1", "CCL21", "CA4", "FCN3", "APLNR"),
                     reduction = "sp", pt.size = 4, raster = TRUE,
                     num_columns = 3, order = TRUE) &
  theme_bw(base_size = 5.5) & pretty_umap & coord_equal() &
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(), 
        plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 5),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(0.5, 0, 0, 0, "mm"))
dev.off()


#### FIGURE 2H/I (MESENCHYMAL SPATIAL EXAMPLES) ----
other_mes_colors <- c(`Activated Fibrotic FBs` = "#c810e8",
                      `Adventitial FBs` = "#ff57b2",
                      `Alveolar FBs` = "blanchedalmond",
                      `Inflammatory FBs` = "black",
                      `Mesothelial` = "#bed00d",
                      `Myofibroblasts` = "#24cdff",
                      `SMCs/Pericytes` = "#ab6600",
                      `Subpleural FBs` = "#00a754")
vuild106 <- subset(xenium, subset = sample == "VUILD106" &
                     final_CT %in% c("Adventitial FBs", "Alveolar FBs", "Activated Fibrotic FBs",
                                     "Inflammatory FBs",
                                     "Mesothelial", "Subpleural FBs", "Myofibroblasts", "SMCs/Pericytes"))
vuild106_dimplot <- DimPlot(vuild106, reduction = "sp", group.by = "final_CT", 
                            raster = TRUE, pt.size = 1.5, na.value = NA, 
                            cols = other_mes_colors) +
  coord_equal() + theme_classic(base_size = 6.5) +
  scale_y_continuous(limits = c(-26650, -20800)) +
  theme_black + pretty_umap +
  guides(color = guide_legend(keywidth = unit(1, "mm"), keyheight = unit(1, "mm"),
                              ncol = 1, byrow = TRUE)) +
  theme(plot.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 5.5),
        plot.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"),
        legend.position = "right",
        legend.justification = "center")

DimPlot(xenium, group.by = "sample", label = TRUE, reduction = "sp") + NoLegend()
vuild106_dimplot

pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2A_vuild106_dimplot_071724.pdf", width = (69.4*0.0393701), height = (41.2*0.0393701))
vuild106_dimplot
dev.off()

vuild106_featplot <- FeaturePlot_scCustom(subset(vuild106, subset = final_lineage == "Mesenchymal"), 
                                          features = c("COL1A1", "CTHRC1", "POSTN", "MFAP5", "PTGDS", "MSLN", "WNT5A", "ACTA2", "CSPG4", "PLIN2"),
                                          reduction = "sp", raster = TRUE, order = TRUE, num_columns = 5, pt.size = 2) & 
  theme_classic(base_size = 6.5) &
  theme_black & coord_equal() & pretty_umap &
  theme(axis.line = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 5.5, hjust = 0.5, face = "bold"),
        legend.background = element_blank(),
        plot.margin = margin(0.5, 0.5, 0, 0, "mm"),
        legend.key.height = unit(2, "mm"),
        legend.spacing = unit(0, "mm"),
        legend.box.background = element_blank(),
        legend.key.width = unit(1, "mm"))
vuild106_featplot

pdf("/scratch/avannan/MANUSCRIPT_figures/NewFig2A_vuild106_featplot_071724.pdf", width = (115*0.0393701), height = (41.2*0.0393701))
vuild106_featplot
dev.off()


## FIGURE 3 (PERCENT PATHOLOGY & ANNOTATION OVERVIEW) ----
#### FIGURE 3B (Percent Pathology x CT Volcano Plots) ----
library(googlesheets4)
gs4_deauth()
path_ct_sheet <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094")
path_ct_sheet <- read_sheet(path_ct_sheet, sheet = 6, skip = 1)
pathscore_ct_df <- as.data.frame(path_ct_sheet)

pathscore_ct_plot <- pathscore_ct_df %>%
  mutate(neglog10.adjP = -log10(adj.P.Val),
         sig = ifelse(adj.P.Val < 0.01, "sig", "n.s."),
         celltype = ifelse(sig == "sig", celltype, "Other"),
         celltype = ifelse(celltype == "Alveolar Macrophages", "Alveolar MΦ", celltype),
         celltype = ifelse(celltype == "Activated Fibrotic FBs", "Activ. Fibr. FBs", celltype),
         label = ifelse(sig == "sig", celltype, NA_character_)) %>%
  ggplot(aes(x = coefficient, y = neglog10.adjP, fill = celltype, color = celltype)) +
  geom_hline(yintercept = -log10(0.01), lty = "longdash", color = "grey40", linewidth = 0.25) +
  geom_vline(xintercept = 0, lty = "longdash", color = "grey40", linewidth = 0.25) +
  geom_point(shape = 21, stroke = 0.2, alpha = 0.75, size = 1.5) +
  geom_text_repel(aes(label = label), fontface = "bold", size = 1.7, seed = 715,
                  min.segment.length = unit(0.18, "lines"),
                  box.padding = 0.15, na.rm = TRUE) +
  theme_classic(base_size = 6.5) +
  theme(axis.line = element_line(color = "black"), 
        axis.text = element_text(color = "black", size = 5.5),
        axis.title = element_text(color = "black", size = 6)) +
  scale_x_continuous(limits = c(-0.04, 0.045)) +
  scale_color_manual(values = c(color_list, `Other` = "grey85", `Alveolar MΦ` = "#ffa3a1", `Activ. Fibr. FBs` = "#005a90")) +
  scale_fill_manual(values = c(color_list, `Other` = "grey85", `Alveolar MΦ` = "#ffa3a1", `Activ. Fibr. FBs` = "#005a90")) +
  NoLegend() +
  labs(x = "Rate of Change (Percent Pathology)", y = "-log10(FDR)")
pathscore_ct_plot

# pdf("/scratch/avannan/MANUSCRIPT_figures/oldFig2_volcano_ct.pdf", width = (53*0.0393701), height = (50*0.0393701))
pdf("/scratch/avannan/MANUSCRIPT_figures/oldFig2_volcano_ct.pdf", width = (52*0.0393701), height = (56*0.0393701))
pathscore_ct_plot
dev.off()


#### FIGURE 3C (Annotation x CT Dotplot) ----
# Create dataframe of lineage variables
lineage_df <- xenium@meta.data %>%
  select(final_CT, final_lineage) %>%
  unique()
names(lineage_df) <- c("CT", "Lineage")

new_annotation_order <- c(# Epithelial
  "Normal Alveoli",
  "Minimally Remodeled Alveoli",
  "Hyperplastic AECs",
  "Emphysema",
  "Remodeled Epithelium",
  "Advanced Remodeling",
  "Epithelial Detachment",
  "Remnant Alveoli",
  "Small Airway",
  "Large Airway",
  "Microscopic Honeycombing",
  "Goblet Cell Metaplasia",
  
  # Immune
  "Granuloma",
  "Mixed Inflammation",
  "TLS",
  
  # Endothelial/Mesenchymal
  "Artery",
  "Muscularized Artery",
  "Venule",
  "Interlobular Septum",
  "Airway Smooth Muscle",
  "Fibroblastic Focus",
  "Fibrosis",
  "Severe Fibrosis",
  
  # Other
  "Multinucleated Cell",
  "Giant Cell")

# Annotation Composition
path_ct_prop_df <- as.data.frame((table(final_annotations_df$Annotation_Type, final_annotations_df$final_CT)/
                                    rowSums(table(final_annotations_df$Annotation_Type, final_annotations_df$final_CT))))
names(path_ct_prop_df) <- c("Pathology", "CT", "Proportion")

new_ct_order_4plot <- gsub("Proliferating", "Prolif.", gsub("Activated", "Activ.", gsub("Macrophages", "MΦ", new_ct_order)))
color_list_4plot <- color_list
names(color_list_4plot) <- gsub("Proliferating", "Prolif.", gsub("Activated", "Activ.", gsub("Macrophages", "MΦ", names(color_list_4plot))))
new_annotation_order_4plot <- c("Normal Alveoli", "Min. Remod. Alveoli", "Hyperplastic AECs",
                                "Advanced Remodeling", "Epithelial Detachment", "Micro. Honeycombing",
                                "Goblet Cell Metaplasia", "Fibroblastic Focus", "Severe Fibrosis",
                                "Granuloma", "TLS")
path_ct_prop_plot <- path_ct_prop_df %>%
  left_join(lineage_df) %>%
  # Only keep some annotations for this plot
  filter(Pathology %in% c("Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                          "Advanced Remodeling", "Epithelial Detachment", "Microscopic Honeycombing",
                          "Goblet Cell Metaplasia", "Granuloma", "TLS",
                          "Fibroblastic Focus", "Severe Fibrosis")) %>%
  mutate(Pathology = gsub("Minimally Remodeled", "Min. Remod.", gsub("Microscopic", "Micro.", Pathology)),
         Pathology = ordered(Pathology, levels = new_annotation_order_4plot)) %>%
  mutate(CT = gsub("Proliferating", "Prolif.", gsub("Activated", "Activ.", gsub("Macrophages", "MΦ", CT))),
         CT = ordered(CT, levels = new_ct_order_4plot),
         Lineage = gsub("Endothelial", "Endo.", Lineage)) %>%
  filter(Proportion >= 0.01) %>%
  ggplot(aes(x = Pathology,
             y = reorder(CT, dplyr::desc(CT)),
             size = Proportion, fill = CT)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6.5) +
  scale_fill_manual(values = color_list_4plot) +
  theme(axis.text = element_text(color = "black", size = 5),
        strip.text = element_text(color = "white", size = 5),
        axis.text.y = element_text(size = 5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5),
        panel.grid.major.y = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, size = 6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.background = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.background = element_blank(),
        plot.margin = margin(5.5, 5.5, 8, 8, "pt"),
        legend.margin = margin(1, 1, 1, 1, "pt"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm")) +
  theme_angle +
  labs(title = "Annotation Composition (Cells)") +
  guides(fill = "none",
         size = guide_legend(title.theme = element_text(size = 5),
                             size = 2,
                             nrow = 1, ncol = 6, title.hjust = 0.4,
                             label.position = "bottom",
                             title.position = "top",
                             label.theme = element_text(size = 5))) +
  scale_size_continuous(limits = c(0, 1), range = c(0, 4), breaks = seq(0.2, 1, 0.2)) +
  facet_grid(Lineage~., space = "free", scales = "free")
path_ct_prop_plot

# Change facet colors
path_ct_prop_plot_table <- ggplot_gtable(ggplot_build(path_ct_prop_plot))
striprt <- which(grepl("strip-r", path_ct_prop_plot_table$layout$name) | 
                   grepl("strip-t", path_ct_prop_plot_table$layout$name))
fills <- c("chocolate2", "chartreuse4", "maroon3", "cornflowerblue") # Lineages
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", path_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  path_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(path_ct_prop_plot_table)

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig2_anno_ct_partial.pdf", width = (54*0.0393701), height = (115*0.0393701))
grid::grid.draw(path_ct_prop_plot_table)
dev.off()


#### FIGURE 3D (Example Annotations) ----
vuild91lf <- subset(xenium, subset = sample == "VUILD91LF" & CNiche == "C3")
plt <- DimPlot(vuild91lf, reduction = "sp", group.by = "CNiche", cols = nuclei_niche_color_list) + coord_equal() + pretty_umap + NoLegend() +
  theme(plot.title = element_blank())
pdf("/scratch/avannan/MANUSCRIPT_figures/VUILD91LF_C3niche.pdf", width = (200*0.0393701), height = (200*0.0393701))
plt
dev.off()


vuild91mf <- subset(xenium, subset = sample == "VUILD91MF")
plt <- DimPlot(vuild91mf, reduction = "sp", group.by = "CNiche", cols = nuclei_niche_color_list) + coord_equal() + pretty_umap + NoLegend() +
  theme(plot.title = element_blank())
pdf("/scratch/avannan/MANUSCRIPT_figures/test.pdf", width = (200*0.0393701), height = (200*0.0393701))
plt
dev.off()

vuild91mf <- subset(xenium, subset = sample == "VUILD91MF")
tmp <- vuild91mf@meta.data %>% filter(CNiche == "C3")
write.csv(tmp, "/scratch/avannan/MANUSCRIPT_figures/VUILD91MF_seuratmetadata_v2.csv")

vuild107mf <- subset(xenium, subset = sample == "VUILD107MF")
tmp <- vuild107mf@meta.data %>% filter(CNiche == "C3")
write.csv(tmp, "/scratch/avannan/MANUSCRIPT_figures/vuild107mf_seuratmetadata_v2.csv")


plt <- DimPlot(vuild107mf, reduction = "sp", group.by = "CNiche", cols = nuclei_niche_color_list) + coord_equal() + pretty_umap + NoLegend() +
  theme(plot.title = element_blank())
plt
pdf("/scratch/avannan/MANUSCRIPT_figures/test.pdf", width = (200*0.0393701), height = (200*0.0393701))
plt
dev.off()


## FIGURE 4 (OVERVIEW OF NICHE ANALYSIS) ----
lineage_df <- xenium@meta.data %>%
  select(final_CT, final_lineage) %>%
  unique()
names(lineage_df) <- c("CT", "Lineage")

#### FIGURE 4A (Example CNiche DimPlots) ----
vuhd113 <- subset(xenium, subset = sample == "VUHD113")
vuild107mf <- subset(xenium, subset = sample == "VUILD107MF")

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig3_VUHD113_CNiches.pdf", width = (150*0.0393701), height = (150*0.0393701))
DimPlot(vuhd113, reduction = "sp", group.by = "CNiche", cols = nuclei_niche_color_list,
        pt.size = 0.01) +
  coord_equal() + pretty_umap +
  theme(plot.title = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank()) +
  NoLegend()
dev.off()

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig3_VUILD107mf_CNiches.pdf", width = (150*0.0393701), height = (150*0.0393701))
DimPlot(vuild107mf, reduction = "sp", group.by = "CNiche", cols = nuclei_niche_color_list,
        pt.size = 0.01) +
  coord_equal() + pretty_umap +
  theme(plot.title = element_blank(), 
        axis.line = element_blank(),
        axis.title = element_blank()) +
  NoLegend()
dev.off()


#### FIGURE 4B (Niche x CT Dotplots - Cell Assignment to Niches) ----
# Transcript Niches
ct_tniche_prop_df <- as.data.frame((table(xenium$final_CT, xenium$TNiche)/
                                      rowSums(table(xenium$final_CT, xenium$TNiche))))
names(ct_tniche_prop_df) <- c("CT", "Niche", "Proportion")
ct_tniche_prop_plot <- ct_tniche_prop_df %>%
  filter(Niche != "TNA", Niche != "NA", !is.na(Niche)) %>%
  full_join(lineage_df) %>%
  mutate(CT = ordered(CT, levels = new_ct_order),
         Niche = ordered(Niche, levels = paste0("T", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  ggplot(aes(y = reorder(Niche, dplyr::desc(Niche)), x = CT,
             size = Proportion, color = as.factor(Niche),
             alpha = Proportion_Visible,
             fill = as.factor(Niche))) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw(base_size = 6) +
  scale_color_manual(values = c("black", "grey80", rep("black", 10)), guide = "none") +
  scale_fill_manual(values = transcript_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2.5), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 6),
        strip.text = element_text(color = "white"),
        legend.title = element_text(size = 6),
        legend.text = element_text(margin = margin(0, 0, 0, 0, "pt")),
        legend.spacing = unit(0, "pt"),
        legend.position = "bottom", legend.title.align = 0.5,
        legend.key.width = unit(3.5, "mm"),
        legend.key.height = unit(2, "mm")) +
  guides(size = guide_legend(title.position = "top", label.position = "bottom")) +
  # ggh4x::force_panelsizes(cols = c(2.1, 3.8, 4.1, 3.2)) +
  ggh4x::force_panelsizes(cols = c(1.5, 4.2, 8, 3.5)) +
  labs(title = "Cell Assignment to Niches") +
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  NoLegend()
ct_tniche_prop_plot

# Change facet colors
ct_tniche_prop_df_table <- ggplot_gtable(ggplot_build(ct_tniche_prop_plot))
striprt <- which(grepl("strip-r", ct_tniche_prop_df_table$layout$name) |
                   grepl("strip-t", ct_tniche_prop_df_table$layout$name))
fills <- rep(c("chocolate2", "chartreuse4", "maroon3", "cornflowerblue"), 2)
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", ct_tniche_prop_df_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  ct_tniche_prop_df_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

# Cell Niches
ct_niche_prop_df <- as.data.frame((table(xenium$final_CT, xenium$CNiche)/
                                     rowSums(table(xenium$final_CT, xenium$CNiche))))
names(ct_niche_prop_df) <- c("CT", "Niche", "Proportion")

ct_niche_prop_plot_with_legend <- ct_niche_prop_df %>%
  full_join(lineage_df) %>%
  mutate(CT = ordered(CT, levels = new_ct_order),
         Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  filter(!is.na(Niche)) %>%
  ggplot(aes(y = reorder(Niche, dplyr::desc(Niche)), x = CT, 
             size = Proportion, color = as.factor(Niche), 
             alpha = Proportion_Visible,
             fill = as.factor(Niche))) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw(base_size = 6) +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_color_manual(values = c(rep("black", 7), "grey80", rep("black", 4)), guide = "none") +
  scale_size_continuous(range = c(0, 2.5), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust = 1, size = 5),
        axis.title = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3), 
        strip.text = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        strip.background = element_blank(),
        axis.text.y = element_text(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        plot.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom") + 
  # ggh4x::force_panelsizes(cols = c(2.1, 3.8, 4.1, 3.2)) +
  ggh4x::force_panelsizes(cols = c(1.5, 4.2, 8, 3.5)) +
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  guides(size = guide_legend(keywidth = unit(1, "mm"), 
                             title.position = "top",
                             title.hjust = 0.5))
ct_niche_prop_plot <- ct_niche_prop_plot_with_legend + NoLegend()
ct_niche_prop_plot_with_legend

# Get legend
leg <- get_legend(ct_niche_prop_plot_with_legend)
a <- cowplot::plot_grid(ct_tniche_prop_df_table, ct_niche_prop_plot, ncol = 1, align = "hv")
cowplot::plot_grid(a, leg, ncol = 1, align = "hv")

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig3_cellassign2niches.pdf.pdf", width = (115*0.0393701), height = (220*0.0393701))
cowplot::plot_grid(a, leg, ncol = 1, align = "hv")
dev.off()



#### FIGURE 4C (Niche by Sample Type; Cells) ----
sample_type_tranx_niche_plot <- xenium@meta.data %>%
  select(TNiche, sample_affect) %>%
  mutate(TNiche = ordered(TNiche, levels = paste0("T", 1:12)),
         sample_affect = ordered(sample_affect, levels = c("Unaffected", "Less Affected", "More Affected"))) %>%
  filter(!is.na(TNiche)) %>%
  ggplot(aes(x = sample_affect, fill = TNiche)) +
  geom_bar(position = position_fill()) +
  scale_fill_manual(values = transcript_niche_color_list) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of Cells") +
  theme_classic(base_size = 5.5) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.key.size = unit(3, "mm")) +
  guides(fill = guide_legend(title.position = "top", ncol = 6,
                             title = "Transcript Niches", title.hjust = 0.5,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) +
  NoLegend()
sample_type_tranx_niche_plot

sample_type_cell_niche_plot <- xenium@meta.data %>%
  select(CNiche, sample_affect) %>%
  mutate(CNiche = ordered(CNiche, levels = paste0("C", 1:12)),
         sample_affect = ordered(sample_affect, levels = c("Unaffected", "Less Affected", "More Affected"))) %>%
  ggplot(aes(x = sample_affect, fill = CNiche)) +
  geom_bar(position = position_fill()) +
  scale_fill_manual(values = nuclei_niche_color_list) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of Cells") +
  theme_classic(base_size = 5.5) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.key.size = unit(3, "mm"),
        axis.title.y = element_text(color = NA)) +
  guides(fill = guide_legend(title.position = "top", ncol = 6,
                             title = "Cell Niches", title.hjust = 0.5,
                             label.theme = element_text(size = 5),
                             by_row = TRUE)) +
  NoLegend()
sample_type_cell_niche_plot

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig3_niche_overall_prop_sample_type.pdf",
    width = 67*0.0393701, height = 30*0.0393701)
ggarrange(sample_type_tranx_niche_plot, sample_type_cell_niche_plot, ncol = 2)
dev.off()


#### FIGURE 4D (Annotations x Niche Dotplots) ----
new_annotation_order <- c(# Epithelial
  "Normal Alveoli",
  "Minimally Remodeled Alveoli",
  "Hyperplastic AECs",
  "Emphysema",
  "Remodeled Epithelium",
  "Advanced Remodeling",
  "Epithelial Detachment",
  "Remnant Alveoli",
  "Small Airway",
  "Large Airway",
  "Microscopic Honeycombing",
  "Goblet Cell Metaplasia",
  
  # Immune
  "Granuloma",
  "Mixed Inflammation",
  "TLS",
  
  # Endothelial/Mesenchymal
  "Artery",
  "Muscularized Artery",
  "Venule",
  "Interlobular Septum",
  "Airway Smooth Muscle",
  "Fibroblastic Focus",
  "Fibrosis",
  "Severe Fibrosis",
  
  # Other
  "Multinucleated Cell",
  "Giant Cell")

# Relabel some things
new_annotation_order_4plot <- c("Normal Alveoli", "Min. Remod. Alveoli", "Hyperplastic AECs",
                                "Advanced Remodeling", "Epithelial Detachment", "Micro. Honeycombing",
                                "Goblet Cell Metaplasia", "Fibroblastic Focus", "Severe Fibrosis",
                                "Granuloma", "TLS")

# B1: Transcript Niches (Annotation Composition)
tniche_path_prop_df <- as.data.frame((table(final_annotations_df$Annotation_Type, final_annotations_df$TNiche)/
                                        rowSums(table(final_annotations_df$Annotation_Type, final_annotations_df$TNiche))))
names(tniche_path_prop_df) <- c("Pathology", "Niche", "Proportion")
tniche_path_prop_plot <- tniche_path_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("T", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  filter(Pathology %in% c("Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                          "Advanced Remodeling", "Epithelial Detachment", "Microscopic Honeycombing",
                          "Goblet Cell Metaplasia", "Granuloma", "TLS",
                          "Fibroblastic Focus", "Severe Fibrosis")) %>%
  mutate(Pathology = gsub("Minimally Remodeled", "Min. Remod.", gsub("Microscopic", "Micro.", Pathology)),
         Pathology = ordered(Pathology, levels = new_annotation_order_4plot)) %>%
  ggplot(aes(x = Niche, y = reorder(Pathology, dplyr::desc(Pathology)), 
             size = Proportion, alpha = Proportion_Visible, fill = Niche)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6.5) +
  scale_fill_manual(values = transcript_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2.5), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5), 
        axis.ticks.y = element_line(color = "black"),
        plot.margin = margin(0, 1, 0, 0, "pt"),
        legend.key = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0, "pt"),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.title = element_text(size = 6),
        legend.text = element_text(margin = margin(0, 0, 0, 0, "pt")),
        legend.spacing = unit(0, "pt"),
        legend.position = "bottom", legend.title.align = 0.5,
        legend.key.width = unit(3.5, "mm"),
        legend.key.height = unit(2, "mm")) +
  guides(size = guide_legend(title.position = "top", label.position = "bottom")) +
  scale_alpha_manual(values = c(0, 1), guide = "none")
tniche_path_prop_plot

# B2: Cell Niches (Annotation Composition)
cniche_path_prop_df <- as.data.frame((table(final_annotations_df$Annotation_Type, final_annotations_df$CNiche)/
                                        rowSums(table(final_annotations_df$Annotation_Type, final_annotations_df$CNiche))))
names(cniche_path_prop_df) <- c("Pathology", "Niche", "Proportion")
cniche_path_prop_plot <- cniche_path_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  filter(Pathology %in% c("Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                          "Advanced Remodeling", "Epithelial Detachment", "Microscopic Honeycombing",
                          "Goblet Cell Metaplasia", "Granuloma", "TLS",
                          "Fibroblastic Focus", "Severe Fibrosis")) %>%
  mutate(Pathology = gsub("Minimally Remodeled", "Min. Remod.", gsub("Microscopic", "Micro.", Pathology)),
         Pathology = ordered(Pathology, levels = new_annotation_order_4plot)) %>%
  ggplot(aes(x = Niche, y = reorder(Pathology, dplyr::desc(Pathology)), 
             size = Proportion, alpha = Proportion_Visible, fill = Niche)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6.5) +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2.5), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(0, 0, 0, 1, "pt"),
        legend.key = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0, "pt"),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.title = element_text(size = 6),
        legend.text = element_text(margin = margin(0, 0, 0, 0, "pt")),
        legend.spacing = unit(0, "pt"),
        legend.position = "bottom", legend.title.align = 0.5,
        legend.key.width = unit(3.5, "mm"),
        legend.key.height = unit(2, "mm")) +
  guides(size = guide_legend(title.position = "top", label.position = "bottom")) +
  scale_alpha_manual(values = c(0, 1), guide = "none")
cniche_path_prop_plot

# Combine and save
pdf("/scratch/avannan/MANUSCRIPT_figures/Fig3_niches_assign2anno.pdf", width = (63.5*0.0393701), height = (60*0.0393701))
(tniche_path_prop_plot + cniche_path_prop_plot) + 
  patchwork::plot_layout(widths = unit(c(21.5, 21.5), "mm"),
                         heights = unit(c(33.5, 33.5), "mm"),
                         ncol = 2, nrow = 1)
dev.off()


## FIGURE 5 (EPITHELIAL DETACHMENT) ----
#### FIGURE 5H (Niche x Sample Type Boxplots - T3/C7) ----
# T3 niche
tniche_prop_df <- (table(xenium$sample, xenium$TNiche)/rowSums(table(xenium$sample, xenium$TNiche))) %>%
  as.data.frame() %>%
  dplyr::rename(sample = "Var1", TNiche = "Var2", Prop = "Freq") %>%
  full_join(xenium@meta.data %>% select(sample, sample_affect) %>% unique()) %>%
  mutate(sample_affect = case_when(sample_affect == "Less Affected" ~ "Less",
                                   sample_affect == "More Affected" ~ "More",
                                   TRUE ~ sample_affect)) %>%
  mutate(sample_affect = ordered(sample_affect, levels = c("Unaffected", "Less", "More")))

t3_boxplot <- tniche_prop_df %>%
  filter(TNiche == "T3") %>%
  ggplot(aes(x = sample_affect, y = Prop, fill = sample_affect)) +
  geom_boxplot(width = 0.8, position = position_dodge(width = 1),
               outlier.shape = 21, color = "black", size = 0.1, 
               outlier.size = 0.1) +
  geom_point(size = 0.1) +
  scale_fill_manual(values = c(list(Unaffected = "white", `Less` = "grey85", `More` = "grey50"))) +
  labs(y = "Proportion of Cells") +
  theme_bw(base_size = 6) +
  facet_wrap(~TNiche) +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        # axis.text.x = element_text(size = 5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "white", size = 6),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey90"),
        plot.title = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        strip.background = element_rect(fill = "#8F7C00"),
        axis.title.y = element_text(size = 6),
        plot.margin = margin(0.5,2.5,0.5,0, "pt"),
        legend.title = element_blank()) +
  NoLegend()
t3_boxplot

# C3 niche
cniche_prop_df <- (table(xenium$sample, xenium$CNiche)/rowSums(table(xenium$sample, xenium$CNiche))) %>%
  as.data.frame() %>%
  dplyr::rename(sample = "Var1", CNiche = "Var2", Prop = "Freq") %>%
  full_join(xenium@meta.data %>% select(sample, sample_affect) %>% unique()) %>%
  mutate(sample_affect = case_when(sample_affect == "Less Affected" ~ "Less",
                                   sample_affect == "More Affected" ~ "More",
                                   sample_affect == "Unaffected" ~ "Unaff.")) %>%
  mutate(sample_affect = ordered(sample_affect, levels = c("Unaff.", "Less", "More")))

c3_boxplot <- cniche_prop_df %>%
  filter(CNiche == "C3") %>%
  ggplot(aes(x = sample_affect, y = Prop, fill = sample_affect)) +
  geom_boxplot(width = 0.8, position = position_dodge(width = 1),
               outlier.shape = 21, color = "black", size = 0.1, 
               outlier.size = 0.1) +
  geom_point(size = 0.1) +
  scale_fill_manual(values = c(list(Unaff. = "white", `Less` = "grey85", `More` = "grey50"))) +
  labs(y = "Proportion of Cells") +
  theme_bw(base_size = 6) +
  facet_wrap(~CNiche) +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "white", size = 6),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey90"),
        plot.title = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        strip.background = element_rect(fill = "#d10040"),
        axis.title.y = element_text(size = 6),
        plot.margin = margin(0.5,0,0.5,2.5, "pt"),
        legend.title = element_blank()) +
  NoLegend()
c3_boxplot

plt <- (t3_boxplot + c3_boxplot) + patchwork::plot_layout(ncol = 1)
pdf("/scratch/avannan/MANUSCRIPT_figures/Fig4_t3_c3_boxplots.pdf", width = 28.5*0.0393701, height = 51*0.0393701)
plt
dev.off()


#### FIGURE 5I: Cell Proximity by Niche; KRT5-/KRT17+ vs. Activated Fibrotic Fibroblasts ----
library(data.table)
library(dplyr)
library(EnrichedHeatmap)
library(circlize)
library(RColorBrewer)

'%!in%' <- function(x,y)!('%in%'(x,y))
filter <- dplyr::filter
select <- dplyr::select
#### load data and define code paramters
workDir <- '/scratch/avannan/late_IPF_spatial/xenium/REVISION_prox_analysis_av/Plots/'


logreg_cniche <- read.csv(file.path(workDir, 'LogReg_ProximityResults_cNiche.csv')) %>%
  mutate(Significance = ifelse(sig == "fdr < 0.05", "FDR < 0.05", "n.s."))
logreg_cniche$niche <- factor(logreg_cniche$niche, levels = paste0('C', 1:12))

p1_withlegend <- logreg_cniche %>%
  filter(target == 'Activated Fibrotic FBs') %>%
  filter(source == 'KRT5-/KRT17+') %>%
  ggplot(aes(x = niche, y = logOR, color = Significance)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = logOR_lo, ymax = logOR_hi), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme_black +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.title = element_blank(),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 6),
        legend.position = "bottom", legend.title = element_blank()) +
  labs(title = 'Proximity of KRT5-/KRT17+ to Activated Fibrotic FBs',
       x = NULL,  y = 'log(OR)') +
  scale_color_manual(values = c('purple', 'black')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2)
p1 <- p1_withlegend + theme(legend.position = "none")

p2 <- logreg_cniche %>%
  filter(target == 'KRT5-/KRT17+') %>%
  filter(source == 'Activated Fibrotic FBs') %>%
  ggplot(aes(x = niche, y = logOR, color = Significance)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = logOR_lo, ymax = logOR_hi), width = 0.5, linewidth = 0.3) +
  theme_bw(base_size = 6.5) +
  theme_black +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.margin = margin(0,0,0,0, "mm"),
        plot.title = element_blank(),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.title.y = element_text(size = 6),
        legend.position = "none", legend.title = element_blank()) +
  labs(title = 'Proximity of Activated Fibrotic FBs to KRT5-/KRT17+',
       x = NULL,  y = 'log(OR)')+
  scale_color_manual(values = c('purple', 'black')) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.2)

a <- p1 +
  patchwork::plot_annotation(title = "KRT5-/KRT17+ to Activated Fibrotic FBs",
                             theme = theme(plot.title = element_text(hjust = 0.5, size = 6,
                                                                     margin = margin(0,0,1,0, "mm"))))
b <- p2 +
  patchwork::plot_annotation(title = "Activated Fibrotic FBs to KRT5-/KRT17+",
                             theme = theme(plot.title = element_text(hjust = 0.5, size = 6,
                                                                     margin = margin(0,0,1,0, "mm"))))
leg <- cowplot::get_legend(p1_withlegend)
plots_withoutlegend <- cowplot::plot_grid(a, b, ncol = 1,  align = "hv", axis = "rl")
plots_withoutlegend

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig4_detachment_cell_proximity_KRT5_ActFB.pdf",
    width = 60.5*0.0393701, height = 45*0.0393701)
cowplot::plot_grid(plots_withoutlegend, leg, ncol = 1, nrow = 2, rel_heights = c(1, 0.1))
dev.off()


#### FIGURE 5J: Dotplot Heatmap of Epithelial vs. KRT5-/KRT17+ Expression ----
# Separate out epithelial cells
epithelial_cells <- subset(xenium, subset = final_lineage == "Epithelial")
epithelial_cells <- SetIdent(epithelial_cells, value = "final_CT")

# Add another variable for KRT5-/KRT17+ vs. all other epithelial
epithelial_cells$tmp_ct <- "Other Epithelial"
epithelial_cells$tmp_ct[epithelial_cells$final_CT == "KRT5-/KRT17+"] <- "KRT5-/KRT17+"
epithelial_cells <- SetIdent(epithelial_cells, value = "tmp_ct")

# Get markers for KRT5-/KRT17+ cells
krt5_markers <- FindMarkers(epithelial_cells, ident.1 = "KRT5-/KRT17+", 
                            ident.2 = "Other Epithelial", only.pos = TRUE) %>%
  mutate(pct.diff = abs(pct.1 - pct.2),
         fdr = p.adjust(p_val, method = "fdr")) %>%
  rownames_to_column(var = "Gene") %>%
  filter(fdr < 0.05)
krt5_markers_conservative <- krt5_markers %>% filter(avg_log2FC > 0.5, pct.diff > 0.1) %>% pull(Gene)
krt5_markers_conservative <- krt5_markers %>% filter(avg_log2FC > 0.5, pct.diff > 0.1) %>% pull(Gene)
krt5_genes <- krt5_markers %>% pull(Gene)

# Pick out genes for heatmap
krt5_genes_4heat <- c("COL1A1", "COL1A2", "COL3A1", "COL4A3", "ITGB1", "MMP7", "POSTN", "SOX9", "TGFB1", "TGFB2", # GO:0030198: extracellular matrix organization
                      "COL1A2", "COL3A1", "ITGB6", "TGFB1", "TGFB2", # GO:0007179: transforming growth factor beta receptor signaling pathway
                      "CTNNB1", "FN1", "ITGA3", "ITGAV", "ITGB1", "ITGB6", "MSLN", # GO:0007160: cell-matrix adhesion
                      "CEACAM6", "CTNNB1", "ICAM1", "ITGA3", "ITGAV", "ITGB1", "ITGB6", "KRT18", "SOX9", "TGFB2", # GO:0098609: cell-cell adhesion
                      "AXL", "CEACAM6", "CTNNB1", "COL3A1", "COL4A3", "FN1", "ICAM1", "ITGA3", "ITGAV", "ITGB1", 
                      "ITGB6", "KRT18", "MSLN", "POSTN", "SOX9", "TGFB2", # GO:0007155: cell adhesion
                      "AXL", "BAX", "CCL18", "CTNNB1", "COL3A1", "FN1", "ICAM1", "ITGA3", "ITGAV", "ITGB1", "ITGB6", "KDR", "TGFB1", "TGFB2", # GO:0016477: cell migration & GO:0048870: cell motility
                      "COL3A1", "FN1", "ITGA3", "ITGAV", "ITGB1", "ITGB6", # GO:0007229: integrin-mediated signaling pathway
                      "KRT14", "KRT18", "VIM", # GO:0045104: intermediate filament cytoskeleton organization
                      "BAX", "BCL2L1", "CTNNB1", "FAS", "ITGAV", "KRT18", "TGFB1", "TGFB2", # GO:0097190: apoptotic signaling pathway
                      "CTNNB1", "KDR", "SOX9", "SPINK1", "TGFB1", # GO:0050679: positive regulation of epithelial cell proliferation
                      "COL1A1", "SOX4", "TGFB1", # GO:0030177: positive regulation of canonical Wnt singaling pathway
                      "POSTN", "SOX4", "TGFB1", "TGFB2", "SOX4", # GO:0001666: response to hypoxia
                      "SOX4", "SOX9") # PTHR10270: SOX transcription factor
krt5_genes_4heat %in% krt5_markers_conservative
krt5_genes_4heat <- sort(unique(krt5_genes_4heat))

# Cell categories
extracellular_mat <- c("COL1A1", "COL1A2", "COL3A1", "COL4A3", "ITGB1", "MMP7", "POSTN", "SOX9", "TGFB1", "TGFB2") # GO:0030198: extracellular matrix organization
tgfb_sig <- c("COL1A2", "COL3A1", "ITGB6", "TGFB1", "TGFB2") # GO:0007179: transforming growth factor beta receptor signaling pathway
cell_adhesion <- c("CTNNB1", "FN1", "ITGA3", "ITGAV", "ITGB1", "ITGB6", "MSLN", # GO:0007160: cell-matrix adhesion
                   "CEACAM6", "CTNNB1", "ICAM1", "ITGA3", "ITGAV", "ITGB1", "ITGB6", "KRT18", "SOX9", "TGFB2", # GO:0098609: cell-cell adhesion
                   "AXL", "CEACAM6", "CTNNB1", "COL3A1", "COL4A3", "FN1", "ICAM1", "ITGA3", "ITGAV", "ITGB1", 
                   "ITGB6", "KRT18", "MSLN", "POSTN", "SOX9", "TGFB2", # GO:0007155: cell adhesion
                   "AXL", "BAX", "CCL18", "CTNNB1", "COL3A1", "FN1", "ICAM1", "ITGA3", "ITGAV", "ITGB1", "ITGB6", "KDR", "TGFB1", "TGFB2") # GO:0016477: cell migration & GO:0048870: cell motility
integrin_sig <- c("COL3A1", "FN1", "ITGA3", "ITGAV", "ITGB1", "ITGB6") # GO:0007229: integrin-mediated signaling pathway
cytoskeleton <- c("KRT14", "KRT18", "VIM") # GO:0045104: intermediate filament cytoskeleton organization
apoptosis <- c("BAX", "BCL2L1", "CTNNB1", "FAS", "ITGAV", "KRT18", "TGFB1", "TGFB2") # GO:0097190: apoptotic signaling pathway
hypoxia <- c("POSTN", "SOX4", "TGFB1", "TGFB2", "SOX4") # GO:0001666: response to hypoxia
sox <- c("SOX4", "SOX9") # PTHR10270: SOX transcription factor
epi_prolif <- c("CTNNB1", "KDR", "SOX9", "SPINK1", "TGFB1") # GO:0050679: positive regulation of epithelial cell proliferation

# Make sure all genes have at least 1 category
krt5_genes_4heat %in% unique(c(extracellular_mat, tgfb_sig, cell_adhesion, integrin_sig, cytoskeleton, apoptosis, hypoxia, sox, epi_prolif))

# Log-normalizing and scaling all features in the RNA assay
# Scaling so that all features can be visualized using the same color scale
epithelial_cells <- ScaleData(epithelial_cells)

# Create dotplot base
p <- DotPlot(epithelial_cells, features = sort(unique(krt5_genes_4heat)), 
             group.by = "final_CT", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p

# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  tidyr::pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# Replace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 2 will be mapped to bright red
col_fun = circlize::colorRamp2(c(-1, 0, 2), colorspace::diverge_hsv(3))

# Adding annotation colors
library(scales) 
epi_color_list <- color_list[c("AT1", "Transitional AT2", "AT2",
                               "Proliferating AT2", "KRT5-/KRT17+",
                               "Basal", "Multiciliated", "Goblet", 
                               "RASC", "Secretory", "PNEC",
                               "Proliferating Airway")]
col_for_plot <- epi_color_list
col_for_plot_heatmap <- list(cluster = unlist(epi_color_list))
exp_mat <- exp_mat[, names(epi_color_list)]
percent_mat <- percent_mat[, names(epi_color_list)]

# Creating a layer to add to plot
layer_fun = function(j, i, x, y, w, h, fill){
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex(exp_mat, i, j))), # t(exp_mat)
              size = pindex(percent_mat, i, j)/100 * unit(2, "mm"), # t(percent_mat)
              pch = 19)
}

# List of legends
lgd_list = list(
  Legend(labels = c(0,0.25,0.5,0.75,1), title = "Proportion", type = "points", 
         pch = 19, size = c(0,0.25,0.5,0.75,1) * unit(2, "mm"),
         legend_gp = gpar(col = "black"), direction = "horizontal", ncol = 5,
         title_position = "topcenter", background = NA, grid_height = unit(1, "mm"),
         grid_width = unit(1, "mm"),
         title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5),
         legend_height = unit(20, "mm")))

# Create heatmap
row_ha = rowAnnotation("Adhesion/Motility" = krt5_genes_4heat %in% cell_adhesion,
                       "Extracellular Matrix" = krt5_genes_4heat %in% extracellular_mat,
                       "TGF-beta Signaling" = krt5_genes_4heat %in% tgfb_sig,
                       "Integrin Signaling" = krt5_genes_4heat %in% integrin_sig,
                       "Cytoskeleton" = krt5_genes_4heat %in% cytoskeleton,
                       "Apoptosis" = krt5_genes_4heat %in% apoptosis,
                       "Hypoxia" = krt5_genes_4heat %in% hypoxia,
                       "Epithelial Proliferation" = krt5_genes_4heat %in% epi_prolif,
                       "SOX TFs" = krt5_genes_4heat %in% sox,
                       col = list("Extracellular Matrix" = c(`TRUE` = "black", `FALSE` = "white"),
                                  "TGF-beta Signaling" = c(`TRUE` = "black", `FALSE` = "white"),
                                  "Adhesion/Motility" = c(`TRUE` = "black", `FALSE` = "white"),
                                  "Integrin Signaling" = c(`TRUE` = "black", `FALSE` = "white"),
                                  "Cytoskeleton" = c(`TRUE` = "black", `FALSE` = "white"),
                                  "Apoptosis" = c(`TRUE` = "black", `FALSE` = "white"),
                                  "Hypoxia" = c(`TRUE` = "black", `FALSE` = "white"),
                                  "SOX TFs" = c(`TRUE` = "black", `FALSE` = "white"),
                                  "Epithelial Proliferation" = c(`TRUE` = "black", `FALSE` = "white")),
                       simple_anno_size = unit(1.9, "mm"),
                       annotation_name_gp = gpar(fontsize = 5),
                       annotation_name_rot = 45,
                       show_legend = FALSE)

colnames(exp_mat) <- c("AT1", "Transitional AT2", "AT2",
                       "Proliferating AT2", "KRT5-/KRT17+",
                       "Basal", "Multiciliated", "Goblet", 
                       "RASC", "Secretory", "PNEC",
                       "Proliferating Airway")
hp <- Heatmap(exp_mat, # t(exp_mat)
              name = "Scaled Expression",
              heatmap_legend_param = list(title = "Scaled",
                                          labels = c("-1", "0", "1", "2+"),
                                          legend_direction = "horizontal",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_width = unit(1, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(20, "mm")),
              column_title = " ", 
              row_title = " ",
              row_title_gp = gpar(fontsize = 6),
              column_title_gp = gpar(fontsize = 5),
              row_title_side = "right",
              column_title_side = "top",
              column_names_rot = 45,
              # row_names_rot = 45,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 5),
              row_names_side = "left",
              column_names_gp = gpar(fontsize = 5),
              cluster_rows = FALSE,  cluster_columns = FALSE,
              column_dend_height = unit(1, "mm"),
              border = "black",
              
              left_annotation = row_ha)

# Add annotations to heatmap
ht <- draw(hp, 
           heatmap_legend_list = lgd_list,
           heatmap_legend_side = "bottom",
           align_heatmap_legend = "global_center")
ht

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig4_krt5_expr.pdf", width = 59*0.0393701, height = 69*0.0393701)
ht <- draw(hp, 
           legend_gap = unit(2, "mm"),
           padding = margin(0,0,0,0, "mm"),
           heatmap_legend_list = lgd_list,
           heatmap_legend_side = "bottom",
           align_heatmap_legend = "global_center")
dev.off()



# FIGURE 6 (Macrophage Accumulation) ----
#### FIGURE 6A (Niche x Sample Type Boxplots - C11) ----
# C11 niche
cniche_prop_df <- (table(xenium$sample, xenium$CNiche)/rowSums(table(xenium$sample, xenium$CNiche))) %>%
  as.data.frame() %>%
  dplyr::rename(sample = "Var1", CNiche = "Var2", Prop = "Freq") %>%
  full_join(xenium@meta.data %>% select(sample, sample_affect) %>% unique()) %>%
  mutate(sample_affect = case_when(sample_affect == "Less Affected" ~ "Less",
                                   sample_affect == "More Affected" ~ "More",
                                   TRUE ~ sample_affect)) %>%
  mutate(sample_affect = ordered(sample_affect, levels = c("Unaffected", "Less", "More")))

c11_boxplot <- cniche_prop_df %>%
  filter(CNiche == "C11") %>%
  ggplot(aes(x = sample_affect, y = Prop, fill = sample_affect)) +
  geom_boxplot(width = 0.8, position = position_dodge(width = 1.5),
               outlier.shape = 21, color = "black", size = 0.3, 
               outlier.size = 0.1) +
  geom_point(size = 0.1) +
  scale_fill_manual(values = c(list(Unaffected = "white", `Less` = "grey85", `More` = "grey50"))) +
  labs(y = "Proportion of Cells") +
  theme_bw(base_size = 6) +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 5),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "white", size = 6.5),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey90"),
        plot.title = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        strip.background = element_rect(fill = "#72c7ff"),
        axis.title.y = element_text(size = 6),
        legend.title = element_blank()) +
  facet_wrap(~CNiche) +
  NoLegend()
c11_boxplot

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig5_c11_boxplot.pdf", width = 38*0.0393701, height = 27*0.0393701)
c11_boxplot
dev.off()


#### FIGURE 6B (Example Sample - VUILD102LF) ----
vuild102lf <- subset(xenium, subset = sample == "VUILD102LF")
plt <- DimPlot(vuild102lf, group.by = "CNiche", cols = nuclei_niche_color_list, reduction = "sp", pt.size = 0.001) + coord_equal() + 
  pretty_umap + theme(axis.line = element_blank(), axis.title = element_blank(), plot.title = element_blank()) + NoLegend()
FeaturePlot(vuild102lf, features = c("CTHRC1", "FAP", "POSTN"), reduction = "sp", order = TRUE, ncol = 3) & coord_equal()

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig5_vuild102lf_cniches.pdf", width = 140*0.0393701, height = 140*0.0393701)
plt
dev.off()


#### FIGURE 6F (Proportion of Macrophages) ----
mac_color_list <- color_list[which(names(color_list) %in% c("Alveolar Macrophages", "Interstitial Macrophages",
                                                            "SPP1+ Macrophages", "Macrophages - IFN-activated",
                                                            "Monocytes/MDMs", "Proliferating Myeloid"))]
names(mac_color_list) <- gsub("Macrophages", "MΦ", names(mac_color_list))
macro_prop_plot <- xenium@meta.data %>%
  filter(final_CT %in% c("Alveolar Macrophages", "Interstitial Macrophages",
                         "SPP1+ Macrophages", "Macrophages - IFN-activated",
                         "Monocytes/MDMs", "Proliferating Myeloid")) %>%
  mutate(final_CT = gsub("Macrophages", "MΦ", final_CT),
         final_CT = ordered(final_CT, levels = names(mac_color_list))) %>%
  mutate(sample_affect = ordered(case_when(sample_affect == "Less Affected" ~ "Less",
                                           sample_affect == "More Affected" ~ "More",
                                           TRUE ~ sample_affect),
                                 levels = c("Unaffected", "Less", "More"))) %>%
  ggplot(aes(x = sample_affect, fill = final_CT)) +
  geom_bar(position = position_fill()) +
  scale_fill_manual(values = mac_color_list) +
  theme_classic(base_size = 6.5) +
  theme_black +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of MΦ") +
  theme(axis.title.x = element_blank(),
        legend.key.size = unit(2, "mm"),
        axis.title.y = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(0,0,0,0, "mm"),
        legend.box.spacing = unit(0, "mm"))
macro_prop_plot

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig5_macro_prop.pdf", width = 52*0.0393701, height = 27*0.0393701)
macro_prop_plot
dev.off()


## FIGURE 7 (LUMEN ANALYSIS) ----
#### FIGURE 7B (Heatmap of Lumens - Partial, w/ Gene Line Plots) ----
# Possible labels to add to heatmap
trans_tniche_lumen_metadata <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/trans_tniche_lumen_metadata_080724.csv") %>%
  left_join(xenium@meta.data %>% select(sample, sample_affect, percent_pathology) %>% unique()) %>%
  arrange(pseudotime_rank_trans_t4_tie_avg)
pseudotime_value <- trans_tniche_lumen_metadata$pseudotime_rank_trans_t4_tie_avg
sample_affect_value <- trans_tniche_lumen_metadata$sample_affect
path_scores <- trans_tniche_lumen_metadata$percent_pathology


ct_props <- xenium@meta.data %>%
  filter(!is.na(lumen_id)) %>%
  group_by(lumen_id) %>%
  mutate(cell_count = length(final_CT)) %>%
  group_by(lumen_id, final_CT) %>%
  mutate(celltype_count = length(final_CT),
         celltype_prop = celltype_count/cell_count) %>%
  select(sample, lumen_id, cell_id, final_CT, celltype_prop) %>%
  select(lumen_id, final_CT, celltype_prop) %>%
  unique() %>%
  pivot_wider(names_from = "final_CT", values_from = "celltype_prop") %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) 
ct_mat <- as.matrix(ct_props[, -1])
rownames(ct_mat) <- ct_props$lumen_id
ct_mat <- ct_mat[trans_tniche_lumen_metadata$lumen_id, ]


top_annotation =
  HeatmapAnnotation(
    AT1 = ct_mat[, "AT1"],
    Capillary = ct_mat[, "Capillary"],
    #`Inflammatory FBs` = ct_mat[, "Inflammatory FBs"],
    #`Monocytes/MDMs` = ct_mat[, "Monocytes/MDMs"],
    #`NK/NKT` = ct_mat[, "NK/NKT"],
    #Neutrophils = ct_mat[, "Neutrophils"],
    #Arteriole = ct_mat[, "Arteriole"],
    
    #`CD8+ T-cells` = ct_mat[, "CD8+ T-cells"],
    #Tregs = ct_mat[, "Tregs"],
    #`Macrophages - IFN-activated` = ct_mat[, "Macrophages - IFN-activated"],
    #`CD4+ T-cells` = ct_mat[, "CD4+ T-cells"],
    #Plasma = ct_mat[, "Plasma"],
    `Alveolar Mϕ` = ct_mat[, "Alveolar Macrophages"],
    #`Migratory DCs` = ct_mat[, "Migratory DCs"],
    #`Interstitial Macrophages` = ct_mat[, "Interstitial Macrophages"],
    #Mast = ct_mat[, "Mast"],
    #`Proliferating AT2` = ct_mat[, "Proliferating AT2"],
    `Activ. Fibrotic FBs` = ct_mat[, "Activated Fibrotic FBs"],
    #Venous = ct_mat[, "Venous"],
    #cDCs = ct_mat[, "cDCs"],
    #Multiciliated = ct_mat[, "Multiciliated"],
    #Myofibroblasts = ct_mat[, "Myofibroblasts"],
    
    `SPP1+ Mϕ` = ct_mat[, "SPP1+ Macrophages"],
    `Transitional AT2` = ct_mat[, "Transitional AT2"],
    #RASC = ct_mat[, "RASC"],
    #Secretory = ct_mat[, "Secretory"],
    #Basal = ct_mat[, "Basal"],
    `KRT5-/KRT17+` = ct_mat[, "KRT5-/KRT17+"],
    
    C8 = trans_tniche_lumen_metadata$C8,
    #C5 = trans_tniche_lumen_metadata$C5,
    #C12 = trans_tniche_lumen_metadata$C12,
    C3 = trans_tniche_lumen_metadata$C3,
    C4 = trans_tniche_lumen_metadata$C4,
    C11 = trans_tniche_lumen_metadata$C11,
    C2 = trans_tniche_lumen_metadata$C2,
    
    T4 = trans_tniche_lumen_metadata$T4,
    #T10 = trans_tniche_lumen_metadata$T10,
    #T2 = trans_tniche_lumen_metadata$T2,
    #T9 = trans_tniche_lumen_metadata$T9,
    T11 = trans_tniche_lumen_metadata$T11,
    #T5 = trans_tniche_lumen_metadata$T5,
    T6 = trans_tniche_lumen_metadata$T6,
    #T12 = trans_tniche_lumen_metadata$T12,
    T8 = trans_tniche_lumen_metadata$T8,
    T3 = trans_tniche_lumen_metadata$T3,
    
    `Disease Severity` = sample_affect_value,
    `% Pathology` = path_scores,
    Pseudotime = pseudotime_value,
    
    # lumen_yes = trans_tniche_lumen_metadata$lumen_yes,
    
    simple_anno_size = unit(1.6, "mm"),
    height = unit(1.6, "mm"),
    annotation_name_gp = gpar(fontsize = 5),
    gap = unit(0, "mm"),
    annotation_name_side = "left",
    
    col = list(
      AT1 = colorRamp2(c(0, max(ct_mat[, "AT1"])), c("white", "chartreuse4")),
      Capillary = colorRamp2(c(0, max(ct_mat[, "Capillary"])), c("white", "chocolate2")),
      `Inflammatory FBs` = colorRamp2(c(0, max(ct_mat[, "Inflammatory FBs"])), c("white", "cornflowerblue")),
      `Monocytes/MDMs` = colorRamp2(c(0, max(ct_mat[, "Monocytes/MDMs"])), c("white", "deeppink3")),
      `NK/NKT` = colorRamp2(c(0, max(ct_mat[, "NK/NKT"])), c("white", "darkmagenta")),
      Neutrophils = colorRamp2(c(0, max(ct_mat[, "Neutrophils"])), c("white", "deeppink3")),
      Arteriole = colorRamp2(c(0, max(ct_mat[, "Arteriole"])), c("white", "chocolate2")),
      
      `CD8+ T-cells` = colorRamp2(c(0, max(ct_mat[, "CD8+ T-cells"])), c("white", "darkmagenta")),
      Tregs = colorRamp2(c(0, max(ct_mat[, "Tregs"])), c("white", "darkmagenta")),
      `Macrophages - IFN-activated` = colorRamp2(c(0, max(ct_mat[, "Macrophages - IFN-activated"])), c("white", "deeppink3")),
      `CD4+ T-cells` = colorRamp2(c(0, max(ct_mat[, "CD4+ T-cells"])), c("white", "darkmagenta")),
      Plasma = colorRamp2(c(0, max(ct_mat[, "Plasma"])), c("white", "darkmagenta")),
      `Alveolar Mϕ` = colorRamp2(c(0, max(ct_mat[, "Alveolar Macrophages"])), c("white", "deeppink3")),
      `Migratory DCs` = colorRamp2(c(0, max(ct_mat[, "Migratory DCs"])), c("white", "deeppink3")),
      `Interstitial Macrophages` = colorRamp2(c(0, max(ct_mat[, "Interstitial Macrophages"])), c("white", "deeppink3")),
      Mast = colorRamp2(c(0, max(ct_mat[, "Mast"])), c("white", "deeppink3")),
      `Proliferating AT2` = colorRamp2(c(0, max(ct_mat[, "Proliferating AT2"])), c("white", "chartreuse4")),
      `Activ. Fibrotic FBs` = colorRamp2(c(0, max(ct_mat[, "Activated Fibrotic FBs"])), c("white", "cornflowerblue")),
      Venous = colorRamp2(c(0, max(ct_mat[, "Venous"])), c("white", "chocolate2")),
      cDCs = colorRamp2(c(0, max(ct_mat[, "cDCs"])), c("white", "deeppink3")),
      Multiciliated = colorRamp2(c(0, max(ct_mat[, "Multiciliated"])), c("white", "chartreuse4")),
      Myofibroblasts = colorRamp2(c(0, max(ct_mat[, "Myofibroblasts"])), c("white", "cornflowerblue")),
      
      `SPP1+ Mϕ` = colorRamp2(c(0, max(ct_mat[, "SPP1+ Macrophages"])), c("white", "deeppink3")),
      `Transitional AT2` = colorRamp2(c(0, max(ct_mat[, "Transitional AT2"])), c("white", "chartreuse4")),
      RASC = colorRamp2(c(0, max(ct_mat[, "RASC"])), c("white", "chartreuse4")),
      Secretory = colorRamp2(c(0, max(ct_mat[, "Secretory"])), c("white", "chartreuse4")),
      Basal = colorRamp2(c(0, max(ct_mat[, "Basal"])), c("white", "chartreuse4")),
      `KRT5-/KRT17+` = colorRamp2(c(0, max(ct_mat[, "KRT5-/KRT17+"])), c("white", "chartreuse4")),
      
      # lumen_yes = c(`1` = "red", `0` = "white"),
      
      C1 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C1)), c("white", "#1036b5")),
      C2 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C2)), c("white", "#6f774b")),
      C3 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C3)), c("white", "#d10040")),
      C4 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C4)), c("white", "#84cd5d")),
      C5 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C5)), c("white", "#8E76DB")),
      C6 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C6)), c("white", "#71c8a5")),
      C7 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C7)), c("white", "#ffa26d")),
      C8 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C8)), c("white", "#924373")),
      C9 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C9)), c("white", "#a6513c")),
      C10 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C10)), c("white", "black")),
      C11 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C11)), c("white", "#72c7ff")),
      C12 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$C12)), c("white", "#cabd00")),
      
      T1 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T1)), c("white", "#94FFB5")),
      T2 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T2)), c("white", "#191919")),
      T3 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T3)), c("white", "#8F7C00")),
      T4 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T4)), c("white", "#FFCC99")),
      T5 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T5)), c("white", "#2BCE48")),
      T6 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T6)), c("white", "#993F00")),
      T7 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T7)), c("white", "#005C31")),
      T8 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T8)), c("white", "#0075DC")),
      T9 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T9)), c("white", "#F0A0FF")),
      T10 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T10)), c("white", "#C20088")),
      T11 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T11)), c("white", "#4C005C")),
      T12 = colorRamp2(c(0, max(trans_tniche_lumen_metadata$T12)), c("white", "#9DCC00")),
      
      `Disease Severity` = c(`Unaffected` = "#88CCEE", `Less Affected` = "pink", `More Affected` = "#ac3e6f"),
      `% Pathology` = colorRamp2(c(min(path_scores), max(path_scores)), c("white", "black")),
      Pseudotime = colorRamp2(c(min(pseudotime_value), max(pseudotime_value)), c("white", "black"))),
    
    show_legend = FALSE)

ps_ht <- Heatmap(matrix(nc = nrow(trans_tniche_lumen_metadata), nr = 0),
                 top_annotation = top_annotation)
ps_ht

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig6b_top_heatmap_anno.pdf", width = (52.1*0.0393701), height = (40*0.0393701))
ps_ht
dev.off()

# Gene line plots
gene_order_df <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/specc_cluster_offset_total_tx_tma_corrected_change_time_in_chunk_4c_line_plot_df.csv")
gene_line_plot <- gene_order_df %>%
  mutate(sc_clusters = ordered(case_when(sc_clusters == "4" ~ "Homeostasis",
                                         sc_clusters == "3" ~ "Early\nRemodeling",
                                         sc_clusters == "1" ~ "Intermediate\nRemodeling",
                                         sc_clusters == "2" ~ "Late\nRemodeling"), 
                               levels = c("Homeostasis", "Early\nRemodeling", "Intermediate\nRemodeling", "Late\nRemodeling"))) %>%
  ggplot(aes(x = time_order, y = value, group = gene)) +
  geom_line(alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  geom_smooth(aes(group = sc_clusters, color = sc_clusters), linewidth = 1) +
  theme_bw(base_size = 6.5) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 5.5),
        axis.title.y = element_text(vjust = -19),
        strip.text = element_text(size = 5.5),
        strip.switch.pad.grid = unit(3.8, "mm"),
        strip.placement = "outside"
  ) +
  theme_black +
  labs(x = "Pseudotime", y = "Gene Expression Z-Score") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("#04846E", "#FFC107", "#FF6A0F", "#D81B60")) +
  facet_grid(sc_clusters~., scales = "free_y", switch = "y") +
  NoLegend()
gene_line_plot

# Change facet colors
gene_line_plot_table <- ggplot_gtable(ggplot_build(gene_line_plot))
striprt <- which(grepl("strip-l", gene_line_plot_table$layout$name) | 
                   grepl("strip-t", gene_line_plot_table$layout$name))
# fills <- c("#D81B60", "#FF6A0F", "#FFC107", "#04846E")
# font_colors <- c("white", "white", "black", "white")
fills <- c("#04846E", "#FFC107", "#FF6A0F", "#D81B60")
font_colors <- c("white", "black", "white", "white")
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", gene_line_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  gene_line_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  gene_line_plot_table$grobs[[i]]$grobs[[1]]$children[[j+1]]$children[[1]]$gp$col <- font_colors[k]
  k <- k+1
}
grid::grid.draw(gene_line_plot_table)

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig6b_bottom_gene_lineplot.pdf", width = (52.2*0.0393701), height = (63.4*0.0393701))
grid::grid.draw(gene_line_plot_table)
dev.off()


#### FIGURE 7C (Heatmap of Homeostatic and Early Remodeling Expression) ----
# Load in lumen genes and their order
gene_mode_and_sc_clusters <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/gene_mode_and_sc_clusters.csv")

# Load in significant terms x CT and filter
library(googlesheets4)
gs4_deauth()
sig_terms_perCT_sheet <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094")
sig_terms_perCT_sheet <- read_sheet(sig_terms_perCT_sheet, sheet = 17, skip = 1)
sig_terms_perCT <- as.data.frame(sig_terms_perCT_sheet) %>%
  filter(padj < 0.05) %>%
  full_join(gene_mode_and_sc_clusters, by = c("gene" = "genes")) %>%
  mutate(total_sig = length(gene)) %>%
  group_by(celltype) %>%
  mutate(prop_sig_ct = length(gene)/total_sig) %>%
  ungroup()

# Load in lumen metadata and list of final lumens
lumen_nuclei_RNA_sce <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/lumen_nuclei_RNA_sce_t4_ordered.rds")
lumen_metadata <- xenium@meta.data %>%
  filter(lumen_id %in% colnames(lumen_nuclei_RNA_sce))

# Subset object to contain only cells found in final lumens
xenium_lumens_cells_only <- subset(xenium, cells = rownames(lumen_metadata))

# Pick out genes for heatmap
heatmap_genes <- sig_terms_perCT %>%
  filter(sc_clusters %in% c(4, 3)) %>%
  arrange(desc(sc_clusters), gene) %>%
  pull(gene) %>%
  unique()
heatmap4_genes <- sig_terms_perCT %>%
  filter(sc_clusters == 4) %>%
  arrange(desc(sc_clusters), gene) %>%
  pull(gene) %>%
  unique()
heatmap3_genes <- sig_terms_perCT %>%
  filter(sc_clusters == 3) %>%
  arrange(desc(sc_clusters), gene) %>%
  pull(gene) %>%
  unique()

# Log-normalizing and scaling all features in the RNA assay
# Scaling so that all features can be visualized using the same color scale
xenium_lumens_cells_only <- ScaleData(xenium_lumens_cells_only)

# Create dotplot base
p <- DotPlot(xenium_lumens_cells_only, features = heatmap_genes, 
             group.by = "final_CT", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p

# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()
dim(exp_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# Replace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Create dataframe with significance per lineage, as 1/0; use this as row annotation
# TRUE = Significant across pseudotime for at least one CT in that lineage
gene_lineage_lumen_heatmap_df_og <- sig_terms_perCT %>%
  dplyr::rename(final_CT = "celltype") %>%
  mutate(final_CT = gsub("\\.", " ", final_CT),
         final_CT = case_when(final_CT == "KRT5  KRT17 " ~ "KRT5-/KRT17+",
                              final_CT == "SPP1  Macrophages" ~ "SPP1+ Macrophages",
                              final_CT == "Monocytes MDMs" ~ "Monocytes/MDMs" ,
                              final_CT == "Macrophages   IFN activated" ~ "Macrophages - IFN-activated",
                              final_CT == "SMCs Pericytes" ~ "SMCs/Pericytes",
                              final_CT == "NK NKT" ~ "NK/NKT",
                              final_CT == "CD4  T cells" ~ "CD4+ T-cells",
                              final_CT == "CD8  T cells" ~ "CD8+ T-cells",
                              TRUE ~ final_CT)) %>%
  left_join(xenium@meta.data %>% select(final_CT, final_lineage, final_sublineage) %>% unique()) %>%
  mutate(lineage = ifelse(final_lineage == "Epithelial", "Epithelial", final_sublineage)) %>%
  mutate(sig = 1) %>%
  filter(!is.na(final_CT)) %>%
  filter(!is.na(gene)) %>%
  unique()

# Get proportion/percentage of significant tests
prop_gene_lumen_lineage_df <- gene_lineage_lumen_heatmap_df_og %>%
  group_by(sc_clusters) %>%
  mutate(num_sig_tests = sum(sig)) %>%
  ungroup() %>%
  group_by(lineage, sc_clusters) %>%
  mutate(num_sig_lineage = sum(sig),
         prop_sig_lineage = sum(sig)/num_sig_tests) %>%
  select(final_lineage, prop_sig_lineage, num_sig_lineage, num_sig_tests) %>%
  unique() %>%
  mutate(title = paste0(num_sig_lineage, "/", num_sig_tests, " (", round(prop_sig_lineage*100, 1), "%)"))
prop_gene_lumen_lineage_df

# Prep dataframe for heatmap
gene_lineage_lumen_heatmap_df <- gene_lineage_lumen_heatmap_df_og %>%
  select(gene, sc_clusters, lineage, sig) %>%
  unique() %>% 
  pivot_wider(names_from = "lineage", values_from = "sig") %>%
  arrange(desc(sc_clusters), gene) %>% 
  filter(!is.na(gene))
# mutate(across(everything(), ~ ifelse(.x == 0, NA, .x)))
gene_lineage_lumen_heatmap_df2 <- gene_lineage_lumen_heatmap_df %>%
  select(Endothelial, Epithelial, Lymphoid, Myeloid, Mesenchymal) %>%
  as.data.frame()
# mutate(across(everything(), ~ ifelse(.x == 0, NA, .x)))
rownames(gene_lineage_lumen_heatmap_df2) <- gene_lineage_lumen_heatmap_df$gene
# Including genes not significant for any lineage
nonsig_genes <- heatmap_genes[!(heatmap_genes %in% gene_lineage_lumen_heatmap_df$gene)]
nonsig_gene_lineage_lumen_heatmap_df <- data.frame(matrix(NA, nrow = length(nonsig_genes), ncol = 5))
rownames(nonsig_gene_lineage_lumen_heatmap_df) <- nonsig_genes
colnames(nonsig_gene_lineage_lumen_heatmap_df) <- colnames(gene_lineage_lumen_heatmap_df2)

# Join together
gene_lineage_lumen_heatmap_df3 <- rbind(gene_lineage_lumen_heatmap_df2, 
                                        nonsig_gene_lineage_lumen_heatmap_df)[rownames(exp_mat), ]

# Create heatmap
# New CT order and color list
new_ct_order3 <- c(
  # Endothelial 
  "Arteriole", "Capillary", "Venous", "Lymphatic",
  # Epithelial - Alveolar
  "AT1", "Transitional AT2", "AT2", "Proliferating AT2", "KRT5-/KRT17+",              
  # Epithelial - Airway
  "Basal", "Multiciliated", "Goblet", "RASC", "Secretory", "PNEC", "Proliferating Airway",
  # Immune - Lymphoid
  "B cells", 
  # "Proliferating B cells", # Not tested
  "NK/NKT", "Proliferating NK/NKT", "Tregs", 
  "CD4+ T-cells", "CD8+ T-cells", "Proliferating T-cells", "Plasma", "pDCs",
  # Immune - Myeloid
  "Basophils", "cDCs", "Migratory DCs", "Langerhans cells", "Mast", 
  "Alveolar Macrophages", "Interstitial Macrophages", "SPP1+ Macrophages",
  "Macrophages - IFN-activated", "Monocytes/MDMs", "Neutrophils", "Proliferating Myeloid",
  # Mesenchymal
  "Activated Fibrotic FBs", "Adventitial FBs", "Alveolar FBs", "Inflammatory FBs", 
  "Myofibroblasts", "Subpleural FBs", "Proliferating FBs", "SMCs/Pericytes", "Mesothelial")
new_color_list3 <- color_list[new_ct_order3]

# Get row split variable (clusters)
sc_cluster_order_genes <- sig_terms_perCT %>%
  filter(gene %in% rownames(exp_mat)) %>%
  select(gene, sc_clusters) %>%
  unique() %>%
  arrange(desc(sc_clusters), gene) %>%
  as.vector()
all.equal(sc_cluster_order_genes$gene, rownames(exp_mat)) # Must be TRUE
rowsplit <- ordered(as.factor(sc_cluster_order_genes$sc_clusters), levels = c("4", "3"))

ht_opt$TITLE_PADDING = unit(1.2, "mm")
hp <- Heatmap(exp_mat[, new_ct_order3], # t(exp_mat)
              name = "Scaled Expression",
              width = unit(49, "mm"),
              heatmap_legend_param = list(title = "Scaled Expression",
                                          legend_direction = "horizontal",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_width = unit(2, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),
              column_title_side = "top",
              row_title = " ",
              column_title_gp = gpar(fontsize = 5.5, col = c("chocolate2", "chartreuse4",
                                                             "darkmagenta", "deeppink3", "cornflowerblue")),
              column_names_side = "top",
              show_column_names = FALSE,
              column_names_rot = 90,
              row_names_side = "left",
              column_split = c(rep(" Endo. ", 4), rep(" Epithelial ", 12),
                               rep(" Lymphoid ", 9), rep(" Myeloid ", 12), rep("Mes.", 9)),
              row_split = rowsplit,
              row_names_gp = gpar(fontsize = 4),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_order = new_ct_order3,
              row_dend_width = unit(3, "mm"),
              row_dend_side = "right",
              border = "black",
              top_annotation = 
                columnAnnotation(`Cell Type` = colnames(exp_mat[, new_ct_order3]),
                                 col = list(`Cell Type` = unlist(new_color_list3)),
                                 show_legend = FALSE,
                                 annotation_name_gp = gpar(fontsize = 0),
                                 height = unit(1.2, "mm"),
                                 simple_anno_size = unit(1.2, "mm")),
              
              right_annotation =
                rowAnnotation(
                  `- 6.5%` = anno_points(gene_lineage_lumen_heatmap_df3$Mesenchymal, which = "row",
                                         ylim = c(0.99, 1.01), extend = 0, 
                                         size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                         pch = 15,
                                         axis = FALSE, border = FALSE, gp = gpar(col = "cornflowerblue")),
                  `- 29%` = anno_points(gene_lineage_lumen_heatmap_df3$Myeloid, which = "row",
                                        ylim = c(0.99, 1.01), extend = 0,
                                        size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                        pch = 15,
                                        axis = FALSE, border = FALSE, gp = gpar(col = "deeppink3")),
                  `- 6.5% ` = anno_points(gene_lineage_lumen_heatmap_df3$Lymphoid, which = "row",
                                          ylim = c(0.99, 1.01), extend = 0,
                                          size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                          pch = 15,
                                          axis = FALSE, border = FALSE, gp = gpar(col = "darkmagenta")),
                  `- 32.3%` = anno_points(gene_lineage_lumen_heatmap_df3$Epithelial, which = "row",
                                          ylim = c(0.99, 1.01), extend = 0,
                                          size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                          pch = 15,
                                          axis = FALSE, border = FALSE, gp = gpar(col = "chartreuse4")),
                  `- 25.8%` = anno_points(gene_lineage_lumen_heatmap_df3$Endothelial, which = "row",
                                          ylim = c(0.99, 1.01), extend = 0,
                                          size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                          pch = 15,
                                          axis = FALSE, border = FALSE, gp = gpar(col = "chocolate2")),
                  annotation_name_gp = gpar(fontsize = 5),
                  show_annotation_name = TRUE,
                  annotation_name_rot = -90,
                  annotation_name_side = "bottom",
                  gap = unit(1.4, "mm"),
                  "     " = rowsplit,
                  col = list("     " = c("4" = "#04846E", "3" = "#FFC107")),
                  show_legend = FALSE,
                  simple_anno_size = unit(2, "mm")
                ),
              gap = unit(1, "mm")
)

# Add annotations to heatmap
ht <- draw(hp, 
           heatmap_legend_side = "bottom",
           align_heatmap_legend = "heatmap_center")
ht

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig6c_heatmap_homeostasis_early.pdf", width = (79.5*0.0393701), height = (180.6*0.0393701))
ht
dev.off()
52*ht_opt(RESET = TRUE)


#### FIGURE 7D (SHOW MACROPHAGE LUMENS) ----
trans_tniche_lumen_metadata <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/trans_tniche_lumen_metadata_080724.csv") %>%
  left_join(xenium@meta.data %>% select(sample, sample_affect, percent_pathology) %>% unique()) %>%
  arrange(pseudotime_rank_trans_t4_tie_avg) %>%
  mutate(lumen_yes = ifelse(lumen_id %in% c("VUILD115_90", "VUILD78MF_27"), 1, 0))

top_annotation =
  HeatmapAnnotation(
    lumen_yes = trans_tniche_lumen_metadata$lumen_yes,
    `Alveolar Mϕ` = ct_mat[, "Alveolar Macrophages"],
    `SPP1+ Mϕ` = ct_mat[, "SPP1+ Macrophages"],
    Pseudotime = pseudotime_value,
    
    
    simple_anno_size = unit(2, "mm"),
    # height = unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = 5),
    gap = unit(0, "mm"),
    annotation_name_side = "left",
    
    col = list(
      lumen_yes = c(`1` = "red", `0` = "white"),
      `Alveolar Mϕ` = colorRamp2(c(0, max(ct_mat[, "Alveolar Macrophages"])), c("white", "darkorange3")),
      `SPP1+ Mϕ` = colorRamp2(c(0, max(ct_mat[, "SPP1+ Macrophages"])), c("white", "cornflowerblue")),
      Pseudotime = colorRamp2(c(min(pseudotime_value), max(pseudotime_value)), c("white", "black"))),
    
    show_legend = FALSE)
mac_lumens <- Heatmap(matrix(nc = nrow(trans_tniche_lumen_metadata), nr = 0),
                      top_annotation = top_annotation)
mac_lumens

pdf("/scratch/avannan/MANUSCRIPT_figures/Fig6d_mac_pseudotime.pdf", width = (75*0.0393701), height = (16*0.0393701))
mac_lumens
dev.off()


## SUPPLEMENTARY FIGURE 1 (Demographic Information) ----
# Demographic information
library(googlesheets4)
gs4_deauth()
ipf_sheet <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?usp=sharing")
ipf_sheet <- read_sheet(ipf_sheet, sheet = 1, skip = 1)

nSamples <- 45
nDonors <- nrow(ipf_sheet)
demo_plot <- ipf_sheet %>%
  select(-Age) %>%
  replace(is.na(.), "N/A") %>%
  mutate(`Clinical Diagnosis` = case_when(Clinical_Diagnosis == "IPF" ~ "IPF",
                                          Clinical_Diagnosis == "Control" ~ "Unaffected",
                                          Clinical_Diagnosis == "PH" ~ "PH",
                                          TRUE ~ "Other ILD"),
         Sample_Affect_Pairing2 = case_when(Sample_Affect_Pairing == "Duplicate_Unaffected" ~ "Rep. Unaffected",
                                            Sample_Affect_Pairing %in% c("More_Affected", "Less_Affected", "Unaffected") ~ "No Replicates",
                                            Sample_Affect_Pairing %in% c("Duplicate_Less_Affected", "Duplicate_More_Affected",
                                                                         "Paired_Less_and_More_Affected",
                                                                         "Paired_Less_and_More_Affected_with_Duplicate_More") ~ "Rep. Disease"),
         Punch_Size_Shape = case_when(Punch_Size == "3mm" & Punch_Shape == "Circular" ~ "3mm Round",
                                      Punch_Size == "5mm" & Punch_Shape == "Circular" ~ "5mm Round",
                                      Punch_Size == "3mm" & Punch_Shape == "Square" ~ "3mm Square",
                                      Punch_Size == "3mm" & Punch_Shape == "Circular & Square" ~ "3mm Round & Square")) %>%
  select(-Clinical_Diagnosis) %>% 
  pivot_longer(cols = c("Gender", "Ethnicity", "Tobacco", "Clinical Diagnosis", "Sample_Affect_Pairing2", "Punch_Size_Shape"),
               names_to = "Demographic_Feature", values_to = "Value") %>% 
  group_by(Demographic_Feature, Value) %>%
  mutate(percentage = nrow(Value)) %>%
  summarize(prop = length(Value)/nDonors) %>% 
  mutate(Demographic_Feature = ordered(Demographic_Feature, levels = c("Clinical Diagnosis", "Sample_Affect_Pairing2", "Ethnicity", "Gender", "Tobacco", "Punch_Size_Shape")),
         # Value = ordered(Value, levels = c("N/A", "PH", "Other ILD", "IPF", "Unaffected", "No Replicates", "Paired LF & MF", "Rep. Unaffected",
         #                                   "African American", "American Indian", "Hispanic", "European", "M", "F",  "N", "Y",
         #                                   "3mm Round", "5mm Round", "3mm Square"))) %>% 
         Value = ordered(Value, levels = c("N/A", "Other ILD", "IPF", "Unaffected", 
                                           "Rep. Disease",  "Rep. Unaffected",  "No Replicates", 
                                           "African American", "American Indian", "Hispanic", "European", "M", "F",  "N", "Y",
                                           "3mm Round", "5mm Round", "3mm Round & Square", "3mm Square"))) %>%
  filter(Demographic_Feature != "Gender") %>%
  ggplot(aes(x = Demographic_Feature, y = prop*100, fill = Value)) +
  geom_col(position = "stack", color = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = c("Diagnosis", "Replicates", "Ethnicity", "Ever Smoker", "Tissue Punch")) +
  scale_fill_manual(values = c("grey80", "#ac3e6f", "#d76160", "#88CCEE",
                               "#C25068","#88CCEE", "grey80", 
                               "#a265c2", "#009E73", "#FE6100", "#F0E442",
                               "#005AB5", "#DC3220",
                               "cornsilk1", "burlywood2", "bisque4", "chocolate3")) +
  labs(y = "% of Donors") +
  theme_bw(base_size = 6.5) + 
  theme(axis.text = element_text(color = "black"),
        axis.title.y = element_text(size = 7),
        axis.ticks = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        panel.grid.minor.y = element_blank()) +
  geom_text(aes(label = Value), position = position_stack(vjust = 0.5), size = 2)
demo_plot

table(xenium$final_CT, xenium$disease_status)["SPP1+ Macrophages", ]

png("/scratch/avannan/MANUSCRIPT_figures/SuppFig1_demographic_plot.png", width = 180*0.0393701, height = 120*0.0393701, units = "in", res = 450)
demo_plot
dev.off()


## SUPPLEMENTARY FIGURE 2 (Cell Types, Transcript Counts, etc.) ----
#### SUPPLEMENTARY FIGURE 2A (Lineage Proportions vs. Other Lung Studies) ----
# Add category to metadata
xenium$disease_status <- "Disease"
xenium$disease_status[xenium$sample_type == "Unaffected"] <- "Control"

# Load in and tidy data for comparison
lung_ct_compare <- readxl::read_excel("/scratch/avannan/MANUSCRIPT_figures/CT Recovery Spatial vs SC_06-24-24.xlsx",
                                      sheet = 1, skip = 3) %>%
  dplyr::rename(Overall_Study = `...1`) %>%
  # mutate(across(3:ncol(.), function(x) as.numeric(x))) %>%
  .[-c(16:21), -c(15:17)] %>%
  pivot_longer(3:ncol(.), names_to = "Disease_Status", values_to = "Percentage") %>%
  mutate(Disease_Status = ordered(case_when(grepl("Disease", Disease_Status) ~ "Disease",
                                            grepl("Control", Disease_Status) ~ "Control",
                                            grepl("Total", Disease_Status) ~ "Total"),
                                  levels = c("Control", "Disease"))) %>%
  filter(!is.na(Disease_Status), !is.na(Overall_Study)) %>%
  cbind(Lineage = c(rep("Epithelial", 2), rep("Endothelial", 2), rep("Immune", 2), rep("Mesenchymal", 2))) %>%
  mutate(Lineage = ordered(Lineage, levels = c("Endothelial", "Epithelial", "Immune", "Mesenchymal")))

# Get dataframe of sample types for each sample
sample_type_df <- xenium@meta.data %>%
  select(sample, disease_status) %>%
  unique() %>%
  dplyr::rename(Sample = "sample", Disease_Status = "disease_status")
rownames(sample_type_df) <- NULL

# Dataframe with cell lineage proportions for each sample for spatial data
spatial_df <- (table(xenium$sample, xenium$final_lineage)/
                 rowSums(table(xenium$sample, xenium$final_lineage))) %>%
  as.data.frame() %>%
  dplyr::rename(Sample = "Var1", Lineage = "Var2", Percentage = "Freq") %>%
  mutate(Percentage = Percentage*100, Overall_Study = "Spatial",
         Lineage = ordered(Lineage, levels = c("Endothelial", "Epithelial", "Immune", "Mesenchymal"))) %>%
  full_join(sample_type_df, relationship = "many-to-many")

# Compare to eQTL study as well
# eqtl.ref <- readRDS("/scratch/avannan/eQTL_rds_files/ILD_all_celltypes_Seurat.rds")
# eqtl <- eqtl.ref@meta.data
# write.csv(eqtl, "/scratch/avannan/eQTL_rds_files/SE227136_ILD_all_celltypes_seurat_meta.csv")
eqtl <- read.csv("/scratch/avannan/eQTL_rds_files/SE227136_ILD_all_celltypes_seurat_meta.csv") %>%
  select(Sample_Name, Status, lineage) %>%
  mutate(Lineage = case_when(lineage %in% c("MyoFB", "MyoFB - Activated", "PLIN2+ FB") ~ "Mesenchymal",
                             lineage == "PNEC" ~ "Epithelial",
                             lineage == "Inflamed" ~ "Endothelial",
                             TRUE ~ lineage),
         Dataset = Sample_Name) %>%
  filter(grepl("T", .$Dataset)) %>% # Only retain TGen samples
  select(-lineage, -Sample_Name) %>%
  dplyr::rename(Disease_Status = "Status")

# Create final eQTL df
eqtl_df <- (table(eqtl$Dataset, eqtl$Lineage)/rowSums(table(eqtl$Dataset, eqtl$Lineage))) %>%
  as.data.frame() %>%
  dplyr::rename(Dataset = "Var1", Lineage = "Var2", Percentage = "Freq") %>%
  mutate(Percentage = Percentage*100, Overall_Study = "Natri et al. (2023)",
         Lineage = ordered(Lineage, levels = c("Endothelial", "Epithelial", "Immune", "Mesenchymal"))) %>%
  full_join(eqtl %>% select(Dataset, Disease_Status) %>% unique())

# Create dataframe that compares spatial data to other datasets
lung_ct_compare2 <- spatial_df %>%
  dplyr::rename(Dataset = "Sample") %>%
  full_join(lung_ct_compare) %>%
  full_join(eqtl_df) %>%
  mutate(Proportion = Percentage/100,
         Overall_Study = ordered(Overall_Study,
                                 levels = c("Spatial", "Natri et al. (2023)", "Human Lung Cell Atlas")))

# Create dataframe for plotting
overall_prop_table <- table(xenium$disease_status, xenium$final_lineage)/rowSums(table(xenium$disease_status, xenium$final_lineage))
lung_ct_compare_4plot <- lung_ct_compare2 %>%
  mutate(yintercept = case_when(Lineage == "Epithelial" & Disease_Status == "Control" ~ overall_prop_table["Control", "Epithelial"],
                                Lineage == "Endothelial" & Disease_Status == "Control" ~ overall_prop_table["Control", "Endothelial"],
                                Lineage == "Immune" & Disease_Status == "Control" ~ overall_prop_table["Control", "Immune"],
                                Lineage == "Mesenchymal" & Disease_Status == "Control" ~overall_prop_table["Control", "Mesenchymal"],
                                Lineage == "Epithelial" & Disease_Status == "Disease" ~ overall_prop_table["Disease", "Epithelial"],
                                Lineage == "Endothelial" & Disease_Status == "Disease" ~ overall_prop_table["Disease", "Endothelial"],
                                Lineage == "Immune" & Disease_Status == "Disease" ~ overall_prop_table["Disease", "Immune"],
                                Lineage == "Mesenchymal" & Disease_Status == "Disease" ~ overall_prop_table["Disease", "Mesenchymal"]))

# Plot
lung_compare_plot <- lung_ct_compare_4plot %>% 
  ggplot(aes(x = Overall_Study, y = Proportion, fill = Overall_Study)) +
  geom_blank(data = lung_ct_compare_4plot) +
  geom_hline(aes(yintercept = yintercept, color = Lineage), lty = "longdash",
             show.legend = FALSE) +
  geom_boxplot(width = 0.65, show.legend = FALSE, outlier.shape = 21, 
               size = 0.25, outlier.size = 0.75, outlier.stroke = 0.25,
               color = "black", na.rm = TRUE) +
  theme_bw(base_size = 6) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.background = element_rect(linewidth = 50),
        axis.text = element_text(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        strip.text = element_text(color = "white", size = 6),
        axis.title = element_text(color = "black"),
        plot.margin = margin(5.5, 3, 3, 3, "pt"),
        axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.035, 0)) +
  scale_fill_manual(values = c("white", "grey80", "grey35", "black")) +
  scale_color_manual(values = c("chocolate2", "chartreuse4", "deeppink3", "cornflowerblue")) +
  scale_x_discrete(labels = c("Spatial", "Natri et al.\n(2023)", "HLCA")) +
  # facet_grid(Lineage~Disease_Status, scales = "free", space = "free")
  facet_grid(Disease_Status~Lineage, scales = "free", space = "free")
lung_compare_plot

# Change facet colors
lung_compare_plot_table <- ggplot_gtable(ggplot_build(lung_compare_plot))
striprt <- which(grepl("strip-r", lung_compare_plot_table$layout$name) | 
                   grepl("strip-t", lung_compare_plot_table$layout$name))
fills <- c("chocolate2", "chartreuse4", "maroon3", "cornflowerblue", NA, NA) # Lineages & Sample Types
colors <- c(rep("black", 4), rep(NA, 2))
font_colors <- c(rep("white", 4), rep("black", 2))
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", lung_compare_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  lung_compare_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  lung_compare_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- colors[k]
  lung_compare_plot_table$grobs[[i]]$grobs[[1]]$children[[j+1]]$children[[1]]$gp$col <- font_colors[k]
  k <- k+1
}
grid::grid.draw(lung_compare_plot_table)

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig2B_ct_prop_and_compare_062824.pdf", width = (123*0.0393701), height = (60*0.0393701))
plot(lung_compare_plot_table)
dev.off()


#### SUPPLEMENTARY FIGURE 2B: (Transcripts by Sample and Transcripts by Cell by Sample) ----
# Load in transcript counts and format dataframe
transcript_counts <- readxl::read_excel("/scratch/avannan/MANUSCRIPT_figures/sample_node_counts_TIDY_062624.xlsx") %>%
  separate(col = Filename, into = c("sample", "other"), sep = "_", remove = FALSE) %>%
  separate(col = other, into = c("tma", "other", sep = ".", remove = FALSE)) %>%
  select(Filename, sample, tma, "AfterQC/InputToGraph", "Clustered/Special_AfterRemoveExtraTranx") %>%
  mutate(tma = ifelse(tma == "detected", "TMA5", tma),
         sample = ifelse(sample == "TILD117MF" & tma == "TMA5", "TILD117MFB", sample)) %>%
  filter(Filename != "VUILD105LF_TMA2.SOMEREMOVED_node_meta.csv") %>%
  mutate(sample_size = ordered(case_when(tma == "TMA5" ~ "3mm, Square",
                                         tma == "TMA3" ~ "5mm, Round",
                                         TRUE ~ "3mm, Round"),
                               levels = c("3mm, Round", "3mm, Square", "5mm, Round")))
colnames(transcript_counts)

transcript_counts_plot <- transcript_counts %>%
  left_join(xenium@meta.data %>% select(sample, new_sample_name, tma, sample_type) %>% unique()) %>%
  mutate(sample_type = ordered(sample_type, levels = c("Unaffected", "LF", "MF", "INT"))) %>%
  full_join(xenium@meta.data %>% dplyr::select(sample, percent_pathology, sample_affect) %>% unique) %>%
  mutate(sample_affect = ordered(sample_affect, levels = c("Unaffected", "Less Affected", "More Affected"))) %>%
  ggplot(aes(x = reorder(new_sample_name, percent_pathology), y = `Clustered/Special_AfterRemoveExtraTranx`,
             fill = sample_affect)) +
  geom_col() +
  geom_col(aes(lty = sample_size), color = "black", show.legend = FALSE) +
  geom_line(aes(lty = sample_size)) +
  scale_fill_manual(values = sample_affect_color_list) +
  scale_linetype_manual(values = c("solid", "dashed", NA)) +
  theme_classic(base_size = 6) +
  scale_shape_manual(values = c(21, 22)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Transcript Count") +
  theme(legend.text = element_text(color = "black"),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        # legend.position = c(0.18, 0.8),
        legend.position = "right",
        plot.margin = margin(10,5.5,5.5,5.5, "pt")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, keywidth = 1, keyheight = 1),
         lty = guide_legend(title.hjust = 0.5, keywidth = 2, keyheight = 0.2))
transcript_counts_plot

# Get data per cell and plot
transcript_counts_plot2 <- xenium@meta.data %>%
  mutate(sample_size = ordered(case_when(tma == "TMA5" ~ "3mm, Square",
                                         tma == "TMA3" ~ "5mm, Round",
                                         TRUE ~ "3mm, Round"),
                               levels = c("3mm, Round", "3mm, Square", "5mm, Round"))) %>%
  ggplot(aes(x = reorder(new_sample_name, percent_pathology), y = nCount_RNA,
             fill = sample_affect, lty = sample_size)) +
  geom_violin() +
  labs(y = "Transcripts per Cell") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = sample_affect_color_list) +
  scale_linetype_manual(values = c("solid", "dashed", NA)) +
  theme_classic(base_size = 6) +
  theme(legend.text = element_text(color = "black"),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        strip.text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.margin = margin(5.5,5.5,5.5,5.5, "pt")) +
  NoLegend()
transcript_counts_plot2

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig2_transcript_counts_plot2.pdf", width = 207*0.0393701, height = 90*0.0393701)
(transcript_counts_plot + transcript_counts_plot2) + patchwork::plot_layout(ncol = 1)
dev.off()


#### SUPPLEMENTARY FIGURE 2C: (Transcripts by Cell by Sample Type) ----
transcript_counts_plot3 <- xenium@meta.data %>%
  mutate(sample_size = ordered(case_when(tma == "TMA5" ~ "3mm, Square",
                                         tma == "TMA3" ~ "5mm, Round",
                                         TRUE ~ "3mm, Round"),
                               levels = c("3mm, Round", "3mm, Square", "5mm, Round")),
         sample_affect = ordered(sample_affect, levels = c("Unaffected", "Less Affected", "More Affected"))) %>%
  ggplot(aes(x = sample_affect, y = nCount_RNA, fill = sample_affect)) +
  geom_violin() +
  geom_boxplot(fill = "white", width = 0.3, outlier.alpha = 0.5, outlier.size = 1) +
  labs(y = "Transcripts per Cell") +
  scale_fill_manual(values = sample_affect_color_list) +
  theme_classic(base_size = 6) +
  scale_y_continuous(limits = c(0, 810), expand = c(0, 0)) +
  theme(legend.text = element_text(color = "black"),
        axis.text.x = element_text(size = 5),
        strip.text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.position = c(0.18, 0.8),
        plot.margin = margin(10,5.5,5.5,5.5, "pt")) +
  NoLegend()
transcript_counts_plot3

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig2_transcript_counts_plot3.pdf", width = 52*0.0393701, height = 46*0.0393701)
transcript_counts_plot3
dev.off()


## SUPPLEMENTARY FIGURE 3 (Cell Type Expression Heatmap) ----
# Set up genes for heatmaps
epi_heat <-  c("AGER", "AGR3", "BPIFA1", "C20orf85", "CALCA", "CCN2", "CCNA1", "CD44",    
               "CDH26", "CDK1", "CDKN2A", "CEACAM5", "CEACAM6", "CFTR", "CHGB", "COL1A1",  
               "COL4A3", "DMBT1", "DUOX1", "EPCAM", "ERN2", "FN1", "FOXI1", "FOXJ1",   
               "GKN2", "ICAM1", "IFIT2", "IFIT3", "ITGAV", "ITGB1", "ITGB6", "KRT14",   
               "KRT15", "KRT17", "KRT18", "KRT5", "KRT6A", "KRT8", "LAMP3", "LTF", "MKI67",   
               "MMP10", "MMP7", "MUC5AC", "MUC5B", "NAPSA", "NKX2-1", "OAS2", "PDIA4",   
               "PGC", "RTKN2", "S100A2", "SAA2", "SCG2", "SCGB1A1", "SCGB3A2", "SEC11C",  
               "SFTA2", "SFTPC", "SFTPD", "SOX2", "SOX4", "SOX9", "SPINK1", "TGFB2",   
               "TOP2A", "TP63", "TP73", "VEGFA", "WFDC2")
imm_heat <-  c("AIF1", "AXL", "BANK1", "BCL2L11", "CCL18", "CCL2", "CCL22", "CCL5",     
               "CCR7", "CD14", "CD19", "CD1A", "CD1C", "CD2", "CD247", "CD27",     
               "CD28", "CD3D", "CD3E", "CD4", "CD68", "CD69", "CD79A", "CD79B",    
               "CD86", "CD8A", "CD8B", "CDK1", "CPA3", "CST3", "CTLA4", "CXCR4", "FABP4",    
               "FASLG", "FCER1A", "FCER1G", "FCGBP", "FCGR3A", "FCN1", "FGFBP2", "FOXP3",    
               "GNLY", "GPR183", "GZMA", "GZMB", "GZMK", "HAVCR2", "HLA-DQA1", "HLA-DRA",  
               "IDH1", "IFIT1", "IFIT3", "IFNG", "IL2RA", "IL7R", "IRF7", "ITGAE",    
               "ITGAM", "ITGAX", "ITM2C", "JCHAIN", "KIT", "KLRB1", "KLRG1", "LAG3",     
               "LCK", "LEF1", "LILRA4", "LTB", "LYZ", "MARCO", "MCEMP1", "MKI67",    
               "MS4A1", "MS4A7", "NFKB1", "NKG7", "OAS2", "PIM2", "PKM", "PLIN2",    
               "PPARG", "PTPRC", "S100A12", "S100A8", "S100A9", "SLC1A3", "SLC25A37", 
               "SOD2", "SPP1", "TNFRSF13C", "TOP2A", "TPSAB1", "TRAC")
endo_heat <-   c("ACKR1", "APLN", "APLNR", "BMPR2", "CA4", "CCL21", "CD34",
                 "CLDN5", "EPAS1", "FCN3", "GNG11", "HEY1", "KDR", "PECAM1", 
                 "PLVAP", "RAMP2", "TBXA2R")
mes_heat <-  c("ACTA2", "CDK1", "COL1A1", "COL1A2", "COL3A1", "CSPG4", "CTHRC1", 
               "DCN", "ELN", "EPAS1", "FABP4", "FAP", "FGF10", "FGF2", "FGF7", 
               "FN1", "HAS1", "HAS2", "HIF1A", "KRT18", "KRT8", "LGR5", "LGR6",
               "LUM", "MEG3", "MFAP5", "MKI67", "MSLN", "PDGFRA", "PDGFRB", "PI16",
               "PLIN2", "POSTN", "PPARG", "PTGDS", "SCGB3A2", "SFRP2", "SFRP4", "SFTPC", 
               "SPARCL1", "SPRY2", "TGFB3", "TOP2A", "VIM", "WNT5A", "YAP1")

all_genes_for_heat <- c(endo_heat, mes_heat, epi_heat, imm_heat)
genes_4heat <- unique(all_genes_for_heat)

# Simplify list of genes
genes_4heat <- genes_4heat[-which(genes_4heat %in% 
                                    c("CCN2", "FGF2", "LGR6", "SPRY2", "FOXI1", "AGR3", 
                                      "KRT18", "IRF7", "IFIT1", "FGFBP2", "CDKN2A", "CDK1",
                                      "BCL2L11", "CD19", "CD79A", "CD79B", "LTB", "LAG3",
                                      "CD27", "ITGAE", "IFNG", "FASLG", "GZMK", "CCL5",
                                      "CD247", "LCK", "KLRG1", "CD69", "CXCR4", "SLC1A3",
                                      "ITGAM", "HAVCR2", "AIF1", "MMP10", "S100A2",
                                      "KRT6A", "KRT15", "KRT14", "LTF", "BPIFA1",
                                      "SAA2", "CCNA1", "ERN2", "LEF1", "SOX2", "KLRB1",
                                      "IL2RA", "CD28", "GKN2", "CDH26", "SOX9", "CEACAM5",
                                      "DMBT1", "SFTPD", "ICAM1", "SEC11C", "PDIA4",
                                      "CD8B", "CD86", "VIM", "EPAS1", "APLNR",
                                      "RAMP2", "TGFB3", "FGF10", "SPARCL1", "ELN",
                                      "MEG3", "SFRP4", "TBXA2R", "PDGFRB", "CD44",
                                      "PKM", "IDH1", "NFKB1", "MUC5AC", "GNG11",
                                      "KDR", "S100A12", "COL1A2", "COL3A1", "HAS2",
                                      "FN1", "TGFB2", "NAPSA", "SFTA2", "TP73",
                                      "WFDC2", "SCG2", "LGR5", "IFIT2", "TRAC", 
                                      "ITGAV", "HLA-DRA", "ITGB1"))]
genes_4heat <- unique(c(genes_4heat, "MUC5AC", "PDGFRB", "KRT6A", "WT1", "COL3A1", "FGF2", "FN1", "TREM2"))

# Set up objects for heatmaps!
# Load in main object and separate by lineage
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_070924.rds")
# all_epi <- subset(xenium, subset = final_lineage == "Epithelial")
# all_imm <- subset(xenium, subset = final_lineage == "Immune")
# all_endo_mes <- subset(xenium, subset = final_lineage %in% c("Endothelial", "Mesenchymal"))

# Scale
xenium <- ScaleData(xenium)
# all_epi <- ScaleData(all_epi)
# all_imm <- ScaleData(all_imm)
# all_endo_mes <- ScaleData(all_endo_mes)

# Set cell type order
new_ct_order <- c(
  # Endothelial 
  "Arteriole", "Capillary", "Venous", "Lymphatic",
  # Epithelial - Alveolar
  "AT1", "Transitional AT2", "AT2", "Proliferating AT2", "KRT5-/KRT17+",              
  # Epithelial - Airway
  "Basal", "Multiciliated", "Goblet", "RASC", "Secretory", "PNEC", "Proliferating Airway",
  # Immune - Lymphoid
  "B cells", "Proliferating B cells", "NK/NKT", "Proliferating NK/NKT", "Tregs", 
  "CD4+ T-cells", "CD8+ T-cells", "Proliferating T-cells", "Plasma", "pDCs",
  # Immune - Myeloid
  "Basophils", "cDCs", "Migratory DCs", "Langerhans cells", "Mast", 
  "Alveolar Macrophages", "Interstitial Macrophages", "SPP1+ Macrophages",
  "Macrophages - IFN-activated", "Monocytes/MDMs", "Neutrophils", "Proliferating Myeloid",
  # Mesenchymal
  "Activated Fibrotic FBs", "Adventitial FBs", "Alveolar FBs", "Inflammatory FBs", 
  "Myofibroblasts", "Subpleural FBs", "Proliferating FBs", "SMCs/Pericytes", "Mesothelial")

# Create heatmap!
# Create dotplot base
p <- DotPlot(xenium, features = genes_4heat, group.by = "final_CT", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
p

# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() %>%
  dplyr::select(features.plot, new_ct_order)

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()

# The percentage of cells expressing a feature
percent_mat <- df %>% 
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%  
  tidyr::pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() %>%
  dplyr::select(features.plot, new_ct_order)

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

dim(exp_mat)
dim(percent_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# Replace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Any value that is greater than 3 will be mapped to bright red
col_fun = circlize::colorRamp2(c(-1, 0, 2), colorspace::diverge_hsv(3))

# Creating a layer to add to plot
layer_fun = function(j, i, x, y, w, h, fill){
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex(exp_mat, i, j))), # t(exp_mat)
              size = pindex(percent_mat, i, j)/100 * unit(1.4, "mm"), # t(percent_mat)
              pch = 19)
}

# List of legends
lgd_list = list(
  Legend(labels = c(0,0.25,0.5,0.75,1), title = "Proportion", type = "points", 
         pch = 19, size = c(0,0.25,0.5,0.75,1) * unit(1.4, "mm"),
         legend_gp = gpar(col = "black"), direction = "vertical", ncol = 1,
         title_position = "topcenter", background = NA, 
         grid_height = unit(2, "mm"),
         grid_width = unit(2, "mm"),
         title_gp = gpar(fontsize = 5), labels_gp = gpar(fontsize = 5),
         legend_height = unit(5, "mm")))

# Create heatmap
rename_cts <- colnames(exp_mat)
rename_cts <- gsub("Proliferating", "Prolif.", rename_cts)
rename_cts <- gsub("Macrophages", "MΦ", rename_cts)
rename_cts <- gsub("Activated", "Activ.", rename_cts)
rename_cts <- gsub("Fibrotic", "Fibr.", rename_cts)
rename_cts <- gsub("activated", "activ.", rename_cts)
rename_cts <- gsub("Inflammatory", "Inflamm.", rename_cts)
rename_cts <- gsub("Monocytes", "Mono.", rename_cts)
rename_cts <- gsub("Langerhans cells", "Langerhans", rename_cts)
colnames(exp_mat) <- rename_cts
ct_col_order <- c(rep("   Endothelial   ", 4),
                  rep("  Epithelial  ", 12), 
                  rep("  Lymphoid  ", 10),
                  rep("  Myeloid  ", 12),
                  rep(" Mesenchymal ", 9))
ht_opt$TITLE_PADDING = unit(1, "mm")
ht_opt$ANNOTATION_LEGEND_PADDING = unit(0, "mm")
ht_opt$HEATMAP_LEGEND_PADDING = unit(0, "mm")
ht_opt$DIMNAME_PADDING = unit(0.5, "mm")
hp <- Heatmap(exp_mat, # t(exp_mat)
              name = "Scaled Expression",
              height = unit(170, "mm"),
              width = unit(150, "mm"),
              heatmap_legend_param = list(title = "Scaled\n",
                                          labels = c("-1", "0", "1", "2+"),
                                          legend_direction = "vertical",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          background = NA,
                                          grid_width = unit(2, "mm"),
                                          grid_height = unit(6, "mm"),
                                          legend_height = unit(10, "mm"),
                                          legend_width = unit(2, "mm")),
              row_title = " ",
              row_title_gp = gpar(fontsize = 5),
              column_title_gp = gpar(fontsize = 5),
              row_title_side = "right",
              column_title_side = "top",
              column_names_rot = 35,
              column_split = ct_col_order,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 4),
              row_names_side = "left",
              column_names_gp = gpar(fontsize = 5),
              cluster_rows = TRUE, cluster_columns = FALSE,
              cluster_column_slices = FALSE,
              row_dend_width = unit(3, "mm"),
              border = "black",
              gap = unit(2, "mm"))

ht <- draw(hp, 
           heatmap_legend_list = lgd_list,
           heatmap_legend_side = "right",
           annotation_legend_side = "right",
           align_heatmap_legend = "global_center",
           padding = unit(c(0, 0, 0, 0), "mm"))
ht

png("/scratch/avannan/MANUSCRIPT_figures/SuppFig3_heatmap_all_celltypes_071924.png",
    width = 180*0.0393701, height = 183*0.0393701, units = "in", res = 450)
ht
dev.off()
ht_opt(RESET = TRUE)

# # Get rowname order
# rownames(ht@ht_list[["Scaled Expression"]]@matrix[unlist(ht@ht_list[["Scaled Expression"]]@row_order_list), ]) %>%
#   as.data.frame()


## SUPPLEMENTARY FIGURE 9 (Cell Proximity Heatmap Overall) ----
library(data.table)
library(dplyr)
library(EnrichedHeatmap)
library(circlize)
library(RColorBrewer)

calc_heatmap_w <- function(df){width <- ncol(df) * 2.5; return(width)}
calc_heatmap_h <- function(df){height <- nrow(df) * 2.5; return(height)}

'%!in%' <- function(x,y)!('%in%'(x,y))
filter <- dplyr::filter
select <- dplyr::select

# load data and define code paramters
workDir <- '/scratch/avannan/late_IPF_spatial/xenium/REVISION_prox_analysis_av/Plots/'

celltype_logit_res <- read.csv(file.path(workDir, 'Cell_ProxinmityEnrichment_wOR_forAV3.csv'))
colnames(celltype_logit_res) <- c('target', 'prob', 'coefficient', 'std.err', 'z.value', 'p.value', 'OR', 'ci2.5', 'ci97.5', 'source')

celltype_logit_res <- celltype_logit_res %>%
  group_by(source) %>%
  mutate(fdr = p.adjust(p.value, method = 'fdr'))
celltype_logit_res$logOR <- log(celltype_logit_res$OR)

prob <- reshape2::dcast(source ~ target, value.var = 'prob', data = celltype_logit_res)
pval <- reshape2::dcast(source ~ target, value.var = 'fdr', data = celltype_logit_res)
or <- reshape2::dcast(source ~ target, value.var = 'logOR', data = celltype_logit_res)

rownames(prob) <- prob$source; prob$source <- NULL
rownames(pval) <- pval$source; pval$source <- NULL
rownames(or) <- or$source; or$source <- NULL

# convert NA to 0s
ct_prob <- prob
# ct_prob[is.na(ct_prob)] = 0

ct_pval<- pval
# ct_pval[is.na(ct_pval)] = 1

ct_or <- or
# ct_or[is.na(ct_or)] = 0

# check all tables are in the same order
mean(rownames(ct_prob) == rownames(ct_pval))
mean(colnames(ct_prob) == colnames(ct_pval))
mean(rownames(ct_or) == rownames(ct_pval))
mean(colnames(ct_or) == colnames(ct_pval))

# convert none-significant results to 0s
ct_prob[ct_pval > 0.05] = 0
ct_or[ct_pval > 0.05] = 0

# plot odds ratio
min01 <- quantile(as.matrix(ct_or), 0.01, na.rm = TRUE)
max99 <- quantile(as.matrix(ct_or), 0.99, na.rm = TRUE)
col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))
ct_or_ht <- Heatmap(as.matrix(ct_or), name = 'log(OR)',
                    heatmap_legend_param = list(legend_direction = "vertical",
                                                title_position = "topcenter",
                                                title_gp = gpar(fontsize = 6),
                                                labels_gp = gpar(fontsize = 6),
                                                background = NA,
                                                grid_width = unit(2, "mm"),
                                                grid_height = unit(6, "mm"),
                                                legend_height = unit(10, "mm"),
                                                legend_width = unit(2, "mm")),
                    clustering_distance_columns = , clustering_distance_rows = ,
                    clustering_method_columns = , clustering_method_rows = ,
                    row_title = 'Starting Cell', row_title_side = 'left', row_title_gp = gpar(fontsize = 7),
                    row_names_gp = gpar(fontsize = 6),
                    row_dend_side = "right", row_names_side = "left",
                    column_title = 'Nearest Neighbor Cell', column_names_side = 'top', column_dend_side = 'bottom',  column_title_gp = gpar(fontsize = 7),
                    column_names_gp = gpar(fontsize = 6),
                    rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
                    width = unit(128.1, "mm"), height = unit(128.1, "mm"),
                    na_col = 'black', col = col_fun, column_names_rot = 45, row_dend_width = unit(0.5, 'cm'), column_dend_height = unit(0.5, 'cm'))
draw(ct_or_ht)

png("/scratch/avannan/MANUSCRIPT_figures/SuppFig_all_cell_prox_results.png",
    width = (180*0.0393701), height = (165*0.0393701), unit = "in", res = 450)
draw(ct_or_ht)
dev.off()



## SUPPLEMENTARY FIGURE 10 (Absolute Proximity of KRT5-/KRT17+ to Other Epithelial Cells) ----
onion = read.csv("/scratch/avannan/MANUSCRIPT_figures/cell_proximity_by_sample_KRTs.csv")
onion = onion[grep("HD", onion$sample, invert = T),]
onion = onion[,grep("distal", colnames(onion), invert = T)]
onion = onion[,grep("ct", colnames(onion), invert = T)]
onion = onion[onion$prox.cell.type != "KRT5-/KRT17+",]
pathology = xenium@meta.data %>% select(sample, sample_affect) %>% unique()
onion = base::merge(onion, pathology, all.x = T)

epi_ct_order <- c(
  # Epithelial - Alveolar
  "AT1", "Transitional AT2", "AT2", "Proliferating AT2",            
  # Epithelial - Airway
  "Basal", "Multiciliated", "Goblet", "RASC", "Secretory", "PNEC", "Proliferating Airway")
epi_ct_colors <- color_list[epi_ct_order]
epi_ct_order <- gsub("Proliferating", "Prolif.", epi_ct_order)

krt5_proximity_epi_plot <- onion %>%
  mutate(sample_affect = "All Disease") %>%
  rbind(onion) %>%
  mutate(prox.cell.type = gsub("Proliferating", "Prolif.", prox.cell.type),
         prox.cell.type = ordered(prox.cell.type, epi_ct_order),
         sample_affect = ordered(sample_affect, levels = c("Less Affected", "More Affected", "All Disease"))) %>%
  ggplot(aes(x = prox.cell.type, y=N.proximal, fill = sample_affect, color = sample_affect)) + 
  geom_boxplot(outlier.shape = 21, width = 1, position = position_dodge(width = 1.2)) +
  theme_classic(base_size = 6.5) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 130)) +
  scale_fill_manual(values = c("pink", "#ac3e6f", "#5A012D")) +
  scale_color_manual(values = c("black", "black", "grey40")) +
  theme_black +
  labs(x = "Nearest Neighbor Cell Type", y = "Number of Proximal Cells") +
  theme(legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(0.92, 0.6)) +
  facet_grid(~prox.cell.type, scales = "free")
krt5_proximity_epi_plot

# Change facet colors
krt5_proximity_epi_plot_table <- ggplot_gtable(ggplot_build(krt5_proximity_epi_plot))
striprt <- which(grepl("strip-r", krt5_proximity_epi_plot_table$layout$name) | 
                   grepl("strip-t", krt5_proximity_epi_plot_table$layout$name))
fills <- unlist(epi_ct_colors)
font_colors <- c(rep("black", 2), rep("white", 2), rep("black", 6), "white")
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", krt5_proximity_epi_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  krt5_proximity_epi_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  krt5_proximity_epi_plot_table$grobs[[i]]$grobs[[1]]$children[[j+1]]$children[[1]]$gp$col <- font_colors[k]
  k <- k+1
}
grid::grid.draw(krt5_proximity_epi_plot_table)

png("/scratch/avannan/MANUSCRIPT_figures/SuppFig_absolute_prox_KRT5.png", width = (180*0.0393701), height = (48*0.0393701),
    units = "in", res = 450)
grid::grid.draw(krt5_proximity_epi_plot_table)
dev.off()

onion %>%
  group_by(prox.cell.type) %>%
  mutate(sum_all = sum(N.proximal)) %>%
  select(prox.cell.type, sum_all) %>%
  unique() %>%
  arrange(sum_all)


## SUPPLEMENTARY FIGURE 11 (Cell Type Proportions by Sample & Sample Type) ----
# Cell type proportions
a <- xenium@meta.data %>%
  group_by(final_sublineage) %>%
  mutate(count = length(final_sublineage),
         sample_affect = ordered(sample_affect, levels = c("Unaffected", "Less Affected", "More Affected"))) %>%
  ggplot(aes(x = sample_affect, fill = final_sublineage)) +
  geom_bar(position = position_fill(), width = 1, color = "black", linewidth = 0.25) +
  labs(y = "Cell Type Proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("chartreuse4", "chartreuse3", 
                               "darkmagenta", "deeppink3", 
                               "chocolate2", "cornflowerblue")) +
  theme_classic(base_size = 6.5) +
  theme(legend.text = element_text(color = "black"),
        legend.title = element_blank(),
        strip.text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), 
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = NA, color = "black"),
        legend.key.size = unit(5, "mm"),
        legend.position = "bottom",
        plot.margin = margin(10,5.5,5.5,5.5, "pt"))
a

# Cell type proportions by sample type
b <- xenium@meta.data %>%
  group_by(final_sublineage) %>%
  mutate(count = length(final_sublineage)) %>%
  ggplot(aes(x = reorder(new_sample_name, percent_pathology), fill = final_sublineage)) +
  geom_bar(position = position_fill(), width = 1, color = "black", size = 0.25) +
  labs(y = "Cell Type Proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("chartreuse4", "chartreuse3", 
                               "darkmagenta", "deeppink3", 
                               "chocolate2", "cornflowerblue")) +
  theme_grey(base_size = 6.5) +
  theme(legend.text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.ticks = element_line(color = "black"), 
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(10,5.5,5.5,5.5, "pt")) +
  NoLegend()
b

ggsave("/scratch/avannan/MANUSCRIPT_figures/SuppFig5A_ct_sample_type.pdf", width = 50*0.0393701, height = 85*0.0393701, device = "pdf", units = "in")
a
dev.off()

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig5B_ct_sample.pdf", width = 130*0.0393701, height = 80*0.0393701)
b
dev.off()


## SUPPLEMENTARY FIGURE 12 (Cell Type by Sample Type Boxplots) ----
# Load dataframes
library(googlesheets4)
gs4_deauth()
propeller_res <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094")
propeller_res_sheet <- read_sheet(propeller_res, sheet = 3, skip = 1)
propeller_res_df <- propeller_res_sheet %>%
  as.data.frame() %>%
  select(-1)
propeller_res2 <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094")
propeller_res_sheet2 <- read_sheet(propeller_res2, sheet = 4, skip = 1)
propeller_res_df2 <- propeller_res_sheet2 %>%
  as.data.frame()

# Prepare dataframes
celltypeprop_df <- xenium@meta.data %>%
  group_by(final_CT, sample, sample_affect) %>%
  summarise(celltype_count = n()) %>%
  group_by(sample, sample_affect) %>%
  mutate(sample_total_cells = sum(celltype_count)) %>%
  mutate(celltype_prop_in_sample = celltype_count/sample_total_cells) %>%
  ungroup() %>%
  left_join(propeller_res_df, by = c("final_CT")) %>%
  mutate(final_CT2 = ifelse(!is.na(`Significance level`), paste0(final_CT, " ", `Significance level`), final_CT),
         dummy_x = case_when(sample_affect == "Unaffected" ~ 0,
                             sample_affect == "Less Affected" ~ 1,
                             sample_affect == "More Affected" ~ 2))

# # Set significance labels and locations
# propeller_res_df3 <- propeller_res_df2 %>%
#   left_join(celltypeprop_df, by = c("celltype" = "final_CT")) %>%
#   unique() %>%
#   group_by(celltype) %>%
#   summarise(max_pro = max(celltype_prop_in_sample)+1.2*sqrt(var(celltype_prop_in_sample)),
#             max_pro2 = max(celltype_prop_in_sample)+0.9*sqrt(var(celltype_prop_in_sample)),
#             max_pro3 = max(celltype_prop_in_sample)+0.7*sqrt(var(celltype_prop_in_sample))) %>%
#   ungroup() %>%
#   left_join(propeller_res_df2) %>%
#   mutate(`Significance Level` = case_when(is.na(`Significance Level`) ~ "", 
#                                           `Significance Level` == "**" ~ " ** ",
#                                           `Significance Level` == "*" ~ " * ")) %>%
#   dplyr::rename(final_CT = celltype) %>%
#   select(final_CT, max_pro, max_pro2, max_pro3, Comparison, `Significance Level`) %>%
#   unique()
# propeller_res_df3_a <- propeller_res_df3 %>%
#   filter(Comparison == "Less Affected vs. Unaffected")
# propeller_res_df3_b <- propeller_res_df3 %>%
#   filter(Comparison == "Less Affected vs. More Affected")

# Plot
sample_type_color_list <- list(`Unaffected` = "#88CCEE",
                               `Less Affected` = "pink",
                               `More Affected` = "#ac3e6f")
sample_affect_color_list <- list(`Unaffected` = "#88CCEE",
                                 `Less Affected` = "pink",
                                 `More Affected` = "#ac3e6f")
var_width = 6000000
celltypeprop_df <- mutate(celltypeprop_df, final_CT = str_wrap(final_CT, width = var_width))
ct_labels <- unique(celltypeprop_df$final_CT2)
names(ct_labels) <- unique(celltypeprop_df$final_CT)

celltype_boxplot <- celltypeprop_df %>%
  mutate(sample_affect = factor(sample_affect, level = c("Unaffected", "Less Affected", "More Affected"))) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = sample_affect, y = celltype_prop_in_sample, fill = sample_affect),
               outlier.shape = 21, size = 0.4, outlier.size = 1.25, width = 0.7,
               position = position_dodge2(preserve = "single")) +
  ylab("Proportion") +
  # scale_x_continuous(limits = c(-0.5, 2.5)) +
  theme_classic(base_size = 6.5) +
  facet_wrap(.~final_CT, scales = "free", ncol = 6, labeller = labeller(.cols = ct_labels)) +
  guides(fill = guide_legend(label.position = "top",
                             keywidth = 1.7,
                             keyheight = 1.7)) +
  theme_black +
  theme(strip.background.x = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black", size = 5.5),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = c(0.92, 0.05),
        legend.direction = "horizontal",
        strip.background = element_blank()) +
  scale_fill_manual(values = sample_type_color_list, labels = c("Unaffected", "Less\nAffected", "More\nAffected"))
celltype_boxplot

ggsave("/scratch/avannan/MANUSCRIPT_figures/Supp4_celltype_boxplot.pdf", width = 180, height = 170, units = "mm", device = "pdf")
celltype_boxplot
dev.off()


## SUPPLEMENTARY FIGURE 13 (CT Proportion Correlation with Pathology Score) ----
# Load dataframes
library(googlesheets4)
gs4_deauth()
path_corr <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094")
path_corr_sheet <- read_sheet(path_corr, sheet = 6, skip = 1)
path_corr_df <- path_corr_sheet %>%
  as.data.frame() %>%
  select(-1)
ctprop_path_df <- (table(xenium$sample, xenium$final_CT)/rowSums(table(xenium$sample, xenium$final_CT))) %>%
  as.data.frame() %>%
  dplyr::rename(sample = "Var1", final_CT = "Var2", Prop = "Freq") %>%
  full_join(xenium@meta.data %>% 
              select(sample, sample_affect, percent_pathology) %>% 
              unique()) %>%
  left_join(path_corr_df, by = c("final_CT" = "celltype")) %>%
  mutate(final_CT2 = ifelse(!is.na(`Significance level`), paste0(final_CT, " ", `Significance level`), final_CT))

# Create plot
ctprop_path_plot <- ctprop_path_df %>%
  ggplot(aes(x = percent_pathology, y = Prop)) +
  geom_point(aes(fill = final_CT), shape = 21) +
  geom_smooth(aes(color = final_CT), method = "lm") +
  theme_classic(base_size = 7) +
  theme_black +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list) +
  labs(x = "Percent Pathology", y = "Proportion") +
  facet_wrap(~final_CT2, scales = "free_y", ncol = 6) +
  NoLegend()
ctprop_path_plot

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ctprop_pathology_corr.pdf",
    width = 180*0.0393701, height = 170*0.0393701)
ctprop_path_plot
dev.off()


## SUPPLEMENTARY FIGURE 16 (Annotation Counts) ----
anno_colors <- distinctColorPalette(25)

# Number of Nuclei
num_nuclei_df <- final_annotations_df %>%
  group_by(Annotation_Type) %>%
  mutate(num_cells_annotation_type = length(full_cell_id)) %>%
  mutate(Annotation_Type = ordered(Annotation_Type, levels = new_annotation_order)) %>%
  select(Annotation_Type, num_cells_annotation_type) %>%
  unique()
tmp1 <- num_nuclei_df %>% 
  select(Annotation_Type, num_cells_annotation_type) %>% 
  mutate(Pos = ifelse(num_cells_annotation_type < 100, num_cells_annotation_type+6000,
                      ifelse(num_cells_annotation_type > 2000, num_cells_annotation_type+11000,
                             num_cells_annotation_type+8000)),
         num_cells_annotation_type) %>%
  # mutate(Pos = case_when(num_cells_annotation_type >= 5000 ~ num_cells_annotation_type+5000
  #                        num_cells_annotation_type < 5000 ~ num_cells_annotation_type+2500,
  #                        num_cells_annotation_type < 1000 ~ num_cells_annotation_type+1500)) %>% 
  unique()

a1 <- num_nuclei_df %>%
  full_join(tmp1) %>%
  arrange(desc(num_cells_annotation_type)) %>%
  ggplot(aes(y = reorder(Annotation_Type, dplyr::desc(Annotation_Type)),
             x = num_cells_annotation_type, fill = Annotation_Type)) +
  geom_col() +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), position = "bottom", limits = c(0, 220000),
                     labels = c("0", "50k", "100k", "150k", "200k")) +
  scale_fill_manual(values = anno_colors) +
  scale_y_discrete(position = "left") +
  theme_classic(base_size = 6.5) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 6),
        axis.line.x.top = element_line(color = "black"),
        panel.grid.major.x = element_line(color = "grey65"),
        plot.margin = margin(0,1,0,2, "mm")) +
  labs(x = "Number of Nuclei") +
  geom_text(aes(x = Pos, label = num_cells_annotation_type), size = 2) +
  NoLegend()
a1

# Number of Annotations
anno_num_df <- final_annotations_df %>%
  select(sample, Annotation_Type, annotation_type_instance) %>%
  unique() %>%
  group_by(Annotation_Type) %>%
  mutate(num_annotation_type = length(Annotation_Type)) %>%
  mutate(Annotation_Type = ordered(Annotation_Type, levels = new_annotation_order)) %>%
  ungroup() %>%
  select(Annotation_Type, num_annotation_type) %>%
  unique()
tmp2 <- anno_num_df %>% 
  select(Annotation_Type, num_annotation_type) %>% 
  mutate(Pos = ifelse(num_annotation_type < 50, 
                      num_annotation_type+2, num_annotation_type+3)) %>% 
  unique()

a2 <- anno_num_df %>%
  full_join(tmp2) %>%
  ggplot(aes(y = reorder(Annotation_Type, dplyr::desc(Annotation_Type)),
             x = num_annotation_type, fill = Annotation_Type)) +
  geom_col() +
  theme_classic() +
  # scale_x_continuous(expand = c(0, 0), position = "bottom", limits = c(0, 62)) +
  scale_fill_manual(values = anno_colors) +
  scale_y_discrete(position = "left") +
  theme_classic(base_size = 6.5) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 6),
        axis.line.x.top = element_line(color = "black"),
        panel.grid.major.x = element_line(color = "grey65"),
        plot.margin = margin(0,2,0,2, "mm")) +
  labs(x = "Number of Annotations") +
  geom_text(aes(x = Pos, label = num_annotation_type), size = 2) +
  NoLegend()

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_A_num_anno.pdf",
    width = 180*0.0393701, height = 85*0.0393701)
a1 + a2 + patchwork::plot_layout(ncol = 2)
dev.off()


# Per sample info
sample_order <- xenium@meta.data %>%
  select(sample, new_sample_name, percent_pathology) %>%
  unique() %>%
  arrange(percent_pathology, new_sample_name) %>%
  pull(new_sample_name)
b1 <- final_annotations_df %>%
  select(sample, annotation_type_instance, Annotation_Type) %>%
  full_join(xenium@meta.data %>% select(sample, new_sample_name, percent_pathology) %>% unique()) %>%
  unique() %>%
  filter(!is.na(Annotation_Type)) %>%
  mutate(Annotation_Type = ordered(Annotation_Type, levels = new_annotation_order),
         percent_pathology = as.numeric(percent_pathology)) %>%
  arrange(percent_pathology, new_sample_name) %>%
  ggplot(aes(x = reorder(new_sample_name, percent_pathology), fill = Annotation_Type)) +
  scale_fill_manual(values = anno_colors) +
  geom_bar() +
  theme_classic(base_size = 6.5) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = c(0.48, 0.7),
        plot.margin = margin(5.5,7,5.5,5.5, "pt")) +
  labs(y = "Number of Annotations") +
  guides(fill = guide_legend(ncol = 2, keywidth = unit(3, "mm"),
                             keyheight = unit(3, "mm"))) +
  NoLegend()
b1

b2 <- final_annotations_df %>%
  select(sample, annotation_type_instance, Annotation_Type) %>%
  full_join(xenium@meta.data %>% select(sample, new_sample_name, percent_pathology) %>% unique()) %>%
  unique() %>%
  mutate(Annotation_Type = ordered(Annotation_Type, levels = new_annotation_order),
         percent_pathology = as.numeric(percent_pathology)) %>%
  arrange(percent_pathology, new_sample_name) %>%
  filter(!is.na(Annotation_Type)) %>%
  ggplot(aes(x = reorder(new_sample_name, percent_pathology), fill = Annotation_Type)) +
  scale_fill_manual(values = anno_colors) +
  geom_bar(position = position_fill(), color = "black", width = 1) +
  theme_classic(base_size = 6.5) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title = element_blank(),
        plot.margin = margin(5.5,5.5,5.5,7, "pt")) +
  labs(y = "Proportion of Annotations") +
  NoLegend()
b2

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_B_num_anno.pdf",  width = 180*0.0393701, height = 85*0.0393701)
b1+b2
dev.off()


## SUPPLEMENTARY FIGURE 20 (Niche x Sample) ----
# TNiches
tniche_sample_prop_df <- as.data.frame((table(xenium$new_sample_name, xenium$TNiche)/
                                          rowSums(table(xenium$new_sample_name, xenium$TNiche))))
names(tniche_sample_prop_df) <- c("Sample", "Niche", "Proportion")
a <- tniche_sample_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("T", 1:12)),
         Sample = ordered(Sample, levels = sample_order)) %>%
  ggplot(aes(x = Sample, y = Proportion, fill = Niche)) +
  geom_col(color = "black", width = 1, linewidth = 0.2) +
  scale_fill_manual(values = transcript_niche_color_list) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 6.5) +
  theme_angle +
  theme_black +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm")) +
  labs(y = "Proportion of Cells")

# CNiches
cniche_sample_prop_df <- as.data.frame((table(xenium$new_sample_name, xenium$CNiche)/
                                          rowSums(table(xenium$new_sample_name, xenium$CNiche))))
names(cniche_sample_prop_df) <- c("Sample", "Niche", "Proportion")
b <- cniche_sample_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Sample = ordered(Sample, levels = sample_order)) %>%
  ggplot(aes(x = Sample, y = Proportion, fill = Niche)) +
  geom_col(color = "black", width = 1, linewidth = 0.2) +
  scale_fill_manual(values = nuclei_niche_color_list) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 6.5) +
  theme_angle +
  theme_black +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(4, "mm")) +
  labs(y = "Proportion of Cells")
b

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig12_niche_by_sample.pdf", width = (180*0.0393701), height = (168*0.0393701))
ggarrange(a, b, ncol = 1)
dev.off()


## SUPPLEMENTARY FIGURE 22 (Niche x CT Dotplots - Niche Composition) ----
lineage_df <- xenium@meta.data %>%
  select(final_CT, final_lineage) %>%
  unique()
names(lineage_df) <- c("CT", "Lineage")

# Transcript Niches
tniche_ct_prop_df <- as.data.frame((table(xenium$TNiche, xenium$final_CT)/
                                      rowSums(table(xenium$TNiche, xenium$final_CT))))
names(tniche_ct_prop_df) <- c("Niche", "CT", "Proportion")
tniche_ct_prop_plot <- tniche_ct_prop_df %>%
  filter(Niche != "TNA", Niche != "NA", !is.na(Niche)) %>%
  full_join(lineage_df) %>%
  mutate(Niche = ordered(Niche, levels = paste0("T", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  ggplot(aes(y = reorder(Niche, dplyr::desc(Niche)), x = CT,
             size = Proportion, color = as.factor(Niche),
             alpha = Proportion_Visible,
             fill = as.factor(Niche))) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw(base_size = 6) +
  scale_color_manual(values = c("black", "grey80", rep("black", 10)), guide = "none") +
  scale_fill_manual(values = transcript_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 4), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 6),
        strip.text = element_text(color = "white"),
        legend.position = "bottom") + 
  # ggh4x::force_panelsizes(cols = c(2.1, 3.8, 4.1, 3.2)) +
  ggh4x::force_panelsizes(cols = c(1.5, 4.2, 8, 3.5)) +
  labs(title = "Niche Composition") +
  scale_alpha_manual(values = c(0, 1), guide = "none")
tniche_ct_prop_plot

# Change facet colors
tniche_ct_prop_plot_table <- ggplot_gtable(ggplot_build(tniche_ct_prop_plot))
striprt <- which(grepl("strip-r", tniche_ct_prop_plot_table$layout$name) |
                   grepl("strip-t", tniche_ct_prop_plot_table$layout$name))
fills <- rep(c("chocolate2", "chartreuse4", "maroon3", "cornflowerblue"), 2)
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", tniche_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  tniche_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

# Cell Niches
cniche_ct_prop_df <- as.data.frame((table(xenium$CNiche, xenium$final_CT)/
                                      rowSums(table(xenium$CNiche, xenium$final_CT))))
names(cniche_ct_prop_df) <- c("Niche", "CT", "Proportion")

cniche_ct_prop_plot <- cniche_ct_prop_df %>%
  full_join(lineage_df) %>%
  mutate(Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  filter(!is.na(Niche)) %>%
  ggplot(aes(y = reorder(Niche, dplyr::desc(Niche)), x = CT, 
             size = Proportion, color = as.factor(Niche), 
             alpha = Proportion_Visible,
             fill = as.factor(Niche))) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw(base_size = 6) +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_color_manual(values = c(rep("black", 7), "grey80", rep("black", 4)), guide = "none") +
  scale_size_continuous(range = c(0, 4), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust = 1, size = 5),
        axis.title = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3), 
        strip.text = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        strip.background = element_blank(),
        axis.text.y = element_text(color = "black"), 
        axis.ticks.y = element_line(color = "black"),
        plot.title = element_blank(),
        legend.position = "bottom") + 
  # ggh4x::force_panelsizes(cols = c(2.1, 3.8, 4.1, 3.2)) +
  ggh4x::force_panelsizes(cols = c(1.5, 4.2, 8, 3.5)) +
  scale_alpha_manual(values = c(0, 1), guide = "none")
cniche_ct_prop_plot

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig14_nichecomposition.pdf", width = (160*0.0393701), height = (140*0.0393701))
cowplot::plot_grid(tniche_ct_prop_plot_table, cniche_ct_prop_plot, ncol = 1, align = "hv")
dev.off()


## SUPPLEMENTARY FIGURE 23 (Annotation x CT and Niche Dotplots) ----
# Create dataframe of lineage variables
lineage_df <- xenium@meta.data %>%
  select(final_CT, final_lineage) %>%
  unique()
names(lineage_df) <- c("CT", "Lineage")

new_annotation_order <- c(# Epithelial
  "Normal Alveoli",
  "Minimally Remodeled Alveoli",
  "Hyperplastic AECs",
  "Emphysema",
  "Remodeled Epithelium",
  "Advanced Remodeling",
  "Epithelial Detachment",
  "Remnant Alveoli",
  "Small Airway",
  "Large Airway",
  "Microscopic Honeycombing",
  "Goblet Cell Metaplasia",
  
  # Immune
  "Granuloma",
  "Mixed Inflammation",
  "TLS",
  
  # Endothelial/Mesenchymal
  "Artery",
  "Muscularized Artery",
  "Venule",
  "Interlobular Septum",
  "Airway Smooth Muscle",
  "Fibroblastic Focus",
  "Fibrosis",
  "Severe Fibrosis",
  
  # Other
  "Multinucleated Cell",
  "Giant Cell")

# A: Annotation Composition
path_ct_prop_df <- as.data.frame((table(final_annotations_df$Annotation_Type, final_annotations_df$final_CT)/
                                    rowSums(table(final_annotations_df$Annotation_Type, final_annotations_df$final_CT))))
names(path_ct_prop_df) <- c("Pathology", "CT", "Proportion")

full_path_ct_prop_plot <- path_ct_prop_df %>%
  left_join(lineage_df) %>%
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  mutate(CT = ordered(CT, levels = new_ct_order)) %>%
  filter(Proportion >= 0.01) %>%
  ggplot(aes(x = Pathology,
             y = reorder(CT, dplyr::desc(CT)),
             size = Proportion, fill = CT)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6.5) +
  scale_fill_manual(values = color_list) +
  scale_size_continuous(range = c(0, 3), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text = element_text(color = "black", size = 5),
        strip.text = element_text(color = "white", size = 5),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5),
        panel.grid.major.y = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, size = 6),
        legend.position = c(-0.3, -0.155),
        legend.box.background = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.background = element_blank(),
        plot.margin = margin(5.5, 5.5, 8, 8, "pt"),
        legend.margin = margin(1.5, 1.5, 1.5, 1.5, "pt")) +
  theme_angle +
  labs(title = "Annotation Composition (Cells)") +
  guides(fill = "none", 
         size = guide_legend(title.theme = element_text(size = 5),
                             size = 2,
                             nrow = 6, ncol = 1, title.hjust = 0.5,
                             label.position = "right",
                             label.theme = element_text(size = 3.5))) +
  facet_grid(Lineage~., space = "free", scales = "free")
full_path_ct_prop_plot

# Change facet colors
full_path_ct_prop_plot_table <- ggplot_gtable(ggplot_build(full_path_ct_prop_plot))
striprt <- which(grepl("strip-r", full_path_ct_prop_plot_table$layout$name) | 
                   grepl("strip-t", full_path_ct_prop_plot_table$layout$name))
fills <- c("chocolate2", "chartreuse4", "maroon3", "cornflowerblue") # Sample Types
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", full_path_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  full_path_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(full_path_ct_prop_plot_table)

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFigA_full_ct_anno.pdf", width = (92*0.0393701), height = (114*0.0393701))
grid::grid.draw(full_path_ct_prop_plot_table)
dev.off()


# B1: Transcript Niches (Annotation Composition)
tniche_path_prop_df <- as.data.frame((table(final_annotations_df$Annotation_Type, final_annotations_df$TNiche)/
                                        rowSums(table(final_annotations_df$Annotation_Type, final_annotations_df$TNiche))))
names(tniche_path_prop_df) <- c("Pathology", "Niche", "Proportion")
tniche_path_prop_plot <- tniche_path_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("T", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  ggplot(aes(x = Niche, y = reorder(Pathology, dplyr::desc(Pathology)), 
             size = Proportion, alpha = Proportion_Visible, fill = Niche)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6.5) +
  scale_fill_manual(values = transcript_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 3), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 5), axis.ticks.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 2, "pt"),
        legend.key = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0, "pt"),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(0, 0, 0, 0, "pt")),
        legend.spacing = unit(0, "pt"),
        legend.position = "bottom", legend.title.align = 0.5) +
  guides(size = guide_legend(title.position = "top", label.position = "bottom")) +
  scale_alpha_manual(values = c(0, 1), guide = "none")
tniche_path_prop_plot

# B2: Cell Niches (Annotation Composition)
cniche_path_prop_df <- as.data.frame((table(final_annotations_df$Annotation_Type, final_annotations_df$CNiche)/
                                        rowSums(table(final_annotations_df$Annotation_Type, final_annotations_df$CNiche))))
names(cniche_path_prop_df) <- c("Pathology", "Niche", "Proportion")
cniche_path_prop_plot <- cniche_path_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  ggplot(aes(x = Niche, y = reorder(Pathology, dplyr::desc(Pathology)), 
             size = Proportion, alpha = Proportion_Visible, fill = Niche)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6.5) +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 3), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 2, "pt"),
        legend.key = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0, "pt"),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(0, 0, 0, 0, "pt")),
        legend.spacing = unit(0, "pt"),
        legend.position = "bottom", legend.title.align = 0.5,
        legend.key.width = unit(3.5, "mm"),
        legend.key.height = unit(2, "mm")) +
  guides(size = guide_legend(title.position = "top", label.position = "bottom")) +
  scale_alpha_manual(values = c(0, 1), guide = "none")
cniche_path_prop_plot

# Combine and save
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFigB_full_niche_anno.pdf", width = (90*0.0393701), height = (125*0.0393701))
(tniche_path_prop_plot + cniche_path_prop_plot) + 
  patchwork::plot_layout(widths = unit(c(30, 30), "mm"),
                         heights = unit(c(100, 100), "mm"),
                         ncol = 2, nrow = 1)
dev.off()



## SUPPLEMENTARY FIGURE 27 (Cell Proximity Heatmap for C3) ----
logreg_cniche <- read.csv(file.path(workDir, 'LogReg_ProximityResults_cNiche.csv')) %>%
  mutate(Significance = ifelse(sig == "fdr < 0.05", "FDR < 0.05", "n.s."))
logreg_cniche$niche <- factor(logreg_cniche$niche, levels = paste0('C', 1:12))

## Cell niches
prob <- reshape2::dcast(niche + source ~ target, value.var = 'prob', data = logreg_cniche) %>% split(., .$niche)
pval <- reshape2::dcast(niche + source ~ target, value.var = 'fdr', data = logreg_cniche) %>% split(., .$niche)
or <- reshape2::dcast(niche + source ~ target, value.var = 'logOR', data = logreg_cniche) %>% split(., .$niche)

prob <- lapply(prob, function(x) {
  zz <- x
  rownames(zz) <- zz$source
  zz$source <- NULL
  zz$niche <- NULL
  return(zz)
})

pval <- lapply(pval, function(x) {
  zz <- x
  rownames(zz) <- zz$source
  zz$source <- NULL
  zz$niche <- NULL
  return(zz)
})

or <- lapply(or, function(x) {
  zz <- x
  rownames(zz) <- zz$source
  zz$source <- NULL
  zz$niche <- NULL
  return(zz)
})

ht_opt$HEATMAP_LEGEND_PADDING = unit(5, "mm")
c_ht_list <- list()
for(i in 3){
  
  niche_id <- paste0('C', i)
  table_name <- paste0(paste0('C', i), '\n','log(OR)')
  ct_prob<- prob[[niche_id]] %>% as.matrix()
  ct_prob[is.na(ct_prob)] = 0
  
  ct_pval <- pval[[niche_id]] %>% as.matrix()
  ct_pval[is.na(ct_pval)] = 0
  
  ct_or <- or[[niche_id]] %>% as.matrix()
  ct_or[is.na(ct_or)] = 0
  
  mean(rownames(ct_prob) == rownames(ct_pval))
  mean(colnames(ct_prob) == colnames(ct_pval))
  mean(rownames(ct_or) == rownames(ct_pval))
  mean(colnames(ct_or) == colnames(ct_pval))
  
  ct_prob[ct_pval > 0.05] = 0
  ct_or[ct_pval > 0.05] = 0
  
  # check if pvalue rownames and colnames match to those in or table before masking p>0.05 values
  mean(rownames(ct_or) == rownames(ct_pval))
  mean(colnames(ct_or) == colnames(ct_pval))
  # mask OR for values that are p-val > 0.05
  ct_or[ct_pval > 0.05] = 0
  
  # Remove any columns/rows where all values are 0 or NA
  remove_cols <- which(colSums(ct_or, na.rm = TRUE) == 0)
  remove_rows <- which(rowSums(ct_or, na.rm = TRUE) == 0)
  ct_or <- ct_or[-remove_rows, -remove_cols]
  
  # plot odds ratio
  min01 <- quantile(as.matrix(ct_or), 0.01, na.rm = TRUE)
  max99 <- quantile(as.matrix(ct_or), 0.99, na.rm = TRUE)
  col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))
  
  ht <- Heatmap(as.matrix(ct_or), name = table_name,
                      heatmap_legend_param = list(legend_direction = "vertical",
                                                  title_position = "topcenter",
                                                  title_gp = gpar(fontsize = 6),
                                                  labels_gp = gpar(fontsize = 6),
                                                  background = NA,
                                                  grid_width = unit(2, "mm"),
                                                  grid_height = unit(6, "mm"),
                                                  legend_height = unit(10, "mm"),
                                                  legend_width = unit(2, "mm")),
                      clustering_distance_columns = , clustering_distance_rows = ,
                      clustering_method_columns = , clustering_method_rows = ,
                      row_title = 'Starting Cell', row_title_side = 'left', row_title_gp = gpar(fontsize = 7),
                      row_names_gp = gpar(fontsize = 6),
                      row_dend_side = "right", row_names_side = "left",
                      column_title = 'Nearest Neighbor Cell', column_names_side = 'top', column_dend_side = 'bottom',  column_title_gp = gpar(fontsize = 7),
                      column_names_gp = gpar(fontsize = 6),
                      rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
                      width = unit(128.1, "mm"), height = unit(128.1, "mm"),
                      na_col = 'black', col = col_fun, column_names_rot = 45, row_dend_width = unit(0.5, 'cm'), column_dend_height = unit(0.5, 'cm'))
  c_ht_list[[table_name]] <- ComplexHeatmap::draw(ht)
}
draw(c_ht_list[[1]])
ht_opt(RESET=TRUE)

png("/scratch/avannan/MANUSCRIPT_figures/SuppFig_C3_prox_results.png",
    width = (180*0.0393701), height = (165*0.0393701), unit = "in", res = 450)
draw(c_ht_list[[1]])
dev.off()


## SUPPLEMENTARY FIGURE 28 (Cell Proximity Heatmap for T3) ----
logreg_tniche <- read.csv(file.path(workDir, 'LogReg_ProximityResults_tNiche.csv')) %>%
  mutate(Significance = ifelse(sig == "fdr < 0.05", "FDR < 0.05", "n.s."))
logreg_tniche$niche <- factor(logreg_tniche$niche, levels = paste0('T', 1:12))

## Cell niches
prob <- reshape2::dcast(niche + source ~ target, value.var = 'prob', data = logreg_tniche) %>% split(., .$niche)
pval <- reshape2::dcast(niche + source ~ target, value.var = 'fdr', data = logreg_tniche) %>% split(., .$niche)
or <- reshape2::dcast(niche + source ~ target, value.var = 'logOR', data = logreg_tniche) %>% split(., .$niche)

prob <- lapply(prob, function(x) {
  zz <- x
  rownames(zz) <- zz$source
  zz$source <- NULL
  zz$niche <- NULL
  return(zz)
})

pval <- lapply(pval, function(x) {
  zz <- x
  rownames(zz) <- zz$source
  zz$source <- NULL
  zz$niche <- NULL
  return(zz)
})

or <- lapply(or, function(x) {
  zz <- x
  rownames(zz) <- zz$source
  zz$source <- NULL
  zz$niche <- NULL
  return(zz)
})

ht_opt$HEATMAP_LEGEND_PADDING = unit(5, "mm")
t_ht_list <- list()
for(i in 3){
  
  niche_id <- paste0('T', i)
  table_name <- paste0(paste0('T', i), '\n','log(OR)')
  ct_prob<- prob[[niche_id]] %>% as.matrix()
  ct_prob[is.na(ct_prob)] = 0
  
  ct_pval <- pval[[niche_id]] %>% as.matrix()
  ct_pval[is.na(ct_pval)] = 0
  
  ct_or <- or[[niche_id]] %>% as.matrix()
  ct_or[is.na(ct_or)] = 0
  
  mean(rownames(ct_prob) == rownames(ct_pval))
  mean(colnames(ct_prob) == colnames(ct_pval))
  mean(rownames(ct_or) == rownames(ct_pval))
  mean(colnames(ct_or) == colnames(ct_pval))
  
  ct_prob[ct_pval > 0.05] = 0
  ct_or[ct_pval > 0.05] = 0
  
  # check if pvalue rownames and colnames match to those in or table before masking p>0.05 values
  mean(rownames(ct_or) == rownames(ct_pval))
  mean(colnames(ct_or) == colnames(ct_pval))
  # mask OR for values that are p-val > 0.05
  ct_or[ct_pval > 0.05] = 0
  
  # Remove any columns/rows where all values are 0 or NA
  remove_cols <- which(colSums(ct_or, na.rm = TRUE) == 0)
  remove_rows <- which(rowSums(ct_or, na.rm = TRUE) == 0)
  ct_or <- ct_or[-remove_rows, -remove_cols]
  
  # plot odds ratio
  min01 <- quantile(as.matrix(ct_or), 0.01, na.rm = TRUE)
  max99 <- quantile(as.matrix(ct_or), 0.99, na.rm = TRUE)
  col_fun = colorRamp2(c(min01, 0, max99), c("#377EB8", "white", "#E41A1C"))
  
  ht <- Heatmap(as.matrix(ct_or), name = table_name,
                heatmap_legend_param = list(legend_direction = "vertical",
                                            title_position = "topcenter",
                                            title_gp = gpar(fontsize = 6),
                                            labels_gp = gpar(fontsize = 6),
                                            background = NA,
                                            grid_width = unit(2, "mm"),
                                            grid_height = unit(6, "mm"),
                                            legend_height = unit(10, "mm"),
                                            legend_width = unit(2, "mm")),
                clustering_distance_columns = , clustering_distance_rows = ,
                clustering_method_columns = , clustering_method_rows = ,
                row_title = 'Starting Cell', row_title_side = 'left', row_title_gp = gpar(fontsize = 7),
                row_names_gp = gpar(fontsize = 6),
                row_dend_side = "right", row_names_side = "left",
                column_title = 'Nearest Neighbor Cell', column_names_side = 'top', column_dend_side = 'bottom',  column_title_gp = gpar(fontsize = 7),
                column_names_gp = gpar(fontsize = 6),
                rect_gp = gpar(col = "#b8b7b6", lwd = 0.5),
                width = unit(128.1, "mm"), height = unit(128.1, "mm"),
                na_col = 'black', col = col_fun, column_names_rot = 45, row_dend_width = unit(0.5, 'cm'), column_dend_height = unit(0.5, 'cm'))
  t_ht_list[[table_name]] <- ComplexHeatmap::draw(ht)
}
draw(t_ht_list[[1]])
ht_opt(RESET=TRUE)

png("/scratch/avannan/MANUSCRIPT_figures/SuppFig_T3_prox_results.png",
    width = (180*0.0393701), height = (165*0.0393701), unit = "in", res = 450)
draw(t_ht_list[[1]])
dev.off()


## SUPPLEMENTARY FIGURE 31 (Heatmap of Intermediate/Late Remodeling Expression) ----
# Load in lumen genes and their order
gene_mode_and_sc_clusters <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/gene_mode_and_sc_clusters.csv")

# Load in significant terms x CT and filter
library(googlesheets4)
gs4_deauth()
sig_terms_perCT_sheet <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094")
sig_terms_perCT_sheet <- read_sheet(sig_terms_perCT_sheet, sheet = 17, skip = 1)
sig_terms_perCT <- as.data.frame(sig_terms_perCT_sheet) %>%
  filter(padj < 0.05) %>%
  full_join(gene_mode_and_sc_clusters, by = c("gene" = "genes")) %>%
  mutate(total_sig = length(gene)) %>%
  group_by(celltype) %>%
  mutate(prop_sig_ct = length(gene)/total_sig) %>%
  ungroup()

# Load in lumen metadata and list of final lumens
lumen_nuclei_RNA_sce <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_lumen_order_t4/lumen_nuclei_RNA_sce_t4_ordered.rds")
lumen_metadata <- xenium@meta.data %>%
  filter(lumen_id %in% colnames(lumen_nuclei_RNA_sce))

# Subset object to contain only cells found in final lumens
xenium_lumens_cells_only <- subset(xenium, cells = rownames(lumen_metadata))

# Pick out genes for heatmap
heatmap_genes <- sig_terms_perCT %>%
  filter(sc_clusters %in% c(1, 2)) %>%
  arrange(desc(sc_clusters), gene) %>%
  pull(gene) %>%
  unique()
heatmap1_genes <- sig_terms_perCT %>% # Intermediate
  filter(sc_clusters == 1) %>%
  arrange(desc(sc_clusters), gene) %>%
  pull(gene) %>%
  unique()
heatmap2_genes <- sig_terms_perCT %>% # Late
  filter(sc_clusters == 2) %>%
  arrange(desc(sc_clusters), gene) %>%
  pull(gene) %>%
  unique()

# Log-normalizing and scaling all features in the RNA assay
# Scaling so that all features can be visualized using the same color scale
xenium_lumens_cells_only <- ScaleData(xenium_lumens_cells_only)

# Create dotplot base
p <- DotPlot(xenium_lumens_cells_only, features = heatmap_genes, 
             group.by = "final_CT", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p

# Reproducing the dotplot using ComplexHeatmap
df <- p$data
head(df)

# (Scaled) expression levels
exp_mat <- df %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  tidyr::pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()
dim(exp_mat)

# Get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99), na.rm=T)

# Replace NA with 0 in matrix
exp_mat[!is.finite(exp_mat)] <- 0
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Create dataframe with significance per lineage, as 1/0; use this as row annotation
# TRUE = Significant across pseudotime for at least one CT in that lineage
gene_lineage_lumen_heatmap_df_og <- sig_terms_perCT %>%
  dplyr::rename(final_CT = "celltype") %>%
  mutate(final_CT = gsub("\\.", " ", final_CT),
         final_CT = case_when(final_CT == "KRT5  KRT17 " ~ "KRT5-/KRT17+",
                              final_CT == "SPP1  Macrophages" ~ "SPP1+ Macrophages",
                              final_CT == "Monocytes MDMs" ~ "Monocytes/MDMs" ,
                              final_CT == "Macrophages   IFN activated" ~ "Macrophages - IFN-activated",
                              final_CT == "SMCs Pericytes" ~ "SMCs/Pericytes",
                              final_CT == "NK NKT" ~ "NK/NKT",
                              final_CT == "CD4  T cells" ~ "CD4+ T-cells",
                              final_CT == "CD8  T cells" ~ "CD8+ T-cells",
                              TRUE ~ final_CT)) %>%
  left_join(xenium@meta.data %>% select(final_CT, final_lineage, final_sublineage) %>% unique()) %>%
  mutate(lineage = ifelse(final_lineage == "Epithelial", "Epithelial", final_sublineage)) %>%
  mutate(sig = 1) %>%
  filter(!is.na(final_CT)) %>%
  filter(!is.na(gene)) %>%
  unique()

# Get proportion/percentage of significant tests
prop_gene_lumen_lineage_df <- gene_lineage_lumen_heatmap_df_og %>%
  group_by(sc_clusters) %>%
  mutate(num_sig_tests = sum(sig)) %>%
  ungroup() %>%
  group_by(lineage, sc_clusters) %>%
  mutate(num_sig_lineage = sum(sig),
         prop_sig_lineage = sum(sig)/num_sig_tests) %>%
  select(final_lineage, prop_sig_lineage, num_sig_lineage, num_sig_tests) %>%
  unique() %>%
  mutate(title = paste0(num_sig_lineage, "/", num_sig_tests, " (", round(prop_sig_lineage*100, 1), "%)"))
prop_gene_lumen_lineage_df

# Prep dataframe for heatmap
gene_lineage_lumen_heatmap_df <- gene_lineage_lumen_heatmap_df_og %>%
  select(gene, sc_clusters, lineage, sig) %>%
  unique() %>% 
  pivot_wider(names_from = "lineage", values_from = "sig") %>%
  arrange(desc(sc_clusters), gene) %>% 
  filter(!is.na(gene))
# mutate(across(everything(), ~ ifelse(.x == 0, NA, .x)))
gene_lineage_lumen_heatmap_df2 <- gene_lineage_lumen_heatmap_df %>%
  select(Endothelial, Epithelial, Lymphoid, Myeloid, Mesenchymal) %>%
  as.data.frame()
# mutate(across(everything(), ~ ifelse(.x == 0, NA, .x)))
rownames(gene_lineage_lumen_heatmap_df2) <- gene_lineage_lumen_heatmap_df$gene
# Including genes not significant for any lineage
nonsig_genes <- heatmap_genes[!(heatmap_genes %in% gene_lineage_lumen_heatmap_df$gene)]
nonsig_gene_lineage_lumen_heatmap_df <- data.frame(matrix(NA, nrow = length(nonsig_genes), ncol = 5))
rownames(nonsig_gene_lineage_lumen_heatmap_df) <- nonsig_genes
colnames(nonsig_gene_lineage_lumen_heatmap_df) <- colnames(gene_lineage_lumen_heatmap_df2)

# Join together
gene_lineage_lumen_heatmap_df3 <- rbind(gene_lineage_lumen_heatmap_df2, 
                                        nonsig_gene_lineage_lumen_heatmap_df)[rownames(exp_mat), ]

# Create heatmap
# New CT order and color list
new_ct_order3 <- c(
  # Endothelial 
  "Arteriole", "Capillary", "Venous", "Lymphatic",
  # Epithelial - Alveolar
  "AT1", "Transitional AT2", "AT2", "Proliferating AT2", "KRT5-/KRT17+",              
  # Epithelial - Airway
  "Basal", "Multiciliated", "Goblet", "RASC", "Secretory", "PNEC", "Proliferating Airway",
  # Immune - Lymphoid
  "B cells", 
  # "Proliferating B cells", # Not tested
  "NK/NKT", "Proliferating NK/NKT", "Tregs", 
  "CD4+ T-cells", "CD8+ T-cells", "Proliferating T-cells", "Plasma", "pDCs",
  # Immune - Myeloid
  "Basophils", "cDCs", "Migratory DCs", "Langerhans cells", "Mast", 
  "Alveolar Macrophages", "Interstitial Macrophages", "SPP1+ Macrophages",
  "Macrophages - IFN-activated", "Monocytes/MDMs", "Neutrophils", "Proliferating Myeloid",
  # Mesenchymal
  "Activated Fibrotic FBs", "Adventitial FBs", "Alveolar FBs", "Inflammatory FBs", 
  "Myofibroblasts", "Subpleural FBs", "Proliferating FBs", "SMCs/Pericytes", "Mesothelial")
new_color_list3 <- color_list[new_ct_order3]

# Get row split variable (clusters)
sc_cluster_order_genes <- sig_terms_perCT %>%
  filter(gene %in% rownames(exp_mat)) %>%
  select(gene, sc_clusters) %>%
  unique() %>%
  arrange(desc(sc_clusters), gene) %>%
  as.vector()
all.equal(sc_cluster_order_genes$gene, rownames(exp_mat)) # Must be TRUE
rowsplit <- ordered(as.factor(sc_cluster_order_genes$sc_clusters), levels = c("1", "2"))

ht_opt$TITLE_PADDING = unit(1.2, "mm")
hp <- Heatmap(exp_mat[, new_ct_order3], # t(exp_mat)
              name = "Scaled Expression",
              width = unit(50, "mm"),
              heatmap_legend_param = list(title = "Scaled Expression",
                                          legend_direction = "horizontal",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_width = unit(2, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),
              column_title_side = "top",
              row_title = " ",
              column_title_gp = gpar(fontsize = 5.5, col = c("chocolate2", "chartreuse4",
                                                             "darkmagenta", "deeppink3", "cornflowerblue")),
              column_names_side = "top",
              show_column_names = FALSE,
              column_names_rot = 90,
              row_names_side = "left",
              column_split = c(rep(" Endo. ", 4), rep(" Epithelial ", 12),
                               rep(" Lymphoid ", 9), rep(" Myeloid ", 12), rep("Mes.", 9)),
              row_split = rowsplit,
              row_names_gp = gpar(fontsize = 4),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_order = new_ct_order3,
              row_dend_width = unit(3, "mm"),
              row_dend_side = "right",
              border = "black",
              top_annotation = 
                columnAnnotation(`Cell Type` = colnames(exp_mat[, new_ct_order3]),
                                 col = list(`Cell Type` = unlist(new_color_list3)),
                                 show_legend = FALSE,
                                 annotation_name_gp = gpar(fontsize = 0),
                                 height = unit(1.2, "mm"),
                                 simple_anno_size = unit(1.2, "mm")),
              
              right_annotation =
                rowAnnotation(
                  `- 13.4%` = anno_points(gene_lineage_lumen_heatmap_df3$Mesenchymal, which = "row",
                                          ylim = c(0.99, 1.01), extend = 0, 
                                          size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                          pch = 15,
                                          axis = FALSE, border = FALSE, gp = gpar(col = "cornflowerblue")),
                  `- 37.8%` = anno_points(gene_lineage_lumen_heatmap_df3$Myeloid, which = "row",
                                          ylim = c(0.99, 1.01), extend = 0,
                                          size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                          pch = 15,
                                          axis = FALSE, border = FALSE, gp = gpar(col = "deeppink3")),
                  `- 1.7% ` = anno_points(gene_lineage_lumen_heatmap_df3$Lymphoid, which = "row",
                                          ylim = c(0.99, 1.01), extend = 0,
                                          size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                          pch = 15,
                                          axis = FALSE, border = FALSE, gp = gpar(col = "darkmagenta")),
                  `- 34.5%` = anno_points(gene_lineage_lumen_heatmap_df3$Epithelial, which = "row",
                                          ylim = c(0.99, 1.01), extend = 0,
                                          size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                          pch = 15,
                                          axis = FALSE, border = FALSE, gp = gpar(col = "chartreuse4")),
                  `- 12.6%` = anno_points(gene_lineage_lumen_heatmap_df3$Endothelial, which = "row",
                                          ylim = c(0.99, 1.01), extend = 0,
                                          size = unit(1.4, "mm"), height = unit(0.8, "mm"), width = unit(0.24, "mm"),
                                          pch = 15,
                                          axis = FALSE, border = FALSE, gp = gpar(col = "chocolate2")),
                  annotation_name_gp = gpar(fontsize = 5),
                  show_annotation_name = TRUE,
                  annotation_name_rot = -90,
                  annotation_name_side = "bottom",
                  gap = unit(1.4, "mm"),
                  "     " = rowsplit,
                  col = list("     " = c("1" = "#FF6A0F", "2" = "#D81B60")),
                  show_legend = FALSE,
                  simple_anno_size = unit(2, "mm")
                ),
              gap = unit(1, "mm")
)

# Add annotations to heatmap
ht <- draw(hp, 
           heatmap_legend_side = "bottom",
           align_heatmap_legend = "heatmap_center")
ht

pdf("/scratch/avannan/MANUSCRIPT_figures/intermediate_late_remodeling_heatmap_supp.pdf", width = (80*0.0393701), height = (180.6*0.0393701))
ht
dev.off()
#ht_opt(RESET = TRUE)


# col = list(Trajectory = c(" Normal " = "#04846E", 
#                           "Early Transition" = "#FFC107", 
#                           "#FF6A0F",
#                           "Late Remodeling" = "#D81B60")



## RENAME SAMPLES FOR FINAL OBJECT ----
rename_samples_meta <- xenium@meta.data %>%
  # Initial rename
  mutate(new_sample_name = case_when(sample_affect == "Less Affected" ~ paste0(patient, "LA"),
                                     sample_affect == "More Affected"  ~ paste0(patient, "MA"),
                                     TRUE ~ sample)) %>%
  # Fix additional names
  mutate(new_sample_name = case_when(sample == "TILD117MF" ~ "TILD117MA1",
                                     sample == "TILD117MFB" ~ "TILD117MA2",
                                     sample == "VUILD104LF" ~ "VUILD104MA1",
                                     sample == "VUILD104MF" ~ "VUILD104MA2",
                                     sample == "VUILD105LF" ~ "VUILD105MA1",
                                     sample == "VUILD105MF" ~ "VUILD105MA2",
                                     sample == "VUILD48LF" ~ "VUILD48LA1",
                                     sample == "VUILD48MF" ~ "VUILD48LA2",
                                     TRUE ~ new_sample_name)) %>%
  dplyr::rename("old_sample_name" = "sample", "sample" = "new_sample_name") %>%
  # Fix cell ids
  mutate(full_cell_id = paste0(sample, "_", cell_id))
rownames(rename_samples_meta) <- rename_samples_meta$full_cell_id


# Update metadata
xenium@meta.data <- rename_samples_meta %>%
  separate(lumen_id, into = c("tmp", "lumen_id")) %>%
  mutate(lumen_id = paste0(sample, "_", lumen_id)) %>%
  select(sample, patient, cell_id, full_cell_id,
         sample_type, sample_affect, disease_status, percent_pathology,
         tma, run, 
         final_CT, final_lineage,
         CNiche, TNiche, lumen_id, lumen_rank,
         x_centroid, y_centroid, adj_x_centroid, adj_y_centroid,
         super_adj_x_centroid, super_adj_y_centroid,
         nCount_RNA, nFeature_RNA,
         perc_negcontrolprobe, perc_negcontrolcodeword,
         perc_unassigned, perc_negcontrolorunassigned)

# Re-save object (final - no annotations - RENAMED)
# saveRDS(xenium, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_filtered_lumens_RENAMED_082024.rds")
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_filtered_lumens_RENAMED_082024.rds")



