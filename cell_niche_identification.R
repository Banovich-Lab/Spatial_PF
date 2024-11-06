############################################
# Cell Niches
# Author: Annika Vannan (avannan@tgen.org)
# Date: 06/21/2024
# Description: Cell niches with Seurat v5
############################################

## SET ENVIRONMENT AND LOAD DATA ----
library(Seurat)
library(ggplot2)
library(tidyverse)
library(Cairo)
set.seed(0309)

# Better UMAP function
pretty_umap <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     axis.title = element_text(hjust = 1))

# Load data
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/all_celltypes_final_RAPIDS_clustered_062024.rds")

# Add spatial dimension reduction object separately
position_xy <- cbind(xenium$adj_x_centroid, xenium$adj_y_centroid)
row.names(position_xy) <- row.names(xenium@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
xenium[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_", assay = DefaultAssay(xenium))

# Set directories for saving plots
scratch_dir <- "/scratch/avannan/MANUSCRIPT_figures/cell_niches/"

## CREATE COLOR LISTS ----
color_list <- list(# Endothelial
  `Arteriole` = "#d59445",
  `Capillary` = "#cdbf5e",
  `Lymphatic` = "#975f47",
  `Venous` = "#e46300",
  
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
  `cDCs/Langerhans cells` =  "#cb547b",
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

nuclei_niche_color_list <- list(`1` = "#ffa26d",
                                `2` = "#1036b5",
                                `3` = "#cabd00",
                                `4` = "#a6513c",
                                `5` = "#84cd5d",
                                `6` = "#72c7ff",
                                `7` = "#d10040",
                                `8` = "black",
                                `9` = "#71c8a5",
                                `10` = "#6f774b",
                                `11` = "#8E76DB",
                                `12` = "#924373")

nuclei_niche_color_list <- list(`7` = "#ffa26d",
                                `1` = "#1036b5",
                                `12` = "#cabd00",
                                `9` = "#a6513c",
                                `4` = "#84cd5d",
                                `11` = "#72c7ff",
                                `3` = "#d10040",
                                `10` = "black",
                                `6` = "#71c8a5",
                                `2` = "#6f774b",
                                `5` = "#8E76DB",
                                `8` = "#924373")


## ADJUST COORDINATES ----
DimPlot(xenium, group.by = "sample", label = TRUE, reduction = "sp") + NoLegend() + coord_equal()
# Set spatial rows and columns
xenium$spatial_row <- ""
xenium$spatial_row[xenium$sample %in% c("VUHD116A", "VUILD102LF", "VUHD095", "VUILD105MF", "VUHD069", "TILD167LF", "VUILD142", "TILD130MF")] <- "R1"
xenium$spatial_row[xenium$sample %in% c("VUILD96MF", "VUILD102MF", "VUHD116B", "VUILD104LF", "VUHD113", "VUILD105LF", "VUILD58", "TILD117MFB", "VUILD141")] <- "R2"
xenium$spatial_row[xenium$sample %in% c("VUILD96LF", "VUILD107MF", "VUILD48MF", "VUILD104MF", "VUILD48LF", "TILD080MF", "VUILD49", "TILD113MF")] <- "R3"
xenium$spatial_row[xenium$sample %in% c("VUHD090", "TILD111LF", "VUHD049")] <- "R4"
xenium$spatial_row[xenium$sample %in% c("VUILD115", "VUILD110", "THD0011", "TILD117MF", "TILD315LF", "TILD299LF", "TILD049LF")] <- "R5"
xenium$spatial_row[xenium$sample %in% c("VUILD78MF", "TILD117LF", "VUILD91MF", "TILD028MF", "VUHD038")] <- "R6"
xenium$spatial_row[xenium$sample %in% c("VUILD106", "THD0008", "VUILD91LF", "VUILD78LF", "TILD175")] <- "R7"
xenium$spatial_column <- ""
xenium$spatial_column[xenium$sample %in% c("VUILD96MF", "VUILD115", "VUILD106")] <- "C1"
xenium$spatial_column[xenium$sample %in% c("VUHD116A", "VUILD102MF", "VUILD96LF")] <- "C2"
xenium$spatial_column[xenium$sample %in% c("VUILD102LF", "VUHD116B", "VUILD107MF", "VUILD110", "THD0008")] <- "C3"
xenium$spatial_column[xenium$sample %in% c("VUHD095", "VUILD104LF", "VUILD48MF", "VUILD78MF", "VUILD91LF")] <- "C4"
xenium$spatial_column[xenium$sample %in% c("VUILD105MF", "VUHD113", "VUILD104MF", "THD0011", "TILD117LF", "VUILD78LF")] <- "C5"
xenium$spatial_column[xenium$sample %in% c("VUHD069", "VUILD105LF", "VUILD48LF", "TILD117MF", "VUILD91MF", "TILD175")] <- "C6"
xenium$spatial_column[xenium$sample %in% c("TILD167LF", "VUILD58", "TILD080MF", "VUHD090", "TILD315LF", "TILD028MF")] <- "C7"
xenium$spatial_column[xenium$sample %in% c("VUILD142", "TILD117MFB", "VUILD49", "TILD111LF", "TILD299LF", "VUHD038")] <- "C8"
xenium$spatial_column[xenium$sample %in% c("TILD130MF", "VUILD141", "TILD113MF", "VUHD049", "TILD049LF")] <- "C9"

# Adjust coordinates based on row and column positions
xenium$super_adj_y_centroid <- xenium$adj_y_centroid
xenium$super_adj_y_centroid[xenium$spatial_row == "R1"] <- xenium$adj_y_centroid[xenium$spatial_row == "R1"] + 5000
xenium$super_adj_y_centroid[xenium$spatial_row == "R2"] <- xenium$adj_y_centroid[xenium$spatial_row == "R2"] + 2500
xenium$super_adj_y_centroid[xenium$spatial_row == "R4"] <- xenium$adj_y_centroid[xenium$spatial_row == "R4"] - 2000
xenium$super_adj_y_centroid[xenium$spatial_row == "R5" & xenium$tma != "TMA5"] <- xenium$adj_y_centroid[xenium$spatial_row == "R5" & xenium$tma != "TMA5"] - 3500
xenium$super_adj_y_centroid[xenium$spatial_row == "R5" & xenium$tma == "TMA5"] <- xenium$adj_y_centroid[xenium$spatial_row == "R5" & xenium$tma == "TMA5"] - 4200
xenium$super_adj_y_centroid[xenium$spatial_row == "R6" & xenium$tma != "TMA5"] <- xenium$adj_y_centroid[xenium$spatial_row == "R6" & xenium$tma != "TMA5"] - 5000
xenium$super_adj_y_centroid[xenium$spatial_row == "R6" & xenium$tma == "TMA5"] <- xenium$adj_y_centroid[xenium$spatial_row == "R6" & xenium$tma == "TMA5"] - 6500
xenium$super_adj_y_centroid[xenium$spatial_row == "R7"] <- xenium$adj_y_centroid[xenium$spatial_row == "R7"] - 7000
xenium$super_adj_x_centroid <- xenium$adj_x_centroid
xenium$super_adj_x_centroid[xenium$spatial_column == "C1"] <- xenium$adj_x_centroid[xenium$spatial_column == "C1"] - 4200
xenium$super_adj_x_centroid[xenium$spatial_column == "C2"] <- xenium$adj_x_centroid[xenium$spatial_column == "C2"] - 3000
xenium$super_adj_x_centroid[xenium$spatial_column == "C3"] <- xenium$adj_x_centroid[xenium$spatial_column == "C3"] - 1000
xenium$super_adj_x_centroid[xenium$spatial_column == "C5"] <- xenium$adj_x_centroid[xenium$spatial_column == "C5"] + 2500
xenium$super_adj_x_centroid[xenium$spatial_column == "C6"] <- xenium$adj_x_centroid[xenium$spatial_column == "C6"] + 5000
xenium$super_adj_x_centroid[xenium$spatial_column == "C7"] <- xenium$adj_x_centroid[xenium$spatial_column == "C7"] + 7500
xenium$super_adj_x_centroid[xenium$spatial_column == "C8"] <- xenium$adj_x_centroid[xenium$spatial_column == "C8"] + 10000
xenium$super_adj_x_centroid[xenium$spatial_column == "C9"] <- xenium$adj_x_centroid[xenium$spatial_column == "C9"] + 12500
# Add "super" spatial dimension reduction object
position_xy <- cbind(xenium$super_adj_x_centroid, xenium$super_adj_y_centroid)
row.names(position_xy) <- row.names(xenium@meta.data)
colnames(position_xy) <- c("Super_SP_1", "Super_SP_2")
xenium[["super_sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "Super_SP_", assay = DefaultAssay(xenium))
DimPlot(xenium, group.by = "sample", label = TRUE, reduction = "super_sp") + NoLegend() + coord_equal()

# Create new CT var
xenium$CT_secondpass3 <- xenium$CT_secondpass
xenium$CT_secondpass3[xenium$CT_firstpass_full %in% c("cDCs", "cDC2 (Air Obj)", "Langerhans DCs")] <- "cDCs"
xenium$CT_secondpass3[xenium$CT_firstpass_full == "Transitional AT2 (Myeloid Obj)"] <- "Langerhans cells"

## SET VARIABLE NAMES ----
# Broad, 25, 8-12 niches
ct_var = xenium$CT_secondpass3
ct_name_var = "CT_secondpass3"
num_neighbors_list = 25
num_niches_list = c(12, 10, 11)


## RUN NICHE ANALYSIS ----
for (j in 1:length(num_neighbors_list)) {
  num_neighbors = num_neighbors_list[j]
  print(paste0("Neighbors K: ", num_neighbors))
  
  for (i in num_niches_list) {
    num_niches = i
    print(paste0("Niches K: ", num_niches))
    
    # Set directory names for saving files
    niche_name_var <-  paste0(ct_name_var, "_nichek", num_niches, "_neighborsk",
                              num_neighbors)
    plot_save_dir <- paste0(scratch_dir, niche_name_var, "/")
    
    print(paste0("Saving plots to ", plot_save_dir))
    
    ## CREATE FOV AND FIND CELL NICHES ----
    # Create FOV
    coords <- xenium@meta.data[, c("super_adj_x_centroid", "super_adj_y_centroid")]
    xenium[["fov"]]  <- CreateFOV(coords, assay = "RNA", type = "centroids")

    Sys.time()
    xenium <- BuildNicheAssay(object = xenium, fov = "fov", group.by = ct_name_var,
                              niches.k = num_niches, neighbors.k = num_neighbors)
    Sys.time()

    # Set more variables
    niche_var <- xenium$niches
    xenium@meta.data[[niche_name_var]] <- xenium$niches
    
    ## GENERAL SUMMARY FIGURES ----
    # Niche by Cell Type
    CairoPNG(paste0(plot_save_dir, niche_name_var, "_by_ct.png"), width = 10, height = 6, units = "in", dpi = 300)
    a <- table(ct_var, niche_var) %>%
      as.data.frame() %>%
      ggplot(aes(x = niche_var, y = Freq, fill = ct_var)) +
      geom_col(position = position_fill()) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = color_list) +
      labs(x = "Cell Niche", y = "Cell Type Proportion", fill = "", title = niche_name_var) +
      theme_bw() +
      theme(axis.text = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5))
    print(a)
    dev.off()
    
    # Cell Type by Niche
    CairoPNG(paste0(plot_save_dir, "ct_by_", niche_name_var, ".png"), width = 10, height = 6, units = "in", dpi = 300)
    a <- table(ct_var, niche_var) %>%
      as.data.frame() %>%
      ggplot(aes(x = ct_var, y = Freq, fill = niche_var)) +
      geom_col(position = position_fill(), width = 0.8) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = nuclei_niche_color_list) +
      labs(y = "Cell Niche Proportion", fill = "", title = niche_name_var) +
      theme_bw(base_size = 12) +
      theme(axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5))
    print(a)
    dev.off()
    
    # Niche by Sample Type
    CairoPNG(paste0(plot_save_dir, niche_name_var, "_by_sample_type.png"), width = 5, height = 5, units = "in", dpi = 300)
    a <- table(xenium$sample_type, niche_var) %>% 
      as.data.frame() %>%
      mutate(Var1 = ordered(Var1, levels = c("Unaffected", "LF", "MF", "INT"))) %>%
      ggplot(aes(x = Var1, y = Freq, fill = niche_var)) +
      geom_col(position = position_fill()) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = nuclei_niche_color_list) +
      theme_classic(base_size = 14) +
      labs(x = "", y = "Niche Proportions", title = niche_name_var) +
      theme(legend.title = element_blank(), 
            axis.text = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5))
    print(a)
    dev.off()
    
    # Niche by Sample
    CairoPNG(paste0(plot_save_dir, niche_name_var, "_by_sample.png"), width = 8, height = 5, units = "in", dpi = 300)
    a <- table(xenium$sample, niche_var) %>%
      as.data.frame() %>%
      # mutate(Var1 = ordered(Var1, levels = c("VUHD116A", "VUHD116B", "VUHD095", "VUHD113", 
      #                                        "VUHD069", "THD0011", "THD0008", "TILD117LF",
      #                                        "VUILD96LF", "VUILD91LF", "VUILD78LF", "VUILD105LF",
      #                                        "VUILD48LF", "VUILD102LF", "VUILD104LF", "VUILD78MF",
      #                                        "VUILD96MF", "VUILD48MF", "VUILD104MF", "TILD117MF", 
      #                                        "VUILD102MF", "VUILD107MF", "VUILD91MF",  "VUILD110", 
      #                                        "VUILD105MF", "TILD175", "VUILD115", "VUILD106"))) %>%
      ggplot(aes(x = Var1, y = Freq, fill = niche_var)) +
      geom_col(position = position_fill()) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = nuclei_niche_color_list) +
      theme_classic(base_size = 14) +
      labs(x = "", y = "Niche Proportions", title = niche_name_var) +
      theme(legend.title = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
            plot.title = element_text(hjust = 0.5))
    print(a)
    dev.off()
    
    ## SPLIT SAMPLE PLOTS ----
    # Split samples
    sample_ids <- unique(xenium$sample)
    sample_list <- lapply(sample_ids, function(XX) { subset(xenium, subset = sample == XX) })
    names(sample_list) <- sample_ids
    
    # Create DimPlots for each sample
    dimplot_list <- lapply(sample_ids, function(XX) {
      DimPlot(sample_list[[XX]], reduction = "sp", 
              group.by = niche_name_var, cols = nuclei_niche_color_list, 
              raster = FALSE) + pretty_umap + 
        ggtitle(paste0(XX, "\n", niche_name_var)) +
        coord_equal()
    })
    names(dimplot_list) <- sample_ids
    
    # Create PNGs for each figure
    for (sm in sample_ids) {
      CairoPNG(paste0(plot_save_dir, sm, "_", niche_name_var, ".png"), width = 7.5, height = 7.5, units = "in", dpi = 300)
      print(dimplot_list[[sm]])
      dev.off()
    }
    # Save object
    saveRDS(xenium, paste0(scratch_dir, "xenium_", niche_name_var, ".rds"))
  }
}

saveRDS(xenium, paste0(scratch_dir, "xenium_", paste0(ct_name_var, "_neighborsk", num_neighbors_list), ".rds"))
