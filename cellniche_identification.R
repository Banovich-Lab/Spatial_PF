############################################
# Cell Niches
# Author: Annika Vannan (avannan@tgen.org)
# Date: 12/18/2023
# Description: Cell niches with Seurat v5
############################################


## SETTING ENVIRONMENT ----
# # Set library paths so they are consistent in RStudio Server and command line R
# .libPaths(c("/home/avannan/R/x86_64-pc-linux-gnu-library/4.2",
#             "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

# .libPaths(c("/home/avannan/R/rstudio-4.3.0-3-with_modules.sif/libs",
#             "/home/avannan/miniconda3/lib/R/library",
#             "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

library(Seurat)
library(SeuratObject)
library(tidyverse)
library(gplots)
library(tibble)


## CREATE COLOR LISTS ----
color_list <- list(`aCap` = "#aa4e26",
                   `Arteriole` = "#d59445",
                   `gCap` = "#cdbf5e",
                   `DCs/Lymphatic` = "#e46300",
                   `Venous` = "#975f47",
                   `Proliferating Endothelial` = "#dacf3c",  
                   
                   `Activated FBs (CTHRC1+/FAP+)` = "#765fa5",        
                   `Activated FBs` = "#765fa5",  
                   `Adventitial FBs` = "#01b0ea",                     
                   `Alveolar FBs` = "#6d4bce",        
                   `Fibroblasts` = "#7cb3f5",                        
                   `HAS1+ Fibroblasts` = "#0065f5",  
                   `Lipofibroblasts` = "#c6adf2",    
                   `Mesothelial` = "#374ca3",    
                   `Peribronchial FBs` = "#9da7ff",
                   `PLIN2+ Fibroblasts` = "#4c4e8c",
                   `SMCs` = "#0192fe",    
                   `Proliferating Mesenchymal` = "#7895ff",
                   
                   `B cells` = "#d84265",
                   `NK cells` = "#b16f6f",
                   `T-cells` = "#e62f4c",  
                   `Plasma` = "#df4738",
                   `pDCs` = "#7b39c7",
                   `Proliferating Lymphoid` = "#993432",
                   
                   `cDCs` =  "#7a4ca1",
                   `DCs` = "#ac67dc",
                   `FABP4+ Macrophages` = "#cc85d9",
                   `Macrophages` = "#dc48bd",
                   `Mast` = "#bf4ce3",
                   `Monocyte-derived Macrophages` = "#963095",
                   `Monocyte-derived Macrophages (FCN1+)`  = "#d497ba",                          
                   `SPP1+ Macrophages` = "#9e3c6d",
                   `Proliferating Myeloid` = "#e25292",  
                   
                   `AT1` = "#6cd66d",                                 
                   `AT2` = "#34491d",                                                        
                   `Basal` = "#acde47",                                                     
                   `Ciliated` = "#b9b700",                          
                   `Differentiating Ciliated` = "#01eaaa",           
                   `KRT5-/KRT17+` = "#a0cd9e",                                        
                   `MUC5B+ Secretory` = "#31a200",      
                   `SCGB3A2+` = "#ccd87a",
                   `SCGB3A2+/SCGB1A1+` = "#006e3e",  
                   `Transitional AT2` = "#778d39",   
                   `Proliferating Epithelial` = "#586551")

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
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/rds_files/xenium_full_CT_051923.rds")

# Set directories for saving plots
scratch_dir <- "/Users/avannan/Documents/nuclei_niches/"


## CREATE COLOR LISTS ----
color_list <- list(`aCap` = "#aa4e26",
                   `Arteriole` = "#d59445",
                   `gCap` = "#cdbf5e",
                   `DCs/Lymphatic` = "#e46300",
                   `Venous` = "#975f47",
                   `Proliferating Endothelial` = "#dacf3c",  
                   
                   `Activated FBs (CTHRC1+/FAP+)` = "#765fa5",        
                   `Activated FBs` = "#765fa5",  
                   `Adventitial FBs` = "#01b0ea",                     
                   `Alveolar FBs` = "#6d4bce",        
                   `Fibroblasts` = "#7cb3f5",                        
                   `HAS1+ Fibroblasts` = "#0065f5",  
                   `Lipofibroblasts` = "#c6adf2",    
                   `Mesothelial` = "#374ca3",    
                   `Peribronchial FBs` = "#9da7ff",
                   `PLIN2+ Fibroblasts` = "#4c4e8c",
                   `SMCs` = "#0192fe",    
                   `Proliferating Mesenchymal` = "#7895ff",
                   
                   `B cells` = "#d84265",
                   `NK cells` = "#b16f6f",
                   `T-cells` = "#e62f4c",  
                   `Plasma` = "#df4738",
                   `pDCs` = "#7b39c7",
                   `Proliferating Lymphoid` = "#993432",
                   
                   `cDCs` =  "#7a4ca1",
                   `DCs` = "#ac67dc",
                   `FABP4+ Macrophages` = "#cc85d9",
                   `Macrophages` = "#dc48bd",
                   `Mast` = "#bf4ce3",
                   `Monocyte-derived Macrophages` = "#963095",
                   `Monocyte-derived Macrophages (FCN1+)`  = "#d497ba",                          
                   `SPP1+ Macrophages` = "#9e3c6d",
                   `Proliferating Myeloid` = "#e25292",  
                   
                   `AT1` = "#6cd66d",                                 
                   `AT2` = "#34491d",                                                        
                   `Basal` = "#acde47",                                                     
                   `Ciliated` = "#b9b700",                          
                   `Differentiating Ciliated` = "#01eaaa",           
                   `KRT5-/KRT17+` = "#a0cd9e",                                        
                   `MUC5B+ Secretory` = "#31a200",      
                   `SCGB3A2+` = "#ccd87a",
                   `SCGB3A2+/SCGB1A1+` = "#006e3e",  
                   `Transitional AT2` = "#778d39",   
                   `Proliferating Epithelial` = "#586551")

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


## FIX COORDINATES ----
# Set spatial rows and columns
xenium$spatial_row <- ""
xenium$spatial_row[xenium$sample %in% c("VUHD116A", "VUILD102LF", "VUHD095", "VUILD105MF", "VUHD069")] <- "R1"
xenium$spatial_row[xenium$sample %in% c("VUILD96MF", "VUILD102MF", "VUHD116B", "VUILD104LF", "VUHD113", "VUILD105LF")] <- "R2"
xenium$spatial_row[xenium$sample %in% c("VUILD96LF", "VUILD107MF", "VUILD48MF", "VUILD104MF", "VUILD48LF")] <- "R3"
xenium$spatial_row[xenium$sample %in% c("VUILD115", "VUILD110", "THD0011", "TILD117MF")] <- "R4"
xenium$spatial_row[xenium$sample %in% c("VUILD78MF", "TILD117LF", "VUILD91MF")] <- "R5"
xenium$spatial_row[xenium$sample %in% c("VUILD106", "THD0008", "VUILD91LF", "VUILD78LF", "TILD175")] <- "R6"
xenium$spatial_column <- ""
xenium$spatial_column[xenium$sample %in% c("VUILD96MF", "VUILD115", "VUILD106")] <- "C1"
xenium$spatial_column[xenium$sample %in% c("VUHD116A", "VUILD102MF", "VUILD96LF")] <- "C2"
xenium$spatial_column[xenium$sample %in% c("VUILD102LF", "VUHD116B", "VUILD107MF", "VUILD110", "THD0008")] <- "C3"
xenium$spatial_column[xenium$sample %in% c("VUHD095", "VUILD104LF", "VUILD48MF", "VUILD78MF", "VUILD91LF")] <- "C4"
xenium$spatial_column[xenium$sample %in% c("VUILD105MF", "VUHD113", "VUILD104MF", "THD0011", "TILD117LF", "VUILD78LF")] <- "C5"
xenium$spatial_column[xenium$sample %in% c("VUHD069", "VUILD105LF", "VUILD48LF", "TILD117MF", "VUILD91MF", "TILD175")] <- "C6"

# Adjust coordinates based on row and column positions
xenium$super_adj_y_centroid <- xenium$adj_y_centroid
xenium$super_adj_y_centroid[xenium$spatial_row == "R1"] <- xenium$adj_y_centroid[xenium$spatial_row == "R1"] + 4400
xenium$super_adj_y_centroid[xenium$spatial_row == "R2"] <- xenium$adj_y_centroid[xenium$spatial_row == "R2"] + 2200
xenium$super_adj_y_centroid[xenium$spatial_row == "R4"] <- xenium$adj_y_centroid[xenium$spatial_row == "R4"] -1000
xenium$super_adj_y_centroid[xenium$spatial_row == "R5"] <- xenium$adj_y_centroid[xenium$spatial_row == "R5"] - 2200
xenium$super_adj_y_centroid[xenium$spatial_row == "R6"] <- xenium$adj_y_centroid[xenium$spatial_row == "R6"] - 4400
xenium$super_adj_x_centroid <- xenium$adj_x_centroid
xenium$super_adj_x_centroid[xenium$spatial_column == "C1"] <- xenium$adj_x_centroid[xenium$spatial_column == "C1"] - 3000
xenium$super_adj_x_centroid[xenium$spatial_column == "C2"] <- xenium$adj_x_centroid[xenium$spatial_column == "C2"] - 2000
xenium$super_adj_x_centroid[xenium$spatial_column == "C3"] <- xenium$adj_x_centroid[xenium$spatial_column == "C3"] - 1000
xenium$super_adj_x_centroid[xenium$spatial_column == "C5"] <- xenium$adj_x_centroid[xenium$spatial_column == "C5"] + 2200
xenium$super_adj_x_centroid[xenium$spatial_column == "C6"] <- xenium$adj_x_centroid[xenium$spatial_column == "C6"] + 4400

# Add "super" spatial dimension reduction object
position_xy <- cbind(xenium$super_adj_x_centroid,
                     xenium$super_adj_y_centroid)
row.names(position_xy) <- row.names(xenium@meta.data)
colnames(position_xy) <- c("Super_SP_1", "Super_SP_2")
xenium[["super_sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "Super_SP_", assay = DefaultAssay(xenium))



# SET VARIABLE NAMES ----
# Broad, 25, 8-12 niches
ct_var = xenium$broad_CT5
ct_name_var = "broad_CT5"
num_neighbors_list = 25
num_niches_list = c(8, 9, 10, 11, 12)


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
      labs(x = "Nuclei Niche", y = "Cell Type Proportion", fill = "", title = niche_name_var) +
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
      labs(y = "Nuclei Niche Proportion", fill = "", title = niche_name_var) +
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
      mutate(Var1 = ordered(Var1, levels = c("Unaffected", "LF", "MF", "ILD"))) %>%
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
    CairoPNG(paste0(plot_save_dir, niche_name_var, "_by_sample.png"), width = 5, height = 5, units = "in", dpi = 300)
    a <- table(xenium$sample, niche_var) %>%
      as.data.frame() %>%
      mutate(Var1 = ordered(Var1, levels = c("VUHD116A", "VUHD116B", "VUHD095", "VUHD113", 
                                             "VUHD069", "THD0011", "THD0008", "TILD117LF",
                                             "VUILD96LF", "VUILD91LF", "VUILD78LF", "VUILD105LF",
                                             "VUILD48LF", "VUILD102LF", "VUILD104LF", "VUILD78MF",
                                             "VUILD96MF", "VUILD48MF", "VUILD104MF", "TILD117MF", 
                                             "VUILD102MF", "VUILD107MF", "VUILD91MF",  "VUILD110", 
                                             "VUILD105MF", "TILD175", "VUILD115", "VUILD106"))) %>%
      ggplot(aes(x = Var1, y = Freq, fill = niche_var)) +
      geom_col(position = position_fill()) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = nuclei_niche_color_list) +
      theme_classic(base_size = 14) +
      labs(x = "", y = "Niche Proportions", title = niche_name_var) +
      theme(legend.title = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
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
              group.by = "niches", cols = nuclei_niche_color_list, 
              raster = FALSE) + pretty_umap + 
        ggtitle(paste0(XX, "\n", niche_name_var))
    })
    names(dimplot_list) <- sample_ids
    
    # Create PNGs for each figure
    for (sm in sample_ids) {
      CairoPNG(paste0(plot_save_dir, sm, "_", niche_name_var, ".png"), width = 7.5, height = 7.5, units = "in", dpi = 300)
      print(dimplot_list[[sm]])
      dev.off()
    }
    # ave object
    saveRDS(xenium, paste0(scratch_dir, "xenium_", niche_name_var, ".rds"))
  }
}

saveRDS(xenium, paste0(scratch_dir, "xenium_", paste0(ct_name_var, "_neighborsk", num_neighbors_list), ".rds"))




