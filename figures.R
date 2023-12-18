############################################
# Figures for Late Spatial PF Project
# Author: Annika Vannan (avannan@tgen.org)
# Date: 12/18/2023
# Description: Code for creating select
#              figures
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


## SET COLORS ----
color_list <- list(`aCap` = "#aa4e26",
                   `Arteriole` = "#d59445",
                   `gCap` = "#cdbf5e",
                   `DCs/Lymphatic` = "#975f47",
                   `Lymphatic/DCs` = "#975f47",
                   `Venous` = "#e46300",
                   `Proliferating Endothelial` = "#ffd741",  
                   
                   `Activated FBs (CTHRC1+/FAP+)` = "#005a90",        
                   `Activated FBs` = "#005a90",      
                   `Adventitial FBs` = "#00cace",                     
                   `Alveolar FBs` = "#6d4bce",        
                   `Fibroblasts` = "#7cb3f5",                        
                   `HAS1+ Fibroblasts` = "#0065f5",  
                   `Lipofibroblasts` = "#c6adf2",    
                   `Mesothelial` = "#374ca3",    
                   `Peribronchial FBs` = "#4c4e8c",
                   `PLIN2+ Fibroblasts` = "#7895ff",
                   `SMCs` = "#00aec9",
                   `Proliferating Mesenchymal` = "#7dfff6",
                   
                   `B cells` = "#d84265",
                   `NK cells` = "#b16f6f",
                   `T-cells` = "#e62f4c",  
                   `Plasma` = "#df4738",
                   `pDCs` = "#7b39c7",
                   `Proliferating Lymphoid` = "#993432",
                   
                   `cDCs` =  "#7a4ca1",
                   `DCs` = "#ac67dc",
                   `FABP4+ Macrophages` = "#ffa3a1",
                   `FABP4+ Mϕs` = "#ffa3a1",
                   `Macrophages` = "#dc48bd",
                   `Mast` = "#bf4ce3",
                   `Interstitial Macrophages` = "#963095",
                   `Interstitial Macrophages (FCN1+)`  = "#f2294e",            
                   `SPP1+ Macrophages` = "#9e3c6d",
                   `SPP1+ Mϕs` = "#9e3c6d",
                   `Proliferating Myeloid` = "#e25292",  
                   
                   `AT1` = "#6cd66d",
                   `AT2` = "#34491d",
                   `Basal` = "#31a200",
                   `Ciliated` = "#b9b700",
                   `Differentiating Ciliated` = "#006e3e",
                   `KRT5-/KRT17+` = "#a0cd9e",
                   `MUC5B+` = "#acde47",
                   `SCGB3A2+` = "#ccd87a",
                   `SCGB3A2+/SCGB1A1+` = "#01eaaa",
                   `Transitional AT2` = "#778d39",
                   `Proliferating Epithelial` = "#5e6b00")

nuclei_niche_color_list <- list(`C1` = "#ffa26d",
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



## LOAD AND ORGANIZE DATA ----
# xenium <- readRDS("/Volumes/dback_scratch/avannan/full_xenium_07-13-23.rds")
# sort(table(xenium$sample))
# 

lumen_nuclei_RNA_sce <- readRDS(file = "/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/lumen_analysis/lumen_analysis_Nov15/RDS/lumen_nuclei_sce_slingshot_Nov15.rds")
lumen_nuclei_RNA_sce_ts_knot14 <- readRDS(file = "/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/lumen_analysis/lumen_analysis_Nov15/RDS/lumen_nuclei_sce_tradeseq_knot14_Nov15.rds")
lumen_cell_types_sce <- readRDS("/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/lumen_analysis/lumen_analysis_Nov15/RDS/lumen_cell_types_filter_sce.rds")
lumen_cell_types_sce_sling <- readRDS("/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/lumen_analysis/lumen_analysis_Nov15/RDS/lumen_cell_types_sce_filter_slingshot_Nov15.rds")
lumen_cell_types_sce_ts <- readRDS("/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/lumen_analysis/lumen_analysis_Nov15/RDS/lumen_cell_types_sce_filter_tradeseq_Nov15.rds")
lumen_tranxniche_sce <- readRDS(file = "/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/lumen_analysis/lumen_analysis_Nov15/RDS/lumen_tranx_niche_sce_filter_15Nov.rds" )
lumen_cellniche_sce_ts <- readRDS(file = "/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/lumen_analysis/lumen_analysis_Nov15/RDS/lumen_cellniche_sce_filter_15Nov.rds" )
annotation_df <- read.csv("/Users/avannan/Downloads/annotation_df_final_Dec112023.csv", skip = 1)


#### ADD ANNOTATION DATA ----
annotation_df <- read.csv("/Users/avannan/Downloads/annotation_df_final_Dec112023.csv", skip = 1)
xenium$Annotation_Type_Instance <- annotation_df$Annotation_Type_Instance[match(colnames(xenium), annotation_df$cell_id)]
xenium$Annotation_Type <- annotation_df$Annotation_Type[match(colnames(xenium), annotation_df$cell_id)]
xenium$Annotation_Instance <- annotation_df$Annotation_Instance[match(colnames(xenium), annotation_df$cell_id)]
xenium$Num_Annotation_Type <- annotation_df$Num_Annotation_Type[match(colnames(xenium), annotation_df$cell_id)]
xenium$CNiche <- paste0("C", xenium$broad_CT5_nichek12_neighborsk25)
xenium$TNiche <- paste0("T", xenium$gmm12_5k_trained_hex_gmm)
xenium@meta.data <- xenium@meta.data %>%
  select(orig.ident, cell_id, x_centroid, y_centroid, adj_x_centroid, adj_y_centroid,
         super_adj_x_centroid, super_adj_y_centroid, nucleus_area, nCount_RNA,
         nFeature_RNA, nCount_cell_RNA, nFeature_cell_RNA, percent.blank, sample,
         patient, sample_type, tma, run, lineage, finest_CT1, fine_CT2, fine_CT3,
         fine_CT4, broad_CT5, broad_CT6, TNiche, CNiche,
         Annotation_Type, Annotation_Type_Instance)
xenium$Annotation_Type[xenium$Annotation_Type == "smaller_airway"] <- NA


#### REMOVE SELECT NUCLEI ----
# Nuclei to remove where appropriate
remove_nuclei_df <- read.csv("/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/remove_nuclei_VUILD105LF.csv") %>%
  dplyr::rename(cell_id = "Cell.ID") %>%
  mutate(cell_id = paste0("VUILD105LF_", cell_id))
remove_nuclei <- remove_nuclei_df$cell_id

# Remove nuclei and create separate object
xenium <- subset(xenium, cells = remove_nuclei, invert = TRUE)

# Save Seurat object
# saveRDS(xenium, "/Volumes/dback_scratch/avannan/full_xenium_12-11-23_rm_nuc.rds")
xenium <- readRDS("/Volumes/dback_scratch/avannan/full_xenium_12-11-23_rm_nuc.rds")

# Reorder cell types
ct_order <- xenium@meta.data %>% 
  select(broad_CT5, lineage) %>% 
  unique() %>% arrange(lineage, broad_CT5) %>% 
  pull(broad_CT5)
new_ct_order <- ct_order[c(5, 6, 1, 2, 4, 3,
                           7:13, 15:17, 14,
                           18, 19, 24, 29, 25, 30, 26,
                           22, 21, 23, 20, 28, 27,
                           31:37, 39, 38)]

# Add pathology score information
library(googlesheets4)
gs4_deauth()
path_sheet <- gs4_get("https://docs.google.com/spreadsheets/d/1Bdcu3R-ptuEkZv_19rX6pcCojOaI6v2aRHac7evkeS0/edit?usp=sharing")
path_sheet <- read_sheet(path_sheet, sheet = 5, skip = 0)
path_sheet <- path_sheet %>%
  filter(Feature == "Adjusted_path_score") %>%
  select(-1) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "sample", values_to = "adjusted_pathology_score")
xenium$adjusted_pathology_score <- path_sheet$adjusted_pathology_score[match(xenium$sample, path_sheet$sample)]


#### RENAME CELLS TO SHORTENED NAMES ---
# Named vector for renaming cells on plots 
rename_cells <- unique(xenium$broad_CT5)
names(rename_cells) <- unique(xenium$broad_CT5)
rename_cells[which(rename_cells == "DCs/Lymphatic")] <- "Lymphatic/DCs"
rename_cells[which(rename_cells %in% c("Proliferating Endothelial",
                                       "Proliferating Epithelial",
                                       "Proliferating Mesenchymal"))] <- "Proliferating"


#### CREATE FACTOR VARIABLES ----
xenium$sample_type <- ordered(xenium$sample_type,
                              levels = c("Unaffected", "LF", "MF", "ILD"))


## FIGURE 1 (UMAP; Part of Biorender Figure) ----
# This panel shows the data processing and was made primarily in Biorender. 
# It also includes the UMAP below.
CairoPNG("/Users/avannan/Documents/spatial_figures/manuscript_submission/xenium_late_UMAP.png", width = 3.5, height = 3.3, units = "in", dpi = 300)
DimPlot(xenium, group.by = "broad_CT5", cols = color_list, raster = FALSE) + NoLegend() + 
  pretty_umap + theme(plot.title = element_blank()) + labs(x = "UMAP 1", y = "UMAP 2")
dev.off()


## FIGURE 2 (Pathology Score and Pathology Annotations) ----
#### FIGURE 2B (Pathology Score x CT Volcano Plots) ----
pathscore_ct_df <- read_csv("/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/adjusted_pathology_scores_2023-11-07/celltype_prop_all_changes_adjpathology_scores.csv")
pathscore_ct_plot <- pathscore_ct_df %>%
  mutate(neglog10.adjP = -log10(adj.P.Val),
         sig = ifelse(adj.P.Val < 0.01, "sig", "n.s."),
         celltype = ifelse(sig == "sig", celltype, "Other"),
         celltype = ifelse(celltype == "Monocyte-derived Macrophages", "Interstitial Macrophages", celltype),
         celltype = ifelse(celltype == "FABP4+ Macrophages", "FABP4+\nMacrophages", celltype),
         celltype = ifelse(celltype == "Proliferating Lymphoid", "Proliferating\nLymphoid", celltype),
         celltype = ifelse(celltype == "Proliferating Endothelial", "Proliferating\nEndothelial", celltype),
         label = ifelse(sig == "sig", celltype, NA_character_)) %>%
  ggplot(aes(x = logFC, y = neglog10.adjP, fill = celltype, color = celltype)) +
  geom_hline(yintercept = -log10(0.01), lty = "longdash", color = "grey40", linewidth = 0.25) +
  geom_vline(xintercept = 0, lty = "longdash", color = "grey40", linewidth = 0.25) +
  geom_point(shape = 21, stroke = 0.2, alpha = 0.75, size = 1.25) +
  geom_text_repel(aes(label = label), fontface = "bold", size = 1.7, seed = 720,
                  min.segment.length = unit(0.2, "lines"),
                  box.padding = 0.1, na.rm = TRUE) +
  theme_classic() +
  theme(axis.line = element_line(color = "black"), 
        axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(color = "black", size = 6),
        text = element_text(family = "Arial")) +
  scale_color_manual(values = c(color_list, `Other` = "grey", `Interstitial Macrophages` = "#963095", `FABP4+\nMacrophages` = "#ffa3a1", 
                               `Proliferating\nEndothelial` = "#ffd741", `Proliferating\nLymphoid` = "#993432")) +
  scale_fill_manual(values = c(color_list, `Other` = "grey", `Interstitial Macrophages` = "#963095", `FABP4+\nMacrophages` = "#ffa3a1", 
                               `Proliferating\nEndothelial` = "#ffd741", `Proliferating\nLymphoid` = "#993432")) + 
  NoLegend() +
  labs(x = "Rate of Change (Pathology Score)", y = "-log10(FDR)")

CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/Fig3_volcano_ct.pdf", width = (50*0.0393701), height = (50*0.0393701))
pathscore_ct_plot
dev.off()


#### FIGURE 2C (Annotation x CT Dotplot) ----
# Create dataframe of lineage variables
lineage_df <- xenium@meta.data %>%
  select(broad_CT5, lineage) %>%
  unique()
names(lineage_df) <- c("CT", "Lineage")

new_annotation_order <- c("artery",
                          "muscularized_artery",
                          "vein",
                          "venule",
                          "airway_smooth_muscle",
                          "interlobular_septum",
                          "fibroblastic_focus",
                          "severe_fibrosis",
                          "normal_alveoli",
                          "minimally_remodeled_alveoli",
                          "hyperplastic_aec",
                          "epithelial_detachment",
                          "multinucleated_cell",
                          "remodeled_epithelium",
                          "small_airway",
                          "microscopic_honeycombing",
                          "goblet_cell_metaplasia",
                          "giant_cell",
                          "granuloma",
                          "mixed_inflammation",
                          "TLS")


# Annotation Composition
path_ct_prop_df <- as.data.frame((table(xenium$Annotation_Type, xenium$broad_CT5)/
                                    rowSums(table(xenium$Annotation_Type, xenium$broad_CT5))))
names(path_ct_prop_df) <- c("Pathology", "CT", "Proportion")

# Only keep some annotations for figure
path_ct_prop_plot <- path_ct_prop_df %>%
  left_join(lineage_df) %>%
  filter(Pathology != "remnant_alveoli") %>% # Remove annotations with only 1 example
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  filter(Pathology != "small_airway") %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   Pathology %in% c("smaller_airway", "small_airway") ~ "Small Airway",
                                   # Pathology == "microscopic_honeycombing" ~ "Micro. Honeycombing",
                                   # Pathology == "minimally_remodeled_alveoli" ~ "Min. Remodeled Alveoli",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Microscopic Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS")),
         CT = ordered(CT, levels = new_ct_order)) %>%
  # filter(Proportion >= 0.01) %>%
  # Only keep some annotations for figure
  filter(Path2 %in% c("Fibroblastic Focus", "Severe Fibrosis", "Normal Alveoli", "Minimally Remodeled Alveoli",
                      "Hyperplastic AECs", "Epithelial Detachment", "Smaller Airway", "Microscopic Honeycombing",
                      "Goblet Cell Metaplasia", "Granuloma", "TLS")) %>%
  droplevels() %>%
  ggplot(aes(x = Path2,
             y = reorder(CT, dplyr::desc(CT)),
             size = Proportion, fill = CT)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6) +
  scale_fill_manual(values = color_list) +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0.01, 1)) +
  scale_y_discrete(labels = rename_cells) +
  theme(axis.text = element_text(color = "black", size = 4.5),
        strip.text = element_text(color = "white", size = 4.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.y = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, size = 6),
        legend.position = c(-0.3, -0.155),
        legend.box.background = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.background = element_blank(),
        plot.margin = margin(5.5, 5.5, 8, 5.5, "pt"),
        legend.margin = margin(1.5, 1.5, 1.5, 1.5, "pt")) +
  labs(title = "Annotation Composition (Cells)") +
  guides(fill = "none", 
         size = guide_legend(title.theme = element_text(size = 5),
                             size = 2,
                              nrow = 6, ncol = 1, title.hjust = 0.5,
                             label.position = "right",
                             label.theme = element_text(size = 3.5))) +
  facet_grid(Lineage~., space = "free", scales = "free") +
  NoLegend()
path_ct_prop_plot

# Change facet colors
path_ct_prop_plot_table <- ggplot_gtable(ggplot_build(path_ct_prop_plot))
striprt <- which(grepl("strip-r", path_ct_prop_plot_table$layout$name) | 
                   grepl("strip-t", path_ct_prop_plot_table$layout$name))
fills <- c("chocolate2", "chartreuse4", "maroon3", "cornflowerblue") # Sample Types
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", path_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  path_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(path_ct_prop_plot_table)

CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/Fig3_anno_ct_comp.pdf", width = (54*0.0393701), height = (92*0.0393701))
grid::grid.draw(path_ct_prop_plot_table)
dev.off()


#### FIGURE 2D (Examples of Annotations and Cell Types) ----
vuild96mf <- subset(xenium, subset = sample == "VUILD96MF" & Annotation_Type == "granuloma")
plt <- vuild96mf@meta.data %>%
  mutate(new_adj_x_centroid = adj_x_centroid/0.2125, 
         new_adj_y_centroid = adj_y_centroid/0.2125) %>%
  ggplot(aes(x = new_adj_x_centroid, y = new_adj_y_centroid, color = broad_CT5)) +
  geom_point(size = 2) +
  theme_classic() +
  coord_equal() +
  scale_color_manual(values = color_list) +
  scale_x_continuous(limits = c(-7500, -5000)) +
  scale_y_continuous(limits = c(-24500, -23000)) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  NoLegend()
pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig2_granuloma_example.pdf", width = (100*0.0393701), height = (100*0.0393701))
plt
dev.off()



vuild78mf <- subset(xenium, subset = sample == "VUILD78MF" & Annotation_Type == "microscopic_honeycombing")
plt <- vuild78mf@meta.data %>%
  mutate(new_adj_x_centroid = adj_x_centroid/0.2125, 
         new_adj_y_centroid = adj_y_centroid/0.2125) %>%
  ggplot(aes(x = new_adj_x_centroid, y = new_adj_y_centroid, color = broad_CT5)) +
  geom_point(size = 1) +
  theme_classic() +
  coord_equal() +
  scale_color_manual(values = color_list) +
  scale_y_continuous(limits = c(-100000, -95000)) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  NoLegend()
plt
pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig2_microhoneycombing_example.pdf", width = (120*0.0393701), height = (120*0.0393701))
plt
dev.off()



vuild105lf <- subset(xenium, subset = sample == "VUILD105LF" & Annotation_Type == "fibroblastic_focus")
plt <- vuild105lf@meta.data %>%
  mutate(new_adj_x_centroid = adj_x_centroid/0.2125, 
         new_adj_y_centroid = adj_y_centroid/0.2125) %>%
  ggplot(aes(x = new_adj_x_centroid, y = new_adj_y_centroid, color = broad_CT5)) +
  geom_point(size = 1.5) +
  theme_classic() +
  coord_equal() +
  scale_color_manual(values = color_list) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  NoLegend()
plt
pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig2_fibroblastic_focus.pdf", width = (80*0.0393701), height = (80*0.0393701))
plt
dev.off()



vuild107mf <- subset(xenium, subset = sample == "VUILD107MF" & Annotation_Type == "hyperplastic_aec")
plt <- vuild107mf@meta.data %>%
  mutate(new_adj_x_centroid = adj_x_centroid/0.2125, 
         new_adj_y_centroid = adj_y_centroid/0.2125) %>%
  ggplot(aes(x = new_adj_x_centroid, y = new_adj_y_centroid, color = broad_CT5)) +
  geom_point(size = 2) +
  theme_classic() +
  coord_equal() +
  scale_color_manual(values = color_list) +
  scale_x_continuous(limits = c(29000, 30000)) +
  scale_y_continuous(limits = c(-42000, -40500)) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  NoLegend()
plt
pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig2_hyperplastic_aecs.pdf", width = (80*0.0393701), height = (80*0.0393701))
plt
dev.off()


vuild104mf <- subset(xenium, subset = sample == "VUILD104MF" & Annotation_Type == "goblet_cell_metaplasia")
plt <- vuild104mf@meta.data %>%
  mutate(new_adj_x_centroid = adj_x_centroid/0.2125, 
         new_adj_y_centroid = adj_y_centroid/0.2125) %>%
  ggplot(aes(x = new_adj_x_centroid, y = new_adj_y_centroid, color = broad_CT5)) +
  geom_point(size = 0.5) +
  theme_classic() +
  coord_equal() +
  scale_color_manual(values = color_list) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  NoLegend()
plt
pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig2_goblet_metaplasia.pdf", width = (80*0.0393701), height = (80*0.0393701))
plt
dev.off()

vuild110 <- subset(xenium, subset = sample == "VUILD110" & Annotation_Type == "TLS")
plt <- vuild110@meta.data %>%
  mutate(new_adj_x_centroid = adj_x_centroid/0.2125, 
         new_adj_y_centroid = adj_y_centroid/0.2125) %>%
  ggplot(aes(x = new_adj_x_centroid, y = new_adj_y_centroid, color = broad_CT5)) +
  geom_point(size = 2) +
  theme_classic() +
  coord_equal() +
  scale_color_manual(values = color_list) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  NoLegend()
plt
pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig2_TLS.pdf", width = (80*0.0393701), height = (80*0.0393701))
plt
dev.off()


## FIGURE 3 (Niches) ----
#### FIGURE 3B (Niche x CT Dotplots - Cell Assignment to Niches) ----
# Transcript Niches
ct_tniche_prop_df <- as.data.frame((table(xenium$broad_CT5, xenium$TNiche)/
                                      rowSums(table(xenium$broad_CT5, xenium$TNiche))))
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
  scale_x_discrete(labels = rename_cells) +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw(base_size = 6) +
  scale_color_manual(values = c("black", "grey80", rep("black", 10)), guide = "none") +
  scale_fill_manual(values = transcript_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 6),
        strip.text = element_text(color = "white")) +
  ggh4x::force_panelsizes(cols = c(2.1, 3.8, 4.1, 3.2)) + 
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
ct_niche_prop_df <- as.data.frame((table(xenium$broad_CT5, xenium$CNiche)/
                                     rowSums(table(xenium$broad_CT5, xenium$CNiche))))
names(ct_niche_prop_df) <- c("CT", "Niche", "Proportion")

ct_niche_prop_plot <- ct_niche_prop_df %>%
  full_join(lineage_df) %>%
  mutate(CT = ordered(CT, levels = new_ct_order),
         Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  ggplot(aes(y = reorder(Niche, dplyr::desc(Niche)), x = CT, 
             size = Proportion, color = as.factor(Niche), 
             alpha = Proportion_Visible,
             fill = as.factor(Niche))) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw(base_size = 6) +
  scale_x_discrete(labels = rename_cells) +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_color_manual(values = c(rep("black", 7), "grey80", rep("black", 4)), guide = "none") +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
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
        plot.title = element_blank()) + 
  ggh4x::force_panelsizes(cols = c(2.1, 3.8, 4.1, 3.2)) + 
  scale_alpha_manual(values = c(0, 1), guide = "none") + 
  NoLegend()
ct_niche_prop_plot

CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/Fig3_cellassign2niches.pdf", width = (110*0.0393701), height = (100*0.0393701))
cowplot::plot_grid(ct_tniche_prop_df_table, ct_niche_prop_plot, ncol = 1, align = "hv")
dev.off()


#### FIGURE 3C (Niche by Sample Type; Cells) ----
sample_type_tranx_niche_plot <- xenium@meta.data %>%
  filter(sample_type != "ILD") %>%
  mutate(TNiche = ordered(TNiche, levels = paste0("T", 1:12))) %>%
  ggplot(aes(x = sample_type, fill = TNiche)) +
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
  filter(sample_type != "ILD") %>%
  mutate(CNiche = ordered(CNiche, levels = paste0("C", 1:12))) %>%
  ggplot(aes(x = sample_type, fill = CNiche)) +
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

CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/Fig3_niche_overall_prop_sample_type.pdf", width = 67*0.0393701, height = 30*0.0393701)
ggarrange(sample_type_tranx_niche_plot, sample_type_cell_niche_plot, ncol = 2)
dev.off()


#### FIGURE 3D (Annotation Composition - Niches) ----
# Transcript Niches (Annotation Composition)
tniche_path_prop_df <- as.data.frame((table(xenium$Annotation_Type, xenium$TNiche)/
                                        rowSums(table(xenium$Annotation_Type, xenium$TNiche))))
names(tniche_path_prop_df) <- c("Pathology", "Niche", "Proportion")
tniche_path_prop_df <- tniche_path_prop_df %>%
  filter(Niche != "TNA", Niche != "NA", !is.na(Niche))
tniche_path_prop_plot <- tniche_path_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("T", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  filter(Pathology != "remnant_alveoli") %>% # Remove annotations with only 1 example
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   Pathology %in% c("small_airway", "smaller_airway") ~ "Small Airway",
                                   Pathology == "microscopic_honeycombing" ~ "Micro. Honeycombing",
                                   Pathology == "minimally_remodeled_alveoli" ~ "Min. Remodeled Alveoli",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Min. Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Micro. Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS"))) %>%
  # Only keep some annotations for figure
  filter(Path2 %in% c("Fibroblastic Focus", "Severe Fibrosis", "Normal Alveoli", "Min. Remodeled Alveoli",
                      "Hyperplastic AECs", "Epithelial Detachment", "Small Airway", "Micro. Honeycombing",
                      "Goblet Cell Metaplasia", "Granuloma", "TLS")) %>%
  ggplot(aes(x = Niche, y = reorder(Path2, dplyr::desc(Path2)), 
             size = Proportion, alpha = Proportion_Visible, fill = Niche)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6) +
  scale_fill_manual(values = transcript_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(), 
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(5.5, 2, 5.5, 5.5, "pt")) + 
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  NoLegend()
tniche_path_prop_plot

# Cell Niches (Annotation Composition)
cniche_path_prop_df <- as.data.frame((table(xenium$Annotation_Type, xenium$CNiche)/
                                        rowSums(table(xenium$Annotation_Type, xenium$CNiche))))
names(cniche_path_prop_df) <- c("Pathology", "Niche", "Proportion")
cniche_path_prop_plot_with_legend <- cniche_path_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  filter(Pathology != "remnant_alveoli") %>% # Remove annotations with only 1 example
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   Pathology %in% c("small_airway", "smaller_airway") ~ "Small Airway",
                                   Pathology == "microscopic_honeycombing" ~ "Micro. Honeycombing",
                                   Pathology == "minimally_remodeled_alveoli" ~ "Min. Remodeled Alveoli",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Min. Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Micro. Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS"))) %>%
  # Only keep some annotations for figure
  filter(Path2 %in% c("Fibroblastic Focus", "Severe Fibrosis", "Normal Alveoli", "Min. Remodeled Alveoli",
                      "Hyperplastic AECs", "Epithelial Detachment", "Small Airway", "Micro. Honeycombing",
                      "Goblet Cell Metaplasia", "Granuloma", "TLS")) %>%
  ggplot(aes(x = Niche, y = reorder(Path2, dplyr::desc(Path2)), 
             size = Proportion, alpha = Proportion_Visible, fill = Niche)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 5) +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
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
        legend.position = "bottom", legend.title.align = 0.5
        # legend.key.width = unit(3.5, "mm"),
        # legend.key.height = unit(2, "mm")
        ) +
  guides(size = guide_legend(title.position = "top", label.position = "bottom")) +
  scale_alpha_manual(values = c(0, 1), guide = "none")
cniche_path_prop_plot <- cniche_path_prop_plot_with_legend + NoLegend()
cniche_path_prop_plot_with_legend


# Combine and save
CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/Fig3_path_niche_dotplot.pdf", width = (100*0.0393701), height = (70*0.0393701))
(tniche_path_prop_plot + cniche_path_prop_plot) + 
  patchwork::plot_layout(widths = unit(c(21.5, 21.5), "mm"), 
                         heights = unit(c(32, 32), "mm"), 
                         guides = "keep", ncol = 2, nrow = 1)
dev.off()


CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/Fig3_dotplot_LEGENDS.pdf", width = (100*0.0393701), height = (70*0.0393701))
cniche_path_prop_plot_with_legend
dev.off()



# FIGURE 4 (KRT5-/KRT17+ Epithelial Detachment) ----
#### FIGURE 4H (Niche x Sample Type Boxplots - T3/C7) ----
# T3 niche
tniche_prop_df <- (table(xenium$sample, xenium$TNiche)/rowSums(table(xenium$sample, xenium$TNiche))) %>%
  as.data.frame() %>%
  rename(sample = "Var1", TNiche = "Var2", Prop = "Freq") %>%
  full_join(xenium@meta.data %>% select(sample, sample_type) %>% unique()) %>%
  mutate(sample_type = ordered(sample_type, levels = c("Unaffected", "LF", "MF", "ILD")))

t3_boxplot <-tniche_prop_df %>%
  filter(TNiche == "T3") %>%
  ggplot(aes(x = sample_type, y = Prop, fill = sample_type)) +
  geom_boxplot(width = 0.8, position = position_dodge(width = 1),
               outlier.shape = 21, color = "black", size = 0.1, 
               outlier.size = 0.1) +
  geom_point(size = 0.1) +
  scale_fill_manual(values = c(list(Unaffected = "white", LF = "grey85", MF = "grey50", ILD = "grey30"))) +
  labs(y = "Proportion of Cells") +
  theme_bw(base_size = 6) +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(color = "white", size = 5.5),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "white", size = 6.5),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey90"),
        plot.title = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        strip.background = element_rect(fill = "#8F7C00"),
        axis.title.y = element_text(size = 6),
        plot.margin = margin(0,0,0.5,0, "pt"),
        legend.title = element_blank()) +
  facet_wrap(~TNiche) +
  NoLegend()
t3_boxplot

# C7 niche
cniche_prop_df <- (table(xenium$sample, xenium$CNiche)/rowSums(table(xenium$sample, xenium$CNiche))) %>%
  as.data.frame() %>%
  rename(sample = "Var1", CNiche = "Var2", Prop = "Freq") %>%
  full_join(xenium@meta.data %>% select(sample, sample_type) %>% unique()) %>%
  mutate(sample_type = ordered(sample_type, levels = c("Unaffected", "LF", "MF", "ILD")))

c7_boxplot <- cniche_prop_df %>%
  filter(CNiche == "C7") %>%
  ggplot(aes(x = sample_type, y = Prop, fill = sample_type)) +
  geom_boxplot(width = 0.8, position = position_dodge(width = 1),
               outlier.shape = 21, color = "black", size = 0.1, 
               outlier.size = 0.1) +
  geom_point(size = 0.1) +
  scale_fill_manual(values = c(list(Unaffected = "white", LF = "grey85", MF = "grey50", ILD = "grey30"))) +
  labs(y = "Proportion of Cells") +
  theme_bw(base_size = 6) +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 5.5),
        axis.title.x = element_blank(),
        strip.text = element_text(color = "white", size = 6.5),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "grey65"),
        panel.grid.minor.y = element_line(color = "grey90"),
        plot.title = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        strip.background = element_rect(fill = "#d10040"),
        axis.title.y = element_text(size = 6),
        plot.margin = margin(0.5,0,0,0, "pt"),
        legend.title = element_blank()) +
  facet_wrap(~CNiche) +
  NoLegend()
c7_boxplot


plt <- (t3_boxplot + c7_boxplot) + patchwork::plot_layout(ncol = 1)
ggsave(plot = plt, "/Users/avannan/Documents/spatial_figures/manuscript/Fig4_t3_c7_boxplots.pdf", width = 53*0.0393701, height = 60*0.0393701, device = "pdf", units = "in")
dev.off()


# FIGURE 5 (Macrophage Accumulation) ----
#### FIGURE 5B (Niche x Sample Type Boxplots - C6) ----
# C6 niche
cniche_prop_df <- (table(xenium$sample, xenium$CNiche)/rowSums(table(xenium$sample, xenium$CNiche))) %>%
  as.data.frame() %>%
  rename(sample = "Var1", CNiche = "Var2", Prop = "Freq") %>%
  full_join(xenium@meta.data %>% select(sample, sample_type) %>% unique()) %>%
  mutate(sample_type = ordered(sample_type, levels = c("Unaffected", "LF", "MF", "ILD")))

c6_boxplot <- cniche_prop_df %>%
  filter(CNiche == "C6") %>%
  ggplot(aes(x = sample_type, y = Prop, fill = sample_type)) +
  geom_boxplot(width = 0.8, position = position_dodge(width = 1),
               outlier.shape = 21, color = "black", size = 0.3, 
               outlier.size = 0.1) +
  geom_point(size = 0.1) +
  scale_fill_manual(values = c(list(Unaffected = "white", LF = "grey85", MF = "grey50", ILD = "grey30"))) +
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
c6_boxplot

ggsave(plot = c6_boxplot, "/Users/avannan/Documents/spatial_figures/manuscript/Fig5_c6_boxplot.pdf", 
       width = 36*0.0393701, height = 34*0.0393701, device = "pdf", units = "in")
dev.off()


# FIGURE 6 (AIRSPACE ANALYSIS) ----
#### FIGURE 6B (PCA and Pseudotime Plots) ----
pca_data <- lumen_cell_types_sce_sling@int_colData@listData$reducedDims[["broad_CT5_prop_PCA_t3"]][, 1:2] %>%
  as.data.frame() %>%
  rownames_to_column(var = "lumen_id") %>%
  mutate(sample_type = ordered(case_when(grepl("LF", lumen_id) == TRUE ~ "LF",
                                         grepl("MF", lumen_id) == TRUE ~ "MF",
                                         grepl("HD", lumen_id) == TRUE ~ "Unaffected",
                                         TRUE ~ "ILD"),
                               levels = rev(c("Unaffected", "LF", "MF", "ILD"))))
curves <- slingCurves(lumen_cell_types_sce_sling, as.df = TRUE)
pseudotime_df <- lumen_cell_types_sce_sling$slingPseudotime_1 %>% as.data.frame() %>% dplyr::rename("pseudotime" = 1)
rownames(pseudotime_df) <- colnames(lumen_cell_types_sce_sling)
pseudotime_df <- pseudotime_df %>%
  rownames_to_column(var = "lumen_id") %>%
  full_join(pca_data)

a <- ggplot(data = pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = sample_type, color = sample_type), size = 0.25) +
  theme_classic() +
  scale_color_manual(values = c(list(Unaffected = "#88CCEE", LF = "pink", MF ="#ac3e6f", ILD = "#d76160"))) +
  scale_shape_manual(values = c(16, 17, 15, 1)) +
  theme(axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.205),
        legend.key.size = unit(2, "mm"),
        legend.background = element_blank()) +
  geom_path(data = curves %>% arrange(Order), aes(x = PC1, y = PC2, group = Lineage), size = 0.5) +
  guides(fill = guide_legend(title = NA)) +
  coord_equal()
a

b <- ggplot(data = pseudotime_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = pseudotime), size = 0.01) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(color = "black", size = 6),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title = element_text(size = 5.5),
        legend.text = element_text(size = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.86, 0.15),
        legend.key.height = unit(1, "mm"),
        legend.key.width = unit(2.7, "mm"),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100),
                        breaks = c(0, 0.625, 1.25), labels = c(0, 0.625, 1.25)) +
  geom_path(data = curves %>% arrange(Order), aes(x = PC1, y = PC2, group = Lineage), size = 0.5) +
  guides(color = guide_colorbar(title = "Pseudotime", title.position = "top", title.hjust = 0.5)) +
  coord_equal()
b

pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig6b_pca_pseudotime.pdf", width = 100*0.0393701, height = 46*0.0393701)
a+b
dev.off()


#### FIGURE 6C (Heatmap of Lumens - Partial) ----
library(data.table)

# Load significant genes
genes_lumen_ordered_df <- read_csv("/Users/avannan/Downloads/lumen_gene_order_12122023.csv")
de_genes_knot14 <- genes_lumen_ordered_df$gene_order_smooth
length(de_genes_knot14)

# Predict expression and normalize
predicted_gene_expr_lumen <- data.frame(predictCells(models = lumen_nuclei_RNA_sce_ts_knot14,
                                                     gene = de_genes_knot14))

predicted_gene_expr_lumen_size_normed <- apply(predicted_gene_expr_lumen,1, function(x){
  log(x)-log(lumen_nuclei_RNA_sce[,colnames(predicted_gene_expr_lumen)]$lumen_size)
})
predicted_gene_expr_lumen_size_normed_zscore <- scale(predicted_gene_expr_lumen_size_normed)

# Get order of genes
predicted_gene_expr_lumen_size_normed_zscore_time_order <- 
  predicted_gene_expr_lumen_size_normed_zscore[order(lumen_nuclei_RNA_sce_ts_knot14$crv$pseudotime),]
gene_order <- de_genes_knot14
set.seed(100)

test <- lumen_nuclei_RNA_sce_ts_knot14$crv$pseudotime[rownames(predicted_gene_expr_lumen_size_normed_zscore_time_order)]
test2 <- lumen_nuclei_RNA_sce_ts_knot14$crv$pseudotime[order(lumen_nuclei_RNA_sce_ts_knot14$crv$pseudotime)]

# Set up data for heatmap
gene_by_lumen <- t(predicted_gene_expr_lumen_size_normed_zscore_time_order)
all_lumens_in_order <- names(test2)
lumens_to_plot <- all_lumens_in_order

# Other info for heatmap
ct_mat <- lumen_cell_types_sce@assays@data@listData[["broad_CT5_prop"]][, lumens_to_plot]
cniche_mat <- lumen_cellniche_sce_ts@assays@data@listData[["bCT5_nichek12_prop"]][, lumens_to_plot]
tniche_mat <- lumen_tranxniche_sce@assays@data@listData[["tranx_nichek12_prop"]][, lumens_to_plot]
rownames(tniche_mat) <- as.numeric(rownames(tniche_mat)) + 1

# Possible labels to add to heatmap
gmm_label = lumen_nuclei_RNA_sce[,rownames(predicted_gene_expr_lumen_size_normed_zscore_time_order)]$GMM_broad_CT5_prop_PCA
sample_type_label = lumen_nuclei_RNA_sce[,rownames(predicted_gene_expr_lumen_size_normed_zscore_time_order)]$sample_type
pseudotime_value = lumen_nuclei_RNA_sce_ts_knot14$crv$pseudotime[order(lumen_nuclei_RNA_sce_ts_knot14$crv$pseudotime)]
path_scores = lumen_nuclei_RNA_sce[,rownames(predicted_gene_expr_lumen_size_normed_zscore_time_order)]$path_scores

# Separate unimodal and bimodal genes
bimodal <- genes_lumen_ordered_df %>% filter(gene_mode == "bimodal") %>% pull(gene_order_smooth)
unimodal <- genes_lumen_ordered_df %>% filter(gene_mode == "unimodal") %>% pull(gene_order_smooth)

# Set up final dataframe and get trajectories
df <- gene_by_lumen[unimodal, lumens_to_plot]
broad_trajectories <- genes_lumen_ordered_df %>% filter(gene_mode == "unimodal") %>% pull(trajectory)
broad_trajectories[broad_trajectories == "Normal"] <- " Normal "

# Get genes in each set
normal_genes <- genes_lumen_ordered_df %>% filter(gene_order_smooth %in% unimodal, gene_order <= 84) %>% pull(gene_order_smooth)
late_genes <- genes_lumen_ordered_df %>% filter(gene_order_smooth %in% unimodal, gene_order >= 180) %>% pull(gene_order_smooth)
early_trans_genes <- genes_lumen_ordered_df %>% filter(gene_order_smooth %in% unimodal, !(gene_order_smooth %in% c(normal_genes, late_genes))) %>% pull(gene_order_smooth)

pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig6c_heatmap_pseudotime_lumens_partial.pdf", width = 60*0.0393701, height = 120*0.0393701)
ht <- Heatmap(df,
              cluster_columns = FALSE, cluster_rows = FALSE,
              row_names_gp = gpar(fontsize = 1),
              column_title_gp = gpar(fontsize = 0),
              show_column_names = FALSE,
              border = TRUE,
              height = unit(57, "mm"),
              width = unit(31, "mm"),
              row_split = broad_trajectories,
              row_title_gp = gpar(fontsize = 5),
              row_title_side = "left",
              
              heatmap_legend_param = list(title = "Z-Score",
                                          legend_direction = "horizontal",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          labels_side = "left",
                                          grid_width = unit(1, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),
              
              left_annotation = rowAnnotation(
                Trajectory = broad_trajectories,
                annotation_name_gp = gpar(fontsize = 0),
                simple_anno_size = unit(0.2, "cm"),
                show_legend = FALSE,
                
                col = list(Trajectory = c(" Normal " = "#04846E", 
                                          "Early Transition" = "#FFC107", 
                                          "Late Remodeling" = "#D81B60")),
                na_col = "white"
              ),
              
              top_annotation =
                HeatmapAnnotation(
                  # Cell types
                  `Interstitial MPs` = ct_mat["Monocyte-derived Macrophages", ],
                  aCap = ct_mat["aCap", ],
                  AT1 = ct_mat["AT1", ],
                  AT2 = ct_mat["AT2", ],
                  `Activated FBs` = ct_mat["Activated FBs", ],
                  `Transitional AT2` = ct_mat["Transitional AT2", ],
                  `FABP4+ MPs` = ct_mat["FABP4+ Macrophages", ],
                  `SPP1+ MPs` = ct_mat["SPP1+ Macrophages", ],
                  `Prolif. Myeloid` = ct_mat["Proliferating Myeloid", ],
                  
                  
                  # CNiche
                  C12 = cniche_mat["12", ],
                  C10 = cniche_mat["10", ],
                  C7 = cniche_mat["7", ],
                  C6 = cniche_mat["6", ],
                  
                  # TNiche
                  T4 = tniche_mat["4", ],
                  T3 = tniche_mat["3", ],
                  T6 = tniche_mat["6", ],
                  T8 = tniche_mat["8", ],
                  
                  # Pseudotime
                  Pseudotime = pseudotime_value,
                  
                  simple_anno_size = unit(0.2, "cm"),
                  simple_anno_size_adjust = TRUE,
                  annotation_name_gp = gpar(fontsize = 5),
                  
                  col = list(
                    `Interstitial MPs` = colorRamp2(c(0, max(ct_mat["Monocyte-derived Macrophages", ])), c("white", "deeppink3")),
                    aCap = colorRamp2(c(0, max(ct_mat["aCap", ])), c("white", "chocolate2")),
                    AT1 = colorRamp2(c(0, max(ct_mat["AT1", ])), c("white", "chartreuse3")),
                    AT2 = colorRamp2(c(0, max(ct_mat["AT2", ])), c("white", "chartreuse3")),
                    `Activated FBs` = colorRamp2(c(0, max(ct_mat["Activated FBs", ])), c("white", "cornflowerblue")),
                    `Transitional AT2` = colorRamp2(c(0, max(ct_mat["Activated FBs", ])), c("white", "chartreuse3")),
                    `FABP4+ MPs` = colorRamp2(c(0, max(ct_mat["FABP4+ Macrophages", ])), c("white", "deeppink3")),
                    `SPP1+ MPs` = colorRamp2(c(0, max(ct_mat["SPP1+ Macrophages", ])), c("white", "deeppink3")),
                    `Prolif. Myeloid` = colorRamp2(c(0, max(ct_mat["Proliferating Myeloid", ])), c("white", "deeppink3")),
                    
                    C12 = colorRamp2(c(0, max(cniche_mat["12", ])), c("white", "#924373")),
                    C10 = colorRamp2(c(0, max(cniche_mat["10", ])), c("white", "#6f774b")),
                    C7 = colorRamp2(c(0, max(cniche_mat["7", ])), c("white", "#d10040")),
                    C6 = colorRamp2(c(0, max(cniche_mat["6", ])), c("white", "#72c7ff")),
                    
                    T4 = colorRamp2(c(0, max(tniche_mat["4", ])), c("white", "#FFCC99")),
                    T3 = colorRamp2(c(0, max(tniche_mat["3", ])), c("white", "#8F7C00")),
                    T6 = colorRamp2(c(0, max(tniche_mat["5", ])), c("white", "#993F00")),
                    T8 = colorRamp2(c(0, max(tniche_mat["8", ])), c("white", "#0075DC")),
                    
                    Pseudotime = colorRamp2(c(min(pseudotime_value), max(pseudotime_value)), c("white", "black"))),
                  
                  annotation_name_side = "left", show_legend = FALSE, gap = unit(0.01, "mm")))
draw(ht, heatmap_legend_side = "bottom", align_heatmap_legend = "heatmap_center")
dev.off()


#### FIGURE 6D (Heatmap of Early Transition Expression) ----
# Load in significant terms x CT and filter to only the middle tranche of genes
sig_terms_perCT <- read_csv("/Users/avannan/Downloads/sig_genes_perCT/perCT/all_perCT_sig_terms_s15.csv") %>%
  full_join(genes_lumen_ordered_df, by = c("gene" = "gene_order_smooth")) %>%
  mutate(total_sig = length(gene)) %>%
  group_by(celltype) %>%
  mutate(prop_sig_ct = length(gene)/total_sig) %>%
  ungroup()

# Load in lumen metadata
lumen_metadata <- readRDS("/Volumes/dback_scratch/avannan/lumen_pca/corrected_lumen_data.rds") %>%
  filter(lumen_id %in% colnames(lumen_nuclei_RNA_sce))

# Subset object to contain only cells found in final lumens
xenium_lumens_cells_only <- subset(xenium, cells = rownames(lumen_metadata))

# Pick out genes for heatmap
heatmap_genes <- early_trans_genes

# Log-normalizing and scaling all features in the RNA assay
# Scaling so that all features can be visualized using the same color scale
xenium_lumens_cells_only <- ScaleData(xenium_lumens_cells_only)

# Create dotplot base
p <- DotPlot(xenium_lumens_cells_only, features = heatmap_genes, 
             group.by = "broad_CT5", scale = TRUE) + 
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
  select(gene, celltype) %>%
  unique() %>%
  filter(gene %in% early_trans_genes, !(gene %in% bimodal)) %>%
  dplyr::rename(broad_CT5 = "celltype") %>%
  mutate(broad_CT5 = case_when(broad_CT5 == "T.cells" ~ "T-cells",
                               broad_CT5 == "KRT5..KRT17." ~ "KRT5-/KRT17+",
                               broad_CT5 == "FABP4..Macrophages" ~ "FABP4+ Macrophages",
                               broad_CT5 == "SPP1..Macrophages" ~ "SPP1+ Macrophages",
                               broad_CT5 == "Monocyte.derived.Macrophages" ~ "Interstitial Macrophages",
                               broad_CT5 == "Transitional.AT2" ~ "Transitional AT2",
                               broad_CT5 == "Alveolar.FBs" ~ "Alveolar FBs",
                               TRUE ~ broad_CT5)) %>%
  full_join(xenium@meta.data %>% select(broad_CT5, lineage) %>% unique()) %>% 
  mutate(broad_CT5 = case_when(broad_CT5 == "Proliferating Lymphoid" ~ "Prolif. Lymphoid",
                               broad_CT5 == "Proliferating Myeloid" ~ "Prolif. Myeloid",
                               broad_CT5 == "Proliferating Mesenchymal" ~ "Prolif. Mesenchymal",
                               broad_CT5 == "Proliferating Epithelial" ~ "Prolif. Epithelial",
                               broad_CT5 == "Proliferating Endothelial" ~ "Prolif. Endothelial",
                               TRUE ~ broad_CT5),
         lineage = case_when(broad_CT5 %in% c("Mast", "FABP4+ Macrophages", "Macrophages", "SPP1+ Macrophages",
                                              "Prolif. Myeloid", "Interstitial Macrophages") ~ "Myeloid",
                             lineage == "Immune" ~ "Lymphoid",
                             TRUE ~ lineage)) %>%
  mutate(sig = 1)

# Get proportion/percentage of significant tests
prop_gene_lumen_lineage_df <- gene_lineage_lumen_heatmap_df_og %>% 
  filter(!is.na(broad_CT5)) %>%
  mutate(num_sig_tests = sum(sig)) %>%
  group_by(lineage) %>%
  mutate(num_sig_lineage = sum(sig),
         prop_sig_lineage = sum(sig)/num_sig_tests) %>%
  select(lineage, prop_sig_lineage, num_sig_lineage, num_sig_tests) %>%
  unique() %>%
  mutate(title = paste0(num_sig_lineage, "/", num_sig_tests, " (", round(prop_sig_lineage*100, 1), "%)"))

# Prep dataframe for heatmap
gene_lineage_lumen_heatmap_df <- gene_lineage_lumen_heatmap_df_og %>%
  select(-broad_CT5) %>%
  unique() %>% 
  pivot_wider(names_from = "lineage", values_from = "sig") %>%
  select(-`NA`) %>%
  filter(!is.na(gene)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))
gene_lineage_lumen_heatmap_df2 <- gene_lineage_lumen_heatmap_df %>%
  select(Endothelial, Epithelial, Lymphoid, Myeloid, Mesenchymal) %>%
  as.data.frame()
rownames(gene_lineage_lumen_heatmap_df2) <- gene_lineage_lumen_heatmap_df$gene

# Put genes in order
genes_in_order <- genes_lumen_ordered_df %>%
  filter(gene_mode %in% "unimodal", trajectory == "Early Transition") %>% 
  pull(gene_order_smooth)
gene_lineage_lumen_heatmap_df2 <- gene_lineage_lumen_heatmap_df2[genes_in_order, ]

# Fix labels
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Lymphoid"] <- "Prolif. Lymphoid"
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Myeloid"] <- "Prolif. Myeloid"
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Mesenchymal"] <- "Prolif. Mesenchymal"
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Endothelial"] <- "Prolif. Endothelial"
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Epithelial"] <- "Prolif. Epithelial"

new_ct_order2 <- new_ct_order
new_ct_order2[new_ct_order2 == "Proliferating Lymphoid"] <- "Prolif. Lymphoid"
new_ct_order2[new_ct_order2 == "Proliferating Myeloid"] <- "Prolif. Myeloid"
new_ct_order2[new_ct_order2 == "Proliferating Mesenchymal"] <- "Prolif. Mesenchymal"
new_ct_order2[new_ct_order2 == "Proliferating Endothelial"] <- "Prolif. Endothelial"
new_ct_order2[new_ct_order2 == "Proliferating Epithelial"] <- "Prolif. Epithelial"

# Fix colors
new_color_list <- color_list[new_ct_order]
names(new_color_list) <- new_ct_order2

# Make heatmap
hp <- Heatmap(exp_mat[, new_ct_order2], # t(exp_mat)
              name = "Scaled Expression",
              heatmap_legend_param = list(title = "Scaled Expression",
                                          legend_direction = "horizontal",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          grid_width = unit(1, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),
              column_title_side = "top",
              row_title = " ",
              column_title_gp = gpar(fontsize = 5, col = c("chocolate2", "chartreuse4",
                                                           "darkmagenta", "deeppink3", "cornflowerblue")),
              column_names_side = "top",
              show_column_names = FALSE,
              column_names_rot = 90,
              row_names_side = "right",
              column_split = c(rep(" Endothelial  ", 6), rep(" Epithelial ", 11),
                               rep(" Lymphoid ", 7), rep(" Myeloid ", 6), rep("Mesenchymal", 9)),
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize = 5),
              cluster_rows = FALSE,  cluster_columns = FALSE,
              column_order = new_ct_order2,
              border = "black",
              top_annotation = 
                columnAnnotation(`Cell Type` = colnames(exp_mat[, new_ct_order2]),
                                 col = list(`Cell Type` = unlist(new_color_list)),
                                 show_legend = FALSE,
                                 annotation_name_gp = gpar(fontsize = 0),
                                 simple_anno_size = unit(1, "mm")),
              left_annotation =
                rowAnnotation(`13.2% -`= anno_points(gene_lineage_lumen_heatmap_df2$Endothelial, which = "row",
                                                   ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                   pch = 15,
                                                   axis = FALSE, border = FALSE, gp = gpar(col = "chocolate2"),
                                                   size = unit(1, "mm")),
                              `36.8% -` = anno_points(gene_lineage_lumen_heatmap_df2$Epithelial, which = "row",
                                                      ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                      pch = 15,
                                                      axis = FALSE, border = FALSE, gp = gpar(col = "chartreuse4"),
                                                      size = unit(1, "mm")),
                              `13.2% -` = anno_points(gene_lineage_lumen_heatmap_df2$Lymphoid, which = "row",
                                                      ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                      pch = 15,
                                                      axis = FALSE, border = FALSE, gp = gpar(col = "darkmagenta"),
                                                      size = unit(1, "mm")),
                              `18.4% -` = anno_points(gene_lineage_lumen_heatmap_df2$Myeloid, which = "row",
                                                      ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                      pch = 15,
                                                      axis = FALSE, border = FALSE, gp = gpar(col = "deeppink3"),
                                                      size = unit(1, "mm")),
                              ` 18.4% -` = anno_points(gene_lineage_lumen_heatmap_df2$Mesenchymal, which = "row",
                                                       ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                       pch = 15,
                                                       axis = FALSE, border = FALSE, gp = gpar(col = "cornflowerblue"),
                                                       size = unit(1, "mm")),
                              annotation_name_gp = gpar(fontsize = 5),
                              # show_annotation_name = FALSE,
                              annotation_name_rot = 90,
                              annotation_name_side = "bottom",
                              gap = unit(0.5, "mm")
                ),
              gap = unit(10, "mm")
)

# Add annotations to heatmap
ht <- draw(hp, 
           heatmap_legend_side = "bottom",
           align_heatmap_legend = "heatmap_center")
ht

pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig6D_heatmap_early_transition.pdf", width = (85*0.0393701), height = (100*0.0393701))
ht
dev.off()


#### FIGURE 6E (Legend for Example Lumens) ----
example_lumens <- lumens_to_plot
example_lumens[example_lumens %in% c("VUILD115_L1184", "VUILD102MF_L463")] <- "KEEP"
example_lumens[example_lumens != "KEEP"] <- NA

pdf("/Users/avannan/Documents/spatial_figures/manuscript/Fig6_example_lumen_legend.pdf",
    width = (50*0.0393701), height = (20*0.0393701))
plot(HeatmapAnnotation(
  lumens = example_lumens,
  Pseudotime = pseudotime_value,
  `FABP4+ Mϕs` = ct_mat["FABP4+ Macrophages", ],
  `SPP1+ Mϕs` = ct_mat["SPP1+ Macrophages", ],
  simple_anno_size = unit(0.2, "cm"),
  simple_anno_size_adjust = TRUE,
  annotation_name_gp = gpar(fontsize = 5),
  col = list(Pseudotime = colorRamp2(c(min(pseudotime_value), max(pseudotime_value)), c("white", "black")),
             `FABP4+ Mϕs` = colorRamp2(c(0, max(ct_mat["FABP4+ Macrophages", ])), c("white", "#FF8C37")),
             `SPP1+ Mϕs` = colorRamp2(c(0, max(ct_mat["SPP1+ Macrophages", ])), c("white", "#3DB5DF"))),
  annotation_name_side = "left", show_legend = FALSE, gap = unit(0.01, "mm"), na_col = "white", width = unit(35, "mm")))
dev.off()




## SUPPLEMENTARY FIGURE 1 (Demographic Information) ----
# Demographic information
library(googlesheets4)
gs4_deauth()
ipf_sheet <- gs4_get("https://docs.google.com/spreadsheets/d/1Bdcu3R-ptuEkZv_19rX6pcCojOaI6v2aRHac7evkeS0/edit?usp=sharing")
ipf_sheet <- read_sheet(ipf_sheet, sheet = 1, skip = 1)

nSamples <- 28
nDonors <- nrow(ipf_sheet)
demo_plot <- ipf_sheet %>%
  select(-Age) %>%
  replace(is.na(.), "N/A") %>%
  mutate(`Clinical Diagnosis` = case_when(Clinical_Diagnosis == "IPF" ~ "IPF",
                                          Clinical_Diagnosis == "Control" ~ "Unaffected",
                                          TRUE ~ "Other ILD"),
         Sample_Type = case_when(Sample_Type == "Duplicate_control" ~ "Rep. Unaffected",
                                 Sample_Type == "More_and_less_fibrotic" ~ "Paired LF & MF",
                                 Sample_Type %in% c("ILD", "More_fibrotic_only", "Control") ~ "No Replicates")) %>%
  select(-Clinical_Diagnosis) %>%
  pivot_longer(cols = c("Gender", "Ethnicity", "Tobacco", "Clinical Diagnosis", "Sample_Type"),
               names_to = "Demographic_Feature", values_to = "Value") %>%
  group_by(Demographic_Feature, Value) %>%
  mutate(percentage = nrow(Value)) %>%
  summarize(prop = length(Value)/nDonors) %>% 
  mutate(Demographic_Feature = ordered(Demographic_Feature, levels = c("Clinical Diagnosis", "Sample_Type", "Ethnicity", "Gender", "Tobacco")),
         Value = ordered(Value, levels = c("N/A", "Other ILD", "IPF", "Unaffected", "No Replicates", "Paired LF & MF", "Rep. Unaffected",
                                           "African American", "American Indian", "European", "M", "F",  "N", "Y"))) %>% 
  filter(Demographic_Feature != "Gender") %>%
  ggplot(aes(x = Demographic_Feature, y = prop*100, fill = Value)) +
  geom_col(position = "stack", color = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = c("Diagnosis", "Replicates", "Ethnicity", "Ever Smoker")) +
  scale_fill_manual(values = c("grey80", 
                               "#ac3e6f", "#d76160", "#88CCEE",
                               "grey80", "#C95565", "#88CCEE",
                               "#a265c2", "#009E73", "#F0E442",
                               "#005AB5", "#DC3220")) +
  labs(y = "% of Donors") +
  theme_bw(base_size = 7) + 
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

CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig1_demographic_plot.pdf", width = 90*0.0393701, height = 85*0.0393701)
demo_plot
dev.off()

## SUPPLEMENTARY FIGURE 2 (Cell Types) ----
#### SUPPLEMENTARY FIGURE 2A (Cell Type Proportions Overall) ----
# Cell type proportions
lineage_color_list <- list(Endothelial = "chocolate2",
                           Epithelial = "chartreuse4",
                           Immune = "deeppink3",
                           Mesenchymal = "cornflowerblue")

ct_lineage_plot_list <- lapply(c("Endothelial", "Epithelial", "Immune", "Mesenchymal"), function(XX) {
  # Get strip color
  strip_color = lineage_color_list[[XX]]
  
  plt <- xenium@meta.data %>%
    group_by(broad_CT5, lineage) %>%
    summarize(count = length(fine_CT4)) %>%
    mutate(round_count = ifelse(count >= 1000,
                                paste0(format(round(count/1000, 1), n.small = 2, decimal.mark = "."), "k"),
                                as.character(count))) %>%
    filter(lineage == XX) %>%
    ggplot(aes(x = reorder(broad_CT5, -count), y = count, fill = broad_CT5)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~lineage, scales = "free") +
    scale_fill_manual(values = color_list) +
    # scale_x_discrete(labels = rename_cells) +
    scale_y_continuous(limits = c(0, 170000), expand = expansion(mult = c(0, 0)),
                       breaks = seq(0, 160000, 40000), 
                       labels = c("0", "40k", "80k", "120k", "160k")) +
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
    geom_text(aes(label = round_count, vjust = -0.4), size = 1.55)
  
  # Fix y axis label
  if (XX == "Endothelial" | XX == "Epithelial") {
    plt <- plt + 
      labs(y = "Number of Cells") +
      theme(plot.margin = margin(5.5, 0, 0, 3, "pt"))
  } else {
    plt <- plt + theme(axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       plot.margin = margin(2, 0, 0, 4, "pt"))
  }
  
  plt
})

a <- wrap_plots(ct_lineage_plot_list[[1]], ct_lineage_plot_list[[3]], ncol = 2, widths = c(0.55, 1))
b <- wrap_plots(ct_lineage_plot_list[[2]], ct_lineage_plot_list[[4]], ncol = 2, widths = c(1, 0.78))
wrap_plots(a, b, ncol = 1)

CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/Fig2A_ct_prop_and_compare.pdf", width = (110*0.0393701), height = (80*0.0393701))
wrap_plots(a, b, ncol = 1)
dev.off()


#### SUPPLEMENTARY FIGURE 2B (Cell Type Proportions vs. Other Lung Studies) ----
# Add category to metadata
xenium$disease_status <- "Disease"
xenium$disease_status[xenium$sample_type == "Unaffected"] <- "Control"

# Load in and tidy data for comparison
lung_ct_compare <- readxl::read_excel("/Users/avannan/Downloads/CT Recovery Spatial vs SC_10-10-23.xlsx", sheet = 1, skip = 2) %>%
  dplyr::rename(Overall_Study = `...1`) %>%
  mutate(across(3:ncol(.), function(x) as.numeric(x))) %>%
  pivot_longer(3:ncol(.), names_to = "Disease_Status", values_to = "Percentage") %>%
  filter(!grepl("%...15|%...16|%...17", Disease_Status)) %>%
  mutate(Lineage = ordered(case_when(grepl("%...3|%...4|%...5", Disease_Status) ~ "Epithelial",
                                     grepl("%...6|%...7|%...8", Disease_Status) ~ "Endothelial",
                                     grepl("%...9|%...10|%...11", Disease_Status) ~ "Immune",
                                     grepl("%...12|%...13|%...14", Disease_Status) ~ "Mesenchymal"),
                           levels = c("Endothelial", "Epithelial", "Immune", "Mesenchymal")),
         Disease_Status = ordered(case_when(grepl("Disease", Disease_Status) ~ "Disease",
                                            grepl("Control", Disease_Status) ~ "Control",
                                            grepl("Total", Disease_Status) ~ "Total"),
                                  levels = c("Control", "Disease"))) %>%
  filter(!is.na(Overall_Study), Overall_Study != "IPF Cell Atlas", Disease_Status != "Total")

# Get dataframe of sample types for each sample
sample_type_df <- xenium@meta.data %>%
  select(sample, disease_status) %>%
  unique() %>%
  dplyr::rename(Sample = "sample", Disease_Status = "disease_status")
rownames(sample_type_df) <- NULL

# Dataframe with cell lineage proportions for each sample for spatial data
spatial_df <- (table(xenium$sample, xenium$lineage)/
                 rowSums(table(xenium$sample, xenium$lineage))) %>%
  as.data.frame() %>%
  dplyr::rename(Sample = "Var1", Lineage = "Var2", Percentage = "Freq") %>%
  mutate(Percentage = Percentage*100, Overall_Study = "Spatial",
         Lineage = ordered(Lineage, levels = c("Endothelial", "Epithelial", "Immune", "Mesenchymal"))) %>%
  full_join(sample_type_df, relationship = "many-to-many")

# Compare to eQTL study as well
eqtl <- read.csv("/Users/avannan/Documents/GSE227136_ILD_all_celltypes_seurat_meta.csv") %>%
  select(Sample_Name, Status, lineage) %>%
  mutate(Lineage = case_when(lineage %in% c("MyoFB", "MyoFB - Activated", "PLIN2+ FB") ~ "Mesenchymal",
                             lineage == "PNEC" ~ "Epithelial",
                             TRUE ~ lineage),
         Dataset = Sample_Name) %>%
  filter(Lineage != "Inflamed") %>% # Lineage is not clear
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
overall_prop_table <- table(xenium$disease_status, xenium$lineage)/rowSums(table(xenium$disease_status, xenium$lineage))
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
  facet_grid(Lineage~Disease_Status, scales = "free", space = "free")
lung_compare_plot

# Change facet colors
lung_compare_plot_table <- ggplot_gtable(ggplot_build(lung_compare_plot))
striprt <- which(grepl("strip-r", lung_compare_plot_table$layout$name) | 
                   grepl("strip-t", lung_compare_plot_table$layout$name))
fills <- c(NA, NA, "chocolate2", "chartreuse4", "maroon3", "cornflowerblue") # Lineages & Sample Types
colors <- c(rep(NA, 2), rep("black", 4))
font_colors <- c(rep("black", 2), rep("white", 4))
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", lung_compare_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  lung_compare_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  lung_compare_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- colors[k]
  lung_compare_plot_table$grobs[[i]]$grobs[[1]]$children[[j+1]]$children[[1]]$gp$col <- font_colors[k]
  k <- k+1
}
grid::grid.draw(lung_compare_plot_table)

CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/Fig2B_ct_prop_and_compare.pdf", width = (68*0.0393701), height = (77*0.0393701))
plot(lung_compare_plot_table)
dev.off()


#### SUPPLEMENTARY FIGURE 3: Dotplot Heatmap Cell Type Markers ----
# Pick out genes for heatmap
genes_4heat <- c("ACKR1", "AGER", "APLN", "BANK1", "BMPR2", 
                 "C20orf85", "CA4", "CCL2", "CD4", "CD14", "CD19", "CD1C", 
                 "CD34", "CD3E", "CD68", "CD79A", "CD8A",
                 "CFTR", "COL1A1", 
                 "COL1A2", "COL3A1", "CPA3", "CTHRC1",
                 "DCN", "EPCAM", "FABP4", "FCER1G", 
                 "FCGR3A", "FN1", "FOXJ1", "GNLY", "GZMA", "GZMB",
                 "HES1", "HEY1", "HLA-DQA1", "HLA-DRA", "IL7R",
                 "ITGAV", "KRT17", "KRT5",
                 "KRT8", "LILRA4", "LUM", "MARCO",
                 "MKI67", "MMP7", "MS4A1", "MSLN", "MUC5B",
                 "NKG7", "NKX2-1", "PECAM1", "PI16", 
                 "PLVAP", "PPARG", "PTPRC", "RTKN2", "S100A8", "S100A9", 
                 "SCGB1A1", "SCGB3A2", "SFTPC", "SOX4",
                 "SPP1", "TOP2A", "TP63", "TP73", 
                 "TPSAB1", "WNT5A")

# Log-normalizing and scaling all features in the RNA assay
# Scaling so that all features can be visualized using the same color scale
xenium <- ScaleData(xenium)

# Create dotplot base
p <- DotPlot(xenium, features = genes_4heat, 
             group.by = "broad_CT5", scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
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

# Any value that is greater than 2 will be mapped to bright red
col_fun = circlize::colorRamp2(c(-1, 0, 2), colorspace::diverge_hsv(3))

# Creating a layer to add to plot
layer_fun = function(j, i, x, y, w, h, fill){
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex(exp_mat, i, j))), # t(exp_mat)
              size = pindex(percent_mat, i, j)/100 * unit(2.25, "mm"), # t(percent_mat)
              pch = 19)
}

# List of legends
lgd_list = list(
  Legend(labels = c(0,0.25,0.5,0.75,1), title = "Proportion", type = "points", 
         pch = 19, size = c(0,0.25,0.5,0.75,1) * unit(2.25, "mm"),
         legend_gp = gpar(col = "black"), direction = "vertical", ncol = 1,
         title_position = "topcenter", background = NA, grid_height = unit(1, "mm"),
         grid_width = unit(1, "mm"),
         title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6),
         legend_height = unit(12, "mm")))


# Create heatmap
test <- c(rep(" Endothelial ", 6), rep(" Epithelial ", 11), " Lymphoid ", 
          " Myeloid ", rep(" Lymphoid ", 5), rep(" Myeloid ", 6), rep("Mesenchymal", 9))
hp <- Heatmap(exp_mat, # t(exp_mat)
              name = "Scaled Expression",
              height = unit(143, "mm"),
              heatmap_legend_param = list(title = "Scaled",
                                          legend_direction = "vertical",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 6),
                                          labels_gp = gpar(fontsize = 6),
                                          grid_width = unit(1, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),
              row_title = " ",
              row_title_gp = gpar(fontsize = 5),
              column_title_gp = gpar(fontsize = 7),
              row_title_side = "right",
              column_title_side = "top",
              column_names_rot = 45,
              column_split = test,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 6),
              row_names_side = "left",
              column_names_gp = gpar(fontsize = 6),
              cluster_rows = TRUE,  cluster_columns = FALSE,
              column_dend_height = unit(1, "mm"),
              border = "black",
              gap = unit(1, "mm"))

pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_marker_gene_expr.pdf", width = 180*0.0393701, height = 170*0.0393701)
ht <- draw(hp, 
           heatmap_legend_list = lgd_list,
           heatmap_legend_side = "right",
           annotation_legend_side = "right",
           align_heatmap_legend = "global_center")
dev.off()


## SUPPLEMENTARY FIGURE 5 (Cell Type Proportions by Sample & Sample Type) ----
# Cell type proportions
a <- xenium@meta.data %>%
  mutate(new_category = ordered(case_when(broad_CT5 %in% c("AT1", "AT2", "Transitional AT2", "KRT5-/KRT17+") ~ "Alveolar",
                                          lineage == "Epithelial" & !(broad_CT5 %in% c("AT1", "AT2", "Transitional AT2", "KRT5-/KRT17+")) ~ "Airway",
                                          broad_CT5 %in% c("T-cells", "NK cells", "B cells", "Plasma", "pDCs") ~ "Lymphoid",
                                          lineage == "Immune" & !(broad_CT5 %in% c("T-cells", "NK cells", "B cells", "Plasma", "pDCs")) ~ "Myeloid",
                                          TRUE ~ lineage), 
                                levels = c("Airway", "Alveolar", "Lymphoid", "Myeloid", "Endothelial", "Mesenchymal"))) %>%
  group_by(new_category) %>%
  mutate(count = length(new_category),
         sample_type = ordered(sample_type, levels = c("Unaffected", "LF", "MF", "ILD"))) %>%
  ggplot(aes(x = sample_type, fill = new_category)) +
  geom_bar(position = position_fill(), width = 1, color = "black", size = 0.25) +
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
        legend.position = "bottom",
        plot.margin = margin(10,5.5,5.5,5.5, "pt"))
a

# Cell type proportions by sample type
b <- xenium@meta.data %>%
  mutate(new_category = ordered(case_when(broad_CT5 %in% c("AT1", "AT2", "Transitional AT2", "KRT5-/KRT17+") ~ "Alveolar Epithelium",
                                          lineage == "Epithelial" & !(broad_CT5 %in% c("AT1", "AT2", "Transitional AT2", "KRT5-/KRT17+")) ~ "Airway Epithelium",
                                          broad_CT5 %in% c("T-cells", "NK cells", "B cells", "Plasma", "pDCs") ~ "Lymphoid Immune",
                                          lineage == "Immune" & !(broad_CT5 %in% c("T-cells", "NK cells", "B cells", "Plasma", "pDCs")) ~ "Myeloid Immune",
                                          TRUE ~ lineage), 
                                levels = c("Airway Epithelium", "Alveolar Epithelium", "Lymphoid Immune", "Myeloid Immune", "Endothelial", "Mesenchymal"))) %>%
  group_by(new_category) %>%
  mutate(count = length(new_category)) %>%
  ggplot(aes(x = reorder(sample, adjusted_pathology_score), fill = new_category)) +
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

ggsave("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig3A_ct_sample_type.pdf", width = 67.5*0.0393701, height = 85*0.0393701, device = "pdf", units = "in")
a
dev.off()

ggsave("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig3B_ct_sample.pdf", width = 112.5*0.0393701, height = 80*0.0393701, device = "pdf", units = "in")
b
dev.off()



## SUPPLEMENTARY FIGURE 6 (Cell Type by Sample Type Boxplots) ----
# Load dataframes
celltype_counts_threegs <- read.csv(file ="/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/cell_type_cell_niche_prop_test_Nov14/output/savedFiles/celltype_counts_threeGs_fixed105LF.csv")
propeller_res <- read.csv(file ="/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/manuscript_submission/cell_type_cell_niche_prop_test_Nov14/output/savedFiles/celltype_counts_threeGs_propeller_res_fixed105LF.csv") %>%
  select(-X)

# Prepare dataframes
celltypeprop_df <- celltype_counts_threegs %>%
  group_by(broad_CT5, sample) %>% 
  summarise(celltype_count = n()) %>% 
  group_by(sample) %>%
  mutate(sample_total_cells = sum(celltype_count)) %>%
  mutate(celltype_prop_in_sample = celltype_count/sample_total_cells) %>%
  dplyr::left_join(unique(celltype_counts_threegs[,c("sample","sample_type")])) %>%
  filter(sample_type!="ILD") %>%
  mutate(broad_CT5 = case_when(broad_CT5 == "MUC5B+ Secretory" ~ "MUC5B+",
                               broad_CT5 == "DCs/Lymphatic" ~ "Lymphatic/DCs",
                               broad_CT5 == "Monocyte-derived Macrophages" ~ "Interstitial Macrophages",
                               TRUE ~ broad_CT5))
propeller_res <- propeller_res %>%
  mutate(broad_CT5 = case_when(broad_CT5 == "MUC5B+ Secretory" ~ "MUC5B+",
                               broad_CT5 == "DCs/Lymphatic" ~ "Lymphatic/DCs",
                               broad_CT5 == "Monocyte-derived Macrophages" ~ "Interstitial Macrophages",
                               TRUE ~ broad_CT5))
propeller_res2 <- propeller_res %>%
  right_join(celltypeprop_df, by = c("broad_CT5")) %>%
  group_by(broad_CT5) %>%
  summarise(max_pro = max(celltype_prop_in_sample)+1.1*sqrt(var(celltype_prop_in_sample)),
            max_pro2 = max(celltype_prop_in_sample)+0.45*sqrt(var(celltype_prop_in_sample))) %>%
  ungroup() %>%
  left_join(propeller_res) %>%
  mutate(label_symbol = ifelse(is.na(label_symbol), "", label_symbol))

# Plot
sample_type_color_list <- list(`Unaffected` = "#88CCEE",
                               `LF` = "pink",
                               `MF` = "#ac3e6f",
                               `ILD` = "#d76160")
celltype_boxplot <- celltypeprop_df %>% 
  mutate(sample_type = factor(sample_type, level =c("Unaffected","LF","MF"))) %>%
  ggplot()+
  geom_blank(data = propeller_res2, mapping = aes(x = broad_CT5, y = max_pro))+
  geom_text(data = propeller_res2,
            mapping = aes(x = broad_CT5, y = max_pro2, label = label_symbol), size=4) +
  geom_boxplot(mapping = aes(x = broad_CT5,y = celltype_prop_in_sample, fill = sample_type),
               outlier.shape = 21, size = 0.4, outlier.size = 1.25,
               position = position_dodge2(preserve = "single"))+
  ylab("Proportion") +
  theme_classic(base_size = 6.5)+
  facet_wrap(.~broad_CT5, scales = "free", ncol = 5)+
  theme(strip.background.x = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black", size = 5.5),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.06),
        legend.direction = "horizontal")+
  scale_fill_manual(values = sample_type_color_list) +
  guides(fill = guide_legend(label.position = "top", override.aes = list(size = 2)))
celltype_boxplot

ggsave("/Users/avannan/Documents/spatial_figures/manuscript/Supp4_celltype_boxplot.pdf", width = 180, height = 170, units = "mm", device = "pdf")
celltype_boxplot
dev.off()


## SUPPLEMENTARY FIGURE 7 (CT Proportion Correlation with Pathology Score) ----
ctprop_path_df <- (table(xenium$sample, xenium$broad_CT5)/rowSums(table(xenium$sample, xenium$broad_CT5))) %>%
  as.data.frame() %>%
  dplyr::rename(sample = "Var1", broad_CT5 = "Var2", Prop = "Freq") %>%
  full_join(xenium@meta.data %>% 
              select(sample, sample_type, adjusted_pathology_score) %>% 
              unique())
ctprop_path_plot <- ctprop_path_df %>%
  ggplot(aes(x = adjusted_pathology_score, y = Prop)) +
  geom_point(aes(fill = broad_CT5), shape = 21) +
  geom_smooth(aes(color = broad_CT5), method = "lm") +
  theme_classic(base_size = 6) +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  scale_color_manual(values = color_list) +
  scale_fill_manual(values = color_list) +
  labs(x = "Pathology Score", y = "Proportion") +
  facet_wrap(~broad_CT5, scales = "free_y", ncol = 5) +
  NoLegend()
ctprop_path_plot

ggsave("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_ctprop_pathology_corr.pdf", 
       width = 180*0.0393701, height = 170*0.0393701, device = "pdf", units = "in")
ctprop_path_plot
dev.off()

## SUPPLEMENTARY FIGURE 10 (Annotation Counts) ----
anno_colors <- c("#7CA2DE", "#C1A456", "#E1895D", "#B4E952", "#D9BDDE", "#70E194", 
                 "#907298", "#91CCE3", "#62918D", "#E48CA5", "#7C3FE0", "#E8DC5C",
                 "#DC90E2", "#D1E3A5", "#68E95C", "#D4B39F", "#7BE5D5", "#DC4CA6",
                 "#D5E5DD", "#D14746", "#776FD8", "#D648E0")

# Number of Nuclei
num_nuclei_df <- table(xenium$Annotation_Type) %>%
  as.data.frame() %>%
  dplyr::rename("Pathology" = "Var1") %>%
  full_join(annotation_df %>% dplyr::rename(Pathology = "Annotation_Type")) %>%
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   Pathology == "smaller_airway" ~ "Small Airway",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Microscopic Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS"))) %>%
  filter(!is.na(Pathology), !is.na(Path2)) %>%
  select(Path2, Freq) %>%
  unique()
tmp <- num_nuclei_df %>% 
  select(Path2, Freq) %>% 
  mutate(Pos = ifelse(Freq < 1000, Freq+170, Freq+220)) %>% 
  unique()
a1 <- num_nuclei_df %>%
  ggplot(aes(y = reorder(Path2, dplyr::desc(Path2)), x = Freq, fill = Path2)) +
  geom_col() +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), position = "bottom", limits = c(0, 6000)) +
  scale_fill_manual(values = anno_colors) +
  scale_y_discrete(position = "left") +
  theme_classic(base_size = 6) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.line.x.top = element_line(color = "black"),
        panel.grid.major.x = element_line(color = "grey65"),
        plot.margin = margin(5.5,2,8,5.5, "pt")) +
  labs(x = "Number of Nuclei") +
  geom_text(data = tmp, aes(x = Pos, label = Freq), size = 2) +
  NoLegend()

# Number of Annotations
anno_num_df <- annotation_df %>%
  select(Sample, Annotation_Type, Annotation_Type_Instance) %>%
  unique()
a2 <- table(anno_num_df$Annotation_Type) %>%
  as.data.frame() %>%
  dplyr::rename(Pathology = "Var1") %>%
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   Pathology == "smaller_airway" ~ "Small Airway",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Microscopic Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS"))) %>%
  filter(!is.na(Pathology), !is.na(Path2)) %>%
  ggplot(aes(y = reorder(Path2, dplyr::desc(Path2)), x = Freq, fill = Path2)) +
  geom_col() +
  theme_classic(base_size = 6) +
  scale_x_continuous(expand = c(0, 0), position = "bottom", breaks = seq(0, 36, 2)) +
  scale_y_discrete(position = "left") +
  scale_fill_manual(values = anno_colors) +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.line.x.top = element_line(color = "black"),
        panel.grid.major.x = element_line(color = "grey65"),
        plot.margin = margin(5.5,5.5,8,2, "pt")) +
  labs(x = "Number of Annotations") +
  NoLegend()


# Per sample info
sample_order <- xenium@meta.data %>%
  select(sample, adjusted_pathology_score) %>%
  unique() %>%
  arrange(adjusted_pathology_score, sample) %>%
  pull(sample)
b1 <- annotation_df %>%
  select(Sample, Annotation_Type_Instance, Annotation_Type) %>%
  full_join(xenium@meta.data %>% select(sample, adjusted_pathology_score) %>% unique() %>% dplyr::rename(Sample = "sample")) %>%
  unique() %>%
  mutate(Pathology = ordered(Annotation_Type, levels = new_annotation_order)) %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   Pathology == "smaller_airway" ~ "Small Airway",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Microscopic Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS")),
         adjusted_pathology_score = as.numeric(adjusted_pathology_score)) %>%
  arrange(adjusted_pathology_score, Sample) %>%
  filter(!is.na(Pathology), !is.na(Path2)) %>%
  ggplot(aes(x = reorder(Sample, adjusted_pathology_score), fill = Path2)) +
  scale_fill_manual(values = anno_colors) +
  geom_bar() +
  theme_classic(base_size = 6) +
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
                             keyheight = unit(3, "mm")))

b2 <- annotation_df %>%
  select(Sample, Annotation_Type_Instance, Annotation_Type) %>%
  full_join(xenium@meta.data %>% select(sample, adjusted_pathology_score) %>% unique() %>% dplyr::rename(Sample = "sample")) %>%
  unique() %>%
  mutate(Pathology = ordered(Annotation_Type, levels = new_annotation_order)) %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   Pathology == "smaller_airway" ~ "Small Airway",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Microscopic Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS")),
         adjusted_pathology_score = as.numeric(adjusted_pathology_score)) %>%
  arrange(adjusted_pathology_score, Sample) %>%
  filter(!is.na(Pathology), !is.na(Path2)) %>%
  ggplot(aes(x = reorder(Sample, adjusted_pathology_score), fill = Path2)) +
  scale_fill_manual(values = anno_colors) +
  geom_bar(position = position_fill(), color = "black", width = 1) +
  theme_classic(base_size = 6) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title = element_blank(),
        plot.margin = margin(5.5,5.5,5.5,7, "pt")) +
  labs(y = "Proportion of Annotations") +
  NoLegend()

pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_A_num_anno.pdf",  width = 180*0.0393701, height = 85*0.0393701)
a1+a2
dev.off()

pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_B_num_anno.pdf",  width = 180*0.0393701, height = 85*0.0393701)
b1+b2
dev.off()



## SUPPLEMENTARY FIGURE 12 (Niche by Sample) ----
sample_tranx_niche_plot <- xenium@meta.data %>%
  filter(sample_type != "ILD") %>%
  mutate(TNiche = ordered(TNiche, levels = paste0("T", 1:12))) %>%
  ggplot(aes(x = reorder(sample, adjusted_pathology_score), fill = TNiche)) +
  geom_bar(position = position_fill(), width = 1, color = "black", size = 0.25) +
  scale_fill_manual(values = transcript_niche_color_list) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of Cells") +
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
        plot.margin = margin(10,5.5,5.5,5.5, "pt"),
        legend.key.size = unit(4, "mm"))
sample_tranx_niche_plot

sample_cell_niche_plot <- xenium@meta.data %>%
  filter(sample_type != "ILD") %>%
  mutate(CNiche = ordered(CNiche, levels = paste0("C", 1:12))) %>%
  ggplot(aes(x = reorder(sample, adjusted_pathology_score), fill = CNiche)) +
  geom_bar(position = position_fill(), width = 1, color = "black", size = 0.25) +
  scale_fill_manual(values = nuclei_niche_color_list) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of Cells") +
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
        plot.margin = margin(10,5.5,5.5,5.5, "pt"),
        legend.key.size = unit(4, "mm"))
sample_cell_niche_plot

ggsave("/Users/avannan/Documents/spatial_figures/manuscript/test.pdf", 
       width = 180*0.0393701, height = 82*0.0393701, device = "pdf", units = "in")
sample_tranx_niche_plot + sample_cell_niche_plot
dev.off()


## SUPPLEMENTARY FIGURE 14 (Niche x CT Dotplots - Niche Composition) ----
# Transcript Niches
tniche_ct_prop_df <- as.data.frame((table(xenium$TNiche, xenium$broad_CT5)/
                                      rowSums(table(xenium$TNiche, xenium$broad_CT5))))
names(tniche_ct_prop_df) <- c("Niche", "CT", "Proportion")
tniche_ct_prop_plot <- tniche_ct_prop_df %>%
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
  scale_x_discrete(labels = rename_cells) +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw(base_size = 6) +
  scale_color_manual(values = c("black", "grey80", rep("black", 10)), guide = "none") +
  scale_fill_manual(values = transcript_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 6),
        strip.text = element_text(color = "white")) +
  ggh4x::force_panelsizes(cols = c(2.1, 3.8, 4.1, 3.2)) + 
  labs(title = "Niche Composition") +
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  NoLegend()
tniche_ct_prop_plot

# Change facet colors
tniche_ct_prop_df_table <- ggplot_gtable(ggplot_build(tniche_ct_prop_plot))
striprt <- which(grepl("strip-r", tniche_ct_prop_df_table$layout$name) |
                   grepl("strip-t", tniche_ct_prop_df_table$layout$name))
fills <- rep(c("chocolate2", "chartreuse4", "maroon3", "cornflowerblue"), 2)
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", tniche_ct_prop_df_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  tniche_ct_prop_df_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

# Cell Niches
niche_ct_prop_df <- as.data.frame((table(xenium$CNiche, xenium$broad_CT5)/
                                     rowSums(table(xenium$CNiche, xenium$broad_CT5))))
names(niche_ct_prop_df) <- c("Niche", "CT", "Proportion")

niche_ct_prop_plot <- niche_ct_prop_df %>%
  full_join(lineage_df) %>%
  mutate(CT = ordered(CT, levels = new_ct_order),
         Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  ggplot(aes(y = reorder(Niche, dplyr::desc(Niche)), x = CT, 
             size = Proportion, color = as.factor(Niche), 
             alpha = Proportion_Visible,
             fill = as.factor(Niche))) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  facet_grid(~Lineage, scales = "free", space = "free") +
  theme_bw(base_size = 6) +
  scale_x_discrete(labels = rename_cells) +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_color_manual(values = c(rep("black", 7), "grey80", rep("black", 4)), guide = "none") +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
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
        legend.direction = "horizontal",
        legend.position = "bottom") + 
  ggh4x::force_panelsizes(cols = c(2.1, 3.8, 4.1, 3.2)) + 
  scale_alpha_manual(values = c(0, 1), guide = "none")
niche_ct_prop_plot

CairoPDF("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_niche_composition_celltypes.pdf", width = (180*0.0393701), height = (100*0.0393701))
cowplot::plot_grid(tniche_ct_prop_df_table, niche_ct_prop_plot, ncol = 1, align = "hv")
dev.off()


## SUPPLEMENTARY FIGURE 15 (Annotations, Cells, and Niches) ----
#### SUPPLEMENTARY FIGURE 15A (Annotation x CT Dotplot - FULL) ----
# Create dataframe of lineage variables
lineage_df <- xenium@meta.data %>%
  select(broad_CT5, lineage) %>%
  unique()
names(lineage_df) <- c("CT", "Lineage")

new_annotation_order <- c("artery",
                          "muscularized_artery",
                          "vein",
                          "venule",
                          "airway_smooth_muscle",
                          "interlobular_septum",
                          "fibroblastic_focus",
                          "severe_fibrosis",
                          "normal_alveoli",
                          "minimally_remodeled_alveoli",
                          "hyperplastic_aec",
                          "epithelial_detachment",
                          "multinucleated_cell",
                          "remodeled_epithelium",
                          "small_airway",
                          "microscopic_honeycombing",
                          "goblet_cell_metaplasia",
                          "giant_cell",
                          "granuloma",
                          "mixed_inflammation",
                          "TLS")


# Annotation Composition
path_ct_prop_df <- as.data.frame((table(xenium$Annotation_Type, xenium$broad_CT5)/
                                    rowSums(table(xenium$Annotation_Type, xenium$broad_CT5))))
names(path_ct_prop_df) <- c("Pathology", "CT", "Proportion")

path_ct_prop_plot <- path_ct_prop_df %>%
  left_join(lineage_df) %>%
  filter(Pathology != "remnant_alveoli") %>% # Remove annotations with only 1 example
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  filter(Pathology != "small_airway") %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   Pathology %in% c("smaller_airway", "small_airway") ~ "Small Airway",
                                   # Pathology == "microscopic_honeycombing" ~ "Micro. Honeycombing",
                                   # Pathology == "minimally_remodeled_alveoli" ~ "Min. Remodeled Alveoli",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Microscopic Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS")),
         CT = ordered(CT, levels = new_ct_order)) %>%
  droplevels() %>%
  ggplot(aes(x = Path2,
             y = reorder(CT, dplyr::desc(CT)),
             size = Proportion, fill = CT)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6) +
  scale_fill_manual(values = color_list) +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0.01, 1)) +
  scale_y_discrete(labels = rename_cells) +
  theme(axis.text = element_text(color = "black", size = 4.5),
        strip.text = element_text(color = "white", size = 4.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.y = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, size = 6),
        legend.box.background = element_blank(),
        legend.key.size = unit(3, "mm"),
        legend.background = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(5.5, 5.5, 5.5, 10, "pt"),
        legend.margin = margin(1.5, 1.5, 1.5, 1.5, "pt")) +
  labs(title = "Annotation Composition (Cells)") +
  guides(fill = "none", 
         size = guide_legend(title.theme = element_text(size = 5),
                             size = 2,
                             nrow = 1, title.hjust = 0.5,
                             label.position = "bottom",
                             label.theme = element_text(size = 3.5))) +
  facet_grid(Lineage~., space = "free", scales = "free")
path_ct_prop_plot

# Change facet colors
path_ct_prop_plot_table <- ggplot_gtable(ggplot_build(path_ct_prop_plot))
striprt <- which(grepl("strip-r", path_ct_prop_plot_table$layout$name) | 
                   grepl("strip-t", path_ct_prop_plot_table$layout$name))
fills <- c("chocolate2", "chartreuse4", "maroon3", "cornflowerblue") # Sample Types
k <- 1
for (i in striprt) {
  j <- which(grepl("rect", path_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$childrenOrder))
  path_ct_prop_plot_table$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(path_ct_prop_plot_table)

pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_anno_ct_comp_FULL.pdf",
    width = (90*0.0393701), height = (100*0.0393701))
grid::grid.draw(path_ct_prop_plot_table)
dev.off()


#### SUPPLEMENTARY FIGURE 15B (Annotation Composition - Niches - FULL) ----
# Transcript Niches (Annotation Composition)
tniche_path_prop_df <- as.data.frame((table(xenium$Annotation_Type, xenium$TNiche)/
                                        rowSums(table(xenium$Annotation_Type, xenium$TNiche))))
names(tniche_path_prop_df) <- c("Pathology", "Niche", "Proportion")
tniche_path_prop_df <- tniche_path_prop_df %>%
  filter(Niche != "TNA", Niche != "NA", !is.na(Niche))
tniche_path_prop_plot <- tniche_path_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("T", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Microscopic Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS"))) %>%
  filter(!is.na(Path2)) %>%
  ggplot(aes(x = Niche, y = reorder(Path2, dplyr::desc(Path2)), 
             size = Proportion, alpha = Proportion_Visible, fill = Niche)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6) +
  scale_fill_manual(values = transcript_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_text(size = 5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(), 
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(5.5, 2, 5.5, 5.5, "pt")) + 
  scale_alpha_manual(values = c(0, 1), guide = "none") +
  NoLegend()
tniche_path_prop_plot

# Cell Niches (Annotation Composition)
cniche_path_prop_df <- as.data.frame((table(xenium$Annotation_Type, xenium$CNiche)/
                                        rowSums(table(xenium$Annotation_Type, xenium$CNiche))))
names(cniche_path_prop_df) <- c("Pathology", "Niche", "Proportion")
cniche_path_prop_plot_with_legend <- cniche_path_prop_df %>%
  mutate(Niche = ordered(Niche, levels = paste0("C", 1:12)),
         Proportion_Visible = ifelse(Proportion >= 0.01, TRUE, FALSE)) %>%
  mutate(Pathology = ordered(Pathology, levels = new_annotation_order)) %>%
  mutate(Path2 = ordered(case_when(Pathology == "TLS" ~ "TLS",
                                   Pathology == "hyperplastic_aec" ~ "Hyperplastic AECs",
                                   TRUE ~ str_to_title(gsub("_", " ", Pathology))),
                         levels = c("Artery", "Muscularized Artery", "Vein", "Venule",
                                    "Airway Smooth Muscle", "Interlobular Septum",
                                    "Fibroblastic Focus", "Severe Fibrosis",
                                    "Normal Alveoli", "Minimally Remodeled Alveoli", "Hyperplastic AECs",
                                    "Epithelial Detachment", "Multinucleated Cell",
                                    "Remodeled Epithelium", "Microscopic Honeycombing",
                                    "Small Airway", "Goblet Cell Metaplasia", "Giant Cell",
                                    "Granuloma", "Mixed Inflammation", "TLS"))) %>%
  filter(!is.na(Path2)) %>%
  ggplot(aes(x = Niche, y = reorder(Path2, dplyr::desc(Path2)), 
             size = Proportion, alpha = Proportion_Visible, fill = Niche)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6) +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_size_continuous(range = c(0, 2), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
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
        legend.position = "bottom", legend.title.align = 0.5
        # legend.key.width = unit(3.5, "mm"),
        # legend.key.height = unit(2, "mm")
  ) +
  guides(size = guide_legend(title.position = "top", label.position = "bottom")) +
  scale_alpha_manual(values = c(0, 1), guide = "none")
cniche_path_prop_plot <- cniche_path_prop_plot_with_legend + NoLegend()
cniche_path_prop_plot_with_legend


# Combine and save
pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_niche_anno_prop_FULL.pdf", width = (90*0.0393701), height = (100*0.0393701))
(tniche_path_prop_plot + cniche_path_prop_plot) + 
  patchwork::plot_layout(widths = unit(c(31, 31), "mm"), 
                         heights = unit(c(81, 81), "mm"), 
                         guides = "keep", ncol = 2, nrow = 1)
dev.off()


pdf("/Users/avannan/Documents/spatial_figures/manuscript/test.pdf", width = (90*0.0393701), height = (100*0.0393701))
cniche_path_prop_plot_with_legend
dev.off()




## SUPPLEMENTARY FIGURE 16 (Heatmap of Lumens - Full) ----
# Set up final dataframe and get trajectories
df <- gene_by_lumen[c(unimodal, bimodal), lumens_to_plot]
full_gene_order <- rownames(df)
broad_trajectories <- genes_lumen_ordered_df %>% filter(gene_order_smooth %in% full_gene_order) %>% arrange(desc(gene_mode)) %>% pull(trajectory)
broad_trajectories[broad_trajectories == "Normal"] <- "  Normal  "
broad_trajectories[broad_trajectories == "Late Remodeling"] <- " Late Remodeling "
broad_trajectories[broad_trajectories == "Early Transition"] <- " Early Transition "
broad_trajectories[is.na(broad_trajectories)] <- "Bimodal"

sample_names <- data.frame(sample = gsub("_L[[:digit:]]*", "", lumens_to_plot))
pathology_scores <- left_join(sample_names, xenium@meta.data %>% select(sample, adjusted_pathology_score) %>% unique()) %>% pull(adjusted_pathology_score)

pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_heatmap_pseudotime_lumens_FULL.pdf", width = 90*0.0393701, height = 175*0.0393701)
ht <- Heatmap(df,
              cluster_columns = FALSE, cluster_rows = FALSE,
              show_row_names = FALSE,
              column_title_gp = gpar(fontsize = 0),
              show_column_names = FALSE,
              border = TRUE,
              height = unit(82, "mm"),
              width = unit(35, "mm"),
              row_split = broad_trajectories,
              row_title_gp = gpar(fontsize = 5),
              row_title_side = "left",
              
              heatmap_legend_param = list(title = "Z-Score",
                                          legend_direction = "vertical",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          labels_side = "left",
                                          grid_width = unit(1, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),
              
              left_annotation = rowAnnotation(
                Trajectory = broad_trajectories,
                annotation_name_gp = gpar(fontsize = 0),
                simple_anno_size = unit(0.15, "cm"),
                show_legend = FALSE,
                
                col = list(Trajectory = c("  Normal  " = "#04846E", 
                                          " Early Transition " = "#FFC107", 
                                          " Late Remodeling " = "#D81B60",
                                          "Bimodal" = "grey70")),
                na_col = "white"
              ),
              
              top_annotation =
                HeatmapAnnotation(
                  # Cell types
                  `Interstitial Macrophages` = ct_mat["Monocyte-derived Macrophages", ],
                  `NK cells` = ct_mat["NK cells", ],
                  `Proliferating Endothelial` = ct_mat["Proliferating Endothelial", ],
                  aCap = ct_mat["aCap", ],
                  `Alveolar FBs` = ct_mat["Alveolar FBs", ],
                  Arteriole = ct_mat["Arteriole", ],
                  AT1 = ct_mat["AT1", ],
                  `Lymphatic/DCs` = ct_mat["DCs/Lymphatic", ],
                  gCap = ct_mat["gCap", ],
                  SMCs = ct_mat["SMCs", ],
                  `T-cells` = ct_mat["T-cells", ],
                  `Adventitial FBs` = ct_mat["Adventitial FBs", ],
                  Fibroblasts = ct_mat["Fibroblasts", ],
                  `Peribronchial FBs` = ct_mat["Peribronchial FBs", ],
                  AT2 = ct_mat["AT2", ],
                  Plasma = ct_mat["Plasma", ],
                  Mast = ct_mat["Mast", ],
                  Lipofibroblasts = ct_mat["Lipofibroblasts", ],
                  `Activated FBs` = ct_mat["Activated FBs", ],
                  `Transitional AT2` = ct_mat["Transitional AT2", ],
                  Venous = ct_mat["Venous", ],
                  `SCGB3A2+` = ct_mat["SCGB3A2+", ],
                  `FABP4+ Macrophages` = ct_mat["FABP4+ Macrophages", ],
                  `Proliferating Epithelial` = ct_mat["Proliferating Epithelial", ],
                  Macrophages = ct_mat["Macrophages", ],
                  DCs = ct_mat["DCs", ],
                  `SPP1+ Macrophages` = ct_mat["SPP1+ Macrophages", ],
                  `Proliferating Mesenchymal` = ct_mat["Proliferating Mesenchymal", ],
                  Ciliated = ct_mat["Ciliated", ],
                  `Proliferating Myeloid` = ct_mat["Proliferating Myeloid", ],
                  `KRT5-/KRT17+` = ct_mat["KRT5-/KRT17+", ],
                  
                  
                  # CNiche
                  C12 = cniche_mat["12", ],
                  C3 = cniche_mat["3", ],
                  C8 = cniche_mat["8", ],
                  C10 = cniche_mat["10", ],
                  C5 = cniche_mat["5", ],
                  C7 = cniche_mat["7", ],
                  C4 = cniche_mat["4", ],
                  C6 = cniche_mat["6", ],
                  
                  # TNiche
                  T5 = tniche_mat["5", ],
                  T4 = tniche_mat["4", ],
                  T2 = tniche_mat["2", ],
                  T10 = tniche_mat["10", ],
                  T9 = tniche_mat["9", ],
                  T11 = tniche_mat["11", ],
                  T12 = tniche_mat["12", ],
                  T3 = tniche_mat["3", ],
                  T6 = tniche_mat["6", ],
                  T8 = tniche_mat["8", ],
                  
                  # Pathology Score
                  `Pathology Score` = pathology_scores,
                  
                  # Pseudotime
                  Pseudotime = pseudotime_value,
                  
                  simple_anno_size = unit(0.16, "cm"),
                  simple_anno_size_adjust = TRUE,
                  annotation_name_gp = gpar(fontsize = 5),
                  
                  col = list(
                    `Interstitial Macrophages` = colorRamp2(c(0, max(ct_mat["Monocyte-derived Macrophages", ])), c("white", "deeppink3")),
                    `NK cells` = colorRamp2(c(0, max(ct_mat["NK cells", ])), c("white", "darkmagenta")),
                    `Proliferating Endothelial` = colorRamp2(c(0, max(ct_mat["Proliferating Endothelial", ])), c("white", "chocolate2")),
                    aCap = colorRamp2(c(0, max(ct_mat["aCap", ])), c("white", "chocolate2")),
                    `Alveolar FBs` = colorRamp2(c(0, max(ct_mat["Alveolar FBs", ])), c("white", "cornflowerblue")),
                    Arteriole = colorRamp2(c(0, max(ct_mat["Arteriole", ])), c("white", "chocolate2")),
                    AT1 = colorRamp2(c(0, max(ct_mat["AT1", ])), c("white", "chartreuse3")),
                    `Lymphatic/DCs` = colorRamp2(c(0, max(ct_mat["DCs/Lymphatic", ])), c("white", "chocolate2")),
                    gCap = colorRamp2(c(0, max(ct_mat["gCap", ])), c("white", "chocolate2")),
                    SMCs = colorRamp2(c(0, max(ct_mat["SMCs", ])), c("white", "cornflowerblue")),
                    `T-cells` = colorRamp2(c(0, max(ct_mat["T-cells", ])), c("white", "darkmagenta")),
                    `Adventitial FBs` = colorRamp2(c(0, max(ct_mat["Adventitial FBs", ])), c("white", "cornflowerblue")),
                    Fibroblasts = colorRamp2(c(0, max(ct_mat["Fibroblasts", ])), c("white", "cornflowerblue")),
                    `Peribronchial FBs` = colorRamp2(c(0, max(ct_mat["Peribronchial FBs", ])), c("white", "cornflowerblue")),
                    AT2 = colorRamp2(c(0, max(ct_mat["AT2", ])), c("white", "chartreuse3")),
                    Plasma = colorRamp2(c(0, max(ct_mat["Plasma", ])), c("white", "darkmagenta")),
                    Mast = colorRamp2(c(0, max(ct_mat["Mast", ])), c("white", "deeppink3")),
                    Lipofibroblasts = colorRamp2(c(0, max(ct_mat["Lipofibroblasts", ])), c("white", "cornflowerblue")),
                    `Activated FBs` = colorRamp2(c(0, max(ct_mat["Activated FBs", ])), c("white", "cornflowerblue")),
                    `Transitional AT2` = colorRamp2(c(0, max(ct_mat["Activated FBs", ])), c("white", "chartreuse3")),
                    Venous = colorRamp2(c(0, max(ct_mat["Venous", ])), c("white", "chocolate2")),
                    `SCGB3A2+` = colorRamp2(c(0, max(ct_mat["SCGB3A2+", ])), c("white", "chartreuse3")),
                    `FABP4+ Macrophages` = colorRamp2(c(0, max(ct_mat["FABP4+ Macrophages", ])), c("white", "deeppink3")),
                    `Proliferating Epithelial` = colorRamp2(c(0, max(ct_mat["Proliferating Epithelial", ])), c("white", "chartreuse3")),
                    Macrophages = colorRamp2(c(0, max(ct_mat["Macrophages", ])), c("white", "deeppink3")),
                    DCs = colorRamp2(c(0, max(ct_mat["DCs", ])), c("white", "deeppink3")),
                    `SPP1+ Macrophages` = colorRamp2(c(0, max(ct_mat["SPP1+ Macrophages", ])), c("white", "deeppink3")),
                    `Proliferating Mesenchymal` = colorRamp2(c(0, max(ct_mat["Proliferating Mesenchymal", ])), c("white", "cornflowerblue")),
                    Ciliated = colorRamp2(c(0, max(ct_mat["Ciliated", ])), c("white", "chartreuse3")),
                    `Proliferating Myeloid` = colorRamp2(c(0, max(ct_mat["Proliferating Myeloid", ])), c("white", "deeppink3")),
                    `KRT5-/KRT17+` = colorRamp2(c(0, max(ct_mat["KRT5-/KRT17+", ])), c("white", "chartreuse3")),
                    
                    C12 = colorRamp2(c(0, max(cniche_mat["12", ])), c("white", "#924373")),
                    C3 = colorRamp2(c(0, max(cniche_mat["3", ])), c("white", "#cabd00")),
                    C8 = colorRamp2(c(0, max(cniche_mat["8", ])), c("white", "black")),
                    C10 = colorRamp2(c(0, max(cniche_mat["10", ])), c("white", "#6f774b")),
                    C5 = colorRamp2(c(0, max(cniche_mat["5", ])), c("white", "#84cd5d")),
                    C7 = colorRamp2(c(0, max(cniche_mat["7", ])), c("white", "#d10040")),
                    C4 = colorRamp2(c(0, max(cniche_mat["4", ])), c("white", "#a6513c")),
                    C6 = colorRamp2(c(0, max(cniche_mat["6", ])), c("white", "#72c7ff")),
                    
                    
                    T5 = colorRamp2(c(0, max(tniche_mat["5", ])), c("white", "#2BCE48")),
                    T4 = colorRamp2(c(0, max(tniche_mat["4", ])), c("white", "#FFCC99")),
                    T2 = colorRamp2(c(0, max(tniche_mat["2", ])), c("white", "#191919")),
                    T10 = colorRamp2(c(0, max(tniche_mat["10", ])), c("white", "#C20088")),
                    T9 = colorRamp2(c(0, max(tniche_mat["3", ])), c("white", "#4C005C")),
                    T11 = colorRamp2(c(0, max(tniche_mat["11", ])), c("white", "#4C005C")),
                    T12 = colorRamp2(c(0, max(tniche_mat["12", ])), c("white", "#9DCC00")),
                    T3 = colorRamp2(c(0, max(tniche_mat["3", ])), c("white", "#8F7C00")),
                    T6 = colorRamp2(c(0, max(tniche_mat["6", ])), c("white", "#993F00")),
                    T8 = colorRamp2(c(0, max(tniche_mat["8", ])), c("white", "#0075DC")),
                    
                    `Pathology Score` = colorRamp2(c(min(pathology_scores), max(pathology_scores)), c("white", "black")),
                    Pseudotime = colorRamp2(c(min(pseudotime_value), max(pseudotime_value)), c("white", "black"))),
                  
                  annotation_name_side = "left", show_legend = FALSE, gap = unit(0.01, "mm")))
draw(ht, heatmap_legend_side = "left", align_heatmap_legend = "heatmap_center")
dev.off()

# Second panel
df2 <- df[broad_trajectories %in% c("  Normal  ", " Early Transition "), ]
broad_trajectories2 <- broad_trajectories[which(broad_trajectories %in% c("  Normal  ", " Early Transition "))]
pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_heatmap_pseudotime_lumens_FULL_2ndpanel.pdf", width = 90*0.0393701, height = 170*0.0393701)
ht <- Heatmap(df2,
              cluster_columns = FALSE, cluster_rows = FALSE,
              row_names_gp = gpar(fontsize = 3),
              column_title_gp = gpar(fontsize = 0),
              show_column_names = FALSE,
              border = TRUE,
              height = unit(165, "mm"),
              width = unit(35, "mm"),
              row_split = broad_trajectories2,
              row_title_gp = gpar(fontsize = 5),
              row_title_side = "left",
              
              heatmap_legend_param = list(title = "Z-Score",
                                          legend_direction = "vertical",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          labels_side = "left",
                                          grid_width = unit(1, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),

              left_annotation = rowAnnotation(
                Trajectory = broad_trajectories2,
                annotation_name_gp = gpar(fontsize = 0),
                simple_anno_size = unit(0.15, "cm"),
                show_legend = FALSE,

              col = list(Trajectory = c("  Normal  " = "#04846E", 
                                        " Early Transition " = "#FFC107"))),
              na_col = "white",
              gap = unit(1, "mm"))
draw(ht, heatmap_legend_side = "left", align_heatmap_legend = "heatmap_center")
dev.off()


# Third panel
df3 <- df[broad_trajectories %in% c(" Late Remodeling ", "Bimodal"), ]
broad_trajectories3 <- broad_trajectories[which(broad_trajectories %in% c(" Late Remodeling ", "Bimodal"))]
pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_heatmap_pseudotime_lumens_FULL_3rdpanel.pdf", width = 90*0.0393701, height = 170*0.0393701)
ht <- Heatmap(df3,
              cluster_columns = FALSE, cluster_rows = FALSE,
              row_names_gp = gpar(fontsize = 4),
              column_title_gp = gpar(fontsize = 0),
              show_column_names = FALSE,
              border = TRUE,
              height = unit(165, "mm"),
              width = unit(35, "mm"),
              row_split = broad_trajectories3,
              row_title_gp = gpar(fontsize = 5),
              row_title_side = "left",
              
              heatmap_legend_param = list(title = "Z-Score",
                                          legend_direction = "vertical",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 5),
                                          labels_gp = gpar(fontsize = 5),
                                          labels_side = "left",
                                          grid_width = unit(1, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),
              
              left_annotation = rowAnnotation(
                Trajectory = broad_trajectories3,
                annotation_name_gp = gpar(fontsize = 0),
                simple_anno_size = unit(0.15, "cm"),
                show_legend = FALSE,
                
                col = list(Trajectory = c(" Late Remodeling " = "#D81B60", 
                                          "Bimodal" = "grey70"))),
              na_col = "white",
              gap = unit(1, "mm"))
draw(ht, heatmap_legend_side = "left", align_heatmap_legend = "heatmap_center")
dev.off()


## SUPPLEMENTARY FIGURE 19 (Heatmap of Late Remodeling Expression) ----
# Load in significant terms x CT and filter to only the middle tranche of genes
sig_terms_perCT <- read_csv("/Users/avannan/Downloads/sig_genes_perCT/perCT/all_perCT_sig_terms_s15.csv") %>%
  full_join(genes_lumen_ordered_df, by = c("gene" = "gene_order_smooth")) %>%
  mutate(total_sig = length(gene)) %>%
  group_by(celltype) %>%
  mutate(prop_sig_ct = length(gene)/total_sig) %>%
  ungroup()

# Subset object to contain only cells found in final lumens
xenium_lumens_cells_only <- subset(xenium, cells = rownames(lumen_metadata))

# Pick out genes for heatmap
heatmap_genes <- late_genes

# Log-normalizing and scaling all features in the RNA assay
# Scaling so that all features can be visualized using the same color scale
xenium_lumens_cells_only <- ScaleData(xenium_lumens_cells_only)

# Create dotplot base
p <- DotPlot(xenium_lumens_cells_only, features = heatmap_genes, 
             group.by = "broad_CT5", scale = TRUE) + 
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
  select(gene, celltype) %>%
  unique() %>%
  filter(gene %in% late_genes, !(gene %in% bimodal)) %>%
  dplyr::rename(broad_CT5 = "celltype") %>%
  mutate(broad_CT5 = case_when(broad_CT5 == "T.cells" ~ "T-cells",
                               broad_CT5 == "KRT5..KRT17." ~ "KRT5-/KRT17+",
                               broad_CT5 == "FABP4..Macrophages" ~ "FABP4+ Macrophages",
                               broad_CT5 == "SPP1..Macrophages" ~ "SPP1+ Macrophages",
                               broad_CT5 == "Monocyte.derived.Macrophages" ~ "Interstitial Macrophages",
                               broad_CT5 == "Transitional.AT2" ~ "Transitional AT2",
                               broad_CT5 == "Alveolar.FBs" ~ "Alveolar FBs",
                               TRUE ~ broad_CT5)) %>%
  full_join(xenium@meta.data %>% select(broad_CT5, lineage) %>% unique()) %>% 
  mutate(broad_CT5 = case_when(broad_CT5 == "Proliferating Lymphoid" ~ "Prolif. Lymphoid",
                               broad_CT5 == "Proliferating Myeloid" ~ "Prolif. Myeloid",
                               broad_CT5 == "Proliferating Mesenchymal" ~ "Prolif. Mesenchymal",
                               broad_CT5 == "Proliferating Epithelial" ~ "Prolif. Epithelial",
                               broad_CT5 == "Proliferating Endothelial" ~ "Prolif. Endothelial",
                               TRUE ~ broad_CT5),
         lineage = case_when(broad_CT5 %in% c("Mast", "FABP4+ Macrophages", "Macrophages", "SPP1+ Macrophages",
                                              "Prolif. Myeloid", "Interstitial Macrophages") ~ "Myeloid",
                             lineage == "Immune" ~ "Lymphoid",
                             TRUE ~ lineage)) %>%
  mutate(sig = 1)

# Get proportion/percentage of significant tests
prop_gene_lumen_lineage_df <- gene_lineage_lumen_heatmap_df_og %>% 
  filter(!is.na(broad_CT5)) %>%
  mutate(num_sig_tests = sum(sig)) %>%
  group_by(lineage) %>%
  mutate(num_sig_lineage = sum(sig),
         prop_sig_lineage = sum(sig)/num_sig_tests) %>%
  select(lineage, prop_sig_lineage, num_sig_lineage, num_sig_tests) %>%
  unique() %>%
  mutate(title = paste0(num_sig_lineage, "/", num_sig_tests, " (", round(prop_sig_lineage*100, 1), "%)"))

# Prep dataframe for heatmap
gene_lineage_lumen_heatmap_df <- gene_lineage_lumen_heatmap_df_og %>%
  select(-broad_CT5) %>%
  unique() %>% 
  pivot_wider(names_from = "lineage", values_from = "sig") %>%
  select(-`NA`) %>%
  filter(!is.na(gene)) %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))
gene_lineage_lumen_heatmap_df2 <- gene_lineage_lumen_heatmap_df %>%
  select(Endothelial, Epithelial, Lymphoid, Myeloid, Mesenchymal) %>%
  as.data.frame()
rownames(gene_lineage_lumen_heatmap_df2) <- gene_lineage_lumen_heatmap_df$gene

# Put genes in order
genes_in_order <- genes_lumen_ordered_df %>%
  filter(gene_mode %in% "unimodal", trajectory == "Late Remodeling") %>% 
  pull(gene_order_smooth)
gene_lineage_lumen_heatmap_df2 <- gene_lineage_lumen_heatmap_df2[genes_in_order, ]

# Fix labels
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Lymphoid"] <- "Prolif. Lymphoid"
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Myeloid"] <- "Prolif. Myeloid"
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Mesenchymal"] <- "Prolif. Mesenchymal"
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Endothelial"] <- "Prolif. Endothelial"
colnames(exp_mat)[colnames(exp_mat) == "Proliferating Epithelial"] <- "Prolif. Epithelial"

new_ct_order2 <- new_ct_order
new_ct_order2[new_ct_order2 == "Proliferating Lymphoid"] <- "Prolif. Lymphoid"
new_ct_order2[new_ct_order2 == "Proliferating Myeloid"] <- "Prolif. Myeloid"
new_ct_order2[new_ct_order2 == "Proliferating Mesenchymal"] <- "Prolif. Mesenchymal"
new_ct_order2[new_ct_order2 == "Proliferating Endothelial"] <- "Prolif. Endothelial"
new_ct_order2[new_ct_order2 == "Proliferating Epithelial"] <- "Prolif. Epithelial"

# Fix colors
new_color_list <- color_list[new_ct_order]
names(new_color_list) <- new_ct_order2

# Make heatmap
hp <- Heatmap(exp_mat[, new_ct_order2], # t(exp_mat)
              name = "Scaled Expression",
              heatmap_legend_param = list(title = "Scaled Expression",
                                          legend_direction = "horizontal",
                                          title_position = "topcenter",
                                          title_gp = gpar(fontsize = 6),
                                          labels_gp = gpar(fontsize = 6),
                                          grid_width = unit(1, "mm"),
                                          grid_height = unit(1, "mm"),
                                          legend_height = unit(15, "mm")),
              column_title_side = "top",
              row_title = " ",
              column_title_gp = gpar(fontsize = 6, col = c("chocolate2", "chartreuse4",
                                                           "darkmagenta", "deeppink3", "cornflowerblue")),
              column_names_side = "top",
              show_column_names = FALSE,
              column_names_rot = 90,
              row_names_side = "right",
              column_split = c(rep(" Endothelial  ", 6), rep(" Epithelial ", 11),
                               rep(" Lymphoid ", 7), rep(" Myeloid ", 6), rep("Mesenchymal", 9)),
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 7),
              cluster_rows = FALSE,  cluster_columns = FALSE,
              column_order = new_ct_order2,
              border = "black",
              top_annotation = 
                columnAnnotation(`Cell Type` = colnames(exp_mat[, new_ct_order2]),
                                 col = list(`Cell Type` = unlist(new_color_list)),
                                 show_legend = FALSE,
                                 annotation_name_gp = gpar(fontsize = 0),
                                 simple_anno_size = unit(1, "mm")),
              left_annotation =
                rowAnnotation(`11% -`= anno_points(gene_lineage_lumen_heatmap_df2$Endothelial, which = "row",
                                                   ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                   pch = 15,
                                                   axis = FALSE, border = FALSE, gp = gpar(col = "chocolate2"),
                                                   size = unit(1, "mm")),
                              `20.9% -` = anno_points(gene_lineage_lumen_heatmap_df2$Epithelial, which = "row",
                                                      ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                      pch = 15,
                                                      axis = FALSE, border = FALSE, gp = gpar(col = "chartreuse4"),
                                                      size = unit(1, "mm")),
                              `9.9% -` = anno_points(gene_lineage_lumen_heatmap_df2$Lymphoid, which = "row",
                                                     ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                     pch = 15,
                                                     axis = FALSE, border = FALSE, gp = gpar(col = "darkmagenta"),
                                                     size = unit(1, "mm")),
                              `37.4% -` = anno_points(gene_lineage_lumen_heatmap_df2$Myeloid, which = "row",
                                                      ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                      pch = 15,
                                                      axis = FALSE, border = FALSE, gp = gpar(col = "deeppink3"),
                                                      size = unit(1, "mm")),
                              ` 20.9% -` = anno_points(gene_lineage_lumen_heatmap_df2$Mesenchymal, which = "row",
                                                       ylim = c(0.99, 1.01), extend = 0, width = unit(1, "mm"),
                                                       pch = 15,
                                                       axis = FALSE, border = FALSE, gp = gpar(col = "cornflowerblue"),
                                                       size = unit(1, "mm")),
                              annotation_name_gp = gpar(fontsize = 5),
                              # show_annotation_name = FALSE,
                              annotation_name_rot = 90,
                              annotation_name_side = "bottom",
                              gap = unit(0.5, "mm")
                ),
              gap = unit(10, "mm")
)

# Add annotations to heatmap
ht <- draw(hp, 
           heatmap_legend_side = "bottom",
           align_heatmap_legend = "heatmap_center")
ht

pdf("/Users/avannan/Documents/spatial_figures/manuscript/SuppFig_heatmap_late_transition.pdf",
    width = (90*0.0393701), height = (170*0.0393701))
ht
dev.off()



saveRDS(xenium, "/Volumes/dback_scratch/avannan/late_IPF_spatial/xenium/Seurat_GSE250346.rds")

