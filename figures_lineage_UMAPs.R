## SETTING ENVIRONMENT ----
# library(rlang, lib.loc = "/usr/local/lib/R/site-library") # 1.1.1
library(scCustomize)
library(rlang) # 1.1.3
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


## LOAD IN FILES ----
xenium <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_final_CT_plus_niches_070924.rds")

# Save files
# endothelial <- subset(xenium, subset = final_lineage == "Endothelial")
# epithelial <- subset(xenium, subset = final_lineage == "Epithelial")
# airway <- subset(xenium, subset = final_sublineage == "Airway")
# alveolar <- subset(xenium, subset = final_sublineage == "Alveolar")
# immune <- subset(xenium, subset = final_lineage == "Immune")
# lymphoid <- subset(xenium, subset = final_sublineage == "Lymphoid")
# myeloid <- subset(xenium, subset = final_sublineage == "Myeloid")
# mesenchymal <- subset(xenium, subset = final_lineage == "Mesenchymal")
# saveRDS(endothelial, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/ENDOTHELIAL_xenium_final_CT_plus_niches_070224.rds")
# saveRDS(epithelial, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/EPITHELIAL_xenium_final_CT_plus_niches_070224.rds")
# saveRDS(airway, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/AIRWAY_xenium_final_CT_plus_niches_070224.rds")
# saveRDS(alveolar, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/ALVEOLAR_xenium_final_CT_plus_niches_070224.rds")
# saveRDS(immune, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/IMMUNE_xenium_final_CT_plus_niches_070224.rds")
# saveRDS(lymphoid, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/LYMPHOID_xenium_final_CT_plus_niches_070224.rds")
# saveRDS(myeloid, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/MYELOID_xenium_final_CT_plus_niches_070224.rds")
# saveRDS(mesenchymal, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/MESENCHYMAL_xenium_final_CT_plus_niches_070224.rds")

# Load files
# Lineages
endothelial <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/ENDOTHELIAL_xenium_final_CT_plus_niches_RAPIDS_clustered_070224.rds")
epithelial <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/EPITHELIAL_xenium_final_CT_plus_niches_RAPIDS_clustered_070524.rds")
immune <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/IMMUNE_xenium_final_CT_plus_niches_RAPIDS_clustered_070524.rds")
mesenchymal <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/MESENCHYMAL_xenium_final_CT_plus_niches_RAPIDS_clustered_070224.rds")
# Sub-lineages
airway <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/AIRWAY_xenium_final_CT_plus_niches_RAPIDS_clustered_071624.rds")
alveolar <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/ALVEOLAR_xenium_final_CT_plus_niches_RAPIDS_clustered_071624.rds")
lymphoid <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/LYMPHOID_xenium_final_CT_plus_niches_RAPIDS_clustered_071624.rds")
myeloid <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/MYELOID_xenium_final_CT_plus_niches_RAPIDS_clustered_071624.rds")


## SET COLORS ----
# Endothelial
endo_colors <- list(`Arteriole` = "#e46300",
                    `Capillary` = "#cdbf5e",
                    `Lymphatic` = "#975f47",
                    `Venous` = "#d59445")
endo_colors2 <- list(`Arteriole` = "#e22039",
                     # `Capillary` = "#c8c0f2",
                     `Capillary` = "#AFA6DE",
                     `Lymphatic` = "#0272db",
                     `Venous` = "#ccc100")

# Epithelial
alv_colors <- list(`AT1` = "#6cd66d",
                   `AT2` = "#34491d",
                   `Transitional AT2` = "#4fa381",
                   `Proliferating AT2` = "#5e6b00",
                   `KRT5-/KRT17+` = "#a0cd9e")
air_colors <- list(`Basal` = "#31a200",
                   `Multiciliated` = "#b9b700",
                   `Goblet` = "#acde47",
                   `RASC` = "#ccd87a",
                   `Secretory` = "#01eaaa",
                   `PNEC` = "#aef359",
                   `Proliferating Airway` = "#006e3e")
epi_colors2 <- list(`AT1` = "#d70043",
                   `AT2` = "#bea73f",
                   `Transitional AT2` = "#e491ff",
                   `Proliferating AT2` = "#575604",
                   `KRT5-/KRT17+` = "#0059bf",
                   `Basal` = "#ff8039",
                   `Multiciliated` = "#64c5ff",
                   `Goblet` = "#8e1e8c",
                   `RASC` = "#dfc389",
                   `Secretory` = "#db0073",
                   `PNEC` = "#bfd035",
                   `Proliferating Airway` = "#ff7dd3")

# Immune
lym_colors <- list(`CD4+ T-cells` = "#e62f4c",
                   `CD8+ T-cells` = "#dc6fd5",
                   `Tregs` = "#9e548d",
                   `B cells` = "#d84265",
                   `NK/NKT` = "#b16f6f",
                   `Plasma` = "#df4738",
                   `pDCs` = "#842b8c",
                   `Proliferating T-cells` = "#993432",
                   `Proliferating NK/NKT` = "#ab3b99",
                   `Proliferating B cells` = "#724173")
mye_colors <- list(`Alveolar Macrophages` = "#ffa3a1",
                   `Basophils` = "#d93e77",
                   `cDCs` =  "#cb547b",
                   `Langerhans cells` = "#e5849a",
                   `Interstitial Macrophages` = "#dc48bd",
                   `Macrophages - IFN-activated` = "#963095",
                   `Mast` = "#bf4ce3",
                   `Migratory DCs` = "#a94b59",
                   `Monocytes/MDMs`  = "#f2294e",
                   `Neutrophils` = "#a51219",
                   `SPP1+ Macrophages` = "#9e3c6d",
                   `Proliferating Myeloid` = "#e25292")
imm_colors2 <- list(`CD4+ T-cells` = "#fcab11",
                    `CD8+ T-cells` = "#ff507b",
                    `Tregs` = "#ae54bb",
                    `B cells` = "#ee39b0", 
                    `NK/NKT` = "#5d9c55",
                    `Plasma` = "#5753a7",
                    `pDCs` = "#93320a",
                    `Proliferating T-cells` = "#db6440",
                    `Proliferating NK/NKT` = "#ba9800",
                    `Proliferating B cells` = "#724173",
                    `Alveolar Macrophages` = "#005bbd",
                    `Basophils` = "#d7c84a",
                    `cDCs` =  "#bc547b",
                    `Langerhans cells` = "#769500",
                    `Interstitial Macrophages` = "#d19cff",
                    `Macrophages - IFN-activated` = "#61d751",
                    `Mast` = "#6f4275",
                    `Migratory DCs` = "#a80045",
                    `Monocytes/MDMs`  = "#018977",
                    `Neutrophils` = "#7291c8",
                    `SPP1+ Macrophages` = "#ffa0aa",
                    `Proliferating Myeloid` = "#a17637")

# Mesenchymal
mes_colors <- list(`Activated Fibrotic FBs` = "#005a90",
                   `Adventitial FBs` = "#00cace",
                   `Alveolar FBs` = "#6d4bce",
                   `Inflammatory FBs` = "#7895ff",
                   `Mesothelial` = "#374ca3",
                   `Myofibroblasts` = "#4c4e8c",
                   `SMCs/Pericytes` = "#00aec9",
                   `Subpleural FBs` = "#0065f5",
                   `Proliferating FBs` = "#7dfff6")
mes_colors2 <- c(`Activated Fibrotic FBs` = "#c810e8",
                 `Adventitial FBs` = "#ff57b2",
                 `Alveolar FBs` = "blanchedalmond",
                 `Inflammatory FBs` = "black",
                 `Mesothelial` = "#bed00d",
                 `Myofibroblasts` = "#24cdff",
                 `SMCs/Pericytes` = "#ab6600",
                 `Subpleural FBs` = "#00a754",
                 `Proliferating FBs` = "darkred")


## ENDOTHELIAL ----
# A: DimPlot of cell types
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ENDOTHELIAL_DimPlot.pdf", 
    width = 60*0.0393701, height = 35*0.0393701)
DimPlot(endothelial, group.by = "final_CT", cols = unlist(endo_colors2)) +
  theme_bw(base_size = 6.5) +
  theme_black +
  pretty_umap + coord_equal() + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 5),
        plot.margin = margin(0, 0, 0, 0, "mm"), panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 5.5),
        legend.key.height = unit(0.8, "mm"),
        legend.key.width = unit(0.8, "mm"))
dev.off()

# B: Proportion of endothelial cells
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ENDOTHELIAL_propcells.pdf",
    width = 44*0.0393701, height = 29.7*0.0393701)
xenium@meta.data %>%
  filter(final_lineage == "Endothelial") %>%
  mutate(sample_affect = ordered(case_when(sample_affect == "Less Affected" ~ "Less",
                                           sample_affect == "More Affected" ~ "More",
                                           TRUE ~ sample_affect),
                                 levels = c("Unaffected", "Less", "More"))) %>%
  ggplot(aes(x = sample_affect, fill = final_CT)) +
  geom_bar(position = position_fill()) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic(base_size = 6.5) +
  theme_black +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 5.5),
        legend.title = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.margin = margin(0,0,0,0, "mm")) +
  labs(y = "Prop. of Endothelial Cells") +
  scale_fill_manual(values = endo_colors2)
dev.off()

# C: Overall FeaturePlots of marker genes
endo_feats <- FeaturePlot_scCustom(endothelial, 
                     features = c("PECAM1", "BMPR2", "CD34", # General endothelial
                                  "HEY1", "VEGFA", # Arteriole
                                  "CA4", # Capillary 1
                                  "CCL21", "FABP4", # Lymphatic
                                  "FCN3", # Capillary 2
                                  "ACKR1", "PLVAP", # Venous
                                  "APLN" # Capillary 3
                     ),
                     num_columns = 3, order = TRUE) &
  theme_bw(base_size = 5.5) & pretty_umap & coord_equal() &
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(), 
        plot.title = element_text(size = 6, hjust = 0.5, face = "bold",
                                  margin = margin(0,0,0.4,0, "mm")),
        legend.text = element_text(size = 5),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        legend.box.margin = margin(0,0,0.5,0, "mm"),
        plot.margin = margin(0.5, 0, 0, 0, "mm"))
endo_feats
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ENDOTHELIAL_FeaturePlots.pdf", 
    width = 71.2*0.0393701, height = 64*0.0393701)
endo_feats
dev.off()


# D: Split FeaturePlots with select genes that are harder to see on main FeaturePlot
endo_genes_to_plot <- c("HEY1", "ACKR1", "CA4", "FCN3")
endo_plots <- lapply(1:length(endo_genes_to_plot), function(XX) { 
  if (XX == 1) {
    FeaturePlot_scCustom(endothelial, features = endo_genes_to_plot[XX], order = TRUE, 
                         split.by = "final_CT", pt.size = 1) &
      theme_bw(base_size = 6.5) & pretty_umap & coord_equal() &
      theme(panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 6),
            axis.title = element_blank(),
            axis.line.y.right = element_blank(),
            axis.title.y.right = element_blank(),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_text(face = "bold", size = 5),
            plot.margin = margin(0, 0, 0, 0, "mm"),
            legend.key.width = unit(2, "mm"),
            legend.key.height = unit(2, "mm"),
            legend.box.margin = margin(0,0,0,0, "mm"),
            legend.justification = "left",
            legend.box.background = element_blank(),
            legend.box.just = "top")
  } else {
    FeaturePlot_scCustom(endothelial, features = endo_genes_to_plot[XX], order = TRUE, 
                         split.by = "final_CT", pt.size = 1) &
      theme_bw(base_size = 6.5) & pretty_umap & coord_equal() &
      theme(panel.background = element_blank(),
            axis.line.y.right = element_blank(),
            plot.title = element_blank(),
            axis.title = element_blank(),
            axis.title.y.right = element_blank(),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_text(face = "bold", size = 5),
            plot.margin = margin(0, 0, 0, 0, "mm"),
            legend.key.width = unit(2, "mm"),
            legend.key.height = unit(2, "mm"),
            legend.box.margin = margin(0,0,0,0, "mm"),
            legend.justification = "left",
            legend.box.background = element_blank(),
            legend.box.just = "top")
  }
})

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ENDOTHELIAL_splitFeaturePlots.pdf",
    width = 68.5*0.0393701, height = 65*0.0393701)
cowplot::plot_grid(endo_plots[[1]], endo_plots[[2]], endo_plots[[3]], endo_plots[[4]], ncol = 1, align = "hv")
dev.off()

# E: Spatial example (DimPlot)
endo_colors2 <- list(`Arteriole` = "#e22039",
                     `Capillary` = "#9D93D2",
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
        cols = c(unlist(endo_colors2), `\nSMCs/Pericytes\n(Mesenchymal)` = "grey90")) + 
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

# F: Spatial example (FeaturePlots)
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ENDOTHELIAL_VUILD107MF_FeaturePlot.pdf",
    width = 138*0.0393701, height = 46.5*0.0393701)
FeaturePlot_scCustom(subset(vuild107mf, subset = final_lineage == "Endothelial"), 
                     features = c("PECAM1", "BMPR2", "HEY1", "VEGFA", "ACKR1",
                                  "CCL21", "FABP4", "CA4", "FCN3", "PLVAP"),
                     reduction = "sp", pt.size = 4, raster = TRUE,
                     num_columns = 5, order = TRUE) &
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



## EPITHELIAL ----
# A: Overall DimPlot of cell types (all epithelial)
epi_ct_order <- c(
  # Epithelial - Alveolar
  "AT1", "Transitional AT2", "AT2", "Proliferating AT2", "KRT5-/KRT17+",
  # Epithelial - Airway
  "Basal", "Multiciliated", "Goblet", "RASC", "Secretory", "PNEC", "Proliferating Airway")
epithelial$final_CT <- ordered(epithelial$final_CT, levels = epi_ct_order)
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_EPITHELIAL_DimPlot.pdf",
    width = 86*0.0393701, height = 44*0.0393701)
DimPlot(epithelial, group.by = "final_CT", cols = epi_colors2) +
  theme_bw(base_size = 6.5) +
  theme_black +
  pretty_umap + coord_equal() + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 5),
        plot.margin = margin(0, 0, 0, 0, "mm"), panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "left",
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(3, "mm")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.5)))
dev.off()

# B: FeaturePlot of alveolar vs. airway
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_EPITHELIAL_FeaturePlot.pdf",
    width = 106*0.0393701, height = 40*0.0393701)
FeaturePlot_scCustom(epithelial, 
                     features = c("ITGB6", "AGER", # Alveolar 1
                                  "KRT8", # Several
                                  "KRT5", "FOXJ1", # Airway 1
                                  "MKI67", # Proliferating
                                  "SFTPC",  # Alveolar 2
                                  "SCGB3A2", # Several
                                  "MUC5B", "SCGB1A1"), # Airway 2
                     num_columns = 5, order = TRUE) &
  theme_bw(base_size = 5.5) & pretty_umap & coord_equal() &
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(), 
        plot.title = element_text(size = 5.5, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 5),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.spacing = unit(0, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(0.5, 0, 0, 0, "mm"))
dev.off()

# C: Alveolar DimPlot
alveolar$final_CT <- ordered(alveolar$final_CT, levels = epi_ct_order)
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ALVEOLAR_DimPlot.pdf",
    width = 52*0.0393701, height = 46*0.0393701)
DimPlot(alveolar, group.by = "final_CT", cols = unlist(epi_colors2),
        order = c("AT1", "KRT5-/KRT17+", "Transitional AT2", "Proliferating AT2", "AT2")) +
  theme_bw(base_size = 6.5) +
  theme_black +
  pretty_umap + coord_equal() + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 5),
        plot.margin = margin(0, 0, 0, 0, "mm"), panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.height = unit(0.5, "mm"),
        legend.key.width = unit(0.5, "mm")) +
  NoLegend()
dev.off()

# D: Alveolar feature plot
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_ALVEOLAR_FeaturePlot.pdf",
    width = 127.5*0.0393701, height = 45*0.0393701)
FeaturePlot_scCustom(alveolar, 
                     features = c("EPCAM", "PECAM1", "COL3A1", "PTPRC", # Lineage markers
                                  "AGER", "RTKN2", "COL4A3", "VEGFA", # AT1
                                  "SFTPC", "LAMP3", "PGC", # AT2 & Transitional AT2
                                  "CEACAM6", # AT1, Transitional AT2
                                  "SCGB3A2", # Transitional AT2
                                  "KRT8", # AT2, Transitional AT2, & KRT5-/KRT17+
                                  "DUOX1", # AT1, AT2, Transitional AT2
                                  "KRT5", "KRT17", "COL1A1", "SOX4", "MMP7", # KRT5-/KRT17+
                                  "MKI67"), # Proliferation
                     num_columns = 7, order = TRUE, raster = TRUE) &
  theme_bw(base_size = 5.5) & pretty_umap & coord_equal() &
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(), 
        plot.title = element_text(size = 6, hjust = 0.5, face = "bold",
                                  margin = margin(0,0,1.6,0, "mm")),
        legend.text = element_text(size = 5),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.spacing = unit(0, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(0.5, 0, 0, 0, "mm"))
dev.off()

# E: Airway DimPlot
airway$final_CT <- ordered(airway$final_CT, levels = epi_ct_order)
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_AIRWAY_DimPlot.pdf",
    width = 52*0.0393701, height = 62*0.0393701)
DimPlot(airway, group.by = "final_CT", cols = unlist(epi_colors2)) +
  theme_bw(base_size = 6.5) +
  theme_black +
  pretty_umap + coord_equal() + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 5),
        plot.margin = margin(0, 0, 0, 0, "mm"), panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.height = unit(0.5, "mm"),
        legend.key.width = unit(0.5, "mm")) +
  NoLegend()
dev.off()

# F: Airway feature plot
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_AIRWAY_FeaturePlot.pdf",
    width = 127.2*0.0393701, height = 62*0.0393701)
FeaturePlot_scCustom(airway, 
                     features = c("EPCAM", "PECAM1", "COL3A1", "PTPRC", # Lineage markers
                                  "SOX4", # General Airway
                                  "MMP7", # General Airway (+); Multiciliated (-)
                                  "KRT5", "KRT6A", "TP63", "KRT17", # Basal
                                  "FOXJ1", "C20orf85", "TP73", # Multiciliated (+); Secretory & RASC (-)
                                  "MUC5B", "MUC5AC", # Goblet (+); Secretory & RASC (-)
                                  "SCGB1A1", # Secretory (+) & RASC (-)
                                  "SCGB3A2", "CEACAM6", # RASC
                                  "CALCA", "CHGB", # PNEC
                                  "MKI67"), # Proliferating
                     num_columns = 7, order = TRUE, raster = TRUE) &
  theme_bw(base_size = 5.5) & pretty_umap & coord_equal() &
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(), 
        plot.title = element_text(size = 6, hjust = 0.5, face = "bold",
                                  margin = margin(0,0,2.5,0, "mm")),
        legend.text = element_text(size = 5),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.spacing = unit(0, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(0, 0, 1.5, 0, "mm"))
dev.off()


## IMMUNE ----
# A: Overall DimPlot of cell types (all immune)
imm_ct_order <- c(
  # Immune - Lymphoid
  "B cells", "Proliferating B cells", "NK/NKT", "Proliferating NK/NKT", "Tregs", 
  "CD4+ T-cells", "CD8+ T-cells", "Proliferating T-cells", "Plasma", "pDCs",
  # Immune - Myeloid
  "Basophils", "cDCs", "Migratory DCs", "Langerhans cells", "Mast", 
  "Alveolar Macrophages", "Interstitial Macrophages", "SPP1+ Macrophages",
  "Macrophages - IFN-activated", "Monocytes/MDMs", "Neutrophils", "Proliferating Myeloid")
immune$final_CT <- ordered(immune$final_CT, levels = imm_ct_order)
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_IMMUNE_DimPlot.pdf",
    width = 86*0.0393701, height = 44*0.0393701)
DimPlot(immune, group.by = "final_CT",
        cols = unlist(imm_colors2)) +
  theme_bw(base_size = 6.5) +
  theme_black +
  pretty_umap + coord_equal() + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 5),
        plot.margin = margin(0, 0, 0, 0, "mm"), panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "left",
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.box.spacing = unit(0, "mm"),
        legend.key.height = unit(1.1, "mm"),
        legend.key.width = unit(1.1, "mm")) +
  guides(color = guide_legend(ncol = 1))
dev.off()

# B: FeaturePlot of myeloid vs. lymphoid
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_IMMUNE_FeaturePlot.pdf",
    width = 88*0.0393701, height = 45*0.0393701)
FeaturePlot_scCustom(immune, 
                     features = c("MS4A1", "GNLY", "CD3E","JCHAIN", # Lymphoid
                                  "MKI67", # Proliferating
                                  "CPA3", "LYZ", "FCER1G", "CD68", "CD14", 
                                  "MARCO", "PPARG", "S100A8", # Myeloid
                                  "CD1A", "CD1C"), # cDC markers
                     num_columns = 5, order = TRUE) &
  theme_bw(base_size = 5.5) & pretty_umap & coord_equal() &
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(), 
        plot.title = element_text(size = 5.5, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 5),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.spacing = unit(0, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(0.5, 0, 0, 0, "mm"))
dev.off()

# C: Lymphoid DimPlot
lymphoid$final_CT <- ordered(lymphoid$final_CT, levels = imm_ct_order)
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_LYMPHOID_DimPlot.pdf",
    width = 55.3*0.0393701, height = 52.5*0.0393701)
DimPlot(lymphoid, group.by = "final_CT",
        cols = unlist(imm_colors2)) +
  theme_bw(base_size = 6.5) +
  theme_black +
  pretty_umap + coord_equal() + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 5),
        plot.margin = margin(0, 0, 0, 0, "mm"), panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.height = unit(0.5, "mm"),
        legend.key.width = unit(0.5, "mm")) +
  NoLegend()
dev.off()

# D: Lymphoid feature plot
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_LYMPHOID_FeaturePlot.pdf",
    width = 123*0.0393701, height = 48*0.0393701)
FeaturePlot_scCustom(lymphoid, 
                     features = c("JCHAIN", "ITM2C", "PIM2", # pDCs/Plasma
                                  "LILRA4", # pDCs
                                  "BANK1", "MS4A1", "TNFRSF13C", # B cells
                                  "CD3E", "CD3D", # T-cells
                                  "CD8A", "CD8B", # CD8+ T-cells
                                  "CD4", # CD4+ T-cells
                                  "FOXP3", "CTLA4", "RTKN2", # Tregs
                                  "GNLY", "NKG7", # NK/NKT
                                  "MKI67"), # Proliferation
                     num_columns = 6, order = TRUE) &
  theme_bw(base_size = 5.5) & pretty_umap & coord_equal() &
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(), 
        plot.title = element_text(size = 6, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 5),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.spacing = unit(0, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(0.5, 0, 0, 0, "mm"))
dev.off()

# E: Myeloid DimPlot
myeloid$final_CT <- ordered(myeloid$final_CT, levels = imm_ct_order)
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_MYELOID_DimPlot.pdf",
    width = 55.3*0.0393701, height = 52*0.0393701)
DimPlot(myeloid, group.by = "final_CT",
        cols = unlist(imm_colors2),
        order = "Langerhans cells") +
  # DimPlot(myeloid, group.by = "final_CT",
  #         cols = unlist(imm_colors2),
  #         order = "Langerhans cells") +
  theme_bw(base_size = 6.5) +
  theme_black +
  pretty_umap + coord_equal() + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 5),
        plot.margin = margin(0, 0, 0, 0, "mm"), panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.height = unit(0.5, "mm"),
        legend.key.width = unit(0.5, "mm")) +
  NoLegend()
dev.off()

# F: Myeloid feature plot
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_MYELOID_FeaturePlot.pdf",
    width = 123*0.0393701, height = 100*0.0393701)
FeaturePlot_scCustom(myeloid, 
                     features = c("LYZ", # Myeloid overall
                                  "CPA3", "TPSAB1", # Mast/Basophils
                                  "KIT", # Mast
                                  "LAMP3", "CCR7", 
                                  
                                  "CCL22", # Migratory DCs
                                  "CD1C", "FCER1A", # cDCs
                                  "CD1A", "FCGBP", # Langerhans DCs
                                  "CD68", # Macrophages/Monocytes
                                  
                                  "MRC1", "CD14", # Macrophages/Monocytes
                                  "FCN1", # Mono/MDMs
                                  "SLC25A37", # Neutrophils
                                  "S100A8", "S100A9", # Neutrophils & Mono/MDMS
                                  
                                  "AGER", "RTKN2", # Placeholder genes
                                  "MARCO", # Alveolar/SPP1+ Macrophages
                                  "FABP4", # Alveolar Macrophages
                                  "SPP1", # SPP1+ Macrophages
                                  "S100A12", # Neutrophils & Mono/MDMS
                                  
                                  "SFTPC", "FOXJ1", # Placeholder genes
                                  "PPARG", # Alveolar/SPP1+ Macrophages
                                  "OAS3", "IFIT1", # IFN-activated
                                  "MKI67"), # Proliferating
                     num_columns = 6, order = TRUE) &
  theme_bw(base_size = 5.5) & pretty_umap & coord_equal() &
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(), 
        plot.title = element_text(size = 6, hjust = 0.5, face = "bold",
                                  margin = margin(0,0,1.23,0, "mm")),
        legend.text = element_text(size = 5),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.box.spacing = unit(0, "mm"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm"),
        plot.margin = margin(0, 0, 0, 0, "mm"))
dev.off()



## MESENCHYMAL ----
# A: DimPlot of cell types
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_MESENCHYMAL_DimPlot.pdf", 
    width = 82*0.0393701, height = 44*0.0393701)
DimPlot(mesenchymal, group.by = "final_CT", cols = unlist(mes_colors2),
        # label.color = c("white", "black", rep("white", 3), rep("black", 3), "white"),
        label.color = c(rep("black", 6), "white", "white", "black"),
        order = c("Myofibroblasts", "Proliferating FBs", "Inflammatory FBs",
                  "Alveolar FBs", "Activated Fibrotic FBs", "Subpleural FBs",
                  "Adventitial FBs", "SMCs/Pericytes", "Mesothelial")) + 
  theme_bw(base_size = 6.5) +
  theme_black +
  pretty_umap + coord_equal() + 
  theme(plot.title = element_blank(), axis.title = element_text(size = 5),
        plot.margin = margin(0, 0, 0, 0, "mm"), panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.key.height = unit(0.8, "mm"),
        legend.key.width = unit(0.8, "mm"))
dev.off()

# B: Proportion of Fibroblasts by group; Alveolar FBs highest in Unaffected
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_MESENCHYMAL_prop_fibroblasts.pdf", 
    width = 30*0.0393701, height = 29*0.0393701)
xenium@meta.data %>%
  filter(final_lineage == "Mesenchymal",
         !(final_CT %in% c("SMCs/Pericytes", "Mesothelial"))) %>%
  mutate(sample_affect = case_when(sample_affect == "Less Affected" ~ "Less",
                                   sample_affect == "More Affected" ~ "More",
                                   TRUE ~ sample_affect),
         sample_affect = ordered(sample_affect, levels = c("Unaffected", "Less", "More"))) %>%
  group_by(final_CT) %>%
  mutate(num_cells_type = length(final_CT)) %>%
  ungroup() %>%
  ggplot(aes(x = sample_affect, fill = reorder(final_CT, num_cells_type))) +
  geom_bar(position = position_fill()) +
  scale_fill_manual(values = unlist(mes_colors2)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 6.5) +
  theme_black +
  labs(y = "Proportion of Fibroblasts") +
  theme(axis.title.x = element_blank(), legend.title = element_blank(),
        axis.title.y = element_text(size = 6)) +
  NoLegend()
dev.off()


# C: Overall FeaturePlots of marker genes
pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_MESENCHYMAL_FeaturePlots.pdf", 
    width = 125*0.0393701, height = 80*0.0393701)
FeaturePlot_scCustom(mesenchymal, 
                     features = c("COL1A1", "FN1", # General FB markers
                                  "CTHRC1", "FAP", "POSTN", # Activated Fibrotic FBs
                                  "PI16", "MFAP5", # Adventitial FBs
                                  # Alveolar FBs - no specific gene markers; use spatial
                                  "CCL2", "PTGDS", "SOD2", # Inflammatory FBs
                                  "WNT5A", # Myofibroblasts
                                  "PLIN2", "HAS1", "WT1", # Subpleural FBs
                                  "MKI67", # Proliferating FBs (1)
                                  "ACTA2", "CSPG4", # SMCs/Pericytes
                                  "MSLN", "KRT8", # Mesothelial
                                  "TOP2A"), # Proliferating FBs (2)
                     num_columns = 5, order = TRUE) &
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

# D: Split FeaturePlots with select genes that are harder to see on main FeaturePlot
mes_genes_to_plot <- c("WNT5A", "HAS1", "POSTN", "CCL2")
mes_plots <- lapply(1:length(mes_genes_to_plot), function(XX) { 
  if (XX == 1) {
    FeaturePlot_scCustom(mesenchymal, features = mes_genes_to_plot[XX], order = TRUE, 
                         split.by = "final_CT", pt.size = 1.5) &
      theme_bw(base_size = 6.5) & pretty_umap & coord_equal() &
      theme(panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 6),
            axis.title = element_text(color = "white"),
            axis.line.y.right = element_blank(),
            axis.title.y.right = element_blank(),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_text(face = "bold", size = 5),
            plot.margin = margin(0.1, 0, 0, 0.5, "mm"),
            legend.key.width = unit(2, "mm"),
            legend.key.height = unit(2, "mm"),
            legend.justification = "left")
  } else {
    FeaturePlot_scCustom(mesenchymal, features = mes_genes_to_plot[XX], order = TRUE, 
                         split.by = "final_CT", pt.size = 1.5) &
      theme_bw(base_size = 6.5) & pretty_umap & coord_equal() &
      theme(panel.background = element_blank(),
            axis.line.y.right = element_blank(),
            plot.title = element_blank(),
            axis.title = element_text(color = "white"),
            axis.title.y.right = element_blank(),
            axis.line = element_line(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_text(face = "bold", size = 5),
            plot.margin = margin(0.1, 0, 0, 0.5, "mm"),
            legend.key.width = unit(2, "mm"),
            legend.key.height = unit(2, "mm"),
            legend.justification = "left")
  }
})

pdf("/scratch/avannan/MANUSCRIPT_figures/SuppFig_MESENCHYMAL_splitFeaturePlots.pdf", 
    width = 179*0.0393701, height = 75*0.0393701)
cowplot::plot_grid(mes_plots[[1]], mes_plots[[2]], mes_plots[[3]], mes_plots[[4]], ncol = 1, align = "hv")
dev.off()






















# a <- FeaturePlot_scCustom(mesenchymal_split[["Mesothelial"]], features = "MSLN", split.by = "final_CT") & 
#   theme_bw(base_size = 6.5) & pretty_umap & coord_equal() &
#   theme(axis.line = element_blank(), axis.title.y.right = element_blank(), 
#         plot.title = element_text(hjust = 0.5),
#         panel.grid = element_blank())
# FeaturePlot_scCustom(mesenchymal_split[["Activated Fibrotic FBs"]], features = c("CTHRC1", "FAP", "POSTN"), split.by = "final_CT") & 
#   theme_bw(base_size = 6.5) & pretty_umap & coord_equal() &
#   theme(axis.line = element_blank(), axis.title.y.right = element_blank(), 
#         plot.title = element_text(hjust = 0.5),
#         panel.grid = element_blank())
# ggarrange(a,b)
# FeaturePlot_scCustom(mesenchymal_split[["Subpleural FBs"]], features = c("PLIN2", "HAS1"), split.by = "final_CT") & 
#   theme_bw(base_size = 6.5) & pretty_umap & coord_equal() &
#   theme(axis.line = element_blank(), axis.title.y.right = element_blank(), 
#         plot.title = element_text(hjust = 0.5),
#         panel.grid = element_blank())


unique(xenium$patient) %in% "TILD270"

FeaturePlot_scCustom(mesenchymal, features = c("ACTA2", "CSPG4", "MFAP5", "WNT5A", "MSLN",
                                               "PLIN2", "CTHRC1", "FAP", "CCL2", "PTGDS", "TOP2A"), 
                     split.by = "final_CT", order = TRUE) &
  theme_bw(base_size = 6.5) & pretty_umap & coord_equal() & 
  theme(axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y.left = element_blank(),
        axis.title.y.right = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())


library(scCustomize)
test <- FeaturePlot_scCustom(mesenchymal, features = "COL1A1", split.by = "final_CT", order = TRUE, combine = FALSE)
test2 <- lapply(1:(length(test)-3), function(XX) { test[[XX]] + NoLegend() })
test3 <- list(test2, test[length(test)-2], length(test)-1, length(test))
ggarrange(test3, common.legend = TRUE)
ggarrange(test)

FeaturePlot_scCustom(mesenchymal, features = "COL1A1", split.by = "final_CT", order = TRUE, common.legend = TRUE)

# c("COL1A1", "ACTA2", "CSPG4", "MFAP5", "WNT5A", "MSLN", "PLIN2", "CTHRC1", "SOD2", "PTGDS", "YAP1")
FeaturePlot(mesenchymal, features = c("COL1A1", "ACTA2"), order = TRUE,
            split.by = "final_CT", keep.scale = all) &
  theme_bw(base_size = 6.5) & pretty_umap & coord_equal() & theme(axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y.left = element_blank(),
                                                                  axis.title.y.right = element_text(hjust = 0.5))

FeaturePlot(mesenchymal, features = c("COL1A1", "ACTA2", "CSPG4", "MFAP5", "WNT5A", "MSLN", "PLIN2", "CTHRC1", "SOD2", "PTGDS", "YAP1"), order = TRUE) &
  pretty_umap & coord_equal() & theme(axis.line = element_blank(), axis.title = element_blank()) & NoLegend()



mesenchymal$label_CT <- ""

mesenchymal$label_CT <- unname(unlist(test[mesenchymal$final_CT]))

test <- list(`Activated Fibrotic FBs` = 1,
             `Adventitial FBs` = 2,
             `Alveolar FBs` = 3,
             `Inflammatory FBs` = 4,
             `Mesothelial` = 5,
             `Myofibroblasts` = 6,
             `SMCs/Pericytes` = 7,
             `Subpleural FBs` = 8,
             `Proliferating FBs` = 9)
DimPlot(mesenchymal, group.by = "final_CT", cols = unname(unlist(mes_colors)), label = TRUE, label.box = TRUE, repel = TRUE) + 
  pretty_umap + coord_equal() + theme(plot.title = element_blank())
FeaturePlot(mesenchymal, features = c("ACTA2", "CSPG4", "CTHRC1", "FAP", "MSLN", "HAS1", "PLIN2", "MKI67"), 
            order = TRUE) & coord_equal() & pretty_umap

DimPlot(airway, group.by = "final_CT", cols = unlist(air_colors), label = TRUE, repel = TRUE) + pretty_umap + coord_equal()
FeaturePlot(airway, features = c("KRT5", "MUC5B", "C20orf85", "CALCA", "SCGB3A2", "SCGB1A1", "MKI67")) & coord_equal() & pretty_umap



vuild96mf <- subset(airway, subset = sample == "VUILD96MF")
DimPlot(vuild96mf, reduction = "sp", group.by = "final_CT", cols = air_colors) + coord_equal() + pretty_umap + 
  theme(plot.title = element_blank())
DimPlot(vuild96mf, reduction = "sp", group.by = "final_CT", cols = other_air_colors) + coord_equal() + pretty_umap + 
  theme(plot.title = element_blank())

DimPlot(vuild96mf, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
  theme(plot.title = element_blank()) +
  scale_y_continuous(limits = c(-5700, -3700)) +
  scale_x_continuous(limits = c(-3800, -950))

vuild110 <- subset(airway, subset = sample == "VUILD110")
DimPlot(vuild110, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
  theme(plot.title = element_blank()) +
  scale_y_continuous(limits = c(-17750, -15800), expand = c(0, 0)) +
  scale_x_continuous(limits = c(2100, 2800), expand = c(0, 0)) +
  pretty_umap


vuhd049 <- subset(airway, subset = sample == "VUHD049")
DimPlot(vuhd049, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
  theme(plot.title = element_blank())

# vuild106 <- subset(airway, subset = sample == "VUILD106")
# DimPlot(vuild106, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
#   theme(plot.title = element_blank())

# vuild58 <- subset(airway, subset = sample == "VUILD58")
# DimPlot(vuild58, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
#   theme(plot.title = element_blank())

vuild107mf <- subset(airway, subset = sample == "VUILD107MF")
DimPlot(vuild107mf, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
  theme(plot.title = element_blank()) +
  scale_y_continuous(limits = c(-11450, -10400), expand = c(0, 0)) +
  scale_x_continuous(limits = c(4550, 6200), expand = c(0, 0))

# tild117mfb <- subset(airway, subset = sample == "TILD117MFB")
# DimPlot(tild117mfb, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
#   theme(plot.title = element_blank()) +
#   scale_y_continuous(limits = c(-5900, -5000), expand = c(0, 0)) +
#   scale_x_continuous(limits = c(26000, 27500), expand = c(0, 0))

vuild104mf <- subset(airway, subset = sample == "VUILD104MF")
DimPlot(vuild104mf, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
  theme(plot.title = element_blank()) +
  scale_y_continuous(limits = c(-9400, -8500), expand = c(0, 0)) +
  scale_x_continuous(limits = c(15200, 16700), expand = c(0, 0))
# scale_y_continuous(limits = c(-11490, -10800), expand = c(0, 0)) +
# scale_x_continuous(limits = c(15540, 16200), expand = c(0, 0))
# scale_y_continuous(limits = c(-11500, -9800), expand = c(0, 0)) +
# scale_x_continuous(limits = c(14300, 16200), expand = c(0, 0))


sort(table(airway$final_CT, airway$sample)["RASC",])

# tild117mf <- subset(airway, subset = sample == "TILD117MF")
# DimPlot(tild117mf, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
#   theme(plot.title = element_blank()) +
#   scale_y_continuous(limits = c(-17480, -16300), expand = c(0, 0)) +
#   scale_x_continuous(limits = c(17500, 18150), expand = c(0, 0))
other_air_colors <- c("#c30039",
                      "#b2d350",
                      "#498aff",
                      "#175f26",
                      "#ff316a",
                      "#714270",
                      "#fead15")
tild299lf <- subset(xenium, subset = sample == "TILD299LF" & (final_sublineage == "Airway"))
DimPlot(tild299lf, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.8) + coord_equal() +
  theme(plot.title = element_blank()) +
  scale_y_continuous(limits = c(-15250, -14250), expand = c(0, 0)) +
  scale_x_continuous(limits = c(25450, 26500), expand = c(0, 0)) +
  theme(legend.position = c(0.6, 0.75)) +
  pretty_umap
FeaturePlot(tild299lf, reduction = "sp", features = c("TP63", "MUC5B", "C20orf85", "SCGB3A2", "SCGB1A1", "CALCA"),
            pt.size = 0.1, order = TRUE, ncol = 3) &
  pretty_umap &
  theme(axis.line = element_blank(), axis.title = element_blank()) &
  coord_equal() &
  scale_y_continuous(limits = c(-15250, -14250), expand = c(0, 0)) &
  scale_x_continuous(limits = c(25450, 26500), expand = c(0, 0))


# tild049lf <- subset(airway, subset = sample == "TILD049LF")
# DimPlot(tild049lf, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.2) + coord_equal() +
#   theme(plot.title = element_blank())
# vuhd095 <- subset(airway, subset = sample == "VUHD095")
# DimPlot(vuhd095, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
#   theme(plot.title = element_blank())
# vuhd116a <- subset(airway, subset = sample == "VUHD116A")
# DimPlot(vuhd116a, reduction = "sp", group.by = "final_CT", cols = other_air_colors, pt.size = 0.5) + coord_equal() +
#   theme(plot.title = element_blank())




DimPlot(alveolar, reduction = "sp", group.by = "sample", label = TRUE) + coord_equal() + NoLegend()

DimPlot(endothelial, reduction = "sp", group.by = "sample", label = TRUE) + coord_equal() + NoLegend()
DimPlot(endothelial, reduction = "sp", group.by = "final_CT", cols = c("red", "blue", "green", "orange")) + coord_equal()
DimPlot(mesenchymal, reduction = "sp", group.by = "final_CT") + coord_equal()




other_mes_colors <- c("#c810e8",
                      "#ff57b2",
                      "grey90",
                      "#bed00d",
                      "#24cdff",
                      "#ab6600",
                      "#00a754")
vuild106 <- subset(xenium, subset = sample == "VUILD106" &
                     final_CT %in% c("Adventitial FBs", "Alveolar FBs", "Activated Fibrotic FBs",
                                     "Mesothelial", "Subpleural FBs", "Myofibroblasts", "SMCs/Pericytes"))
DimPlot(vuild106, reduction = "sp", group.by = "final_CT", cols = other_mes_colors, pt.size = 0.25) + coord_equal() +
  theme(plot.title = element_blank()) + pretty_umap
FeaturePlot(vuild106, features = c("COL1A1", "ACTA2", "CSPG4", "MFAP5", "WNT5A", "MSLN", "PLIN2", "CTHRC1"), reduction = "sp",
            pt.size = 0.0001, order = TRUE, ncol = 4) & 
  coord_equal() & pretty_umap & theme(axis.line = element_blank(), axis.title = element_blank())


# scale_y_continuous(limits = c(-15250, -14250), expand = c(0, 0)) +
# scale_x_continuous(limits = c(25450, 26500), expand = c(0, 0)) +
# theme(legend.position = c(0.6, 0.75)) +
# pretty_umap

tmp <- xenium@meta.data %>%
  filter(sample == "VUHD113") %>%
  select(cell_id, final_CT) %>%
  dplyr::rename(group = final_CT)
rownames(tmp) <- NULL
write.csv(tmp, "/scratch/avannan/VUHD113_celltypes.csv")
other_alv_colors <- c("#ff37c0",
                      "#01bc8a",
                      "black",
                      "#647c00",
                      "#c8c0f2",
                      "#ff8941",
                      "#93002c")
vuhd113_a <- subset(xenium, subset = sample == "VUHD113" & final_CT %in%
                      c("AT1", "AT2", "Alveolar Macrophages", "Interstitial Macrophages", "Capillary", "Alveolar FBs", "Neutrophils"))
DimPlot(vuhd113_a, reduction = "sp", group.by = "final_CT", cols = other_alv_colors, pt.size = 0.25) + coord_equal() +
  theme(plot.title = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(-7650, -4000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(13500, 16400)) +
  pretty_umap
FeaturePlot(vuhd113_a, reduction = "sp", features = c("AGER", "SFTPC", "CA4"), pt.size = 0.25) & coord_equal() &
  scale_y_continuous(expand = c(0, 0), limits = c(-7650, -4000)) &
  scale_x_continuous(expand = c(0, 0), limits = c(13500, 16400)) &
  pretty_umap & theme(axis.line = element_blank(), axis.title = element_blank())



DimPlot(xenium, group.by = "sample", label = TRUE, reduction = "sp") + coord_equal() + NoLegend()
DimPlot(xenium, group.by = "final_lineage", reduction = "sp", pt.size = 1, cols = c("")) + coord_equal()
vuild110 <- subset(xenium, subset = sample == "VUILD110" & final_CT %in%
                     c("AT1", "AT2", "Alveolar Macrophages", "Interstitial Macrophages", "Capillary", "Alveolar FBs", "Neutrophils"))
DimPlot(vuild110, reduction = "sp", group.by = "final_CT", cols = other_alv_colors, pt.size = 0.25) + coord_equal() +
  theme(plot.title = element_blank()) +
  # scale_y_continuous(expand = c(0, 0), limits = c(-7650, -4000)) +
  scale_x_continuous(expand = c(0, 0), limits = c(1700, 7000))
pretty_umap
FeaturePlot(vuild110, reduction = "sp", features = c("AGER", "SFTPC", "CA4"), pt.size = 0.25, order = TRUE) & coord_equal() &
  pretty_umap & theme(axis.line = element_blank(), axis.title = element_blank())

vuhd113_b <- subset(xenium, subset = sample == "VUHD113")
vuhd113_b_range <- subset(vuhd113_b, subset = adj_y_centroid >= -5650 & adj_y_centroid <= -5275 &
                            adj_x_centroid >= 13691 & adj_x_centroid <= 13850)
DimPlot(vuhd113_b, reduction = "sp", group.by = "final_CT", cols = unlist(color_list[vuhd113_b_celltypes]), pt.size = 3) & coord_equal() &
  theme(plot.title = element_blank()) &
  scale_y_continuous(expand = c(0, 0), limits = c(-5650, -5275)) &
  scale_x_continuous(expand = c(0, 0), limits = c(13691, 13850)) &
  pretty_umap
FeaturePlot(vuhd113_b, reduction = "sp", features = c("EPCAM", "AGER", "RTKN2", "SFTPC", "LAMP3",
                                                      "PECAM1", "CA4", "FCN3",
                                                      "PTPRC", "S100A8"), ncol = 5, pt.size = 1) & coord_equal() &
  scale_y_continuous(expand = c(0, 0), limits = c(-5650, -5275)) &
  scale_x_continuous(expand = c(0, 0), limits = c(13691, 13850)) &
  pretty_umap & theme(axis.line = element_blank(), axis.title = element_blank())

FeaturePlot(vuhd113_b_range, reduction = "sp", features = c("AGER", "RTKN2"), blend = TRUE, pt.size = 1) & coord_equal()
FeaturePlot(vuhd113_b_range, reduction = "sp", features = c("SFTPC", "LAMP3"), blend = TRUE, pt.size = 1) & coord_equal()




# vuild106_a <- subset(xenium, subset = sample == "VUILD106" & (final_lineage == "Endothelial" | final_CT == "SMCs/Pericytes"))
# DimPlot(vuild106_a, reduction = "sp", group.by = "final_CT", cols = other_endo_colors, pt.size = 0.25) + coord_equal() +
#   theme(plot.title = element_blank())
# vuild106_b <- subset(xenium, subset = sample == "VUILD106" & final_lineage == "Endothelial")
# DimPlot(vuild106_b, reduction = "sp", group.by = "final_CT", cols = c("#e22039",
#                                                                       "#c8c0f2",
#                                                                       "#0272db",
#                                                                       "#ccc100"), pt.size = 0.25) + coord_equal() +
#   theme(plot.title = element_blank()) & pretty_umap
# FeaturePlot(vuild106, features = c("HEY1", "ACKR1", "CCL21", "CA4", "FCN3", "APLN", "APLNR", "IL7R"), reduction = "sp",
#             pt.size = 0.01, order = TRUE) & 
#   coord_equal() & pretty_umap & theme(axis.line = element_blank(), axis.title = element_blank())



https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094
library(googlesheets4)
gs4_deauth()
path_sheet <- gs4_get("https://docs.google.com/spreadsheets/d/1SQLg5dvS69YwEaXyUevn0j04MLGRIKkNI9zC0kECvVE/edit?gid=1329421094#gid=1329421094")
path_sheet <- read_sheet(path_sheet, sheet = 5, skip = 0)

path_df <- path_sheet %>%
  filter(Feature %in% c("Pathology_score", "Percent pathology", "Adjustment_factor", "Adjusted_path_score")) %>%
  t() %>%
  .[-1, ] %>%
  as.data.frame() %>%
  mutate(Pathology_score = as.numeric(V1),
         Percent_pathology = as.numeric(V2),
         Adjustment_factor = as.numeric(V3),
         Adjusted_path_score = as.numeric(V4)) %>%
  select(Pathology_score, Percent_pathology, Adjustment_factor, Adjusted_path_score) %>%
  rownames_to_column(var = "sample")

sample_colors <- distinctColorPalette(45)

path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, tma) %>% unique()) %>%
  filter(tma != "TMA5") %>%
  ggplot(aes(x = Pathology_score, y = Percent_pathology)) +
  geom_point(aes(color = sample)) +
  geom_smooth(method = "lm", alpha = 0.2) +
  scale_color_manual(values = sample_colors)
path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, tma) %>% unique()) %>%
  filter(tma != "TMA5") %>%
  ggplot(aes(x = Pathology_score, y = Adjusted_path_score)) +
  geom_point(aes(color = sample)) +
  geom_smooth(method = "lm", alpha = 0.2) +
  scale_color_manual(values = sample_colors)
path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, tma) %>% unique()) %>%
  filter(tma != "TMA5") %>%
  ggplot(aes(x = Percent_pathology, y = Adjusted_path_score)) +
  geom_point(aes(color = sample)) +
  geom_smooth(method = "lm", alpha = 0.2) +
  scale_color_manual(values = sample_colors)

path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, tma) %>% unique()) %>%
  mutate(sample_type = ifelse(is.na(sample_type), "MF", sample_type)) %>%
  ggplot(aes(x = Percent_pathology, fill = sample_type)) +
  geom_histogram() +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0))

path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, tma) %>% unique()) %>%
  ggplot(aes(x = Percent_pathology, fill = sample)) +
  geom_histogram() +
  scale_fill_manual(values = sample_colors) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0))



path_df <- path_sheet %>%
  filter(Feature %in% c("Pathology_score", "Percent pathology", "Adjustment_factor", "Adjusted_path_score")) %>%
  t() %>%
  .[-1, ] %>%
  as.data.frame() %>%
  mutate(Pathology_score = as.numeric(V1),
         Percent_pathology = as.numeric(V2),
         Adjustment_factor = as.numeric(V3),
         Adjusted_path_score = as.numeric(V4)) %>%
  select(Pathology_score, Percent_pathology, Adjustment_factor, Adjusted_path_score) %>%
  rownames_to_column(var = "sample")

path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, patient, tma) %>% unique()) %>%
  mutate(sample_affect = case_when(sample_type == "Unaffected" ~ "Unaffected",
                                   Percent_pathology > 75 ~ "More Affected",
                                   Percent_pathology <= 75 ~ "Less Affected")) %>%
  mutate(patient = ifelse(is.na(patient), "VUILD090", patient)) %>%
  
  
  path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, tma) %>% unique()) %>%
  filter(sample_type != "Unaffected") %>%
  filter(Percent_pathology > 75) %>%
  dim()

path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, patient, tma) %>% unique()) %>%
  mutate(sample_affect = case_when(sample_type == "Unaffected" ~ "Unaffected",
                                   Percent_pathology > 75 ~ "More Affected",
                                   Percent_pathology <= 75 ~ "Less Affected")) %>%
  mutate(patient = ifelse(is.na(patient), "VUILD090", patient)) %>%
  ggplot(aes(x = sample_affect, fill = patient)) +
  geom_bar() +
  scale_fill_manual(values = distinctColorPalette(45)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0))


xenium@meta.data %>%
  
  
  test <- path_df %>%
  left_join(xenium@meta.data %>% select(sample_type, sample, patient, tma) %>% unique()) %>%
  mutate(sample_affect = case_when(sample_type == "Unaffected" ~ "Unaffected",
                                   Percent_pathology > 75 ~ "More Affected",
                                   Percent_pathology <= 75 ~ "Less Affected")) %>%
  mutate(patient = ifelse(is.na(patient), "VUILD090", patient))
table(test$patient, test$sample_affect)

7, 13, 25 (50 as cutoff)
7, 18, 20 (70 as cutoff)

# Loupe broswer ideas
# Endothelial genes
# TLSs? Granulomas? I don't think these are present on sample of interest













