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

