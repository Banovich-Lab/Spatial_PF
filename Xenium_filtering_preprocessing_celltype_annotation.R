############################################
# Cell Type Annotation and Quality Filtering
# Author: Annika Vannan (avannan@tgen.org)
# Date: 12/18/2023
# Description: Initial quality filter and
#              cell type annotation for 
#              Xenium spatial PF project
############################################

## SETTING ENVIRONMENT ----
# Set library paths so they are consistent in RStudio Server and command line R
.libPaths(c("/home/avannan/R/x86_64-pc-linux-gnu-library/4.2",
            "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))


library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(tidyverse)
library(gplots)
library(tibble)
library(ggpubr)
library(ggrepel)


# Set seed
set.seed(0317)
work_dir <- "/scratch/avannan"
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)
filter <- dplyr::filter
select <- dplyr::select
pretty_umap <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     axis.title = element_text(hjust = 1))


## SET LISTS OF LINEAGE MARKERS TO EXPLORE ----
epithelial_features <- c("EGFR", "DUOX1", "NKX2-1", "AGER", "RTKN2", "NAPSA", "PGC", "SFTA2",
                         "SFTPC", "SFTPD", "KRT14", "KRT15", "KRT5", "KRT6A", "S100A2",
                         "TP63", "KRT17", "AGR3", "C20orf85", "FOXJ1", "GCLM", "DMBT1",
                         "EPCAM", "KRT18", "MGST1", "MMP7", "FOXI1", "MUC5B", "SCGB1A1", 
                         "SCGB3A2", "WFDC2", "ATF3", "KRT8", "SOX9", "SPINK1", "GKN2",
                         "MMP10", "SOX2", "CDH26", "TP73", "CFTR", "HES1", "PKM", "SOX4",
                         "NUCB2", "RNASE1", "SAA2", "AKR1C1", "AKR1C2", "BPIFA1", "CEACAM5",
                         "ERN2", "FCGBP", "GSR", "LTF", "CCNA1", "ICAM1", "ITGB1")
endothelial_features <- c("APLN", "CA4", "HEY1", "BMPR2", "CD34", "EPAS1", "FCN3", 
                          "GNG11", "PECAM1", "APLNR", "COL15A1", "PLVAP", "ACKR1",
                          "POSTN", "CLDN5", "RAMP2", "ZEB1", "HAS1", "KDR", "CDKN2A")
immune_features <- c("PPARG", "BANK1", "CD19", "CD79A", "LTB", "MS4A1", "TNFRSF13C", 
                     "CD86", "GZMB", "HLA-DRA", "CCR7", "CXCR4", "PTPRC", "TCL1A", 
                     "CD69", "CD4", "CD8A", "CD8B", "CD2", "CD28", "CD3D", "CD3E", 
                     "CD3G", "FOXP3", "GZMK", "TRAC", "ITM2C", "CD27", "CCL5", "LCK", 
                     "FABP4", "MARCO", "MCEMP1", "SPP1", "FCN1", "S100A12", "S100A8", 
                     "S100A9", "CCL22", "ITGAM", "NFKB1",  "IFIT2", "FGFBP2", "GNLY", 
                     "KLRB1", "KLRC1", "NKG7", "LILRA4", "BCL2", "CD79B", "CXCR5", 
                     "CXCL9", "GPR183", "HLA-DQA1", "KLRG1", "BCL2L11", "CD52", 
                     "SLC1A3", "TNFRSF9", "CTLA4",  "IL2RA", "LAG3", "PDCD1", "PDIA6", 
                     "PIM2", "IL7R", "LEF1", "FASLG", "HAVCR2", "ISG20", "CPA3", 
                     "KIT", "TPSAB1", "C1QC", "CD68", "MS4A7", "AIF1", "CD14", 
                     "FCGR3A", "FCER1G", "SLC25A37", "CD247", "GZMA", "IRF7")
mesenchymal_features <- c("MFAP5", "PI16", "SFRP2", "ELN", "FAP", "AXL", "LGR6", 
                          "COL1A1", "COL1A2", "COL3A1", "DCN", "FN1", "HAS2", 
                          "LUM", "MEG3", "SPARCL1", "CTHRC1")


## SUBDIRECTORY NAMES ----
# Get subdirectory names that correspond with sample IDs
work_dir <- "/scratch/avannan/xenium_files/xenium_ouput/"
id_list <- c(
  # TMA1
  VUHD116A = "output-XETG00048__0003817__VUHD116A__20230308__003730",
  VUHD116B = "output-XETG00048__0003817__VUHD116B__20230308__003731",
  VUILD102LF = "output-XETG00048__0003817__VUILD102LF__20230308__003731",
  VUILD102MF = "output-XETG00048__0003817__VUILD102MF__20230308__003730",
  VUILD107MF = "output-XETG00048__0003817__VUILD107MF__20230308__003731",
  VUILD96LF = "output-XETG00048__0003817__VUILD96LF__20230308__003730",
  VUILD96MF = "output-XETG00048__0003817__VUILD96MF__20230308__003730",
  # # TMA2
  VUHD069 = "output-XETG00048__0003789__VUHD069__20230308__003731",
  VUHD095 = "output-XETG00048__0003789__VUHD095__20230308__003731",
  VUHD113 = "output-XETG00048__0003789__VUHD113__20230308__003731",
  VUILD104LF = "output-XETG00048__0003789__VUILD104LF__20230308__003731",
  VUILD104MFVUILD48LFVUILD105LF = "output-XETG00048__0003789__VUILD104MFVUILD48LFVUILD105LF__20230308__003731",
  VUILD105MF = "output-XETG00048__0003789__VUILD105MF__20230308__003731",
  VUILD48MF = "output-XETG00048__0003789__VUILD48MF__20230308__003731",
  # TMA3
  THD0008 = "output-XETG00048__0003392__THD0008__20230313__191400",
  VUILD106 = "output-XETG00048__0003392__VUILD106__20230313__191400",
  VUILD110 = "output-XETG00048__0003392__VUILD110__20230313__191400",
  VUILD115 = "output-XETG00048__0003392__VUILD115__20230313__191400",
  # TMA4
  THD0011 = "output-XETG00048__0003400__THD0011__20230313__191400",
  TILD117LF = "output-XETG00048__0003400__TILD117LF__20230313__191400",
  TILD117MF = "output-XETG00048__0003400__TILD117MF__20230313__191400",
  TILD175 = "output-XETG00048__0003400__TILD175__20230313__191400",
  VUILD78LF = "output-XETG00048__0003400__VUILD78LF__20230313__191400",
  VUILD78MF = "output-XETG00048__0003400__VUILD78MF__20230313__191400",
  VUILD91LF = "output-XETG00048__0003400__VUILD91LF__20230313__191400",
  VUILD91MF = "output-XETG00048__0003400__VUILD91MF__20230313__191400")

# Get subdirectory names for obtaining file paths
subdirs <- unname(id_list)

# Get TMA identity
tmas <- unlist(lapply(str_split(id_list, "__"), function(XX) { XX[2]}))
names(tmas) <- names(id_list)
tmas <- gsub("0003817", "TMA1", tmas)
tmas <- gsub("0003789", "TMA2", tmas)
tmas <- gsub("0003392", "TMA3", tmas)
tmas <- gsub("0003400", "TMA4", tmas)

# Get run identity
run_ids <- unlist(lapply(str_split(id_list, "__"), function(XX) { XX[4]}))
names(run_ids) <- names(id_list)
run_ids <- gsub("20230308", "Run1", run_ids)
run_ids <- gsub("20230308", "Run1", run_ids)
run_ids <- gsub("20230313", "Run2", run_ids)
run_ids <- gsub("20230313", "Run2", run_ids)


## LOAD IN TRANSCRIPT AND METADATA FILES ----
# Get count and metadata files
all_files <- list.files(file.path(work_dir, subdirs), full.names = TRUE)
h5_files <- all_files[grep(".h5", all_files)]
transcript_files <- all_files[grep("transcripts.csv.gz", all_files)]
meta_files <- all_files[grep("cells.csv.gz", all_files)]

# Get sample IDs
sample_ids <- names(id_list)

# Read in files
counts <- lapply(h5_files, Read10X_h5)

transcripts <- lapply(transcript_files, function(XX) {
  read_csv(XX, col_types = c(transcript_id = "c", cell_id = "c")) })

metadata <- lapply(meta_files, function(XX) {
  tmp_meta <- read.delim(XX, sep = ",", colClasses = c(cell_id = "character"))
  rownames(tmp_meta) <- tmp_meta$cell_id
  tmp_meta })

# Rename files in lists
sample_ids <- unlist(lapply(str_split(meta_files, "__"), function(XX) { XX[3] }))
names(counts) <- sample_ids
names(transcripts) <- sample_ids
names(metadata) <- sample_ids

# Get transcripts that only overlap the nucleus and create cell x gene matrix
# Also count the number of blanks per cell
all_transcripts <- list()
nuc_transcripts <- list()
updated_metadata <- list()
for (sm in sample_ids) {
  message(paste("Getting nuclei counts for sample", sm))
  
  # Filter out low quality transcripts 
  all_transcripts[[sm]] <- transcripts[[sm]][transcripts[[sm]]$qv > 20, ]
  
  # Find transcripts that overlap a nucleus
  nuc_transcripts[[sm]] <- transcripts[[sm]][transcripts[[sm]]$overlaps_nucleus == "1", ]
  
  # Create cell x gene dataframe
  nuc_transcripts[[sm]] <- as.data.frame(table(nuc_transcripts[[sm]]$cell_id, 
                                               nuc_transcripts[[sm]]$feature_name))
  names(nuc_transcripts[[sm]]) <- c("cell_id", "feature_name", "Count")
  nuc_transcripts[[sm]] <- nuc_transcripts[[sm]] %>% 
    pivot_wider(names_from = "feature_name", values_from = "Count")
  
  # Get blanks count per nucleus
  blank_nuc_ids <- nuc_transcripts[[sm]]$cell_id
  blank_nuc_mat <- nuc_transcripts[[sm]][, grep("BLANK", 
                                                colnames(nuc_transcripts[[sm]]))]
  blank_nuc_counts <- as.data.frame(rowSums(blank_nuc_mat))
  blank_nuc_counts$cell_id <- blank_nuc_ids
  
  # Remove negative controls and convert to cell x gene matrix
  nuc_transcripts[[sm]] <- nuc_transcripts[[sm]][, grep("NegControl", 
                                                        colnames(nuc_transcripts[[sm]]), 
                                                        invert = TRUE)]
  nuc_transcripts[[sm]] <- nuc_transcripts[[sm]][, grep("BLANK", 
                                                        colnames(nuc_transcripts[[sm]]), 
                                                        invert = TRUE)]
  keep_cells <- nuc_transcripts[[sm]]$cell_id
  nuc_transcripts[[sm]] <- as.data.frame(nuc_transcripts[[sm]])
  rownames(nuc_transcripts[[sm]]) <- keep_cells
  nuc_transcripts[[sm]] <- nuc_transcripts[[sm]][, -1]
  nuc_transcripts[[sm]] <- as.matrix(t(nuc_transcripts[[sm]]))
  
  # Subset nuclear metadata to "cells" with transcripts that overlap nuclei
  updated_metadata[[sm]] <- metadata[[sm]][metadata[[sm]]$cell_id %in% keep_cells, ]
  
  # Add blank counts to metadata
  updated_metadata[[sm]] <- full_join(updated_metadata[[sm]], blank_nuc_counts,
                                      by = "cell_id")
  updated_metadata[[sm]] <- updated_metadata[[sm]] %>%
    rename(num.blank = `rowSums(blank_nuc_mat)`)
  rownames(updated_metadata[[sm]]) <- updated_metadata[[sm]]$cell_id
}


## COORDINATE ADJUSTMENT LIST ----
# Create coordinate adjustment list
coord_adjust_list <- list(
  # TMA1 
  VUHD095 = c(0, 400),
  VUHD069 = c(6300, 0),
  VUILD104LF = c(0, 3400),
  VUHD113 = c(3000, 3700),
  VUILD105MF = c(2800, 0),
  VUILD48MF = c(-200, 7800),
  VUILD104MFVUILD48LFVUILD105LF = c(0, 6000),
  # TMA2
  VUHD116A = c(0, 0),
  VUILD102LF = c(3800, -800),
  VUILD96MF = c(-4000, 3000),
  VUILD102MF= c(0, 3400),
  VUHD116B = c(4000, 4000),
  VUILD96LF = c(0, 7000),
  VUILD107MF = c(4000, 7800),
  # TMA3
  VUILD115 = c(0, 0),
  VUILD110 = c(6000, 500),
  VUILD106 = c(0, 7000),
  THD0008 = c(6000, 7500),
  # TMA4
  THD0011 = c(0, 0),
  TILD117MF = c(4000, 0),
  VUILD78MF = c(-4000, 3600),
  TILD117LF = c(-500, 3600),
  VUILD91MF = c(4000, 4200),
  VUILD91LF = c(-3500, 7500),
  VUILD78LF = c(-500, 7000),
  TILD175 = c(4000, 7200))


## CREATE SEURAT OBJECTS ----
obj_list <- list()
obj_list <- sapply(sample_ids, function(XX) {
  # Create a Seurat object containing the RNA adata
  sobj <- CreateSeuratObject(counts = nuc_transcripts[[XX]], 
                             assay = "RNA")
  
  # Add metadata
  sobj <- AddMetaData(sobj, metadata = updated_metadata[[XX]])
  sobj$sample <- XX
  sobj$tma <- tmas[[XX]]
  sobj$run <- run_ids[[XX]]
  
  # Calculate percent blank
  sobj$percent.blank <- sobj$num.blank/(sobj$nCount_RNA + sobj$num.blank)*100
  
  # Remove cells with 0 nCount_RNA
  sobj <- subset(sobj, subset = nCount_RNA != 0)
  
  # Rename cells to add sample ID as prefix
  if (XX != "VUILD104MFVUILD48LFVUILD105LF") {
    sobj <- RenameCells(sobj, add.cell.id = XX)
  }
  
  # Adjust coordinates
  if (tmas[[XX]] == "TMA1") {
    sobj$adj_x_centroid <- sobj$x_centroid + coord_adjust_list[[XX]][1]
    sobj$adj_y_centroid <- (sobj$y_centroid + coord_adjust_list[[XX]][2])*-1
  } else if (tmas[[XX]] == "TMA2") {
    sobj$adj_x_centroid <- (sobj$x_centroid + coord_adjust_list[[XX]][1])+10500
    sobj$adj_y_centroid <- (sobj$y_centroid + coord_adjust_list[[XX]][2])*-1
  } else if (tmas[[XX]] == "TMA3") {
    sobj$adj_x_centroid <- (sobj$x_centroid + coord_adjust_list[[XX]][1])-4600
    sobj$adj_y_centroid <- ((sobj$y_centroid + coord_adjust_list[[XX]][2])*-1)-13000
  } else if (tmas[[XX]] == "TMA4") {
    sobj$adj_x_centroid <- (sobj$x_centroid + coord_adjust_list[[XX]][1])+12900
    sobj$adj_y_centroid <- ((sobj$y_centroid + coord_adjust_list[[XX]][2])*-1)-15100
  }
  
  # Add spatial coordinates as dimension reduction objects
  position_xy <- cbind(sobj$adj_x_centroid, sobj$adj_y_centroid)
  row.names(position_xy) <- row.names(sobj@meta.data)
  colnames(position_xy) <- c("SP_1", "SP_2")
  sobj[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                       assay = DefaultAssay(sobj))
  obj_list[[XX]] <- sobj
})


sample_list[["VUILD78MF"]]@meta.data %>%
  ggplot(aes(x = x_centroid, y = y_centroid)) +
  geom_point()


## FIX MERGED SAMPLES AND MERGE ALL XENIUM DATA ----
# Split merged samples
obj_list[["VUILD104MF"]] <- subset(obj_list[["VUILD104MFVUILD48LFVUILD105LF"]],
                                   subset = y_centroid < 3800)
obj_list[["VUILD48LF"]] <- subset(obj_list[["VUILD104MFVUILD48LFVUILD105LF"]],
                                  subset = x_centroid < 3950 & y_centroid > 3800)
obj_list[["VUILD105LF"]] <- subset(obj_list[["VUILD104MFVUILD48LFVUILD105LF"]],
                                   subset = x_centroid > 3950 & y_centroid > 3800)

# Adjust coordinates for merged samples
obj_list[["VUILD104MF"]]$adj_x_centroid <- obj_list[["VUILD104MF"]]$adj_x_centroid+3000
obj_list[["VUILD104MF"]]$adj_y_centroid <- obj_list[["VUILD104MF"]]$adj_y_centroid-2200
obj_list[["VUILD48LF"]]$adj_x_centroid <- obj_list[["VUILD48LF"]]$adj_x_centroid+5700
obj_list[["VUILD48LF"]]$adj_y_centroid <- obj_list[["VUILD48LF"]]$adj_y_centroid+1800
obj_list[["VUILD105LF"]]$adj_x_centroid <- obj_list[["VUILD105LF"]]$adj_x_centroid+2300
obj_list[["VUILD105LF"]]$adj_y_centroid <- obj_list[["VUILD105LF"]]$adj_y_centroid+6200

# Change sample names and rename cells
obj_list[["VUILD104MF"]]$sample <- "VUILD104MF"
obj_list[["VUILD48LF"]]$sample <- "VUILD48LF"
obj_list[["VUILD105LF"]]$sample <- "VUILD105LF"

### ERROR HERE
obj_list[["VUILD104MF"]] <- RenameCells(obj_list[["VUILD104MF"]], 
                                        add.cell.id = paste0("VUILD104MF", "_"))
obj_list[["VUILD48LF"]] <- RenameCells(obj_list[["VUILD48LF"]], 
                                       add.cell.id = paste0("VUILD48LF", "_"))
obj_list[["VUILD105LF"]] <- RenameCells(obj_list[["VUILD105LF"]], 
                                        add.cell.id = paste0("VUILD105LF", "_"))

# Remove merged samples
obj_list <- obj_list[-which(names(obj_list) == "VUILD104MFVUILD48LFVUILD105LF")]

# Get sample IDs again
sample_ids <- names(obj_list)

# Label patients
for (sm in names(obj_list)) {
  obj_list[[sm]]$patient <- sm
  obj_list[[sm]]$patient <- gsub("A", "", obj_list[[sm]]$patient)
  obj_list[[sm]]$patient <- gsub("B", "", obj_list[[sm]]$patient)
  obj_list[[sm]]$patient <- gsub("LF", "", obj_list[[sm]]$patient)
  obj_list[[sm]]$patient <- gsub("MF", "", obj_list[[sm]]$patient)
}

# Label which samples are replicates of the MERFISH TMAs
for (sm in names(obj_list)) {
  obj_list[[sm]]$merfish_replicate <- ""
  if (unique(obj_list[[sm]]$sample) == "THD0008" | 
      unique(obj_list[[sm]]$tma) == "TMA4") {
    obj_list[[sm]]$merfish_replicate <- "Xenium Only"
  } else {
    obj_list[[sm]]$merfish_replicate <- "MERFISH + Xenium"
  }
}

# Merge objects (cannot do spatial DimPlots for this)
merged_spatial_unfiltered <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

# Add spatial dimension reduction object separately
position_xy <- cbind(merged_spatial_unfiltered$adj_x_centroid,
                     merged_spatial_unfiltered$adj_y_centroid)
row.names(position_xy) <- row.names(merged_spatial_unfiltered@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
merged_spatial_unfiltered[["sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(merged_spatial_unfiltered))

# View TMAs!
DimPlot(merged_spatial_unfiltered, reduction = "sp", group.by = "sample", label = TRUE)
DimPlot(merged_spatial_unfiltered, reduction = "sp", group.by = "tma", label = TRUE)
DimPlot(merged_spatial_unfiltered, reduction = "sp", group.by = "run", label = TRUE)

# saveRDS(merged_spatial_unfiltered, "/scratch/avannan/xenium_files/xenium_spatial_nuclei_only_unfiltered.rds")
merged_spatial_unfiltered <- readRDS("/scratch/avannan/xenium_files/xenium_spatial_nuclei_only_unfiltered.rds")


## ADD CELL-LEVEL COUNT DATA ----
# Get sample IDs
sample_ids <- names(id_list)

cell_obj_list <- list()
cell_obj_list <- sapply(sample_ids, function(XX) {
  message(paste("Creating cell Seurat object for sample", XX))
  
  # Create a Seurat object containing the RNA cell information
  sobj <- CreateSeuratObject(counts = counts[[XX]]$`Gene Expression`,
                             assay = "RNA")
  rownames(metadata[[XX]]) <- metadata[[XX]]$cell_id
  sobj <- AddMetaData(sobj, metadata = metadata[[XX]])
  
  # Rename cells to add sample ID as prefix
  if (XX != "VUILD104MFVUILD48LFVUILD105LF") {
    sobj <- RenameCells(sobj, add.cell.id = XX)
  }
  
  cell_obj_list[[XX]] <- sobj
})

# Split merged samples
cell_obj_list[["VUILD104MF"]] <- subset(cell_obj_list[["VUILD104MFVUILD48LFVUILD105LF"]],
                                        subset = y_centroid < 3800)
cell_obj_list[["VUILD48LF"]] <- subset(cell_obj_list[["VUILD104MFVUILD48LFVUILD105LF"]],
                                       subset = x_centroid < 3950 & y_centroid > 3800)
cell_obj_list[["VUILD105LF"]] <- subset(cell_obj_list[["VUILD104MFVUILD48LFVUILD105LF"]],
                                        subset = x_centroid > 3950 & y_centroid > 3800)

# Change sample names and rename cells
cell_obj_list[["VUILD104MF"]]$sample <- "VUILD104MF"
cell_obj_list[["VUILD48LF"]]$sample <- "VUILD48LF"
cell_obj_list[["VUILD105LF"]]$sample <- "VUILD105LF"
cell_obj_list[["VUILD104MF"]] <- RenameCells(cell_obj_list[["VUILD104MF"]],
                                             add.cell.id = paste0("VUILD104MF", "_"))
cell_obj_list[["VUILD48LF"]] <- RenameCells(cell_obj_list[["VUILD48LF"]],
                                            add.cell.id = paste0("VUILD104MF", "_"))
cell_obj_list[["VUILD105LF"]] <- RenameCells(cell_obj_list[["VUILD105LF"]],
                                             add.cell.id = paste0("VUILD104MF", "_"))

# Merge cell information
cell_merged <- merge(cell_obj_list[[1]], y = cell_obj_list[2:length(cell_obj_list)])

# Add cell information to nuclei object
cell_count_matrix <- cell_merged@assays$RNA@counts
keep_cells <- colnames(merged_spatial_unfiltered)
cell_count_matrix <- cell_count_matrix[, keep_cells]
merged_spatial_unfiltered[["cell_RNA"]] <- CreateAssayObject(counts = cell_count_matrix)


## QUALITY CONTROL ----
# Number of cells per sample before filtering
summary(as.factor(merged_spatial_unfiltered$sample))
merged_spatial_unfiltered@meta.data %>%
  ggplot(aes(y = sample, fill = run)) +
  geom_bar()

# Percent.blank
merged_spatial_unfiltered@meta.data %>%
  ggplot(aes(x = percent.blank, fill = sample)) +
  geom_histogram(bins = 50, show.legend = FALSE, color = "black") +
  theme_classic() +
  theme(title = element_text(color = "black"), 
        axis.text = element_text(color = "black")) +
  facet_wrap(~sample, scales = "free")
# nCount_RNA
merged_spatial@meta.data %>%
  ggplot(aes(x = nCount_RNA, fill = sample)) +
  geom_histogram(bins = 50, show.legend = FALSE, color = "black") +
  theme_classic() +
  theme(title = element_text(color = "black"), 
        axis.text = element_text(color = "black")) +
  facet_wrap(~sample, scales = "free")
# nucleus_area
merged_spatial_unfiltered@meta.data %>%
  ggplot(aes(x = nucleus_area, fill = sample)) +
  geom_histogram(bins = 50, show.legend = FALSE, color = "black") +
  theme_classic() +
  theme(title = element_text(color = "black"), 
        axis.text = element_text(color = "black")) +
  facet_wrap(~sample, scales = "free")

# Create function for better VlnPlot
BetterVlnPlot <- function(data, features, ylim = NA){
  VlnPlot(data, pt.size = 0, features = features, 
          group.by = "sample", y.max = ylim) + labs(x = "") + NoLegend()
}

BetterVlnPlot(merged_spatial_unfiltered, features = "percent.blank")
BetterVlnPlot(merged_spatial_unfiltered, features = "nCount_RNA")
BetterVlnPlot(merged_spatial_unfiltered, features = "nFeature_RNA")
BetterVlnPlot(merged_spatial_unfiltered, features = "nucleus_area")

# nCount_RNA vs. percent.blank
smoothScatter(merged_spatial_unfiltered@meta.data$percent.blank,
              log(merged_spatial_unfiltered@meta.data$nCount_RNA),
              cex = 0.5, pch = 16)
abline(v = 4, h = log(12), lty = "dashed", col = "black")
text(5, 5, col = "black", adj = c(0, -.1),
     "nCount_RNA >= 12 & percent.blank <= 4")

# nFeature_RNA vs. percent.blank
smoothScatter(merged_spatial_unfiltered@meta.data$percent.blank,
              log(merged_spatial_unfiltered@meta.data$nFeature_RNA),
              cex = 0.5, pch = 16)
abline(v = 4, h = log(10), lty = "dashed", col = "black")
text(5, 4, col = "black", adj = c(0, -.1),
     "nFeature_RNA >= 10 & percent.blank <= 4")

# nCount_RNA vs. nFeature_RNA
smoothScatter(log(merged_spatial_unfiltered$nCount_RNA),
              log(merged_spatial_unfiltered$nFeature_RNA),
              cex = 0.5, pch = 16)
abline(v = log(12), h = log(10), lty = "dashed", col = "black")
text(0.3, 4.6, col = "black", adj = c(0, -.1),
     "nCount_RNA >= 12 & nFeature_RNA >= 10")

# nCount RNA vs. nucleus_area
smoothScatter(merged_spatial_unfiltered$nucleus_area,
              log(merged_spatial_unfiltered$nCount_RNA),
              cex = 0.5, pch = 16)
abline(v = c(6, 80), h = log(12), lty = "dashed", col = "black")
text(120, 0.7, col = "black", adj = c(0, -.1),
     "nCount_RNA >= 12 & nucleus_area between 6-80")

# nFeature RNA vs. nucleus_area
smoothScatter(merged_spatial_unfiltered$nucleus_area,
              log(merged_spatial_unfiltered$nFeature_RNA),
              cex = 0.5, pch = 16)
abline(v = c(6, 80), h = log(10), lty = "dashed", col = "black")
text(120, 0.4, col = "black", adj = c(0, -.1),
     "nFeature_RNA >= 10 & & nucleus_area between 6-80")

min(merged_spatial_unfiltered$nucleus_area)
max(merged_spatial_unfiltered$nucleus_area)

# Filter merged and individual data
merged_spatial2 <- subset(merged_spatial_unfiltered,
                          subset = nCount_RNA >= 12 & nFeature_RNA >= 10 &
                            percent.blank <= 5 & 
                            nucleus_area >= 6 & nucleus_area <= 80)

# Number of nuclei before and after filtering
bf_cells <- table(merged_spatial_unfiltered$sample)
aft_cells <- table(merged_spatial$sample)
diff_cells <- bf_cells - aft_cells
prop_kept_cells <- round(aft_cells/bf_cells*100, 2)
prop_kept_cells

# DimPlots of before and after for each sample
DimPlotCompare <- function(sm){
  bf_cells <- ncol(subset(merged_spatial_unfiltered, subset = sample == sm))
  a <- DimPlot(subset(merged_spatial_unfiltered, subset = sample == sm),
               reduction = "sp") + NoLegend() +
    labs(title = paste0(sm, ", Unfiltered, ", bf_cells, " nuclei"))
  
  aft_cells <- ncol(subset(merged_spatial, subset = sample == sm))
  b <- DimPlot(subset(merged_spatial, subset = sample == sm),
               reduction = "sp") + NoLegend() +
    labs(title = paste0(sm, ", Filtered, ", aft_cells, " nuclei"))
  ggarrange(a,b)
}

# TMA 1
DimPlotCompare("VUHD116A")
DimPlotCompare("VUHD116B")
DimPlotCompare("VUILD102LF")
DimPlotCompare("VUILD102MF")
DimPlotCompare("VUILD96LF")
DimPlotCompare("VUILD96MF")
DimPlotCompare("VUILD107LF")
DimPlotCompare("VUILD107MF")
# TMA 2
DimPlotCompare("VUHD069")
DimPlotCompare("VUHD095")
DimPlotCompare("VUHD113")
DimPlotCompare("VUILD104LF")
DimPlotCompare("VUILD104MF")
DimPlotCompare("VUILD48LF")
DimPlotCompare("VUILD48MF")
DimPlotCompare("VUILD105LF")
DimPlotCompare("VUILD105MF")
# TMA 3
DimPlotCompare("THD0008")
DimPlotCompare("VUILD106")
DimPlotCompare("VUILD115")
DimPlotCompare("VUILD110")
# TMA 4
DimPlotCompare("THD0011")
DimPlotCompare("TILD117LF")
DimPlotCompare("TILD117MF")
DimPlotCompare("TILD175")
DimPlotCompare("VUILD78LF")
DimPlotCompare("VUILD78MF")
DimPlotCompare("VUILD91LF")
DimPlotCompare("VUILD91MF")

# Label disease
merged_spatial$sample_type <- ""
merged_spatial$sample_type[merged_spatial$sample %in% 
                             c("VUHD116A", "VUHD116B", "VUHD069", "VUHD095", 
                               "VUHD113", "THD0008", "THD0011")] <- "Unaffected"
merged_spatial$sample_type[merged_spatial$sample %in% 
                             c("VUILD102LF", "VUILD107LF", "VUILD96LF", 
                               "VUILD104LF", "VUILD105LF", "VUILD48LF", 
                               "TILD117LF", "VUILD78LF", "VUILD91LF")] <- "LF"
merged_spatial$sample_type[merged_spatial$sample %in% 
                             c("VUILD102MF", "VUILD107MF", "VUILD96MF", 
                               "VUILD104MF", "VUILD105MF", "VUILD48MF", 
                               "TILD117MF", "VUILD78MF", "VUILD91MF")] <- "MF"
merged_spatial$sample_type[merged_spatial$sample %in% 
                             c("VUILD106", "VUILD110", "VUILD115", "TILD175")] <- "ILD"


## UMAP AND CLUSTERING ----
# Normalize and scale data and run PCA
DefaultAssay(merged_spatial) <- "RNA"
merged_spatial <- NormalizeData(merged_spatial)
merged_spatial <- ScaleData(merged_spatial, features = rownames(merged_spatial))
merged_spatial <- RunPCA(merged_spatial, features = rownames(merged_spatial))

# Function for getting PCs
get_pcs <- function(seurat_obj, reduc = "pca", var = 0.05) {
  # Determine percent of variation associated with each PC
  pct <- seurat_obj[[reduc]]@stdev / sum(seurat_obj[[reduc]]@stdev)*100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % 
  # variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  # Last point where change of % of variation is more than 0.05%
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > var), 
              decreasing = T)[1] + 1
  c(co1, co2)
}

# Get PCs and make elbow plot
npcs <- min(get_pcs(merged_spatial))
ElbowPlot(merged_spatial, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 22

# Find neighbors, cluster and UMAP
merged_spatial <- FindNeighbors(merged_spatial, dims = 1:npcs) # 10 min
merged_spatial <- RunUMAP(test, dims = 1:npcs, n.neighbors = 100) # 40 min
merged_spatial <- FindClusters(merged_spatial, resolution = c(0.5)) # 15 min per resolution
DimPlot(merged_spatial, group.by = "RNA_snn_res.0.3", label = TRUE) + NoLegend()

# better_merged_spatial <- RunHarmony(merged_spatial, group.by.vars = "sample", max.iter.harmony = 50, theta = 1, dims.use = 1:npcs)
# better_merged_spatial <- FindNeighbors(better_merged_spatial, reduction = "harmony", dims = 1:npcs)
# better_merged_spatial <- FindClusters(better_merged_spatial, resolution = c(0.2, 0.3))
# better_merged_spatial <- RunUMAP(better_merged_spatial, dims = 1:npcs, reduction = "harmony")
# DimPlot(better_merged_spatial, group.by = "sample")
# saveRDS(better_merged_spatial, "/scratch/avannan/better_xenium_spatial_maybe_04142023.rds")
# better_merged_spatial <- readRDS("/scratch/avannan/better_xenium_spatial_maybe_04142023.rds")

DimPlot(merged_spatial, group.by = "sample")
DimPlot(merged_spatial, group.by = "sample", ncol = 4, split.by = "sample") + NoLegend()
DimPlot(merged_spatial, group.by = "tma")
DimPlot(merged_spatial, group.by = "tma", split.by = "tma")
DimPlot(merged_spatial, group.by = "run", split.by = "run")
DimPlot(merged_spatial, group.by = "merfish_replicate")
DimPlot(merged_spatial, group.by = "sample_type", split.by = "sample_type")

# Cluster quality
merged_spatial <- SetIdent(merged_spatial, value = "RNA_snn_res.0.3")
VlnPlot(merged_spatial, features = c("nCount_RNA", "nFeature_RNA", "percent.blank"), group.by = "RNA_snn_res.0.3", pt.size = 0)
DimPlot(merged_spatial, group.by = "RNA_snn_res.0.3", split.by = "RNA_snn_res.0.3", ncol = 5) + NoLegend()

# Find Markers
merged_spatial <- SetIdent(merged_spatial, value = "RNA_snn_res.0.3")
top_markers03 <- FindAllMarkers(merged_spatial, only.pos = TRUE)
top_markers03 <- top_markers03 %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

merged_spatial <- SetIdent(merged_spatial, value = "RNA_snn_res.0.4")
top_markers04 <- FindAllMarkers(merged_spatial, only.pos = TRUE)
top_markers04 <- top_markers04 %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

merged_spatial <- SetIdent(merged_spatial, value = "RNA_snn_res.0.5")
top_markers05 <- FindAllMarkers(merged_spatial, only.pos = TRUE)
top_markers05 <- top_markers05 %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")


DimPlot(merged_spatial, group.by = "RNA_snn_res.0.4", label = TRUE)
merged_spatial <- SetIdent(merged_spatial, value = "RNA_snn_res.0.4")
FeaturePlot(merged_spatial, features = c("EPCAM", "PECAM1", "PTPRC", "LUM", "DCN", "ACTA2"))

merged_spatial <- SetIdent(merged_spatial, value = "RNA_snn_res.0.5")
DotPlot(merged_spatial, features = c("EPCAM", "PECAM1", "PTPRC", "LUM", "DCN", 
                                     "ACTA2", "ELN", "KIT", "MS4A1", "CCL21", 
                                     "KDR", "CCL5", "MSLN", "PDGFRB", "VIM",
                                     "ITGB6", "COL3A1"))
VlnPlot(merged_spatial, features = c("EPCAM", "PECAM1", "PTPRC", "LUM", "DCN", 
                                     "ACTA2", "ELN", "CCL21", "KDR", "PDGFRA", "PDGFRB",
                                     "VIM"),
        pt.size = 0)

subset22 <- subset(merged_spatial, subset = RNA_snn_res.0.5 == "12")
subset22_counts <- rowMeans(subset22@assays$RNA@counts)
sort(subset22_counts)


## EPITHELIAL ----
# Subset any clusters that have EPCAM expression and EPCAM as a marker gene, even if there are other markers
epithelial <- subset(merged_spatial, subset = RNA_snn_res.0.5 %in% c(3, 8, 10, 17, 21))

DefaultAssay(epithelial) <- "RNA"
epithelial <- NormalizeData(epithelial)
epithelial <- ScaleData(epithelial, features = rownames(epithelial))
epithelial <- RunPCA(epithelial, features = rownames(epithelial))

# Get PCs and make elbow plot
npcs <- min(get_pcs(epithelial))
ElbowPlot(epithelial, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 13

# Find neighbors, cluster and UMAP
epithelial <- FindNeighbors(epithelial, dims = 1:npcs)
epithelial <- RunUMAP(epithelial, dims = 1:npcs, n.neighbors = 100)
epithelial <- FindClusters(epithelial, resolution = c(0.2))
DimPlot(epithelial, group.by = "RNA_snn_res.0.2", label = TRUE)

# Look at epithelial cell type markers, as well as general lineage markers
DotPlot(epithelial, features = c("EPCAM", "SFTPC", "LAMP3",
                                 "AGER", "RTKN2", "CEACAM6", "KRT8",
                                 "SCGB3A2", "SCGB1A1", "MUC5B", "MUC5AC",
                                 "KRT17", "KRT5", "TP63", "TP73", "FOXJ1",
                                 "CCL21", "CD3E", "TRAC", "GZMA"))
DotPlot(epithelial, features = c("EPCAM", "PECAM1", "PTPRC", "LUM", "ACTA2", "COL1A1"))

# Look at other lineage markers
# Cluster 5 appears to be a mixed of different cell types, with primarily immune cells
FeaturePlot(epithelial, features = c("EPCAM", "PTPRC"), split.by = "RNA_snn_res.0.2", blend = TRUE)
DotPlot(epithelial, features = c("EPCAM", "PECAM1", "PTPRC"))

epithelial <- SetIdent(epithelial, value = "RNA_snn_res.0.2")
epi_markers <- FindAllMarkers(epithelial, only.pos = TRUE) %>%
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# 5 has immune cells
FeaturePlot(epithelial, features = c("CD3E", "MS4A7", "CCL21"), split.by = "RNA_snn_res.0.2")
epithelial <- FindSubCluster(epithelial, cluster = "5", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub5")
DimPlot(epithelial, group.by = "sub5", label = TRUE)
epithelial <- SetIdent(epithelial, value = "sub5")
immune_epi_markers <- FindMarkers(epithelial, ident.1 = "5_0", ident.2 = "5_1")

# Start with general labels
epithelial$CT_firstpass <- ""
epithelial$CT_firstpass[epithelial$RNA_snn_res.0.2 == 5] <- "Immune"
epithelial$CT_firstpass[epithelial$RNA_snn_res.0.2 %in% c(1, 3)] <- "Airway Epithelium"
epithelial$CT_firstpass[epithelial$RNA_snn_res.0.2 %in% c(0, 2, 4)] <- "Alveolar Epithelium"

epithelial <- FindSubCluster(epithelial, cluster = "1", graph.name = "RNA_snn", resolution = 0.18)
subset1 <- subset(epithelial, subset = RNA_snn_res.0.2 == "1")
DimPlot(subset1, group.by = "sub.cluster", label = TRUE)

DotPlot(subset1, group.by = "sub.cluster",
        features = c("EPCAM", "KRT8", "SCGB3A2", "SCGB1A1", "MUC5B", "MUC5AC", 
                     "KRT17", "KRT5", "TP63", "TP73", "FOXJ1", "SOX4", "SOX9",
                     "MKI67"))

FeaturePlot(subset1, features = c("SFTPC", "KRT8"), split.by = "sub.cluster" )# Lots of KRT8 throughout
FeaturePlot(subset1, features = c("TP73", "FOXJ1"), split.by = "sub.cluster") # 1_0 and 1_4 = Ciliated
FeaturePlot(subset1, features = c("MUC5B", "MUC5AC", "JCHAIN"), split.by = "sub.cluster") # 1_3 = MUC5B+, plasma, and mixed cells

# 1_2 does not have strong distinguishing markers

FeaturePlot(subset1, features = c("SCGB3A2", "SCGB1A1"), split.by = "sub.cluster") # 1_1 is SCGB secretory
FeaturePlot(subset1, features = c("SCGB3A2", "SCGB1A1"), split.by = "sub.cluster", blend = TRUE) # Contains both SCGB3A2+/SCGB1A1+ and SCGB3A2+ cells

subset1 <- SetIdent(subset1, value = "sub.cluster")
epi1_markers <- FindAllMarkers(subset1, only.pos = TRUE)
epi1_markers <- epi1_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# Based on marker genes, 1_1 shares some similarities with alveolar epithelium
# 1_0 and 1_4 = Ciliated
# 1_1 is SCGB secretory (need to split)
# 1_2 does not have strong distinguishing markers; has MUC5B, FOXJ1, SCGB, etc. expression
# 1_3 = MUC5B+, plasma, and mixed cells (need to split)

# Look at SCGB markers - nuclear and cellular
scgb1a1_nuc_expr <- names(which((alv_air@assays$RNA@counts["SCGB1A1", ] > 0), arr.ind = TRUE))
scgb3a2_nuc_expr <- names(which((alv_air@assays$RNA@counts["SCGB3A2", ] > 0), arr.ind = TRUE))

alv_air$scgb_nuc_expr <- "SCGB3A2-/SCGB1A1-"
alv_air$scgb_nuc_expr[colnames(alv_air) %in% scgb1a1_nuc_expr] <- "SCGB3A2-/SCGB1A1+"
alv_air$scgb_nuc_expr[colnames(alv_air) %in% scgb3a2_nuc_expr] <- "SCGB3A2+/SCGB1A1-"
alv_air$scgb_nuc_expr[colnames(alv_air) %in% scgb1a1_nuc_expr & colnames(alv_air) %in% scgb3a2_nuc_expr] <- "SCGB3A2+/SCGB1A1+"
DimPlot(alv_air, group.by = "scgb_nuc_expr")

scgb1a1_cell_expr <- names(which((alv_air@assays$cell_RNA@counts["SCGB1A1", ] > 0), arr.ind = TRUE))
scgb3a2_cell_expr <- names(which((alv_air@assays$cell_RNA@counts["SCGB3A2", ] > 0), arr.ind = TRUE))
alv_air$scgb_cell_expr <- "SCGB3A2-/SCGB1A1-"
alv_air$scgb_cell_expr[colnames(alv_air) %in% scgb1a1_cell_expr] <- "SCGB3A2-/SCGB1A1+"
alv_air$scgb_cell_expr[colnames(alv_air) %in% scgb3a2_cell_expr] <- "SCGB3A2+/SCGB1A1-"
alv_air$scgb_cell_expr[colnames(alv_air) %in% scgb1a1_cell_expr & colnames(alv_air) %in% scgb3a2_cell_expr] <- "SCGB3A2+/SCGB1A1+"
DimPlot(alv_air, group.by = "scgb_cell_expr", cols = c("grey", "red", "blue", "green"), split.by = "sub2")
DimPlot(alv_air, group.by = "scgb_nuc_expr", cols = c("grey", "red", "blue", "green"), split.by = "sub2")

table(alv_air$sub2, alv_air$scgb_nuc_expr)
table(alv_air$sub2, alv_air$scgb_cell_expr)

# Label cells - "airway"
epithelial$airway_firstpass <- ""
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 0))] <- "Ciliated"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 1))] <- "KRT5-low Basal"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 2))] <- "SCGB+ Secretory"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 3))] <- "MUC5B+ Secretory"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = sub.cluster == "4_0"))] <- "Basal"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = sub.cluster == "4_1"))] <- "Proliferating Epithelial"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 5))] <- "Immune"
epithelial$airway_firstpass[epithelial$RNA_snn_res.0.2 == 5] <- "Immune"
epithelial$airway_firstpass[epithelial$RNA_snn_res.0.2 %in% c(0, 2, 4)] <- "Alveolar Epithelium"
airway$CT_firstpass <- epithelial$airway_firstpass
DimPlot(airway, group.by = "CT_firstpass", label = TRUE, raster = FALSE)



#### AIRWAY ----
airway <- subset(epithelial, subset = RNA_snn_res.0.2 %in% c("1", "3"))
DefaultAssay(airway) <- "RNA"
airway <- NormalizeData(airway)
airway <- ScaleData(airway, features = rownames(airway))
airway <- RunPCA(airway, features = rownames(airway))

# Get PCs and make elbow plot
npcs <- min(get_pcs(airway, var = 0.01))
ElbowPlot(airway, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 12 or 23

airway <- FindNeighbors(airway, dims = 1:npcs)
airway <- RunUMAP(airway, dims = 1:npcs)
airway <- FindClusters(airway, resolution = c(0.2))
DimPlot(airway)

airway <- SetIdent(airway, value = "RNA_snn_res.0.2")
air_markers <- FindAllMarkers(airway, only.pos = TRUE)
air_markers <- air_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Look at markers
DotPlot(airway, features = c("EPCAM", "PTPRC", "FOXJ1", "KRT5", "KRT17", "TP63", 
                             "SCGB3A2", "SCGB1A1", "MMP7", "MUC5B", "CPA3"))
FeaturePlot(airway, features = c("FOXJ1", "KRT5", "KRT17", "TP63"), pt.size = 0.0001)
FeaturePlot(airway, features = c("SCGB3A2", "SCGB1A1", "MMP7", "SFTPC"), pt.size = 0.0001)
FeaturePlot(airway, features = c("MUC5B", "MUC5AC", "CPA3"), pt.size = 0.0001)
FeaturePlot(airway, features = c("FOXJ1", "KRT5", "KRT17", "TP63"), pt.size = 0.0001, 
            split.by = "RNA_snn_res.0.2")

hist(subset(airway, subset = RNA_snn_res.0.2 == "1")@assays$RNA@counts["KRT5", ], breaks = 100)
hist(subset(airway, subset = RNA_snn_res.0.2 == "4")@assays$RNA@counts["KRT5", ], breaks = 100)
hist(subset(airway, subset = RNA_snn_res.0.2 == "1")@assays$RNA@data["KRT5", ], breaks = 100)
hist(subset(airway, subset = RNA_snn_res.0.2 == "4")@assays$RNA@data["KRT5", ], breaks = 100)

# Can't subcluster "mast" population
airway <- FindSubCluster(airway, cluster = "5", graph.name = "RNA_snn", resolution = 0.1)
DimPlot(airway, group.by = "sub.cluster", label = TRUE)
DotPlot(airway, features = c("EPCAM", "PTPRC", "FOXJ1", "KRT5", "KRT17", "TP63", 
                             "SCGB3A2", "SCGB1A1", "MMP7", "MUC5B", "CPA3"),
        group.by = "sub.cluster")
airway <- SetIdent(airway, value = "sub.cluster")
air_markers2 <- FindMarkers(airway, ident.1 = "5_0", ident.2 = "5_1")
air_markers2 <- air_markers2 %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2") %>% rownames_to_column()

FeaturePlot(airway, features = c("MKI67"), pt.size = 0.0001)

# Subcluster 1
airway <- FindSubCluster(airway, cluster = "1", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub1")
DimPlot(airway, group.by = "sub1", label = TRUE)
sub1 <- subset(airway, subset = RNA_snn_res.0.2 == "1")
DimPlot(sub1, reduction = "sp", group.by = "sub1")
DotPlot(airway, features = c(epithelial_features, mesenchymal_features), group.by = "sub1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(airway, features = c("EPCAM", "PTPRC", "FOXJ1", "KRT5", "KRT17", "TP63", 
                             "SCGB3A2", "SCGB1A1", "MMP7", "MUC5B", "CPA3"),
        group.by = "sub.cluster")
airway <- SetIdent(airway, value = "sub.cluster")
air_markers2 <- FindMarkers(airway, ident.1 = "5_0", ident.2 = "5_1")
air_markers2 <- air_markers2 %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2") %>% rownames_to_column()

FeaturePlot(airway, features = c("MKI67"), pt.size = 0.0001)

# Subcluster proliferation marker
airway <- SetIdent(airway, value = "RNA_snn_res.0.2")
airway <- FindSubCluster(airway, cluster = "4", graph.name = "RNA_snn", resolution = 0.1)
DotPlot(airway, features = c("FOXJ1", "TP63", "KRT5", "KRT17", "MKI67"), group.by = "sub.cluster")
FeaturePlot(airway, features = "MKI67", split.by = "sub.cluster", pt.size = 0.0001) # 4_1 is proliferating


DotPlot(airway, features = c("FOXJ1", "TP63", "KRT5", "KRT17", "MKI67", "SOX9", "SOX4"), group.by = "sub.cluster")

# Check difference between basal clusters
airway <- SetIdent(airway, value = "sub.cluster")
air_markers3 <- FindMarkers(airway, ident.1 = "1", ident.2 = "4_0")
air_markers3 <- air_markers3 %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2") %>% rownames_to_column()

hist(subset(airway, subset = sub.cluster == "1")@assays$RNA@counts["KRT5", ], breaks = 100)
hist(subset(airway, subset = sub.cluster == "4_0")@assays$RNA@counts["KRT5", ], breaks = 100)
hist(subset(airway, subset = sub.cluster == "1")@assays$RNA@data["KRT5", ], breaks = 100)
hist(subset(airway, subset = sub.cluster == "4_0")@assays$RNA@data["KRT5", ], breaks = 100)

hist(subset(airway, subset = sub.cluster == "1")@assays$RNA@counts["SOX4", ], breaks = 100)
hist(subset(airway, subset = sub.cluster == "4_0")@assays$RNA@counts["SOX4", ], breaks = 100)
hist(subset(airway, subset = sub.cluster == "1")@assays$RNA@data["SOX4", ], breaks = 100)
hist(subset(airway, subset = sub.cluster == "4_0")@assays$RNA@data["SOX4", ], breaks = 100)

# Look again at markers
DotPlot(airway, features = c("EPCAM", "PTPRC", "FOXJ1", "KRT5", "KRT17", "TP63", 
                             "SCGB3A2", "SCGB1A1", "MMP7", "MUC5B", "CPA3", "MKI67"))

# Label cells - epthelial overall
epithelial$airway_firstpass <- ""
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 0))] <- "Ciliated"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 1))] <- "KRT5-low Basal"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 2))] <- "SCGB+ Secretory"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 3))] <- "MUC5B+ Secretory"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = sub.cluster == "4_0"))] <- "Basal"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = sub.cluster == "4_1"))] <- "Proliferating Epithelial"
epithelial$airway_firstpass[colnames(epithelial) %in% colnames(subset(airway, subset = RNA_snn_res.0.2 == 5))] <- "Immune"
epithelial$airway_firstpass[epithelial$RNA_snn_res.0.2 == 5] <- "Immune"
epithelial$airway_firstpass[epithelial$RNA_snn_res.0.2 %in% c(0, 2, 4)] <- "Alveolar Epithelium"

DimPlot(epithelial, group.by = "CT_firstpass", split.by = "CT_firstpass")
DimPlot(epithelial, group.by = "airway_firstpass", label = TRUE, repel = TRUE, raster = FALSE)

DimPlot(epithelial, group.by = "CT_firstpass", label = TRUE)


#### ALVEOLAR ----
alv <- subset(epithelial, subset = RNA_snn_res.0.2 %in% c(0, 2, 4))
DefaultAssay(alv) <- "RNA"
alv <- NormalizeData(alv)
alv <- ScaleData(alv, features = rownames(alv))
alv <- RunPCA(alv, features = rownames(alv))

# Get PCs and make elbow plot
npcs <- min(get_pcs(alv, var = 0.01))
ElbowPlot(alv, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 12 or 23

alv <- FindNeighbors(alv, dims = 1:npcs)
alv <- RunUMAP(alv, dims = 1:npcs)
alv <- FindClusters(alv, resolution = c(0.1, 0.2))
DimPlot(alv)

alv <- SetIdent(alv, value = "RNA_snn_res.0.1")
alv_markers <- FindAllMarkers(alv, only.pos = TRUE)
alv_markers <- alv_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

DotPlot(alv, features = c("SFTPC", "GKN2", "NAPSA", "PGC", "SFTA2", "SFTPD",
                          "RTKN2", "AGER", "CEACAM6", "KRT8", "SOX9", "SPINK1",
                          "SCGB3A2", "SCGB1A1", "KRT5", "KRT17", "TP63", "CPA3"))

FeaturePlot(alv, features = c("SFTPC", "CPA3"), blend = TRUE)
FeaturePlot(alv, features = c("SFTPC", "RTKN2"), split.by = "RNA_snn_res.0.1")
FeaturePlot(airway, features = c("SFTPC"), split.by = "RNA_snn_res.0.2")

FeaturePlot(alv, features = c("MKI67"), split.by = "RNA_snn_res.0.1")
FeaturePlot(alv, features = c("SFTPC", "AGER", "CEACAM6"), split.by = "RNA_snn_res.0.1")

# 0 = AT2
# 1 = AT1
# 2 = Airway? or transitional? (may need to split)
# 3 = Mast/SFTPC (may need to split)
# 4 = Proliferating!

# Subset airway-like clusters to get finer detail
alv <- FindSubCluster(alv, cluster = "2", graph.name = "RNA_snn", resolution = 0.05)
DimPlot(alv, group.by = "sub.cluster", label = TRUE)
FeaturePlot(alv, features = c("SFTPC", "AGER", "KRT8"), split.by = "sub.cluster")

DotPlot(alv, features = c("SFTPC", "GKN2", "NAPSA", "PGC", "SFTA2", "SFTPD",
                          "RTKN2", "AGER", "CEACAM6", "KRT8", "SOX9", "SPINK1",
                          "SCGB3A2", "SCGB1A1", "KRT5", "KRT17", "TP63", "CPA3", 
                          "MRC1", "MKI67"),
        group.by = "sub.cluster")

DimPlot(alv, split.by = "RNA_snn_res.0.1")
FeaturePlot(alv, features = c("SCGB3A2", "KRT5", "KRT17", "FOXJ1", "TP73", "MMP7"))

FeaturePlot(alv, features = "MKI67")

alv_air <- subset(alv, subset = RNA_snn_res.0.1 == "2")
alv_air <- SetIdent(alv_air, value = "sub.cluster")
alv_air_markers <- FindAllMarkers(alv_air, only.pos = TRUE)
alv_air_markers <- alv_air_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

FeaturePlot(alv_air, features = c("SCGB3A2", "SCGB1A1", "KRT5", "KRT17", "FOXJ1"), split.by = "sub.cluster")
FeaturePlot(alv_air, features = c("MUC5B", "FOXJ1"), split.by = "sub.cluster")

alv <- SetIdent(alv, value = "sub.cluster")
alv <- FindSubCluster(alv, cluster = "2_0", graph.name = "RNA_snn", 
                      resolution = 0.3, subcluster.name = "sub2")
DimPlot(alv, group.by = "sub2", label = TRUE)
alv_air <- subset(alv, subset = RNA_snn_res.0.1 == "2")
FeaturePlot(alv_air, features = c("MUC5B", "FOXJ1", "SCGB3A2", "SCGB1A1"), split.by = "sub2")
FeaturePlot(alv_air, features = c("MUC5B", "FOXJ1", "SCGB3A2", "SCGB1A1"))
VlnPlot(alv_air, features = c("MUC5B", "FOXJ1", "SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA"), group.by = "sub2", pt.size = 0, ncol = 2)
DotPlot(alv_air, features = c("MUC5B", "FOXJ1", "SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA",
                              "KRT5", "KRT17", "SPINK1", "SOX9", "TP63"), group.by = "sub2")
DimPlot(alv_air, group.by = "sub2", label = TRUE)
FeaturePlot(alv_air, features = c("KRT5", "KRT17", "SPINK1", "SOX9", "KRT8")) # Basal/KRT5-/KRT17+
FeaturePlot(alv_air, features = c("KRT17", "CCL18"), split.by = "sub2", blend = TRUE)

# Lots of JCHAIN in the ciliated cluster - can't separate
DotPlot(alv, features = c("JCHAIN", "FOXJ1"), group.by = "sub2")
FeaturePlot(alv_air, features = c("JCHAIN", "FOXJ1"), blend = TRUE)
alv <- SetIdent(alv, value = "sub2")
alv <- FindSubCluster(alv, cluster = "2_0_3", graph.name = "RNA_snn", 
                      resolution = 0.3, subcluster.name = "sub_ciliated")
DimPlot(alv, group.by = "sub_ciliated", label = TRUE)
DotPlot(alv, features = c("JCHAIN", "FOXJ1"), group.by = "sub_ciliated")

alv <- SetIdent(alv, value = "sub2")
alv_air <- SetIdent(alv_air, value = "sub2")
alv_air_markers2 <- FindAllMarkers(alv_air, only.pos = TRUE)
alv_air_markers2 <- alv_air_markers2 %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

alv_air2 <- SetIdent(alv_air, value = "sub2")
alv_air_markers2 <- FindAllMarkers(alv_air2, only.pos = TRUE)
alv_air_markers2 <- alv_air_markers2 %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(alv, features = c("PECAM1", "EPCAM", "SPARCL1"), group.by = "sub2")
DotPlot(alv, features = c("SCGB3A2", "SFTPD", "SCGB1A1"), group.by = "sub2")

# Split mast-like cluster
alv <- FindSubCluster(alv, cluster = "3", graph.name = "RNA_snn", 
                      resolution = 0.05, subcluster.name = "mast_sub")
DimPlot(alv, group.by = "mast_sub", label = TRUE)
FeaturePlot(alv, features = "CPA3", split.by = "mast_sub")
DotPlot(alv, features = c("CPA3", "SFTPC", "SCGB3A2", "KRT8", "FOXJ1"), group.by = "mast_sub")

alv <- SetIdent(alv, value = "mast_sub")
alv_mast_markers <- FindMarkers(alv, ident.1 = "3_0", ident.2 = "3_1")
alv_mast_markers <- alv_mast_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

DimPlot(alv, group.by = "CT_firstpass", label = TRUE)
DimPlot(alv, group.by = "sub2", label = TRUE)

DotPlot(alv, features = c("SFTPC", "GKN2", "NAPSA", "PGC", "SFTA2", "SFTPD",
                          "RTKN2", "AGER", "CEACAM6", "KRT8", "SOX9", "SPINK1",
                          "SCGB3A2", "SCGB1A1", "KRT5", "KRT17", "TP63", "CPA3", 
                          "MRC1", "MKI67"), group.by = "sub2")

FeaturePlot(alv, features = c("CPA3", "KIT", "CD34", "FCN3", "MKI67"))
VlnPlot(alv, features = c("CPA3", "CD34", "FCN3", "MKI67"), pt.size = 0, group.by = "sub2", ncol = 2)

FeaturePlot(alv, features = c("SFTPC", "NAPSA", "AGER", "RTKN2", "SCGB3A2", "SCGB3A2"))
VlnPlot(alv, features = c("SFTPC", "NAPSA", "AGER", "RTKN2", "SCGB3A2", "SCGB3A2"), group.by = "sub2", pt.size = 0, ncol = 2)
FeaturePlot(alv, features = c("SCGB3A2", "SCGB3A2", "TP63", "TP73", "FOXJ1"))
FeaturePlot(alv, features = c("TP63", "TP73", "FOXJ1", "MUC5B", "MUC5AC"))

# Keep how it is without mast split
# Label cell types
epithelial$CT_firstpass <- epithelial$airway_firstpass
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 == 0))] <- "AT2"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 == 1))] <- "AT1"

epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 %in% c("2_0_0", "2_0_2")))] <- "SCGB+ Secretory"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 == "2_0_1"))] <- "AT2"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 == "2_0_3"))] <- "Differentiating Ciliated"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 == "2_0_4"))] <- "MUC5B+ Secretory"

epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 == "2_1"))] <- "KRT5-/KRT17+"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 == "2_2"))] <- "Endothelial"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 %in% c("3_0", "3_1")))] <- "Immune"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = sub2 == 4))] <- "Proliferating Epithelial"
DimPlot(epithelial, group.by = "CT_firstpass", label = TRUE, raster = FALSE)
DimPlot(epithelial, group.by = "RNA_snn_res.0.2", label = TRUE)

# Look at lineage markers
VlnPlot(epithelial, features = c("PECAM1", "EPCAM", "PTPRC", "LUM"),
        pt.size = 0, group.by = "CT_firstpass")

# Look at typical epithelial markers and other top markers
VlnPlot(epithelial, features = c("SFTPC", "NAPSA", "SFTPD", "MKI67", "RTKN2", 
                                 "AGER", "KRT8", "CPA3"), # MRC1 has strange expression
        pt.size = 0, group.by = "CT_firstpass")
VlnPlot(epithelial, features = c("SOX4", "SOX9", "SCGB3A2", "MUC5B",
                                 "SCGB1A1", "KRT5", "KRT17", "TP63"),
        pt.size = 0, group.by = "CT_firstpass")
VlnPlot(epithelial, features = c("MRC1", "KRT5", "KRT17", "TP63", "EPCAM", "PTPRC", "LYZ"),
        pt.size = 0, group.by = "CT_firstpass")
VlnPlot(epithelial, features = c("TP73", "FOXJ1", "MUC5B", "SCGB1A1", "SCGB3A2"),
        pt.size = 0, group.by = "CT_firstpass")

# Look at SCGB markers - nuclear and cellular
scgb1a1_nuc_expr <- names(which((alv_air@assays$RNA@counts["SCGB1A1", ] > 0), arr.ind = TRUE))
scgb3a2_nuc_expr <- names(which((alv_air@assays$RNA@counts["SCGB3A2", ] > 0), arr.ind = TRUE))

alv_air$scgb_nuc_expr <- "SCGB3A2-/SCGB1A1-"
alv_air$scgb_nuc_expr[colnames(alv_air) %in% scgb1a1_nuc_expr] <- "SCGB3A2-/SCGB1A1+"
alv_air$scgb_nuc_expr[colnames(alv_air) %in% scgb3a2_nuc_expr] <- "SCGB3A2+/SCGB1A1-"
alv_air$scgb_nuc_expr[colnames(alv_air) %in% scgb1a1_nuc_expr & colnames(alv_air) %in% scgb3a2_nuc_expr] <- "SCGB3A2+/SCGB1A1+"
DimPlot(alv_air, group.by = "scgb_nuc_expr")

scgb1a1_cell_expr <- names(which((alv_air@assays$cell_RNA@counts["SCGB1A1", ] > 0), arr.ind = TRUE))
scgb3a2_cell_expr <- names(which((alv_air@assays$cell_RNA@counts["SCGB3A2", ] > 0), arr.ind = TRUE))
alv_air$scgb_cell_expr <- "SCGB3A2-/SCGB1A1-"
alv_air$scgb_cell_expr[colnames(alv_air) %in% scgb1a1_cell_expr] <- "SCGB3A2-/SCGB1A1+"
alv_air$scgb_cell_expr[colnames(alv_air) %in% scgb3a2_cell_expr] <- "SCGB3A2+/SCGB1A1-"
alv_air$scgb_cell_expr[colnames(alv_air) %in% scgb1a1_cell_expr & colnames(alv_air) %in% scgb3a2_cell_expr] <- "SCGB3A2+/SCGB1A1+"
DimPlot(alv_air, group.by = "scgb_cell_expr", cols = c("grey", "red", "blue", "green"), split.by = "sub2")
DimPlot(alv_air, group.by = "scgb_nuc_expr", cols = c("grey", "red", "blue", "green"), split.by = "sub2")

# Look at KRT5/KRT17 markers - nuclear and cellular
krt5_nuc_expr <- names(which((alv_air@assays$RNA@counts["KRT5", ] > 0), arr.ind = TRUE))
krt17_nuc_expr <- names(which((alv_air@assays$RNA@counts["KRT17", ] > 0), arr.ind = TRUE))

alv_air$krt_nuc_expr <- "KRT5-/KRT17-"
alv_air$krt_nuc_expr[colnames(alv_air) %in% krt5_nuc_expr] <- "KRT5-/KRT17+"
alv_air$krt_nuc_expr[colnames(alv_air) %in% krt17_nuc_expr] <- "KRT5-/KRT17+"
alv_air$krt_nuc_expr[colnames(alv_air) %in% krt5_nuc_expr & colnames(alv_air) %in% krt17_nuc_expr] <-"KRT5+/KRT17+"
DimPlot(alv_air, group.by = "krt_nuc_expr")

krt5_cell_expr <- names(which((alv_air@assays$cell_RNA@counts["KRT5", ] > 0), arr.ind = TRUE))
krt17_cell_expr <- names(which((alv_air@assays$cell_RNA@counts["KRT17", ] > 0), arr.ind = TRUE))
alv_air$krt_cell_expr <- "KRT5-/KRT17-"
alv_air$krt_cell_expr[colnames(alv_air) %in% krt5_cell_expr] <- "KRT5-/KRT17+"
alv_air$krt_cell_expr[colnames(alv_air) %in% krt17_cell_expr] <- "KRT5-/KRT17+"
alv_air$krt_cell_expr[colnames(alv_air) %in% krt5_cell_expr & colnames(alv_air) %in% krt17_cell_expr] <-"KRT5+/KRT17+"
DimPlot(alv_air, group.by = "krt_cell_expr")

table(alv_air$sub2, alv_air$krt_nuc_expr)
table(alv_air$sub2, alv_air$krt_cell_expr)

# Label cells - "alveolar"
alv$CT_firstpass <- ""
alv$CT_firstpass[alv$sub2 %in% c("0", "2_0_1")] <- "AT2"
alv$CT_firstpass[alv$sub2 == "1"] <- "AT1"
alv$CT_firstpass[alv$sub2 %in% c("2_0_0", "2_0_2")] <- "SCGB+ Secretory"
alv$CT_firstpass[alv$sub2 == "2_0_3"] <- "Differentiating Ciliated"
alv$CT_firstpass[alv$sub2 == "2_0_4"] <- "MUC5B+ Secretory"
alv$CT_firstpass[alv$sub2 == "2_1"] <- "KRT5-/KRT17+"
alv$CT_firstpass[alv$sub2 == "2_2"] <- "Endothelial"
alv$CT_firstpass[alv$sub2 %in% c("3_0", "3_1")] <- "Immune"
alv$CT_firstpass[alv$sub2 == "4"] <- "Proliferating Epithelial"
DimPlot(alv, group.by = "CT_firstpass", label = TRUE)

# Label cells - epithelial
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "AT2"))] <- "AT2"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "AT1"))] <- "AT1"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "SCGB+ Secretory"))] <- "SCGB+ Secretory"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "Differentiating Ciliated"))] <- "Differentiating Ciliated"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "MUC5B+ Secretory"))] <- "MUC5B+ Secretory"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "KRT5-/KRT17+"))] <- "KRT5-/KRT17+"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "Endothelial"))] <- "Endothelial"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "Immune"))] <- "Immune"
epithelial$CT_firstpass[colnames(epithelial) %in% colnames(subset(alv, subset = CT_firstpass == "Proliferating Epithelial"))] <- "Proliferating Epithelial"
DimPlot(epithelial, group.by = "CT_firstpass", label = TRUE)

#saveRDS(airway, "/scratch/avannan/xenium_files/airway_firstpass.rds")
airway <- readRDS("/scratch/avannan/xenium_files/airway_firstpass.rds")
#saveRDS(alv, "/scratch/avannan/xenium_files/alveolar_firstpass.rds")
#saveRDS(epithelial, "/scratch/avannan/xenium_files/epithelial_firstpass.rds")
test <- readRDS("/scratch/avannan/xenium_spatial.rds")
DimPlot(xenium, group.by = "CT_firstpass")
DimPlot(epithelial_test, group.by = "sub.cluster", label = TRUE)


DimPlot(epithelial, group.by = "CT_firstpass", label = TRUE, repel = TRUE, raster = FALSE)
DimPlot(epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, repel = TRUE, raster = FALSE)


## ENDOTHELIAL ----
# Subset endothelial clusters from starting annotation and from epithelial clusters
endothelial <- subset(merged_spatial, subset = RNA_snn_res.0.5 %in% c(4, 7, 16))

DefaultAssay(endothelial) <- "RNA"
endothelial <- NormalizeData(endothelial)
endothelial <- ScaleData(endothelial, features = rownames(endothelial))
endothelial <- RunPCA(endothelial, features = rownames(endothelial))

# Get PCs and make elbow plot
npcs <- min(get_pcs(endothelial))
ElbowPlot(endothelial, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 11

# Find neighbors, cluster and UMAP
endothelial <- FindNeighbors(endothelial, dims = 1:npcs)
endothelial <- RunUMAP(endothelial, dims = 1:npcs)
endothelial <- FindClusters(endothelial, resolution = c(0.1, 0.2))
DimPlot(endothelial, group.by = "RNA_snn_res.0.2", label = TRUE)

endothelial <- SetIdent(endothelial, value = "RNA_snn_res.0.2")
endo_markers <- FindAllMarkers(endothelial, only.pos = TRUE)
endo_markers <- endo_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Look at lineage markers
DotPlot(endothelial, features = c("PECAM1", "EPCAM", "PTPRC", "LUM", "ACTA2", "DCN"))
FeaturePlot(endothelial, features = c("PECAM1", "EPCAM", "PTPRC", "LUM", "ACTA2", "DCN"))

# Look at some of the top markers and lineage markers
FeaturePlot(endothelial, features = c("POSTN", "ACKR1", "PLVAP", "SPARCL1", "COL1A1", "MYC", "LEF1"))
DotPlot(endothelial, features = c("POSTN", "ACKR1", "PLVAP", "SPARCL1", "COL1A1", 
                                  "PDGFRA", "S100A9", "S100A8", "FCN1", "FCER1G", "MS4A7"))
VlnPlot(endothelial, features = c("PECAM1", "EPCAM", "PTPRC", "LUM", "ACTA2", "DCN",
                                  "GZMB", "HLA-DRA", "XBP1", "PPARG", "CD79A", "S100A9"),
        pt.size = 0, group.by = "RNA_snn_res.0.2")

# Look at typical endothelial markers
VlnPlot(endothelial, features = c("CA4", "APLN", "APLNR", "CCL21", "ACKR1", 
                                  "HEY1", "PLVAP", "COL15A1", "CCL21"),
        pt.size = 0, group.by = "RNA_snn_res.0.2")
FeaturePlot(endothelial, features = c("CA4", "APLN", "APLNR", "CCL21", "ACKR1", 
                                      "HEY1", "PLVAP", "COL15A1", "CCL21"))
DotPlot(endothelial, features = c("CA4", "APLN", "APLNR", "CCL21", "ACKR1", 
                                  "HEY1", "PLVAP", "COL15A1"))

# Immune top markers
VlnPlot(endothelial, features = c("PDGFRB", "CD14", "FCER1G", "PPARG",
                                  "S100A8", "S100A9"), pt.size = 0)

# HLA-DRA is high - check if other DC markers are high. They aren't
VlnPlot(endothelial, features = c("CD1A", "CD1C", "CXCL9", "GPR183", "GZMB", "HLA-DQA1", 
                                  "HLA-DQB1", "HLA-DRA", "MMP12"), pt.size = 0)

# Check markers of cluster 2 and other mesenchymal genes
VlnPlot(endothelial, features = c("CCL21", "FABP4", "ACTA2"), pt.size = 0)
FeaturePlot(endothelial, features = c("ACTA2", "LUM", "DCN", "COL3A1"))

# Split 2 further to see if there are endothelial cells
endothelial <- FindSubCluster(endothelial, cluster = "2", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub2")
DimPlot(endothelial, group.by = "sub.cluster", label = TRUE)
endothelial <- SetIdent(endothelial, value = "sub.cluster")
endo2_markers <- FindMarkers(endothelial, ident.1 = "2_0", ident.2 = "2_1") %>%
  mutate(pct.diff = abs(pct.1 - pct.2))
endothelial <- SetIdent(endothelial, value = "RNA_snn_res.0.2")
og_endothelial <- endothelial
DimPlot(og_endothelial, group.by = "CT_firstpass", label = TRUE)


#### TRUE ENDOTHELIAL ----
# Recluster true endothelial cells
# Re-subset endothelial clusters from starting annotation and from epithelial clusters
endothelial <- subset(endothelial, subset = sub2 %in% c("0", "2_0", "1", "3"))

DefaultAssay(endothelial) <- "RNA"
endothelial <- NormalizeData(endothelial)
endothelial <- ScaleData(endothelial, features = rownames(endothelial))
endothelial <- RunPCA(endothelial, features = rownames(endothelial))

# Get PCs and make elbow plot
npcs <- min(get_pcs(endothelial))
ElbowPlot(endothelial, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 12

# Find neighbors, cluster and UMAP
endothelial <- FindNeighbors(endothelial, dims = 1:npcs)
endothelial <- RunUMAP(endothelial, dims = 1:npcs)
endothelial <- FindClusters(endothelial, resolution = c(0.2, 0.3))
DimPlot(endothelial, group.by = "RNA_snn_res.0.2", label = TRUE)

endothelial <- SetIdent(endothelial, value = "RNA_snn_res.0.2")
endo_markers <- FindAllMarkers(endothelial, only.pos = TRUE)
endo_markers <- endo_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Look at lineage markers
VlnPlot(endothelial, features = c("PECAM1", "LUM", "ACTA2", "FABP4", "KDR", "GNG11"), pt.size = 0)
FeaturePlot(endothelial, features = c("PECAM1", "EPCAM", "PTPRC", "LUM", "ACTA2", "DCN", "FABP4", "KDR", "GNG11"))

# Look at general endothelial markers
FeaturePlot(endothelial, features = c("CA4", "APLN", "APLNR", "CCL21", "ACKR1", "HEY1", "PLVAP", "COL15A1", "RSPO3", "FCN3"))
FeaturePlot(endothelial, features = c("CA4", "APLN", "APLNR", "FCN3"), ncol = 4)
FeaturePlot(endothelial, features = c("APLN", "FCN3"), blend = TRUE)
DotPlot(endothelial, features = c("CA4", "APLN", "APLNR", "CCL21", "ACKR1", "HEY1", "PLVAP", "COL15A1", "RSPO3", "FCN3",
                                  "SFTPC"))
DimPlot(endothelial, label = TRUE)

# Subcluster 1
endothelial <- FindSubCluster(endothelial, cluster = "1", graph.name = "RNA_snn", resolution = 0.3)
DimPlot(endothelial, group.by = "sub.cluster", label = TRUE)
VlnPlot(endothelial, features = c("ACKR1", "PLVAP", "COL15A1"), pt.size = 0, group.by = "sub.cluster")

# Subcluster 0 (capillary)
endothelial <- FindSubCluster(endothelial, cluster = "0", graph.name = "RNA_snn", resolution = 0.3)
DimPlot(endothelial, group.by = "sub.cluster", label = TRUE)
VlnPlot(endothelial, features = c("CA4", "APLN", "APLNR", "TBXA2R", "VEGFA", "RSPO3"), pt.size = 0, group.by = "sub.cluster")
DotPlot(endothelial, features = c("CA4", "APLN", "APLNR", "TBXA2R", "VEGFA", "RSPO3"), group.by = "sub.cluster")
FeaturePlot(endothelial, features = c("CA4", "APLN", "APLNR", "TBXA2R", "VEGFA", "RSPO3", "PLVAP"))

# Label cell types
DimPlot(endothelial)
endothelial$CT_firstpass <- ""
endothelial$CT_firstpass[endothelial$RNA_snn_res.0.2 == 0] <- "Capillary"
endothelial$CT_firstpass[endothelial$RNA_snn_res.0.2 == 1] <- "Venous"
endothelial$CT_firstpass[endothelial$RNA_snn_res.0.2 == 2] <- "Arteriole"
DimPlot(endothelial, group.by = "CT_firstpass", label = TRUE, repel = TRUE)

#saveRDS(endothelial, "/scratch/avannan/xenium_files/endothelial_firstpass_04072023.rds")
endothelial <- readRDS("/scratch/avannan/xenium_files/endothelial_firstpass_04072023.rds")
DimPlot(endothelial, group.by = "CT_firstpass", label = TRUE)

# Lable endothelial as a whole
og_endothelial$CT_firstpass <- endothelial$CT_firstpass
og_endothelial$CT_firstpass[og_endothelial$RNA_snn_res.0.2 == "2"] <- "Mesenchymal"
og_endothelial$CT_firstpass[og_endothelial$RNA_snn_res.0.2 == "4"] <- "Immune"
DimPlot(og_endothelial, group.by = "CT_firstpass", label = TRUE, raster = FALSE)
#saveRDS(og_endothelial, "/scratch/avannan/xenium_files/og_endothelial_firstpass_04142023.rds")
og_endothelial <- readRDS("/scratch/avannan/xenium_files/og_endothelial_firstpass_04142023.rds")


## IMMUNE ----
immune <- subset(merged_spatial, subset = RNA_snn_res.0.5 %in% c(1, 6, 11, 13, 14, 18))

DefaultAssay(immune) <- "RNA"
immune <- NormalizeData(immune)
immune <- ScaleData(immune, features = rownames(immune))
immune <- RunPCA(immune, features = rownames(immune))

# Get PCs and make elbow plot
npcs <- min(get_pcs(immune))
ElbowPlot(immune, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 16

# Find neighbors, cluster and UMAP
immune <- FindNeighbors(immune, dims = 1:npcs)
immune <- RunUMAP(immune, dims = 1:npcs)
immune <- FindClusters(immune, resolution = c(0.1, 0.2, 0.3))
DimPlot(immune, group.by = "RNA_snn_res.0.2", label = TRUE)

# General lineage markers
VlnPlot(immune, features = c("PECAM1", "EPCAM", "PTPRC", "DCN"), pt.size = 0, ncol = 2)
FeaturePlot(immune, features = c("PECAM1", "EPCAM", "PTPRC", "DCN"))
DotPlot(immune, features = c("PECAM1", "EPCAM", "PTPRC", "DCN"))

immune <- SetIdent(immune, value = "RNA_snn_res.0.2")
imm_markers <- FindAllMarkers(immune, only.pos = TRUE)
imm_markers <- imm_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Look at some of the top markers
VlnPlot(immune, features = c("CPA3", "KIT", "MS4A1", "BANK1", "CD3E", "CD4", 
                             "CD8A", "GNLY", "NKG7", "MKI67", "S100A8", "CD14",
                             "FCGR3A", "LILRA4", "CCL22", "HLA-DQA1"), pt.size = 0, ncol = 4)

FeaturePlot(immune, features = c("CPA3", "KIT", "MS4A1", "BANK1")) # Mast and B
VlnPlot(immune, features = c("CPA3", "KIT", "MS4A1", "BANK1"), pt.size = 0, ncol = 4) # Mast and B
FeaturePlot(immune, features = c("CCL22", "HLA-DQA1", "COL4A3")) # Dendritic cells - CCL22+, with high COL4A3
FeaturePlot(immune, features = c("MKI67")) # Proliferating - throughout
VlnPlot(immune, features = c("MKI67"), pt.size = 0) # Proliferating - throughout, but mostly in 6
FeaturePlot(immune, features = c("CD3E", "CD4", "CD8A", "FOXP3")) # T-cells
DotPlot(immune, features = c("CD3E", "CD4", "CD8A", "FOXP3")) # T-cells
VlnPlot(immune, features = c("CD3E", "CD4", "CD8A", "FOXP3"), pt.size = 0, ncol = 4) # T-cells
FeaturePlot(immune, features = c("GNLY", "NKG7")) # NK cells
DotPlot(immune, features = c("GNLY", "NKG7")) # NK cells
VlnPlot(immune, features = c("GNLY", "NKG7"), pt.size = 0) # NK cells
FeaturePlot(immune, features = c("S100A8", "S100A9", "CD14", "FCGR3A")) # Myeloid
DotPlot(immune, features = c("S100A8", "S100A9", "CD14", "FCGR3A")) # Myeloid
VlnPlot(immune, features = c("S100A8", "S100A9", "CD14", "FCGR3A"), pt.size = 0, ncol = 2) # Myeloid

VlnPlot(immune, features = c("IRF7", "LILRA4"), pt.size = 0) # pDCs
FeaturePlot(immune, features = c("IRF7", "LILRA4")) # pDCs

# Plasma markers and pDCs
FeaturePlot(immune, features = c("EHMT1", "FKBP11", "LMAN1", "PDIA4", "PIM2",
                                 "SSR3", "UBE2J1", "LILRA4", "IRF7")) # Last 2 are pDC markers

# Plasma/pDC/DC
DotPlot(immune, features = c("CCL22", "HLA-DQA1", "COL4A3", # Collagen + DC
                             "DERL3", "LMAN1", "PDIA4", "PIM2",
                             "SPCS3", "SSR3", "UBE2J1",
                             "LILRA4", "IRF7")) # Last 2 are pDC markers
# Cluster 8 is plasma/pDC

# There are ACTA2+ cells that will need to be separated from mast and myeloid
FeaturePlot(immune, features = c("S100A8", "ACTA2"), blend = TRUE)
FeaturePlot(immune, features = c("CPA3", "ACTA2"), blend = TRUE)


###### LYMPHOID ----
# Re-UMAP lymphoid cells
lymphoid <- subset(immune, subset = RNA_snn_res.0.2 %in% c(0, 2, 5, 6, 8))
DefaultAssay(lymphoid) <- "RNA"
lymphoid <- NormalizeData(lymphoid)
lymphoid <- ScaleData(lymphoid, features = rownames(lymphoid))
lymphoid <- RunPCA(lymphoid, features = rownames(lymphoid))

# Get PCs and make elbow plot
npcs <- min(get_pcs(lymphoid))
ElbowPlot(lymphoid, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs

# Find neighbors, cluster and UMAP
lymphoid <- FindNeighbors(lymphoid, dims = 1:npcs)
lymphoid <- RunUMAP(lymphoid, dims = 1:npcs)
lymphoid <- FindClusters(lymphoid, resolution = c(0.2, 0.25, 0.3))
DimPlot(lymphoid, group.by = "RNA_snn_res.0.25", label = TRUE)
DimPlot(lymphoid, group.by = "tma")
DimPlot(lymphoid, group.by = "sample")

lymphoid <- SetIdent(lymphoid, value = "RNA_snn_res.0.25")
FeaturePlot(lymphoid, features = c("CD3E", "CD8A", "CD4", "FOXP3", "IL2RA", "LEF1", "NKG7", "GNLY",
                                   "MS4A1", "LILRA4", "XBP1", "MKI67", "CD247"))
DotPlot(lymphoid, features = c("CD3E", "CD8A", "CD4", "FOXP3", "IL2RA", "LEF1", "NKG7", "GNLY",
                               "MS4A1", "LILRA4", "XBP1", "MKI67", "SFTPC", "PDGFRB", "MEG3",
                               "FCER1G", "CCL21", "CD68", "COL1A1", "COL3A1", "FN1"), group.by = "RNA_snn_res.0.3")
VlnPlot(lymphoid, features = c("CD3E", "CD8A", "CD4", "NKG7", "GNLY", "LEF1", "KLRB1",
                               "MS4A1", "LILRA4"), pt.size = 0)

FeaturePlot(lymphoid, features = c("CD4", "CD8A"), blend = TRUE, split.by = "RNA_snn_res.0.25")
FeaturePlot(lymphoid, features = c("CD4", "CD8A", "CD3E", "FOXP3", "NKG7"), split.by = "RNA_snn_res.0.25")
FeaturePlot(lymphoid, features = c("FN1", "MEG3", "COL1A1"), split.by = "RNA_snn_res.0.25")
FeaturePlot(lymphoid, features = c("SFTPC", "NAPSA", "PGC"), ncol = 3)

# Look at top markers
lymphoid <- SetIdent(lymphoid, value = "RNA_snn_res.0.25")
lym_markers <- FindAllMarkers(lymphoid, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
VlnPlot(lymphoid, features = c("CD3E", "GNLY", "COL3A1", "COL1A1", "FN1", 
                               "S100A8", "S100A9", "LEF1", "FOXP3", "PDGFRB", 
                               "ACTA2", "MS4A7", "PLIN2", "CD14", "MKI67", 
                               "TOP2A"), pt.size = 0)
DotPlot(lymphoid, features = c("CD3E", "GNLY", "COL3A1", "FN1", "S100A8", 
                               "S100A9", "PDGFRB", "MS4A7", "PLIN2", "CD14",
                               "MKI67", "TOP2A"))

# Check cluster 6 for overlap of CD3E and SFTPC
FeaturePlot(subset(lymphoid, subset = RNA_snn_res.0.25 == "6"), features = c("SFTPC", "CD3E"), blend = TRUE, pt.size = 0.0001)

# Subcluster 6
lymphoid <- FindSubCluster(lymphoid, cluster = "6", subcluster.name = "sub6",
                           graph.name = "RNA_snn", resolution = 0.15)
sub6 <- subset(lymphoid, subset = RNA_snn_res.0.25 == "6")
DimPlot(lymphoid, group.by = "sub6", label = TRUE)
VlnPlot(lymphoid, features = c("SFTPC", "NAPSA", "CD3E"), group.by = "sub6", pt.size = 0)

# Subcluster 4 (NK + something else)
lymphoid <- FindSubCluster(lymphoid, cluster = "4", subcluster.name = "sub4",
                           graph.name = "RNA_snn", resolution = 0.1)
DimPlot(lymphoid, group.by = "sub4", label = TRUE)
lymphoid <- SetIdent(lymphoid, value = "sub4")
sub4_markers <- FindMarkers(lymphoid, ident.1 = "4_0", ident.2 = "4_1")
VlnPlot(lymphoid, features = c("CD3E", "GNLY", "NKG7", "COL3A1", "COL1A1", "FN1",
                               "PECAM1", "CD34", "FCN3", "KDR"), pt.size = 0)
DotPlot(lymphoid, features = c("CD3E", "CD8A", "CD4", "FOXP3", "IL2RA", "LEF1", "NKG7", "GNLY",
                               "MS4A1", "LILRA4", "XBP1", "MKI67", "SFTPC", "PDGFRB", "MEG3",
                               "FCER1G", "CCL21", "CD68", "COL1A1", "COL3A1", "FN1",
                               "PECAM1", "CD34", "FCN3", "KDR", "CA4", "HEY1", "ACKR1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(lymphoid, features = c("PDGFRB", "CCL21", "COL1A1", "COL3A1", "FN1", "DCN", "LUM", "IL7R", "ITGB6"), pt.size = 0)
FeaturePlot(lymphoid, features = c("FCN3", "PECAM1", "CD34", "CA4"))

# Subcluster 0
lymphoid <- FindSubCluster(lymphoid, cluster = "0", subcluster.name = "sub0",
                           graph.name = "RNA_snn", resolution = 0.1)
DimPlot(lymphoid, group.by = "sub0", label = TRUE)
DotPlot(lymphoid, features = c("CD3E", "CD3D", "CD4", "CD8A", "CD8B", "FOXP3"), group.by = "sub0")
lymphoid <- SetIdent(lymphoid, value = "sub1")
sub1 <- subset(lymphoid, subset = sub1 %in% c("1_0", "1_1", "1_2"))
sub1_markers <- FindAllMarkers(sub1, only.pos = TRUE) %>% mutate(pct.diff = abs(pct.1 - pct.2))

# Subcluster 1
lymphoid <- FindSubCluster(lymphoid, cluster = "1", subcluster.name = "sub1",
                           graph.name = "RNA_snn", resolution = 0.2)
DimPlot(lymphoid, group.by = "sub1", label = TRUE)
lymphoid <- SetIdent(lymphoid, value = "sub1")
sub1 <- subset(lymphoid, subset = sub1 %in% c("1_0", "1_1", "1_2"))
sub1_markers <- FindAllMarkers(sub1, only.pos = TRUE) %>% mutate(pct.diff = abs(pct.1 - pct.2))


# Look at fibroblast markers
FeaturePlot(lymphoid, features = c("COL1A1", "COL3A1", "FN1"))
DotPlot(lymphoid, features = c("AXL", "LGR6", "SPARCL1", "CCL2", "COL1A1", "COL1A2", 
                               "COL3A1", "CXCL14", "DCN", "ELN", "FAP", "FGF10", 
                               "FGF2", "FGF7", "FN1", "HAS2", "IL11", "ITGAV", "LUM", 
                               "MAL", "MEG3", "PDGFRA", "RSPO3", "SNAI2", "TGFB3", 
                               "UGDH", "WNT2", "WNT5A")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(lymphoid, features = c("CD247", "FGFBP2", "GNLY", "GZMA", "KLRB1", 
                               "KLRC1", "NKG7", "CD3E", "CD3D", "CD4", "CD8A", "FOXP3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Label cells
lymphoid$CT_firstpass <- ""
lymphoid$CT_firstpass[lymphoid$sub0 == "0_1"] <- "CD4+ Treg"
lymphoid$CT_firstpass[lymphoid$sub0 == "0_0"] <- "CD4+ T-cell"
lymphoid$CT_firstpass[lymphoid$sub1 == "1_0"] <- "T-cell"
lymphoid$CT_firstpass[lymphoid$sub1 == "1_1"] <- "Mesenchymal"
lymphoid$CT_firstpass[lymphoid$sub1 == "1_2"] <- "Mesenchymal"
lymphoid$CT_firstpass[lymphoid$RNA_snn_res.0.25 == "2"] <- "CD8+ T-cell"
lymphoid$CT_firstpass[lymphoid$RNA_snn_res.0.25 == "3"] <- "B cell"
lymphoid$CT_firstpass[lymphoid$sub4 == "4_0"] <- "NK cell"
lymphoid$CT_firstpass[lymphoid$sub4 == "4_1"] <- "Endothelial"
lymphoid$CT_firstpass[lymphoid$RNA_snn_res.0.25 == "5"] <- "Proliferating T-cell"
lymphoid$CT_firstpass[lymphoid$sub6 == "6_0"] <- "T-cell"
lymphoid$CT_firstpass[lymphoid$sub6 == "6_1"] <- "Epithelial"
lymphoid$CT_firstpass[lymphoid$RNA_snn_res.0.25 == "7"] <- "Plasma/pDC"
DimPlot(lymphoid, group.by = "CT_firstpass", label = TRUE, raster = FALSE)

# Label cells on immune object
immune_cell_labels <- merge(myeloid, lymphoid)
immune$CT_firstpass <- immune_cell_labels$CT_firstpass
DimPlot(immune, group.by = "CT_firstpass", label = TRUE)

# saveRDS(immune, "/scratch/avannan/xenium_files/immune_firstpass_04172023.rds")
# saveRDS(lymphoid, "/scratch/avannan/xenium_files/lymphoid_firstpass_04172023.rds")
# lymphoid <- readRDS("/scratch/avannan/xenium_files/lymphoid_firstpass_04172023.rds")


###### MYELOID ----
# Re-UMAP myeloid cells
# Also include cells that seemed myeloid-like from the lymphoid group
myeloid <- subset(immune, subset = RNA_snn_res.0.2 %in% c(1, 3, 4, 7))
DefaultAssay(myeloid) <- "RNA"
myeloid <- NormalizeData(myeloid)
myeloid <- ScaleData(myeloid, features = rownames(myeloid))
myeloid <- RunPCA(myeloid, features = rownames(myeloid))

# Get PCs and make elbow plot
npcs <- min(get_pcs(myeloid))
ElbowPlot(myeloid, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 12

# Find neighbors, cluster and UMAP
myeloid <- FindNeighbors(myeloid, dims = 1:npcs)
myeloid <- RunUMAP(myeloid, dims = 1:npcs)
myeloid <- FindClusters(myeloid, resolution = c(0.3))
DimPlot(myeloid, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(myeloid, group.by = "sample")
DimPlot(myeloid, group.by = "sample", split.by = "sample", ncol = 6) + NoLegend()
DimPlot(myeloid, group.by = "sample_type")

myeloid <- SetIdent(myeloid, value = "RNA_snn_res.0.3")
FeaturePlot(myeloid, features = c("PECAM1", "EPCAM", "PTPRC", "CD14", "FCGR3A", "LYZ"))
DotPlot(myeloid, features = c("PECAM1", "EPCAM", "PTPRC", "CD14", "FCGR3A", "LYZ"))
VlnPlot(myeloid, features = c("PECAM1", "EPCAM", "PTPRC"), pt.size = 0)

# Look at top markers
myeloid <- SetIdent(myeloid, value = "RNA_snn_res.0.3")
my_markers <- FindAllMarkers(myeloid, only.pos = TRUE)
my_markers <- my_markers %>% mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# General myeloid markers, with possible epithelial contamination
FeaturePlot(myeloid, features = c("S100A8", "S100A9", "SLC25A37", "S100A12"))
VlnPlot(myeloid, features = c("S100A8", "S100A9", "SLC25A37", "S100A12"), pt.size = 0, ncol = 2)
FeaturePlot(myeloid, features = c("EPCAM", "SFTPC", "SCGB3A2", "SCGB1A1"))
VlnPlot(myeloid, features = c("EPCAM", "SFTPC", "SCGB3A2", "SCGB1A1"), pt.size = 0, ncol = 2)

# Lymphoid markers - none
VlnPlot(myeloid, features = c("NKG7", "GNLY", "CD3E", "CD4", "CD8A"), pt.size = 0)

# Myeloid markers
VlnPlot(myeloid, features = c("CCL18", "CCL22", "FCER1A", "FCER1G", "IL1B", 
                              "ITGAM", "ITGAX", "LYZ", "MRC1", "NFKB1", 
                              "PLIN2", "RHOA", "TREM2", "UQCRHL"), pt.size = 0)
# DC - CCL22+
FeaturePlot(myeloid, features = c("CD86", "CD1A", "CD1C", "CXCL9", "GPR183",
                                  "GZMB", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", 
                                  "MMP12"))
FeaturePlot(myeloid, features = "CCL22")

# Macrophage markers
FeaturePlot(myeloid, features = c("C1QC", "CD68", "FABP4", "IL1A", "MARCO", 
                                  "MCEMP1", "MS4A7", "SPP1"))
VlnPlot(myeloid, features = c("CD68", "MCEMP1", "MS4A7", "ITGAM", "ITGAX", 
                              "PPARG", "MRC1", "PLIN2"), pt.size = 0)

# Epithelial markers (and LYZ)
DotPlot(myeloid, features = c("SFTPD", "NAPSA", "PGC", "SFTPC", "LYZ"))
FeaturePlot(myeloid, features = c("SFTPD", "NAPSA", "PGC", "SFTPC", "LYZ"))

FeaturePlot(myeloid, features = c("S100A8", "S100A9", "PGC", "SFTPC", "LYZ"))

# Subcluster 0
myeloid <- FindSubCluster(myeloid, cluster = "0", subcluster.name = "sub0",
                          graph.name = "RNA_snn", resolution = 0.3)
myeloid <- SetIdent(myeloid, value = "sub0")
DimPlot(myeloid, group.by = "sub0", label = TRUE)
DotPlot(myeloid, features = c("SFTPD", "NAPSA", "PGC", "SFTPC", "LYZ", "S100A8", "S100A9", "PPARG"), group.by = "sub0")

# Subcluster 2
myeloid <- FindSubCluster(myeloid, cluster = "2", subcluster.name = "sub2",
                          graph.name = "RNA_snn", resolution = 0.1)
DimPlot(myeloid, group.by = "sub2", label = TRUE)
DotPlot(myeloid, features = c("SFTPD", "NAPSA", "PGC", "SFTPC", "LYZ", "S100A8", "S100A9", "PPARG"), group.by = "sub3")
sub3 <- subset(myeloid, subset = RNA_snn_res.0.3 == "3")
FeaturePlot(sub3, features = c("SFTPC", "S100A8"), blend = TRUE, pt.size = 0.0001, split.by = "sub3")

# Subcluster 3
myeloid <- FindSubCluster(myeloid, cluster = "3", subcluster.name = "sub3",
                           graph.name = "RNA_snn", resolution = 0.3)
myeloid <- SetIdent(myeloid, value = "sub3")
DimPlot(myeloid, group.by = "sub3", label = TRUE)
DotPlot(myeloid, features = c("SFTPD", "NAPSA", "PGC", "SFTPC", "LYZ", "S100A8", "S100A9", "PPARG"), group.by = "sub3")
sub3 <- subset(myeloid, subset = RNA_snn_res.0.3 == "3")
FeaturePlot(sub3, features = c("SFTPC", "S100A8"), blend = TRUE, pt.size = 0.0001, split.by = "sub3")

# Subcluster 4 - doesn't work well; 4 is mesenchymal/mast-like
myeloid <- FindSubCluster(myeloid, cluster = "4", subcluster.name = "sub4",
                          graph.name = "RNA_snn", resolution = 0.2)
myeloid <- SetIdent(myeloid, value = "sub4")
DimPlot(myeloid, group.by = "sub4", label = TRUE)
DotPlot(myeloid, features = c("KIT", "CPA3", "CD44", "CD69", "DDIT3", "TPSAB1",
                              "MEG3", "COL1A1", "FN1", "CCL2", "LUM", "DCN", "ELN",
                              "HAS2", "ITGAV"), group.by = "sub4")
sub4 <- subset(myeloid, subset = RNA_snn_res.0.3 == "4")


myeloid <- FindSubCluster(myeloid, cluster = "6", subcluster.name = "sub6",
                          graph.name = "RNA_snn", resolution = 0.2)
myeloid <- SetIdent(myeloid, value = "sub6")
DimPlot(myeloid, group.by = "sub6", label = TRUE)
myeloid <- FindSubCluster(myeloid, cluster = "6_1", subcluster.name = "sub6",
                          graph.name = "RNA_snn", resolution = 0.2)
DimPlot(myeloid, group.by = "sub6", label = TRUE)

# Subcluster 7
myeloid <- FindSubCluster(myeloid, cluster = "7", subcluster.name = "sub7",
                          graph.name = "RNA_snn", resolution = 0.2)
myeloid <- SetIdent(myeloid, value = "sub7")
DimPlot(myeloid, group.by = "sub7", label = TRUE)
DotPlot(myeloid, features = c("MUC5B", "SCGB1A1", "SCGB3A2", "MMP7", "LYZ",
                              "S100A8", "S100A9", "PPARG"), group.by = "sub7")

DimPlot(myeloid, group.by = "sub7", label = TRUE)

# Label cells
myeloid$CT_firstpass <- ""
myeloid$CT_firstpass[myeloid$sub0 == "0_2"] <- "Endothelial"
myeloid$CT_firstpass[myeloid$sub0 %in% c("0_0", "0_1", "0_3", "0_4")] <- "Myeloid"
myeloid$CT_firstpass[myeloid$RNA_snn_res.0.3 == "1"] <- "Macrophage"
myeloid$CT_firstpass[myeloid$RNA_snn_res.0.3 %in% c("2", "4")] <- "Mast"
myeloid$CT_firstpass[myeloid$RNA_snn_res.0.3 == "3"] <- "Epithelial"
myeloid$CT_firstpass[myeloid$RNA_snn_res.0.3 == "5"] <- "DC - CCL22+"
myeloid$CT_firstpass[myeloid$RNA_snn_res.0.3 == "6"] <- "Monocyte"
myeloid$CT_firstpass[myeloid$sub7 == "7_0"] <- "Myeloid"
myeloid$CT_firstpass[myeloid$sub7 == "7_1"] <- "Epithelial"
myeloid$CT_firstpass[myeloid$RNA_snn_res.0.3 == "8"] <- "Monocyte"
DimPlot(myeloid, group.by = "CT_firstpass", label = TRUE, raster = FALSE)



## MESENCHYMAL ----
mesenchymal <- subset(merged_spatial, subset = RNA_snn_res.0.5 %in% c(0, 2, 5, 9, 12, 15, 19, 20, 22))

DefaultAssay(mesenchymal) <- "RNA"
mesenchymal <- NormalizeData(mesenchymal)
mesenchymal <- ScaleData(mesenchymal, features = rownames(mesenchymal))
mesenchymal <- RunPCA(mesenchymal, features = rownames(mesenchymal))

# Get PCs and make elbow plot
npcs <- min(get_pcs(mesenchymal))
ElbowPlot(mesenchymal, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 16

# Find neighbors, cluster and UMAP
mesenchymal <- FindNeighbors(mesenchymal, dims = 1:npcs)
mesenchymal <- RunUMAP(mesenchymal, dims = 1:npcs)
mesenchymal <- FindClusters(mesenchymal, resolution = c(0.1, 0.2, 0.3))
DimPlot(mesenchymal, group.by = "RNA_snn_res.0.3", label = TRUE)

mesenchymal <- SetIdent(mesenchymal, value = "RNA_snn_res.0.3")
mes_markers <- FindAllMarkers(mesenchymal, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

mesenchymal <- SetIdent(mesenchymal, value = "RNA_snn_res.0.3")
FeaturePlot(mesenchymal, features = c("PECAM1", "EPCAM", "PTPRC", "LUM", "DCN", "ACTA2"))
DotPlot(mesenchymal, features = c("PECAM1", "EPCAM", "PTPRC", "LUM", "DCN", "ACTA2"))
VlnPlot(mesenchymal, features = c("PECAM1", "EPCAM", "PTPRC", "LUM", "DCN", "ACTA2"), pt.size = 0)

# Looking at lymphoid markers
VlnPlot(mesenchymal, features = c("CD3E", "CD8A", "GNLY", "NKG7"), pt.size = 0, ncol = 2)
FeaturePlot(mesenchymal, features = c("CD3E", "CD8A", "GNLY", "NKG7"))

# Epithelial markers
FeaturePlot(mesenchymal, features = c("SFTPC", "SCGB3A2", "SCGB1A1", "MUC5B", "FOXJ1", "NAPSA"))
VlnPlot(mesenchymal, features = c("SFTPC", "SCGB3A2", "SCGB1A1", "MUC5B", "FOXJ1", "NAPSA"), pt.size = 0)

# Myeloid markers
FeaturePlot(mesenchymal, features = c("MKI67", "S100A8", "S100A9", "CD14", "FCGR3A", "FCN3"))

# SMC/FB markers (plus immune markers that are in the ACTA2+ cluster)
FeaturePlot(mesenchymal, features = c("ACTA2", "AXL", "LGR6", "SPARCL1"))
DotPlot(mesenchymal, features = c("ACTA2", "AXL", "LGR6", "SPARCL1", "MKI67", 
                                  "CD14", "S100A8", "S100A9", "CSPG4", "PDGFRB"))

# Other mesenchymal markers
FeaturePlot(mesenchymal, features = c("MFAP5", "PI16", "SFRP2", "SFRP4")) # Adventitial FB
FeaturePlot(mesenchymal, features = c("CSPG4", "PDGFRB")) # Pericytes

FeaturePlot(mesenchymal, features = c("CCL2", "COL1A1", "COL1A2", "COL3A1", 
                                      "CXCL14", "DCN", "ELN", "FAP", "FGF2", 
                                      "FGF7", "FN1", "HAS2", "ITGAV", "LUM", 
                                      "MEG3", "PDGFRA", "SNAI2", "TGFB3", "UGDH"))
DotPlot(mesenchymal, features = c("CCL2", "COL1A1", "COL1A2", "COL3A1", 
                                  "CXCL14", "DCN", "ELN", "FAP", "FGF2", 
                                  "FGF7", "FN1", "HAS2", "ITGAV", "LUM", 
                                  "MEG3", "PDGFRA", "SNAI2", "TGFB3", "UGDH"))

# Subcluster 5
mesenchymal <- FindSubCluster(mesenchymal, cluster = "5", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub5")
mesenchymal <- SetIdent(mesenchymal, value = "sub5")
DimPlot(mesenchymal, group.by = "sub5", label = TRUE)
mes_sub5 <- subset(mesenchymal, subset = RNA_snn_res.0.3 == "5")
DimPlot(mes_sub5, group.by = "sub5", label = TRUE)
mes_sub5 <- SetIdent(mes_sub5, value = "sub5")
mes_sub5_markers <- FindAllMarkers(mes_sub5, only.pos = TRUE) %>%
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(mesenchymal, features = c("FABP4", "MCEMP1",  "CD52", "SFTPD",
                                  "CCL2", "SPP1", "GPR183", "HMOX1", "HIF1A",
                                  "MKI67", "TOP2A",
                                  "ACTA2", "CSPG4"), group.by = "sub5") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(mesenchymal, features = c("FN1", "ELN", "COL1A1", "COL1A2", 
                                  "FABP4", "MCEMP1",  "CD52", "SFTPD",
                                  "CCL2", "SPP1", "GPR183", "HMOX1", "HIF1A",
                                  "MKI67", "TOP2A",
                                  "ACTA2"), pt.size = 0, group.by = "sub5")
# 5 is all immune, despite high ACTA2 expression

# Subcluster 7
mesenchymal <- FindSubCluster(mesenchymal, cluster = "7", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub7")
DimPlot(mesenchymal, group.by = "sub7", label = TRUE)
mesenchymal <- SetIdent(mesenchymal, value = "sub7")
mes_sub7 <- subset(mesenchymal, subset = RNA_snn_res.0.3 == "7")
mes_sub7 <- SetIdent(mes_sub7, value = "sub7")
mes_sub7_markers <- FindAllMarkers(mes_sub7, only.pos = TRUE) %>%
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(mesenchymal, features = c("PECAM1", "EPAS1", "CD34", "SPARCL1", "KDR", 
                                  "GNG11", "PLVAP", "ACKR1", "FCN1", "S100A8", 
                                  "S100A9", "MCEMP1", "ITGAM", "PTPRC", 
                                  "FCER1G", "COL1A1", "COL3A1", "FN1", "LUM", "MEG3"), group.by = "sub7") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(mesenchymal, features = c("PECAM1", "EPAS1", "CD34", "SPARCL1", "KDR", 
                                  "GNG11", "S100A8", "S100A9", "MCEMP1", "ITGAM", "PTPRC", 
                                  "FCER1G", "COL1A1", "COL3A1", "FN1", "LUM", "MEG3"),
        pt.size = 0, group.by = "sub7")

# Subcluster 0
mesenchymal <- FindSubCluster(mesenchymal, cluster = "0", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub0")
DimPlot(mesenchymal, group.by = "sub0", label = TRUE)
mes_sub0 <- subset(mesenchymal, subset = RNA_snn_res.0.3 == "0")
mes_sub0 <- SetIdent(mes_sub0, value = "sub0")
mes_sub0_markers <- FindAllMarkers(mes_sub0, only.pos = TRUE) %>%
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(mesenchymal, features = c("SFTPC", "NAPSA", "PGC", "SFTPD", "COL1A1", 
                                  "COL1A2", "COL3A1", "FN1", "ELN", "ACTA2", 
                                  "CD14", "FCER1G", "IFIT1", "IFIT3", "OAS2", 
                                  "OAS3", "ITGAX", "SFRP2", "STAT1"), group.by = "sub0") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(mesenchymal, features = c("SFTPC", "NAPSA", "PGC",  "COL1A1", 
                                  "COL1A2", "COL3A1", "FN1", "ELN", "ACTA2", 
                                  "CD14", "FCER1G", "IFIT1", "IFIT3", 
                                  "OAS3", "ITGAX", "SFRP2"),
        pt.size = 0, group.by = "sub0")

# General fibroblast genes and other markers
DotPlot(mesenchymal, features = c("PECAM1", "PTPRC", "EPCAM", "CCL2", "COL1A1",
                                  "COL1A2", "COL3A1", "DCN", "ELN", "FGF2", 
                                  "FGF7", "FN1", "HAS2", "ITGAV", "LUM", "MEG3",
                                  "SFTPC", "NAPSA", "SFTPD", "KRT5", "KRT17", 
                                  "RTKN2", "S100A8", "S100A9", "ACKR1", "PLVAP"))

# Label cell types
mesenchymal$CT_firstpass <- ""
mesenchymal$CT_firstpass[mesenchymal$sub0 %in% c("0_0", "0_1", "0_2")] <- "Immune"
mesenchymal$CT_firstpass[mesenchymal$sub0 == "0_3"] <- "Epithelial"
mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 == "1"] <- "Mesenchymal"
mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 == "2"] <- "Adventitial FB"
mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 == "3"] <- "Mesenchymal"
mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 == "4"] <- "Immune"
mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 == "5"] <- "Immune"
mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 == "6"] <- "Epithelial"
mesenchymal$CT_firstpass[mesenchymal$sub7 == "7_0"] <- "Endothelial"
mesenchymal$CT_firstpass[mesenchymal$sub7 == "7_1"] <- "Immune"
mesenchymal$CT_firstpass[mesenchymal$sub7 == "7_2"] <- "Mesenchymal"
mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 == "8"] <- "Immune"
mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 == "9"] <- "Mesothelial"
DimPlot(mesenchymal, group.by = "CT_firstpass", label = TRUE, raster = FALSE)

mesenchymal$CT_firstpass <- og_mesenchymal$CT_firstpass





###### TRUE MESENCHYMAL ----
# mesenchymal <- subset(mesenchymal, subset = CT_firstpass %in% 
#                         c("Epithelial", "Immune", "Endothelial"), invert = TRUE)
# 
# DefaultAssay(mesenchymal) <- "RNA"
# mesenchymal <- NormalizeData(mesenchymal)
# mesenchymal <- ScaleData(mesenchymal, features = rownames(mesenchymal))
# mesenchymal <- RunPCA(mesenchymal, features = rownames(mesenchymal))
# 
# # Get PCs and make elbow plot
# npcs <- min(get_pcs(mesenchymal))
# ElbowPlot(mesenchymal, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
# npcs # 12
# 
# # Find neighbors, cluster and UMAP
# mesenchymal <- FindNeighbors(mesenchymal, dims = 1:npcs)
# mesenchymal <- RunUMAP(mesenchymal, dims = 1:npcs)
# mesenchymal <- FindClusters(mesenchymal, resolution = c(0.2, 0.3))
# DimPlot(mesenchymal, group.by = "RNA_snn_res.0.3", label = TRUE)
# 
# # FindMarkers
# mesenchymal <- SetIdent(mesenchymal, value = "RNA_snn_res.0.3")
# mes_markers <- FindAllMarkers(mesenchymal, only.pos = TRUE) %>% 
#   mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# 
# # Lineage markers
# FeaturePlot(mesenchymal, features = c("PECAM1", "EPCAM", "PTPRC", "LUM", "ACTA2", "MSLN"))
# 
# # Possible non-mesenchymal genes/clusters
# DotPlot(mesenchymal, features = c("COL1A2", "SPARCL1", "EPAS1", "ELN", "ITGAX", "TGFB3", "LGR6",
#                                   "POSTN", "CTHRC1", "COL1A1", "LUM", "COL3A1", "SFRP2",
#                                   "PTGDS", "DCN", "MEG3", "TGFB2", "ITGA3", "HAS2",
#                                   "ACTA2", "MS4A7", "PDGFRA", "HLA-DQB1", "CD14", 
#                                   "FCER1G", "CCL21", "PDGFRB", "IFIT1", "IFIT3", 
#                                   "OAS2", "OAS3", "STAT1", "GDF15", 
#                                   "LGR5", "WNT5A", "SCGB3A2", "SCGB1A1",
#                                   "CDK1", "TOP2A", "MKI67", "CENPF")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# FeaturePlot(mesenchymal, features = c("COL1A2", "SPARCL1", "ELN", "ITGAX", "TGFB3",
#                                       "CTHRC1", "FAP", "COL1A1", "LUM", "SFRP2",
#                                       "PTGDS", "DCN", "MEG3", "TGFB2", "ITGA3", 
#                                       "ACTA2", "MS4A7", "PDGFRA", "FCER1G", 
#                                       "CCL21", "PDGFRB", "IFIT1", "OAS2", 
#                                       "STAT1", "GDF15", "WNT5A", "SCGB3A2", "MKI67"), ncol = 5)
# VlnPlot(mesenchymal, features = c("COL1A2", "SPARCL1", "EPAS1", "ELN", "ITGAX", "TGFB3",
#                                   "POSTN", "CTHRC1", "COL1A1", "LUM", "COL3A1", "SFRP2",
#                                   "PTGDS", "DCN", "MEG3", "TGFB2", "ITGA3",
#                                   "ACTA2", "MS4A7", "FCER1G", "CCL21", "PDGFRB", 
#                                   "IFIT1", "IFIT3", "STAT1", "GDF15", "WNT5A", 
#                                   "SCGB3A2","TOP2A", "MKI67"), pt.size = 0, ncol = 6)
# 
# # Subcluster 3
# sub3 <- subset(mesenchymal, subset = RNA_snn_res.0.3 == "3")
# FeaturePlot(sub3, features = c("ACTA2", "MS4A7", "CD14", "C1QC"))
# FeaturePlot(sub3, features = c("ACTA2", "MS4A7"), blend = TRUE)
# mesenchymal <- FindSubCluster(mesenchymal, cluster = "3", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub3")
# DimPlot(mesenchymal, group.by = "sub3", label = TRUE)
# mes_sub3 <- subset(mesenchymal, subset = RNA_snn_res.0.3 == "3")
# DimPlot(mes_sub3, group.by = "sub3", label = TRUE)
# mes_sub3 <- SetIdent(mes_sub3, value = "sub3")
# mes_sub3_markers <- FindAllMarkers(mes_sub3, only.pos = TRUE) %>%
#   mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# VlnPlot(mes_sub3, features = c("COL1A2", "SPARCL1", "EPAS1", "ELN", "ITGAX", "TGFB3",
#                                "CTHRC1", "COL1A1", "LUM", "COL3A1", "PTGDS", 
#                                "DCN", "ACTA2", "MS4A7", "FCER1G", "CCL21", "PDGFRB",
#                                "PDGFRA"), pt.size = 0, ncol = 6)
# 
# # SMC/FB/pericyte
# FeaturePlot(mesenchymal, features = c("ACTA2", "AXL", "LGR6", "SPARCL1", "PDGFRA"))
# 
# # Label cell types
# mesenchymal$CT_firstpass <- ""
# mesenchymal$CT_firstpass[mesenchymal$RNA_snn_res.0.3 %in% c("0", "1", "2", "4", "5", "6")] <- "Mesenchymal"
# mesenchymal$CT_firstpass[mesenchymal$sub3 == "3_0"] <- "Mesenchymal"
# mesenchymal$CT_firstpass[mesenchymal$sub3 == "3_1"] <- "Epithelial"
# DimPlot(mesenchymal, group.by = "CT_firstpass", label = TRUE, raster = FALSE)
# 
# # Label cell types
# og_mesenchymal$CT_firstpass <- mesenchymal$CT_firstpass
# og_mesenchymal$CT_firstpass[og_mesenchymal$sub0 == "0_2"] <- "Immune"
# og_mesenchymal$CT_firstpass[og_mesenchymal$sub0 == "0_3"] <- "Epithelial"
# og_mesenchymal$CT_firstpass[og_mesenchymal$RNA_snn_res.0.3 == "4"] <- "Immune"
# og_mesenchymal$CT_firstpass[og_mesenchymal$RNA_snn_res.0.3 == "5"] <- "Immune"
# og_mesenchymal$CT_firstpass[og_mesenchymal$RNA_snn_res.0.3 == "6"] <- "Epithelial"
# og_mesenchymal$CT_firstpass[og_mesenchymal$sub7 == "7_0"] <- "Endothelial"
# og_mesenchymal$CT_firstpass[og_mesenchymal$sub7 == "7_1"] <- "Immune"
# og_mesenchymal$CT_firstpass[og_mesenchymal$RNA_snn_res.0.3 == "8"] <- "Immune"
# DimPlot(og_mesenchymal, group.by = "CT_firstpass", label = TRUE, raster = FALSE, repel = TRUE)
# 
# # saveRDS(mesenchymal, "/scratch/avannan/xenium_files/mesenchymal_firstpass_04182023.rds")
# mesenchymal <- readRDS("/scratch/avannan/xenium_files/mesenchymal_firstpass_04182023.rds")
# # saveRDS(og_mesenchymal, "/scratch/avannan/xenium_files/og_mesenchymal_firstpass_04182023.rds")


## LABEL FIRSTPASS CELL TYPES ON MAIN OBJECT ----
full_cell_labels <- merge(epithelial, y = list(og_endothelial, immune, og_mesenchymal))
merged_spatial$CT_firstpass <- ""
merged_spatial$CT_firstpass <- full_cell_labels$CT_firstpass
DimPlot(merged_spatial, group.by = "CT_firstpass", label = TRUE, repel = TRUE)


## SECOND PASS ANNOTATIONS ----
#### EPITHELIAL ----
second_epi <- subset(merged_spatial, 
                     subset = CT_firstpass %in% c("AT1", "AT2", "Basal", "Ciliated", "Differentiating Ciliated",
                                                  "Epithelial", "KRT5-/KRT17+", "KRT5-low Basal", "MUC5B+ Secretory",
                                                  "Proliferating Epithelial", "SCGB+ Secretory"))

DefaultAssay(second_epi) <- "RNA"
second_epi <- NormalizeData(second_epi)
second_epi <- ScaleData(second_epi, features = rownames(second_epi))
second_epi <- RunPCA(second_epi, features = rownames(second_epi))

# Get PCs and make elbow plot
npcs <- min(get_pcs(second_epi))
ElbowPlot(second_epi, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 14

# Find neighbors, cluster and UMAP
second_epi <- FindNeighbors(second_epi, dims = 1:npcs)
second_epi <- RunUMAP(second_epi, dims = 1:npcs)
second_epi <- FindClusters(second_epi, resolution = c(0.3, 0.4))
DimPlot(second_epi, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(second_epi, group.by = "CT_firstpass", label = TRUE)

second_epi <- SetIdent(second_epi, value = "RNA_snn_res.0.3")
epi_markers <- FindAllMarkers(second_epi, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(second_epi, features = epithelial_features, group.by = "CT_firstpass") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Look at cell distribution
round(table(second_epi$sample, second_epi$RNA_snn_res.0.3)/rowSums(table(second_epi$sample, second_epi$RNA_snn_res.0.3)), 3)*100

# Subcluster 0
sub0 <- subset(second_epi, subset = RNA_snn_res.0.3 == "0")
DimPlot(sub0, reduction = "sp")
second_epi <- FindSubCluster(second_epi, cluster = "0", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub0")
DimPlot(second_epi, group.by = "sub0", label = TRUE)
second_epi <- SetIdent(second_epi, value = "sub0")
second_epi <- FindSubCluster(second_epi, cluster = "0_1", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "sub0_1")
DimPlot(second_epi, group.by = "sub0_1", label = TRUE)
DotPlot(second_epi, features = c("SFTPC", "KRT8", "KRT17", "KRT5"), group.by = "sub0")
DotPlot(sub0, features = c("SFTPC", "KRT8", "KRT17", "KRT5"), group.by = "sub0")
VlnPlot(sub0, features = c("SFTPC", "KRT8", "KRT17", "KRT5"), group.by = "sub0", pt.size = 0)
sub0 <- subset(second_epi, subset = RNA_snn_res.0.3 == "0")
DimPlot(sub0, reduction = "sp", group.by = "sub0")
VlnPlot(second_epi, features = c("SFTPC", "KRT8", "KRT17", "KRT5"), group.by = "sub0_1", pt.size = 0, ncol = 2)
FeaturePlot(sub0, features = c("KRT5", "KRT17", "KRT8", "SFTPC"), pt.size = 0.0001, split.by = "sub0_1")

# Subcluster 2
sub2 <- subset(second_epi, subset = RNA_snn_res.0.3 == "2")
DimPlot(sub2, reduction = "sp")

# Subcluster 3
second_epi <- FindSubCluster(second_epi, cluster = "3", graph.name = "RNA_snn", resolution = 0.25, subcluster.name = "sub3")
DimPlot(second_epi, group.by = "sub3", label = TRUE)
sub3 <- subset(second_epi, subset = RNA_snn_res.0.3 == "3")
FeaturePlot(sub3, features = c("C20orf85", "SCGB1A1"), blend = TRUE, split.by = "sub3")
VlnPlot(sub3, features = c("C20orf85", "SCGB1A1", "MUC5B", "SCGB3A2"), group.by = "sub3", pt.size = 0, ncol = 2)
DotPlot(sub3, features = epithelial_features, group.by = "sub3") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
sub3 <- SetIdent(sub3, value = "sub3")
epi3_markers <- FindAllMarkers(sub3, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(second_epi, features = epithelial_features, group.by = "CT_firstpass") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Find proliferating cells
# Subcluster 7
second_epi <- FindSubCluster(second_epi, cluster = "7", graph.name = "RNA_snn", resolution = 0.15, subcluster.name = "sub7")
DimPlot(second_epi, group.by = "sub7", label = TRUE)
DotPlot(second_epi, features = c("TOP2A", "MKI67"), group.by = "sub7")
VlnPlot(second_epi, features = c("TOP2A", "MKI67"), group.by = "sub7", pt.size = 0)
sub7 <- subset(second_epi, subset = RNA_snn_res.0.3 == "7")
FeaturePlot(sub7, features = "TOP2A", split.by = "sub7")

# Subcluster immune-like/AT2 (4)
FeaturePlot(second_epi, features = "PTPRC")
DotPlot(second_epi, features = c("CD3E", "CD4", "CD8A", "CD8B", "IL7R", "GNLY", "NKG7"))
second_epi <- FindSubCluster(second_epi, cluster = "4", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub4")
DimPlot(second_epi, group.by = "sub4", label = TRUE)
sub4 <- subset(second_epi, subset = RNA_snn_res.0.3 == "4")
sub4 <- SetIdent(sub4, value = "sub4")
epi4_markers <- FindAllMarkers(sub4, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
VlnPlot(second_epi, features = c("CD3E", "MS4A7", "S100A8"), group.by = "sub4", pt.size = 0)

# Label cells
second_epi$CT_secondpassA <- ""
second_epi$CT_secondpassA[second_epi$RNA_snn_res.0.3 == "0"] <- "Transitional AT2"
# second_epi$CT_secondpassA[second_epi$sub0_1 == "0_1_0"] <- "KRT8+/SFTPC+"
second_epi$CT_secondpassA[second_epi$sub0_1 == "0_1_3"] <- "KRT5-/KRT17+"
second_epi$CT_secondpassA[second_epi$RNA_snn_res.0.3 == "1"] <- "Ciliated"
second_epi$CT_secondpassA[second_epi$RNA_snn_res.0.3 == "2"] <- "Alveolar FBs/Fibrotic AT2"
second_epi$CT_secondpassA[second_epi$sub3 %in% c("3_0", "3_3")] <- "MUC5B+ Secretory"
second_epi$CT_secondpassA[second_epi$sub3 %in% c("3_1", "3_2")] <- "SCGB+ Secretory"
second_epi$CT_secondpassA[second_epi$sub3 == "3_4"] <- "Differentiating Ciliated"
second_epi$CT_secondpassA[second_epi$sub4 == "4_0"] <- "AT2"
second_epi$CT_secondpassA[second_epi$sub4 == "4_1"] <- "Macrophages"
second_epi$CT_secondpassA[second_epi$RNA_snn_res.0.3 == "6"] <- "AT1"
second_epi$CT_secondpassA[second_epi$RNA_snn_res.0.3 %in% c("5", "7")] <- "Basal-like"
second_epi$CT_secondpassA[second_epi$sub7 == "7_3"] <- "Proliferating Epithelial"
second_epi$CT_secondpassA[second_epi$RNA_snn_res.0.3 == "8"] <- "Proliferating Epithelial"
second_epi$CT_secondpassA[second_epi$RNA_snn_res.0.3 == "9"] <- "Interstitial Macrophages"
second_epi$CT_secondpassA[second_epi$RNA_snn_res.0.3 == "10"] <- "T-cells"
DimPlot(second_epi, group.by = "CT_secondpassA", label = TRUE)
DimPlot(second_epi, group.by = "CT_firstpass", label = TRUE)
DotPlot(second_epi, features = c(epithelial_features, "COL1A1", "COL3A1", "FN1", "MKI67", "TOP2A", 
                                 "CD3E", "CD8A", "MS4A7", "MARCO", "S100A8", "S100A9", "S100A12"), 
        group.by = "CT_secondpassA") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(second_epi$CT_firstpass, second_epi$CT_secondpassA)
#saveRDS(second_epi, "/scratch/avannan/xenium_files/second_epi_secondpass_04212023.rds")

# Look at potential alveolar FBs and compare to AT1/2 cells
at_cells <- subset(second_epi, subset = CT_secondpassA %in% c("AT1", "Transitional AT2", "AT2", "Alveolar FBs/Fibrotic AT2"))
VlnPlot(at_cells, features = c("SFTPC", "NAPSA", "AGER", "RTKN2", "COL1A1", "COL3A1", "FN1"), pt.size = 0)
# Not sure if these are AT2 fibrotic or alveolar FBs


###### AIRWAY ----
second_air <- subset(second_epi, 
                     subset = CT_secondpassA %in% c("Ciliated", "MUC5B+ Secretory", 
                                                    "SCGB+ Secretory", "Differentiating Ciliated",
                                                    "Basal-like", "KRT5-/KRT17+"))
DefaultAssay(second_air) <- "RNA"
second_air <- NormalizeData(second_air)
second_air <- ScaleData(second_air, features = rownames(second_air))
second_air <- RunPCA(second_air, features = rownames(second_air))

# Get PCs and make elbow plot
npcs <- min(get_pcs(second_air))
ElbowPlot(second_air, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 13

# Find neighbors, cluster and UMAP
second_air <- FindNeighbors(second_air, dims = 1:npcs)
second_air <- RunUMAP(second_air, dims = 1:npcs)
second_air <- FindClusters(second_air, resolution = c(0.2, 0.3))
DimPlot(second_air, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(second_epi, group.by = "sample")
DimPlot(second_air, group.by = "CT_firstpass", label = TRUE)
DimPlot(second_air, group.by = "CT_secondpassA", label = TRUE)

air_markers <- FindAllMarkers(second_air, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
table(second_air$RNA_snn_res.0.3, second_air$CT_secondpassA)
table(second_air$RNA_snn_res.0.3, second_air$CT_firstpass)

second_air <- SetIdent(second_air, value = "RNA_snn_res.0.3")
DotPlot(second_air, features = c(epithelial_features, "COL1A1", "COL3A1", "FN1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(second_air, features = c("KRT5", "KRT8", "KRT17", "TP63", "SFTPC", 
                                 "NAPSA", "AGER", "RTKN2", "FOXJ1", "C20orf85",
                                 "MUC5B", "MUC5AC", "SCGB1A1", "SCGB3A2", "SPINK1",
                                 "SOX4", "MMP7", "IFIT1", "IFIT3", "OAS2", "OAS3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Look at ciliated/MUC5B/SCGB genes
FeaturePlot(second_air, features = c("FOXJ1", "MUC5B", "MUC5AC", "SCGB3A2", "SCGB1A1", "KRT5"))
VlnPlot(second_air, features = c("FOXJ1", "MUC5B", "SCGB3A2", "SCGB1A1"), pt.size = 0, ncol = 2)

# Subcluster 2 - separate MUC5B+, SCGB+, and Ciliated
second_air <- FindSubCluster(second_air, cluster = "2", graph.name = "RNA_snn", resolution = 0.3, subcluster.name = "sub2")
DimPlot(second_air, group.by = "sub2", label = TRUE)
sub2 <- subset(second_air, subset = RNA_snn_res.0.3 == "2")
FeaturePlot(sub2, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "MUC5B"), split.by = "sub2")
VlnPlot(sub2, features = c("SCGB3A2", "SCGB1A1", "MUC5B", "C20orf85", "FOXJ1"), group.by = "sub2", pt.size = 0)
second_air <- SetIdent(second_air, value = "sub2")
second_air <- FindSubCluster(second_air, cluster = "2_0", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub2_0")
DimPlot(second_air, group.by = "sub2_0", label = TRUE)
sub2 <- subset(second_air, subset = RNA_snn_res.0.3 == "2")
FeaturePlot(sub2, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "MUC5B", "C20orf85"), split.by = "sub2_0")
VlnPlot(sub2, features = c("SCGB3A2", "SCGB1A1", "MUC5B", "C20orf85", "FOXJ1", "SFTPC"), group.by = "sub2_0", pt.size = 0)
sub2 <- SetIdent(sub2, value = "sub2_0")
sub2_markers <- FindAllMarkers(sub2, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")


# Subcluster 4 (Ciliated and SCGB+)
second_air <- SetIdent(second_air, value = "RNA_snn_res.0.3")
second_air <- FindSubCluster(second_air, cluster = "4", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub4")
VlnPlot(second_air, features = c("C20orf85", "SCGB1A1", "MUC5B", "SCGB3A2"), group.by = "sub4", pt.size = 0, ncol = 2)
VlnPlot(sub4, features = c("C20orf85", "SCGB1A1", "MUC5B", "SCGB3A2"), group.by = "sub4", pt.size = 0, ncol = 2)


second_air$CT_secondpassB <- ""
second_air$CT_secondpassB[second_air$RNA_snn_res.0.3 == "0"] <- "Basal"
second_air$CT_secondpassB[second_air$RNA_snn_res.0.3 == "1"] <- "Ciliated"
second_air$CT_secondpassB[second_air$sub2_0 == "2_0_0"] <- "SCGB+ Secretory"
second_air$CT_secondpassB[second_air$sub2_0 %in% c("2_1", "2_2", "2_4")] <- "MUC5B+ Secretory"
second_air$CT_secondpassB[second_air$sub2_0 %in% c("2_0_1", "2_3")] <- "Differentiating Ciliated"
second_air$CT_secondpassB[second_air$RNA_snn_res.0.3 == "3"] <- "Ciliated"
second_air$CT_secondpassB[second_air$RNA_snn_res.0.3 == "4"] <- "SCGB+ Secretory"
second_air$CT_secondpassB[second_air$RNA_snn_res.0.3 == "5"] <- "Basal (Interferon High)"
second_air$CT_secondpassB[second_air$RNA_snn_res.0.3 == "6"] <- "KRT5-/KRT17+"
second_air$CT_secondpassB[second_air$RNA_snn_res.0.3 == "7"] <- "MUC5B+ Secretory (Interferon High)"
DimPlot(second_air, group.by = "CT_secondpassB", label = TRUE)
DimPlot(second_air, group.by = "CT_secondpassA", label = TRUE)
DimPlot(second_air, group.by = "CT_firstpass", label = TRUE)
table(second_air$CT_secondpassB, second_air$CT_secondpassA)
DotPlot(second_air, features = epithelial_features, group.by = "CT_secondpassB") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###### ALVEOLAR ----
second_alv <- subset(second_epi,  subset = CT_secondpassA %in% c("Alveolar FBs/Fibrotic AT2", "AT1", 
                                                                 "Transitional AT2", "AT2"))

DefaultAssay(second_alv) <- "RNA"
second_alv <- NormalizeData(second_alv)
second_alv <- ScaleData(second_alv, features = rownames(second_alv))
second_alv <- RunPCA(second_alv, features = rownames(second_alv))

# Get PCs and make elbow plot
npcs <- min(get_pcs(second_alv))
ElbowPlot(second_alv, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs #

# Find neighbors, cluster and UMAP
second_alv <- FindNeighbors(second_alv, dims = 1:npcs)
second_alv <- RunUMAP(second_alv, dims = 1:npcs)
second_alv <- FindClusters(second_alv, resolution = c(0.2, 0.3))
DimPlot(second_alv, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(second_alv, group.by = "CT_firstpass", label = TRUE)
DimPlot(second_alv, group.by = "sample")
DimPlot(second_alv, group.by = "sample_type")

second_alv <- SetIdent(second_alv, value = "RNA_snn_res.0.3")
alv_markers <- FindAllMarkers(second_alv, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

DotPlot(second_alv, features = epithelial_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(second_alv, features = mesenchymal_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(second_alv, features = immune_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
FeaturePlot(second_alv, features = c("SFTPC", "NAPSA", "AGER", "RTKN2", "CEACAM6", 
                                     "KRT8", "KRT15", "SCGB3A2", "SCGB1A1", "KRT5", "KRT17"))

# Subcluster 4
second_alv <- FindSubCluster(second_alv, cluster = "4", graph.name = "RNA_snn", resolution = 0.3, subcluster.name = "sub4")
DimPlot(second_alv, group.by = "sub4", label = TRUE)
sub4 <- subset(second_alv, subset = RNA_snn_res.0.3 == "4")
FeaturePlot(sub4, features = c("SFTPC", "SCGB3A2"), blend = TRUE, split.by = "sub4")
FeaturePlot(sub4, features = c("SFTPC", "KRT15"), blend = TRUE, split.by = "sub4")

FeaturePlot(second_alv, features = c("MUC5B", "MUC5AC"))
FeaturePlot(second_alv, features = mesenchymal_features)
FeaturePlot(second_alv, features = c("SOX4", "SOX9"))

# AT2, AT1, and Transitional AT2 markers
FeaturePlot(second_alv, features = c("SFTPC", "SCGB3A2", "RTKN2", "CEACAM5", "KRT8", "SOX4"))

# Markers for 1 and 6
FeaturePlot(second_alv, features = c("TGFB2", "ITGA3", "PTGDS", "LUM", "TGFB3", "ELN"))

# Look at AT2 and fibrotic clusters
second_alv <- FindSubCluster(second_alv, cluster = "1", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub1")
DimPlot(second_alv, group.by = "sub1", label = TRUE)
at2_fb_sub <- subset(second_alv, subset = RNA_snn_res.0.3 %in% c("3", "1", "6"))
DimPlot(at2_fb_sub, reduction = "sp", group.by = "sub1")
sub1 <- subset(at2_fb_sub, subset = RNA_snn_res.0.3 == "1")
DimPlot(subset(sub1, subset = sub1 == "1_0"), reduction = "sp", group.by = "sub1") # 1_1 appears to be alveolar fibroblasts
sub1 <- SetIdent(sub1, value = "sub1")
sub1_markers <- FindAllMarkers(sub1, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
VlnPlot(second_alv, group.by = "sub1", pt.size = 0, features = c("SFTPC", "ELN", "FN1", "NAPSA", "PGC", "SCGB3A2"))
DimPlot(subset(at2_fb_sub, subset = RNA_snn_res.0.3 %in% c("1", "3") & sample == "VUILD115"), reduction = "sp", group.by = "sub1", cols = c("red", "green", "blue"))
FeaturePlot(subset(at2_fb_sub, subset = RNA_snn_res.0.3 %in% c("1", "3") & sample == "VUILD115"), reduction = "sp", features = c("TGFB2", "SCGB3A2"))
sub6 <- subset(at2_fb_sub, subset = RNA_snn_res.0.3 == "6")
DimPlot(sub6, reduction = "sp")
# Difference between 1 and 6
sub1_6 <- subset(second_alv, subset = RNA_snn_res.0.3 %in% c("1", "6"))
sub1_6 <- SetIdent(sub1_6, value = "sub1")
sub1_6_markers <- FindAllMarkers(sub1_6, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Examine FB markers
FeaturePlot(second_alv, features = c("CTHRC1", "FAP", "HAS1", "PLIN2", "COL1A2", "PI16", "SFRP2"))

# Look at healthy vs. disease AT2
at2_sub <- subset(second_alv, subset = RNA_snn_res.0.3 %in% c("0", "3"))
DimPlot(at2_sub, reduction = "sp", split.by = "RNA_snn_res.0.3")
at2_markers <- FindMarkers(second_alv, ident.1 = "0", ident.2 = "3") %>%
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Subcluster 2
second_alv <- FindSubCluster(second_alv, cluster = "2", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub2")
DimPlot(second_alv, group.by = "sub2", label = TRUE) 
DotPlot(second_alv, features = immune_features, group.by = "sub2") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(second_alv, features = c("S100A8", "RTKN2"), group.by = "sub2", pt.size = 0)
sub2 <- subset(second_alv, subset = RNA_snn_res.0.3 == "2")
DimPlot(subset(sub2, subset = sample == "THD0008"), reduction = "sp", group.by = "sub2", cols = c("grey", "red", "grey", "grey"))

# Label cells
second_alv$CT_secondpassB <- ""
second_alv$CT_secondpassB[second_alv$RNA_snn_res.0.3 == "0"] <- "AT2 (Healthy)"
second_alv$CT_secondpassB[second_alv$sub1 == "1_0"] <- "Alveolar FBs (Disease)"
second_alv$CT_secondpassB[second_alv$sub1 == "1_1"] <- "Alveolar FBs (Healthy)"
second_alv$CT_secondpassB[second_alv$RNA_snn_res.0.3 == "2"] <- "AT1"
second_alv$CT_secondpassB[second_alv$sub2 == "2_1"] <- "Interstitial Macrophages"
second_alv$CT_secondpassB[second_alv$RNA_snn_res.0.3 == "3"] <- "AT2 (Disease)"
second_alv$CT_secondpassB[second_alv$RNA_snn_res.0.3 == "4"] <- "Transitional AT2"
second_alv$CT_secondpassB[second_alv$RNA_snn_res.0.3 == "5"] <- "AT2 (Interferon High)"
second_alv$CT_secondpassB[second_alv$RNA_snn_res.0.3 == "6"] <- "Fibroblasts"
DimPlot(second_alv, group.by = "CT_firstpass", label = TRUE, raster = FALSE)
DimPlot(second_alv, group.by = "CT_secondpassA", label = TRUE, raster = FALSE)
DimPlot(second_alv, group.by = "CT_secondpassB", label = TRUE, raster = FALSE)
table(second_alv$CT_secondpassA, second_alv$CT_secondpassB)
#saveRDS(second_alv, "/scratch/avannan/xenium_files/second_alv_secondpass_04252023.rds")

# Add annotations to epithelial object
second_alv$CT_secondpass_full <- second_alv$CT_secondpassB
second_air$CT_secondpass_full <- second_air$CT_secondpassB
second_epi$CT_secondpass_full <- second_epi$CT_secondpassA
second_epi_sub <- subset(second_epi, cells = c(colnames(second_air), colnames(second_alv)), invert = TRUE)
DimPlot(second_epi_sub, group.by = "CT_secondpassA")
second_epi_sub$CT_secondpass_full <- second_epi_sub$CT_secondpassA
full_epi_labels <- merge(second_air, y = list(second_alv, second_epi_sub))
second_epi$CT_secondpass_full <- full_epi_labels$CT_secondpass_full
second_epi$CT_secondpass_full[second_epi$CT_secondpass_full == "Interstitial Macrophages"] <- "Interstitial Macrophages"
DimPlot(second_epi, group.by = "CT_secondpass_full", repel = TRUE, label = TRUE)


#### ENDOTHELIAL ----
second_endo <- subset(merged_spatial, 
                      subset = CT_firstpass %in% c("Arteriole", "Capillary", "Endothelial", "Venous"))

DefaultAssay(second_endo) <- "RNA"
second_endo <- NormalizeData(second_endo)
second_endo <- ScaleData(second_endo, features = rownames(second_endo))
second_endo <- RunPCA(second_endo, features = rownames(second_endo))

# Get PCs and make elbow plot
npcs <- min(get_pcs(second_endo))
ElbowPlot(second_endo, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 11

# Find neighbors, cluster and UMAP
second_endo <- FindNeighbors(second_endo, dims = 1:npcs)
second_endo <- RunUMAP(second_endo, dims = 1:npcs)
second_endo <- FindClusters(second_endo, resolution = c(0.1, 0.2, 0.3, 0.4))
DimPlot(second_endo, group.by = "RNA_snn_res.0.2", label = TRUE)
DimPlot(second_endo, group.by = "CT_firstpass", label = TRUE)

second_endo <- SetIdent(second_endo, value = "RNA_snn_res.0.2")
endo_markers <- FindAllMarkers(second_endo, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

FeaturePlot(second_endo, features = c("ACKR1", "GNG11", "APLN", "CA4", "HEY1", "FCN3", "CCL21", "KDR", "COL15A1", "PLVAP"))
DotPlot(og_endothelial, features = c("ACKR1", "GNG11", "APLN", "CA4", "HEY1", "FCN3", "CCL21", "KDR", "COL15A1", "PLVAP", "MKI67")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Subcluster 0
second_endo <- FindSubCluster(second_endo, cluster = "0", graph.name = "RNA_snn", resolution = 0.25, subcluster.name = "sub0")
DimPlot(second_endo, group.by = "sub0", label = TRUE)
endo_sub0 <- subset(second_endo, subset = RNA_snn_res.0.2 == "0")
DotPlot(second_endo, features = c("CA4", "FCN3", "KDR", "GNG11", "APLN", "CCL21"), group.by = "sub0")
VlnPlot(second_endo, features = c("CA4", "FCN3", "KDR", "GNG11", "APLN", "CCL21"), group.by = "sub0", pt.size = 0)
VlnPlot(endo_sub0, features = c("CA4", "FCN3", "KDR", "GNG11", "APLN", "CCL21"), group.by = "sub0", pt.size = 0)
FeaturePlot(endo_sub0, features = c("APLN"), split.by = "sub0")
         
# Label cells
second_endo$CT_secondpass <- ""
second_endo$CT_secondpass[second_endo$RNA_snn_res.0.2 == "0"] <- "aCap"
second_endo$CT_secondpass[second_endo$sub0 == "0_4"] <- "gCap"
second_endo$CT_secondpass[second_endo$RNA_snn_res.0.2 == "1"] <- "Venous"
second_endo$CT_secondpass[second_endo$RNA_snn_res.0.2 == "2"] <- "Arteriole"
second_endo$CT_secondpass[second_endo$RNA_snn_res.0.2 == "3"] <- "aCap"
second_endo$CT_secondpass[second_endo$RNA_snn_res.0.2 == "4"] <- "Proliferating Endothelial"
DimPlot(second_endo, group.by = "CT_secondpass", label = TRUE, raster = FALSE)
table(second_endo$CT_firstpass, second_endo$CT_secondpass)
second_endo$CT_secondpass_full <- second_endo$CT_secondpass
#saveRDS(second_endo, "/scratch/avannan/second_endo_secondpass_042423.rds")


#### IMMUNE ----
second_imm <- subset(merged_spatial, 
                     subset = CT_firstpass %in% c("Immune", "B cell", "CD4+ T-cell", "CD8+ T-cell",
                                                  "DC - CCL22+", "Macrophage", "Mast", "Monocyte",
                                                  "Myeloid", "T-cell", "Proliferating T-cell",
                                                  "NK cell", "Plasma/pDC", "CD4+ Treg"))
DefaultAssay(second_imm) <- "RNA"
second_imm <- NormalizeData(second_imm)
second_imm <- ScaleData(second_imm, features = rownames(second_imm))
second_imm <- RunPCA(second_imm, features = rownames(second_imm))

# Get PCs and make elbow plot
npcs <- min(get_pcs(second_imm))
ElbowPlot(second_imm, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 16

# Find neighbors, cluster and UMAP
second_imm <- FindNeighbors(second_imm, dims = 1:npcs)
second_imm <- RunUMAP(second_imm, dims = 1:npcs)
second_imm <- FindClusters(second_imm, resolution = c(0.1, 0.2))
DimPlot(second_imm, group.by = "RNA_snn_res.0.2", label = TRUE)
DimPlot(second_imm, group.by = "CT_firstpass", label = TRUE)
DimPlot(second_imm, group.by = "sample")

second_imm <- SetIdent(second_imm, value = "RNA_snn_res.0.2")
imm_markers <- FindAllMarkers(second_imm, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

DotPlot(second_imm, features = epithelial_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(second_imm, features = immune_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Find pDCs in Plasma (4)
FeaturePlot(second_imm, features = c("LILRA4", "IRF7"))
second_imm <- FindSubCluster(second_imm, cluster = "4", graph.name = "RNA_snn", resolution = 0.08, subcluster.name = "sub4")
DimPlot(second_imm, group.by = "sub6", label = TRUE)
VlnPlot(second_imm, features = "TOP2A", group.by = "sub6", pt.size = 0)
sub6 <- subset(second_imm, subset = RNA_snn_res.0.2 == "6")
FeaturePlot(sub6, features = "TOP2A", split.by = "sub6")

# Subcluster B-cells for proliferation markers (6)
second_imm <- FindSubCluster(second_imm, cluster = "6", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub6")
DimPlot(second_imm, group.by = "sub6", label = TRUE)
sub6 <- subset(second_imm, subset = RNA_snn_res.0.2 == "6")
FeaturePlot(sub6, features = "TOP2A", split.by = "sub6")

# Subcluster Mast cells for proliferation markers (7)
second_imm <- FindSubCluster(second_imm, cluster = "7", graph.name = "RNA_snn", resolution = 0.05, subcluster.name = "sub7")
DimPlot(second_imm, group.by = "sub7", label = TRUE)
sub7 <- subset(second_imm, subset = RNA_snn_res.0.2 == "7")
FeaturePlot(sub7, features = c("TOP2A", "MKI67"), split.by = "sub7")

# Subcluster 5
second_imm <- FindSubCluster(second_imm, cluster = "5", graph.name = "RNA_snn", resolution = 0.08, subcluster.name = "sub5")
DimPlot(second_imm, group.by = "sub5", label = TRUE)
sub5 <- subset(second_imm, subset = RNA_snn_res.0.2 == "5")
sub5 <- SetIdent(sub5, value = "sub5")
sub5_markers <- FindAllMarkers(sub5, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(sub5, features = c("MKI67"), split.by = "sub5")
FeaturePlot(sub5, features = c("CTHRC1"), split.by = "sub5")

# Subcluster 0 (T-cells)
second_imm <- FindSubCluster(second_imm, cluster = "0", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub0")
DimPlot(second_imm, group.by = "sub0", label = TRUE)
sub0 <- subset(second_imm, subset = RNA_snn_res.0.2 == "0")
FeaturePlot(sub0, features = "CCL21", split.by = "sub0")
DotPlot(sub0, features = "CCL21", group.by = "sub0")
sub0 <- SetIdent(sub0, value = "sub0")
sub0_markers <- FindAllMarkers(sub0, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(sub0, features = c("CCL22", "SFTPC"), split.by = "sub0")
second_imm <- SetIdent(second_imm, value = "sub0")
second_imm <- FindSubCluster(second_imm, cluster = "0_4", graph.name = "RNA_snn", resolution = 0.25, subcluster.name = "sub0_4")
sub0_4 <- subset(second_imm, subset = sub0 == "0_4")
FeaturePlot(sub0_4, features = c("CCL21", "SFTPC"), split.by = "sub0_4")
DotPlot(sub0_4, features = c("CCL21", "SFTPC"), group.by = "sub0_4")
DimPlot(second_imm, group.by = "sub0_4", label = TRUE)
VlnPlot(sub0, features = c("CD3E", "CD3D", "CD4", "CD8A", "CD8B", "FOXP3"), group.by = "sub0", pt.size = 0)
DotPlot(second_imm, features = c("CD3E", "CD3D", "CD4", "IL7R", "CD8A", "CD8B", "FOXP3", "LEF1", "CD69"), group.by = "sub0")

# Subcluster 9 (Epithelial cells and other things)
second_imm <- FindSubCluster(second_imm, cluster = "9", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub9")
DimPlot(second_imm, group.by = "sub9", label = TRUE)
sub9 <- subset(second_imm, subset = RNA_snn_res.0.2 == "9")
sub9 <- SetIdent(sub9, value = "sub9")
sub9_markers <- FindAllMarkers(sub9, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(sub9, features = "MS4A7", reduction = "sp")
DimPlot(subset(sub9, subset = sample == "VUILD96MF"), reduction = "sp", group.by = "sub9")
FeaturePlot(sub9, features = c("SCGB3A2", "KRT5"), blend = TRUE, split.by = "sub9")
DotPlot(sub9, features = c(epithelial_features, "MS4A7", "CPA3")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(sub9, features = c("SCGB3A2", "SCGB1A1", "EPCAM", "MMP7",
                           "KRT5", "KRT8", "KRT17", "CPA3", "TPSAB1"), pt.size = 0)
FeaturePlot(sub9, features = c("RTKN2", "SCGB3A2"), blend = TRUE, split.by = "sub9", pt.size = 0.0001)
FeaturePlot(sub9, features = c("RTKN2", "CPA3"), blend = TRUE, split.by = "sub9", pt.size = 0.0001)
second_imm <- SetIdent(second_imm, value = "sub9")
second_imm <- FindSubCluster(second_imm, cluster = "9_1", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub9_1")
DimPlot(second_imm, group.by = "sub9_1", label = TRUE)
sub9_1 <- subset(second_imm, subset = sub9 == "9_1")
sub9_1 <- SetIdent(sub9_1, value = "sub9_1")
DimPlot(sub9_1, group.by = "sub9_1", label = TRUE)
sub9_markers <- FindAllMarkers(sub9_1, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
VlnPlot(sub9_1, features = c("SCGB3A2", "SCGB1A1", "EPCAM", "MMP7",
                             "KRT5", "KRT8", "KRT17", "CPA3", "TPSAB1"), pt.size = 0)
DimPlot(subset(sub9_1, subset = sample == "VUILD78MF"), reduction = "sp")

# View fibroblast markers in mesenchymal cluster
sub5 <- subset(second_imm, subset = RNA_snn_res.0.2 == "5")
DimPlot(sub5, group.by = "RNA_snn_res.0.2")
FeaturePlot(sub5, features = mesenchymal_features)

# Subcluster 8
second_imm <- FindSubCluster(second_imm, cluster = "8", graph.name = "RNA_snn", resolution = 0.15, subcluster.name = "sub8")
DimPlot(second_imm, group.by = "sub8", label = TRUE)
FeaturePlot(second_imm, features = c("KDR", "HLA-DRA", "GNG11"))
sub8 <- subset(second_imm, subset = RNA_snn_res.0.2 == "8")
DimPlot(sub8, reduction = "sp")
DimPlot(subset(sub8, subset = sample == "TILD175"), reduction = "sp", group.by = "sub8")
FeaturePlot(subset(sub8, subset = sample == "TILD175"), reduction = "sp", features = c("GNG11", "HLA-DRA"))
DimPlot(second_endo, reduction = "sp", group.by = "CT_secondpass", raster = FALSE, cols = c("red", "blue", "orange", "grey", "grey", "green"))
DotPlot(sub9, features = immune_features, group.by = "sub9") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Macrophages
second_imm <- SetIdent(second_imm, value = "RNA_snn_res.0.2")
FeaturePlot(second_imm, features = c("MS4A7", "SPP1", "FABP4", "CD14", "FCGR3A", "EREG", "PPARG", "LYZ"))
DotPlot(second_imm, features = c("MS4A7", "SPP1", "FABP4", "CD14", "FCGR3A", "EREG", "PPARG", "LYZ"))
macro_sub <- subset(second_imm, subset = RNA_snn_res.0.2 %in% c(1, 2, 3))
DimPlot(macro_sub, reduction = "sp", raster = FALSE)
DimPlot(subset(macro_sub, subset = RNA_snn_res.0.2 != "3"), reduction = "sp", raster = FALSE)

# Not sure if plasma cells are actually PPARG+ alveolar macrophages
FeaturePlot(second_imm, features = c("PPARG", "LILRA4", "XBP1", "HERPUD1"))

# Label cells
second_imm$CT_secondpassA <- ""
second_imm$CT_secondpassA[second_imm$sub0 == "0_0"] <- "CD8+ T-cells"
second_imm$CT_secondpassA[second_imm$sub0 == "0_1"] <- "Tregs"
second_imm$CT_secondpassA[second_imm$sub0 == "0_2"] <- "T-cells"
second_imm$CT_secondpassA[second_imm$sub0 == "0_3"] <- "T-cells (Interferon High)"
second_imm$CT_secondpassA[second_imm$sub0 == "0_4"] <- "DCs - CCL22+"
second_imm$CT_secondpassA[second_imm$RNA_snn_res.0.2 == "1"] <- "Macrophages"
second_imm$CT_secondpassA[second_imm$RNA_snn_res.0.2 == "2"] <- "Macrophages"
second_imm$CT_secondpassA[second_imm$RNA_snn_res.0.2 == "3"] <- "FABP4+/SPP1+ Macrophages"
second_imm$CT_secondpassA[second_imm$sub4 == "4_0"] <- "Plasma"
second_imm$CT_secondpassA[second_imm$sub4 == "4_1"] <- "Proliferating Plasma"
second_imm$CT_secondpassA[second_imm$sub4 == "4_2"] <- "pDCs"
second_imm$CT_secondpassA[second_imm$sub5 == "5_0"] <- "SMCs"
second_imm$CT_secondpassA[second_imm$sub5 == "5_1"] <- "Fibroblasts"
second_imm$CT_secondpassA[second_imm$sub5 == "5_2"] <- "Proliferating Mesenchymal"
second_imm$CT_secondpassA[second_imm$sub6 == "6_1"] <- "Proliferating B cells"
second_imm$CT_secondpassA[second_imm$sub6 == "6_0"] <- "B cells"
second_imm$CT_secondpassA[second_imm$sub7 %in% c("7_0", "7_1")] <- "Mast"
second_imm$CT_secondpassA[second_imm$sub7 == "7_2"] <- "Proliferating Mast"
second_imm$CT_secondpassA[second_imm$RNA_snn_res.0.2 == "8"] <- "DC/Lymphatic"
second_imm$CT_secondpassA[second_imm$sub9 == "9_0"] <- "Macrophages"
second_imm$CT_secondpassA[second_imm$sub9 == "9_1"] <- "Mast"
second_imm$CT_secondpassA[second_imm$sub9 == "9_2"] <- "Transitional AT2"
second_imm$CT_secondpassA[second_imm$sub9 == "9_3"] <- "CD8+ T-cells"
second_imm$CT_secondpassA[second_imm$sub9 == "9_4"] <- "cDCs"
second_imm$CT_secondpassA[second_imm$RNA_snn_res.0.2 == "10"] <- "NK cells"
second_imm$CT_secondpassA[second_imm$RNA_snn_res.0.2 == "11"] <- "Proliferating T-cells"
DimPlot(second_imm, group.by = "CT_secondpassA", label = TRUE, raster = FALSE)
DotPlot(second_imm, features = immune_features, group.by = "CT_secondpassA") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(second_imm$CT_firstpass, second_imm$CT_secondpassA)
# saveRDS(second_imm, "/scratch/avannan/xenium_files/second_imm_secondpass_04242023.rds")


###### LYMPHOID ----
second_lym <- subset(second_imm, 
                     subset = CT_secondpassA %in% c("T-cells", "Proliferating T-cells", "CD8+ T-cells", "NK cells",
                                                    "T-cells (Interferon High)", "Plasma", "pDCs", "Proliferating Plasma",
                                                    "B cells", "Proliferating B cells", "Tregs"))
                     
DefaultAssay(second_lym) <- "RNA"
second_lym <- NormalizeData(second_lym)
second_lym <- ScaleData(second_lym, features = rownames(second_lym))
second_lym <- RunPCA(second_lym, features = rownames(second_lym))

# Get PCs and make elbow plot
npcs <- min(get_pcs(second_lym))
ElbowPlot(second_lym, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 12

# Find neighbors, cluster and UMAP
second_lym <- FindNeighbors(second_lym, dims = 1:npcs)
second_lym <- RunUMAP(second_lym, dims = 1:npcs)
second_lym <- FindClusters(second_lym, resolution = c(0.1, 0.2))
DimPlot(second_lym, group.by = "RNA_snn_res.0.2", label = TRUE)
DimPlot(second_lym, group.by = "CT_firstpass", label = TRUE)
DimPlot(second_lym, group.by = "CT_secondpassA", label = TRUE)
DimPlot(second_lym, group.by = "sample")

second_lym <- SetIdent(second_lym, value = "RNA_snn_res.0.2")
lym_markers <- FindAllMarkers(second_lym, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(second_lym, features = immune_features[1:20])
FeaturePlot(second_lym, features = immune_features[21:40])
FeaturePlot(second_lym, features = immune_features[41:60])
FeaturePlot(second_lym, features = immune_features[61:84])
DotPlot(second_lym, features = immune_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Proliferating cells
FeaturePlot(second_lym, features = c("TOP2A", "MKI67"))

# Subset 3 (B cells)
second_lym <- FindSubCluster(second_lym, cluster = "3", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub3")
DimPlot(second_lym, group.by = "sub3", label = TRUE)

# Subset 5 (Proliferating cluster)
second_lym <- FindSubCluster(second_lym, cluster = "5", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub5")
DimPlot(second_lym, group.by = "sub5", label = TRUE)
DotPlot(second_lym, features = c("MKI67", "TOP2A", "SEC11C", "CD3E", "MS4A7"), group.by = "sub5")

# Subset 4 (NK cells)
second_lym <- FindSubCluster(second_lym, cluster = "4", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub4")
DimPlot(second_lym, group.by = "sub4", label = TRUE)
DotPlot(second_lym, features = c("MKI67", "TOP2A"), group.by = "sub4")
sub4 <- subset(second_lym, subset = RNA_snn_res.0.2 == "4")
FeaturePlot(sub4, features = c("GNLY", "TOP2A"), split.by = "sub4")

# Subset 0 (majority of T-cells)
second_lym <- FindSubCluster(second_lym, cluster = "0", graph.name = "RNA_snn", resolution = 0.15, subcluster.name = "sub0")
DimPlot(second_lym, group.by = "sub0", label = TRUE)
second_lym <- SetIdent(second_lym, value = "RNA_snn_res.0.2")
sub0 <- subset(second_lym, subset = RNA_snn_res.0.2 == "0")
sub0 <- SetIdent(sub0, value = "sub0")
lym0_markers <- FindAllMarkers(sub0, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(sub0, features = immune_features, group.by = "sub0") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(sub0, features = endothelial_features, group.by = "sub0") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DimPlot(sub0, group.by = "sub0", reduction = "sp", cols = c("red", "blue", "green", "orange"))
DimPlot(subset(sub0, subset = sample == "VUILD96MF"), group.by = "sub0", 
        reduction = "sp", cols = c("red", "blue", "green", "orange"))
DimPlot(subset(sub0, subset = sample == "VUILD115"), group.by = "sub0", 
        reduction = "sp", cols = c("red", "blue", "green", "orange"))
FeaturePlot(sub0, features = "HLA-DQA1", split.by = "sub0")
FeaturePlot(subset(sub0, subset = sample == "VUILD96MF"), features = c("FABP4", "SPP1", "CD4", "CD8A"), reduction = "sp", split.by = "sub0")
FeaturePlot(subset(sub0, subset = sample == "VUILD110"), features = c("FABP4", "SPP1", "CD4", "CD8A"), reduction = "sp", split.by = "sub0")
FeaturePlot(subset(sub0, subset = sample == "VUILD115"), features = c("CD3E", "CD4", "CD8A"), reduction = "sp", split.by = "sub0")
FeaturePlot(sub0, features = c("CD4", "CD8A", "IL7R"))
second_lym <- SetIdent(second_lym, value = "sub0")
second_lym <- FindSubCluster(second_lym, cluster = "0_1", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub0_1")
DimPlot(second_lym, group.by = "sub0_1", label = TRUE)
sub0 <- subset(second_lym, subset = RNA_snn_res.0.2 == "0")
DotPlot(sub0, features = c("CD4", "CD8A", "IL7R"), group.by = "sub0_1")
VlnPlot(sub0, features = c("CD4", "CD8A", "IL7R"), group.by = "sub0_1", pt.size = 0)
VlnPlot(sub0, pt.size = 0, group.by = "sub0", features = c("CD4", "IL7R", "SPP1", "MS4A7", "HLA-DQA1", "CD3E"), ncol = 2)

# Subset 6 (Epithelial contamination)
second_lym <- FindSubCluster(second_lym, cluster = "6", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub6")
DimPlot(second_lym, group.by = "sub6", label = TRUE)
sub6 <- subset(second_lym, subset = RNA_snn_res.0.2 == "6")
DimPlot(sub6, group.by = "sub6", reduction = "sp")
DotPlot(second_lym, features = epithelial_features, group.by = "sub6") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(second_lym, features = immune_features, group.by = "sub6") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(second_lym, features = c("CD3E", "CD4", "CD8A", "NKG7", "SCGB3A2", "SCGB1A1"), group.by = "sub6", pt.size = 0, ncol = 2)
DimPlot(subset(sub6, subset = sample == "VUILD96MF"), reduction = "sp", group.by = "sub6")
DimPlot(subset(sub6, subset = sample == "VUILD110"), reduction = "sp", group.by = "sub6")

# Label cells
second_lym$CT_secondpassB <- ""
second_lym$CT_secondpassB[second_lym$sub0 == "0_0"] <- "CD4+ T-cells"
second_lym$CT_secondpassB[second_lym$sub0 == "0_1"] <- "CD8+ T-cells"
second_lym$CT_secondpassB[second_lym$sub0 == "0_2"] <- "CD4+ Tregs"
second_lym$CT_secondpassB[second_lym$sub0 == "0_3"] <- "DCs"
second_lym$CT_secondpassB[second_lym$RNA_snn_res.0.2 == "1"] <- "Plasma"
second_lym$CT_secondpassB[second_lym$RNA_snn_res.0.2 == "2"] <- "CD8+ T-cells"
second_lym$CT_secondpassB[second_lym$sub3 == "3_0"] <- "B cells"
second_lym$CT_secondpassB[second_lym$sub3 == "3_1"] <- "Proliferating B cells"
second_lym$CT_secondpassB[second_lym$sub4 == "4_0"] <- "NK cells"
second_lym$CT_secondpassB[second_lym$sub4 == "4_1"] <- "Proliferating NK cells"
second_lym$CT_secondpassB[second_lym$sub5 == "5_0"] <- "Proliferating T-cells"
second_lym$CT_secondpassB[second_lym$sub5 == "5_1"] <- "Proliferating Macrophages"
second_lym$CT_secondpassB[second_lym$sub5 == "5_2"] <- "Proliferating Plasma"
second_lym$CT_secondpassB[second_lym$sub6 == "6_0"] <- "CD8+ T-cells"
second_lym$CT_secondpassB[second_lym$sub6 == "6_1"] <- "pDCs"
DimPlot(second_lym, group.by = "CT_secondpassA", label = TRUE, raster = FALSE)
DimPlot(second_lym, group.by = "CT_secondpassB", label = TRUE, raster = FALSE, repel = TRUE)
DotPlot(second_lym, features = c(immune_features, "MKI67", "TOP2A"), group.by = "CT_secondpassB") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(second_lym$CT_secondpassA, second_lym$CT_secondpassB)
# saveRDS(second_lym, "/scratch/avannan/second_lym_secondpass042423.rds")



###### MYELOID ----
second_my <- subset(second_imm, 
                    subset = CT_secondpassA %in% c("Macrophages", "FABP4+/SPP1+ Macrophages",
                                                   "Mast", "Proliferating Mast",
                                                   "cDCs", "DC/Lymphatic", "DCs - CCL22+"))

DefaultAssay(second_my) <- "RNA"
second_my <- NormalizeData(second_my)
second_my <- ScaleData(second_my, features = rownames(second_my))
second_my <- RunPCA(second_my, features = rownames(second_my))

# Get PCs and make elbow plot
npcs <- min(get_pcs(second_my))
ElbowPlot(second_my, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 17

# Find neighbors, cluster and UMAP
second_my <- FindNeighbors(second_my, dims = 1:npcs)
second_my <- RunUMAP(second_my, dims = 1:npcs)
second_my <- FindClusters(second_my, resolution = 0.2)
DimPlot(second_my, group.by = "RNA_snn_res.0.2", label = TRUE)
DimPlot(second_my, group.by = "CT_firstpass", label = TRUE)
DimPlot(second_my, group.by = "CT_secondpassA", label = TRUE)
DimPlot(second_my, group.by = "sample")

second_my <- SetIdent(second_my, value = "RNA_snn_res.0.2")
my_markers <- FindAllMarkers(second_my, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(second_my, features = c(immune_features, "MKI67", "TOP2A")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Look at general markers
FeaturePlot(second_my, features = immune_features[1:20])
FeaturePlot(second_my, features = immune_features[21:40])
FeaturePlot(second_my, features = immune_features[41:60])
FeaturePlot(second_my, features = immune_features[61:84])

# Proliferating
FeaturePlot(second_my, features = c("TOP2A", "MKI67"))

# FABP4 and SPP1 - both in 1
FeaturePlot(second_my, features = c("FABP4", "SPP1"))
second_my <- FindSubCluster(second_my, cluster = "1", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub1")
DimPlot(second_my, group.by = "sub1", label = TRUE)
VlnPlot(second_my, features = c("FABP4", "SPP1"), pt.size = 0, group.by = "sub1")

# Mast cells
FeaturePlot(second_my, features = c("CPA3", "KIT"))
second_my <- FindSubCluster(second_my, cluster = "3", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub3")
DimPlot(second_my, group.by = "sub3", label = TRUE)
sub3 <- subset(second_my, subset = RNA_snn_res.0.2 == "3")
FeaturePlot(sub3, features = c("CPA3", "MKI67"), split.by = "sub3")

# CCL22+ DCs
sub10 <- subset(second_my, subset = RNA_snn_res.0.2 == "10")
DimPlot(sub10, reduction = "sp")

# Other DCs
FeaturePlot(second_my, features = c("FCER1A", "CD1C", "CD11"))

# Subcluster proliferating cells from 7
second_my <- FindSubCluster(second_my, cluster = "7", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub7")
DimPlot(second_my, group.by = "sub7", label = TRUE)
sub7 <- subset(second_my, subset = RNA_snn_res.0.2 == "7")
DimPlot(sub7, group.by = "sub7", reduction = "sp")
DimPlot(sub7, group.by = "sub7", label = TRUE)
FeaturePlot(sub7, features = c("ACKR1", "HEY1", "APLN", "PLVAP", "COL15A1", "CA4", "MKI67", "TOP2A"), pt.size = 0.0001)
VlnPlot(sub7, features = c("ACKR1", "HEY1", "APLN", "PLVAP", "COL15A1", "CA4", "MKI67", "TOP2A"), group.by = "sub7", pt.size = 0, ncol = 4)

# Look at two populations of monocye-derived macrophages
mono_sub <- subset(second_my, subset = RNA_snn_res.0.2 %in% c("2", "4"))
DimPlot(mono_sub, reduction = "sp")
FeaturePlot(mono_sub, reduction = "sp", features = c("SFTPC", "SCGB1A1", "RTKN2"))

# Epithelial/endothelial contamination
FeaturePlot(second_my, features = c("SFTPC", "SCGB1A1", "RTKN2"))
FeaturePlot(second_my, features = c("ACKR1", "PLVAP", "HEY1", "COL15A1"))

# Subcluster T-cells (9)
second_my <- FindSubCluster(second_my, cluster = "9", graph.name = "RNA_snn", resolution = 0.3, subcluster.name = "sub9")
DimPlot(second_my, group.by = "sub9", label = TRUE)
FeaturePlot(second_my, features = c("CD3E", "CD4", "FOXP3", "CD8A"))
DotPlot(second_my, features = c("CD3E", "CD4", "FOXP3", "CD8A", "IL7R", "GNLY"), group.by = "sub9")

# Monocyte markers
FeaturePlot(second_my, features = c("EREG", "STAT6", "FCN1", "S100A12", "S100A8",
                                    "S100A9", "AIF1", "CD14", "FCGR3A"))
VlnPlot(second_my, features = c("EREG", "STAT6", "FCN1", "S100A12", "S100A8",
                                "S100A9", "AIF1", "CD14", "FCGR3A")), pt.size = 0)
second_my <- FindSubCluster(second_my, cluster = "2", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub2")
DimPlot(second_my, group.by = "sub2", label = TRUE)
DimPlot(subset(macro_sub, subset = RNA_snn_res.0.2 %in% c("2", "4")), group.by = "sub2", reduction = "sp", raster = FALSE)
# Split macro
macro_sub <- subset(second_my, subset = sub1 %in% c("1_0", "1_1", "0", "2", "4"))
DimPlot(macro_sub, group.by = "sub1", reduction = "sp", raster = FALSE)
macro_sub <- SetIdent(macro_sub, value = "sub1")
macro_markers <- FindAllMarkers(macro_sub, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DimPlot(subset(macro_sub, subset = RNA_snn_res.0.2 %in% c("2", "4")), group.by = "sub1", reduction = "sp", raster = FALSE)

# Split cluster 0
second_my <- FindSubCluster(second_my, cluster = "0", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub0")
DimPlot(second_my, group.by = "sub0", label = TRUE)
FeaturePlot(second_my, features = c("OAS2", "OAS3", "IFIT1", "IFIT3"))
sub0 <- subset(second_my, subset = RNA_snn_res.0.2 == "0")
DimPlot(sub0, reduction = "sp", group.by = "sub0")
FeaturePlot(sub0, reduction = "sp", features = "IFIT1", split.by = "sub0")

# Split cluster 5
second_my <- FindSubCluster(second_my, cluster = "5", graph.name = "RNA_snn", resolution = 0.15, subcluster.name = "sub5")
DimPlot(second_my, group.by = "sub5", label = TRUE)
DotPlot(second_my, features = c("OAS2", "OAS3", "IFIT1", "IFIT3"), group.by = "sub5")

# Look further at 6
FeaturePlot(second_my, features = c("CPA3", "TPSAB1", "SCGB3A2", "SCGB1A1")) # These are mast cells with epithelial contamination
second_my <- FindSubCluster(second_my, cluster = "6", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub6")
DimPlot(second_my, group.by = "sub6", label = TRUE)
sub6 <- subset(second_my, subset = RNA_snn_res.0.2 == "6")
sub6 <- SetIdent(sub6, value = "sub6")
sub6_markers <- FindMarkers(sub6, ident.1 = "6_0", ident.2 = "6_1") %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DimPlot(subset(sub6, subset = sample == "VUILD96MF"), group.by = "sub6", reduction = "sp")
FeaturePlot(second_my, features = c("CST3", "FCER1A", "CD1A", "CD1C", "HLA-DQB1", "CXCL9", "GPR183", "HLA-DQA1"))

# Label cells
second_my$CT_secondpassB <- ""
second_my$CT_secondpassB[second_my$sub0 == "0_0"] <- "Macrophages"
second_my$CT_secondpassB[second_my$sub0 == "0_1"] <- "Macrophages (Interferon High)"
second_my$CT_secondpassB[second_my$sub1 == "1_0"] <- "SPP1+ Macrophages"
second_my$CT_secondpassB[second_my$sub1 == "1_1"] <- "FABP4+ Macrophages"
second_my$CT_secondpassB[second_my$sub2 == "2_0"] <- "Interstitial Macrophages"
second_my$CT_secondpassB[second_my$sub2 == "2_1"] <- "Interstitial Macrophages (Interferon High)"
second_my$CT_secondpassB[second_my$RNA_snn_res.0.2 == "3"] <- "Mast"
second_my$CT_secondpassB[second_my$RNA_snn_res.0.2 == "4"] <- "Interstitial Macrophages (FCN1+)"
second_my$CT_secondpassB[second_my$RNA_snn_res.0.2 == "5"] <- "DCs/Lymphatic"
second_my$CT_secondpassB[second_my$sub6 == "6_0"] <- "cDCs"
second_my$CT_secondpassB[second_my$sub6 == "6_1"] <- "Mast"
second_my$CT_secondpassB[second_my$sub7 == "7_0"] <- "Venous"
second_my$CT_secondpassB[second_my$sub7 == "7_1"] <- "Proliferating Endothelial"
second_my$CT_secondpassB[second_my$RNA_snn_res.0.2 == "8"] <- "cDCs"
second_my$CT_secondpassB[second_my$RNA_snn_res.0.2 == "9"] <- "T-cells"
second_my$CT_secondpassB[second_my$RNA_snn_res.0.2 == "10"] <- "DCs - CCL22+"
DimPlot(second_my, group.by = "CT_secondpassA", label = TRUE, raster = FALSE)
DimPlot(second_my, group.by = "CT_secondpassB", label = TRUE, raster = FALSE, repel = TRUE)
DotPlot(second_my, features = c(immune_features[-c(2:6)], "CD1C", "CST3", "FCER1A", "IFIT1", "IFIT3",
                                "SCGB3A2", "SCGB1A1", "ACKR1", "PLVAP", "EPAS1", "MKI67", "TOP2A"), 
        group.by = "CT_secondpassB") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
table(second_my$CT_secondpassA, second_my$CT_secondpassB)
#saveRDS(second_my, "/scratch/avannan/xenium_files/second_my_secondpass_042523.rds")

# Add annotations to immune object
full_imm_labels <- merge(second_my, y = second_lym)
second_imm$CT_secondpassB <- second_imm$CT_secondpassA
second_imm$CT_secondpassB <- full_imm_labels$CT_secondpassB
second_imm$CT_secondpassB[is.na(second_imm$CT_secondpassB)] <- second_imm$CT_secondpassB
second_imm$CT_secondpassB <- full_imm_labels$CT_secondpassB
DimPlot(second_imm, group.by = "CT_secondpassA", label = TRUE)
DimPlot(second_imm, group.by = "CT_secondpassB", label = TRUE)
test <- subset(second_imm, subset = CT_secondpassA %in% c("Fibroblasts", "Proliferating Mesenchymal", "Transitional AT2", "SMCs"))
DimPlot(test, group.by = "CT_secondpassA")

# Add annotations to immune object
full_imm_labels <- merge(second_my, y = second_lym)
second_imm$CT_secondpass_full <- full_imm_labels$CT_secondpassB
second_imm$CT_secondpass_full[second_imm$CT_secondpassA == "Fibroblasts"] <- "Fibroblasts"
second_imm$CT_secondpass_full[second_imm$CT_secondpassA == "SMCs"] <- "SMCs"
second_imm$CT_secondpass_full[second_imm$CT_secondpassA == "Proliferating Mesenchymal"] <- "Proliferating Mesenchymal"
second_imm$CT_secondpass_full[second_imm$CT_secondpassA == "Transitional AT2"] <- "Transitional AT2"
DimPlot(second_imm, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE)


#### MESENCHYMAL ----
second_mes <- subset(merged_spatial, subset = CT_firstpass %in% c("Mesenchymal", "Adventitial FB"))

DefaultAssay(second_mes) <- "RNA"
second_mes <- NormalizeData(second_mes)
second_mes <- ScaleData(second_mes, features = rownames(second_mes))
second_mes <- RunPCA(second_mes, features = rownames(second_mes))

# Get PCs and make elbow plot
npcs <- min(get_pcs(second_mes))
ElbowPlot(second_mes, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 16

# Find neighbors, cluster and UMAP
second_mes <- FindNeighbors(second_mes, dims = 1:npcs)
second_mes <- RunUMAP(second_mes, dims = 1:npcs)
second_mes <- FindClusters(second_mes, resolution = c(0.2, 0.3))
DimPlot(second_mes, group.by = "RNA_snn_res.0.3", label = TRUE)
DimPlot(second_mes, group.by = "CT_firstpass", label = TRUE)

second_mes <- SetIdent(second_mes, value = "RNA_snn_res.0.3")
mes_markers <- FindAllMarkers(second_mes, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(second_mes, features = c("FAP", "CTHRC1"))

# Look at fibroblast and other markers
DotPlot(second_mes, features = c("COL1A2", "SPARCL1", "EPAS1", "ELN", "ITGAX", "TGFB3", "LGR6",
                                 "POSTN", "CTHRC1", "COL1A1", "LUM", "COL3A1", "SFRP2",
                                 "PTGDS", "DCN", "MEG3", "TGFB2", "ITGA3", "HAS2",
                                 "ACTA2", "MS4A7", "SPP1", "FABP4", "KDR", "GNG11", 
                                 "PDGFRA", "HLA-DQB1", "CD14", 
                                 "FCER1G", "CCL21", "PDGFRB", "IFIT1", "IFIT3", 
                                 "OAS2", "OAS3", "STAT1", "GDF15", 
                                 "LGR5", "WNT5A", "SCGB3A2", "SCGB1A1",
                                 "CDK1", "TOP2A", "MKI67", "CENPF",
                                 "MSLN", "ELANE", "HAS1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(second_mes, features = c(immune_features, epithelial_features)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = ))
FeaturePlot(second_mes, features = mesenchymal_features)

# Airway markers around CTHRC1+/FAP+ FBs?
second_mes <- FindSubCluster(second_mes, cluster = "3", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub3")
DimPlot(second_mes, group.by = "sub3", label = TRUE)
DotPlot(second_mes, features = c("SCGB3A2", "SCGB1A1", "CTHRC1", "FAP", "LGR5"), group.by = "sub3")
FeaturePlot(second_mes, features = c("SCGB3A2", "SCGB1A1", "CTHRC1", "FAP", "LGR5"))
sub3_1 <- subset(second_mes, subset = sub3 == "3_1")
DimPlot(sub3_1, reduction = "sp")
DimPlot(subset(sub3_1, subset = sample == "VUILD96LF"), reduction = "sp")

# Subcluster 8 - proliferating
second_mes <- FindSubCluster(second_mes, cluster = "8", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub8")
DimPlot(second_mes, group.by = "sub8", label = TRUE)
second_mes <- SetIdent(second_mes, value = "sub8")

# Adventitial
FeaturePlot(second_mes, features = c("MFAP5", "PI16", "HAS1"))
DotPlot(second_mes, features = c("MFAP5", "PI16", "HAS1"))
VlnPlot(second_mes, features = c("MFAP5", "PI16"), pt.size = 0)

# Interferon response
FeaturePlot(second_mes, features = c("IFIT1", "IFIT3", "OAS2", "OAS3"))

# Cluster 6 markers
FeaturePlot(second_mes, features = c("FABP4", "HLA-DRA", "MUC5AC", "SPRY2"))

# Subcluster 0
second_mes <- FindSubCluster(second_mes, cluster = "0", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub0")
DimPlot(second_mes, group.by = "sub0", label = TRUE)
second_mes <- SetIdent(second_mes, value = "sub0")
second_mes <- FindSubCluster(second_mes, cluster = "0_2", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "sub0_2")
DimPlot(second_mes, group.by = "sub0_2", label = TRUE)
DimPlot(sub0, group.by = "sub0_2", label = TRUE, repel = TRUE)
sub0 <- subset(second_mes, subset = RNA_snn_res.0.3 == "0")
sub0_2 <- subset(second_mes, subset = sub0 == "0_2")
DimPlot(sub0, reduction = "sp", group.by = "sub0_2")
DimPlot(sub0_2, reduction = "sp", group.by = "sub0_2")
sub0 <- SetIdent(sub0, value = "sub0_2")
mes0_markers <- FindAllMarkers(sub0, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(second_mes, features = c("MFAP5", "PI16", "HAS1", "PLIN2"), group.by = "sub0_2")
DotPlot(subset(second_mes, subset = RNA_snn_res.0.3 != "9"), features = c("MFAP5", "PI16", "HAS1", "PLIN2"), group.by = "sub0_2")
DotPlot(sub0, features = c("MFAP5", "PI16", "HAS1", "PLIN2"), group.by = "sub0_2")
FeaturePlot(sub0, features = "HAS1", split.by = "sub0_2", pt.size = 0.001)

# Subcluster 1 - SMCs
second_mes <- FindSubCluster(second_mes, cluster = "1", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub1")
DimPlot(second_mes, group.by = "sub1", label = TRUE)
sub1 <- subset(second_mes, subset = RNA_snn_res.0.3 == "1")
DimPlot(sub1, reduction = "sp", group.by = "sub1")
FeaturePlot(second_mes, features = c("AXL", "ACTA2", "LGR6", "COL1A2"))

# Proliferating
FeaturePlot(second_mes, features = c("MKI67", "TOP2A", "FN1", "MS4A7"))

# Lipofibroblast markers (6)
sub6 <- subset(second_mes, subset = RNA_snn_res.0.3 == "6")
FeaturePlot(sub6, reduction = "sp", features = "FABP4")
DotPlot(second_mes, features = c("PPARG", "MRC1", "CD44", "FABP4"), group.by = "RNA_snn_res.0.3")
FeaturePlot(subset(sub6, subset = sample == "VUILD104MF"), reduction = "sp", features = "FABP4")
FeaturePlot(subset(sub6, subset = sample == "VUILD106"), reduction = "sp", features = "FABP4")
second_mes <- FindSubCluster(second_mes, cluster = "6", graph.name = "RNA_snn", resolution = 0.3, subcluster.name = "sub6")
sub6 <- subset(second_mes, subset = RNA_snn_res.0.3 == "6")
DimPlot(subset(sub6, subset = sample == "VUILD106"), reduction = "sp", group.by = "sub6")
DimPlot(subset(sub6, subset = sub6 == "6_0"), reduction = "sp", group.by = "sub6")
DimPlot(second_mes, group.by = "sub6", label = TRUE)
# Can't subset

# Label cell types
second_mes$CT_secondpass <- ""
second_mes$CT_secondpass[second_mes$RNA_snn_res.0.3 == "0"] <- "Fibroblasts"
second_mes$CT_secondpass[second_mes$sub0_2 == "0_1"] <- "PLIN2+ Fibroblasts"
second_mes$CT_secondpass[second_mes$sub0_2 == "0_2_2"] <- "HAS1+ Fibroblasts"
second_mes$CT_secondpass[second_mes$sub1 %in% c("1_0", "1_1")] <- "SMCs"
second_mes$CT_secondpass[second_mes$sub1 == "1_2"] <- "Fibroblasts"
second_mes$CT_secondpass[second_mes$RNA_snn_res.0.3 == "2"] <- "Fibroblasts"
second_mes$CT_secondpass[second_mes$sub3 == "3_0"] <- "Activated FBs (CTHRC1+/FAP+)"
second_mes$CT_secondpass[second_mes$sub3 == "3_1"] <- "Peribronchial FBs"
second_mes$CT_secondpass[second_mes$RNA_snn_res.0.3 == "4"] <- "Fibroblasts (Interferon High)"
second_mes$CT_secondpass[second_mes$RNA_snn_res.0.3 == "5"] <- "T-cells"
second_mes$CT_secondpass[second_mes$RNA_snn_res.0.3 == "6"] <- "Lipofibroblasts"
second_mes$CT_secondpass[second_mes$RNA_snn_res.0.3 == "7"] <- "Adventitial FBs"
second_mes$CT_secondpass[second_mes$sub8 %in% c("8_0", "8_2")] <- "Proliferating Mesenchymal"
second_mes$CT_secondpass[second_mes$sub8 == "8_1"] <- "Proliferating Macrophages"
second_mes$CT_secondpass[second_mes$RNA_snn_res.0.3 == "9"] <- "Mesothelial"
second_mes$CT_secondpass_full <- second_mes$CT_secondpass

DimPlot(second_mes, group.by = "CT_secondpass", label = TRUE, repel = TRUE)
DotPlot(second_mes, features = c("MSLN", "HAS1", "PLIN2", "COL1A1", 
                                 "COL1A2", "COL3A1", "FN1", "VIM",
                                 "MFAP5", "PI16", "SFRP2", "LGR6", "AXL",
                                 "MKI67", "TOP2A", "CD3E", "CD3D", 
                                 "CD4", "CD8A", "CD8B",  "NKG7", "GNLY",
                                 "MS4A7", "SCGB3A2", "SCGB1A1", "LGR5",
                                 "IFIT1", "IFIT3", "OAS2", "OAS3", "CTHRC1", "FAP",
                                 "MRC1", "FABP4", "CD44", "PPARG"),
        group.by = "CT_secondpass") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(second_mes, features = immune_features, group.by = "CT_secondpass") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(second_mes$CT_secondpass, second_mes$CT_firstpass)


## FINISH ANNOTATIONS ----
# saveRDS(merged_spatial, "/scratch/avannan/xenium_files/xenium_firstpass_04202023.rds")
# saveRDS(epithelial, "/scratch/avannan/xenium_files/epithelial_firstpass_04202023.rds")
# saveRDS(airway, "/scratch/avannan/xenium_files/airway_firstpass_04202023.rds")
# saveRDS(alv, "/scratch/avannan/xenium_files/alv_firstpass_04202023.rds")
# saveRDS(og_endothelial, "/scratch/avannan/xenium_files/og_endothelial_firstpass_04202023.rds")
# saveRDS(endothelial, "/scratch/avannan/xenium_files/endothelial_firstpass_04202023.rds")
# saveRDS(immune, "/scratch/avannan/xenium_files/immune_firstpass_04202023.rds")
# saveRDS(lymphoid, "/scratch/avannan/xenium_files/lymphoid_firstpass_04202023.rds")
# saveRDS(myeloid, "/scratch/avannan/xenium_files/myeloid_firstpass_04202023.rds")
# saveRDS(og_mesenchymal, "/scratch/avannan/xenium_files/og_mesenchymal_firstpass_04202023.rds")

#saveRDS(merged_spatial, "/scratch/avannan/xenium_files/xenium_secondpass_04262023.rds")

test_epi <- readRDS("/scratch/avannan/xenium_files/epithelial_firstpass_04202023.rds")

merged_spatial <- readRDS("/scratch/avannan/xenium_files/xenium_secondpass_04262023.rds")
second_epi <- readRDS("/scratch/avannan/xenium_files/second_epi_secondpass_04212023.rds")
second_air <- readRDS("/scratch/avannan/xenium_files/second_airway_secondpass_0421023.rds")
second_alv <- readRDS("/scratch/avannan/xenium_files/second_alv_secondpass_04252023.rds")
second_endo <- readRDS("/scratch/avannan/second_endo_secondpass_042423.rds")
second_imm <- readRDS("/scratch/avannan/xenium_files/second_imm_secondpass_04242023.rds")
second_lym <- readRDS("/scratch/avannan/second_lym_secondpass042423.rds")
second_my <- readRDS("/scratch/avannan/xenium_files/second_my_secondpass_042523.rds")
second_mes <- readRDS("/scratch/avannan/xenium_files/second_mes_secondpass_04212023.rds")

second_endo$CT_secondpass_full <- second_endo$CT_secondpass
second_mes$CT_secondpass_full <- second_mes$CT_secondpass

# Fix T-cell labels on epithelial, mesenchymal, and immune
t_mes <- subset(second_mes, subset = CT_secondpass_full == "T-cells")
t_mes$CT_secondpass_full <- "T-cells (Mes)"
t_epi <- subset(second_epi, subset = CT_secondpass_full == "T-cells")
t_epi$CT_secondpass_full <- "T-cells (Epi)"
tcell_test <- merge(subset(second_imm, 
                           subset = CT_secondpass_full %in% c("T-cells", "CD4+ Tregs", 
                                                          "CD4+ T-cells", "CD8+ T-cells")),
                    y = list(t_mes, t_epi))
DotPlot(tcell_test, features = immune_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(tcell_test, features = c(mesenchymal_features, "CD3E"), group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DimPlot(t_epi, reduction = "sp", group.by = "CT_secondpass_full")
DimPlot(t_mes, reduction = "sp", group.by = "CT_secondpass_full")
FeaturePlot(t_mes, reduction = "sp", features = c("COL1A1", "CD3E"))
# Change T-cells (Mes) to T-cells
# Change T-cells (Epi) to T-cells

# Fix fibroblast labels on epithelial and immune
fb_epi <- subset(second_epi, subset = CT_secondpass_full == "Fibroblasts")
fb_epi$CT_secondpass_full <- "Fibroblasts (Epi)"
fb_imm <- subset(second_imm, subset = CT_secondpass_full == "Fibroblasts")
fb_imm$CT_secondpass_full <- "Fibroblasts (Imm)"
fb_test <- merge(second_mes, y = list(fb_epi, fb_imm))
DotPlot(fb_test, features = mesenchymal_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(fb_test, features = epithelial_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Change Fibroblasts (Epi) to Alevolar FBs (Healthy)
# Change Fibroblasts (Imm) to Fibroblasts

# Fix myeloid labels on epithelial
my_epi <- subset(second_epi, subset = CT_secondpass_full %in% c("Macrophages",
                                                                "Interstitial Macrophages",
                                                                "Interstitial Macrophages"))
my_epi$CT_secondpass_full[my_epi$CT_secondpassB == "Macrophages"] <- "Macrophages (Epi)"
my_epi$CT_secondpass_full[my_epi$CT_secondpassB == "Interstitial Macrophages"] <- "Interstitial Macrophages (Epi)"
my_epi$CT_secondpass_full[my_epi$CT_secondpassB == "Interstitial Macrophages"] <- "Interstitial Macrophages (Epi)"
my_test <- merge(second_imm, y = my_epi)
DotPlot(my_test, features = immune_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(my_test, features = epithelial_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
table(my_test$CT_secondpassB, my_test$sample_type)
DimPlot(my_epi, reduction = "sp", group.by = "CT_secondpass_full")
# Change Interstitial/Derived Macrophages (Epi) to Interstitial Macrophages
# Change Macrophages (Epi) to Macrophages

# Fix labels
full_cell_labels <- merge(second_epi, y = list(second_endo, second_imm, second_mes))
merged_spatial$CT_secondpass_full <- full_cell_labels$CT_secondpass_full
merged_spatial$CT_secondpass_full[merged_spatial$CT_secondpass_full == "Interstitial Macrophages"] <- "Interstitial Macrophages"
epi_alv_fb_cells <- colnames(subset(second_epi, subset = CT_secondpassB == "Fibroblasts"))
merged_spatial$CT_secondpass_full[colnames(merged_spatial) %in% epi_alv_fb_cells] <- "Alveolar FBs (Healthy)"
merged_spatial$CT_secondpass_full[merged_spatial$CT_secondpass_full == "pDcs"] <- "pDCs"
merged_spatial$CT_secondpass_full[merged_spatial$CT_secondpass_full == "Interferon Response"] <- "Fibroblasts (Interferon High)"
merged_spatial$CT_secondpass_full[merged_spatial$CT_secondpass_full == "Lymphoid"] <- "aCap"
merged_spatial$CT_secondpass_full[merged_spatial$CT_secondpass_full == "Proliferating Immune"] <- "Proliferating Macrophages"
DimPlot(merged_spatial, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE) + NoLegend()

DotPlot(merged_spatial, features = epithelial_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(merged_spatial, features = endothelial_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(merged_spatial, features = immune_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(merged_spatial, features = mesenchymal_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



###### FINE-TUNE SCGB3A2+/SFTPC+/KRT17+ ----
# Get the cells of interest from the airway and alveolar objects
fine_epi <- subset(merged_spatial, 
                   subset = CT_secondpass_full %in% c("SCGB+ Secretory", "KRT5-/KRT17+",
                                                      "AT2 (Healthy)", "AT2 (Disease)", 
                                                      "AT2 (Interferon High)", "AT1",
                                                      "Transitional AT2"))
DefaultAssay(fine_epi) <- "RNA"
fine_epi <- NormalizeData(fine_epi)
fine_epi <- ScaleData(fine_epi, features = rownames(fine_epi))
fine_epi <- RunPCA(fine_epi, features = rownames(fine_epi))

# Get PCs and make elbow plot
npcs <- min(get_pcs(fine_epi))
ElbowPlot(fine_epi, ndims = 30) + geom_vline(xintercept = npcs, lty = "dashed")
npcs # 13

# Find neighbors, cluster and UMAP
fine_epi <- FindNeighbors(fine_epi, dims = 1:npcs)
fine_epi <- RunUMAP(fine_epi, dims = 1:npcs)
fine_epi <- FindClusters(fine_epi, resolution = 0.4)
DimPlot(fine_epi, group.by = "RNA_snn_res.0.4", label = TRUE)
DimPlot(fine_epi, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE)
DimPlot(fine_epi, group.by = "sample")
DimPlot(fine_epi, group.by = "sample_type")

fine_epi <- SetIdent(fine_epi, value = "RNA_snn_res.0.4")
fine_markers <- FindAllMarkers(fine_epi, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Look at secretory things
FeaturePlot(fine_epi, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA", "KRT5", 
                                   "KRT17", "KRT8", "SOX4", "SOX9", "RTKN2", "AGER"))
VlnPlot(fine_epi, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA", "KRT5", 
                               "KRT17", "KRT8", "SOX4", "RTKN2", "AGER",
                               "MS4A7", "CD3E", "S100A8"), pt.size = 0)
DotPlot(fine_epi, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA", "KRT5", 
                               "KRT17", "KRT8", "SOX4", "SOX9", "RTKN2", "AGER")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Look at where they are located
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "0"), group.by = "RNA_snn_res.0.4", reduction = "sp") # Airway/cryptic alveoli?
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "1"), group.by = "RNA_snn_res.0.4", reduction = "sp") # Healthy
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "2"), group.by = "RNA_snn_res.0.4", reduction = "sp") # Airway
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "3"), group.by = "RNA_snn_res.0.4", reduction = "sp") # Airway
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "4"), group.by = "RNA_snn_res.0.4", reduction = "sp") # Healthy-ish
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "5"), group.by = "RNA_snn_res.0.4", reduction = "sp") # AT2 Interferon
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "6"), group.by = "RNA_snn_res.0.4", reduction = "sp") # AT2 cells surrounded by macrophages
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "7"), group.by = "RNA_snn_res.0.4", reduction = "sp") # THD0008 (Proliferating)
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "8"), group.by = "RNA_snn_res.0.4", reduction = "sp") # Slight T-cell contamination
DimPlot(subset(fine_epi, subset = RNA_snn_res.0.4 == "9"), group.by = "RNA_snn_res.0.4", reduction = "sp") # THD0008 (AT2 cells surrounded by interstitial macrophages)

# AT1: 4
# AT2-like: 0, 1, 2, 5, 7, 8, 9
# KRT5-/KRT17+ & Secretory-like: 3
table(fine_epi$CT_secondpass_full, fine_epi$RNA_snn_res.0.4)
table(fine_epi$CT_secondpass_full, fine_epi$sub2)

# Subcluster 2 (KRT5-/KRT17+)
fine_epi <- FindSubCluster(fine_epi, cluster = "2", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub2")
DimPlot(fine_epi, group.by = "sub2", label = TRUE)
VlnPlot(fine_epi, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA", "KRT5", 
                           "KRT17", "KRT8", "SPINK1", "SOX4"), 
        pt.size = 0, group.by = "sub2")
sub2 <- subset(fine_epi, subset = RNA_snn_res.0.4 == "2")
VlnPlot(sub2, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA", "KRT5", 
                           "KRT17", "KRT8", "SPINK1", "SOX4"), 
        pt.size = 0, group.by = "sub2")
DotPlot(fine_epi, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA", "KRT5", 
                               "KRT17", "KRT8", "SPINK1", "SOX4"), group.by = "sub2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sub2 <- SetIdent(sub2, value = "sub2")
sub2_markers <- FindMarkers(sub2, ident.1 = "2_0", ident.2 = "2_1") %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
fine_epi <- SetIdent(fine_epi, value = "sub2")
fine_epi <- FindSubCluster(fine_epi, cluster = "2_1", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "sub2_1")
sub2_1 <- subset(fine_epi, subset = sub2 == "2_1")
sub2_1 <- SetIdent(sub2_1, value = "sub2_1")
sub2_1_markers <- FindAllMarkers(sub2_1, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
sub2 <- subset(fine_epi, subset = RNA_snn_res.0.4 == "2")
DimPlot(sub2, reduction = "sp", group.by = "sub2_1")
DimPlot(sub2,  group.by = "sub2_1")
FeaturePlot(sub2,  features = c("KRT5", "KRT17", "KRT8", "SPINK1", "TP63"), split.by = "sub2_1")
DimPlot(sub2, reduction = "sp", group.by = "sub2")
VlnPlot(sub2_1, features = c("KRT17", "KRT8", "SPINK1", "KRT5", "SOX4", "MMP7", "TP63"), pt.size = 0)
DotPlot(fine_epi, features = c("KRT17", "KRT8", "SPINK1", "KRT5", "SOX4", "MMP7", 
                               "TP63", "SCGB3A2", "SFTPC", "TP73", "FOXJ1", "CEACAM5", 
                               "COL4A3", "AGER", "RTKN2"), group.by = "sub2_1")
DotPlot(fine_epi, features = c("KRT17", "KRT8", "SPINK1", "KRT5", "SOX4", "MMP7", 
                               "TP63", "SCGB3A2", "SFTPC", "TP73", "FOXJ1", "CEACAM5", 
                               "COL4A3", "AGER", "RTKN2", "LAMP3"), group.by = "sub0")
VlnPlot(sub2_1, features = c("SFTPC", "NAPSA", "KRT17", "KRT8", "SOX4", "MMP7"), pt.size = 0, group.by = "sub2")

# Subcluster 3 (SCGB+)
fine_epi <- FindSubCluster(fine_epi, cluster = "3", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "sub3")
DimPlot(fine_epi, group.by = "sub3", label = TRUE)
DotPlot(fine_epi, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA", "KRT5", 
                               "KRT17", "KRT8", "SPINK1", "SOX4"), group.by = "sub3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sub3 <- subset(fine_epi, subset = RNA_snn_res.0.4 == "3")
VlnPlot(sub3, features = c("SCGB3A2", "SCGB1A1", "SFTPC", "NAPSA", "KRT5", 
                           "KRT17", "KRT8", "SPINK1", "SOX4"), 
        pt.size = 0, group.by = "sub3")
DimPlot(sub3, group.by = "sub3", reduction = "sp")
FeaturePlot(merged_spatial, features = c("SCGB3A2", "SCGB1A1"), reduction = "sp", pt.size = 0.0001)

# Subcluster AT2 Disease (0)?
fine_epi <- FindSubCluster(fine_epi, cluster = "0", graph.name = "RNA_snn", resolution = 0.15, subcluster.name = "sub0")
DimPlot(fine_epi, group.by = "sub0", label = TRUE)
sub0 <- subset(fine_epi, subset = RNA_snn_res.0.4 == "0")
DimPlot(sub0, reduction = "sp", split.by = "sub0")

# Subcluster AT2 Healthy (1)
fine_epi <- FindSubCluster(fine_epi, cluster = "1", graph.name = "RNA_snn", resolution = 0.15, subcluster.name = "sub1")
DimPlot(fine_epi, group.by = "sub1", label = TRUE)
sub1 <- subset(fine_epi, subset = RNA_snn_res.0.4 == "1")
DimPlot(sub1, reduction = "sp", split.by = "sub1")

# Look further at 6
VlnPlot(fine_epi, features = c("S100A8", "MS4A7"), group.by = "RNA_snn_res.0.4", pt.size = 0)
DotPlot(fine_epi, features = c("SPP1", "FABP4"))
sub6 <- subset(fine_epi, subset = RNA_snn_res.0.4 == "6")
FeaturePlot(sub6, features = c("SPP1", "FABP4", "S100A8"), pt.size = 0.0001)

# Label cells
fine_epi$fine_CT <- ""
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "0"] <- "AT2 0"
fine_epi$fine_CT[fine_epi$sub1 == "1_0"] <- "AT2 (Healthy)"
fine_epi$fine_CT[fine_epi$sub1 == "1_1"] <- "AT2 (Disease)"
fine_epi$fine_CT[fine_epi$sub2 == "2_0"] <- "Transitional AT2"
fine_epi$fine_CT[fine_epi$sub2_1 == "2_1_0"] <- "KRT5-/KRT17+"
fine_epi$fine_CT[fine_epi$sub2_1 == "2_1_1"] <- "KRT5-/KRT17+"
fine_epi$fine_CT[fine_epi$sub3 == "3_0"] <- "SCGB3A2+/SCGB1A1+"
fine_epi$fine_CT[fine_epi$sub3 == "3_1"] <- "SCGB3A2+"
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "4"] <- "AT1"
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "5"] <- "AT2 (Interferon High)"
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "6"] <- "Macrophages"
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "7"] <- "Proliferating Epithelial"
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "8"] <- "AT2 8"
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "9"] <- "AT2 (Healthy)"
DimPlot(fine_epi, group.by = "fine_CT", label = TRUE)
DimPlot(subset(xenium, subset = sample == "VUILD104MF"), reduction = "sp")
group.by = "fine_CT", cols = c(rep("grey", 5), "red", "blue", rep("grey", 5)))
DimPlot(fine_epi, group.by = "CT_secondpassA", label = TRUE)
DimPlot(fine_epi, group.by = "CT_firstpass", label = TRUE)
table(fine_epi$CT_secondpass_full, fine_epi$RNA_snn_res.0.4)
table(fine_epi$fine_CT, fine_epi$sample_type)
table(fine_epi$sub1, fine_epi$sample_type)
DotPlot(fine_epi, features = c("TGFB1", "TGFB2", "TGFB3", "SPINK1", "HES1", "SMAD4", 
                               "NOX4", "SOD2", "GCLM", "DNAJB9", "DUOX1", "HMOX1", "PRDX4"), group.by = "fine_CT") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(fine_epi, features = c(epithelial_features, "MKI67", "TOP2A"), group.by = "fine_CT") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(fine_epi, features = immune_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(fine_epi, features = mesenchymal_features, group.by = "fine_CT") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(fine_epi, features = endothelial_features, group.by = "fine_CT") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Markers for AT2 cells
at2_sub <- subset(fine_epi, subset = fine_CT %in% c("AT2 0", "AT2 (Healthy)",
                                                    "AT2 (Disease)", "AT2 8", "Transitional AT2"))
at2_sub <- SetIdent(at2_sub, value = "fine_CT")
at2_markers <- FindAllMarkers(at2_sub, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(at2_sub, features = c(epithelial_features, "MKI67", "TOP2A"), group.by = "fine_CT") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(at2_sub, features = c("SFTPC", "SCGB3A2", "AGER", "RTKN2", "NAPSA", "KRT8", "SOX4"), pt.size = 0)

# Re-label cells
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "8"] <- "AT2"
fine_epi$fine_CT[fine_epi$RNA_snn_res.0.4 == "0"] <- "AT2"
DimPlot(fine_epi, group.by = "fine_CT", label = TRUE, raster = FALSE, repel = TRUE)
DotPlot(fine_epi, features = c(epithelial_features, "MKI67", "TOP2A"), group.by = "fine_CT") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# saveRDS(fine_epi, "/scratch/avannan/xenium_files/fine_epi_04272023.rds")
fine_epi <- readRDS("/scratch/avannan/xenium_files/fine_epi_04272023.rds")

DimPlot(fine_epi, group.by = "fine_CT", label = TRUE)


## RELABEL CELLS IN MAIN OBJECT ----
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "Transitional AT2"))] <- "Transitional AT2"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "SCGB3A2+/SCGB1A1+"))] <- "SCGB3A2+/SCGB1A1+"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "SCGB3A2+"))] <- "SCGB3A2+"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "Proliferating Epithelial"))] <- "Proliferating Epithelial"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "Macrophages"))] <- "Macrophages"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "KRT5-/KRT17+"))] <- "KRT5-/KRT17+"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "AT2 (Interferon High)"))] <- "AT2 (Interferon High)"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "AT2 (Healthy)"))] <- "AT2 (Healthy)"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "AT2 (Disease)"))] <- "AT2 (Disease)"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "AT2"))] <- "AT2"
merged_spatial$CT_secondpass_full[colnames(subset(fine_epi, subset = fine_CT == "AT1"))] <- "AT1"
DimPlot(merged_spatial, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE) + NoLegend()

# Dot plots
DotPlot(merged_spatial, features = epithelial_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(merged_spatial, features = endothelial_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(merged_spatial, features = immune_features, group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(merged_spatial, features = c(mesenchymal_features, "HAS1"), group.by = "CT_secondpass_full") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Label lineages
merged_spatial$lineage <- ""
merged_spatial$lineage[
  merged_spatial$CT_secondpass_full %in% c("AT1", "AT2", "AT2 (Healthy)", "AT2 (Disease)",
                                           "AT2 (Interferon High)", "Basal", "Basal (Interferon High)",
                                           "Ciliated", "Differentiating Ciliated", "KRT5-/KRT17+",
                                           "MUC5B+ Secretory", "MUC5B+ Secretory (Interferon High)",
                                           "Proliferating Epithelial", "SCGB3A2+", "SCGB3A2+/SCGB1A1+",
                                           "Transitional AT2")] <- "Epithelial"
merged_spatial$lineage[
  merged_spatial$CT_secondpass_full %in% c("aCap", "Arteriole", "gCap", "Proliferating Endothelial", "Venous")] <- "Endothelial"
merged_spatial$lineage[
  merged_spatial$CT_secondpass_full %in% c("B cells", "CD4+ T-cells", "CD8+ T-cells", "CD4+ Tregs",
                                           "cDCs", "DCs", "DCs - CCL22+", "FABP4+ Macrophages",
                                           "Macrophages", "Macrophages (Interferon High)", "Mast",
                                           "Interstitial Macrophages", "Interstitial Macrophages (FCN1+)",
                                           "Interstitial Macrophages (Interferon High)", "NK cells",
                                           "pDCs", "Plasma", "Proliferating B cells", "Proliferating Macrophages",
                                           "Proliferating T-cells", "Proliferating NK cells", "SPP1+ Macrophages",
                                           "T-cells", "Proliferating Plasma")] <- "Immune"
merged_spatial$lineage[
  merged_spatial$CT_secondpass_full %in% c("Activated FBs (CTHRC1+/FAP+)", "Adventitial FBs",
                                           "Alveolar FBs (Healthy)", "Alveolar FBs (Disease)",
                                           "Fibroblasts", "Fibroblasts (Interferon High)",
                                           "HAS1+ Fibroblasts", "Lipofibroblasts", "Mesothelial",
                                           "Peribronchial FBs", "PLIN2+ Fibroblasts", 
                                           "Proliferating Mesenchymal", "SMCs")] <- "Mesenchymal"
merged_spatial$lineage[merged_spatial$CT_secondpass_full == "DCs/Lymphatic"] <- "Endothelial" # This is really Immune/Endothelial
DimPlot(merged_spatial, group.by = "lineage", label = TRUE, repel = TRUE)


# saveRDS(fine_epi, "/scratch/avannan/xenium_files/second_mes_secondpass_04272023.rds")
# saveRDS(merged_spatial, "/scratch/avannan/xenium_files/xenium_secondpass_full_04272023.rds")
xenium <- readRDS("/scratch/avannan/xenium_files/xenium_secondpass_full_04272023.rds")


## REMAKE UMAPS FOR EACH LINEAGE ----
epi_new <- subset(xenium, subset = lineage == "Epithelial")
endo_new <- subset(xenium, subset = lineage == "Endothelial")
imm_new <- subset(xenium, subset = lineage == "Immune")
mes_new <- subset(xenium, subset = lineage == "Mesenchymal")

DefaultAssay(epi_new) <- "RNA"
epi_new <- NormalizeData(epi_new)
epi_new <- ScaleData(epi_new, features = rownames(epi_new))
epi_new <- RunPCA(epi_new, features = rownames(epi_new))
npcs <- min(get_pcs(epi_new))
epi_new <- FindNeighbors(epi_new, dims = 1:npcs)
epi_new <- RunUMAP(epi_new, dims = 1:npcs)
DimPlot(epi_new, group.by = "CT_secondpass_full")

DefaultAssay(endo_new) <- "RNA"
endo_new <- NormalizeData(endo_new)
endo_new <- ScaleData(endo_new, features = rownames(endo_new))
endo_new <- RunPCA(endo_new, features = rownames(endo_new))
npcs <- min(get_pcs(endo_new))
endo_new <- FindNeighbors(endo_new, dims = 1:npcs)
endo_new <- RunUMAP(endo_new, dims = 1:npcs)
DimPlot(endo_new, group.by = "CT_secondpass_full", label = TRUE)

DefaultAssay(imm_new) <- "RNA"
imm_new <- NormalizeData(imm_new)
imm_new <- ScaleData(imm_new, features = rownames(imm_new))
imm_new <- RunPCA(imm_new, features = rownames(imm_new))
npcs <- min(get_pcs(imm_new))
imm_new <- FindNeighbors(imm_new, dims = 1:npcs)
imm_new <- RunUMAP(imm_new, dims = 1:npcs)
DimPlot(imm_new, group.by = "CT_secondpass_full", label = TRUE)

DefaultAssay(mes_new) <- "RNA"
mes_new <- NormalizeData(mes_new)
mes_new <- ScaleData(mes_new, features = rownames(mes_new))
mes_new <- RunPCA(mes_new, features = rownames(mes_new))
npcs <- min(get_pcs(mes_new))
mes_new <- FindNeighbors(mes_new, dims = 1:npcs)
mes_new <- RunUMAP(mes_new, dims = 1:npcs)
DimPlot(mes_new, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE) + pretty_umap + ggtitle("")

xenium@meta.data %>%
  filter(sample_type != "ILD") %>%
  ggplot(aes(x = sample_type, fill = lineage)) +
  geom_bar(position = position_fill())
mes_new@meta.data %>%
  filter(sample_type != "ILD") %>%
  ggplot(aes(x = sample_type, fill = CT_secondpass_full)) +
  geom_bar(position = position_fill())

# saveRDS(epi_new, "/scratch/avannan/xenium_files/epi_new_umap_050123.rds")
# saveRDS(endo_new, "/scratch/avannan/xenium_files/endo_new_umap_050123.rds")
# saveRDS(imm_new, "/scratch/avannan/xenium_files/imm_new_umap_050123.rds")
# saveRDS(mes_new, "/scratch/avannan/xenium_files/mes_new_umap_050123.rds")
mes_new <- readRDS("/scratch/avannan/xenium_files/mes_new_umap_050123.rds")

DimPlot(epi_new, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE) + ggtitle("")
DimPlot(endo_new, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE) + ggtitle("")
DimPlot(imm_new, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE) + ggtitle("")
DimPlot(mes_new, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE) + ggtitle("")
DimPlot(xenium, group.by = "CT_secondpass_full", label = TRUE, repel = TRUE) + ggtitle("")
DimPlot(xenium, group.by = "lineage", label = TRUE) + ggtitle("") + NoLegend()


## FIX SOME CELL TYPE LABELS ----
aCap_cells <- colnames(xenium)[xenium$broad_CT5 == "aCap"]
gCap_cells <- colnames(xenium)[xenium$broad_CT5 == "gCap"]

xenium$broad_CT5[colnames(xenium) %in% aCap_cells] <- "gCap"
xenium$fine_CT4[colnames(xenium) %in% aCap_cells] <- "gCap"
xenium$fine_CT3[colnames(xenium) %in% aCap_cells] <- "gCap"
xenium$fine_CT2[colnames(xenium) %in% aCap_cells] <- "gCap"
xenium$finest_CT1[colnames(xenium) %in% aCap_cells] <- "gCap"

xenium$broad_CT5[colnames(xenium) %in% gCap_cells] <- "aCap"
xenium$fine_CT4[colnames(xenium) %in% gCap_cells] <- "aCap"
xenium$fine_CT3[colnames(xenium) %in% gCap_cells] <- "aCap"
xenium$fine_CT2[colnames(xenium) %in% gCap_cells] <- "aCap"
xenium$finest_CT1[colnames(xenium) %in% gCap_cells] <- "aCap"

xenium$broad_CT5[xenium$broad_CT5 == "MUC5B+ Secretory"] <- "MUC5B+"

xenium$broad_CT5[xenium$broad_CT5 == "Monocyte-derived Macrophages"] <- "Interstitial Macrophages"
xenium$fine_CT3[xenium$fine_CT3 == "Monocyte-derived Macrophages"] <- "Interstitial Macrophages"
xenium$fine_CT3[xenium$fine_CT3 == "Monocyte-derived Macrophages (FCN1+)"] <- "Interstitial Macrophages (FCN1+)"
xenium$fine_CT2[xenium$fine_CT2 == "Monocyte-derived Macrophages"] <- "Interstitial Macrophages"
xenium$fine_CT2[xenium$fine_CT2 == "Monocyte-derived Macrophages (FCN1+)"] <- "Interstitial Macrophages (FCN1+)"
xenium$finest_CT1[xenium$finest_CT1 == "Monocyte-derived Macrophages"] <- "Interstitial Macrophages"
xenium$finest_CT1[xenium$finest_CT1 == "Monocyte-derived Macrophages (FCN1+)"] <- "Interstitial Macrophages (FCN1+)"
xenium$finest_CT1[xenium$finest_CT1 == "Monocyte-derived Macrophages (Interferon High)"] <- "Interstitial Macrophages (Interferon High)"
xenium$broad_CT5[which(xenium$broad_CT5 == "DCs/Lymphatic")] <- "Lymphatic/DCs"


## SAVE OBJECTS ----
# saveRDS(epithelial, "/scratch/avannan/xenium_files/epithelial_firstpass_04072023.rds")
# saveRDS(endothelial, "/scratch/avannan/xenium_files/endothelial_firstpass_04072023.rds")
# saveRDS(immune, "/scratch/avannan/xenium_files/immune_firstpass_04072023.rds")
# saveRDS(lymphoid, "/scratch/avannan/xenium_files/lymphoid_firstpass_04072023.rds")
# saveRDS(myeloid, "/scratch/avannan/xenium_files/myeloid_firstpass_04072023.rds")
# saveRDS(mesenchymal, "/scratch/avannan/xenium_files/mesenchymal_firstpass_04072023.rds")
epithelial <- readRDS("/scratch/avannan/xenium_files/epithelial_firstpass_04072023.rds")
endothelial <- readRDS("/scratch/avannan/xenium_files/endothelial_firstpass_04072023.rds")
immune <- readRDS("/scratch/avannan/xenium_files/immune_firstpass_04072023.rds")
lymphoid <- readRDS("/scratch/avannan/xenium_files/lymphoid_firstpass_04072023.rds")
myeloid <- readRDS("/scratch/avannan/xenium_files/myeloid_firstpass_04072023.rds")
mesenchymal <- readRDS("/scratch/avannan/xenium_files/mesenchymal_firstpass_04072023.rds")

