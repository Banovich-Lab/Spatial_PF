## SETTING ENVIRONMENT ----
.libPaths(c("/home/avannan/R/x86_64-pc-linux-gnu-library/4.2", 
            "/home/avannan/R/rstudio-4.3.0-4-with_modules.sif/libs",
            "/usr/local/lib/R/site-library",                        
            "/usr/local/lib/R/library"))

# unloadNamespace("pillar")
# unloadNamespace("vctrs")
library(rlang, lib.loc = "/usr/local/lib/R/site-library")
library(tidyverse, lib.loc = "/usr/local/lib/R/site-library")
library(ggplot2)
library(Matrix)
library(ggpubr)
library(Seurat)
library(SeuratObject)
library(gplots)
library(tibble)
library(randomcoloR)
library(patchwork)

# Set seed
set.seed(0309)
work_dir <- "/scratch/avannan"
options(scipen = 99999)
options(ggrepel.max.overlaps = Inf)
'%!in%' <- function(x,y)!('%in%'(x,y))
filter <- dplyr::filter
select <- dplyr::select
pretty_umap <- theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     axis.title = element_text(hjust = 1))
theme_angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

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

## SET LISTS OF LINEAGE MARKERS TO EXPLORE ----
endothelial_features <- c("APLN", "CA4", "HEY1", "BMPR2", "CD34", "EPAS1", "FCN3", 
                          "GNG11", "PECAM1", "APLNR", "PLVAP", "ACKR1", "CLDN5", 
                          "RAMP2", "KDR")
epithelial_features <- c("DUOX1", "NKX2-1", "AGER", "RTKN2", "NAPSA", "PGC", "SFTA2",
                         "SFTPC", "SFTPD", "KRT14", "KRT15", "KRT5", "KRT6A", "S100A2",
                         "TP63", "KRT17", "AGR3", "C20orf85", "FOXJ1", "DMBT1", "EPCAM", 
                         "KRT18", "MMP7", "FOXI1", "MUC5B", "SCGB1A1", "SCGB3A2", 
                         "WFDC2", "KRT8", "SOX9", "SPINK1", "GKN2", "MMP10", "SOX2",
                         "CDH26", "TP73", "CFTR", "SOX4", "SAA2", "BPIFA1", "CEACAM5",
                         "CEACAM6", "ERN2", "LTF", "CCNA1", "ICAM1", "ITGB1", "TGFB2", "CDKN2A")
immune_features <- c("PPARG", "BANK1", "CD19", "CD79A", "LTB", "MS4A1", "TNFRSF13C", 
                     "CD86", "GZMB", "HLA-DRA", "CCR7", "CXCR4", "PTPRC", "TCL1A", 
                     "CD69", "CD4", "CD8A", "CD8B", "CD2", "CD28", "CD3D", "CD3E", 
                     "CD3G", "FOXP3", "GZMK", "TRAC", "ITM2C", "CD27", "CCL5", "LCK", 
                     "FABP4", "MARCO", "MCEMP1", "SPP1", "FCN1", "S100A12", "S100A8", 
                     "S100A9", "CCL22", "ITGAM", "NFKB1",  "IFIT2", "FGFBP2", "GNLY", 
                     "KLRB1", "KLRC1", "NKG7", "LILRA4", "CD79B", "CXCR5", 
                     "CXCL9", "GPR183", "HLA-DQA1", "KLRG1", "BCL2L11", "CD52", 
                     "SLC1A3", "TNFRSF9", "CTLA4",  "IL2RA", "LAG3", "PDCD1",
                     "PIM2", "IL7R", "LEF1", "FASLG", "HAVCR2", "ISG20", "CPA3", 
                     "KIT", "TPSAB1", "C1QC", "CD68", "MS4A7", "CD14", "FCGR3A", 
                     "FCER1G",  "CD247", "GZMA", "IRF7", "JCHAIN", "LYZ", "CD1C", "CD1A")
# Note that CD4 & ITM2C also mark pericytes, + IL7R and CD14 also mark some endothelial
mesenchymal_features <- c("MFAP5", "PI16", "SFRP2", "ELN", "FAP", "COL1A1", 
                          "COL1A2", "COL3A1", "DCN", "FN1", "HAS1", "HAS2", "LUM", 
                          "MEG3", "SPARCL1", "CTHRC1", "TGFB3", "POSTN", "ACTA2")
# Note that SPARCL1 also marks some endothelial

# random_sample_colors <- distinctColorPalette(45)
random_sample_colors <- c("#E055D4", "#DCEC45", "#9CA45E", "#E587B8", "#98739C", "#689CE2", "#E94B55",
                          "#77B5B9", "#587CE2", "#A3E460", "#9A6CE8", "#62B492", "#F0E198", "#64E6EC",
                          "#5CC2E8", "#DF4D99", "#E28642", "#CFE9AF", "#B2D8EE", "#E2DBE1", "#E5EF86",
                          "#4EE6C7", "#95E9B0", "#A796E5", "#DA39ED", "#E3B5CF", "#DA8BE1", "#A0A996",
                          "#A7B5DC", "#CCF1E3", "#D6A769", "#4A6A83", "#61EC6E", "#68B355", "#9AE4CF",
                          "#E2A99A", "#5AEB9A", "#BD6A70", "#7E3EE3", "#D8B1EB", "#84F533", "#A7EA93",
                          "#DCC140", "#903898", "#E8DDBD")

# Load eQTL file for comparison
eqtl.ref <- readRDS("/scratch/avannan/eQTL_rds_files/ILD_all_celltypes_Seurat.rds")
eqtl.ref <- SetIdent(eqtl.ref, value = "manual_annotation_1")
DotPlot(eqtl.ref, features = endothelial_features, group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(eqtl.ref, features = epithelial_features, group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(eqtl.ref, features = immune_features, group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(eqtl.ref, features = mesenchymal_features, group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


## SUBDIRECTORY NAMES ----
# Get subdirectory names that correspond with sample IDs
work_dir <- "/scratch/avannan/late_IPF_spatial/xenium/REVISION/REVISION_xenium_output/"
id_list <- c(
  # TMA1
  VUHD116A = "relabel_output-XETG00048__0003817__VUHD116A__20230308__003730",
  VUHD116B = "relabel_output-XETG00048__0003817__VUHD116B__20230308__003731",
  VUILD102LF = "relabel_output-XETG00048__0003817__VUILD102LF__20230308__003731",
  VUILD102MF = "relabel_output-XETG00048__0003817__VUILD102MF__20230308__003730",
  VUILD107MF = "relabel_output-XETG00048__0003817__VUILD107MF__20230308__003731",
  VUILD96LF = "relabel_output-XETG00048__0003817__VUILD96LF__20230308__003730",
  VUILD96MF = "relabel_output-XETG00048__0003817__VUILD96MF__20230308__003730",
  # # TMA2
  VUHD069 = "relabel_output-XETG00048__0003789__VUHD069__20230308__003731",
  VUHD095 = "relabel_output-XETG00048__0003789__VUHD095__20230308__003731",
  VUHD113 = "relabel_output-XETG00048__0003789__VUHD113__20230308__003731",
  VUILD104LF = "relabel_output-XETG00048__0003789__VUILD104LF__20230308__003731",
  VUILD104MFVUILD48LFVUILD105LF = "relabel_output-XETG00048__0003789__VUILD104MFVUILD48LFVUILD105LF__20230308__003731",
  VUILD105MF = "relabel_output-XETG00048__0003789__VUILD105MF__20230308__003731",
  VUILD48MF = "relabel_output-XETG00048__0003789__VUILD48MF__20230308__003731",
  # TMA3
  THD0008 = "relabel_output-XETG00048__0003392__THD0008__20230313__191400",
  VUILD106 = "relabel_output-XETG00048__0003392__VUILD106__20230313__191400",
  VUILD110 = "relabel_output-XETG00048__0003392__VUILD110__20230313__191400",
  VUILD115 = "relabel_output-XETG00048__0003392__VUILD115__20230313__191400",
  # TMA4
  THD0011 = "relabel_output-XETG00048__0003400__THD0011__20230313__191400",
  TILD117LF = "relabel_output-XETG00048__0003400__TILD117LF__20230313__191400",
  TILD117MF = "relabel_output-XETG00048__0003400__TILD117MF__20230313__191400",
  TILD175 = "relabel_output-XETG00048__0003400__TILD175__20230313__191400",
  VUILD78LF = "relabel_output-XETG00048__0003400__VUILD78LF__20230313__191400",
  VUILD78MF = "relabel_output-XETG00048__0003400__VUILD78MF__20230313__191400",
  VUILD91LF = "relabel_output-XETG00048__0003400__VUILD91LF__20230313__191400",
  VUILD91MF = "relabel_output-XETG00048__0003400__VUILD91MF__20230313__191400")

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
all_files <- list.files(file.path(work_dir, subdirs, "outs"), full.names = TRUE)
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
  all_transcripts[[sm]] <- transcripts[[sm]][transcripts[[sm]]$qv >= 20, ]
  
  # Find transcripts that overlap a nucleus
  nuc_transcripts[[sm]] <- transcripts[[sm]][transcripts[[sm]]$overlaps_nucleus == "1", ]
  
  # Create cell x gene dataframe
  nuc_transcripts[[sm]] <- as.data.frame(table(nuc_transcripts[[sm]]$cell_id, 
                                               nuc_transcripts[[sm]]$feature_name))
  names(nuc_transcripts[[sm]]) <- c("cell_id", "feature_name", "Count")
  nuc_transcripts[[sm]] <- nuc_transcripts[[sm]] %>% 
    pivot_wider(names_from = "feature_name", values_from = "Count")
  
  # Get counts of negative control probes, negative control codewords, and
  # unassigned codewords count per nucleus
  # Negative Control Probes:
  negcontrolprobe_nuc_ids <- nuc_transcripts[[sm]]$cell_id
  negcontrolprobe_nuc_mat <- nuc_transcripts[[sm]][, grep("NegControlProbe", 
                                                          colnames(nuc_transcripts[[sm]]))]
  negcontrolprobe_nuc_counts <- as.data.frame(rowSums(negcontrolprobe_nuc_mat))
  negcontrolprobe_nuc_counts$cell_id <- negcontrolprobe_nuc_ids
  
  # Negative Control Codewords:
  negcontrolcodeword_nuc_ids <- nuc_transcripts[[sm]]$cell_id
  negcontrolcodeword_nuc_mat <- nuc_transcripts[[sm]][, grep("NegControlCodeword", 
                                                             colnames(nuc_transcripts[[sm]]))]
  negcontrolcodeword_nuc_counts <- as.data.frame(rowSums(negcontrolcodeword_nuc_mat))
  negcontrolcodeword_nuc_counts$cell_id <- negcontrolcodeword_nuc_ids
  
  # Unassigned Codewords:
  unassigned_nuc_ids <- nuc_transcripts[[sm]]$cell_id
  unassigned_nuc_mat <- nuc_transcripts[[sm]][, grep("UnassignedCodeword_", 
                                                     colnames(nuc_transcripts[[sm]]))]
  unassigned_nuc_counts <- as.data.frame(rowSums(unassigned_nuc_mat))
  unassigned_nuc_counts$cell_id <- unassigned_nuc_ids
  
  # Remove negative controls and unassigned, and convert to cell x gene matrix
  nuc_transcripts[[sm]] <- nuc_transcripts[[sm]][, grep("NegControl",
                                                        colnames(nuc_transcripts[[sm]]),
                                                        invert = TRUE)]
  nuc_transcripts[[sm]] <- nuc_transcripts[[sm]][, grep("UnassignedCodeword_",
                                                        colnames(nuc_transcripts[[sm]]),
                                                        invert = TRUE)]
  keep_cells <- nuc_transcripts[[sm]]$cell_id
  nuc_transcripts[[sm]] <- as.data.frame(nuc_transcripts[[sm]])
  rownames(nuc_transcripts[[sm]]) <- keep_cells
  nuc_transcripts[[sm]] <- nuc_transcripts[[sm]][, -1]
  nuc_transcripts[[sm]] <- as.matrix(t(nuc_transcripts[[sm]]))
  
  # Subset nuclear metadata to "cells" with transcripts that overlap nuclei
  updated_metadata[[sm]] <- metadata[[sm]][metadata[[sm]]$cell_id %in% keep_cells, ]
  
  # Add unassigned codeword counts to metadata
  updated_metadata[[sm]] <- full_join(updated_metadata[[sm]], 
                                      negcontrolprobe_nuc_counts,
                                      by = "cell_id") %>%
    full_join(negcontrolcodeword_nuc_counts, by = "cell_id") %>%
    full_join(unassigned_nuc_counts, by = "cell_id")
  updated_metadata[[sm]] <- updated_metadata[[sm]] %>%
    dplyr::rename(num.negcontrolprobe = `rowSums(negcontrolprobe_nuc_mat)`,
                  num.negcontrolcodeword = `rowSums(negcontrolcodeword_nuc_mat)`,
                  num.unassigned = `rowSums(unassigned_nuc_mat)`)
  rownames(updated_metadata[[sm]]) <- updated_metadata[[sm]]$cell_id
}


#### COORDINATE ADJUSTMENT LIST ----
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


smpls <- c("TILD028MF", "TILD049LF", "TILD080MF",  "TILD111LF", "TILD113MF",
           "TILD117MF", "TILD130MF", "TILD167LF", "TILD299LF", "TILD315LF", "VUHD038",
           "VUHD049", "VUHD090", "VUILD141", "VUILD142", "VUILD49", "VUILD58")




#### CREATE SEURAT OBJECTS ----
obj_list <- list()
obj_list <- sapply(sample_ids, function(XX) {
  message(paste("Creating Seurat object for sample", XX))
  
  # Create a Seurat object containing the RNA info 
  sobj <- CreateSeuratObject(counts = nuc_transcripts[[XX]], 
                             assay = "RNA")
  
  # Add metadata
  sobj <- AddMetaData(sobj, metadata = updated_metadata[[XX]])
  sobj$sample <- XX
  sobj$tma <- tmas[[XX]]
  sobj$run <- run_ids[[XX]]
  
  # Calculate quality metrics, including total number of transcripts + negative
  # control probes + codewords + unassigned codewords
  sobj$total_counts_ALL <- sobj$nCount_RNA + sobj$num.negcontrolprobe + 
    sobj$num.negcontrolcodeword + sobj$num.unassigned
  sobj$perc_negcontrolprobe <- (sobj$num.negcontrolprobe/sobj$total_counts_ALL)*100
  sobj$perc_negcontrolcodeword <- (sobj$num.negcontrolcodeword/sobj$total_counts_ALL)*100
  sobj$perc_unassigned <- (sobj$num.unassigned/sobj$total_counts_ALL)*100
  
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


#### FIX MERGED SAMPLES AND MERGE ALL XENIUM DATA ----
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

# Fix cell ids
obj_list[["VUILD104MF"]] <- RenameCells(obj_list[["VUILD104MF"]], add.cell.id = "VUILD104MF")
obj_list[["VUILD48LF"]] <- RenameCells(obj_list[["VUILD48LF"]], add.cell.id = "VUILD48LF")
obj_list[["VUILD105LF"]] <- RenameCells(obj_list[["VUILD105LF"]], add.cell.id = "VUILD105LF")

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

# Merge objects (cannot do spatial DimPlots for this)
merged_spatial_unfiltered <- merge(x = obj_list[[1]], y = obj_list[2:length(obj_list)])

# Add spatial dimension reduction object separately
position_xy <- cbind(merged_spatial_unfiltered$adj_x_centroid,
                     merged_spatial_unfiltered$adj_y_centroid)
row.names(position_xy) <- row.names(merged_spatial_unfiltered@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
merged_spatial_unfiltered[["sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(merged_spatial_unfiltered))


#### REMOVE CELL IDs ----
# Remove cell IDs where the sample they originate from is unclear
# First load in files with those cell ids
rm_nuc_VUILD105LF <- read_csv("/scratch/avannan/late_IPF_spatial/xenium/REVISION/REVISION_remove_nuclei_transcripts/remove_nuclei_VUILD105LF.csv") %>%
  mutate(cell_id = paste0("VUILD105LF_", cell_id)) %>%
  pull(cell_id)
keep_cells <- colnames(merged_spatial_unfiltered)[colnames(merged_spatial_unfiltered) %!in% rm_nuc_VUILD105LF]

# Remove cells
merged_spatial_unfiltered <- subset(merged_spatial_unfiltered, cells = keep_cells)

# View TMAs!
DimPlot(merged_spatial_unfiltered, reduction = "sp", group.by = "sample", label = TRUE)
DimPlot(merged_spatial_unfiltered, reduction = "sp", group.by = "tma", label = TRUE)
DimPlot(merged_spatial_unfiltered, reduction = "sp", group.by = "run", label = TRUE)

# JoinLayers since this is Seurat v5
merged_spatial_unfiltered <- JoinLayers(merged_spatial_unfiltered)

# Save and load file
#saveRDS(merged_spatial_unfiltered, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_spatial_nuclei_only_unfiltered.rds")
merged_spatial_unfiltered <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/xenium_spatial_nuclei_only_unfiltered.rds")




## CREATE OBJECT FOR NEW SLIDE, TMA5 ----
#### LOAD IN DATA ----
path <- "/tgen_labs/banovich/xenium_run_folders/20240508__191935__NB_ProbeDiagnostic/output-XETG00048__0024909__CustomPanel__20240508__193012"
transcript_file <- list.files(path, pattern = "transcripts.csv.gz", full.name = TRUE)
metadata_file <- list.files(path, pattern = "cells.csv.gz", full.name = TRUE)

transcripts <- read_csv(transcript_file, col_types = c(transcript_id = "c", cell_id = "c"))

metadata <- read.delim(metadata_file, sep = ",", colClasses = c(cell_id = "character"))
rownames(metadata) <- metadata$cell_id

# Get transcripts that only overlap the nucleus and create cell x gene matrix
# Also count the number of blanks per cell

# Filter out low quality transcripts 
all_transcripts <- transcripts[transcripts$qv >= 20, ]

# Find transcripts that overlap a nucleus
nuc_transcripts <- transcripts[transcripts$overlaps_nucleus == "1", ]

# Create cell x gene dataframe
nuc_transcripts <- as.data.frame(table(nuc_transcripts$cell_id, 
                                       nuc_transcripts$feature_name))
names(nuc_transcripts) <- c("cell_id", "feature_name", "Count")
nuc_transcripts <- nuc_transcripts %>% 
  pivot_wider(names_from = "feature_name", values_from = "Count")

# Get counts of negative control probes, negative control codewords, and
# unassigned codewords count per nucleus
# Negative Control Probes:
negcontrolprobe_nuc_ids <- nuc_transcripts$cell_id
negcontrolprobe_nuc_mat <- nuc_transcripts[, grep("NegControlProbe", 
                                                  colnames(nuc_transcripts))]
negcontrolprobe_nuc_counts <- as.data.frame(rowSums(negcontrolprobe_nuc_mat))
negcontrolprobe_nuc_counts$cell_id <- negcontrolprobe_nuc_ids

# Negative Control Codewords:
negcontrolcodeword_nuc_ids <- nuc_transcripts$cell_id
negcontrolcodeword_nuc_mat <- nuc_transcripts[, grep("NegControlCodeword", 
                                                     colnames(nuc_transcripts))]
negcontrolcodeword_nuc_counts <- as.data.frame(rowSums(negcontrolcodeword_nuc_mat))
negcontrolcodeword_nuc_counts$cell_id <- negcontrolcodeword_nuc_ids

# Unassigned Codewords:
unassigned_nuc_ids <- nuc_transcripts$cell_id
unassigned_nuc_mat <- nuc_transcripts[, grep("UnassignedCodeword_", 
                                             colnames(nuc_transcripts))]
unassigned_nuc_counts <- as.data.frame(rowSums(unassigned_nuc_mat))
unassigned_nuc_counts$cell_id <- unassigned_nuc_ids

# Remove negative controls and unassigned, and convert to cell x gene matrix
nuc_transcripts <- nuc_transcripts[, grep("NegControl",
                                          colnames(nuc_transcripts),
                                          invert = TRUE)]
nuc_transcripts <- nuc_transcripts[, grep("UnassignedCodeword_",
                                          colnames(nuc_transcripts),
                                          invert = TRUE)]
keep_cells <- nuc_transcripts$cell_id
nuc_transcripts <- as.data.frame(nuc_transcripts)
rownames(nuc_transcripts) <- keep_cells
nuc_transcripts <- nuc_transcripts[, -1]
nuc_transcripts <- as.matrix(t(nuc_transcripts))

# Subset nuclear metadata to "cells" with transcripts that overlap nuclei
updated_metadata <- metadata[metadata$cell_id %in% keep_cells, ]

# Add unassigned codeword counts to metadata
updated_metadata <- full_join(updated_metadata, 
                              negcontrolprobe_nuc_counts,
                              by = "cell_id") %>%
  full_join(negcontrolcodeword_nuc_counts, by = "cell_id") %>%
  full_join(unassigned_nuc_counts, by = "cell_id")
updated_metadata <- updated_metadata %>%
  dplyr::rename(num.negcontrolprobe = `rowSums(negcontrolprobe_nuc_mat)`,
                num.negcontrolcodeword = `rowSums(negcontrolcodeword_nuc_mat)`,
                num.unassigned = `rowSums(unassigned_nuc_mat)`)
rownames(updated_metadata) <- updated_metadata$cell_id


#### CREATE SEURAT OBJECT ----
# Create a Seurat object containing the RNA info 
sobj <- CreateSeuratObject(counts = nuc_transcripts, assay = "RNA")

# Add metadata
sobj <- AddMetaData(sobj, metadata = updated_metadata)
sobj$tma <- "TMA5"
sobj$run <- "Run3"

# Calculate quality metrics, including total number of transcripts + negative
# control probes + codewords + unassigned codewords
sobj$total_counts_ALL <- sobj$nCount_RNA + sobj$num.negcontrolprobe + 
  sobj$num.negcontrolcodeword + sobj$num.unassigned
sobj$perc_negcontrolprobe <- (sobj$num.negcontrolprobe/sobj$total_counts_ALL)*100
sobj$perc_negcontrolcodeword <- (sobj$num.negcontrolcodeword/sobj$total_counts_ALL)*100
sobj$perc_unassigned <- (sobj$num.unassigned/sobj$total_counts_ALL)*100

# Get which cells correspond to which sample and label
path2 <- "/scratch/avannan/late_IPF_spatial/xenium/REVISION/REVISION_remove_nuclei_transcripts/keep_cells/"

smpls <- c("TILD028MF", "TILD049LF", "TILD080MF",  "TILD111LF", "TILD113MF",
           "TILD117MF", "TILD130MF", "TILD167LF", "TILD299LF", "TILD315LF", "VUHD038",
           "VUHD049", "VUHD090", "VUILD141", "VUILD142", "VUILD49", "VUILD58")

smpl_cell_ids_all <- c()
for (smpl in smpls) {
  print(smpl)
  smpl_cell_ids <- read.csv(paste(path2, smpl, "_cells_stats.csv", sep = ""),
                            header = T, comment.char = "#")
  # Add sample 
  smpl_cell_ids$Sample <- smpl
  smpl_cell_ids_sub <- smpl_cell_ids %>% 
    dplyr::select(Cell.ID, Sample)
  colnames(smpl_cell_ids_sub) <- c("cell_id", "sample")
  smpl_cell_ids_all <- rbind(smpl_cell_ids_all, smpl_cell_ids_sub)
}

# Join the two files
cell_meta_merge <- merge(updated_metadata, smpl_cell_ids_all)

# Filter out cells not assigned to a sample
# Since some of the samples overlap they will not be included
sobj@meta.data$cell_id <- rownames(sobj@meta.data)
cells_to_keep <- cell_meta_merge$cell_id
sobj_sub <- subset(sobj, subset = cell_id %in% cells_to_keep)

# Add metadata to filtered object
sobj_meta <- sobj_sub@meta.data
sobj_meta_cell_meta <- full_join(sobj_meta, cell_meta_merge) %>% filter(cell_id %in% cells_to_keep)
rownames(sobj_meta_cell_meta) <- rownames(sobj_sub@meta.data)
sobj_sub@meta.data <- sobj_meta_cell_meta
# Total number of cells: 600224

# Remove cells with 0 nCount_RNA
sobj_sub <- subset(sobj_sub, subset = nCount_RNA != 0)

# Rename cells to add sample ID as prefix and fix other metadata
sobj_sub$sample[sobj_sub$sample == "TILD117MF"] <- "TILD117MFB"
sobj_sub$patient <- gsub("LF|MF}MFB", "", sobj_sub$sample)
sobj_sub <- RenameCells(sobj_sub, add.cell.id = sobj_sub)
sobj_sub$adj_x_centroid <- sobj_sub$x_centroid + 20000
sobj_sub$adj_y_centroid <- sobj_sub$y_centroid*-1

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(sobj_sub$adj_x_centroid, sobj_sub$adj_y_centroid)
row.names(position_xy) <- row.names(sobj_sub@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
sobj_sub[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                         assay = DefaultAssay(sobj_sub))


## FILTERING AND MERGING ----
# Calculate new metric of percent negcontrolprobe+negcontrolcodeword+unassignedcodeword
sobj_sub$perc_negcontrolorunassigned <- ((sobj_sub$num.negcontrolprobe + sobj_sub$num.negcontrolcodeword + sobj_sub$num.unassigned)/
                                           (sobj_sub$num.negcontrolprobe + sobj_sub$num.negcontrolcodeword + sobj_sub$num.unassigned + sobj_sub$nCount_RNA)*100)

merged_spatial_unfiltered$perc_negcontrolorunassigned <- ((merged_spatial_unfiltered$num.negcontrolprobe +
                                                             merged_spatial_unfiltered$num.negcontrolcodeword + merged_spatial_unfiltered$num.unassigned)/merged_spatial_unfiltered$total_counts_ALL)*100


sobj_sub_filtered <- subset(sobj_sub,
                            subset = nCount_RNA >= 12 & nFeature_RNA >= 10 &
                              perc_negcontrolprobe <= 5 & 
                              perc_negcontrolcodeword <= 5 &
                              perc_unassigned <= 5 &
                              perc_negcontrolorunassigned <= 5 & 
                              nucleus_area >= 6 & nucleus_area <= 80)
# Number of nuclei before and after filtering
bf_cells <- table(sobj_sub$sample)
aft_cells <- table(sobj_sub_filtered$sample)
diff_cells <- bf_cells - aft_cells
prop_kept_cells <- round(aft_cells/bf_cells*100, 2)
prop_kept_cells

merged_spatial <- subset(merged_spatial_unfiltered,
                         subset = nCount_RNA >= 12 & nFeature_RNA >= 10 &
                           perc_negcontrolprobe <= 5 & 
                           perc_negcontrolcodeword <= 5 &
                           perc_unassigned <= 5 &
                           perc_negcontrolorunassigned <= 5 & 
                           nucleus_area >= 6 & nucleus_area <= 80)
# Number of nuclei before and after filtering
bf_cells <- table(xenium.obj_sub$sample)
aft_cells <- table(xenium.obj_sub_filtered$sample)
diff_cells <- bf_cells - aft_cells
prop_kept_cells <- round(aft_cells/bf_cells*100, 2)
prop_kept_cells

# More prep for merging
meta_to_keep <- colnames(merged_spatial@meta.data)
meta_to_keep <- meta_to_keep[-which(meta_to_keep == "total_counts_ALL")]
sobj_sub_filtered@meta.data <- sobj_sub_filtered@meta.data %>% select(meta_to_keep)
merged_spatial@meta.data <- merged_spatial@meta.data %>% select(meta_to_keep)

# Merge
full_merged_nucleus_nucleus <- merge(sobj_sub_filtered, merged_spatial)
#saveRDS(full_merged_nucleus_nucleus, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/full_merged_nucleus_nucleus_052824.rds")
ncol(full_merged_nucleus_nucleus) # 1650233



## CELL TYPE ANNOTATION (FIRSTPASS) ----
# Here, cells are given a firstpass annotation and lineage so that they can be further split by lineage and re-annotated more accurately.
clustered_obj <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION/REVISION_rds_files/full_merged_nucleus_nucleus_RAPIDS_clustered_052824.rds")

# Initial clustering and quality plots
DimPlot(clustered_obj, reduction = "umap")
#random_cluster_colors <- distinctColorPalette(23)
random_cluster_colors <-  c("#C3E9C6", "#DCC6A8", "#8966D5", "#DA9658", "#77BDE0", "#D7C047", "#85E851",
                            "#65DFB9", "#828FD5", "#D68FDB", "#5FE990", "#DB49D9", "#DC549F", "#C17480",
                            "#DDB7D8", "#D8EB54", "#EB5659", "#E0E398", "#779285", "#98D884", "#CEDDE4",
                            "#8840E5", "#73E2E3")
DimPlot(clustered_obj, reduction = "umap", group.by = "leiden_res0.5", cols = random_cluster_colors) + coord_equal()
DimPlot(clustered_obj, reduction = "umap", split.by = "leiden_res0.5", ncol = 5) + coord_equal() + NoLegend()
FeaturePlot(clustered_obj, features = c("nCount_RNA", "nFeature_RNA", "perc_negcontrolprobe",
                                        "perc_negcontrolcodeword", "perc_unassigned", 
                                        "perc_negcontrolorunassigned"))  + NoLegend() & coord_equal()
VlnPlot(clustered_obj, features = c("nCount_RNA", "nFeature_RNA", "perc_negcontrolprobe",
                                    "perc_negcontrolcodeword", "perc_unassigned", 
                                    "perc_negcontrolorunassigned"), pt.size = 0, group.by = "leiden_res0.5")
# Fix sample names
clustered_obj$sample[clustered_obj$sample == "TILD117MF" & clustered_obj$tma == "TMA5"] <- "TILD117MFB"
DimPlot(clustered_obj, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(clustered_obj, group.by = "sample", ncol = 8, split.by = "sample", cols = random_sample_colors) + NoLegend() + coord_equal()
DimPlot(clustered_obj, group.by = "tma") + coord_equal()
DimPlot(clustered_obj, group.by = "tma", split.by = "tma") + coord_equal()
DimPlot(clustered_obj, group.by = "run", split.by = "run") + coord_equal()
clustered_obj@meta.data %>%
  ggplot(aes(x = as.factor(leiden_res0.5), fill = tma)) +
  geom_bar(position = position_fill())

# Create sample type variable
clustered_obj$sample_type <- "INT"
clustered_obj$sample_type[grepl("HD", clustered_obj$sample)] <- "Unaffected"
clustered_obj$sample_type[grepl("LF", clustered_obj$sample)] <- "LF"
clustered_obj$sample_type[grepl("MF", clustered_obj$sample)] <- "MF"
clustered_obj$sample_type <- ordered(clustered_obj$sample_type, levels = c("Unaffected", "LF", "MF", "INT"))
DimPlot(clustered_obj, group.by = "sample_type", split.by = "sample_type") + coord_equal()
clustered_obj@meta.data %>%
  ggplot(aes(x = as.factor(leiden_res0.5), fill = sample_type)) +
  geom_bar(position = position_fill())

# Find Markers
clustered_obj <- SetIdent(clustered_obj, value = "leiden_res0.5")
top_markers05 <- FindAllMarkers(clustered_obj, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Look at major lineage markers
FeaturePlot(clustered_obj, features = c("EPCAM", "PECAM1", "PTPRC", "LUM", "DCN", "ACTA2")) & coord_equal()
DotPlot(clustered_obj, features = c("EPCAM", "PECAM1", "PTPRC", "LUM", "DCN", "ACTA2", 
                                    "MSLN", "HAS1", "CPA3", "KIT"), group.by = "leiden_res0.5")
VlnPlot(clustered_obj, features = c("EPCAM", "PECAM1", "PTPRC", "LUM", "DCN", "ACTA2"),
        group.by = "leiden_res0.4", cols = random_cluster_colors, pt.size = 0, ncol = 2)
DotPlot(clustered_obj, features = c("MSLN", "HAS1", # 0: Mesothelial & HAS1+ FBs
                                    "S100A8", "S100A9", "ITGAX", "S100A12", # 1: S100A+ Macrophages
                                    "PLVAP", "PECAM1", "ACKR1", "CD34", # 2: Venous
                                    "CDK1", "TOP2A", "MKI67", # 3: Proliferating mix
                                    "COL3A1", "COL1A2", "COL1A1", "LUM", "DCN", "MEG3", # 4,16: Fibroblasts
                                    "CPA3", "TPSAB1", "KIT", # 5: Mast
                                    "AGR3", "SCGB1A1", "C20orf85", "FOXJ1", "MUC5B", "SCGB3A2", # 6: Airway epithelium
                                    "KRT15", "KRT17", "KRT5", "NKX2-1", "TP63", # 7: Basal
                                    "NAPSA", "CEACAM6", "SFTPD", "SFTPC", "LAMP3", "AGER", "RTKN2", # 8: Alveolar epithelium
                                    "MS4A7", "LYZ", "FCER1G", "CD68", "MRC1", # 9,12: Macrophages
                                    "ACTA2", "PDGFRB", "SPARCL1", "CSPG4", # 10: SMCs/Pericytes
                                    "EPAS1", "FCN3", "HEY1", # 11: Capillary/Arteriole
                                    "MS4A1", "TNFRSF13C", "BANK1", # 13: B cells
                                    "PIM2", "JCHAIN", "HERPUD1", "SEC11C", # 14: Plasma/pDCs
                                    "IFIT1", "IFIT3", "OAS2", # 15: Interferon-high Epithelial
                                    # 17: Proliferating Mesenchymal
                                    # 18: AT2
                                    # 19: Proliferating Plasma/pDCs/B cells
                                    "TRAC", "CD3E", "CD3D", # 20: T-cells
                                    "NKG7", "GNLY", "GZMA", "GZMB" # 21: NK cells
                                    # 22: Proliferating T-cells
), group.by = "leiden_res0.5") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

# Check subclustering for cluster 3
sub3 <- subset(clustered_obj, subset = leiden_res0.5 == "3")
DimPlot(sub3, group.by = "leiden_res0.5_sub3_res0.1", cols = distinctColorPalette(2)) + coord_equal()
DimPlot(sub3, group.by = "leiden_res0.5_sub3_res0.2", cols = distinctColorPalette(5)) + coord_equal()
DimPlot(sub3, group.by = "leiden_res0.5_sub3_res0.3", cols = distinctColorPalette(5)) + coord_equal()
DimPlot(sub3, group.by = "leiden_res0.5_sub3_res0.4", cols = distinctColorPalette(7)) + coord_equal()
FeaturePlot(sub3, features = c("EPCAM", "PECAM1", "PTPRC", "LUM", "DCN", "MKI67", "TOP2A", "CENPF", "CDK1"))
sub3 <- SetIdent(sub3, value = "leiden_res0.5_sub3_res0.3")
sub3_markers <- FindAllMarkers(sub3, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(sub3, features = c("COL1A2", "COL3A1", "COL1A1", "FN1", "LUM", "MEG3", # 3_0: Proliferating Mesenchymal
                           "SFTPD", "CEACAM6", "NAPSA", "EPCAM", "SCGB3A2", "C20orf85", # 3_1: Proliferating Epithelial
                           "HIST1H1C", "LYZ", "HLA-DQA1", "HLA-DRA", "MS4A7", "MRC1", "ITGAX", # 3_2, 3_3: Proliferating Macrophages 
                           "IL7R", "TRAC", "CD3E" # 3_4: Proliferating T-cells
), group.by = "leiden_res0.5_sub3_res0.3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
VlnPlot(sub3, group.by = "leiden_res0.5_sub3_res0.3", pt.size = 0, features = c("MKI67", "TOP2A", "CENPF", "CDK1"), ncol = 2, layer = "counts")

# Create initial CT variable
clustered_obj$initial_CT <- ""
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "0"] <- "Mesothelial/HAS1+ FBs"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "1"] <- "S100A+ Macrophages"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "2"] <- "Venous"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5_sub3_res0.3 == "3,0" |
                           clustered_obj$leiden_res0.5 == "17"] <- "Proliferating Mesenchymal"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5_sub3_res0.3 == "3,1"] <- "Proliferating Epithelial"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5_sub3_res0.3 %in% c("3,2", "3,3")] <- "Proliferating Macrophages"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5_sub3_res0.3 == "3,4" | 
                           clustered_obj$leiden_res0.5 == "22"] <- "Proliferating T-cells"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 %in% c("4", "16")] <- "Fibroblasts"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "5"] <- "Mast"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "6"] <- "Airway Epithelium"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "7"] <- "Basal"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "8"] <- "Alveolar Epithelium"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 %in% c("9", "12")] <- "Macrophages"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "10"] <- "SMCs/Pericytes"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "11"] <- "Capillary/Arteriole"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "13"] <- "B cells"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "14"] <- "Plasma/pDCs"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "15"] <- "Interferon-high Epithelial"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "18"] <- "AT2"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "19"] <- "Proliferating Plasma/pDCs/B cells"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "20"] <- "T-cells"
clustered_obj$initial_CT[clustered_obj$leiden_res0.5 == "21"] <- "NK cells"
# random_ct_colors <- distinctColorPalette(23)
random_ct_colors <- c("#92D260", "#7DEE42", "#E6D8BE", "#CAD9E5", "#DDE34B", "#736C91", "#5EEA90",
                      "#61E0BC", "#9E3FE8", "#D94BC0", "#AFE99F", "#DFB6D6", "#C0ECD2", "#DE607D",
                      "#8997E6", "#E1924D", "#8B9E83", "#7CB7DE", "#6ADCE3", "#C9907D", "#E1D98C",
                      "#D68BD5", "#8866D8")
DimPlot(clustered_obj, group.by = "initial_CT", cols = random_ct_colors) + coord_equal()
DotPlot(clustered_obj, features = c("MSLN", "HAS1", # 0: Mesothelial & HAS1+ FBs
                                    "S100A8", "S100A9", "ITGAX", "S100A12", # 1: S100A+ Macrophages
                                    "PLVAP", "PECAM1", "ACKR1", "CD34", # 2: Venous
                                    "CDK1", "TOP2A", "MKI67", # 3: Proliferating mix
                                    "COL3A1", "COL1A2", "COL1A1", "LUM", "DCN", "MEG3", # 4,16: Fibroblasts
                                    "CPA3", "TPSAB1", "KIT", # 5: Mast
                                    "AGR3", "SCGB1A1", "C20orf85", "FOXJ1", "MUC5B", "SCGB3A2", # 6: Airway epithelium
                                    "KRT15", "KRT17", "KRT5", "NKX2-1", "TP63", # 7: Basal
                                    "NAPSA", "CEACAM6", "SFTPD", "SFTPC", "LAMP3", "AGER", "RTKN2", # 8: Alveolar epithelium
                                    "MS4A7", "LYZ", "FCER1G", "CD68", "MRC1", # 9,12: Macrophages
                                    "ACTA2", "PDGFRB", "SPARCL1", "CSPG4", # 10: SMCs/Pericytes
                                    "EPAS1", "FCN3", "HEY1", # 11: Capillary/Arteriole
                                    "MS4A1", "TNFRSF13C", "BANK1", # 13: B cells
                                    "PIM2", "JCHAIN", "HERPUD1", "SEC11C", # 14: Plasma/pDCs
                                    "IFIT1", "IFIT3", "OAS2", # 15: Interferon-high Epithelial
                                    # 17: Proliferating Mesenchymal
                                    # 18: AT2
                                    # 19: Proliferating Plasma/pDCs/B cells
                                    "TRAC", "CD3E", "CD3D", # 20: T-cells
                                    "NKG7", "GNLY", "GZMA", "GZMB" # 21: NK cells
                                    # 22: Proliferating T-cells
), group.by = "initial_CT") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

# Create intial lineage variable
clustered_obj$initial_lineage <- ""
clustered_obj$initial_lineage[clustered_obj$initial_CT %in% c("Capillary/Arteriole", "Venous")] <- "Endothelial"
clustered_obj$initial_lineage[clustered_obj$initial_CT %in% c("Airway Epithelium", "Alveolar Epithelium", "AT2", "Basal", 
                                                              "Interferon-high Epithelial", "Proliferating Epithelial")] <- "Epithelial"
clustered_obj$initial_lineage[clustered_obj$initial_CT %in% c("Fibroblasts", "Mesothelial/HAS1+ FBs",
                                                              "Proliferating Mesenchymal", "SMCs/Pericytes")] <- "Mesenchymal"
clustered_obj$initial_lineage[clustered_obj$initial_CT %in% c("B cells", "Macrophages", "Mast", "NK cells", "Plasma/pDCs", 
                                                              "Proliferating Macrophages", "Proliferating Plasma/pDCs/B cells",
                                                              "Proliferating T-cells", "S100A+ Macrophages", "T-cells")] <- "Immune"
DimPlot(clustered_obj, group.by = "initial_lineage", cols = c("orange", "green", "blue", "red")) + coord_equal()


# Fix cell names
clustered_obj <- RenameCells(clustered_obj, new.names = paste0(clustered_obj$sample, "_", clustered_obj$cell_id))

# Save object
# saveRDS(clustered_obj, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/full_merged_nucleus_nucleus_RAPIDS_clustered_annotated_052824.rds")
clustered_annotated_obj <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/full_merged_nucleus_nucleus_RAPIDS_clustered_annotated_052824.rds")

# Split into lineage-specific objects to cluster
endo <- subset(clustered_obj, subset = initial_lineage == "Endothelial")
epi <- subset(clustered_obj, subset = initial_lineage == "Epithelial")
imm <- subset(clustered_obj, subset = initial_lineage == "Immune")
mes <- subset(clustered_obj, subset = initial_lineage == "Mesenchymal")
ncol(endo) + ncol(epi) + ncol(imm) + ncol(mes) == ncol(clustered_obj) # Must be TRUE
# saveRDS(endo, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/endothelial_lineage_split_nucnuc_052824.rds")
# saveRDS(epi, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_lineage_split_nucnuc_052824.rds")
# saveRDS(imm, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lineage_split_nucnuc_052824.rds")
# saveRDS(mes, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/mesenchymal_lineage_split_nucnuc_052824.rds")

ncol(endo)
ncol(endo1)

ncol(mes)
ncol(mes1)

ncol(air) + ncol(alv)
ncol(epi1) - (ncol(air) + ncol(alv))


ncol(lym) + ncol(mye)
ncol(imm1) - (ncol(lym) + ncol(mye))
2288 + 2223



#### ENDOTHELIAL ----
endo <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/endothelial_lineage_split_nucnuc_RAPIDS_clustered_052824.rds")
DimPlot(endo, group.by = "leiden_res0.4", label = TRUE) + coord_equal()
DimPlot(endo, group.by = "leiden_res0.4", split.by = "leiden_res0.4") + coord_equal()
DimPlot(endo, group.by = "initial_CT", label = TRUE) + coord_equal()

# Add spatial dimension reduction object separately
position_xy <- cbind(endo$adj_x_centroid, endo$adj_y_centroid)
row.names(position_xy) <- row.names(endo@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
endo[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_", assay = DefaultAssay(endo))

# Basic DimPlots and FeaturePlots
DimPlot(endo, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(endo, group.by = "tma") + coord_equal()
DimPlot(endo, group.by = "run", split.by = "run") + coord_equal()
DimPlot(endo, group.by = "sample_type", split.by = "sample_type") + coord_equal()
FeaturePlot(endo, features = c("nCount_RNA", "nFeature_RNA")) & coord_equal()

# Find Markers
endo <- SetIdent(endo, value = "leiden_res0.4")
endo_markers <- FindAllMarkers(endo, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Create and view module scores
endo$endo_score <- colSums(endo@assays$RNA@counts[endothelial_features, ])/colSums(endo@assays$RNA@counts)
endo$epi_score <- colSums(endo@assays$RNA@counts[epithelial_features, ])/colSums(endo@assays$RNA@counts)
endo$imm_score <- colSums(endo@assays$RNA@counts[immune_features, ])/colSums(endo@assays$RNA@counts)
endo$mes_score <- colSums(endo@assays$RNA@counts[mesenchymal_features, ])/colSums(endo@assays$RNA@counts)
FeaturePlot(endo, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(endo, features = c("endo_score", "epi_score", "imm_score", "mes_score"), pt.size = 0, ncol = 2)

# Look at lineage and top markers
FeaturePlot(endo, features = c("EPCAM", "PECAM1", "PTPRC", "DCN", "LUM")) & coord_equal()
FeaturePlot(endo, features = c("PLVAP", "ACKR1", "POSTN", "COL15A1", # 1
                               "CCL21", "FABP4", "MRC1", "GNG11", # 0
                               "PTPRC", "TRAC", "CD3E", # 2
                               "IL7R", "FCN3", # 3
                               # 4 No markers for 4
                               "CA4", "ITGA3", "KDR", # 5
                               "HEY1", "VEGFA", # 6
                               "APLN", "APLNR" # aCap/gCap
)) & coord_equal()
DotPlot(endo, features = c("PLVAP", "ACKR1", "POSTN", "COL15A1", # 1
                           "CCL21", "FABP4", "MRC1", "GNG11", # 0
                           "PTPRC", "TRAC", "CD3E", # 2
                           "IL7R", "FCN3", # 3
                           # 4 No markers for 4
                           "CA4", "ITGA3", "KDR", # 5
                           "HEY1", "VEGFA", "HERPUD1", # 6
                           "APLN", "APLNR" # aCap/gCap
)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
FeaturePlot(endo, features = c("CCL21", "FABP4", "MRC1", "GNG11")) & coord_equal() # 0: Lymphatic
FeaturePlot(endo, features = c("PLVAP", "ACKR1", "POSTN", "COL15A1")) & coord_equal() # 1: Venous
FeaturePlot(endo, features = c("PLVAP", "ACKR1", "POSTN", "COL15A1"), max.cutoff = 1) & coord_equal() # 1: Venous
FeaturePlot(endo, features = c("PTPRC", "TRAC", "CD3E")) & coord_equal() # 2: T-cells
FeaturePlot(endo, features = c("PTPRC", "TRAC", "CD3E"), split.by = "leiden_res0.4") & coord_equal() # 2: T-cells
FeaturePlot(endo, features = c("IL7R", "FCN3", "CA4", "CLDN5", "APLN", "APLNR")) & coord_equal() # 3, 4, & 5: Capillary
FeaturePlot(endo, features = c("IL7R", "FCN3", "CA4", "CLDN5", "APLN", "APLNR"), max.cutoff = 1) & coord_equal() # 3, 4, & 5: Capillary
FeaturePlot(endo, features = c("IL7R", "FCN3", "CA4", "APLN", "APLNR"), split.by = "leiden_res0.4", pt.size = 0.001) & coord_equal() # 3, 4, & 5: Capillary
FeaturePlot(endo, features = c("HEY1", "VEGFA", "HERPUD1")) & coord_equal() # 6: Arteriole
VlnPlot(endo, features = c("IL7R", "FCN3", # 3
                           "CA4", "ITGA3", # 5
                           "APLN", "APLNR" # aCap/gCap
), pt.size = 0)

# Check all genes
DotPlot(endo, features = rownames(endo)[1:115]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(endo, features = rownames(endo)[116:231]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(endo, features = rownames(endo)[232:343]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# Spot-check genes
FeaturePlot(endo, features = c("SFTPC", "AGER", "S100A8", "S100A9",
                               "COL1A1", "COL1A2", "COL3A1")) & coord_equal() # 3, 4, & 5: Capillary

###### VIEW SUBCLUSTERING RESULTS ----
# Subcluster 1 into peribronchial (top) and venous (bottom)
DimPlot(endo, group.by = "leiden_res0.4_sub1_res0.4", label = TRUE) + coord_equal()
VlnPlot(endo, features = c("ACKR1", "COL15A1"), group.by = "leiden_res0.4_sub1_res0.4", pt.size = 0)
FeaturePlot(endo, features = c("ACKR1", "COL15A1"), split.by = "leiden_res0.4_sub1_res0.4", pt.size = 0.001) & coord_equal()
DotPlot(endo, features = c("ACKR1", "COL15A1", "PLVAP"), group.by = "leiden_res0.4_sub1_res0.4")
# Compare to eQTL
eqtl.endo <- subset(eqtl.ref, subset = lineage == "Endothelial")
DotPlot(eqtl.endo, features = c("ACKR1", "COL15A1", "PLVAP"), group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Subcluster 2 to get T-cells
sub2 <- subset(endo, subset = leiden_res0.4 == "2")
DimPlot(sub2, group.by = "leiden_res0.4_sub2_res0.4", label = TRUE) + coord_equal()
FeaturePlot(sub2, features = c("CD3E", "TRAC"), split.by = "leiden_res0.4_sub2_res0.4", pt.size = 0.001) & coord_equal()
FeaturePlot(sub2, features = c("endo_score", "epi_score", "imm_score", "mes_score"), split.by = "leiden_res0.4_sub2_res0.4", pt.size = 0.001) & coord_equal()
VlnPlot(sub2, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        same.y.lims = TRUE, group.by = "leiden_res0.4_sub2_res0.4", ncol = 2, pt.size = 0)
DotPlot(sub2, features = rownames(endo)[1:115], group.by = "leiden_res0.4_sub2_res0.4") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(sub2, features = rownames(endo)[116:231], group.by = "leiden_res0.4_sub2_res0.4") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(sub2, features = rownames(endo)[232:343], group.by = "leiden_res0.4_sub2_res0.4") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
# These cells are not strongly marked by T-cell markers in comparison to endothelial
# What cell type in 2?
FeaturePlot(endo, features = c("PLVAP", "ACKR1", "COL15A1", "CD3E"), split.by = "leiden_res0.4", pt.size = 0.001) & coord_equal()

# Subcluster 3 to get fibroblasts
DimPlot(endo, group.by = "leiden_res0.4_sub3_res0.3", label = TRUE) + coord_equal()
sub3 <- subset(endo, subset = leiden_res0.4 == "3")
VlnPlot(sub3, features = c("endo_score", "epi_score", "imm_score", "mes_score"), group.by = "leiden_res0.4_sub3_res0.3", ncol = 2, pt.size = 0)
DimPlot(sub3, group.by = "leiden_res0.4_sub3_res0.3", label = TRUE) + coord_equal()
FeaturePlot(sub3, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(sub3, features = c("endo_score", "epi_score", "imm_score", "mes_score"), group.by = "leiden_res0.4_sub3_res0.2", ncol = 2, pt.size = 0)
FeaturePlot(sub3, features = "SPARCL1") & coord_equal()
DotPlot(sub3, group.by = "leiden_res0.4_sub3_res0.3", features = c(endothelial_features, mesenchymal_features)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
VlnPlot(sub3, group.by = "leiden_res0.4_sub3_res0.3", features = c(endothelial_features, mesenchymal_features), pt.size = 0, ncol = 6)
vuild91lf <- subset(sub3, subset = sample == "VUILD91LF")
DimPlot(vuild91lf, group.by = "leiden_res0.4_sub3_res0.3", reduction = "sp", cols = distinctColorPalette(10), pt.size = 1)
# 3,3 = Fibroblasts

# Subcluster 5 to get immune and epithelial clusters out
DimPlot(endo, group.by = "leiden_res0.4_sub5_res0.4_sub4_res0.3", label = TRUE) + coord_equal()
sub5 <- subset(endo, subset = leiden_res0.4 == "5")
DimPlot(sub5, group.by = "leiden_res0.4_sub5_res0.4_sub4_res0.3", label = TRUE) + coord_equal()
FeaturePlot(sub5, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(sub5, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        same.y.lims = TRUE, group.by = "leiden_res0.4_sub5_res0.4_sub4_res0.3", ncol = 2, pt.size = 0)
FeaturePlot(sub5, features = c("PECAM1", "SFTPC"), blend = TRUE) & coord_equal()
FeaturePlot(sub5, features = c("PECAM1", "SFTPC"), split.by = "leiden_res0.4_sub5_res0.4_sub4_res0.3", blend = TRUE, pt.size = 0.0001) & coord_equal()
DotPlot(sub5, group.by = "leiden_res0.4_sub5_res0.4_sub4_res0.3", features = c("SFTPC", "SFTPD", "EPCAM", "CA4", "PECAM1", "S100A8", "S100A9"))
DotPlot(sub5, group.by = "leiden_res0.4_sub5_res0.4_sub4_res0.3", features = c("S100A8", "S100A9"))
sub5 <- SetIdent(sub5, value = "leiden_res0.4_sub5_res0.4_sub4_res0.3")
sub5_markers <- FindAllMarkers(sub5, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(sub5, group.by = "leiden_res0.4_sub5_res0.4_sub4_res0.3", features = endothelial_features) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(sub5, group.by = "leiden_res0.4_sub5_res0.4_sub4_res0.3", 
        features = c(endothelial_features, epithelial_features, immune_features)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# Subcluster 4 to split endothelial cell types
DimPlot(endo, group.by = "leiden_res0.4_sub4_res0.2", label = TRUE) + coord_equal()
sub4 <- subset(endo, subset = leiden_res0.4 == "4")
DimPlot(sub4, group.by = "leiden_res0.4_sub4_res0.4", label = TRUE) + coord_equal()
VlnPlot(sub4, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.4_sub4_res0.2", ncol = 2, pt.size = 0)
DotPlot(endo, features = c("PLVAP", "ACKR1", "POSTN", "COL15A1", # 1
                           "CCL21", "FABP4", "MRC1", "GNG11", # 0
                           "PTPRC", "TRAC", "CD3E", # 2
                           "IL7R", "FCN3", # 3
                           "CA4", "ITGA3", "KDR", # 5
                           "HEY1", "VEGFA", "HERPUD1", # 6
                           "APLN", "APLNR" # aCap/gCap
), group.by = "leiden_res0.4_sub4_res0.2") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
sub4 <- SetIdent(sub4, value = "leiden_res0.4_sub4_res0.2")
sub4_markers <- FindAllMarkers(sub4, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(sub4, features = c("PLVAP", "ACKR1", "POSTN", "COL15A1", # 1
                               "CCL21", "FABP4", "MRC1", "GNG11", # 0
                               "PTPRC", "TRAC", "CD3E", # 2
                               "IL7R", "FCN3", # 3
                               # 4 No markers for 4
                               "CA4", "ITGA3", "KDR", # 5
                               "HEY1", "VEGFA", # 6
                               "APLN", "APLNR" # aCap/gCap
), pt.size = 0.001) & coord_equal()
DotPlot(sub4, features = c("PLVAP", "ACKR1", "POSTN", "COL15A1", # 1
                           "CCL21", "FABP4", "MRC1", "GNG11", # 0
                           "PTPRC", "TRAC", "CD3E", # 2
                           "IL7R", "FCN3", # 3
                           "CA4", "ITGA3", "KDR", # 5
                           "HEY1", "VEGFA", "HERPUD1", # 6
                           "APLN", "APLNR" # aCap/gCap
), group.by = "leiden_res0.4_sub4_res0.2") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Create firstpass CT variable
endo$CT_firstpass <- ""
endo$CT_firstpass[endo$leiden_res0.4 == "0"] <- "Lymphatic"
endo$CT_firstpass[endo$leiden_res0.4 == "1" |
                    endo$leiden_res0.4_sub4_res0.2 == "4,0"] <- "Venous"
endo$CT_firstpass[endo$leiden_res0.4 %in% c("3", "5") |
                    endo$leiden_res0.4_sub4_res0.2 == "4,1"] <- "Capillary"
endo$CT_firstpass[endo$leiden_res0.4_sub3_res0.3 == "3,3"] <- "Fibroblasts (Endo Obj)"
endo$CT_firstpass[endo$leiden_res0.4 == "5"] <- "Arteriole"




endo$CT_firstpass <- ""
endo$CT_firstpass[endo$leiden_res0.4 == "0"] <- "Lymphatic"
endo$CT_firstpass[endo$leiden_res0.4 == "1" |
                    endo$leiden_res0.4_sub2_res0.4_sub3_res0.3 %in% c("2,0", "2,1") |
                    endo$leiden_res0.4_sub4_res0.2 == "4,0"] <- "Venous"
endo$CT_firstpass[endo$leiden_res0.4_sub2_res0.4_sub3_res0.3 == "2,3,2"] <- "T-cells (Endo Obj)"
endo$CT_firstpass[endo$leiden_res0.4_sub2_res0.4_sub3_res0.3 %in% c("2,2", "2,3,0", "2,3,1", "2,4", "2,5") |
                    endo$leiden_res0.4 %in% c("3", "5") |
                    endo$leiden_res0.4_sub4_res0.2 == "4,1"] <- "Capillary"
endo$CT_firstpass[endo$leiden_res0.4_sub3_res0.3 == "3,3"] <- "Fibroblasts (Endo Obj)"
endo$CT_firstpass[endo$leiden_res0.4_sub2_res0.4_sub3_res0.3 == "2,6" | 
                    endo$leiden_res0.4 == "6"] <- "Arteriole"
table(endo$CT_firstpass, endo$leiden_res0.4_sub2_res0.4_sub3_res0.3)
DimPlot(endo, group.by = "initial_CT", label = TRUE) + coord_equal()
DimPlot(endo, group.by = "CT_firstpass", label = TRUE) + coord_equal()

# Save object
#saveRDS(endo, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/endothelial_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_061324.rds")


#### EPITHELIAL ----
epi <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_lineage_split_nucnuc_RAPIDS_clustered_052824.rds")
DimPlot(epi, group.by = "leiden_res0.3", label = TRUE) + coord_equal()
DimPlot(epi, group.by = "initial_CT", label = TRUE, repel = TRUE) + coord_equal()

# Basic DimPlots and FeaturePlots
DimPlot(epi, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(epi, group.by = "tma") + coord_equal()
DimPlot(epi, group.by = "run", split.by = "run") + coord_equal()
DimPlot(epi, group.by = "sample_type", split.by = "sample_type") + coord_equal()
FeaturePlot(epi, features = c("nCount_RNA", "nFeature_RNA")) & coord_equal()
DimPlot(epi, group.by = "leiden_res0.3", split.by = "leiden_res0.3") + coord_equal()

# Create and view module scores
epi$endo_score <- colSums(epi@assays$RNA@counts[endothelial_features, ])/colSums(epi@assays$RNA@counts)
epi$epi_score <- colSums(epi@assays$RNA@counts[epithelial_features, ])/colSums(epi@assays$RNA@counts)
epi$imm_score <- colSums(epi@assays$RNA@counts[immune_features, ])/colSums(epi@assays$RNA@counts)
epi$mes_score <- colSums(epi@assays$RNA@counts[mesenchymal_features, ])/colSums(epi@assays$RNA@counts)
FeaturePlot(epi, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(epi, features = c("endo_score", "epi_score", "imm_score", "mes_score"), pt.size = 0, ncol = 2)

# Find Markers
epi <- SetIdent(epi, value = "leiden_res0.3")
epi_markers <- FindAllMarkers(epi, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(epi$adj_x_centroid, epi$adj_y_centroid)
row.names(position_xy) <- row.names(epi@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
epi[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                    assay = DefaultAssay(epi))

# Look at lineage and top markers
FeaturePlot(epi, features = c("EPCAM", "PECAM1", "PTPRC", "DCN", "LUM")) & coord_equal()
FeaturePlot(epi, features = c("SFTPC", "SFTPD", "AGER", "RTKN2", "KRT8", "CEACAM5", "CEACAM6",
                              "SCGB3A2", "SCGB1A1", "MUC5B", "C20orf85",
                              "KRT5", "KRT17", "MKI67", "TOP2A"
), ncol = 6) & coord_equal()
FeaturePlot(epi, features = c("RTKN2", "AGER", "COL4A3", # 0: AT1
                              "SCG2", "CHGB", "CALCA", "ASCL1", # 1: PNEC + ?
                              "KRT5", "KRT15", "KRT17", "TP63", # 2: Basal
                              "C20orf85", "FOXJ1", "SCGB1A1", "TP73", # 3: Ciliated
                              "GZMA", "CD3E", "TRAC", # 4: T-cells + ?
                              "MUC5B", "WFDC2", "MMP7", "SCGB3A2", # 5: MUC5B+ and SCGB+
                              "SFTPC", "SFTPD", "PGC", "NAPSA", "LAMP3", # 6: AT2 and Transitional AT2
                              "CPA3", "TPSAB1", "KIT", # 7: Mast
                              "TOP2A", "MKI67", "CDK1", # 8: Proliferating
                              "NKX2-1", "BMP4", "SPINK1" # 9: ?
)) & coord_equal()
DotPlot(epi, features = c("RTKN2", "AGER", "COL4A3", # 0: AT1
                          "SCG2", "CHGB", "CALCA", "ASCL1", # 1: PNEC + ?
                          "KRT5", "KRT15", "KRT17", "TP63", # 2: Basal
                          "C20orf85", "FOXJ1", "SCGB1A1", "TP73", # 3: Ciliated
                          "GZMA", "CD3E", "TRAC", # 4: T-cells + ?
                          "MUC5B", "WFDC2", "MMP7", "SCGB3A2", # 5: MUC5B+ and SCGB+
                          "SFTPC", "SFTPD", "PGC", "NAPSA", "LAMP3", # 6: AT2 and Transitional AT2
                          "CPA3", "TPSAB1", "KIT", # 7: Mast
                          "TOP2A", "MKI67", "CDK1", # 8: Proliferating
                          "NKX2-1", "BMP4", "SPINK1" # 9: ?
)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
FeaturePlot(epi, features = c("RTKN2", "AGER", "COL4A3")) & coord_equal() # 0: AT1
FeaturePlot(epi, features = c("SCG2", "CHGB", "CALCA", "SFTPC", "C20orf85", "MUC5B")) & coord_equal() # 1: PNEC + Other markers
FeaturePlot(epi, features = c("KRT5", "KRT15", "KRT17", "TP63")) & coord_equal() # 2: Basal

DotPlot(epi, features = rownames(epi)[1:115]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(epi, features = rownames(epi)[116:231]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(epi, features = rownames(epi)[232:343]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))


###### VIEW SUBCLUSTERING RESULTS ----
# Subcluster 1 to get PNECs split from other cell types
sub1 <- subset(epi, subset = leiden_res0.3 == "1")
DimPlot(epi, group.by = "leiden_res0.3_sub1_res0.13", label = TRUE) + coord_equal()
DimPlot(sub1, group.by = "leiden_res0.3_sub1_res0.13", label = TRUE) + coord_equal()
FeaturePlot(sub1, features = c("CALCA", "CHGB", "SCG2", "C20orf85", "SCGB1A1", "FOXJ1", "MUC5B", "NAPSA", "SFTPC",
                               "SFTPD", "SFTPD", "AGER", "COL4A3")) & coord_equal()
VlnPlot(sub1, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.3_sub1_res0.13", pt.size = 0, ncol = 2)
DotPlot(epi, group.by = "leiden_res0.3_sub1_res0.13", features = epithelial_features) +
  theme_angle + theme(axis.text.x = element_text(size = 8))
sub1 <- SetIdent(sub1, value = "leiden_res0.3_sub1_res0.13")
sub1_markers <- FindAllMarkers(sub1, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# 1,0: PNEC
# 1,1: Ciliated
# 1,2: AT2

# Subcluster 7 to possibly split out mast cells
sub7 <- subset(epi, subset = leiden_res0.3 == "7")
DimPlot(epi, group.by = "leiden_res0.3_sub7_res0.3", label = TRUE) + coord_equal()
DimPlot(sub7, group.by = "leiden_res0.3_sub7_res0.3", label = TRUE) + coord_equal()
VlnPlot(sub7, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.3_sub7_res0.3", pt.size = 0, ncol = 2)
DotPlot(sub7, group.by = "leiden_res0.3_sub7_res0.3", 
        features = c(epithelial_features, "CPA3", "TPSAB1", "KIT")) +
  theme_angle + theme(axis.text.x = element_text(size = 8))
DotPlot(sub7, group.by = "leiden_res0.3_sub7_res0.3", features = endothelial_features) + theme_angle
sub7 <- SetIdent(sub7, value = "leiden_res0.3_sub7_res0.3")
sub7_markers <- FindAllMarkers(sub7, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# 7,0: Ciliated
# 7,1 + 7,4: AT2
# 7,2 + 7,3: Secretory
# 7,5: Mast
# 7,6: Venous (Epithelial Obj)
# 7,7: Basal

# Subcluster 4 to examine low-quality cluster
sub4 <- subset(epi, subset = leiden_res0.3 == "4")
DimPlot(epi, group.by = "leiden_res0.3_sub4_res0.2", label = TRUE) + coord_equal()
DimPlot(sub4, group.by = "leiden_res0.3_sub4_res0.2", label = TRUE) + coord_equal()
FeaturePlot(sub4, features = c("RTKN2", "AGER", "COL4A3", "KRT5", "KRT15", "KRT17", "TP63",
                               "C20orf85", "FOXJ1", "SCGB1A1", "TP73", "MUC5B", "MMP7",
                               "SCGB3A2", "SFTPC", "SFTPD", "NAPSA", "LAMP3"), pt.size = 0.001) & coord_equal()
VlnPlot(sub4, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.3_sub4_res0.3", pt.size = 0, ncol = 2)
DotPlot(sub4, group.by = "leiden_res0.3_sub4_res0.2", features = c(epithelial_features, endothelial_features)) + theme_angle
sub4 <- SetIdent(sub4, value = "leiden_res0.3_sub4_res0.2")
sub4_markers <- FindAllMarkers(sub4, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# 4,0: Secretory
# 4,1 + 4,2: AT2
# 4,3: Differentiating Ciliated
# 4,4: Venous (Epi Obj)

# Create firstpass variable
epi$CT_firstpass <- ""
epi$CT_firstpass[epi$leiden_res0.3 == "0"] <- "AT1"
epi$CT_firstpass[epi$leiden_res0.3_sub1_res0.13 == "1,0"] <- "PNEC"
epi$CT_firstpass[epi$leiden_res0.3 == "2" |
                   epi$leiden_res0.3_sub7_res0.3 == "7,7"] <- "Basal"
epi$CT_firstpass[epi$leiden_res0.3_sub1_res0.13 == "1,1" |
                   epi$leiden_res0.3 == "3" |
                   epi$leiden_res0.3_sub7_res0.3 == "7,0"] <- "Ciliated"
epi$CT_firstpass[epi$leiden_res0.3_sub4_res0.2 == "4,0" |
                   epi$leiden_res0.3_sub7_res0.3 %in% c("7,2", "7,3")] <- "Secretory"
epi$CT_firstpass[epi$leiden_res0.3_sub1_res0.13 == "1,2" |
                   epi$leiden_res0.3_sub4_res0.2 %in% c("4,1", "4,2") |
                   epi$leiden_res0.3 == "6" |
                   epi$leiden_res0.3_sub7_res0.3 %in% c("7,1", "7,4")] <- "AT2"
epi$CT_firstpass[epi$leiden_res0.3_sub4_res0.2 == "4,3"] <- "Differentiating Ciliated"
epi$CT_firstpass[epi$leiden_res0.3_sub4_res0.2 == "4,4" |
                   epi$leiden_res0.3_sub7_res0.3 == "7,6"] <- "Venous (Epi Obj)"
epi$CT_firstpass[epi$leiden_res0.3 == "5"] <- "Mucous/Secretory"
epi$CT_firstpass[epi$leiden_res0.3_sub7_res0.3 == "7,5"] <- "Mast (Epi Obj)"
epi$CT_firstpass[epi$leiden_res0.3 == "8"] <- "Proliferating Epithelial"
epi$CT_firstpass[epi$leiden_res0.3 == "9"] <- "Transitional AT2"

# Split into airway and alveolar objects and save
air <- subset(epi, subset = CT_firstpass %in% c("PNEC", "Basal", "Ciliated", "Secretory", "Differentiating Ciliated", "Mucous/Secretory"))
alv <- subset(epi, subset = CT_firstpass %in% c("AT1", "AT2", "Proliferating Epithelial", "Transitional AT2"))
other <- subset(epi, subset = CT_firstpass %in% c("Mast (Epi Obj)", "Venous (Epi Obj)"))
ncol(air) + ncol(alv) + ncol(other) == ncol(epi) # Must be TRUE
# saveRDS(epi, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_052824.rds")
# saveRDS(other, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_052824_OTHER_LINEAGE.rds")
# saveRDS(air, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_airway_lineage_split_nucnuc_053024.rds")
# saveRDS(alv, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_alveolar_lineage_split_nucnuc_053024.rds")


###### AIRWAY ----
air <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_airway_lineage_split_nucnuc_RAPIDS_clustered_053024.rds")
DimPlot(air, group.by = "leiden_res0.4", label = TRUE) + coord_equal()
DimPlot(air, group.by = "initial_CT", label = TRUE, repel = TRUE) + coord_equal()

# Basic DimPlots and FeaturePlots
DimPlot(air, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(air, group.by = "tma") + coord_equal()
DimPlot(air, group.by = "run", split.by = "run") + coord_equal()
DimPlot(air, group.by = "sample_type", split.by = "sample_type") + coord_equal()
FeaturePlot(air, features = c("nCount_RNA", "nFeature_RNA")) & coord_equal()
DimPlot(air, group.by = "leiden_res0.4", split.by = "leiden_res0.4") + coord_equal()

# Create and view module scores
air$endo_score <- colSums(air@assays$RNA@counts[endothelial_features, ])/colSums(air@assays$RNA@counts)
air$epi_score <- colSums(air@assays$RNA@counts[epithelial_features, ])/colSums(air@assays$RNA@counts)
air$imm_score <- colSums(air@assays$RNA@counts[immune_features, ])/colSums(air@assays$RNA@counts)
air$mes_score <- colSums(air@assays$RNA@counts[mesenchymal_features, ])/colSums(air@assays$RNA@counts)
FeaturePlot(air, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(air, features = c("endo_score", "epi_score", "imm_score", "mes_score"), pt.size = 0, ncol = 2)

# Find Markers
air <- SetIdent(air, value = "leiden_res0.4")
air_markers <- FindAllMarkers(air, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(air$adj_x_centroid, air$adj_y_centroid)
row.names(position_xy) <- row.names(air@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
air[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                    assay = DefaultAssay(air))

# Look at lineage and top markers
FeaturePlot(air, features = c("EPCAM", "PECAM1", "PTPRC", "DCN", "LUM")) & coord_equal()
FeaturePlot(air, features = c("SFTPC", "SFTPD", "AGER", "RTKN2", "KRT8", "CEACAM5", "CEACAM6",
                              "SCGB3A2", "SCGB1A1", "MUC5B", "C20orf85", "KRT5", "KRT17", 
                              "MKI67", "TOP2A", "CALCA"
), ncol = 6) & coord_equal()
FeaturePlot(air, features = c("MS4A7", "HLA-DQA1", "LYZ", "PTPRC")) & coord_equal() # 0: Macrophages?
FeaturePlot(air, features = c("IFIT1", "IFIT3", "OAS2", "OAS3", "KRT5", "KRT17", "KRT15", "TP63")) & coord_equal() # 2: Interferon-high Basal; 4: Basal
FeaturePlot(air, features = c("C20orf85", "FOXJ1", "TP73", "MUC5B", "SCGB3A2", "SCGB1A1")) & coord_equal() # Secretory, Mucous, and Ciliated markers
DotPlot(air, features = c(epithelial_features, "CALCA")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Investigate specific clusters
# Investigate 1 and 8 (Mucous/Secretory)
sub1 <- subset(air, subset = leiden_res0.4 == "1")
FeaturePlot(sub1, features = c("SCGB3A2", "SCGB1A1"), blend = TRUE) & coord_equal()
sub1_8 <- subset(air, subset = leiden_res0.4 %in% c("1", "8"))
FeaturePlot(sub1_8, features = c("SCGB3A2", "SCGB1A1"), split.by = "leiden_res0.4") & coord_equal()
# Investigate cluster 0
sub0 <- subset(air, subset = leiden_res0.4 == "0")
FeaturePlot(sub0, features = c("MS4A7", "HLA-DQA1", "LYZ", "PTPRC", "RNASE1", "CEACAM6", "MMP7", "SCGB3A2",
                               "NAPSA", "C20orf85", "FOXJ1", "TP73", "KRT5", "KRT17", "KRT15", "TP63",
                               "SCGB1A1", "MUC5B"), pt.size = 0.001) & coord_equal() # 0: Macrophages? Immune contamination?


###### VIEW SUBCLUSTERING RESULTS ----
# Subcluster 1 to split SCGB3A2+/SCGB1A1+
VlnPlot(air, features = c("SCGB3A2", "SCGB1A1"), group.by = "leiden_res0.4_sub1_res0.3", pt.size = 0)
DimPlot(air, group.by = "leiden_res0.4_sub1_res0.3", label = TRUE) + coord_equal()
DimPlot(sub1, group.by = "leiden_res0.4_sub1_res0.3", label = TRUE) + coord_equal()
sub1 <- SetIdent(sub1, value = "leiden_res0.4_sub1_res0.3")
sub1_markers <- FindAllMarkers(sub1, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# 1,0: ? - Need to split further
# 1,1: MUC5B+
# 1,2: SCGB3A2+
# 1,3: SCGB3A2+/SCGB1A1+

# Further subcluster 1,0 to split out mucous/secretory
DimPlot(sub1, group.by = "leiden_res0.4_sub1_res0.3_sub0_res0.2", label = TRUE) + coord_equal()
FeaturePlot(sub1, features = c("SCGB3A2", "SCGB1A1", "MUC5B", "C20orf85", "KRT8"),
            split.by = "leiden_res0.4_sub1_res0.3_sub0_res0.2", pt.size = 0.01) & coord_equal()
DotPlot(air, features = c("SCGB3A2", "SCGB1A1", "MUC5B", "C20orf85", "KRT8"), 
        group.by = "leiden_res0.4_sub1_res0.3_sub0_res0.2") + theme_angle
# 1,0,1 + 1,2: SCGB3A2+
# 1,1: MUC5B+
# 1,0,0 + 1,3: SCGB3A2+/SCGB1A1+

# Subcluster 6 to get PNEC
DimPlot(air, group.by = "leiden_res0.4_sub6_res0.2", label = TRUE) + coord_equal()
DotPlot(air, features = c(epithelial_features, "CALCA", "CHGB"), group.by = "leiden_res0.4_sub6_res0.2") + theme_angle
sub6 <- subset(air, subset = leiden_res0.4 == "6")
DimPlot(sub6, group.by = "leiden_res0.4_sub6_res0.2", label = TRUE) + coord_equal()
sub6 <- SetIdent(sub6, value = "leiden_res0.4_sub6_res0.2")
sub6_markers <- FindAllMarkers(sub6, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# 6,0: PNEC
# 6,1: Differentiating Ciliated
# 6,2: Basal
# 6,3: Ciliated

# Subcluster 0 - immune contamination?
sub0 <- subset(air, subset = leiden_res0.4 == "0")
DimPlot(air, group.by = "leiden_res0.4_sub0_res0.2", label = TRUE) + coord_equal()
DimPlot(sub0, group.by = "leiden_res0.4_sub0_res0.15", label = TRUE) + coord_equal()
DimPlot(sub0, group.by = "leiden_res0.4_sub0_res0.2", label = TRUE) + coord_equal()
VlnPlot(sub0, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.4_sub0_res0.2", pt.size = 0, ncol = 2)
DotPlot(sub0, group.by = "leiden_res0.4_sub0_res0.2", features = c(epithelial_features, immune_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 6))
sub0 <- SetIdent(sub0, value = "leiden_res0.4_sub0_res0.2")
sub0_markers <- FindAllMarkers(sub0, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(sub0, features = c("KRT5", "TP63", "IFIT1", "IFIT2", "CD1A", "CD1C", "CCL22", "SFTPC", "NAPSA"), pt.size = 0.001) & coord_equal()
# 0,0: AT2
# 0,1: Interferon-high Basal
# 0,2 + 0,4: cDC2
# 0,3: Basal

# Subcluster 5 - this cluster has a lot of alveolar markers
DimPlot(air, group.by = "leiden_res0.4_sub5_res0.2") + coord_equal()
sub5 <- subset(air, subset = leiden_res0.4 == "5")
DimPlot(sub5, group.by = "leiden_res0.4_sub5_res0.2") + coord_equal()
VlnPlot(sub5, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.4_sub5_res0.2", pt.size = 0, ncol = 2)
DotPlot(sub5, group.by = "leiden_res0.4_sub5_res0.2", 
        features = c(epithelial_features, immune_features, "MKI67", "TOP2A", "FN1", "VIM")) +
  theme_angle + theme(axis.text.x = element_text(size = 5))
sub5 <- SetIdent(sub5, value = "leiden_res0.4_sub5_res0.2")
sub5_markers <- FindAllMarkers(sub5, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# # Further subcluster 5,3: Investigate immune contamination and possibly splitting epithelial further
FeaturePlot(sub5, split.by = "leiden_res0.4_sub5_res0.2", features = c("HLA-DRA", "SCGB3A2"), blend = TRUE) & coord_equal()
FeaturePlot(sub5, split.by = "leiden_res0.4_sub5_res0.2", features = c("SCGB3A2", "SCGB1A1", "MUC5B", "KRT5", "TP63", "C20orf85")) & coord_equal()
VlnPlot(sub5, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.4_sub5_res0.2", pt.size = 0, ncol = 4)
sub5_3 <- subset(sub5, subset = leiden_res0.4_sub5_res0.2 == "5,3")
DimPlot(air, group.by = "leiden_res0.4_sub5_res0.2_sub3_res0.3") + coord_equal()
DimPlot(sub5, group.by = "leiden_res0.4_sub5_res0.2_sub3_res0.3") + coord_equal()
DotPlot(sub5_3, group.by = "leiden_res0.4_sub5_res0.2_sub3_res0.3", 
        features = c(epithelial_features, immune_features, "MKI67", "TOP2A", "FN1", "VIM")) +
  theme_angle + theme(axis.text.x = element_text(size = 5))
# Can't split 5,3; look at 5,0
sub5_0 <- subset(sub5, subset = leiden_res0.4_sub5_res0.2 == "5,0")
DimPlot(sub5, group.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", label = TRUE) + coord_equal()
VlnPlot(sub5, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", pt.size = 0, ncol = 2)
DotPlot(sub5_0, group.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", 
        features = c(epithelial_features, immune_features, "MKI67", "TOP2A", "FN1", "VIM")) +
  theme_angle + theme(axis.text.x = element_text(size = 5))
DotPlot(sub5_0, group.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", features = epithelial_features) + theme_angle
DimPlot(sub5_0, group.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", label = TRUE) + coord_equal()
DimPlot(sub5_0, group.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", split.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4",
        ncol = 2, pt.size = 0.001) + coord_equal()
FeaturePlot(sub5_0, split.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", pt.size = 0.001,
            features = c("SCGB3A2", "SCGB1A1", "CD1C", "CD1A", "CCL22", "CD4")) & coord_equal()
# Cluster 5 is not airway-specific; investigate specific samples
DimPlot(air, reduction = "sp", group.by = "leiden_res0.4") + coord_equal()
DimPlot(air, reduction = "sp", group.by = "leiden_res0.4_sub5_res0.2", 
        cols = list(`2` = "black", `7` = "black", `8` = "black", `5,0` = "blue", `5,3` = "red")) + coord_equal()
vuild141_vuild142 <- subset(air, subset = sample %in% c("VUILD141", "VUILD142"))
DimPlot(vuild141_vuild142, reduction = "sp", group.by = "leiden_res0.4_sub5_res0.2", 
        cols = list(`2` = "black", `7` = "black", `8` = "black", `5,0` = "blue", `5,3` = "red")) + coord_equal()
DimPlot(vuild141_vuild142, reduction = "sp", group.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", 
        cols = list(`2` = "black", `7` = "black", `8` = "black", `5,0,0` = "blue", `5,0,2` = "red"), pt.size = 1) + coord_equal()
vuild141_vuild142_sub5 <- subset(air, subset = sample %in% c("VUILD141", "VUILD142") & leiden_res0.4_sub5_res0.2_sub0_res0.4 %in% c("5,0,0", "5,0,2"))
FeaturePlot(vuild141_vuild142_sub5, reduction = "sp", features = c("SFTPC", "RTKN2", "SCGB3A2", "SCGB1A1", "KRT8"),
            split.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", pt.size = 0.5) & coord_equal()
FeaturePlot(vuild141_vuild142_sub5, reduction = "sp", features = c("RTKN2", "SCGB3A2", "KRT8"), max.cutoff = 1, 
            split.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", pt.size = 0.5) & coord_equal()
FeaturePlot(vuild141_vuild142_sub5, reduction = "sp", features = c("SCGB3A2", "SCGB1A1"), blend = TRUE,
            split.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", pt.size = 0.5) & coord_equal()
FeaturePlot(vuild141_vuild142_sub5, reduction = "sp", features = c("RTKN2", "KRT8"), blend = TRUE,
            split.by = "leiden_res0.4_sub5_res0.2_sub0_res0.4", pt.size = 0.5) & coord_equal()
# 5,0,0: AT1 (Air Obj)
# 5,0,1: cDC1s
# 5,0,2: Transitional AT2 (Air Obj)
# 5,0,3 + 5,3: Differentiating Ciliated

# Look at all genes
DotPlot(air, features = rownames(air)[1:115]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(air, features = rownames(air)[116:231]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(air, features = rownames(air)[232:343]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# Create new celltype variable
air$CT_firstpassB <- ""
air$CT_firstpassB[air$leiden_res0.4_sub0_res0.2 == "0,0" |
                    air$leiden_res0.4_sub5_res0.2 == "5,1"] <- "AT2 (Air Obj)"
air$CT_firstpassB[air$leiden_res0.4_sub0_res0.2 %in% c("0,2", "0,4")] <- "cDC2 (Air Obj)"
air$CT_firstpassB[air$leiden_res0.4_sub0_res0.2 == "0,1" |
                    air$leiden_res0.4 == "2" ] <- "Interferon-high Basal"
air$CT_firstpassB[air$leiden_res0.4_sub0_res0.2 == "0,3" |
                    air$leiden_res0.4 == "4" |
                    air$leiden_res0.4_sub6_res0.2 == "6,2"] <- "Basal"
air$CT_firstpassB[air$leiden_res0.4_sub5_res0.2_sub0_res0.4 == "5,0,0"] <- "AT1 (Air Obj)"  # *****
air$CT_firstpassB[air$leiden_res0.4_sub5_res0.2_sub0_res0.4 == "5,0,2"] <- "Transitional AT2 (Air Obj)" # *****
air$CT_firstpassB[air$leiden_res0.4_sub1_res0.3_sub0_res0.2 %in% c("1,0,0", "1,3")] <- "SCGB3A2+/SCGB1A1+"
air$CT_firstpassB[air$leiden_res0.4_sub1_res0.3_sub0_res0.2 %in% c("1,0,1", "1,2")] <- "SCGB3A2+"
air$CT_firstpassB[air$leiden_res0.4_sub5_res0.2 == "5,2"] <- "Proliferating Fibroblasts (Air Obj)"
air$CT_firstpassB[air$leiden_res0.4_sub6_res0.2 == "6,0"] <- "PNEC"
air$CT_firstpassB[air$leiden_res0.4_sub1_res0.3_sub0_res0.2 == "1,1" |
                    air$leiden_res0.4 == "8"] <- "MUC5B+"
air$CT_firstpassB[air$leiden_res0.4_sub5_res0.2_sub0_res0.4 == "5,0,3" | 
                    air$leiden_res0.4_sub5_res0.2 == "5,3" |
                    air$leiden_res0.4_sub6_res0.2 == "6,1"] <- "Differentating Ciliated"
air$CT_firstpassB[air$leiden_res0.4 %in% c("3", "7") |
                    air$leiden_res0.4_sub6_res0.2 == "6,3"] <- "Ciliated"
air$CT_firstpassB[air$leiden_res0.4_sub5_res0.2_sub0_res0.4 == "5,0,1"] <- "cDC1 (Air Obj)"

# Plot
DimPlot(air, group.by = "CT_firstpassB", cols = distinctColorPalette(14), label = TRUE, repel = TRUE) + coord_equal()
VlnPlot(air, features = c("endo_score", "epi_score", "imm_score", "mes_score"), group.by = "CT_firstpassB", pt.size = 0, ncol = 2)
DotPlot(air, group.by = "CT_firstpassB", features = c(epithelial_features, "CALCA", "CD1C", "CCL22", "CD4", "VIM", "MKI67", "TOP2A")) + theme_angle
air_only_air <- subset(air, subset = CT_firstpassB %in% c("Basal", "Ciliated", "Differentiating Ciliated", "Interferon-high Basal",
                                                          "MUC5B+", "PNEC", "SCGB3A2+", "SCGB3A2+/SCGB1A1+"))
tmp_colors <- distinctColorPalette(8)
DimPlot(air_only_air, reduction = "sp", group.by = "CT_firstpassB", cols = tmp_colors, raster = FALSE) + coord_equal()

# Spot-check some samples
DimPlot(air, group.by = "sample", reduction = "sp", label = TRUE) + coord_equal() + NoLegend()
vuild96mf <- subset(air_only_air, subset = sample == "VUILD96MF")
DimPlot(vuild96mf, reduction = "sp", group.by = "CT_firstpassB", cols = c("red", "blue", "red", "green", "grey90", "purple", "orange")) + coord_equal()
vuild106 <- subset(air_only_air, subset = sample == "VUILD106")
DimPlot(vuild106, reduction = "sp", group.by = "CT_firstpassB", cols = c("red", "blue", "red", "green", "grey90", "purple", "orange")) + coord_equal()
vuild110 <- subset(air_only_air, subset = sample == "VUILD110")
DimPlot(vuild110, reduction = "sp", group.by = "CT_firstpassB", cols = c("red", "blue", "red", "green", "grey90", "purple", "orange")) + coord_equal()
FeaturePlot(vuild110, reduction = "sp", features = c("KRT5", "KRT17"), ncol = 1) & coord_equal()
vuild141 <- subset(air_only_air, subset = sample == "VUILD141")
DimPlot(vuild141, reduction = "sp", group.by = "CT_firstpassB", cols = c("red", "blue", "red", "green", "grey90", "purple", "orange")) + coord_equal()
FeaturePlot(vuild141, reduction = "sp", features = c("KRT5", "C20orf85", "SCGB3A2", "SCGB1A1", "MUC5B")) & coord_equal()
vuild78mf <- subset(air_only_air, subset = sample == "VUILD78MF")
DimPlot(vuild78mf, reduction = "sp", group.by = "CT_firstpassB", cols = c("red", "blue", "red", "green", "grey90", "purple", "orange")) + coord_equal()
FeaturePlot(vuild78mf, reduction = "sp", features = c("C20orf85", "SCGB3A2", "SCGB1A1"), ncol = 1) & coord_equal()
vuild91lf <- subset(air_only_air, subset = sample == "VUILD91LF")
DimPlot(vuild91lf, reduction = "sp", group.by = "CT_firstpassB", cols = c("red", "blue", "red", "green", "grey90", "purple", "orange")) + coord_equal()
FeaturePlot(vuild91lf, reduction = "sp", features = c("KRT5", "KRT17"), ncol = 1) & coord_equal()

# Save object
#saveRDS(air, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_airway_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS-B_LABELED_060424.rds")
air <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_airway_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS-B_LABELED_060424.rds")
table(air$CT_firstpassB)



# VUILD141, VUILD142 - AT1 (Air Obj) + Transitional AT2 (Air Obj)





###### ALVEOLAR ----
alv <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_alveolar_lineage_split_nucnuc_RAPIDS_clustered_060624.rds")
DimPlot(alv, group.by = "leiden_res0.5", label = TRUE) + coord_equal()
DimPlot(alv, group.by = "initial_CT", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(alv, group.by = "CT_firstpass", label = TRUE, repel = TRUE) + coord_equal()
FeaturePlot(alv, features = c("KRT5", "KRT17", "KRT8", "COL1A1")) & coord_equal() # Look for KRT5-/KRT17+ cells
table(alv$CT_firstpass, alv$leiden_res0.5)
# KRT5-/KRT17+ cells are in cluster 3 in leiden_res0.5

# Basic DimPlots and FeaturePlots
DimPlot(alv, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(alv, group.by = "tma") + coord_equal()
DimPlot(alv, group.by = "run", split.by = "run") + coord_equal()
DimPlot(alv, group.by = "sample_type", split.by = "sample_type") + coord_equal()
FeaturePlot(alv, features = c("nCount_RNA", "nFeature_RNA")) & coord_equal()
DimPlot(alv, group.by = "leiden_res0.5", split.by = "leiden_res0.5") + coord_equal()

# Create and view module scores
alv$endo_score <- colSums(alv@assays$RNA@counts[endothelial_features, ])/colSums(alv@assays$RNA@counts)
alv$epi_score <- colSums(alv@assays$RNA@counts[epithelial_features, ])/colSums(alv@assays$RNA@counts)
alv$imm_score <- colSums(alv@assays$RNA@counts[immune_features, ])/colSums(alv@assays$RNA@counts)
alv$mes_score <- colSums(alv@assays$RNA@counts[mesenchymal_features, ])/colSums(alv@assays$RNA@counts)
FeaturePlot(alv, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(alv, features = c("endo_score", "epi_score", "imm_score", "mes_score"), pt.size = 0, ncol = 2)

# Find Markers
alv <- SetIdent(alv, value = "leiden_res0.5")
alv_markers <- FindAllMarkers(alv, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(alv$adj_x_centroid, alv$adj_y_centroid)
row.names(position_xy) <- row.names(alv@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
alv[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                    assay = DefaultAssay(alv))

# Look at lineage and top markers
DotPlot(alv, features = c(epithelial_features, "MKI67", "TOP2A", "CDK1")) + theme_angle
FeaturePlot(alv, features = c("EPCAM", "PECAM1", "PTPRC", "DCN", "LUM")) & coord_equal()
FeaturePlot(alv, features = c("IFIT1", "IFIT3", "OAS2", "OAS3", # 0: Interferon markers
                              "LAMP3", "DUOX1", "SFTPC", "SFTPD", # 1: AT2
                              "MKI67", "TOP2A", "CDK1", "CENPF", # 2: Proliferating markers
                              "CEACAM6", "CEACAM5", "MMP7", "KRT8", # 3: Transitional markers
                              "SCGB3A2", # 4: Secretory marker
                              "LYZ", "CD68", "MS4A7", "FCER1G", # 5: Myeloid markers
                              "CALCA", "COL3A1", "COL1A2", # 6: Miscellaneous low expression
                              "AGER", "RTKN2", "COL4A3", "VEGFA" # 7: AT1
), ncol = 6) & coord_equal()
FeaturePlot(alv, features = epithelial_features, ncol = 8) & coord_equal()
FeaturePlot(alv, features = immune_features, ncol = 8) & coord_equal()
DotPlot(alv, features = c("IFIT1", "IFIT3", "OAS2", "OAS3", # 0: Interferon markers
                          "LAMP3", "DUOX1", "SFTPC", "SFTPD", # 1: AT2
                          "MKI67", "TOP2A", "CDK1", "CENPF", # 2: Proliferating markers
                          "CEACAM6", "CEACAM5", "MMP7", "KRT8", # 3: Transitional markers
                          "SCGB3A2", # 4: Secretory marker
                          "LYZ", "CD68", "MS4A7", "FCER1G", # 5: Myeloid markers
                          "CALCA", "COL3A1", "COL1A2", # 6: Miscellaneous low expression
                          "AGER", "RTKN2", "COL4A3", "VEGFA" # 7: AT1
)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
FeaturePlot(alv, features = c("AGER", "RTKN2", "COL4A3", "VEGFA", "LAMP3", 
                              "DUOX1", "SFTPC", "SFTPD", "SCGB3A2", "CEACAM6",
                              "MMP7", "KRT8", "SOX9", "SOX4")) & coord_equal() # AT1, AT2, and transitional markers
DotPlot(alv, features = c("AGER", "RTKN2", "COL4A3", "VEGFA", "LAMP3", "DUOX1", "SFTPC",
                          "SFTPD", "SCGB3A2", "CEACAM6", "MMP7", "KRT8", "SOX9", "SOX4")) + theme_angle # AT1, AT2, and transitional markers

DotPlot(alv, features = rownames(alv)[1:115]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(alv, features = rownames(alv)[116:231]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
DotPlot(alv, features = rownames(alv)[232:343]) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))


######## VIEW SUBCLUSTERING RESULTS ----
# Further investigate possible immune contamination in 5
sub5 <- subset(alv, subset = leiden_res0.5 == "5")
DimPlot(sub5, group.by = "leiden_res0.5_sub5_res0.3", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(sub5, group.by = "leiden_res0.5_sub5_res0.3", split.by = "leiden_res0.5_sub5_res0.3") + coord_equal()
FeaturePlot(sub5, features = c("LYZ", "CD68", "MS4A7", "FCER1G", "SFTPC", "LAMP3", "HLA-DRA", "PGC"))
FeaturePlot(sub5, features = c("HLA-DRA", "SFTPC"), blend = TRUE) & coord_equal()
DotPlot(sub5, group.by = "leiden_res0.5_sub5_res0.3", features = c(epithelial_features, immune_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 7))
VlnPlot(sub5, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "leiden_res0.5_sub5_res0.3",
        pt.size = 0, ncol = 2, same.y.lims = TRUE)
sub5 <- SetIdent(sub5, value = "leiden_res0.5_sub5_res0.3")
sub5_markers <- FindAllMarkers(sub5, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# Look at immune contamination with blend plots
FeaturePlot(sub5, features = c("LYZ", "SFTPC"), blend = TRUE) & coord_equal()
FeaturePlot(sub5, features = c("MARCO", "SFTPC"), blend = TRUE) & coord_equal()
FeaturePlot(sub5, features = c("HLA-DRA", "SFTPC"), blend = TRUE) & coord_equal()
FeaturePlot(sub5, features = c("XBP1", "SFTPC"), blend = TRUE) & coord_equal()
FeaturePlot(sub5, features = c("S100A8", "SFTPC"), blend = TRUE) & coord_equal()
FeaturePlot(sub5, features = c("SFTPC", "HLA-DRA"), blend = TRUE, split.by = "leiden_res0.5_sub5_res0.3", max.cutoff = 3, pt.size = 0.01) & coord_equal()
# Only HLA-DRA and LYZ are somewhat problematic
DimPlot(sub5, group.by = "leiden_res0.5_sub5_res0.3", reduction = "sp") + coord_equal()
sort(table(alv$sample, alv$leiden_res0.5)[, "5"]) # Highest in THD0008 (squished sample)
DotPlot(sub5, group.by = "leiden_res0.5_sub5_res0.3", features = c(epithelial_features, immune_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 7))
thd0008 <- subset(sub5, subset = sample == "THD0008")
FeaturePlot(thd0008, features = c("SFTPC", "SCGB3A2", "KRT8", "CEACAM6", "LYZ", "HLA-DRA"), reduction = "sp", keep.scale = "all") & coord_equal()
FeaturePlot(thd0008, features = c("SFTPC", "HLA-DRA"), reduction = "sp", blend = TRUE) & coord_equal()
FeaturePlot(thd0008, features = c("SFTPC", "HLA-DRA"), reduction = "sp", blend = TRUE, split.by = "leiden_res0.5_sub5_res0.3") & coord_equal()# 5,1 = AT2
# 5,0 = AT2 with immune contamination
VlnPlot(sub5, features = c("SFTPC", "SCGB3A2", "KRT8", "CEACAM6", "LYZ", "HLA-DRA"),
        group.by = "leiden_res0.5_sub5_res0.3",
        pt.size = 0, ncol = 2)
DotPlot(alv, group.by = "leiden_res0.5_sub5_res0.3", features = c(epithelial_features, immune_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 7))

# Further investigate possible mesenchymal contamination in 6
sub6 <- subset(alv, subset = leiden_res0.5 == "6")
FeaturePlot(sub6, features = c("COL3A1", "COL1A2", "SPARCL1", "SFTPC", "LAMP3"))
DimPlot(sub6, group.by = "leiden_res0.5_sub6_res0.3", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(sub6, group.by = "leiden_res0.5_sub6_res0.3", split.by = "leiden_res0.5_sub6_res0.3") + coord_equal()
DotPlot(sub6, group.by = "leiden_res0.5_sub6_res0.3", features = c(epithelial_features, mesenchymal_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 7))
DotPlot(sub6, group.by = "leiden_res0.5_sub6_res0.3", features = c(epithelial_features, endothelial_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 7))
DotPlot(alv, group.by = "leiden_res0.5_sub6_res0.3", features = c(epithelial_features, mesenchymal_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 7))
FeaturePlot(sub6, features = c("SFTPC", "FN1"), blend = TRUE, split.by = "leiden_res0.5_sub6_res0.3") & coord_equal()
sub6 <- SetIdent(sub6, value = "leiden_res0.5_sub6_res0.3")
sub6_markers <- FindAllMarkers(sub6, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DimPlot(sub6, reduction = "sp", group.by = "leiden_res0.5_sub6_res0.3", cols = distinctColorPalette(6))
vuild141 <- subset(epi, subset = sample == "VUILD141")
FeaturePlot(vuild141, reduction = "sp", features = c("SCGB3A2", "SCGB1A1", "MUC5B", "C20orf85", "SFTPC", "KRT5")) & coord_equal()
vuild141 <- subset(sub6, subset = sample == "VUILD141")
DimPlot(vuild141, reduction = "sp", group.by = "leiden_res0.5_sub6_res0.3", cols = distinctColorPalette(6)) + coord_equal()
FeaturePlot(vuild141, reduction = "sp", features = c("SCGB3A2", "SCGB1A1", "MUC5B", "C20orf85", "SFTPC", "KRT5")) & coord_equal()
air_sub6 <- merge(air, sub6)
# Add spatial coordinates as dimension reduction object
position_xy <- cbind(air_sub6$adj_x_centroid, air_sub6$adj_y_centroid)
row.names(position_xy) <- row.names(air_sub6@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
air_sub6[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                         assay = DefaultAssay(air_sub6))
air_sub6$leiden_res0.5_sub6_res0.3[is.na(air_sub6$leiden_res0.5_sub6_res0.3)] <- "NA"
air_sub6$leiden_res0.5_sub6_res0.3_sub2_res0.2[is.na(air_sub6$leiden_res0.5_sub6_res0.3_sub2_res0.2)] <- "NA"
table(air_sub6$leiden_res0.5_sub6_res0.3_sub2_res0.2, air_sub6$sample)
vuild58 <- subset(air_sub6, subset = sample == "VUILD58")
DimPlot(vuild58, group.by = "leiden_res0.5_sub6_res0.3", cols = c(distinctColorPalette(6), "grey90"), reduction = "sp", pt.size = 1) + coord_equal()
DimPlot(vuild58, group.by = "leiden_res0.5_sub6_res0.3_sub2_res0.2", cols = c(distinctColorPalette(8), "grey95"), reduction = "sp", pt.size = 1) + coord_equal()
DimPlot(vuild58, group.by = "leiden_res0.5_sub6_res0.3", cols = c("blue", "green", "red", "green", "blue", "blue", "grey90"), reduction = "sp") + coord_equal()
FeaturePlot(vuild58, features = c("SFTPC", "SCGB3A2", "MUC5B", "SCGB1A1"), reduction = "sp") & coord_equal()
vuild78mf <- subset(air_sub6, subset = sample == "VUILD78MF")
DimPlot(vuild78mf, group.by = "leiden_res0.5_sub6_res0.3", cols = c(distinctColorPalette(6), "grey90"), reduction = "sp") + coord_equal()
DimPlot(vuild78mf, group.by = "leiden_res0.5_sub6_res0.3_sub2_res0.2", cols = c(distinctColorPalette(8), "grey90"), reduction = "sp", pt.size = 1) + coord_equal()
FeaturePlot(vuild78mf, features = c("SFTPC", "SCGB3A2", "MUC5B", "SCGB1A1"), reduction = "sp") & coord_equal()
vuild141 <- subset(air_sub6, subset = sample == "VUILD141")
DimPlot(vuild141, group.by = "leiden_res0.5_sub6_res0.3", cols = c("red", distinctColorPalette(5), "grey90"), reduction = "sp") + coord_equal()
DimPlot(vuild141, group.by = "leiden_res0.5_sub6_res0.3_sub2_res0.2", cols = c(distinctColorPalette(8), "grey90"), reduction = "sp", pt.size = 1) + coord_equal()
FeaturePlot(vuild141, features = c("SFTPC", "SCGB3A2", "MUC5B", "SCGB1A1"), reduction = "sp") & coord_equal()
vuild142 <- subset(air_sub6, subset = sample == "VUILD142")
DimPlot(vuild142, group.by = "leiden_res0.5_sub6_res0.3", cols = c("red", distinctColorPalette(5), "grey90"), reduction = "sp") + coord_equal()
vuild106 <- subset(air_sub6, subset = sample == "VUILD106")
DimPlot(vuild106, group.by = "leiden_res0.5_sub6_res0.3", cols = c("red", distinctColorPalette(5), "grey90"), reduction = "sp") + coord_equal()
VlnPlot(sub6, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.5_sub6_res0.3", pt.size = 0, ncol = 2)
FeaturePlot(sub6, features = c("SCGB3A2", "MUC5B", "C20orf85", "SFTPC", "KRT5"),
            split.by = "leiden_res0.5_sub6_res0.3", pt.size = 0.001) & coord_equal()
FeaturePlot(sub6, features = c("SFTPC", "PECAM1"), blend = TRUE, 
            split.by = "leiden_res0.5_sub6_res0.3", pt.size = 0.001) & coord_equal()
vuild141 <- subset(air_sub6, subset = sample == "VUILD141")
DimPlot(vuild141, group.by = "leiden_res0.5_sub6_res0.3_sub0_res0.2", cols = c("red", "blue", rep("grey90", 6)), reduction = "sp")
vuild142 <- subset(air_sub6, subset = sample == "VUILD142")
DimPlot(vuild142, group.by = "leiden_res0.5_sub6_res0.3_sub0_res0.2", cols = c("red", "blue", rep("grey90", 6)), reduction = "sp")

# Further subcluster 6,0 + 6,1 + 6,2 + 6,5 and combine labels into 1 new variable
sub6$new_cluster6_var <- sub6$leiden_res0.5_sub6_res0.3
sub6$new_cluster6_var[sub6$leiden_res0.5_sub6_res0.3 == "6,0"] <- sub6$leiden_res0.5_sub6_res0.3_sub0_res0.2[sub6$leiden_res0.5_sub6_res0.3 == "6,0"]
sub6$new_cluster6_var[sub6$leiden_res0.5_sub6_res0.3 == "6,1"] <- sub6$leiden_res0.5_sub6_res0.3_sub1_res0.2[sub6$leiden_res0.5_sub6_res0.3 == "6,1"]
sub6$new_cluster6_var[sub6$leiden_res0.5_sub6_res0.3 == "6,2"] <- sub6$leiden_res0.5_sub6_res0.3_sub2_res0.2[sub6$leiden_res0.5_sub6_res0.3 == "6,2"]
sub6$new_cluster6_var[sub6$leiden_res0.5_sub6_res0.3 == "6,3"] <- sub6$leiden_res0.5_sub6_res0.3_sub3_res0.2[sub6$leiden_res0.5_sub6_res0.3 == "6,3"]
sub6$new_cluster6_var[sub6$leiden_res0.5_sub6_res0.3 == "6,4"] <- sub6$leiden_res0.5_sub6_res0.3_sub4_res0.3[sub6$leiden_res0.5_sub6_res0.3 == "6,4"]
sub6$new_cluster6_var[sub6$leiden_res0.5_sub6_res0.3 == "6,5"] <- sub6$leiden_res0.5_sub6_res0.3_sub5_res0.3[sub6$leiden_res0.5_sub6_res0.3 == "6,5"]
DimPlot(sub6, group.by = "new_cluster6_var", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(sub6, group.by = "new_cluster6_var", split.by = "new_cluster6_var", ncol = 5) & coord_equal()
sub6 <- SetIdent(sub6, value = "new_cluster6_var")
sub6_markers2 <- FindAllMarkers(sub6, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(sub6, features = epithelial_features, group.by = "new_cluster6_var") + theme_angle
DotPlot(sub6, features = endothelial_features, group.by = "new_cluster6_var") + theme_angle
DotPlot(sub6, features = immune_features, group.by = "new_cluster6_var") + theme_angle
DotPlot(sub6, features = mesenchymal_features, group.by = "new_cluster6_var") + theme_angle
DotPlot(sub6, features = c(epithelial_features, "EPAS1", "FCN3", "PECAM1", "S100A12", "S100A8", "S100A9", 
                           "ACTA2", "COL1A1", "COL1A2", "LUM"), group.by = "new_cluster6_var") + theme_angle
VlnPlot(sub6, group.by = "new_cluster6_var", features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        pt.size = 0, ncol = 2, same.y.lims = TRUE)
VlnPlot(sub6, group.by = "new_cluster6_var", features = c("SFTPC", "LAMP3", "AGER", "RTKN2"), pt.size = 0, ncol = 2)
VlnPlot(sub6, group.by = "leiden_res0.5_sub6_res0.3", features = c("SFTPC", "LAMP3", "AGER", "RTKN2"), pt.size = 0, ncol = 2)
vuild78mf <- subset(sub6, subset = sample == "VUILD78MF" & leiden_res0.5_sub6_res0.3 == "6,3" )
DimPlot(vuild78mf, group.by = "leiden_res0.5_sub6_res0.3_sub3_res0.2", cols = distinctColorPalette(5), reduction = "sp")
VlnPlot(sub6, features = c("AGER", "RTKN2", "SFTPC", "SFTPD", "LAMP3", "KRT5", "KRT17", "TP63", 
                           "SCGB3A2", "SCGB1A1", "KRT8", "COL1A1", "CEACAM5", "CEACAM6", "SOX9", "SOX4"), 
        group.by = "leiden_res0.5_sub6_res0.3", pt.size = 0)
DotPlot(sub6, features = c("AGER", "RTKN2", "SFTPC", "SFTPD", "LAMP3", "KRT5", "KRT17", "TP63", 
                           "SCGB3A2", "SCGB1A1", "KRT8", "COL1A1", "CEACAM5", "CEACAM6", "SOX9", "SOX4"), 
        group.by = "leiden_res0.5_sub6_res0.3") + theme_angle
sub1_4_6_7 <- subset(alv, subset = leiden_res0.5_sub6_res0.3 %in% c("1", "4", "6,0", "6,4", "6,5", "7")) # Compare to other AT1 and AT2 cell populations
DotPlot(sub1_4_6_7, features = c("AGER", "RTKN2", "SFTPC", "SFTPD", "LAMP3", "KRT5", "KRT17", "TP63", 
                                 "SCGB3A2", "SCGB1A1", "KRT8", "COL1A1", "CEACAM5", "CEACAM6", "SOX9", "SOX4"), 
        group.by = "leiden_res0.5_sub6_res0.3") + theme_angle
sub6_0 <- subset(alv, subset = leiden_res0.5_sub6_res0.3 == "6,0")
sort(rowMeans(sub6_0@assays$RNA@counts))

DimPlot(sub6, reduction = "sp", group.by = "leiden_res0.5_sub6_res0.3", 
        cols = c("blue", "red", "red", "red", "blue", "blue")) + coord_equal()
# 6,1 + 6,3 are airway cell types
# Look at AT2/possible transitional AT2 populations
# 6,2 could be airway or transitional alveolar
sort(table(sub6$leiden_res0.5_sub6_res0.3, sub6$sample)["6,0", ])
sort(table(sub6$leiden_res0.5_sub6_res0.3, sub6$sample)["6,4", ])
sort(table(sub6$leiden_res0.5_sub6_res0.3, sub6$sample)["6,5", ])
vuild91mf <- subset(sub6, subset = sample == "VUILD91MF")
DimPlot(vuild91mf, reduction = "sp", group.by = "leiden_res0.5_sub6_res0.3", 
        cols = c("blue", "grey90", "grey90", "grey90", "red", "green")) + coord_equal()
FeaturePlot(vuild91mf, features = c("SCGB3A2", "SFTPC", "KRT8", "AGER"), reduction = "sp") & coord_equal()
# Not many cells have direct co-expression of SFTPC and AGER
FeaturePlot(sub1_4_6_7, features = c("SFTPC", "AGER"), blend = TRUE, split.by = "leiden_res0.5") & coord_equal()
FeaturePlot(sub1_4_6_7, features = c("AGER", "SFTPC"), blend = TRUE, split.by = "leiden_res0.5") & coord_equal()
FeaturePlot(sub1_4_6_7, features = c("AGER", "SFTPC"), blend = TRUE, split.by = "leiden_res0.5", max.cutoff = 3) & coord_equal()
# Annotations:
# 6,0 + 6,4 + 6,5 = AT2
# 6,1 = Ciliated
# 6,2,0: SCGB3A2+
# 6,2,1: Transitional AT2
# 6,2,2: MUC5B+
# 6,3 = Basal

# Subcluster 7 to split AT1 and Transitional AT2, if necessary
sub7 <- subset(alv, subset = leiden_res0.5 == "7")
DimPlot(alv, group.by = "leiden_res0.5_sub7_res0.2", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(sub7, group.by = "leiden_res0.5_sub7_res0.15", label = TRUE, repel = TRUE) + coord_equal()
FeaturePlot(sub7, features = c("SFTPC", "RTKN2"), split.by = "leiden_res0.5_sub7_res0.15", blend = TRUE, pt.size = 0.001) & coord_equal()
FeaturePlot(sub7, features = c("SFTPC", "SFTPD", "LAMP3", "AGER", "RTKN2"), split.by = "leiden_res0.5_sub7_res0.15", pt.size = 0.001) & coord_equal()
DotPlot(sub7, features = epithelial_features, group.by = "leiden_res0.5_sub7_res0.4") + theme_angle
table(sub7$leiden_res0.5_sub7_res0.35, sub7$sample_type)
VlnPlot(sub7, group.by = "leiden_res0.5_sub7_res0.4", pt.size = 0, ncol = 2,
        features = c("SFTPC",  "SFTPD", "LAMP3", "AGER", "RTKN2", "CEACAM5", "CEACAM6", "KRT8", "SCGB3A2", "SOX9"))
table(sub7$leiden_res0.5_sub7_res0.4, sub7$sample_type)
# 7 is AT1

# Subcluster 3 to find KRT5-/KRT17+ cells
sub3 <- subset(alv, subset = leiden_res0.5 == "3")
DimPlot(alv, group.by = "leiden_res0.5_sub3_res0.4", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(sub3, group.by = "leiden_res0.5_sub3_res0.4", label = TRUE, repel = TRUE) + coord_equal()
vuild91lf <- subset(alv, subset = sample == "VUILD91LF" & leiden_res0.5 == "3")
DimPlot(vuild91lf, group.by = "leiden_res0.5_sub3_res0.3", reduction = "sp", cols = c(`3,0` = "grey90", `3,1` = "grey90", `3,2` = "red")) + coord_equal()
DimPlot(vuild91lf, group.by = "leiden_res0.5_sub3_res0.4", reduction = "sp", cols = c(`3,0` = "grey90", `3,1` = "grey90", `3,2` = "blue", `3,3` = "red")) + coord_equal()
FeaturePlot(vuild91lf, reduction = "sp", features = c("KRT5", "KRT17")) & coord_equal()
vuild91mf <- subset(alv, subset = sample == "VUILD91MF" & leiden_res0.5 == "3")
DimPlot(vuild91mf, group.by = "leiden_res0.5_sub3_res0.4", reduction = "sp") + coord_equal()
DotPlot(sub3, features = c(epithelial_features, "COL1A1"), group.by = "leiden_res0.5_sub3_res0.4") + theme_angle
VlnPlot(sub3, features = c("AGER", "RTKN2", "SFTPC", "SFTPD", "LAMP3", "KRT5", "KRT17", "TP63", 
                           "SCGB3A2", "SCGB1A1", "KRT8", "COL1A1", "CEACAM5", "CEACAM6", "SOX9", "SOX4"), 
        group.by = "leiden_res0.5_sub3_res0.4", pt.size = 0)
sub3 <- SetIdent(sub3, value = "leiden_res0.5_sub3_res0.4")
sub3_markers <- FindAllMarkers(sub3, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(sub3, features = c("AGER", "RTKN2", "COL4A3", "VEGFA", "LAMP3", "DUOX1", "SFTPC",
                           "SFTPD", "SCGB3A2", "CEACAM6", "MMP7", "KRT8", "SOX9", "SOX4")) + theme_angle # AT1, AT2, and transitional markers
DotPlot(alv, features = c("AGER", "RTKN2", "COL4A3", "VEGFA", "LAMP3", "DUOX1", "SFTPC",
                          "SFTPD", "SCGB3A2", "CEACAM6", "MMP7", "KRT8", "SOX9", "SOX4"),
        group.by = "leiden_res0.5_sub3_res0.4") + theme_angle # AT1, AT2, and transitional markers
FeaturePlot(sub6, features = c("AGER", "RTKN2", "COL4A3", "VEGFA", "LAMP3", "DUOX1", "SFTPC",
                               "SFTPD", "SCGB3A2", "CEACAM6", "MMP7", "KRT8", "SOX9", "SOX4")) + coord_equal() # AT1, AT2, and transitional markers
# 3,0 + 3,1 = Transitional AT2
# 3,2 = Basal
# 3,3 = KRT5-/KRT17+

# Subcluster 2 to split proliferating cell types
sub2 <- subset(alv, subset = leiden_res0.5 == "2")
DimPlot(alv, group.by = "leiden_res0.5_sub2_res0.5_sub9_res0.3", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(sub2, group.by = "leiden_res0.5_sub2_res0.5_sub9_res0.3", label = TRUE, repel = TRUE) + coord_equal()
DotPlot(sub2, group.by = "leiden_res0.5_sub2_res0.5_sub9_res0.3", features = c(epithelial_features, "MKI67", "TOP2A", "CDK1")) + theme_angle
sub2_9 <- subset(alv, subset = leiden_res0.5_sub2_res0.5 == "2,9")
FeaturePlot(sub2_9, features = c("SCGB3A2", "C20orf85", "KRT5", "MUC5B"),
            split.by = "leiden_res0.5_sub2_res0.5", pt.size = 0.01) & coord_equal() & patchwork::plot_layout(ncol = 4, nrow = 1)
FeaturePlot(sub2_9, features = c("SCGB3A2", "C20orf85", "KRT5", "MUC5B"),
            split.by = "leiden_res0.5_sub2_res0.5_sub9_res0.3", pt.size = 0.01) & coord_equal()
FeaturePlot(sub2_9, features = c("SCGB3A2", "C20orf85", "KRT5", "MUC5B"), pt.size = 0.01) & coord_equal()
FeaturePlot(sub2_9, features = c("SCGB3A2", "C20orf85"), blend = TRUE, pt.size = 0.01, max.cutoff = 1) & coord_equal()
# 2,0 = Proliferating AT2
# 2,1 + 2,2 + 2,4 + 2,5 + 2,6 + 2,7 + 2,9,0 = Proliferating Transitional AT2
# 2,3 + 2,9,2 = Proliferating Basal
# 2,8 = Proliferating SCGB3A2+/SCGB1A1+
# 2,9,1 = Proliferating Differentiating Ciliated

# Subcluster 0 to understand interferon-high cell types
sub0 <- subset(alv, subset = leiden_res0.5 == "0")
DimPlot(alv, group.by = "leiden_res0.5_sub0_res0.5", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(sub0, group.by = "leiden_res0.5_sub0_res0.5", label = TRUE, repel = TRUE) + coord_equal()
DotPlot(sub0, group.by = "leiden_res0.5_sub0_res0.5", features = c(epithelial_features, "CEACAM6", immune_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 7))
table(sub0$sample, sub0$leiden_res0.5_sub0_res0.5) # All clusters are nearly exclusive to VUILD110
DimPlot(alv, reduction = "sp", group.by = "leiden_res0.5_sub0_res0.5")
DimPlot(sub0, reduction = "sp", group.by = "leiden_res0.5_sub0_res0.5", cols = distinctColorPalette(8))
# Subclusters 0,0 + 0,3 + 0,4 + 0,6 + 0,7 are also present in VUILD48LF and VUILD48MF, which do not contain obvious airways
# From the DimPlot we can see that 0,1 + 0,2 appear to be airway cell types
# From the DotPlot it looks like 0,5 is airway as well
# 0,3 and 0,6 may be transitional alveolar or airway
FeaturePlot(sub0, features = c("SCGB3A2", "KRT5"), blend = TRUE)
sub0_0 <- subset(sub0, subset = leiden_res0.5_sub0_res0.5 == "0,0")
FeaturePlot(sub0_0, reduction = "sp", features = c("KRT5", "KRT17", "KRT8", "SPINK1",
                                                   "SFTPC", "SFTPD", "DUOX1")) & coord_equal()
DimPlot(sub0_0, reduction = "sp", group.by = "sample", label = TRUE) + coord_equal()
DimPlot(sub0_0, reduction = "sp", group.by = "sample", label = TRUE) + coord_equal()
FeaturePlot(sub0, split.by = "leiden_res0.5_sub0_res0.5", features = c("SCGB3A2", "C20orf85"), pt.size = 0.01) & coord_equal()
FeaturePlot(sub0, split.by = "leiden_res0.5_sub0_res0.5", features = c("SCGB3A2", "MUC5B"), pt.size = 0.01) & coord_equal()
FeaturePlot(sub0, split.by = "leiden_res0.5_sub0_res0.5", features = c("AGER", "SFTPC"), pt.size = 0.01) & coord_equal()
VlnPlot(sub0, group.by = "leiden_res0.5_sub0_res0.5", features = c("AGER", "SFTPC"), pt.size = 0)
sub0_1 <- subset(sub0, subset = leiden_res0.5_sub0_res0.5 == "0,1")
FeaturePlot(sub0_1, split.by = "leiden_res0.5_sub0_res0.5", features = c("SCGB3A2", "C20orf85"), blend = TRUE)
sub0_2 <- subset(sub0, subset = leiden_res0.5_sub0_res0.5 == "0,2")
FeaturePlot(sub0_2, features = c("SCGB3A2", "MUC5B"), blend = TRUE, pt.size = 0.01) & coord_equal()
sub0_7 <- subset(sub0, subset = leiden_res0.5_sub0_res0.5 == "0,7")
DimPlot(sub0_7, group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", label = TRUE) + coord_equal()
DimPlot(sub0_7, group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", split.by = "leiden_res0.5_sub0_res0.5_sub7_res0.5") + coord_equal()
FeaturePlot(sub0_7, split.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", features = c("AGER", "SFTPC")) & coord_equal()
FeaturePlot(sub0_7, split.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", features = c("AGER", "SFTPC"), blend = TRUE) & coord_equal()
DotPlot(sub0_7, group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", features = epithelial_features) + theme_angle
VlnPlot(sub0_7, group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", pt.size = 0,
        features = c("AGER", "RTKN2", "SFTPC", "SFTPD", "SCGB3A2", "KRT8", "CEACAM6"))
table(sub0_7$leiden_res0.5_sub0_res0.5_sub7_res0.3)
# 0,7 = Transitional AT2
DimPlot(sub0_7, group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", label = TRUE) + coord_equal()
DimPlot(sub0_7, group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", split.by = "leiden_res0.5_sub0_res0.5_sub7_res0.5") + coord_equal()
FeaturePlot(sub0_7, split.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", features = c("AGER", "SFTPC")) & coord_equal()
FeaturePlot(sub0_7, split.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", features = c("AGER", "SFTPC"), blend = TRUE) & coord_equal()
DotPlot(sub0_7, group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", features = epithelial_features) + theme_angle
VlnPlot(sub0_7, group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3", pt.size = 0,
        features = c("AGER", "RTKN2", "SFTPC", "SFTPD", "SCGB3A2", "KRT8", "CEACAM6"))
VlnPlot(sub0_7, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "leiden_res0.5_sub0_res0.5_sub7_res0.3",
        pt.size = 0, ncol = 4, same.y.lims = TRUE)
table(sub0_7$leiden_res0.5_sub0_res0.5_sub7_res0.3)
# 0,7 = Transitional AT2
sub0_2 <- subset(sub0, subset = leiden_res0.5_sub0_res0.5 == "0,2")
DimPlot(sub0_2, group.by = "leiden_res0.5_sub0_res0.5_sub2_res0.5", label = TRUE) + coord_equal()
DimPlot(sub0_2, group.by = "leiden_res0.5_sub0_res0.5_sub2_res0.5", split.by = "leiden_res0.5_sub0_res0.5_sub2_res0.5") + coord_equal()
FeaturePlot(sub0_2, split.by = "leiden_res0.5_sub0_res0.5_sub2_res0.5", features = c("MUC5B", "SCGB3A2"), blend = TRUE) & coord_equal()
FeaturePlot(sub0_2, split.by = "leiden_res0.5_sub0_res0.5_sub2_res0.5", features = c("MUC5B", "SCGB3A2"), blend = TRUE) & coord_equal()
VlnPlot(sub0_2, group.by = "leiden_res0.5_sub0_res0.5_sub2_res0.5", pt.size = 0, features = c("MUC5B", "SCGB3A2", "SCGB1A1", "C20orf85", "KRT5"))
FeaturePlot(sub0_2, features = c("MUC5B", "SCGB3A2", "SCGB1A1", "C20orf85", "KRT5")) & coord_equal()
DotPlot(sub0_2, features = c(epithelial_features, "CEACAM6"), group.by = "leiden_res0.5_sub0_res0.5_sub2_res0.5") + theme_angle
DotPlot(eqtl.ref, features = c(epithelial_features, "CEACAM6")) + theme_angle
VlnPlot(sub0_2, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "leiden_res0.5_sub0_res0.5_sub2_res0.5",
        pt.size = 0, ncol = 2, same.y.lims = TRUE)
# 0,0 = KRT5-/KRT17+
# 0,1 = Ciliated
# 0,2,0 + 0,2,2: MUC5B+
# 0,2,1 + 0,5 = Basal (light SCGB3A2 expression)
# 0,3 + 0,4 + 0,6 + 0,7 = Transitional AT2 (SCGB3A2, CEACAM5/6, and KRT8 expression, with some AGER)

# Create new celltype variable
alv$CT_firstpassB <- ""
alv$CT_firstpassB[alv$leiden_res0.5_sub0_res0.5_sub2_res0.5 %in% c("0,2,0", "0,2,2") | 
                    alv$leiden_res0.5_sub6_res0.3_sub2_res0.2 == "6,2,2"] <- "MUC5B+ (Alv Obj)"
alv$CT_firstpassB[alv$leiden_res0.5 %in% c("1", "4", "5") |
                    alv$leiden_res0.5_sub6_res0.3 %in% c("6,0", "6,4", "6,5")] <- "AT2"
alv$CT_firstpassB[alv$leiden_res0.5_sub2_res0.5 == "2,0"] <- "Proliferating AT2"
alv$CT_firstpassB[alv$leiden_res0.5_sub2_res0.5_sub9_res0.3 %in% c("2,1", "2,2", "2,4", "2,5",
                                                                   "2,6", "2,7", "2,9,0")] <- "Proliferating Transitional AT2"
alv$CT_firstpassB[alv$leiden_res0.5_sub2_res0.5_sub9_res0.3 %in% c("2,3", "2,9,2")] <- "Proliferating Basal (Alv Obj)"
alv$CT_firstpassB[alv$leiden_res0.5_sub2_res0.5 == "2,8"] <- "Proliferating SCGB3A2+/SCGB1A1+ (Alv Obj)"
alv$CT_firstpassB[alv$leiden_res0.5_sub2_res0.5_sub9_res0.3 == "2,9,1"] <- "Proliferating Differentiating Ciliated (Alv Obj)"
alv$CT_firstpassB[alv$leiden_res0.5_sub0_res0.5 %in% c("0,3", "0,4", "0,6", "0,7") |
                    alv$leiden_res0.5_sub3_res0.4 %in% c("3,0", "3,1") | 
                    alv$leiden_res0.5_sub6_res0.3_sub2_res0.2 == "6,2,1"] <- "Transitional AT2"
alv$CT_firstpassB[alv$leiden_res0.5_sub0_res0.5_sub2_res0.5 %in% c("0,2,1", "0,5") |
                    alv$leiden_res0.5_sub3_res0.4 == "3,2" |
                    alv$leiden_res0.5_sub6_res0.3 == "6,3"] <- "Basal (Alv Obj)"
alv$CT_firstpassB[alv$leiden_res0.5_sub0_res0.5 == "0,0" |
                    alv$leiden_res0.5_sub3_res0.4 == "3,3"] <- "KRT5-/KRT17+"
alv$CT_firstpassB[alv$leiden_res0.5_sub0_res0.5 == "0,1" |
                    alv$leiden_res0.5_sub6_res0.3 == "6,1"] <- "Ciliated (Alv Obj)"
alv$CT_firstpassB[alv$leiden_res0.5_sub6_res0.3_sub2_res0.2 == "6,2,0"] <- "SCGB3A2+ (Alv Obj)"
alv$CT_firstpassB[alv$leiden_res0.5 == "7"] <- "AT1"


table(alv$CT_firstpassB, alv$leiden_res0.5_sub6_res0.3_sub2_res0.2)
DimPlot(alv, group.by = "CT_firstpassB", label = TRUE, repel = TRUE, cols = distinctColorPalette(39)) + coord_equal()
DotPlot(alv, group.by = "CT_firstpassB", features = c(epithelial_features, "CEACAM6", "MKI67", "TOP2A", "CDK1")) + theme_angle

# Save object
#saveRDS(alv, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_alveolar_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS-B_LABELED_061324.rds")


#### IMMUNE ----
imm <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lineage_split_nucnuc_RAPIDS_clustered_052824.rds")
DimPlot(imm, group.by = "leiden_res0.4", label = TRUE, cols = distinctColorPalette(12)) + coord_equal()
DimPlot(imm, group.by = "initial_CT", label = TRUE, cols = distinctColorPalette(10)) + coord_equal()
table(imm$leiden_res0.4, imm$initial_CT)

# Basic DimPlots and FeaturePlots
DimPlot(imm, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(imm, group.by = "tma") + coord_equal()
DimPlot(imm, group.by = "run", split.by = "run") + coord_equal()
DimPlot(imm, group.by = "sample_type", split.by = "sample_type") + coord_equal()
FeaturePlot(imm, features = c("nCount_RNA", "nFeature_RNA")) & coord_equal()
DimPlot(imm, group.by = "leiden_res0.4", split.by = "leiden_res0.4", ncol = 5) + coord_equal()

# Create and view module scores
imm$endo_score <- colSums(imm@assays$RNA@counts[endothelial_features, ])/colSums(imm@assays$RNA@counts)
imm$epi_score <- colSums(imm@assays$RNA@counts[epithelial_features, ])/colSums(imm@assays$RNA@counts)
imm$imm_score <- colSums(imm@assays$RNA@counts[immune_features, ])/colSums(imm@assays$RNA@counts)
imm$mes_score <- colSums(imm@assays$RNA@counts[mesenchymal_features, ])/colSums(imm@assays$RNA@counts)
FeaturePlot(imm, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(imm, features = c("endo_score", "epi_score", "imm_score", "mes_score"), group.by = "leiden_res0.4", pt.size = 0, ncol = 4)

# General mesenchymal contamination throughout - check mesenchymal genes
FeaturePlot(imm, features = c(mesenchymal_features, "MSLN")) & coord_equal()
# Contamination comes from general mesenchymal lineage markers

# Find Markers
imm <- SetIdent(imm, value = "leiden_res0.4")
imm_markers <- FindAllMarkers(imm, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(imm$adj_x_centroid, imm$adj_y_centroid)
rownames(position_xy) <- rownames(imm@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
imm[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                    assay = DefaultAssay(imm))

# Look at lineage and specific cell type markers
FeaturePlot(imm, features = c("GNLY", "CD3E", "FOXP3", "CD4", "CD8A")) & coord_equal() # T-cells and NK cells
FeaturePlot(imm, features = c("FABP4", "SPP1", "PPARG", "MARCO", "FCGR3A", "CCL22", "CD1C")) & coord_equal() # Macrophages and DCs
FeaturePlot(imm, features = c("LILRA4", "CCL22", "CD1C", "CD1A", "ITGAM", "MRC1", "FCER1A")) & coord_equal() # cDCs, pDCs, and moDCs
FeaturePlot(imm, features = c("MS4A7", "LYZ", "FCER1G", "CD14", # 0: Macrophages
                              "CPA3", "TPSAB1", "KIT", # 1: Mast
                              "CD68", "MARCO", "PPARG", # 2: Macrophages
                              "SPP1", # 3: Macrophages; not SPP1+ specifically
                              "NKG7", "GNLY", "CD247", "GZMA", # 4: NK cells
                              "XBP1", "PIM2", "JCHAIN", # 5: Plasma
                              "S100A8", "S100A9", "SLC25A37", "S100A12", # 6: S100A+ Macrophages
                              "MKI67", "TOP2A", # 7: Proliferating mix
                              "FCER1A", "CD1C", "HLA-DQA1", # 8: cDCs
                              "TRAC", "CD3E", "IL7R", "CD4", "CD8A", # 9: T-cells
                              "MS4A1", "TNFRSF13C", "BANK1", # 10: B cells
                              "LAMP3", "CCL22", "CCR7" # 11: cDCs
), ncol = 6) & coord_equal()
DotPlot(imm, features = c("MS4A7", "LYZ", "FCER1G", "CD14", # 0: Macrophages
                          "CPA3", "TPSAB1", "KIT", # 1: Mast
                          "CD68", "MARCO", "PPARG", # 2: Macrophages
                          "SPP1", # 3: Macrophages; not SPP1+ specifically
                          "NKG7", "GNLY", "CD247", "GZMA", # 4: NK cells
                          "XBP1", "PIM2", "JCHAIN", # 5: Plasma
                          "S100A8", "S100A9", "SLC25A37", "S100A12", # 6: S100A+ Macrophages
                          "MKI67", "TOP2A", # 7: Proliferating ?
                          "FCER1A", "CD1C", "HLA-DQA1", # 8: cDC1s
                          "TRAC", "CD3E", "IL7R", "CD4", "CD8A", # 9: T-cells
                          "MS4A1", "TNFRSF13C", "BANK1", # 10: B cells
                          "LAMP3", "CCL22", "CCR7" # 11: cDC2s
)) + theme_angle
# Compare to eQTL celltypes
DotPlot(eqtl.ref, features = c("MS4A7", "LYZ", "FCER1G", "CD14", # 0: Macrophages
                               "CPA3", "TPSAB1", "KIT", # 1: Mast
                               "CD68", "MARCO", "PPARG", # 2: Macrophages
                               "SPP1", # 3: Macrophages; not SPP1+ specifically
                               "NKG7", "GNLY", "CD247", "GZMA", # 4: NK cells
                               "XBP1", "PIM2", "JCHAIN", # 5: Plasma
                               "S100A8", "S100A9", "SLC25A37", "S100A12", # 6: S100A+ Macrophages
                               "MKI67", "TOP2A", # 7: Proliferating ?
                               "FCER1A", "CD1C", "HLA-DQA1", # 8: cDC1s
                               "TRAC", "CD3E", "IL7R", "CD4", "CD8A", # 9: T-cells
                               "MS4A1", "TNFRSF13C", "BANK1", # 10: B cells
                               "LAMP3", "CCL22", "CCR7" # 11: cDC2s
), group.by = "manual_annotation_1") + theme_angle
DotPlot(eqtl.ref, features = immune_features, group.by = "manual_annotation_1") + theme_angle


###### VIEW SUBCLUSTERING RESULTS ----
# Subcluster 7 to split out different types fo proliferating cells
sub7 <- subset(imm, subset = leiden_res0.4 == "7")
DimPlot(imm, group.by = "leiden_res0.4_sub7_res0.2", label = TRUE) + coord_equal()
DimPlot(sub7, group.by = "leiden_res0.4_sub7_res0.2", label = TRUE) + coord_equal()
VlnPlot(sub7, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.4_sub7_res0.2", pt.size = 0, ncol = 2)
DotPlot(sub7, group.by = "leiden_res0.4_sub7_res0.2",
        features = c(immune_features, epithelial_features)) +
  theme_angle + theme(axis.text.x = element_text(size = 6))
sub7 <- SetIdent(sub7, value = "leiden_res0.4_sub7_res0.2")
sub7_markers <- FindAllMarkers(sub7, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(sub7, group.by = "leiden_res0.4_sub7_res0.2", features = epithelial_features) + theme_angle
FeaturePlot(sub7, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
DotPlot(imm, group.by = "leiden_res0.4_sub7_res0.2", features = c("MKI67", "TOP2A", "CDK1")) + theme_angle
# 7,0: Proliferating Epithelial
# 7,1 + 7,2: Proliferating Macrophages
# 7,3: pDCs
# 7,4: Proliferating T-cells

# Create new celltype variable
table(imm$initial_CT, imm$leiden_res0.4)
imm$CT_firstpass <- ""
imm$CT_firstpass[imm$leiden_res0.4 %in% c("0", "2", "3")] <- "Macrophages"
imm$CT_firstpass[imm$leiden_res0.4 == "1"] <- "Mast"
imm$CT_firstpass[imm$leiden_res0.4 == "4"] <- "NK cells"
imm$CT_firstpass[imm$leiden_res0.4 == "5"] <- "Plasma"
imm$CT_firstpass[imm$leiden_res0.4 == "6"] <- "S100A+ Macrophages"
imm$CT_firstpass[imm$leiden_res0.4_sub7_res0.2 == "7,0"] <- "Proliferating Epithelial"
imm$CT_firstpass[imm$leiden_res0.4_sub7_res0.2 %in% c("7,1", "7,2")] <- "Proliferating Macrophages"
imm$CT_firstpass[imm$leiden_res0.4_sub7_res0.2 == "7,3"] <- "pDCs"
imm$CT_firstpass[imm$leiden_res0.4_sub7_res0.2 == "7,4"] <- "Proliferating T-cells"
imm$CT_firstpass[imm$leiden_res0.4 == "8"] <- "cDC1s"
imm$CT_firstpass[imm$leiden_res0.4 == "9"] <- "T-cells"
imm$CT_firstpass[imm$leiden_res0.4 == "10"] <- "B cells"
imm$CT_firstpass[imm$leiden_res0.4 == "11"] <- "cDC2s"

# Plot
DimPlot(imm, group.by = "CT_firstpass", cols = distinctColorPalette(14), label = TRUE, repel = TRUE) + coord_equal()
VlnPlot(imm, features = c("endo_score", "epi_score", "imm_score", "mes_score"), group.by = "CT_firstpass", pt.size = 0, ncol = 2)
DotPlot(imm, group.by = "CT_firstpass", features = c(immune_features, "SFTPC", "EPCAM", "KRT8", "MKI67", "TOP2A")) + theme_angle

# These celltypes will be broken down further at the lymphoid/myeloid level
lym <- subset(imm, subset = CT_firstpass %in% c("NK cells", "Plasma", "pDCs", "Proliferating T-cells",
                                                "T-cells", "B cells"))
mye <- subset(imm, subset = CT_firstpass %in% c("Macrophages", "Mast", "S100A+ Macrophages", 
                                                "Proliferating Macrophages", "cDC1s", "cDC2s"))
other <- subset(imm, subset = CT_firstpass == "Proliferating Epithelial")
ncol(lym) + ncol(mye) + ncol(other) == ncol(imm)
#saveRDS(imm, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_052824.rds")
#saveRDS(other, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_052824_OTHER_LINEAGE.rds")
#saveRDS(lym, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lymphoid_lineage_split_nucnuc_060524.rds")
#saveRDS(mye, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_myeloid_lineage_split_nucnuc_060524.rds")


###### LYMPHOID ----
lym <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lymphoid_lineage_split_nucnuc_RAPIDS_clustered_060524.rds")
DimPlot(lym, group.by = "leiden_res0.3", label = TRUE, cols = distinctColorPalette(9)) + coord_equal()
DimPlot(lym, group.by = "initial_CT", label = TRUE, repel = TRUE, cols = distinctColorPalette(10)) + coord_equal()
DimPlot(lym, group.by = "CT_firstpass", label = TRUE, cols = distinctColorPalette(10)) + coord_equal()
table(lym$leiden_res0.3, lym$initial_CT)
table(lym$leiden_res0.3, lym$CT_firstpass)

# Basic DimPlots and FeaturePlots
DimPlot(lym, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(lym, group.by = "tma") + coord_equal()
DimPlot(lym, group.by = "run", split.by = "run") + coord_equal()
DimPlot(lym, group.by = "sample_type", split.by = "sample_type") + coord_equal()
FeaturePlot(lym, features = c("nCount_RNA", "nFeature_RNA")) & coord_equal()
DimPlot(lym, group.by = "leiden_res0.3", split.by = "leiden_res0.3", ncol = 5) + coord_equal()

# Create and view module scores
lym$endo_score <- colSums(lym@assays$RNA@counts[endothelial_features, ])/colSums(lym@assays$RNA@counts)
lym$epi_score <- colSums(lym@assays$RNA@counts[epithelial_features, ])/colSums(lym@assays$RNA@counts)
lym$imm_score <- colSums(lym@assays$RNA@counts[immune_features, ])/colSums(lym@assays$RNA@counts)
lym$mes_score <- colSums(lym@assays$RNA@counts[mesenchymal_features, ])/colSums(lym@assays$RNA@counts)
FeaturePlot(lym, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(lym, features = c("endo_score", "epi_score", "imm_score", "mes_score"), group.by = "leiden_res0.4", pt.size = 0, ncol = 4)

# General mesenchymal contamination throughout - check mesenchymal genes
FeaturePlot(lym, features = c(mesenchymal_features, "MSLN")) & coord_equal()
# Contamination comes from general mesenchymal lineage markers

# Find Markers
lym <- SetIdent(lym, value = "leiden_res0.3")
lym_markers <- FindAllMarkers(lym, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(lym$adj_x_centroid, lym$adj_y_centroid)
rownames(position_xy) <- rownames(lym@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
lym[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                    assay = DefaultAssay(lym))

# View cell type markers
FeaturePlot(lym, features = immune_features[1:40], ncol = 6) & coord_equal() # All immune features, part 1
FeaturePlot(lym, features = immune_features[41:81], ncol = 6) & coord_equal() # All immune features, part 2
FeaturePlot(lym, features = c("BANK1", "CD19", "CD79A", "LTB", "MS4A1", "TNFRSF13C", "HLA-DRA",
                              "CCR7", "CXCR4", "PTPRC", "CD69", "CD4", "CD8A", "CD2", "CD28", "CD3D",
                              "CD3E", "FOXP3", "GZMK", "TRAC", "ITM2C", "CD27", "CCL5", "LCK",
                              "FGFBP2", "GNLY", "NKG7", "HLA-DQA1", "KLRG1", "BCL2L11", "CD52",
                              "CTLA4", "PIM2", "IL7R", "LEF1", "CPA3", "KIT", "FCGR3A", "CD247", 
                              "GZMA", "JCHAIN"), ncol = 6) & coord_equal()
FeaturePlot(lym, features = c("LILRA4", "GZMB", "GPR183", "JCHAIN")) & coord_equal() # 0: pDCs
FeaturePlot(lym, features = c("TRAC", "IL7R", "CD3E", "CD2", "CD28", "CD3D")) & coord_equal() # 1: T-cells
FeaturePlot(lym, features = c("CD4", "CD8A", "CD8B", "FOXP3")) & coord_equal() # T-cell markers
FeaturePlot(lym, features = c("MKI67", "TOP2A", "CDK1")) & coord_equal() # 2: Proliferating markers
FeaturePlot(lym, features = c("LYZ", "CCL18", "MRC1", "CD14")) & coord_equal() # 3: Monocyte/macrophage markers
FeaturePlot(lym, features = c("GNLY", "NKG7", "EPAS1", "GZMB", "CD247", "GZMA")) & coord_equal() # 4: NK cells
FeaturePlot(lym, features = c("XBP1", "HERPUD1", "SEC11C", "PIM2", "JCHAIN",
                              "CPA3", "KIT", "TPSAB1"), ncol = 2) & coord_equal() # 5: Plasma; 7: Mast cells?
FeaturePlot(lym, features = c("GZMA", "GZMK", "CCL5", "CD8A", "CD3E")) & coord_equal() # 6: CD8+ T-cells
FeaturePlot(lym, features = c("MS4A1", "TNFRSF13C", "BANK1")) & coord_equal() # 8: B cells
DotPlot(lym, features = immune_features, group.by = "leiden_res0.3") + theme_angle

# Is cluster 3 CD4+ or CD8+?
# Remove all other T-cell clusters and CD4/CD8+ clusters and view dot plot
# Mix of both!
temp_lym <- subset(lym, subset = leiden_res0.3 %!in% c("0", "1", "2", "4", "6"))
DotPlot(temp_lym, features = c("CD4", "CD8A", "CD8B", "FOXP3", "GZMA", "GZMK", "CCL5", "CD3E"), group.by = "leiden_res0.3")


###### VIEW SUBCLUSTERING RESULTS ----
# Subcluster 1 to split CD4+ T-cells and Tregs
sub1 <- subset(lym, subset = leiden_res0.3 == "1")
DimPlot(lym, group.by = "leiden_res0.3_sub1_res0.2", label = TRUE) + coord_equal()
DimPlot(sub1, group.by = "leiden_res0.3_sub1_res0.2", label = TRUE) + coord_equal()
DimPlot(sub1, group.by = "leiden_res0.3_sub1_res0.2", split.by = "leiden_res0.3_sub1_res0.2") + coord_equal()
FeaturePlot(sub1, features = c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A", "CD8B")) & coord_equal()
FeaturePlot(sub1, features = c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A", "CD8B"), split.by = "leiden_res0.3_sub1_res0.2") & coord_equal()
DotPlot(sub1, features = c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A", "CD8B"), group.by = "leiden_res0.3_sub1_res0.2")
sub1 <- SetIdent(sub1, value = "leiden_res0.3_sub1_res0.2")
sub1_markers <- FindAllMarkers(sub1, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(lym, features = c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A", "CD8B"), group.by = "leiden_res0.3_sub1_res0.2")

# Subcluster 2 to split proliferating cell types
sub2 <- subset(lym, subset = leiden_res0.3 == "2")
DimPlot(lym, group.by = "leiden_res0.3_sub2_res0.2", label = TRUE) + coord_equal()
DimPlot(sub2, group.by = "leiden_res0.3_sub2_res0.2", label = TRUE) + coord_equal()
VlnPlot(sub2, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.3_sub2_res0.2", pt.size = 0, ncol = 4)
FeaturePlot(sub2, features = c("GZMA", "CD3E")) & coord_equal() # NK cells and T-cells are separate
FeaturePlot(sub2, features = c("GZMA", "GNLY", "CD3E", "MS4A1", "JCHAIN", "MKI67", "TOP2A", "CDK1", "CENPF")) & coord_equal()
sub2 <- SetIdent(sub2, value = "leiden_res0.3_sub2_res0.2")
sub2_markers <- FindAllMarkers(sub2, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(sub2, group.by = "leiden_res0.3_sub2_res0.2", features = c(immune_features, "MKI67", "TOP2A", "CDK1")) + theme_angle
VlnPlot(sub2, group.by = "leiden_res0.3_sub2_res0.2", pt.size = 0, 
        features = c("FABP4", "IL7R", "TRAC", "VIM", "MS4A1", "TNFRSF13C", "BANK1", "JCHAIN", "FCER1G", "CD247", "MKI67", "TOP2A"))

# Further subcluster 2,1 to split proliferating T-cells and NK cells
sub2 <- subset(lym, subset = leiden_res0.3 == "2")
DimPlot(sub2, group.by = "leiden_res0.3_sub2_res0.2_sub1_res0.4", label = TRUE) + coord_equal()
VlnPlot(sub2, group.by = "leiden_res0.3_sub2_res0.2_sub1_res0.4", pt.size = 0, 
        features = c("GZMA", "GNLY", "CD3E", "MKI67", "TOP2A"), ncol = 2)

# Subcluster 3 to split CD4+ and CD8+ T-cells
sub3 <- subset(lym, subset = leiden_res0.3 == "3")
FeaturePlot(sub3, features = c("CD3E", "MS4A1"))
DimPlot(lym, group.by = "leiden_res0.3_sub3_res0.4", label = TRUE) + coord_equal()
DimPlot(sub3, group.by = "leiden_res0.3_sub3_res0.4", label = TRUE) + coord_equal()
DimPlot(sub3, group.by = "leiden_res0.3_sub3_res0.4", split.by = "leiden_res0.3_sub3_res0.4", ncol = 3) + coord_equal()
FeaturePlot(sub3, features =  c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A", "CD8B"), pt.size = 0.01) & coord_equal()
FeaturePlot(sub3, features =  c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A", "CD8B"), pt.size = 0.01, max.cutoff = 1) & coord_equal()
FeaturePlot(sub3, features =  c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A"), 
            split.by = "leiden_res0.3_sub3_res0.2", pt.size = 0.001) & coord_equal()
FeaturePlot(sub3, features =  c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A"), 
            split.by = "leiden_res0.3_sub3_res0.4", pt.size = 0.001, max.cutoff = 1) & coord_equal()
DotPlot(sub3, features =  c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A"), group.by = "leiden_res0.3_sub3_res0.4")
DotPlot(lym, features =  c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A"), group.by = "leiden_res0.3_sub3_res0.4")
sub3 <- SetIdent(sub3, value = "leiden_res0.3_sub3_res0.4")
sub3_markers <- FindAllMarkers(sub3, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
DotPlot(lym, features = c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A", "HLA-DRA", "LYZ"), 
        group.by = "leiden_res0.3_sub3_res0.4") # 3,4 markers: LYZ, HLA-DRA
DotPlot(sub3, features = c("FOXP3", "CTLA4", "CD3E", "CD4", "CD8A", "HLA-DRA", "LYZ"), 
        group.by = "leiden_res0.3_sub3_res0.4") # 3,4 markers: LYZ, HLA-DRA
FeaturePlot(sub3, features ="LYZ", split.by = "leiden_res0.3_sub3_res0.4", pt.size = 0.001, max.cutoff = 1)
DotPlot(sub3, features = immune_features, group.by = "leiden_res0.3_sub3_res0.4") + theme_angle

# Check spatial organization of 3,4 cluster
DimPlot(sub3, reduction = "sp", group.by = "leiden_res0.3_sub3_res0.4", cols = list(`3,1` = "red", `3,4` = "blue")) + coord_equal()
vuild96mf <- subset(imm, subset = sample == "VUILD96MF")
FeaturePlot(vuild96mf, features = c("LYZ", "MARCO", "BANK1", "CD3E", "CD4", "CD8A"), reduction = "sp") & coord_equal()
vuild96mf_sub3 <- subset(sub3, subset = sample == "VUILD96MF")
DimPlot(vuild96mf_sub3, reduction = "sp", group.by = "leiden_res0.3_sub3_res0.4", cols = list(`3,1` = "red", `3,4` = "blue")) + coord_equal()
FeaturePlot(vuild96mf_sub3, features = c("LYZ", "MARCO", "BANK1", "CD3E", "CD4", "CD8A"), reduction = "sp", pt.size = 0.01) & coord_equal()
FeaturePlot(vuild96mf_sub3, features = c("LYZ", "MARCO", "BANK1", "CD3E", "CD4", "CD8A"), reduction = "sp", pt.size = 0.01, max.cutoff = 1) & coord_equal()
# 3,4: Granulomatous T-cells/macrophages

# Further subcluster 3,3 to get B cells
sub3 <- subset(lym, subset = leiden_res0.3 == "3")
sub3_3 <- subset(lym, subset = leiden_res0.3_sub3_res0.4 == "3,3")
DimPlot(sub3_3, group.by = "leiden_res0.3_sub3_res0.4_sub3_res0.3", label = TRUE) + coord_equal()
FeaturePlot(sub3_3, features = c("MS4A1", "CD3E", "JCHAIN"), split.by = "leiden_res0.3_sub3_res0.4_sub3_res0.3", pt.size = 0.001) & coord_equal()
sub3_3 <- SetIdent(sub3_3, value = "leiden_res0.3_sub3_res0.4_sub3_res0.3")
sub3_3_markers <- FindAllMarkers(sub3_3, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
VlnPlot(sub3_3, features = c("MS4A1", "CD3E", "JCHAIN", "CD4", "CD8A"), group.by = "leiden_res0.3_sub3_res0.4_sub3_res0.3", pt.size = 0)
DotPlot(sub3, features = c("MS4A1", "CD3E", "JCHAIN", "CD4", "CD8A"), group.by = "leiden_res0.3_sub3_res0.4_sub3_res0.3")
# 3,3,0 + 3,3,1 = T-cells
# 3,3,2 = Plasma
# 3,3,3: B cells

# Subcluster 7 to split B cells and plasma
sub7 <- subset(lym, subset = leiden_res0.3 == "7")
DimPlot(lym, group.by = "leiden_res0.3_sub7_res0.4", label = TRUE) + coord_equal()
DimPlot(sub7, group.by = "leiden_res0.3_sub7_res0.4", label = TRUE) + coord_equal()
VlnPlot(sub7, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.3_sub7_res0.4", pt.size = 0, ncol = 4)
DotPlot(sub7, group.by = "leiden_res0.3_sub7_res0.4",
        features = c("XBP1", "HERPUD1", "SEC11C", "PIM2", "JCHAIN", "CPA3", "KIT", "TPSAB1",
                     "CD3E", "CD4", "CD8A", "FOXP3")) + theme_angle
FeaturePlot(sub7, features = c("XBP1", "HERPUD1", "SEC11C", "PIM2", "JCHAIN", "CPA3", "KIT", "TPSAB1",
                               "CD3E", "CD4", "CD8A", "FOXP3"), pt.size = 0.01) & coord_equal()
sub7 <- SetIdent(sub7, value = "leiden_res0.3_sub7_res0.4")
sub7_markers <- FindAllMarkers(sub7, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(sub7, features = c("CPA3", "KIT", "TPSAB1", "CD3E", "CD4", "CD8A", "FOXP3"), 
            pt.size = 0.01, split.by = "leiden_res0.3_sub7_res0.4") & coord_equal()
FeaturePlot(sub7, features = c("CPA3", "CD3E"), group.by = "leiden_res0.3_sub7_res0.4", blend = TRUE) & coord_equal()

# Create celltype variable
lym$CT_firstpassB <- ""
lym$CT_firstpassB[lym$leiden_res0.3 == "0"] <- "pDCs"
lym$CT_firstpassB[lym$leiden_res0.3_sub1_res0.2 %in% c("1,0", "1,1") |
                    lym$leiden_res0.3_sub2_res0.2 == "2,0" |
                    lym$leiden_res0.3_sub3_res0.4 %in% c("3,0", "3,2") | 
                    lym$leiden_res0.3_sub3_res0.4_sub3_res0.3 %in% c("3,3,0", "3,3,1") |
                    lym$leiden_res0.3_sub7_res0.4 == "7,0"] <- "T-cells"
lym$CT_firstpassB[lym$leiden_res0.3_sub1_res0.2 == "1,2" |
                    lym$leiden_res0.3_sub7_res0.4 == "7,1"] <- "Tregs"
lym$CT_firstpassB[lym$leiden_res0.3_sub1_res0.2 == "1,3"] <- "CD4+ T-cells"
lym$CT_firstpassB[lym$leiden_res0.3_sub2_res0.2_sub1_res0.4 == "2,1,0"] <- "Proliferating NK cells"
lym$CT_firstpassB[lym$leiden_res0.3_sub2_res0.2_sub1_res0.4 %in% c("2,1,1", "2,1,2", "2,1,3", "2,1,4")] <- "Proliferating T-cells"
lym$CT_firstpassB[lym$leiden_res0.3_sub2_res0.2 == "2,2"] <- "Proliferating B cells"
lym$CT_firstpassB[lym$leiden_res0.3_sub3_res0.4 == "3,1" | 
                    lym$leiden_res0.3 == "6"] <- "CD8+ T-cells"
lym$CT_firstpassB[lym$leiden_res0.3_sub3_res0.4 == "3,4"] <- "Granulomatous T-cells"
lym$CT_firstpassB[lym$leiden_res0.3 == "4"] <- "NK cells"
lym$CT_firstpassB[lym$leiden_res0.3_sub3_res0.4_sub3_res0.3 == "3,3,2" |
                    lym$leiden_res0.3 == "5" | 
                    lym$leiden_res0.3_sub7_res0.4 %in% c("7,2", "7,3", "7,5")] <- "Plasma"
lym$CT_firstpassB[lym$leiden_res0.3_sub7_res0.4 == "7,4"] <- "Mast (Lymphoid Obj)"
lym$CT_firstpassB[lym$leiden_res0.3_sub3_res0.4_sub3_res0.3 == "3,3,3" | 
                    lym$leiden_res0.3 == "8"] <- "B cells"

# Plot
#tmp_colors <- distinctColorPalette(13)
DimPlot(lym, group.by = "initial_CT", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(lym, group.by = "CT_firstpass", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(lym, group.by = "CT_firstpassB", cols = tmp_colors, label = TRUE, repel = TRUE) + coord_equal()
DimPlot(lym, group.by = "CT_firstpassB", reduction = "sp", cols = tmp_colors, raster = FALSE, pt.size = 0.001) + coord_equal()
DotPlot(lym, group.by = "CT_firstpassB", features = c(immune_features, "MKI67", "TOP2A", "CDK1")) + theme_angle + theme(axis.text.x = element_text(size = 9))

# Spot-check samples
DimPlot(lym, group.by = "sample", reduction = "sp", label = TRUE) + coord_equal() + NoLegend()
vuild96mf <- subset(lym, subset = sample == "VUILD96MF")
DimPlot(vuild96mf, reduction = "sp", group.by = "CT_firstpassB", cols = tmp_colors) + coord_equal()
vuild106 <- subset(lym, subset = sample == "VUILD106")
DimPlot(vuild106, reduction = "sp", group.by = "CT_firstpassB", cols = tmp_colors) + coord_equal()
vuild110_thd0008 <- subset(lym, subset = sample %in% c("VUILD110", "THD0008"))
DimPlot(vuild110_thd0008, reduction = "sp", group.by = "CT_firstpassB", cols = tmp_colors) + coord_equal()
FeaturePlot(vuild110_thd0008, reduction = "sp", features = "GNLY") + coord_equal()
vuild141 <- subset(lym, subset = sample == "VUILD141")
DimPlot(vuild141, reduction = "sp", group.by = "CT_firstpassB", cols = tmp_colors) + coord_equal()
FeaturePlot(vuild141, reduction = "sp", features = c("BANK1", "JCHAIN"), ncol = 1) & coord_equal()

# Save object
#saveRDS(lym, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lymphoid_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS-B_LABELED_060524.rds")


###### MYELOID ----
mye <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_myeloid_lineage_split_nucnuc_RAPIDS_clustered_060724.rds")
DimPlot(mye, group.by = "leiden_res0.8", label = TRUE) + coord_equal()
DimPlot(mye, group.by = "leiden_res0.8", split.by = "leiden_res0.8", ncol = 5) & coord_equal()
DimPlot(mye, group.by = "initial_CT", label = TRUE) + coord_equal()
DimPlot(mye, group.by = "CT_firstpass", label = TRUE) + coord_equal()
table(mye$CT_firstpass, mye$leiden_res0.8)
table(mye$leiden_res0.8, mye$sample)

FeaturePlot(mye, features = c("CCL22", "FCER1A", "FCER1G", "CD2", "CD1C", "CD1A"))
DimPlot(mye, group.by = "leiden_res0.8", reduction = "sp", cols = distinctColorPalette(712)) + coord_equal()
DimPlot(mye, group.by = "leiden_res0.8", reduction = "sp", cols = list(`12`= "red")) + coord_equal() # Interesting sample: TILD117MFB
tild117mfb <- subset(mye, subset = sample == "TILD117MF" & tma == "TMA5") 
DimPlot(tild117mfb, group.by = "leiden_res0.8", reduction = "sp", cols = list(`12`= "red")) + coord_equal()
FeaturePlot(tild117mfb, features = c("MSLN", "HAS1", "FCN1", "S100A8"), reduction = "sp")

# Basic DimPlots and FeaturePlots
DimPlot(mye, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(mye, group.by = "tma") + coord_equal()
DimPlot(mye, group.by = "run", split.by = "run") + coord_equal()
DimPlot(mye, group.by = "sample_type", split.by = "sample_type") + coord_equal()
FeaturePlot(mye, features = c("nCount_RNA", "nFeature_RNA")) & coord_equal()
DimPlot(mye, group.by = "leiden_res0.3", split.by = "leiden_res0.3") + coord_equal()

# Create and view module scores
mye$endo_score <- colSums(mye@assays$RNA@counts[endothelial_features, ])/colSums(mye@assays$RNA@counts)
mye$epi_score <- colSums(mye@assays$RNA@counts[epithelial_features, ])/colSums(mye@assays$RNA@counts)
mye$imm_score <- colSums(mye@assays$RNA@counts[immune_features, ])/colSums(mye@assays$RNA@counts)
mye$mes_score <- colSums(mye@assays$RNA@counts[mesenchymal_features, ])/colSums(mye@assays$RNA@counts)
FeaturePlot(mye, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(mye, features = c("endo_score", "epi_score", "imm_score", "mes_score"), group.by = "leiden_res0.8", pt.size = 0, ncol = 4)

# Find Markers
mye <- SetIdent(mye, value = "leiden_res0.8")
mye_markers <- FindAllMarkers(mye, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(mye$adj_x_centroid, mye$adj_y_centroid)
row.names(position_xy) <- row.names(mye@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
mye[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                    assay = DefaultAssay(mye))

# Compare markers in eQTL dataset
eqtl.imm <- subset(eqtl.ref, subset = lineage == "Immune")
DotPlot(eqtl.imm, features = c(immune_features, "CD44", "ELANE", "MMP7", "LAMP3", "FCER1A", "RNASE1"), group.by = "manual_annotation_1") + theme_angle
eqtl.mye <- subset(eqtl.ref, subset = manual_annotation_1 %in% 
                     c("Monocyte-derived macrophage", "Monocyte", "moDC", "Mast", "Interstitial macrophage",
                       "cDC2", "cDC1", "Alveolar macrophage", "Inflammatory monocyte"))
DotPlot(eqtl.mye, features = c(immune_features, "CD44", "ELANE", "MMP7", "LAMP3", "FCER1A", "RNASE1"), group.by = "manual_annotation_1") + theme_angle
eqtl.mye <- SetIdent(eqtl.mye, value = "manual_annotation_1")
eqtl_mye_markers <- FindAllMarkers(eqtl.mye, features = rownames(mye)[-which(rownames(mye) == "CCN2")], 
                                   only.pos = TRUE) %>%
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
eqtl.mye2 <- subset(eqtl.ref, subset = manual_annotation_1 %in% 
                      c("Monocyte-derived macrophage", "Monocyte", "Interstitial macrophage", "Alveolar macrophage", "Inflammatory monocyte"))
DotPlot(eqtl.mye2, features = c(immune_features, "CD44", "ELANE", "MMP7", "LAMP3", "FCER1A", "RNASE1"), group.by = "manual_annotation_1") + theme_angle
eqtl.mye2 <- SetIdent(eqtl.mye2, value = "manual_annotation_1")
eqtl_mye2_markers <- FindAllMarkers(eqtl.mye2, features = rownames(mye)[-which(rownames(mye) == "CCN2")], 
                                    only.pos = TRUE) %>%
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Look at lineage and specific cell type markers
FeaturePlot(mye, features = c("EPCAM", "PECAM1", "PTPRC", "DCN", "LUM")) & coord_equal() # Lineage markers
DotPlot(mye, group.by = "leiden_res0.8", features = c(immune_features, "CD44", "ELANE", "MMP7", "LAMP3", "FCER1A", "RNASE1")) + theme_angle
FeaturePlot(mye, features = immune_features[1:40], ncol = 8) & coord_equal()
FeaturePlot(mye, features = immune_features[41:83], ncol = 9) & coord_equal()
FeaturePlot(mye, features = c("IFIT1", "IFIT3", "OAS2", "OAS3")) & coord_equal() # 0: Interferon response markers
FeaturePlot(mye, features = c("PIM2", "JCHAIN", "CD1A")) & coord_equal() # 1: pDC markers
FeaturePlot(mye, features = c("FABP4", "SPP1", "S100A8", "S100A9", "S100A12", "FCN1", "MARCO", "LYZ")) & coord_equal() # Macrophage population markers
FeaturePlot(mye, features = c("FABP4", "SPP1", "S100A8", "FCN1"), split.by = "leiden_res0.8") & coord_equal() # Macrophage population markers

# Look at spatial plots + specific samples to assess cell types
table(mye$leiden_res0.8, mye$sample)["2", ] %>% sort()
c3_13 <- subset(mye, subset = leiden_res0.8 %in% c("3", "13"))
DimPlot(c3_13, reduction = "sp", group.by = "leiden_res0.8") + coord_equal()
FeaturePlot(c3_13, reduction = "sp", features = c("FABP4", "SPP1"), max.cutoff = 1) + coord_equal()
vuild110 <- subset(mye, subset = sample == "VUILD110")
DimPlot(vuild110, reduction = "sp", group.by = "leiden_res0.8", cols = list(`0` = "red")) + coord_equal()
vuild48 <- subset(mye, subset = sample %in% c("VUILD48LF", "VUILD48MF"))
DimPlot(vuild48, reduction = "sp", group.by = "leiden_res0.8", cols = list(`0` = "red")) + coord_equal()
FeaturePlot(vuild48, reduction = "sp", features = c("IFIT1", "IFIT3")) + coord_equal()


######### VIEW SUBCLUSTERING RESULTS ----
# Subcluster 4 for cDC2s
sub4 <- subset(mye, subset = leiden_res0.8 == "4")
DimPlot(sub4, group.by = "leiden_res0.8_sub4_res0.1", label = TRUE) + coord_equal()
DimPlot(sub4, group.by = "leiden_res0.8_sub4_res0.1_sub0_res0.5", label = TRUE) + coord_equal()
DimPlot(sub4, group.by = "leiden_res0.8_sub4_res0.1", split.by = "leiden_res0.8_sub4_res0.1") + coord_equal()
sub4 <- SetIdent(sub4, value = "leiden_res0.8_sub4_res0.2")
sub4_markers <- FindAllMarkers(sub4, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
FeaturePlot(sub4, features =  c("CD1A", "CD1C", "HLA-DQA1"), pt.size = 0.001) & coord_equal()
FeaturePlot(sub4, features =  c("CD1A", "CD1C", "HLA-DQA1"), 
            split.by = "leiden_res0.8_sub4_res0.1_sub0_res0.5", pt.size = 0.001) & coord_equal()
FeaturePlot(sub4, features =  c("CD1A", "CD1C", "HLA-DQA1", "LTB", "CD4"), max.cutoff = 2,
            split.by = "leiden_res0.8_sub4_res0.5", pt.size = 0.001) & coord_equal()
DotPlot(sub4, features = c("CD1A", "CD1C", "HLA-DQA1", "BANK1", "CD86", "ITGAM", 
                           "CD2", "CD14", "MRC1", "FCER1A", "LTB", "CD4", "VIM"),
        group.by = "leiden_res0.8_sub4_res0.5") + theme_angle

# Subcluster 4,0 further
DimPlot(sub4, group.by = "leiden_res0.8_sub4_res0.1_sub0_res0.5", label = TRUE) + coord_equal()
FeaturePlot(sub4, features =  c("CD1A", "CD1C", "HLA-DQA1", "LTB", "CD4"), max.cutoff = 2,
            split.by = "leiden_res0.8_sub4_res0.1_sub0_res0.4", pt.size = 0.001) & coord_equal()
DotPlot(sub4, features = c("CD1A", "CD1C", "HLA-DQA1", "BANK1", "CD86", "ITGAM", 
                           "CD2", "CD14", "MRC1", "FCER1A", "LTB", "CD4", "VIM"),
        group.by = "leiden_res0.8_sub4_res0.5") + theme_angle
FeaturePlot(sub4, features =  c("CD1A", "CD1C", "HLA-DQA1", "VIM"), order = TRUE,
            split.by = "leiden_res0.8_sub4_res0.1_sub0_res0.5", pt.size = 0.001) & coord_equal()
DotPlot(sub4, features = immune_features, group.by = "leiden_res0.8_sub4_res0.5") + theme_angle
DotPlot(mye, features = immune_features, group.by = "leiden_res0.8_sub4_res0.1_sub0_res0.4") + theme_angle
sub4 <- SetIdent(sub4, value = "leiden_res0.8_sub4_res0.1_sub0_res0.4")
sub4_markers <- FindAllMarkers(sub4, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
VlnPlot(sub4, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.8_sub4_res0.1_sub0_res0.4", pt.size = 0, ncol = 4, same.y.lims = TRUE)
DotPlot(mye, features = c(immune_features, epithelial_features), group.by = "leiden_res0.8_sub4_res0.1_sub0_res0.4") + 
  theme_angle + theme(axis.text.x = element_text(size = 6))
DimPlot(sub4, group.by = "leiden_res0.8_sub4_res0.1_sub0_res0.4", reduction = "sp", 
        cols = c(distinctColorPalette(5), "red", distinctColorPalette(2))) + coord_equal()
FeaturePlot(sub4, features =  c("SFTPC", "SCGB3A2", "KRT8"), order = TRUE,
            split.by = "leiden_res0.8_sub4_res0.1_sub0_res0.5", pt.size = 0.001) & coord_equal()
table(sub4$leiden_res0.8_sub4_res0.1_sub0_res0.5)
# 4,0,0 + 4,0,4: Langerhans DC
# 4,0,1 + 4,0,2 + 4,0,5: cDC2s 
# 4,0,3: CCL22+ cDC1s
# 4,0,6: pDCs
# 4,1: Transitional AT2

# Subcluster 14 to split T-cells and NK cells
sub14 <- subset(mye, subset = leiden_res0.8 == "14")
FeaturePlot(sub14, features = c("CD3E", "CD4", "CD8A", "GNLY", "GZMB"), pt.size = 0.01, order = TRUE,
            split.by = "leiden_res0.8_sub14_res0.4") & coord_equal()
DotPlot(sub14, group.by = "leiden_res0.8_sub14_res0.4", features = c(immune_features, epithelial_features)) + theme_angle
VlnPlot(sub14, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.8_sub14_res0.4", pt.size = 0, ncol = 4, same.y.lims = TRUE)
# Can't split this cluster

# Subcluster all clusters without clear separation into FCN1+, S100A+, SPP1+, or FABP4+
# Subcluster 2
sub2 <- subset(mye, subset = leiden_res0.8 == "2")
DimPlot(sub2, group.by = "leiden_res0.8_sub2_res0.5") + coord_equal()
DotPlot(sub2, group.by = "leiden_res0.8_sub2_res0.5", features = immune_features) + theme_angle
DotPlot(sub2, group.by = "leiden_res0.8_sub2_res0.5", features = c("FABP4", "SPP1", "S100A8", "S100A9", "S100A12", "FCN1", "MARCO", "LYZ")) + theme_angle
VlnPlot(sub2, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.8_sub2_res0.5", pt.size = 0, ncol = 4, same.y.lims = TRUE)
# Subcluster 6
sub6 <- subset(mye, subset = leiden_res0.8 == "6")
DimPlot(sub6, group.by = "leiden_res0.8_sub6_res0.5") + coord_equal()
DotPlot(sub6, group.by = "leiden_res0.8_sub6_res0.5", features = immune_features) + theme_angle
DotPlot(sub6, group.by = "leiden_res0.8_sub6_res0.5", features = c("FABP4", "SPP1", "S100A8", "S100A9", "S100A12", "FCN1", "MARCO", "LYZ")) + theme_angle
VlnPlot(sub6, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.8_sub6_res0.5", pt.size = 0, ncol = 4, same.y.lims = TRUE)
# Subcluster 7
sub7 <- subset(mye, subset = leiden_res0.8 == "7")
DimPlot(sub7, group.by = "leiden_res0.8_sub7_res0.2") + coord_equal()
DotPlot(sub7, group.by = "leiden_res0.8_sub7_res0.2", features = immune_features) + theme_angle
DotPlot(sub7, group.by = "leiden_res0.8_sub7_res0.2", features = c("FABP4", "SPP1", "S100A8", "S100A9", "S100A12", "FCN1", "MARCO", "LYZ")) + theme_angle
VlnPlot(sub7, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.8_sub7_res0.2", pt.size = 0, ncol = 4, same.y.lims = TRUE)
# Subcluster 9 - does not add any meaning
sub9 <- subset(mye, subset = leiden_res0.8 == "9")
DimPlot(sub9, group.by = "leiden_res0.8_sub9_res0.5") + coord_equal()
DotPlot(sub9, group.by = "leiden_res0.8_sub9_res0.5", features = immune_features) + theme_angle
DotPlot(sub9, group.by = "leiden_res0.8_sub9_res0.5", features = c("FABP4", "SPP1", "S100A8", "S100A9", "S100A12", "FCN1", "MARCO", "LYZ")) + theme_angle
VlnPlot(sub9, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.8_sub9_res0.5", pt.size = 0, ncol = 4, same.y.lims = TRUE)
# Subcluster 15
sub15 <- subset(mye, subset = leiden_res0.8 == "15")
DimPlot(sub15, group.by = "leiden_res0.8_sub15_res0.1") + coord_equal()
DotPlot(sub15, group.by = "leiden_res0.8_sub15_res0.1", features = immune_features) + theme_angle
DotPlot(sub15, group.by = "leiden_res0.8_sub15_res0.1", features = c("FABP4", "SPP1", "S100A8", "S100A9", "S100A12", "FCN1", "MARCO", "LYZ")) + theme_angle
VlnPlot(sub15, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.8_sub15_res0.1", pt.size = 0, ncol = 4, same.y.lims = TRUE)
# Subcluster 16
sub16 <- subset(mye, subset = leiden_res0.8 == "16")
DimPlot(sub16, group.by = "leiden_res0.8_sub16_res0.2") + coord_equal()
DotPlot(sub16, group.by = "leiden_res0.8_sub16_res0.2", features = immune_features) + theme_angle
DotPlot(sub16, group.by = "leiden_res0.8_sub16_res0.2", features = c("FABP4", "SPP1", "S100A8", "S100A9", "S100A12", "FCN1", "MARCO", "LYZ")) + theme_angle
VlnPlot(sub16, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "leiden_res0.8_sub16_res0.2", pt.size = 0, ncol = 4, same.y.lims = TRUE)
# Create merged variable to view in DotPlots
mye$merged_subcluster_var <- mye$leiden_res0.8
mye$merged_subcluster_var[mye$leiden_res0.8 == "2"] <- mye$leiden_res0.8_sub2_res0.5[mye$leiden_res0.8 == "2"]
mye$merged_subcluster_var[mye$leiden_res0.8 == "6"] <- mye$leiden_res0.8_sub6_res0.5[mye$leiden_res0.8 == "6"]
mye$merged_subcluster_var[mye$leiden_res0.8 == "7"] <- mye$leiden_res0.8_sub7_res0.2[mye$leiden_res0.8 == "7"]
mye$merged_subcluster_var[mye$leiden_res0.8 == "15"] <- mye$leiden_res0.8_sub15_res0.1[mye$leiden_res0.8 == "15"]
mye$merged_subcluster_var[mye$leiden_res0.8 == "16"] <- mye$leiden_res0.8_sub16_res0.2[mye$leiden_res0.8 == "16"]
DotPlot(mye, group.by = "merged_subcluster_var", features = c("PPARG", "FABP4", "SPP1", "S100A8", "S100A9", "S100A12", "FCN1", "MARCO", "LYZ")) + theme_angle
DotPlot(mye, group.by = "merged_subcluster_var", features = immune_features) + theme_angle
DotPlot(mye, group.by = "merged_subcluster_var", features = endothelial_features) + theme_angle
DotPlot(mye, group.by = "merged_subcluster_var", features = epithelial_features) + theme_angle
DotPlot(mye, group.by = "merged_subcluster_var", features = mesenchymal_features) + theme_angle
DotPlot(mye, group.by = "merged_subcluster_var", features = c(immune_features, epithelial_features, "MFAP5", "MEG3", "COL1A2")) + 
  theme_angle + theme(axis.text.x = element_text(size = 6))
VlnPlot(mye, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        group.by = "merged_subcluster_var", pt.size = 0, ncol = 2, same.y.lims = TRUE)
sort(table(mye$merged_subcluster_var))
table(mye$merged_subcluster_var, mye$sample)
# S100A+ Macrophages: 15,3 + 2,0 + 7,1
# SPP1+ Macrophages: 16,2 + 6,7
# Plasma: 6,4
# All others: Macrophages

# Note that cDC1s and cDC2s were switched around on CT_firstpass; this will be fixed here
# Create new celltype variable
mye$CT_firstpassB <- ""
mye$CT_firstpassB[mye$leiden_res0.8 %in% c("2", "6", "7", "8", "9", "15", "16")] <- "Macrophages"
mye$CT_firstpassB[mye$leiden_res0.8 == "0"] <- "Interferon-high Macrophages"
mye$CT_firstpassB[mye$leiden_res0.8 == "1" |
                    mye$leiden_res0.8_sub4_res0.1_sub0_res0.5 == "4,0,6"] <- "pDCs (Myeloid Obj)"
mye$CT_firstpassB[mye$leiden_res0.8 == "3"] <- "FABP4+ Macrophages"
mye$CT_firstpassB[mye$leiden_res0.8_sub4_res0.1_sub0_res0.5 %in% c("4,0,0", "4,0,4")] <- "Langerhans DCs"
mye$CT_firstpassB[mye$leiden_res0.8_sub4_res0.1_sub0_res0.5 %in% c("4,0,1", "4,0,2", "4,0,5")] <- "cDC2s"
mye$CT_firstpassB[mye$leiden_res0.8_sub4_res0.1_sub0_res0.5 == "4,1"] <- "Transitional AT2 (Myeloid Obj)"
mye$CT_firstpassB[mye$leiden_res0.8 == "5"] <- "Mast"
mye$CT_firstpassB[mye$leiden_res0.8_sub4_res0.1_sub0_res0.5 == "4,0,3" |
                    mye$leiden_res0.8 == "8"] <- "CCL22+ cDC1s"
mye$CT_firstpassB[mye$leiden_res0.8 == "10"] <- "FCN1+/S100A+ Macrophages"
mye$CT_firstpassB[mye$leiden_res0.8 == "11"] <- "Proliferating Myeloid"
mye$CT_firstpassB[mye$leiden_res0.8 == "12"] <- "S100A+ Macrophages"
mye$CT_firstpassB[mye$leiden_res0.8 == "13"] <- "SPP1+ Macrophages"
mye$CT_firstpassB[mye$leiden_res0.8 == "14"] <- "T-cells (Myeloid Obj)"
mye$CT_firstpassB[mye$merged_subcluster_var == "6,4"] <- "Plasma (Myeloid Obj)"
mye$CT_firstpassB[mye$merged_subcluster_var %in% c("16,2", "6,7")] <- "SPP1+ Macrophages"
mye$CT_firstpassB[mye$merged_subcluster_var %in% c("2,0", "7,1", "15,3")] <- "S100A+ Macrophages"
DimPlot(mye, group.by = "CT_firstpassB", cols = c("grey90", distinctColorPalette(715))) + coord_equal()
DotPlot(mye, features = c(immune_features, epithelial_features), group.by = "CT_firstpassB") + theme_angle + theme(axis.text.x = element_text(size = 6))

# Save object
#saveRDS(mye, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_myeloid_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_061324.rds")


#### MESENCHYMAL ----
mes <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/mesenchymal_lineage_split_nucnuc_RAPIDS_clustered_052824.rds")
DimPlot(mye, group.by = "leiden_res0.3", label = TRUE) + coord_equal()
DimPlot(mye, group.by = "initial_CT", label = TRUE) + coord_equal()

# Basic DimPlots and FeaturePlots
DimPlot(mye, group.by = "sample", cols = random_sample_colors) + coord_equal()
DimPlot(mye, group.by = "tma") + coord_equal()
DimPlot(mye, group.by = "run", split.by = "run") + coord_equal()
DimPlot(mye, group.by = "sample_type", split.by = "sample_type") + coord_equal()
FeaturePlot(mye, features = c("nCount_RNA", "nFeature_RNA")) & coord_equal()
DimPlot(mye, group.by = "leiden_res0.3", split.by = "leiden_res0.3") + coord_equal()

# Create and view module scores
mes$endo_score <- colSums(mes@assays$RNA@counts[endothelial_features, ])/colSums(mes@assays$RNA@counts)
mes$epi_score <- colSums(mes@assays$RNA@counts[epithelial_features, ])/colSums(mes@assays$RNA@counts)
mes$imm_score <- colSums(mes@assays$RNA@counts[immune_features, ])/colSums(mes@assays$RNA@counts)
mes$mes_score <- colSums(mes@assays$RNA@counts[mesenchymal_features, ])/colSums(mes@assays$RNA@counts)
FeaturePlot(mes, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(mes, features = c("endo_score", "epi_score", "imm_score", "mes_score"), group.by = "leiden_res0.3", pt.size = 0, ncol = 2)
# Only cluster 0 does not have strong mesenchymal identity

# Find Markers
mes <- SetIdent(mes, value = "leiden_res0.3")
mes_markers <- FindAllMarkers(mes, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(mes$adj_x_centroid, mes$adj_y_centroid)
row.names(position_xy) <- row.names(mes@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
mes[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                    assay = DefaultAssay(mes))

# Look at lineage and top markers
# Trying to find: activated, peribronchial, alveolar, lipo, adventitial, & general fibroblasts
# Also looking for SMCs/pericytes and mesothelial
FeaturePlot(mes, features = c("EPCAM", "PECAM1", "PTPRC", "DCN", "LUM")) & coord_equal()
FeaturePlot(mes, features = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "FN1",
                              "MSLN", "CTHRC1", "FAP", "CSPG4", "ACTA2",
                              "WNT5A", "SCGB3A2", "FABP4", "SFTPC", "COL15A1", "PI16", 
                              "SFRP2", "SFRP4", "MKI67")) & coord_equal()
FeaturePlot(mes, features = c("MSLN", "HAS1", "KR1T8", "KRT8", # 0: Mesothelial
                              "HLA-DRA", "LYZ", # 1: ?
                              "TOP2A", "CDK1", "CENPF", "MKI67" # 5: Proliferating Fibroblasts
)) & coord_equal()
FeaturePlot(mes, features = c("MSLN", "HAS1", "KR1T8", "KRT8")) & coord_equal() # 0: Mesothelial
DotPlot(mes, features = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1", "FN1",
                          "MSLN", "HAS1", "KRT18", "KRT8",
                          "CTHRC1", "FAP", "CSPG4", "ACTA2", "PLIN2", 
                          "WNT5A", "SCGB3A2", "FABP4", "SFTPC", "PI16", 
                          "SFRP2", "SFRP4", "MKI67", "TOP2A", "CDK1",
                          "FGF2", "FGF10", "FGF7")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(mes, features = mesenchymal_features) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 0: Mesothelial
# 1,2,3: General fibroblasts
# 4,6: SMCs/Pericytes, and 6 contains lipofibroblasts
# 5: Proliferating Fibroblasts
# 7: Contains adventitial and activated fibroblasts, with some peribronchial and HAS1+ FBs
# 8: Contains peribronchial and alveolar fibroblasts

# View in spatial context
DimPlot(mes, reduction = "sp", group.by ="leiden_res0.3", cols = distinctColorPalette(9)) + coord_equal()
table(mes$leiden_res0.3, mes$sample)
## Specific example 1: VUILD106
# 4 clearly marks SMCs, and 0 marks mesothelial (pleural edge)
# Peribronchial FBs are not clearly marked by cluster 8 in this sample
vuild106 <- subset(mes, subset = sample == "VUILD106")
DimPlot(vuild106, reduction = "sp", group.by = "leiden_res0.3", cols = distinctColorPalette(9)) + coord_equal()
FeaturePlot(vuild106, features = c("WNT5A", "PI16", "CTHRC1", "FABP4"), reduction = "sp") & coord_equal()
## Specific example 2: VUILD115
# 7 marks adventitial FBs, 4 clearly marks SMCs, and 8 marks some alveolar FBs
vuild115 <- subset(mes, subset = sample == "VUILD115")
DimPlot(vuild115, reduction = "sp", group.by = "leiden_res0.3", cols = distinctColorPalette(9)) + coord_equal()
DimPlot(vuild115, reduction = "sp", group.by = "leiden_res0.3", cols = list(`7` = "blue", `8` = "red")) + coord_equal()
FeaturePlot(vuild115, features = c("ACTA2", "PI16", "CTHRC1", "MSLN", "HAS1"), reduction = "sp") & coord_equal()
## Specific example 3: VUHD116B
# In this sample, 8 seems to mark alveolar FBs
vuhd116b <- subset(mes, subset = sample == "VUHD116B")
DimPlot(vuhd116b, reduction = "sp", group.by = "leiden_res0.3", cols = distinctColorPalette(9)) + coord_equal()
FeaturePlot(vuhd116b, features = c("ACTA2", "PI16", "CTHRC1", "SFTPC"), reduction = "sp") & coord_equal()
# Specific example 4: TILD117MFB
# 0 clearly marks mesothelial in this sample; 7 contains some HAS1+ FBs
tild117mf <- subset(mes, subset = sample == "TILD117MF" & tma == "TMA5")
DimPlot(tild117mf, reduction = "sp", group.by = "leiden_res0.3", cols = distinctColorPalette(9)) + coord_equal()
FeaturePlot(tild117mf, features = c("ACTA2", "PI16", "MSLN", "HAS1"), reduction = "sp") & coord_equal()
## Specific example 5: VUILD110
# 4 marks SMCs, including airway SMCs; 7 marks some peribronchial
vuild110 <- subset(mes, subset = sample == "VUILD110")
DimPlot(vuild110, reduction = "sp", group.by = "leiden_res0.3", cols = distinctColorPalette(9)) + coord_equal()
DimPlot(vuild110, reduction = "sp", group.by = "leiden_res0.3", cols = list(`4` = "blue", `7` = "red", `8` = "cyan")) + coord_equal()
FeaturePlot(vuild110, features = c("ACTA2", "CSPG4", "WNT5A", "SFTPC"), reduction = "sp") & coord_equal()


###### VIEW SUBCLUSTERING RESULTS ----
# Subcluster 7 to split adventitial, activated, peribronchial, and HAS1+ FBs
DimPlot(mes, group.by = "leiden_res0.3_sub7_res0.4", label = TRUE) + coord_equal()
FeaturePlot(mes, features = c("ACTA2", "PDGFRA", "CTHRC1", "FAP")) & coord_equal()
sub7 <- subset(mes, subset = leiden_res0.3 == "7")
table(sub7$sample, sub7$leiden_res0.3_sub7_res0.4)
DimPlot(sub7, group.by = "leiden_res0.3_sub7_res0.4", label = TRUE, repel = TRUE) + coord_equal()
DimPlot(sub7, group.by = "leiden_res0.3_sub7_res0.4", split.by = "leiden_res0.3_sub7_res0.4") & coord_equal()
FeaturePlot(sub7, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(sub7, features = c("endo_score", "epi_score", "imm_score", "mes_score"), 
        same.y.lims = TRUE, group.by = "leiden_res0.3_sub7_res0.4", ncol = 2, pt.size = 0)
FeaturePlot(sub7, features = c("PI16", "SFRP4", "CTHRC1", "FAP", "WNT5A", "HAS1")) & coord_equal()
DotPlot(sub7, features = c("ACTA2", "PDGFRA", "SFRP2", "PI16", "SFRP4", "MFAP5", "CTHRC1", "FAP",
                           "WNT5A", "SCGB3A2", "SCGB1A1", "HAS1", "POSTN", "EPAS1", "HIF1A", "YAP1", 
                           "JCHAIN", "XBP1", "TRAC", "IL7R",
                           "FGF7", "TGFB3", "UBE2J1", "ITGAX", "HLA-DRA"), 
        group.by = "leiden_res0.3_sub7_res0.4") + theme_angle
FeaturePlot(sub7, features = c("PI16", "SFRP4", "CTHRC1", "FAP", "WNT5A", "HAS1"), split.by = "leiden_res0.3_sub7_res0.4") & coord_equal()
DimPlot(vuild110, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4", cols = list(`4` = "blue", `7,2` = "orange", `7,5` = "cyan", `7,7` = "red")) + coord_equal()
DimPlot(vuild106, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4", cols = list(`4` = "blue", `7,2` = "orange", `7,5` = "cyan", `7,7` = "red")) + coord_equal()
DimPlot(vuild106, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4", cols = list(`4` = "blue", `7,5` = "cyan", `7,0` = "red")) + coord_equal()
FeaturePlot(vuild106, features = c("ACTA2", "CTHRC1", "FAP"), reduction = "sp", ncol = 3, max.cutoff = 1) & coord_equal()
DimPlot(tild117mf, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4", cols = list(`4` = "blue", `7,5` = "cyan", `7,4` = "red")) + coord_equal()
FeaturePlot(tild117mf, features = c("ACTA2", "PDGFRA"), reduction = "sp", max.cutoff = 1) & coord_equal() 
DimPlot(tild117mf, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4", cols = list(`4` = "blue", `7,2` = "orange", `7,5` = "cyan", `7,7` = "red")) + coord_equal()
vuild91lf <- subset(mes, subset = sample == "VUILD91LF")
DimPlot(vuild91lf, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4", cols = list(`4` = "blue", `7,2` = "orange", `7,5` = "cyan", `7,0` = "red")) + coord_equal()
FeaturePlot(vuild91lf, features = c("ACTA2", "CTHRC1", "FAP"), reduction = "sp", ncol = 3, max.cutoff = 1) & coord_equal()
vuild91mf <- subset(mes, subset = sample == "VUILD91MF")
DimPlot(vuild91mf, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4", cols = list(`4` = "blue", `7,2` = "orange", `7,5` = "cyan", `7,0` = "red")) + coord_equal()
FeaturePlot(vuild91mf, features = c("ACTA2", "CTHRC1", "FAP"), reduction = "sp", ncol = 3, max.cutoff = 1) & coord_equal()
sub7 <- SetIdent(sub7, value = "leiden_res0.3_sub7_res0.4")
sub7_markers <- FindAllMarkers(sub7, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# 7,1 + 7,3 + 7,6: General fibroblasts
# 7,5: Activated MyoFBs
# 7,7: Peribronchial FBs
# 7,4: Hypoxic FBs

# Split 7_0 further to split out activated fibroblasts
sub7_0 <- subset(mes, subset = leiden_res0.3_sub7_res0.4 == "7,0")
DimPlot(sub7_0, group.by = "leiden_res0.3_sub7_res0.4_sub0_res0.15", label = TRUE) + coord_equal()
FeaturePlot(sub7_0, features = c("CTHRC1", "FAP")) & coord_equal()
FeaturePlot(sub7_0, features = c("CTHRC1", "FAP", "ACTA2"), split.by = "leiden_res0.3_sub7_res0.4_sub0_res0.15") & coord_equal()
DotPlot(sub7_0, features = c("CTHRC1", "FAP", "ACTA2"), group.by = "leiden_res0.3_sub7_res0.4_sub0_res0.15")
DotPlot(sub7_0, features = c("CTHRC1", "FAP", "ACTA2", "FGF7", "SFRP2", "TGFB3", "UBE2J1", "ITGAX"), group.by = "leiden_res0.3_sub7_res0.4_sub0_res0.15")
# Can't split further
# 7,0: Activated MyoFBs

# Split 7_2 further into adventitial and HAS1+ FBs
sub7_2 <- subset(mes, subset = leiden_res0.3_sub7_res0.4 == "7,2")
DimPlot(sub7_2, group.by = "leiden_res0.3_sub7_res0.4_sub2_res0.4", label = TRUE) + coord_equal()
FeaturePlot(sub7_2, features = c("PI16", "MFAP5", "SFRP4", "FABP4", "HAS1", "PLIN2")) & coord_equal()
FeaturePlot(sub7_2, features = c("PI16", "MFAP5", "SFRP4", "FABP4", "HAS1", "PLIN2"), 
            split.by = "leiden_res0.3_sub7_res0.4_sub2_res0.4", pt.size = 0.01) & coord_equal()
DotPlot(sub7_2, features = c("PI16", "MFAP5", "SFRP4", "FABP4", "HAS1"), group.by = "leiden_res0.3_sub7_res0.4_sub2_res0.4")
VlnPlot(sub7_2, features = c("PI16", "MFAP5", "SFRP4", "FABP4", "HAS1", "PLIN2"), group.by = "leiden_res0.3_sub7_res0.4_sub2_res0.4", pt.size = 0, ncol = 2)
VlnPlot(sub7_2, features = c("nCount_RNA", "nFeature_RNA"), group.by = "leiden_res0.3_sub7_res0.4_sub2_res0.4", pt.size = 0, ncol = 2)
DimPlot(sub7_2, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4_sub2_res0.4") + coord_equal()
DimPlot(subset(sub7_2, subset = tma == "TMA3"), reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4_sub2_res0.4") + coord_equal()
vuild106 <- subset(sub7_2, subset = sample == "VUILD106")
DimPlot(vuild106, reduction = "sp", group.by = "leiden_res0.3_sub7_res0.4_sub2_res0.4", cols = distinctColorPalette(4)) + coord_equal()
FeaturePlot(vuild106, features = c("PI16", "MFAP5", "SFRP4", "FABP4", "HAS1", "PLIN2"), reduction = "sp", max.cutoff = 1, pt.size = 0.01) & coord_equal()
sub7_2 <- SetIdent(sub7_2, value = "leiden_res0.3_sub7_res0.4_sub2_res0.4")
sub7_2_markers <- FindAllMarkers(sub7_2, only.pos = TRUE) %>% 
  mutate(pct.diff = abs(pct.1 - pct.2), .after = "pct.2")
# 7,2,0 = General Fibroblasts
# 7,2,1 and 7,2,3 = Adventitial FBs
# 7,2,2 = HAS1+/PLIN2+ FBs

# Compare to eQTL dataset
DotPlot(eqtl.ref, features = mesenchymal_features, group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(eqtl.ref, features = immune_features, group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
eqtl.mes <- subset(eqtl.ref, subset = lineage %in% c("MyoFB", "PLIN2+ FB", "Mesenchymal", "MyoFB - Activated"))
DotPlot(eqtl.mes, features = rownames(mes)[1:150], group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
                                                                                             axis.text.y = element_text(size = 5))
DotPlot(eqtl.mes, features = rownames(mes)[151:343], group.by = "manual_annotation_1") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
                                                                                               axis.text.y = element_text(size = 5))
# Subcluster 4 to split SMCS/pericytes
DimPlot(mes, group.by = "leiden_res0.3_sub7_res0.4", label = TRUE) + coord_equal()
FeaturePlot(mes, features = c("ACTA2", "PDGFRA", "CTHRC1", "FAP")) & coord_equal()
sub4 <- subset(mes, subset = leiden_res0.3 == "4")
table(sub4$sample, sub4$leiden_res0.3_sub4_res0.2)
DimPlot(sub4, group.by = "leiden_res0.3_sub4_res0.2", label = TRUE, repel = TRUE) + coord_equal()
FeaturePlot(sub4, features = c("CSPG4", "ACTA2"), blend = TRUE) & coord_equal()
DimPlot(sub4, group.by = "leiden_res0.3_sub4_res0.2", split.by = "leiden_res0.3_sub4_res0.2") & coord_equal()
FeaturePlot(sub4, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
VlnPlot(sub4, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        same.y.lims = TRUE, group.by = "leiden_res0.3_sub7_res0.4", ncol = 2, pt.size = 0)
FeaturePlot(sub4, features = c("PI16", "SFRP4", "CTHRC1", "FAP", "WNT5A", "HAS1")) & coord_equal()
DotPlot(sub4, features = c("ACTA2", "PDGFRB", "CSPG4"),
        group.by = "leiden_res0.3_sub4_res0.3") + theme_angle
DimPlot(sub4, reduction = "sp", group.by = "")
vuild106 <- subset(sub4, subset = sample == "VUILD106")
DimPlot(vuild106, reduction = "sp", group.by = "leiden_res0.3_sub4_res0.2", cols = c("red", "blue", "cyan")) + coord_equal()
FeaturePlot(vuild106, features = c("ACTA2", "PDGFRB", "CSPG4"), reduction = "sp") & coord_equal()
# Can't split SMCs vs. Pericytes

# Subcluster 6 to split SMCs/pericytes and lipofibroblasts
sub6 <- subset(mes, subset = leiden_res0.3 == "6")
DimPlot(sub6, group.by = "leiden_res0.3_sub6_res0.15", label = TRUE, repel = TRUE) + coord_equal()
FeaturePlot(sub6, features = c("CSPG4", "ACTA2"), blend = TRUE) & coord_equal()
FeaturePlot(sub6, features = c("CSPG4", "ACTA2", "FABP4")) & coord_equal()
FeaturePlot(sub6, features = c("CSPG4", "ACTA2", "FABP4"), split.by = "leiden_res0.3_sub6_res0.15") & coord_equal()
vuild106 <- subset(sub6, subset = sample == "VUILD106")
DimPlot(vuild106, reduction = "sp", group.by = "leiden_res0.3_sub6_res0.2", cols = c("yellow", "blue", "cyan", "red")) + coord_equal()
FeaturePlot(vuild106, features = c("ACTA2", "PDGFRB", "CSPG4", "FABP4"), reduction = "sp") & coord_equal()
DotPlot(sub6, features = c("ACTA2", "PDGFRB", "CSPG4", "FABP4"), group.by = "leiden_res0.3_sub6_res0.15")
# Split 6,0 further to get lipofibroblasts
# Can't split SMCs/pericytes
sub6_0 <- subset(mes, subset = leiden_res0.3_sub6_res0.15 == "6,0")
DimPlot(sub6_0, group.by = "leiden_res0.3_sub6_res0.15_sub0_res0.3", label = TRUE, repel = TRUE) + coord_equal()
DotPlot(sub6_0, features = c("ACTA2", "PDGFRB", "CSPG4", "FABP4"), group.by = "leiden_res0.3_sub6_res0.15_sub0_res0.3")
DotPlot(sub6, features = c("ACTA2", "PDGFRB", "CSPG4", "FABP4"), group.by = "leiden_res0.3_sub6_res0.15_sub0_res0.3")
DotPlot(mes, features = c("ACTA2", "PDGFRB", "CSPG4", "FABP4"), group.by = "leiden_res0.3_sub6_res0.15_sub0_res0.3")
VlnPlot(sub6_0, features = c("ACTA2", "PDGFRB", "CSPG4", "FABP4", "CCL21", "PPARG"), group.by = "leiden_res0.3_sub6_res0.15_sub0_res0.3", pt.size = 0)
FeaturePlot(sub6_0, features = c("ACTA2", "FABP4"), blend = TRUE) & coord_equal()
vuild106 <- subset(mes, subset = sample == "VUILD106")
DimPlot(vuild106, reduction = "sp", group.by = "leiden_res0.3_sub6_res0.15_sub0_res0.3", 
        cols = list(`4` = "grey90", `6,0,0` = "red", `6,0,1` = "pink", `6,0,2` = "cyan", `6,0,3` = "blue")) + coord_equal()
FeaturePlot(vuild106, features = c("FABP4", "PPARG"), max.cutoff = 1, reduction = "sp") & coord_equal()
# 6,0,0 + 6,0,3 = Lipofibroblasts
# 6,0,1 + 6,0,2 + 6,1 = General FBs
# 6,2 = SMCs/Pericytes

# Subcluster 8 to split peribronchial and alveolar fibroblasts
sub8 <- subset(mes, subset = leiden_res0.3 == "8")
table(sub8$leiden_res0.3_sub8_res0.2, sub8$sample)
DimPlot(sub8, group.by = "leiden_res0.3_sub8_res0.2", label = TRUE, repel = TRUE) + coord_equal()
FeaturePlot(sub8, features = c("SCGB3A2", "SCGB1A1", "WNT5A", "AXIN2", "SFTPC", "PDGFRA")) & coord_equal()
FeaturePlot(sub8, features = c("endo_score", "epi_score", "imm_score", "mes_score"), keep.scale = "all") & coord_equal()
DotPlot(sub8, features = c("SCGB3A2", "SCGB1A1", "WNT5A", "AXIN2", "SFTPC", "PDGFRA"), group.by = "leiden_res0.3_sub8_res0.2")
VlnPlot(sub8, features = c("SCGB3A2", "SCGB1A1", "WNT5A", "AXIN2", "SFTPC", "PDGFRA"), group.by = "leiden_res0.3_sub8_res0.2", pt.size = 0)
vuild110 <- subset(mes, subset = sample == "VUILD110")
DimPlot(vuild110, reduction = "sp", group.by = "leiden_res0.3_sub8_res0.2", 
        cols = list(`7` = "grey90", `8,0` = "blue", `8,1` = "orange", `8,2` = "red")) + coord_equal()
# 8,0: Peribronchial FBs
# 8,1 + 8,2 = Alveolar FBs

# Create cell type variable
mes$CT_firstpass <- ""
mes$CT_firstpass[mes$leiden_res0.3 == "0"] <- "Mesothelial"
mes$CT_firstpass[mes$leiden_res0.3 %in% c("1", "2", "3") |
                   mes$leiden_res0.3_sub6_res0.15_sub0_res0.3 %in% c("6,0,1", "6,0,2", "6,1") |
                   mes$leiden_res0.3_sub7_res0.4 %in% c("7,1", "7,3", "7,6") |
                   mes$leiden_res0.3_sub7_res0.4_sub2_res0.4 == "7,2,0"] <- "Fibroblasts"
mes$CT_firstpass[mes$leiden_res0.3 == "4" | 
                   mes$leiden_res0.3_sub6_res0.15 == "6,2"] <- "SMCs/Pericytes"
mes$CT_firstpass[mes$leiden_res0.3 == "5"] <- "Proliferating Fibroblasts"
mes$CT_firstpass[mes$leiden_res0.3_sub6_res0.15_sub0_res0.3 %in% c("6,0,0", "6,0,3")] <- "Lipofibroblasts"
mes$CT_firstpass[mes$leiden_res0.3_sub7_res0.4 %in% c("7,0", "7,5")] <- "Activated MyoFBs"
mes$CT_firstpass[mes$leiden_res0.3_sub7_res0.4 == "7,4"] <- "Hypoxic FBs"
mes$CT_firstpass[mes$leiden_res0.3_sub7_res0.4_sub2_res0.4 %in% c("7,2,1", "7,2,3")] <- "Adventitial FBs"
mes$CT_firstpass[mes$leiden_res0.3_sub7_res0.4_sub2_res0.4 == "7,2,2"] <- "HAS1+/PLIN2+ FBs"
mes$CT_firstpass[mes$leiden_res0.3_sub7_res0.4 == "7,7" |
                   mes$leiden_res0.3_sub8_res0.2 == "8,0"] <- "Peribronchial FBs"
mes$CT_firstpass[mes$leiden_res0.3_sub8_res0.2 %in% c("8,1", "8,2")] <- "Alveolar FBs"

# Plot
DimPlot(mes, group.by = "CT_firstpass", cols = distinctColorPalette(11)) + coord_equal()
FeaturePlot(mes, features = c("CTHRC1", "FAP")) + coord_equal()
DimPlot(mes, group.by = "leiden_res0.3", cols = distinctColorPalette(10), label = TRUE)
table(mes$CT_firstpass, mes$leiden_res0.3)
VlnPlot(mes, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        same.y.lims = TRUE, group.by = "CT_firstpass", ncol = 2, pt.size = 0)
DotPlot(mes, features = c(mesenchymal_features, "CSPG4", "HIF1A", "EPAS1", "YAP1", "FABP4", "PPARG",
                          "SCGB3A2", "SCGB1A1", "SFTPC", "AGER", "RTKN2", "AXIN2", "PLIN2", "MSLN", 
                          "KRT8", "MKI67", "TOP2A", "CDK1"), group.by = "CT_firstpass") + theme_angle

# Save object
#saveRDS(mes, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/mesenchymal_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_060424.rds")



## CELL TYPE ANNOTATION (SECONDPASS) ----
# Load original object
full_merged_nucleus_nucleus <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/full_merged_nucleus_nucleus_RAPIDS_clustered_annotated_052824.rds")

# Load in all firstpass objects
endo <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/endothelial_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_061324.rds")
air <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_airway_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS-B_LABELED_060424.rds")
alv <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_alveolar_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS-B_LABELED_061324.rds")
lym <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lymphoid_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS-B_LABELED_060524.rds")
mye <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_myeloid_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_061324.rds")
mes <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/mesenchymal_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_060424.rds")

# Load in stray clusters from epithelial and immune objects
strays_from_epi <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/epithelial_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_052824_OTHER_LINEAGE.rds")
strays_from_imm <- readRDS("/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/immune_lineage_split_nucnuc_RAPIDS_clustered_FIRSTPASS_LABELED_052824_OTHER_LINEAGE.rds")

# Merge all together
all_merged <- merge(endo, y = list(air, alv, lym, mye, mes, strays_from_epi, strays_from_imm))
all_merged$CT_firstpass_full <- all_merged$CT_firstpassB
all_merged$CT_firstpass_full[is.na(all_merged$CT_firstpassB)] <- all_merged$CT_firstpass[is.na(all_merged$CT_firstpassB)]

# Check that all original cells are present
ncol(endo) + ncol(air) + ncol(alv) + ncol(lym) + ncol(mye) + ncol(mes) +
  ncol(strays_from_epi) + ncol(strays_from_imm) == ncol(full_merged_nucleus_nucleus)
ncol(all_merged) == ncol(full_merged_nucleus_nucleus)

# Create module scores
all_merged$endo_score <- colSums(all_merged@assays$RNA@counts[endothelial_features, ])/colSums(all_merged@assays$RNA@counts)
all_merged$epi_score <- colSums(all_merged@assays$RNA@counts[epithelial_features, ])/colSums(all_merged@assays$RNA@counts)
all_merged$imm_score <- colSums(all_merged@assays$RNA@counts[immune_features, ])/colSums(all_merged@assays$RNA@counts)
all_merged$mes_score <- colSums(all_merged@assays$RNA@counts[mesenchymal_features, ])/colSums(all_merged@assays$RNA@counts)

# Add spatial coordinates as dimension reduction objects
position_xy <- cbind(all_merged$adj_x_centroid, all_merged$adj_y_centroid)
row.names(position_xy) <- row.names(all_merged@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
all_merged[["sp"]] <- CreateDimReducObject(embeddings = position_xy, key = "SP_",
                                           assay = DefaultAssay(all_merged))

# Fix select cell type labels
all_merged$CT_firstpass_full[all_merged$CT_firstpass_full == "Proliferating Epithelial"] <- "Proliferating Epithelial (Imm Obj)"
all_merged$CT_firstpass_full[all_merged$CT_firstpass_full == "Differentating Ciliated"] <- "Differentiating Ciliated"


# Separate by lineage/sub-lineage
all_endo <- subset(all_merged, subset = CT_firstpass_full %in% c("Venous", "Capillary", "Lymphatic", "Arteriole", "Venous (Epi Obj)"))
all_air <- subset(all_merged, subset = CT_firstpass_full %in% c("SCGB3A2+", "SCGB3A2+/SCGB1A1+", "Ciliated", "PNEC",
                                                                "Differentiating Ciliated", "Basal", "MUC5B+", "Interferon-high Basal",
                                                                "SCGB3A2+ (Alv Obj)", "Ciliated (Alv Obj)", "Basal (Alv Obj)",
                                                                "Proliferating Basal (Alv Obj)", "Proliferating Differentiating Ciliated (Alv Obj)",
                                                                "MUC5B+ (Alv Obj)", "Proliferating SCGB3A2+/SCGB1A1+ (Alv Obj)"))
all_alv <- subset(all_merged, subset = CT_firstpass_full %in% c("AT2", "Transitional AT2", "AT1", "KRT5-/KRT17+",
                                                                "Proliferating Transitional AT2", "Proliferating AT2",
                                                                "Transitional AT2 (Air Obj)", "AT2 (Air Obj)",
                                                                "Transitional AT2 (Myeloid Obj)", "Proliferating Epithelial (Imm Obj)"))
all_lym <- subset(all_merged, subset = CT_firstpass_full %in% c("Tregs", "CD4+ T-cells", "NK cells", "CD8+ T-cells", "T-cells",
                                                                "Proliferating T-cells", "B cells", "Granulomatous T-cells",
                                                                "Plasma", "pDCs", "Proliferating NK cells", "Proliferating B cells",
                                                                "T-cells (Endo Obj)", "T-cells (Myeloid Obj)", 
                                                                "pDCs (Myeloid Obj)", "Plasma (Myeloid Obj)"))
all_mye <- subset(all_merged, subset = CT_firstpass_full %in% c("Macrophages", "FABP4+ Macrophages", "SPP1+ Macrophages",
                                                                "S100A+ Macrophages", "Proliferating Myeloid", "CCL22+ cDC1s",
                                                                "Interferon-high Macrophages", "Mast", "FCN1+/S100A+ Macrophages",
                                                                "cDC2s", "Langerhans DCs", "cDC1 (Air Obj)", "cDC2 (Air Obj)",
                                                                "Mast (Lymphoid Obj)", "Mast (Epi Obj)"))
all_mes <- subset(all_merged, subset = CT_firstpass_full %in% c("SMCs/Pericytes", "Fibroblasts", "Proliferating Fibroblasts",
                                                                "Alveolar FBs", "Hypoxic FBs", "Activated MyoFBs", 
                                                                "Peribronchial FBs", "HAS1+/PLIN2+ FBs", "Lipofibroblasts",
                                                                "Adventitial FBs", "Mesothelial", "Fibroblasts (Endo Obj)",
                                                                "Proliferating Fibroblasts (Air Obj)"))
# Check that all cells were split
ncol(all_endo) + ncol(all_air) + ncol(all_alv) + 
  ncol(all_lym) + ncol(all_mye) + ncol(all_mes) == ncol(full_merged_nucleus_nucleus)

# Check module scores for each cell type
VlnPlot(all_merged, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "CT_firstpass_full", pt.size = 0, ncol = 2, same.y.lims = TRUE) &
  theme(axis.text = element_text(size = 7), axis.title.x = element_blank())
VlnPlot(all_endo, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "CT_firstpass_full", pt.size = 0, ncol = 2, same.y.lims = TRUE) & theme(axis.title.x = element_blank())
VlnPlot(all_air, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "CT_firstpass_full", pt.size = 0, ncol = 2, same.y.lims = TRUE) & theme(axis.title.x = element_blank())
VlnPlot(all_alv, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "CT_firstpass_full", pt.size = 0, ncol = 2, same.y.lims = TRUE) & theme(axis.title.x = element_blank())
VlnPlot(all_lym, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "CT_firstpass_full", pt.size = 0, ncol = 2, same.y.lims = TRUE) & theme(axis.title.x = element_blank())
VlnPlot(all_mye, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "CT_firstpass_full", pt.size = 0, ncol = 2, same.y.lims = TRUE) & theme(axis.title.x = element_blank())
VlnPlot(all_mes, features = c("endo_score", "epi_score", "imm_score", "mes_score"),
        group.by = "CT_firstpass_full", pt.size = 0, ncol = 2, same.y.lims = TRUE) & theme(axis.title.x = element_blank())


## RENAME/RECATEGORIZE CELL TYPES ----
all_epi <- merge(all_air, y = all_alv)
all_imm <- merge(all_lym, y = all_mye)
all_endo_mes <- merge(all_endo, y = all_mes)
all_celltypes <- merge(all_epi, y = list(all_imm, all_endo_mes))
all_celltypes <- ScaleData(all_celltypes)

# Rename/recategorize:
all_celltypes$CT_secondpass <- all_celltypes$CT_firstpass_full
# Endothelial:
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Venous (Epi Obj)"] <- "Venous" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Fibroblasts (Endo Obj)"] <- "Capillary" # Confirmed
# Mesenchymal:
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Activated MyoFBs"] <- "Activated Fibrotic FBs" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Peribronchial FBs"] <- "Myofibroblasts" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "HAS1+/PLIN2+ FBs"] <- "Subpleural FBs" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Hypoxic FBs"] <- "Inflammatory FBs" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("Fibroblasts", "Alveolar FBs")] <- "Alveolar FBs" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Lipofibroblasts"] <- "SMCs/Pericytes" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Proliferating Fibroblasts"] <- "Proliferating FBs" # Confirmed
# Epithelial (Airway):
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("Basal (Alv Obj)", "Interferon-high Basal")] <- "Basal" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("Ciliated (Alv Obj)", "Ciliated")] <- "Multiciliated" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("MUC5B+ (Alv Obj)", "MUC5B+")] <- "Goblet" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("Proliferating Basal (Alv Obj)",
                                                                   "Proliferating Differentiating Ciliated (Alv Obj)",
                                                                   "Proliferating SCGB3A2+/SCGB1A1+ (Alv Obj)")] <- "Proliferating Airway" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "SCGB3A2+"] <- "RASC" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "SCGB3A2+/SCGB1A1+"] <- "Secretory" # Confirmed
# Epithelial (Alveolar):
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "AT2 (Air Obj)"] <- "AT2" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Proliferating Transitional AT2"] <- "Proliferating AT2" # Confirmed
# Immune (Lymphoid):
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "T-cells (Myeloid Obj)"] <- "CD8+ T-cells" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Plasma (Myeloid Obj)"] <- "Plasma" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "pDCs (Myeloid Obj)"] <- "pDCs" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("T-cells", "Granulomatous T-cells")] <- "CD4+ T-cells" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("T-cells (Endo Obj)", "NK cells")] <- "NK/NKT" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Proliferating NK cells"] <- "Proliferating NK/NKT" # Confirmed
# Immune (Myeloid)
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Mast (Lymphoid Obj)"] <- "Mast" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("cDC2s", "cDC2 (Air Obj)", "Langerhans DCs")] <- "cDCs" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "CCL22+ cDC1s"] <- "Migratory DCs" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Transitional AT2 (Myeloid Obj)"] <- "Langerhans cells" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Mast (Epi Obj)"] <- "Basophils" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "S100A+ Macrophages"] <- "Neutrophils" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "FCN1+/S100A+ Macrophages"] <- "Monocytes/MDMs" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Interferon-high Macrophages"] <- "Macrophages - IFN-activated" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "FABP4+ Macrophages"] <- "Alveolar Macrophages" # Confirmed
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full == "Macrophages"] <- "Interstitial Macrophages" # Confirmed

# Set low-quality clusters/subclusters to remove
all_celltypes$CT_secondpass[all_celltypes$CT_firstpass_full %in% c("Proliferating Fibroblasts (Air Obj)",
                                                                   "SCGB3A2+ (Alv Obj)",
                                                                   "Transitional AT2 (Air Obj)",
                                                                   "Proliferating Epithelial (Imm Obj)",
                                                                   "Differentiating Ciliated",
                                                                   "cDC1 (Air Obj)")] <- "REMOVE"
all_celltypes <- subset(all_celltypes, subset = CT_secondpass != "REMOVE")

# Set lineage/sublineage variables
all_celltypes$final_lineage <- ""
all_celltypes$final_sublineage <- ""
all_celltypes$final_sublineage[all_celltypes$CT_secondpass %in% c("RASC", "Secretory", "Multiciliated",
                                                                  "PNEC", "Basal", "Goblet", "Proliferating Airway")] <- "Airway"
all_celltypes$final_sublineage[all_celltypes$CT_secondpass %in% c("AT1", "Transitional AT2", "AT2",
                                                                  "Proliferating AT2", "KRT5-/KRT17+")] <- "Alveolar"
all_celltypes$final_lineage[all_celltypes$final_sublineage %in% c("Airway", "Alveolar")] <- "Epithelial"
all_celltypes$final_sublineage[all_celltypes$CT_secondpass %in% c("NK/NKT", "Tregs", "CD4+ T-cells",
                                                                  "CD8+ T-cells", "Proliferating T-cells",
                                                                  "B cells", "Plasma", "pDCs",
                                                                  "Proliferating NK cells", "Proliferating B cells")] <- "Lymphoid"
all_celltypes$final_sublineage[all_celltypes$CT_secondpass %in% c("cDCs", "Migratory DCs", "Interstitial Macrophages",
                                                                  "Alveolar Macrophages", "Mast", "Proliferating Myeloid",
                                                                  "SPP1+ Macrophages", "Neutrophils",
                                                                  "Macrophages - IFN-activated", "Monocytes/MDMs",
                                                                  "Langerhans cells", "Basophils")] <- "Myeloid"
all_celltypes$final_lineage[all_celltypes$final_sublineage %in% c("Lymphoid", "Myeloid")] <- "Immune"
all_celltypes$final_lineage[all_celltypes$CT_secondpass %in% c("Venous", "Capillary", "Lymphatic", "Arteriole")] <- "Endothelial"
all_celltypes$final_sublineage[all_celltypes$final_lineage == "Endothelial"] <- "Endothelial"
all_celltypes$final_lineage[all_celltypes$CT_secondpass %in% c("SMCs/Pericytes", "Alveolar FBs", "Proliferating FBs",
                                                               "Inflammatory FBs", "Activated Fibrotic FBs",
                                                               "Myofibroblasts", "Subpleural FBs", "Adventitial FBs",
                                                               "Mesothelial")] <- "Mesenchymal"
all_celltypes$final_sublineage[all_celltypes$final_lineage == "Mesenchymal"] <- "Mesenchymal"

# Fix sample name and cell IDs
xenium$sample[xenium$sample == "TILD117MF" & xenium$tma == "TMA5"] <- "TILD117MFB"
old_colnames <- colnames(xenium)
new_colnames <- paste0(xenium$sample, "_", xenium$cell_id)

# Save object
# saveRDS(xenium, "/scratch/avannan/late_IPF_spatial/xenium/REVISION_rds_files/all_celltypes_final_RAPIDS_clustered_062024.rds")





