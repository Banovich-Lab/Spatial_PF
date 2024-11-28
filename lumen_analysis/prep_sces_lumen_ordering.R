## prepare SCE objects for lumen ordering analysis

library(Seurat)
library(ggplot2)
library(dplyr)
library(scater)

cell_objects <- readRDS(file = "../../../Dataset/LFST_2022/REVISION_lumen_pca/xenium_final_CT_plus_niches_lumens_071624.rds")

dim(cell_objects)
#sample_name_to_type <- unique(cell_objects@meta.data[,c("sample","sample_type")])
sample_name_to_type_affect <- unique(cell_objects@meta.data[,c("sample","sample_type", "tma",
                                                               "sample_affect","percent_pathology","disease_status")])
write.csv(sample_name_to_type_affect,file = "../../../Dataset/LFST_2022/REVISION_lumen_pca/sample_type_affect_071624.csv")

View(head(cell_objects@meta.data))
colnames(head(cell_objects@meta.data))
cell_lumen_assignment <- cell_objects@meta.data[,c("sample","sample_type",
                                                   "cell_id","lumen_id",
                                                   "final_CT","CNiche",
                                                   "TNiche","percent_pathology","disease_status",
                                                   "sample_affect","final_lineage","tma")] 
dim(cell_lumen_assignment)
cell_lumen_assignment <- cell_lumen_assignment %>% filter(!is.na(lumen_id))
dim(cell_lumen_assignment)

length(unique(cell_lumen_assignment$lumen_id))
dim(cell_lumen_meta)

## based on the input lumen_pca..
filtered_lumen <- unique(cell_lumen_assignment$lumen_id)[!unique(cell_lumen_assignment$lumen_id) %in% rownames(lumen_pca$x)]
length(filtered_lumen)

length(unique(cell_lumen_assignment$lumen_id))
dim(lumen_pca$x)

## keep the 1747 lumen ids 
cell_lumen_assignment <- cell_lumen_assignment[cell_lumen_assignment$lumen_id %in% rownames(lumen_pca$x),]
dim(cell_lumen_assignment)
length(unique(cell_lumen_assignment$lumen_id))

write.csv(cell_lumen_assignment,"../../../Dataset/LFST_2022/REVISION_lumen_pca/cell_lumen_assignment_071624.csv")

### final CT count and prop
cell_lumen_assignment <- read.csv("../../../Dataset/LFST_2022/REVISION_lumen_pca/cell_lumen_assignment_071624.csv")
cell_lumen_assignment$final_CT
dim(cell_lumen_assignment)
length(unique(cell_lumen_assignment$lumen_id))

## get lumen cell type counts
lumen_cell_types <- cell_lumen_assignment %>% group_by(lumen_id, final_CT) %>% 
  summarise(celltype_count=n()) %>% 
  tidyr::pivot_wider(names_from = "final_CT",
                     values_from = celltype_count) 

lumen_cell_types[is.na(lumen_cell_types)] <- 0 

dim(lumen_cell_types)
#[1] 1747   47

lumen_cell_types_count_m <- lumen_cell_types[,2:ncol(lumen_cell_types)]
lumen_cell_types_count_m <- as.matrix(lumen_cell_types_count_m)
## 46 distinct cell types 
dim(lumen_cell_types_count_m)

rownames(lumen_cell_types_count_m) <- lumen_cell_types$lumen_id
lumen_cell_types_sce <- SingleCellExperiment(assays= list("final_CT_count" = t(lumen_cell_types_count_m)))
dim(lumen_cell_types_sce)

colData(lumen_cell_types_sce) <- cbind(colData(lumen_cell_types_sce),
                                       cell_lumen_meta[match(colnames(lumen_cell_types_sce),
                                                             cell_lumen_meta$lumen_id),])
rownames(sample_name_to_type_affect) <- sample_name_to_type_affect$sample

lumen_cell_types_sce$sample_type <- 
  sample_name_to_type_affect[match(lumen_cell_types_sce$sample,sample_name_to_type_affect$sample),"sample_type"]
lumen_cell_types_sce$percent_pathology <- 
  sample_name_to_type_affect[match(lumen_cell_types_sce$sample,sample_name_to_type_affect$sample),"percent_pathology"]
lumen_cell_types_sce$disease_status <- 
  sample_name_to_type_affect[match(lumen_cell_types_sce$sample,sample_name_to_type_affect$sample),"disease_status"]
lumen_cell_types_sce$tma <- 
  sample_name_to_type_affect[match(lumen_cell_types_sce$sample,sample_name_to_type_affect$sample),"tma"]


lumen_cell_types_sce <- addPerCellQC(lumen_cell_types_sce, 
                                     assay.type = "final_CT_count")

dim(lumen_cell_types_sce)

plotColData(lumen_cell_types_sce,x="sample",
            y=c("total"),other_fields = "sample_type")+
  facet_wrap(.~sample_type,scale="free_x")+
  plotColData(lumen_cell_types_sce,x="sample",y=c("detected"),
              other_fields = "sample_type")+
  facet_wrap(.~sample_type,scale="free_x")

dim(lumen_cell_types_sce)
## add cell type prop assay
assay(lumen_cell_types_sce,"final_CT_prop") <- apply(lumen_cell_types_sce@assays@data$final_CT_count,
                                                      2,function(x){x/sum(x)})
max(assay(lumen_cell_types_sce,"final_CT_prop"))
min(assay(lumen_cell_types_sce,"final_CT_prop"))
assay(lumen_cell_types_sce,"final_CT_logit") <- apply(lumen_cell_types_sce@assays@data$final_CT_prop,
                                                            2,function(x){ 
                                                              res  = log(x/(1-x))
                                                              res[is.infinite(res)] <- log(1e-6)
                                                              res
                                                            })
assay(lumen_cell_types_sce,"final_CT_prop_asin") <- apply(lumen_cell_types_sce@assays@data$final_CT_prop,
                                                      2,function(x){ 
                                                        res = asin(sqrt(x))
                                                        res
                                                      })


table(colSums(assay(lumen_cell_types_sce,"final_CT_prop")))
final_CT_prop_offset <- (lumen_cell_types_sce@assays@data$final_CT_count + 0.5)

final_CT_prop_offset <- apply(final_CT_prop_offset,
                               2,function(x){x/sum(x)})
final_CT_prop_offset_logit <- apply(final_CT_prop_offset,
                                     2,function(x){log(x) - log(1-x)})

assay(lumen_cell_types_sce,"final_CT_prop_offset_logit") <- final_CT_prop_offset_logit

##
#dir.create("output/savedRDS/revisions2024/")
saveRDS(lumen_cell_types_sce,file="output/savedRDS/revisions2024/lumen_cell_types_sce.all.ct.rds")

## remove low count cell type 
keep_celltypes <- colSums(lumen_cell_types_count_m >=3 ) >= 20

sum(keep_celltypes)
(keep_celltypes[!keep_celltypes])
## Mesothelial were removed

lumen_cell_types_filter_sce <- lumen_cell_types_sce[names(keep_celltypes[keep_celltypes]),]

##

assay(lumen_cell_types_filter_sce,"final_CT_prop") <- apply(lumen_cell_types_filter_sce@assays@data$final_CT_count,
                                                     2,function(x){x/sum(x)})

max(assay(lumen_cell_types_filter_sce,"final_CT_prop"))
min(assay(lumen_cell_types_filter_sce,"final_CT_prop"))
assay(lumen_cell_types_filter_sce,"final_CT_logit") <- apply(lumen_cell_types_filter_sce@assays@data$final_CT_prop,
                                                      2,function(x){ 
                                                        res  = log(x/(1-x))
                                                        res[is.infinite(res)] <- log(1e-6)
                                                        res
                                                      })
assay(lumen_cell_types_filter_sce,"final_CT_prop_asin") <- apply(lumen_cell_types_filter_sce@assays@data$final_CT_prop,
                                                             2,function(x){ 
                                                               res = asin(sqrt(x))
                                                               res
                                                             })
table(colSums(assay(lumen_cell_types_filter_sce,"final_CT_prop")))

final_CT_prop_offset <- (lumen_cell_types_filter_sce@assays@data$final_CT_count + 0.5)

final_CT_prop_offset <- apply(final_CT_prop_offset,
                              2,function(x){x/sum(x)})

final_CT_prop_offset_logit <- apply(final_CT_prop_offset,
                                    2,function(x){log(x) - log(1-x)})

assay(lumen_cell_types_filter_sce,"final_CT_prop_offset_logit") <- final_CT_prop_offset_logit

saveRDS(lumen_cell_types_filter_sce,file="output/savedRDS/revisions2024/lumen_cell_types_filter_sce.rds")


### cell niche 12 count and prop

lumen_cellniche <- cell_lumen_assignment %>% group_by(lumen_id, CNiche) %>% 
  summarise(cellniche_count=n()) %>% 
  tidyr::pivot_wider(names_from = "CNiche",values_from = cellniche_count)

lumen_cellniche[is.na(lumen_cellniche)] <- 0
dim(lumen_cellniche)

lumen_cellniche_count_m <- lumen_cellniche[,2:ncol(lumen_cellniche)]
lumen_cellniche_count_m <- as.matrix(lumen_cellniche_count_m)

dim(lumen_cellniche_count_m)
rownames(lumen_cellniche_count_m) <- lumen_cellniche$lumen_id

lumen_cellniche_sce <- SingleCellExperiment(assays= list("CNiche_count" = t(lumen_cellniche_count_m)))
dim(lumen_cellniche_sce)

colData(lumen_cellniche_sce) <- cbind(colData(lumen_cellniche_sce),
                                      colData(lumen_cell_types_sce[,colnames(lumen_cellniche_sce)])[,1:49])


dim(lumen_cellniche_sce)
## add prop assay
assay(lumen_cellniche_sce,"CNiche_prop") <- apply(lumen_cellniche_sce@assays@data$CNiche_count,
                                                         2,function(x){x/sum(x)})
## checks
## sum(assay(lumen_cellniche_sce,"CNiche_prop")["C12",] == lumen_cellniche_sce$C12)
max(lumen_cellniche_sce@assays@data$CNiche_prop)
min(lumen_cellniche_sce@assays@data$CNiche_prop)
table(colSums(assay(lumen_cellniche_sce,"CNiche_prop")))

assay(lumen_cellniche_sce,"CNiche_prop_logit") <- apply(lumen_cellniche_sce@assays@data$CNiche_prop,
                                                               2,function(x){ 
                                                                 res  = log(x/(1-x))
                                                                 res[is.infinite(res)] <- log(1e-6)
                                                                 res
                                                               })


CNiche_prop_offset <- (lumen_cellniche_sce@assays@data$CNiche_count + 0.5)

CNiche_prop_offset <- apply(CNiche_prop_offset,
                                   2,function(x){x/sum(x)})
CNiche_prop_offset_logit <- apply(CNiche_prop_offset,
                                         2,function(x){log(x) - log(1-x)})

assay(lumen_cellniche_sce,"CNiche_prop_offset_logit") <- CNiche_prop_offset_logit
#lumen_cellniche_sce$total

saveRDS(lumen_cellniche_sce,file = "output/savedRDS/revisions2024/lumen_cellniche_sce.rds")

## transcript based niches

### trans niche 12 (5k,dmax3.0) count and prop
lumen_tranx_niche <- cell_lumen_assignment %>% group_by(lumen_id, TNiche) %>% 
  summarise(tniche_count=n()) %>% 
  tidyr::pivot_wider(names_from = "TNiche",values_from = tniche_count) %>%
  mutate_all( ~coalesce(.,0))

dim(lumen_tranx_niche)

lumen_tranx_niche_count_m <- lumen_tranx_niche[,2:ncol(lumen_tranx_niche)]
lumen_tranx_niche_count_m <- as.matrix(lumen_tranx_niche_count_m)
dim(lumen_tranx_niche_count_m)
rownames(lumen_tranx_niche_count_m) <- lumen_tranx_niche$lumen_id

lumen_tranx_niche_sce <- SingleCellExperiment(assays = 
                                                list("tranx_nichek12_count" = t(lumen_tranx_niche_count_m)))
dim(lumen_tranx_niche_sce)
colData(lumen_tranx_niche_sce) <- cbind(colData(lumen_tranx_niche_sce),
                                        colData(lumen_cell_types_sce[,colnames(lumen_tranx_niche_sce)])[,1:49])

lumen_tranx_niche_sce <- addPerCellQC(lumen_tranx_niche_sce, 
                                      assay.type = "tranx_nichek12_count")
dim(lumen_tranx_niche_sce)
## add prop assay
assay(lumen_tranx_niche_sce,"tranx_nichek12_prop") <- apply(lumen_tranx_niche_sce@assays@data$tranx_nichek12_count,
                                                            2,function(x){x/sum(x)})
max(assay(lumen_tranx_niche_sce,"tranx_nichek12_prop"))
min(assay(lumen_tranx_niche_sce,"tranx_nichek12_prop"))
assay(lumen_tranx_niche_sce,"tranx_nichek12_prop_logit") <- apply(lumen_tranx_niche_sce@assays@data$tranx_nichek12_prop,
                                                                  2,function(x){ 
                                                                    res  = log(x/(1-x))
                                                                    res[is.infinite(res)] <- log(1e-6)
                                                                    res
                                                                  })
table(colSums(assay(lumen_tranx_niche_sce,"tranx_nichek12_prop")))
tranx_nichek12_prop_offset <- (lumen_tranx_niche_sce@assays@data$tranx_nichek12_count + 0.5)
tranx_nichek12_prop_offset <- apply(tranx_nichek12_prop_offset,
                                    2,function(x){x/sum(x)})
tranx_nichek12_prop_offset_logit <- apply(tranx_nichek12_prop_offset,
                                          2,function(x){log(x) - log(1-x)})

assay(lumen_tranx_niche_sce,"tranx_nichek12_prop_offset_logit") <- tranx_nichek12_prop_offset_logit
rowSums(lumen_tranx_niche_sce@assays@data$tranx_nichek12_count > 1)
saveRDS(lumen_tranx_niche_sce,file = "output/savedRDS/revisions2024/lumen_tranx_niche_sce.rds")


## filter low niche count

keep_cniche <- rowSums(lumen_cellniche_sce@assays@data$CNiche_count>0) > 10
keep_cniche
sum(keep_cniche)


keep_tniche <- rowSums(lumen_tranx_niche_sce@assays@data$tranx_nichek12_count>0) > 10
keep_tniche

## prepare sce for gene expression based psedudomtime analysis
## add trans count assay
dim(cell_lumen_assignment)
dim(cell_objects)
nuclei_RNA <-  t(as.matrix(cell_objects@assays$RNA@counts[,rownames(cell_lumen_assignment)]))
nuclei_RNA[1:5,1:5]
lumen_nuclei_RNA <- t(sapply(by(nuclei_RNA,cell_lumen_assignment$lumen_id,colSums),identity))
dim(lumen_nuclei_RNA)
lumen_nuclei_RNA[colnames(lumen_cell_types_sce),]

lumen_nuclei_sce <- SingleCellExperiment(assays = 
                                           list("trans_count" =
                                                  t(lumen_nuclei_RNA[colnames(lumen_cell_types_sce),])))

## add trans 
# trans prop 
assay(lumen_nuclei_sce,"trans_prop") <- t(apply(assay(lumen_nuclei_sce,"trans_count"),
                                                1,function(x){x/sum(x)}))
assay(lumen_nuclei_sce,"trans_prop_offset") <- t(apply(assay(lumen_nuclei_sce,"trans_count"),
                                                       1,function(x){
                                                         (x + 0.5)/sum(x +0.5)
                                                       }))
table(rowSums(assay(lumen_nuclei_sce,"trans_prop_offset")))
table(rowSums(assay(lumen_nuclei_sce,"trans_prop")))

assay(lumen_nuclei_sce,"trans_prop_offset_logit") <- t(apply(assay(lumen_nuclei_sce,
                                                                   "trans_prop_offset"), 1, 
                                                             function(x){log(x/(1-x))}))

#sum(lumen_nuclei_sce$total == lumen_cell_types_sce$total)
colData(lumen_nuclei_sce) <- colData(lumen_cell_types_sce)[rownames(colData(lumen_nuclei_sce)),1:49]
## lumen size is the number of cells 
lumen_nuclei_sce$lumen_size <- lumen_cell_types_sce$sum

dim(lumen_nuclei_sce)
saveRDS(lumen_nuclei_sce,file = "output/savedRDS/revisions2024/lumen_nuclei_sce.rds")


## per ct 

## add trans count per cell type 
nuclei_RNA <-  t(as.matrix(cell_objects@assays$RNA@counts[,rownames(cell_lumen_assignment)]))
dim(nuclei_RNA)

lumen_nuclei_RNA_per_celltype <- t(sapply(by(nuclei_RNA,paste0(cell_lumen_assignment$final_CT,"_",
                                                               cell_lumen_assignment$lumen_id),colSums),
                                          identity))
dim(lumen_nuclei_RNA_per_celltype)
lumen_nuclei_RNA_per_celltype_sce <- SingleCellExperiment(assays= 
                                                            list("nuclei_trans_count_perCT" =
                                                                   t(lumen_nuclei_RNA_per_celltype)))

lumen_nuclei_RNA_per_celltype_sce$final_CT <- sapply(strsplit(colnames(lumen_nuclei_RNA_per_celltype_sce),
                                                              "_"),`[[`,1)

lumen_nuclei_RNA_per_celltype_sce$sample_name <- sapply(strsplit(colnames(lumen_nuclei_RNA_per_celltype_sce),
                                                                 "_"),`[[`,2)
lumen_nuclei_RNA_per_celltype_sce$lumen_id <- sapply(strsplit(colnames(lumen_nuclei_RNA_per_celltype_sce),
                                                              "_"),`[[`,3)
lumen_nuclei_RNA_per_celltype_sce$lumen_id <- paste0(lumen_nuclei_RNA_per_celltype_sce$sample_name,
                                                     "_",
                                                     lumen_nuclei_RNA_per_celltype_sce$lumen_id )

colnames(lumen_nuclei_RNA_per_celltype_sce)
n_cells <- cell_lumen_assignment %>% group_by(lumen_id) %>% summarise(n_cells = length(unique(final_CT)))
sum(n_cells$n_cells)
dim(lumen_nuclei_RNA_per_celltype_sce)
# rm(cell_obj)

## test for these cell types only
# [1] "Activated Fibrotic FBs"      "Alveolar FBs"                "Alveolar Macrophages"        "Arteriole"                  
# [5] "AT1"                         "AT2"                         "B cells"                     "Basal"                      
# [9] "Capillary"                   "CD4+ T-cells"                "CD8+ T-cells"                "cDCs"                       
# [13] "Inflammatory FBs"            "Interstitial Macrophages"    "KRT5-/KRT17+"                "Lymphatic"                  
# [17] "Macrophages - IFN-activated" "Mast"                        "Migratory DCs"               "Monocytes/MDMs"             
# [21] "Multiciliated"               "Myofibroblasts"              "Neutrophils"                 "NK/NKT"                     
# [25] "Plasma"                      "Proliferating AT2"           "Proliferating FBs"           "Proliferating Myeloid"      
# [29] "Proliferating T-cells"       "RASC"                        "SMCs/Pericytes"              "SPP1+ Macrophages"          
# [33] "Transitional AT2"            "Tregs"                       "Venous"  
##  The above cell types appeared in > 100 lumens
(keep_celltypes <- 
    names(table(lumen_nuclei_RNA_per_celltype_sce$final_CT))[table(lumen_nuclei_RNA_per_celltype_sce$final_CT) >
                                                               100])


lumen_nuclei_RNA_per_celltype_sce_filter <- lumen_nuclei_RNA_per_celltype_sce[,lumen_nuclei_RNA_per_celltype_sce$final_CT %in%
                                                                                keep_celltypes]

##
saveRDS(lumen_nuclei_RNA_per_celltype_sce,
        file="output/savedRDS/revisions2024/lumen_nuclei_RNA_per_celltype_sce.rds")

saveRDS(lumen_nuclei_RNA_per_celltype_sce_filter,
        file="output/savedRDS/revisions2024/lumen_nuclei_RNA_per_celltype_sce_filter.rds")

