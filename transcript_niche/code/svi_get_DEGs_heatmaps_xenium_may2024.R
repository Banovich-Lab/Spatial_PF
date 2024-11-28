
### generate hex plot for xenium sample

.libPaths("/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/Rlibs/4.1.2")

suppressPackageStartupMessages({
  library(ggplot2)
  library(speckle)
  library(dplyr)
  library(RColorBrewer)
  library(limma)
  library(ComplexHeatmap)
  library(Seurat)
  })

args <- (commandArgs(trailingOnly = TRUE))

for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

print(k_dir)
print(k)
if(grepl("gmm",k)) {
  k = gsub("gmm","",strsplit(k,"_")[[1]][1])
}
k <- as.integer(k)
print(k)
print(node_meta_dir)
print(out_png_topFDR)
print(out_png_topPropRatio)
print(out_png_topFDR_tmean)
print(out_png_topPropRatio_tmean)
print(out_png_wide_dot)
print(out_trans_counts_rds)
print(top_n)
top_n = as.integer(top_n)
print(out_niche_prop_bar_plot)
print(out_niche_trans_count_plot)
# print(gmm_id_dir)
samples = c("VUILD96LF","VUHD116A","VUILD115","VUILD91LF","VUILD104MF",
            "VUILD48MF","VUILD105LF","VUILD107MF","VUHD116B","VUILD102LF",
            "VUHD095","VUILD102MF","VUILD48LF","VUILD106","VUILD91MF",
            "VUHD113","VUILD104LF","VUILD78MF","THD0011","VUILD96MF",
            "VUHD069","VUILD78LF","TILD117MF","TILD117LF","TILD175",
            "VUILD110","THD0008","VUILD105MF")

TMAs = c("TMA1","TMA1","TMA3","TMA4","TMA2",
         "TMA2","TMA2","TMA1","TMA1","TMA1",
         "TMA2","TMA1","TMA2","TMA3","TMA4",
         "TMA2","TMA2","TMA4","TMA4","TMA1",
         "TMA2","TMA4","TMA4","TMA4","TMA4",
         "TMA3","TMA3","TMA2")

##. remove vuild110, 106 temporarily 
## 
# samples = c("VUILD96LF","VUHD116A","VUILD115","VUILD91LF","VUILD104MF",
#             "VUILD48MF","VUILD105LF","VUILD107MF","VUHD116B","VUILD102LF",
#             "VUHD095","VUILD102MF","VUILD48LF","VUILD91MF",
#             "VUHD113","VUILD104LF","VUILD78MF","THD0011","VUILD96MF",
#             "VUHD069","VUILD78LF","TILD117MF","TILD117LF","TILD175",
#             "THD0008","VUILD105MF")

# TMAs = c("TMA1","TMA1","TMA3","TMA4","TMA2",
#          "TMA2","TMA2","TMA1","TMA1","TMA1",
#          "TMA2","TMA1","TMA2","TMA4",
#          "TMA2","TMA2","TMA4","TMA4","TMA1",
#          "TMA2","TMA4","TMA4","TMA4","TMA4",
#          "TMA3","TMA2")


sample_levels <- c(  "VUHD116A","VUHD116B","VUHD095","VUHD113","VUHD069","THD0011","THD0008",
                     "TILD117LF","VUILD96LF","VUILD91LF", "VUILD78LF", "VUILD105LF","VUILD48LF",
                     "VUILD102LF", "VUILD104LF","VUILD78MF","VUILD96MF", "VUILD48MF",
                     "VUILD104MF","TILD117MF", "VUILD102MF","VUILD107MF","VUILD91MF","VUILD105MF",
                     "VUILD110","TILD175","VUILD115","VUILD106")

getTransformedPropForK <- function(k_dir = "gmm10_xenium_fullpanel_30ktrained_3NB_aug"){
  transformed_prop <- list()
  
  for(s_n in paste0(samples,"_", TMAs,"_node_meta.csv")) {
    
    node_meta <- readr::read_csv(paste0(node_meta_dir,
                                        s_n),num_threads = 4,
                                        col_types = c("transcript_id"="c","cell_id"="c"),
                                 show_col_types = FALSE)
    node_meta <- node_meta[,c("transcript_id" ,"x_location", "y_location",
                              "z_location","feature_name")]
    region_name <- strsplit(s_n ,"_")[[1]][1]
    tma <- strsplit(s_n ,"_")[[1]][2]
    tma <- gsub(".csv","",tma)
    node_meta$region_name <- region_name
    
    gmm_id <- read.table(file =paste0(node_meta_dir,"/",k_dir,"/",
                                      region_name,"_noMeans_",tma,".txt"))
    
    node_meta$gmm <- gmm_id$V1
    ## remove extra overlapping region from VUILD105LF
      ## remove the unresolvable overlapping transcripts in two regions
    # if(region_name == "VUILD105LF"){
    #   remove_trans_coords <- readr::read_csv("../../../Dataset/LFST_2022/Xenium/remove_nuclei_transcripts_091223/remove_coordinates_VUILD105LF.csv")
    #   remove_trans_coords <- remove_trans_coords[1:4,]
    #   node_meta_removes <- node_meta %>% filter(x_location >  as.numeric(remove_trans_coords[1,"X"]),
    #                                         x_location <  as.numeric(remove_trans_coords[2,"X"]),
    #                                         y_location <  as.numeric(remove_trans_coords[4,"Y"]),
    #                                         y_location >  as.numeric(remove_trans_coords[1,"Y"]))
          
    #   print("dim node_meta_removes: ")
    #   print(dim(node_meta_removes))
    #   print(dim(node_meta))
    #   node_meta <- node_meta %>% filter(!transcript_id %in% node_meta_removes$transcript_id)
    #   print("after removing overlapping trans")
    #   print(dim(node_meta))
    # }
    transformed_prop[[s_n]] = speckle::getTransformedProps(clusters = node_meta$feature_name, 
                                                           sample =node_meta$gmm)
  }
  return(transformed_prop)
}

# if(file.exists(out_trans_counts_rds)){
#   transformed_prop_k = readRDS(out_trans_counts_rds)
# } else {
#   transformed_prop_k = getTransformedPropForK(k_dir=k_dir)
  
# }
transformed_prop_k = getTransformedPropForK(k_dir=k_dir)

# kmeans13_colors <- c(
#   "#2f4f4f","#000000","#008000","#4b0082","#ff0000","#ffd700", '#b5bd61',"#00ffff",
#   "#0000ff","#ff69b4","#1e90ff","#ffdab9","#ff00ff")
# kmeans14_colors <- c(
#   "#2f4f4f","#000000","#008000","#4b0082","#ff0000","#ffd700", '#b5bd61',"#00ffff",
#   "#0000ff","#ff69b4","#1e90ff","#ffdab9","#ff00ff","#ff88ff")


# k_colors = DiscretePalette(n = k)
# names(k_colors) = as.character(0:(k-1))
# if(length(k_colors) > 14) {
#   kmeans14_colors = k_colors
# }


k_15_colors <- c("#F0A0FF", "#0075DC", "#993F00", "#4C005C", "#191919",
                 "#005C31", "#2BCE48", 
                 "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", 
                 "#C20088", "#003380", "#FFA405" )
names(k_15_colors) <- paste0("c",0:14)
gmm8_colors <- k_15_colors[paste0("c",0:7)]
gmm8_k <- 0:7
names(gmm8_colors) <- gmm8_k


gmm9_colors <- k_15_colors[c(paste0("c",0:6),"c8","c7")]
gmm9_k <- c(4,5,1,6,2,3,8,0,7)

names(gmm9_colors) <- (gmm9_k)

gmm10_colors <- k_15_colors[c(paste0("c",0:4),"c9","c5","c10","c6","c7")]
gmm10_k <- c(3,7,6,0,1,4,2,5,8,9)
names(gmm10_colors) <- gmm10_k


gmm11_colors <- k_15_colors[c("c0","c12","c2","c3","c4","c9","c5","c10","c6","c7","c11")]
gmm11_k <- c(1,8,6,4,3,7,2,9,10,0,5)
names(gmm11_colors) <- gmm11_k

gmm12_colors <- k_15_colors[c("c0","c12","c2","c3","c4","c9","c5","c1","c10","c6","c7","c11")]
gmm12_k <- c(8,9,5,10,1,0,6,7,2,4,3,11)
names(gmm12_colors) <- gmm12_k


gmm13_colors <- k_15_colors[c("c0","c12","c2","c3","c4","c9","c5","c13","c1","c10","c6","c7","c11")]
gmm13_k <- c(9,4,1,3,7,6,5,8,10,12,2,0,11)
names(gmm13_colors) <- gmm13_k


matched_color_list <- list("gmm8_colors" = gmm8_colors,
                           "gmm9_colors" = gmm9_colors,
                           "gmm10_colors" = gmm10_colors,
                           "gmm11_colors" = gmm11_colors,
                           "gmm12_colors" = gmm12_colors,
                           "gmm13_colors" = gmm13_colors )
#swatch(kmeans13_colors)

# kmeans14_colors <- c(
#   "#2f4f4f","#000000","#008000","#4b0082","#ff0000","#ffd700", '#b5bd61',"#00ffff",
#   "#0000ff","#ff69b4","#1e90ff","#ffdab9","#ff00ff","#ff88ff")
# names(kmeans14_colors) <- 0:13

# k_colors = DiscretePalette(n = length(unique(node_meta$gmm)))

# names(k_colors) = as.character(unique(node_meta$gmm))

# if(length(k_colors) > 14) {
#   kmeans14_colors = k_colors
# }

# k_colors <- matched_color_list[[paste0("gmm",k,"_colors")]]

if(paste0("gmm",k,"_colors") %in% names(matched_color_list)){
  k_colors <- matched_color_list[[paste0("gmm",k,"_colors")]]
} else {
  k_colors = DiscretePalette(n = k)
}


column_annot_colors <- k_colors
print(column_annot_colors)
if(is.null(names(column_annot_colors))){
  names(column_annot_colors) <- as.character(0:(length(k_colors)-1))
}
## cluster prop plot

gene_count_persample_k <- lapply(transformed_prop_k, function(x){x$Counts})

gene_count_samples = list()
for( sample_name in names(gene_count_persample_k)){
  sample_gene_counts = gene_count_persample_k[[sample_name]]
  res = data.frame(t(sample_gene_counts))
  res$trans_niche = res$sample
  res$sample_name <- sample_name
  gene_count_samples[[sample_name]] = res
}
gene_count_samples_df <- do.call(rbind,gene_count_samples)
rownames(gene_count_samples_df) <- NULL

gene_count_samples_df$sample_name_tma <- gsub("_node_meta.csv","",gene_count_samples_df$sample_name)
gene_count_samples_df$sample_name <- sapply(strsplit(gene_count_samples_df$sample_name_tma,"_"),`[[`,1)
gene_count_samples_df$sample_name <- factor(gene_count_samples_df$sample_name,levels = sample_levels)
p1 = gene_count_samples_df %>% group_by(trans_niche,sample_name) %>% summarise(trans_count = sum(Freq)) %>%
  ggplot()+geom_bar(mapping = aes(x = sample_name, y= trans_count,
                                  fill = trans_niche),position = "dodge",stat = "identity")+
  scale_fill_manual(values = k_colors)+facet_wrap(.~sample_name,scales = "free")


p2  = gene_count_samples_df %>% group_by(trans_niche,sample_name) %>% summarise(trans_count = sum(Freq)) %>%
  ggplot()+geom_bar(mapping = aes(x = sample_name, y= trans_count,
                                  fill = trans_niche),position = "fill",stat = "identity")+
  scale_fill_manual(values = k_colors)+theme_bw()+theme(axis.text.x = element_text(angle = 90))
ggsave(p1,filename = out_niche_trans_count_plot,width = 12,height = 10)
ggsave(p2,filename = out_niche_prop_bar_plot,width = 12,height = 10)

##r




saveRDS(transformed_prop_k,file = out_trans_counts_rds)

genes_343 <- rownames(transformed_prop_k$VUILD96MF_TMA1_node_meta.csv$Counts)

prepareProp_list <- function(trans_prop){
  names(trans_prop) <- sapply(strsplit(names(trans_prop),"_"),`[[`,1)
  
  transformed_prop_list = lapply(names(trans_prop), function(x) {
    s = trans_prop[[x]]
    s = s[["TransformedProps"]]
    colnames(s) <- paste0(x,"_",colnames(s))
    s})
  # rownames_removedV <- rownames(transformed_prop_list[[which(sapply(transformed_prop_list,dim)
  #                                                            [1,]==342)]])
  # rownames_removedV <- rownames(transformed_prop_list[[which(sapply(transformed_prop_list,dim)
  #                                                            [1,]==342)]])
  intersect_rownames <- Reduce(intersect, sapply(transformed_prop_list,rownames))
  rownames_removedV <- intersect_rownames
  message(length(rownames_removedV))

  message(length(rownames_removedV))
  message("missing gene ", setdiff(genes_343, rownames_removedV))
  message("from sample ", names(which(sapply(transformed_prop_list,dim)[1,]==342)))
  
  transformed_prop_list = lapply(names(trans_prop), function(x) {
    s = trans_prop[[x]]
    s = s[["TransformedProps"]]
    colnames(s) <- paste0(x,"_",colnames(s))
    s[rownames_removedV,]
  })
  counts_list = lapply(names(trans_prop), function(x) {
    s = trans_prop[[x]]
    s = s[["Counts"]]
    colnames(s) <- paste0(x,"_",colnames(s))
    s[rownames_removedV,]})
  prop_list = lapply(names(trans_prop), function(x) {
    s = trans_prop[[x]]
    s = s[["Proportions"]]
    colnames(s) <- paste0(x,"_",colnames(s))
    s[rownames_removedV,]})
  
  prop_df = do.call(cbind, transformed_prop_list)
  prop_df_list <- list()
  prop_df_list$TransformedProps <- prop_df
  prop_df_list$Proportions <- do.call(cbind,prop_list)
  prop_df_list$Counts <- do.call(cbind,counts_list)
  
  return(prop_df_list)
}

prop_df_list_k = prepareProp_list(trans_prop = transformed_prop_k)

prop_cut <- 5e-4 

remove_sample_cluster_k <- data.frame(mol_counts = colSums(prop_df_list_k$Counts),
                                       cluster_sample = colnames(prop_df_list_k$Counts) ) %>%
  mutate(cluster_id = sapply(strsplit(cluster_sample,
                                      "_"),`[[`,2),
         sample_name = sapply(strsplit(cluster_sample,
                                       "_"),`[[`,1)) %>% 
  group_by(sample_name) %>%
  mutate(mol_sample_total = sum(mol_counts)) %>%
  mutate(mol_prop_insample = mol_counts/mol_sample_total) %>%
  filter(mol_prop_insample < prop_cut)
remove_sample_cluster_k

getDEGResults_otherGroups <- function(clusters=0:10,prop_list){
  test_result_oneVsAll = list()
  
  for(cid in clusters) {
    cluster_group <- sapply(strsplit(colnames(prop_list$TransformedProps),"_"),`[[`, 2)
    # cluster_group <- c("test","other")[ifelse(cluster_group == c,1,2)]
    sample_name <- sapply(strsplit(colnames(prop_list$TransformedProps),"_"),`[[`, 1)
    message("test cluster ",cid)
    
    design <- model.matrix(~ 0 + cluster_group + sample_name)
    #  design <- model.matrix(~ 0 + cluster_group + sample_type)
    contr1 = paste0("cluster_group",cid)
    message("contr1 ", contr1, " \n")
    contr2 = paste0("cluster_group",setdiff(clusters,cid),collapse = "+")
    contr2 = paste0("(",contr2,")/",length(setdiff(clusters,cid)))
    #    message("contr2", contr2)
    mycontr <- makeContrasts(contrasts = paste0(contr1," - ",contr2), levels=design)
    message("mycontr", colnames(mycontr), "\n")
    
    result_test <- propeller.ttest(prop_list, design,
                                   contrasts = mycontr, robust=TRUE,
                                   trend=FALSE, 
                                   sort=TRUE)
    
    result_test$gene <- rownames(result_test)
    result_test$test_cluster <- cid
    rownames(result_test) <- NULL
    
    test_result_oneVsAll[[paste0("c",cid)]] = result_test
  }
  
  one_VS_all <- do.call(rbind, test_result_oneVsAll)
  return(one_VS_all)
}

filter_prop_list <- function(prop_list, rm_sample_cluster){
  prop_list$TransformedProps <- prop_list$TransformedProps[,!colnames(prop_list$TransformedProps) %in% rm_sample_cluster$cluster_sample]
  prop_list$Counts <- prop_list$Counts[,!colnames(prop_list$Counts) %in% rm_sample_cluster$cluster_sample]
  prop_list$Proportions <-  prop_list$Proportions[,!colnames(prop_list$Proportions) %in% 
                                                    rm_sample_cluster$cluster_sample]
  prop_list
}

prop_df_list_filter_k <- filter_prop_list(prop_list = prop_df_list_k, remove_sample_cluster_k)


dim(prop_df_list_filter_k$TransformedProps)

one_VS_all_varg_filter_k <- getDEGResults_otherGroups(0:(k-1),prop_list = prop_df_list_filter_k)

gene_name = data.frame(gene = sort(unique(one_VS_all_varg_filter_k$gene)))
gene_name$group = cut(1:length(unique(gene_name$gene)),4,labels = FALSE)


plot_dot_wide = one_VS_all_varg_filter_k %>% filter(FDR<=0.01) %>% 
  mutate(gene_group = gene_name$group[match(.data$gene, gene_name$gene)]) %>% 
  mutate(test_cluster = factor(as.character(test_cluster),
                               levels = as.character(0:11))) %>%
  ggplot() + geom_point(mapping = aes(x = test_cluster, y = gene,
                                      color  = log2(PropRatio),
                                      size = -log10(FDR))) + 
  scale_color_gradient2(low = "blue",high = "red") + theme_bw() +
  facet_wrap(.~gene_group,scales = "free_y",ncol = 4)+
  theme(  strip.background = element_blank(),
          strip.text.x = element_blank())+xlab("One versus all")
ggsave(plot_dot_wide,filename = out_png_wide_dot,width = 16,height = 10)

#celltypemarkers <- readr::read_csv(file = "../../Datasets/Lung/cell_Seurat_obj/Xenium/May15_obj/XenhLungGeneCT_new.csv")

#Top 5 smallest FDR (up / dn)
# one_VS_all_varg_filter_k %>% mutate(direction = PropRatio > 1) %>%
#   group_by(direction,test_cluster) %>% slice_min(n = top_n, order_by = FDR) %>% 
#   mutate(gene_group = celltypemarkers$Lineage[match(.data$gene, celltypemarkers$Gene)]) %>% 
#   mutate(test_cluster = factor(as.character(test_cluster),
#                                levels = as.character(0:11))) %>%
#   ggplot() + geom_point(mapping = aes(x = test_cluster, y = gene,
#                                       color  = log2(PropRatio),
#                                       size = -log10(FDR))) + 
#   scale_color_gradient2(low = "blue",high = "red") + theme_bw() +
#   theme(  strip.background = element_blank(),
#           strip.text.y = element_text(angle = 0))+xlab("One versus all")

## Select top FDR 
selected_genes_fdr <- one_VS_all_varg_filter_k %>% mutate(direction = PropRatio > 1) %>% 
  group_by(direction,test_cluster) %>% slice_min(n = top_n, order_by = FDR) %>% select("gene")

## Select  FDR < 0.01 and top PropRatio
selected_genes_propRatio <- one_VS_all_varg_filter_k %>% mutate(direction = PropRatio > 1) %>% 
  filter(FDR < 0.01) %>%
  group_by(direction,test_cluster) %>% slice_max(n = top_n, order_by = abs(log2(PropRatio))) %>% 
  select("gene")

row_orders = order(sapply(strsplit(colnames(prop_df_list_filter_k$Proportions),"_"),`[[`, 2))
row_cid = sapply(strsplit(colnames(prop_df_list_filter_k$Proportions),"_"),`[[`, 2)
#ct_group <- celltypemarkers$Lineage[match(unique(selected_genes_fdr$gene), celltypemarkers$Gene)]
# 
# ComplexHeatmap::Heatmap(t(prop_df_list_filter_k$TransformedProps[unique(selected_genes_fdr$gene),]),column_names_gp = gpar(fontsize=8),
#                         cluster_rows = F,row_order = row_orders , row_names_gp = gpar(fontsize = 5),
#                         name = "Gene Proportion \n (logit transformed)")
# ComplexHeatmap::Heatmap(t(prop_df_list_filter_k$TransformedProps[unique(selected_genes_fdr$gene),]),
#                         cluster_rows = F,cluster_columns = F,row_order = row_orders , 
#                         row_names_gp = gpar(fontsize = 5),
#                         name = "Gene Proportion \n (logit transformed)", row_split = row_cid,column_names_gp = gpar(fontsize=8),
#                         column_split = ct_group)
# 
# ComplexHeatmap::Heatmap(t(prop_df_list_filter_k$TransformedProps[unique(selected_genes_fdr$gene),]),
#                         cluster_rows = F,cluster_columns = F,row_order = row_orders , 
#                         row_names_gp = gpar(fontsize = 5),
#                         name = "Gene Proportion \n (logit transformed)", row_split = row_cid,
#                         column_split = ct_group, show_row_names = FALSE)

# ComplexHeatmap::Heatmap(t(prop_df_list_filter_k$TransformedProps[unique(selected_genes_fdr$gene),]),
#                         cluster_rows = F,cluster_columns = T,row_order = row_orders , 
#                         row_names_gp = gpar(fontsize = 5),show_column_dend = FALSE,
#                         name = "Gene Proportion \n (logit transformed)", row_split = row_cid,
#                         show_row_names = FALSE,column_names_gp = gpar(fontsize=8))


ha = HeatmapAnnotation(
  clusters = row_cid,
  border = FALSE,
  col = list(clusters =column_annot_colors),
  which = "column",
  name = "Clusters"
)

ht_k_rowgenes_topFDR <- ComplexHeatmap::Heatmap((prop_df_list_filter_k$TransformedProps[unique(selected_genes_fdr$gene),]),
                                                  cluster_columns   = F,cluster_rows = T,column_order = row_orders , 
                                                  column_names_gp = gpar(fontsize = 5),show_row_dend = FALSE,
                                                  name = "Gene Proportion \n (logit transformed)",column_split = row_cid,
                                                  show_column_names = FALSE, row_names_gp = gpar(fontsize=8),top_annotation = ha)

ht_k_rowgenes_topProp <- ComplexHeatmap::Heatmap((prop_df_list_filter_k$TransformedProps[unique(selected_genes_propRatio$gene),]),
                                                   cluster_columns   = F,cluster_rows = T,column_order = row_orders , 
                                                   column_names_gp = gpar(fontsize = 5),show_row_dend = FALSE,
                                                   name = "Gene Proportion \n (logit transformed)",column_split = row_cid,
                                                   show_column_names = FALSE, row_names_gp = gpar(fontsize=8),top_annotation = ha)
png(filename = out_png_topFDR,width = 10,height = 16,units = "in",res=100)
draw(ht_k_rowgenes_topFDR)
dev.off()
png(filename = out_png_topPropRatio,width = 10,height = 16,units = "in",res=100)

draw(ht_k_rowgenes_topProp)
dev.off()

#### Pool samples per cluster, Select genes by Rank of FDR

selected_genes_fdr <- one_VS_all_varg_filter_k %>% mutate(direction = PropRatio > 1) %>% 
  group_by(direction,test_cluster) %>% slice_min(n = top_n, order_by = FDR) %>% select("gene")



mean_prop <- (one_VS_all_varg_filter_k[match(unique(selected_genes_fdr$gene),
                                               one_VS_all_varg_filter_k$gene),c(1:k)])
rownames(mean_prop) <- unique(selected_genes_fdr$gene)
mean_prop <- apply(mean_prop, 2, function(x)(log(x/(1-x))))
colnames(mean_prop) <- gsub("cluster_group","T",gsub("PropMean\\.","",colnames(mean_prop)))

png(out_png_topFDR_tmean,width = 6,height = 16,units = "in",res=100)
Heatmap(mean_prop,cluster_columns = FALSE,cluster_rows = FALSE,name = "Mean gene proportion \n (logit transformed)")
dev.off()
#Heatmap(mean_prop,cluster_columns = TRUE, name = "Mean gene proportion \n (logit transformed)")
#### Pool samples per cluster, Select genes by change of PropRatio 

selected_genes_prop <- one_VS_all_varg_filter_k %>% mutate(direction = PropRatio > 1) %>% filter(FDR<0.01)%>%
  group_by(direction,test_cluster) %>% slice_max(n = top_n, order_by = abs(log2(PropRatio))) %>% select("gene")

mean_prop <- (one_VS_all_varg_filter_k[match(unique(selected_genes_prop$gene),
                                               one_VS_all_varg_filter_k$gene),c(1:k)])
rownames(mean_prop) <- unique(selected_genes_prop$gene)
#mean_prop <- apply(mean_prop, 2, function(x)(log(x/(1-x))))
mean_prop <- apply(mean_prop, 2, function(x){
  res = log(x/(1-x))
  res[is.infinite(res)] <- log(1e-6)
  res
})

colnames(mean_prop) <- gsub("cluster_group","T",gsub("PropMean\\.","",colnames(mean_prop)))

png(out_png_topPropRatio_tmean,width = 6,height = 16,units = "in",res=100)
Heatmap(mean_prop,cluster_columns = FALSE,cluster_rows = FALSE, name = "Mean gene proportion \n (logit transformed)")
dev.off()
