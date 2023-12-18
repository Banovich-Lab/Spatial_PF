
gene_df <- read_csv("/Users/avannan/Downloads/Gene Categories (Lumen Analysis) - Sheet4.csv") %>%
  mutate(location_annotation = case_when(location_annotation %in% c("ER", "ER_membrane") ~ "ER",
                                         location_annotation %in% c("mitochondrial", "mitochondrial_membrane") ~ "mitochondrial",
                                         TRUE ~ location_annotation),
         specific_annotation = ifelse(specific_annotation == "proease", "protease", specific_annotation),
         trajectory = paste0("G", trajectory),
         new_trajectory = case_when(peak_rank <= 58 ~ "G1",
                                    peak_rank <= 105 ~ "G2",
                                    peak_rank <= 192 ~ "G3",
                                    peak_rank <= 208 ~ "G4",
                                    peak_rank <= 271 ~ "G5")) %>%
  filter(broad_annotation != "uncharacterized")

g1_genes <- gene_df %>% filter(trajectory == "G1") %>% pull(gene)
g2_genes <- gene_df %>% filter(trajectory == "G2") %>% pull(gene)
g3_genes <- gene_df %>% filter(trajectory == "G3") %>% pull(gene)
g4_genes <- gene_df %>% filter(trajectory == "G4") %>% pull(gene)
g5_genes <- gene_df %>% filter(trajectory == "G5") %>% pull(gene)

epithelial_cells <- subset(xenium, subset = lineage == "Epithelial")
DotPlot(epithelial_cells, group.by = "broad_CT5", features = g2_genes) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





#### Dotplot Heatmap  ----




# Pick out genes for heatmap
g2_genes

# Log-normalizing and scaling all features in the RNA assay
# Scaling so that all features can be visualized using the same color scale
epithelial_cells <- ScaleData(epithelial_cells)
xenium <- ScaleData(xenium)

# Create dotplot base
p <- DotPlot(xenium, features = g5_genes, 
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

# Any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis::viridis(20)[c(1,10, 20)])
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis::viridis(20)[c(1,10, 20)])


# Adding annotation colors
library(scales) 
# col_for_plot <- color_list
# col_for_plot_heatmap <- list(cluster = unlist(color_list))
# exp_mat <- exp_mat[, names(color_list)]
# percent_mat <- percent_mat[, names(color_list)]

# Creating a layer to add to plot
layer_fun = function(j, i, x, y, w, h, fill){
  grid.points(x = x, y = y, 
              gp = gpar(col = col_fun(pindex(t(exp_mat), i, j))), # t(exp_mat)
              size = pindex(t(percent_mat), i, j)/100 * unit(5, "mm"), # t(percent_mat)
              pch = 19)
}

# List of legends
lgd_list = list(
  Legend(labels = c(0,0.25,0.5,0.75,1), title = "Proportion     \nExpressing     ", type = "points", 
         pch = 19, size = c(0,0.25,0.5,0.75,1) * unit(5, "mm"),
         legend_gp = gpar(col = "black"), direction = "vertical", ncol = 1,
         title_position = "topcenter", background = NA, grid_height = unit(4.4, "mm")))


# Create heatmap
hp <- Heatmap(t(exp_mat), # t(exp_mat)
              name = "Scaled Expression",
              heatmap_legend_param = list(title = " Scaled   \n Expression   ",
                                          legend_direction = "vertical",
                                          title_position = "topcenter"),
              column_title = " ", 
              row_title = " ",
              row_title_gp = gpar(fontsize = 15),
              column_title_gp = gpar(fontsize = 5),
              row_title_side = "right",
              column_title_side = "top",
              column_names_rot = 45,
              col = col_fun,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun,
              row_names_gp = gpar(fontsize = 10),
              row_names_side = "left",
              column_names_gp = gpar(fontsize = 10),
              cluster_rows = FALSE, 
              row_split = c(rep("Endo", 6), rep("Epi", 11), rep("Imm", 13), rep("Mes", 9)),
              cluster_columns = FALSE,
              column_dend_height = unit(8, "mm"),
              border = "black")

# Add annotations to heatmap
ht <- draw(hp, 
           heatmap_legend_list = lgd_list,
           heatmap_legend_side = "right",
           align_heatmap_legend = "global_center")

# # Add box around KRT5-/KRT17+
# ro <- row_order(ht)
# co <- column_order(ht)
# decorate_heatmap_body("Scaled Expression", row_slice = 1, column_slice = 1, {
#   grid.rect(unit(-0.23, "npc"), unit(0.72, "npc"), # top left
#             width = (length(co[[1]]) + length(co[[2]]))/length(co[[1]]) * unit(0.6105,  "npc") + unit(1, "mm"),
#             height = (length(ro[[1]]) + length(ro[[2]]))/length(ro[[1]]) * unit(0.04, "npc") + unit(0.1, "mm"),
#             gp = gpar(lwd = 1.5, lty = "dashed", fill = NA), just = c("left", "top")
#   )
# })






gene_df %>%
  ggplot(aes(x = as.factor(new_trajectory), fill = broad_annotation)) +
  geom_bar(position = position_fill(), color = "black") +
  scale_fill_manual(values = c("deeppink2", "yellow", "chartreuse3", "grey80")) +
  theme_classic()
gene_df %>%
  ggplot(aes(x = as.factor(new_trajectory), fill = effect_annotation)) +
  geom_bar(position = position_fill(), color = "black") +
  scale_fill_manual(values = c("deeppink2", "chartreuse3", "blue", "grey80")) +
  theme_classic()
gene_df %>%
  ggplot(aes(x = as.factor(new_trajectory), fill = location_annotation)) +
  geom_bar(position = position_fill(), color = "black") +
  # scale_fill_manual(values = c(distinctColorPalette(9), "grey80")) +
  scale_fill_manual(values = c(distinctColorPalette(7), "grey80")) +
  theme_classic()
gene_df %>%
  ggplot(aes(x = as.factor(new_trajectory), fill = type_annotation)) +
  geom_bar(position = position_fill(), color = "black") +
  scale_fill_manual(values = c(distinctColorPalette(13), "grey80")) +
  theme_classic()
gene_df %>%
  ggplot(aes(x = as.factor(new_trajectory), fill = specific_annotation)) +
  geom_bar(position = position_fill(), color = "black") +
  scale_fill_manual(values = c(distinctColorPalette(50), "grey80")) +
  theme_classic()


data_table = table(gene_df$new_trajectory, gene_df$broad_annotation)
data_table
#colSums(data_table)
#rowSums(data_table)
chi_test <- chisq.test(data_table, simulate.p.value = TRUE)

# Extract the standardized residuals
std_residuals <- chi_test$stdres

# Print the standardized residuals
print(std_residuals) 


data_table = table(gene_df$new_trajectory, gene_df$location_annotation)
data_table
#colSums(data_table)
#rowSums(data_table)
chi_test <- chisq.test(data_table, simulate.p.value = TRUE)

# Extract the standardized residuals
std_residuals <- chi_test$stdres

# Print the standardized residuals
print(std_residuals) 
# Cell membrane = enriched G5
# Cytoplasmic = enriched G4, depleted G3
# Nuclear = enriched G2
# Secreted = depleted G2
# Uncharacterized = enriched G1


data_table = table(gene_df$new_trajectory, gene_df$effect_annotation)
data_table
#colSums(data_table)
#rowSums(data_table)
chi_test <- chisq.test(data_table, simulate.p.value = TRUE)

# Extract the standardized residuals
std_residuals <- chi_test$stdres

# Print the standardized residuals
print(std_residuals) 
# Extracellular modulator = depleted G2
# Intracellular modulator = enriched G2
# Transcription factor = enriched G2



data_table = table(gene_df$new_trajectory, gene_df$specific_annotation)
data_table
#colSums(data_table)
#rowSums(data_table)
chi_test <- chisq.test(data_table, simulate.p.value = TRUE)

# Extract the standardized residuals
std_residuals <- chi_test$stdres

# Print the standardized residuals
print(std_residuals) 
# Antimicrobial = enriched G1
# Apolipoprotein = enriched G1
# Channel = enriched G3
# Collagen = enriched G4
# Enzyme = enriched G6
# Intermediate filament = enriched G4
# Isomerase = enriched G6
# Ligand = enriched G1
# Ligase = enriched G5
# Lyase = enriched G1
# Lysozyme = enrichd G3
# Methyltransferase = enrcihed G2
# Receptor = enriched G5
# RNA-binding = enriched G1
# Signaling molecule = enriched G4
# Synthetase = enriched G1
# Transcription factor = enriched G2
# Transferrin = enriched G1
# Translocase = enriched G1



data_table = table(gene_df$new_trajectory, gene_df$type_annotation)
data_table
#colSums(data_table)
#rowSums(data_table)
chi_test <- chisq.test(data_table, simulate.p.value = TRUE)

# Extract the standardized residuals
std_residuals <- chi_test$stdres

# Print the standardized residuals
print(std_residuals) 





