# 2023-11-13_celltype_proprotion_gene_props_regression_pathology_scores
# simplified script for plots
# from 2023-11-13_celltype_proprotion_gene_props_regression_pathology_scores.html
library(ggplot2)
library(dplyr)

# lm result after fitting limma::lmFit,ebayes to logit transformed cell type proportion data
fit_res_logres_all <-  readr::read_csv(file = "output/savedFiles/celltype_prop_all_changes_adjpathology_scores.csv")
fit_res_logres_all %>% ggplot()+geom_point(mapping = aes(x= logFC, y = -log10(adj.P.Val),color= adj.P.Val < 0.01 ))+
  geom_label(mapping = aes(x= logFC, y = -log10(adj.P.Val),label= celltype,color= adj.P.Val < 0.01 ))+
  theme_classic(base_size = 20) +
  theme_classic(base_size = 20) + xlab("Rate of change (logit space)")

fit_res_logres_all %>% ggplot()+
  ggrepel::geom_label_repel(mapping = aes(x= logFC, y = -log10(adj.P.Val),label= celltype,color= adj.P.Val < 0.01 ))+
  geom_vline(aes(xintercept = 0),linetype=3)+
  theme_classic(base_size = 20) +
  theme_classic(base_size = 20) + xlab("Rate of change (logit space)")



# lm result aftering fitting limma lmFit,ebayes to logit transformed gene proportion per sample
all_results <- readr::read_csv( file ="output/savedFiles/all_trans_prop_changes_adjpathology_scores.csv")
prep_plot_df <- all_results %>% mutate(sig = ifelse( adj.P.Val < 0.01 , "sig","nonsig")) %>% 
  mutate(sig_status = ifelse(logFC < 0,"dn","up")) %>% 
  mutate(sig_status = paste0(sig, sig_status)) 

p_pathscore_degs <- prep_plot_df %>%
  ggplot()+
  geom_point(mapping = aes(x= logFC, y = -log10(adj.P.Val),color= sig_status))+
  ggrepel::geom_text_repel(data = prep_plot_df[prep_plot_df$adj.P.Val<0.01,],
                            mapping = aes(x= logFC, y = -log10(adj.P.Val),
                                          label= gene,color = sig_status))+
  geom_vline(aes(xintercept = 0),linetype=3)+
  theme_classic(base_size = 11) + xlab("Rate of change (logit space)")+
  geom_hline(yintercept=2,linetype=2)+geom_vline(xintercept=0,linetype=2)+
  scale_color_manual(values=c("nonsigdn"="grey",
                              "nonsigup"="grey",
                              "sigup"="firebrick4",
                              "sigdn"="#107DAC"))+
  ylab("-log10(FDR)")+
  theme(legend.position = "None")
  #       text = element_text(family = "Courier"))

# p_pathscore_degs2 + p_pathscore_degs
ggsave(filename = "output/figures/manuscript/Dec8/pathology_scores_degs.png",p_pathscore_degs,
       width = 8,height = 8,dpi = 300)
ggsave(filename = "output/figures/manuscript/Dec8/pathology_scores_degs.pdf",p_pathscore_degs,
       device = "pdf",
       width = 8,height = 8,dpi = 300)

# celltype aware pseudobulk aggregation and then testing for DEGs along pathology scores
# results after removing low count genes per celltype
cell_niche_df <- readr::read_csv(file="output/savedFiles/all_trans_prop_changes_per_ct_adjpathology_scores.csv")
cell_niche_df$celltype[cell_niche_df$celltype == "MUC5B+ Secretory"] = "MUC5B+"
cell_niche_df$celltype[cell_niche_df$celltype == "Monocyte-derived Macrophages"] = "Interstitial Macrophages"
# fix the switching of aCap,gCap labels
gCap_cells <- which(cell_niche_df$celltype == "gCap")
aCap_cells <- which(cell_niche_df$celltype == "aCap")
cell_niche_df$celltype[aCap_cells] = "gCap"
cell_niche_df$celltype[gCap_cells] = "aCap"

write.csv(cell_niche_df, "output/savedFiles/DEGS_along_pathology_scores_per_celltype.csv",row.names = FALSE)

cell_niche_df <- cell_niche_df %>% mutate(sig = ifelse( adj.P.Val < 0.01 , "sig","nonsig")) %>% 
  mutate(sig_status = ifelse(logFC < 0,"dn","up")) %>% 
  mutate(sig_status = paste0(sig, sig_status)) 

p_pathscore_degs_perct <- cell_niche_df %>% ggplot()+geom_point(mapping = aes(x= logFC, y = -log10(adj.P.Val),
                                                                              color= sig_status))+
  ggrepel::geom_text_repel(data = cell_niche_df[cell_niche_df$adj.P.Val<0.01,],
                           mapping = aes(x= logFC, y = -log10(adj.P.Val),label = gene,
                                         color = sig_status), max.overlaps =20,
                           size=3,alpha=1)+
  facet_wrap(.~celltype,ncol = 7)+
  geom_vline(aes(xintercept = 0),linetype=3)+
  theme_classic(base_size = 7) + xlab("Rate of change (logit space)")+
  geom_hline(yintercept=2,linetype=3)+
  scale_color_manual(values=c("nonsigdn"="grey",
                              "nonsigup"="grey",
                              "sigup"="firebrick4",
                              "sigdn"="#107DAC"))+
  ylab("-log10(FDR)")+
  theme(legend.position = "None",
        strip.text = element_text(size=8),
        strip.background = element_blank())

p_pathscore_degs_perct <- cell_niche_df %>% ggplot()+geom_point(mapping = aes(x= logFC, y = -log10(adj.P.Val),
                                                    color= sig_status))+
  ggrepel::geom_text_repel(data = cell_niche_df[cell_niche_df$adj.P.Val<0.01,],
                            mapping = aes(x= logFC, y = -log10(adj.P.Val),label = gene,
                                          color = sig_status), max.overlaps = 8,
                           size=3,alpha=1)+
  facet_wrap(.~celltype,ncol = 7,scales = "free_x")+
  geom_vline(aes(xintercept = 0),linetype=3)+
  theme_classic(base_size = 7) + xlab("Rate of change (logit space)")+
  geom_hline(yintercept=2,linetype=3)+
  scale_color_manual(values=c("nonsigdn"="grey",
                              "nonsigup"="grey",
                              "sigup"="firebrick4",
                              "sigdn"="#107DAC"))+
  ylab("-log10(FDR)")+
  theme(legend.position = "None",
        strip.text = element_text(size=8),
        strip.background = element_blank())

ggsave(filename = "output/figures/manuscript/Dec8/pathology_scores_degs_per_ct.png",p_pathscore_degs_perct,
       width = 14,height = 12,dpi = 300)
ggsave(filename = "output/figures/manuscript/Dec8/pathology_scores_degs_per_ct.pdf",p_pathscore_degs_perct,
       device = "pdf",
       width = 14,height = 12,dpi = 300)


### cell niche proportion changes along (adjusted) pathology scores

cell_niche_df <- readr::read_csv(file="output/savedFiles/cellniche_prop_all_changes_adjpathology_scores.csv")
cell_niche_df <- cell_niche_df %>% mutate(sig = ifelse( adj.P.Val < 0.01 , "sig","nonsig")) %>% 
  mutate(sig_status = ifelse(logFC < 0,"dn","up")) %>% 
  mutate(sig_status = paste0(sig, sig_status)) 
cell_niche_df$cell_niche <- paste0("C",cell_niche_df$cell_niche)
write.csv(cell_niche_df,row.names = FALSE,file="output/figures/manuscript/Dec8/pathology_scores_cell_niche.csv")

p_pathscore_cell_niche <- cell_niche_df %>% ggplot()+geom_point(mapping = aes(x= logFC, y = -log10(adj.P.Val),
                                                                              color= sig_status))+
  ggrepel::geom_text_repel(data = cell_niche_df[cell_niche_df$adj.P.Val<0.01,],
                           mapping = aes(x= logFC, y = -log10(adj.P.Val),label = cell_niche,
                                         color = sig_status), max.overlaps =20,
                           size=3,alpha=1)+
  geom_vline(aes(xintercept = 0),linetype=3)+
  theme_classic(base_size = 7) + xlab("Rate of change (logit space)")+
  geom_hline(yintercept=2,linetype=3)+
  scale_color_manual(values=c("nonsigdn"="grey",
                              "nonsigup"="grey",
                              "sigup"="firebrick4",
                              "sigdn"="#107DAC"))+
  ylab("-log10(FDR)")+
  theme(legend.position = "None",
        strip.text = element_text(size=8),
        strip.background = element_blank())


ggsave(filename = "output/figures/manuscript/Dec8/pathology_scores_cell_niche.png",p_pathscore_cell_niche,
       width = 4,height = 4,dpi = 300)
ggsave(filename = "output/figures/manuscript/Dec8/pathology_scores_cell_niche.pdf",p_pathscore_cell_niche,
       device = "pdf",
       width = 4,height = 4,dpi = 300)


### tranx niche proportion changes along (adjusted) pathology scores



tranx_niche_df <- readr::read_csv(file="output/savedFiles/tranxniche_prop_all_changes_adjpathology_scores.csv")
tranx_niche_df <- tranx_niche_df %>% mutate(sig = ifelse( adj.P.Val < 0.01 , "sig","nonsig")) %>% 
  mutate(sig_status = ifelse(logFC < 0,"dn","up")) %>% 
  mutate(sig_status = paste0(sig, sig_status)) 
tranx_niche_df$trans_niche <- as.numeric(tranx_niche_df$trans_niche) + 1
tranx_niche_df$trans_niche <- paste0("T",tranx_niche_df$trans_niche)
write.csv(tranx_niche_df,row.names = FALSE,file="output/figures/manuscript/Dec8/pathology_scores_tranx_niche.csv")

p_pathscore_tranx_niche <- tranx_niche_df %>% ggplot()+geom_point(mapping = aes(x= logFC, y = -log10(adj.P.Val),
                                                                              color= sig_status))+
  ggrepel::geom_text_repel(data = tranx_niche_df[tranx_niche_df$adj.P.Val<0.01,],
                           mapping = aes(x= logFC, y = -log10(adj.P.Val),label = trans_niche,
                                         color = sig_status), max.overlaps =20,
                           size=3,alpha=1)+
  geom_vline(aes(xintercept = 0),linetype=3)+
  theme_classic(base_size = 7) + xlab("Rate of change (logit space)")+
  geom_hline(yintercept=2,linetype=3)+
  scale_color_manual(values=c("nonsigdn"="grey",
                              "nonsigup"="grey",
                              "sigup"="firebrick4",
                              "sigdn"="#107DAC"))+
  ylab("-log10(FDR)")+
  theme(legend.position = "None",
        strip.text = element_text(size=8),
        strip.background = element_blank())
ggsave(filename = "output/figures/manuscript/Dec8/pathology_scores_tranx_niche.png",p_pathscore_tranx_niche,
       width = 4,height = 4,dpi = 300)
ggsave(filename = "output/figures/manuscript/Dec8/pathology_scores_tranx_niche.pdf",p_pathscore_tranx_niche,
       device = "pdf",
       width = 4,height = 4,dpi = 300)

