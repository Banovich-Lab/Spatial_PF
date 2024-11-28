.libPaths("/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/Rlibs/4.1.2")
library(scater)
library(ggplot2)
library(mclust)
library(slingshot)
library(tradeSeq)
library(grDevices)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(Seurat)
library(dplyr)
library(BiocParallel)
library(mgcv)

args <- (commandArgs(trailingOnly = TRUE))

for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

print(gam_list)
print(gam_list_gene)
print(out_sig_t1_sterm_png)
print(out_sig_ct_sterm_png)
print(check_plots_dir)
print(sig_t1_partial_effect_resid_dir)
print(sig_ct_partial_effect_resid_dir)
print(sig_t1_partial_effect_resid_gg_dir)
print(sig_ct_partial_effect_resid_gg_dir)

gam_list <- readRDS(gam_list)
gam_list_gene <- read.table(gam_list_gene)
gam_list_gene
## nterms

gam_list_stables <- lapply(gam_list, 
                           function(md){
                             res <- summary(md)$s.table
                             terms <- rownames(res)
                             res <- data.frame(res)
                             res$term <- terms
                             res
                           })
nterms <- nrow(gam_list_stables[[1]])
gam_list_stables <- do.call(rbind,gam_list_stables)
gam_list_stables_genes <- lapply(gam_list, 
                                 function(md){
                                  md$gene
                                 })
#gam_list_stables$gene <- rep(gam_list_gene$x, each=nterms)
sum(unlist(gam_list_stables_genes) == gam_list_gene$x)
gam_list_stables$gene <- rep(unlist(gam_list_stables_genes), each=nterms) 
gam_list_stables <- gam_list_stables %>% mutate(padj = p.adjust(p.value,"fdr"))

plot_checks <- function(gene_name,gam_list,type="check",pages=1){
  gene_index <- which(gam_list_gene$x==gene_name)
  
  if(type=="check"){
    gam.check(gam_list[[gene_index]])
  } else if(type=="summary"){
    summary(gam_list[[gene_index]])
    
  } else if(type=="plot"){
    plot(gam_list[[gene_index]],pages=pages)
  }
}

sig_t1_genes <- gam_list_stables %>% mutate(padj = p.adjust(p.value,"fdr")) %>%
  filter(term=="s(t1):l1",
         padj < 0.05)

dir.create(sig_t1_partial_effect_resid_dir)
dir.create(check_plots_dir)
dir.create(sig_t1_partial_effect_resid_gg_dir)
dir.create(sig_ct_partial_effect_resid_gg_dir)

for(gene_name in sig_t1_genes$gene){
  png(paste0(sig_t1_partial_effect_resid_dir,"/",gene_name,".png"),
      width = 1800,height = 2000,res = 100,units = "px")
  plot_checks(gene_name,gam_list,type="plot",pages=1)
  dev.off()
}

for(gene_name in sig_t1_genes$gene){
  png(paste0(check_plots_dir,"/",gene_name,".png"),
      width = 1500,height = 1500,res = 200,units = "px")
  par(mfrow=c(2,2))
  plot_checks(gene_name,gam_list,type="check",pages=1)
  dev.off()
}

## give colors for sig smooth terms
term_levels <- rownames(summary(gam_list[[1]])$s.table)
for(gene_name in unique(sig_t1_genes$gene)){
  gene_index <- which(gam_list_gene$x==gene_name)
  png(paste0(sig_t1_partial_effect_resid_dir,"/",gene_name,"_wthResi.png"),
      width = 1800,height = 2000,res = 100,units = "px")
  plot_data =  plot(gam_list[[gene_index]],pages=1,by.resids=TRUE,residuals=TRUE)
  dev.off()
  all_panels = list()
  all_panels <- lapply(seq_len(length(plot_data)), function(s_term){
    p_vals <- gam_list_stables[gam_list_stables$term == 
                                 gsub(",(.*?)[)]",")",plot_data[[s_term]]$ylab) & 
                                 gam_list_stables$gene == gene_name,]
    n_na <- length(plot_data[[s_term]]$raw) - length(plot_data[[s_term]]$x)
    if(n_na < 0){
      data.frame(x = plot_data[[s_term]]$x,
                 fitted = plot_data[[s_term]]$fit,
                 se = plot_data[[s_term]]$se,
                 raw_x = c(plot_data[[s_term]]$raw,rep(NA,-n_na)),
                 p.resid = c(plot_data[[s_term]]$p.resid,rep(NA,-n_na)),
                 p_vals,
                 term = plot_data[[s_term]]$xlab,
                 ylab_label = plot_data[[s_term]]$ylab)
    } else {
      data.frame(x = c(plot_data[[s_term]]$x,rep(NA,n_na)),fitted = c(plot_data[[s_term]]$fit,rep(NA,n_na)),
                 se = c(plot_data[[s_term]]$se,rep(NA,n_na)),
                 raw_x = plot_data[[s_term]]$raw,p.resid = plot_data[[s_term]]$p.resid,
                 p_vals,term = plot_data[[s_term]]$xlab, ylab_label = plot_data[[s_term]]$ylab)
    }
    
  })
  
  all_panels <- do.call(rbind,all_panels)
  p_term = all_panels %>% mutate(term = factor(term,levels = term_levels)) %>%
    ggplot()+geom_path(mapping = aes(x=x,y=fitted,color= padj < 0.05))+
    geom_path(mapping = aes(x=x,y=fitted+se),linetype=2)+
    geom_path(mapping = aes(x=x,y=fitted-se),linetype=2)+
    theme_bw(base_size = 12)+facet_wrap(.~term,scales = "free")+xlab(gene_name)+
    geom_point(aes(x= raw_x,y=p.resid,color= padj < 0.05),alpha=0.1,size=0.1)+
    scale_color_manual(values=c("TRUE" = "lightblue","FALSE"="grey"))
  
  
  ggsave(plot=p_term,filename = paste0(sig_t1_partial_effect_resid_gg_dir,"/",
                                       gene_name,".png"),
         width = 14,height = 14,dpi = 200)
  
}
gam_list_stables %>% filter(gene %in% sig_t1_genes$gene) %>% 
  mutate(term = factor(term,levels = term_levels[nterms:1])) %>% filter(padj <0.05) %>%
  mutate(padj = ifelse(padj <=0, 1e-5,padj)) %>%
  ggplot()+geom_point(mapping = aes(x = gene,y=term,size= -log10(padj)),pch=1)+
  theme_bw(base_size = 15)+theme(axis.text.x = element_text(angle = 90))
ggsave(out_sig_t1_sterm_png,width = 20,height = 10)

## give colors for sig smooth terms
## not sig for t1 but for some cell types

term_levels <- rownames(summary(gam_list[[1]])$s.table)
nonsig_t1_genes <- gam_list_stables %>% filter(term=="s(t1):l1",
                                               padj > 0.05)
sig_ct_genes <- gam_list_stables %>% filter(gene %in% nonsig_t1_genes$gene,
                                            padj < 0.05)
dir.create(sig_ct_partial_effect_resid_dir)
for(gene_name in unique(sig_ct_genes$gene)){
  gene_index <- which(gam_list_gene$x==gene_name)
  png(paste0(sig_ct_partial_effect_resid_dir,"/",gene_name,"_wthResi.png"),
      width = 1800,height = 2000,res = 100,units = "px")
  plot_data =  plot(gam_list[[gene_index]],pages=1,by.resids=TRUE,residuals=TRUE)
  dev.off()
  all_panels = list()
  all_panels <- lapply(seq_len(length(plot_data)), function(s_term){
    p_vals <- gam_list_stables[gam_list_stables$term == 
                                 gsub(",(.*?)[)]",")",plot_data[[s_term]]$ylab) & 
                                 gam_list_stables$gene == gene_name,]
    
    n_na <- length(plot_data[[s_term]]$raw) - length(plot_data[[s_term]]$x)
    if(n_na < 0){
      data.frame(x = plot_data[[s_term]]$x,
                 fitted = plot_data[[s_term]]$fit,
                 se = plot_data[[s_term]]$se,
                 raw_x = c(plot_data[[s_term]]$raw,rep(NA,-n_na)),
                 p.resid = c(plot_data[[s_term]]$p.resid,rep(NA,-n_na)),
                 p_vals,
                 term = plot_data[[s_term]]$xlab,
                 ylab_label = plot_data[[s_term]]$ylab)
    } else {
      data.frame(x = c(plot_data[[s_term]]$x,rep(NA,n_na)),fitted = c(plot_data[[s_term]]$fit,rep(NA,n_na)),
                 se = c(plot_data[[s_term]]$se,rep(NA,n_na)),
                 raw_x = plot_data[[s_term]]$raw,p.resid = plot_data[[s_term]]$p.resid,
                 p_vals,term = plot_data[[s_term]]$xlab, ylab_label = plot_data[[s_term]]$ylab)
    }
  })
  
  all_panels <- do.call(rbind,all_panels)
  p_term = all_panels %>% mutate(term = factor(term,levels = term_levels)) %>%
    ggplot()+geom_path(mapping = aes(x=x,y=fitted,color= padj < 0.05))+
    geom_path(mapping = aes(x=x,y=fitted+se),linetype=2)+
    geom_path(mapping = aes(x=x,y=fitted-se),linetype=2)+
    theme_bw(base_size = 12)+facet_wrap(.~term,scales = "free")+xlab(gene_name)+
    geom_point(aes(x= raw_x,y=p.resid,color= padj < 0.05),alpha=0.1,size=0.1)+
    scale_color_manual(values=c("TRUE" = "lightblue","FALSE"="grey"))
  
  ggsave(plot=p_term,filename = paste0(sig_ct_partial_effect_resid_gg_dir,"/",
                                       gene_name,".png"),
         width = 14,height = 14,dpi = 200)
  
}

gam_list_stables %>% filter(gene %in% sig_ct_genes$gene) %>% 
  mutate(term = factor(term,levels = term_levels[nterms:1])) %>%
  filter(padj <0.05) %>%
  mutate(padj = ifelse(padj <=0, 1e-5,padj)) %>%
  ggplot()+geom_point(mapping = aes(x = gene,y=term,size= -log10(padj)),pch=1)+
  theme_bw(base_size = 15)+theme(axis.text.x = element_text(angle = 90,size=9))
ggsave(out_sig_ct_sterm_png,
       width = 26,height = 10)