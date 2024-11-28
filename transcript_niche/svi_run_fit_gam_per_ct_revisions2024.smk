## fit GAM with time and cell type proportion per CT

celltypes = ["AT1","AT2","Alveolar.FBs","Alveolar.Macrophages","Capillary","Inflammatory.FBs","Interstitial.Macrophages",
            "Mast","Monocytes.MDMs","NK.NKT","Neutrophils","Proliferating.AT2","Proliferating.FBs","SMCs.Pericytes",
            "Activated.Fibrotic.FBs","Arteriole","CD4..T.cells","Proliferating.Myeloid","SPP1..Macrophages",
            "Transitional.AT2","Venous","cDCs","Macrophages...IFN.activated","RASC","Secretory","Multiciliated",
            "pDCs","Adventitial.FBs","CD8..T.cells","Lymphatic","Proliferating.T.cells","Proliferating.B.cells","Subpleural.FBs","Tregs",
            "Goblet","Plasma","Proliferating.NK.NKT","Basal","Myofibroblasts","B.cells","Migratory.DCs","Basophils",
            "KRT5..KRT17.","Proliferating.Airway","Mesothelial","Langerhans.cells","PNEC"]

keep_celltype = ['AT1','AT2','Alveolar.FBs','Alveolar.Macrophages','Capillary','Inflammatory.FBs','Interstitial.Macrophages',
                'Mast','Monocytes.MDMs','NK.NKT','Neutrophils','Proliferating.AT2','Proliferating.FBs','SMCs.Pericytes',
                'Activated.Fibrotic.FBs','Arteriole','CD4..T.cells','Proliferating.Myeloid', 'SPP1..Macrophages',
                  'Transitional.AT2','Venous','cDCs','Macrophages...IFN.activated','RASC','Secretory','Multiciliated',
                  'pDCs','Adventitial.FBs','CD8..T.cells','Lymphatic','Proliferating.T.cells',"Proliferating.B.cells",'Subpleural.FBs','Tregs',
                'Goblet','Plasma','Proliferating.NK.NKT','Basal','Myofibroblasts','B.cells','Migratory.DCs','Basophils',
                'KRT5..KRT17.', 'Proliferating.Airway','Mesothelial','Langerhans.cells','PNEC']

rule all:
        input:
            expand(["output/savedRDS/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/{c}_min{filter_by}_major_{major_only}.gam_list.rds",
            "output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/sig_t1_withCT_gg/{c}_min{filter_by}_major_{major_only}.sig.t1.png"],c=keep_celltype,major_only=["true"],filter_by= [0.5])

rule run_gam:
    input:
        script = "code/revision2024/lumen_orders_gene_prop_gam_fit_per_ct_slurm.R"
    output:
        knots_png = "output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/{c}_min{filter_by}_major_{major_only}.knots.png",
        rds = "output/savedRDS/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/{c}_min{filter_by}_major_{major_only}.gam_list.rds",
        rds_gene_names = "output/savedRDS/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/{c}_min{filter_by}_major_{major_only}.gam_list.rds.txt"
    resources:
        cpus = 5,
        mem = 15000
    params:
        min_lumen = 80
    log:
        "output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/{c}_min{filter_by}_major_{major_only}.fitGAM.log"
    shell:
        """
            /opt/R/4.1.2/bin/R CMD BATCH --no-restore --no-save "--args \
              celltype='{wildcards.c}' outrds='{output.rds}' filter_by='{wildcards.filter_by}' use_major_ct='{wildcards.major_only}' \
               out_knots_png='{output.knots_png}' min_lumen='{params.min_lumen}' " {input.script} {log}
  
        """
rule summary_plots:
    input:
        script = "code/revision2024/lumen_orders_gene_prop_gam_fit_per_ct_plot_slurm.R",
        gam_list = "output/savedRDS/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/{c}_min{filter_by}_major_{major_only}.gam_list.rds",
        gam_list_gene = "output/savedRDS/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/{c}_min{filter_by}_major_{major_only}.gam_list.rds.txt"
    output:
        sig_t1_png = "output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/sig_t1_withCT_gg/{c}_min{filter_by}_major_{major_only}.sig.t1.png",
        sig_ct_png = "output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/nonsig_t1_withCT_gg/{c}_min{filter_by}_major_{major_only}.sig.ct.png"        
    resources:
        cpus = 1,
        mem = 10000
    params:
        sig_t1_partial_effect_resid_dir="output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/sig_t1_withCT/",
        sig_t1_checks="output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/sig_t1_withCT_checks/",
        sig_ct_partial_effect_resid_dir="output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/nonsig_t1_withCT/",
        sig_t1_partial_effect_resid_gg_dir="output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/sig_t1_withCT_gg/",
        sig_ct_partial_effect_resid_gg_dir="output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/nonsig_t1_withCT_gg/",
    log:
        "output/figures/revisions2024/perCT/{c}_min{filter_by}_major_{major_only}/{c}_min{filter_by}_major_{major_only}.gam_plot.log"
    shell:
        """
            /opt/R/4.1.2/bin/R CMD BATCH --no-restore --no-save "--args \
              gam_list='{input.gam_list}' gam_list_gene='{input.gam_list_gene}'  \
               out_sig_t1_sterm_png='{output.sig_t1_png}' out_sig_ct_sterm_png='{output.sig_ct_png}' check_plots_dir='{params.sig_t1_checks}' \
               sig_t1_partial_effect_resid_dir='{params.sig_t1_partial_effect_resid_dir}' sig_ct_partial_effect_resid_dir='{params.sig_ct_partial_effect_resid_dir}' \
                sig_t1_partial_effect_resid_gg_dir='{params.sig_t1_partial_effect_resid_gg_dir}' \
                sig_ct_partial_effect_resid_gg_dir='{params.sig_ct_partial_effect_resid_gg_dir}'  " {input.script} {log}
  
        """
