new_samples = [ "TILD028MF","TILD049LF",
                "TILD080MF","TILD111LF","TILD113MF",
                "TILD117MF","TILD130MF",
                "TILD167LF","TILD299LF",
                "TILD315LF","VUHD038",
                "VUHD049","VUHD090",
                "VUILD141","VUILD142",
                "VUILD49","VUILD58"]
samples = ["VUILD96LF","VUHD116A","VUILD115","VUILD91LF","VUILD104MF",
                    "VUILD48MF","VUILD105LF","VUILD107MF","VUHD116B","VUILD102LF",
                    "VUHD095","VUILD102MF","VUILD48LF","VUILD91MF",
                    "VUHD113","VUILD104LF","VUILD78MF","THD0011","VUILD96MF",
                    "VUHD069","VUILD78LF","TILD117MF","TILD117LF","TILD175",
                    "THD0008","VUILD105MF"]

TMAs = ["TMA1","TMA1","TMA3","TMA4","TMA2",
          "TMA2","TMA2","TMA1","TMA1","TMA1",
          "TMA2","TMA1","TMA2","TMA4",
          "TMA2","TMA2","TMA4","TMA4","TMA1",
          "TMA2","TMA4","TMA4","TMA4","TMA4",
          "TMA3","TMA2"]

def get_mem_mb_sm(wildcards, attempt):
    return attempt * 6000

def get_mem_mb(wildcards, attempt):
    return attempt * 12000

rule all:
    input: 
        expand("output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/{sample}.hex_summary.png",nroots=5000, d = 3.0, sample = new_samples),
        expand(expand("output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/{{sample}}_{{tma}}.hex_summary.png",nroots=5000, d = 3.0), zip,sample = samples,tma=TMAs)

rule plot_hex_new_samples:
    input:
        script = "code/svi_plot_hex_summary_xenium_revision_new_ori_samples_June2024.R",
        node_meta =  "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{d}/revision_new_samples/min10/{sample}_node_meta.csv",
        cid_results = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{d}/revision_new_samples/min10/gmm12_ori_new/{sample}_gmm12.txt"
    output:
        hex_plot = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/{sample}.hex_summary.png"
    params:
        bin_thresh = 10,
        bin_width = 5
    resources:
        cpus = 2,
        mem = get_mem_mb_sm,
        gpus=0
    log:
        "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/{sample}.hex_summary.log"
    shell:
        """
            /opt/R/4.1.2/bin/R CMD BATCH --no-restore --no-save "--args \
               sample='{wildcards.sample}' node_meta_file='{input.node_meta}' \
               outpng='{output.hex_plot}' bin_thresh='{params.bin_thresh}' \
               bin_width='{params.bin_width}' cid_results='{input.cid_results}' " {input.script} {log}
        """

rule plot_hex_ori_samples:
    input:
        script = "code/svi_plot_hex_summary_xenium_revision_new_ori_samples_June2024.R",
        node_meta =  "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{d}/min10/gmm12_ori_new/{sample}_{tma}_node_meta.csv",
        cid_results = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{d}/revision_new_samples/min10/gmm12_ori_new/{sample}_{tma}_gmm12.txt"
    output:
        hex_plot = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/{sample}_{tma}.hex_summary.png"
    params:
        bin_thresh = 10,
        bin_width = 5
    resources:
        cpus = 2,
        mem = get_mem_mb_sm,
        gpus=0
    log:
        "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/{sample}_{tma}.hex_summary.log"
    shell:
        """
            /opt/R/4.1.2/bin/R CMD BATCH --no-restore --no-save "--args \
               sample='{wildcards.sample}' node_meta_file='{input.node_meta}' \
               outpng='{output.hex_plot}' bin_thresh='{params.bin_thresh}' \
               bin_width='{params.bin_width}' cid_results='{input.cid_results}' " {input.script} {log}
        """

rule getDEGs_perk:
    input:
        node_meta = expand(["output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{{nroots}}_aug2023/dmax{{d}}/min10/gmm12_ori_new/{sample}_{tma}_node_meta.csv"],zip, sample = samples, tma=TMAs),
        script = "code/svi_get_DEGs_heatmaps_xenium_may2024.R"
    output:
        out_png_topFDR = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/original_topFDR_ht.png",
        out_png_topPropRatio = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/original_topPropRatio_ht.png",
        out_png_topFDR_tmean = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/original_topFDR_ht_mean.png",
        out_png_topPropRatio_tmean = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/original_topPropRatio_ht_mean.png",
        out_png_wide_dot = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/original_dotplotWide.png",
        out_trans_counts_rds = "output/savedRDS/xenium/{nroots}_dmax{d}/gmm12_ori_new/original_out_trans_counts.rds",
        out_niche_prop_bar_plot = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/original_niche_prop_bar.png",
        out_niche_trans_count_plot = "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/original_niche_count.png"
    params:
        top_n = 5,
        node_meta_dir = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{d}/min10/gmm12_ori_new/"
    resources:
        cpus = 1,
        mem = get_mem_mb
    log:
        "output/figures/xenium/{nroots}_dmax{d}/gmm12_ori_new/generate_DEGs_plots.log"
    shell:
        """
            /opt/R/4.1.2/bin/R CMD BATCH --no-restore --no-save "--args \
               k_dir='gmm12_ori_new' out_png_topFDR='{output.out_png_topFDR}' out_png_topPropRatio='{output.out_png_topPropRatio}' \
               out_png_topFDR_tmean='{output.out_png_topFDR_tmean}' out_png_topPropRatio_tmean='{output.out_png_topPropRatio_tmean}' \
               out_png_wide_dot='{output.out_png_wide_dot}' out_trans_counts_rds='{output.out_trans_counts_rds}' \
               out_niche_prop_bar_plot='{output.out_niche_prop_bar_plot}' \
               top_n='{params.top_n}' node_meta_dir='{params.node_meta_dir}' \
               out_niche_trans_count_plot='{output.out_niche_trans_count_plot}' k='12' " {input.script} {log}
        """       