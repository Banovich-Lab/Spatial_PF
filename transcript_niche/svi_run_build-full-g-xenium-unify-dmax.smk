## run subgraphs for each Xenium csv file with full list of genes

samples = ["VUILD96LF","VUHD116A","VUILD115","VUILD91LF","VUILD104MF","VUILD48MF","VUILD105LF","VUILD107MF","VUHD116B","VUILD102LF",
           "VUHD095","VUILD102MF","VUILD48LF","VUILD106","VUILD91MF","VUHD113","VUILD104LF",
           "VUILD78MF","THD0011","VUILD96MF","VUHD069","VUILD78LF","TILD117MF","TILD117LF","TILD175","VUILD110","THD0008","VUILD105MF"]

TMAs = ["TMA1","TMA1","TMA3","TMA4","TMA2","TMA2","TMA2","TMA1","TMA1","TMA1","TMA2","TMA1","TMA2","TMA3","TMA4","TMA2","TMA2","TMA4","TMA4",
        "TMA1","TMA2","TMA4","TMA4","TMA4","TMA4","TMA3","TMA3","TMA2"]

#scratch_dir= "/data/scratch/projects/punim0741/"
raw_detected_tx_dir = "/mnt/beegfs/mccarthy/backed_up/general/rlyu/Dataset/LFST_2022/Xenium/fullPanelQC/"
scratch_ln = "scratch_ln/"
rule all:
    input:
        expand(expand(["output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/min10/{{vname}}_50_embeddings_{{TMA}}.npy"],nroots=5000,k=3.0),zip,vname=samples, TMA=TMAs)

def get_mem_mb(wildcards, attempt):
    return 20000 + attempt * 40000
def get_mem_mb_merge(wildcards, attempt):
    return 320000 + attempt * 40000

def get_mem_mb_large(wildcards, attempt):
    if wildcards.vname == ("VUILD115") or wildcards.vname == ("VUILD106") or wildcards.vname == ("VUILD110") or wildcards.vname == ("THD0008"): 
        return 160000 + attempt * 40000
    else:
        return 20000 + attempt * 40000

def get_mem_mb_training(wildcards, attempt):
    if wildcards.nroots == ("20000"):
        return 80000 + attempt * 20000
    else:
        return 55000 + attempt * 20000

def get_csv(wildcards):
    file_name = raw_detected_tx_dir + wildcards.vname+"_"+wildcards.TMA+"_detected_transcripts.csv" 
    return file_name


## fullPanelQC/{vname}_{TMA}_detected_transcripts.csv", after filtering out QV < 20

## build the mRNA transcript graph per sample where each node represents one detected transcript, edges are added between two nodes if their spatial distance is < 3.0
rule build_fulllgraph_bydmax:
    input:
        script = "code/svi_buildFullGraph_input_dmax.py",
        fullpn = get_csv,
        xenium_gene_panel = "output/xenium/xenium_gene_panel.csv"
    output:
        fullgraph = temp(scratch_ln+"xenium_fullpanel_graphs_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.gpickle"),
        fullgraph_meta = scratch_ln+"xenium_fullpanel_graphs_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.node_meta.csv",
        n_connected_comps =  scratch_ln+"xenium_fullpanel_graphs_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.gpickle.components.csv"
    log:
        "output/xenium/fullPanel/graphs/logs/xenium_fullpanel_graphs_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}_buildFullGraph.log"
    params:
        text_log = "output/xenium/fullPanel/graphs/logs/xenium_fullpanel_graphs_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}_buildFullGraph.log.txt"
    resources:
        cpus = 1,
        mem = get_mem_mb_large,
        gpus =0 
    conda:
        "/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/mambaForge/mambaforge/envs/graphsageAug"
    shell:
        """
            python -u {input.script} {wildcards.vname} {wildcards.TMA} {input.fullpn} {input.xenium_gene_panel} \
            {output.fullgraph_meta} {output.n_connected_comps} {wildcards.k} {output.fullgraph} {params.text_log} 2>> {log}
        """

rule filter_compoments_and_subgraph:
    input:
        full_graph = scratch_ln+"xenium_fullpanel_graphs_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.gpickle",
        script = "code/svi_run_filter_and_subgraph_3NB_aug2023.py"
#        fullgraph_meta = scratch_ln+"xenium_fullpanel_graphs_aug2023/{vname}_{TMA}.node_meta.csv"
    output:
        subgraph3NB = scratch_ln+"xenium_subgraph3NB{nroots}_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.subgraph3NB.gpickle",
        subgraph3NB_rootnodesID = scratch_ln+"xenium_subgraph3NB{nroots}_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.subgraph3NB.rootNodesID.csv",
    log:
        "output/xenium/fullPanel/graphs/logs/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.subgraph3NB.log"
    params:
        text_log = "output/xenium/fullPanel/graphs/logs/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.subgraph3NB.log.txt",
        remove_comp_min_n = 10
    resources:
        cpus = 1,
        mem = get_mem_mb,
        gpus = 0
    conda:
        "/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/mambaForge/mambaforge/envs/graphsageAug"
    shell:
        """
            python {input.script} {wildcards.vname} {input.full_graph} {params.remove_comp_min_n} {wildcards.nroots} {output.subgraph3NB} {output.subgraph3NB_rootnodesID} {params.text_log} 2>> {log}

        """

rule merge_subgraphs:
    input:
        script="code/svi_run_merge_subgraphs_xenium_aug2023.py",
        gpickles = expand([scratch_ln+"xenium_subgraph3NB{{nroots}}_aug2023/dmax{{k}}/{vname}_{TMA}_dmax{{k}}.subgraph3NB.gpickle"],zip, vname=samples, TMA=TMAs),
    output:
        sg_merge = scratch_ln+"xenium_subgraph3NB{nroots}_aug2023/dmax{k}/merged_int_subgraph_sg_dmax{k}.gpickle",
        rootID_reindex = scratch_ln+"xenium_subgraph3NB{nroots}_aug2023/dmax{k}/merged_int_subgraph_sg_dmax{k}_rootIDs.csv"
    resources:
        cpus = 1,
        mem = get_mem_mb_merge,
        gpus = 0
    benchmark:
        "output/xenium/fullPanel/graphs/benchmarks/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/mergeSubGraphsToSG_dmax{k}.benchmark.txt"
    log:
        "output/xenium/fullPanel/graphs/logs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/mergeSubGraphsToSG_dmax{k}.log"
    conda:
        "/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/mambaForge/mambaforge/envs/graphsageAug"
    params:
        subgraph_dir = scratch_ln+"xenium_subgraph3NB{nroots}_aug2023/dmax{k}/",
        subgraphs = ",".join(expand([scratch_ln+"xenium_subgraph3NB{{nroots}}_aug2023/dmax{{k}}/{vname}_{TMA}_dmax{{k}}.subgraph3NB.gpickle"],zip, vname=samples, TMA=TMAs)),
        rootID_csvs = ",".join(expand([scratch_ln+"xenium_subgraph3NB{{nroots}}_aug2023/dmax{{k}}/{vname}_{TMA}_dmax{{k}}.subgraph3NB.rootNodesID.csv"],zip, vname=samples, TMA=TMAs)),
        rootID_dir = scratch_ln+"xenium_subgraph3NB{nroots}_aug2023/dmax{k}/",
        textlog = "output/xenium/fullPanel/graphs/logs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/mergeSubGraphsToSG.log.txt"
    shell:
        """
            export LD_LIBRARY_PATH="/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/mambaForge/mambaforge/envs/graphsageAug/lib/python3.8/site-packages/tensorrt/:$LD_LIBRARY_PATH"
            python -u {input.script} {output.sg_merge} {params.subgraphs} {params.rootID_csvs} {output.rootID_reindex} {params.textlog} 2>> {log}

        """   
# ## prepared sg graph union from all 
rule trainUsingMergedgraph:
    input:
        gpickle = scratch_ln+"xenium_subgraph3NB{nroots}_aug2023/dmax{k}/merged_int_subgraph_sg_dmax{k}.gpickle",
        rootids = scratch_ln+"xenium_subgraph3NB{nroots}_aug2023/dmax{k}/merged_int_subgraph_sg_dmax{k}_rootIDs.csv",
        script = "code/svi_run_training_subgraphs_xenium_fullpanel.py"
    output:
        trained_model = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/saved_model/model/saved_model.pb",
        trained_emModel = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/saved_model/embmodel/saved_model.pb"
    benchmark:
        "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/benchmarks/model_training_dmax{k}.ben.txt"
    resources:
        cpus = 1,
        mem = get_mem_mb_training,
        gpus = 1
    params:
        text_log = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/model_training_dmax{k}.txt",
        number_of_walks = 5,
        num_length =  2,
        number_of_samples1 = 20,
        number_of_samples2 = 10,
        number_epoch = 10,
        trained_model = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/saved_model/model",
        trained_emModel  = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/saved_model/embmodel"
    log:
        "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/model_training_dmax{k}.log"
    conda:
        "/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/mambaForge/mambaforge/envs/graphsageAug"
    shell:
        """
           export LD_LIBRARY_PATH="/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/mambaForge/mambaforge/envs/graphsageAug/lib/python3.8/site-packages/tensorrt/:$LD_LIBRARY_PATH"

           python -u {input.script} {params.number_of_walks} {params.num_length}  {params.number_of_samples1}  {params.number_of_samples2}  {input.gpickle} {params.number_epoch} \
           {params.trained_model} {params.trained_emModel} {input.rootids} {params.text_log} 2>> {log}

        """

## Get embedding prediction per sample
rule embedd_graph:
    input:
        script = "code/svi_get_emb_xenium_fullpanel_inputg_filterMinComp_aug2023.py",
        full_g = scratch_ln+"xenium_fullpanel_graphs_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.gpickle",
        full_g_csv = scratch_ln+"xenium_fullpanel_graphs_aug2023/dmax{k}/{vname}_{TMA}_dmax{k}.node_meta.csv",
        trained_model = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/saved_model/embmodel/saved_model.pb"
    output:
        npy = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/min10/{vname}_50_embeddings_{TMA}.npy",
        node_meta_min10 = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/min10/{vname}_{TMA}_node_meta.csv"
    resources:
        cpus = 1,
        mem = get_mem_mb_large,
        gpus = 0
    benchmark:
        "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/min10/logs/{vname}_50_embeddings_{TMA}.ben.txt"
    conda:
        "/mnt/beegfs/mccarthy/backed_up/general/rlyu/Software/mambaForge/mambaforge/envs/graphsageAug"
    log:
        "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/min10/logs/{vname}_50_embeddings_{TMA}.log"
    params:
        textlog ="output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/min10/logs/{vname}_{TMA}.get.fullp.logger.txt",
        trained_emModel  = "output/xenium/fullPanel/graphs/3NB/xenium_subgraph3NB{nroots}_aug2023/dmax{k}/saved_model/embmodel/",
        number_of_samples1 = 20,
        number_of_samples2 = 10,
        min_comp = 10
    shell:
        """
           python -u {input.script} {wildcards.vname} {wildcards.TMA} {output.npy} {output.node_meta_min10} \
            {params.trained_emModel} {input.full_g} {input.full_g_csv} {params.number_of_samples1} {params.number_of_samples2} {params.min_comp} {params.textlog} 2>> {log}

        """