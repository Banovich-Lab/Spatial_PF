# Transcript-based analysis

## svi_run_build-full-g-xenium-unify-dmax.smk

**A snakemake file that defines the following rules**: 
1. Build mRNA spatial graphs per sample using radius threshold d
1. Create subgraph per sample from the randomly sample N root nodes
1. Subgraphs from 28 samples were merged into one graph
1. Two-hop GraphSAGE model was trained on the merged graph
1. Embedding prediction using trained model