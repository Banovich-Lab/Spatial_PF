# Transcript-based analysis

## svi_run_build-full-g-xenium-unify-dmax.smk

**A snakemake file that defines the following rules**: 
1. Build mRNA spatial graphs per sample using radius threshold d
1. Create subgraph per sample from the randomly sample N root nodes
1. Subgraphs from 28 samples were merged into one graph
1. Two-hop GraphSAGE model was trained on the merged graph
1. Embedding prediction using trained model

## Main package versions

1. Python 3.8.0
      - networkx==2.8.8
      - stellargraph==1.2.1
      - numba==0.58.1
      - numpy==1.24.4
      - keras==2.11.0
      - nvidia-cublas==11.5.1.101
      - nvidia-cublas-cu11==2022.4.8
      - nvidia-cublas-cu117==11.10.1.25
      - nvidia-cuda-nvrtc==11.1.105
      - nvidia-cuda-runtime==11.3.58
      - nvidia-cudnn==8.2.0.51
      - nvidia-cudnn-cu11==8.6.0.163
      - nvidia-pyindex==1.0.9
      - nvidia-tensorrt==7.2.3.4
      - pandas==1.5.2
      - scikit-learn==1.3.0
      - tensorflow==2.11.0
      - tensorflow-estimator==2.11.0
      - tqdm==4.66.0
