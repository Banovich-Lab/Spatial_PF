## build mRNA spatial graph 
# 
# coding: utf-8
import tensorflow
import tensorrt
import numpy as np
import random as rn

np.random.seed(42)
rn.seed(12345)
from tensorflow.keras import backend as K
tensorflow.random.set_seed(12345)
import networkx as nx
import pandas as pd
import os
import random
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.spatial import cKDTree as KDTree
from tensorflow.keras.utils import to_categorical

import stellargraph as sg
from stellargraph.data import UniformRandomWalk
import tensorflow as tf
from tensorflow import keras
from stellargraph import globalvar
from tensorflow.random import set_seed
import sys
import logging 
import builtins
import os
import shutil

logging.basicConfig(filename=sys.argv[9], 
                        format='%(asctime)s %(message)s', 
                        filemode='w',level=logging.INFO) 

#Let us Create an object for logging
logger=logging.getLogger() 

logger.info('GPU name: '+ str(tensorflow.config.experimental.list_physical_devices("GPU")))
logger.info('Length of list:'+ str(len(sys.argv)))

samplename = sys.argv[1]
## length of walks for sampling positive pairs
tma = sys.argv[2]
raw_csv = sys.argv[3]
xenium_gene_panel = sys.argv[4]
fullgraph_meta_csv = sys.argv[5]
components_len_csv = sys.argv[6]
dmax = float(sys.argv[7])
out_gpickle = sys.argv[8]

logger.info('samplename '+ samplename) 
logger.info('tma '+str(tma))
logger.info('input raw_csv from '+str(raw_csv))
logger.info('full graph to '+str(out_gpickle))
logger.info('full graph meta to '+str(fullgraph_meta_csv))
logger.info('dmax set as '+str(dmax))

logger.info('graph count components to '+str(components_len_csv))
 
samplefile = raw_csv
barcodes_df = pd.read_csv(samplefile, sep = ",", header=0)

logger.info(str(barcodes_df.shape))
logger.info(str(barcodes_df.head(2)))

xenium_gene_panel = pd.read_csv(xenium_gene_panel)
# xenium_gene_panel
tagList_df = pd.DataFrame(np.unique(xenium_gene_panel.x.values),
                         columns=['gene'])
logger.info("gene panel "+ str(tagList_df.head(3)))
logger.info(str(tagList_df.shape))
logger.info("filter nodes by gene panel")
barcodes_df_filter = barcodes_df[barcodes_df["feature_name"].isin(tagList_df["gene"])]

assert(tagList_df.shape[0]==343)
logger.info("# Nodes not in any cells "+str(barcodes_df.cell_id.str.contains("-1").sum()))

logger.info("# Nodes after filtering by gene panel")
logger.info("barcodes_df_filter shape "+str(barcodes_df_filter.shape))

logger.info("Using preset dmax as "+str(dmax))

d_max = dmax
gene_list = tagList_df.gene.values
logger.info("One-hot-encoding for gene list")
one_hot_encoding = dict(zip(gene_list, to_categorical(np.arange(gene_list.shape[0]),num_classes=gene_list.shape[0]).tolist()))

logger.info("building graph and not removing any small components yet")

def buildGraph(barcodes_df, d_th, one_hot_encoding):
    G = nx.Graph()
    features =[] 
    node_removed = []
    barcodes_df.reset_index(drop=True, inplace=True)
    barcodes_df["feature"] = barcodes_df['feature_name'].map(one_hot_encoding).tolist()
    features.append(np.vstack(barcodes_df.feature.values))

    kdT = KDTree(np.array([barcodes_df.x_location.values,barcodes_df.y_location.values]).T)
    res = kdT.query_pairs(d_th)
    res = [(x[0],x[1]) for x in list(res)]

    # Add nodes to graph
    G.add_nodes_from((barcodes_df.index.values), test=False, val=False, label=0)
    # Add node features to graph
    nx.set_node_attributes(G,dict(zip((barcodes_df.index.values), barcodes_df.feature)), 'feature')
    # Add edges to graph
    G.add_edges_from(res)

    return G,barcodes_df

g, graph_meta = buildGraph(barcodes_df_filter, d_max, one_hot_encoding)

logger.info("g.number_of_nodes()" + str(g.number_of_nodes()))

logger.info("graph_meta csv shape" + str(graph_meta.shape))
assert(graph_meta.shape[0] == barcodes_df_filter.shape[0])
logger.info('Finished building graph' + samplename)
logger.info('Graph with nodes: ' + str(g.number_of_nodes()))

logger.info('First 5 node label' + str(list(g.nodes())[0:5]))
logger.info('Writing full graph to ' + out_gpickle)
nx.write_gpickle(g,out_gpickle)

logger.info('Writing graph meta ' + fullgraph_meta_csv)

graph_meta.drop(columns=['feature'],inplace=True)
graph_meta.to_csv(fullgraph_meta_csv)

## finding connected components 
## useful for deciding filtering cut-offs
logger.info("finding connected components ")

components_len = []
for component in tqdm(list(nx.connected_components(g))):
    components_len.append(len(component))
pd.DataFrame(components_len).to_csv(components_len_csv)


