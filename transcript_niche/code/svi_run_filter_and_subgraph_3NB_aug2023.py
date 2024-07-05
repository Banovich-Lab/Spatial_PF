# code/svi_run_filter_and_subgraph_3NB_aug2023.py
# coding: utf-8

## Generate subgraphs

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
from stellargraph.data import EdgeSplitter
from stellargraph.mapper import GraphSAGELinkGenerator
from stellargraph.layer import GraphSAGE, link_classification
from stellargraph.layer.graphsage import AttentionalAggregator
from stellargraph.data import UniformRandomWalk
from stellargraph.data import UnsupervisedSampler
import tensorflow as tf
from tensorflow import keras
from stellargraph import globalvar


from numpy.random import seed
seed(42)
from tensorflow.random import set_seed
set_seed(42)

import sys
import logging 
import builtins
import shutil

logging.basicConfig(filename=sys.argv[7], 
                        format='%(asctime)s %(message)s', 
                        filemode='w',level=logging.INFO) 

logger=logging.getLogger() 

logger.info('GPU name: '+ str(tensorflow.config.experimental.list_physical_devices("GPU")))
logger.info('Length of list:'+ str(len(sys.argv)))

samplename = sys.argv[1]
## length of walks for sampling positive pairs
full_graph_g =  sys.argv[2]
comp_min_n = int(sys.argv[3])
num_root_nodes = int(sys.argv[4])
gpickle = sys.argv[5]
sampledRootsID_csv = sys.argv[6]

logger.info('samplename '+ samplename) 
logger.info('remove comps smaller than '+str(comp_min_n))

logger.info('num sampleRoots '+str(num_root_nodes))
logger.info('out_subgraph gpickle '+str(gpickle))

logger.info('selected rootIDs '+str(sampledRootsID_csv))

def filterGraph(G,min_n):
    node_removed = []
    for component in tqdm(list(nx.connected_components(G))):
        if len(component) < min_n:
            for node in component:
                node_removed.append(node)
                G.remove_node(node)
    logger.info("node removed: " + str(len(node_removed)))
    return G

logger.info('Filtering graph' + samplename)
logger.info('Removed comps smaller than ' + str(comp_min_n))
g = nx.read_gpickle(full_graph_g)
logger.info("full_graph_g number_of_nodes() " + str(g.number_of_nodes()))
logger.info('First 5 node label' + str(list(g.nodes())[0:5]))

g = filterGraph(g, comp_min_n)
logger.info("filtered g  " + str(g.number_of_nodes()))
logger.info('First 5  node label' + str(list(g.nodes())[0:5]))

logger.info("Generate 3NB subgraphs of sampled root nodes")

set_seed(100)
number_roots = num_root_nodes
logger.info('Sampling number of root nodes ' + str(number_roots))

selected_roots = random.sample(g.nodes,number_roots)

three_hop_nbrs = [nx.dfs_edges(g,n,3) for n in selected_roots]

three_hop_nbrs_list = [i for sublist in three_hop_nbrs for i in sublist]
## a list of tupple of size 2
sub_nodes = set([])
for k in three_hop_nbrs_list:
    sub_nodes.add(k[0])
    sub_nodes.add(k[1]) 

logger.info("Keep number of nodes in the subgraph  " + str(len(sub_nodes)))
subgraph = g.subgraph(sub_nodes).copy()
logger.info("Nodes in the subgraph  " + str(subgraph.number_of_nodes()))
logger.info('Relabel nodes in subgraph adding prefix g_ ')
nx.relabel_nodes(subgraph,mapping = lambda x: "g_"+str(x),copy=False)

logger.info('Creating mapping dict')

mapIdInt = dict()
i = 0
for node in subgraph.nodes():
    mapIdInt[node] = i
    i += 1
logger.info("dict length "+str(len(mapIdInt)))
logger.info('Relabel nodes in subgraph with consecutive ints ')

nx.relabel_nodes(subgraph,mapping =mapIdInt,copy=False)

logger.info('Saving the subgraph with node id starting 0 ')

nx.write_gpickle(subgraph,gpickle)

logger.info('Finished saving the subgraph with node id starting to ' + gpickle)

root_nodes_node_id_to_int  = [mapIdInt[("g_"+str(sn))] for sn in selected_roots]

logger.info('Selected root IDs are mapped to the new integer id labels and saved to ' + sampledRootsID_csv)

pd.DataFrame(root_nodes_node_id_to_int).to_csv(sampledRootsID_csv)


