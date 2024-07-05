# code/svi_run_merge_subgraphs_xenium_aug2023.py
# coding: utf-8


import tensorflow
import tensorrt
import numpy as np
import random as rn

np.random.seed(42)

rn.seed(12345)
from tensorflow.keras import backend as K
tensorflow.random.set_seed(1234)

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
from sklearn.model_selection import train_test_split

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

logging.basicConfig(filename=sys.argv[5], 
                        format='%(asctime)s %(message)s', 
                        filemode='w',level=logging.INFO) 

logger=logging.getLogger() 

logger.info('GPU name: '+ str(tensorflow.config.experimental.list_physical_devices("GPU")))

logger.info('Length of list:'+ str(len(sys.argv)))

output_gpickle = sys.argv[1]
## separate by ","
subgraphs =  sys.argv[2]
rootsIDs =  sys.argv[3]

merged_reindex_rootsID = sys.argv[4]

graph_list = []
rootsID_list = []
subgraphs = subgraphs.split(",")
rootsID_csvs = rootsIDs.split(",")

last_graph_n = 0
for subg_g,subg_rootnode in zip(subgraphs,rootsID_csvs):
    logger.info(subg_rootnode)
    g = nx.read_gpickle(subg_g)
    rootsID = pd.read_csv(subg_rootnode)
    nx.relabel_nodes(g, lambda x: x + last_graph_n,copy=False)
    logger.info("last_graph_n "+ str(last_graph_n))
    rootsID["0"] = rootsID["0"] + last_graph_n
    rootsID_list.append(rootsID)
    last_graph_n = last_graph_n + g.number_of_nodes()
    graph_list.append(g)
    logger.info(list(g.nodes())[0:4])

# # Sample nodes from the large graph, also find the two hop neighbours for the root nodes 

logger.info("length graph: " +str(len(graph_list)))
rootsID_df = pd.concat(rootsID_list)
## these nodes have integer IDs of the merged graph
logger.info("length rootsID_df : " + str((rootsID_df).shape))
logger.info((rootsID_df.head(2)))

k = 0 
for i in graph_list:
    k += i.number_of_nodes()

logger.info("total nodes in merged graph "+str(k))
print("total nodes in merged graph  "+str(k))

large_g_roots = nx.union_all(graph_list)

logger.info("number of nodes in joined graph "+str(large_g_roots.number_of_nodes()))

logger.info("Converting to stellargraph class ")
large_g_roots = sg.StellarGraph(large_g_roots, node_features="feature")
logger.info("Converted to sg  graph ")

logger.info(large_g_roots.info())

logger.info("Saving to sg graph " + output_gpickle)

nx.write_gpickle(large_g_roots,output_gpickle)
logger.info("save reindex roots id to " + merged_reindex_rootsID)
rootsID_df.to_csv(merged_reindex_rootsID)