## code/get_emb_xenium_fullpanel_inputg_filterMinComp_aug2023.py
## embed nodes in gpickle using trained model

# coding: utf-8

import tensorflow
print(tensorflow.__version__)
import tensorrt
print(tensorrt.__version__)

import numpy as np
import random as rn

np.random.seed(42)
rn.seed(12345)
from tensorflow.keras import backend as K

tensorflow.random.set_seed(1234)

import networkx as nx
import pandas as pd
import numpy as np
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
from sklearn import preprocessing, feature_extraction, model_selection
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.metrics import accuracy_score

from stellargraph import globalvar

from numpy.random import seed
seed(42)
from tensorflow.random import set_seed
set_seed(42)

import sys
import logging 
import builtins
import os
import shutil
## file.py num_walks, num_root_nodes, number_of_samples1, number_of_samples2,buildGraph: (KNN/radius), output_dir
#now we will Create and configure logger 

logging.basicConfig(filename=sys.argv[11], 
                        format='%(asctime)s %(message)s', 
                        filemode='w',level=logging.INFO) 

#Let us Create an object 
logger=logging.getLogger() 

logger.info('GPU name: '+ str(tensorflow.config.experimental.list_physical_devices("GPU")))

### Read in the graph for each region
logger.info('Length of list:'+ str(len(sys.argv)))

samplename = sys.argv[1]
## length of walks for sampling positive pairs
tma = sys.argv[2]
outnpy =  sys.argv[3]
outnpy_meta =  sys.argv[4]

trained_model =  sys.argv[5]
gpickle = sys.argv[6]
gpickle_meta = sys.argv[7]

number_of_samples1 =  int(sys.argv[8])
number_of_samples2 =  int(sys.argv[9])
min_comp = int(sys.argv[10])

logger.info('samplename '+ samplename) 
logger.info('tma '+str(tma))
logger.info('outnpy '+ outnpy) 
logger.info('outnpy meta csv '+ outnpy_meta) 

logger.info('trained_model '+str(trained_model))
logger.info('the full graph gpickle '+str(gpickle))
logger.info('the full graph meta '+str(gpickle_meta))
logger.info('number_of_samples1 '+str(number_of_samples1))
logger.info('number_of_samples2 '+str(number_of_samples2))


g = nx.read_gpickle(gpickle)
## graph meta saved before. In subsetXeniumGraphFullPanel.py
graph_meta = pd.read_csv(gpickle_meta)

logger.info("g.number_of_nodes()" + str(g.number_of_nodes()))
logger.info("graph meta shape " + str(graph_meta.shape))

logger.info("graph meta laste index " + str(graph_meta.index[-1]))
assert(list(g.nodes())[-1] == (graph_meta.shape[0]-1))


logger.info("loading saved trained model dir" + str(trained_model))

embedding_model = keras.models.load_model(trained_model,
                                          custom_objects={'AttentionalAggregator':AttentionalAggregator})


embedding_model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=1e-3),
    loss=keras.losses.binary_crossentropy,
    metrics=[keras.metrics.binary_accuracy],
)

logger.info(str(embedding_model.summary()))

# ## Get embeddings

from stellargraph.mapper import GraphSAGENodeGenerator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx   
import pandas as pd
import numpy as np
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
from sklearn import preprocessing, feature_extraction, model_selection
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.metrics import accuracy_score

from stellargraph import globalvar

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
#from umap import UMAP
from stellargraph.mapper import GraphSAGENodeGenerator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

logger.info("filter min comp "+str(min_comp))

def filterGraph(G,graph_meta, min_n):
    node_removed = []
    for component in tqdm(list(nx.connected_components(G))):
        if len(component) < min_n:
            for node in component:
                node_removed.append(node)
                G.remove_node(node)
    if(len(node_removed) > 0):
        graph_meta = graph_meta.drop(node_removed,axis=0)
    logger.info("node removed: " + str(len(node_removed)))
    return G,graph_meta


g,graph_meta = filterGraph(g, graph_meta,min_comp)

logger.info("filtered g  " + str(g.number_of_nodes()))
logger.info('First 5  node label' + str(list(g.nodes())[0:5]))
logger.info('Filtered graph meta shape ' + str(graph_meta.shape))
assert(graph_meta.shape[0] == g.number_of_nodes())

graph_meta.to_csv(outnpy_meta)

if isinstance(g,sg.StellarGraph):
    logger.info("Graph is already a sg graph ")
else:
    g = sg.StellarGraph.from_networkx(g,node_features="feature")
    logger.info("Changed to stellargraph class ")

print(g.info())

logger.info(str(g.info()))


batch_size = 100
num_samples = [number_of_samples1, number_of_samples2]

node_ids = g.nodes()
node_gen = GraphSAGENodeGenerator(g, batch_size, num_samples).flow(node_ids)
node_embeddings = embedding_model.predict(node_gen, workers=4, verbose=1)
logger.info("saving output embedding to "+outnpy)

np.save(outnpy,node_embeddings)

logger.info("saving output embedding to "+outnpy+"done")
