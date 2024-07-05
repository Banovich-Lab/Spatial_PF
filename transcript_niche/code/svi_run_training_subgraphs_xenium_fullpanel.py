# code/run_training_subgraphs_xenium_fullpanel.py
#!/usr/bin/env python
# coding: utf-8

# In[3]:



import tensorflow
import tensorrt
import numpy as np
import random as rn

from tensorflow.keras import backend as K

# The below tf.set_random_seed() will make random number generation
# in the TensorFlow backend have a well-defined initial state.
# For further details, see:
# https://www.tensorflow.org/api_docs/python/tf/set_random_seed

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


# The below is necessary for starting Numpy generated random numbers
# in a well-defined initial state.
np.random.seed(42)
rn.seed(12345)

from tensorflow.random import set_seed
set_seed(42)

import sys
import logging 
import builtins
import shutil
## file.py num_walks, num_root_nodes, number_of_samples1, number_of_samples2,buildGraph: (KNN/radius), output_dir
#now we will Create and configure logger 

logging.basicConfig(filename=sys.argv[10], 
                        format='%(asctime)s %(message)s', 
                        filemode='w',level=logging.INFO) 

#Let us Create an object 
logger=logging.getLogger() 

logger.info('GPU name: '+ str(tensorflow.config.experimental.list_physical_devices("GPU")))

# ## Read in the graph for each region

logger.info('Length of list:'+ str(len(sys.argv)))

## length of walks for sampling positive pairs

number_of_walks = int(sys.argv[1])
num_length =  int(sys.argv[2])
number_of_samples1 = int(sys.argv[3])

number_of_samples2 = int(sys.argv[4])
gpickle = sys.argv[5]
number_epoch = int(sys.argv[6])
trained_model = sys.argv[7]
trained_emModel  = sys.argv[8]
sampled_roots_toIntcsv = sys.argv[9]
logger.info(str(sys.argv))

# ### Read in the subgraphs 
print('GPU name: ', tensorflow.config.experimental.list_physical_devices("GPU"))

logger.info("reading gpickle " +gpickle )
joined_g_10k_roots  = nx.read_gpickle(gpickle)
logger.info(str(joined_g_10k_roots.info()))

root_nodes_as_int = pd.read_csv(sampled_roots_toIntcsv)

logger.info(str(root_nodes_as_int["0"][0:4]))

unsupervised_samples = UnsupervisedSampler(joined_g_10k_roots, nodes=root_nodes_as_int["0"], 
                                           length=num_length, number_of_walks=number_of_walks, seed=42)

batch_size = 100
epoch = number_epoch
num_samples = [number_of_samples1, number_of_samples2]
layer_sizes = [50, 50]



generator = GraphSAGELinkGenerator(joined_g_10k_roots, batch_size, num_samples)
train_gen = generator.flow(unsupervised_samples)

# In[39]:
graphsage = GraphSAGE(
    layer_sizes=layer_sizes, generator=generator, bias=True, 
    aggregator = AttentionalAggregator, dropout=0.0, normalize="l2"
)
# In[40]:
tf.keras.initializers.GlorotUniform(
    seed=42
)
# In[41]:
# Build the model and expose input and output sockets of graphsage, for node pair inputs:
x_inp, x_out = graphsage.in_out_tensors()


# In[42]:


prediction = link_classification(
    output_dim=1, output_act="sigmoid", edge_embedding_method="ip"
)(x_out)


# In[43]:


model = keras.Model(inputs=x_inp, outputs=prediction)

model.compile(
    optimizer=keras.optimizers.Adam(learning_rate=1e-3),
    loss=keras.losses.binary_crossentropy,
    metrics=[keras.metrics.binary_accuracy],
)


# In[ ]:


history = model.fit(
    train_gen,
    epochs=number_epoch,
    verbose=1,
    use_multiprocessing=True,
    workers=4,
    shuffle=True
)


# In[50]:

logger.info(str(history.history))


# In[ ]:


# Save the weights
model.save(trained_model)

## change the input/output size before saving the model
x_inp_src = x_inp[0::2]
x_out_src = x_out[0]
embedding_model = keras.Model(inputs=x_inp_src, outputs=x_out_src)

embedding_model.save(trained_emModel)


logger.info(str(embedding_model.summary()))

