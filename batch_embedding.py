import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import kendalltau, spearmanr
import scipy.stats
from sklearn import preprocessing
import matplotlib.pyplot as plt
from pandas import DataFrame
from itertools import chain
import cv2
from PIL import Image
import os
from sklearn.model_selection import train_test_split
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import random
import h5py
import gc
import sys
import tensorflow.compat.v1 as tf
tf.disable_eager_execution()
print(tf.__version__)
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
from keras.models import load_model
from keras import backend as K
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# =====================
# Base directory
# =====================

BASE_DIR = os.getcwd()

# =====================
# Load models (relative path)
# =====================

model_1 = load_model(
    os.path.join(BASE_DIR, 'Model', 'CAE_Encoder.h5')
)
model_2 = load_model(
    os.path.join(BASE_DIR, 'Model', 'CAE.h5')
)

# =====================
# Input arguments
# =====================

filepath = sys.argv[1]      # image directory
res_path = sys.argv[2]      # embedding output directory
batch_size = int(sys.argv[3])
batch_id = int(sys.argv[4])

filelist = os.listdir(filepath)
list_ = []
edges_name = []

# =====================
# Read image names
# =====================

for filename in filelist:
    filename1 = os.path.splitext(filename)[1]
    filename0 = os.path.splitext(filename)[0]
    list_.append(filename0 + filename1)
    edges_name.append(filename0)

# =====================
# Load images
# =====================

x = []
for i in range(batch_id * batch_size, (batch_id + 1) * batch_size):
    if i >= len(list_):
        break
    else:
        path = os.path.join(filepath, list_[i])
        img = Image.open(path)
        image_matrix = np.array(img)
        x.append(image_matrix)

data = np.array(x)
data = np.expand_dims(data, axis=-1)

sess = tf.Session()
data = tf.cast(data, tf.float32) / 255
data = data.eval(session=sess)

# =====================
# Encode
# =====================

pre_test = model_1.predict(data)
X_encoded_reshape = pre_test.reshape(
    pre_test.shape[0],
    pre_test.shape[1] * pre_test.shape[2] * pre_test.shape[3]
)

df = pd.DataFrame(X_encoded_reshape)
df_ = df.T

if (batch_id + 1) * batch_size <= len(edges_name):
    df_.columns = edges_name[
        batch_id * batch_size:(batch_id + 1) * batch_size
    ]
else:
    df_.columns = edges_name[batch_id * batch_size:]

# =====================
# Save embedding
# =====================

os.makedirs(res_path, exist_ok=True)
df_.to_csv(
    os.path.join(res_path, f'embed_{batch_id}.csv')
)
