import pandas as pd
import SpatialDE
import NaiveDE
import clustering
import random_test_rv2
import si_score_rv1
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples

dataset_name = sys.argv[1]

# =====================
# Path settings
# =====================

BASE_DIR = os.getcwd()

DATA_DIR = os.path.join(BASE_DIR, 'data', dataset_name)
RESULT_DIR = os.path.join(BASE_DIR, 'SViP_results', dataset_name)

# Inputs
batch_dir = os.path.join(
    RESULT_DIR, 'step4', 'step4_out', 'embed'
)

df_count = pd.read_csv(
    os.path.join(DATA_DIR, 'count.csv'),
    sep='\t',
    index_col=0
)

loc = pd.read_csv(
    os.path.join(DATA_DIR, 'loc.csv'),
    sep='\t',
    index_col=0
)

df_act = pd.read_csv(
    os.path.join(RESULT_DIR, 'step3', 'step3_out', 'act_new.csv'),
    sep='\t',
    index_col=0
)

act_file = os.path.join(
    RESULT_DIR, 'step3', 'step3_out', 'act_new.csv'
)

domain_file = os.path.join(DATA_DIR, 'domain.xlsx')

# Outputs
STEP5_OUT_DIR = os.path.join(RESULT_DIR, 'step5', 'step5_out')
CLUSTER_DIR = os.path.join(STEP5_OUT_DIR, 'cluster')
PATTERN_DIR = os.path.join(STEP5_OUT_DIR, 'patterns')
RANDOM_SCORE_DIR = os.path.join(STEP5_OUT_DIR, 'random_score')

os.makedirs(CLUSTER_DIR, exist_ok=True)
os.makedirs(PATTERN_DIR, exist_ok=True)
os.makedirs(RANDOM_SCORE_DIR, exist_ok=True)

# =====================
# Clustering
# =====================

clustering.clustering(batch_dir, CLUSTER_DIR, opt='not_auto')

raw_cluster_file = os.path.join(CLUSTER_DIR, 'raw_cluster.csv')
cluster_file = os.path.join(CLUSTER_DIR, 'cluster.csv')
svip_filtered_file = os.path.join(STEP5_OUT_DIR, 'svip_filtered.xlsx')

# =====================
# Sample information
# =====================

sample_info = pd.DataFrame(columns=['', 'x', 'y', 'total_counts'])
sample_info[''] = df_count.index.tolist()
sample_info['x'] = loc['x'].values
sample_info['y'] = loc['y'].values
sample_info['total_counts'] = df_act.sum(0).tolist()
sample_info = sample_info.set_index('')

# =====================
# SpatialDE filtering
# =====================

norm_expr = NaiveDE.stabilize(df_act).T
X = sample_info[['x', 'y']]

results = SpatialDE.run(X, df_act.T)
res = results.sort_values('qval')[['g', 'pval', 'qval']]
res.to_excel(svip_filtered_file, index=False)

res = pd.read_excel(svip_filtered_file, index_col=0)
res = res[res['qval'] < 0.05]
svip = res.index

# =====================
# Cluster selection
# =====================

df_all_svip_clusters = pd.read_csv(raw_cluster_file, index_col=0)
df_clusters = pd.read_csv(cluster_file, index_col=0)
df_clusters.drop_duplicates(inplace=True)

edges_top = df_clusters.index
svip_clustered = [i for i in svip if i in edges_top]

print('Number of clustered SVIPs:', len(svip_clustered))

# =====================
# Pattern visualization
# =====================

for edge in svip_clustered:

    plt.figure(dpi=300, figsize=(6, 4))
    plt.gca().invert_yaxis()
    plt.scatter(
        loc['x'].values,
        loc['y'].values,
        marker='o',
        c=df_act.loc[edge, :],
        s=50
    )
    plt.axis('equal')

    value = df_clusters.loc[edge, 'cluster_label']
    if isinstance(value, pd.Series):
        value = value.iloc[0]
    pattern_id = int(value)

    plt.title(
        'Pattern {} - {} genes'.format(
            pattern_id,
            df_all_svip_clusters['cluster_label'].eq(pattern_id).sum()
        )
    )

    plt.savefig(
        os.path.join(PATTERN_DIR, f'Pattern_{pattern_id}.jpg')
    )
    plt.close()

# =====================
# Silhouette evaluation & random control
# =====================

si_score_rv1.run(
    cluster_file,
    svip_filtered_file,
    act_file,
    domain_file,
    STEP5_OUT_DIR
)

random_test_rv2.random_cal(
    cluster_file,
    svip_filtered_file,
    act_file,
    domain_file,
    RANDOM_SCORE_DIR
)
