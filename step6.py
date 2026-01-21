import pandas as pd
import sys
import numpy as np
import os
from sklearn.metrics import silhouette_samples
import seaborn as sns
import matplotlib.pyplot as plt

dataset_name = sys.argv[1]

random_times = 100
local_option = 'percentage'

# =====================
# Path settings (relative to current directory)
# =====================

BASE_DIR = os.getcwd()

RESULT_DIR = os.path.join(BASE_DIR, 'SViP_results', dataset_name)
DATA_DIR = os.path.join(BASE_DIR, 'data', dataset_name)

# Step 5 outputs
global_score_path = os.path.join(
    RESULT_DIR, 'step5', 'step5_out', 'global_domain_score.csv'
)
local_score_path = os.path.join(
    RESULT_DIR, 'step5', 'step5_out', 'local_domain_score.csv'
)

# Step 3 output
act_path = os.path.join(
    RESULT_DIR, 'step3', 'step3_out', 'act_new.csv'
)

# Data
loc_path = os.path.join(DATA_DIR, 'loc.csv')
domain_path = os.path.join(DATA_DIR, 'domain.xlsx')

# Step 6 output directory
STEP6_OUT_DIR = os.path.join(
    RESULT_DIR, 'step6', 'step6_out'
)
os.makedirs(STEP6_OUT_DIR, exist_ok=True)

# =====================
# Load data
# =====================

global_score = pd.read_csv(global_score_path, index_col=0)
local_score = pd.read_csv(local_score_path, index_col=0)
df_act = pd.read_csv(act_path, sep='\t', index_col=0)
loc = pd.read_csv(loc_path, sep='\t', index_col=0)
domains = pd.read_excel(domain_path, header=0, index_col=0)

domain_labels = np.array(domains.iloc[:, 2].tolist())
domain_label_set = list(set(domain_labels))

df_pval = pd.DataFrame(index=global_score.index, columns=global_score.columns)

# =====================
# Main loop
# =====================

for i in range(len(global_score)):

    global_control_path = os.path.join(
        RESULT_DIR, 'step5', 'step5_out', 'random_score',
        f'global_control_score_{random_times}_domain_{global_score.index[i]}.csv'
    )
    local_control_path = os.path.join(
        RESULT_DIR, 'step5', 'step5_out', 'random_score',
        f'local_control_score_{random_times}_domain_{global_score.index[i]}.csv'
    )

    global_control_score = pd.read_csv(global_control_path, index_col=0).T
    local_control_score = pd.read_csv(local_control_path, index_col=0).T

    in_domain_labels = [x == global_score.index[i] for x in domain_labels]
    out_domain_labels = [not x for x in in_domain_labels]
    domain_n = sum(in_domain_labels)

    print(global_score.index[i])

    max_score = 0
    max_ind = 0
    max_interest_spot = 0

    for j in range(len(global_score.columns)):

        data = df_act.loc[global_score.columns[j], :].values
        mid = (np.max(data) + np.min(data)) * 3 / 4
        interest_spot = [1 if x > mid else 0 for x in data]

        if sum(interest_spot) <= domain_n * 0.2 or sum(interest_spot) >= len(loc) * 0.99:
            score = -1
            control_score = [1] * random_times

        elif sum(interest_spot) <= min(domain_n * 1.5, len(loc) * 0.2):

            if local_option == 'percentage':
                score = sum(
                    1 for k in range(len(interest_spot))
                    if interest_spot[k] == 1 and in_domain_labels[k]
                )
                score = score / sum(interest_spot)
                score = (score - 0.5) * 2
                flag = 'percentage'

            elif local_option == 'local_silhouette_score':
                score = local_score.iloc[i, j]
                flag = 'local_silhouette_score'

            control_score = sorted(local_control_score.iloc[j, :].values, reverse=True)

        elif sum(interest_spot) > len(loc) * 0.2:
            score = global_score.iloc[i, j]
            flag = 'global_silhouette_score'
            control_score = sorted(global_control_score.iloc[j, :].values, reverse=True)

        else:
            score = -1
            control_score = [1] * random_times

        in_domain_exp = np.mean(df_act.loc[global_score.columns[j], df_act.columns[in_domain_labels]])
        out_domain_exp = np.mean(df_act.loc[global_score.columns[j], df_act.columns[out_domain_labels]])

        if in_domain_exp > out_domain_exp:
            k = 0
            while k < random_times and score < control_score[k]:
                k += 1
            df_pval.iloc[i, j] = float(k) / random_times
        else:
            df_pval.iloc[i, j] = 1.0

        if score > max_score and in_domain_exp > out_domain_exp:
            max_score = score
            max_ind = j
            max_interest_spot = sum(interest_spot)
            max_flag = flag
            tmp_in = in_domain_exp
            tmp_out = out_domain_exp

    # =====================
    # Visualization
    # =====================

    if df_pval.iloc[i, max_ind] < 0.05:

        plt_color = df_act.loc[global_score.columns[max_ind], :].values

        plt.figure(dpi=300, figsize=(6, 4))
        plt.gca().invert_yaxis()
        plt.scatter(
            loc['x'].values,
            loc['y'].values,
            marker='o',
            c=plt_color,
            cmap='bwr',
            vmin=min(plt_color) - 0.2,
            s=2
        )

        fig_path = os.path.join(
            STEP6_OUT_DIR,
            f'{global_score.index[i]}_{global_score.columns[max_ind]}.jpg'
        )
        plt.savefig(fig_path)
        plt.close()

# =====================
# Save p-values
# =====================

df_pval.to_csv(
    os.path.join(STEP6_OUT_DIR, 'svip_res.csv')
)
