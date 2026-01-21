import pandas as pd
import encoding
import sys
import os

dataset_name = sys.argv[1]

# =====================
# Base directory
# =====================

BASE_DIR = os.getcwd()

DATA_DIR = os.path.join(BASE_DIR, 'data', dataset_name)
RESULT_DIR = os.path.join(BASE_DIR, 'SViP_results', dataset_name)

# =====================
# Input paths
# =====================

act_csv = os.path.join(
    RESULT_DIR, 'step2', 'step2_out', 'act.csv'
)

act_edge_name_excel = os.path.join(
    RESULT_DIR, 'step1', 'step1_out', 'act_edge_name.xlsx'
)

act_gene_count_csv = os.path.join(
    RESULT_DIR, 'step1', 'step1_out', 'act_gene_count.csv'
)

location_csv = os.path.join(
    DATA_DIR, 'loc.csv'
)

# =====================
# Output paths
# =====================

STEP3_OUT_DIR = os.path.join(
    RESULT_DIR, 'step3', 'step3_out'
)
PIC_DIR = os.path.join(STEP3_OUT_DIR, 'pic')

os.makedirs(STEP3_OUT_DIR, exist_ok=True)
os.makedirs(PIC_DIR, exist_ok=True)

out_act = os.path.join(
    STEP3_OUT_DIR, 'act_new.csv'
)

# =====================
# Run encoding
# =====================

encoding.prep(
    act_csv,
    act_edge_name_excel,
    act_gene_count_csv,
    out_act
)

spot_size = 50
encoding.picture_generate(
    out_act,
    location_csv,
    PIC_DIR,
    spot_size
)
