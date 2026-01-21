import pandas as pd
import find_svg
import find_edges
import sys
import os

dataset_name = sys.argv[1]
use_svg = sys.argv[2]

# =====================
# Base directory
# =====================

BASE_DIR = os.getcwd()

DATA_DIR = os.path.join(BASE_DIR, 'data', dataset_name)
RESULT_DIR = os.path.join(BASE_DIR, 'SViP_results', dataset_name)

STEP1_DIR = os.path.join(RESULT_DIR, 'step1')
STEP1_OUT_DIR = os.path.join(STEP1_DIR, 'step1_out')

os.makedirs(STEP1_DIR, exist_ok=True)
os.makedirs(STEP1_OUT_DIR, exist_ok=True)

# =====================
# Input files
# =====================

df_count = os.path.join(DATA_DIR, 'count.csv')
loc = os.path.join(DATA_DIR, 'loc.csv')

# =====================
# SVG result
# =====================

svg_excel = os.path.join(
    STEP1_DIR,
    f'{use_svg}_res.xlsx'
)

find_svg.find_svg(
    df_count,
    loc,
    svg_excel,
    method=use_svg
)

# =====================
# Interaction filtering
# =====================

species = 'human'

edges_qv_out = os.path.join(
    STEP1_OUT_DIR,
    'edges_qv.xlsx'
)

act_edge_name_out_excel = os.path.join(
    STEP1_OUT_DIR,
    'act_edge_name.xlsx'
)

act_gene_count_out_csv = os.path.join(
    STEP1_OUT_DIR,
    'act_gene_count.csv'
)

find_edges.find_edges(
    svg_excel,
    df_count,
    species,
    edges_qv_out,
    act_edge_name_out_excel,
    act_gene_count_out_csv
)
