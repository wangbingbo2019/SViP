# SViP
# SViP

**SViP (Spatially Variable Interaction Patterns)** is a computational framework for identifying, clustering, and matching spatially variable gene–gene interaction patterns from spatial transcriptomics data.

SViP integrates interaction filtering, spatial correlation analysis, image-based representation learning, clustering, and spatial-domain matching to systematically characterize spatially variable interaction patterns (SVIPs).

---

## Usage Guide

### 1. Running Environment

- Operating system: Linux / macOS (recommended)
- Python ≥ 3.6
- MATLAB (required only for Step 2)
- Main Python dependencies:
  - numpy
  - pandas
  - scipy
  - scikit-learn
  - tensorflow / keras
  - matplotlib
  - pillow

---

### 2. Project Structure

```text
SViP/
├── step1.py
├── step2.m
├── step3.py
├── step4.sh
├── step5.py
├── step6.py
├── batch_embedding.py
├── model_train.ipynb
├── Model/
│   ├── CAE.h5
│   └── CAE_Encoder.h5
├── data/
│   └── <dataset_name>/
│       ├── count.csv
│       ├── loc.csv
│       └── domain.xlsx
└── SViP_results/
```

### 3. Input Data Format

Each dataset should be placed under:

```
./data/<dataset_name>/
```

#### Required files

| File          | Description                                        |
| ------------- | -------------------------------------------------- |
| `count.csv`   | Spot × gene expression matrix (tab-separated)      |
| `loc.csv`     | Spatial coordinates of spots (x, y; tab-separated) |
| `domain.xlsx` | Spatial domain annotations                         |

------

### 4. Workflow  

```
Step 1  Candidate interaction filtering
Step 2  Spatial correlation computation (MATLAB)
Step 3  Interaction-to-image transformation
Step 4  Image embedding via convolutional autoencoder
Step 5  Interaction clustering and spatial pattern filtering
Step 6  SVIP-to-domain matching
```

------

#### Step 1: Candidate Interaction Filtering

##### Command

```
python step1.py <dataset_name> SpatialDE
```

##### Output files

```
SViP_results/<dataset_name>/step1/step1_out/edges_qv.xlsx: filtered gene–gene interactions
SViP_results/<dataset_name>/step1/step1_out/act_edge_name.xlsx: interaction edge names
SViP_results/<dataset_name>/step1/step1_out/act_gene_count.csv: gene statistics for interactions
```

------

#### Step 2: Spatial Correlation Calculation (MATLAB)

##### Command

1. Open `step2.m`
2. Set dataset name:
3. 3.Run the script

##### Output file

```
SViP_results/<dataset_name>/step2/step2_out/act.csv:gene-gene interactions
```

------

#### Step 3: Interaction Image Construction

##### Command

```
python step3.py <dataset_name>
```

##### Output files

```
SViP_results/<dataset_name>/step3/step3_out/act_new.csv: updated interaction matrix
SViP_results/<dataset_name>/step3/step3_out/pic/grey/: grayscale interaction images
SViP_results/<dataset_name>/step3/step3_out/pic/resize/: resized grayscale images
```

------

#### Step 4: Image Embedding Using CAE

##### Command

```
bash step4.sh batch_size total_size <dataset_name>
```

##### Output files

```
SViP_results/<dataset_name>/step4/step4_out/embed/:embeding vetors
```

------

#### Step 5: Interaction Clustering and Spatial Filtering

##### Command

```
python step5.py <dataset_name>
```

##### Output files

```
SViP_results/<dataset_name>/step5/step5_out/cluster/cluster.csv: cluster centers
SViP_results/<dataset_name>/step5/step5_out/raw_cluster.csv: cluster assignment for each interaction
SViP_results/<dataset_name>/step5/step5_out/svip_filtered.xlsx: spatial variability p-values
SViP_results/<dataset_name>/step5/step5_out/local_domain_score.csv: local silhouette scores
SViP_results/<dataset_name>/step5/step5_out/global_domain_score.csv: global silhouette scores
SViP_results/<dataset_name>/step5/step5_out/random_score/: random control scores
SViP_results/<dataset_name>/step5/step5_out/patterns/: spatial visualizations of interaction patterns
```

------

#### Step 6: SViP-to-Domain Matching

##### Command

```
python step6.py <dataset_name>
```

##### Output file

```
SViP_results/<dataset_name>/step6/step6_out/svip_res.csv:Rows correspond to SVIPs and columns correspond to spatial domains.p-values greater than 0.05 indicate no significant SVIP–domain association.
```

## Example: PDAC-A Dataset

This section demonstrates the complete SViP workflow using the **PDAC-A** spatial transcriptomics dataset.

------

### Data Preparation

```
./data/PDAC-A/
├── count.csv
├── loc.csv
└── domain.xlsx
```

------

### Step-by-step Commands

```
python step1.py PDAC-A SpatialDE
dataset_name = 'PDAC-A';
run step2.m
python step3.py PDAC-A
bash step4.sh batch_size total_size PDAC-A
python step5.py PDAC-A
python step6.py PDAC-A
```

------

### Main result Files

```
SViP_results/PDAC-A/
├── step5/step5_out/
│   ├── cluster/
│   ├── svip_filtered.xlsx
│   ├── local_domain_score.csv
│   ├── global_domain_score.csv
│   └── patterns/
└── step6/step6_out/
    └── svip_res.csv
```

- `svip_res.csv` provides p-values for matching SVIPs to spatial domains
- Spatial visualizations are generated for each matched domain

------

## Notes

- MATLAB is required only for Step 2
- All other steps are implemented in Python
- The workflow is fully reproducible and suitable for publication use
