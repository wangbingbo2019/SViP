#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: $0 batch_size total_size dataset_name"
    exit 1
fi

BATCH_SIZE=$1
TOTAL_SIZE=$2
DATASET_NAME=$3
BATCH_NUM=$((TOTAL_SIZE / BATCH_SIZE))

# =====================
# Base directory (current directory)
# =====================

BASE_DIR=$(pwd)

# =====================
# Input & Output paths
# =====================

FILE_PATH="${BASE_DIR}/SViP_results/${DATASET_NAME}/step3/step3_out/pic/resize/"
RES_PATH="${BASE_DIR}/SViP_results/${DATASET_NAME}/step4/step4_out/embed/"

# create output directory if not exists
mkdir -p "${RES_PATH}"

# =====================
# Batch embedding
# =====================

for ((i=0; i<${BATCH_NUM}; i++))
do
    python batch_embedding.py "${FILE_PATH}" "${RES_PATH}" "${BATCH_SIZE}" "${i}"
    echo "Processed batch ${i}"
done

# Process remaining samples
OTHERS=$((TOTAL_SIZE - BATCH_NUM * BATCH_SIZE))
python batch_embedding.py "${FILE_PATH}" "${RES_PATH}" "${BATCH_SIZE}" "${BATCH_NUM}"
