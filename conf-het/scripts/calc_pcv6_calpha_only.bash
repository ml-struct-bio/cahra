#!/usr/bin/bash
METRICS_LABEL="_v6_calpha_only"

LOG_DIR="calc_pcv${METRICS_LABEL}"
OUTPUT_DIR="metric_calculation${METRICS_LABEL}"

TRUTH_ENSEMBLE_DIR="${OUTPUT_DIR}/gt_confs/pdbs"
MAX_RANK=100 # updated from 50 to 100 from v4 and up

mkdir $LOG_DIR

# min-RMSD aligned to GT ensemble PDB locations
lemon1_dir="${OUTPUT_DIR}/pred_confs/lemon1/pdbs"
lemon2_dir="${OUTPUT_DIR}/pred_confs/lemon2/pdbs"
lemon3_dir="${OUTPUT_DIR}/pred_confs/lemon3/pdbs"
lime1_dir="${OUTPUT_DIR}/pred_confs/lime1/pdbs"
lime2_dir="${OUTPUT_DIR}/pred_confs/lime2/pdbs"
orange1_dir="${OUTPUT_DIR}/pred_confs/orange1/pdbs"
orange2_dir="${OUTPUT_DIR}/pred_confs/orange2/pdbs"
kiwi1_dir="${OUTPUT_DIR}/pred_confs/kiwi1/pdbs"
kiwi2_dir="${OUTPUT_DIR}/pred_confs/kiwi2/pdbs"
kiwi3_dir="${OUTPUT_DIR}/pred_confs/kiwi3/pdbs"
kiwi4_dir="${OUTPUT_DIR}/pred_confs/kiwi4/pdbs"
kiwi5_dir="${OUTPUT_DIR}/pred_confs/kiwi5/pdbs"
banana1_dir="${OUTPUT_DIR}/pred_confs/banana1/pdbs"
banana2_dir="${OUTPUT_DIR}/pred_confs/banana2/pdbs"

lemon1_label="lemon1"
lemon2_label="lemon2"
lemon3_label="lemon3"
lime1_label="lime1"
lime2_label="lime2"
orange1_label="orange1"
orange2_label="orange2"
kiwi1_label="kiwi1"
kiwi2_label="kiwi2"
kiwi3_label="kiwi3"
kiwi4_label="kiwi4"
kiwi5_label="kiwi5"
banana1_label="banana1"
banana2_label="banana2"

cahra-metrics calc_pcv -v -r $MAX_RANK -te $TRUTH_ENSEMBLE_DIR \
-pe $lemon1_dir $lemon2_dir $lemon3_dir $lime1_dir $lime2_dir $orange1_dir $orange2_dir $kiwi1_dir $kiwi2_dir $kiwi3_dir $kiwi4_dir $kiwi5_dir $banana1_dir $banana2_dir \
-pl $lemon1_label $lemon2_label $lemon3_label $lime1_label $lime2_label $orange1_label $orange2_label $kiwi1_label $kiwi2_label $kiwi3_label $kiwi4_label $kiwi5_label $banana1_label $banana2_label \
> $LOG_DIR/calc_pcv${METRICS_LABEL}.log 2>&1 &
