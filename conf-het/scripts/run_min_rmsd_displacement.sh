#!/usr/bin/bash

METRICS_LABEL="_v6"
LOG_DIR="stereochem_analysis${METRICS_LABEL}"
RESULTS_DIR="metric_calculation${METRICS_LABEL}"
JSON="all_metrics.json"

mkdir -p $LOG_DIR

lemon1_label="lemon1"
lemon1_json="${RESULTS_DIR}/pred_confs/${lemon1_label}/results/${JSON}"
lemon2_label="lemon2"
lemon2_json="${RESULTS_DIR}/pred_confs/${lemon2_label}/results/${JSON}"
lemon3_label="lemon3"
lemon3_json="${RESULTS_DIR}/pred_confs/${lemon3_label}/results/${JSON}"
lime1_label="lime1"
lime1_json="${RESULTS_DIR}/pred_confs/${lime1_label}/results/${JSON}"
lime2_label="lime2"
lime2_json="${RESULTS_DIR}/pred_confs/${lime2_label}/results/${JSON}"
orange1_label="orange1"
orange1_json="${RESULTS_DIR}/pred_confs/${orange1_label}/results/${JSON}"
orange2_label="orange2"
orange2_json="${RESULTS_DIR}/pred_confs/${orange2_label}/results/${JSON}"
kiwi1_label="kiwi1"
kiwi1_json="${RESULTS_DIR}/pred_confs/${kiwi1_label}/results/${JSON}"
kiwi2_label="kiwi2"
kiwi2_json="${RESULTS_DIR}/pred_confs/${kiwi2_label}/results/${JSON}"
kiwi3_label="kiwi3"
kiwi3_json="${RESULTS_DIR}/pred_confs/${kiwi3_label}/results/${JSON}"
kiwi4_label="kiwi4"
kiwi4_json="${RESULTS_DIR}/pred_confs/${kiwi4_label}/results/${JSON}"
kiwi5_label="kiwi5"
kiwi5_json="${RESULTS_DIR}/pred_confs/${kiwi5_label}/results/${JSON}"
banana1_label="banana1"
banana1_json="${RESULTS_DIR}/pred_confs/${banana1_label}/results/${JSON}"
banana2_label="banana2"
banana2_json="${RESULTS_DIR}/pred_confs/${banana2_label}/results/${JSON}"

cahra-metrics min_rmsd_displacement -pe $lemon1_json \
-pl $lemon1_label \
-v > $LOG_DIR/min_rmsd_displacement${METRICS_LABEL}.log