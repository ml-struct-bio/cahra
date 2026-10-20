#!/usr/bin/bash

METRICS_LABEL="_label"
LOG_DIR="min_rmsd_displacement${METRICS_LABEL}"
RESULTS_DIR="metric_calculation${METRICS_LABEL}"
JSON="all_metrics.json"

mkdir -p $LOG_DIR

poised_bobcat_1_label="poised-bobcat-1"
poised_bobcat_1_json="${RESULTS_DIR}/pred_confs/${poised_bobcat_1_label}/results/${JSON}"
poised_bobcat_2_label="poised-bobcat-2"
poised_bobcat_2_json="${RESULTS_DIR}/pred_confs/${poised_bobcat_2_label}/results/${JSON}"
poised_bobcat_3_label="poised-bobcat-3"
poised_bobcat_3_json="${RESULTS_DIR}/pred_confs/${poised_bobcat_3_label}/results/${JSON}"
exotic_stingray_1_label="exotic-stingray-1"
exotic_stingray_1_json="${RESULTS_DIR}/pred_confs/${exotic_stingray_1_label}/results/${JSON}"
exotic_stingray_2_label="exotic-stingray-2"
exotic_stingray_2_json="${RESULTS_DIR}/pred_confs/${exotic_stingray_2_label}/results/${JSON}"
curious_fulmar_1_label="curious-fulmar-1"
curious_fulmar_1_json="${RESULTS_DIR}/pred_confs/${curious_fulmar_1_label}/results/${JSON}"
curious_fulmar_2_label="curious-fulmar-2"
curious_fulmar_2_json="${RESULTS_DIR}/pred_confs/${curious_fulmar_2_label}/results/${JSON}"
persuasive_penguins_1_label="persuasive-penguins-1"
persuasive_penguins_1_json="${RESULTS_DIR}/pred_confs/${persuasive_penguins_1_label}/results/${JSON}"
persuasive_penguins_2_label="persuasive-penguins-2"
persuasive_penguins_2_json="${RESULTS_DIR}/pred_confs/${persuasive_penguins_2_label}/results/${JSON}"
persuasive_penguins_3_label="persuasive-penguins-3"
persuasive_penguins_3_json="${RESULTS_DIR}/pred_confs/${persuasive_penguins_3_label}/results/${JSON}"
persuasive_penguins_4_label="persuasive-penguins-4"
persuasive_penguins_4_json="${RESULTS_DIR}/pred_confs/${persuasive_penguins_4_label}/results/${JSON}"
persuasive_penguins_5_label="persuasive-penguins-5"
persuasive_penguins_5_json="${RESULTS_DIR}/pred_confs/${persuasive_penguins_5_label}/results/${JSON}"
kiharalab_1_label="kiharalab-1"
kiharalab_1_json="${RESULTS_DIR}/pred_confs/${kiharalab_1_label}/results/${JSON}"
kiharalab_2_label="kiharalab-2"
kiharalab_2_json="${RESULTS_DIR}/pred_confs/${kiharalab_2_label}/results/${JSON}"

cahra-metrics min_rmsd_displacement -pe \
$poised_bobcat_1_json $poised_bobcat_2_json $poised_bobcat_3_json \
$exotic_stingray_1_json $exotic_stingray_2_json \
$curious_fulmar_1_json $curious_fulmar_2_json \
$persuasive_penguins_1_json $persuasive_penguins_2_json $persuasive_penguins_3_json $persuasive_penguins_4_json $persuasive_penguins_5_json $kiharalab_1_json $kiharalab_2_json \
-pl \
$poised_bobcat_1_label $poised_bobcat_2_label $poised_bobcat_3_label \
$exotic_stingray_1_label $exotic_stingray_2_label \
$curious_fulmar_1_label $curious_fulmar_2_label \
$persuasive_penguins_1_label $persuasive_penguins_2_label $persuasive_penguins_3_label $persuasive_penguins_4_label $persuasive_penguins_5_label \
$kiharalab_1_label $kiharalab_2_label \
-v > $LOG_DIR/min_rmsd_displacement${METRICS_LABEL}.log