#!/usr/bin/bash

METRICS_LABEL="_label"
LOG_DIR="ensemble_maps${METRICS_LABEL}"

mkdir -p $LOG_DIR

poised_bobcat_1_label="poised-bobcat-1"
poised_bobcat_2_label="poised-bobcat-2"
poised_bobcat_3_label="poised-bobcat-3"
exotic_stingray_1_label="exotic-stingray-1"
exotic_stingray_2_label="exotic-stingray-2"
curious_fulmar_1_label="curious-fulmar-1"
curious_fulmar_2_label="curious-fulmar-2"
persuasive_penguins_1_label="persuasive-penguins-1"
persuasive_penguins_2_label="persuasive-penguins-2"
persuasive_penguins_3_label="persuasive-penguins-3"
persuasive_penguins_4_label="persuasive-penguins-4"
persuasive_penguins_5_label="persuasive-penguins-5"
kiharalab_1_label="kiharalab-1"
kiharalab_2_label="kiharalab-2"
poised_bobcat_1_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${poised_bobcat_1_label}/pdbs"
poised_bobcat_2_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${poised_bobcat_2_label}/pdbs"
poised_bobcat_3_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${poised_bobcat_3_label}/pdbs"
exotic_stingray_1_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${exotic_stingray_1_label}/pdbs"
exotic_stingray_2_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${exotic_stingray_2_label}/pdbs"
curious_fulmar_1_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${curious_fulmar_1_label}/pdbs"
curious_fulmar_2_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${curious_fulmar_2_label}/pdbs"
persuasive_penguins_1_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${persuasive_penguins_1_label}/pdbs"
persuasive_penguins_2_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${persuasive_penguins_2_label}/pdbs"
persuasive_penguins_3_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${persuasive_penguins_3_label}/pdbs"
persuasive_penguins_4_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${persuasive_penguins_4_label}/pdbs"
persuasive_penguins_5_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${persuasive_penguins_5_label}/pdbs"
kiharalab_1_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${kiharalab_1_label}/pdbs"
kiharalab_2_pred_ensemble_dir="metric_calculation${METRICS_LABEL}/pred_confs/${kiharalab_2_label}/pdbs"

gtlabel="GT"
gt_ensemble_dir="metric_calculation${METRICS_LABEL}/gt_confs/pdbs"

cahra-metrics ensemble_maps \
-pe \
$poised_bobcat_1_pred_ensemble_dir $poised_bobcat_2_pred_ensemble_dir $poised_bobcat_3_pred_ensemble_dir \
$exotic_stingray_1_pred_ensemble_dir $exotic_stingray_2_pred_ensemble_dir \
$curious_fulmar_1_pred_ensemble_dir $curious_fulmar_2_pred_ensemble_dir \
$persuasive_penguins_1_pred_ensemble_dir $persuasive_penguins_2_pred_ensemble_dir $persuasive_penguins_3_pred_ensemble_dir $persuasive_penguins_4_pred_ensemble_dir $persuasive_penguins_5_pred_ensemble_dir \
$kiharalab_1_pred_ensemble_dir $kiharalab_2_pred_ensemble_dir \
-pl \
$poised_bobcat_1_label $poised_bobcat_2_label $poised_bobcat_3_label \
$exotic_stingray_1_label $exotic_stingray_2_label \
$curious_fulmar_1_label $curious_fulmar_2_label \
$persuasive_penguins_1_label $persuasive_penguins_2_label $persuasive_penguins_3_label $persuasive_penguins_4_label $persuasive_penguins_5_label \
$kiharalab_1_label $kiharalab_2_label \
-te $gt_ensemble_dir \
-v > "${LOG_DIR}/ensemble_maps.log"

