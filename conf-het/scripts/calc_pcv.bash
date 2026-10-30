#!/usr/bin/bash
METRICS_LABEL="_label"

LOG_DIR="calc_pcv${METRICS_LABEL}"
OUTPUT_DIR="metric_calculation${METRICS_LABEL}"

TRUTH_ENSEMBLE_DIR="${OUTPUT_DIR}/gt_confs/pdbs"
MAX_RANK=100

mkdir $LOG_DIR

# min-RMSD aligned to GT ensemble PDB locations
poised_bobcat_1_dir="${OUTPUT_DIR}/pred_confs/poised-bobcat-1/pdbs"
poised_bobcat_2_dir="${OUTPUT_DIR}/pred_confs/poised-bobcat-2/pdbs"
poised_bobcat_3_dir="${OUTPUT_DIR}/pred_confs/poised-bobcat-3/pdbs"
exotic_stingray_1_dir="${OUTPUT_DIR}/pred_confs/exotic-stingray-1/pdbs"
exotic_stingray_2_dir="${OUTPUT_DIR}/pred_confs/exotic-stingray-2/pdbs"
curious_fulmar_1_dir="${OUTPUT_DIR}/pred_confs/curious-fulmar-1/pdbs"
curious_fulmar_2_dir="${OUTPUT_DIR}/pred_confs/curious-fulmar-2/pdbs"
persuasive_penguins_1_dir="${OUTPUT_DIR}/pred_confs/persuasive-penguins-1/pdbs"
persuasive_penguins_2_dir="${OUTPUT_DIR}/pred_confs/persuasive-penguins-2/pdbs"
persuasive_penguins_3_dir="${OUTPUT_DIR}/pred_confs/persuasive-penguins-3/pdbs"
persuasive_penguins_4_dir="${OUTPUT_DIR}/pred_confs/persuasive-penguins-4/pdbs"
persuasive_penguins_5_dir="${OUTPUT_DIR}/pred_confs/persuasive-penguins-5/pdbs"
kiharalab_1_dir="${OUTPUT_DIR}/pred_confs/kiharalab-1/pdbs"
kiharalab_2_dir="${OUTPUT_DIR}/pred_confs/kiharalab-2/pdbs"

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

cahra-metrics calc_pcv -v -r $MAX_RANK -te $TRUTH_ENSEMBLE_DIR \
-pe $poised_bobcat_1_dir $poised_bobcat_2_dir $poised_bobcat_3_dir $exotic_stingray_1_dir $exotic_stingray_2_dir $curious_fulmar_1_dir $curious_fulmar_2_dir $persuasive_penguins_1_dir $persuasive_penguins_2_dir $persuasive_penguins_3_dir $persuasive_penguins_4_dir $persuasive_penguins_5_dir $kiharalab_1_dir $kiharalab_2_dir \
-pl $poised_bobcat_1_label $poised_bobcat_2_label $poised_bobcat_3_label $exotic_stingray_1_label $exotic_stingray_2_label $curious_fulmar_1_label $curious_fulmar_2_label $persuasive_penguins_1_label $persuasive_penguins_2_label $persuasive_penguins_3_label $persuasive_penguins_4_label $persuasive_penguins_5_label $kiharalab_1_label $kiharalab_2_label \
> $LOG_DIR/calc_pcv${METRICS_LABEL}.log 2>&1 &
