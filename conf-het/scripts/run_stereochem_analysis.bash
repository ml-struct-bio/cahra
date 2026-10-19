#!/usr/bin/bash

METRICS_LABEL="_label"
LOG_DIR="stereochem_analysis${METRICS_LABEL}"
RESULTS_DIR="metric_calculation${METRICS_LABEL}"

mkdir -p $LOG_DIR

STEREOCHEM_DIR="run_rota_rama_clash${METRICS_LABEL}"

GT_SUBDIR="gt_confs"
ORIG_SUBDIR="orig_6xou"
PRED_SUBDIR="pred_confs"

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

cahra-metrics stereochem_analysis -v \
-o $RESULTS_DIR -s $STEREOCHEM_DIR -gts $GT_SUBDIR -os $ORIG_SUBDIR -ps $PRED_SUBDIR \
-pl $poised_bobcat_1_label $poised_bobcat_2_label $poised_bobcat_3_label $exotic_stingray_1_label $exotic_stingray_2_label $curious_fulmar_1_label $curious_fulmar_2_label $persuasive_penguins_1_label $persuasive_penguins_2_label $persuasive_penguins_3_label $persuasive_penguins_4_label $persuasive_penguins_5_label $kiharalab_1_label $kiharalab_2_label \
> $LOG_DIR/stereochem_analysis${METRICS_LABEL}.log 2>&1 &
